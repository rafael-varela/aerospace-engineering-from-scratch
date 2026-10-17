"""P03.13 Finite Element Basics: reference solution.

A one-dimensional bar finite-element solver written from scratch: linear (two-node) elements, element
stiffness and consistent-mass matrices, Gauss-quadrature load vectors for a distributed body force,
assembly, Dirichlet (fixed-end) constraints by partitioning, element stresses, reactions, strain energy
and the axial natural frequencies. The engineering case is the spar of a helicopter rotor blade under its
own centrifugal load plus a tip mass, for which the exact displacement and stress are closed-form; a
tapered strut under a tip load is the second closed-form case. ``scipy.integrate.quad``,
``scipy.integrate.solve_bvp``, ``scipy.linalg.solve_banded`` and ``scipy.linalg.eigh`` are the oracles.

Governing equation (weak form derived in the lesson):  -(E A u')' = b(x) on 0 < x < L, u(0) = 0,
E A u'(L) = P (tip load).  b is force per unit length [N/m], u displacement [m], sigma = E u' [Pa].

Conventions: nodes x_0 = 0 < x_1 < ... < x_n = L (n elements, n + 1 nodes); element e spans
[x_e, x_{e+1}] with length h_e; the parent coordinate xi runs from -1 to 1; K is the (n+1) x (n+1) dense
global stiffness matrix; dof j is the displacement of node j.

Run:  python solution.py          (prints the verification tables, writes bar_fem.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-03/lesson-13-finite-element-basics
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np


# ---------------------------------------------------------------- the engineering case

# GIVEN
def rotor_blade() -> dict[str, float]:
    """Return the rotor-blade spar case: an aluminium spar modelled as a bar from the hub centre to the tip.

    E [Pa], rho [kg/m^3] for an aluminium alloy; L [m] the rotor radius; A [m^2] the spar cross-section;
    omega [rad/s] the rotor speed (400 rpm); m_tip [kg] a tip mass (balance weight and tip cap). The
    centrifugal body force per unit length is b(x) = rho A omega^2 x and the tip load is P = m_tip omega^2 L.
    The spar is treated as a uniform bar from x = 0 (hub centre, fixed) to x = L (free tip). Illustrative of
    a light helicopter, not a specific aircraft.
    """
    return {"E": 70.0e9, "rho": 2700.0, "L": 5.0, "A": 6.0e-4, "omega": 400.0 * 2.0 * np.pi / 60.0, "m_tip": 2.0}


# ---------------------------------------------------------------- exact solutions

def exact_rotating_bar(x: np.ndarray, E: float, rho: float, omega: float, L: float, A: float,
                       P: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """Return (u, sigma) for a uniform bar fixed at x = 0 spinning at omega about x = 0, with tip load P.

    E A u'' = -rho A omega^2 x with u(0) = 0, E A u'(L) = P gives
    u = (rho omega^2 / (2E)) (L^2 x - x^3 / 3) + P x / (E A),  sigma = rho omega^2 (L^2 - x^2) / 2 + P / A.
    """
    x = np.asarray(x, dtype=float)
    u = rho * omega ** 2 / (2.0 * E) * (L ** 2 * x - x ** 3 / 3.0) + P * x / (E * A)
    sigma = 0.5 * rho * omega ** 2 * (L ** 2 - x ** 2) + P / A
    return u, sigma


def exact_tapered_bar(x: np.ndarray, E: float, A0: float, c: float, L: float, P: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (u, sigma) for a bar with A(x) = A0 (1 - c x / L), fixed at 0, loaded by P at the tip, no body force.

    E A(x) u' = P is constant, so u = (P L / (E A0 c)) ln(1 / (1 - c x / L)) and sigma = P / A(x).
    Requires 0 < c < 1.
    """
    if not 0.0 < c < 1.0:
        raise ValueError("taper ratio c must be in (0, 1)")
    x = np.asarray(x, dtype=float)
    u = P * L / (E * A0 * c) * np.log(1.0 / (1.0 - c * x / L))
    sigma = P / (A0 * (1.0 - c * x / L))
    return u, sigma


# ---------------------------------------------------------------- element level

def shape_functions(xi: np.ndarray) -> np.ndarray:
    """Return the two linear shape functions [N1, N2] = [(1 - xi)/2, (1 + xi)/2] at parent coordinates xi.

    Shape: (len(xi), 2). They sum to one (a rigid translation is represented exactly) and N1 = 1 at
    xi = -1, N2 = 1 at xi = +1 (the nodal values are the coefficients).
    """
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    return np.column_stack([0.5 * (1.0 - xi), 0.5 * (1.0 + xi)])


def element_stiffness(E: float, A: float, h: float) -> np.ndarray:
    """Return k_e = (E A / h) [[1, -1], [-1, 1]] for a two-node element of length h with constant E A.

    From the weak form: k_ij = int E A N_i' N_j' dx with N' = -/+ 1/h. For a linearly varying A use its
    mid-element (mean) value, which integrates the product E A N_i' N_j' exactly. Raise ``ValueError``
    for h <= 0.
    """
    if h <= 0.0:
        raise ValueError("element length must be positive")
    return E * A / h * np.array([[1.0, -1.0], [-1.0, 1.0]])


def element_mass(rho: float, A: float, h: float) -> np.ndarray:
    """Return the consistent mass matrix m_e = (rho A h / 6) [[2, 1], [1, 2]] (int rho A N_i N_j dx)."""
    if h <= 0.0:
        raise ValueError("element length must be positive")
    return rho * A * h / 6.0 * np.array([[2.0, 1.0], [1.0, 2.0]])


def element_load(b: Callable[[np.ndarray], np.ndarray], x1: float, x2: float, n_gauss: int = 3) -> np.ndarray:
    """Return the consistent load vector f_e = int N^T b(x) dx over [x1, x2] by n_gauss-point Gauss-Legendre.

    Map xi in [-1, 1] to x = (x1 + x2)/2 + xi (x2 - x1)/2, dx = (h/2) dxi. Three points integrate b N
    exactly when b is a polynomial of degree 4 or less. The centrifugal load b = rho A omega^2 x is linear,
    so even two points are exact; the quadrature is kept general for Problem 3.
    """
    h = x2 - x1
    if h <= 0.0:
        raise ValueError("element length must be positive")
    xi, w = np.polynomial.legendre.leggauss(n_gauss)
    x = 0.5 * (x1 + x2) + 0.5 * h * xi
    N = shape_functions(xi)
    return 0.5 * h * (N.T @ (w * b(x)))


# ---------------------------------------------------------------- assembly and solution

def assemble(nodes: np.ndarray, E: float, A: Callable[[np.ndarray], np.ndarray] | float,
             b: Callable[[np.ndarray], np.ndarray] | None = None, rho: float | None = None,
             n_gauss: int = 3) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    """Assemble the global stiffness K, load vector f and (if rho is given) consistent mass M.

    ``A`` is a constant or a function of x (evaluated at each element's mid-point); ``b`` is the body force
    per unit length (None for no body force). Each element's 2 x 2 matrices are added into rows/columns
    (e, e+1) of the global arrays: that scatter-add IS the finite element method's assembly.
    Raise ``ValueError`` unless the nodes are strictly increasing with at least two of them.
    """
    nodes = np.asarray(nodes, dtype=float)
    if nodes.size < 2 or np.any(np.diff(nodes) <= 0.0):
        raise ValueError("nodes must be strictly increasing with at least two entries")
    n = nodes.size - 1
    K = np.zeros((n + 1, n + 1))
    f = np.zeros(n + 1)
    M = np.zeros((n + 1, n + 1)) if rho is not None else None
    for e in range(n):
        x1, x2 = nodes[e], nodes[e + 1]
        h = x2 - x1
        A_e = float(A(np.array([0.5 * (x1 + x2)]))[0]) if callable(A) else float(A)
        K[e:e + 2, e:e + 2] += element_stiffness(E, A_e, h)
        if b is not None:
            f[e:e + 2] += element_load(b, x1, x2, n_gauss)
        if M is not None:
            M[e:e + 2, e:e + 2] += element_mass(rho, A_e, h)
    return K, f, M


def apply_point_load(f: np.ndarray, node: int, P: float) -> np.ndarray:
    """Return a copy of f with the point load P [N] added at ``node`` (a natural boundary condition)."""
    out = np.asarray(f, dtype=float).copy()
    out[node] += P
    return out


def solve_bar(K: np.ndarray, f: np.ndarray, fixed: dict[int, float]) -> tuple[np.ndarray, np.ndarray]:
    """Solve K u = f + r with the displacements in ``fixed`` ({node: value}) prescribed; return (u, r).

    Partition the dofs into constrained (c) and free (f): K_ff u_f = f_f - K_fc u_c, then the reactions
    at the constrained dofs are r = K u - f (zero at every free dof). Raise ``ValueError`` if nothing is
    fixed (the bar would have a rigid-body mode and K would be singular).
    """
    K = np.asarray(K, dtype=float)
    f = np.asarray(f, dtype=float)
    n = K.shape[0]
    if not fixed:
        raise ValueError("at least one displacement must be prescribed")
    c = np.array(sorted(fixed), dtype=int)
    free = np.setdiff1d(np.arange(n), c)
    u = np.zeros(n)
    u[c] = [fixed[k] for k in c]
    u[free] = np.linalg.solve(K[np.ix_(free, free)], f[free] - K[np.ix_(free, c)] @ u[c])
    r = K @ u - f
    return u, r


def element_stress(u: np.ndarray, nodes: np.ndarray, E: float) -> np.ndarray:
    """Return the (constant) stress E (u_{e+1} - u_e) / h_e in each element: shape (n,)."""
    u = np.asarray(u, dtype=float)
    nodes = np.asarray(nodes, dtype=float)
    return E * np.diff(u) / np.diff(nodes)


def interpolate(u: np.ndarray, nodes: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Evaluate the piecewise-linear finite-element field u_h(x) = sum_j u_j N_j(x) at the points x."""
    return np.interp(np.asarray(x, dtype=float), np.asarray(nodes, dtype=float), np.asarray(u, dtype=float))


def strain_energy(K: np.ndarray, u: np.ndarray) -> float:
    """Return U = 1/2 u^T K u [J], the elastic energy stored in the discrete bar."""
    u = np.asarray(u, dtype=float)
    return float(0.5 * u @ np.asarray(K, dtype=float) @ u)


def potential_energy(K: np.ndarray, f: np.ndarray, u: np.ndarray) -> float:
    """Return Pi(u) = 1/2 u^T K u - f^T u, which the finite-element solution minimises over all admissible u."""
    u = np.asarray(u, dtype=float)
    return strain_energy(K, u) - float(np.asarray(f, dtype=float) @ u)


# ---------------------------------------------------------------- error measures and convergence

def l2_error(u: np.ndarray, nodes: np.ndarray, u_exact: Callable[[np.ndarray], np.ndarray], n_gauss: int = 4) -> float:
    """Return the L2 norm of u_h - u_exact over the bar, integrated element by element with Gauss quadrature.

    Linear elements give an L2 error of O(h^2) for a smooth exact solution.
    """
    nodes = np.asarray(nodes, dtype=float)
    u = np.asarray(u, dtype=float)
    xi, w = np.polynomial.legendre.leggauss(n_gauss)
    total = 0.0
    for e in range(nodes.size - 1):
        x1, x2 = nodes[e], nodes[e + 1]
        h = x2 - x1
        x = 0.5 * (x1 + x2) + 0.5 * h * xi
        uh = shape_functions(xi) @ u[e:e + 2]
        total += 0.5 * h * np.sum(w * (uh - u_exact(x)) ** 2)
    return float(np.sqrt(total))


def energy_error(u: np.ndarray, nodes: np.ndarray, E: float, sigma_exact: Callable[[np.ndarray], np.ndarray],
                 n_gauss: int = 4) -> float:
    """Return the energy-norm error sqrt( int (sigma_h - sigma_exact)^2 / E dx ): O(h) for linear elements.

    The stress of a linear element is constant, so it can only match a varying exact stress to first
    order; this is why displacements converge one order faster than stresses.
    """
    nodes = np.asarray(nodes, dtype=float)
    sig_h = element_stress(u, nodes, E)
    xi, w = np.polynomial.legendre.leggauss(n_gauss)
    total = 0.0
    for e in range(nodes.size - 1):
        x1, x2 = nodes[e], nodes[e + 1]
        h = x2 - x1
        x = 0.5 * (x1 + x2) + 0.5 * h * xi
        total += 0.5 * h * np.sum(w * (sig_h[e] - sigma_exact(x)) ** 2) / E
    return float(np.sqrt(total))


def observed_order(errors: np.ndarray, ratio: float = 2.0) -> np.ndarray:
    """Return log(e_i / e_{i+1}) / log(ratio) for successive errors obtained by refining by ``ratio``."""
    e = np.asarray(errors, dtype=float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


def solve_rotor_blade(case: dict[str, float], n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Solve the rotor-blade spar on n uniform elements; return (nodes, u, sigma_elements, K, f).

    Body force b = rho A omega^2 x, tip load P = m_tip omega^2 L applied at the last node, u(0) = 0.
    """
    E, rho, L, A, om = case["E"], case["rho"], case["L"], case["A"], case["omega"]
    nodes = np.linspace(0.0, L, n + 1)
    K, f, _ = assemble(nodes, E, A, lambda x: rho * A * om ** 2 * x)
    f = apply_point_load(f, n, case["m_tip"] * om ** 2 * L)
    u, _ = solve_bar(K, f, {0: 0.0})
    return nodes, u, element_stress(u, nodes, E), K, f


def convergence_study(case: dict[str, float], ns: list[int]) -> tuple[np.ndarray, np.ndarray]:
    """Return (L2 displacement errors, energy-norm stress errors) of the rotor-blade solution for each n in ``ns``."""
    E, rho, L, A, om = case["E"], case["rho"], case["L"], case["A"], case["omega"]
    P = case["m_tip"] * om ** 2 * L
    u_ex = lambda x: exact_rotating_bar(x, E, rho, om, L, A, P)[0]      # noqa: E731
    s_ex = lambda x: exact_rotating_bar(x, E, rho, om, L, A, P)[1]      # noqa: E731
    e_l2, e_en = [], []
    for n in ns:
        nodes, u, _, _, _ = solve_rotor_blade(case, n)
        e_l2.append(l2_error(u, nodes, u_ex))
        e_en.append(energy_error(u, nodes, E, s_ex))
    return np.array(e_l2), np.array(e_en)


# ---------------------------------------------------------------- vibration

def axial_frequencies(K: np.ndarray, M: np.ndarray, fixed: list[int]) -> np.ndarray:
    """Return the axial natural frequencies omega_k [rad/s] (ascending) of the discrete bar with the given dofs fixed.

    Generalised eigenproblem K phi = omega^2 M phi on the free dofs. The exact frequencies of a uniform
    fixed-free bar are omega_k = (2k - 1) (pi / 2L) sqrt(E / rho); the consistent-mass finite element
    model overestimates them, with an error that falls as h^2 for the low modes.
    """
    n = K.shape[0]
    free = np.setdiff1d(np.arange(n), np.asarray(fixed, dtype=int))
    Kf, Mf = K[np.ix_(free, free)], M[np.ix_(free, free)]
    # Reduce to a standard symmetric problem with the Cholesky factor of M (P03.08): L^{-1} K L^{-T} y = lambda y.
    Lc = np.linalg.cholesky(Mf)
    Linv = np.linalg.inv(Lc)
    lam = np.linalg.eigvalsh(Linv @ Kf @ Linv.T)
    return np.sqrt(np.maximum(lam, 0.0))


def exact_axial_frequencies(E: float, rho: float, L: float, n_modes: int) -> np.ndarray:
    """Return omega_k = (2k - 1) pi sqrt(E / rho) / (2 L), k = 1 .. n_modes, for a uniform fixed-free bar."""
    k = np.arange(1, n_modes + 1)
    return (2 * k - 1) * np.pi / (2.0 * L) * np.sqrt(E / rho)


# ---------------------------------------------------------------- oracle

# GIVEN
def tapered_bar_bvp_oracle(E: float, A0: float, c: float, L: float, P: float, x: np.ndarray) -> np.ndarray:
    """Solve (E A(x) u')' = 0, u(0) = 0, E A(L) u'(L) = P with scipy's collocation BVP solver; return u(x).

    State y = (u, N) with N = E A u' the axial force: u' = N / (E A), N' = 0.
    """
    from scipy.integrate import solve_bvp

    def rhs(xx, y):
        return np.vstack([y[1] / (E * A0 * (1.0 - c * xx / L)), np.zeros_like(xx)])

    def bc(ya, yb):
        return np.array([ya[0], yb[1] - P])

    xs = np.linspace(0.0, L, 50)
    sol = solve_bvp(rhs, bc, xs, np.vstack([np.zeros_like(xs), P * np.ones_like(xs)]), tol=1e-10)
    return sol.sol(np.asarray(x, dtype=float))[0]


# ---------------------------------------------------------------- main

# GIVEN
def main(out_dir: str | Path = ".") -> None:
    """Solve the rotor-blade spar, verify against the closed form, run the convergence and vibration studies, and plot."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    c = rotor_blade()
    E, rho, L, A, om = c["E"], c["rho"], c["L"], c["A"], c["omega"]
    P = c["m_tip"] * om ** 2 * L
    print(f"rotor-blade spar: L = {L} m, A = {A*1e6:.0f} mm^2, omega = {om:.2f} rad/s ({om*60/(2*np.pi):.0f} rpm),"
          f" tip load P = {P/1e3:.2f} kN, total centrifugal load of the spar {0.5*rho*A*om**2*L**2/1e3:.2f} kN")
    for n in (4, 16, 64):
        nodes, u, sig, K, f = solve_rotor_blade(c, n)
        u_ex, s_ex = exact_rotating_bar(nodes, E, rho, om, L, A, P)
        mid = 0.5 * (nodes[:-1] + nodes[1:])
        s_mid = exact_rotating_bar(mid, E, rho, om, L, A, P)[1]
        _, r = solve_bar(K, f, {0: 0.0})
        print(f"  n = {n:3d}: tip u = {u[-1]*1e3:.4f} mm (exact {u_ex[-1]*1e3:.4f}), max |nodal error| = {np.max(np.abs(u - u_ex)):.1e} m,"
              f" root stress = {sig[0]/1e6:.2f} MPa (exact at root {s_ex[0]/1e6:.2f}, at element centre {s_mid[0]/1e6:.2f}),"
              f" reaction = {r[0]/1e3:.3f} kN, strain energy = {strain_energy(K, u):.2f} J")
    print("  (linear elements are exact at the nodes for a uniform bar, so each element stress is the exact"
          " stress averaged over the element: it differs from the centre value by rho omega^2 h^2 / 24)")

    ns = [4, 8, 16, 32, 64]
    e_l2, e_en = convergence_study(c, ns)
    print(f"\nconvergence (n = {ns}):")
    print(f"  L2 displacement error  {np.array2string(e_l2, formatter={'float_kind': lambda v: f'{v:.2e}'})} orders {np.round(observed_order(e_l2), 2)}")
    print(f"  energy-norm stress err {np.array2string(e_en, formatter={'float_kind': lambda v: f'{v:.2e}'})} orders {np.round(observed_order(e_en), 2)}")

    # Tapered strut: a second closed form and the scipy BVP oracle.
    A0, taper, P2, L2 = 8.0e-4, 0.6, 50.0e3, 1.5
    print(f"\ntapered strut: A from {A0*1e6:.0f} to {A0*(1-taper)*1e6:.0f} mm^2 over {L2} m, P = {P2/1e3:.0f} kN")
    for n in (4, 16, 64):
        nodes = np.linspace(0.0, L2, n + 1)
        K, f, _ = assemble(nodes, E, lambda x: A0 * (1.0 - taper * x / L2))
        f = apply_point_load(f, n, P2)
        u, _ = solve_bar(K, f, {0: 0.0})
        u_ex, s_ex = exact_tapered_bar(nodes, E, A0, taper, L2, P2)
        u_bvp = tapered_bar_bvp_oracle(E, A0, taper, L2, P2, nodes)
        print(f"  n = {n:3d}: tip u = {u[-1]*1e3:.4f} mm (exact {u_ex[-1]*1e3:.4f}, scipy BVP {u_bvp[-1]*1e3:.4f}),"
              f" max nodal error {np.max(np.abs(u - u_ex))*1e6:.3f} um, tip stress {element_stress(u, nodes, E)[-1]/1e6:.1f} MPa (exact {s_ex[-1]/1e6:.1f})")

    # Vibration: axial modes of the spar (not spinning).
    print("\naxial natural frequencies of the spar [Hz], consistent mass:")
    w_ex = exact_axial_frequencies(E, rho, L, 3)
    for n in (4, 16, 64):
        nodes = np.linspace(0.0, L, n + 1)
        K, _, M = assemble(nodes, E, A, None, rho)
        w = axial_frequencies(K, M, [0])[:3]
        print(f"  n = {n:3d}: " + ", ".join(f"{wi/(2*np.pi):8.2f}" for wi in w) + "   (exact " + ", ".join(f"{wi/(2*np.pi):8.2f}" for wi in w_ex) + ")")

    # Figure.
    nodes, u, sig, K, f = solve_rotor_blade(c, 8)
    xs = np.linspace(0.0, L, 400)
    u_ex, s_ex = exact_rotating_bar(xs, E, rho, om, L, A, P)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    ax = axes[0, 0]
    ax.plot(xs, u_ex * 1e3, "k-", lw=2, label="exact")
    ax.plot(nodes, u * 1e3, "o--", label="FE, 8 elements")
    ax.set_xlabel("x [m]"); ax.set_ylabel("u [mm]"); ax.set_title("spar displacement"); ax.legend()
    ax = axes[0, 1]
    ax.plot(xs, s_ex / 1e6, "k-", lw=2, label="exact")
    ax.step(nodes, np.append(sig, sig[-1]) / 1e6, where="post", label="FE element stress")
    ax.set_xlabel("x [m]"); ax.set_ylabel("sigma [MPa]"); ax.set_title("spar stress"); ax.legend()
    ax = axes[1, 0]
    ax.loglog(ns, e_l2, "o-", label="L2 displacement error")
    ax.loglog(ns, e_en, "s-", label="energy-norm stress error")
    ax.loglog(ns, e_l2[0] * (ns[0] / np.array(ns)) ** 2, "k:", label="slope -2")
    ax.loglog(ns, e_en[0] * (ns[0] / np.array(ns)), "k-.", label="slope -1")
    ax.set_xlabel("elements n"); ax.set_title("convergence"); ax.legend(fontsize=8)
    ax = axes[1, 1]
    for n in (4, 8, 16, 32):
        nodes_n = np.linspace(0.0, L, n + 1)
        K_n, _, M_n = assemble(nodes_n, E, A, None, rho)
        w = axial_frequencies(K_n, M_n, [0])[:4]
        ax.plot(np.arange(1, 5), w / exact_axial_frequencies(E, rho, L, 4), "o-", label=f"n = {n}")
    ax.axhline(1.0, color="k", lw=0.5)
    ax.set_xlabel("mode"); ax.set_ylabel("omega_FE / omega_exact"); ax.set_title("axial frequencies, consistent mass"); ax.legend(fontsize=8)
    fig.suptitle("P03.13  One-dimensional bar finite elements")
    fig.tight_layout()
    fig.savefig(out_dir / "bar_fem.png", dpi=150)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'bar_fem.png'}")


if __name__ == "__main__":
    main()
