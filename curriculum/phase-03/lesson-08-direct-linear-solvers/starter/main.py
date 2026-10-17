"""P03.08 Direct Linear Solvers: starter code.

LU factorisation with partial pivoting, Cholesky factorisation, triangular substitutions, the
Thomas algorithm for tridiagonal systems, determinants, condition numbers, backward error and
iterative refinement, all written from scratch in NumPy, with ``scipy.linalg`` as the oracle.

Applications: a planar truss (a cantilevered wing-spar model) whose stiffness matrix is symmetric
positive definite, so Cholesky applies and one factorisation serves every load case; a heated rod
discretised by central differences, whose tridiagonal system the Thomas algorithm solves in O(n)
and whose error shows second-order convergence; and Hilbert matrices, whose conditioning shows
how a tiny backward error becomes a large forward error.

Conventions: P A = L U with L unit lower triangular and P the permutation given by ``perm``
(row i of P A is row perm[i] of A); A = L L^T for Cholesky; norms are 1-norms unless stated.

Run:  python main.py          (prints the verification tables, writes direct_solvers.png)
Test: pytest curriculum/phase-03/lesson-08-direct-linear-solvers
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


# ---------------------------------------------------------------- triangular solves

def forward_substitution(L: np.ndarray, b: np.ndarray, unit_diagonal: bool = False) -> np.ndarray:
    """Solve L y = b for lower-triangular L, row by row: y_i = (b_i - sum_{j<i} L_ij y_j) / L_ii.

    With ``unit_diagonal=True`` the diagonal of L is taken as 1 (the LU convention) and never read.
    Cost n^2 flops. Raises ValueError on a zero diagonal entry.
    """
    raise NotImplementedError


def back_substitution(U: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Solve U x = y for upper-triangular U, from the last row up. Cost n^2 flops."""
    raise NotImplementedError


# ---------------------------------------------------------------- LU with partial pivoting

def lu_decompose(A: np.ndarray, pivot: bool = True) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Gaussian elimination: return (L, U, perm) with A[perm] = L @ U.

    At step k the row with the largest |entry| in column k (rows k..n-1) is swapped into place
    (partial pivoting), the multipliers l_ik = a_ik / a_kk are stored below the diagonal, and the
    trailing block is updated by the rank-1 update a_ij -= l_ik a_kj. Cost 2n^3/3 flops.
    With ``pivot=False`` no swapping is done and a zero pivot raises ValueError; with pivoting a
    pivot below 1e-14 * max|A| raises ValueError (numerically singular).
    """
    raise NotImplementedError


def lu_solve(L: np.ndarray, U: np.ndarray, perm: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Solve A x = b from A[perm] = L U: permute b, forward-substitute, back-substitute. Cost 2n^2."""
    raise NotImplementedError


def solve(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    """LU with partial pivoting followed by the two substitutions; b may have several columns."""
    raise NotImplementedError


def determinant_from_lu(U: np.ndarray, perm: np.ndarray) -> float:
    """det A = sign(perm) * prod(diag U). The sign is (-1)^(number of transpositions)."""
    raise NotImplementedError


def inverse_from_lu(L: np.ndarray, U: np.ndarray, perm: np.ndarray) -> np.ndarray:
    """A^{-1} by solving A x_j = e_j for every column: n solves of 2n^2 each, 2n^3 in all."""
    raise NotImplementedError


def growth_factor(A: np.ndarray) -> float:
    """rho = max |U_ij| / max |A_ij| for partial-pivoting LU: how much elimination can amplify entries.
    Bounded by 2^(n-1), reached by ``wilkinson_matrix``; usually O(1) in practice."""
    raise NotImplementedError


def wilkinson_matrix(n: int) -> np.ndarray:
    """1 on the diagonal, -1 below it, 1 in the last column: growth factor exactly 2^(n-1)."""
    W = -np.tril(np.ones((n, n)), -1) + np.eye(n)
    W[:, -1] = 1.0
    return W


# ---------------------------------------------------------------- Cholesky

def cholesky(A: np.ndarray) -> np.ndarray:
    """Lower-triangular L with A = L L^T for symmetric positive definite A. Cost n^3/3 flops.

    Column by column: L_jj = sqrt(A_jj - sum_k L_jk^2), L_ij = (A_ij - sum_k L_ik L_jk) / L_jj.
    No pivoting is needed: the pivots are the squared diagonal of L and stay positive exactly when
    A is positive definite, so a non-positive pivot *is* the test. Raises ValueError if A is not
    symmetric (to 1e-12 relative) or not positive definite.
    """
    raise NotImplementedError


def cholesky_solve(L: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Solve A x = b from A = L L^T: L y = b, then L^T x = y."""
    raise NotImplementedError


# ---------------------------------------------------------------- tridiagonal

def thomas_solve(lower: np.ndarray, diag: np.ndarray, upper: np.ndarray, d: np.ndarray) -> np.ndarray:
    """Solve a tridiagonal system in O(n): lower[i] x_{i-1} + diag[i] x_i + upper[i] x_{i+1} = d_i.

    ``lower[0]`` and ``upper[-1]`` are ignored. This is LU without pivoting specialised to the
    three diagonals (forward sweep eliminates the sub-diagonal, backward sweep substitutes);
    it is safe when the matrix is diagonally dominant or SPD. Raises ValueError on a zero pivot.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- errors and conditioning

def norm1(A: np.ndarray) -> float:
    """Matrix 1-norm: the largest column sum of |A| (vector 1-norm for a 1-D input)."""
    A = np.asarray(A, dtype=float)
    return float(np.max(np.sum(np.abs(A), axis=0))) if A.ndim == 2 else float(np.sum(np.abs(A)))


def condition_number_1(A: np.ndarray) -> float:
    """kappa_1 = ||A||_1 ||A^{-1}||_1, with the inverse from our own LU."""
    raise NotImplementedError


def backward_error(A: np.ndarray, x: np.ndarray, b: np.ndarray) -> float:
    """Normwise relative backward error ||b - A x||_1 / (||A||_1 ||x||_1 + ||b||_1): the size of the
    smallest perturbation of (A, b) for which x is the exact solution. A stable solver gives O(eps)."""
    raise NotImplementedError


def iterative_refinement(A: np.ndarray, L: np.ndarray, U: np.ndarray, perm: np.ndarray, b: np.ndarray,
                         n_iter: int = 3) -> tuple[np.ndarray, np.ndarray]:
    """x <- x + A^{-1}(b - A x) using the existing factors, with the residual in extended precision
    (np.longdouble). Returns (x, residual norms per iteration). Each iteration costs 2n^2."""
    raise NotImplementedError


def hilbert(n: int) -> np.ndarray:
    """H_ij = 1 / (i + j + 1): symmetric positive definite and catastrophically ill-conditioned."""
    i = np.arange(n)
    return 1.0 / (i[:, None] + i[None, :] + 1.0)


# ---------------------------------------------------------------- problem 1: a planar truss

def spar_truss(n_bays: int = 4, span: float = 2.0, depth: float = 0.3) -> tuple[np.ndarray, np.ndarray]:
    """A cantilevered Warren-type truss: two chords, verticals and diagonals, ``n_bays`` bays.

    Nodes 0..n_bays on the lower chord (y = 0) and n_bays+1..2 n_bays+1 on the upper chord (y = depth),
    root at x = 0 (nodes 0 and n_bays+1 fixed), tip at x = span. Returns (nodes (N, 2) in m,
    elements (M, 2) node indices). Members: lower chord, upper chord, a vertical at every node pair, one diagonal per bay.
    """
    xs = np.linspace(0.0, span, n_bays + 1)
    nodes = np.vstack([np.column_stack([xs, np.zeros_like(xs)]), np.column_stack([xs, np.full_like(xs, depth)])])
    nb = n_bays + 1
    elements = []
    for i in range(n_bays):
        elements.append((i, i + 1))                 # lower chord
        elements.append((nb + i, nb + i + 1))       # upper chord
        elements.append((i, nb + i + 1))            # diagonal (lower-left to upper-right)
    for i in range(nb):
        elements.append((i, nb + i))                # verticals (one at the root too)
    return nodes, np.array(elements)


def truss_stiffness(nodes: np.ndarray, elements: np.ndarray, E: float, A: float) -> np.ndarray:
    """Assemble the global stiffness matrix (2N x 2N) of a pin-jointed planar truss.

    Each bar contributes k = (E A / L) [c, s]^T [c, s] to the 2x2 blocks of its two nodes with
    signs (+, -; -, +), c = cos, s = sin of its orientation. K is symmetric and positive
    semi-definite; it becomes positive definite once enough supports are applied.
    """
    raise NotImplementedError


def apply_supports(K: np.ndarray, fixed_dofs: list[int]) -> tuple[np.ndarray, np.ndarray]:
    """Delete the rows and columns of the fixed degrees of freedom; return (K_free, free_dofs)."""
    raise NotImplementedError


def member_forces(nodes: np.ndarray, elements: np.ndarray, E: float, A: float, u: np.ndarray) -> np.ndarray:
    """Axial force in each bar, N_e = (E A / L) [-c, -s, c, s] . u_e (positive in tension)."""
    raise NotImplementedError


def solve_truss(nodes: np.ndarray, elements: np.ndarray, E: float, A: float, fixed_dofs: list[int],
                loads: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Displacements (2N,) for one or several load vectors (2N,) or (2N, m), by Cholesky of the
    reduced stiffness matrix factored once. Returns (u, reactions) with reactions = K u - f."""
    raise NotImplementedError


# ---------------------------------------------------------------- problem 2: heated rod (tridiagonal)

def heated_rod_system(n: int, length: float, k: float, q0: float, T_wall: float
                      ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Central-difference discretisation of -k T'' = q0 sin(pi x / L), T(0) = T(L) = T_wall, on n
    interior points with spacing h = L/(n+1). Returns (x, lower, diag, upper, rhs) for ``thomas_solve``."""
    h = length / (n + 1)
    x = h * np.arange(1, n + 1)
    lower = np.full(n, -k / h ** 2)
    diag = np.full(n, 2 * k / h ** 2)
    upper = np.full(n, -k / h ** 2)
    rhs = q0 * np.sin(np.pi * x / length)
    rhs[0] += k / h ** 2 * T_wall
    rhs[-1] += k / h ** 2 * T_wall
    return x, lower, diag, upper, rhs


def heated_rod_exact(x: np.ndarray, length: float, k: float, q0: float, T_wall: float) -> np.ndarray:
    """T(x) = T_wall + (q0 L^2 / (k pi^2)) sin(pi x / L)."""
    return T_wall + q0 * length ** 2 / (k * np.pi ** 2) * np.sin(np.pi * x / length)


def heated_rod_errors(ns: list[int], length: float = 0.5, k: float = 15.0, q0: float = 2.0e5,
                      T_wall: float = 300.0) -> np.ndarray:
    """Max-norm error of the Thomas-solved rod for each n: the data for the convergence study."""
    raise NotImplementedError


def convergence_order(hs: np.ndarray, errs: np.ndarray) -> float:
    """Least-squares slope of log(err) against log(h)."""
    hs, errs = np.asarray(hs, dtype=float), np.asarray(errs, dtype=float)
    if np.any(errs <= 0):
        raise ValueError("errors must be positive to take logarithms")
    return float(np.polyfit(np.log(hs), np.log(errs), 1)[0])


def main(out_dir: str | Path = ".") -> None:
    """Run the lesson's experiments, print the tables and save ``direct_solvers.png``."""
    import time
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy import linalg as sla

    out_dir = Path(out_dir)
    rng = np.random.default_rng(308)

    # ---- LU and Cholesky against SciPy
    A = rng.standard_normal((6, 6))
    b = rng.standard_normal(6)
    L, U, perm = lu_decompose(A)
    x = lu_solve(L, U, perm, b)
    print(f"6x6 random: ||A[perm] - L U|| = {np.max(np.abs(A[perm] - L @ U)):.2e}, "
          f"||x - scipy|| = {np.max(np.abs(x - sla.solve(A, b))):.2e}, backward error {backward_error(A, x, b):.2e}, "
          f"det {determinant_from_lu(U, perm):.6f} vs {np.linalg.det(A):.6f}")
    print(f"growth factor: random 6x6 {growth_factor(A):.2f}, Wilkinson 6x6 {growth_factor(wilkinson_matrix(6)):.0f} "
          f"(= 2^5), Wilkinson 30x30 {growth_factor(wilkinson_matrix(30)):.3e} (= 2^29 = {2 ** 29:.3e})")

    # ---- the truss
    E, Ab = 70e9, 2.0e-4
    nodes, elements = spar_truss()
    nb = 5
    fixed = [0, 1, 2 * nb, 2 * nb + 1]
    K = truss_stiffness(nodes, elements, E, Ab)
    Kf, free = apply_supports(K, fixed)
    lamK = np.linalg.eigvalsh(Kf)
    print(f"\ntruss: {nodes.shape[0]} nodes, {elements.shape[0]} members, {Kf.shape[0]} free dofs; K symmetric: "
          f"{np.max(np.abs(Kf - Kf.T)):.1e}; eigenvalues of K_free in [{lamK[0]:.3e}, {lamK[-1]:.3e}] N/m, "
          f"kappa_2 = {lamK[-1] / lamK[0]:.3e}; full K has {np.sum(np.abs(np.linalg.eigvalsh(K)) < 1e-6 * lamK[-1])} zero modes")
    tip = 4
    f_tip = np.zeros(2 * nodes.shape[0]); f_tip[2 * tip + 1] = -4e3
    f_dist = np.zeros(2 * nodes.shape[0])
    for i in range(1, nb):
        f_dist[2 * (nb + i) + 1] = -1e3
    loads = np.column_stack([f_tip, f_dist])
    u, R = solve_truss(nodes, elements, E, Ab, fixed, loads)
    N = member_forces(nodes, elements, E, Ab, u[:, 0])
    print(f"tip load 4 kN down: tip deflection {u[2 * tip + 1, 0] * 1e3:.3f} mm, "
          f"reactions at root: node 0 ({R[0, 0] / 1e3:+.2f}, {R[1, 0] / 1e3:+.2f}) kN, node {nb} ({R[2 * nb, 0] / 1e3:+.2f}, {R[2 * nb + 1, 0] / 1e3:+.2f}) kN; "
          f"sum of vertical reactions {np.sum(R[1::2, 0]) / 1e3:.3f} kN")
    i_max = int(np.argmax(np.abs(N)))
    print(f"member forces: max |N| = {np.abs(N[i_max]) / 1e3:.2f} kN in member {i_max} {tuple(int(v) for v in elements[i_max])} "
          f"(stress {N[i_max] / Ab / 1e6:+.1f} MPa); root lower chord {N[0] / 1e3:+.2f} kN, root upper chord {N[1] / 1e3:+.2f} kN, "
          f"root diagonal {N[2] / 1e3:+.2f} kN")
    print(f"distributed 4 x 1 kN on the upper chord: tip deflection {u[2 * tip + 1, 1] * 1e3:.3f} mm")
    Lc = cholesky(Kf)
    print(f"Cholesky: ||K - L L^T|| / ||K|| = {np.max(np.abs(Kf - Lc @ Lc.T)) / np.max(np.abs(Kf)):.2e}; vs scipy.linalg.cholesky "
          f"{np.max(np.abs(Lc - sla.cholesky(Kf, lower=True))):.2e}; solution vs scipy.linalg.solve {np.max(np.abs(u[free, 0] - sla.solve(Kf, f_tip[free], assume_a='pos'))):.2e} m")
    try:
        cholesky(K)
    except ValueError as err:
        print(f"Cholesky of the unsupported K: {err}")

    # ---- heated rod: Thomas and convergence
    ns = [10, 20, 40, 80, 160]
    errs = heated_rod_errors(ns)
    hs = 0.5 / (np.array(ns) + 1)
    print(f"\nheated rod (L = 0.5 m, k = 15 W/m/K, q0 = 2e5 W/m^3): max error [K] for n = {ns}: {np.array2string(errs, precision=3)}; "
          f"observed order {convergence_order(hs, errs):.2f}; peak temperature rise {2e5 * 0.5 ** 2 / (15 * np.pi ** 2):.2f} K")
    n = 600
    x_, lo, di, up, rhs = heated_rod_system(n, 0.5, 15.0, 2e5, 300.0)
    t0 = time.perf_counter(); T_th = thomas_solve(lo, di, up, rhs); t_th = time.perf_counter() - t0
    Ad = np.diag(di) + np.diag(lo[1:], -1) + np.diag(up[:-1], 1)
    t0 = time.perf_counter(); T_lu = solve(Ad, rhs); t_lu = time.perf_counter() - t0
    print(f"n = {n}: Thomas {t_th * 1e3:.2f} ms, dense LU {t_lu * 1e3:.0f} ms, max difference {np.max(np.abs(T_th - T_lu)):.2e} K")

    # ---- Hilbert: backward vs forward error
    print("\nHilbert matrices, x_true = ones:")
    conds, fwd, bwd = [], [], []
    for n in range(4, 11, 2):
        H = hilbert(n)
        x_true = np.ones(n)
        bb = H @ x_true
        L, U, perm = lu_decompose(H)
        x = lu_solve(L, U, perm, bb)
        kappa = condition_number_1(H)
        conds.append(kappa); fwd.append(norm1(x - x_true) / norm1(x_true)); bwd.append(backward_error(H, x, bb))
        xr, hist = iterative_refinement(H, L, U, perm, bb, 3)
        print(f"  n = {n:2d}: kappa_1 = {kappa:.2e} (numpy {np.linalg.cond(H, 1):.2e}), backward error {bwd[-1]:.1e}, "
              f"forward error {fwd[-1]:.1e}, kappa x eps = {kappa * np.finfo(float).eps:.1e}; after refinement {norm1(xr - x_true) / norm1(x_true):.1e}")
    slope = np.polyfit(np.log(conds), np.log(fwd), 1)[0]
    print(f"  slope of log(forward error) vs log(kappa): {slope:.2f} (theory: 1)")

    # ---- cost scaling
    print("\ncost of our LU (pure NumPy, rank-1 updates):")
    ts = []
    for n in (100, 200, 400):
        M = rng.standard_normal((n, n))
        t0 = time.perf_counter(); lu_decompose(M); t1 = time.perf_counter() - t0
        t0 = time.perf_counter(); sla.lu_factor(M); t2 = time.perf_counter() - t0
        ts.append(t1)
        print(f"  n = {n}: ours {t1 * 1e3:7.1f} ms, scipy.linalg.lu_factor {t2 * 1e3:6.2f} ms, 2n^3/3 = {2 * n ** 3 / 3:.2e} flops")
    print(f"  observed time exponent {np.polyfit(np.log([100, 200, 400]), np.log(ts), 1)[0]:.2f} (theory 3 for large n; Python overhead flattens it)")

    # ---- figure
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    ax = axes[0, 0]
    scale = 20.0
    for (i, j) in elements:
        ax.plot(nodes[[i, j], 0], nodes[[i, j], 1], "k-", lw=0.8)
        ax.plot(nodes[[i, j], 0] + scale * u[[2 * i, 2 * j], 0], nodes[[i, j], 1] + scale * u[[2 * i + 1, 2 * j + 1], 0], "r-", lw=1.2)
    ax.set_aspect("equal"); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
    ax.set_title(f"spar truss, 4 kN at the tip (deflection x{scale:.0f})")
    ax = axes[0, 1]
    ax.bar(range(len(N)), N / 1e3)
    ax.set_xlabel("member"); ax.set_ylabel("axial force [kN] (+ tension)"); ax.set_title("member forces, tip load")
    ax = axes[1, 0]
    ax.loglog(hs, errs, "o-", label=f"Thomas, order {convergence_order(hs, errs):.2f}")
    ax.loglog(hs, errs[0] * (hs / hs[0]) ** 2, "k--", lw=0.8, label="h^2")
    ax.set_xlabel("h [m]"); ax.set_ylabel("max |T - exact| [K]"); ax.set_title("heated rod: convergence"); ax.legend()
    ax = axes[1, 1]
    ax.loglog(conds, fwd, "o-", label="forward error")
    ax.loglog(conds, bwd, "s-", label="backward error")
    ax.loglog(conds, np.array(conds) * np.finfo(float).eps, "k--", lw=0.8, label="kappa x eps")
    ax.set_xlabel("kappa_1(H_n)"); ax.set_ylabel("relative error"); ax.set_title("Hilbert systems: conditioning"); ax.legend()
    fig.suptitle("P03.08  Direct linear solvers")
    fig.tight_layout()
    fig.savefig(out_dir / "direct_solvers.png", dpi=150)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'direct_solvers.png'}")


if __name__ == "__main__":
    main()
