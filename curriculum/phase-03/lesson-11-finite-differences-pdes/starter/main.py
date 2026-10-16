"""P03.11 Finite Differences for PDEs: starter code.

Transient conduction through a titanium bulkhead, T_t = alpha T_xx on 0 < x < L with both faces at
fixed temperatures, marched in time by three finite-difference schemes written from scratch:

    FTCS            (forward time, centred space; explicit, stable only for r <= 1/2),
    BTCS            (backward time, centred space; implicit, unconditionally stable, first order in dt),
    Crank-Nicolson  (trapezoidal in time; implicit, unconditionally stable, second order in dt),

with r = alpha dt / dx^2 the diffusion number. The implicit schemes solve a tridiagonal system by the
Thomas algorithm, also written from scratch. The exact series solution of the step problem, the exact
solution of the *semi-discrete* (method-of-lines) system, von Neumann amplification factors and the
eigenvalues of the update matrix are the analytic references; ``scipy.integrate.solve_ivp``,
``scipy.linalg.solve_banded`` and ``scipy.linalg.expm`` are the oracles.

Conventions: grids have n cells and n + 1 nodes x_j = j dx, j = 0 .. n, with the Dirichlet values held
in u[0] and u[-1]; the interior unknowns are u[1:-1]. Temperatures in kelvin, lengths in metres,
times in seconds.

Run:  python main.py          (prints the verification tables, writes heat_fd_schemes.png)
Test: pytest curriculum/phase-03/lesson-11-finite-differences-pdes
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

SCHEMES = ("ftcs", "btcs", "cn")


# ---------------------------------------------------------------- the engineering case

def titanium_bulkhead() -> dict[str, float]:
    """Return the bulkhead case: a 10 mm Ti-6Al-4V plate between a hot-gas duct and a cooled bay.

    k [W/(m K)], rho [kg/m^3], cp [J/(kg K)] are handbook values for Ti-6Al-4V near 500 K; L [m] is the
    plate thickness; T0 [K] the uniform initial temperature; T_hot [K] the duct-side face temperature
    applied at t = 0 and held; T_cold [K] the coolant-side face temperature; T_limit [K] the allowable
    at the mid-plane. The case is illustrative, not a flight-cleared design.
    """
    return {"k": 7.0, "rho": 4430.0, "cp": 560.0, "L": 0.010, "T0": 300.0, "T_hot": 700.0,
            "T_cold": 300.0, "T_limit": 450.0}


def thermal_diffusivity(k: float, rho: float, cp: float) -> float:
    """Return alpha = k / (rho cp) [m^2/s]. Raise ``ValueError`` unless all three are positive."""
    raise NotImplementedError


# ---------------------------------------------------------------- exact solutions

def heat_exact_mode(x: np.ndarray, t: float, L: float, alpha: float, m: int = 1) -> np.ndarray:
    """Return the single Fourier mode sin(m pi x / L) exp(-alpha (m pi / L)^2 t), zero at both ends.

    Each mode of the Dirichlet problem decays at its own rate alpha (m pi / L)^2; this is the solution
    used for convergence studies because it has no Gibbs transient and a known decay rate.
    """
    raise NotImplementedError


def heat_exact_step(x: np.ndarray, t: float, L: float, alpha: float, T0: float, T_left: float,
                    T_right: float, n_terms: int = 200) -> np.ndarray:
    """Return T(x, t) for a slab initially at T0 whose faces are held at T_left (x = 0) and T_right (x = L).

    T = T_left + (T_right - T_left) x / L + sum_n b_n sin(n pi x / L) exp(-alpha (n pi / L)^2 t), with
    b_n = (2 / (n pi)) [ (T0 - T_left)(1 - (-1)^n) + (T_right - T_left)(-1)^n ]  (P02.11 separation of
    variables, coefficients by hand). At t = 0 the series converges slowly (Gibbs); for t > 0 the
    exponential makes a few hundred terms exact to rounding.
    """
    raise NotImplementedError


def cold_face_flux_exact(t: float, L: float, alpha: float, k_cond: float, T0: float, T_left: float,
                         T_right: float, n_terms: int = 400) -> float:
    """Return the heat flux q = -k dT/dx [W/m^2] leaving the slab through the x = L face at time t.

    Differentiate the series of ``heat_exact_step`` term by term: dT/dx(L) = (T_right - T_left)/L
    + sum_n b_n (n pi / L) (-1)^n exp(-alpha (n pi / L)^2 t). Positive q means heat flows out of the
    slab into the coolant (from hot to cold when T_left > T_right).
    """
    raise NotImplementedError


def semi_discrete_exact(x: np.ndarray, t: float, L: float, alpha: float, m: int = 1) -> np.ndarray:
    """Return the exact solution of the *spatially discretised* system for one mode: the method-of-lines ODE.

    On a uniform grid with spacing dx = x[1] - x[0], the centred second difference turns sin(k x_j) into
    -lambda_h sin(k x_j) with lambda_h = (4 / dx^2) sin^2(k dx / 2), so the semi-discrete solution is
    sin(k x_j) exp(-alpha lambda_h t). Comparing a time-marcher with this isolates the *time* error from
    the spatial one; comparing it with ``heat_exact_mode`` isolates the spatial error.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- the stencils

def second_difference(u: np.ndarray, dx: float) -> np.ndarray:
    """Return the centred second difference (u[j-1] - 2 u[j] + u[j+1]) / dx^2 at the interior nodes 1 .. n-1.

    Its truncation error is (dx^2 / 12) u'''' + O(dx^4): second order. Raise ``ValueError`` for dx <= 0
    or fewer than three nodes.
    """
    raise NotImplementedError


def solve_tridiagonal(lower: np.ndarray, diag: np.ndarray, upper: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """Solve the tridiagonal system with sub-diagonal ``lower`` (length n-1), diagonal ``diag`` (n) and
    super-diagonal ``upper`` (n-1) by the Thomas algorithm: forward elimination then back substitution,
    O(n) operations and no pivoting (safe for the diagonally dominant matrices of BTCS and Crank-Nicolson).

    Raise ``ValueError`` on inconsistent lengths or a zero pivot.
    """
    raise NotImplementedError


def step_ftcs(u: np.ndarray, r: float) -> np.ndarray:
    """Advance one FTCS step: u_j^{n+1} = u_j + r (u_{j-1} - 2 u_j + u_{j+1}), Dirichlet ends unchanged.

    r = alpha dt / dx^2. Explicit; von Neumann analysis gives g = 1 - 4 r sin^2(theta/2), so the
    scheme is stable only for r <= 1/2. Raise ``ValueError`` if r < 0.
    """
    raise NotImplementedError


def step_btcs(u: np.ndarray, r: float) -> np.ndarray:
    """Advance one BTCS (implicit Euler) step: (1 + 2r) u_j^{n+1} - r (u_{j-1}^{n+1} + u_{j+1}^{n+1}) = u_j^n.

    Tridiagonal solve for the interior with the Dirichlet values moved to the right-hand side.
    g = 1 / (1 + 4 r sin^2(theta/2)): stable for every r, first order in dt.
    """
    raise NotImplementedError


def step_crank_nicolson(u: np.ndarray, r: float) -> np.ndarray:
    """Advance one Crank-Nicolson step: the average of FTCS and BTCS right-hand sides.

    (1 + r) u_j^{n+1} - (r/2)(u_{j-1}^{n+1} + u_{j+1}^{n+1}) = (1 - r) u_j^n + (r/2)(u_{j-1}^n + u_{j+1}^n).
    g = (1 - 2 r sin^2(theta/2)) / (1 + 2 r sin^2(theta/2)): |g| <= 1 for every r, second order in dt,
    but g -> -1 for the grid-scale mode as r grows, so sharp initial data ring instead of being damped.
    """
    raise NotImplementedError


def march(u0: np.ndarray, dx: float, dt: float, alpha: float, n_steps: int, scheme: str = "cn") -> np.ndarray:
    """March ``u0`` forward ``n_steps`` steps of ``dt`` with the named scheme ('ftcs', 'btcs' or 'cn').

    The end values u0[0], u0[-1] are Dirichlet data and are held. Raise ``ValueError`` for an unknown
    scheme, and for FTCS when r > 1/2 (an unstable run is an error, not a result).
    """
    raise NotImplementedError


# ---------------------------------------------------------------- stability analysis

def amplification_factor(scheme: str, r: float, theta: np.ndarray) -> np.ndarray:
    """Return the von Neumann amplification factor g(theta) of one step for the Fourier mode e^{i j theta}.

    FTCS: 1 - 4 r s;  BTCS: 1 / (1 + 4 r s);  Crank-Nicolson: (1 - 2 r s) / (1 + 2 r s), s = sin^2(theta/2).
    Stability requires |g| <= 1 for every theta in [0, pi].
    """
    raise NotImplementedError


def max_stable_dt(scheme: str, dx: float, alpha: float) -> float:
    """Return the largest stable time step: dx^2 / (2 alpha) for FTCS, infinity for BTCS and Crank-Nicolson."""
    raise NotImplementedError


def update_matrix(scheme: str, n_interior: int, r: float) -> np.ndarray:
    """Return the dense matrix G such that u^{n+1} = G u^n on the n_interior unknowns (homogeneous Dirichlet).

    With A the tridiagonal (1, -2, 1) matrix: FTCS G = I + r A; BTCS G = (I - r A)^{-1};
    Crank-Nicolson G = (I - r A / 2)^{-1} (I + r A / 2). The eigenvalues of A are
    -4 sin^2(j pi / (2 (n+1))), j = 1 .. n, so the eigenvalues of G are the von Neumann factors at
    theta_j = j pi / (n+1): matrix and Fourier stability analyses agree for this problem.
    """
    raise NotImplementedError


def spectral_radius(scheme: str, n_interior: int, r: float) -> float:
    """Return max |eigenvalue| of the update matrix: the scheme is stable on this grid iff it is <= 1."""
    raise NotImplementedError


# ---------------------------------------------------------------- convergence studies

def observed_order(errors: np.ndarray, ratio: float = 2.0) -> np.ndarray:
    """Return log(e_i / e_{i+1}) / log(ratio) for successive errors obtained by refining by ``ratio``."""
    raise NotImplementedError


def spatial_convergence(scheme: str, ns: list[int], L: float, alpha: float, t_final: float, r: float,
                        m: int = 1) -> np.ndarray:
    """Return the max-norm error against ``heat_exact_mode`` at t_final for each cell count in ``ns``.

    The diffusion number r is held fixed, so dt = r dx^2 / alpha shrinks with dx^2 and the number of
    steps grows as n^2. For FTCS and Crank-Nicolson the error then falls as dx^2 (second order); for
    BTCS the O(dt) = O(dx^2) time error also gives second order at fixed r. t_final is rounded to a
    whole number of steps and the exact solution is evaluated at the time actually reached.
    """
    raise NotImplementedError


def temporal_convergence(scheme: str, dts: list[float], n: int, L: float, alpha: float, t_final: float,
                         m: int = 1) -> np.ndarray:
    """Return the max-norm error against ``semi_discrete_exact`` at t_final for each time step in ``dts``.

    The grid (n cells) is fixed and the reference is the exact solution of the semi-discrete system, so
    only the time-integration error is measured: first order for BTCS, second for Crank-Nicolson (and
    for FTCS where it is stable). Each dt must divide t_final to within rounding.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- the engineering answers

def time_to_limit(case: dict[str, float], scheme: str, n: int, dt: float, x_probe: float | None = None,
                  t_max: float = 200.0) -> tuple[float, int]:
    """March the bulkhead step problem until the probe temperature first exceeds T_limit.

    Return (time [s] by linear interpolation between the bracketing steps, number of steps taken).
    x_probe defaults to the mid-plane L/2 and must be a node. Raise ``RuntimeError`` if the limit is not
    reached by t_max.
    """
    raise NotImplementedError


def time_to_limit_exact(case: dict[str, float], x_probe: float | None = None) -> float:
    """Return the time at which the exact series solution first reaches T_limit at the probe (bisection)."""
    raise NotImplementedError


# ---------------------------------------------------------------- oracles

def method_of_lines_oracle(u0: np.ndarray, dx: float, alpha: float, t_final: float) -> np.ndarray:
    """Integrate the semi-discrete system du/dt = alpha D2 u (interior nodes) with scipy's BDF to rtol 1e-10.

    This is the library way to solve a stiff method-of-lines system; it is the oracle for the time
    marchers (it has only the spatial error of the same grid).
    """
    from scipy.integrate import solve_ivp

    u0 = np.asarray(u0, dtype=float)
    left, right = u0[0], u0[-1]

    def rhs(_t: float, v: np.ndarray) -> np.ndarray:
        full = np.concatenate([[left], v, [right]])
        return alpha * second_difference(full, dx)

    sol = solve_ivp(rhs, (0.0, t_final), u0[1:-1], method="BDF", rtol=1e-10, atol=1e-12)
    return np.concatenate([[left], sol.y[:, -1], [right]])


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Run the bulkhead case with the three schemes, the stability and convergence studies, and plot."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    c = titanium_bulkhead()
    alpha = thermal_diffusivity(c["k"], c["rho"], c["cp"])
    L = c["L"]
    print(f"Ti-6Al-4V bulkhead: alpha = {alpha:.3e} m^2/s, diffusion time L^2/alpha = {L**2 / alpha:.1f} s")
    n = 100
    dx = L / n
    x = np.linspace(0.0, L, n + 1)
    dt_ftcs = max_stable_dt("ftcs", dx, alpha)
    print(f"n = {n} cells, dx = {dx*1e3:.2f} mm: FTCS needs dt <= {dt_ftcs*1e3:.3f} ms;"
          f" implicit schemes have no limit")

    t_exact = time_to_limit_exact(c)
    print(f"\nexact (series) time for the mid-plane to reach {c['T_limit']:.0f} K: {t_exact:.3f} s")
    for scheme, dt in (("ftcs", 0.9 * dt_ftcs), ("btcs", 0.1), ("cn", 0.1), ("btcs", 0.01), ("cn", 0.01)):
        t_hit, steps = time_to_limit(c, scheme, n, dt)
        print(f"  {scheme:5s} dt = {dt*1e3:8.3f} ms  r = {alpha*dt/dx**2:8.3f}: t = {t_hit:.3f} s"
              f" ({(t_hit - t_exact)*1e3:+.1f} ms) in {steps} steps")

    # Profiles at t = 5 s and the Crank-Nicolson ringing on the step initial condition.
    t_show = 5.0
    u0 = np.full(n + 1, c["T0"])
    u0[0], u0[-1] = c["T_hot"], c["T_cold"]
    exact5 = heat_exact_step(x, t_show, L, alpha, c["T0"], c["T_hot"], c["T_cold"])
    print(f"\nprofiles at t = {t_show} s, max |T - series| in K:")
    profiles = {}
    for scheme, dt in (("ftcs", 0.5 * dt_ftcs), ("btcs", 0.5), ("cn", 0.5), ("cn", 0.05)):
        steps = int(round(t_show / dt))
        u = march(u0, dx, dt, alpha, steps, scheme)
        profiles[(scheme, dt)] = u
        over = np.max(u) - c["T_hot"]
        wiggle = np.max(np.abs(np.diff(u, 2)))              # grid-scale ringing shows up in the second difference
        print(f"  {scheme:5s} dt = {dt*1e3:8.2f} ms r = {alpha*dt/dx**2:7.2f}: error {np.max(np.abs(u - exact5)):.3f} K,"
              f" max(T) - T_hot = {over:+.3f} K, max |second difference| = {wiggle:.2f} K,"
              f" cold-face flux {-c['k']*(u[-1]-u[-2])/dx/1e3:.2f} kW/m^2")
    print(f"  exact cold-face flux at {t_show} s: {cold_face_flux_exact(t_show, L, alpha, c['k'], c['T0'], c['T_hot'], c['T_cold'])/1e3:.2f} kW/m^2,"
          f" steady {c['k']*(c['T_hot']-c['T_cold'])/L/1e3:.2f} kW/m^2")
    mol = method_of_lines_oracle(u0, dx, alpha, t_show)
    print(f"  scipy BDF method of lines vs series: {np.max(np.abs(mol - exact5)):.3f} K (spatial error of the grid);"
          f" vs CN dt = 0.05 s: {np.max(np.abs(mol - profiles[('cn', 0.05)])):.4f} K")

    # Stability: von Neumann against the update-matrix eigenvalues.
    print("\nspectral radius of the update matrix (49 interior nodes) vs max |g| from von Neumann:")
    for scheme in SCHEMES:
        for r in (0.25, 0.5, 0.51, 5.0):
            th = np.arange(1, 50) * np.pi / 50
            print(f"  {scheme:5s} r = {r:4.2f}: rho(G) = {spectral_radius(scheme, 49, r):.6f},"
                  f" max|g| = {np.max(np.abs(amplification_factor(scheme, r, th))):.6f}")

    # Convergence orders.
    ns = [10, 20, 40, 80]
    print("\nspatial convergence at fixed r (mode 1, t = L^2/(alpha pi^2)): max error and observed order")
    t_conv = L ** 2 / (alpha * np.pi ** 2)
    for scheme, r in (("ftcs", 0.4), ("ftcs", 1.0 / 6.0), ("btcs", 0.4), ("cn", 0.5), ("cn", 2.0)):
        e = spatial_convergence(scheme, ns, L, alpha, t_conv, r)
        print(f"  {scheme:5s} r = {r:.3f}: errors {np.array2string(e, formatter={'float_kind': lambda v: f'{v:.2e}'})}"
              f" orders {np.round(observed_order(e), 2)}")
    dts = [t_conv / 5, t_conv / 10, t_conv / 20, t_conv / 40]
    print("temporal convergence at fixed grid (n = 200, reference = exact semi-discrete solution):")
    for scheme in ("btcs", "cn"):
        e = temporal_convergence(scheme, dts, 200, L, alpha, t_conv)
        print(f"  {scheme:5s}: errors {np.array2string(e, formatter={'float_kind': lambda v: f'{v:.2e}'})}"
              f" orders {np.round(observed_order(e), 2)}")

    # Figure.
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    ax = axes[0, 0]
    ax.plot(x * 1e3, exact5, "k-", lw=2, label="series (exact)")
    for (scheme, dt), u in profiles.items():
        ax.plot(x * 1e3, u, "--", lw=1, label=f"{scheme.upper()} dt = {dt*1e3:.2f} ms")
    ax.set_xlabel("x [mm]"); ax.set_ylabel("T [K]"); ax.set_title(f"bulkhead profile at t = {t_show} s"); ax.legend(fontsize=8)
    ax = axes[0, 1]
    th = np.linspace(0.0, np.pi, 200)
    for scheme, r in (("ftcs", 0.25), ("ftcs", 0.5), ("ftcs", 0.6), ("btcs", 5.0), ("cn", 5.0)):
        ax.plot(th, amplification_factor(scheme, r, th), label=f"{scheme.upper()} r = {r}")
    ax.axhline(1, color="k", lw=0.5); ax.axhline(-1, color="k", lw=0.5)
    ax.set_xlabel("theta = k dx"); ax.set_ylabel("g"); ax.set_title("von Neumann amplification factor"); ax.legend(fontsize=8)
    ax = axes[1, 0]
    for scheme, r in (("ftcs", 0.4), ("btcs", 0.4), ("cn", 2.0)):
        e = spatial_convergence(scheme, ns, L, alpha, t_conv, r)
        ax.loglog(np.array(ns), e, "o-", label=f"{scheme.upper()} r = {r}")
    ax.loglog(ns, e[0] * (ns[0] / np.array(ns)) ** 2, "k:", label="slope -2")
    ax.set_xlabel("cells n"); ax.set_ylabel("max error"); ax.set_title("spatial convergence (fixed r)"); ax.legend(fontsize=8)
    ax = axes[1, 1]
    for scheme in ("btcs", "cn"):
        e = temporal_convergence(scheme, dts, 200, L, alpha, t_conv)
        ax.loglog(dts, e, "o-", label=scheme.upper())
    ax.loglog(dts, e[0] * (np.array(dts) / dts[0]) ** 2, "k:", label="slope 2")
    ax.loglog(dts, e[0] * 30 * (np.array(dts) / dts[0]), "k-.", label="slope 1")
    ax.set_xlabel("dt [s]"); ax.set_ylabel("max error vs semi-discrete exact"); ax.set_title("temporal convergence (fixed grid)"); ax.legend(fontsize=8)
    fig.suptitle("P03.11  Finite differences for the heat equation")
    fig.tight_layout()
    fig.savefig(out_dir / "heat_fd_schemes.png", dpi=150)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'heat_fd_schemes.png'}")


if __name__ == "__main__":
    main()
