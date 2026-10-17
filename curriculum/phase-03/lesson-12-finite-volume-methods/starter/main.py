"""P03.12 Finite Volume Methods: starter code.

A slug of off-mixture propellant travels down a feed line: linear advection u_t + a u_x = 0 on a periodic
domain of length L, discretised in conservation form on n cells of width dx,

    u_i^{m+1} = u_i^m - (dt/dx) (F_{i+1/2} - F_{i-1/2}),

with the numerical flux F_{i+1/2} at each cell interface written from scratch for four schemes:
upwind (first order, monotone), Lax-Friedrichs (first order, more diffusive), Lax-Wendroff (second order,
oscillatory at discontinuities) and a minmod-limited second-order upwind scheme (TVD). The Courant number
nu = a dt / dx governs stability (|nu| <= 1) and accuracy. The exact solution is a shift of the initial
data, so every error can be measured; the modified equation of the upwind scheme predicts its numerical
diffusion; von Neumann analysis predicts its amplification factors, which the eigenvalues of
``scipy.linalg.circulant`` confirm; ``scipy.integrate.quad`` checks the cell averages and
``scipy.integrate.solve_ivp`` the semi-discrete limit.

Conventions: cell i covers [i dx, (i+1) dx) with centre x_i = (i + 1/2) dx; u_i is the cell AVERAGE, not a
point value; interface i+1/2 lies between cells i and i+1; periodic indices (cell n is cell 0).
Lengths in metres, time in seconds, u dimensionless (a mass or mixture-ratio fraction).

Run:  python main.py          (prints the verification tables, writes fv_advection.png)
Test: pytest curriculum/phase-03/lesson-12-finite-volume-methods
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np

SCHEMES = ("upwind", "lax_friedrichs", "lax_wendroff", "minmod")


# ---------------------------------------------------------------- the engineering case

def feed_line() -> dict[str, float]:
    """Return the feed-line case: a slug of off-mixture propellant convected toward the injector.

    a [m/s] bulk velocity in the line; L [m] length of the periodic model domain; x_pulse [m] centre of the
    slug at t = 0; sigma [m] standard deviation of the Gaussian slug; width [m] of the square (sharp-edged)
    slug; t_final [s] time at which the slug has travelled a L/2 = 1.2 m to the sensor at x_sensor [m].
    The numbers are illustrative of a small liquid engine's feed system, not a specific vehicle.
    """
    return {"a": 12.0, "L": 2.4, "x_pulse": 0.6, "sigma": 0.1, "width": 0.3, "t_final": 0.1, "x_sensor": 1.8}


# ---------------------------------------------------------------- grid, data, exact solution

def cell_centres(L: float, n: int) -> tuple[np.ndarray, float]:
    """Return (x, dx): the n cell centres x_i = (i + 1/2) dx on [0, L) and the cell width dx = L / n."""
    raise NotImplementedError


def cell_averages(f: Callable[[np.ndarray], np.ndarray], L: float, n: int, n_gauss: int = 4) -> np.ndarray:
    """Return the exact cell averages (1/dx) int f dx over each of the n cells by Gauss-Legendre quadrature.

    A finite-volume method evolves cell averages, so the initial data and the starter code must be
    averaged too; sampling f at cell centres is wrong at O(dx^2) and would spoil a second-order study.
    n_gauss points integrate polynomials of degree 2 n_gauss - 1 exactly.
    """
    raise NotImplementedError


def gaussian_pulse(x: np.ndarray, x0: float, sigma: float) -> np.ndarray:
    """Return exp(-(x - x0)^2 / (2 sigma^2)): a smooth slug of unit peak."""
    raise NotImplementedError


def square_pulse(x: np.ndarray, x0: float, width: float) -> np.ndarray:
    """Return 1 where |x - x0| < width / 2 and 0 elsewhere: a sharp-edged slug."""
    raise NotImplementedError


def exact_advection(f: Callable[[np.ndarray], np.ndarray], x: np.ndarray, a: float, t: float, L: float) -> np.ndarray:
    """Return the exact solution f((x - a t) mod L) of u_t + a u_x = 0 with periodic data f on [0, L).

    Linear advection translates the initial profile at speed a without changing its shape: every
    deformation a scheme produces is numerical error.
    """
    raise NotImplementedError


def courant_number(a: float, dt: float, dx: float) -> float:
    """Return nu = a dt / dx, the fraction of a cell the solution moves per step (signed with a)."""
    raise NotImplementedError


# ---------------------------------------------------------------- numerical fluxes and the update

def minmod(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Return the minmod limiter: the smaller-magnitude argument when both have the same sign, else zero."""
    raise NotImplementedError


def numerical_flux(u: np.ndarray, a: float, nu: float, scheme: str) -> np.ndarray:
    """Return the flux F_{i+1/2}, i = 0 .. n-1, through the interface between cells i and i+1 (periodic).

    upwind:          F = a u_i for a > 0 (a u_{i+1} for a < 0): take the value from where the flow comes.
    lax_friedrichs:  F = a (u_i + u_{i+1}) / 2 - (a / (2 nu)) (u_{i+1} - u_i): central flux plus the largest
                     diffusion that keeps the scheme monotone at |nu| <= 1.
    lax_wendroff:    F = a (u_i + u_{i+1}) / 2 - (a nu / 2) (u_{i+1} - u_i): central flux plus exactly the
                     diffusion that cancels the O(dt) error, giving second order and dispersive wiggles.
    minmod:          upwind flux of a piecewise-linear reconstruction with minmod-limited slopes,
                     F = a (u_i + (1 - nu) sigma_i / 2), sigma_i = minmod(u_i - u_{i-1}, u_{i+1} - u_i):
                     second order where the solution is smooth, first order and monotone at extrema and jumps.
    Raise ``ValueError`` for an unknown scheme or for nu = 0 with Lax-Friedrichs.
    """
    raise NotImplementedError


def fv_step(u: np.ndarray, a: float, dt: float, dx: float, scheme: str) -> np.ndarray:
    """Advance the cell averages one step in conservation form: u_i - (dt/dx)(F_{i+1/2} - F_{i-1/2}).

    Because every interface flux is added to one cell and subtracted from its neighbour, the sum of the
    cell averages (the total amount of u) is conserved to rounding whatever the flux formula.
    Raise ``ValueError`` if |nu| > 1: no explicit scheme on this stencil is stable past the CFL limit.
    """
    raise NotImplementedError


def advect(u0: np.ndarray, a: float, dx: float, dt: float, n_steps: int, scheme: str) -> np.ndarray:
    """March the cell averages ``u0`` for ``n_steps`` steps of ``dt`` with the named scheme."""
    raise NotImplementedError


def sensor_history(u0: np.ndarray, a: float, dx: float, dt: float, n_steps: int, scheme: str,
                   i_sensor: int) -> np.ndarray:
    """Return the value in cell ``i_sensor`` after each of ``n_steps`` steps: what a sensor there would log."""
    raise NotImplementedError


def arrival_time(history: np.ndarray, dt: float, threshold: float) -> float:
    """Return the first time at which ``history`` crosses ``threshold`` (linear interpolation), or nan."""
    raise NotImplementedError


# ---------------------------------------------------------------- analysis

def amplification_factor(scheme: str, nu: float, theta: np.ndarray) -> np.ndarray:
    """Return the complex von Neumann factor g(theta) of one step for the mode e^{i j theta} (linear schemes).

    upwind:         g = 1 - nu (1 - e^{-i theta}),      |g|^2 = 1 - 2 nu (1 - nu)(1 - cos theta)
    lax_friedrichs: g = cos theta - i nu sin theta,      |g|^2 = cos^2 theta + nu^2 sin^2 theta
    lax_wendroff:   g = 1 - i nu sin theta - nu^2 (1 - cos theta),  |g|^2 = 1 - nu^2 (1 - nu^2)(1 - cos theta)^2
    All three satisfy |g| <= 1 iff |nu| <= 1. The exact factor is e^{-i nu theta}: unit modulus, phase
    -nu theta. The minmod scheme is nonlinear and has no amplification factor; raise ``ValueError``.
    """
    raise NotImplementedError


def phase_speed_ratio(scheme: str, nu: float, theta: np.ndarray) -> np.ndarray:
    """Return the numerical phase speed over the exact one, -arg g / (nu theta): 1 means no dispersion.

    Lax-Wendroff is slower than exact for every resolved mode (its wiggles trail the front); upwind is
    slower than exact for nu < 1/2 and faster for nu > 1/2 (and exact at nu = 1/2 and nu = 1).
    """
    raise NotImplementedError


def numerical_diffusivity(a: float, dx: float, nu: float) -> float:
    """Return the diffusion coefficient D = a dx (1 - nu) / 2 of the upwind scheme's modified equation.

    Expanding the upwind stencil in Taylor series shows it solves u_t + a u_x = D u_xx + O(dx^2) rather than
    the advection equation: a Gaussian slug of variance sigma^2 therefore spreads to sigma^2 + 2 D t.
    At nu = 1 the diffusion vanishes and the scheme is exact (it shifts the data by one cell).
    """
    raise NotImplementedError


def pulse_variance(x: np.ndarray, u: np.ndarray, dx: float, L: float) -> float:
    """Return the variance of the periodic distribution u(x) about its circular mean (u >= 0, not all zero)."""
    raise NotImplementedError


def total_variation(u: np.ndarray) -> float:
    """Return TV(u) = sum |u_{i+1} - u_i| over the periodic array: new wiggles raise it, a monotone scheme never does."""
    raise NotImplementedError


def l1_error(u: np.ndarray, v: np.ndarray, dx: float) -> float:
    """Return the discrete L1 norm sum |u_i - v_i| dx: the natural norm for conservation laws."""
    raise NotImplementedError


def observed_order(errors: np.ndarray, ratio: float = 2.0) -> np.ndarray:
    """Return log(e_i / e_{i+1}) / log(ratio) for successive errors obtained by refining by ``ratio``."""
    raise NotImplementedError


def convergence_study(scheme: str, ns: list[int], f: Callable[[np.ndarray], np.ndarray], a: float, L: float,
                      t_final: float, nu: float) -> np.ndarray:
    """Return the L1 error against the exact (cell-averaged) solution at t_final for each cell count in ``ns``.

    The Courant number is held fixed, so dt = nu dx / a and the step count grows as n. The initial data and
    the reference are both cell averages (``cell_averages``). t_final is rounded to a whole number of steps
    and the reference is evaluated at the time actually reached.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- oracle

def semi_discrete_upwind_oracle(u0: np.ndarray, a: float, dx: float, t_final: float) -> np.ndarray:
    """Integrate the semi-discrete upwind system du_i/dt = -a (u_i - u_{i-1}) / dx with scipy's RK45 (rtol 1e-10).

    This is the nu -> 0 limit of the fully discrete upwind scheme (method of lines); the fully discrete
    scheme at Courant number nu differs from it by O(nu), which the test checks.
    """
    from scipy.integrate import solve_ivp

    def rhs(_t: float, u: np.ndarray) -> np.ndarray:
        return -a * (u - np.roll(u, 1)) / dx

    sol = solve_ivp(rhs, (0.0, t_final), np.asarray(u0, dtype=float), method="RK45", rtol=1e-10, atol=1e-12)
    return sol.y[:, -1]


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Run the feed-line slug with the four schemes, the sensor study, the stability and convergence studies, and plot."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    c = feed_line()
    a, L, t_final = c["a"], c["L"], c["t_final"]
    n, nu = 120, 0.5
    x, dx = cell_centres(L, n)
    dt = nu * dx / a
    n_steps = int(round(t_final / dt))
    print(f"feed line: a = {a} m/s, L = {L} m, n = {n} cells (dx = {dx*100:.1f} cm), nu = {nu}, dt = {dt*1e3:.3f} ms,"
          f" {n_steps} steps to t = {t_final} s (slug travels {a*t_final:.2f} m)")

    gauss = lambda xx: gaussian_pulse(xx, c["x_pulse"], c["sigma"])          # noqa: E731
    square = lambda xx: square_pulse(xx, c["x_pulse"], c["width"])           # noqa: E731
    results = {}
    for name, f in (("Gaussian", gauss), ("square", square)):
        u0 = cell_averages(f, L, n)
        ref = cell_averages(lambda xx: exact_advection(f, xx, a, n_steps * dt, L), L, n)
        print(f"\n{name} slug after {n_steps} steps: L1 error, peak (exact {np.max(ref):.3f}), min, total variation"
              f" (initial {total_variation(u0):.3f}), mass change")
        for scheme in SCHEMES:
            u = advect(u0, a, dx, dt, n_steps, scheme)
            results[(name, scheme)] = (u, ref)
            print(f"  {scheme:15s} L1 = {l1_error(u, ref, dx):.4f}, peak = {np.max(u):.3f}, min = {np.min(u):+.4f},"
                  f" TV = {total_variation(u):.3f}, sum(u) dx change = {np.sum(u - u0)*dx:+.1e}")

    # Modified equation: the upwind scheme spreads a Gaussian like diffusion with D = a dx (1 - nu) / 2.
    u0 = cell_averages(gauss, L, n)
    u = advect(u0, a, dx, dt, n_steps, "upwind")
    D = numerical_diffusivity(a, dx, nu)
    var0, var1 = pulse_variance(x, u0, dx, L), pulse_variance(x, u, dx, L)
    print(f"\nupwind numerical diffusivity D = a dx (1 - nu)/2 = {D:.4f} m^2/s; variance growth predicted 2 D t = {2*D*n_steps*dt:.5f} m^2,"
          f" measured {var1 - var0:.5f} m^2 (sigma {np.sqrt(var0):.4f} -> {np.sqrt(var1):.4f} m)")
    u1 = advect(u0, a, dx, dx / a, int(round(t_final * a / dx)), "upwind")
    print(f"at nu = 1 the upwind scheme is exact: L1 error {l1_error(u1, cell_averages(lambda xx: exact_advection(gauss, xx, a, t_final, L), L, n), dx):.1e}")

    # Sensor trace for the square slug: what would a mixture-ratio sensor 1.2 m downstream read?
    i_s = int(c["x_sensor"] / dx)
    u0s = cell_averages(square, L, n)
    t_arrive_exact = (c["x_sensor"] - (c["x_pulse"] + 0.5 * c["width"])) / a
    print(f"\nsensor in cell {i_s} (x = {x[i_s]:.3f} m): exact slug front arrives at {t_arrive_exact*1e3:.2f} ms, exact peak 1.000")
    for scheme in SCHEMES:
        h = sensor_history(u0s, a, dx, dt, 2 * n_steps, scheme, i_s)
        print(f"  {scheme:15s} 10 % arrival {arrival_time(h, dt, 0.1)*1e3:6.2f} ms, 50 % arrival {arrival_time(h, dt, 0.5)*1e3:6.2f} ms,"
              f" peak read {np.max(h):.3f}, min read {np.min(h):+.3f}")

    # Stability: von Neumann against circulant eigenvalues, and the CFL limit.
    print("\nmax |g| over theta (von Neumann) for nu = 0.5, 1.0, 1.1:")
    th = np.linspace(0.0, np.pi, 1001)
    for scheme in ("upwind", "lax_friedrichs", "lax_wendroff"):
        print(f"  {scheme:15s} " + ", ".join(f"nu = {v}: {np.max(np.abs(amplification_factor(scheme, v, th))):.4f}" for v in (0.5, 1.0, 1.1)))

    # Convergence orders.
    ns = [240, 480, 960, 1920]
    print(f"\nL1 convergence on the Gaussian slug at nu = {nu} (n = {ns}):")
    for scheme in SCHEMES:
        e = convergence_study(scheme, ns, gauss, a, L, t_final, nu)
        print(f"  {scheme:15s} errors {np.array2string(e, formatter={'float_kind': lambda v: f'{v:.2e}'})} orders {np.round(observed_order(e), 2)}")
    print(f"L1 convergence on the square slug at nu = {nu}:")
    for scheme in ("upwind", "lax_wendroff", "minmod"):
        e = convergence_study(scheme, ns, square, a, L, t_final, nu)
        print(f"  {scheme:15s} errors {np.array2string(e, formatter={'float_kind': lambda v: f'{v:.2e}'})} orders {np.round(observed_order(e), 2)}")

    # Figure.
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    for ax, name in zip(axes[0], ("Gaussian", "square")):
        _, ref = results[(name, "upwind")]
        ax.plot(x, ref, "k-", lw=2, label="exact (cell averages)")
        for scheme in SCHEMES:
            ax.plot(x, results[(name, scheme)][0], "--", lw=1, label=scheme)
        ax.set_xlim(1.2, 2.4); ax.set_xlabel("x [m]"); ax.set_ylabel("u"); ax.set_title(f"{name} slug at t = {t_final} s, n = {n}, nu = {nu}")
        ax.legend(fontsize=8)
    ax = axes[1, 0]
    for scheme in ("upwind", "lax_friedrichs", "lax_wendroff"):
        ax.plot(th, np.abs(amplification_factor(scheme, nu, th)), label=f"|g| {scheme}")
        ax.plot(th[1:], phase_speed_ratio(scheme, nu, th[1:]), ":", label=f"phase speed ratio {scheme}")
    ax.set_ylim(0, 1.3); ax.set_xlabel("theta = k dx"); ax.set_title(f"von Neumann analysis, nu = {nu}"); ax.legend(fontsize=7)
    ax = axes[1, 1]
    for scheme in SCHEMES:
        e = convergence_study(scheme, ns, gauss, a, L, t_final, nu)
        ax.loglog(ns, e, "o-", label=scheme)
    ax.loglog(ns, e[0] * (ns[0] / np.array(ns)) ** 1, "k:", label="slope -1")
    ax.loglog(ns, e[0] * (ns[0] / np.array(ns)) ** 2, "k-.", label="slope -2")
    ax.set_xlabel("cells n"); ax.set_ylabel("L1 error"); ax.set_title("convergence, Gaussian slug"); ax.legend(fontsize=8)
    fig.suptitle("P03.12  Finite-volume advection of a propellant slug")
    fig.tight_layout()
    fig.savefig(out_dir / "fv_advection.png", dpi=150)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'fv_advection.png'}")


if __name__ == "__main__":
    main()
