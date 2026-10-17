"""P02.11 PDEs: Heat, Wave & Laplace Equations: reference solution.

Separation-of-variables series for the three classical linear PDEs, each checked against an
independent method written from scratch: an explicit finite-difference marcher (heat), d'Alembert's
travelling-wave formula and a leapfrog marcher (wave), and Jacobi relaxation on a grid (Laplace).

The three problems:
    * Heat:    T_t = alpha T_xx on 0 < x < L.  A re-entry tile: hot face held at T_s, back face insulated,
               initial temperature T_0. Also the Dirichlet-Dirichlet slab for general initial data.
    * Wave:    u_tt = c^2 u_xx on 0 < x < L, u = 0 at both ends. A tensioned cable plucked into a
               triangle, c = sqrt(T / mu).
    * Laplace: u_xx + u_yy = 0 on a rectangle 0 < x < a, 0 < y < b, u = 0 on three edges and u = f(x)
               on the edge y = b. Steady temperature of a plate with one heated edge.

Conventions: x, y in metres, t in seconds, temperatures in kelvin, alpha in m^2/s, c in m/s. Sine
coefficients b_n multiply sin(n pi x / L), n = 1, 2, ... .

Run:  python solution.py          (prints the tile, cable and plate results, writes pde_solutions.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-02/lesson-11-pdes-heat-wave-laplace
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
from scipy.optimize import brentq


# GIVEN
def reentry_tile() -> dict[str, float]:
    """Return the tile model used in the lesson (plausible values for a rigid silica tile, not a datasheet).

    k = 0.05 W/(m K), rho = 144 kg/m^3, c_p = 700 J/(kg K) give alpha = k / (rho c_p) = 4.96e-7 m^2/s.
    Thickness 50 mm; hot face held at 1300 K from t = 0; initial temperature 300 K; the aluminium
    structure under the tile must stay below 450 K during a 1200 s heating pulse.
    """
    k, rho, cp = 0.05, 144.0, 700.0
    return {"k": k, "rho": rho, "cp": cp, "alpha": k / (rho * cp), "L": 0.050,
            "T0": 300.0, "Ts": 1300.0, "T_limit": 450.0, "t_pulse": 1200.0}


# GIVEN
def sensor_cable() -> dict[str, float]:
    """Return the tensioned-cable model: L = 1.2 m, tension 200 N, mass per length 0.010 kg/m.

    Plucked at a = 0.3 m from one end to a height h = 5 mm and released from rest.
    """
    return {"L": 1.2, "T": 200.0, "mu": 0.010, "a": 0.3, "h": 0.005}


# ---------------------------------------------------------------- Fourier sine coefficients

def sine_coefficients(f: Callable[[np.ndarray], np.ndarray], L: float, n_terms: int, n_quad: int = 8001) -> np.ndarray:
    """Return b_1 .. b_{n_terms} with b_n = (2/L) int_0^L f(x) sin(n pi x / L) dx, by the trapezoidal rule.

    ``n_quad`` points on [0, L]. Raise ``ValueError`` if L <= 0 or n_terms < 1.
    """
    if L <= 0.0 or n_terms < 1:
        raise ValueError("L must be positive and n_terms >= 1")
    x = np.linspace(0.0, L, n_quad)
    fx = np.asarray(f(x), dtype=float)
    n = np.arange(1, n_terms + 1)
    integrand = fx[None, :] * np.sin(np.outer(n, x) * np.pi / L)
    return (2.0 / L) * np.trapezoid(integrand, x, axis=1)


# ---------------------------------------------------------------- the heat equation

def heat_series(x: np.ndarray, t: float, L: float, alpha: float, b: np.ndarray) -> np.ndarray:
    """Return u(x, t) = sum_n b_n sin(n pi x / L) exp(-alpha (n pi / L)^2 t) for u = 0 at both ends.

    Raise ``ValueError`` for t < 0 (the heat equation cannot be run backwards).
    """
    if t < 0.0:
        raise ValueError("the heat equation is solved forward in time only")
    x = np.asarray(x, dtype=float)
    n = np.arange(1, b.size + 1)
    lam = n * np.pi / L
    return np.sin(np.outer(x, lam)) @ (b * np.exp(-alpha * lam ** 2 * t))


def tile_temperature(x: np.ndarray, t: float, L: float, alpha: float, T0: float, Ts: float, n_terms: int = 50) -> np.ndarray:
    """Return T(x, t) in a slab with T(0, t) = Ts (hot face), dT/dx(L, t) = 0 (insulated back), T(x, 0) = T0.

    theta = T - Ts obeys the heat equation with eigenfunctions sin(lambda_n x), lambda_n = (2n+1) pi / (2L),
    and coefficients 4 (T0 - Ts) / ((2n+1) pi):
        T = Ts + (T0 - Ts) sum_{n>=0} 4/((2n+1) pi) sin(lambda_n x) exp(-alpha lambda_n^2 t).
    """
    if t < 0.0:
        raise ValueError("t must be non-negative")
    x = np.asarray(x, dtype=float)
    n = np.arange(n_terms)
    lam = (2 * n + 1) * np.pi / (2.0 * L)
    coeff = 4.0 / ((2 * n + 1) * np.pi) * np.exp(-alpha * lam ** 2 * t)
    return Ts + (T0 - Ts) * (np.sin(np.outer(x, lam)) @ coeff)


def tile_thickness_for_limit(t: float, alpha: float, T0: float, Ts: float, T_limit: float,
                             L_min: float = 1e-3, L_max: float = 1.0) -> float:
    """Return the slab thickness L [m] at which the insulated back face reaches exactly T_limit at time t.

    Root of tile_temperature(L, t, L, ...) - T_limit in L by Brent's method; raise ``ValueError`` if the
    limit is not bracketed between L_min and L_max (for example T_limit outside (T0, Ts)).
    """
    def g(L: float) -> float:
        return float(tile_temperature(np.array([L]), t, L, alpha, T0, Ts)[0]) - T_limit
    if g(L_min) * g(L_max) > 0.0:
        raise ValueError("T_limit is not reached for any thickness in [L_min, L_max]")
    return float(brentq(g, L_min, L_max, xtol=1e-9))


def heat_ftcs(u0: np.ndarray, dx: float, dt: float, alpha: float, n_steps: int,
              left: tuple[str, float] = ("dirichlet", 0.0), right: tuple[str, float] = ("dirichlet", 0.0)) -> np.ndarray:
    """March u_t = alpha u_xx with the explicit forward-time, centred-space scheme and return u after n_steps.

    u[j]^{m+1} = u[j] + r (u[j+1] - 2 u[j] + u[j-1]),  r = alpha dt / dx^2.  Each boundary is
    ("dirichlet", value) (u fixed) or ("neumann", 0.0) (insulated: mirror the interior neighbour, a
    second-order ghost point). Raise ``ValueError`` if r > 0.5 (the scheme is unstable there).
    """
    r = alpha * dt / dx ** 2
    if r > 0.5:
        raise ValueError(f"FTCS is unstable for r = alpha dt / dx^2 = {r:.3f} > 0.5")
    u = np.array(u0, dtype=float)
    for _ in range(n_steps):
        new = u.copy()
        new[1:-1] = u[1:-1] + r * (u[2:] - 2.0 * u[1:-1] + u[:-2])
        for side, (kind, value) in ((0, left), (-1, right)):
            if kind == "dirichlet":
                new[side] = value
            elif kind == "neumann":
                nb = u[1] if side == 0 else u[-2]
                new[side] = u[side] + 2.0 * r * (nb - u[side])
            else:
                raise ValueError(f"unknown boundary type {kind!r}")
        u = new
    return u


# ---------------------------------------------------------------- the wave equation

def string_mode_frequency(n: int, L: float, T: float, mu: float) -> float:
    """Return f_n = (n / 2L) sqrt(T / mu) [Hz], the n-th natural frequency of a string fixed at both ends."""
    if n < 1 or L <= 0.0 or T <= 0.0 or mu <= 0.0:
        raise ValueError("n >= 1 and L, T, mu > 0 are required")
    return n / (2.0 * L) * np.sqrt(T / mu)


def plucked_string_coefficients(L: float, a: float, h: float, n_terms: int) -> np.ndarray:
    """Return the sine coefficients of the triangle f(x) = h x / a (x < a), h (L - x) / (L - a) (x > a).

    b_n = 2 h L^2 / (n^2 pi^2 a (L - a)) sin(n pi a / L).
    """
    if not 0.0 < a < L:
        raise ValueError("the pluck point must lie strictly inside the string")
    n = np.arange(1, n_terms + 1)
    return 2.0 * h * L ** 2 / (n ** 2 * np.pi ** 2 * a * (L - a)) * np.sin(n * np.pi * a / L)


def wave_series(x: np.ndarray, t: float, L: float, c: float, b_disp: np.ndarray, b_vel: np.ndarray | None = None) -> np.ndarray:
    """Return u(x, t) = sum_n [b_n cos(w_n t) + (d_n / w_n) sin(w_n t)] sin(n pi x / L), w_n = n pi c / L.

    b_disp are the sine coefficients of u(x, 0), b_vel those of u_t(x, 0) (zero if None).
    """
    x = np.asarray(x, dtype=float)
    n = np.arange(1, b_disp.size + 1)
    k = n * np.pi / L
    w = k * c
    amp = b_disp * np.cos(w * t)
    if b_vel is not None:
        amp = amp + b_vel / w * np.sin(w * t)
    return np.sin(np.outer(x, k)) @ amp


def odd_periodic_extension(f: Callable[[np.ndarray], np.ndarray], L: float) -> Callable[[np.ndarray], np.ndarray]:
    """Return F with F(x) = f(x) on [0, L], F(-x) = -F(x), F(x + 2L) = F(x): the string's image function."""
    def F(x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        xr = np.mod(x, 2.0 * L)
        out = np.where(xr <= L, f(np.minimum(xr, L)), -f(2.0 * L - np.maximum(xr, L)))
        return out
    return F


def dalembert(x: np.ndarray, t: float, c: float, L: float, f: Callable[[np.ndarray], np.ndarray]) -> np.ndarray:
    """Return u(x, t) = [F(x - ct) + F(x + ct)] / 2 for a string released from rest in shape f.

    F is the odd 2L-periodic extension of f, which makes the two travelling waves satisfy u = 0 at
    x = 0 and x = L automatically (reflection with inversion).
    """
    F = odd_periodic_extension(f, L)
    x = np.asarray(x, dtype=float)
    return 0.5 * (F(x - c * t) + F(x + c * t))


def wave_leapfrog(u0: np.ndarray, v0: np.ndarray, dx: float, dt: float, c: float, n_steps: int) -> tuple[np.ndarray, np.ndarray]:
    """March u_tt = c^2 u_xx with the centred (leapfrog) scheme, u = 0 at both ends; return (u, u_t) after n_steps.

    u[j]^{m+1} = 2 u[j]^m - u[j]^{m-1} + s^2 (u[j+1] - 2 u[j] + u[j-1])^m,  s = c dt / dx (the Courant number).
    The first step uses a Taylor start, u^1 = u^0 + dt v0 + (s^2 / 2) D2 u^0. u_t is returned as the
    centred difference (u^{m+1} - u^{m-1}) / (2 dt) at the final level. Raise ``ValueError`` if s > 1 (the CFL
    condition).
    """
    s = c * dt / dx
    if s > 1.0:
        raise ValueError(f"leapfrog is unstable for Courant number s = {s:.3f} > 1")
    u_prev = np.array(u0, dtype=float)
    u_prev[0] = u_prev[-1] = 0.0
    d2 = np.zeros_like(u_prev)
    d2[1:-1] = u_prev[2:] - 2.0 * u_prev[1:-1] + u_prev[:-2]
    u = u_prev + dt * np.asarray(v0, dtype=float) + 0.5 * s ** 2 * d2
    u[0] = u[-1] = 0.0
    if n_steps == 0:
        return u_prev, np.asarray(v0, dtype=float)
    for _ in range(n_steps - 1):
        d2[1:-1] = u[2:] - 2.0 * u[1:-1] + u[:-2]
        u_next = 2.0 * u - u_prev + s ** 2 * d2
        u_next[0] = u_next[-1] = 0.0
        u_prev, u = u, u_next
    # one more level to centre the velocity estimate at the returned u
    d2[1:-1] = u[2:] - 2.0 * u[1:-1] + u[:-2]
    u_next = 2.0 * u - u_prev + s ** 2 * d2
    u_next[0] = u_next[-1] = 0.0
    ut = (u_next - u_prev) / (2.0 * dt)
    return u, ut


def string_energy(u: np.ndarray, ut: np.ndarray, dx: float, T: float, mu: float) -> float:
    """Return E = (1/2) int (mu u_t^2 + T u_x^2) dx [J] by the trapezoidal rule, u_x by central differences."""
    ux = np.gradient(np.asarray(u, dtype=float), dx)
    return float(0.5 * np.trapezoid(mu * np.asarray(ut) ** 2 + T * ux ** 2, dx=dx))


# ---------------------------------------------------------------- Laplace's equation

def laplace_rectangle(x: np.ndarray, y: np.ndarray, a: float, b: float, b_top: np.ndarray) -> np.ndarray:
    """Return u(x, y) = sum_n b_n sin(n pi x / a) sinh(n pi y / a) / sinh(n pi b / a) on the a x b rectangle.

    u = 0 on x = 0, x = a, y = 0 and u = f(x) on y = b, where b_top are the sine coefficients of f on [0, a].
    x and y are broadcast against each other (pass meshgrid arrays or a column and a row).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = np.arange(1, b_top.size + 1)
    k = n * np.pi / a
    u = np.zeros(np.broadcast(x, y).shape)
    for kn, bn in zip(k, b_top):
        # sinh ratio written with exponentials to avoid overflow for large n
        ratio = np.exp(kn * (y - b)) * (1.0 - np.exp(-2.0 * kn * y)) / (1.0 - np.exp(-2.0 * kn * b))
        u = u + bn * np.sin(kn * x) * ratio
    return u


def laplacian_residual(u: np.ndarray, dx: float, dy: float) -> np.ndarray:
    """Return the five-point Laplacian (u_xx + u_yy) of the interior of a grid u[i, j] = u(x_i, y_j)."""
    u = np.asarray(u, dtype=float)
    return ((u[2:, 1:-1] - 2.0 * u[1:-1, 1:-1] + u[:-2, 1:-1]) / dx ** 2
            + (u[1:-1, 2:] - 2.0 * u[1:-1, 1:-1] + u[1:-1, :-2]) / dy ** 2)


def laplace_jacobi(nx: int, ny: int, a: float, b: float, f_top: Callable[[np.ndarray], np.ndarray],
                   tol: float = 1e-8, max_iter: int = 100000) -> tuple[np.ndarray, int]:
    """Relax u_xx + u_yy = 0 on an (nx x ny)-point grid by Jacobi iteration; return (u, iterations).

    u = 0 on three edges, u = f_top(x) on y = b. Each sweep replaces every interior value by the weighted
    mean of its four neighbours (dx = dy is not assumed). Stops when the largest change in a sweep is
    below ``tol`` times the largest boundary value; raise ``RuntimeError`` at ``max_iter``.
    """
    x = np.linspace(0.0, a, nx)
    y = np.linspace(0.0, b, ny)
    dx, dy = x[1] - x[0], y[1] - y[0]
    u = np.zeros((nx, ny))
    u[:, -1] = f_top(x)
    u[0, -1] = u[-1, -1] = 0.0
    scale = max(np.max(np.abs(u)), 1e-300)
    cx, cy = 1.0 / dx ** 2, 1.0 / dy ** 2
    for it in range(1, max_iter + 1):
        new = u.copy()
        new[1:-1, 1:-1] = (cx * (u[2:, 1:-1] + u[:-2, 1:-1]) + cy * (u[1:-1, 2:] + u[1:-1, :-2])) / (2.0 * (cx + cy))
        change = np.max(np.abs(new - u))
        u = new
        if change < tol * scale:
            return u, it
    raise RuntimeError("Jacobi iteration did not converge")


# ---------------------------------------------------------------- main

# GIVEN
def main(out_dir: str | Path = ".") -> None:
    """Solve the tile, the cable and the plate; compare each series with its independent check; plot."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)

    # --- the re-entry tile ---------------------------------------------------------------
    tile = reentry_tile()
    alpha, L, T0, Ts, t_p = tile["alpha"], tile["L"], tile["T0"], tile["Ts"], tile["t_pulse"]
    print(f"TILE: alpha = {alpha:.3e} m^2/s, L = {L * 1e3:.0f} mm, diffusion time L^2/alpha = {L ** 2 / alpha:.0f} s,"
          f" Fourier number at {t_p:.0f} s = {alpha * t_p / L ** 2:.4f}")
    T_back = tile_temperature(np.array([L]), t_p, L, alpha, T0, Ts)[0]
    one_term = Ts + (T0 - Ts) * 4.0 / np.pi * np.exp(-alpha * (np.pi / (2 * L)) ** 2 * t_p)
    print(f"back-face temperature after {t_p:.0f} s: series {T_back:.1f} K, one-term {one_term:.1f} K,"
          f" limit {tile['T_limit']:.0f} K")
    L_req = tile_thickness_for_limit(t_p, alpha, T0, Ts, tile["T_limit"])
    print(f"thickness that holds the back face at {tile['T_limit']:.0f} K after {t_p:.0f} s: {L_req * 1e3:.1f} mm")
    nx = 51
    x = np.linspace(0.0, L, nx)
    dx = x[1] - x[0]
    dt = 0.4 * dx ** 2 / alpha
    n_steps = int(np.ceil(t_p / dt))
    dt = t_p / n_steps
    u_fd = heat_ftcs(np.full(nx, T0), dx, dt, alpha, n_steps, left=("dirichlet", Ts), right=("neumann", 0.0))
    u_se = tile_temperature(x, t_p, L, alpha, T0, Ts)
    print(f"FTCS with {nx} points, r = {alpha * dt / dx ** 2:.3f}, {n_steps} steps: max |FD - series| = {np.max(np.abs(u_fd - u_se)):.3f} K,"
          f" back face {u_fd[-1]:.1f} K")
    try:
        heat_ftcs(np.full(nx, T0), dx, 1.3 * dx ** 2 / alpha, alpha, 10)
    except ValueError as err:
        print("FTCS refuses r > 1/2:", err)

    # --- the plucked cable ----------------------------------------------------------------
    cab = sensor_cable()
    Lc, Tc, mu, a, h = cab["L"], cab["T"], cab["mu"], cab["a"], cab["h"]
    c = np.sqrt(Tc / mu)
    print(f"\nCABLE: c = sqrt(T/mu) = {c:.2f} m/s, f1 = {string_mode_frequency(1, Lc, Tc, mu):.2f} Hz,"
          f" f2 = {string_mode_frequency(2, Lc, Tc, mu):.2f} Hz, round trip 2L/c = {2 * Lc / c * 1e3:.2f} ms")
    tri = lambda xx: np.where(xx < a, h * xx / a, h * (Lc - xx) / (Lc - a))
    b_tri = plucked_string_coefficients(Lc, a, h, 60)
    b_num = sine_coefficients(tri, Lc, 60)
    print(f"closed-form vs quadrature coefficients: max difference {np.max(np.abs(b_tri - b_num)):.2e} m;"
          f" b1..b4 = {b_tri[:4] * 1e3} mm")
    xs = np.linspace(0.0, Lc, 401)
    for t_ms in (0.0, 1.0, 2.0, 3.0):
        t = t_ms * 1e-3
        err = np.max(np.abs(wave_series(xs, t, Lc, c, b_tri) - dalembert(xs, t, c, Lc, tri)))
        print(f"t = {t_ms:.0f} ms: max |series(60 terms) - d'Alembert| = {err * 1e3:.4f} mm")
    nxw = 241
    xw = np.linspace(0.0, Lc, nxw)
    dxw = xw[1] - xw[0]
    dtw = 0.8 * dxw / c
    E0 = string_energy(tri(xw), np.zeros(nxw), dxw, Tc, mu)
    for t_ms in (2.0, 10.0):
        n_st = int(round(t_ms * 1e-3 / dtw))
        u_lf, ut_lf = wave_leapfrog(tri(xw), np.zeros(nxw), dxw, dtw, c, n_st)
        E = string_energy(u_lf, ut_lf, dxw, Tc, mu)
        print(f"leapfrog s = {c * dtw / dxw:.2f}, {n_st} steps to {n_st * dtw * 1e3:.2f} ms: max |FD - d'Alembert| ="
              f" {np.max(np.abs(u_lf - dalembert(xw, n_st * dtw, c, Lc, tri))) * 1e3:.4f} mm, energy {E * 1e3:.4f} mJ (initial {E0 * 1e3:.4f} mJ)")

    # --- the heated plate -----------------------------------------------------------------
    ap, bp, dT = 0.30, 0.20, 60.0
    f_top = lambda xx: dT * 4.0 * xx * (ap - xx) / ap ** 2          # parabolic edge profile, b_n = 32 dT / (n pi)^3, odd n
    b_top = sine_coefficients(f_top, ap, 20)
    nxp, nyp = 61, 41
    xp, yp = np.linspace(0.0, ap, nxp), np.linspace(0.0, bp, nyp)
    X, Y = np.meshgrid(xp, yp, indexing="ij")
    u_series = laplace_rectangle(X, Y, ap, bp, b_top)
    u_jac, iters = laplace_jacobi(nxp, nyp, ap, bp, f_top, tol=1e-9)
    print(f"\nPLATE {ap} x {bp} m, top edge peaks {dT:.0f} K above the others (parabola): b1..b5 = {b_top[:5]},"
          f" closed form 32 dT/(n pi)^3 = {32 * dT / (np.arange(1, 6) * np.pi) ** 3}")
    print(f"Jacobi on {nxp} x {nyp}: {iters} sweeps; max |Jacobi - series| = {np.max(np.abs(u_jac - u_series)):.4f} K;"
          f" centre temperature series {laplace_rectangle(ap / 2, bp / 2, ap, bp, b_top):.4f} K")
    res = laplacian_residual(u_series, xp[1] - xp[0], yp[1] - yp[0])
    print(f"five-point Laplacian of the series on the grid: max |residual| = {np.max(np.abs(res)):.3e} K/m^2"
          f" (scale {dT / bp ** 2:.0f})")

    # --- figure ---------------------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    ax = axes[0]
    for tt in (60.0, 300.0, 600.0, 1200.0):
        ax.plot(x * 1e3, tile_temperature(x, tt, L, alpha, T0, Ts), label=f"t = {tt:.0f} s")
    ax.plot(x * 1e3, u_fd, "k.", ms=3, label="FTCS, 1200 s")
    ax.axhline(tile["T_limit"], color="r", lw=0.8, ls="--")
    ax.set_xlabel("depth from hot face [mm]"); ax.set_ylabel("T [K]"); ax.set_title("re-entry tile, back face insulated")
    ax.legend(fontsize=7)
    ax = axes[1]
    for t_ms in (0.0, 1.0, 2.0, 3.0):
        ax.plot(xs, wave_series(xs, t_ms * 1e-3, Lc, c, b_tri) * 1e3, label=f"t = {t_ms:.0f} ms")
    ax.plot(xs, dalembert(xs, 2e-3, c, Lc, tri) * 1e3, "k.", ms=2, label="d'Alembert, 2 ms")
    ax.set_xlabel("x [m]"); ax.set_ylabel("u [mm]"); ax.set_title("plucked cable (60 modes)"); ax.legend(fontsize=7)
    ax = axes[2]
    cs = ax.contourf(X, Y, u_series, levels=12)
    ax.contour(X, Y, u_jac, levels=cs.levels, colors="k", linewidths=0.5)
    fig.colorbar(cs, ax=ax, label="T - T_edge [K]")
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.set_title("plate: series (fill), Jacobi (lines)")
    ax.set_aspect("equal")
    fig.suptitle("P02.11  Separation of variables against finite differences")
    fig.tight_layout()
    fig.savefig(out_dir / "pde_solutions.png", dpi=150)
    print(f"\nsaved {out_dir / 'pde_solutions.png'}")


if __name__ == "__main__":
    main()
