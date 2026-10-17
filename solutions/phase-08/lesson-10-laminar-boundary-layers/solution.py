"""P08.10 Boundary Layers: Laminar: reference solution.

How thick is the laminar boundary layer on the light aeroplane's wing, and how much friction drag
does it carry? This module answers with the Blasius similarity solution of Prandtl's boundary-layer
equations for a flat plate at zero pressure gradient, computed two independent ways:

* shooting: integrate f''' + (1/2) f f'' = 0 from the wall with classical RK4 for a guessed wall
  curvature f''(0), and adjust the guess with a secant root finder until f'(eta_max) = 1;
* collocation: hand the same two-point boundary-value problem to ``scipy.integrate.solve_bvp``;

and checked a third way by marching the boundary-layer partial differential equations
u u_x + v u_y = nu u_yy, u_x + v_y = 0 downstream from a uniform profile with an implicit finite-
difference scheme, which never assumes similarity and still collapses onto the Blasius profile.

Conventions: eta = y sqrt(U / (nu x)), u/U = f'(eta), v sqrt(x / (nu U)) = (eta f' - f)/2.
Published anchor (labelled published): f''(0) = 0.332057 (Howarth 1938 and later), giving
c_f sqrt(Re_x) = 0.664, theta sqrt(Re_x)/x = 0.664, delta* sqrt(Re_x)/x = 1.7208, H = 2.591,
delta99 sqrt(Re_x)/x = 4.91.

Units: SI (m, s, kg, Pa).

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes blasius.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-08/lesson-10-laminar-boundary-layers
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp
from scipy.linalg import solve_banded

# ISA constants (shared Phase 08 numbers)
T0, P0, R_AIR, G0, LAPSE = 288.15, 101325.0, 287.05287, 9.80665, -0.0065
FPP0_PUBLISHED = 0.332057          # published Blasius wall curvature f''(0)


# ---------------------------------------------------------------- given helpers

def isa_troposphere(h: float) -> tuple[float, float, float, float]:
    """ISA temperature (K), pressure (Pa), density (kg/m^3) and Sutherland viscosity (Pa s) at a
    geopotential altitude h (m) below 11 km (the lapse-rate layer of P08.01)."""
    T = T0 + LAPSE * h
    p = P0 * (T / T0) ** (-G0 / (R_AIR * LAPSE))
    rho = p / (R_AIR * T)
    mu = 1.458e-6 * T ** 1.5 / (T + 110.4)
    return T, p, rho, mu


def observed_order(errs, ratio: float = 2.0) -> np.ndarray:
    """Observed order of convergence log(e_k / e_{k+1}) / log(ratio) for a sequence refined by ``ratio``."""
    e = np.asarray(errs, dtype=float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


def _cumtrapz(f: np.ndarray, x: np.ndarray) -> np.ndarray:
    return np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(x))])


# ---------------------------------------------------------------- the Blasius problem

def blasius_rhs(eta: float, Y: np.ndarray) -> np.ndarray:
    """Right-hand side of the Blasius equation as a first-order system, Y = (f, f', f''):
    dY/deta = (f', f'', -f f''/2). Works for Y of shape (3,) or (3, m) (vectorised for solve_bvp)."""
    f, fp, fpp = Y
    return np.array([fp, fpp, -0.5 * f * fpp])


def rk4_blasius(fpp0: float, eta_max: float = 10.0, n: int = 1000) -> tuple[np.ndarray, np.ndarray]:
    """Integrate the Blasius system from the wall, Y(0) = (0, 0, fpp0), to eta_max with n classical
    fourth-order Runge-Kutta steps of equal size. Returns (eta, Y) with eta of shape (n+1,) and Y of
    shape (n+1, 3) holding f, f', f''."""
    h = eta_max / n
    eta = np.linspace(0.0, eta_max, n + 1)
    Y = np.empty((n + 1, 3))
    Y[0] = (0.0, 0.0, fpp0)
    y = Y[0].copy()
    for i in range(n):
        t = eta[i]
        k1 = blasius_rhs(t, y)
        k2 = blasius_rhs(t + h / 2, y + h / 2 * k1)
        k3 = blasius_rhs(t + h / 2, y + h / 2 * k2)
        k4 = blasius_rhs(t + h, y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        Y[i + 1] = y
    return eta, Y


def shoot_blasius(eta_max: float = 10.0, n: int = 1000, guesses: tuple[float, float] = (0.2, 0.5),
                  tol: float = 1e-12, maxiter: int = 50) -> tuple[float, np.ndarray, np.ndarray]:
    """Shooting method: find f''(0) such that the RK4 solution satisfies f'(eta_max) = 1, using the
    secant method on the residual r(s) = f'(eta_max; s) - 1 started from the two ``guesses``.
    Returns (fpp0, eta, Y) for the converged wall curvature. Raises RuntimeError if the secant
    iteration does not converge to |r| < tol."""
    def resid(s: float) -> float:
        return rk4_blasius(s, eta_max, n)[1][-1, 1] - 1.0

    s0, s1 = guesses
    r0, r1 = resid(s0), resid(s1)
    for _ in range(maxiter):
        if abs(r1) < tol:
            eta, Y = rk4_blasius(s1, eta_max, n)
            return s1, eta, Y
        s0, s1 = s1, s1 - r1 * (s1 - s0) / (r1 - r0)
        r0, r1 = r1, resid(s1)
    raise RuntimeError("secant shooting did not converge")


def blasius_bvp(eta_max: float = 10.0, tol: float = 1e-8) -> tuple[float, np.ndarray, np.ndarray]:
    """Solve the Blasius boundary-value problem f(0) = f'(0) = 0, f'(eta_max) = 1 by collocation with
    ``scipy.integrate.solve_bvp`` (an independent method). Initial guess f' = 1 - exp(-eta).
    Returns (fpp0, eta, Y) with Y of shape (m, 3) on solve_bvp's final mesh."""
    eta = np.linspace(0.0, eta_max, 101)
    guess = np.vstack([eta - 1 + np.exp(-eta), 1 - np.exp(-eta), np.exp(-eta)])
    bc = lambda ya, yb: np.array([ya[0], ya[1], yb[1] - 1.0])
    sol = solve_bvp(blasius_rhs, bc, eta, guess, tol=tol, max_nodes=100000)
    if not sol.success:
        raise RuntimeError(sol.message)
    return float(sol.y[2, 0]), sol.x, sol.y.T


def blasius_integrals(eta: np.ndarray, Y: np.ndarray) -> dict[str, float]:
    """Boundary-layer thicknesses in similarity units (multiply by x / sqrt(Re_x) for metres):
    delta99 (eta where f' = 0.99, by linear interpolation), delta_star = int (1 - f') deta,
    theta = int f' (1 - f') deta (trapezoid rule on the given points), the shape factor
    H = delta_star / theta, the wall curvature fpp0 = f''(0), and v_edge = (eta f' - f)/2 at
    eta[-1], the outward velocity v sqrt(x / (nu U)) at the edge of the domain."""
    f, fp, fpp = Y[:, 0], Y[:, 1], Y[:, 2]
    i = int(np.argmax(fp >= 0.99))
    d99 = eta[i - 1] + (0.99 - fp[i - 1]) * (eta[i] - eta[i - 1]) / (fp[i] - fp[i - 1])
    ds = float(np.trapezoid(1 - fp, eta))
    th = float(np.trapezoid(fp * (1 - fp), eta))
    return {"delta99": float(d99), "delta_star": ds, "theta": th, "H": ds / th,
            "fpp0": float(fpp[0]), "v_edge": float(0.5 * (eta[-1] * fp[-1] - f[-1]))}


def plate_boundary_layer(x: np.ndarray, U: float, nu: float, rho: float,
                         consts: dict[str, float]) -> dict[str, np.ndarray]:
    """Dimensional laminar boundary layer on a flat plate at stations x (m) from the leading edge,
    from the similarity constants ``consts`` (the output of ``blasius_integrals``):
    Re_x = U x / nu, delta99 = consts['delta99'] x / sqrt(Re_x) (and likewise delta_star and theta),
    the local skin-friction coefficient c_f = tau_w / (rho U^2 / 2) = 2 f''(0) / sqrt(Re_x), and the
    wall shear stress tau_w (Pa). Returns a dict of arrays with keys
    Re_x, delta99, delta_star, theta, cf, tau_w."""
    x = np.asarray(x, dtype=float)
    Re = U * x / nu
    s = x / np.sqrt(Re)
    cf = 2.0 * consts["fpp0"] / np.sqrt(Re)
    return {"Re_x": Re, "delta99": consts["delta99"] * s, "delta_star": consts["delta_star"] * s,
            "theta": consts["theta"] * s, "cf": cf, "tau_w": 0.5 * rho * U ** 2 * cf}


def friction_drag(U: float, L: float, b: float, rho: float, nu: float, fpp0: float,
                  sides: int = 2) -> tuple[float, float]:
    """Laminar friction drag of a flat plate of chord L and span b wetted on ``sides`` sides:
    the integrated coefficient C_F = (1/L) int_0^L c_f dx = 4 f''(0) / sqrt(Re_L) (twice the local
    value at x = L, because c_f ~ x^(-1/2)), and D = (rho U^2 / 2) C_F (sides b L).
    Returns (D in N, C_F)."""
    CF = 4.0 * fpp0 / np.sqrt(U * L / nu)
    return 0.5 * rho * U ** 2 * CF * sides * b * L, float(CF)


def march_boundary_layer(U: float, nu: float, L: float, y_max: float, ny: int = 200, nx: int = 400,
                         x0_frac: float = 1e-4) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """March Prandtl's boundary-layer equations for a flat plate at zero pressure gradient,

        u u_x + v u_y = nu u_yy,   u_x + v_y = 0,   u(x, 0) = v(x, 0) = 0,   u(x, y_max) = U,

    downstream from a uniform profile u = U (u = 0 at the wall) at x0 = x0_frac L to x = L, with no
    similarity assumption. Grid: ny equal intervals in y on [0, y_max]; nx steps in x spaced
    geometrically from x0 to L (the layer grows like sqrt(x), so steps grow with x). Each step is
    implicit (backward Euler in x, central differences in y) with the convecting velocities u, v
    lagged from the previous station, so it is a tridiagonal solve; v then follows from continuity
    by the trapezoid rule. Returns (x, y, u, v) with u, v of shape (nx+1, ny+1)."""
    x = np.geomspace(x0_frac * L, L, nx + 1)
    y = np.linspace(0.0, y_max, ny + 1)
    dy = y[1] - y[0]
    u = np.empty((nx + 1, ny + 1))
    v = np.zeros((nx + 1, ny + 1))
    u[0] = U
    u[0, 0] = 0.0
    m = ny - 1                                   # interior unknowns i = 1..ny-1
    for j in range(nx):
        dx = x[j + 1] - x[j]
        uj, vj = u[j, 1:-1], v[j, 1:-1]
        lower = -vj / (2 * dy) - nu / dy ** 2   # coefficient of u_{i-1}
        diag = uj / dx + 2 * nu / dy ** 2
        upper = vj / (2 * dy) - nu / dy ** 2    # coefficient of u_{i+1}
        rhs = uj * uj / dx
        rhs[-1] -= upper[-1] * U                 # u = U at y_max (u = 0 at the wall adds nothing)
        ab = np.zeros((3, m))
        ab[0, 1:] = upper[:-1]
        ab[1] = diag
        ab[2, :-1] = lower[1:]
        u[j + 1, 1:-1] = solve_banded((1, 1), ab, rhs)
        u[j + 1, 0], u[j + 1, -1] = 0.0, U
        dudx = (u[j + 1] - u[j]) / dx
        v[j + 1] = -_cumtrapz(dudx, y)
    return x, y, u, v


def wall_cf_from_profile(y: np.ndarray, u: np.ndarray, U: float, nu: float) -> float:
    """Skin-friction coefficient from a velocity profile: c_f = 2 nu (du/dy)_wall / U^2, with the
    wall gradient from the second-order one-sided difference (-3u0 + 4u1 - u2)/(2 dy)."""
    dy = y[1] - y[0]
    dudy = (-3 * u[0] + 4 * u[1] - u[2]) / (2 * dy)
    return float(2 * nu * dudy / U ** 2)


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save blasius.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    fpp_s, eta, Y = shoot_blasius()
    fpp_b, eta_b, Y_b = blasius_bvp()
    c = blasius_integrals(eta, Y)
    cb = blasius_integrals(eta_b, Y_b)
    print(f"Blasius f''(0): shooting (RK4, n = 1000, eta_max = 10) {fpp_s:.8f}; solve_bvp {fpp_b:.8f}; "
          f"published {FPP0_PUBLISHED}")
    print(f"  difference shooting - bvp {fpp_s - fpp_b:.1e}")
    print(f"  delta99 {c['delta99']:.4f}, delta* {c['delta_star']:.5f}, theta {c['theta']:.5f}, "
          f"H {c['H']:.4f}, v_edge {c['v_edge']:.4f} (bvp: theta {cb['theta']:.5f}, H {cb['H']:.4f})")

    print("RK4 convergence of f''(0) (eta_max = 10, reference n = 4000):")
    ref = shoot_blasius(n=4000)[0]
    ns = [25, 50, 100, 200]
    errs = [abs(shoot_blasius(n=k)[0] - ref) for k in ns]
    for k, e in zip(ns, errs):
        print(f"  n = {k:4d}: |error| = {e:.3e}")
    print(f"  observed orders {np.array2string(observed_order(errs), precision=2)}")

    print("domain truncation: f''(0) for eta_max = 4, 5, 6, 8, 10, 12 (n = 100 per unit eta):")
    for em in (4.0, 5.0, 6.0, 8.0, 10.0, 12.0):
        print(f"  eta_max {em:4.1f}: {shoot_blasius(em, int(100 * em))[0]:.8f}")

    # the light aeroplane's wing
    h, V, chord, semi = 3000.0, 60.0, 1.5, 5.0
    T, p, rho, mu = isa_troposphere(h)
    nu = mu / rho
    print(f"\nwing at {h:.0f} m: T {T:.2f} K, p {p:.0f} Pa, rho {rho:.4f} kg/m^3, mu {mu:.4e} Pa s, nu {nu:.4e} m^2/s")
    xs = np.array([0.01, 0.1, 0.3, 0.75, 1.5])
    bl = plate_boundary_layer(xs, V, nu, rho, c)
    print("     x (m)     Re_x      delta99 (mm)  delta* (mm)  theta (mm)   c_f        tau_w (Pa)")
    for i, xx in enumerate(xs):
        print(f"  {xx:7.3f}  {bl['Re_x'][i]:10.4g}   {bl['delta99'][i] * 1e3:8.3f}     {bl['delta_star'][i] * 1e3:8.3f}"
              f"    {bl['theta'][i] * 1e3:8.3f}   {bl['cf'][i]:.4e}  {bl['tau_w'][i]:.3f}")
    D, CF = friction_drag(V, chord, 2 * semi, rho, nu, c["fpp0"])
    q = 0.5 * rho * V ** 2
    print(f"all-laminar wing (two sides, {2 * semi * chord:.0f} m^2 planform): Re_c {V * chord / nu:.4g}, "
          f"C_F {CF:.4e}, q {q:.1f} Pa, friction drag {D:.1f} N")

    # independent check: marching the PDEs (on a 1 m plate at 10 m/s in sea-level air, Re_L = 6.8e5)
    Um, num, Lm = 10.0, 1.4607e-5, 1.0
    d99L = 5.0 * Lm / np.sqrt(Um * Lm / num)
    xm, ym, um, vm = march_boundary_layer(Um, num, Lm, 3 * d99L, ny=300, nx=1600)
    print("\nmarching the boundary-layer PDEs from a uniform profile (U = 10 m/s, nu = 1.4607e-5, L = 1 m, ny = 300, nx = 1600):")
    for xf in (0.01, 0.1, 0.5, 1.0):
        j = int(np.argmin(abs(xm - xf)))
        cfm = wall_cf_from_profile(ym, um[j], Um, num)
        print(f"  x = {xm[j]:.3f} m: c_f sqrt(Re_x) = {cfm * np.sqrt(Um * xm[j] / num):.4f} (Blasius {2 * c['fpp0']:.4f})")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].plot(Y[:, 1], eta, label="u/U = f′ (shooting, RK4)")
    ax[0].plot(Y_b[:, 1], eta_b, "--", label="f′ (solve_bvp)")
    ax[0].plot(Y[:, 2] / Y[0, 2], eta, label="f″/f″(0)")
    for xf, mk in ((0.1, "o"), (0.5, "s"), (1.0, "^")):
        j = int(np.argmin(abs(xm - xf)))
        etam = ym * np.sqrt(Um / (num * xm[j]))
        ax[0].plot(um[j, ::6] / Um, etam[::6], mk, ms=3, label=f"marched PDE, x = {xm[j]:.1f} m")
    ax[0].set_ylim(0, 8)
    ax[0].set_xlabel("u/U")
    ax[0].set_ylabel("η = y √(U/(νx))")
    ax[0].set_title("Blasius profile and the marched PDE")
    ax[0].legend(fontsize=7)
    ax[1].loglog(10.0 / np.array(ns), errs, "o-", label="|f″(0) − reference|")
    ax[1].loglog(10.0 / np.array(ns), errs[0] * (ns[0] / np.array(ns)) ** 4, "k--", lw=0.8, label="slope 4")
    ax[1].set_xlabel("RK4 step Δη")
    ax[1].set_ylabel("error in f″(0)")
    ax[1].set_title("RK4 shooting convergence")
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(out_dir / "blasius.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'blasius.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
