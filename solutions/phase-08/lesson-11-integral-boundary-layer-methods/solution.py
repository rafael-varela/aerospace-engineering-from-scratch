"""P08.11 Integral Boundary Layer Methods: reference solution.

How much deceleration can a laminar boundary layer take before it separates? This module answers
with Thwaites's integral method, which turns von Karman's momentum integral equation

    d(theta)/dx + (2 + H) (theta / Ue) dUe/dx = c_f / 2

into one quadrature, theta^2 Ue^6 = theta_0^2 Ue_0^6 + 0.45 nu int Ue^5 dx, with the shape factor H and
the wall-shear parameter l as functions of the pressure-gradient parameter lambda = theta^2 Ue'/nu
(Cebeci and Bradshaw's fits to exact similarity solutions; laminar separation at lambda = -0.09).

Independent checks built here:

* the Falkner-Skan similarity solutions f''' + f f'' + beta (1 - f'^2) = 0 (wedge flows Ue = C x^m,
  beta = 2m/(m+1)) by ``scipy.integrate.solve_bvp``, which include Blasius (beta = 0) and the
  Hiemenz stagnation-point flow (beta = 1), and give exact (lambda, l, H) to compare with Thwaites;
* a marching finite-difference solution of the boundary-layer PDEs with a pressure gradient,
  u u_x + v u_y = Ue Ue' + nu u_yy, which never uses an integral closure and finds separation where
  the wall shear falls to zero (Howarth's linearly decelerating flow: published x_s/L = 0.1199).

Blasius constants (published; derived in P08.10): f''(0) = 0.332057, theta sqrt(Re_x)/x = 0.664115,
delta* sqrt(Re_x)/x = 1.720788, H = 2.591.

Units: SI (m, s); the methods themselves are dimensionless where stated.

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes thwaites.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-08/lesson-11-integral-boundary-layer-methods
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp
from scipy.linalg import solve_banded

BLASIUS = {"fpp0": 0.332057, "theta": 0.664115, "delta_star": 1.720788, "H": 2.591}   # published
HOWARTH_XSEP = 0.1199          # published separation point x_s/L of Ue = U0 (1 - x/L) (exact numerical solutions)
LAMBDA_SEP = -0.09             # Thwaites's separation criterion (l(lambda) = 0)

# ISA constants (shared Phase 08 numbers)
T0, P0, R_AIR, G0, LAPSE = 288.15, 101325.0, 287.05287, 9.80665, -0.0065


# ---------------------------------------------------------------- given helpers

def isa_nu(h: float) -> float:
    """Kinematic viscosity (m^2/s) of the ISA troposphere at geopotential altitude h (m), with
    Sutherland's law (P08.01)."""
    T = T0 + LAPSE * h
    p = P0 * (T / T0) ** (-G0 / (R_AIR * LAPSE))
    rho = p / (R_AIR * T)
    mu = 1.458e-6 * T ** 1.5 / (T + 110.4)
    return mu / rho


def cumtrapz(f: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Cumulative trapezoid integral of samples f(x), starting at 0."""
    return np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(x))])


def observed_order(errs, ratio: float = 2.0) -> np.ndarray:
    """Observed order of convergence log(e_k / e_{k+1}) / log(ratio)."""
    e = np.asarray(errs, dtype=float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


def cylinder_edge_velocity(x: np.ndarray, U: float, R: float) -> tuple[np.ndarray, np.ndarray]:
    """Potential-flow edge velocity round a circular cylinder of radius R in a stream U, with x the
    arc length from the front stagnation point: Ue = 2U sin(x/R), dUe/dx = (2U/R) cos(x/R)."""
    x = np.asarray(x, dtype=float)
    return 2 * U * np.sin(x / R), 2 * U / R * np.cos(x / R)


def model_wing_edge_velocity(x: np.ndarray, U: float, c: float, x_peak: float, peak: float = 1.25,
                             te: float = 0.95, m: float = 0.1) -> tuple[np.ndarray, np.ndarray]:
    """A *model* upper-surface edge velocity for a wing section (not a real aerofoil's): a gently
    accelerating wedge flow Ue = peak U (x/x_peak)^m from the leading edge to the suction peak at
    x_peak, then a linear deceleration to te U at the trailing edge x = c. Returns (Ue, dUe/dx);
    dUe/dx is infinite at x = 0 (the wedge-flow singularity) and returned as inf there."""
    x = np.asarray(x, dtype=float)
    up = x <= x_peak
    with np.errstate(divide="ignore", invalid="ignore"):
        Ue = np.where(up, peak * U * (x / x_peak) ** m, U * (peak - (peak - te) * (x - x_peak) / (c - x_peak)))
        dU = np.where(up, np.where(x > 0, m * peak * U * (x / x_peak) ** m / x, np.inf),
                      -U * (peak - te) / (c - x_peak))
    return Ue, dU


# ---------------------------------------------------------------- Thwaites's method

def thwaites_l_H(lam: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Thwaites's closure functions, in Cebeci and Bradshaw's curve fits to exact solutions:

        0 <= lambda <= 0.1:     l = 0.22 + 1.57 lambda - 1.8 lambda^2,   H = 2.61 - 3.75 lambda + 5.24 lambda^2
        -0.1 <= lambda < 0:     l = 0.22 + 1.402 lambda + 0.018 lambda / (lambda + 0.107),
                                H = 2.088 + 0.0731 / (lambda + 0.14)

    with l = tau_w theta / (mu Ue) and H = delta*/theta. Values of lambda outside [-0.1, 0.1] are
    clipped to that range (the fits are not valid beyond it). Returns (l, H) as arrays."""
    lam = np.clip(np.asarray(lam, dtype=float), -0.1, 0.1)
    pos = lam >= 0
    l = np.where(pos, 0.22 + 1.57 * lam - 1.8 * lam ** 2,
                 0.22 + 1.402 * lam + 0.018 * lam / (lam + 0.107))
    H = np.where(pos, 2.61 - 3.75 * lam + 5.24 * lam ** 2, 2.088 + 0.0731 / (lam + 0.14))
    return l, H


def thwaites(x: np.ndarray, Ue: np.ndarray, dUedx: np.ndarray, nu: float,
             theta0: float = 0.0) -> dict[str, np.ndarray]:
    """Thwaites's method on stations x (m, increasing) with edge velocity Ue (m/s) and its gradient
    dUe/dx (1/s) at those stations:

        theta^2 = (theta0^2 Ue0^6 + 0.45 nu int_x0^x Ue^5 dx) / Ue^6     (trapezoid rule),

    except at a stagnation point (Ue[0] = 0), where the limit theta^2 = 0.075 nu / Ue'(0) is used.
    Then lambda = theta^2 Ue'/nu, (l, H) = thwaites_l_H(lambda), delta* = H theta and
    c_f = 2 nu l / (Ue theta) (infinite at a stagnation point, where Ue = 0: returned as inf).
    Returns a dict with keys theta, lam, l, H, delta_star, cf (arrays like x); lambda is nan at a
    station where theta = 0 and dUe/dx is infinite (the leading edge of a wedge flow with m < 1)."""
    x, Ue, dU = (np.asarray(a, dtype=float) for a in (x, Ue, dUedx))
    I = cumtrapz(Ue ** 5, x)
    with np.errstate(divide="ignore", invalid="ignore"):
        th2 = (theta0 ** 2 * Ue[0] ** 6 + 0.45 * nu * I) / Ue ** 6
    if Ue[0] == 0.0:
        th2[0] = 0.075 * nu / dU[0]
    theta = np.sqrt(th2)
    with np.errstate(invalid="ignore"):
        lam = th2 * dU / nu                      # nan where theta = 0 and Ue' = inf (a wedge's leading edge)
    l, H = thwaites_l_H(lam)
    with np.errstate(divide="ignore"):
        cf = np.where(Ue > 0, 2 * nu * l / (Ue * theta), np.inf)
    return {"theta": theta, "lam": lam, "l": l, "H": H, "delta_star": H * theta, "cf": cf}


def separation_point(x: np.ndarray, lam: np.ndarray, lam_sep: float = LAMBDA_SEP) -> float:
    """First station where lambda falls through lam_sep, by linear interpolation between the two
    bracketing stations; nan if it never does."""
    x, lam = np.asarray(x, dtype=float), np.asarray(lam, dtype=float)
    below = np.nonzero(lam <= lam_sep)[0]
    if below.size == 0:
        return float("nan")
    i = below[0]
    if i == 0:
        return float(x[0])
    return float(x[i - 1] + (lam_sep - lam[i - 1]) * (x[i] - x[i - 1]) / (lam[i] - lam[i - 1]))


def wedge_thwaites(m: float) -> dict[str, float]:
    """Thwaites's method in closed form for a wedge flow Ue = C x^m (m > -1/5): the quadrature gives
    theta^2 = 0.45 nu x / ((5m + 1) Ue), so lambda = 0.45 m / (5m + 1) is constant along the wall.
    Returns lambda, l, H and the similarity-variable momentum thickness theta_eta, defined with
    eta = y sqrt((m + 1) Ue / (2 nu x)) so it can be compared with ``falkner_skan``:
    theta_eta = sqrt(0.45 (m + 1) / (2 (5m + 1)))."""
    lam = 0.45 * m / (5 * m + 1)
    l, H = thwaites_l_H(lam)
    return {"lam": lam, "l": float(l), "H": float(H), "theta_eta": float(np.sqrt(0.45 * (m + 1) / (2 * (5 * m + 1))))}


# ---------------------------------------------------------------- exact similarity solutions

def falkner_skan(beta: float, eta_max: float = 10.0, tol: float = 1e-8) -> dict[str, float]:
    """Falkner-Skan similarity solution f''' + f f'' + beta (1 - f'^2) = 0, f(0) = f'(0) = 0,
    f'(eta_max) = 1, by ``scipy.integrate.solve_bvp`` (eta = y sqrt((m + 1) Ue / (2 nu x)),
    beta = 2m/(m + 1)). beta = 0 is Blasius in this scaling (f''(0) = 0.4696 = 0.332057 sqrt 2),
    beta = 1 is Hiemenz's stagnation-point flow. Initial guess f' = 1 - exp(-eta); for beta < 0 a
    continuation in beta from 0 is used so the attached (upper-branch) solution is found.

    Returns a dict: fpp0, delta_star and theta (in eta units), H, and Thwaites's variables
    m = beta / (2 - beta), lam = theta^2 2m/(m + 1) and l = fpp0 theta (so that they compare
    directly with ``thwaites_l_H``)."""
    eta = np.linspace(0.0, eta_max, 201)
    y = np.vstack([eta - 1 + np.exp(-eta), 1 - np.exp(-eta), np.exp(-eta)])
    path = [beta] if beta >= 0 else list(np.linspace(0.0, beta, 6))
    for b in path:
        rhs = lambda t, Y, b=b: np.vstack([Y[1], Y[2], -Y[0] * Y[2] - b * (1 - Y[1] ** 2)])
        bc = lambda ya, yb: np.array([ya[0], ya[1], yb[1] - 1.0])
        sol = solve_bvp(rhs, bc, eta, y, tol=tol, max_nodes=200000)
        if not sol.success:
            raise RuntimeError(f"Falkner-Skan beta = {b}: {sol.message}")
        eta, y = sol.x, sol.y
    fp = y[1]
    ds = float(np.trapezoid(1 - fp, eta))
    th = float(np.trapezoid(fp * (1 - fp), eta))
    m = beta / (2 - beta)
    return {"fpp0": float(y[2, 0]), "delta_star": ds, "theta": th, "H": ds / th,
            "m": m, "lam": th ** 2 * 2 * m / (m + 1), "l": float(y[2, 0]) * th}


# ---------------------------------------------------------------- independent check: marching the PDEs

def march_boundary_layer(x: np.ndarray, Ue: np.ndarray, dUedx: np.ndarray, nu: float, y_max: float,
                         ny: int = 200) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """March Prandtl's boundary-layer equations with a pressure gradient,

        u u_x + v u_y = Ue dUe/dx + nu u_yy,   u_x + v_y = 0,   u = v = 0 at y = 0,   u = Ue at y_max,

    over the stations x (increasing; x[0] > 0), starting from u = Ue[0] off the wall. Each step is
    implicit in u (backward Euler in x, central differences on ny equal intervals in y) with u, v in
    the convective terms lagged; v follows from continuity by the trapezoid rule. Marching stops at
    the first station where the wall shear is no longer positive (separation: the equations cannot
    be marched into reversed flow) or where the solve fails.

    Returns (x_done, y, tau) with tau = nu (du/dy)_wall (m^2/s^2, i.e. wall shear over density) at
    each completed station, from the second-order one-sided difference."""
    x, Ue, dU = (np.asarray(a, dtype=float) for a in (x, Ue, dUedx))
    y = np.linspace(0.0, y_max, ny + 1)
    dy = y[1] - y[0]
    u = np.full(ny + 1, Ue[0])
    u[0] = 0.0
    v = np.zeros(ny + 1)
    taus = [nu * (-3 * u[0] + 4 * u[1] - u[2]) / (2 * dy)]
    m = ny - 1
    for j in range(len(x) - 1):
        dx = x[j + 1] - x[j]
        ui, vi = u[1:-1], v[1:-1]
        lower = -vi / (2 * dy) - nu / dy ** 2
        diag = ui / dx + 2 * nu / dy ** 2
        upper = vi / (2 * dy) - nu / dy ** 2
        rhs = ui * ui / dx + Ue[j + 1] * dU[j + 1]
        rhs[-1] -= upper[-1] * Ue[j + 1]
        ab = np.zeros((3, m))
        ab[0, 1:] = upper[:-1]
        ab[1] = diag
        ab[2, :-1] = lower[1:]
        unew = np.empty_like(u)
        unew[1:-1] = solve_banded((1, 1), ab, rhs)
        unew[0], unew[-1] = 0.0, Ue[j + 1]
        tau = nu * (-3 * unew[0] + 4 * unew[1] - unew[2]) / (2 * dy)
        if not np.all(np.isfinite(unew)) or tau <= 0.0:
            break
        v = -cumtrapz((unew - u) / dx, y)
        u = unew
        taus.append(tau)
    n = len(taus)
    return x[:n], y, np.array(taus)


def march_separation(x: np.ndarray, tau: np.ndarray) -> float:
    """Separation point from a marched wall-shear history: near separation tau^2 falls linearly in x
    (Goldstein's singularity, tau ~ sqrt(x_s - x)), so fit a straight line to tau^2 over the last 10
    completed stations and return where it reaches zero."""
    p = np.polyfit(x[-10:], tau[-10:] ** 2, 1)
    return float(-p[1] / p[0])


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save thwaites.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    nu = isa_nu(3000.0)
    U, c = 60.0, 1.5

    # flat plate
    x = np.linspace(0.0, c, 3001)
    fp = thwaites(x, np.full_like(x, U), np.zeros_like(x), nu)
    Re = U * x[-1] / nu
    print(f"flat plate (U = {U} m/s, nu = {nu:.4e} m^2/s, x = {c} m, Re_x = {Re:.4g}):")
    print(f"  theta sqrt(Re_x)/x: Thwaites {fp['theta'][-1] * np.sqrt(Re) / c:.4f}, Blasius {BLASIUS['theta']:.4f} "
          f"({(fp['theta'][-1] * np.sqrt(Re) / c / BLASIUS['theta'] - 1) * 100:+.2f} %)")
    print(f"  c_f sqrt(Re_x): Thwaites {fp['cf'][-1] * np.sqrt(Re):.4f}, Blasius {2 * BLASIUS['fpp0']:.4f} "
          f"({(fp['cf'][-1] * np.sqrt(Re) / (2 * BLASIUS['fpp0']) - 1) * 100:+.2f} %)")
    print(f"  H: Thwaites {fp['H'][-1]:.3f}, Blasius {BLASIUS['H']:.3f}")

    # Falkner-Skan against Thwaites
    print("\nwedge flows Ue = C x^m: exact (Falkner-Skan, solve_bvp) against Thwaites")
    print("     m       beta     lam exact  lam Thw    l exact   l Thw     H exact   H Thw    theta_eta err")
    for mm in (-0.0904, -0.06, -0.03, 0.0, 0.1, 0.333333, 1.0):
        beta = 2 * mm / (mm + 1)
        fs = falkner_skan(beta)
        tw = wedge_thwaites(mm)
        print(f"  {mm:7.4f}  {beta:7.4f}   {fs['lam']:8.4f}  {tw['lam']:8.4f}   {fs['l']:7.4f}  {tw['l']:7.4f}   "
              f"{fs['H']:7.3f}  {tw['H']:7.3f}   {(tw['theta_eta'] / fs['theta'] - 1) * 100:+6.2f} %")
    hz = falkner_skan(1.0)
    print(f"Hiemenz stagnation point: f''(0) = {hz['fpp0']:.5f}, theta sqrt(a/nu) = {hz['theta']:.4f} exact, "
          f"{np.sqrt(0.075):.4f} Thwaites; delta* sqrt(a/nu) = {hz['delta_star']:.4f}, H = {hz['H']:.3f}")
    sep = falkner_skan(-0.19883)
    print(f"Falkner-Skan at beta = -0.19883 (m = {sep['m']:.4f}): f''(0) = {sep['fpp0']:.5f}, H = {sep['H']:.3f}, "
          f"lambda = {sep['lam']:.4f}")

    # Howarth's linearly decelerating flow
    L = 1.0
    xh = np.linspace(0.0, 0.2 * L, 4001)
    th = thwaites(xh, U * (1 - xh / L), np.full_like(xh, -U / L), nu)
    xs_t = separation_point(xh, th["lam"])
    xm = np.linspace(1e-4 * L, 0.2 * L, 4001)
    d99 = 5 * np.sqrt(nu * 0.12 * L / U)
    xd, yd, tau = march_boundary_layer(xm, U * (1 - xm / L), np.full_like(xm, -U / L), nu, 4 * d99, ny=400)
    xs_m = march_separation(xd, tau)
    print(f"\nHowarth flow Ue = U (1 - x/L), L = {L} m: separation x_s/L Thwaites {xs_t / L:.4f}, "
          f"marched PDE {xs_m / L:.4f}, published {HOWARTH_XSEP}")
    print(f"  Thwaites: the layer separates after Ue has fallen by {xs_t / L * 100:.1f} %")

    # circular cylinder
    R = 0.05
    xc = np.linspace(0.0, 0.75 * np.pi * R, 6001)
    Uc, dUc = cylinder_edge_velocity(xc, 20.0, R)
    tc = thwaites(xc, Uc, dUc, 1.5e-5)
    xs_c = separation_point(xc, tc["lam"])
    xcm = np.linspace(1e-3 * R, 0.7 * np.pi * R, 8001)
    Ucm, dUcm = cylinder_edge_velocity(xcm, 20.0, R)
    xdc, _, tauc = march_boundary_layer(xcm, Ucm, dUcm, 1.5e-5, 12 * np.sqrt(1.5e-5 * R / 40.0), ny=300)
    xs_cm = march_separation(xdc, tauc)
    print(f"cylinder (potential-flow Ue): separation at {np.degrees(xs_c / R):.1f} deg (Thwaites), "
          f"{np.degrees(xs_cm / R):.1f} deg (marched PDE), from the front stagnation point")

    # the model wing: where does the laminar layer separate behind the suction peak?
    print(f"\nmodel wing upper surface (U = {U} m/s at 3000 m, chord {c} m, peak 1.25 U, trailing edge 0.95 U):")
    xw = np.linspace(0.0, c, 6001)
    for xp in (0.1, 0.3, 0.5, 0.7):
        Uw, dUw = model_wing_edge_velocity(xw, U, c, xp * c)
        tw = thwaites(xw, Uw, dUw, nu)
        xs = separation_point(xw, tw["lam"])
        Us = np.interp(xs, xw, Uw)
        ip = int(np.argmin(abs(xw - xp * c)))
        print(f"  peak at {xp:.1f} c: theta at peak {tw['theta'][ip] * 1e3:.3f} mm; laminar separation at x/c = {xs / c:.3f}, "
              f"{max((xs - xp * c) * 1e3, 0.0):.0f} mm behind the peak, after Ue fell {(1 - Us / (1.25 * U)) * 100:.1f} %")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    lam = np.linspace(-0.09, 0.1, 200)
    l, H = thwaites_l_H(lam)
    ax[0].plot(lam, l, label="l(λ) fit")
    ax[0].plot(lam, H / 10, label="H(λ)/10 fit")
    for bb in (-0.19, -0.14, -0.1, -0.05, 0.0, 0.3, 0.6, 1.0):
        fs = falkner_skan(bb)
        ax[0].plot(fs["lam"], fs["l"], "ko", ms=3)
        ax[0].plot(fs["lam"], fs["H"] / 10, "ks", ms=3)
    ax[0].plot([], [], "ko", ms=3, label="Falkner–Skan exact")
    ax[0].axvline(LAMBDA_SEP, color="r", lw=0.8, ls="--", label="separation λ = −0.09")
    ax[0].set_xlabel("λ = θ² Ue′/ν")
    ax[0].set_title("Thwaites closure against exact wedge flows")
    ax[0].legend(fontsize=8)
    ax[1].plot(xh / L, th["lam"], label="Howarth flow, Thwaites")
    ax[1].plot(np.degrees(xc / R) / 180, tc["lam"], label="cylinder, Thwaites (x = angle/180°)")
    ax[1].axhline(LAMBDA_SEP, color="r", lw=0.8, ls="--")
    ax[1].set_ylim(-0.12, 0.1)
    ax[1].set_xlabel("x/L   or   angle/180°")
    ax[1].set_ylabel("λ")
    ax[1].set_title("pressure-gradient parameter to separation")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "thwaites.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'thwaites.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
