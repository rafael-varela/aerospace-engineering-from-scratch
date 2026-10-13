"""P01.02 Applications of Derivatives: reference solution.

Maximum-range glide from a parabolic drag polar: optimisation by critical points, the second-
derivative test, linear approximation (sensitivity), related rates (best-glide speed during the
descent), Newton's method (best speed in a headwind), and their numerical checks.

Model:
    C_D = C_D0 + K C_L^2,            K = 1 / (pi e AR)
    L = W cos(gamma),  D = W sin(gamma),  tan(gamma) = C_D / C_L       (steady glide, no thrust)
    V = sqrt(2 W cos(gamma) / (rho S C_L)),  sink w = V sin(gamma)
Best glide (max L/D) at C_L* = sqrt(C_D0 / K), (L/D)max = 1 / (2 sqrt(C_D0 K)); min sink at sqrt(3) C_L*.

Run:  python solution.py          (prints the glide tables, writes glide_polar.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-01/lesson-02-applications-of-derivatives
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np

G0 = 9.80665            # m/s^2
RHO0 = 1.225            # kg/m^3, ISA sea level
T0 = 288.15             # K
R_AIR = 287.05287       # J/(kg K)
LAPSE = -0.0065         # K/m, ISA troposphere

# The lesson's training glider.
GLIDER = {"mass": 350.0, "area": 10.5, "aspect_ratio": 20.0, "oswald": 0.9, "cd0": 0.015}


# ---------------------------------------------------------------- drag polar and glide performance

def induced_drag_factor(aspect_ratio: float, oswald: float) -> float:
    """Return K = 1 / (pi e AR), the induced-drag factor of the parabolic polar. Raise ``ValueError`` unless both > 0."""
    if aspect_ratio <= 0.0 or oswald <= 0.0:
        raise ValueError("aspect ratio and Oswald factor must be positive")
    return 1.0 / (np.pi * oswald * aspect_ratio)


def drag_polar(cl: np.ndarray | float, cd0: float, k: float) -> np.ndarray | float:
    """Return C_D = cd0 + k C_L^2 for lift coefficient ``cl`` (scalar or array)."""
    return cd0 + k * np.asarray(cl, dtype=float) ** 2


def lift_to_drag(cl: np.ndarray | float, cd0: float, k: float) -> np.ndarray | float:
    """Return L/D = C_L / C_D along the polar."""
    cl = np.asarray(cl, dtype=float)
    return cl / drag_polar(cl, cd0, k)


def ld_derivative(cl: np.ndarray | float, cd0: float, k: float) -> np.ndarray | float:
    """Return d(L/D)/dC_L = (cd0 - k C_L^2) / (cd0 + k C_L^2)^2 by the quotient rule."""
    cl = np.asarray(cl, dtype=float)
    return (cd0 - k * cl**2) / (cd0 + k * cl**2) ** 2


def ld_second_derivative(cl: np.ndarray | float, cd0: float, k: float) -> np.ndarray | float:
    """Return d^2(L/D)/dC_L^2 = -2 k C_L (3 cd0 - k C_L^2) / (cd0 + k C_L^2)^3."""
    cl = np.asarray(cl, dtype=float)
    return -2.0 * k * cl * (3.0 * cd0 - k * cl**2) / (cd0 + k * cl**2) ** 3


def best_glide_cl(cd0: float, k: float) -> float:
    """Return C_L* = sqrt(cd0 / k), the critical point of L/D (where induced drag equals parasite drag)."""
    return float(np.sqrt(cd0 / k))


def max_lift_to_drag(cd0: float, k: float) -> float:
    """Return (L/D)max = 1 / (2 sqrt(cd0 k))."""
    return float(1.0 / (2.0 * np.sqrt(cd0 * k)))


def min_sink_cl(cd0: float, k: float) -> float:
    """Return the minimum-sink lift coefficient sqrt(3 cd0 / k) (critical point of C_D / C_L^(3/2))."""
    return float(np.sqrt(3.0 * cd0 / k))


def glide_angle(cl: np.ndarray | float, cd0: float, k: float) -> np.ndarray | float:
    """Return the glide angle gamma = arctan(C_D / C_L) in radians (positive downwards), exact for a steady glide."""
    cl = np.asarray(cl, dtype=float)
    return np.arctan(drag_polar(cl, cd0, k) / cl)


def airspeed(cl: np.ndarray | float, weight: float, rho: float, area: float, gamma: np.ndarray | float = 0.0) -> np.ndarray | float:
    """Return V = sqrt(2 W cos(gamma) / (rho S C_L)) in m/s. Raise ``ValueError`` if any C_L <= 0 or rho, S, W <= 0."""
    cl = np.asarray(cl, dtype=float)
    if np.any(cl <= 0.0) or weight <= 0.0 or rho <= 0.0 or area <= 0.0:
        raise ValueError("C_L, weight, density and area must be positive")
    return np.sqrt(2.0 * weight * np.cos(gamma) / (rho * area * cl))


def sink_rate(cl: np.ndarray | float, cd0: float, k: float, weight: float, rho: float, area: float) -> np.ndarray | float:
    """Return the sink rate w = V sin(gamma) in m/s for a steady glide at lift coefficient ``cl`` (exact angles)."""
    gamma = glide_angle(cl, cd0, k)
    return airspeed(cl, weight, rho, area, gamma) * np.sin(gamma)


def sink_rate_of_speed(v: np.ndarray | float, cd0: float, k: float, weight: float, rho: float, area: float) -> np.ndarray | float:
    """Return the small-angle sink rate w(V) = a V^3 + b / V with a = rho S cd0 / (2 W), b = 2 W k / (rho S).

    Follows from w = V C_D / C_L with C_L = 2 W / (rho S V^2) (L = W). Raise ``ValueError`` if any V <= 0.
    """
    v = np.asarray(v, dtype=float)
    if np.any(v <= 0.0):
        raise ValueError("airspeed must be positive")
    a = rho * area * cd0 / (2.0 * weight)
    b = 2.0 * weight * k / (rho * area)
    return a * v**3 + b / v


def glide_range(height: float, cd0: float, k: float) -> float:
    """Return the still-air glide range h (L/D)max in m (small-angle: tan(gamma) ~ sin(gamma))."""
    if height < 0.0:
        raise ValueError("height must be non-negative")
    return height * max_lift_to_drag(cd0, k)


# ---------------------------------------------------------------- Newton's method

def newton(f: Callable[[float], float], dfdx: Callable[[float], float], x0: float,
           tol: float = 1e-12, maxiter: int = 50) -> tuple[float, int]:
    """Solve f(x) = 0 by Newton's method from ``x0``; return (root, iterations used).

    Stops when |x_{n+1} - x_n| <= tol * max(1, |x|). Raise ``ZeroDivisionError`` if f' vanishes and
    ``RuntimeError`` if ``maxiter`` is reached without convergence.
    """
    x = float(x0)
    for n in range(1, maxiter + 1):
        d = dfdx(x)
        if d == 0.0 or not np.isfinite(d):
            raise ZeroDivisionError("derivative vanished in Newton's method")
        step = f(x) / d
        x_new = x - step
        if not np.isfinite(x_new):
            raise RuntimeError("Newton's method diverged")
        if abs(x_new - x) <= tol * max(1.0, abs(x_new)):
            return x_new, n
        x = x_new
    raise RuntimeError(f"Newton's method did not converge in {maxiter} iterations")


def best_glide_cl_newton(cd0: float, k: float, cl0: float = 0.5) -> float:
    """Find the best-glide C_L as the root of ``ld_derivative`` with Newton's method, using ``ld_second_derivative``."""
    root, _ = newton(lambda c: float(ld_derivative(c, cd0, k)), lambda c: float(ld_second_derivative(c, cd0, k)), cl0)
    return root


# ---------------------------------------------------------------- linear approximation (sensitivity)

def ld_max_sensitivity_cd0(cd0: float, k: float) -> float:
    """Return d(L/D)max / d cd0 = -(L/D)max / (2 cd0), the sensitivity of best glide ratio to parasite drag."""
    return -max_lift_to_drag(cd0, k) / (2.0 * cd0)


def linear_approximation(f: Callable[[float], float], dfdx: Callable[[float], float], x0: float, dx: float) -> float:
    """Return the tangent-line estimate f(x0) + f'(x0) dx of f(x0 + dx)."""
    return float(f(x0) + dfdx(x0) * dx)


# ---------------------------------------------------------------- related rates: best-glide speed in the descent

def isa_density(h: np.ndarray | float) -> np.ndarray | float:
    """Return ISA troposphere density rho = rho0 (T/T0)^n, n = -g0/(R L) - 1, for geopotential h in [0, 11000] m."""
    h = np.asarray(h, dtype=float)
    if np.any(h < 0.0) or np.any(h > 11000.0):
        raise ValueError("h must lie in the troposphere [0, 11000] m")
    n = -G0 / (R_AIR * LAPSE) - 1.0
    out = RHO0 * ((T0 + LAPSE * h) / T0) ** n
    return float(out) if out.ndim == 0 else out


def isa_density_gradient(h: np.ndarray | float) -> np.ndarray | float:
    """Return d rho / dh = rho n L / T (kg/m^4) by the chain rule applied to ``isa_density``."""
    h = np.asarray(h, dtype=float)
    n = -G0 / (R_AIR * LAPSE) - 1.0
    out = isa_density(h) * n * LAPSE / (T0 + LAPSE * h)
    return float(out) if np.ndim(out) == 0 else out


def best_glide_speed(h: np.ndarray | float, cd0: float, k: float, weight: float, area: float) -> np.ndarray | float:
    """Return the best-glide airspeed V* = sqrt(2 W / (rho(h) S C_L*)) at altitude ``h`` (small-angle)."""
    return np.sqrt(2.0 * weight / (isa_density(h) * area * best_glide_cl(cd0, k)))


def best_glide_speed_rate(h: float, hdot: float, cd0: float, k: float, weight: float, area: float) -> float:
    """Return dV*/dt = -(V* / (2 rho)) (d rho/dh) (dh/dt): the related rate as the glider descends (hdot < 0)."""
    v = float(best_glide_speed(h, cd0, k, weight, area))
    return -(v / (2.0 * isa_density(h))) * float(isa_density_gradient(h)) * hdot


# ---------------------------------------------------------------- headwind: Newton on a transcendental optimum

def ground_glide_ratio(v: np.ndarray | float, v_wind: float, cd0: float, k: float, weight: float, rho: float, area: float) -> np.ndarray | float:
    """Return (V - V_wind) / w(V): ground distance per metre of height lost, headwind positive."""
    v = np.asarray(v, dtype=float)
    return (v - v_wind) / sink_rate_of_speed(v, cd0, k, weight, rho, area)


def headwind_best_speed(v_wind: float, cd0: float, k: float, weight: float, rho: float, area: float) -> float:
    """Return the airspeed maximising ``ground_glide_ratio`` by Newton's method on its derivative.

    d/dV [(V - Vw)/w] = 0  <=>  F(V) = w(V) - (V - Vw) w'(V) = 0, with F'(V) = -(V - Vw) w''(V).
    For Vw = 0 the root is the still-air best-glide speed V*; the iteration starts from 1.2 V* + max(Vw, 0),
    on the far side of the optimum where F is monotone, so it converges for any headwind below the
    speed at which the glider cannot make headway.
    """
    a = rho * area * cd0 / (2.0 * weight)
    b = 2.0 * weight * k / (rho * area)
    w = lambda v: a * v**3 + b / v
    w1 = lambda v: 3.0 * a * v**2 - b / v**2
    w2 = lambda v: 6.0 * a * v + 2.0 * b / v**3
    f = lambda v: w(v) - (v - v_wind) * w1(v)
    df = lambda v: -(v - v_wind) * w2(v)
    v_still = float((b / a) ** 0.25)
    root, _ = newton(f, df, 1.2 * v_still + max(v_wind, 0.0))
    return root


# ---------------------------------------------------------------- driver

def main(out_dir: str | Path = ".") -> None:
    """Print the glider's performance table and save the sink-rate polar with the best-glide tangents."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    g = GLIDER
    w, s, cd0 = g["mass"] * G0, g["area"], g["cd0"]
    k = induced_drag_factor(g["aspect_ratio"], g["oswald"])
    cl_bg, ld_max, cl_ms = best_glide_cl(cd0, k), max_lift_to_drag(cd0, k), min_sink_cl(cd0, k)
    print(f"W = {w:.1f} N, S = {s} m^2, C_D0 = {cd0}, K = {k:.5f}")
    print(f"best glide:  C_L* = {cl_bg:.4f}, (L/D)max = {ld_max:.2f}, Newton gives C_L = {best_glide_cl_newton(cd0, k):.6f}")
    print(f"             gamma = {np.degrees(glide_angle(cl_bg, cd0, k)):.3f} deg, V = {airspeed(cl_bg, w, RHO0, s, glide_angle(cl_bg, cd0, k)):.2f} m/s,"
          f" sink = {sink_rate(cl_bg, cd0, k, w, RHO0, s):.3f} m/s")
    print(f"min sink:    C_L = {cl_ms:.4f}, V = {airspeed(cl_ms, w, RHO0, s, glide_angle(cl_ms, cd0, k)):.2f} m/s,"
          f" sink = {sink_rate(cl_ms, cd0, k, w, RHO0, s):.3f} m/s, L/D = {lift_to_drag(cl_ms, cd0, k):.2f}")
    print(f"still-air range from 1000 m: {glide_range(1000.0, cd0, k) / 1e3:.2f} km")
    from scipy.optimize import minimize_scalar

    res = minimize_scalar(lambda c: -float(lift_to_drag(c, cd0, k)), bounds=(0.05, 3.0), method="bounded", options={"xatol": 1e-10})
    print(f"scipy minimize_scalar: C_L = {res.x:.8f} ({res.nfev} evaluations), differs from sqrt(C_D0/K) by {res.x - cl_bg:+.1e};"
          f" L/D there differs by {-res.fun - ld_max:+.1e}")

    d = ld_max_sensitivity_cd0(cd0, k)
    for frac in (0.01, 0.10, 0.30):
        dcd0 = frac * cd0
        lin = linear_approximation(lambda c: max_lift_to_drag(c, k), lambda c: ld_max_sensitivity_cd0(c, k), cd0, dcd0)
        exact = max_lift_to_drag(cd0 + dcd0, k)
        print(f"C_D0 + {100 * frac:.0f} %: linear (L/D)max = {lin:.3f}, exact = {exact:.3f}, linear error = {100 * (lin / exact - 1):+.2f} %"
              f"  (d(L/D)max/dC_D0 = {d:.1f} per unit C_D0)")

    h, hdot = 1000.0, -float(sink_rate(cl_bg, cd0, k, w, isa_density(1000.0), s))
    print(f"\nrelated rate at {h:.0f} m descending at {-hdot:.3f} m/s: V* = {best_glide_speed(h, cd0, k, w, s):.3f} m/s,"
          f" dV*/dt = {best_glide_speed_rate(h, hdot, cd0, k, w, s):+.3e} m/s^2")

    print(f"\n{'headwind [m/s]':>15} {'best V [m/s]':>13} {'ground L/D':>11}")
    for vw in (-10.0, -5.0, 0.0, 5.0, 10.0, 15.0):
        vb = headwind_best_speed(vw, cd0, k, w, RHO0, s)
        print(f"{vw:15.1f} {vb:13.2f} {float(ground_glide_ratio(vb, vw, cd0, k, w, RHO0, s)):11.2f}")

    v = np.linspace(14.0, 45.0, 300)
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(v, -sink_rate_of_speed(v, cd0, k, w, RHO0, s), label="sink polar w(V)")
    for vw, style in ((0.0, "k--"), (10.0, "r:")):
        vb = headwind_best_speed(vw, cd0, k, w, RHO0, s)
        wb = float(sink_rate_of_speed(vb, cd0, k, w, RHO0, s))
        ax.plot([vw, vb, 45.0], [0.0, -wb, -wb * (45.0 - vw) / (vb - vw)], style, lw=1, label=f"tangent from V_wind = {vw:.0f} m/s: V = {vb:.1f} m/s")
        ax.plot(vb, -wb, "o", color=style[0])
    ax.set_xlabel("airspeed V [m/s]")
    ax.set_ylabel("sink rate [m/s] (down is negative)")
    ax.set_ylim(-3.0, 0.2)
    ax.set_xlim(0.0, 45.0)
    ax.grid(True, lw=0.3)
    ax.legend(fontsize=8)
    ax.set_title("Glide polar: the tangent from the wind point gives best range")
    fig.tight_layout()
    fig.savefig(out_dir / "glide_polar.png", dpi=150)
    print(f"saved {out_dir / 'glide_polar.png'}")


if __name__ == "__main__":
    main()
