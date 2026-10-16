"""P02.13 Calculus of Variations: starter code.

The brachistochrone (fastest frictionless gravity chute between two points) solved three ways: the
Euler-Lagrange equation and Beltrami identity in SymPy, the closed-form cycloid, and a direct numerical
minimisation (Ritz) of the travel-time functional over a family of trial paths, with the straight ramp
and Snell's invariant as checks. Hamilton's principle is applied to the pendulum of P02.06 and to a
spring-mass system, and the travel time of any path through any speed field (Fermat's principle) is
computed from scratch.

Conventions:
    * x horizontal [m], y measured DOWNWARD from the start [m], so y >= 0 along a gravity chute and the
      speed at depth y is v = sqrt(2 g y) (frictionless, released from rest).
    * A path is a pair of arrays (x, y) of points joined by straight segments.
    * Cycloid through the origin: x = k (theta - sin theta), y = k (1 - cos theta), 0 <= theta <= theta_f.
    * Functional J[y] = integral F(x, y, y') dx; Euler-Lagrange: dF/dy - d/dx (dF/dy') = 0.

Run:  python main.py          (prints the chute comparison and the symbolic results, writes brachistochrone.png)
Test: pytest curriculum/phase-02/lesson-13-calculus-of-variations
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.optimize import brentq, minimize

G0 = 9.80665   # standard gravity [m/s^2]


def chute() -> dict[str, float]:
    """Return the lesson's gravity chute: from (0, 0) to 4.0 m along and 2.0 m down, g = 9.80665 m/s^2.

    A frictionless parts chute between two deck levels of a workshop; the numbers are an example, not
    a specification.
    """
    return {"x_f": 4.0, "y_f": 2.0, "g": G0}


# ---------------------------------------------------------------- symbolic variational calculus

def euler_lagrange(F: sp.Expr, y: sp.Function, x: sp.Symbol) -> sp.Expr:
    """Return the Euler-Lagrange expression dF/dy - d/dx(dF/dy') for the integrand F(x, y(x), y'(x)).

    ``y`` must be an applied SymPy function y(x) and F an expression in x, y(x) and Derivative(y(x), x).
    The extremals of int F dx satisfy (returned expression) == 0.
    """
    raise NotImplementedError


def beltrami_first_integral(F: sp.Expr, y: sp.Function, x: sp.Symbol) -> sp.Expr:
    """Return F - y' dF/dy', which is constant along extremals when F has no explicit x dependence.

    Raise ``ValueError`` if F depends explicitly on x (then the identity does not hold).
    """
    raise NotImplementedError


def lagrangian_equation(T: sp.Expr, V: sp.Expr, q: sp.Function, t: sp.Symbol) -> sp.Expr:
    """Return the equation of motion d/dt(dL/dq') - dL/dq (== 0) from Hamilton's principle, L = T - V.

    ``q`` is an applied function q(t); T and V are expressions in q(t), Derivative(q(t), t) and t.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- travel time of a path

def path_time(x: np.ndarray, y: np.ndarray, speed: Callable[[np.ndarray, np.ndarray], np.ndarray]) -> float:
    """Return the travel time along the polygonal path (x, y) through the speed field speed(x, y).

    Each segment of length ds between speeds v1 and v2 is charged 2 ds / (v1 + v2), which is exact when
    v^2 varies linearly along the segment (the gravity chute, v = sqrt(2 g y) with y linear) and
    second-order accurate otherwise. v1 = 0 at the start is handled without division by zero.
    Raise ``ValueError`` if any speed is negative or both ends of a segment have zero speed.
    """
    raise NotImplementedError


def gravity_speed(g: float) -> Callable[[np.ndarray, np.ndarray], np.ndarray]:
    """Return the speed field v(x, y) = sqrt(2 g y) of a frictionless particle released from rest at y = 0."""
    raise NotImplementedError


def travel_time(x: np.ndarray, y: np.ndarray, g: float = G0) -> float:
    """Return the time to slide along the path (x, y) under gravity from rest at (x[0], 0): path_time with sqrt(2 g y)."""
    raise NotImplementedError


def straight_line_time(x_f: float, y_f: float, g: float = G0) -> float:
    """Return the slide time down a straight ramp from (0, 0) to (x_f, y_f): T = L sqrt(2 / (g y_f)), L = sqrt(x_f^2 + y_f^2).

    Constant acceleration g y_f / L along a ramp of length L from rest. Raise ``ValueError`` if y_f <= 0.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- the brachistochrone

def cycloid_parameter(x_f: float, y_f: float) -> tuple[float, float]:
    """Return (k, theta_f) of the cycloid x = k(theta - sin theta), y = k(1 - cos theta) through (x_f, y_f).

    theta_f solves (theta - sin theta) / (1 - cos theta) = x_f / y_f on (0, 2 pi) by Brent's method; then
    k = y_f / (1 - cos theta_f). Raise ``ValueError`` unless x_f > 0 and y_f > 0.
    """
    raise NotImplementedError


def brachistochrone(x_f: float, y_f: float, n: int = 401) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Return (x, y, k, theta_f): n points of the cycloid from (0, 0) to (x_f, y_f), equally spaced in theta."""
    raise NotImplementedError


def brachistochrone_time(k: float, theta_f: float, g: float = G0) -> float:
    """Return the slide time along the cycloid, T = theta_f sqrt(k / g) (exact: ds / v = sqrt(k / g) d theta)."""
    raise NotImplementedError


def snell_invariant(x: np.ndarray, y: np.ndarray, speed: Callable[[np.ndarray, np.ndarray], np.ndarray]) -> np.ndarray:
    """Return sin(theta) / v at each segment midpoint, theta the angle of the path from the y axis (the gradient direction).

    Fermat's principle in a medium stratified in y makes this constant along a minimum-time path
    (Snell's law in differential form); for the gravity chute it is 1 / sqrt(4 g k) on the cycloid.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- direct minimisation (Ritz)

def ritz_path(coeffs: np.ndarray, x_f: float, y_f: float, n_points: int = 401) -> tuple[np.ndarray, np.ndarray]:
    """Return the trial path y = y_f x / x_f + sum_k a_k sin(k pi x / x_f) on n_points, which meets both end conditions for every a."""
    raise NotImplementedError


def ritz_minimum_time(x_f: float, y_f: float, n_modes: int, g: float = G0, n_points: int = 401) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    """Minimise the travel time over the Ritz family and return (x, y, T, coeffs).

    Nelder-Mead from a = 0 (the straight ramp); a trial path that rises above the start (y < 0) is
    penalised, since the particle could not get there. The minimum must lie above the cycloid time and
    approach it as n_modes grows.
    """
    raise NotImplementedError


def perturbed_time_increase(x: np.ndarray, y: np.ndarray, eps: float, g: float = G0) -> float:
    """Return T[y + eps * eta] - T[y] for eta = sin(pi x / x_f) (vanishing at both ends).

    On an extremal the first variation vanishes, so the increase is O(eps^2): the test checks the exponent.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Derive the brachistochrone symbolically, compare cycloid, ramp and Ritz paths, and plot."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    xs = sp.symbols("x")
    y = sp.Function("y")(xs)
    yp = sp.Derivative(y, xs)
    g = sp.symbols("g", positive=True)

    # Shortest path in the plane.
    F_len = sp.sqrt(1 + yp ** 2)
    print("shortest path: Euler-Lagrange gives", euler_lagrange(F_len, y, xs), "= 0  (so y'' = 0: a straight line)")
    # Brachistochrone integrand and Beltrami identity.
    F_br = sp.sqrt((1 + yp ** 2) / (2 * g * y))
    belt = beltrami_first_integral(F_br, y, xs)
    print("brachistochrone: Beltrami first integral", belt, "= const,")
    print("   so 1/(2 g C^2) =", sp.simplify(1 / (2 * g * belt ** 2)), "= const: y (1 + y'^2) = 2k, the cycloid's differential equation")
    # Hamilton's principle: pendulum and spring-mass.
    t = sp.symbols("t")
    th = sp.Function("theta")(t)
    m, l, k = sp.symbols("m l k", positive=True)
    eom = lagrangian_equation(sp.Rational(1, 2) * m * l ** 2 * sp.Derivative(th, t) ** 2, -m * g * l * sp.cos(th), th, t)
    print("pendulum from Hamilton's principle:", sp.simplify(eom / (m * l ** 2)), "= 0")
    q = sp.Function("q")(t)
    print("spring-mass:", lagrangian_equation(sp.Rational(1, 2) * m * sp.Derivative(q, t) ** 2, sp.Rational(1, 2) * k * q ** 2, q, t), "= 0")

    c = chute()
    x_f, y_f, gg = c["x_f"], c["y_f"], c["g"]
    xb, yb, kk, th_f = brachistochrone(x_f, y_f)
    T_cyc = brachistochrone_time(kk, th_f, gg)
    T_cyc_num = travel_time(xb, yb, gg)
    T_line = straight_line_time(x_f, y_f, gg)
    xl = np.linspace(0.0, x_f, 401)
    T_line_num = travel_time(xl, y_f * xl / x_f, gg)
    print(f"\nCHUTE from (0, 0) to ({x_f}, {y_f}) m: cycloid k = {kk:.4f} m, theta_f = {th_f:.4f} rad ({np.degrees(th_f):.1f} deg),"
          f" lowest point y = {2 * kk:.3f} m at x = {kk * np.pi:.3f} m{' (beyond the end: no dip)' if kk * np.pi > x_f else ''}")
    print(f"cycloid time: closed form {T_cyc:.5f} s, path integral {T_cyc_num:.5f} s")
    print(f"straight ramp: closed form {T_line:.5f} s, path integral {T_line_num:.5f} s ({(T_line / T_cyc - 1) * 100:.1f} % slower)")
    # A circular arc through both points, tangent-vertical at the start.
    rr = (x_f ** 2 + y_f ** 2) / (2 * x_f)
    phi = np.linspace(0.0, np.arccos(1.0 - x_f / rr), 401)
    xc, yc = rr * (1 - np.cos(phi)), rr * np.sin(phi)
    T_arc = travel_time(xc, yc, gg)
    print(f"circular arc (vertical start): {T_arc:.5f} s ({(T_arc / T_cyc - 1) * 100:.1f} % slower)")
    sn = snell_invariant(xb, yb, gravity_speed(gg))
    print(f"Snell invariant sin(theta)/v on the cycloid: mean {sn.mean():.6f}, spread {np.ptp(sn):.1e}, prediction 1/sqrt(4 g k) = {1 / np.sqrt(4 * gg * kk):.6f}")
    for n_modes in (1, 2, 4, 6):
        xr, yr, T_r, a = ritz_minimum_time(x_f, y_f, n_modes, gg)
        print(f"Ritz with {n_modes} sine modes: T = {T_r:.5f} s, excess over the cycloid {(T_r / T_cyc - 1) * 100:.3f} %, a1 = {a[0]:.4f} m")
    for eps in (0.1, 0.05, 0.025):
        print(f"perturb the cycloid by eps = {eps} m: T rises by {perturbed_time_increase(xb, yb, eps, gg):.3e} s")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    ax1.plot(xb, yb, label=f"cycloid, {T_cyc:.3f} s")
    ax1.plot(xl, y_f * xl / x_f, "--", label=f"straight ramp, {T_line:.3f} s")
    ax1.plot(xc, yc, ":", label=f"circular arc, {T_arc:.3f} s")
    ax1.plot(xr, yr, ".", ms=2, label=f"Ritz, 6 modes, {T_r:.3f} s")
    ax1.invert_yaxis(); ax1.set_aspect("equal"); ax1.set_xlabel("x [m]"); ax1.set_ylabel("depth y [m]")
    ax1.set_title("gravity chute: candidate paths"); ax1.legend(fontsize=8); ax1.grid(True, lw=0.4)
    eps = np.array([0.2, 0.1, 0.05, 0.025, 0.0125])
    dT = np.array([perturbed_time_increase(xb, yb, e, gg) for e in eps])
    ax2.loglog(eps, dT, "o-", label="cycloid")
    xl2, yl2 = xl, y_f * xl / x_f
    dTl = np.array([abs(perturbed_time_increase(xl2, yl2, e, gg)) for e in eps])
    ax2.loglog(eps, dTl, "s--", label="straight ramp")
    ax2.loglog(eps, dT[0] * (eps / eps[0]) ** 2, "k:", label="slope 2")
    ax2.set_xlabel("perturbation amplitude eps [m]"); ax2.set_ylabel("|T[y + eps eta] - T[y]| [s]")
    ax2.set_title("first variation vanishes only on the extremal"); ax2.legend(fontsize=8); ax2.grid(True, which="both", lw=0.4)
    fig.suptitle("P02.13  Brachistochrone")
    fig.tight_layout()
    fig.savefig(out_dir / "brachistochrone.png", dpi=150)
    print(f"\nsaved {out_dir / 'brachistochrone.png'}")


if __name__ == "__main__":
    main()
