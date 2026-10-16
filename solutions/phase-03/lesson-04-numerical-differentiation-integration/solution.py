"""P03.04 Numerical Differentiation & Integration: reference solution.

Finite-difference stencils with Richardson extrapolation, and Newton-Cotes (trapezoid, Simpson),
Romberg, Gauss-Legendre and adaptive Simpson quadrature, all written from scratch and applied to
the heating of a ballistic re-entry capsule: the deceleration (a derivative of the velocity
profile) and the total heat load (an integral of the heating rate), both of which have closed
forms in the Allen-Eggers model so that every numerical answer can be checked exactly.
``scipy.integrate`` and ``numpy.polynomial.legendre.leggauss`` are the oracles.

Conventions:
    * Altitude h [m] is the independent variable of the entry model; the capsule descends at a
      constant flight-path angle gamma below the horizontal, so dt = -dh / (V sin gamma).
    * Allen-Eggers (no gravity, no lift, exponential atmosphere rho = rho0 exp(-h/H), constant
      ballistic coefficient beta = m / (C_D A)): V(h) = V_e exp(-B exp(-h/H)), B = rho0 H / (2 beta sin gamma).
    * Convective stagnation-point heating (Sutton-Graves form): qdot = k sqrt(rho / R_n) V^3 [W/m^2].
    * A quadrature with n panels (trapezoid, Simpson) uses n + 1 equally spaced function values;
      Simpson needs n even. Gauss-Legendre with n points is exact for polynomials of degree 2n - 1.
    * Convergence order is the slope of log(error) against log(h) (or against log(1/n)).

Run:  python solution.py          (prints the stencil, quadrature and re-entry tables, writes reentry_heating.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-03/lesson-04-numerical-differentiation-integration
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Callable, NamedTuple

import numpy as np

EPS: float = 2.0 ** -52
G0: float = 9.80665


# ---------------------------------------------------------------- finite differences

def forward_diff(f: Callable[[float], float], x: float, h: float) -> float:
    """Return (f(x + h) - f(x)) / h: first order, truncation error h f''/2."""
    return (f(x + h) - f(x)) / h


def central_diff(f: Callable[[float], float], x: float, h: float) -> float:
    """Return (f(x + h) - f(x - h)) / (2h): second order, truncation error h^2 f'''/6."""
    return (f(x + h) - f(x - h)) / (2.0 * h)


def second_diff(f: Callable[[float], float], x: float, h: float) -> float:
    """Return (f(x + h) - 2 f(x) + f(x - h)) / h^2: second-order estimate of f'', truncation error h^2 f''''/12.

    Rounding error ~ 4 eps |f| / h^2, so the best h is about eps^(1/4) ~ 1e-4 and the best error ~ 1e-8.
    """
    return (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h)


def richardson(estimate: Callable[[float], float], h: float, order: int, ratio: float = 2.0) -> float:
    """Richardson extrapolation: combine estimate(h) and estimate(h / ratio) to cancel the leading error term.

    If estimate(h) = exact + C h^order + O(h^(order+1)), then
        (ratio^order * estimate(h / ratio) - estimate(h)) / (ratio^order - 1) = exact + O(h^(order+1)),
    and O(h^(order+2)) when the error expansion has only even powers (central differences, trapezoid).
    """
    if order < 1 or ratio <= 1.0:
        raise ValueError("need order >= 1 and ratio > 1")
    r = ratio ** order
    return (r * estimate(h / ratio) - estimate(h)) / (r - 1.0)


def optimal_step_central(f_bound: float, f3_bound: float, eps: float = EPS) -> float:
    """Return h* = (3 eps f_bound / f3_bound)^(1/3), the minimiser of h^2 f3/6 + eps f/h for the central difference.

    About 6e-6 for f ~ f''' ~ 1, with best error ~ 1e-11: two thirds of the digits, against half for the forward difference.
    """
    if f_bound <= 0.0 or f3_bound <= 0.0:
        raise ValueError("bounds must be positive")
    return (3.0 * eps * f_bound / f3_bound) ** (1.0 / 3.0)


# ---------------------------------------------------------------- quadrature

def trapezoid(f: Callable[[np.ndarray], np.ndarray], a: float, b: float, n: int) -> float:
    """Composite trapezoid rule with n panels: h [f_0/2 + f_1 + ... + f_{n-1} + f_n/2]. Error -(b - a) h^2 f''/12.

    Exact for linear functions; spectrally accurate for periodic integrands over a full period.
    Raise ``ValueError`` if n < 1.
    """
    if n < 1:
        raise ValueError("need at least one panel")
    x = np.linspace(a, b, n + 1)
    y = np.asarray(f(x), dtype=float)
    h = (b - a) / n
    return float(h * (0.5 * y[0] + np.sum(y[1:-1]) + 0.5 * y[-1]))


def simpson(f: Callable[[np.ndarray], np.ndarray], a: float, b: float, n: int) -> float:
    """Composite Simpson rule with n (even) panels: h/3 [f_0 + 4 f_1 + 2 f_2 + ... + 4 f_{n-1} + f_n]. Error -(b - a) h^4 f''''/180.

    Exact for cubics (a parabola through three points integrates a cubic exactly by symmetry).
    Raise ``ValueError`` if n is odd or < 2.
    """
    if n < 2 or n % 2:
        raise ValueError("Simpson needs an even number of panels >= 2")
    x = np.linspace(a, b, n + 1)
    y = np.asarray(f(x), dtype=float)
    h = (b - a) / n
    return float(h / 3.0 * (y[0] + 4.0 * np.sum(y[1:-1:2]) + 2.0 * np.sum(y[2:-1:2]) + y[-1]))


def romberg(f: Callable[[np.ndarray], np.ndarray], a: float, b: float, levels: int) -> np.ndarray:
    """Return the Romberg table R[i, j] (lower triangle), i = 0..levels-1: R[i, 0] is the trapezoid rule with 2^i panels
    and R[i, j] = (4^j R[i, j-1] - R[i-1, j-1]) / (4^j - 1) is its j-fold Richardson extrapolation.

    R[:, 1] is Simpson's rule, R[:, 2] is Boole's rule; the best estimate is R[-1, -1].
    Raise ``ValueError`` if levels < 1.
    """
    if levels < 1:
        raise ValueError("need at least one level")
    R = np.zeros((levels, levels))
    for i in range(levels):
        R[i, 0] = trapezoid(f, a, b, 2 ** i)
        for j in range(1, i + 1):
            R[i, j] = (4.0 ** j * R[i, j - 1] - R[i - 1, j - 1]) / (4.0 ** j - 1.0)
    return R


def gauss_legendre_nodes(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Return the n Gauss-Legendre nodes and weights on [-1, 1], computed from scratch.

    Nodes are the roots of the Legendre polynomial P_n, found by Newton's method from the
    initial guesses cos(pi (i + 3/4) / (n + 1/2)) using the three-term recurrence
    (k + 1) P_{k+1} = (2k + 1) x P_k - k P_{k-1} for P_n and P_n'; weights are 2 / ((1 - x^2) P_n'(x)^2).
    Equal to ``numpy.polynomial.legendre.leggauss(n)`` to rounding. Raise ``ValueError`` if n < 1.
    """
    if n < 1:
        raise ValueError("need at least one node")
    i = np.arange(n)
    x = np.cos(np.pi * (i + 0.75) / (n + 0.5))
    for _ in range(100):
        p0 = np.ones_like(x)
        p1 = x.copy()
        if n == 1:
            p1, p0 = x.copy(), np.ones_like(x)
        for k in range(1, n):
            p0, p1 = p1, ((2.0 * k + 1.0) * x * p1 - k * p0) / (k + 1.0)
        # p1 = P_n(x), p0 = P_{n-1}(x); P_n' = n (x P_n - P_{n-1}) / (x^2 - 1)
        if n == 1:
            dp = np.ones_like(x)
        else:
            dp = n * (x * p1 - p0) / (x * x - 1.0)
        dx = p1 / dp
        x = x - dx
        if np.max(np.abs(dx)) < 1e-15:
            break
    # Final evaluation of P_n' at the converged nodes for the weights.
    p0 = np.ones_like(x)
    p1 = x.copy()
    for k in range(1, n):
        p0, p1 = p1, ((2.0 * k + 1.0) * x * p1 - k * p0) / (k + 1.0)
    dp = np.ones_like(x) if n == 1 else n * (x * p1 - p0) / (x * x - 1.0)
    w = 2.0 / ((1.0 - x * x) * dp * dp)
    order = np.argsort(x)
    return x[order], w[order]


def gauss_legendre(f: Callable[[np.ndarray], np.ndarray], a: float, b: float, n: int) -> float:
    """n-point Gauss-Legendre quadrature of f on [a, b]: exact for polynomials of degree <= 2n - 1.

    The nodes are mapped by x = (a + b)/2 + (b - a)/2 t and the weights scaled by (b - a)/2.
    """
    t, w = gauss_legendre_nodes(n)
    x = 0.5 * (a + b) + 0.5 * (b - a) * t
    return float(0.5 * (b - a) * np.sum(w * np.asarray(f(x), dtype=float)))


def adaptive_simpson(f: Callable[[float], float], a: float, b: float, tol: float = 1e-10,
                     max_depth: int = 40) -> tuple[float, int]:
    """Adaptive Simpson quadrature: return (integral, number of function evaluations).

    On each interval compare Simpson on the whole with Simpson on the two halves; the difference
    over 15 estimates the error of the finer one (the h^4 error law), and the interval is split
    while that exceeds its share of ``tol``. The returned value includes the Richardson correction
    (S_fine + (S_fine - S_coarse) / 15). Raise ``RuntimeError`` if ``max_depth`` is exceeded.
    """
    count = [0]

    def fe(x: float) -> float:
        count[0] += 1
        return float(f(x))

    def simp(fa: float, fm: float, fb: float, width: float) -> float:
        return width / 6.0 * (fa + 4.0 * fm + fb)

    def recurse(a_: float, b_: float, fa: float, fm: float, fb: float, whole: float, tol_: float, depth: int) -> float:
        m = 0.5 * (a_ + b_)
        lm, rm = 0.5 * (a_ + m), 0.5 * (m + b_)
        flm, frm = fe(lm), fe(rm)
        left = simp(fa, flm, fm, m - a_)
        right = simp(fm, frm, fb, b_ - m)
        if depth > max_depth:
            raise RuntimeError("adaptive Simpson exceeded max_depth")
        if abs(left + right - whole) <= 15.0 * tol_:
            return left + right + (left + right - whole) / 15.0
        return recurse(a_, m, fa, flm, fm, left, tol_ / 2.0, depth + 1) + recurse(m, b_, fm, frm, fb, right, tol_ / 2.0, depth + 1)

    fa, fb = fe(a), fe(b)
    fm = fe(0.5 * (a + b))
    whole = simp(fa, fm, fb, b - a)
    return recurse(a, b, fa, fm, fb, whole, tol, 0), count[0]


def observed_order(hs: np.ndarray, errs: np.ndarray) -> float:
    """Return the least-squares slope of log(err) against log(h). Raise ``ValueError`` for < 2 points or non-positive entries."""
    hs = np.asarray(hs, dtype=float)
    errs = np.asarray(errs, dtype=float)
    if hs.size < 2 or hs.size != errs.size or np.any(hs <= 0.0) or np.any(errs <= 0.0):
        raise ValueError("need >= 2 positive (h, err) pairs")
    return float(np.polyfit(np.log(hs), np.log(errs), 1)[0])


# ---------------------------------------------------------------- the re-entry model (Allen-Eggers)

class EntryParams(NamedTuple):
    """Ballistic entry parameters, SI units."""

    V_e: float = 7800.0        # entry speed, m/s
    gamma: float = math.radians(3.0)   # flight-path angle below horizontal, rad
    beta: float = 300.0        # ballistic coefficient m / (C_D A), kg/m^2
    R_n: float = 1.0           # nose radius, m
    H: float = 7200.0          # atmospheric scale height, m
    rho0: float = 1.225        # sea-level density, kg/m^3
    k_sg: float = 1.7415e-4    # Sutton-Graves constant for Earth air, kg^0.5 / m


def entry_B(p: EntryParams) -> float:
    """Return the dimensionless Allen-Eggers parameter B = rho0 H / (2 beta sin gamma)."""
    return p.rho0 * p.H / (2.0 * p.beta * math.sin(p.gamma))


def density_exp(h: np.ndarray | float, p: EntryParams) -> np.ndarray | float:
    """Return the exponential-atmosphere density rho0 exp(-h / H) [kg/m^3]."""
    return p.rho0 * np.exp(-np.asarray(h, dtype=float) / p.H)


def entry_velocity(h: np.ndarray | float, p: EntryParams) -> np.ndarray | float:
    """Return V(h) = V_e exp(-B exp(-h/H)) [m/s]: the Allen-Eggers velocity profile (no gravity, no lift)."""
    return p.V_e * np.exp(-entry_B(p) * np.exp(-np.asarray(h, dtype=float) / p.H))


def deceleration_analytic(h: np.ndarray | float, p: EntryParams) -> np.ndarray | float:
    """Return the deceleration rho V^2 / (2 beta) [m/s^2] = -dV/dt along the straight-line path."""
    return density_exp(h, p) * entry_velocity(h, p) ** 2 / (2.0 * p.beta)


def deceleration_from_velocity(h: float, p: EntryParams, dh: float = 10.0) -> float:
    """Return -dV/dt = (dV/dh) V sin(gamma) with dV/dh from a central difference of ``entry_velocity`` (step dh in m).

    Checks the differentiation machinery against ``deceleration_analytic``.
    """
    dVdh = central_diff(lambda x: float(entry_velocity(x, p)), h, dh)
    return dVdh * float(entry_velocity(h, p)) * math.sin(p.gamma)


def peak_deceleration(p: EntryParams) -> tuple[float, float]:
    """Return (altitude [m], value [m/s^2]) of the maximum deceleration: h = H ln(2B), a_max = V_e^2 sin(gamma) / (2 e H)."""
    B = entry_B(p)
    return p.H * math.log(2.0 * B), p.V_e ** 2 * math.sin(p.gamma) / (2.0 * math.e * p.H)


def heating_rate(h: np.ndarray | float, p: EntryParams) -> np.ndarray | float:
    """Return the stagnation-point convective heating rate k sqrt(rho / R_n) V^3 [W/m^2] at altitude h."""
    return p.k_sg * np.sqrt(density_exp(h, p) / p.R_n) * entry_velocity(h, p) ** 3


def peak_heating_altitude(p: EntryParams) -> float:
    """Return the altitude of maximum heating rate, h = H ln(6B) (where d ln(qdot)/dh = 0)."""
    return p.H * math.log(6.0 * entry_B(p))


def heat_load_integrand(h: np.ndarray | float, p: EntryParams) -> np.ndarray | float:
    """Return qdot / (V sin gamma) [J/m^2 per m of altitude]: the heat load per unit altitude, since dt = -dh / (V sin gamma)."""
    return heating_rate(h, p) / (entry_velocity(h, p) * math.sin(p.gamma))


def heat_load_numeric(h_top: float, p: EntryParams, method: str = "simpson", n: int = 64) -> float:
    """Return the total heat load [J/m^2] from h_top down to sea level by integrating ``heat_load_integrand``.

    ``method`` is 'trapezoid', 'simpson' or 'gauss' (n panels or n points).
    """
    f = lambda h: heat_load_integrand(h, p)                     # noqa: E731
    if method == "trapezoid":
        return trapezoid(f, 0.0, h_top, n)
    if method == "simpson":
        return simpson(f, 0.0, h_top, n)
    if method == "gauss":
        return gauss_legendre(f, 0.0, h_top, n)
    raise ValueError(f"unknown method {method!r}")


def heat_load_analytic(h_top: float, p: EntryParams) -> float:
    """Return the closed-form heat load [J/m^2] from h_top to 0.

    With u = exp(-h/H): integrand dh = k sqrt(rho0/R_n) V_e^2 H / sin(gamma) * u^(-1/2) exp(-2 B u) du, and
    int u^(-1/2) exp(-2 B u) du = sqrt(pi / (2B)) erf(sqrt(2 B u)), so
    Q = C sqrt(pi / (2B)) [erf(sqrt(2B)) - erf(sqrt(2B exp(-h_top/H)))].
    """
    B = entry_B(p)
    C = p.k_sg * math.sqrt(p.rho0 / p.R_n) * p.V_e ** 2 * p.H / math.sin(p.gamma)
    return C * math.sqrt(math.pi / (2.0 * B)) * (math.erf(math.sqrt(2.0 * B)) - math.erf(math.sqrt(2.0 * B * math.exp(-h_top / p.H))))


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Measure stencil and quadrature orders, compute the entry deceleration and heat load, and plot."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy import integrate

    out_dir = Path(out_dir)

    # --- stencils on sin at x = 1
    print("derivatives of sin at x = 1 (exact cos 1 = 0.540302305868):")
    hs = 10.0 ** -np.arange(1, 9, dtype=float)
    print(f"{'h':>8}{'forward':>14}{'central':>14}{'Richardson(c)':>16}{'second':>14}")
    ef, ec, er, e2 = [], [], [], []
    for h in hs:
        ef.append(abs(forward_diff(math.sin, 1.0, h) - math.cos(1.0)))
        ec.append(abs(central_diff(math.sin, 1.0, h) - math.cos(1.0)))
        er.append(abs(richardson(lambda s: central_diff(math.sin, 1.0, s), h, 2) - math.cos(1.0)))
        e2.append(abs(second_diff(math.sin, 1.0, h) + math.sin(1.0)))
        print(f"{h:>8.0e}{ef[-1]:>14.2e}{ec[-1]:>14.2e}{er[-1]:>16.2e}{e2[-1]:>14.2e}")
    print(f"observed orders (h = 1e-1..1e-3): forward {observed_order(hs[:3], ef[:3]):.2f}, central {observed_order(hs[:3], ec[:3]):.2f}, "
          f"second {observed_order(hs[:3], e2[:3]):.2f}; Richardson of central (h = 1e-1..1e-2, before its rounding floor) {observed_order(hs[:2], er[:2]):.2f}")
    print(f"optimal central step h* = {optimal_step_central(1.0, 1.0):.1e}, best observed central error {min(ec):.1e} at h = {hs[int(np.argmin(ec))]:.0e}")

    # --- quadrature on a smooth integrand
    f = lambda x: np.exp(-x) * np.cos(3.0 * x)                 # noqa: E731
    exact = (1.0 - math.exp(-2.0) * (math.cos(6.0) - 3.0 * math.sin(6.0))) / 10.0
    print(f"\nintegral of exp(-x) cos(3x) on [0, 2] = {exact:.15f}")
    print(f"{'panels':>7}{'trapezoid err':>15}{'Simpson err':>14}{'Gauss n err':>14}")
    ns = [4, 8, 16, 32, 64]
    et, es, eg = [], [], []
    for n in ns:
        et.append(abs(trapezoid(f, 0.0, 2.0, n) - exact))
        es.append(abs(simpson(f, 0.0, 2.0, n) - exact))
        eg.append(abs(gauss_legendre(f, 0.0, 2.0, n) - exact))
        print(f"{n:>7}{et[-1]:>15.2e}{es[-1]:>14.2e}{eg[-1]:>14.2e}")
    print(f"observed orders: trapezoid {observed_order(1 / np.array(ns), et):.2f}, Simpson {observed_order(1 / np.array(ns), es):.2f}; "
          f"Gauss falls faster than any power (spectral)")
    R = romberg(f, 0.0, 2.0, 6)
    print(f"Romberg R[5,5] error {abs(R[5, 5] - exact):.1e} from 33 evaluations; column 1 equals Simpson: {abs(R[3, 1] - simpson(f, 0.0, 2.0, 8)):.1e}")
    q, nev = adaptive_simpson(lambda x: math.exp(-x) * math.cos(3.0 * x), 0.0, 2.0, 1e-10)
    print(f"adaptive Simpson: error {abs(q - exact):.1e} with {nev} evaluations; scipy quad error {abs(integrate.quad(f, 0.0, 2.0)[0] - exact):.1e}")
    x_g, w_g = gauss_legendre_nodes(5)
    print(f"5-point Gauss nodes {np.round(x_g, 10)}, weights sum {w_g.sum():.15f}")

    # --- re-entry
    p = EntryParams()
    B = entry_B(p)
    h_q = peak_heating_altitude(p)
    h_a, a_max = peak_deceleration(p)
    print(f"\nballistic entry: V_e = {p.V_e:.0f} m/s, gamma = {math.degrees(p.gamma):.0f} deg, beta = {p.beta:.0f} kg/m^2, R_n = {p.R_n} m, B = {B:.1f}")
    print(f"peak heating at h = {h_q / 1e3:.1f} km: qdot = {heating_rate(h_q, p) / 1e4:.1f} W/cm^2, V = {entry_velocity(h_q, p):.0f} m/s")
    print(f"peak deceleration at h = {h_a / 1e3:.1f} km: {a_max:.1f} m/s^2 = {a_max / G0:.1f} g; "
          f"from a central difference of V: {deceleration_from_velocity(h_a, p):.1f} m/s^2")
    Q = heat_load_analytic(120e3, p)
    print(f"heat load from 120 km to sea level: analytic {Q / 1e6:.4f} MJ/m^2 (erf closed form)")
    print(f"{'n':>5}{'trapezoid':>14}{'Simpson':>14}{'Gauss':>14}   (relative errors)")
    for n in (8, 16, 32, 64, 128):
        print(f"{n:>5}" + "".join(f"{abs(heat_load_numeric(120e3, p, m, n) - Q) / Q:>14.2e}" for m in ("trapezoid", "simpson", "gauss")))
    hh = np.linspace(0.0, 120e3, 1001)
    q_s = integrate.simpson(heat_load_integrand(hh, p), x=hh)
    print(f"scipy.integrate.simpson on 1001 samples: relative error {abs(q_s - Q) / Q:.1e}; quad: {abs(integrate.quad(lambda h: heat_load_integrand(h, p), 0.0, 120e3)[0] - Q) / Q:.1e}")
    # numeric peak heating altitude: root of the central difference of qdot
    g = lambda h: central_diff(lambda x: float(heating_rate(x, p)), h, 1.0)   # noqa: E731
    lo, hi = 20e3, 80e3
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if g(lo) * g(mid) <= 0.0:
            hi = mid
        else:
            lo = mid
    print(f"peak heating altitude by bisection on the central difference: {0.5 * (lo + hi) / 1e3:.3f} km (analytic {h_q / 1e3:.3f} km)")

    # --- figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    ax = axes[0]
    ax.loglog(hs, ef, "o-", label="forward (p = 1)")
    ax.loglog(hs, ec, "s-", label="central (p = 2)")
    ax.loglog(hs, er, "^-", label="Richardson of central (p = 4)")
    ax.loglog(hs, e2, "v-", label="second derivative (p = 2)")
    ax.set_xlabel("step h")
    ax.set_ylabel("|error|")
    ax.set_title("Stencils on sin at x = 1")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    ax = axes[1]
    ax.loglog(ns, et, "o-", label="trapezoid")
    ax.loglog(ns, es, "s-", label="Simpson")
    ax.semilogy(ns, np.maximum(eg, 1e-17), "^-", label="Gauss-Legendre")
    ax.set_xlabel("panels / points n")
    ax.set_ylabel("|error|")
    ax.set_title("Quadrature of exp(-x) cos 3x on [0, 2]")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    ax = axes[2]
    h_km = hh / 1e3
    ax.plot(h_km, heating_rate(hh, p) / 1e4, label="heating rate [W/cm²]")
    ax.plot(h_km, deceleration_analytic(hh, p) / G0 * 10, label="deceleration [g] × 10")
    ax.axvline(h_q / 1e3, color="gray", ls=":", label=f"peak q at {h_q / 1e3:.1f} km")
    ax.axvline(h_a / 1e3, color="gray", ls="--", label=f"peak a at {h_a / 1e3:.1f} km")
    ax.set_xlabel("altitude [km]")
    ax.set_title("Ballistic entry (Allen-Eggers model)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "reentry_heating.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
