"""P01.03 Integrals & the Fundamental Theorem: reference solution.

Area and centroid of a NACA four-digit airfoil section by Riemann sums, by the Fundamental Theorem
(antiderivative evaluated at the ends) and by SymPy, with a substitution that removes the
square-root singularity at the leading edge.

NACA four-digit half-thickness (x measured from the leading edge, 0 <= x <= c):
    y_t(x) = 5 t c [0.2969 sqrt(x/c) - 0.1260 (x/c) - 0.3516 (x/c)^2 + 0.2843 (x/c)^3 - 0.1015 (x/c)^4]
Camber line (m = max camber, p = its position, as fractions of c):
    y_c = m x (2p - x/c) / p^2                    for x <  p c
    y_c = m (c - x) (1 + x/c - 2p) / (1 - p)^2    for x >= p c
Section area A = 2 int_0^c y_t dx; centroid x_bar = (2/A) int x y_t dx, y_bar = (2/A) int y_c y_t dx
(thickness applied vertically, a stated simplification of the NACA construction).

Run:  python solution.py          (prints the section table, writes airfoil_section.png and riemann_convergence.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-01/lesson-03-integrals-fundamental-theorem
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
import sympy as sp

# NACA four-digit thickness polynomial: coefficients and exponents of (x/c).
NACA_COEFFS = np.array([0.2969, -0.1260, -0.3516, 0.2843, -0.1015])
NACA_COEFFS_CLOSED = np.array([0.2969, -0.1260, -0.3516, 0.2843, -0.1036])
NACA_EXPONENTS = np.array([0.5, 1.0, 2.0, 3.0, 4.0])

RULES = ("left", "right", "midpoint")

# Shared symbols for the symbolic area: thickness ratio t, chord c, chordwise position x.
t_s, c_s, x_s = sp.symbols("t c x", positive=True)


def model_symbols() -> dict[str, sp.Symbol]:
    """Return the symbols used by ``symbolic_section_area``: ``t`` (-), ``c`` (m), ``x`` (m)."""
    return {"t": t_s, "c": c_s, "x": x_s}


# ---------------------------------------------------------------- the section geometry

def naca4_thickness(x: np.ndarray | float, t: float, c: float = 1.0, closed_te: bool = False) -> np.ndarray | float:
    """Return the half-thickness y_t(x) of a NACA four-digit section of thickness ratio ``t`` and chord ``c``.

    ``closed_te`` uses -0.1036 for the last coefficient so the trailing edge has zero thickness.
    Raise ``ValueError`` if any x lies outside [0, c] or if t <= 0 or c <= 0.
    """
    x = np.asarray(x, dtype=float)
    if t <= 0.0 or c <= 0.0:
        raise ValueError("thickness ratio and chord must be positive")
    if np.any(x < 0.0) or np.any(x > c):
        raise ValueError("x must lie in [0, c]")
    xi = x / c
    coeffs = NACA_COEFFS_CLOSED if closed_te else NACA_COEFFS
    out = 5.0 * t * c * sum(a * xi**e for a, e in zip(coeffs, NACA_EXPONENTS))
    return float(out) if np.ndim(out) == 0 else out


def naca4_camber(x: np.ndarray | float, m: float, p: float, c: float = 1.0) -> np.ndarray | float:
    """Return the camber line y_c(x) with maximum camber ``m`` (fraction of c) at ``p`` (fraction of c).

    m = 0 gives a symmetric section (zero everywhere). Raise ``ValueError`` unless 0 < p < 1 when m != 0.
    """
    x = np.asarray(x, dtype=float)
    if m == 0.0:
        out = np.zeros_like(x)
        return float(out) if out.ndim == 0 else out
    if not 0.0 < p < 1.0:
        raise ValueError("camber position p must lie strictly between 0 and 1")
    xi = x / c
    front = m * c * xi * (2.0 * p - xi) / p**2
    back = m * c * (1.0 - xi) * (1.0 + xi - 2.0 * p) / (1.0 - p) ** 2
    out = np.where(xi < p, front, back)
    return float(out) if out.ndim == 0 else out


# ---------------------------------------------------------------- Riemann sums

def riemann_sum(f: Callable[[np.ndarray], np.ndarray], a: float, b: float, n: int, rule: str = "midpoint") -> float:
    """Return the Riemann sum of f on [a, b] with ``n`` equal panels and ``rule`` in {'left', 'right', 'midpoint'}.

    Raise ``ValueError`` for n < 1 or an unknown rule. Works for b < a (the sign flips, as the integral does).
    """
    if n < 1:
        raise ValueError("n must be at least 1")
    if rule not in RULES:
        raise ValueError(f"rule must be one of {RULES}")
    h = (b - a) / n
    edges = a + h * np.arange(n + 1)
    if rule == "left":
        pts = edges[:-1]
    elif rule == "right":
        pts = edges[1:]
    else:
        pts = 0.5 * (edges[:-1] + edges[1:])
    return float(h * np.sum(f(pts)))


def cumulative_integral(f: Callable[[np.ndarray], np.ndarray], a: float, xs: np.ndarray, n: int = 400) -> np.ndarray:
    """Return F(x) = int_a^x f for each x in ``xs``, each by a midpoint sum with ``n`` panels.

    This is the 'accumulation function' of the Fundamental Theorem, part 1: F'(x) = f(x).
    """
    xs = np.atleast_1d(np.asarray(xs, dtype=float))
    return np.array([riemann_sum(f, a, x, n, "midpoint") if x != a else 0.0 for x in xs])


# ---------------------------------------------------------------- closed forms (the Fundamental Theorem, part 2)

def section_area_exact(t: float, c: float = 1.0, closed_te: bool = False) -> float:
    """Return A = 2 int_0^c y_t dx in closed form: 10 t c^2 * sum a_k / (e_k + 1)."""
    coeffs = NACA_COEFFS_CLOSED if closed_te else NACA_COEFFS
    return float(10.0 * t * c**2 * np.sum(coeffs / (NACA_EXPONENTS + 1.0)))


def section_centroid_x_exact(t: float, c: float = 1.0, closed_te: bool = False) -> float:
    """Return x_bar = (2/A) int_0^c x y_t dx in closed form: c * [sum a_k/(e_k+2)] / [sum a_k/(e_k+1)]."""
    coeffs = NACA_COEFFS_CLOSED if closed_te else NACA_COEFFS
    return float(c * np.sum(coeffs / (NACA_EXPONENTS + 2.0)) / np.sum(coeffs / (NACA_EXPONENTS + 1.0)))


def symbolic_section_area(closed_te: bool = False) -> sp.Expr:
    """Return A as a SymPy expression in symbols ``t`` and ``c``, by integrating the thickness polynomial exactly."""
    coeffs = NACA_COEFFS_CLOSED if closed_te else NACA_COEFFS
    yt = 5 * t_s * c_s * sum(sp.Rational(str(a)) * (x_s / c_s) ** sp.Rational(str(e)) for a, e in zip(coeffs, NACA_EXPONENTS))
    return sp.simplify(2 * sp.integrate(yt, (x_s, 0, c_s)))


def definite_from_antiderivative(expr: sp.Expr, x: sp.Symbol, a: float, b: float) -> float:
    """Return F(b) - F(a) with F = sp.integrate(expr, x): the Fundamental Theorem, part 2, as a function."""
    big_f = sp.integrate(expr, x)
    return float(big_f.subs(x, b) - big_f.subs(x, a))


# ---------------------------------------------------------------- numerical area and centroid

def section_area(t: float, c: float = 1.0, n: int = 1000, rule: str = "midpoint", closed_te: bool = False) -> float:
    """Return A = 2 int_0^c y_t dx by a Riemann sum with ``n`` panels."""
    return 2.0 * riemann_sum(lambda x: naca4_thickness(x, t, c, closed_te), 0.0, c, n, rule)


def section_area_substituted(t: float, c: float = 1.0, n: int = 1000, rule: str = "midpoint", closed_te: bool = False) -> float:
    """Return A with the substitution x = c u^2 (dx = 2 c u du), which turns sqrt(x/c) into u and removes the
    infinite slope at the leading edge: A = 2 int_0^1 y_t(c u^2) 2 c u du, a polynomial integrand in u."""
    return 2.0 * riemann_sum(lambda u: naca4_thickness(c * u**2, t, c, closed_te) * 2.0 * c * u, 0.0, 1.0, n, rule)


def section_centroid(t: float, c: float = 1.0, m: float = 0.0, p: float = 0.4, n: int = 1000, closed_te: bool = False) -> tuple[float, float]:
    """Return (x_bar, y_bar) of the section with thickness applied vertically about the camber line.

    Upper surface y_c + y_t, lower y_c - y_t: A = int 2 y_t dx, x_bar = (1/A) int x 2 y_t dx,
    y_bar = (1/A) int y_c 2 y_t dx. Midpoint sums with ``n`` panels.
    """
    yt = lambda x: naca4_thickness(x, t, c, closed_te)
    area = 2.0 * riemann_sum(yt, 0.0, c, n)
    sx = 2.0 * riemann_sum(lambda x: x * yt(x), 0.0, c, n)
    sy = 2.0 * riemann_sum(lambda x: naca4_camber(x, m, p, c) * yt(x), 0.0, c, n)
    return sx / area, sy / area


# ---------------------------------------------------------------- net change

def net_change(rate: Callable[[np.ndarray], np.ndarray], a: float, b: float, n: int = 1000) -> float:
    """Return int_a^b rate dt: the net change of the quantity whose rate is ``rate`` (Net Change Theorem)."""
    return riemann_sum(rate, a, b, n, "midpoint")


def total_variation(rate: Callable[[np.ndarray], np.ndarray], a: float, b: float, n: int = 1000) -> float:
    """Return int_a^b |rate| dt: the total distance travelled, as opposed to the displacement."""
    return riemann_sum(lambda s: np.abs(rate(s)), a, b, n, "midpoint")


# ---------------------------------------------------------------- convergence

def convergence_study(f: Callable[[np.ndarray], np.ndarray], exact: float, a: float, b: float,
                      ns: list[int] | np.ndarray, rule: str = "midpoint") -> np.ndarray:
    """Return |riemann_sum - exact| for each panel count in ``ns``."""
    return np.array([abs(riemann_sum(f, a, b, int(n), rule) - exact) for n in ns])


def observed_order(ns: np.ndarray, errs: np.ndarray) -> float:
    """Return the order of convergence p in err ~ n^-p from a least-squares fit of log(err) against log(n).

    Raise ``ValueError`` for fewer than two points or a non-positive error.
    """
    ns, errs = np.asarray(ns, dtype=float), np.asarray(errs, dtype=float)
    if ns.size < 2 or np.any(errs <= 0.0):
        raise ValueError("need at least two points with positive errors")
    slope, _ = np.polyfit(np.log(ns), np.log(errs), 1)
    return float(-slope)


# ---------------------------------------------------------------- driver

def main(out_dir: str | Path = ".") -> None:
    """Print the section areas and centroids, the convergence table, and save the two figures."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    t, c = 0.12, 1.0
    a_exact = section_area_exact(t, c)
    print(f"NACA 0012, c = {c} m: exact area {a_exact:.6f} m^2 = {a_exact / (t * c**2):.4f} t c^2,"
          f" x_bar = {section_centroid_x_exact(t, c):.4f} c")
    print(f"closed trailing edge: area {section_area_exact(t, c, True):.6f} m^2, x_bar = {section_centroid_x_exact(t, c, True):.4f} c")
    print("SymPy: A =", symbolic_section_area(), "->", float(symbolic_section_area().subs({t_s: t, c_s: c})))
    xb, yb = section_centroid(t, c, 0.02, 0.4, n=4000)
    print(f"NACA 2412 (thickness applied vertically): x_bar = {xb:.4f} c, y_bar = {yb:.5f} c")

    ns = np.array([10, 20, 40, 80, 160, 320, 640, 1280])
    f = lambda x: 2.0 * naca4_thickness(x, t, c)
    g = lambda u: 2.0 * naca4_thickness(c * u**2, t, c) * 2.0 * c * u
    print(f"\n{'n':>6} {'left':>11} {'right':>11} {'midpoint':>11} {'midpoint, u-subst':>18}")
    errs = {}
    for rule in RULES:
        errs[rule] = convergence_study(f, a_exact, 0.0, c, ns, rule)
    errs["subst"] = convergence_study(g, a_exact, 0.0, 1.0, ns, "midpoint")
    for i, n in enumerate(ns):
        print(f"{n:6d} {errs['left'][i]:11.2e} {errs['right'][i]:11.2e} {errs['midpoint'][i]:11.2e} {errs['subst'][i]:18.2e}")
    for key in ("left", "right", "midpoint", "subst"):
        print(f"observed order, {key}: {observed_order(ns[2:], errs[key][2:]):.2f}")

    # FTC part 1: the derivative of the accumulated area is the local thickness.
    xs = np.linspace(0.1, 0.9, 9)
    h = 1e-3
    big_f = lambda x: cumulative_integral(f, 0.0, x, n=2000)
    dfdx = (big_f(xs + h) - big_f(xs - h)) / (2.0 * h)
    print(f"\nFTC1 check: max |dF/dx - 2 y_t| on x in [0.1, 0.9] = {np.max(np.abs(dfdx - f(xs))):.2e}")

    # Net change: a glider in a thermal, vertical speed v(t) = 2 sin(pi t / 60) - 0.8 m/s over 0..120 s.
    v = lambda s: 2.0 * np.sin(np.pi * s / 60.0) - 0.8
    print(f"net height change over 120 s: {net_change(v, 0.0, 120.0):+.2f} m; total vertical distance: {total_variation(v, 0.0, 120.0):.2f} m")

    x = np.linspace(0.0, c, 400)
    yt, yc = naca4_thickness(x, t, c), naca4_camber(x, 0.02, 0.4, c)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    ax = axes[0]
    ax.fill_between(x, yc - yt, yc + yt, alpha=0.3, label="NACA 2412 (vertical thickness)")
    ax.plot(x, yt, "k--", lw=0.8, label="NACA 0012 upper")
    ax.plot(x, -yt, "k--", lw=0.8)
    ax.plot(xb, yb, "ro", label=f"centroid ({xb:.3f}, {yb:.4f}) c")
    ax.set_aspect("equal")
    ax.set_xlabel("x / c")
    ax.set_ylabel("y / c")
    ax.legend(fontsize=7)
    ax.grid(True, lw=0.3)
    ax = axes[1]
    for key, marker in zip(("left", "right", "midpoint", "subst"), ("v", "^", "o", "s")):
        ax.loglog(ns, errs[key], marker=marker, label=key)
    ax.loglog(ns, errs["midpoint"][0] * (ns / ns[0]) ** -1.5, "k:", lw=0.8, label=r"$n^{-1.5}$")
    ax.loglog(ns, errs["subst"][0] * (ns / ns[0]) ** -2.0, "k--", lw=0.8, label=r"$n^{-2}$")
    ax.set_xlabel("panels n")
    ax.set_ylabel("|area error| [m$^2$]")
    ax.grid(True, which="both", lw=0.3)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out_dir / "airfoil_section.png", dpi=150)
    fig2, ax2 = plt.subplots(figsize=(5, 3.5))
    for key, marker in zip(("left", "right", "midpoint", "subst"), ("v", "^", "o", "s")):
        ax2.loglog(ns, errs[key], marker=marker, label=key)
    ax2.set_xlabel("panels n")
    ax2.set_ylabel("|area error| [m$^2$]")
    ax2.grid(True, which="both", lw=0.3)
    ax2.legend(fontsize=8)
    fig2.tight_layout()
    fig2.savefig(out_dir / "riemann_convergence.png", dpi=150)
    print(f"saved {out_dir / 'airfoil_section.png'} and {out_dir / 'riemann_convergence.png'}")


if __name__ == "__main__":
    main()
