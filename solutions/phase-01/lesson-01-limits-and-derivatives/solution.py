"""P01.01 Limits, Continuity & Derivatives: reference solution.

Numerical derivatives (forward, backward, central, five-point, Richardson) compared with symbolic
derivatives from SymPy, applied to the lift-curve slope of an airfoil and to the kink in the ISA
temperature profile at the tropopause.

Lift model (thin-airfoil slope softened towards stall):
    C_L(alpha) = C_Lmax * tanh(a0 * alpha / C_Lmax),      dC_L/dalpha = a0 / cosh^2(a0 alpha / C_Lmax)

Run:  python solution.py          (prints the comparison tables, writes derivative_error_vs_h.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-01/lesson-01-limits-and-derivatives
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
import sympy as sp

A0_THIN_AIRFOIL = 2.0 * np.pi   # lift-curve slope of a thin airfoil, per radian
CL_MAX = 1.4                    # maximum lift coefficient of the model section
LAPSE_TROPOSPHERE = -0.0065     # K/m, ISA layer 1
H_TROPOPAUSE = 11000.0          # m (geopotential)
T_SEA_LEVEL = 288.15            # K
EPS = np.finfo(float).eps       # machine epsilon of float64, about 2.2e-16

# Shared symbols for the symbolic lift-curve slope.
alpha_s, a0_s, clmax_s = sp.symbols("alpha a0 C_Lmax", positive=True)


def model_symbols() -> dict[str, sp.Symbol]:
    """Return the symbols of the lift model by name: ``alpha`` (rad), ``a0`` (1/rad), ``clmax`` (-)."""
    return {"alpha": alpha_s, "a0": a0_s, "clmax": clmax_s}


# ---------------------------------------------------------------- the lift model

def lift_coefficient(alpha: np.ndarray | float, a0: float = A0_THIN_AIRFOIL, cl_max: float = CL_MAX) -> np.ndarray | float:
    """Return C_L = cl_max * tanh(a0 * alpha / cl_max) for angle of attack ``alpha`` in radians.

    Slope a0 at alpha = 0 (thin-airfoil theory gives 2 pi per radian), bending smoothly towards
    the plateau cl_max. Accepts a float or an array of any shape.
    """
    return cl_max * np.tanh(a0 * np.asarray(alpha, dtype=float) / cl_max)


def lift_curve_slope_exact(alpha: np.ndarray | float, a0: float = A0_THIN_AIRFOIL, cl_max: float = CL_MAX) -> np.ndarray | float:
    """Return dC_L/dalpha = a0 / cosh^2(a0 alpha / cl_max), the derivative of ``lift_coefficient`` by hand."""
    return a0 / np.cosh(a0 * np.asarray(alpha, dtype=float) / cl_max) ** 2


def symbolic_lift_slope() -> sp.Expr:
    """Return dC_L/dalpha as a SymPy expression in ``alpha``, ``a0``, ``C_Lmax`` using ``sp.diff``."""
    cl = clmax_s * sp.tanh(a0_s * alpha_s / clmax_s)
    return sp.diff(cl, alpha_s)


def lambdify_lift_slope() -> Callable[..., np.ndarray]:
    """Return a NumPy function ``slope(alpha, a0, cl_max)`` built from ``symbolic_lift_slope``."""
    return sp.lambdify((alpha_s, a0_s, clmax_s), symbolic_lift_slope(), modules="numpy")


# ---------------------------------------------------------------- limits

def limit_sequence(f: Callable[[np.ndarray], np.ndarray], a: float, hs: np.ndarray, side: str = "right") -> np.ndarray:
    """Return f(a + h) for ``side='right'`` or f(a - h) for ``side='left'`` at each positive step in ``hs``.

    This is the numerical face of ``lim_{x -> a+} f(x)``: watch the values settle as h shrinks.
    Raise ``ValueError`` for a side other than ``'left'``/``'right'`` or a non-positive h.
    """
    hs = np.asarray(hs, dtype=float)
    if np.any(hs <= 0.0):
        raise ValueError("every step in hs must be positive")
    if side == "right":
        return f(a + hs)
    if side == "left":
        return f(a - hs)
    raise ValueError("side must be 'left' or 'right'")


def one_sided_limits(f: Callable[[np.ndarray], np.ndarray], a: float, h0: float = 1e-2, n: int = 12) -> tuple[float, float]:
    """Estimate (lim_{x->a-} f, lim_{x->a+} f) from the values at a -/+ h0 * 2**-k, k = 0..n-1.

    The estimate is the value at the smallest step; it is not extrapolated, so a function whose
    limit is approached slowly (like |x|^0.1) will be off. Returns plain floats.
    """
    hs = h0 * 0.5 ** np.arange(n)
    left = limit_sequence(f, a, hs, "left")
    right = limit_sequence(f, a, hs, "right")
    return float(left[-1]), float(right[-1])


def is_continuous_at(f: Callable[[np.ndarray], np.ndarray], a: float, tol: float = 1e-6) -> bool:
    """Return True if the left limit, the right limit and f(a) agree to within ``tol`` (absolute)."""
    left, right = one_sided_limits(f, a)
    fa = float(f(np.asarray(a, dtype=float)))
    return abs(left - fa) <= tol and abs(right - fa) <= tol


def difference_quotients(f: Callable[[np.ndarray], np.ndarray], x: float, hs: np.ndarray) -> np.ndarray:
    """Return the forward difference quotients (f(x + h) - f(x)) / h for each h in ``hs``.

    The derivative is the limit of this sequence as h -> 0: the definition, made numerical.
    """
    hs = np.asarray(hs, dtype=float)
    if np.any(hs <= 0.0):
        raise ValueError("every step in hs must be positive")
    return (f(x + hs) - f(x)) / hs


# ---------------------------------------------------------------- finite differences

def _check_h(h: float) -> None:
    if not np.isfinite(h) or h <= 0.0:
        raise ValueError("h must be a positive finite number")


def forward_difference(f: Callable[[np.ndarray], np.ndarray], x: np.ndarray | float, h: float) -> np.ndarray | float:
    """Return (f(x + h) - f(x)) / h, first-order accurate: error ~ (h/2) f''."""
    _check_h(h)
    x = np.asarray(x, dtype=float)
    return (f(x + h) - f(x)) / h


def backward_difference(f: Callable[[np.ndarray], np.ndarray], x: np.ndarray | float, h: float) -> np.ndarray | float:
    """Return (f(x) - f(x - h)) / h, first-order accurate: error ~ -(h/2) f''."""
    _check_h(h)
    x = np.asarray(x, dtype=float)
    return (f(x) - f(x - h)) / h


def central_difference(f: Callable[[np.ndarray], np.ndarray], x: np.ndarray | float, h: float) -> np.ndarray | float:
    """Return (f(x + h) - f(x - h)) / (2h), second-order accurate: error ~ (h^2/6) f'''."""
    _check_h(h)
    x = np.asarray(x, dtype=float)
    return (f(x + h) - f(x - h)) / (2.0 * h)


def five_point_difference(f: Callable[[np.ndarray], np.ndarray], x: np.ndarray | float, h: float) -> np.ndarray | float:
    """Return (-f(x+2h) + 8 f(x+h) - 8 f(x-h) + f(x-2h)) / (12 h), fourth-order: error ~ (h^4/30) f^(5)."""
    _check_h(h)
    x = np.asarray(x, dtype=float)
    return (-f(x + 2.0 * h) + 8.0 * f(x + h) - 8.0 * f(x - h) + f(x - 2.0 * h)) / (12.0 * h)


def richardson_extrapolation(f: Callable[[np.ndarray], np.ndarray], x: np.ndarray | float, h: float) -> np.ndarray | float:
    """Return (4 D(h/2) - D(h)) / 3 where D is the central difference: the h^2 error term cancels.

    Algebraically this is the five-point formula with step h/2; the test checks that identity.
    """
    _check_h(h)
    return (4.0 * central_difference(f, x, 0.5 * h) - central_difference(f, x, h)) / 3.0


def one_sided_derivatives(f: Callable[[np.ndarray], np.ndarray], a: float, h: float) -> tuple[float, float]:
    """Return (left derivative, right derivative) at ``a`` from backward and forward differences.

    Each one-sided estimate uses Richardson extrapolation on the first-order formula,
    (2 D(h/2) - D(h)), so its error is O(h^2). They differ at a kink, where f'(a) does not exist.
    """
    _check_h(h)
    a = float(a)
    left = 2.0 * backward_difference(f, a, 0.5 * h) - backward_difference(f, a, h)
    right = 2.0 * forward_difference(f, a, 0.5 * h) - forward_difference(f, a, h)
    return float(left), float(right)


def isa_temperature_two_layer(h: np.ndarray | float) -> np.ndarray | float:
    """Return the ISA temperature in K for geopotential altitude ``h`` (m) in the first two layers.

    T = 288.15 - 0.0065 h below 11 km, 216.65 K from 11 km to 20 km: continuous at 11 km with a
    kink (the derivative jumps from -0.0065 K/m to 0). Raise ``ValueError`` outside [0, 20000] m.
    """
    h = np.asarray(h, dtype=float)
    if np.any(h < 0.0) or np.any(h > 20000.0):
        raise ValueError("h must lie in [0, 20000] m for the two-layer model")
    t_trop = T_SEA_LEVEL + LAPSE_TROPOSPHERE * h
    t_tropopause = T_SEA_LEVEL + LAPSE_TROPOSPHERE * H_TROPOPAUSE
    out = np.where(h < H_TROPOPAUSE, t_trop, t_tropopause)
    return float(out) if out.ndim == 0 else out


def gust_lift_increment(slope: float, rho: float, v: float, s: float, u_gust: float) -> float:
    """Return the lift increment dL = (1/2) rho V^2 S * slope * (u_gust / V) in N from a vertical gust.

    A vertical gust u_gust (m/s) met at airspeed V (m/s) changes the angle of attack by about
    u_gust / V (small angle); the lift-curve slope (1/rad) turns that into a lift change. This is the
    derivative used as a sensitivity: dL = (dL/dalpha) * dalpha. Raise ``ValueError`` if V <= 0.
    """
    if v <= 0.0:
        raise ValueError("airspeed must be positive")
    return 0.5 * rho * v**2 * s * slope * (u_gust / v)


# ---------------------------------------------------------------- error analysis

METHODS: dict[str, Callable[..., np.ndarray | float]] = {
    "forward": forward_difference,
    "central": central_difference,
    "five_point": five_point_difference,
    "richardson": richardson_extrapolation,
}


def error_vs_step(f: Callable[[np.ndarray], np.ndarray], dfdx: Callable[[np.ndarray], np.ndarray],
                  x: float, hs: np.ndarray) -> dict[str, np.ndarray]:
    """Return {method: |numerical - exact| at each h in hs} for forward, central, five_point, richardson."""
    hs = np.asarray(hs, dtype=float)
    exact = float(dfdx(np.asarray(x, dtype=float)))
    return {name: np.array([abs(float(method(f, x, h)) - exact) for h in hs]) for name, method in METHODS.items()}


def observed_order(hs: np.ndarray, errs: np.ndarray) -> float:
    """Return the slope of log(err) against log(h) by least squares: the observed order of accuracy.

    Use only points in the truncation regime (errors falling as h shrinks); rounding noise at tiny h
    will drag the slope towards -1. Raise ``ValueError`` if any error is zero or fewer than two points.
    """
    hs, errs = np.asarray(hs, dtype=float), np.asarray(errs, dtype=float)
    if hs.size < 2 or np.any(errs <= 0.0) or np.any(hs <= 0.0):
        raise ValueError("need at least two points with positive h and positive error")
    slope, _ = np.polyfit(np.log(hs), np.log(errs), 1)
    return float(slope)


def optimal_step_central(f_scale: float, f3_scale: float, eps: float = EPS) -> float:
    """Return h_opt = (3 eps |f| / |f'''|)^(1/3), where truncation (h^2 f'''/6) and rounding (eps |f|/h) balance.

    ``f_scale`` is the size of f near the point and ``f3_scale`` the size of its third derivative.
    Raise ``ValueError`` if either scale is not positive.
    """
    if f_scale <= 0.0 or f3_scale <= 0.0:
        raise ValueError("scales must be positive")
    return float((3.0 * eps * f_scale / f3_scale) ** (1.0 / 3.0))


# ---------------------------------------------------------------- driver

def main(out_dir: str | Path = ".") -> None:
    """Print the symbolic/numerical comparison and save the error-versus-step plot."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    print("symbolic dC_L/dalpha =", symbolic_lift_slope())
    slope_fn = lambdify_lift_slope()
    print(f"{'alpha [deg]':>12} {'C_L':>8} {'slope exact':>12} {'central h=1e-3':>15} {'|diff|':>10}")
    for deg in (0.0, 4.0, 8.0, 12.0, 16.0):
        a = np.radians(deg)
        cl = lift_coefficient(a)
        ex = float(slope_fn(a, A0_THIN_AIRFOIL, CL_MAX))
        num = float(central_difference(lift_coefficient, a, 1e-3))
        print(f"{deg:12.1f} {cl:8.4f} {ex:12.5f} {num:15.5f} {abs(ex - num):10.2e}")

    # Gust response: a 5 m/s vertical gust on a light aircraft (S = 15 m^2, W = 10 kN) at 50 m/s.
    rho, v, s_wing, u_gust = 1.225, 50.0, 15.0, 5.0
    q_s = 0.5 * rho * v**2 * s_wing
    for deg in (0.0, 8.0):
        a = np.radians(deg)
        dl = gust_lift_increment(float(slope_fn(a, A0_THIN_AIRFOIL, CL_MAX)), rho, v, s_wing, u_gust)
        dl_exact = q_s * (lift_coefficient(a + u_gust / v) - lift_coefficient(a))
        print(f"gust of {u_gust} m/s at {v} m/s (d alpha = {np.degrees(u_gust / v):.1f} deg) from alpha = {deg:.0f} deg:"
              f" linear dL = {dl / 1e3:.2f} kN ({dl / 10e3:.2f} g on 10 kN), exact {dl_exact / 1e3:.2f} kN,"
              f" linear over by {100 * (dl / dl_exact - 1):.0f} %")

    # Limits: the derivative of sin at 0 is the limit of sin(h)/h.
    hs = 10.0 ** -np.arange(1, 9, dtype=float)
    print("\nsin(h)/h for h = 1e-1 .. 1e-8:", np.array2string(limit_sequence(lambda x: np.sin(x) / x, 0.0, hs), precision=10))

    # The ISA kink.
    left, right = one_sided_derivatives(isa_temperature_two_layer, H_TROPOPAUSE, 10.0)
    print(f"continuous at 11 km: {is_continuous_at(isa_temperature_two_layer, H_TROPOPAUSE)};"
          f" left derivative {left:.5f} K/m, right derivative {right:.5f} K/m")

    # Error against step size at alpha = 8 deg.
    a8 = np.radians(8.0)
    hs = 10.0 ** -np.arange(1, 13, dtype=float)
    errs = error_vs_step(lift_coefficient, lift_curve_slope_exact, a8, hs)
    print(f"\n{'h':>8} " + " ".join(f"{k:>12}" for k in errs))
    for i, h in enumerate(hs):
        print(f"{h:8.0e} " + " ".join(f"{errs[k][i]:12.2e}" for k in errs))
    # Fit the order only where truncation dominates: the five-point error reaches rounding by h = 1e-4.
    for name, window in (("forward", slice(1, 4)), ("central", slice(1, 4)), ("five_point", slice(1, 3))):
        print(f"observed order of {name} on h = {hs[window][0]:.0e}..{hs[window][-1]:.0e}: "
              f"{observed_order(hs[window], errs[name][window]):.2f}")
    cl_expr = clmax_s * sp.tanh(a0_s * alpha_s / clmax_s)
    f3 = abs(float(sp.diff(cl_expr, alpha_s, 3).subs({a0_s: A0_THIN_AIRFOIL, clmax_s: CL_MAX, alpha_s: a8})))
    print(f"|C_L'''| at 8 deg = {f3:.2f} 1/rad^3;  h_opt (central) ~ {optimal_step_central(abs(lift_coefficient(a8)), f3):.1e} rad")

    fig, ax = plt.subplots(figsize=(6, 4))
    for name, marker in zip(errs, ("o", "s", "^", "x")):
        ax.loglog(hs, errs[name], marker=marker, label=name)
    ax.loglog(hs, EPS * abs(lift_coefficient(a8)) / hs, "k:", label=r"rounding $\epsilon |f| / h$")
    ax.set_xlabel("step h [rad]")
    ax.set_ylabel(r"|error| in $dC_L/d\alpha$ [1/rad]")
    ax.set_title(r"Finite-difference error at $\alpha$ = 8°")
    ax.grid(True, which="both", lw=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "derivative_error_vs_h.png", dpi=150)
    print(f"saved {out_dir / 'derivative_error_vs_h.png'}")


if __name__ == "__main__":
    main()
