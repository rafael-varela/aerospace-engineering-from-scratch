"""P03.03 Interpolation & Splines: starter code.

Atmospheric table interpolator: the ISA (International Standard Atmosphere) is tabulated at a
coarse altitude spacing and recovered between the rows by polynomial (Lagrange, Newton),
piecewise-linear, log-linear and cubic-spline interpolation, all written from scratch, with
``scipy.interpolate`` as the oracle and the analytic ISA as the truth.

Complete every function marked NotImplementedError; ``main()`` already produces the lesson output once they work.

Conventions:
    * Altitudes are geopotential metres h; the ISA layers run to 86 km (seven lapse-rate layers).
    * ``isa_atmosphere(h)`` is the analytic reference (the "truth" the tables are sampled from);
      it is given, having been built in P00.03.
    * Nodes x_0 < x_1 < ... < x_n are strictly increasing; interpolants are evaluated only inside
      [x_0, x_n] and raise ``ValueError`` outside (no silent extrapolation).
    * A cubic spline is stored as the knots, the values and the second derivatives M_i at the
      knots (``Spline``); boundary conditions 'natural' (M_0 = M_n = 0), 'clamped' (given end
      slopes) or 'not-a-knot' (third derivative continuous at x_1 and x_{n-1}).
    * Errors are max |interpolant - truth| over a fine grid; convergence order is the slope of
      log(error) against log(spacing).

Run:  python main.py          (prints the interpolation tables, writes atmosphere_interpolation.png)
Test: pytest curriculum/phase-03/lesson-03-interpolation-splines
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Callable, NamedTuple

import numpy as np

# ---------------------------------------------------------------- the ISA (given)

G0: float = 9.80665           # m/s^2
R_AIR: float = 287.05287      # J/(kg K)
T0: float = 288.15            # K
P0: float = 101325.0          # Pa
# Layer base altitudes (geopotential m) and lapse rates (K/m) of the layer starting there.
ISA_BASES: np.ndarray = np.array([0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0, 84852.0])
ISA_LAPSE: np.ndarray = np.array([-6.5e-3, 0.0, 1.0e-3, 2.8e-3, 0.0, -2.8e-3, -2.0e-3])


def _isa_layer_bases() -> tuple[np.ndarray, np.ndarray]:
    """Return the temperature and pressure at each layer base, built layer by layer."""
    T_b = np.empty(ISA_BASES.size)
    p_b = np.empty(ISA_BASES.size)
    T_b[0], p_b[0] = T0, P0
    for i in range(ISA_BASES.size - 1):
        dh = ISA_BASES[i + 1] - ISA_BASES[i]
        L = ISA_LAPSE[i]
        T_b[i + 1] = T_b[i] + L * dh
        if L == 0.0:
            p_b[i + 1] = p_b[i] * math.exp(-G0 * dh / (R_AIR * T_b[i]))
        else:
            p_b[i + 1] = p_b[i] * (T_b[i + 1] / T_b[i]) ** (-G0 / (R_AIR * L))
    return T_b, p_b


def isa_atmosphere(h: np.ndarray | float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (T [K], p [Pa], rho [kg/m^3]) of the ISA at geopotential altitude h [m], 0 <= h <= 86 km.

    Given (from P00.03). Piecewise: T linear in each layer; p by the hydrostatic integral
    (power law in a lapse layer, exponential in an isothermal one); rho = p / (R T).
    Raise ``ValueError`` outside [0, 86000] m.
    """
    h = np.asarray(h, dtype=float)
    if np.any(h < 0.0) or np.any(h > 86000.0):
        raise ValueError("ISA is defined here for 0 <= h <= 86 km geopotential")
    T_b, p_b = _isa_layer_bases()
    idx = np.clip(np.searchsorted(ISA_BASES, h, side="right") - 1, 0, ISA_LAPSE.size - 1)
    L = ISA_LAPSE[idx]
    dh = h - ISA_BASES[idx]
    T = T_b[idx] + L * dh
    with np.errstate(divide="ignore", invalid="ignore"):
        p_lapse = p_b[idx] * (T / T_b[idx]) ** np.where(L == 0.0, 0.0, -G0 / (R_AIR * np.where(L == 0.0, 1.0, L)))
    p_iso = p_b[idx] * np.exp(-G0 * dh / (R_AIR * T_b[idx]))
    p = np.where(L == 0.0, p_iso, p_lapse)
    rho = p / (R_AIR * T)
    if T.ndim == 0:
        return float(T), float(p), float(rho)
    return T, p, rho


def isa_table(h_step: float, h_max: float = 86000.0) -> np.ndarray:
    """Return the ISA sampled every ``h_step`` m from 0 to ``h_max``: array of shape (n, 4) with columns h, T, p, rho.

    Raise ``ValueError`` if h_step <= 0 or h_max > 86000.
    """
    if h_step <= 0.0 or h_max > 86000.0 or h_max <= 0.0:
        raise ValueError("need 0 < h_step and 0 < h_max <= 86000")
    h = np.arange(0.0, h_max + 0.5 * h_step, h_step)
    h = h[h <= h_max + 1e-9]
    T, p, rho = isa_atmosphere(h)
    return np.column_stack([h, T, p, rho])


# ---------------------------------------------------------------- polynomial interpolation

def _check_nodes(xs: np.ndarray, ys: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    raise NotImplementedError


def _check_inside(xs: np.ndarray, x: np.ndarray) -> None:
    raise NotImplementedError


def lagrange_interpolate(xs: np.ndarray, ys: np.ndarray, x: np.ndarray | float) -> np.ndarray | float:
    """Evaluate the degree-n interpolating polynomial through (xs, ys) at x in the Lagrange form.

    p(x) = sum_i y_i l_i(x), l_i(x) = prod_{j != i} (x - x_j) / (x_i - x_j). O(n^2) per point.
    Raise ``ValueError`` for repeated nodes or x outside [x_0, x_n].
    """
    raise NotImplementedError


def newton_divided_differences(xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
    """Return the Newton coefficients c_k = f[x_0, ..., x_k] of the interpolating polynomial.

    p(x) = c_0 + c_1 (x - x_0) + c_2 (x - x_0)(x - x_1) + ... Built by the divided-difference
    table in O(n^2); adding a node adds one coefficient without recomputing the others.
    """
    raise NotImplementedError


def newton_evaluate(xs: np.ndarray, coeffs: np.ndarray, x: np.ndarray | float) -> np.ndarray | float:
    """Evaluate the Newton-form polynomial by nested multiplication (Horner on the centred factors)."""
    raise NotImplementedError


def chebyshev_nodes(n: int, a: float, b: float) -> np.ndarray:
    """Return n Chebyshev nodes of the first kind on [a, b], sorted increasing.

    x_k = (a + b)/2 + (b - a)/2 cos((2k + 1) pi / (2n)), k = 0..n-1. They cluster near the ends and
    make the node polynomial prod (x - x_k) as small as possible, which defeats the Runge phenomenon.
    """
    raise NotImplementedError


def runge(x: np.ndarray | float) -> np.ndarray | float:
    """Return 1 / (1 + 25 x^2), the classic function on which equispaced interpolation diverges."""
    raise NotImplementedError


# ---------------------------------------------------------------- piecewise interpolation

def piecewise_linear(xs: np.ndarray, ys: np.ndarray, x: np.ndarray | float) -> np.ndarray | float:
    """Linear interpolation between neighbouring nodes (what ``numpy.interp`` does). Error <= h^2 max|f''| / 8."""
    raise NotImplementedError


def log_linear_interpolate(xs: np.ndarray, ys: np.ndarray, x: np.ndarray | float) -> np.ndarray | float:
    """Interpolate ln(y) linearly and exponentiate: exact for y = A exp(B x), so for pressure in an isothermal layer.

    Raise ``ValueError`` if any y <= 0.
    """
    raise NotImplementedError


class Spline(NamedTuple):
    """A cubic spline: knots xs, values ys and second derivatives M at the knots."""

    xs: np.ndarray
    ys: np.ndarray
    M: np.ndarray


def cubic_spline(xs: np.ndarray, ys: np.ndarray, bc: str = "natural",
                 end_slopes: tuple[float, float] = (0.0, 0.0)) -> Spline:
    """Build the C^2 cubic spline through (xs, ys) by solving the linear system for the knot second derivatives M_i.

    Interior equations (continuity of the first derivative at each interior knot):
        h_{i-1} M_{i-1} + 2 (h_{i-1} + h_i) M_i + h_i M_{i+1} = 6 (d_i - d_{i-1}),
    with h_i = x_{i+1} - x_i and d_i = (y_{i+1} - y_i) / h_i. Boundary rows:
        'natural'    M_0 = M_n = 0;
        'clamped'    2 h_0 M_0 + h_0 M_1 = 6 (d_0 - s_0), h_{n-1} M_{n-1} + 2 h_{n-1} M_n = 6 (s_n - d_{n-1})
                     for prescribed end slopes (s_0, s_n) = ``end_slopes``;
        'not-a-knot' third derivative continuous at x_1 and x_{n-1} (needs >= 4 knots).
    The system is solved with ``numpy.linalg.solve`` (tridiagonal solvers come in P03.08).
    """
    raise NotImplementedError


def spline_evaluate(spline: Spline, x: np.ndarray | float, deriv: int = 0) -> np.ndarray | float:
    """Evaluate the spline (deriv = 0), its first (1) or second (2) derivative at x inside the knots.

    On [x_i, x_{i+1}] with h = x_{i+1} - x_i, a = x_{i+1} - x, b = x - x_i:
    S = M_i a^3 / (6h) + M_{i+1} b^3 / (6h) + (y_i - M_i h^2 / 6) a / h + (y_{i+1} - M_{i+1} h^2 / 6) b / h.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- error models and studies

def linear_error_bound(h: float, f2_max: float) -> float:
    """Return h^2 f2_max / 8: the maximum error of piecewise-linear interpolation on spacing h when |f''| <= f2_max."""
    raise NotImplementedError


def spline_error_bound(h: float, f4_max: float) -> float:
    """Return 5 h^4 f4_max / 384: the maximum error of the clamped (or not-a-knot, in the interior) cubic spline when |f''''| <= f4_max."""
    raise NotImplementedError


def observed_order(hs: np.ndarray, errs: np.ndarray) -> float:
    """Return the least-squares slope of log(err) against log(h). Raise ``ValueError`` for < 2 points or non-positive entries."""
    raise NotImplementedError


def max_interpolation_error(f: Callable[[np.ndarray], np.ndarray], a: float, b: float, n_intervals: int,
                            method: str, n_fine: int = 2001) -> float:
    """Return max |interpolant - f| on a fine grid of [a, b] for equispaced nodes with ``n_intervals`` intervals.

    ``method`` is 'linear', 'loglinear', 'natural', 'clamped', 'not-a-knot' or 'lagrange'. For
    'clamped' the end slopes are estimated by one-sided second-order differences of f with step (b - a) 1e-5.
    """
    raise NotImplementedError


def interpolate_isa(table: np.ndarray, h: np.ndarray | float, column: str = "rho", method: str = "natural") -> np.ndarray | float:
    """Interpolate the ISA ``table`` (from ``isa_table``) at altitude h for ``column`` 'T', 'p' or 'rho' with ``method``.

    ``method`` is 'linear', 'loglinear', 'natural' or 'not-a-knot'. For 'p' and 'rho' the log-linear
    form is the physically motivated one (both fall nearly exponentially with altitude).
    """
    raise NotImplementedError


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Interpolate a coarse ISA table four ways, measure errors and convergence orders, show Runge, and plot."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy import interpolate as si

    out_dir = Path(out_dir)
    T11, p11, rho11 = isa_atmosphere(11000.0)
    print(f"ISA check: 11 km T = {T11:.2f} K (216.65), p = {p11:.0f} Pa (22632), rho = {rho11:.4f} kg/m^3 (0.3639)")

    # --- a 2 km table, evaluated at 7.5 km and across the whole range
    table = isa_table(2000.0)
    h_fine = np.linspace(0.0, 86000.0, 8601)
    _, _, rho_true = isa_atmosphere(h_fine)
    print(f"\n2 km table ({table.shape[0]} rows): density at 7.5 km, truth {isa_atmosphere(7500.0)[2]:.6f} kg/m^3")
    print(f"{'method':<12}{'rho(7.5 km)':>14}{'rel. error':>12}{'max rel. error 0-86 km':>26}{'max rel. error 0-11 km':>26}")
    for method in ("linear", "loglinear", "natural", "not-a-knot"):
        r75 = interpolate_isa(table, 7500.0, "rho", method)
        ri = interpolate_isa(table, h_fine, "rho", method)
        rel = np.abs(ri - rho_true) / rho_true
        print(f"{method:<12}{r75:>14.6f}{abs(r75 - isa_atmosphere(7500.0)[2]) / isa_atmosphere(7500.0)[2]:>12.2e}"
              f"{np.max(rel):>26.2e}{np.max(rel[h_fine <= 11000.0]):>26.2e}")
    cs = si.CubicSpline(table[:, 0], table[:, 3], bc_type="natural")
    ours = spline_evaluate(cubic_spline(table[:, 0], table[:, 3], "natural"), h_fine)
    print(f"natural spline against scipy.interpolate.CubicSpline: max |diff| = {np.max(np.abs(ours - cs(h_fine))):.2e}")

    # --- convergence study in the troposphere (smooth: no layer boundary)
    f_rho = lambda h: isa_atmosphere(h)[2]                   # noqa: E731
    steps = np.array([2000.0, 1000.0, 500.0, 250.0])
    print("\nconvergence of the density interpolant on 0-11 km (max abs error, kg/m^3):")
    print(f"{'step m':>8}" + "".join(f"{m:>14}" for m in ("linear", "loglinear", "natural", "not-a-knot", "clamped")))
    errs = {m: [] for m in ("linear", "loglinear", "natural", "not-a-knot", "clamped")}
    for hstep in steps:
        row = f"{hstep:>8.0f}"
        for m in errs:
            e = max_interpolation_error(f_rho, 0.0, 11000.0, int(round(11000.0 / hstep)), m)
            errs[m].append(e)
            row += f"{e:>14.2e}"
        print(row)
    print("observed orders: " + ", ".join(f"{m} {observed_order(steps, errs[m]):.2f}" for m in errs))
    print("theory: linear 2, loglinear 2, natural 2 (end layers) / 4 interior, not-a-knot 4, clamped 4")

    # --- the tropopause kink
    table1 = isa_table(1000.0)
    hk = np.linspace(9000.0, 13000.0, 401)
    T_true = isa_atmosphere(hk)[0]
    T_sp = interpolate_isa(table1, hk, "T", "natural")
    T_li = interpolate_isa(table1, hk, "T", "linear")
    print(f"\ntemperature across the tropopause from a 1 km table: max error spline {np.max(np.abs(T_sp - T_true)):.3f} K, "
          f"linear {np.max(np.abs(T_li - T_true)):.1e} K (the table nodes include 11 km, and T is piecewise linear)")

    # --- Runge
    xf = np.linspace(-1.0, 1.0, 1001)
    print("\nRunge's function 1/(1+25x^2) on [-1, 1]: max error of the interpolating polynomial")
    print(f"{'n nodes':>8}{'equispaced':>14}{'Chebyshev':>14}")
    for n in (5, 9, 13, 17, 21):
        xe = np.linspace(-1.0, 1.0, n)
        xc = chebyshev_nodes(n, -1.0, 1.0)
        ee = np.max(np.abs(lagrange_interpolate(xe, runge(xe), xf) - runge(xf)))
        xcf = xf[(xf >= xc[0]) & (xf <= xc[-1])]
        ec = np.max(np.abs(lagrange_interpolate(xc, runge(xc), xcf) - runge(xcf)))
        print(f"{n:>8}{ee:>14.3e}{ec:>14.3e}")

    # --- figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    ax = axes[0]
    for method in ("linear", "loglinear", "natural"):
        ri = interpolate_isa(table, h_fine, "rho", method)
        ax.semilogy(h_fine / 1e3, np.maximum(np.abs(ri - rho_true) / rho_true, 1e-12), label=method)
    for hb in ISA_BASES[1:-1]:
        ax.axvline(hb / 1e3, color="gray", ls=":", lw=0.8)
    ax.set_xlabel("altitude [km]")
    ax.set_ylabel("relative density error")
    ax.set_title("2 km ISA table: interpolation error")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    ax = axes[1]
    for m in ("linear", "natural", "not-a-knot", "clamped"):
        ax.loglog(steps, errs[m], "o-", label=f"{m} (p = {observed_order(steps, errs[m]):.2f})")
    ax.set_xlabel("table spacing [m]")
    ax.set_ylabel("max |error| in density, 0-11 km")
    ax.set_title("Convergence with table spacing")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    ax = axes[2]
    xe = np.linspace(-1.0, 1.0, 13)
    xc = chebyshev_nodes(13, -1.0, 1.0)
    ax.plot(xf, runge(xf), "k", label="1/(1+25x²)")
    ax.plot(xf, lagrange_interpolate(xe, runge(xe), xf), label="13 equispaced nodes")
    xcf = xf[(xf >= xc[0]) & (xf <= xc[-1])]
    ax.plot(xcf, lagrange_interpolate(xc, runge(xc), xcf), label="13 Chebyshev nodes")
    ax.set_ylim(-1.5, 2.0)
    ax.set_title("Runge phenomenon")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "atmosphere_interpolation.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
