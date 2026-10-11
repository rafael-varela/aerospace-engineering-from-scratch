"""P00.14 Measurement, Calibration & Uncertainty: reference solution.

An uncertainty toolkit following the GUM (JCGM 100:2008) and NIST TN 1297:

* Type A standard uncertainty from repeated readings,
* Type B standard uncertainty from a resolution or a stated half-width,
* combined standard uncertainty through any function by a finite-difference Jacobian
  (the law of propagation of uncertainty),
* expanded uncertainty with a coverage factor k,
* rounding and reporting to the right number of significant figures,
* a straight-line calibration fit (slope, offset, residuals and their uncertainties).

Run:  python solution.py [output.png]
Test: AEFS_IMPL=solution pytest curriculum/phase-00/lesson-14-measurement-calibration-uncertainty
"""
from __future__ import annotations

import math
import sys
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ArrayLike = Sequence[float] | np.ndarray

_DIVISORS: dict[str, float] = {"rectangular": math.sqrt(3.0), "triangular": math.sqrt(6.0), "normal": 2.0}
"""Half-width divisors: a/sqrt(3) for a uniform distribution, a/sqrt(6) for a triangular one,
and a/2 when the half-width is a stated 95 % (k = 2) bound."""


def type_a(readings: ArrayLike) -> tuple[float, float]:
    """Mean and Type A standard uncertainty of the mean from repeated readings.

    Returns ``(mean, s / sqrt(n))`` with ``s`` the sample standard deviation (n - 1 in the
    denominator). Requires at least two readings.
    """
    x = np.asarray(readings, dtype=float).ravel()
    n = x.size
    if n < 2:
        raise ValueError("Type A evaluation needs at least two readings")
    return float(x.mean()), float(x.std(ddof=1) / math.sqrt(n))


def type_b_half_width(a: float, distribution: str = "rectangular") -> float:
    """Type B standard uncertainty from a half-width ``a``.

    ``distribution`` is ``"rectangular"`` (any value in ±a equally likely: a/sqrt(3)),
    ``"triangular"`` (values near the centre more likely: a/sqrt(6)), or ``"normal"``
    (``a`` is a quoted 95 % bound, k = 2: a/2).
    """
    if a < 0:
        raise ValueError("half-width must be non-negative")
    try:
        return a / _DIVISORS[distribution]
    except KeyError as exc:
        raise ValueError(f"unknown distribution {distribution!r}") from exc


def type_b_resolution(resolution: float) -> float:
    """Type B standard uncertainty of a digital display with the given step (resolution).

    The true value lies anywhere within ±resolution/2 of the displayed value with equal
    probability, so u = (resolution / 2) / sqrt(3) = resolution / sqrt(12).
    """
    return type_b_half_width(resolution / 2.0, "rectangular")


def combine(*components: float) -> float:
    """Root-sum-square of independent standard uncertainties."""
    return math.sqrt(sum(float(c) ** 2 for c in components))


def jacobian_fd(f: Callable[[np.ndarray], float], x: ArrayLike, rel_step: float = 1e-6) -> np.ndarray:
    """Partial derivatives of scalar ``f`` at ``x`` by central finite differences.

    The step in component ``i`` is ``rel_step * max(|x_i|, 1)``. Returns an array of the
    same length as ``x``. Central differences are exact for functions that are linear in
    each input and second-order accurate otherwise.
    """
    x0 = np.asarray(x, dtype=float).ravel()
    jac = np.empty_like(x0)
    for i in range(x0.size):
        h = rel_step * max(abs(x0[i]), 1.0)
        xp, xm = x0.copy(), x0.copy()
        xp[i] += h
        xm[i] -= h
        jac[i] = (float(f(xp)) - float(f(xm))) / (2.0 * h)
    return jac


@dataclass(frozen=True)
class Propagated:
    """Result of propagating uncertainty through a function."""

    value: float
    u: float
    sensitivities: np.ndarray
    """The partial derivatives dy/dx_i evaluated at the inputs (the sensitivity coefficients)."""
    contributions: np.ndarray
    """Each input's contribution |c_i| * u_i to the combined uncertainty, in units of y."""


def propagate(
    f: Callable[[np.ndarray], float],
    x: ArrayLike,
    u: ArrayLike,
    corr: np.ndarray | None = None,
) -> Propagated:
    """Combined standard uncertainty of ``y = f(x)`` by the law of propagation of uncertainty.

    ``x`` are the input estimates, ``u`` their standard uncertainties. Without ``corr`` the
    inputs are independent and ``u_y^2 = sum_i (c_i u_i)^2`` with ``c_i = df/dx_i``. With a
    correlation matrix ``corr`` (ones on the diagonal), ``u_y^2 = (c*u)^T corr (c*u)``.
    """
    x0 = np.asarray(x, dtype=float).ravel()
    u0 = np.asarray(u, dtype=float).ravel()
    if x0.shape != u0.shape:
        raise ValueError("x and u must have the same length")
    if np.any(u0 < 0):
        raise ValueError("standard uncertainties must be non-negative")
    c = jacobian_fd(f, x0)
    cu = c * u0
    if corr is None:
        var = float(cu @ cu)
    else:
        r = np.asarray(corr, dtype=float)
        var = float(cu @ r @ cu)
    return Propagated(value=float(f(x0)), u=math.sqrt(max(var, 0.0)), sensitivities=c, contributions=np.abs(cu))


def expanded(u: float, k: float = 2.0) -> float:
    """Expanded uncertainty U = k * u. k = 2 is about 95 % coverage for a normal distribution."""
    return k * u


def round_to_uncertainty(value: float, u: float, sig: int = 2) -> tuple[float, float, int]:
    """Round ``u`` to ``sig`` significant figures and ``value`` to the same decimal place.

    Returns ``(rounded_value, rounded_u, decimals)`` where ``decimals`` is the number of
    digits after the decimal point (negative means rounding to tens, hundreds, ...).
    """
    if u <= 0 or not math.isfinite(u):
        raise ValueError("uncertainty must be positive and finite")
    decimals = sig - 1 - math.floor(math.log10(u))
    u_r = round(u, decimals)
    # rounding can carry u to the next decade (0.0996 -> 0.10): recompute once
    if u_r > 0 and math.floor(math.log10(u_r)) != math.floor(math.log10(u)):
        decimals = sig - 1 - math.floor(math.log10(u_r))
        u_r = round(u, decimals)
    return round(value, decimals), u_r, decimals


def format_measurement(value: float, u: float, unit: str = "", sig: int = 2) -> str:
    """Report ``value ± u unit`` with ``u`` at ``sig`` significant figures.

    Examples: ``format_measurement(9.8123, 0.1234, "m/s^2")`` is ``"9.81 ± 0.12 m/s^2"``.
    """
    v_r, u_r, decimals = round_to_uncertainty(value, u, sig)
    if decimals > 0:
        body = f"{v_r:.{decimals}f} ± {u_r:.{decimals}f}"
    else:
        body = f"{int(round(v_r))} ± {int(round(u_r))}"
    return f"{body} {unit}".rstrip()


@dataclass(frozen=True)
class CalibrationFit:
    """Straight-line calibration ``reading = slope * reference + offset``."""

    slope: float
    offset: float
    u_slope: float
    u_offset: float
    residuals: np.ndarray
    """``reading - (slope * reference + offset)`` for each point, in reading units."""
    s_residual: float
    """Standard deviation of the residuals about the line (n - 2 degrees of freedom)."""


def calibration_fit(reference: ArrayLike, reading: ArrayLike) -> CalibrationFit:
    """Least-squares straight line through (reference, reading) pairs.

    Needs at least three points so that the residual scatter can be estimated.
    The standard uncertainties of slope and offset follow the usual ordinary-least-squares
    formulas with the residual standard deviation as the per-point uncertainty.
    """
    xr = np.asarray(reference, dtype=float).ravel()
    yr = np.asarray(reading, dtype=float).ravel()
    n = xr.size
    if n < 3 or yr.size != n:
        raise ValueError("need at least three (reference, reading) pairs of equal length")
    xbar, ybar = xr.mean(), yr.mean()
    sxx = float(np.sum((xr - xbar) ** 2))
    if sxx == 0.0:
        raise ValueError("reference values must not all be equal")
    slope = float(np.sum((xr - xbar) * (yr - ybar)) / sxx)
    offset = float(ybar - slope * xbar)
    residuals = yr - (slope * xr + offset)
    s = math.sqrt(float(np.sum(residuals**2)) / (n - 2))
    u_slope = s / math.sqrt(sxx)
    u_offset = s * math.sqrt(1.0 / n + xbar**2 / sxx)
    return CalibrationFit(slope, offset, u_slope, u_offset, residuals, s)


def apply_calibration(reading: ArrayLike | float, fit: CalibrationFit) -> np.ndarray | float:
    """Convert raw reading(s) to corrected value(s): ``(reading - offset) / slope``."""
    r = np.asarray(reading, dtype=float)
    out = (r - fit.offset) / fit.slope
    return float(out) if out.ndim == 0 else out


def corrected_value(reading: float, u_reading: float, fit: CalibrationFit) -> Propagated:
    """Corrected value and its uncertainty, propagating reading, slope and offset uncertainties.

    Slope and offset are treated as independent, which slightly overstates the uncertainty
    away from the centre of the calibration range; good enough for a first report.
    """
    return propagate(lambda p: (p[0] - p[2]) / p[1], [reading, fit.slope, fit.offset],
                     [u_reading, fit.u_slope, fit.u_offset])


# Illustrative example numbers for main(): they stand in for your own data/calibration.csv and are
# NOT a measurement. Reference masses are nominal coin masses; readings are made-up scale displays.
EXAMPLE_REFERENCE_G = np.array([0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0])
EXAMPLE_READING_G = np.array([0.0, 4.9, 9.9, 15.0, 19.8, 24.8, 29.7])
EXAMPLE_RESOLUTION_G = 0.1


def main() -> None:
    """Run the toolkit on the illustrative example and save a residual plot."""
    out = sys.argv[1] if len(sys.argv) > 1 else "calibration_residuals.png"
    print("Illustrative example (replace with your own data from projects/p00-calibration/)")

    repeats = np.array([9.9, 10.0, 9.9, 10.0, 9.8])
    mean, u_a = type_a(repeats)
    u_b = type_b_resolution(EXAMPLE_RESOLUTION_G)
    u_c = combine(u_a, u_b)
    print(f"repeated readings of a 10 g reference: mean {mean:.3f} g, u_A {u_a:.3f} g, "
          f"u_B(resolution) {u_b:.3f} g, combined {u_c:.3f} g, U(k=2) {expanded(u_c):.3f} g")

    fit = calibration_fit(EXAMPLE_REFERENCE_G, EXAMPLE_READING_G)
    print(f"calibration: slope {format_measurement(fit.slope, fit.u_slope)}, "
          f"offset {format_measurement(fit.offset, fit.u_offset, 'g')}, "
          f"residual s {fit.s_residual:.3f} g")

    raw = 17.3
    corr = corrected_value(raw, u_c, fit)
    print(f"raw reading {raw} g -> corrected {format_measurement(corr.value, corr.u, 'g')} "
          f"(contributions: reading {corr.contributions[0]:.3f}, slope {corr.contributions[1]:.3f}, "
          f"offset {corr.contributions[2]:.3f} g)")

    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.axhline(0.0, color="0.5", lw=0.8)
    ax.errorbar(EXAMPLE_REFERENCE_G, fit.residuals, yerr=u_c, fmt="o", capsize=3)
    ax.set_xlabel("reference mass (g)")
    ax.set_ylabel("reading minus fit (g)")
    ax.set_title("Calibration residuals (illustrative example)")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
