"""P03.01 Floating Point & Numerical Error: starter code.

Catastrophic cancellation demo: every function here either exposes how IEEE 754 doubles are
stored and rounded, or computes the same quantity two ways (a textbook formula that cancels and
a rearranged one that does not) so the loss of digits can be measured, not guessed.

Complete every function marked NotImplementedError; ``main()`` already produces the lesson output once they work.

Conventions:
    * ``eps`` is machine epsilon, the spacing between 1.0 and the next double, 2^-52 = 2.22e-16.
      The unit roundoff u = eps / 2 bounds the relative error of one rounded operation:
      fl(x op y) = (x op y)(1 + delta) with |delta| <= u.
    * ``ulp(x)`` is the spacing of doubles at x: 2^(e - 52) where e = floor(log2 |x|).
    * Relative error is |approx - exact| / |exact|; "correct digits" is -log10 of that.
    * Every "cancelling" and "stable" pair computes the same mathematical quantity; the stable
      form is algebraically identical and differs only in which rounding errors it exposes.
    * Any dtype argument accepts ``numpy.float32`` or ``numpy.float64``; float32 is used only to
      make effects visible that float64 hides in the 9th-16th digit.

Run:  python main.py          (prints the cancellation tables, writes floating_point_error.png)
Test: pytest curriculum/phase-03/lesson-01-floating-point-numerical-error
"""
from __future__ import annotations

import math
import struct
from fractions import Fraction
from pathlib import Path
from typing import Callable

import numpy as np

EPS64: float = 2.0 ** -52      # machine epsilon of a double
EPS32: float = 2.0 ** -23      # machine epsilon of a single


# ---------------------------------------------------------------- the representation

def float_components(x: float) -> tuple[int, int, int]:
    """Return (sign, biased_exponent, fraction) of the IEEE 754 double ``x``.

    sign is 0 or 1, biased_exponent is the 11-bit field (0..2047; the true exponent is this
    minus 1023 for normal numbers), fraction is the 52-bit field (the mantissa without its
    implicit leading 1). Decoded by packing the double as big-endian bytes.
    """
    raise NotImplementedError


def machine_epsilon(dtype: type = np.float64) -> float:
    """Return machine epsilon for ``dtype`` by halving: the smallest power of two e with 1 + e != 1.

    Returns 2^-52 for float64 and 2^-23 for float32, equal to ``numpy.finfo(dtype).eps``.
    Raise ``TypeError`` for a non-floating dtype.
    """
    raise NotImplementedError


def ulp(x: float) -> float:
    """Return the spacing of doubles at ``x``: 2^(e - 52) with e = floor(log2 |x|), as ``numpy.spacing``.

    For x = 0 return the smallest subnormal, 2^-1074. Raise ``ValueError`` for inf or nan.
    """
    raise NotImplementedError


def relative_error(approx: np.ndarray | float, exact: np.ndarray | float) -> np.ndarray | float:
    """Return |approx - exact| / |exact| elementwise. Raise ``ValueError`` if any exact value is 0."""
    raise NotImplementedError


def correct_digits(approx: np.ndarray | float, exact: np.ndarray | float) -> np.ndarray | float:
    """Return -log10(relative error), capped at 16 when the relative error is zero."""
    raise NotImplementedError


def exact_rounding_error(a: float, b: float, op: str) -> float:
    """Return delta such that fl(a op b) = (a op b)(1 + delta), computed exactly with ``Fraction``.

    ``op`` is '+', '-', '*' or '/'. This checks the rounding model: |delta| <= u = eps / 2 always.
    Raise ``ValueError`` for an unknown op or when the exact result is zero.
    """
    raise NotImplementedError


def fixed_point_truncation_error(value: float, bits: int) -> float:
    """Return value - trunc(value) where trunc keeps ``bits`` binary digits after the point.

    This is the per-tick error of a clock that stores a decimal tick length (such as 0.1 s) in a
    binary fixed-point register. Exact arithmetic via ``Fraction``. Raise ``ValueError`` if bits < 1.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- cancellation, two ways each

def one_minus_cos_naive(theta: np.ndarray | float) -> np.ndarray | float:
    """Return 1 - cos(theta) exactly as written: cancels for small theta."""
    raise NotImplementedError


def one_minus_cos_stable(theta: np.ndarray | float) -> np.ndarray | float:
    """Return 1 - cos(theta) as 2 sin^2(theta / 2): no cancellation, relative error ~ eps for all theta."""
    raise NotImplementedError


def cw_position_coefficient(n: float, t: float, stable: bool = True) -> float:
    """Return 2 (1 - cos(n t)) / n [s], the Clohessy-Wiltshire coefficient of the along-track velocity in x(t).

    From P02.05: x(t) = (4 - 3 cos nt) x0 + (sin nt / n) xdot0 + (2 (1 - cos nt) / n) ydot0. For a 1 s
    step at n ~ 1.1e-3 rad/s the argument is 1e-3 rad and the naive form loses about six digits.
    """
    raise NotImplementedError


def quadratic_roots_naive(a: float, b: float, c: float) -> tuple[float, float]:
    """Return the two real roots of a x^2 + b x + c = 0 by the textbook formula (-b -+ sqrt(D)) / 2a.

    Returned as (smaller magnitude, larger magnitude). Raise ``ValueError`` if a == 0 or D < 0.
    One root cancels when |b| >> |a c|.
    """
    raise NotImplementedError


def quadratic_roots_stable(a: float, b: float, c: float) -> tuple[float, float]:
    """Return the roots with no cancellation: q = -(b + sign(b) sqrt(D)) / 2, x1 = q / a, x2 = c / q.

    The large root is formed by adding like signs; the small one by Vieta's product x1 x2 = c / a.
    Same ordering and errors as ``quadratic_roots_naive``.
    """
    raise NotImplementedError


def exp_series(x: float, max_terms: int = 200) -> float:
    """Return e^x by summing the Taylor series until a term is below eps times the running sum.

    Accurate for x >= 0. For x << 0 the terms alternate and grow to e^|x| before shrinking, so the
    sum of large cancelling terms carries an absolute error ~ eps e^|x| that swamps the tiny answer.
    """
    raise NotImplementedError


def exp_stable(x: float) -> float:
    """Return e^x by the series for |x|, inverting when x < 0: e^x = 1 / e^|x|. Relative error ~ eps."""
    raise NotImplementedError


def subtract_positions(r_station: np.ndarray, separation: np.ndarray, dtype: type = np.float64) -> np.ndarray:
    """Return the relative position recovered by storing both absolute positions in ``dtype`` and subtracting.

    r_chaser = r_station + separation is rounded to ``dtype`` first, as a navigation filter that
    stores absolute Earth-centred positions would do. The error against ``separation`` is the ulp
    of the absolute position, not of the separation.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- summation and variance

def naive_sum(x: np.ndarray, dtype: type = np.float64) -> float:
    """Return the running left-to-right sum of ``x`` with every partial sum rounded to ``dtype``.

    Worst-case error ~ n u sum|x|; for random rounding errors the typical growth is ~ sqrt(n).
    """
    raise NotImplementedError


def kahan_sum(x: np.ndarray, dtype: type = np.float64) -> float:
    """Return the compensated (Kahan) sum of ``x`` in ``dtype``: the rounding error of each add is
    carried to the next. Error bound ~ 2 u sum|x| independent of n."""
    raise NotImplementedError


def pairwise_sum(x: np.ndarray, dtype: type = np.float64) -> float:
    """Return the sum of ``x`` by recursive halving (pairwise summation). Error ~ u log2(n) sum|x|.

    This is the strategy inside ``numpy.sum`` (blocked), which is why ``numpy.sum`` beats a loop.
    """
    raise NotImplementedError


def variance_one_pass(x: np.ndarray, dtype: type = np.float64) -> float:
    """Return the population variance by the textbook one-pass formula E[x^2] - E[x]^2 in ``dtype``.

    Cancels catastrophically when the standard deviation is small compared with the mean (an
    altimeter reading 400 km with centimetre scatter). Can even return a negative number.
    """
    raise NotImplementedError


def variance_two_pass(x: np.ndarray, dtype: type = np.float64) -> float:
    """Return the population variance by the two-pass formula mean((x - mean)^2) in ``dtype``.

    Shifting by the mean first removes the cancellation; the relative error is ~ eps.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- truncation versus rounding

def forward_difference(f: Callable[[float], float], x: float, h: float) -> float:
    """Return (f(x + h) - f(x)) / h, the first-order forward-difference estimate of f'(x)."""
    raise NotImplementedError


def fd_total_error_bound(h: np.ndarray | float, f2_bound: float, f_bound: float, eps: float = EPS64) -> np.ndarray | float:
    """Return the error model of the forward difference: h f2_bound / 2 + 2 eps f_bound / h.

    The first term is the Taylor truncation error (|f''| <= f2_bound on the interval); the second
    is the rounding error of the two function values and of x + h, divided by h, where f_bound
    bounds both |f| and |x f'| near x. The sum is minimised at h* = 2 sqrt(eps f_bound / f2_bound).
    """
    raise NotImplementedError


def fd_optimal_step(f2_bound: float, f_bound: float, eps: float = EPS64) -> float:
    """Return h* = 2 sqrt(eps f_bound / f2_bound), the step that minimises ``fd_total_error_bound``.

    For f ~ f'' ~ 1 this is about 3e-8, and the best achievable error is about 3e-8 as well: the
    forward difference can never deliver more than half the digits of the machine.
    """
    raise NotImplementedError


def relative_condition_number(f: Callable[[float], float], dfdx: Callable[[float], float], x: float) -> float:
    """Return kappa = |x f'(x) / f(x)|, the factor by which a relative error in x is amplified in f(x).

    kappa ~ 1 is well conditioned (sqrt has 1/2, exp has |x|); kappa >> 1 means the *problem*, not
    the algorithm, loses digits: no rearrangement can recover them. Raise ``ValueError`` if f(x) = 0.
    """
    raise NotImplementedError


def observed_order(hs: np.ndarray, errs: np.ndarray) -> float:
    """Return the least-squares slope of log(err) against log(h): the observed order of an error law err ~ C h^p.

    Positive for truncation-dominated errors, negative (about -1) where rounding dominates a
    difference quotient. Raise ``ValueError`` for fewer than two points or a non-positive entry.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the representation facts and the cancellation tables; write floating_point_error.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    eps = machine_epsilon()
    print(f"machine epsilon (double) {eps:.3e} = 2^-52; (single) {machine_epsilon(np.float32):.3e} = 2^-23")
    s, e, f = float_components(0.1)
    print(f"0.1 stored as sign={s} exponent={e} (2^{e - 1023}) fraction=0x{f:013x}; exactly {Fraction(0.1)}")
    print(f"ulp(1.0) = {ulp(1.0):.3e}, ulp(6.778e6 m) = {ulp(6.778e6):.3e} m, "
          f"float32 ulp at 6.778e6 m = {float(np.spacing(np.float32(6.778e6))):.3f} m")

    # --- separation of two spacecraft from absolute positions
    r_st = np.array([6.778e6, 0.0, 0.0])
    sep = np.array([-50.0, -1000.0, 0.0]) + np.array([0.0123, 0.0456, 0.0])
    for dtype in (np.float64, np.float32):
        rec = subtract_positions(r_st, sep, dtype)
        print(f"separation recovered in {dtype.__name__}: error {np.max(np.abs(rec - sep)):.3e} m")

    # --- 1 - cos theta
    thetas = 10.0 ** -np.arange(1, 9)
    exact = 2.0 * np.sin(thetas / 2.0) ** 2                  # the stable form is the reference (error ~ eps)
    print("\n   theta     naive rel. error   stable rel. error   digits lost (naive)")
    for th, ex in zip(thetas, exact):
        rn = relative_error(one_minus_cos_naive(th), ex)
        rs = abs(one_minus_cos_stable(th) - ex) / ex
        print(f"  {th:8.0e}   {rn:16.3e}   {rs:17.3e}   {max(0.0, math.log10(max(rn, 1e-17) / EPS64)):6.1f}")
    n_cw = 1.1314e-3
    print(f"CW coefficient 2(1 - cos nt)/n at t = 1 s: naive {cw_position_coefficient(n_cw, 1.0, False):.12e}, "
          f"stable {cw_position_coefficient(n_cw, 1.0, True):.12e} s")

    # --- quadratic roots: closing time with tiny thrust
    a, v, s0 = 0.5e-6, 1.0, 10.0                              # ½ a t² + v t - s = 0
    tn = quadratic_roots_naive(a, v, -s0)[0]
    ts = quadratic_roots_stable(a, v, -s0)[0]
    print(f"\nclosing time: naive {tn:.12f} s, stable {ts:.12f} s, residual naive {a*tn*tn + v*tn - s0:.2e}, stable {a*ts*ts + v*ts - s0:.2e}")
    rn, rs = quadratic_roots_naive(1.0, 1e8, 1.0)[0], quadratic_roots_stable(1.0, 1e8, 1.0)[0]
    print(f"x^2 + 1e8 x + 1 = 0 small root: naive {rn:.6e}, stable {rs:.6e} (exact -1.0000000000000001e-08)")

    # --- e^x by series
    print("\n     x    series       stable       exact        rel.err series")
    for x in (1.0, -5.0, -20.0, -30.0):
        print(f"  {x:5.0f}  {exp_series(x):.6e} {exp_stable(x):.6e} {math.exp(x):.6e}  {relative_error(exp_series(x), math.exp(x)):.2e}")

    # --- summation and variance in float32
    rng = np.random.default_rng(301)
    x = rng.uniform(0.0, 1.0, 1_000_000).astype(np.float32)
    exact_sum = math.fsum(x.astype(np.float64).tolist())
    print(f"\nsum of 1e6 uniforms in float32: naive err {abs(naive_sum(x, np.float32) - exact_sum):.3f}, "
          f"pairwise err {abs(pairwise_sum(x, np.float32) - exact_sum):.4f}, kahan err {abs(kahan_sum(x, np.float32) - exact_sum):.5f}  (exact {exact_sum:.3f})")
    alt = 400_000.0 + 0.01 * rng.standard_normal(10_000)
    v_exact = float(np.var(alt))
    print(f"variance of altimeter readings (mean 400 km, sigma 1 cm): one-pass float64 {variance_one_pass(alt):.4e}, "
          f"two-pass {variance_two_pass(alt):.4e}, numpy {v_exact:.4e}; one-pass float32 {variance_one_pass(alt, np.float32):.4e}")

    # --- forward difference: truncation versus rounding
    hs = 10.0 ** -np.arange(1, 17, dtype=float)
    errs = np.array([abs(forward_difference(math.sin, 1.0, h) - math.cos(1.0)) for h in hs])
    bound = fd_total_error_bound(hs, 1.0, 1.0)
    h_star = fd_optimal_step(1.0, 1.0)
    i_best = int(np.argmin(errs))
    print(f"\nforward difference of sin at 1: best h = 1e-{i_best + 1} with error {errs[i_best]:.2e}; model h* = {h_star:.2e}, "
          f"model error there {fd_total_error_bound(h_star, 1.0, 1.0):.2e}")
    print(f"observed order, truncation regime (h = 1e-2..1e-4): {observed_order(hs[1:4], errs[1:4]):.3f}; "
          f"rounding regime (h = 1e-10..1e-15): {observed_order(hs[9:15], errs[9:15]):.3f}")

    # --- fixed point clock
    per_tick = fixed_point_truncation_error(0.1, 24)
    drift = per_tick * 10 * 3600 * 100
    print(f"\n0.1 s tick in 24-bit fixed point: error {per_tick:.3e} s per tick, {drift:.3f} s after 100 h, "
          f"{drift * 1676:.0f} m at 1676 m/s")

    # --- figure
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    ax.loglog(thetas, [relative_error(one_minus_cos_naive(t), e_) for t, e_ in zip(thetas, exact)], "o-", label="1 - cos θ (naive)")
    ax.loglog(thetas, np.maximum([abs(one_minus_cos_stable(t) - e_) / e_ for t, e_ in zip(thetas, exact)], 1e-17), "s-", label="2 sin²(θ/2) (stable)")
    ax.loglog(thetas, EPS64 / (thetas ** 2 / 2), "--", color="gray", label="eps / (θ²/2) model")
    ax.set_xlabel("θ [rad]")
    ax.set_ylabel("relative error")
    ax.set_title("Cancellation in 1 - cos θ")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    ax = axes[1]
    ax.loglog(hs, errs, "o-", label="|FD - cos 1| (observed)")
    ax.loglog(hs, bound, "--", label="h/2 + 2 eps/h (model)")
    ax.axvline(h_star, color="gray", ls=":", label=f"h* = {h_star:.1e}")
    ax.set_xlabel("step h")
    ax.set_ylabel("absolute error")
    ax.set_title("Forward difference: truncation vs rounding")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "floating_point_error.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
