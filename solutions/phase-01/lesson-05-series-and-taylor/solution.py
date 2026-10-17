"""P01.05 Sequences, Series & Taylor Expansions: reference solution.

Partial sums of geometric, Taylor and binomial series, the remainder bounds that say how far a
partial sum is from the limit, and the error of the small-angle approximations used throughout
flight dynamics (sin a ~ a, cos a ~ 1, tan a ~ a), the Prandtl-Glauert compressibility factor
1/sqrt(1 - M^2) as a binomial series, and the period of a large-amplitude pendulum (the swing test
that measures an aircraft's moment of inertia) as a power series in sin(theta0/2).

Run:  python solution.py          (prints the tables, writes small_angle_error.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-01/lesson-05-series-and-taylor
"""
from __future__ import annotations

from math import factorial
from pathlib import Path

import numpy as np
import sympy as sp

x_s = sp.Symbol("x", real=True)

APPROXIMATIONS = ("sin", "cos", "cos2", "tan")   # sin ~ x, cos ~ 1, cos ~ 1 - x^2/2, tan ~ x


def model_symbols() -> dict[str, sp.Symbol]:
    """Return the symbol used by the symbolic helpers: ``x`` (real)."""
    return {"x": x_s}


# ---------------------------------------------------------------- sequences and series

def partial_sums(terms: np.ndarray) -> np.ndarray:
    """Return the sequence of partial sums S_k = a_0 + ... + a_k of the 1-D array ``terms``."""
    return np.cumsum(np.asarray(terms, dtype=float))


def geometric_terms(r: float, n: int) -> np.ndarray:
    """Return the first ``n`` terms r^0, r^1, ..., r^(n-1) of the geometric series with ratio ``r``."""
    if n < 1:
        raise ValueError("n must be at least 1")
    return float(r) ** np.arange(n, dtype=float)


def geometric_sum(r: float) -> float:
    """Return the sum 1/(1 - r) of the geometric series, or raise ValueError when |r| >= 1 (divergent)."""
    if abs(r) >= 1.0:
        raise ValueError("the geometric series diverges for |r| >= 1")
    return 1.0 / (1.0 - r)


def harmonic_partial_sums(n: int) -> np.ndarray:
    """Return H_k = 1 + 1/2 + ... + 1/k for k = 1..n: the textbook divergent series that looks convergent."""
    if n < 1:
        raise ValueError("n must be at least 1")
    return np.cumsum(1.0 / np.arange(1, n + 1, dtype=float))


def ratio_test_estimate(terms: np.ndarray) -> float:
    """Return |a_{n}/a_{n-1}| for the last two non-zero entries of ``terms``: the ratio-test limit, estimated.

    Below 1 the series converges absolutely, above 1 it diverges, at 1 the test is silent.
    """
    t = np.abs(np.asarray(terms, dtype=float))
    t = t[t > 0.0]
    if t.size < 2:
        raise ValueError("need at least two non-zero terms")
    return float(t[-1] / t[-2])


def radius_of_convergence(coeffs: np.ndarray) -> float:
    """Return R estimated from the last two non-zero power-series coefficients, R = |c_j / c_k|^(1/(k-j)), j < k.

    ``coeffs[k]`` multiplies x^k. For a series in x^2 (every other coefficient zero) the exponent 1/2
    is what the ratio test gives. The estimate tends to the true radius as the series is extended.
    """
    c = np.abs(np.asarray(coeffs, dtype=float))
    nz = np.flatnonzero(c > 0.0)
    if nz.size < 2:
        raise ValueError("need at least two non-zero coefficients")
    j, k = nz[-2], nz[-1]
    return float((c[j] / c[k]) ** (1.0 / (k - j)))


def alternating_series_bound(terms: np.ndarray, k: int) -> float:
    """Return |a_{k+1}|, the alternating-series bound on |S - S_k| for an alternating series of decreasing terms."""
    t = np.asarray(terms, dtype=float)
    if k + 1 >= t.size:
        raise ValueError("no term after index k to bound the remainder with")
    return float(abs(t[k + 1]))


# ---------------------------------------------------------------- Taylor polynomials

def _centre(x0: float | sp.Expr) -> sp.Expr:
    """Turn a float centre into an exact SymPy number when it is a simple multiple of pi, else keep it as a Float."""
    if isinstance(x0, sp.Expr):
        return x0
    guess = sp.nsimplify(x0, [sp.pi], tolerance=1e-12, rational=False)
    return guess if abs(float(guess) - float(x0)) < 1e-12 and len(str(guess)) < 12 else sp.Float(x0)


def taylor_polynomial(expr: sp.Expr, x0: float | sp.Expr, n: int) -> sp.Expr:
    """Return the degree-``n`` Taylor polynomial of the SymPy expression ``expr`` in ``x`` about ``x0``.

    Built from derivatives, sum_{k=0}^{n} f^(k)(x0) (x - x0)^k / k!, left in powers of (x - x0).
    A float centre that is a simple multiple of pi (pi/3, pi/4, ...) is recognised and kept exact.
    """
    if n < 0:
        raise ValueError("degree must be non-negative")
    c = _centre(x0)
    poly = sp.Integer(0)
    d = expr
    for k in range(n + 1):
        poly += sp.simplify(d.subs(x_s, c)) * (x_s - c) ** k / sp.factorial(k)
        d = sp.diff(d, x_s)
    return poly


def taylor_coefficients(expr: sp.Expr, x0: float | sp.Expr, n: int) -> np.ndarray:
    """Return the numeric Taylor coefficients c_k = f^(k)(x0)/k!, k = 0..n, of a SymPy expression in ``x`` about ``x0``."""
    if n < 0:
        raise ValueError("degree must be non-negative")
    c = _centre(x0)
    out = np.empty(n + 1)
    d = expr
    for k in range(n + 1):
        out[k] = float(d.subs(x_s, c)) / factorial(k)
        d = sp.diff(d, x_s)
    return out


def series_terms(x: np.ndarray | float, kind: str, n_terms: int) -> np.ndarray:
    """Return the first ``n_terms`` Maclaurin terms of ``kind`` in {"sin", "cos", "exp"} at ``x``, shape (n_terms, *x.shape).

    Terms are the non-zero ones: sin -> x, -x^3/3!, ...; cos -> 1, -x^2/2!, ...; exp -> 1, x, x^2/2!, ...
    Each term is built from the previous by a ratio, which avoids overflow in the factorials.
    """
    if n_terms < 1:
        raise ValueError("n_terms must be at least 1")
    x = np.asarray(x, dtype=float)
    out = np.empty((n_terms,) + x.shape)
    if kind == "sin":
        term = x.copy()
        for k in range(n_terms):
            out[k] = term
            term = -term * x * x / ((2 * k + 2) * (2 * k + 3))
    elif kind == "cos":
        term = np.ones_like(x)
        for k in range(n_terms):
            out[k] = term
            term = -term * x * x / ((2 * k + 1) * (2 * k + 2))
    elif kind == "exp":
        term = np.ones_like(x)
        for k in range(n_terms):
            out[k] = term
            term = term * x / (k + 1)
    else:
        raise ValueError(f"unknown series kind {kind!r}")
    return out


def maclaurin_partial_sum(x: np.ndarray | float, kind: str, n_terms: int) -> np.ndarray | float:
    """Return the sum of the first ``n_terms`` non-zero Maclaurin terms of sin, cos or exp at ``x``."""
    return np.sum(series_terms(x, kind, n_terms), axis=0)


def lagrange_remainder_bound(m_next: float, x: float, x0: float, n: int) -> float:
    """Return M |x - x0|^(n+1) / (n+1)!, the Lagrange bound on |f(x) - T_n(x)| when |f^(n+1)| <= M between x0 and x."""
    if n < 0 or m_next < 0.0:
        raise ValueError("degree must be non-negative and the derivative bound non-negative")
    return float(m_next * abs(x - x0) ** (n + 1) / factorial(n + 1))


def terms_needed(x: float, kind: str, tol: float, max_terms: int = 200) -> int:
    """Return the smallest number of non-zero Maclaurin terms of ``kind`` whose remainder bound at ``x`` is below ``tol``.

    For sin and cos the series alternates with decreasing terms once |x|^2 < (2k+2)(2k+3), so the
    alternating bound |next term| applies; for exp the Lagrange bound e^|x| |x|^(n)/n! is used.
    """
    if tol <= 0.0:
        raise ValueError("tol must be positive")
    x = float(x)
    t = series_terms(x, kind, max_terms + 1)
    for k in range(1, max_terms + 1):
        nxt = abs(t[k])
        if kind == "exp":
            nxt = np.exp(abs(x)) * abs(x) ** k / factorial(k)
        elif abs(t[k - 1]) < nxt:
            continue                      # terms still growing: alternating bound not yet valid
        if nxt < tol:
            return k
    raise RuntimeError("remainder bound did not fall below tol within max_terms")


# ---------------------------------------------------------------- small-angle approximations

def small_angle_error(theta: np.ndarray | float, approx: str) -> np.ndarray | float:
    """Return the relative error |approximation - exact| / |exact| of a small-angle approximation at ``theta`` (rad).

    ``approx``: "sin" (sin ~ theta), "cos" (cos ~ 1), "cos2" (cos ~ 1 - theta^2/2), "tan" (tan ~ theta).
    """
    th = np.asarray(theta, dtype=float)
    if approx == "sin":
        exact, appr = np.sin(th), th
    elif approx == "cos":
        exact, appr = np.cos(th), np.ones_like(th)
    elif approx == "cos2":
        exact, appr = np.cos(th), 1.0 - th * th / 2.0
    elif approx == "tan":
        exact, appr = np.tan(th), th
    else:
        raise ValueError(f"unknown approximation {approx!r}; choose from {APPROXIMATIONS}")
    with np.errstate(divide="ignore", invalid="ignore"):
        err = np.abs(appr - exact) / np.abs(exact)
    return np.where(th == 0.0, 0.0, err) if err.ndim else (0.0 if th == 0.0 else float(err))


def small_angle_error_leading(theta: np.ndarray | float, approx: str) -> np.ndarray | float:
    """Return the leading-order relative error from the next Taylor term: theta^2/6, theta^2/2, theta^4/24, theta^2/3."""
    th = np.asarray(theta, dtype=float)
    lead = {"sin": th**2 / 6.0, "cos": th**2 / 2.0, "cos2": th**4 / 24.0, "tan": th**2 / 3.0}
    if approx not in lead:
        raise ValueError(f"unknown approximation {approx!r}; choose from {APPROXIMATIONS}")
    return lead[approx]


def max_angle_for_error(approx: str, tol: float) -> float:
    """Return the largest angle (rad, in (0, pi/2)) at which ``approx`` is still within relative error ``tol``, by bisection."""
    if not 0.0 < tol < 1.0:
        raise ValueError("tol must lie in (0, 1)")
    lo, hi = 0.0, 0.5 * np.pi - 1e-9
    if small_angle_error(hi, approx) < tol:
        return hi
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if small_angle_error(mid, approx) < tol:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------- binomial series: Prandtl-Glauert

def binomial_coefficients(alpha: float, n: int) -> np.ndarray:
    """Return the generalised binomial coefficients C(alpha, k) = alpha (alpha-1) ... (alpha-k+1) / k! for k = 0..n."""
    if n < 0:
        raise ValueError("n must be non-negative")
    c = np.empty(n + 1)
    c[0] = 1.0
    for k in range(1, n + 1):
        c[k] = c[k - 1] * (alpha - k + 1) / k
    return c


def prandtl_glauert(mach: np.ndarray | float) -> np.ndarray | float:
    """Return the Prandtl-Glauert factor 1/sqrt(1 - M^2) for subsonic Mach number ``mach`` in [0, 1)."""
    m = np.asarray(mach, dtype=float)
    if np.any(m < 0.0) or np.any(m >= 1.0):
        raise ValueError("Prandtl-Glauert is defined for 0 <= M < 1")
    return 1.0 / np.sqrt(1.0 - m * m)


def prandtl_glauert_series(mach: np.ndarray | float, n_terms: int) -> np.ndarray | float:
    """Return the binomial-series partial sum sum_{k<n_terms} C(-1/2, k) (-M^2)^k = 1 + M^2/2 + 3M^4/8 + ..."""
    m = np.asarray(mach, dtype=float)
    c = binomial_coefficients(-0.5, n_terms - 1)
    u = -m * m
    return sum(c[k] * u**k for k in range(n_terms))


# ---------------------------------------------------------------- pendulum period: the inertia swing test

def pendulum_period_ratio(theta0: np.ndarray | float, n_terms: int) -> np.ndarray | float:
    """Return T/T0 for a pendulum released from amplitude ``theta0`` (rad), as the series in k^2 = sin^2(theta0/2).

    T/T0 = sum_{n>=0} [ (2n)! / (2^(2n) n!^2) ]^2 k^(2n) = 1 + k^2/4 + 9 k^4/64 + ...  (T0 = 2 pi sqrt(L/g)).
    The small-angle approximation of P00 is the first term.
    """
    if n_terms < 1:
        raise ValueError("n_terms must be at least 1")
    k2 = np.sin(0.5 * np.asarray(theta0, dtype=float)) ** 2
    coef = binomial_coefficients(-0.5, n_terms - 1) ** 2        # [C(-1/2, n)]^2 = [(2n)!/(4^n n!^2)]^2
    return sum(coef[n] * k2**n for n in range(n_terms))


def pendulum_period_ratio_exact(theta0: float, n: int = 2000) -> float:
    """Return T/T0 = (2/pi) K(k), k = sin(theta0/2), by the midpoint rule on the complete elliptic integral."""
    k2 = np.sin(0.5 * float(theta0)) ** 2
    phi = (np.arange(n) + 0.5) * (0.5 * np.pi / n)
    return float((2.0 / np.pi) * np.sum(1.0 / np.sqrt(1.0 - k2 * np.sin(phi) ** 2)) * (0.5 * np.pi / n))


# ---------------------------------------------------------------- Euler's formula

def euler_formula_residual(theta: float, n_terms: int) -> float:
    """Return |exp-series(i theta) - (cos-series + i sin-series)| using ``n_terms`` of each real series and 2 n_terms of exp.

    With matching truncations the residual is rounding only: the exp series splits into cos (even powers)
    and sin (odd powers) term by term, which is Euler's formula.
    """
    z = complex(0.0, float(theta))
    term, total = 1.0 + 0.0j, 0.0 + 0.0j
    for k in range(2 * n_terms):
        total += term
        term = term * z / (k + 1)
    c = maclaurin_partial_sum(theta, "cos", n_terms)
    s = maclaurin_partial_sum(theta, "sin", n_terms)
    return float(abs(total - complex(c, s)))


# ---------------------------------------------------------------- driver

def main(out_dir: str | Path = ".") -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    print("geometric series, r = 0.5: partial sums", np.round(partial_sums(geometric_terms(0.5, 8)), 5), "limit", geometric_sum(0.5))
    print("harmonic series: H_10 = %.4f, H_1000 = %.4f, H_100000 = %.4f (no limit)" % tuple(harmonic_partial_sums(100000)[[9, 999, 99999]]))
    print("Taylor of sin about 0, degree 7:", taylor_polynomial(sp.sin(x_s), 0.0, 7))
    print("Taylor of cos about pi/3, degree 2:", taylor_polynomial(sp.cos(x_s), np.pi / 3, 2))
    print("radius of convergence estimated from 1/sqrt(1-x) coefficients:", radius_of_convergence(binomial_coefficients(-0.5, 60)), "(exact 1)")

    print("\nsmall-angle approximations: angle (deg) at which the relative error reaches 1 % and 0.1 %")
    for a, name in zip(APPROXIMATIONS, ("sin ~ a", "cos ~ 1", "cos ~ 1 - a^2/2", "tan ~ a")):
        print(f"  {name:18s} 1 %: {np.degrees(max_angle_for_error(a, 1e-2)):6.2f} deg   0.1 %: {np.degrees(max_angle_for_error(a, 1e-3)):6.2f} deg")
    for deg in (5.0, 10.0, 15.0, 30.0):
        th = np.radians(deg)
        print(f"  at {deg:4.0f} deg: sin {small_angle_error(th, 'sin'):.2e}  cos {small_angle_error(th, 'cos'):.2e}  "
              f"cos2 {small_angle_error(th, 'cos2'):.2e}  tan {small_angle_error(th, 'tan'):.2e}")
    print("terms of the sin series for 1e-10 at x = 0.5, 2, 10:", [terms_needed(x, "sin", 1e-10) for x in (0.5, 2.0, 10.0)])

    print("\nPrandtl-Glauert 1/sqrt(1 - M^2) against its binomial series")
    for m in (0.3, 0.6, 0.8, 0.95):
        ex = prandtl_glauert(m)
        print(f"  M = {m:4.2f}: exact {ex:.5f}  2 terms {prandtl_glauert_series(m, 2):.5f}  4 terms {prandtl_glauert_series(m, 4):.5f}  "
              f"10 terms {prandtl_glauert_series(m, 10):.5f}  ratio of last terms {binomial_coefficients(-0.5, 10)[10] / binomial_coefficients(-0.5, 10)[9] * (-m * m):.4f}")

    print("\npendulum swing test: T/T0 against amplitude")
    for deg in (5.0, 20.0, 45.0, 90.0):
        th = np.radians(deg)
        print(f"  {deg:3.0f} deg: 1 term 1.00000  2 terms {pendulum_period_ratio(th, 2):.5f}  3 terms {pendulum_period_ratio(th, 3):.5f}  exact {pendulum_period_ratio_exact(th):.5f}")
    print("Euler residual at theta = 1 with 10 terms:", euler_formula_residual(1.0, 10))

    th = np.radians(np.linspace(0.5, 60.0, 400))
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for a, name in zip(APPROXIMATIONS, ("sin a ~ a", "cos a ~ 1", "cos a ~ 1 - a^2/2", "tan a ~ a")):
        ax.loglog(np.degrees(th), small_angle_error(th, a), label=name)
        ax.loglog(np.degrees(th), small_angle_error_leading(th, a), "k:", lw=0.8)
    ax.axhline(1e-2, color="grey", lw=0.8)
    ax.set_xlabel("angle (deg)")
    ax.set_ylabel("relative error")
    ax.set_title("Small-angle approximations: error (solid) and leading Taylor term (dotted)")
    ax.legend()
    ax.grid(True, which="both", lw=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "small_angle_error.png", dpi=120)
    print(f"\nwrote {out_dir / 'small_angle_error.png'}")


if __name__ == "__main__":
    main()
