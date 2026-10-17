"""P02.18 Regression & Curve Fitting: starter code.

Thrust curve fit with confidence bands. Simple linear regression (slope, intercept, r, r^2,
standard errors, confidence and prediction bands, outlier flags, the significance test of r),
general linear least squares by QR with the parameter covariance (polynomials of any degree,
leave-one-out cross-validation to choose the degree), a damped Gauss-Newton solver for a
nonlinear thrust model with delta-method and Monte Carlo uncertainty on the total impulse, and the
log-transform exponential fit. ``scipy.stats.linregress``, ``numpy.linalg.lstsq`` and
``scipy.optimize.curve_fit`` are the oracles in the tests.

All thrust and pressure data in this module are SYNTHETIC: they are generated from the stated
models plus Gaussian noise with stated seeds, not measured on any motor.

Conventions:
    * x is treated as exact; y carries independent errors of equal variance (ordinary least squares).
    * s_e is the residual standard deviation sqrt(SS_res / (n - p)) with p fitted parameters.
    * Bands are returned as half-widths; add and subtract them from the fitted curve.
    * Thrust in N, time in s, chamber pressure in MPa, total impulse in N s.

Complete every function marked NotImplementedError (the synthetic-data generators and ``thrust_model`` are given).

Run:  python main.py          (prints the fit tables, writes thrust_curve_fit.png)
Test: pytest curriculum/phase-02/lesson-18-regression-curve-fitting
"""
from __future__ import annotations

import math
from collections.abc import Callable
from pathlib import Path

import numpy as np
from scipy import stats

# ---------------------------------------------------------------- synthetic data

def synthetic_pressure_thrust(seed: int = 18, n: int = 12, noise: float = 25.0,
                              with_outlier: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """Return SYNTHETIC (chamber pressure [MPa], thrust [N]) pairs from a static-fire series.

    Model: F = C_F A_t p_c with C_F = 1.45 and A_t = 2.0e-4 m^2 (a 16 mm throat), so the true slope is
    290 N/MPa and the true intercept 0; Gaussian noise of standard deviation ``noise`` is added to F.
    With ``with_outlier`` the 9th point is raised by 120 N (a data-logging fault).
    """
    rng = np.random.default_rng(seed)
    p = np.linspace(2.0, 6.0, n)
    F = 1.45 * 2.0e-4 * (p * 1e6) + rng.normal(0.0, noise, size=n)    # C_F A_t p_c [Pa] = 290 N/MPa x p
    if with_outlier:
        F = F.copy()
        F[8] += 120.0
    return p, F


def thrust_model(t: np.ndarray, params: np.ndarray) -> np.ndarray:
    """Return F(t) = a t^b exp(-c t): a three-parameter empirical thrust-curve shape (rise then tail-off).

    params = (a, b, c) with a > 0 [N s^-b], b > 0, c > 0 [1/s]. Peak at t = b/c; total impulse a Gamma(b+1)/c^(b+1).
    """
    a, b, c = params
    t = np.asarray(t, dtype=float)
    return a * np.power(np.clip(t, 0.0, None), b) * np.exp(-c * t)


def thrust_model_jacobian(t: np.ndarray, params: np.ndarray) -> np.ndarray:
    """Return the (n, 3) matrix of partial derivatives of ``thrust_model`` with respect to (a, b, c)."""
    raise NotImplementedError


TRUE_THRUST_PARAMS = np.array([1550.0, 0.8, 1.6])


def synthetic_thrust_curve(seed: int = 18, n: int = 60, noise: float = 10.0) -> tuple[np.ndarray, np.ndarray]:
    """Return SYNTHETIC (time [s], thrust [N]) samples of ``thrust_model`` with TRUE_THRUST_PARAMS plus N(0, noise^2)."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0.05, 4.0, n)
    return t, thrust_model(t, TRUE_THRUST_PARAMS) + rng.normal(0.0, noise, size=n)


# ---------------------------------------------------------------- simple linear regression

def linear_regression(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    """Fit y = intercept + slope x by least squares. Return a dict with
    'slope', 'intercept', 'r', 'r2', 's_e', 'se_slope', 'se_intercept', 'n', 'df', 'x_mean', 'sxx'.

    slope = S_xy / S_xx, intercept = y_mean - slope x_mean, r = S_xy / sqrt(S_xx S_yy), r^2 = 1 - SS_res/SS_tot,
    s_e = sqrt(SS_res / (n - 2)), se_slope = s_e / sqrt(S_xx), se_intercept = s_e sqrt(1/n + x_mean^2 / S_xx).
    """
    raise NotImplementedError


def predict_line(fit: dict[str, float], x_new: np.ndarray) -> np.ndarray:
    """Return intercept + slope x_new."""
    raise NotImplementedError


def residuals(x: np.ndarray, y: np.ndarray, fit: dict[str, float]) -> np.ndarray:
    """Return y - yhat."""
    raise NotImplementedError


def correlation_coefficient(x: np.ndarray, y: np.ndarray) -> float:
    """Return Pearson's r = S_xy / sqrt(S_xx S_yy)."""
    raise NotImplementedError


def correlation_t_test(r: float, n: int) -> tuple[float, float]:
    """Test H0: the population correlation is zero. Return (t, two-sided p) with t = r sqrt(n-2) / sqrt(1-r^2), df = n-2."""
    raise NotImplementedError


def slope_confidence_interval(fit: dict[str, float], conf: float = 0.95) -> tuple[float, float]:
    """Return slope -/+ t*_{n-2} se_slope."""
    raise NotImplementedError


def confidence_band(fit: dict[str, float], x_new: np.ndarray, conf: float = 0.95) -> np.ndarray:
    """Return the half-width of the ``conf`` confidence band for the *mean* response at x_new:
    t* s_e sqrt(1/n + (x_new - x_mean)^2 / S_xx). Narrowest at x_mean, where it equals t* s_e / sqrt(n)."""
    raise NotImplementedError


def prediction_band(fit: dict[str, float], x_new: np.ndarray, conf: float = 0.95) -> np.ndarray:
    """Return the half-width of the ``conf`` prediction band for a *new observation* at x_new:
    t* s_e sqrt(1 + 1/n + (x_new - x_mean)^2 / S_xx). Always wider than the confidence band by the extra 1."""
    raise NotImplementedError


def outlier_flags(x: np.ndarray, y: np.ndarray, fit: dict[str, float], k: float = 2.0) -> np.ndarray:
    """Return a boolean mask of points whose residual exceeds k s_e in magnitude (the textbook rule uses about 2 s)."""
    raise NotImplementedError


# ---------------------------------------------------------------- general linear least squares

def design_matrix_poly(x: np.ndarray, degree: int) -> np.ndarray:
    """Return the (n, degree+1) Vandermonde matrix [1, x, x^2, ..., x^degree]."""
    raise NotImplementedError


def linear_least_squares(A: np.ndarray, y: np.ndarray) -> dict[str, np.ndarray | float]:
    """Solve min |A beta - y|^2 by QR. Return {'beta', 'cov', 's_e', 'df', 'residuals'}.

    A = Q R, beta = R^{-1} Q^T y (no normal equations: cond(R) = cond(A), not cond(A)^2). The parameter
    covariance is s_e^2 (A^T A)^{-1} = s_e^2 R^{-1} R^{-T} with s_e^2 = SS_res / (n - p). Raise ``ValueError``
    if n <= p or A is rank deficient.
    """
    raise NotImplementedError


def general_band(A_new: np.ndarray, cov: np.ndarray, s_e: float, df: float, conf: float = 0.95,
                 prediction: bool = False) -> np.ndarray:
    """Return band half-widths for a linear model at the rows of A_new.

    Confidence: t* sqrt(a^T cov a); prediction: t* sqrt(a^T cov a + s_e^2). Reduces to ``confidence_band`` /
    ``prediction_band`` for the straight line.
    """
    raise NotImplementedError


def polynomial_fit(x: np.ndarray, y: np.ndarray, degree: int) -> dict[str, np.ndarray | float]:
    """Fit a polynomial of the given degree by ``linear_least_squares`` on ``design_matrix_poly``."""
    raise NotImplementedError


def loo_cv_error(x: np.ndarray, y: np.ndarray, degree: int) -> float:
    """Return the leave-one-out cross-validation RMS prediction error of a polynomial of the given degree.

    For each i, fit to all points but i and predict y_i; RMS of the n prediction errors. Unlike the residual
    s_e, this penalises over-fitting, so its minimum over degrees picks the model that predicts best.
    """
    raise NotImplementedError


def choose_degree(x: np.ndarray, y: np.ndarray, degrees: list[int]) -> tuple[int, np.ndarray]:
    """Return (degree with the smallest LOO error, array of LOO errors in the order of ``degrees``)."""
    raise NotImplementedError


# ---------------------------------------------------------------- nonlinear least squares

def gauss_newton(model: Callable[[np.ndarray, np.ndarray], np.ndarray],
                 jac: Callable[[np.ndarray, np.ndarray], np.ndarray],
                 x: np.ndarray, y: np.ndarray, p0: np.ndarray, tol: float = 1e-10,
                 max_iter: int = 100) -> dict[str, np.ndarray | float | int]:
    """Minimise |y - model(x, p)|^2 by damped Gauss-Newton. Return {'params', 'cov', 's_e', 'df', 'iterations', 'ss_res'}.

    Each step solves the linearised least-squares problem J dp = r by QR (``linear_least_squares`` on the
    Jacobian) and halves the step until the sum of squares decreases. Stops when the relative step is below
    ``tol``. The covariance is s_e^2 (J^T J)^{-1} at the solution (the linearisation that also gives the
    bands). Raise ``RuntimeError`` if ``max_iter`` is reached.
    """
    raise NotImplementedError


def nonlinear_band(jac: Callable[[np.ndarray, np.ndarray], np.ndarray], x_new: np.ndarray, params: np.ndarray,
                   cov: np.ndarray, s_e: float, df: float, conf: float = 0.95, prediction: bool = False) -> np.ndarray:
    """Return band half-widths for a nonlinear fit by the delta method: ``general_band`` on the Jacobian rows at x_new."""
    raise NotImplementedError


# ---------------------------------------------------------------- total impulse

def total_impulse_trapezoid(t: np.ndarray, F: np.ndarray) -> float:
    """Return the trapezoid-rule integral of F over t [N s], written out (no library call)."""
    raise NotImplementedError


def total_impulse_model(params: np.ndarray) -> float:
    """Return the closed-form integral of ``thrust_model`` over [0, inf): a Gamma(b + 1) / c^(b + 1)."""
    raise NotImplementedError


def impulse_gradient(params: np.ndarray) -> np.ndarray:
    """Return d(total impulse)/d(a, b, c) analytically: I/a, I (psi(b+1) - ln c), -I (b+1)/c with psi the digamma function."""
    raise NotImplementedError


def impulse_uncertainty_delta(params: np.ndarray, cov: np.ndarray) -> float:
    """Return the delta-method standard deviation of the total impulse: sqrt(g^T cov g) with g = ``impulse_gradient``."""
    raise NotImplementedError


def impulse_uncertainty_mc(params: np.ndarray, cov: np.ndarray, n_samples: int = 20000, seed: int = 0) -> float:
    """Return the standard deviation of the total impulse over ``n_samples`` SYNTHETIC parameter draws from N(params, cov)."""
    raise NotImplementedError


# ---------------------------------------------------------------- exponential tail by log transform

def exponential_fit_by_log(t: np.ndarray, F: np.ndarray) -> dict[str, float]:
    """Fit F = F0 exp(-t / tau) by a straight line through (t, ln F). Return {'F0', 'tau', 'r2_log'}.

    Requires F > 0. Minimises squared errors in ln F, which weights small F heavily: compare with a direct fit.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- the lesson's numbers

def main(out_dir: str | Path = ".") -> None:
    """Print the fit tables and write thrust_curve_fit.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit

    out_dir = Path(out_dir)
    np.set_printoptions(precision=4, suppress=True)

    print("== 1. thrust against chamber pressure (SYNTHETIC static-fire series) ==")
    p, F = synthetic_pressure_thrust(seed=18)
    fit = linear_regression(p, F)
    print(f"n = {fit['n']:.0f}: slope {fit['slope']:.2f} +/- {fit['se_slope']:.2f} N/MPa (true 290.00), "
          f"intercept {fit['intercept']:.1f} +/- {fit['se_intercept']:.1f} N (true 0), r = {fit['r']:.5f}, r^2 = {fit['r2']:.5f}, s_e = {fit['s_e']:.1f} N")
    lo, hi = slope_confidence_interval(fit, 0.95)
    print(f"95 % CI for the slope: ({lo:.2f}, {hi:.2f}) N/MPa; C_F A_t = slope x 1e-6 = {fit['slope'] * 1e-6:.3e} m^2")
    t_r, p_r = correlation_t_test(fit["r"], int(fit["n"]))
    print(f"significance of r: t = {t_r:.1f}, p = {p_r:.1e}")
    lr = stats.linregress(p, F)
    print(f"scipy.stats.linregress: slope {lr.slope:.2f}, intercept {lr.intercept:.1f}, r {lr.rvalue:.5f}, stderr {lr.stderr:.2f}, p {lr.pvalue:.1e}")
    for pc in (2.0, 4.0, 6.0, 8.0):
        print(f"  at p_c = {pc:.0f} MPa: F = {float(predict_line(fit, pc)):.0f} N, 95 % confidence +/- {float(confidence_band(fit, pc)):.1f} N, "
              f"95 % prediction +/- {float(prediction_band(fit, pc)):.1f} N{'   (extrapolation)' if pc > 6.0 else ''}")
    p_o, F_o = synthetic_pressure_thrust(seed=18, with_outlier=True)
    fit_o = linear_regression(p_o, F_o)
    flags = outlier_flags(p_o, F_o, fit_o, k=2.0)
    print(f"with a +120 N logging fault at point 9: slope {fit_o['slope']:.2f}, s_e {fit_o['s_e']:.1f} N, flagged points {np.flatnonzero(flags)}, "
          f"residual there {residuals(p_o, F_o, fit_o)[8]:.1f} N = {residuals(p_o, F_o, fit_o)[8] / fit_o['s_e']:.2f} s_e")
    fit_c = linear_regression(p_o[~flags], F_o[~flags])
    print(f"refit without it: slope {fit_c['slope']:.2f}, s_e {fit_c['s_e']:.1f} N")

    print("\n== 2. the thrust-time curve (SYNTHETIC, 60 samples, noise 10 N) ==")
    t, Ft = synthetic_thrust_curve(seed=18)
    degrees = [1, 2, 3, 4, 5, 6, 8, 10]
    best, errs = choose_degree(t, Ft, degrees)
    print("polynomial degree :", degrees)
    print("residual s_e [N]  :", np.array([polynomial_fit(t, Ft, d)["s_e"] for d in degrees]))
    print("LOO RMS error [N] :", errs, f"-> degree {best}")
    gn = gauss_newton(thrust_model, thrust_model_jacobian, t, Ft, np.array([1000.0, 0.5, 1.0]))
    a, b, c = gn["params"]
    sd = np.sqrt(np.diag(gn["cov"]))
    print(f"Gauss-Newton from (1000, 0.5, 1.0): {gn['iterations']} iterations, a = {a:.0f} +/- {sd[0]:.0f}, b = {b:.3f} +/- {sd[1]:.3f}, "
          f"c = {c:.3f} +/- {sd[2]:.3f} 1/s, s_e = {gn['s_e']:.2f} N (true 1550, 0.8, 1.6; noise 10)")
    popt, pcov = curve_fit(lambda tt, aa, bb, cc: thrust_model(tt, (aa, bb, cc)), t, Ft, p0=[1000.0, 0.5, 1.0])
    print(f"scipy curve_fit: {popt}, parameter sd {np.sqrt(np.diag(pcov))}")
    print(f"peak thrust at t = b/c = {b / c:.3f} s: {float(thrust_model(b / c, gn['params'])):.1f} N")
    I_trap = total_impulse_trapezoid(t, Ft)
    I_model = total_impulse_model(gn["params"])
    I_true = total_impulse_model(TRUE_THRUST_PARAMS)
    sd_delta = impulse_uncertainty_delta(gn["params"], gn["cov"])
    sd_mc = impulse_uncertainty_mc(gn["params"], gn["cov"], 20000, seed=18)
    print(f"total impulse: trapezoid of the data {I_trap:.1f} N s (misses the first 0.05 s and the tail beyond 4 s), "
          f"model {I_model:.1f} +/- {sd_delta:.1f} N s (delta method), +/- {sd_mc:.1f} (Monte Carlo); true {I_true:.1f} N s")
    window = (t >= 1.0) & (t <= 2.5)
    ex = exponential_fit_by_log(t[window], Ft[window])
    print(f"tail-off (1.0 <= t <= 2.5 s, thrust 60-310 N) by log transform: F0 = {ex['F0']:.0f} N, tau = {ex['tau']:.3f} s, "
          f"r^2 in log space {ex['r2_log']:.3f}; the model's local decay time 1/(c - b/t) is {1 / (c - b / 1.0):.3f} s at 1.0 s "
          f"and {1 / (c - b / 2.5):.3f} s at 2.5 s")
    late = t > 2.5
    ex_late = exponential_fit_by_log(t[late], np.clip(Ft[late], 1.0, None))
    print(f"same fit on t > 2.5 s (thrust below 60 N, noise 10 N, values clipped at 1 N): tau = {ex_late['tau']:.3f} s, "
          f"r^2 in log space {ex_late['r2_log']:.3f}: the log transform amplifies noise near zero")

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.2))
    pp = np.linspace(1.5, 8.0, 200)
    yhat = predict_line(fit, pp)
    axes[0].plot(p, F, "o", label="synthetic static fires")
    axes[0].plot(pp, yhat, "C1", label=f"fit: {fit['slope']:.0f} p + {fit['intercept']:.0f}")
    axes[0].fill_between(pp, yhat - confidence_band(fit, pp), yhat + confidence_band(fit, pp), color="C1", alpha=0.3, label="95 % confidence band")
    axes[0].fill_between(pp, yhat - prediction_band(fit, pp), yhat + prediction_band(fit, pp), color="C1", alpha=0.12, label="95 % prediction band")
    axes[0].axvline(6.0, color="gray", ls=":")
    axes[0].set_xlabel("chamber pressure [MPa]")
    axes[0].set_ylabel("thrust [N]")
    axes[0].legend(fontsize=7)
    axes[0].set_title("linear regression (synthetic)")
    tt = np.linspace(0.0, 4.5, 400)
    Fhat = thrust_model(tt, gn["params"])
    cb = nonlinear_band(thrust_model_jacobian, tt, gn["params"], gn["cov"], gn["s_e"], gn["df"])
    pb = nonlinear_band(thrust_model_jacobian, tt, gn["params"], gn["cov"], gn["s_e"], gn["df"], prediction=True)
    axes[1].plot(t, Ft, ".", ms=4, label="synthetic thrust samples")
    axes[1].plot(tt, Fhat, "C1", label="Gauss-Newton fit")
    axes[1].fill_between(tt, Fhat - cb, Fhat + cb, color="C1", alpha=0.35, label="95 % confidence band")
    axes[1].fill_between(tt, Fhat - pb, Fhat + pb, color="C1", alpha=0.12, label="95 % prediction band")
    axes[1].plot(tt, thrust_model(tt, TRUE_THRUST_PARAMS), "k--", lw=0.8, label="true curve")
    axes[1].set_xlabel("time [s]")
    axes[1].set_ylabel("thrust [N]")
    axes[1].legend(fontsize=7)
    axes[1].set_title("thrust curve fit with bands (synthetic)")
    axes[2].semilogy(degrees, [polynomial_fit(t, Ft, d)["s_e"] for d in degrees], "s-", label="residual $s_e$")
    axes[2].semilogy(degrees, errs, "o-", label="leave-one-out RMS error")
    axes[2].axhline(10.0, color="gray", ls=":", label="true noise 10 N")
    axes[2].set_xlabel("polynomial degree")
    axes[2].set_ylabel("error [N]")
    axes[2].legend(fontsize=8)
    axes[2].set_title("choosing the degree")
    fig.tight_layout()
    fig.savefig(out_dir / "thrust_curve_fit.png", dpi=150)
    print(f"saved {out_dir / 'thrust_curve_fit.png'}")


if __name__ == "__main__":
    main()
