"""Public tests for P02.18 Regression & Curve Fitting. They assert properties, with scipy.stats.linregress,
numpy.linalg.lstsq and scipy.optimize.curve_fit as the oracles.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-02/parameterized-cases/.

Statistical tolerances: every simulation uses a fixed seed, so results are deterministic; coverage tolerances
are set at four standard errors (two-sided miss probability 6e-5 for a different valid generator stream) and
justified where they appear.
"""
import math

import numpy as np
import numpy.testing as npt
import pytest
from scipy import stats
from scipy.optimize import curve_fit

from conftest import load_impl

impl = load_impl(__file__)


def pressure_data():
    return impl.synthetic_pressure_thrust(seed=18)


# ---------------------------------------------------------------- simple linear regression

def test_linear_regression_matches_scipy_linregress():
    p, F = pressure_data()
    fit = impl.linear_regression(p, F)
    lr = stats.linregress(p, F)
    assert fit["slope"] == pytest.approx(lr.slope, rel=1e-12)
    assert fit["intercept"] == pytest.approx(lr.intercept, rel=1e-12)
    assert fit["r"] == pytest.approx(lr.rvalue, rel=1e-12)
    assert fit["se_slope"] == pytest.approx(lr.stderr, rel=1e-12)
    assert fit["se_intercept"] == pytest.approx(lr.intercept_stderr, rel=1e-12)
    t, pv = impl.correlation_t_test(fit["r"], int(fit["n"]))
    assert pv == pytest.approx(lr.pvalue, rel=1e-8)
    assert impl.correlation_coefficient(p, F) == pytest.approx(lr.rvalue, rel=1e-12)
    assert fit["slope"] == pytest.approx(290.0, abs=3 * fit["se_slope"])         # recovers the synthetic truth
    with pytest.raises(ValueError):
        impl.linear_regression(np.ones(5), np.arange(5.0))                      # no spread in x


def test_residual_properties_and_r_squared():
    p, F = pressure_data()
    fit = impl.linear_regression(p, F)
    res = impl.residuals(p, F, fit)
    assert np.sum(res) == pytest.approx(0.0, abs=1e-9)                           # normal equation 1
    assert np.dot(res, p) == pytest.approx(0.0, abs=1e-9)                        # normal equation 2
    ss_res, ss_tot = float(res @ res), float(np.sum((F - F.mean()) ** 2))
    assert fit["r2"] == pytest.approx(1 - ss_res / ss_tot) and fit["r2"] == pytest.approx(fit["r"] ** 2)
    assert fit["s_e"] == pytest.approx(math.sqrt(ss_res / (len(p) - 2)))
    # An exact line: r = 1, zero residual scatter, zero standard errors.
    x = np.arange(1.0, 8.0)
    exact = impl.linear_regression(x, 3.0 * x - 2.0)
    assert exact["slope"] == pytest.approx(3.0) and exact["intercept"] == pytest.approx(-2.0)
    assert exact["r"] == pytest.approx(1.0) and exact["s_e"] == pytest.approx(0.0, abs=1e-12)
    assert impl.linear_regression(x, -3.0 * x)["r"] == pytest.approx(-1.0)


def test_slope_is_invariant_to_shifts_and_scales_correctly():
    p, F = pressure_data()
    f0 = impl.linear_regression(p, F)
    f1 = impl.linear_regression(p + 10.0, F)
    assert f1["slope"] == pytest.approx(f0["slope"]) and f1["intercept"] == pytest.approx(f0["intercept"] - 10.0 * f0["slope"])
    f2 = impl.linear_regression(p, 2.0 * F)
    assert f2["slope"] == pytest.approx(2.0 * f0["slope"]) and f2["s_e"] == pytest.approx(2.0 * f0["s_e"])
    assert f2["r"] == pytest.approx(f0["r"])                                      # r is dimensionless
    f3 = impl.linear_regression(1e6 * p, F)                                        # pressure in Pa
    assert f3["slope"] == pytest.approx(f0["slope"] * 1e-6)


def test_slope_interval_coverage_over_synthetic_datasets():
    # 500 datasets: SE of the coverage = sqrt(0.95 * 0.05 / 500) = 0.0097; tolerance 4 SE = 0.04.
    hits = 0
    for seed in range(500):
        p, F = impl.synthetic_pressure_thrust(seed=1000 + seed)
        lo, hi = impl.slope_confidence_interval(impl.linear_regression(p, F), 0.95)
        hits += lo <= 290.0 <= hi
    assert hits / 500 == pytest.approx(0.95, abs=0.04)


def test_bands_have_the_right_shape_and_coverage():
    p, F = pressure_data()
    fit = impl.linear_regression(p, F)
    xs = np.linspace(1.0, 8.0, 50)
    cb, pb = impl.confidence_band(fit, xs), impl.prediction_band(fit, xs)
    assert np.all(pb > cb)                                                       # a new observation is less certain
    assert np.argmin(cb) == np.argmin(np.abs(xs - fit["x_mean"]))                # narrowest at x_mean
    tcrit = stats.t.ppf(0.975, fit["df"])
    assert float(impl.confidence_band(fit, fit["x_mean"])) == pytest.approx(tcrit * fit["s_e"] / math.sqrt(fit["n"]))
    assert float(impl.prediction_band(fit, fit["x_mean"])) == pytest.approx(tcrit * fit["s_e"] * math.sqrt(1 + 1 / fit["n"]))
    # Symmetric about x_mean, and growing without bound.
    assert float(impl.confidence_band(fit, fit["x_mean"] + 2.0)) == pytest.approx(float(impl.confidence_band(fit, fit["x_mean"] - 2.0)))
    assert float(impl.confidence_band(fit, 100.0)) > 10 * float(impl.confidence_band(fit, fit["x_mean"]))
    # Coverage at x0 = 5 MPa over 400 synthetic datasets: SE 0.011, tolerance 4 SE = 0.044.
    x0, true_mean = 5.0, 290.0 * 5.0
    rng = np.random.default_rng(77)
    hit_c = hit_p = 0
    for seed in range(400):
        pp, FF = impl.synthetic_pressure_thrust(seed=5000 + seed)
        f = impl.linear_regression(pp, FF)
        yhat = float(impl.predict_line(f, x0))
        hit_c += abs(yhat - true_mean) <= float(impl.confidence_band(f, x0))
        hit_p += abs(yhat - (true_mean + rng.normal(0.0, 25.0))) <= float(impl.prediction_band(f, x0))
    assert hit_c / 400 == pytest.approx(0.95, abs=0.044)
    assert hit_p / 400 == pytest.approx(0.95, abs=0.044)


def test_outlier_flags_catch_the_injected_fault_only():
    p, F = impl.synthetic_pressure_thrust(seed=18, with_outlier=True)
    fit = impl.linear_regression(p, F)
    flags = impl.outlier_flags(p, F, fit, k=2.0)
    assert flags.dtype == bool and flags.shape == p.shape
    assert np.flatnonzero(flags).tolist() == [8]
    assert not np.any(impl.outlier_flags(p, F, fit, k=3.0))
    clean_fit = impl.linear_regression(*pressure_data())
    assert not np.any(impl.outlier_flags(*pressure_data(), clean_fit, k=2.5))
    assert clean_fit["s_e"] < fit["s_e"]                                          # the fault inflates the scatter


# ---------------------------------------------------------------- general linear least squares

def test_polynomial_fit_matches_lstsq_and_recovers_exact_polynomials():
    t, Ft = impl.synthetic_thrust_curve(seed=18)
    A = impl.design_matrix_poly(t, 4)
    assert A.shape == (60, 5)
    npt.assert_allclose(A[:, 0], 1.0) and npt.assert_allclose(A[:, 3], t ** 3)
    fit = impl.linear_least_squares(A, Ft)
    beta_np = np.linalg.lstsq(A, Ft, rcond=None)[0]
    npt.assert_allclose(fit["beta"], beta_np, rtol=1e-8, atol=1e-8)
    resid = Ft - A @ fit["beta"]
    npt.assert_allclose(fit["residuals"], resid, atol=1e-9)
    assert fit["s_e"] == pytest.approx(math.sqrt(resid @ resid / (60 - 5)))
    npt.assert_allclose(fit["cov"], fit["s_e"] ** 2 * np.linalg.inv(A.T @ A), rtol=1e-8)
    npt.assert_allclose(A.T @ resid, 0.0, atol=1e-6)                             # residual orthogonal to every column
    # Degree 1 is the straight line.
    p, F = pressure_data()
    lin = impl.linear_regression(p, F)
    f1 = impl.polynomial_fit(p, F, 1)
    npt.assert_allclose(f1["beta"], [lin["intercept"], lin["slope"]], rtol=1e-10)
    assert math.sqrt(f1["cov"][1, 1]) == pytest.approx(lin["se_slope"], rel=1e-10)
    # An exact cubic is recovered.
    x = np.linspace(-1, 2, 15)
    c_true = np.array([1.0, -2.0, 0.5, 3.0])
    npt.assert_allclose(impl.polynomial_fit(x, impl.design_matrix_poly(x, 3) @ c_true, 3)["beta"], c_true, atol=1e-10)
    with pytest.raises(ValueError):
        impl.linear_least_squares(np.ones((3, 3)), np.ones(3))                   # n <= p
    with pytest.raises(ValueError):
        impl.linear_least_squares(np.column_stack([np.ones(6), np.ones(6)]), np.arange(6.0))   # rank deficient


def test_general_band_reduces_to_the_straight_line_bands():
    p, F = pressure_data()
    lin = impl.linear_regression(p, F)
    f1 = impl.polynomial_fit(p, F, 1)
    xs = np.linspace(1.0, 8.0, 20)
    A_new = impl.design_matrix_poly(xs, 1)
    npt.assert_allclose(impl.general_band(A_new, f1["cov"], f1["s_e"], f1["df"]), impl.confidence_band(lin, xs), rtol=1e-10)
    npt.assert_allclose(impl.general_band(A_new, f1["cov"], f1["s_e"], f1["df"], prediction=True), impl.prediction_band(lin, xs), rtol=1e-10)


def test_loo_cross_validation_picks_the_true_degree_and_exceeds_the_residual_rms():
    rng = np.random.default_rng(1818)
    x = np.linspace(0, 3, 40)
    y = 2.0 - 1.5 * x + 0.8 * x ** 2 - 0.3 * x ** 3 + rng.normal(0.0, 0.05, size=40)
    best, errs = impl.choose_degree(x, y, [1, 2, 3, 4, 5, 6])
    assert best == 3
    assert errs.shape == (6,) and np.all(errs > 0)
    # The LOO residual e_i / (1 - h_ii) is never smaller in magnitude than the ordinary residual.
    for d in (1, 3, 5):
        fit = impl.polynomial_fit(x, y, d)
        assert impl.loo_cv_error(x, y, d) >= math.sqrt(np.mean(fit["residuals"] ** 2)) - 1e-12
    with pytest.raises(ValueError):
        impl.loo_cv_error(x[:5], y[:5], 4)


# ---------------------------------------------------------------- nonlinear least squares

def test_thrust_model_shape_and_jacobian_against_finite_differences():
    params = np.array([1550.0, 0.8, 1.6])
    t = np.linspace(0.05, 4.0, 30)
    F = impl.thrust_model(t, params)
    assert np.all(F >= 0.0) and float(impl.thrust_model(0.0, params)) == 0.0
    assert t[np.argmax(F)] == pytest.approx(0.8 / 1.6, abs=t[1] - t[0])          # peak near b/c
    J = impl.thrust_model_jacobian(t, params)
    assert J.shape == (30, 3)
    for k in range(3):
        h = 1e-6 * max(1.0, abs(params[k]))
        dp = np.zeros(3)
        dp[k] = h
        fd = (impl.thrust_model(t, params + dp) - impl.thrust_model(t, params - dp)) / (2 * h)
        npt.assert_allclose(J[:, k], fd, rtol=1e-5, atol=1e-6)


def test_gauss_newton_matches_curve_fit_and_handles_the_linear_case():
    t, Ft = impl.synthetic_thrust_curve(seed=18)
    gn = impl.gauss_newton(impl.thrust_model, impl.thrust_model_jacobian, t, Ft, np.array([1000.0, 0.5, 1.0]))
    popt, pcov = curve_fit(lambda tt, a, b, c: impl.thrust_model(tt, (a, b, c)), t, Ft, p0=[1000.0, 0.5, 1.0])
    npt.assert_allclose(gn["params"], popt, rtol=1e-6)
    npt.assert_allclose(gn["cov"], pcov, rtol=1e-3)
    assert gn["ss_res"] == pytest.approx(float(np.sum((Ft - impl.thrust_model(t, popt)) ** 2)), rel=1e-8)
    assert gn["df"] == 57.0 and gn["iterations"] < 50
    # The fit recovers the synthetic truth within 3 standard errors on every parameter.
    sd = np.sqrt(np.diag(gn["cov"]))
    assert np.all(np.abs(gn["params"] - impl.TRUE_THRUST_PARAMS) < 3 * sd)
    # A linear model converges in one step to the linear least-squares answer.
    p, F = pressure_data()
    lin = impl.linear_least_squares(impl.design_matrix_poly(p, 1), F)
    gl = impl.gauss_newton(lambda x, q: q[0] + q[1] * x, lambda x, q: impl.design_matrix_poly(x, 1), p, F, np.array([0.0, 0.0]))
    npt.assert_allclose(gl["params"], lin["beta"], rtol=1e-10)
    npt.assert_allclose(gl["cov"], lin["cov"], rtol=1e-8)
    assert gl["iterations"] <= 2
    with pytest.raises(RuntimeError):
        impl.gauss_newton(impl.thrust_model, impl.thrust_model_jacobian, t, Ft, np.array([1000.0, 0.5, 1.0]), max_iter=0)


def test_nonlinear_bands_are_ordered_and_consistent_with_general_band():
    t, Ft = impl.synthetic_thrust_curve(seed=18)
    gn = impl.gauss_newton(impl.thrust_model, impl.thrust_model_jacobian, t, Ft, np.array([1000.0, 0.5, 1.0]))
    tt = np.linspace(0.1, 4.0, 25)
    cb = impl.nonlinear_band(impl.thrust_model_jacobian, tt, gn["params"], gn["cov"], gn["s_e"], gn["df"])
    pb = impl.nonlinear_band(impl.thrust_model_jacobian, tt, gn["params"], gn["cov"], gn["s_e"], gn["df"], prediction=True)
    assert np.all(pb > cb) and np.all(cb > 0)
    J = impl.thrust_model_jacobian(tt, gn["params"])
    npt.assert_allclose(cb, impl.general_band(J, gn["cov"], gn["s_e"], gn["df"]), rtol=1e-12)
    # The true curve lies inside the 95 % confidence band at most points (pointwise 95 %).
    inside = np.abs(impl.thrust_model(tt, impl.TRUE_THRUST_PARAMS) - impl.thrust_model(tt, gn["params"])) <= cb
    assert inside.mean() >= 0.8


# ---------------------------------------------------------------- total impulse

def test_total_impulse_closed_form_quadrature_gradient_and_uncertainty():
    params = np.array([1550.0, 0.8, 1.6])
    I = impl.total_impulse_model(params)
    assert I == pytest.approx(1550.0 * math.gamma(1.8) / 1.6 ** 1.8)
    assert impl.total_impulse_model(np.array([100.0, 0.0, 2.0])) == pytest.approx(50.0)     # b = 0: a / c
    t = np.linspace(0.0, 15.0, 20001)
    assert impl.total_impulse_trapezoid(t, impl.thrust_model(t, params)) == pytest.approx(I, rel=1e-4)
    assert impl.total_impulse_trapezoid(np.array([0.0, 1.0, 2.0]), np.array([0.0, 2.0, 0.0])) == pytest.approx(2.0)
    g = impl.impulse_gradient(params)
    for k in range(3):
        h = 1e-6 * max(1.0, abs(params[k]))
        dp = np.zeros(3)
        dp[k] = h
        fd = (impl.total_impulse_model(params + dp) - impl.total_impulse_model(params - dp)) / (2 * h)
        assert g[k] == pytest.approx(fd, rel=1e-5)
    tt, Ft = impl.synthetic_thrust_curve(seed=18)
    gn = impl.gauss_newton(impl.thrust_model, impl.thrust_model_jacobian, tt, Ft, np.array([1000.0, 0.5, 1.0]))
    sd_delta = impl.impulse_uncertainty_delta(gn["params"], gn["cov"])
    sd_mc = impl.impulse_uncertainty_mc(gn["params"], gn["cov"], 20000, seed=18)
    assert sd_delta == pytest.approx(sd_mc, rel=0.15)                             # linearisation holds at this noise
    assert abs(impl.total_impulse_model(gn["params"]) - I) < 3 * sd_delta
    with pytest.raises(ValueError):
        impl.total_impulse_model(np.array([1.0, 0.5, -1.0]))


def test_exponential_fit_by_log_recovers_a_noiseless_decay_and_rejects_bad_input():
    t = np.linspace(0.0, 3.0, 25)
    ex = impl.exponential_fit_by_log(t, 250.0 * np.exp(-t / 0.7))
    assert ex["F0"] == pytest.approx(250.0, rel=1e-10) and ex["tau"] == pytest.approx(0.7, rel=1e-10)
    assert ex["r2_log"] == pytest.approx(1.0)
    with pytest.raises(ValueError):
        impl.exponential_fit_by_log(t, np.exp(-t) - 0.5)                          # non-positive values
    with pytest.raises(ValueError):
        impl.exponential_fit_by_log(t, np.exp(+t))                                # growing, tau < 0


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "thrust_curve_fit.png").stat().st_size > 1000
