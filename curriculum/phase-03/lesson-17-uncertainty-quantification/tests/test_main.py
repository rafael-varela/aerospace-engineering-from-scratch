"""Public tests for P03.17 Uncertainty Quantification. They assert properties of the samplers, of the
Monte Carlo estimators (with seeded generators and tolerances justified by the estimators' own standard
errors), of the Hermite polynomials and Gauss-Hermite rule against numpy.polynomial.hermite_e and
scipy.stats, and of the polynomial-chaos expansion (exactness on polynomials, the closed-form
coefficients of exp(a xi), moments and Sobol indices), plus the observed N^(-1/2) order of Monte Carlo.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-03/parameterized-cases/.
"""
import math

import numpy as np
import numpy.testing as npt
import pytest
from numpy.polynomial import hermite_e as H
from scipy import stats

from conftest import load_impl

impl = load_impl(__file__)

P = impl.RangeInputs()
F = lambda xi: impl.range_model(xi, P)
D = 3


# ---------------------------------------------------------------- the model

def test_range_model_is_breguet_and_has_the_right_monotonicities():
    R0 = impl.range_model(np.zeros(3), P)
    expected = (P.V / (impl.G0 * P.c_T_mean)) * P.LD_mean * np.log((P.W_zf + P.W_fuel_mean) / P.W_zf)
    assert R0 == pytest.approx(expected)
    assert 1.4e6 < R0 < 1.6e6                                     # about 1500 km: a light twin's range
    assert impl.range_model(np.array([1.0, 0.0, 0.0]), P) < R0    # worse fuel consumption, less range
    assert impl.range_model(np.array([0.0, 1.0, 0.0]), P) > R0    # better L/D, more range
    assert impl.range_model(np.array([0.0, 0.0, 1.0]), P) > R0    # more fuel, more range
    out = impl.range_model(np.zeros((5, 3)), P)
    assert out.shape == (5,) and np.all(out == R0)
    # Linear in L/D: equal steps in xi_2 give equal steps in R.
    r = [impl.range_model(np.array([0.0, k, 0.0]), P) for k in (-1.0, 0.0, 1.0)]
    assert (r[2] - r[1]) == pytest.approx(r[1] - r[0], rel=1e-12)
    # Convex in c_T (1/x): the chord lies above the curve, so the mean shifts up.
    r = [impl.range_model(np.array([k, 0.0, 0.0]), P) for k in (-1.0, 0.0, 1.0)]
    assert 0.5 * (r[0] + r[2]) > r[1]


# ---------------------------------------------------------------- sampling

def test_box_muller_produces_standard_normals():
    n = 200_000
    z = impl.standard_normals_box_muller(n, seed=1)
    assert z.shape == (n,)
    # Tolerances: the sample mean has sd 1/sqrt(n) = 0.0022, the variance 1.4/sqrt(n) = 0.0032,
    # skewness and excess kurtosis about sqrt(6/n) = 0.0055 and sqrt(24/n) = 0.011; we allow 5 sd.
    assert abs(np.mean(z)) < 0.011
    assert abs(np.var(z) - 1.0) < 0.016
    assert abs(stats.skew(z)) < 0.03
    assert abs(stats.kurtosis(z)) < 0.06
    ks = stats.kstest(z, "norm")
    assert ks.pvalue > 1e-3                                       # not rejected at any sane level
    npt.assert_array_equal(z, impl.standard_normals_box_muller(n, seed=1))       # reproducible
    assert impl.standard_normals_box_muller(7, seed=2).shape == (7,)             # odd n


def test_inversion_method_reproduces_an_exponential_and_a_normal():
    rng = np.random.default_rng(2)
    u = rng.random(100_000)
    x = impl.inverse_cdf_sample(u, lambda uu: -np.log(1.0 - uu) / 2.0)        # exponential, rate 2
    assert np.mean(x) == pytest.approx(0.5, abs=5 * 0.5 / math.sqrt(u.size))  # mean 1/rate, sd 1/rate
    assert stats.kstest(x, "expon", args=(0, 0.5)).pvalue > 1e-3
    z = impl.inverse_cdf_sample(u, stats.norm.ppf)
    assert stats.kstest(z, "norm").pvalue > 1e-3
    with pytest.raises(ValueError):
        impl.inverse_cdf_sample(np.array([0.0, 0.5]), stats.norm.ppf)


# ---------------------------------------------------------------- Monte Carlo

def test_monte_carlo_estimator_and_its_standard_error():
    mc = impl.monte_carlo(F, D, 20_000, seed=17)
    assert mc["samples"].shape == (20_000,)
    assert mc["se_mean"] == pytest.approx(mc["std"] / math.sqrt(20_000))
    c, idx = impl.pce_projection(F, D, 4, 5)
    mean_ref, var_ref = impl.pce_moments(c, idx)
    # The MC mean lies within 4 standard errors of the (quadrature-exact) PCE mean: a 6e-5 false-alarm rate.
    assert abs(mc["mean"] - mean_ref) < 4 * mc["se_mean"]
    # The sample std is within 5 % of the PCE std (its own sd is about std / sqrt(2 n) = 0.5 %).
    assert mc["std"] == pytest.approx(math.sqrt(var_ref), rel=0.05)
    # Same seed, same numbers; different seed, different numbers.
    assert impl.monte_carlo(F, D, 500, seed=17)["mean"] == impl.monte_carlo(F, D, 500, seed=17)["mean"]
    assert impl.monte_carlo(F, D, 500, seed=17)["mean"] != impl.monte_carlo(F, D, 500, seed=18)["mean"]
    # Tail probability and its binomial standard error, against the normal approximation of a linear model.
    pr, se = impl.mc_probability(mc["samples"], P.R_required)
    assert se == pytest.approx(math.sqrt(pr * (1 - pr) / 20_000))
    lin = impl.linear_propagation(F, D)
    assert abs(pr - stats.norm.cdf((P.R_required - lin["mean1"]) / lin["std"])) < 0.01    # ~4 se at p ~ 0.05
    assert impl.mc_probability(np.array([1.0, 2.0, 3.0, 4.0]), 2.5) == (0.5, pytest.approx(0.25))


def test_bootstrap_interval_covers_the_reference_and_narrows_with_n():
    mc = impl.monte_carlo(F, D, 4_000, seed=17)
    lo, hi = impl.bootstrap_ci(mc["samples"], np.mean, n_boot=1000, seed=1)
    assert lo < mc["mean"] < hi
    # For the mean the percentile bootstrap matches the normal interval mean +- 1.96 se to about 10 %.
    half = 0.5 * (hi - lo)
    assert half == pytest.approx(1.96 * mc["se_mean"], rel=0.10)
    c, idx = impl.pce_projection(F, D, 4, 5)
    assert lo < impl.pce_moments(c, idx)[0] < hi
    # A statistic without a textbook standard error: the shortfall probability.
    lo_p, hi_p = impl.bootstrap_ci(mc["samples"], lambda y: float(np.mean(y < P.R_required)), n_boot=1000, seed=2)
    pr, se = impl.mc_probability(mc["samples"], P.R_required)
    assert lo_p < pr < hi_p
    assert 0.5 * (hi_p - lo_p) == pytest.approx(1.96 * se, rel=0.15)
    # More data, narrower interval (roughly 1/sqrt(2) for twice the data).
    mc2 = impl.monte_carlo(F, D, 16_000, seed=17)
    lo2, hi2 = impl.bootstrap_ci(mc2["samples"], np.mean, n_boot=1000, seed=1)
    assert (hi2 - lo2) < 0.6 * (hi - lo)


def test_antithetic_sampling_reduces_the_variance_of_a_nearly_linear_model():
    plain = impl.monte_carlo(F, D, 10_000, seed=17)
    anti = impl.antithetic_monte_carlo(F, D, 10_000, seed=17)
    assert anti["se_mean"] < 0.2 * plain["se_mean"]               # the linear part cancels exactly
    c, idx = impl.pce_projection(F, D, 4, 5)
    mean_ref = impl.pce_moments(c, idx)[0]
    assert abs(anti["mean"] - mean_ref) < 4 * anti["se_mean"]
    # For a purely odd function antithetic pairs give the exact mean (zero) with zero spread.
    odd = lambda xi: xi[:, 0] ** 3 + 2.0 * xi[:, 1]
    res = impl.antithetic_monte_carlo(odd, 2, 1000, seed=3)
    assert abs(res["mean"]) < 1e-12 and res["se_mean"] < 1e-12
    with pytest.raises(ValueError):
        impl.antithetic_monte_carlo(F, D, 1001, seed=3)


def test_monte_carlo_error_falls_as_one_over_sqrt_n():
    """Observed order test: the RMS error of the MC mean over 100 replications at N = 100, 400, 1600 has
    log-log slope -0.5. With 100 replications the RMS estimate has a relative sd of about 7 %, which over
    the 16-fold range of N perturbs the fitted slope by about 0.05; we allow 0.12."""
    c, idx = impl.pce_projection(F, D, 4, 5)
    mean_ref = impl.pce_moments(c, idx)[0]
    ns = np.array([100, 400, 1600])
    rms = impl.mc_error_study(F, D, ns, n_rep=100, seed=3, exact_mean=mean_ref)
    assert np.all(np.diff(rms) < 0)
    slope = np.polyfit(np.log(ns), np.log(rms), 1)[0]
    assert slope == pytest.approx(-0.5, abs=0.12)
    # And the absolute level is the predicted std / sqrt(N) to 20 %.
    npt.assert_allclose(rms, math.sqrt(impl.pce_moments(c, idx)[1]) / np.sqrt(ns), rtol=0.2)


# ---------------------------------------------------------------- linear propagation

def test_linear_propagation_matches_the_hand_sensitivities():
    lin = impl.linear_propagation(F, D)
    R0 = impl.range_model(np.zeros(3), P)
    # dR/dxi: -s1 R0 (1/x), +s2 R0 (linear), s3 R0 f / ((1 + f) ln(1 + f)) with f = W_fuel / W_zf.
    f = P.W_fuel_mean / P.W_zf
    expected = R0 * np.array([-P.c_T_rel_std, P.LD_rel_std, P.W_fuel_rel_std * f / ((1 + f) * math.log(1 + f))])
    npt.assert_allclose(lin["grad"], expected, rtol=1e-6)
    assert lin["std"] == pytest.approx(np.linalg.norm(expected), rel=1e-6)
    assert lin["contributions"].sum() == pytest.approx(1.0)
    assert np.argmax(lin["contributions"]) == 1                   # L/D dominates
    assert lin["mean2"] > lin["mean1"]                            # convexity in c_T shifts the mean up
    # Exact on a linear function; exact second-order mean on a quadratic.
    lin_q = impl.linear_propagation(lambda xi: 3.0 + 2.0 * np.atleast_2d(xi)[:, 0] + np.atleast_2d(xi)[:, 1] ** 2, 2)
    assert lin_q["std"] == pytest.approx(2.0, rel=1e-6) and lin_q["mean2"] == pytest.approx(4.0, rel=1e-5)


# ---------------------------------------------------------------- Hermite polynomials and quadrature

def test_hermite_polynomials_match_numpy_and_are_orthogonal():
    x = np.linspace(-3, 3, 25)
    for n in range(7):
        npt.assert_allclose(impl.hermite_e(n, x), H.hermeval(x, [0] * n + [1]), rtol=1e-12, atol=1e-12)
    npt.assert_allclose(impl.hermite_e(3, x), x ** 3 - 3 * x)
    nodes, w = impl.gauss_hermite_e(10)
    for m in range(6):
        for n in range(6):
            inner = np.sum(w * impl.hermite_e(m, nodes) * impl.hermite_e(n, nodes))
            assert inner == pytest.approx(math.factorial(n) if m == n else 0.0, abs=1e-9)
    with pytest.raises(ValueError):
        impl.hermite_e(-1, x)


def test_gauss_hermite_rule_matches_numpy_and_integrates_moments_exactly():
    for n in (1, 2, 5, 8):
        nodes, w = impl.gauss_hermite_e(n)
        x_ref, w_ref = H.hermegauss(n)                            # weight exp(-x^2/2): divide by sqrt(2 pi)
        npt.assert_allclose(np.sort(nodes), np.sort(x_ref), atol=1e-10)
        npt.assert_allclose(w[np.argsort(nodes)], w_ref[np.argsort(x_ref)] / math.sqrt(2 * math.pi), atol=1e-12)
        assert np.sum(w) == pytest.approx(1.0)
        # E[xi^(2k)] = (2k - 1)!! exactly for 2k <= 2n - 1, and odd moments vanish.
        for k in range(0, n):
            assert np.sum(w * nodes ** (2 * k)) == pytest.approx(math.prod(range(1, 2 * k, 2)) if k else 1.0, rel=1e-10, abs=1e-10)
        for k in range(0, n):
            assert np.sum(w * nodes ** (2 * k + 1)) == pytest.approx(0.0, abs=1e-9)
    # The rule with n nodes is NOT exact for degree 2n: E[xi^4] = 3 needs at least 3 nodes.
    nodes, w = impl.gauss_hermite_e(2)
    assert np.sum(w * nodes ** 4) != pytest.approx(3.0)
    with pytest.raises(ValueError):
        impl.gauss_hermite_e(0)


# ---------------------------------------------------------------- polynomial chaos

def test_pce_basis_norms_and_multi_indices():
    idx = impl.total_degree_multi_index(3, 2)
    assert idx.shape == (10, 3) and np.all(idx.sum(axis=1) <= 2)
    npt.assert_array_equal(idx[0], [0, 0, 0])
    norms = impl.pce_norms(np.array([[0, 0], [1, 0], [2, 1], [3, 3]]))
    npt.assert_array_equal(norms, [1, 1, 2, 36])
    xi = np.array([[0.5, -1.0]])
    Psi = impl.pce_basis(np.array([[0, 0], [1, 1], [2, 0]]), xi)
    npt.assert_allclose(Psi[0], [1.0, 0.5 * -1.0, 0.5 ** 2 - 1.0])


def test_pce_is_exact_for_polynomials_and_for_the_exponential_coefficients():
    # A quadratic in two variables is reproduced exactly by a degree-2 expansion with 3 nodes per axis.
    g = lambda xi: 1.0 + 2.0 * xi[:, 0] - xi[:, 1] + 0.5 * xi[:, 0] ** 2 + 3.0 * xi[:, 0] * xi[:, 1]
    c, idx = impl.pce_projection(g, 2, 2, 3)
    rng = np.random.default_rng(4)
    xi = rng.standard_normal((50, 2))
    npt.assert_allclose(impl.pce_evaluate(c, idx, xi), g(xi), rtol=1e-12, atol=1e-12)
    m, v = impl.pce_moments(c, idx)
    assert m == pytest.approx(1.5)                                # E[g] = 1 + 0.5 E[xi^2] = 1.5
    # Var g by hand: 2^2 + 1 + 0.5^2 Var(xi^2 = He2 + 1) = 0.25 * 2 + 3^2 * 1 = 4 + 1 + 0.5 + 9 = 14.5
    assert v == pytest.approx(14.5)
    # exp(a xi) = e^{a^2/2} sum_n a^n / n! He_n(xi): closed-form coefficients.
    a = 0.7
    c, idx = impl.pce_projection(lambda xi: np.exp(a * xi[:, 0]), 1, 6, 12)
    expected = math.exp(a ** 2 / 2) * a ** np.arange(7) / np.array([math.factorial(n) for n in range(7)])
    npt.assert_allclose(c, expected, rtol=1e-8, atol=1e-10)
    m, v = impl.pce_moments(c, idx)
    assert m == pytest.approx(math.exp(a ** 2 / 2), rel=1e-9)               # lognormal mean
    assert v == pytest.approx(math.exp(a ** 2) * (math.exp(a ** 2) - 1), rel=1e-5)   # truncation at degree 6


def test_pce_moments_converge_with_degree_and_regression_agrees_with_projection():
    stds = []
    for deg in (1, 2, 3, 4):
        c, idx = impl.pce_projection(F, D, deg, deg + 1)
        stds.append(math.sqrt(impl.pce_moments(c, idx)[1]))
    diffs = np.abs(np.diff(stds))
    assert diffs[0] > diffs[1] > diffs[2]                        # fast (spectral) convergence
    assert diffs[2] < 1e-6 * stds[-1]
    c_proj, idx = impl.pce_projection(F, D, 2, 3)
    c_reg, idx2 = impl.pce_regression(F, D, 2, 600, seed=11)
    npt.assert_array_equal(idx, idx2)
    npt.assert_allclose(c_reg, c_proj, atol=2e-4 * abs(c_proj[0]))  # sampling error of a 600-point fit
    with pytest.raises(ValueError):
        impl.pce_regression(F, D, 2, 15, seed=11)
    # The PCE mean reproduces the second-order linear mean (the model is nearly quadratic) ...
    lin = impl.linear_propagation(F, D)
    m, v = impl.pce_moments(c_proj, idx)
    assert m == pytest.approx(lin["mean2"], rel=1e-5)
    # ... and the PCE std is a little above the linear std (the curvature adds variance).
    assert math.sqrt(v) > lin["std"] and math.sqrt(v) == pytest.approx(lin["std"], rel=0.01)


def test_sobol_indices_from_pce_match_the_variance_shares():
    c, idx = impl.pce_projection(F, D, 4, 5)
    S, ST = impl.sobol_indices(c, idx)
    assert np.all(S >= 0) and np.all(ST >= S - 1e-12)
    assert 0.99 < S.sum() <= 1.0 + 1e-12                          # nearly additive model
    assert S[1] > S[0] > S[2]                                     # L/D, then c_T, then fuel
    lin = impl.linear_propagation(F, D)
    npt.assert_allclose(S, lin["contributions"], atol=0.01)
    # An additive test function: S = ST exactly; a pure interaction: S = 0, ST = 1.
    c, idx = impl.pce_projection(lambda xi: xi[:, 0] + 3.0 * xi[:, 1] ** 2, 2, 3, 4)
    S, ST = impl.sobol_indices(c, idx)
    npt.assert_allclose(S, ST, atol=1e-12)
    npt.assert_allclose(S, [1.0 / 19.0, 18.0 / 19.0], atol=1e-10)            # Var = 1 + 9 * 2
    c, idx = impl.pce_projection(lambda xi: xi[:, 0] * xi[:, 1], 2, 2, 3)
    S, ST = impl.sobol_indices(c, idx)
    npt.assert_allclose(S, [0.0, 0.0], atol=1e-12)
    npt.assert_allclose(ST, [1.0, 1.0], atol=1e-12)


def test_surrogate_tail_probability_agrees_with_monte_carlo():
    c, idx = impl.pce_projection(F, D, 3, 4)
    rng = np.random.default_rng(5)
    xi = rng.standard_normal((200_000, D))
    p_sur = np.mean(impl.pce_evaluate(c, idx, xi) < P.R_required)
    mc = impl.monte_carlo(F, D, 50_000, seed=17)
    p_mc, se = impl.mc_probability(mc["samples"], P.R_required)
    assert 0.01 < p_mc < 0.2                                      # a tail, not a coin toss
    assert abs(p_sur - p_mc) < 4 * se
    # The surrogate is the model to 0.1 % everywhere in the +-3 sigma box.
    npt.assert_allclose(impl.pce_evaluate(c, idx, xi[:2000]), F(xi[:2000]), rtol=1e-3)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "range_uq.png").stat().st_size > 1000
