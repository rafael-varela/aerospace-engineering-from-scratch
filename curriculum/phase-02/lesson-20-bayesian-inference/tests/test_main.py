"""Public tests for P02.20 Bayesian Inference. They assert properties, with scipy.stats as the oracle.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-02/parameterized-cases/.

Every random quantity is seeded. Statistical tolerances are derived in a comment next to each one
from the number of (effective) samples, at four or more standard deviations.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import quad
from scipy.stats import beta as beta_dist, binom, chi2, gamma as gamma_dist, poisson

from conftest import load_impl

impl = load_impl(__file__)

M = 1e6
A0, B0 = 2.5, 1.2e6                 # heritage prior: Jeffreys + 2 failures in 1.2e6 h
KF, TF = 1, 6.0e4                   # new fleet: 1 failure in 60 000 h
LAM_REQ = -np.log(0.95) / (3.0 * 8766.0)


def heritage_grid():
    grid = np.linspace(0.0, 15.0 / M, 6001)
    pdf = impl.grid_posterior(grid, impl.log_gamma_pdf(grid, A0, B0), impl.log_poisson_likelihood(KF, TF, grid))
    return grid, pdf


# ---------------------------------------------------------------- discrete Bayes

def test_discrete_update_matches_hand_calculation_and_sums_to_one():
    prior = np.array([0.9, 0.1])
    like = np.exp(impl.log_binomial_likelihood(1, 20, np.array([0.01, 0.10])))
    post, ev = impl.bayes_discrete(prior, like)
    assert post.sum() == pytest.approx(1.0)
    # By hand: P(1 of 20 | 0.01) = 20 * 0.01 * 0.99^19, P(1 of 20 | 0.1) = 20 * 0.1 * 0.9^19.
    l_good, l_bad = 20 * 0.01 * 0.99**19, 20 * 0.1 * 0.9**19
    assert ev == pytest.approx(0.9 * l_good + 0.1 * l_bad)
    assert post[0] == pytest.approx(0.9 * l_good / (0.9 * l_good + 0.1 * l_bad))
    # A likelihood equal under all hypotheses leaves the prior unchanged; a zero likelihood kills a hypothesis.
    npt.assert_allclose(impl.bayes_discrete(prior, np.array([0.3, 0.3]))[0], prior)
    npt.assert_allclose(impl.bayes_discrete(prior, np.array([0.0, 0.3]))[0], [0.0, 1.0])
    # Updating twice with independent data equals updating once with the product of likelihoods.
    p1, _ = impl.bayes_discrete(prior, like)
    p12, _ = impl.bayes_discrete(p1, like)
    npt.assert_allclose(p12, impl.bayes_discrete(prior, like**2)[0])
    for bad_prior, bad_like in ((np.array([0.5, 0.4]), like), (np.array([-0.1, 1.1]), like), (prior, np.zeros(2))):
        with pytest.raises(ValueError):
            impl.bayes_discrete(bad_prior, bad_like)


def test_odds_form_agrees_with_the_probability_form():
    prior = np.array([0.9, 0.1])
    like = np.array([0.1652, 0.2702])
    post, _ = impl.bayes_discrete(prior, like)
    post_odds = impl.odds_update(prior[0] / prior[1], like[0] / like[1])
    assert post_odds == pytest.approx(post[0] / post[1])
    assert impl.odds_update(1.0, 1.0) == 1.0
    with pytest.raises(ValueError):
        impl.odds_update(-1.0, 2.0)


# ---------------------------------------------------------------- likelihoods and densities

def test_log_likelihoods_match_scipy_logpmf_including_edges():
    p = np.array([0.0, 1e-4, 0.01, 0.3, 0.999, 1.0])
    for k, n in ((0, 40), (1, 20), (7, 10), (10, 10)):
        ours = impl.log_binomial_likelihood(k, n, p)
        ref = binom.logpmf(k, n, p)
        npt.assert_allclose(ours, ref, rtol=1e-12, atol=1e-12)
    lam = np.array([0.0, 1e-9, 1e-6, 1e-4])
    for k, T in ((0, 6e4), (1, 6e4), (3, 1.2e6)):
        npt.assert_allclose(impl.log_poisson_likelihood(k, T, lam), poisson.logpmf(k, lam * T), rtol=1e-12, atol=1e-12)
    with pytest.raises(ValueError):
        impl.log_binomial_likelihood(5, 4, p)
    with pytest.raises(ValueError):
        impl.log_poisson_likelihood(1, 0.0, lam)


def test_gamma_and_beta_log_densities_match_scipy_and_integrate_to_one():
    lam = np.linspace(1e-8, 20.0 / M, 7)
    for a, b in ((A0, B0), (3.5, 1.26e6), (1.0, 5e4), (0.5, 1e5)):
        npt.assert_allclose(impl.log_gamma_pdf(lam, a, b), gamma_dist.logpdf(lam, a, scale=1.0 / b), rtol=1e-10)
        total, _ = quad(lambda x: np.exp(impl.log_gamma_pdf(x, a, b)), 0.0, 60.0 / b)    # 60/beta holds all but e^-60
        assert total == pytest.approx(1.0, rel=1e-6)
    assert impl.log_gamma_pdf(-1e-9, A0, B0) == -np.inf
    # Improper priors: Jeffreys is lambda^(-1/2), flat is constant (zero log density).
    x = np.array([1e-7, 1e-6, 1e-5])
    npt.assert_allclose(impl.log_gamma_pdf(x, 0.5, 0.0), -0.5 * np.log(x))
    npt.assert_allclose(impl.log_gamma_pdf(x, 1.0, 0.0), 0.0, atol=1e-15)
    p = np.array([0.0, 0.01, 0.05, 0.5, 0.99, 1.0])
    for a, b in ((1.0, 41.0), (0.5, 40.5), (3.0, 7.0)):
        npt.assert_allclose(impl.log_beta_pdf(p[1:-1], a, b), beta_dist.logpdf(p[1:-1], a, b), rtol=1e-10)
        total, _ = quad(lambda x: np.exp(impl.log_beta_pdf(x, a, b)), 0.0, 1.0)
        assert total == pytest.approx(1.0, rel=1e-6)
    assert impl.log_beta_pdf(1.5, 2.0, 2.0) == -np.inf
    with pytest.raises(ValueError):
        impl.log_gamma_pdf(lam, 0.0, 1.0)
    with pytest.raises(ValueError):
        impl.log_beta_pdf(p, 1.0, 0.0)


# ---------------------------------------------------------------- conjugate updates against the grid

def test_grid_posterior_reproduces_the_conjugate_gamma_pointwise():
    grid, pdf = heritage_grid()
    a, b = impl.gamma_poisson_update(A0, B0, KF, TF)
    assert (a, b) == (A0 + KF, B0 + TF)
    ref = gamma_dist.pdf(grid, a, scale=1.0 / b)
    # The grid normalisation uses the trapezoid rule on 6001 points: agreement to 1e-5 relative where the density matters.
    mask = ref > 1e-3 * ref.max()
    npt.assert_allclose(pdf[mask], ref[mask], rtol=2e-5)
    assert np.trapezoid(pdf, grid) == pytest.approx(1.0, rel=1e-10)
    assert impl.grid_mean(grid, pdf) == pytest.approx(a / b, rel=1e-4)
    # Improper priors give proper posteriors once there is data.
    for a0 in (0.5, 1.0):
        pdf0 = impl.grid_posterior(np.linspace(0.0, 200.0 / M, 20001), impl.log_gamma_pdf(np.linspace(0.0, 200.0 / M, 20001), a0, 0.0),
                                   impl.log_poisson_likelihood(KF, TF, np.linspace(0.0, 200.0 / M, 20001)))
        assert impl.grid_mean(np.linspace(0.0, 200.0 / M, 20001), pdf0) == pytest.approx((a0 + KF) / TF, rel=1e-3)
    with pytest.raises(ValueError):
        impl.grid_posterior(grid, np.full_like(grid, -np.inf), np.zeros_like(grid))
    with pytest.raises(ValueError):
        impl.grid_posterior(grid[::-1], np.zeros_like(grid), np.zeros_like(grid))


def test_posterior_mean_is_a_weighted_average_and_data_eventually_dominate():
    a, b = impl.gamma_poisson_update(A0, B0, KF, TF)
    w_prior = B0 / (B0 + TF)
    assert a / b == pytest.approx(w_prior * (A0 / B0) + (1.0 - w_prior) * (KF / TF))
    # With a long record at a true rate of 10 per 1e6 h the posterior mean converges to k/T and its SD shrinks as 1/sqrt(T).
    lam_true = 10.0 / M
    sds = []
    for T in (1e7, 1e8, 1e9):
        k = int(round(lam_true * T))
        aa, bb = impl.gamma_poisson_update(A0, B0, k, T)
        assert abs(aa / bb - lam_true) < 5.0 * A0 / T                       # prior pull is O(alpha0 / T)
        sds.append(np.sqrt(aa) / bb)
    # Once T >> beta0 = 1.2e6 h the posterior SD falls as 1/sqrt(T): the ratio per decade tends to sqrt(0.1).
    assert sds[1] / sds[0] == pytest.approx(np.sqrt(0.1), rel=0.10) and sds[2] / sds[1] == pytest.approx(np.sqrt(0.1), rel=0.02)
    # Sequential updating is order-independent and equals one batch update.
    events = [(0, 2.0e4), (1, 1.5e4), (0, 2.5e4)]
    assert impl.sequential_update(A0, B0, events) == impl.sequential_update(A0, B0, events[::-1])
    assert impl.sequential_update(A0, B0, events) == impl.gamma_poisson_update(A0, B0, 1, 6.0e4)
    assert impl.beta_binomial_update(1.0, 1.0, 0, 40) == (1.0, 41.0)
    assert impl.beta_binomial_update(2.0, 3.0, 2, 5) == (4.0, 6.0)
    with pytest.raises(ValueError):
        impl.gamma_poisson_update(A0, B0, -1, TF)
    with pytest.raises(ValueError):
        impl.beta_binomial_update(1.0, 1.0, 5, 4)


# ---------------------------------------------------------------- credible intervals on the grid

def test_grid_quantiles_and_intervals_match_scipy_and_hpd_is_shortest():
    grid, pdf = heritage_grid()
    a, b = A0 + KF, B0 + TF
    cdf = impl.grid_cdf(grid, pdf)
    assert cdf[0] == 0.0 and cdf[-1] == pytest.approx(1.0) and np.all(np.diff(cdf) >= 0.0)
    npt.assert_allclose(cdf[::500], gamma_dist.cdf(grid[::500], a, scale=1.0 / b), atol=2e-5)
    for q in (0.05, 0.5, 0.95):
        assert impl.grid_quantile(grid, pdf, q) == pytest.approx(gamma_dist.ppf(q, a, scale=1.0 / b), rel=2e-3)
    lo, hi = impl.credible_interval(grid, pdf, 0.9)
    npt.assert_allclose([lo, hi], gamma_dist.ppf([0.05, 0.95], a, scale=1.0 / b), rtol=2e-3)
    h_lo, h_hi = impl.hpd_interval(grid, pdf, 0.9)
    # The HPD interval holds 90 % (to the grid spacing), is shorter than the equal-tailed one, and contains the mode.
    mass = np.interp(h_hi, grid, cdf) - np.interp(h_lo, grid, cdf)
    assert mass == pytest.approx(0.9, abs=2e-3)
    assert h_hi - h_lo < hi - lo
    mode = (a - 1.0) / b
    assert h_lo <= mode <= h_hi
    assert h_lo < lo and h_hi < hi                                             # a right-skewed density: HPD sits to the left
    # Equal-tailed widths grow with the level.
    assert np.diff(impl.credible_interval(grid, pdf, 0.5))[0] < np.diff(impl.credible_interval(grid, pdf, 0.95))[0]
    with pytest.raises(ValueError):
        impl.credible_interval(grid, pdf, 1.0)
    with pytest.raises(ValueError):
        impl.grid_quantile(grid, pdf, 1.5)


# ---------------------------------------------------------------- Metropolis

def test_metropolis_samples_the_conjugate_posterior():
    a, b = A0 + KF, B0 + TF

    def log_target(x):
        return float(impl.log_gamma_pdf(x, A0, B0) + impl.log_poisson_likelihood(KF, TF, x)) if x > 0.0 else -np.inf

    samples, acc = impl.metropolis(log_target, a / b, 40_000, 1.5 * np.sqrt(a) / b, np.random.default_rng(20))
    assert samples.shape == (40_000,) and np.all(samples > 0.0)
    assert 0.2 < acc < 0.8
    tau = impl.autocorrelation_time(samples)
    assert tau >= 1.0
    n_eff = samples.size / tau
    assert n_eff > 2000
    # Mean: SD of the sample mean is sigma/sqrt(n_eff) = 1.5e-6/sqrt(6000) = 2e-8, i.e. 0.7 % of the mean; 5 % is 7 sigma.
    assert samples.mean() == pytest.approx(a / b, rel=0.05)
    # SD: relative scatter 1/sqrt(2 n_eff) = 0.9 %; 10 % is 11 sigma.
    assert samples.std() == pytest.approx(np.sqrt(a) / b, rel=0.10)
    # Tail probability: P(lambda <= requirement) = 0.33; binomial SD with n_eff = 6000 is 0.6 %; 0.03 is 5 sigma.
    assert np.mean(samples <= LAM_REQ) == pytest.approx(impl.prob_rate_below(a, b, LAM_REQ), abs=0.03)
    # Seeded: the same call gives the same chain; another seed gives another chain.
    s2, _ = impl.metropolis(log_target, a / b, 40_000, 1.5 * np.sqrt(a) / b, np.random.default_rng(20))
    npt.assert_array_equal(s2, samples)
    s3, _ = impl.metropolis(log_target, a / b, 2000, 1.5 * np.sqrt(a) / b, np.random.default_rng(21))
    assert not np.array_equal(s3, samples[:2000])
    # A tiny step accepts almost everything but mixes slowly; a huge step rejects almost everything.
    _, acc_small = impl.metropolis(log_target, a / b, 2000, 0.01 * np.sqrt(a) / b, np.random.default_rng(1))
    _, acc_big = impl.metropolis(log_target, a / b, 2000, 100.0 * np.sqrt(a) / b, np.random.default_rng(1))
    assert acc_small > 0.9 > 0.2 > acc_big
    with pytest.raises(ValueError):
        impl.metropolis(log_target, -1.0, 100, 1e-6, np.random.default_rng(0))
    with pytest.raises(ValueError):
        impl.metropolis(log_target, a / b, 0, 1e-6, np.random.default_rng(0))


# ---------------------------------------------------------------- mission quantities

def test_probability_below_requirement_matches_the_grid_cdf_and_is_monotone():
    grid, pdf = heritage_grid()
    a, b = A0 + KF, B0 + TF
    cdf = impl.grid_cdf(grid, pdf)
    for lam in (1.0 / M, LAM_REQ, 5.0 / M):
        assert impl.prob_rate_below(a, b, lam) == pytest.approx(np.interp(lam, grid, cdf), abs=1e-4)
    assert impl.prob_rate_below(a, b, 0.0) == 0.0
    assert impl.prob_rate_below(a, b, 1.0) == pytest.approx(1.0)
    lams = np.linspace(0.0, 20.0 / M, 50)
    assert np.all(np.diff([impl.prob_rate_below(a, b, x) for x in lams]) >= 0.0)
    # The lesson's decision: with the heritage prior the requirement is far from demonstrated at 90 %.
    assert impl.prob_rate_below(a, b, LAM_REQ) < 0.5
    with pytest.raises(ValueError):
        impl.prob_rate_below(a, 0.0, LAM_REQ)


def test_posterior_predictive_reliability_matches_monte_carlo_and_jensen():
    a, b = A0 + KF, B0 + TF
    t = 3.0 * 8766.0
    R = impl.posterior_predictive_reliability(a, b, t)
    # Monte Carlo over 200 000 seeded gamma draws: SD of the mean of e^(-lambda t) is below 0.05/sqrt(2e5) = 1e-4; 1e-3 is > 8 sigma.
    lam = gamma_dist.rvs(a, scale=1.0 / b, size=200_000, random_state=np.random.default_rng(3))
    assert R == pytest.approx(np.mean(np.exp(-lam * t)), abs=1e-3)
    assert R > np.exp(-a / b * t)                                               # Jensen: convex e^(-lambda t)
    assert impl.posterior_predictive_reliability(a, b, 0.0) == 1.0
    ts = np.linspace(0.0, 1e6, 20)
    assert np.all(np.diff([impl.posterior_predictive_reliability(a, b, x) for x in ts]) < 0.0)
    # Known rate (alpha, beta -> infinity at fixed ratio): the predictive tends to e^(-lambda t).
    assert impl.posterior_predictive_reliability(1e6 * a, 1e6 * b, t) == pytest.approx(np.exp(-a / b * t), rel=1e-5)


def test_hours_to_demonstrate_closes_the_loop_and_is_zero_when_already_met():
    a, b = A0 + KF, B0 + TF
    T_add = impl.hours_to_demonstrate(a, b, LAM_REQ, 0.9)
    assert T_add > 0.0
    assert impl.prob_rate_below(a, b + T_add, LAM_REQ) == pytest.approx(0.9, abs=1e-9)
    assert impl.hours_to_demonstrate(a, b, LAM_REQ, 0.99) > T_add                # a higher confidence costs more hours
    assert impl.hours_to_demonstrate(a, b, 100.0 / M, 0.9) == 0.0               # a lax requirement is already met
    with pytest.raises(ValueError):
        impl.hours_to_demonstrate(a, b, LAM_REQ, 1.0)


def test_frequentist_bound_equals_flat_prior_credible_bound():
    for k, T, conf in ((1, 6e4, 0.95), (0, 1e5, 0.9), (5, 3e6, 0.99)):
        ub = impl.frequentist_upper_bound(k, T, conf)
        assert ub == pytest.approx(chi2.ppf(conf, 2 * k + 2) / (2.0 * T))
        assert ub == pytest.approx(gamma_dist.ppf(conf, k + 1, scale=1.0 / T), rel=1e-10)
        assert impl.prob_rate_below(k + 1, T, ub) == pytest.approx(conf)
    assert impl.frequentist_upper_bound(0, 1e5, 0.95) == pytest.approx(-np.log(0.05) / 1e5)   # zero failures: 3/T rule
    with pytest.raises(ValueError):
        impl.frequentist_upper_bound(-1, 1e5, 0.95)


def test_zero_failure_demonstration_and_beta_cdf():
    n = impl.zero_failure_tests_required(0.05, 0.95)
    assert n == 58
    assert 1.0 - 0.95 ** (n + 1) >= 0.95 > 1.0 - 0.95 ** n
    assert impl.zero_failure_tests_required(0.01, 0.95) == 298
    for a, b, pm in ((1.0, 41.0, 0.05), (0.5, 40.5, 0.05), (3.0, 7.0, 0.4)):
        assert impl.beta_cdf_grid(a, b, pm) == pytest.approx(beta_dist.cdf(pm, a, b), abs=2e-4)
    assert impl.beta_cdf_grid(1.0, 41.0, 0.05) == pytest.approx(1.0 - 0.95**41, abs=2e-4)
    with pytest.raises(ValueError):
        impl.zero_failure_tests_required(0.0, 0.95)


# ---------------------------------------------------------------- synthetic fleets and calibration

def test_simulated_failures_are_seeded_poisson_counts():
    rng = np.random.default_rng(9)
    counts = np.array([impl.simulate_fleet_failures(2.0 / M, 1.0e6, rng) for _ in range(4000)])
    assert counts.dtype.kind == "i" and np.all(counts >= 0)
    # Poisson mean 2: the sample mean of 4000 draws has SD sqrt(2/4000) = 0.022; 0.1 is 4.5 sigma.
    assert counts.mean() == pytest.approx(2.0, abs=0.1)
    assert counts.var() == pytest.approx(2.0, abs=0.3)                          # variance = mean for Poisson
    rng2 = np.random.default_rng(9)
    assert [impl.simulate_fleet_failures(2.0 / M, 1.0e6, rng2) for _ in range(10)] == counts[:10].tolist()
    assert impl.simulate_fleet_failures(0.0, 1e6, rng) == 0
    with pytest.raises(ValueError):
        impl.simulate_fleet_failures(-1.0, 1e6, rng)


def test_credible_intervals_are_calibrated_under_the_prior():
    # Binomial SD of a coverage estimate from 4000 trials: sqrt(0.9 * 0.1 / 4000) = 0.0047 (0.02 is 4 sigma);
    # for level 0.5, sqrt(0.25 / 4000) = 0.0079 (0.035 is 4.4 sigma).
    assert impl.credible_interval_coverage(A0, B0, TF, 4000, 0.9, seed=5) == pytest.approx(0.9, abs=0.02)
    assert impl.credible_interval_coverage(A0, B0, 1.0e6, 4000, 0.5, seed=6) == pytest.approx(0.5, abs=0.035)
    assert impl.credible_interval_coverage(A0, B0, TF, 500, 0.9, seed=5) == impl.credible_interval_coverage(A0, B0, TF, 500, 0.9, seed=5)
    with pytest.raises(ValueError):
        impl.credible_interval_coverage(A0, B0, TF, 0, 0.9, seed=5)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "failure_rate_posterior.png").stat().st_size > 1000
