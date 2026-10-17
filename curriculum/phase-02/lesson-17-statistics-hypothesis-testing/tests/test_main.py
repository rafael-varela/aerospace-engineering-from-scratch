"""Public tests for P02.17 Statistics & Hypothesis Testing. They assert properties, with scipy.stats as the oracle.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-02/parameterized-cases/.

Statistical tolerances: every simulation uses a fixed seed, so results are deterministic; tolerances are
nevertheless set at four standard errors (two-sided miss probability 6e-5 for a different valid generator
stream) and each is justified in a comment where it appears.
"""
import math

import numpy as np
import pytest
from scipy import stats

from conftest import load_impl

impl = load_impl(__file__)

BATCH = np.array([465.22, 456.06, 445.52, 436.88, 429.26, 452.22, 442.27, 441.53])   # the lesson's synthetic batch
REF = np.array([478.1, 462.3, 471.9, 455.0, 480.2, 468.7, 474.4, 459.9, 483.3, 467.8])


# ---------------------------------------------------------------- statistics and intervals

def test_sample_statistics_match_numpy():
    assert impl.sample_mean(BATCH) == pytest.approx(BATCH.mean())
    assert impl.sample_variance(BATCH) == pytest.approx(BATCH.var(ddof=1))
    assert impl.standard_error(BATCH) == pytest.approx(BATCH.std(ddof=1) / math.sqrt(BATCH.size))
    assert impl.sample_variance(np.array([3.0, 3.0, 3.0])) == 0.0
    with pytest.raises(ValueError):
        impl.sample_variance(np.array([1.0]))
    b = impl.synthetic_batch(8, 452.0, 12.0, seed=17)
    assert b.shape == (8,) and np.allclose(b, impl.synthetic_batch(8, 452.0, 12.0, seed=17))   # seeded


def test_t_critical_values_and_interval_shape():
    assert impl.t_critical(0.95, 7) == pytest.approx(2.364624, rel=1e-6)
    assert impl.t_critical(0.95, 1e6) == pytest.approx(1.959964, rel=1e-5)        # t -> normal
    assert impl.t_critical(0.99, 7) > impl.t_critical(0.95, 7) > impl.t_critical(0.90, 7)
    lo, hi = impl.confidence_interval_mean(BATCH, 0.95)
    m, se = BATCH.mean(), BATCH.std(ddof=1) / math.sqrt(8)
    assert 0.5 * (lo + hi) == pytest.approx(m)
    assert 0.5 * (hi - lo) == pytest.approx(stats.t.ppf(0.975, 7) * se)
    lo99, hi99 = impl.confidence_interval_mean(BATCH, 0.99)
    assert lo99 < lo and hi99 > hi                                                 # more confidence, wider interval
    # Half-width scales as t*/sqrt(n) for fixed s: tile the data to quadruple n.
    big = np.tile(BATCH - BATCH.mean(), 4) * math.sqrt(31 / 28) + BATCH.mean()     # same s with n = 32
    assert np.std(big, ddof=1) == pytest.approx(BATCH.std(ddof=1))
    lo4, hi4 = impl.confidence_interval_mean(big, 0.95)
    assert (hi4 - lo4) / (hi - lo) == pytest.approx(0.5 * stats.t.ppf(0.975, 31) / stats.t.ppf(0.975, 7))
    with pytest.raises(ValueError):
        impl.t_critical(1.0, 7)


def test_interval_coverage_is_nominal_over_synthetic_batches():
    # 2000 batches: SE of the coverage = sqrt(0.95 * 0.05 / 2000) = 0.0049; tolerance 4 SE = 0.02.
    rng = np.random.default_rng(173)
    hits = 0
    for _ in range(2000):
        lo, hi = impl.confidence_interval_mean(rng.normal(452.0, 12.0, size=8), 0.95)
        hits += lo <= 452.0 <= hi
    assert hits / 2000 == pytest.approx(0.95, abs=0.02)


# ---------------------------------------------------------------- one-sample test

@pytest.mark.parametrize("alternative", ["two-sided", "greater", "less"])
def test_one_sample_t_matches_scipy_for_every_alternative(alternative):
    for mu0 in (440.0, 446.0, 470.0):
        r = impl.one_sample_t_test(BATCH, mu0, alternative)
        ref = stats.ttest_1samp(BATCH, mu0, alternative=alternative)
        assert r["t"] == pytest.approx(ref.statistic, rel=1e-10)
        assert r["p"] == pytest.approx(ref.pvalue, rel=1e-10)
        assert r["df"] == 7.0
        assert 0.0 <= r["p"] <= 1.0


def test_one_sample_t_properties_and_duality():
    r = impl.one_sample_t_test(BATCH, BATCH.mean(), "two-sided")
    assert r["t"] == pytest.approx(0.0, abs=1e-12) and r["p"] == pytest.approx(1.0)
    # Two-sided p is twice the smaller one-sided p.
    g, l = impl.one_sample_t_test(BATCH, 440.0, "greater")["p"], impl.one_sample_t_test(BATCH, 440.0, "less")["p"]
    assert g + l == pytest.approx(1.0)
    assert impl.one_sample_t_test(BATCH, 440.0, "two-sided")["p"] == pytest.approx(2 * min(g, l))
    # Shift and scale invariance of t.
    t0 = impl.one_sample_t_test(BATCH, 440.0, "greater")["t"]
    assert impl.one_sample_t_test(BATCH + 100.0, 540.0, "greater")["t"] == pytest.approx(t0)
    assert impl.one_sample_t_test(BATCH * 1e-3, 0.440, "greater")["t"] == pytest.approx(t0)
    # Test-interval duality: one-sided test at alpha rejects iff the one-sided lower bound exceeds mu0.
    alpha = 0.05
    lower_bound = BATCH.mean() - stats.t.ppf(1 - alpha, 7) * BATCH.std(ddof=1) / math.sqrt(8)
    for mu0 in np.linspace(430.0, 450.0, 41):
        rejects = impl.one_sample_t_test(BATCH, mu0, "greater")["p"] < alpha
        assert rejects == (lower_bound > mu0)
    assert impl.decide(0.01) == "reject H0" and impl.decide(0.2) == "fail to reject H0"
    with pytest.raises(ValueError):
        impl.one_sample_t_test(BATCH, 440.0, "bigger")
    with pytest.raises(ValueError):
        impl.decide(0.5, alpha=1.5)


# ---------------------------------------------------------------- two samples

def test_welch_matches_scipy_and_df_bounds():
    for alt in ("two-sided", "greater", "less"):
        r = impl.welch_t_test(BATCH, REF, alt)
        ref = stats.ttest_ind(BATCH, REF, equal_var=False, alternative=alt)
        assert r["t"] == pytest.approx(ref.statistic, rel=1e-10)
        assert r["p"] == pytest.approx(ref.pvalue, rel=1e-10)
    r = impl.welch_t_test(BATCH, REF)
    assert min(BATCH.size, REF.size) - 1 <= r["df"] <= BATCH.size + REF.size - 2
    assert r["diff"] == pytest.approx(BATCH.mean() - REF.mean())
    assert r["se"] == pytest.approx(math.sqrt(BATCH.var(ddof=1) / 8 + REF.var(ddof=1) / 10))
    # Antisymmetry and the identical-sample case.
    assert impl.welch_t_test(REF, BATCH)["t"] == pytest.approx(-r["t"])
    assert impl.welch_t_test(BATCH, BATCH + 0.0)["p"] == pytest.approx(1.0)


def test_permutation_test_agrees_with_welch_on_normal_data_and_has_finite_resolution():
    p_perm = impl.permutation_test(BATCH, REF, n_perm=5000, seed=17)
    p_welch = impl.welch_t_test(BATCH, REF)["p"]
    assert p_perm >= 1 / 5001
    assert p_perm < 0.05 and p_welch < 0.05                                        # same verdict
    assert p_perm / p_welch < 5.0                                                  # same order of magnitude
    # Samples from the same distribution give a large p.
    rng = np.random.default_rng(3)
    x, y = rng.normal(size=30), rng.normal(size=30)
    assert impl.permutation_test(x, y, 2000, seed=4) > 0.05
    assert impl.permutation_test(x, y, 2000, seed=4) == impl.permutation_test(x, y, 2000, seed=4)   # seeded
    assert impl.permutation_test(x, x, 500, seed=1) == pytest.approx(1.0)          # identical groups: every shuffle ties


# ---------------------------------------------------------------- proportions

def test_proportion_z_test_against_normal_and_exact_binomial():
    r = impl.proportion_z_test(5, 120, 0.02, "greater")
    assert r["p_hat"] == pytest.approx(5 / 120)
    assert r["z"] == pytest.approx((5 / 120 - 0.02) / math.sqrt(0.02 * 0.98 / 120))
    assert r["p"] == pytest.approx(stats.norm.sf(r["z"]))
    # Small n p0: the approximation is far from exact (the lesson's verdict flip).
    assert r["p"] < 0.05 < stats.binomtest(5, 120, 0.02, alternative="greater").pvalue
    # Large n p0: approximation within 0.06 of the exact test (continuity correction omitted).
    r2 = impl.proportion_z_test(60, 1000, 0.05, "greater")
    assert abs(r2["p"] - stats.binomtest(60, 1000, 0.05, alternative="greater").pvalue) < 0.06
    assert impl.proportion_z_test(50, 1000, 0.05, "two-sided")["p"] == pytest.approx(1.0)
    with pytest.raises(ValueError):
        impl.proportion_z_test(5, 4, 0.5)


# ---------------------------------------------------------------- errors, power, sample size

def test_type_i_error_rate_equals_alpha():
    # 4000 batches: SE = sqrt(0.05 * 0.95 / 4000) = 0.0034; tolerance 4 SE = 0.014.
    for n in (3, 8, 50):
        rate = impl.type_i_error_rate_mc(440.0, 12.0, n, 0.05, 4000, seed=170 + n, alternative="greater")
        assert rate == pytest.approx(0.05, abs=0.014)
    rate2 = impl.type_i_error_rate_mc(0.0, 1.0, 8, 0.10, 4000, seed=9, alternative="two-sided")
    assert rate2 == pytest.approx(0.10, abs=0.019)


def test_power_monte_carlo_matches_noncentral_t_and_is_monotone():
    for delta in (5.0, 10.0, 15.0, 20.0):
        pa = impl.power_analytic(delta, 12.0, 8, 0.05, "greater")
        pm = impl.power_monte_carlo(delta, 12.0, 8, 0.05, 4000, seed=171, alternative="greater")
        se = math.sqrt(pa * (1 - pa) / 4000)                                       # 4 SE tolerance
        assert abs(pm - pa) < 4 * se
    assert impl.power_analytic(0.0, 12.0, 8, 0.05, "greater") == pytest.approx(0.05)
    assert impl.power_analytic(0.0, 12.0, 8, 0.05, "two-sided") == pytest.approx(0.05)
    assert impl.power_analytic(100.0, 12.0, 8, 0.05, "greater") == pytest.approx(1.0, abs=1e-9)
    assert impl.power_analytic(-10.0, 12.0, 8, 0.05, "less") == pytest.approx(impl.power_analytic(10.0, 12.0, 8, 0.05, "greater"))
    powers_d = [impl.power_analytic(d, 12.0, 8, 0.05, "greater") for d in (2, 5, 10, 15)]
    powers_n = [impl.power_analytic(10.0, 12.0, n, 0.05, "greater") for n in (4, 8, 16, 32)]
    assert np.all(np.diff(powers_d) > 0) and np.all(np.diff(powers_n) > 0)
    assert impl.power_analytic(10.0, 12.0, 8, 0.05, "greater") == pytest.approx(0.6833, abs=5e-4)
    # Power depends on delta / sigma only.
    assert impl.power_analytic(20.0, 24.0, 8) == pytest.approx(impl.power_analytic(10.0, 12.0, 8))
    with pytest.raises(ValueError):
        impl.power_monte_carlo(10.0, 12.0, 1)


def test_required_sample_size_is_minimal():
    for target in (0.8, 0.9, 0.95):
        n = impl.required_sample_size(10.0, 12.0, 0.05, target, "greater")
        assert impl.power_analytic(10.0, 12.0, n, 0.05, "greater") >= target
        assert impl.power_analytic(10.0, 12.0, n - 1, 0.05, "greater") < target
    assert impl.required_sample_size(10.0, 12.0, 0.05, 0.9, "greater") == 14
    assert impl.required_sample_size(5.0, 12.0, 0.05, 0.9, "greater") > 4 * 14 * 0.8   # roughly 4x for half the margin
    with pytest.raises(ValueError):
        impl.required_sample_size(0.0, 12.0)


# ---------------------------------------------------------------- bootstrap

def test_bootstrap_interval_agrees_with_t_interval_for_large_normal_samples():
    rng = np.random.default_rng(5)
    x = rng.normal(450.0, 12.0, size=200)
    lo_b, hi_b = impl.bootstrap_ci_mean(x, n_boot=4000, conf=0.95, seed=6)
    lo_t, hi_t = impl.confidence_interval_mean(x, 0.95)
    assert (hi_b - lo_b) == pytest.approx(hi_t - lo_t, rel=0.10)
    assert abs(0.5 * (lo_b + hi_b) - x.mean()) < 0.1 * (hi_t - lo_t)
    assert lo_b < x.mean() < hi_b
    assert impl.bootstrap_ci_mean(x, 500, 0.95, seed=1) == impl.bootstrap_ci_mean(x, 500, 0.95, seed=1)
    # At n = 8 the percentile bootstrap is a little narrower than t (it cannot see unseen tails).
    lo8, hi8 = impl.bootstrap_ci_mean(BATCH, 5000, 0.95, seed=17)
    lo8t, hi8t = impl.confidence_interval_mean(BATCH, 0.95)
    assert hi8 - lo8 < hi8t - lo8t
    with pytest.raises(ValueError):
        impl.bootstrap_ci_mean(BATCH, 100, 1.0)


# ---------------------------------------------------------------- tolerance bounds

def test_tolerance_factor_exceeds_z_and_converges_and_covers():
    zp = stats.norm.ppf(0.90)
    ks = [impl.one_sided_tolerance_factor(n, 0.90, 0.95) for n in (3, 8, 30, 1000, 100000)]
    assert np.all(np.diff(ks) < 0.0) and ks[-1] > zp
    assert ks[-1] == pytest.approx(zp, abs=0.01)
    assert impl.one_sided_tolerance_factor(8, 0.90, 0.95) == pytest.approx(2.582, abs=2e-3)
    assert impl.one_sided_tolerance_factor(8, 0.99, 0.95) > impl.one_sided_tolerance_factor(8, 0.90, 0.95)
    assert impl.one_sided_tolerance_factor(8, 0.90, 0.99) > impl.one_sided_tolerance_factor(8, 0.90, 0.95)
    lb = impl.lower_tolerance_bound(BATCH, 0.90, 0.95)
    assert lb == pytest.approx(BATCH.mean() - ks[1] * BATCH.std(ddof=1))
    assert lb < impl.confidence_interval_mean(BATCH, 0.95)[0]                    # a tail bound sits below a mean bound
    # Coverage: 4000 trials, SE = sqrt(0.95 * 0.05 / 4000) = 0.0034, tolerance 4 SE.
    assert impl.tolerance_bound_coverage_mc(8, 0.90, 0.95, 4000, seed=172) == pytest.approx(0.95, abs=0.014)
    with pytest.raises(ValueError):
        impl.one_sided_tolerance_factor(1, 0.90, 0.95)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "acceptance_testing.png").stat().st_size > 1000
