"""Public tests for P02.15 Probability Fundamentals. They assert properties, with exact counting and
scipy.special / scipy.stats as the oracles.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-02/parameterized-cases/.

Statistical tolerances: every Monte Carlo assertion uses a fixed seed, so the result is deterministic;
the tolerance is nevertheless set at four standard errors so that a correct implementation with a
different (but valid) generator stream would still pass with probability 1 - 6e-5.
"""
import math

import numpy as np
import numpy.testing as npt
import pytest
from scipy.special import comb
from scipy.stats import binom

from conftest import load_impl

impl = load_impl(__file__)


# ---------------------------------------------------------------- counting

def test_binomial_coefficient_matches_scipy_and_pascals_rule():
    for n in range(0, 31):
        for k in range(0, n + 1):
            assert impl.n_choose_k(n, k) == comb(n, k, exact=True)
        assert sum(impl.n_choose_k(n, k) for k in range(n + 1)) == 2 ** n          # binomial theorem at a = b = 1
    for n in range(1, 25):
        for k in range(1, n):
            assert impl.n_choose_k(n, k) == impl.n_choose_k(n - 1, k - 1) + impl.n_choose_k(n - 1, k)
    assert impl.n_choose_k(5, -1) == 0 and impl.n_choose_k(5, 6) == 0
    assert impl.n_choose_k(1000, 500) == comb(1000, 500, exact=True)              # exact integers, no overflow
    assert impl.permutations_count(10, 3) == 720 and impl.permutations_count(4, 4) == 24
    with pytest.raises(ValueError):
        impl.n_choose_k(-1, 0)


def test_two_dice_sample_space_and_hand_counts():
    outcomes = impl.enumerate_dice(2)
    assert outcomes.shape == (36, 2)
    assert len({tuple(o) for o in outcomes}) == 36                                 # every outcome once
    assert impl.event_probability(outcomes, lambda o: o.sum(axis=1) == 7) == pytest.approx(6 / 36)
    assert impl.event_probability(outcomes, lambda o: o.sum(axis=1) >= 10) == pytest.approx(6 / 36)
    assert impl.event_probability(outcomes, lambda o: o[:, 0] == o[:, 1]) == pytest.approx(6 / 36)
    assert impl.event_probability(outcomes, lambda o: o.sum(axis=1) == 2) == pytest.approx(1 / 36)
    four = impl.enumerate_dice(4)
    assert impl.event_probability(four, lambda o: np.any(o == 6, axis=1)) == pytest.approx(1 - (5 / 6) ** 4)
    with pytest.raises(ValueError):
        impl.event_probability(outcomes, lambda o: np.ones(3, dtype=bool))


def test_dice_sum_pmf_sums_to_one_and_matches_enumeration():
    for n in (1, 2, 3):
        sums, probs = impl.dice_sum_pmf(n)
        assert sums[0] == n and sums[-1] == 6 * n
        assert probs.sum() == pytest.approx(1.0)
        npt.assert_allclose(probs, probs[::-1])                                      # symmetric about the centre
        outcomes = impl.enumerate_dice(n)
        for s, p in zip(sums, probs):
            assert p == pytest.approx(impl.event_probability(outcomes, lambda o, s=s: o.sum(axis=1) == s))
    sums, probs = impl.dice_sum_pmf(2, sides=4)
    assert sums.tolist() == list(range(2, 9)) and probs[sums == 5] == pytest.approx(4 / 16)


# ---------------------------------------------------------------- the rules

def test_addition_complement_and_conditional_rules():
    assert impl.complement(0.3) == pytest.approx(0.7)
    assert impl.addition_rule(0.3, 0.5, 0.0) == pytest.approx(0.8)                 # disjoint: axiom (iii)
    assert impl.addition_rule(0.3, 0.5, 0.15) == pytest.approx(0.65)
    assert impl.addition_rule(0.3, 0.3, 0.3) == pytest.approx(0.3)                 # A = B
    assert impl.conditional_probability(0.2, 0.2) == pytest.approx(1.0)            # P(A | A) = 1
    assert impl.conditional_probability(0.15, 0.5) == pytest.approx(0.3)
    for bad in (impl.complement, ):
        with pytest.raises(ValueError):
            bad(1.3)
    with pytest.raises(ValueError):
        impl.conditional_probability(0.1, 0.0)
    with pytest.raises(ValueError):
        impl.addition_rule(0.2, 0.3, 0.4)                                          # P(A and B) > P(A)


def test_independence_on_dice_including_the_classic_trap():
    outcomes = impl.enumerate_dice(2)
    P = lambda f: impl.event_probability(outcomes, f)
    first1 = lambda o: o[:, 0] == 1
    second3 = lambda o: o[:, 1] == 3
    assert impl.are_independent(P(first1), P(second3), P(lambda o: first1(o) & second3(o)))
    sum7 = lambda o: o.sum(axis=1) == 7
    sum6 = lambda o: o.sum(axis=1) == 6
    assert impl.are_independent(P(sum7), P(first1), P(lambda o: sum7(o) & first1(o)))       # 1/36 = 1/6 * 1/6
    assert not impl.are_independent(P(sum6), P(first1), P(lambda o: sum6(o) & first1(o)))   # 1/36 != 5/216
    # Mutually exclusive events with positive probability are never independent.
    double = lambda o: o[:, 0] == o[:, 1]
    assert P(lambda o: sum7(o) & double(o)) == 0.0
    assert not impl.are_independent(P(sum7), P(double), 0.0)


def test_total_probability_and_bayes_posterior_properties():
    priors, lik = [0.02, 0.98], [0.90, 0.05]
    assert impl.total_probability(priors, lik) == pytest.approx(0.067)
    post = impl.bayes_posterior(priors, lik)
    assert post.sum() == pytest.approx(1.0)
    assert post[0] == pytest.approx(0.018 / 0.067)
    # An uninformative test (same likelihood on every branch) leaves the prior unchanged.
    npt.assert_allclose(impl.bayes_posterior(priors, [0.3, 0.3]), priors)
    # A perfectly specific test (no false alarms) makes a flag conclusive.
    assert impl.bayes_posterior(priors, [0.9, 0.0])[0] == pytest.approx(1.0)
    # Posterior is monotone in the prior.
    posts = [impl.bayes_posterior([q, 1 - q], lik)[0] for q in (0.001, 0.01, 0.1, 0.5)]
    assert np.all(np.diff(posts) > 0.0)
    with pytest.raises(ValueError):
        impl.total_probability([0.5, 0.4], [0.1, 0.2])                              # not a partition
    with pytest.raises(ValueError):
        impl.bayes_posterior(priors, [0.0, 0.0])                                    # impossible event


def test_inspection_tree_leaves_sum_to_one_and_match_bayes():
    tree = impl.inspection_tree(0.02, 0.90, 0.05)
    leaves = [tree[k] for k in ("defect_and_flag", "defect_and_pass", "good_and_flag", "good_and_pass")]
    assert sum(leaves) == pytest.approx(1.0)
    assert tree["p_flag"] == pytest.approx(tree["defect_and_flag"] + tree["good_and_flag"])
    assert tree["p_defect_given_flag"] == pytest.approx(impl.bayes_posterior([0.02, 0.98], [0.90, 0.05])[0])
    assert tree["p_defect_given_pass"] == pytest.approx(impl.bayes_posterior([0.02, 0.98], [0.10, 0.95])[0])
    assert tree["p_defect_given_pass"] < 0.02 < tree["p_defect_given_flag"]        # a flag raises, a pass lowers
    # Total probability of "defect" is conserved: P(D and F) + P(D and pass) = prior.
    assert tree["defect_and_flag"] + tree["defect_and_pass"] == pytest.approx(0.02)


# ---------------------------------------------------------------- reliability

def test_series_and_parallel_limits():
    ps = [0.995, 0.998, 0.990, 0.997, 0.999]
    s, p = impl.series_reliability(ps), impl.parallel_reliability(ps)
    assert s == pytest.approx(np.prod(ps)) and s <= min(ps)
    assert p == pytest.approx(1 - np.prod(1 - np.array(ps))) and p >= max(ps)
    assert impl.series_reliability([0.7]) == impl.parallel_reliability([0.7]) == pytest.approx(0.7)
    assert impl.series_reliability([1.0, 1.0]) == 1.0 and impl.parallel_reliability([0.0, 0.0]) == 0.0
    with pytest.raises(ValueError):
        impl.series_reliability([0.5, 1.2])


@pytest.mark.parametrize("n,p", [(3, 0.98), (5, 0.9), (7, 0.6)])
def test_k_of_n_matches_binomial_tail_and_reduces_to_series_and_parallel(n, p):
    for k in range(1, n + 1):
        assert impl.k_of_n_reliability(n, k, p) == pytest.approx(binom.sf(k - 1, n, p), rel=1e-12)
    assert impl.k_of_n_reliability(n, 1, p) == pytest.approx(impl.parallel_reliability([p] * n))
    assert impl.k_of_n_reliability(n, n, p) == pytest.approx(impl.series_reliability([p] * n))
    rel = [impl.k_of_n_reliability(n, k, p) for k in range(1, n + 1)]
    assert np.all(np.diff(rel) <= 1e-15)                                           # demanding more units can only hurt
    rel_p = [impl.k_of_n_reliability(n, 2, q) for q in (0.5, 0.7, 0.9, 0.99)]
    assert np.all(np.diff(rel_p) > 0.0)                                            # better units, better system
    assert sum(impl.binomial_probability(k, n, p) for k in range(n + 1)) == pytest.approx(1.0)
    with pytest.raises(ValueError):
        impl.k_of_n_reliability(n, 0, p)


def test_two_of_three_voting_value_and_crossover():
    assert impl.k_of_n_reliability(3, 2, 0.98) == pytest.approx(0.998816, abs=1e-6)
    # Voting beats a single unit only above p = 1/2 (3p^2 - 2p^3 > p).
    assert impl.k_of_n_reliability(3, 2, 0.8) > 0.8
    assert impl.k_of_n_reliability(3, 2, 0.4) < 0.4
    assert impl.k_of_n_reliability(3, 2, 0.5) == pytest.approx(0.5)


# ---------------------------------------------------------------- Monte Carlo

def test_proportion_standard_error_shape():
    assert impl.proportion_standard_error(0.0, 100) == 0.0 and impl.proportion_standard_error(1.0, 100) == 0.0
    assert impl.proportion_standard_error(0.5, 100) == pytest.approx(0.05)
    assert impl.proportion_standard_error(0.5, 400) == pytest.approx(0.025)        # halves when N quadruples
    with pytest.raises(ValueError):
        impl.proportion_standard_error(0.5, 0)


def test_monte_carlo_dice_sum_within_four_standard_errors():
    # Tolerance: 4 SE. With a fixed seed the outcome is deterministic; 4 SE keeps the test honest for any
    # valid generator stream (two-sided 4-sigma miss probability 6e-5).
    p_hat, se = impl.monte_carlo_dice_sum(2, 7, 100_000, seed=15)
    assert se == pytest.approx(math.sqrt(p_hat * (1 - p_hat) / 100_000))
    assert abs(p_hat - 1 / 6) < 4 * se
    assert se < 0.0015
    # Generic engine: P(at least one six in four rolls).
    p_hat2, se2 = impl.monte_carlo_probability(
        lambda rng, n: rng.integers(1, 7, size=(n, 4)), lambda r: np.any(r == 6, axis=1), 50_000, seed=3)
    assert abs(p_hat2 - (1 - (5 / 6) ** 4)) < 4 * se2
    # Same seed, same answer: the generator is seeded inside the function.
    assert impl.monte_carlo_dice_sum(2, 7, 1000, seed=7) == impl.monte_carlo_dice_sum(2, 7, 1000, seed=7)
    with pytest.raises(ValueError):
        impl.monte_carlo_probability(lambda rng, n: rng.integers(1, 7, size=n), lambda r: r[:3] == 1, 10, seed=0)


def test_monte_carlo_error_falls_as_one_over_sqrt_n():
    ns = [100, 400, 1600, 6400, 25600]
    rms = impl.monte_carlo_error_study(ns, n_repeats=20, seed=15)
    assert rms.shape == (5,) and np.all(rms > 0.0)
    slope = np.polyfit(np.log(ns), np.log(rms), 1)[0]
    # Tolerance 0.15: each RMS value is estimated from 20 repeats, so it carries about 16 % scatter
    # (SE of an SD ~ 1/sqrt(2*20)), which propagates to roughly 0.08 in a five-point log-log slope.
    assert slope == pytest.approx(-0.5, abs=0.15)
    predicted = np.sqrt((1 / 6) * (5 / 6) / np.array(ns))
    assert np.all(rms < 2.0 * predicted) and np.all(rms > 0.5 * predicted)


# ---------------------------------------------------------------- dice to dispersion

def test_walk_endpoint_pmf_is_binomial_and_symmetric():
    for n in (1, 5, 20):
        pos, probs = impl.walk_endpoint_pmf(n, 15.0)
        assert probs.sum() == pytest.approx(1.0)
        npt.assert_allclose(probs, binom.pmf(np.arange(n + 1), n, 0.5), rtol=1e-12)
        npt.assert_allclose(pos, 15.0 * (2 * np.arange(n + 1) - n))
        npt.assert_allclose(probs, probs[::-1])
        assert np.sum(pos * probs) == pytest.approx(0.0, abs=1e-12)                 # mean zero
        assert np.sum(pos ** 2 * probs) == pytest.approx(n * 15.0 ** 2)             # variance n d^2


def test_simulated_landings_have_binomial_moments():
    n, d, N = 20, 15.0, 20_000
    landings = impl.random_walk_landings(n, d, N, seed=2015)
    assert landings.shape == (N, 2)
    # Every offset is on the grid d (2h - n).
    assert np.all(np.isclose(((landings / d) + n) % 2, 0.0))
    sigma = d * math.sqrt(n)
    # Mean within 4 SE of zero (SE of a mean = sigma / sqrt(N)); variance within 5 % of n d^2
    # (SE of a variance estimate is about sigma^2 sqrt(2/N) = 1 %).
    assert np.all(np.abs(landings.mean(axis=0)) < 4 * sigma / math.sqrt(N))
    npt.assert_allclose(landings.var(axis=0, ddof=1), sigma ** 2, rtol=0.05)
    # The two axes are independent: correlation within 4 / sqrt(N).
    assert abs(np.corrcoef(landings.T)[0, 1]) < 4 / math.sqrt(N)
    with pytest.raises(ValueError):
        impl.random_walk_landings(0, d, 10)


def test_exact_recovery_probability_against_simulation_and_limits():
    n, d = 20, 15.0
    landings = impl.random_walk_landings(n, d, 20_000, seed=2015)
    for R in (50.0, 100.0, 150.0):
        p_exact = impl.probability_within_radius_exact(n, d, R)
        p_mc, se = impl.probability_within_radius_mc(landings, R)
        assert abs(p_mc - p_exact) < 4 * se                                        # 4 SE, see module docstring
    radii = np.linspace(0.0, 500.0, 51)
    ps = [impl.probability_within_radius_exact(n, d, R) for R in radii]
    assert np.all(np.diff(ps) >= -1e-15)                                           # non-decreasing in R
    assert ps[-1] == pytest.approx(1.0)                                             # R >= sqrt(2) n d catches everything
    assert impl.probability_within_radius_exact(n, d, 0.0) == pytest.approx(binom.pmf(n // 2, n, 0.5) ** 2)
    assert impl.probability_within_radius_exact(21, d, 0.5 * d) == 0.0             # odd n never lands on the origin
    # A grid point exactly on the circle counts as inside.
    assert impl.probability_within_radius_exact(2, 1.0, 2.0) == pytest.approx(
        impl.probability_within_radius_exact(2, 1.0, 2.0 + 1e-9))


def test_gaussian_formula_converges_to_the_exact_count_as_layers_grow():
    d = 15.0
    gaps = []
    for n in (4, 10, 20, 50, 100):
        sigma = d * math.sqrt(n)
        exact = impl.probability_within_radius_exact(n, d, 1.5 * sigma)
        gauss = impl.gaussian_radius_probability(n, d, 1.5 * sigma)
        assert gauss == pytest.approx(1 - math.exp(-1.125))
        gaps.append(exact - gauss)
    assert np.all(np.array(gaps) > 0.0)                                            # binomial tails are lighter
    assert gaps[-1] < 0.1 * gaps[0]                                                # and the gap closes


def test_radius_for_probability_is_consistent_with_the_exact_count():
    n, d = 20, 15.0
    for prob in (0.5, 0.9, 0.99):
        R = impl.radius_for_probability_exact(n, d, prob)
        assert impl.probability_within_radius_exact(n, d, R) >= prob - 1e-12
        # The next grid radius inwards would not be enough.
        assert impl.probability_within_radius_exact(n, d, R - 1e-6) < prob
    assert impl.radius_for_probability_exact(n, d, 0.9) < impl.radius_for_probability_exact(n, d, 0.99)
    with pytest.raises(ValueError):
        impl.radius_for_probability_exact(n, d, 1.5)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "probability_fundamentals.png").stat().st_size > 1000
