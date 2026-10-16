"""Public tests for P03.01 Floating Point & Numerical Error. They assert the error model, not magic numbers.

Oracles: ``numpy.finfo``, ``numpy.spacing``, ``math.fsum`` (exactly rounded sum), ``fractions.Fraction``
(exact rational arithmetic) and ``numpy.var``. Run `pytest <lesson dir>` against your starter/main.py,
or `AEFS_IMPL=solution pytest <lesson dir>` against the reference solution.
"""
import math
from fractions import Fraction

import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(301)


# ---------------------------------------------------------------- representation

def test_machine_epsilon_matches_finfo():
    assert impl.machine_epsilon(np.float64) == np.finfo(np.float64).eps == 2.0 ** -52
    assert impl.machine_epsilon(np.float32) == np.finfo(np.float32).eps == 2.0 ** -23
    with pytest.raises(TypeError):
        impl.machine_epsilon(np.int64)


def test_float_components_decode_known_doubles():
    assert impl.float_components(1.0) == (0, 1023, 0)
    assert impl.float_components(-2.0) == (1, 1024, 0)
    s, e, f = impl.float_components(0.1)
    assert (s, e) == (0, 1019)                     # 0.1 = 1.6 x 2^-4
    # Rebuild the value from its fields and compare: (1 + f / 2^52) 2^(e - 1023).
    assert (1.0 + f / 2.0 ** 52) * 2.0 ** (e - 1023) == 0.1
    assert impl.float_components(0.0) == (0, 0, 0)


@pytest.mark.parametrize("x", [1.0, 0.1, 6.778e6, 1e300, 3e-310, -2.5])
def test_ulp_matches_numpy_spacing(x):
    assert impl.ulp(x) == np.spacing(abs(x))
    assert impl.ulp(0.0) == 2.0 ** -1074
    with pytest.raises(ValueError):
        impl.ulp(float("inf"))


def test_rounding_model_bounds_every_operation_by_unit_roundoff():
    u = 2.0 ** -53
    a = RNG.uniform(-10, 10, 200)
    b = RNG.uniform(-10, 10, 200)
    for op in "+-*/":
        deltas = np.array([impl.exact_rounding_error(x, y, op) for x, y in zip(a, b)])
        assert np.all(np.abs(deltas) <= u), op
    # And the model is tight: 0.1 + 0.2 is not 0.3, with a known exact rounding error.
    d = impl.exact_rounding_error(0.1, 0.2, "+")
    assert d == float((Fraction(0.1 + 0.2) - (Fraction(0.1) + Fraction(0.2))) / (Fraction(0.1) + Fraction(0.2)))
    assert 0 < abs(d) <= u
    with pytest.raises(ValueError):
        impl.exact_rounding_error(1.0, 1.0, "^")


def test_relative_error_and_digits():
    assert impl.relative_error(1.01, 1.0) == pytest.approx(0.01)
    assert impl.correct_digits(1.0 + 1e-6, 1.0) == pytest.approx(6.0, abs=1e-6)
    assert impl.correct_digits(1.0, 1.0) == 16.0
    npt.assert_allclose(impl.relative_error([2.0, 4.0], [1.0, 2.0]), [1.0, 1.0])
    with pytest.raises(ValueError):
        impl.relative_error(1.0, 0.0)


def test_fixed_point_truncation_is_exact_and_below_one_lsb():
    for bits in (8, 16, 24, 32):
        err = impl.fixed_point_truncation_error(0.1, bits)
        assert 0.0 < err < 2.0 ** -bits                         # truncation loses less than one LSB
        assert err == float(Fraction(1, 10) - Fraction(math.floor(Fraction(1, 10) * 2 ** bits), 2 ** bits))
    assert impl.fixed_point_truncation_error(0.5, 4) == 0.0       # 0.5 = 0.1 binary is exact
    assert impl.fixed_point_truncation_error(0.1, 24) > impl.fixed_point_truncation_error(0.1, 32) * 10
    with pytest.raises(ValueError):
        impl.fixed_point_truncation_error(0.1, 0)


# ---------------------------------------------------------------- cancellation pairs

def test_one_minus_cos_stable_form_keeps_eps_accuracy_while_naive_loses_digits():
    thetas = 10.0 ** -np.arange(1, 8, dtype=float)
    exact = np.array([float(2 * Fraction(math.sin(t / 2)) ** 2) for t in thetas])   # stable form in exact arithmetic
    stable = impl.one_minus_cos_stable(thetas)
    naive = impl.one_minus_cos_naive(thetas)
    assert np.all(impl.relative_error(stable, exact) < 4 * 2.0 ** -52)
    # The naive error follows eps / (theta^2 / 2) within a factor 10 and grows as theta shrinks.
    model = 2.0 ** -52 / (thetas ** 2 / 2.0)
    rel_naive = impl.relative_error(naive, exact)
    assert np.all(rel_naive < 10 * model)
    assert rel_naive[-1] > 1e-4 and rel_naive[0] < 1e-13
    # The naive form has lost every digit by theta = 1e-8 (1 - cos rounds to exactly 0).
    assert impl.one_minus_cos_naive(1e-8) == 0.0 and impl.one_minus_cos_stable(1e-8) == pytest.approx(5e-17, rel=1e-10)
    # Both forms agree where there is no cancellation.
    assert impl.one_minus_cos_naive(1.0) == pytest.approx(impl.one_minus_cos_stable(1.0), rel=1e-15)


def test_cw_coefficient_matches_small_angle_series_only_in_stable_form():
    n, t = 1.1314e-3, 1.0
    series = n * t * t * (1.0 - (n * t) ** 2 / 12.0 + (n * t) ** 4 / 360.0)    # 2(1 - cos nt)/n = n t^2 (1 - (nt)^2/12 + ...)
    assert impl.cw_position_coefficient(n, t, stable=True) == pytest.approx(series, rel=1e-14)
    assert impl.relative_error(impl.cw_position_coefficient(n, t, stable=False), series) > 1e-12
    # For a long step there is no cancellation and both agree.
    assert impl.cw_position_coefficient(n, 2000.0, False) == pytest.approx(impl.cw_position_coefficient(n, 2000.0, True), rel=1e-13)


def test_quadratic_roots_stable_form_satisfies_vieta_and_the_naive_one_does_not():
    a, b, c = 1.0, 1e8, 1.0
    small_s, big_s = impl.quadratic_roots_stable(a, b, c)
    small_n, big_n = impl.quadratic_roots_naive(a, b, c)
    assert small_s * big_s == pytest.approx(c / a, rel=1e-15)           # Vieta's product
    assert small_s + big_s == pytest.approx(-b / a, rel=1e-15)
    assert impl.relative_error(small_n, -1e-8) > 0.1                      # the textbook formula is 25 % off
    assert impl.relative_error(small_s, -1e-8) < 1e-15
    assert big_n == big_s
    # Ordinary coefficients: both agree with numpy.roots.
    for _ in range(20):
        p, q, r = RNG.uniform(-3, 3, 3)
        if q * q - 4 * p * r < 0 or abs(p) < 1e-3:
            continue
        ref = sorted(np.roots([p, q, r]).real, key=abs)
        npt.assert_allclose(impl.quadratic_roots_stable(p, q, r), ref, rtol=1e-10, atol=1e-12)
        npt.assert_allclose(impl.quadratic_roots_naive(p, q, r), ref, rtol=1e-8, atol=1e-10)
    # b = 0 and degenerate inputs.
    assert impl.quadratic_roots_stable(1.0, 0.0, -4.0) == (-2.0, 2.0)
    with pytest.raises(ValueError):
        impl.quadratic_roots_stable(0.0, 1.0, 1.0)
    with pytest.raises(ValueError):
        impl.quadratic_roots_naive(1.0, 0.0, 1.0)


def test_exp_series_is_fine_for_positive_x_and_fails_for_large_negative_x():
    for x in (0.0, 0.5, 1.0, 5.0, 20.0):
        assert impl.exp_series(x) == pytest.approx(math.exp(x), rel=1e-14)
        assert impl.exp_stable(-x) == pytest.approx(math.exp(-x), rel=1e-14)
    # Alternating series for e^-20: the largest term is ~ 20^20/20! ~ 4e7, so the absolute error is
    # ~ eps 4e7 ~ 1e-8 against an answer of 2e-9: zero correct digits.
    assert impl.relative_error(impl.exp_series(-20.0), math.exp(-20.0)) > 0.1
    assert impl.relative_error(impl.exp_series(-30.0), math.exp(-30.0)) > 1e3


def test_subtracting_absolute_positions_costs_the_ulp_of_the_large_number():
    r_station = np.array([6.778e6, 0.0, 0.0])
    sep = np.array([-50.0123, -1000.0456, 0.0])
    err64 = np.max(np.abs(impl.subtract_positions(r_station, sep, np.float64) - sep))
    err32 = np.max(np.abs(impl.subtract_positions(r_station, sep, np.float32) - sep))
    assert err64 <= impl.ulp(6.778e6)                     # about 1e-9 m: harmless
    assert 0.001 < err32 <= float(np.spacing(np.float32(6.778e6)))   # up to 0.5 m: a centimetre requirement is dead
    npt.assert_allclose(impl.subtract_positions(np.zeros(3), sep, np.float64), sep)


# ---------------------------------------------------------------- summation and variance

def test_sums_match_fsum_in_double_and_the_error_ordering_shows_in_single():
    x = RNG.uniform(0.0, 1.0, 20_000)
    exact = math.fsum(x.tolist())
    for fn in (impl.naive_sum, impl.kahan_sum, impl.pairwise_sum):
        assert fn(x) == pytest.approx(exact, rel=1e-12)
    x32 = x.astype(np.float32)
    exact32 = math.fsum(x32.astype(np.float64).tolist())
    e_naive = abs(impl.naive_sum(x32, np.float32) - exact32)
    e_pair = abs(impl.pairwise_sum(x32, np.float32) - exact32)
    e_kahan = abs(impl.kahan_sum(x32, np.float32) - exact32)
    assert e_kahan <= e_pair * 1.5 + 1e-3 and e_pair < e_naive / 5     # compensated <= pairwise << naive
    assert e_naive > 1e-3                                                # naive single: visibly wrong
    assert e_kahan <= 2 * 2.0 ** -24 * np.sum(np.abs(x32)) + np.spacing(np.float32(exact32))   # the Kahan bound
    assert impl.pairwise_sum(np.array([])) == 0.0 and impl.pairwise_sum(np.array([3.0])) == 3.0


def test_naive_sum_error_grows_with_n_while_kahan_does_not():
    errs_naive, errs_kahan = [], []
    for n in (2_000, 20_000, 200_000):
        x = RNG.uniform(0.0, 1.0, n).astype(np.float32)
        exact = math.fsum(x.astype(np.float64).tolist())
        errs_naive.append(abs(impl.naive_sum(x, np.float32) - exact))
        errs_kahan.append(abs(impl.kahan_sum(x, np.float32) - exact) / exact)
    # Naive error grows at least linearly in n over a factor 100 (the model is between sqrt(n) and n
    # growth in error for a sum whose magnitude also grows as n: the ratio here must exceed 10).
    assert errs_naive[2] / errs_naive[0] > 10
    assert max(errs_kahan) < 1e-6                                      # relative error stays ~ eps32


def test_variance_two_pass_matches_numpy_and_one_pass_collapses_for_a_large_mean():
    alt = 400_000.0 + 0.01 * RNG.standard_normal(5_000)
    ref = float(np.var(alt))
    assert impl.variance_two_pass(alt) == pytest.approx(ref, rel=1e-10)
    # One pass in double: E[x^2] ~ 1.6e11 with ulp ~ 3e-5 against a variance of 1e-4: no digits.
    assert impl.relative_error(impl.variance_one_pass(alt), ref) > 0.1
    # One pass in single precision is garbage by many orders of magnitude.
    assert abs(impl.variance_one_pass(alt, np.float32)) > 1e3 * ref
    # Shift invariance: subtracting the mean does not change the true variance, and both formulas agree then.
    centred = alt - alt.mean()
    assert impl.variance_one_pass(centred) == pytest.approx(impl.variance_two_pass(centred), rel=1e-10)


# ---------------------------------------------------------------- truncation versus rounding

def test_forward_difference_observed_order_is_one_then_minus_one():
    hs = 10.0 ** -np.arange(1, 17, dtype=float)
    errs = np.array([abs(impl.forward_difference(math.sin, 1.0, h) - math.cos(1.0)) for h in hs])
    # Truncation regime h = 1e-2 .. 1e-4: error = h |sin 1| / 2 + O(h^2), observed order 1.
    assert impl.observed_order(hs[1:4], errs[1:4]) == pytest.approx(1.0, abs=0.05)
    # Rounding regime h = 1e-10 .. 1e-15: error ~ eps / h, observed order -1 (rounding is noisy: +-0.4).
    assert impl.observed_order(hs[9:15], errs[9:15]) == pytest.approx(-1.0, abs=0.4)
    # The model bounds the observed error everywhere and predicts the best step to within a decade.
    bound = impl.fd_total_error_bound(hs, 1.0, 1.0)
    assert np.all(errs <= bound)
    h_star = impl.fd_optimal_step(1.0, 1.0)
    assert h_star == pytest.approx(2.0 * math.sqrt(2.0 ** -52))
    best = hs[np.argmin(errs)]
    assert 0.1 * h_star <= best <= 10 * h_star
    assert np.min(errs) < 1e-7                                          # about half the digits survive
    with pytest.raises(ValueError):
        impl.fd_optimal_step(0.0, 1.0)


def test_observed_order_recovers_a_known_power_law():
    hs = np.array([0.1, 0.05, 0.025, 0.0125])
    assert impl.observed_order(hs, 3.0 * hs ** 2) == pytest.approx(2.0, abs=1e-10)
    assert impl.observed_order(hs, 0.5 / hs) == pytest.approx(-1.0, abs=1e-10)
    with pytest.raises(ValueError):
        impl.observed_order([0.1], [1.0])
    with pytest.raises(ValueError):
        impl.observed_order([0.1, 0.0], [1.0, 1.0])


def test_condition_numbers_distinguish_problem_from_algorithm():
    assert impl.relative_condition_number(math.sqrt, lambda x: 0.5 / math.sqrt(x), 4.0) == pytest.approx(0.5)
    assert impl.relative_condition_number(math.exp, math.exp, 3.0) == pytest.approx(3.0)
    # 1 - cos x is well conditioned near 0 (kappa -> 2): the digits lost by the naive formula are the
    # algorithm's fault, not the problem's; log x near 1 is ill conditioned (kappa -> infinity).
    assert impl.relative_condition_number(lambda x: 1 - math.cos(x), math.sin, 1e-3) == pytest.approx(2.0, rel=1e-6)
    assert impl.relative_condition_number(math.log, lambda x: 1 / x, 1.0 + 1e-6) > 1e5
    with pytest.raises(ValueError):
        impl.relative_condition_number(math.sin, math.cos, 0.0)


def test_main_runs_and_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "floating_point_error.png").exists()
