"""Public tests for P00.14. They assert the GUM formulas on analytic cases, not magic numbers.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Nothing is written outside tmp_path.
"""
import math

import numpy as np
import pytest

from conftest import load_impl

impl = load_impl(__file__)


# ---------- Type A ----------

def test_type_a_mean_and_standard_error():
    mean, u = impl.type_a([1.0, 2.0, 3.0])
    assert mean == pytest.approx(2.0)
    assert u == pytest.approx(1.0 / math.sqrt(3.0))      # s = 1, n = 3


def test_type_a_of_identical_readings_is_zero():
    mean, u = impl.type_a(np.full(7, 4.2))
    assert mean == pytest.approx(4.2) and u == 0.0


def test_type_a_shrinks_as_one_over_sqrt_n():
    rng = np.random.default_rng(0)
    big = rng.normal(10.0, 0.5, 40_000)
    _, u_small = impl.type_a(big[:100])
    _, u_large = impl.type_a(big)
    # expected ratio sqrt(400) = 20; allow sampling scatter on the 100-point sample
    assert 14 < u_small / u_large < 28


def test_type_a_rejects_single_reading():
    with pytest.raises(ValueError):
        impl.type_a([1.0])


# ---------- Type B ----------

def test_type_b_resolution_is_step_over_sqrt_12():
    assert impl.type_b_resolution(1.0) == pytest.approx(1.0 / math.sqrt(12.0))
    assert impl.type_b_resolution(0.1) == pytest.approx(impl.type_b_half_width(0.05))


@pytest.mark.parametrize("dist, divisor", [("rectangular", math.sqrt(3)), ("triangular", math.sqrt(6)), ("normal", 2.0)])
def test_type_b_half_width_divisors(dist, divisor):
    assert impl.type_b_half_width(0.6, dist) == pytest.approx(0.6 / divisor)


def test_type_b_rejects_bad_input():
    with pytest.raises(ValueError):
        impl.type_b_half_width(-1.0)
    with pytest.raises(ValueError):
        impl.type_b_half_width(1.0, "cauchy")


def test_combine_is_root_sum_square():
    assert impl.combine(3.0, 4.0) == pytest.approx(5.0)
    assert impl.combine(0.0, 0.0, 2.5) == pytest.approx(2.5)


# ---------- Jacobian and propagation ----------

def test_jacobian_matches_analytic_partials():
    f = lambda p: p[0] ** 2 + 3.0 * p[1] - math.sin(p[2])
    x = np.array([2.0, 1.0, 0.5])
    expected = np.array([4.0, 3.0, -math.cos(0.5)])
    assert np.allclose(impl.jacobian_fd(f, x), expected, rtol=1e-6, atol=1e-9)


def test_jacobian_error_decreases_with_step():
    f = lambda p: math.exp(p[0])
    coarse = abs(impl.jacobian_fd(f, [1.0], rel_step=1e-2)[0] - math.e)
    fine = abs(impl.jacobian_fd(f, [1.0], rel_step=1e-3)[0] - math.e)
    assert fine < coarse / 50          # central difference is second order: 100x per decade


def test_propagation_sum_adds_in_quadrature():
    r = impl.propagate(lambda p: p[0] + p[1], [3.0, 5.0], [0.3, 0.4])
    assert r.value == pytest.approx(8.0)
    assert r.u == pytest.approx(0.5)


def test_propagation_difference_same_as_sum():
    r = impl.propagate(lambda p: p[0] - p[1], [3.0, 5.0], [0.3, 0.4])
    assert r.u == pytest.approx(0.5)


def test_propagation_product_adds_relative_uncertainties():
    a, b, ua, ub = 4.0, 2.5, 0.2, 0.1
    r = impl.propagate(lambda p: p[0] * p[1], [a, b], [ua, ub])
    rel = math.hypot(ua / a, ub / b)
    assert r.value == pytest.approx(a * b)
    assert r.u / r.value == pytest.approx(rel, rel=1e-6)


def test_propagation_quotient_adds_relative_uncertainties():
    a, b, ua, ub = 9.0, 3.0, 0.3, 0.06
    r = impl.propagate(lambda p: p[0] / p[1], [a, b], [ua, ub])
    rel = math.hypot(ua / a, ub / b)
    assert r.value == pytest.approx(3.0)
    assert r.u / r.value == pytest.approx(rel, rel=1e-6)


def test_propagation_scales_linearly_with_constant():
    r = impl.propagate(lambda p: 7.0 * p[0], [1.3], [0.2])
    assert r.u == pytest.approx(1.4)
    assert r.sensitivities[0] == pytest.approx(7.0)
    assert r.contributions[0] == pytest.approx(1.4)


def test_propagation_power_law_fall_time():
    # t = sqrt(2h/g): relative u_t = 0.5 * sqrt((u_h/h)^2 + (u_g/g)^2)
    h, g, uh, ug = 1.2, 9.81, 0.01, 0.02
    r = impl.propagate(lambda p: math.sqrt(2 * p[0] / p[1]), [h, g], [uh, ug])
    assert r.u / r.value == pytest.approx(0.5 * math.hypot(uh / h, ug / g), rel=1e-6)


def test_zero_input_uncertainty_gives_zero_output_uncertainty():
    r = impl.propagate(lambda p: p[0] ** 3 / p[1], [2.0, 4.0], [0.0, 0.0])
    assert r.u == 0.0


def test_fully_correlated_difference_cancels():
    corr = np.array([[1.0, 1.0], [1.0, 1.0]])
    r = impl.propagate(lambda p: p[0] - p[1], [10.0, 10.0], [0.3, 0.3], corr=corr)
    assert r.u == pytest.approx(0.0, abs=1e-9)
    r2 = impl.propagate(lambda p: p[0] + p[1], [10.0, 10.0], [0.3, 0.3], corr=corr)
    assert r2.u == pytest.approx(0.6)


def test_propagate_rejects_mismatched_or_negative():
    with pytest.raises(ValueError):
        impl.propagate(lambda p: p[0], [1.0, 2.0], [0.1])
    with pytest.raises(ValueError):
        impl.propagate(lambda p: p[0], [1.0], [-0.1])


def test_expanded_uncertainty():
    assert impl.expanded(0.25) == pytest.approx(0.5)
    assert impl.expanded(0.25, k=3.0) == pytest.approx(0.75)


# ---------- Reporting ----------

@pytest.mark.parametrize(
    "value, u, unit, expected",
    [
        (9.8123, 0.1234, "m/s^2", "9.81 ± 0.12 m/s^2"),
        (1234.4, 25.0, "", "1234 ± 25"),
        (0.012343, 0.00047, "kg", "0.01234 ± 0.00047 kg"),
        (3.14159, 0.0996, "", "3.14 ± 0.10"),
    ],
)
def test_format_measurement(value, u, unit, expected):
    assert impl.format_measurement(value, u, unit) == expected


def test_round_to_uncertainty_keeps_value_and_u_at_same_place():
    v, u, d = impl.round_to_uncertainty(9.8123, 0.1234)
    assert (v, u, d) == (9.81, 0.12, 2)
    with pytest.raises(ValueError):
        impl.round_to_uncertainty(1.0, 0.0)


# ---------- Calibration ----------

def test_calibration_fit_recovers_exact_line():
    ref = np.array([0.0, 5.0, 10.0, 20.0, 50.0])
    fit = impl.calibration_fit(ref, 1.02 * ref - 0.3)
    assert fit.slope == pytest.approx(1.02, rel=1e-9)
    assert fit.offset == pytest.approx(-0.3, abs=1e-9)
    assert np.allclose(fit.residuals, 0.0, atol=1e-9)
    assert fit.s_residual == pytest.approx(0.0, abs=1e-9)


def test_calibration_fit_with_noise_is_consistent():
    rng = np.random.default_rng(42)
    ref = np.linspace(0.0, 100.0, 11)
    reading = 0.98 * ref + 0.5 + rng.normal(0.0, 0.2, ref.size)
    fit = impl.calibration_fit(ref, reading)
    assert abs(fit.slope - 0.98) < 4 * fit.u_slope
    assert abs(fit.offset - 0.5) < 4 * fit.u_offset
    assert abs(fit.residuals.mean()) < 1e-9          # least squares with an intercept: residuals sum to zero
    assert 0.1 < fit.s_residual < 0.4


def test_apply_calibration_inverts_the_line():
    ref = np.array([0.0, 10.0, 20.0, 30.0])
    fit = impl.calibration_fit(ref, 1.05 * ref + 2.0)
    assert np.allclose(impl.apply_calibration(1.05 * ref + 2.0, fit), ref)
    assert impl.apply_calibration(12.5, fit) == pytest.approx((12.5 - 2.0) / 1.05)


def test_calibration_fit_rejects_too_few_points():
    with pytest.raises(ValueError):
        impl.calibration_fit([0.0, 1.0], [0.0, 1.0])
    with pytest.raises(ValueError):
        impl.calibration_fit([2.0, 2.0, 2.0], [1.0, 2.0, 3.0])


def test_corrected_value_uncertainty_is_at_least_the_reading_uncertainty():
    rng = np.random.default_rng(1)
    ref = np.linspace(0.0, 50.0, 6)
    fit = impl.calibration_fit(ref, 1.01 * ref + rng.normal(0.0, 0.05, ref.size))
    res = impl.corrected_value(20.0, 0.1, fit)
    assert res.value == pytest.approx((20.0 - fit.offset) / fit.slope)
    assert res.u >= 0.1 / fit.slope * 0.999
    assert res.contributions.shape == (3,)
