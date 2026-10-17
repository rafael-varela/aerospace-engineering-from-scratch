"""Public tests for P01.01 Limits, Continuity & Derivatives. They assert properties, not magic numbers.

Run `pytest curriculum/phase-01/lesson-01-limits-and-derivatives` against your starter/main.py, or with
`AEFS_IMPL=solution` against the reference solution. Hidden validation cases live in
assessments/phase-01/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
import sympy as sp

from conftest import load_impl

impl = load_impl(__file__)

A0, CLMAX = 2.0 * np.pi, 1.4
S = impl.model_symbols()
alpha, a0, clmax = S["alpha"], S["a0"], S["clmax"]


# ---------------------------------------------------------------- the lift model

def test_lift_model_limits_and_symmetry():
    assert impl.lift_coefficient(0.0) == 0.0
    al = np.radians(np.array([2.0, 6.0, 10.0, 20.0, 40.0]))
    cl = impl.lift_coefficient(al)
    assert cl.shape == al.shape
    assert np.all(np.diff(cl) > 0.0) and np.all(cl < CLMAX)                  # monotone, below the plateau
    npt.assert_allclose(impl.lift_coefficient(-al), -cl, rtol=1e-14)         # symmetric section: odd in alpha
    assert impl.lift_coefficient(np.radians(1.0)) == pytest.approx(A0 * np.radians(1.0), rel=5e-3)  # thin-airfoil at small alpha
    assert impl.lift_coefficient(np.radians(60.0), A0, 1.2) == pytest.approx(1.2, rel=1e-3)        # far past stall: cl_max


def test_exact_slope_is_the_derivative_of_the_lift_model():
    assert impl.lift_curve_slope_exact(0.0) == pytest.approx(A0, rel=1e-14)
    al = np.radians(np.array([0.0, 5.0, 10.0, 15.0]))
    slope = impl.lift_curve_slope_exact(al)
    assert np.all(np.diff(slope) < 0.0)                                      # slope decays towards stall
    npt.assert_allclose(slope, impl.central_difference(impl.lift_coefficient, al, 1e-5), rtol=1e-8)


def test_symbolic_slope_matches_hand_derivative_and_numerics():
    expr = impl.symbolic_lift_slope()
    hand = a0 / sp.cosh(a0 * alpha / clmax) ** 2
    assert sp.simplify(expr - hand) == 0
    fn = impl.lambdify_lift_slope()
    al = np.radians(np.linspace(0.0, 20.0, 9))
    npt.assert_allclose(fn(al, A0, CLMAX), impl.lift_curve_slope_exact(al, A0, CLMAX), rtol=1e-13)
    npt.assert_allclose(fn(al, A0, CLMAX), impl.central_difference(lambda x: impl.lift_coefficient(x, A0, CLMAX), al, 1e-4), rtol=1e-6)


# ---------------------------------------------------------------- limits and continuity

def test_limit_sequence_of_sin_x_over_x_tends_to_one_from_both_sides():
    hs = 10.0 ** -np.arange(1, 8, dtype=float)
    f = lambda x: np.sin(x) / x
    right = impl.limit_sequence(f, 0.0, hs, "right")
    left = impl.limit_sequence(f, 0.0, hs, "left")
    assert np.all(np.diff(np.abs(right - 1.0)) < 0.0)                        # error shrinks monotonically
    assert abs(right[-1] - 1.0) < 1e-12 and abs(left[-1] - 1.0) < 1e-12
    npt.assert_allclose(np.abs(right[:3] - 1.0), hs[:3] ** 2 / 6.0, rtol=1e-2)  # the error is h^2/6: sin x = x - x^3/6 + ...
    with pytest.raises(ValueError):
        impl.limit_sequence(f, 0.0, hs, "up")
    with pytest.raises(ValueError):
        impl.limit_sequence(f, 0.0, np.array([0.0]), "right")


def test_one_minus_cos_over_x_squared_tends_to_one_half():
    hs = 10.0 ** -np.arange(1, 5, dtype=float)
    vals = impl.limit_sequence(lambda x: (1.0 - np.cos(x)) / x**2, 0.0, hs)
    assert abs(vals[-1] - 0.5) < 1e-7


def test_one_sided_limits_and_continuity_detector():
    step = lambda x: np.where(np.asarray(x) < 0.0, -1.0, 1.0)                # jump at 0
    left, right = impl.one_sided_limits(step, 0.0)
    assert left == pytest.approx(-1.0) and right == pytest.approx(1.0)
    assert not impl.is_continuous_at(step, 0.0)
    assert impl.is_continuous_at(np.abs, 0.0, tol=1e-5)                      # |x|: continuous (tol above the last probe step) but not differentiable
    assert impl.is_continuous_at(impl.isa_temperature_two_layer, 11000.0)


def test_isa_two_layer_temperature_is_continuous_with_a_kink():
    t = impl.isa_temperature_two_layer(np.array([0.0, 11000.0, 15000.0]))
    assert t[0] == pytest.approx(288.15) and t[1] == pytest.approx(216.65) and t[2] == pytest.approx(216.65)
    left, right = impl.one_sided_derivatives(impl.isa_temperature_two_layer, 11000.0, 10.0)
    assert left == pytest.approx(-0.0065, abs=1e-9) and right == pytest.approx(0.0, abs=1e-9)
    smooth = impl.one_sided_derivatives(impl.isa_temperature_two_layer, 5000.0, 10.0)
    assert smooth[0] == pytest.approx(smooth[1], abs=1e-9)                   # no kink inside a layer
    with pytest.raises(ValueError):
        impl.isa_temperature_two_layer(25000.0)


def test_difference_quotients_converge_to_the_derivative_at_first_order():
    hs = 10.0 ** -np.arange(1, 7, dtype=float)
    q = impl.difference_quotients(np.exp, 0.5, hs)
    err = np.abs(q - np.exp(0.5))
    assert np.all(np.diff(err) < 0.0) and err[-1] < 1e-5
    rates = np.log10(err[:-1] / err[1:])
    assert np.all(rates > 0.9) and np.all(rates < 1.1)                       # error ~ h: one decade per decade


# ---------------------------------------------------------------- finite differences

@pytest.mark.parametrize("method", ["forward_difference", "backward_difference", "central_difference",
                                    "five_point_difference", "richardson_extrapolation"])
def test_every_method_rejects_bad_steps_and_preserves_shape(method):
    fn = getattr(impl, method)
    for bad in (0.0, -1e-3, float("nan")):
        with pytest.raises(ValueError):
            fn(np.sin, 1.0, bad)
    x = np.linspace(0.0, 2.0, 12).reshape(3, 4)
    out = fn(np.sin, x, 1e-3)
    assert np.shape(out) == x.shape
    npt.assert_allclose(out, np.cos(x), atol=2e-3)


def test_central_difference_is_exact_for_quadratics_and_five_point_for_quartics():
    quad = lambda x: 3.0 * x**2 - 2.0 * x + 1.0
    assert abs(impl.central_difference(quad, 1.3, 0.1) - (6.0 * 1.3 - 2.0)) < 1e-12
    quart = lambda x: x**4 - x**3 + 0.5 * x
    assert abs(impl.five_point_difference(quart, 1.3, 0.1) - (4.0 * 1.3**3 - 3.0 * 1.3**2 + 0.5)) < 1e-12
    assert abs(impl.central_difference(quart, 1.3, 0.1) - (4.0 * 1.3**3 - 3.0 * 1.3**2 + 0.5)) > 1e-3  # but not central


def test_forward_and_backward_errors_are_opposite_and_first_order():
    f, df = np.exp, np.exp
    for h in (1e-2, 1e-3):
        ef = impl.forward_difference(f, 0.3, h) - df(0.3)
        eb = impl.backward_difference(f, 0.3, h) - df(0.3)
        assert ef > 0.0 and eb < 0.0                                         # f'' > 0: forward overshoots, backward undershoots
        assert ef == pytest.approx(0.5 * h * np.exp(0.3), rel=0.02)          # leading term (h/2) f''
        assert (ef + eb) == pytest.approx(0.0, abs=h**2)                     # they cancel to O(h^2)


def test_richardson_equals_five_point_with_half_the_step():
    x = np.array([0.2, 0.7, 1.9])
    npt.assert_allclose(impl.richardson_extrapolation(np.sin, x, 0.1), impl.five_point_difference(np.sin, x, 0.05), rtol=1e-12)


def test_observed_orders_forward_one_central_two_five_point_four():
    hs = np.array([2e-2, 1e-2, 5e-3])                   # halving steps, truncation-dominated for exp near 1
    errs = impl.error_vs_step(np.exp, np.exp, 1.0, hs)
    assert set(errs) >= {"forward", "central", "five_point", "richardson"}
    assert impl.observed_order(hs, errs["forward"]) == pytest.approx(1.0, abs=0.1)
    assert impl.observed_order(hs, errs["central"]) == pytest.approx(2.0, abs=0.1)
    assert impl.observed_order(hs, errs["five_point"]) == pytest.approx(4.0, abs=0.25)
    assert np.all(errs["five_point"] <= errs["central"]) and np.all(errs["central"] <= errs["forward"])


def test_observed_order_recovers_a_synthetic_slope_and_rejects_zero_errors():
    hs = np.array([0.1, 0.05, 0.025, 0.0125])
    assert impl.observed_order(hs, 3.7 * hs**2) == pytest.approx(2.0, abs=1e-9)
    with pytest.raises(ValueError):
        impl.observed_order(hs, np.array([1.0, 0.0, 1.0, 1.0]))


def test_rounding_takes_over_below_the_optimal_step():
    hs = 10.0 ** -np.arange(1, 13, dtype=float)
    err = impl.error_vs_step(np.sin, np.cos, 1.0, hs)["central"]
    h_opt = impl.optimal_step_central(np.sin(1.0), np.cos(1.0))
    assert 1e-6 < h_opt < 1e-4
    i_min = int(np.argmin(err))
    assert 1e-7 <= hs[i_min] <= 1e-4                                         # the measured minimum sits near h_opt
    assert err[-1] > 100.0 * err[i_min]                                      # h = 1e-12 is far worse than the optimum
    assert err[0] > 100.0 * err[i_min]                                       # and so is h = 1e-1
    with pytest.raises(ValueError):
        impl.optimal_step_central(0.0, 1.0)


def test_lift_slope_by_every_method_agrees_with_symbolic_in_its_sweet_spot():
    a8 = np.radians(8.0)
    exact = float(impl.lambdify_lift_slope()(a8, A0, CLMAX))
    assert abs(impl.central_difference(impl.lift_coefficient, a8, 1e-4) - exact) < 1e-6
    assert abs(impl.five_point_difference(impl.lift_coefficient, a8, 1e-2) - exact) < 1e-4
    assert abs(impl.forward_difference(impl.lift_coefficient, a8, 1e-4) - exact) < 2e-3
    assert abs(impl.forward_difference(impl.lift_coefficient, a8, 1e-4) - exact) > 1e-4   # forward is visibly worse


def test_gust_increment_is_the_derivative_used_as_a_sensitivity():
    rho, v, s, u = 1.225, 50.0, 15.0, 5.0
    a8 = np.radians(8.0)
    slope = impl.lift_curve_slope_exact(a8)
    dl = impl.gust_lift_increment(slope, rho, v, s, u)
    assert dl > 0.0 and impl.gust_lift_increment(slope, rho, v, s, 2.0 * u) == pytest.approx(2.0 * dl)   # linear in the gust
    # The derivative is a local tool: for a 1 m/s gust (0.02 rad) the linear estimate is within 10 % of the
    # true lift change, for a 5 m/s gust (0.1 rad, towards stall) it overshoots by more than 20 %.
    q = 0.5 * rho * v**2
    exact_small = q * s * (impl.lift_coefficient(a8 + 1.0 / v) - impl.lift_coefficient(a8))
    assert impl.gust_lift_increment(slope, rho, v, s, 1.0) == pytest.approx(exact_small, rel=0.10)
    exact_big = q * s * (impl.lift_coefficient(a8 + u / v) - impl.lift_coefficient(a8))
    assert dl > 1.2 * exact_big
    # The same gust from alpha = 0 produces more lift: the slope there is larger.
    assert impl.gust_lift_increment(impl.lift_curve_slope_exact(0.0), rho, v, s, u) > dl
    with pytest.raises(ValueError):
        impl.gust_lift_increment(slope, rho, 0.0, s, u)
