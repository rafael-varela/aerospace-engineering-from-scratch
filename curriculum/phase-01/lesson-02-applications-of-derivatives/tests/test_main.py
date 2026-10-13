"""Public tests for P01.02 Applications of Derivatives. They assert optimality, limits and agreement, not magic numbers.

Run `pytest curriculum/phase-01/lesson-02-applications-of-derivatives` against your starter/main.py, or with
`AEFS_IMPL=solution` against the reference solution. Hidden validation cases live in
assessments/phase-01/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

G0, RHO0 = 9.80665, 1.225
W, S, CD0 = 350.0 * G0, 10.5, 0.015
K = 1.0 / (np.pi * 0.9 * 20.0)


def central(f, x, h=1e-6):
    return (f(x + h) - f(x - h)) / (2.0 * h)


# ---------------------------------------------------------------- the polar and its derivatives

def test_induced_drag_factor_and_polar_limits():
    assert impl.induced_drag_factor(20.0, 0.9) == pytest.approx(K)
    assert impl.induced_drag_factor(40.0, 0.9) == pytest.approx(K / 2.0)        # doubling AR halves K
    with pytest.raises(ValueError):
        impl.induced_drag_factor(0.0, 0.9)
    assert impl.drag_polar(0.0, CD0, K) == CD0                                   # no lift, only parasite drag
    cl = np.linspace(0.1, 1.5, 15)
    assert np.all(np.diff(impl.drag_polar(cl, CD0, K)) > 0.0)                    # drag grows with lift
    assert impl.lift_to_drag(0.0, CD0, K) == 0.0


def test_ld_derivatives_match_finite_differences():
    cl = np.array([0.3, 0.7, 0.9, 1.2])
    npt.assert_allclose(impl.ld_derivative(cl, CD0, K), central(lambda c: impl.lift_to_drag(c, CD0, K), cl), rtol=1e-7)
    npt.assert_allclose(impl.ld_second_derivative(cl, CD0, K), central(lambda c: impl.ld_derivative(c, CD0, K), cl), rtol=1e-6)


def test_best_glide_is_a_maximum_of_l_over_d():
    cl_star = impl.best_glide_cl(CD0, K)
    assert cl_star == pytest.approx(np.sqrt(CD0 / K))
    assert impl.ld_derivative(cl_star, CD0, K) == pytest.approx(0.0, abs=1e-12)  # critical point
    assert impl.ld_second_derivative(cl_star, CD0, K) < 0.0                       # second-derivative test: maximum
    ld_max = impl.max_lift_to_drag(CD0, K)
    assert ld_max == pytest.approx(impl.lift_to_drag(cl_star, CD0, K))
    cl = np.linspace(0.05, 2.5, 500)
    assert np.all(impl.lift_to_drag(cl, CD0, K) <= ld_max + 1e-12)              # nothing beats the optimum
    # At best glide, induced drag equals parasite drag.
    assert K * cl_star**2 == pytest.approx(CD0)


def test_newton_finds_the_critical_point_inside_its_basin_and_fails_outside():
    for cl0 in (0.3, 0.6, 1.2):
        assert impl.best_glide_cl_newton(CD0, K, cl0) == pytest.approx(impl.best_glide_cl(CD0, K), rel=1e-10)
    # The basin ends at the inflection point of L/D, C_L = sqrt(3 cd0/k) (= the min-sink C_L), where
    # (L/D)'' = 0 and Newton's step blows up: from beyond it the iteration diverges.
    assert impl.ld_second_derivative(impl.min_sink_cl(CD0, K), CD0, K) == pytest.approx(0.0, abs=1e-12)
    with pytest.raises(RuntimeError):
        impl.best_glide_cl_newton(CD0, K, 2.0)


def test_min_sink_is_sqrt3_times_best_glide_and_minimises_sink():
    cl_ms = impl.min_sink_cl(CD0, K)
    assert cl_ms == pytest.approx(np.sqrt(3.0) * impl.best_glide_cl(CD0, K))
    cl = np.linspace(0.3, 2.5, 800)
    w = impl.sink_rate(cl, CD0, K, W, RHO0, S)
    assert np.all(w >= impl.sink_rate(cl_ms, CD0, K, W, RHO0, S) * (1.0 - 1e-3))   # exact-angle sink is minimal there to 0.1 %
    assert impl.sink_rate(cl_ms, CD0, K, W, RHO0, S) < impl.sink_rate(impl.best_glide_cl(CD0, K), CD0, K, W, RHO0, S)


def test_glide_geometry_and_scaling():
    cl_star = impl.best_glide_cl(CD0, K)
    gamma = impl.glide_angle(cl_star, CD0, K)
    assert np.tan(gamma) == pytest.approx(1.0 / impl.max_lift_to_drag(CD0, K))
    v = impl.airspeed(cl_star, W, RHO0, S, gamma)
    assert impl.airspeed(cl_star, 4.0 * W, RHO0, S, gamma) == pytest.approx(2.0 * v)     # V ~ sqrt(W)
    assert impl.airspeed(cl_star, W, RHO0 / 4.0, S, gamma) == pytest.approx(2.0 * v)     # V ~ 1/sqrt(rho)
    assert impl.airspeed(cl_star, W, RHO0, S) > v                                         # cos(gamma) < 1
    with pytest.raises(ValueError):
        impl.airspeed(0.0, W, RHO0, S)
    with pytest.raises(ValueError):
        impl.airspeed(-0.5, W, RHO0, S)


def test_small_angle_sink_agrees_with_exact_for_a_good_glider():
    cl = np.linspace(0.4, 1.6, 7)
    v = impl.airspeed(cl, W, RHO0, S)
    exact = impl.sink_rate(cl, CD0, K, W, RHO0, S)
    approx = impl.sink_rate_of_speed(v, CD0, K, W, RHO0, S)
    npt.assert_allclose(approx, exact, rtol=2e-3)                                 # gamma < 3 deg: cos ~ 1 to 0.15 %
    assert np.all(approx > exact)                                                 # small-angle slightly overestimates sink
    with pytest.raises(ValueError):
        impl.sink_rate_of_speed(0.0, CD0, K, W, RHO0, S)


def test_glide_range_is_linear_in_height_and_in_l_over_d():
    r = impl.glide_range(1000.0, CD0, K)
    assert r == pytest.approx(1000.0 * impl.max_lift_to_drag(CD0, K))
    assert impl.glide_range(2000.0, CD0, K) == pytest.approx(2.0 * r)
    assert impl.glide_range(1000.0, 4.0 * CD0, K) == pytest.approx(r / 2.0)     # (L/D)max ~ 1/sqrt(cd0)
    assert impl.glide_range(0.0, CD0, K) == 0.0
    with pytest.raises(ValueError):
        impl.glide_range(-1.0, CD0, K)


# ---------------------------------------------------------------- Newton's method itself

def test_newton_converges_quadratically_and_reports_failure():
    root, n = impl.newton(lambda x: x**2 - 2.0, lambda x: 2.0 * x, 1.0)
    assert root == pytest.approx(np.sqrt(2.0), rel=1e-14) and n <= 7           # quadratic: ~doubling digits each step
    with pytest.raises(RuntimeError):
        impl.newton(lambda x: x**2 + 1.0, lambda x: 2.0 * x, 0.5, maxiter=30)   # no real root
    with pytest.raises(ZeroDivisionError):
        impl.newton(lambda x: x**2 - 2.0, lambda x: 0.0, 1.0)


# ---------------------------------------------------------------- linear approximation

def test_linear_approximation_error_is_second_order_in_the_step():
    f = lambda c: impl.max_lift_to_drag(c, K)
    df = lambda c: impl.ld_max_sensitivity_cd0(c, K)
    assert df(CD0) == pytest.approx(central(f, CD0, 1e-7), rel=1e-6)
    assert df(CD0) < 0.0                                                          # more parasite drag, less glide
    errs = np.array([abs(impl.linear_approximation(f, df, CD0, d * CD0) - f(CD0 + d * CD0)) for d in (0.04, 0.02, 0.01)])
    rates = np.log2(errs[:-1] / errs[1:])
    assert np.all(rates > 1.9) and np.all(rates < 2.1)
    assert impl.linear_approximation(np.sin, np.cos, 0.0, 0.1) == pytest.approx(0.1)


# ---------------------------------------------------------------- related rates

def test_isa_density_and_its_gradient():
    assert impl.isa_density(0.0) == pytest.approx(RHO0)
    h = np.array([100.0, 2000.0, 5000.0, 9000.0])
    rho = impl.isa_density(h)
    assert np.all(np.diff(rho) < 0.0)
    npt.assert_allclose(impl.isa_density_gradient(h), central(impl.isa_density, h, 1.0), rtol=1e-7)
    assert impl.isa_density_gradient(0.0) == pytest.approx(-RHO0 * 4.2559 * 0.0065 / 288.15, rel=1e-4)   # rho n L / T at sea level
    assert np.all(impl.isa_density_gradient(h) < 0.0)
    with pytest.raises(ValueError):
        impl.isa_density(12000.0)


def test_best_glide_speed_rate_matches_a_finite_difference_in_time():
    hdot = -0.8
    h0, dt = 3000.0, 1.0
    rate = impl.best_glide_speed_rate(h0, hdot, CD0, K, W, S)
    fd = (impl.best_glide_speed(h0 + hdot * dt, CD0, K, W, S) - impl.best_glide_speed(h0 - hdot * dt, CD0, K, W, S)) / (2.0 * dt)
    assert rate == pytest.approx(fd, rel=1e-6)
    assert rate < 0.0                                                             # descending into denser air: slow down
    assert impl.best_glide_speed_rate(h0, 0.0, CD0, K, W, S) == 0.0               # level: no change
    assert impl.best_glide_speed(0.0, CD0, K, W, S) == pytest.approx(impl.airspeed(impl.best_glide_cl(CD0, K), W, RHO0, S))


# ---------------------------------------------------------------- headwind optimum

def test_headwind_optimum_reduces_to_best_glide_in_still_air_and_shifts_correctly():
    v_still = impl.headwind_best_speed(0.0, CD0, K, W, RHO0, S)
    assert v_still == pytest.approx(impl.airspeed(impl.best_glide_cl(CD0, K), W, RHO0, S), rel=1e-9)
    v_head = impl.headwind_best_speed(10.0, CD0, K, W, RHO0, S)
    v_tail = impl.headwind_best_speed(-10.0, CD0, K, W, RHO0, S)
    assert v_head > v_still > v_tail                                              # speed up into wind, slow down with it
    v_min_sink = impl.airspeed(impl.min_sink_cl(CD0, K), W, RHO0, S)
    assert v_tail > v_min_sink                                                    # never slower than minimum sink
    for vw, vb in ((10.0, v_head), (-10.0, v_tail), (0.0, v_still)):
        g_best = impl.ground_glide_ratio(vb, vw, CD0, K, W, RHO0, S)
        v = np.linspace(0.6 * vb, 1.6 * vb, 400)
        assert np.all(impl.ground_glide_ratio(v, vw, CD0, K, W, RHO0, S) <= g_best + 1e-9)   # it is the maximum
    assert impl.ground_glide_ratio(v_still, 0.0, CD0, K, W, RHO0, S) == pytest.approx(impl.max_lift_to_drag(CD0, K), rel=1e-9)
