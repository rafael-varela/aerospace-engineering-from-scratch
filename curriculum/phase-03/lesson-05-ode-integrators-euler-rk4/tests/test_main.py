"""Public tests for P03.05 ODE Integrators I: Euler to RK4. They assert the Butcher order conditions, exactness
classes, local and global convergence orders, linear-stability limits and the two-body orbit results, with
scipy.integrate.solve_ivp and the Kepler closed form as oracles.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import math

import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import solve_ivp

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(305)
METHODS = ["euler", "midpoint", "heun", "kutta3", "rk4"]
ORDERS = {"euler": 1, "midpoint": 2, "heun": 2, "kutta3": 3, "rk4": 4}


def f_harm(t, y):
    return np.array([y[1], -y[0]])


def ex_harm(t):
    return np.array([math.cos(t), -math.sin(t)])


def f_gauss(t, y):
    return -2.0 * t * y


def ex_gauss(t):
    return np.array([math.exp(-t * t)])


# ---------------------------------------------------------------- tableaus and steps

def test_tableaus_satisfy_the_order_conditions():
    for name, tab in impl.TABLEAUS.items():
        s = tab.b.size
        assert tab.a.shape == (s, s) and tab.c.shape == (s,)
        assert np.allclose(np.triu(tab.a), 0.0), f"{name}: tableau must be explicit"
        npt.assert_allclose(tab.c, tab.a.sum(axis=1), atol=1e-15)
        assert tab.b.sum() == pytest.approx(1.0, abs=1e-15)                      # order 1
        if tab.order >= 2:
            assert tab.b @ tab.c == pytest.approx(0.5, abs=1e-15)                # order 2
        if tab.order >= 3:
            assert tab.b @ tab.c ** 2 == pytest.approx(1.0 / 3.0, abs=1e-15)
            assert tab.b @ (tab.a @ tab.c) == pytest.approx(1.0 / 6.0, abs=1e-15)
        if tab.order >= 4:
            assert tab.b @ tab.c ** 3 == pytest.approx(0.25, abs=1e-15)
            assert tab.b @ (tab.c * (tab.a @ tab.c)) == pytest.approx(1.0 / 8.0, abs=1e-15)
            assert tab.b @ (tab.a @ tab.c ** 2) == pytest.approx(1.0 / 12.0, abs=1e-15)
            assert tab.b @ (tab.a @ (tab.a @ tab.c)) == pytest.approx(1.0 / 24.0, abs=1e-15)
    assert [impl.stages(m) for m in METHODS] == [1, 2, 2, 3, 4]


def test_hand_written_steps_equal_the_generic_tableau_step():
    y = RNG.normal(size=4)
    t, h = 0.7, 0.13
    f = lambda t, y: np.array([y[1], -math.sin(y[0]) + 0.1 * t, y[3], -y[2] ** 3])
    for name, step in [("euler", impl.euler_step), ("midpoint", impl.midpoint_step), ("heun", impl.heun_step), ("rk4", impl.rk4_step)]:
        npt.assert_allclose(step(f, t, y, h), impl.explicit_rk_step(f, t, y, h, impl.TABLEAUS[name]), rtol=1e-14, atol=1e-15)


def test_one_step_on_y_prime_equals_y_is_the_truncated_exponential():
    f = lambda t, y: y
    h = 0.5
    y0 = np.array([1.0])
    for m in METHODS:
        p = ORDERS[m]
        expected = sum(h ** k / math.factorial(k) for k in range(p + 1))
        assert impl.STEPPERS[m](f, 0.0, y0, h)[0] == pytest.approx(expected, rel=1e-15)
    # the local error of RK4 is then exactly the remainder of the series: about h^5 / 120
    err = impl.rk4_step(f, 0.0, y0, h)[0] - math.exp(h)
    assert abs(err) == pytest.approx(h ** 5 / 120.0, rel=0.15)


def test_rk4_is_exact_for_quartic_polynomials_in_t_and_euler_for_linear_ones():
    # y' = g(t) with polynomial g: RK4 reduces to Simpson's rule, exact for cubic g (quartic y)
    for coeffs in [np.array([1.0, -2.0, 3.0, 0.5]), RNG.normal(size=4)]:
        g = lambda t, y, c=coeffs: np.array([np.polyval(c, t)])
        G = lambda t, c=coeffs: np.polyval(np.polyint(c), t)
        y1 = impl.rk4_step(g, 0.3, np.array([G(0.3)]), 0.8)[0]
        assert y1 == pytest.approx(G(1.1), rel=1e-13)
    g1 = lambda t, y: np.array([3.0])
    assert impl.euler_step(g1, 0.0, np.array([1.0]), 0.4)[0] == pytest.approx(2.2, rel=1e-15)
    # ... and the quintic breaks RK4
    g5 = lambda t, y: np.array([5.0 * t ** 4])
    assert impl.rk4_step(g5, 0.0, np.array([0.0]), 1.0)[0] != pytest.approx(1.0, rel=1e-6)


# ---------------------------------------------------------------- integrate

def test_integrate_lands_on_t1_with_the_right_shapes_and_adjusts_h():
    t, Y = impl.integrate(f_harm, (0.0, 1.0), ex_harm(0.0), 0.3, "rk4")
    assert t.shape == (5,) and Y.shape == (5, 2)
    assert t[-1] == 1.0 and t[0] == 0.0
    npt.assert_allclose(np.diff(t), 0.25, rtol=1e-14)
    npt.assert_allclose(Y[0], ex_harm(0.0))
    t2, Y2 = impl.integrate(f_harm, (0.0, 1.0), ex_harm(0.0), 0.25, "rk4")
    npt.assert_allclose(Y2, Y, rtol=1e-14)
    for bad in [dict(h=0.0), dict(h=-1.0), dict(method="rk5")]:
        with pytest.raises(ValueError):
            impl.integrate(f_harm, (0.0, 1.0), ex_harm(0.0), bad.get("h", 0.1), bad.get("method", "rk4"))
    with pytest.raises(ValueError):
        impl.integrate(f_harm, (1.0, 0.0), ex_harm(0.0), 0.1)


def test_rk4_reproduces_the_exponential_and_matches_solve_ivp():
    t, Y = impl.integrate(lambda t, y: y, (0.0, 2.0), np.array([1.0]), 0.01, "rk4")
    assert abs(Y[-1, 0] - math.exp(2.0)) < 1e-8
    # a nonlinear problem without a closed form: the Van der Pol oscillator against DOP853 at tight tolerance
    vdp = lambda t, y: np.array([y[1], 2.0 * (1.0 - y[0] ** 2) * y[1] - y[0]])
    y0 = np.array([2.0, 0.0])
    ref = solve_ivp(vdp, (0.0, 5.0), y0, method="DOP853", rtol=1e-12, atol=1e-12).y[:, -1]
    _, Y = impl.integrate(vdp, (0.0, 5.0), y0, 1e-3, "rk4")
    npt.assert_allclose(Y[-1], ref, rtol=1e-8, atol=1e-8)


# ---------------------------------------------------------------- orders

def test_local_error_order_is_p_plus_one():
    hs = np.array([0.4, 0.2, 0.1])
    for m in METHODS:
        errs = [impl.local_error(f_harm, ex_harm, 0.3, h, m) for h in hs]
        assert impl.observed_order(hs, errs) == pytest.approx(ORDERS[m] + 1, abs=0.1), m
    # on y' = -2 t y the RK4 leading error term changes sign near t = 0.4, so its apparent order there is
    # meaningless; the lower methods still show p + 1 to within 0.3
    for m in ["euler", "midpoint", "heun", "kutta3"]:
        errs = [impl.local_error(f_gauss, ex_gauss, 0.4, h, m) for h in hs]
        assert impl.observed_order(hs, errs) == pytest.approx(ORDERS[m] + 1, abs=0.35), m


def test_global_error_order_is_p_on_the_harmonic_oscillator():
    hs = 2.0 * math.pi / np.array([100, 200, 400])
    for m in METHODS:
        errs = [impl.global_error(f_harm, ex_harm, (0.0, 2.0 * math.pi), h, m) for h in hs]
        assert impl.observed_order(hs, errs) == pytest.approx(ORDERS[m], abs=0.12), m
    # and the methods rank as their orders say at h = 2 pi / 400
    e = [impl.global_error(f_harm, ex_harm, (0.0, 2.0 * math.pi), hs[-1], m) for m in ["euler", "heun", "kutta3", "rk4"]]
    assert e[0] > e[1] > e[2] > e[3]
    assert e[3] < 1e-8 and e[0] > 1e-2


def test_observed_order_and_step_for_tolerance():
    hs = np.array([0.1, 0.05, 0.025])
    assert impl.observed_order(hs, 3.0 * hs ** 4) == pytest.approx(4.0, abs=1e-10)
    with pytest.raises(ValueError):
        impl.observed_order(hs[:1], hs[:1])
    with pytest.raises(ValueError):
        impl.observed_order(hs, np.array([1.0, 0.0, 1.0]))
    assert impl.step_for_tolerance(0.1, 1e-4, 4, 1e-8) == pytest.approx(0.01)
    assert impl.step_for_tolerance(0.1, 1e-4, 1, 1e-8) == pytest.approx(1e-5)
    with pytest.raises(ValueError):
        impl.step_for_tolerance(0.1, 1e-4, 0, 1e-8)
    # the prediction works on a real integration (RK4, harmonic oscillator): within 30 % of the target error
    h0 = 2.0 * math.pi / 100
    e0 = impl.global_error(f_harm, ex_harm, (0.0, 2.0 * math.pi), h0, "rk4")
    h1 = impl.step_for_tolerance(h0, e0, 4, 1e-8)
    e1 = impl.global_error(f_harm, ex_harm, (0.0, 2.0 * math.pi), h1, "rk4")
    assert e1 == pytest.approx(1e-8, rel=0.3)


# ---------------------------------------------------------------- stability

def test_stability_function_and_real_axis_limits():
    z = -0.3 + 0.4j
    assert impl.stability_function(z, 1) == pytest.approx(1 + z)
    assert impl.stability_function(z, 4) == pytest.approx(1 + z + z ** 2 / 2 + z ** 3 / 6 + z ** 4 / 24)
    zs = np.array([-1.0, -2.0, 0.5j])
    npt.assert_allclose(impl.stability_function(zs, 2), 1 + zs + zs ** 2 / 2)
    assert abs(impl.stability_function(-2.0, 1)) == pytest.approx(1.0)
    assert abs(impl.stability_function(-2.0, 2)) == pytest.approx(1.0)
    limits = [impl.stability_limit_real(p) for p in range(1, 5)]
    npt.assert_allclose(limits, [2.0, 2.0, 2.5127, 2.7853], atol=1e-3)
    for p in range(1, 5):
        assert abs(impl.stability_function(-limits[p - 1] * 1.01, p)) > 1.0
        assert abs(impl.stability_function(-limits[p - 1] * 0.99, p)) < 1.0
    # RK4 is stable on the imaginary axis up to |h lambda| = 2 sqrt 2, Euler nowhere
    assert abs(impl.stability_function(2.8j, 4)) < 1.0 < abs(impl.stability_function(2.85j, 4))
    assert abs(impl.stability_function(0.1j, 1)) > 1.0
    with pytest.raises(ValueError):
        impl.stability_function(z, 5)


def test_the_stability_limit_decides_blow_up_on_a_stiff_decay():
    lam = -1000.0
    f = lambda t, y: lam * y
    for m, p in [("euler", 1), ("heun", 2), ("rk4", 4)]:
        hmax = impl.max_stable_step(p, lam)
        _, Y_ok = impl.integrate(f, (0.0, 0.2), np.array([1.0]), 0.95 * hmax, m)
        _, Y_bad = impl.integrate(f, (0.0, 0.2), np.array([1.0]), 1.05 * hmax, m)
        assert abs(Y_ok[-1, 0]) < 1.0, m
        assert abs(Y_bad[-1, 0]) > 10.0, m
    with pytest.raises(ValueError):
        impl.max_stable_step(4, 1.0)


# ---------------------------------------------------------------- two-body problem

def test_two_body_rhs_and_the_kepler_reference():
    mu = impl.MU_EARTH
    r = 7.0e6
    v = math.sqrt(mu / r)
    rhs = impl.two_body_rhs(0.0, np.array([r, 0.0, 0.0, v]), mu)
    npt.assert_allclose(rhs, [0.0, v, -mu / r ** 2, 0.0], rtol=1e-15, atol=1e-12)
    a, e = 7.5e6, 0.1
    y0 = impl.kepler_state(0.0, a, e, mu)
    npt.assert_allclose(y0[:2], [a * (1 - e), 0.0], atol=1e-6)
    assert y0[3] == pytest.approx(math.sqrt(mu * (1 + e) / (a * (1 - e))), rel=1e-13)
    T = impl.orbit_period(a, mu)
    assert T == pytest.approx(2 * math.pi * math.sqrt(a ** 3 / mu), rel=1e-15)
    npt.assert_allclose(impl.kepler_state(T, a, e, mu), y0, rtol=1e-9, atol=1e-6)
    ts = T * RNG.uniform(0.0, 3.0, size=20)
    Ys = np.array([impl.kepler_state(t, a, e, mu) for t in ts])
    eps = impl.specific_energy(Ys, mu)
    npt.assert_allclose(eps, -mu / (2 * a), rtol=1e-12)
    npt.assert_allclose(impl.angular_momentum(Ys), math.sqrt(mu * a * (1 - e * e)), rtol=1e-12)
    # the Kepler state is the solution of the ODE: compare with DOP853 at tight tolerance after half an orbit
    sol = solve_ivp(lambda t, y: impl.two_body_rhs(t, y, mu), (0.0, 0.5 * T), y0, method="DOP853", rtol=1e-13, atol=1e-6)
    npt.assert_allclose(sol.y[:, -1], impl.kepler_state(0.5 * T, a, e, mu), rtol=1e-9, atol=1e-3)
    with pytest.raises(ValueError):
        impl.solve_kepler(1.0, 1.2)
    assert impl.solve_kepler(0.0, 0.5) == pytest.approx(0.0, abs=1e-15)
    E = impl.solve_kepler(2.0, 0.7)
    assert E - 0.7 * math.sin(E) == pytest.approx(2.0, abs=1e-13)


def test_rk4_orbit_error_converges_at_fourth_order_and_beats_the_lower_methods():
    a, e = 7.5e6, 0.1
    T = impl.orbit_period(a)
    spo = np.array([200, 400, 800])
    errs = {m: [impl.orbit_position_error(a, e, 1.0, int(n), m) for n in spo] for m in ["heun", "kutta3", "rk4"]}
    hs = T / spo
    assert impl.observed_order(hs, errs["rk4"]) == pytest.approx(4.0, abs=0.3)
    assert impl.observed_order(hs, errs["kutta3"]) == pytest.approx(3.0, abs=0.3)
    assert impl.observed_order(hs, errs["heun"]) == pytest.approx(2.0, abs=0.3)
    assert errs["rk4"][-1] < 1e-2 < 1.0 < errs["kutta3"][-1] < errs["heun"][-1]
    assert impl.orbit_position_error(a, e, 1.0, 400, "euler") > 1e6


def test_energy_and_angular_momentum_drift_follow_the_method():
    a, e = 7.5e6, 0.1
    t, Y = impl.propagate_orbit(a, e, 3.0, 400, "rk4")
    assert t[-1] == pytest.approx(3.0 * impl.orbit_period(a), rel=1e-14)
    de = impl.energy_drift(Y)
    assert de[0] == 0.0
    assert abs(de[-1]) < 1e-8
    dh = impl.angular_momentum(Y) / impl.angular_momentum(Y[0]) - 1.0
    assert np.max(np.abs(dh)) < 1e-8
    # RK4 on an orbit is dissipative: the energy error after three orbits is negative and about three times
    # the error after one (linear growth of a per-step bias)
    i1 = np.searchsorted(t, impl.orbit_period(a))
    assert de[-1] < 0.0
    assert de[-1] / de[i1] == pytest.approx(3.0, rel=0.25)
    # Euler spirals outward: energy increases monotonically on average
    _, Ye = impl.propagate_orbit(a, e, 1.0, 400, "euler")
    dee = impl.energy_drift(Ye)
    assert dee[-1] > 1e-3 and np.all(np.diff(dee)[::50] > 0.0)


def test_cost_for_a_tolerance_predicts_the_rk4_step():
    a, e = 7.5e6, 0.1
    T = impl.orbit_period(a)
    e400 = impl.orbit_position_error(a, e, 1.0, 400, "rk4")
    h_req = impl.step_for_tolerance(T / 400, e400, 4, 1.0)
    n_req = int(math.ceil(T / h_req))
    assert 150 < n_req < 300
    err = impl.orbit_position_error(a, e, 1.0, n_req, "rk4")
    assert 0.5 < err < 2.0
    # same evaluation budget, lower order: Heun with 2 n_req steps is worse by more than a factor 1000
    assert impl.orbit_position_error(a, e, 1.0, 2 * n_req, "heun") > 1e3 * err


# ---------------------------------------------------------------- main

def test_main_runs_and_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "rk4_orbit.png").exists()
