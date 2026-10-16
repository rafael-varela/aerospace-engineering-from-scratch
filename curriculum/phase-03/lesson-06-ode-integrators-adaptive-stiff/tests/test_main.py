"""Public tests for P03.06 ODE Integrators II: Adaptive & Stiff. They assert properties, with
``scipy.integrate.solve_ivp`` and ``scipy.linalg.expm`` as the oracles.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-03/parameterized-cases/.
"""
import math

import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import solve_ivp

from conftest import load_impl

impl = load_impl(__file__)


def oscillator(t, y):
    return np.array([y[1], -y[0]])


def oscillator_exact(t):
    return np.array([np.cos(t), -np.sin(t)])


def decay(t, y):
    return -y


# ---------------------------------------------------------------- one step of each method

def test_rk45_step_is_fifth_order_accurate_on_the_exponential():
    f = lambda t, y: y
    y0 = np.array([1.0])
    for h in (0.2, 0.1):
        y5, err = impl.rk45_step(f, 0.0, y0, h)
        # Fifth order: one-step error O(h^6); the embedded estimate is the fourth-order error O(h^5).
        assert abs(y5[0] - np.exp(h)) < 2.0 * h ** 6
        assert abs(err[0]) < 2.0 * h ** 5
        assert abs(err[0]) > abs(y5[0] - np.exp(h))        # the estimate is conservative


def test_embedded_error_estimate_scales_as_h_to_the_fifth():
    y0 = np.array([1.0, 0.0])
    est = [np.linalg.norm(impl.rk45_step(oscillator, 0.0, y0, h)[1]) for h in (0.2, 0.1, 0.05)]
    ratios = np.array(est[:-1]) / np.array(est[1:])
    npt.assert_allclose(ratios, 32.0, rtol=0.15)          # 2^5 = 32 per halving


def test_implicit_euler_step_is_exact_backward_euler_on_a_linear_problem():
    lam, h = -3.0, 0.7
    f = lambda t, y: lam * y
    jac = lambda t, y: np.array([[lam]])
    z, n_it = impl.implicit_euler_step(f, 0.0, np.array([2.0]), h, jac=jac)
    assert z[0] == pytest.approx(2.0 / (1.0 - h * lam), rel=1e-12)     # y_{n+1} = y_n / (1 - h lambda)
    assert n_it <= 2                                                   # Newton is exact in one step on a linear problem
    # Without an analytic Jacobian the finite-difference one gives the same answer.
    z_fd, _ = impl.implicit_euler_step(f, 0.0, np.array([2.0]), h)
    assert z_fd[0] == pytest.approx(z[0], rel=1e-8)


def test_implicit_euler_raises_when_newton_cannot_converge():
    f = lambda t, y: y ** 2
    with pytest.raises(RuntimeError):
        impl.implicit_euler_step(f, 0.0, np.array([1.0]), 0.5, max_iter=1, tol=1e-14)


def test_fd_jacobian_matches_the_analytic_thermal_jacobian():
    p = impl.SkinStructure()
    y = np.array([350.0, 310.0])
    J_fd = impl.fd_jacobian(lambda t, y: impl.thermal_rhs(t, y, p), 0.0, y)
    J = impl.thermal_jac(0.0, y, p)
    npt.assert_allclose(J_fd, J, rtol=1e-5, atol=1e-9)
    assert J.shape == (2, 2) and J[0, 0] < 0 and J[1, 1] < 0


# ---------------------------------------------------------------- observed convergence orders (fixed step)

def test_fixed_step_rk45_is_fifth_order_and_rk4_fourth_on_the_oscillator():
    y0 = np.array([1.0, 0.0])
    ns = [20, 40, 80, 160]
    hs = 2 * np.pi / np.array(ns)
    e45 = impl.fixed_step_errors(impl.rk45_step, oscillator, oscillator_exact, (0.0, 2 * np.pi), y0, ns)
    e4 = impl.fixed_step_errors(impl.rk4_step, oscillator, oscillator_exact, (0.0, 2 * np.pi), y0, ns)
    assert 4.7 < impl.convergence_order(hs, e45) < 5.3, e45
    assert 3.8 < impl.convergence_order(hs, e4) < 4.2, e4
    assert np.all(e45 < e4)


def test_fixed_step_implicit_euler_is_first_order():
    ns = [50, 100, 200]
    hs = 2.0 / np.array(ns)
    errs = impl.fixed_step_errors(impl.implicit_euler_step, decay, lambda t: np.array([np.exp(-t)]), (0.0, 2.0),
                                  np.array([1.0]), ns, jac=lambda t, y: np.array([[-1.0]]))
    assert 0.9 < impl.convergence_order(hs, errs) < 1.1, errs
    # Backward Euler lags a decaying exponential: 1/(1 + h) > e^{-h}, so it stays above the truth.
    res = impl.solve_fixed(impl.implicit_euler_step, decay, (0.0, 2.0), np.array([1.0]), 50)
    assert np.all(res.y[0, 1:] > np.exp(-res.t[1:]))


def test_convergence_order_helper():
    hs = np.array([0.4, 0.2, 0.1])
    assert impl.convergence_order(hs, 3.0 * hs ** 2) == pytest.approx(2.0, abs=1e-10)
    with pytest.raises(ValueError):
        impl.convergence_order(hs, np.array([1.0, 0.0, -1.0]))
    with pytest.raises(ValueError):
        impl.solve_fixed(impl.rk4_step, decay, (0.0, 1.0), np.array([1.0]), 0)


# ---------------------------------------------------------------- linear stability

def test_stability_functions_match_the_hand_formulas():
    z = complex(-0.7, 0.4)
    assert impl.stability_function("explicit_euler", z) == pytest.approx(1 + z)
    assert impl.stability_function("rk4", z) == pytest.approx(1 + z + z ** 2 / 2 + z ** 3 / 6 + z ** 4 / 24)
    assert impl.stability_function("implicit_euler", z) == pytest.approx(1 / (1 - z))
    # The fifth-order DP weights reproduce e^z through the z^5 term.
    r = impl.stability_function("rk45", z)
    taylor5 = sum(z ** k / math.factorial(k) for k in range(6))
    assert abs(r - taylor5) < abs(z) ** 6 / 100 and abs(r - np.exp(z)) < 1e-3   # agree through z^5
    with pytest.raises(ValueError):
        impl.stability_function("leapfrog", z)


def test_real_stability_intervals():
    assert impl.stability_interval("explicit_euler") == pytest.approx(2.0, abs=1e-6)
    assert impl.stability_interval("rk4") == pytest.approx(2.7853, abs=2e-3)
    assert 3.2 < impl.stability_interval("rk45") < 3.4
    assert impl.stability_interval("implicit_euler") == np.inf


def test_stiffness_ratio_and_explicit_euler_blow_up():
    assert impl.stiffness_ratio(np.diag([-1.0, -1000.0])) == pytest.approx(1000.0)
    assert impl.stiffness_ratio(np.array([[0.0, 1.0], [-1.0, 0.0]])) == np.inf
    p = impl.SkinStructure()
    T0 = np.array([300.0, 300.0])
    J = impl.thermal_jac(0.0, T0, p)
    assert 3e4 < impl.stiffness_ratio(J) < 5e4
    lam_fast = np.max(np.abs(np.linalg.eigvals(J).real))
    f = lambda t, y: impl.thermal_rhs(t, y, p)
    h_ok, h_bad = 0.8 * 2 / lam_fast, 1.2 * 2 / lam_fast
    stable = impl.solve_fixed(impl.explicit_euler_step, f, (0.0, 200 * h_ok), T0, 200)
    assert np.all(np.isfinite(stable.y)) and np.max(np.abs(stable.y - 300.0)) < 50.0
    with np.errstate(over="ignore", invalid="ignore"):
        unstable = impl.solve_fixed(impl.explicit_euler_step, f, (0.0, 200 * h_bad), T0, 200)
    assert (not np.all(np.isfinite(unstable.y))) or np.max(np.abs(unstable.y)) > 1e4


# ---------------------------------------------------------------- the adaptive driver

def test_solve_rk45_reaches_t1_exactly_and_bookkeeps_its_steps():
    res = impl.solve_rk45(oscillator, (0.0, 10.0), np.array([1.0, 0.0]), rtol=1e-8, atol=1e-10)
    assert res.t[0] == 0.0 and res.t[-1] == pytest.approx(10.0, abs=1e-12)
    assert np.all(np.diff(res.t) > 0)
    assert res.y.shape == (2, res.t.size)
    assert res.h.size == res.n_accepted == res.t.size - 1
    npt.assert_allclose(res.h, np.diff(res.t), rtol=1e-12)
    assert res.nfev >= 7 * (res.n_accepted + res.n_rejected)
    npt.assert_allclose(res.y[:, -1], oscillator_exact(10.0), atol=1e-6)


def test_solve_rk45_error_falls_and_work_rises_with_the_tolerance():
    y0 = np.array([1.0])
    exact = np.exp(-5.0)
    errs, steps = [], []
    for rtol in (1e-4, 1e-6, 1e-8, 1e-10):
        r = impl.solve_rk45(decay, (0.0, 5.0), y0, rtol=rtol, atol=1e-14)
        errs.append(abs(r.y[0, -1] - exact)); steps.append(r.n_accepted)
    errs, steps = np.array(errs), np.array(steps)
    assert np.all(np.diff(errs) < 0) and np.all(np.diff(steps) > 0)
    assert np.all(errs < 100 * np.array([1e-4, 1e-6, 1e-8, 1e-10]))
    # Local error O(h^5) at a fixed tolerance means the step count grows like rtol^(-1/5):
    # a factor 100 in rtol should cost about 100^(1/5) = 2.5 in steps.
    growth = np.log(steps[1:] / steps[:-1]) / np.log(100.0)
    assert np.all((0.12 < growth) & (growth < 0.30)), (steps, growth)


def test_solve_rk45_rejects_steps_and_stops_on_request():
    cap = impl.Capsule()
    y0 = impl.entry_state()
    res = impl.solve_rk45(lambda t, y: impl.reentry_rhs(t, y, cap), (0.0, 1500.0), y0, rtol=1e-6, atol=1e-3,
                          stop=lambda t, y: y[0] < 20e3)
    assert res.n_rejected > 0                                   # the peak forces at least one rejection
    assert res.y[0, -1] < 20e3 and res.y[0, -2] >= 20e3        # stopped right after crossing 20 km
    assert res.t[-1] < 1500.0
    with pytest.raises(RuntimeError):
        impl.solve_rk45(oscillator, (0.0, 100.0), np.array([1.0, 0.0]), rtol=1e-10, atol=1e-12, max_steps=10)


def test_rk45_matches_solve_ivp_on_the_re_entry_and_concentrates_steps_at_the_peak():
    cap = impl.Capsule()
    y0 = impl.entry_state()
    f = lambda t, y: impl.reentry_rhs(t, y, cap)
    res = impl.solve_rk45(f, (0.0, 1500.0), y0, rtol=1e-8, atol=1e-3, stop=lambda t, y: y[0] < 20e3)
    ref = solve_ivp(f, (0.0, res.t[-1]), y0, method="DOP853", rtol=1e-12, atol=1e-9, t_eval=res.t)
    assert np.max(np.abs(res.y[1] - ref.y[1])) < 1e-2             # m/s, over a 7.8 km/s entry
    assert np.max(np.abs(res.y[0] - ref.y[0])) < 1.0              # m
    a = impl.deceleration(res, cap)
    i = int(np.argmax(a))
    assert 17.5 < a[i] / cap.g0 < 18.5                            # with gravity: a little above Allen-Eggers
    assert 40e3 < res.y[0, i] < 44e3
    # After the start-up steps, the smallest step is taken near the deceleration peak and the
    # controller has stretched the step by more than a factor 3 elsewhere.
    h_run, t_run = res.h[4:], res.t[5:]
    assert abs(t_run[np.argmin(h_run)] - res.t[i]) < 40.0, "smallest step should sit near the peak"
    assert np.max(h_run) > 3.0 * np.min(h_run)


def test_drag_only_entry_reproduces_allen_eggers():
    cap = impl.Capsule()
    y0 = impl.entry_state()
    res = impl.solve_rk45(lambda t, y: impl.reentry_rhs(t, y, cap, gravity=False), (0.0, 1500.0), y0,
                          rtol=1e-9, atol=1e-3, stop=lambda t, y: y[0] < 20e3)
    v_ae = impl.allen_eggers_velocity(res.y[0], y0[0], y0[1], y0[2], cap)
    assert np.max(np.abs(res.y[1] - v_ae)) < 1e-2
    a_max, h_star = impl.allen_eggers_peak(y0[1], y0[2], cap)
    a = impl.deceleration(res, cap)
    assert np.max(a) == pytest.approx(a_max, rel=5e-3)
    assert abs(res.y[0, np.argmax(a)] - h_star) < 1.5e3          # within one coarse step of the predicted altitude
    # a_max is independent of the ballistic coefficient; h* is not.
    heavy = impl.Capsule(mass=2 * cap.mass)
    a2, h2 = impl.allen_eggers_peak(y0[1], y0[2], heavy)
    assert a2 == pytest.approx(a_max) and h2 < h_star
    assert heavy.beta == pytest.approx(2 * cap.beta)


# ---------------------------------------------------------------- the stiff thermal problem

def test_implicit_euler_matches_the_matrix_exponential_on_the_linear_model():
    p = impl.SkinStructure(eps=0.0)
    T0 = np.array([300.0, 300.0])
    A, c = impl.thermal_linear_system(p)
    npt.assert_allclose(A @ T0 + c, impl.thermal_rhs(0.0, T0, p), rtol=1e-12, atol=1e-12)
    npt.assert_allclose(A, impl.thermal_jac(0.0, T0, p), rtol=1e-12)
    exact = impl.thermal_exact_linear(3600.0, T0, p)[:, 0]
    res = impl.solve_implicit_euler(lambda t, y: impl.thermal_rhs(t, y, p), (0.0, 3600.0), T0, 600,
                                    jac=lambda t, y: impl.thermal_jac(t, y, p))
    assert np.max(np.abs(res.y[:, -1] - exact)) < 0.02         # K, with h = 6 s = 30 x the fast time constant
    # Equilibrium: skin warmer than structure by q/k_c, structure above the sink by q/k_b.
    y_eq = -np.linalg.solve(A, c)
    assert y_eq[1] == pytest.approx(p.T_inf + p.q_in / p.k_b)
    assert y_eq[0] - y_eq[1] == pytest.approx(p.q_in / p.k_c)


def test_implicit_euler_is_first_order_and_accurate_on_the_stiff_nonlinear_model():
    p = impl.SkinStructure()
    T0 = np.array([300.0, 300.0])
    f = lambda t, y: impl.thermal_rhs(t, y, p)
    jac = lambda t, y: impl.thermal_jac(t, y, p)
    ref = solve_ivp(f, (0.0, 3600.0), T0, method="Radau", jac=jac, rtol=1e-12, atol=1e-10).y[:, -1]
    errs = []
    for n in (60, 120, 240):
        res = impl.solve_implicit_euler(f, (0.0, 3600.0), T0, n, jac=jac)
        errs.append(np.max(np.abs(res.y[:, -1] - ref)))
        assert res.n_newton >= n and res.nfev >= res.n_newton
    assert errs[0] < 0.2                                        # K after an hour with 60 s steps
    assert 0.9 < impl.convergence_order(3600.0 / np.array([60, 120, 240]), errs) < 1.1, errs


def test_adaptive_rk45_is_stability_limited_on_the_stiff_problem():
    p = impl.SkinStructure()
    T0 = np.array([300.0, 300.0])
    f = lambda t, y: impl.thermal_rhs(t, y, p)
    r45 = impl.solve_rk45(f, (0.0, 1800.0), T0, rtol=1e-6, atol=1e-6)
    rie = impl.solve_implicit_euler(f, (0.0, 1800.0), T0, 30, jac=lambda t, y: impl.thermal_jac(t, y, p))
    lam_fast = np.max(np.abs(np.linalg.eigvals(impl.thermal_jac(0.0, T0, p)).real))
    assert np.max(r45.h) < 1.2 * impl.stability_interval("rk45") / lam_fast   # steps capped by stability, not accuracy
    assert r45.n_accepted > 50 * rie.n_accepted
    assert r45.nfev > 50 * rie.nfev
    assert abs(r45.y[1, -1] - rie.y[1, -1]) < 0.1                # both are right; one is 100 times dearer


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "adaptive_stiff.png").stat().st_size > 1000
