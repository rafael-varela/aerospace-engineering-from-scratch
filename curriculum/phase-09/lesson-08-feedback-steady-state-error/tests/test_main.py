"""Public tests for P09.08 Feedback & Steady-State Error. They assert properties, not magic numbers: the
closed loop, the sensitivity and the disturbance path equal L/(1+L), 1/(1+L) and P/(1+CP) pointwise and
S + T = 1; the system type and the error constants follow from the integrators; the steady-state
errors to step, ramp and parabola follow from the constants; the final-value theorem agrees with a
long from-scratch simulation for steps, ramps and disturbances (the independent check) and refuses to
answer when the closed loop is unstable, where the simulation diverges; the critical proportional gain
agrees with the Routh condition for a cubic; and the limits K -> 0 and an integrator behave as derived.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy import signal

from conftest import load_impl

impl = load_impl(__file__)

S_PTS = np.array([0.3j, 1.0 + 2.0j, -0.2 + 0.7j, 5.0j])
NUM_P, DEN_P = impl.cruise_plant()
C_AERO = 2 * impl.D0 / impl.V0


# ---------------------------------------------------------------- loop algebra

def test_given_simulator_agrees_with_scipy_lsim():
    t = np.arange(0.0, 60.0, 0.05)
    u = np.sin(0.3 * t) + (t > 5)
    _, ys, _ = signal.lsim((NUM_P, DEN_P), u, t, interp=False)
    npt.assert_allclose(impl.simulate_tf(NUM_P, DEN_P, t, u), ys, atol=1e-10)


def test_closed_loop_sensitivity_and_disturbance_paths_pointwise():
    nc, dc = impl.pi_controller(400.0, 40.0)
    nl, dl = impl.open_loop(nc, dc, NUM_P, DEN_P)
    L = impl.tf_eval(nc, dc, S_PTS) * impl.tf_eval(NUM_P, DEN_P, S_PTS)
    npt.assert_allclose(impl.tf_eval(nl, dl, S_PTS), L, rtol=1e-12)
    npt.assert_allclose(impl.tf_eval(*impl.closed_loop(nl, dl), S_PTS), L / (1 + L), rtol=1e-12)
    npt.assert_allclose(impl.tf_eval(*impl.sensitivity(nl, dl), S_PTS), 1 / (1 + L), rtol=1e-12)
    P = impl.tf_eval(NUM_P, DEN_P, S_PTS)
    npt.assert_allclose(impl.tf_eval(*impl.disturbance_to_output(nc, dc, NUM_P, DEN_P), S_PTS), P / (1 + L), rtol=1e-12)


def test_sensitivity_plus_complementary_sensitivity_is_one():
    nl, dl = impl.open_loop(*impl.p_controller(1000.0), NUM_P, DEN_P)
    w = 1j * np.logspace(-3, 2, 200)
    S = impl.tf_eval(*impl.sensitivity(nl, dl), w)
    T = impl.tf_eval(*impl.closed_loop(nl, dl), w)
    npt.assert_allclose(S + T, 1.0, atol=1e-12)
    assert abs(S[0]) == pytest.approx(1 / (1 + 1000.0 / C_AERO), rel=1e-3)   # |S| at low frequency = step error
    assert abs(S[-1]) == pytest.approx(1.0, abs=1e-4)                         # no feedback action at high frequency


# ---------------------------------------------------------------- type and error constants

@pytest.mark.parametrize("ctrl,typ", [((np.array([400.0]), np.array([1.0])), 0),
                                      ((np.array([400.0, 40.0]), np.array([1.0, 0.0])), 1),
                                      ((np.array([1.0, 1.0]), np.array([1.0, 0.0, 0.0])), 2)])
def test_system_type_counts_the_integrators(ctrl, typ):
    nl, dl = impl.open_loop(*ctrl, NUM_P, DEN_P)
    assert impl.system_type(nl, dl) == typ


def test_error_constants_of_types_0_1_2():
    K0 = impl.error_constants(*impl.open_loop(*impl.p_controller(400.0), NUM_P, DEN_P))
    assert K0["Kp"] == pytest.approx(400.0 / C_AERO) and K0["Kv"] == 0.0 and K0["Ka"] == 0.0
    K1 = impl.error_constants(*impl.open_loop(*impl.pi_controller(400.0, 40.0), NUM_P, DEN_P))
    assert np.isinf(K1["Kp"]) and K1["Kv"] == pytest.approx(40.0 / C_AERO) and K1["Ka"] == 0.0
    # a satellite attitude loop: PD on 1/(J s^2) is type 2, Ka = K / J
    K2 = impl.error_constants(np.array([2.0, 0.5]), np.array([10.0, 0.0, 0.0]))
    assert np.isinf(K2["Kp"]) and np.isinf(K2["Kv"]) and K2["Ka"] == pytest.approx(0.05)


def test_steady_state_errors_follow_from_the_constants():
    nl, dl = impl.open_loop(*impl.p_controller(400.0), NUM_P, DEN_P)
    assert impl.steady_state_error(nl, dl, "step", 2.0) == pytest.approx(2.0 / (1 + 400.0 / C_AERO))
    assert np.isinf(impl.steady_state_error(nl, dl, "ramp"))
    nl, dl = impl.open_loop(*impl.pi_controller(400.0, 40.0), NUM_P, DEN_P)
    assert impl.steady_state_error(nl, dl, "step") == 0.0
    assert impl.steady_state_error(nl, dl, "ramp", 0.5) == pytest.approx(0.5 * C_AERO / 40.0)
    assert np.isinf(impl.steady_state_error(nl, dl, "parabola"))
    with pytest.raises(ValueError):                                   # unstable loop: no steady state
        impl.steady_state_error(*impl.open_loop(*impl.p_controller(1e4), NUM_P, DEN_P), "step")


def test_gain_limits_k_to_zero_and_large_k():
    errs = [impl.steady_state_error(*impl.open_loop(*impl.p_controller(k), NUM_P, DEN_P)) for k in (1e-6, 10.0, 100.0, 1000.0, 4000.0)]
    assert errs[0] == pytest.approx(1.0, abs=1e-6)                     # no feedback: the whole step is error
    assert all(e1 > e2 for e1, e2 in zip(errs, errs[1:]))              # more gain, less error ...
    assert errs[-1] > 0.0                                              # ... but never zero for type 0


# ---------------------------------------------------------------- final-value theorem: the independent check

@pytest.mark.parametrize("case", [("P", "step", 1.0), ("P", "disturbance", 200.0), ("PI", "step", 1.0),
                                  ("PI", "disturbance", 200.0), ("PI", "ramp", 0.5)])
def test_final_value_theorem_agrees_with_a_long_simulation(case):
    ctrl, kind, A = case
    c = impl.p_controller(400.0) if ctrl == "P" else impl.pi_controller(400.0, 40.0)
    pred, sim = impl.fvt_vs_simulation(*c, kind, t_end=600.0, dt=0.05, amplitude=A)
    # a held (zero-order-hold) ramp lags the true ramp by dt/2 on average, so the simulated ramp error
    # exceeds the theorem's by A dt / 2; steps are exact
    extra = A * 0.05 / 2 if kind == "ramp" else 0.0
    assert sim == pytest.approx(pred + extra, abs=1e-6)


def test_final_value_formula_matches_the_error_constants():
    nl, dl = impl.open_loop(*impl.p_controller(400.0), NUM_P, DEN_P)
    ns, ds = impl.sensitivity(nl, dl)
    assert impl.final_value(ns, np.convolve(ds, [1.0, 0.0])) == pytest.approx(impl.steady_state_error(nl, dl, "step"))
    # disturbance: P control leaves d / (c + Kp); open loop would lose d / c
    nd, dd = impl.disturbance_to_output(*impl.p_controller(400.0), NUM_P, DEN_P)
    assert impl.final_value(200.0 * nd, np.convolve(dd, [1.0, 0.0])) == pytest.approx(200.0 / (C_AERO + 400.0))


def test_final_value_theorem_refuses_an_unstable_loop_whose_simulation_diverges():
    nl, dl = impl.open_loop(*impl.p_controller(6000.0), NUM_P, DEN_P)
    ns, ds = impl.sensitivity(nl, dl)
    with pytest.raises(ValueError):
        impl.final_value(ns, np.convolve(ds, [1.0, 0.0]))
    t = np.arange(0.0, 60.0, 0.05)
    y = impl.simulate_tf(*impl.closed_loop(nl, dl), t, np.ones_like(t))
    assert abs(1 - y[-1]) > 1.0                                        # the formula's 0.65 % would be a lie
    with pytest.raises(ValueError):                                    # an undamped oscillator has no final value
        impl.final_value(np.array([1.0]), np.array([1.0, 0.0, 4.0]))


def test_critical_gain_agrees_with_routh_for_the_cubic():
    a3, a2, a1, a0 = DEN_P
    # Routh for a3 s^3 + a2 s^2 + a1 s + (a0 + K): stable while a2 a1 > a3 (a0 + K)
    assert impl.critical_gain(NUM_P, DEN_P) == pytest.approx(a2 * a1 / a3 - a0, rel=1e-8)
    K = impl.critical_gain(NUM_P, DEN_P)
    assert impl.is_stable(impl.poly_add(DEN_P, 0.99 * K * NUM_P))
    assert not impl.is_stable(impl.poly_add(DEN_P, 1.01 * K * NUM_P))
    assert np.isinf(impl.critical_gain(np.array([1.0]), np.array([1.0, 3.0, 2.0])))   # second order: never unstable


def test_zero_input_gives_zero_error_and_responses_superpose():
    nl, dl = impl.open_loop(*impl.pi_controller(400.0, 40.0), NUM_P, DEN_P)
    nt, dt_ = impl.closed_loop(nl, dl)
    t = np.arange(0.0, 30.0, 0.05)
    npt.assert_allclose(impl.simulate_tf(nt, dt_, t, np.zeros_like(t)), 0.0, atol=0.0)
    r1, r2 = np.ones_like(t), 0.5 * t
    npt.assert_allclose(impl.simulate_tf(nt, dt_, t, r1 + r2),
                        impl.simulate_tf(nt, dt_, t, r1) + impl.simulate_tf(nt, dt_, t, r2), atol=1e-10)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "cruise_control.png").exists()
