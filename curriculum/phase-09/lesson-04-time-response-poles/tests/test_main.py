"""Public tests for P09.04 Time Response & Poles.

The independent checks: step-response metrics measured from a sampled response against the zeta-omega_n
closed-form formulas; the partial-fraction (residue) step response against an RK4 simulation and against
``scipy.signal.step`` (oracle only). Also: limiting cases (zeta = 0, zeta = 1, overdamped, a far third
pole, a far zero), first-order metrics, invariance of the metrics under time and amplitude scaling, the
fixed real part of a proportional servo loop, and the undershoot of a right-half-plane zero.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
from __future__ import annotations

import numpy as np
import numpy.testing as npt
import pytest
from scipy import signal

from conftest import load_impl

impl = load_impl(__file__)


# ---------------------------------------------------------------- metrics against the formulas

@pytest.mark.parametrize("zeta", [0.2, 0.45, 0.7])
def test_measured_metrics_match_the_second_order_formulas(zeta):
    wn = 3.0
    t = np.linspace(0.0, 12.0, 24001)                 # dt = 0.5 ms
    m = impl.step_metrics(t, impl.second_order_step(t, zeta, wn), 1.0)
    assert m["overshoot"] == pytest.approx(impl.overshoot_formula(zeta), abs=1e-4)
    assert m["overshoot"] == pytest.approx(100 * np.exp(-np.pi * zeta / np.sqrt(1 - zeta ** 2)), abs=1e-4)
    assert m["peak_time"] == pytest.approx(np.pi / (wn * np.sqrt(1 - zeta ** 2)), abs=6e-4)
    assert m["rise_time"] == pytest.approx(impl.rise_time_exact(zeta, wn), abs=1e-6)
    # the envelope bound is an upper bound; 4/(zeta wn) is within about 25 %
    assert m["settling_time"] <= impl.settling_time_envelope(zeta, wn) + 1e-9
    assert m["settling_time"] == pytest.approx(4 / (zeta * wn), rel=0.25)
    assert m["undershoot"] < 1e-9


def test_first_order_metrics():
    tau = 0.05
    t = np.linspace(0.0, 0.5, 50001)
    m = impl.step_metrics(t, 1 - np.exp(-t / tau), 1.0)
    assert m["rise_time"] == pytest.approx(tau * np.log(9.0), rel=1e-5)
    assert m["settling_time"] == pytest.approx(tau * np.log(50.0), rel=1e-5)
    assert m["overshoot"] == 0.0
    assert m["steady_state"] == 1.0 and impl.step_metrics(t, 1 - np.exp(-t / tau))["steady_state"] == pytest.approx(1 - np.exp(-10.0))


def test_metrics_scale_with_time_and_ignore_amplitude():
    zeta = 0.4
    t = np.linspace(0.0, 20.0, 40001)
    m1 = impl.step_metrics(t, impl.second_order_step(t, zeta, 1.0), 1.0)
    m5 = impl.step_metrics(t / 5.0, impl.second_order_step(t / 5.0, zeta, 5.0), 1.0)
    for k in ("rise_time", "peak_time", "settling_time"):
        assert m5[k] == pytest.approx(m1[k] / 5.0, rel=1e-6)
    mA = impl.step_metrics(t, -3.0 * impl.second_order_step(t, zeta, 1.0), -3.0)     # a negative step
    assert mA["overshoot"] == pytest.approx(m1["overshoot"], rel=1e-9)
    assert mA["settling_time"] == pytest.approx(m1["settling_time"], rel=1e-9)


# ---------------------------------------------------------------- the closed forms

@pytest.mark.parametrize("zeta", [0.0, 0.3, 1.0, 2.5])
def test_second_order_closed_form_matches_scipy(zeta):
    wn = 2.0
    t = np.linspace(0.0, 10.0, 1001)
    _, y = signal.step(signal.lti([wn ** 2], [1.0, 2 * zeta * wn, wn ** 2]), T=t)
    npt.assert_allclose(impl.second_order_step(t, zeta, wn), y, atol=1e-9)


def test_damping_limits():
    t = np.linspace(0.0, 20.0, 20001)
    y0 = impl.second_order_step(t, 0.0, 1.0)
    npt.assert_allclose(y0, 1 - np.cos(t), atol=1e-12)                 # undamped: between 0 and 2 for ever
    for z in (1 - 1e-6, 1 + 1e-6):                                      # continuous through critical damping
        npt.assert_allclose(impl.second_order_step(t, z, 1.0), impl.second_order_step(t, 1.0, 1.0), atol=1e-5)
    assert impl.overshoot_formula(1.0) == 0.0 and impl.overshoot_formula(0.0) == pytest.approx(100.0)
    assert impl.peak_time_formula(1.2, 1.0) == np.inf
    for os_ in (1.0, 5.0, 10.0, 25.0, 60.0):
        assert impl.overshoot_formula(impl.zeta_from_overshoot(os_)) == pytest.approx(os_, rel=1e-12)
    z, wn = impl.second_order_params(np.array([0.05, 1.0, 14.0]))
    assert wn == pytest.approx(np.sqrt(280.0)) and z == pytest.approx(1 / (2 * np.sqrt(14 * 0.05)))


# ---------------------------------------------------------------- partial fractions

def test_partial_fraction_step_matches_simulation_and_scipy_on_the_pitch_loop():
    G = impl.pitch_closed_loop()
    t = np.linspace(0.0, 6.0, 1201)
    y_pf = impl.residue_step_response(G, t)
    npt.assert_allclose(y_pf, impl.simulate_step(G, t), atol=1e-8)
    _, y_sp = signal.step(signal.lti(*G), T=t)
    npt.assert_allclose(y_pf, y_sp, atol=1e-9)


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_partial_fractions_on_random_stable_systems(seed):
    rng = np.random.default_rng(seed)
    p = np.concatenate([-rng.uniform(0.3, 4.0, 2), [complex(-rng.uniform(0.2, 2), rng.uniform(0.5, 5))]])
    p = np.concatenate([p, [np.conj(p[-1])]])
    den = np.real(np.poly(p))
    num = rng.standard_normal(3)
    G = impl.tf(num, den)
    t = np.linspace(0.0, 15.0, 1501)
    _, y_sp = signal.step(signal.lti(num, den), T=t)
    npt.assert_allclose(impl.residue_step_response(G, t), y_sp, atol=1e-9 * max(1, np.abs(y_sp).max()))


def test_residues_give_the_dc_gain_and_start_from_rest():
    G = impl.pitch_closed_loop()
    p, r = impl.step_residues(G)
    i0 = int(np.argmin(np.abs(p)))
    assert r[i0] == pytest.approx(1.0, rel=1e-10)                      # residue at s = 0 is G(0)
    assert abs(np.sum(r)) < 1e-9                                        # y(0) = 0 for a strictly proper G
    with pytest.raises(ValueError):
        impl.residue_step_response(impl.second_order_tf(1.0, 2.0), np.linspace(0, 1, 5))   # double pole


# ---------------------------------------------------------------- the servo, dominance and zeros

@pytest.mark.parametrize("os_target", [5.0, 10.0, 20.0])
def test_servo_gain_gives_the_target_overshoot_and_a_fixed_real_part(os_target):
    tau = 0.05
    K = impl.servo_gain_for_overshoot(os_target, tau)
    G = impl.servo_closed_loop(K, tau)
    t = np.linspace(0.0, 1.0, 2001)
    m = impl.step_metrics(t, impl.simulate_step(G, t), 1.0)
    assert m["overshoot"] == pytest.approx(os_target, abs=0.02)
    npt.assert_allclose(np.roots(G[1]).real, -1 / (2 * tau), rtol=1e-9)    # underdamped: Re s = -1/(2 tau) whatever K
    assert G[0][-1] / G[1][-1] == pytest.approx(1.0)                        # type 1: unit DC gain


def test_a_far_third_pole_leaves_second_order_behaviour_and_a_near_one_does_not():
    zeta, wn = 0.5, 1.0
    t = np.linspace(0.0, 30.0, 6001)
    base = impl.step_metrics(t, impl.second_order_step(t, zeta, wn), 1.0)
    far = impl.with_third_pole(zeta, wn, 20 * zeta * wn)
    dom, ratio = impl.dominant_poles(far)
    assert ratio == pytest.approx(20.0) and len(dom) == 2
    m_far = impl.step_metrics(t, impl.residue_step_response(far, t), 1.0)
    assert m_far["overshoot"] == pytest.approx(base["overshoot"], rel=0.05)
    near = impl.with_third_pole(zeta, wn, 0.5 * zeta * wn)
    dom_n, _ = impl.dominant_poles(near)
    assert len(dom_n) == 1 and dom_n[0].real == pytest.approx(-0.25)
    assert impl.step_metrics(t, impl.residue_step_response(near, t), 1.0)["overshoot"] < 0.5


def test_zero_closed_form_lhp_overshoot_and_rhp_undershoot():
    zeta, wn = 0.5, 1.0
    t = np.linspace(0.0, 25.0, 5001)
    over = []
    for z in (-5.0, -2.0, -1.0, -0.5):
        y = impl.zero_step_closed_form(t, zeta, wn, z)
        npt.assert_allclose(y, impl.residue_step_response(impl.with_zero(zeta, wn, z), t), atol=1e-10)
        over.append(impl.step_metrics(t, y, 1.0)["overshoot"])
    assert np.all(np.diff(over) > 0)                                     # a zero nearer the origin: more overshoot
    y_rhp = impl.zero_step_closed_form(t, zeta, wn, 2.0)
    _, y_sp = signal.step(signal.lti(*impl.with_zero(zeta, wn, 2.0)), T=t)
    npt.assert_allclose(y_rhp, y_sp, atol=1e-9)
    assert impl.step_metrics(t, y_rhp, 1.0)["undershoot"] > 5.0          # starts the wrong way
    # far zero -> plain second order
    npt.assert_allclose(impl.zero_step_closed_form(t, zeta, wn, -1e6), impl.second_order_step(t, zeta, wn), atol=1e-5)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "time_response.png").stat().st_size > 1000
