"""Public tests for P09.03 Transfer Functions & Block Diagrams.

The independent checks: every algebraic reduction is compared pointwise with the block formula evaluated
at random complex s (G1 G2, G1 + G2, G / (1 + G H)); the reduced pitch-attitude transfer function is
compared with a time simulation of the UNREDUCED block diagram and with ``scipy.signal.lti`` / ``step``
(oracle only). Also: limiting cases of feedback (K -> 0, K -> infinity), the state-space realisation
reproduces the TF, poles and zeros, DC gain, pole-zero cancellation, and stability of the loop.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
from __future__ import annotations

import numpy as np
import numpy.testing as npt
import pytest
from scipy import signal

pytestmark = pytest.mark.filterwarnings("ignore::scipy.signal.BadCoefficients")

from conftest import load_impl

impl = load_impl(__file__)
RNG = np.random.default_rng(20260903)
S_TEST = np.array([0.3 + 1.7j, -0.8 + 0.4j, 2.5 - 3.1j, 0.05j + 0.01, -4.0 + 9.0j])


def G_at(num, den, s):
    return np.polyval(num, s) / np.polyval(den, s)


def random_tf(nz: int, npole: int):
    num = RNG.standard_normal(nz + 1)
    den = np.concatenate([[1.0], RNG.uniform(0.5, 3.0, npole)])
    return impl.tf(num, den)


# ---------------------------------------------------------------- algebra

def test_tf_normalises_and_rejects_bad_input():
    num, den = impl.tf([0.0, 4.0, 8.0], [2.0, 6.0, 4.0])
    npt.assert_allclose(den, [1.0, 3.0, 2.0])
    npt.assert_allclose(num, [2.0, 4.0])
    with pytest.raises(ValueError):
        impl.tf([1.0, 0.0, 0.0], [1.0, 1.0])            # improper
    with pytest.raises(ValueError):
        impl.tf([1.0], [0.0, 0.0])
    npt.assert_allclose(impl.evaluate((num, den), S_TEST), G_at([4.0, 8.0], [2.0, 6.0, 4.0], S_TEST))


@pytest.mark.parametrize("orders", [((0, 1), (1, 2)), ((1, 3), (0, 1)), ((2, 2), (1, 3))])
def test_series_parallel_feedback_match_the_block_formulas_pointwise(orders):
    G1, G2 = random_tf(*orders[0]), random_tf(*orders[1])
    g1, g2 = G_at(*G1, S_TEST), G_at(*G2, S_TEST)
    npt.assert_allclose(impl.evaluate(impl.series(G1, G2), S_TEST), g1 * g2, rtol=1e-10)
    npt.assert_allclose(impl.evaluate(impl.parallel(G1, G2), S_TEST), g1 + g2, rtol=1e-10)
    npt.assert_allclose(impl.evaluate(impl.feedback(G1, G2), S_TEST), g1 / (1 + g1 * g2), rtol=1e-9)
    npt.assert_allclose(impl.evaluate(impl.feedback(G1, G2, sign=+1), S_TEST), g1 / (1 - g1 * g2), rtol=1e-9)
    npt.assert_allclose(impl.evaluate(impl.feedback(G1), S_TEST), g1 / (1 + g1), rtol=1e-9)


def test_feedback_limits_small_and_large_loop_gain():
    G = impl.tf([1.0], [0.05, 1.0, 0.0])                   # servo K/(s (tau s + 1)) with K = 1
    H = impl.tf([1.0], [0.01, 1.0])
    g, h = G_at(*G, S_TEST), G_at(*H, S_TEST)
    small = impl.feedback(impl.series(impl.tf([1e-9], [1.0]), G), H)
    npt.assert_allclose(impl.evaluate(small, S_TEST), 1e-9 * g, rtol=1e-6)        # K -> 0: open loop
    big = impl.feedback(impl.series(impl.tf([1e9], [1.0]), G), H)
    npt.assert_allclose(impl.evaluate(big, S_TEST), 1.0 / h, rtol=1e-6)            # K -> inf: 1/H


def test_poles_of_a_feedback_loop_are_the_roots_of_one_plus_GH():
    G, H = random_tf(1, 3), random_tf(0, 1)
    cl = impl.feedback(G, H)
    p = impl.poles(cl)
    assert len(p) == 4
    npt.assert_allclose(np.abs(1 + G_at(*G, p) * G_at(*H, p)), 0.0, atol=1e-7)      # characteristic equation
    # series: poles are the union
    Gs = impl.series(impl.tf([1.0], [1.0, 2.0]), impl.tf([1.0], [1.0, 2.0, 5.0]))
    npt.assert_allclose(np.sort_complex(impl.poles(Gs)), np.sort_complex(np.array([-2.0, -1 + 2j, -1 - 2j])), atol=1e-10)
    assert np.all(np.diff(impl.poles(Gs).real) >= 0)
    npt.assert_allclose(impl.zeros(impl.tf([2.0, 6.0], [1.0, 1.0, 1.0])), [-3.0])
    assert impl.zeros(impl.tf([5.0], [1.0, 1.0])).size == 0


def test_dc_gain_including_integrators():
    assert impl.dc_gain(impl.tf([3.0], [2.0, 1.0])) == pytest.approx(3.0)
    assert impl.dc_gain(impl.tf([1.0], [1.0, 0.0])) == np.inf
    assert impl.dc_gain(impl.parallel(impl.tf([3.0], [1.0, 1.0]), impl.tf([2.0], [1.0, 4.0]))) == pytest.approx(3.5)
    # a type-1 loop with unity feedback has unit DC gain whatever the gain
    for K in (0.1, 1.0, 50.0):
        assert impl.dc_gain(impl.feedback(impl.tf([K], [1.0, 1.0, 0.0]))) == pytest.approx(1.0)


def test_minreal_cancels_coincident_pole_zero_pairs_only():
    G = impl.series(impl.tf([1.0, 2.0], [1.0, 1.0]), impl.tf([3.0], [1.0, 2.0]))     # (s + 2)/(s + 1) * 3/(s + 2)
    assert len(G[1]) == 3                                                           # algebra does not cancel
    R = impl.minreal(G)
    npt.assert_allclose(R[1], [1.0, 1.0], atol=1e-9)
    npt.assert_allclose(R[0], [3.0], atol=1e-9)
    npt.assert_allclose(impl.evaluate(R, S_TEST), impl.evaluate(G, S_TEST), rtol=1e-8)
    keep = impl.tf([1.0, 2.1], [1.0, 3.0, 2.0])
    assert len(impl.minreal(keep)[1]) == 3                                          # -2.1 is not -2


# ---------------------------------------------------------------- realisation and simulation

@pytest.mark.parametrize("G", [([1.0], [1.0, 1.0]), ([2.0, 1.0], [1.0, 3.0, 5.0]),
                               ([1.0, 0.5, 4.0], [1.0, 2.0, 3.0]), ([240.0, 360.0], [1.0, 23.2, 80.0, 320.0])])
def test_state_space_realisation_reproduces_the_tf(G):
    G = impl.tf(*G)
    A, B, C, D = impl.tf_to_ss(G)
    n = len(G[1]) - 1
    assert A.shape == (n, n) and B.shape == (n, 1) and C.shape == (1, n) and D.shape == (1, 1)
    for s in S_TEST:
        val = (C @ np.linalg.solve(s * np.eye(n) - A, B) + D)[0, 0]
        assert val == pytest.approx(impl.evaluate(G, s), rel=1e-9)
    npt.assert_allclose(np.sort_complex(np.linalg.eigvals(A)), np.sort_complex(impl.poles(G)), atol=1e-8)
    oracle = signal.lti(A, B, C, D).to_tf()                                         # oracle
    num_sp, den_sp = np.atleast_1d(np.squeeze(oracle.num)), oracle.den
    npt.assert_allclose(np.polyval(num_sp, S_TEST) / np.polyval(den_sp, S_TEST), impl.evaluate(G, S_TEST), rtol=1e-8)


@pytest.mark.parametrize("G", [([4.0], [1.0, 2.0, 4.0]), ([1.0, 2.0], [1.0, 3.0, 2.0]), ([2.0, 1.0], [1.0, 1.0])])
def test_simulated_step_matches_scipy(G):
    G = impl.tf(*G)
    t = np.linspace(0.0, 8.0, 801)
    y = impl.simulate_tf(G, lambda tt: 1.0, t)
    _, y_sp = signal.step(signal.lti(*G), T=t)
    npt.assert_allclose(y, y_sp, atol=1e-8)


# ---------------------------------------------------------------- the pitch-attitude loop

def hand_reduction(p):
    """The pitch loop reduced by hand (written out independently of the impl's reduce_pitch_loop)."""
    s = S_TEST
    servo = 1 / (p["tau_s"] * s + 1)
    sp = p["M_delta"] * (s + p["inv_T_theta2"]) / (s ** 2 + 2 * p["zeta_sp"] * p["w_sp"] * s + p["w_sp"] ** 2)
    gyro = 1 / (p["tau_g"] * s + 1)
    inner = servo * sp / (1 + p["K_q"] * gyro * servo * sp)
    L = p["K_theta"] * inner / s
    return L / (1 + L)


def test_reduced_pitch_loop_equals_the_hand_reduction():
    p = impl.pitch_model()
    red = impl.reduce_pitch_loop(p)
    npt.assert_allclose(impl.evaluate(red["closed"], S_TEST), hand_reduction(p), rtol=1e-9)
    assert len(red["closed"][1]) - 1 == 5                                           # servo + 2 + gyro + integrator
    assert impl.dc_gain(red["closed"]) == pytest.approx(1.0, rel=1e-12)            # integrator in the forward path
    z = impl.zeros(red["closed"])
    npt.assert_allclose(np.sort(z.real), [-1.0 / p["tau_g"], -p["inv_T_theta2"]], rtol=1e-9)   # feedback-path pole -> zero


def test_unreduced_diagram_simulation_matches_the_reduced_tf_and_scipy():
    p = impl.pitch_model()
    t = np.linspace(0.0, 6.0, 601)
    hist = impl.simulate_diagram(impl.pitch_blocks(p), impl.pitch_connect(p), lambda tt: 1.0, t)
    red = impl.reduce_pitch_loop(p)["closed"]
    _, y_sp = signal.step(signal.lti(*red), T=t)
    npt.assert_allclose(hist["integrator"], y_sp, atol=1e-7)
    npt.assert_allclose(impl.simulate_tf(red, lambda tt: 1.0, t), hist["integrator"], atol=1e-9)
    # the measured rate converges to zero and the attitude to the command
    assert abs(hist["gyro"][-1]) < 5e-3 and hist["integrator"][-1] == pytest.approx(1.0, abs=5e-3)


def test_pitch_damper_stabilises_the_attitude_loop():
    p = impl.pitch_model()
    assert np.all(impl.poles(impl.reduce_pitch_loop(p)["closed"]).real < 0)
    bad = dict(p, K_q=0.0, K_theta=4.0)
    assert impl.poles(impl.reduce_pitch_loop(bad)["closed"]).real.max() > 0
    t = np.linspace(0.0, 60.0, 3001)
    hist = impl.simulate_diagram(impl.pitch_blocks(bad), impl.pitch_connect(bad), lambda tt: 1.0, t)
    assert np.max(np.abs(hist["integrator"][-500:] - 1)) > 2 * np.max(np.abs(hist["integrator"][1000:1500] - 1))


def test_servo_position_loop_diagram_equals_its_reduction():
    K, tau = 20.0, 0.05                                     # second order: wn = 20 rad/s, zeta = 0.5
    blocks = {"motor": impl.tf([K], [tau, 1.0, 0.0])}
    connect = lambda y, r: {"motor": r - y["motor"]}
    t = np.linspace(0.0, 1.0, 501)
    hist = impl.simulate_diagram(blocks, connect, lambda tt: 1.0, t)
    cl = impl.feedback(blocks["motor"])
    npt.assert_allclose(cl[1], [1.0, 20.0, 400.0])
    _, y_sp = signal.step(signal.lti(*cl), T=t)
    npt.assert_allclose(hist["motor"], y_sp, atol=1e-7)
    with pytest.raises(ValueError):
        impl.simulate_diagram({"pd": impl.tf([1.0, 1.0], [0.1, 1.0])}, lambda y, r: {"pd": r}, lambda tt: 1.0, t)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "block_diagram.png").stat().st_size > 1000
