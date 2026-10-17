"""Public tests for P09.15 Pole Placement & Observers.

They assert design properties, not magic numbers: eig(A - BK) equals the target poles and K equals the
unique single-input gain from ``scipy.signal.place_poles``; the observer gain by duality gives
eig(A - LC) equal to its targets and matches ``place_poles`` on (A^T, C^T); the closed loop with the
observer has exactly the union of the two pole sets (separation); the observer error follows
e^{(A - LC) t} e0 whatever the reference and the regulator, and decays at the slowest placed observer
pole (fitted decay rate, the lesson's independent check); the reference gain removes the steady error;
uncontrollable and unobservable cases raise; and SYNTHETIC tracker noise hurts more with a faster observer.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy import signal
from scipy.linalg import expm

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(91501)
SAT = impl.satellite_axis(3.0)
PITCH = impl.pitch_model()


def random_controllable(n):
    A = RNG.standard_normal((n, n))
    B = RNG.standard_normal((n, 1))
    return A, B


CASES = [
    (SAT[0], SAT[1], impl.second_order_poles(0.707, 0.1)),
    (SAT[0], SAT[1], np.array([-0.2, -0.5])),
    (PITCH[0], PITCH[1], np.array([-4 + 4j, -4 - 4j, -6.0, -20.0])),
    (*random_controllable(3), np.array([-1.0, -2.0 + 1.0j, -2.0 - 1.0j])),
    (*random_controllable(5), np.array([-1.0, -1.5, -2.0, -3.0, -4.0])),
]


def sorted_c(x):
    return np.sort_complex(np.asarray(x, dtype=complex))


# ---------------------------------------------------------------- regulator

@pytest.mark.parametrize("A,B,poles", CASES)
def test_ackermann_places_every_pole(A, B, poles):
    K = impl.ackermann(A, B, poles)
    assert K.shape == (1, A.shape[0])
    npt.assert_allclose(sorted_c(np.linalg.eigvals(A - B @ K)), sorted_c(poles), atol=1e-6 * max(1.0, np.abs(poles).max()))


@pytest.mark.parametrize("A,B,poles", CASES[:4])
def test_ackermann_matches_scipy_place_poles(A, B, poles):
    K = impl.ackermann(A, B, poles)
    ref = signal.place_poles(A, B, poles).gain_matrix          # single input: the gain is unique
    npt.assert_allclose(K, ref, rtol=1e-7, atol=1e-9)


def test_satellite_gain_in_closed_form_and_open_loop_poles_need_no_gain():
    J, wn, zeta = 3.0, 0.1, 0.707
    A, B, C = SAT
    K = impl.ackermann(A, B, impl.second_order_poles(zeta, wn))
    npt.assert_allclose(K, [[J * wn ** 2, J * 2 * zeta * wn]], rtol=1e-10)
    # targets equal to the open-loop eigenvalues: K = 0
    Ap = PITCH[0]
    npt.assert_allclose(impl.ackermann(Ap, PITCH[1], np.linalg.eigvals(Ap)), 0.0, atol=1e-9)


def test_gain_transforms_with_the_coordinates():
    A, B = random_controllable(3)
    poles = np.array([-1.0, -2.0, -3.0])
    T = RNG.standard_normal((3, 3)) + 2 * np.eye(3)
    K = impl.ackermann(A, B, poles)
    Kz = impl.ackermann(np.linalg.solve(T, A @ T), np.linalg.solve(T, B), poles)
    npt.assert_allclose(Kz, K @ T, rtol=1e-8, atol=1e-10)                  # u = -K x = -(K T) z


def test_uncontrollable_and_unobservable_cases_raise():
    with pytest.raises(ValueError):
        impl.ackermann(np.diag([-1.0, -2.0]), np.array([[1.0], [0.0]]), np.array([-3.0, -4.0]))
    with pytest.raises(ValueError):
        impl.observer_gain(SAT[0], np.array([[0.0, 1.0]]), np.array([-1.0, -2.0]))   # gyro only: angle unobservable
    with pytest.raises(ValueError):
        impl.char_poly(np.array([-1.0 + 1.0j, -2.0]))                                 # pair incomplete


# ---------------------------------------------------------------- observer and separation

@pytest.mark.parametrize("A,C,poles", [(SAT[0], SAT[2], np.array([-0.3, -0.6])),
                                       (PITCH[0], PITCH[2], np.array([-8.0, -9.0, -10.0, -25.0]))])
def test_observer_gain_places_the_observer_poles_and_matches_place_poles(A, C, poles):
    L = impl.observer_gain(A, C, poles)
    assert L.shape == (A.shape[0], 1)
    npt.assert_allclose(sorted_c(np.linalg.eigvals(A - L @ C)), sorted_c(poles), atol=1e-6 * np.abs(poles).max())
    ref = signal.place_poles(A.T, C.T, poles).gain_matrix.T
    npt.assert_allclose(L, ref, rtol=1e-6, atol=1e-9)


def test_separation_principle_poles_are_the_union():
    A, B, C = PITCH
    reg = np.array([-4 + 4j, -4 - 4j, -6.0, -20.0])
    obs = np.array([-8.0, -9.0, -10.0, -25.0])
    K = impl.ackermann(A, B, reg)
    L = impl.observer_gain(A, C, obs)
    N = impl.reference_gain(A, B, C, K)
    Acl, Bcl = impl.observer_closed_loop(A, B, C, K, L, N)
    assert Acl.shape == (8, 8) and Bcl.shape == (8, 2)
    npt.assert_allclose(sorted_c(np.linalg.eigvals(Acl)), sorted_c(np.concatenate([reg, obs])), atol=1e-5)


def test_reference_gain_gives_zero_steady_state_error():
    A, B, C = SAT
    K = impl.ackermann(A, B, impl.second_order_poles(0.707, 0.1))
    L = impl.observer_gain(A, C, np.array([-0.3, -0.6]))
    N = impl.reference_gain(A, B, C, K)
    assert N == pytest.approx(K[0, 0])                                    # theta is the first state
    t = np.arange(0.0, 200.0, 0.1)
    run = impl.simulate_observer_loop(A, B, C, K, L, N, np.zeros(2), np.zeros(2), np.full(len(t), 0.2), np.zeros(len(t)), 0.1)
    assert run["x"][-1, 0] == pytest.approx(0.2, rel=1e-5)


def test_observer_error_follows_its_own_dynamics_whatever_the_loop_does():
    A, B, C = SAT
    L = impl.observer_gain(A, C, np.array([-0.3, -0.6]))
    dt = 0.1
    t = np.arange(0.0, 40.0, dt)
    x0, xh0 = np.array([0.01, 0.002]), np.array([-0.02, 0.0])
    runs = []
    for wn, r in ((0.1, 0.17), (0.3, -0.5)):
        K = impl.ackermann(A, B, impl.second_order_poles(0.707, wn))
        N = impl.reference_gain(A, B, C, K)
        runs.append(impl.simulate_observer_loop(A, B, C, K, L, N, x0, xh0, np.full(len(t), r), np.zeros(len(t)), dt))
    npt.assert_allclose(runs[0]["e"], runs[1]["e"], atol=1e-12)          # independent of K and r
    e_exact = np.array([expm((A - L @ C) * ti) @ (x0 - xh0) for ti in t[::50]])
    npt.assert_allclose(runs[0]["e"][::50], e_exact, atol=1e-12)


def test_perfect_initial_estimate_reproduces_state_feedback():
    A, B, C = SAT
    K = impl.ackermann(A, B, impl.second_order_poles(0.707, 0.1))
    L = impl.observer_gain(A, C, np.array([-0.3, -0.6]))
    N = impl.reference_gain(A, B, C, K)
    dt = 0.1
    t = np.arange(0.0, 60.0, dt)
    x0 = np.array([0.0, 0.001])
    run = impl.simulate_observer_loop(A, B, C, K, L, N, x0, x0, np.full(len(t), 0.1), np.zeros(len(t)), dt)
    npt.assert_allclose(run["e"], 0.0, atol=1e-14)
    Ad = expm(np.block([[A - B @ K, B * N], [np.zeros((1, 3))]]) * dt)   # state feedback alone, oracle
    z = np.concatenate([x0, [0.1]])
    for k in range(len(t)):
        assert run["x"][k] == pytest.approx(z[:2], abs=1e-12)
        z = Ad @ z


# ---------------------------------------------------------------- the independent check: fitted decay rate

@pytest.mark.parametrize("obs", [np.array([-0.3, -0.6]), np.array([-0.5, -1.5]), np.array([-0.2, -0.45])])
def test_observer_error_decays_at_the_slowest_placed_pole(obs):
    A, B, C = SAT
    K = impl.ackermann(A, B, impl.second_order_poles(0.707, 0.1))
    L = impl.observer_gain(A, C, obs)
    N = impl.reference_gain(A, B, C, K)
    dt = 0.05
    slow, fast = -obs.max(), -obs.min()
    # window: start when the fast mode is e^-12 below the slow one, stop before |e| reaches ~1e-12 of the state
    t0 = 12.0 / (fast - slow)
    t1 = t0 + 10.0 / slow
    t = np.arange(0.0, t1 + dt, dt)
    run = impl.simulate_observer_loop(A, B, C, K, L, N, np.array([0.0, 0.001]), np.zeros(2),
                                      np.full(len(t), 0.1), np.zeros(len(t)), dt)
    assert impl.fit_decay_rate(t, run["e"], t0, t1) == pytest.approx(slow, rel=1e-4)


def test_fit_decay_rate_on_a_known_exponential():
    t = np.linspace(0.0, 10.0, 201)
    assert impl.fit_decay_rate(t, 3.0 * np.exp(-0.7 * t), 0.0) == pytest.approx(0.7, rel=1e-12)
    e2 = np.column_stack([np.exp(-0.4 * t), 2 * np.exp(-0.4 * t)])
    assert impl.fit_decay_rate(t, e2, 2.0, 8.0) == pytest.approx(0.4, rel=1e-12)


# ---------------------------------------------------------------- trade-offs

def test_faster_regulator_needs_more_torque():
    A, B, C = SAT
    L = impl.observer_gain(A, C, np.array([-0.3, -0.6]))
    t = np.arange(0.0, 100.0, 0.1)
    peaks = []
    for wn in (0.05, 0.1, 0.2):
        K = impl.ackermann(A, B, impl.second_order_poles(0.707, wn))
        N = impl.reference_gain(A, B, C, K)
        run = impl.simulate_observer_loop(A, B, C, K, L, N, np.zeros(2), np.zeros(2), np.full(len(t), np.radians(10.0)), np.zeros(len(t)), 0.1)
        peaks.append(np.abs(run["u"]).max())
        assert peaks[-1] == pytest.approx(3.0 * wn ** 2 * np.radians(10.0), rel=1e-9)   # peak = N r at t = 0
    assert peaks[0] < peaks[1] < peaks[2]


def test_faster_observer_passes_more_synthetic_tracker_noise():
    A, B, C = SAT
    K = impl.ackermann(A, B, impl.second_order_poles(0.707, 0.1))
    N = impl.reference_gain(A, B, C, K)
    dt = 0.1
    t = np.arange(0.0, 300.0, dt)
    v = np.radians(10.0 / 3600.0) * np.random.default_rng(15).standard_normal(len(t))   # SYNTHETIC, seeded
    errs = []
    for wo in (0.15, 0.6):
        L = impl.observer_gain(A, C, np.array([-wo, -2 * wo]))
        run = impl.simulate_observer_loop(A, B, C, K, L, N, np.zeros(2), np.zeros(2), np.zeros(len(t)), v, dt)
        errs.append(np.std(run["e"][t >= 40.0, 1]))
    # theory: rms rate error grows roughly as wo^1.5 (factor ~8 for 4x); 2600 samples give a few-% scatter,
    # so a factor of 3 is a safe threshold that still fails if the observer ignores the noise or inverts the trend
    assert errs[1] > 3.0 * errs[0]


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "pole_placement.png").stat().st_size > 1000
