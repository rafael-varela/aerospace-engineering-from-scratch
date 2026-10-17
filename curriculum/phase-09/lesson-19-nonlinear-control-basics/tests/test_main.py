"""Tests for P09.19 Nonlinear Control Basics.

Run against the starter (default) or the reference solution (AEFS_IMPL=solution).
The independent checks: the Lyapunov-equation residual and ``scipy.linalg.solve_continuous_lyapunov``; the
describing-function limit-cycle amplitude and frequency against a simulation of the nonlinear loop.
"""
from __future__ import annotations

import numpy as np
import numpy.testing as npt
import pytest
from conftest import load_impl
from scipy.integrate import solve_ivp
from scipy.linalg import expm, solve_continuous_lyapunov

impl = load_impl(__file__)

RNG = np.random.default_rng(20250919)
P_PEND = impl.pendulum_params()
K_PEND = impl.pendulum_gains(P_PEND)
NG, DG = impl.servo_loop_plant()


# ---------------------------------------------------------------- Lyapunov equation

@pytest.mark.parametrize("n", [2, 3, 6])
def test_lyapunov_solve_residual_and_scipy(n):
    A = RNG.standard_normal((n, n))
    A -= (np.max(np.linalg.eigvals(A).real) + 0.5) * np.eye(n)
    C = RNG.standard_normal((n, n))
    Q = C @ C.T + 0.1 * np.eye(n)
    P = impl.lyapunov_solve(A, Q)
    npt.assert_allclose(A.T @ P + P @ A + Q, 0.0, atol=1e-10 * np.abs(Q).max())
    npt.assert_allclose(P, solve_continuous_lyapunov(A.T, -Q), rtol=1e-9, atol=1e-12)
    assert np.all(np.linalg.eigvalsh(P) > 0)                         # Hurwitz A, Q > 0  =>  P > 0


def test_lyapunov_singular_and_unstable_cases():
    with pytest.raises(ValueError):
        impl.lyapunov_solve(np.array([[0.0, 1.0], [-1.0, 0.0]]), np.eye(2))      # eigenvalues +-j sum to 0
    P = impl.lyapunov_solve(np.array([[0.5, 0.0], [0.0, -1.0]]), np.eye(2))     # unstable A: P not positive definite
    assert np.min(np.linalg.eigvalsh(P)) < 0


def test_v_decreases_at_rate_xQx_along_linear_trajectories():
    A = np.array([[0.0, 1.0], [-16.0, -5.6]])
    Q = np.diag([1.0, 0.1])
    P = impl.lyapunov_solve(A, Q)
    x0 = np.array([0.4, -1.0])
    h = 1e-5
    for t in (0.0, 0.3, 1.1):
        V = lambda tt: (expm(A * tt) @ x0) @ P @ (expm(A * tt) @ x0)
        x = expm(A * t) @ x0
        npt.assert_allclose((V(t + h) - V(t - h)) / (2 * h), -x @ Q @ x, rtol=1e-6)
        npt.assert_allclose(impl.lyapunov_derivative(lambda xx: A @ xx, P, x), -x @ Q @ x, rtol=1e-12)


# ---------------------------------------------------------------- linearisation and the region of attraction

@pytest.mark.parametrize("theta0", [0.0, 0.7, np.pi])
def test_jacobian_matches_the_analytic_linearisation(theta0):
    ml2 = P_PEND["m"] * P_PEND["l"] ** 2
    J = impl.jacobian_fd(lambda x: impl.pendulum_rhs(x, 0.3, P_PEND), np.array([theta0, 0.2]))
    ref = np.array([[0.0, 1.0], [P_PEND["g"] / P_PEND["l"] * np.cos(theta0), -P_PEND["b"] / ml2]])
    npt.assert_allclose(J, ref, rtol=1e-7, atol=1e-8)


def test_upright_is_a_saddle_and_the_controller_places_the_poles():
    J0 = impl.jacobian_fd(lambda x: impl.pendulum_rhs(x, 0.0, P_PEND), np.zeros(2))
    lam = np.linalg.eigvals(J0)
    assert lam.real.max() > 0 > lam.real.min()
    Acl = impl.jacobian_fd(impl.saturated_feedback_rhs(P_PEND, K_PEND), np.zeros(2))
    npt.assert_allclose(np.sort_complex(np.linalg.eigvals(Acl)), np.sort_complex(np.roots([1.0, 5.6, 16.0])), atol=1e-6)


@pytest.fixture(scope="module")
def roa():
    f = impl.saturated_feedback_rhs(P_PEND, K_PEND)
    Acl = impl.jacobian_fd(f, np.zeros(2))
    P = impl.lyapunov_solve(Acl, np.diag([1.0, 0.1]))
    return f, P, impl.roa_estimate(f, P)


def test_roa_estimate_has_decreasing_v_inside(roa):
    f, P, c = roa
    assert c > 0
    L = np.linalg.cholesky(np.linalg.inv(P))
    for phi in np.linspace(0, 2 * np.pi, 72, endpoint=False):
        for frac in (0.1, 0.5, 0.95):
            x = L @ (np.sqrt(frac * c) * np.array([np.cos(phi), np.sin(phi)]))
            assert impl.lyapunov_derivative(f, P, x) < 0


def test_starts_inside_the_estimate_converge_and_the_estimate_is_conservative(roa):
    f, P, c = roa
    L = np.linalg.cholesky(np.linalg.inv(P))
    for _ in range(25):
        z = RNG.standard_normal(2)
        x0 = L @ (z / np.linalg.norm(z) * np.sqrt(c) * 0.99 * RNG.uniform() ** 0.5)
        sol = solve_ivp(lambda t, x: f(x), (0.0, 8.0), x0, rtol=1e-8, atol=1e-10, max_step=0.05)
        assert np.linalg.norm(sol.y[:, -1]) < 1e-3
    # at rest the estimate reaches less than the static limit asin(u_max / m g l): it is inside the true region
    theta_rest = np.sqrt(c / P[0, 0])
    assert theta_rest < np.arcsin(P_PEND["u_max"] / (P_PEND["m"] * P_PEND["g"] * P_PEND["l"]))


def test_more_torque_gives_a_larger_region():
    cs = []
    for u_max in (3.0, 6.0, 12.0):
        p = dict(P_PEND, u_max=u_max)
        f = impl.saturated_feedback_rhs(p, K_PEND)
        P = impl.lyapunov_solve(impl.jacobian_fd(f, np.zeros(2)), np.diag([1.0, 0.1]))
        cs.append(impl.roa_estimate(f, P, n_dir=180, n_r=300))
    assert cs[0] < cs[1] < cs[2]


# ---------------------------------------------------------------- describing functions

def _saturation(delta, k):
    return lambda th, A: k * np.clip(A * np.sin(th), -delta, delta)


def _relay_hyst(M, h):
    def y(th, A):
        on = np.arcsin(h / A)                       # rising input reaches +h: switch to +M; falling reaches -h: to -M
        return np.where((th >= on) & (th < np.pi + on), M, -M)
    return y


@pytest.mark.parametrize("A", [0.5, 1.0, 1.5, 3.0, 10.0])
def test_saturation_describing_function_against_fourier_coefficient(A):
    delta, k = 1.0, 2.5
    N = impl.df_saturation(A, delta, k)[0]
    npt.assert_allclose(N, impl.first_harmonic_gain(_saturation(delta, k), A).real, rtol=1e-6)
    assert N <= k + 1e-15
    if A <= delta:
        assert N == k


def test_saturation_df_limits_and_monotonicity():
    A = np.logspace(0, 3, 50)
    N = impl.df_saturation(A, 1.0, 1.0)
    assert np.all(np.diff(N) < 0)
    npt.assert_allclose(N[-1], 4.0 / (np.pi * A[-1]), rtol=1e-6)    # large A: saturation looks like a relay of height k delta


@pytest.mark.parametrize("h", [0.0, 0.1, 0.3])
def test_relay_describing_function_against_fourier_coefficient(h):
    for A in (0.5, 1.0, 2.0):
        N = impl.df_relay(A, 2.0, h)[0]
        npt.assert_allclose(N, impl.first_harmonic_gain(_relay_hyst(2.0, h), A, n=200000), rtol=2e-4, atol=1e-5)
    with pytest.raises(ValueError):
        impl.df_relay(0.05, 1.0, 0.1)


def test_phase_crossover_and_measurement_helpers():
    npt.assert_allclose(impl.phase_crossover(NG, DG), 1.0 / np.sqrt(0.04 * 0.01), rtol=1e-9)
    t = np.linspace(0.0, 5.0, 50001)
    amp, w = impl.oscillation_amplitude_frequency(t, 0.7 * np.sin(13.0 * t + 0.3) + 0.2, 1.0)
    npt.assert_allclose([amp, w], [0.7, 13.0], rtol=1e-5)


@pytest.fixture(scope="module")
def relay_sim():
    t, e = impl.simulate_relay_loop(NG, DG, 1.0, t_end=2.5)
    th, eh = impl.simulate_relay_loop(NG, DG, 1.0, hysteresis=0.02, t_end=2.5)
    return impl.oscillation_amplitude_frequency(t, e, 1.8), impl.oscillation_amplitude_frequency(th, eh, 1.8)


def test_relay_limit_cycle_prediction_matches_simulation(relay_sim):
    """Harmonic balance neglects the higher harmonics of the square wave; for this low-pass plant the
    prediction is within 6 % in amplitude and 4 % in frequency."""
    A_p, w_p = impl.predict_relay_limit_cycle(NG, DG, 1.0)
    npt.assert_allclose(A_p, 4.0 * 0.24 / np.pi, rtol=1e-6)            # |G(j50)| = 30 * 0.0004 / 0.05
    (A_s, w_s), _ = relay_sim
    assert abs(A_p / A_s - 1) < 0.06
    assert abs(w_p / w_s - 1) < 0.04


def test_relay_with_hysteresis_prediction_matches_simulation(relay_sim):
    A_p, w_p = impl.predict_relay_hysteresis_limit_cycle(NG, DG, 1.0, 0.02)
    G = np.polyval(NG, 1j * w_p) / np.polyval(DG, 1j * w_p)
    npt.assert_allclose(G, -1.0 / impl.df_relay(A_p, 1.0, 0.02)[0], rtol=1e-8)   # harmonic balance holds exactly
    _, (A_s, w_s) = relay_sim
    assert abs(A_p / A_s - 1) < 0.06
    assert abs(w_p / w_s - 1) < 0.04
    A0, w0 = impl.predict_relay_limit_cycle(NG, DG, 1.0)
    assert A_p > A0 and w_p < w0                                        # hysteresis: bigger, slower cycle


def test_saturating_amplifier_limit_cycle_prediction_matches_simulation():
    k_crit = 1.0 / 0.24
    A_p, w_p = impl.predict_saturation_limit_cycle(NG, DG, 2.0 * k_crit, 1.0)
    npt.assert_allclose(impl.df_saturation(A_p, 1.0 / (2.0 * k_crit), 2.0 * k_crit)[0] * 0.24, 1.0, rtol=1e-9)
    t, e = impl.simulate_saturation_loop(NG, DG, 2.0 * k_crit, 1.0, t_end=3.0)
    A_s, w_s = impl.oscillation_amplitude_frequency(t, e, 2.2)
    assert abs(A_p / A_s - 1) < 0.05
    assert abs(w_p / w_s - 1) < 0.04


def test_below_critical_gain_no_limit_cycle_is_predicted_and_none_appears():
    with pytest.raises(ValueError):
        impl.predict_saturation_limit_cycle(NG, DG, 0.8 / 0.24, 1.0)
    t, e = impl.simulate_saturation_loop(NG, DG, 0.5 / 0.24, 1.0, t_end=2.0, y0=0.2)   # gain margin 2
    assert np.abs(e[t > 1.5]).max() < 2e-3 * 0.2


# ---------------------------------------------------------------- feedback linearisation

def test_feedback_linearising_torque_cancels_the_nonlinearity():
    for _ in range(10):
        x = RNG.uniform([-2.0, -5.0], [2.0, 5.0])
        v = RNG.uniform(-10.0, 10.0)
        u = impl.feedback_linearising_torque(x, v, P_PEND)
        npt.assert_allclose(impl.pendulum_rhs(x, u, P_PEND)[1], v, rtol=1e-12, atol=1e-12)


def test_feedback_linearised_loop_is_exactly_linear_at_large_angle():
    wn, zeta = 4.0, 0.7
    ctrl = lambda x: impl.feedback_linearising_torque(x, -(wn ** 2) * x[0] - 2 * zeta * wn * x[1], P_PEND)
    x0 = np.array([np.radians(60.0), 0.5])
    t, x, u = impl.simulate_pendulum(x0, ctrl, P_PEND, 3.0, n_out=61)
    Alin = np.array([[0.0, 1.0], [-wn ** 2, -2 * zeta * wn]])
    ref = np.array([expm(Alin * ti) @ x0 for ti in t])
    npt.assert_allclose(x, ref, atol=1e-8)
    # the linear (gravity-blind) controller from the same start does not follow the linear model
    t, xl, _ = impl.simulate_pendulum(x0, lambda xx: float(-K_PEND @ xx), P_PEND, 3.0, n_out=61)
    assert np.abs(xl - ref).max() > 1e-2


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "nonlinear.png").stat().st_size > 1000
