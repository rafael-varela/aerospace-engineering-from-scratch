"""Public tests for P04.06 Oscillations & Waves. They assert the closed forms against solve_ivp, the
limiting cases (zero damping, small angles, Courant 1), energy conservation of the oscillator and of the
finite-difference string, the pendulum period against the elliptic-integral result, the resonance
curve, the CFL limit, and the reduction of a SYNTHETIC pendulum dataset (labelled synthetic: it is
generated, not measured, and exists so the pipeline can be checked against a known truth).

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. The only file written goes to tmp_path.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import solve_ivp

from conftest import load_impl

impl = load_impl(__file__)

G = impl.G_STANDARD


# ---------------------------------------------------------------- SHM and damping

def test_shm_satisfies_the_equation_of_motion_and_conserves_energy():
    w = 3.0
    t = np.linspace(0.0, 5.0, 501)
    x = impl.shm(t, 0.1, -0.4, w)
    assert x[0] == pytest.approx(0.1)
    h = 1e-5
    v = (impl.shm(t + h, 0.1, -0.4, w) - impl.shm(t - h, 0.1, -0.4, w)) / (2 * h)
    acc = (impl.shm(t + h, 0.1, -0.4, w) - 2 * x + impl.shm(t - h, 0.1, -0.4, w)) / h ** 2
    assert v[0] == pytest.approx(-0.4, rel=1e-6)
    npt.assert_allclose(acc, -w ** 2 * x, atol=1e-4)
    E = impl.shm_energy(x, v, 2.0, 2.0 * w ** 2)                      # m = 2, k = m w^2
    npt.assert_allclose(E, E[0], rtol=1e-6)


def test_damped_response_matches_solve_ivp_in_all_three_regimes():
    t = np.linspace(0.0, 3.0, 301)
    for zeta in (0.0, 0.1, 1.0, 2.5):
        sol = solve_ivp(impl.damped_oscillator_rhs, (0.0, 3.0), [0.5, 1.5], args=(4.0, zeta), t_eval=t,
                        method="DOP853", rtol=1e-11, atol=1e-13)
        npt.assert_allclose(impl.damped_response(t, 0.5, 1.5, 4.0, zeta), sol.y[0], atol=1e-8)


def test_zero_damping_recovers_shm_and_damping_removes_energy():
    t = np.linspace(0.0, 10.0, 1001)
    npt.assert_allclose(impl.damped_response(t, 0.3, 0.2, 2.0, 0.0), impl.shm(t, 0.3, 0.2, 2.0), atol=1e-12)
    x = impl.damped_response(t, 0.3, 0.0, 2.0, 0.05)
    peaks = x[np.r_[False, (x[1:-1] > x[:-2]) & (x[1:-1] > x[2:]), False]]
    assert np.all(np.diff(peaks) < 0)
    delta = np.log(peaks[0] / peaks[1])
    assert impl.log_decrement_to_zeta(delta) == pytest.approx(0.05, rel=1e-3)
    assert impl.log_decrement_to_zeta(0.0) == 0.0


# ---------------------------------------------------------------- the pendulum

def test_small_angle_limit_of_the_pendulum_period():
    L = 0.8
    T0 = impl.pendulum_period_small(L)
    assert T0 == pytest.approx(2 * np.pi * np.sqrt(L / G))
    assert impl.pendulum_period_exact(L, 1e-4) == pytest.approx(T0, rel=1e-8)
    assert impl.pendulum_period_series(L, 1e-4) == pytest.approx(T0, rel=1e-8)
    assert impl.amplitude_correction(0.0) == pytest.approx(1.0)
    assert impl.pendulum_period_small(4 * L) == pytest.approx(2 * T0)   # T ~ sqrt(L)


def test_exact_series_and_simulated_periods_agree():
    L = 1.0
    for deg in (10.0, 30.0, 60.0):
        th = np.radians(deg)
        Te = impl.pendulum_period_exact(L, th)
        Ts = impl.pendulum_period_series(L, th)
        Tsim = impl.pendulum_period_simulated(L, th)
        assert Tsim == pytest.approx(Te, rel=1e-7)                       # oracle: solve_ivp with an event
        assert Ts == pytest.approx(Te, rel=2e-5 if deg <= 30 else 5e-4)  # series truncated at theta0^4
    assert impl.pendulum_period_exact(L, np.radians(30.0)) / impl.pendulum_period_small(L) - 1 == pytest.approx(0.0174, abs=2e-4)
    assert impl.pendulum_period_exact(L, np.radians(90.0)) / impl.pendulum_period_small(L) == pytest.approx(1.1803, abs=1e-3)


def test_period_grows_monotonically_with_amplitude():
    ths = np.radians(np.linspace(1.0, 150.0, 30))
    T = np.array([impl.pendulum_period_exact(1.0, th) for th in ths])
    assert np.all(np.diff(T) > 0)
    assert impl.pendulum_period_exact(1.0, np.radians(179.9)) > 5 * impl.pendulum_period_small(1.0)   # diverges toward 180 deg


def test_pendulum_rhs_small_angle_limit_and_damping_sign():
    y = np.array([0.01, 0.0])
    d = impl.pendulum_rhs(0.0, y, 0.5)
    assert d[0] == 0.0 and d[1] == pytest.approx(-(G / 0.5) * 0.01, rel=1e-4)
    assert impl.pendulum_rhs(0.0, np.array([0.0, 2.0]), 0.5, G, 0.3)[1] == pytest.approx(-0.6)


# ---------------------------------------------------------------- resonance

def test_resonance_curve_limits_and_peak():
    w0, z, f0 = 10.0, 0.05, 2.0
    assert impl.steady_state_amplitude(0.0, w0, z, f0) == pytest.approx(f0 / w0 ** 2)            # static deflection
    assert impl.steady_state_amplitude(w0, w0, z, f0) == pytest.approx(impl.quality_factor(z) * f0 / w0 ** 2)
    assert impl.quality_factor(z) == pytest.approx(10.0)
    w = np.linspace(0.5 * w0, 1.5 * w0, 20001)
    X = impl.steady_state_amplitude(w, w0, z, f0)
    assert w[np.argmax(X)] == pytest.approx(impl.resonance_peak_frequency(w0, z), rel=1e-3)
    assert impl.resonance_peak_frequency(w0, z) < w0
    assert impl.resonance_peak_frequency(w0, 0.9) == 0.0                                           # no peak when zeta >= 1/sqrt2
    # half-power points: X = X_peak / sqrt(2), separated by about 2 zeta w0
    above = w[X >= X.max() / np.sqrt(2.0)]
    assert above[-1] - above[0] == pytest.approx(impl.half_power_bandwidth(w0, z), rel=0.02)
    # high-frequency roll-off ~ 1/w^2 and zero damping gives an infinite peak
    assert impl.steady_state_amplitude(100 * w0, w0, z, f0) == pytest.approx(f0 / (100 * w0) ** 2, rel=1e-3)
    assert np.isinf(impl.steady_state_amplitude(w0, w0, 0.0, f0))


def test_phase_lag_runs_from_zero_through_90_degrees_to_180():
    w0, z = 10.0, 0.1
    assert impl.steady_state_phase(0.0, w0, z) == pytest.approx(0.0)
    assert impl.steady_state_phase(w0, w0, z) == pytest.approx(np.pi / 2)
    assert impl.steady_state_phase(1e3 * w0, w0, z) == pytest.approx(np.pi, abs=1e-3)
    w = np.linspace(0.0, 5 * w0, 500)
    assert np.all(np.diff(impl.steady_state_phase(w, w0, z)) > 0)


# ---------------------------------------------------------------- waves

def test_string_wave_speed_and_modes():
    c = impl.string_wave_speed(80.0, 5e-3)
    assert c == pytest.approx(np.sqrt(80.0 / 5e-3))
    assert impl.string_wave_speed(320.0, 5e-3) == pytest.approx(2 * c)   # four times the tension, twice the speed
    f = impl.standing_wave_frequencies(0.5, c, 5)
    assert f.shape == (5,)
    npt.assert_allclose(f / f[0], np.arange(1, 6))                         # harmonics
    assert f[0] == pytest.approx(c / 1.0)


def test_standing_wave_satisfies_the_wave_equation_and_the_boundary_conditions():
    L, c, n = 0.7, 120.0, 2
    x = np.linspace(0.0, L, 141)
    h = 1e-4
    u = impl.standing_wave(x, 0.3e-3, n, L, c)
    assert u[0] == pytest.approx(0.0, abs=1e-12) and u[-1] == pytest.approx(0.0, abs=1e-12)
    utt = (impl.standing_wave(x, 0.3e-3 + h, n, L, c) - 2 * u + impl.standing_wave(x, 0.3e-3 - h, n, L, c)) / h ** 2
    dx = x[1] - x[0]
    uxx = (u[2:] - 2 * u[1:-1] + u[:-2]) / dx ** 2
    npt.assert_allclose(utt[1:-1], c ** 2 * uxx, rtol=2e-3, atol=1e-3 * np.max(np.abs(utt)))


def test_fd_wave_solver_is_exact_at_courant_one_and_conserves_energy():
    L, c = 1.0, 2.0
    pulse = lambda x: np.exp(-((x - 0.5 * L) / 0.05) ** 2)
    x, t, U = impl.wave_equation_fd(c, L, 200, 2.0 * L / c, 1.0, pulse)   # one round trip
    assert U.shape == (len(t), 201)
    assert t[-1] == pytest.approx(2.0 * L / c)
    ref = np.array([impl.dalembert(x, ti, c, pulse, L) for ti in t])
    npt.assert_allclose(U, ref, atol=1e-12)
    npt.assert_allclose(U[-1], U[0], atol=1e-12)                           # back where it started
    half = np.argmin(np.abs(t - L / c))
    npt.assert_allclose(U[half], -U[0], atol=1e-12)                        # inverted after reflecting at both ends
    E = impl.wave_energy(U, t[1] - t[0], x[1] - x[0], c)
    assert (E.max() - E.min()) / E.mean() < 1e-12


def test_dalembert_splits_the_pulse_into_two_halves_moving_at_c():
    L, c = 2.0, 1.0
    pulse = lambda x: np.exp(-((x - 1.0) / 0.05) ** 2)
    x = np.linspace(0.0, L, 2001)
    u = impl.dalembert(x, 0.3, c, pulse, L)
    assert u.max() == pytest.approx(0.5, rel=1e-3)                        # each half has half the height
    peaks = x[np.r_[False, (u[1:-1] > u[:-2]) & (u[1:-1] > u[2:]), False]]
    assert peaks.size == 2
    npt.assert_allclose(np.sort(peaks), [1.0 - 0.3, 1.0 + 0.3], atol=2e-3)
    F = impl.odd_periodic_extension(pulse, L)
    assert F(np.array([-0.4]))[0] == pytest.approx(-pulse(np.array([0.4]))[0])   # odd about 0
    assert F(np.array([2.0 * L + 0.4]))[0] == pytest.approx(pulse(np.array([0.4]))[0])   # period 2L


def test_fd_solver_converges_at_second_order_below_courant_one_and_blows_up_above():
    L, c = 1.0, 1.0
    errs, hs = [], []
    for nx in (50, 100, 200, 400):
        x, t, U = impl.wave_equation_fd(c, L, nx, 0.6, 0.5, lambda xx: impl.standing_wave(xx, 0.0, 1, L, c))
        ref = impl.standing_wave(x, t[-1], 1, L, c)
        errs.append(np.max(np.abs(U[-1] - ref)))
        hs.append(L / nx)
    slope = np.polyfit(np.log(hs), np.log(errs), 1)[0]
    assert 1.9 < slope < 2.1
    assert impl.is_cfl_stable(1.0) and not impl.is_cfl_stable(1.01)
    x, t, Ubad = impl.wave_equation_fd(c, L, 100, 1.0, 1.1, lambda xx: np.exp(-((xx - 0.5) / 0.05) ** 2))
    assert np.max(np.abs(Ubad)) > 1e3


# ---------------------------------------------------------------- reducing pendulum data

def test_period_reduction_type_a_and_type_b_terms():
    times = np.array([20.10, 20.07, 20.13, 20.09, 20.11])
    T, u = impl.period_from_timings(times, 10)
    assert T == pytest.approx(times.mean() / 10)
    assert u == pytest.approx(times.std(ddof=1) / 10 / np.sqrt(5))
    # the Type B term is divided by the swing count but does not shrink with repeats
    T50, u50 = impl.period_from_timings(np.full(50, 20.1), 10, u_single=0.0136)
    assert u50 == pytest.approx(0.0136 / 10)
    assert T50 == pytest.approx(2.01)


def test_g_from_pendulum_and_its_uncertainty_against_finite_differences():
    L, T, uL, uT = 1.0, 2.006, 0.002, 0.002
    g, ug = impl.g_from_pendulum(L, T, uL, uT)
    assert g == pytest.approx(4 * np.pi ** 2 * L / T ** 2)
    f = lambda l, t: 4 * np.pi ** 2 * l / t ** 2
    eps = 1e-7
    dL = (f(L + eps, T) - f(L - eps, T)) / (2 * eps)
    dT = (f(L, T + eps) - f(L, T - eps)) / (2 * eps)
    assert ug == pytest.approx(np.hypot(dL * uL, dT * uT), rel=1e-6)
    assert impl.g_from_pendulum(L, T, 0.0, 0.0)[1] == 0.0


def test_synthetic_dataset_is_reduced_to_g_within_its_uncertainty():
    # SYNTHETIC data (generated at g = 9.81, theta0 = 0.15 rad): the reduction must recover g within the
    # stated uncertainty once the amplitude correction is applied, for several seeds. Frame quantisation
    # at 30 fps gives u_T ~ 1.5 ms per period over 10 swings; the tolerance of 3 u_g covers 99.7 % of draws.
    u_frame = np.sqrt(2.0) * 0.5 / (np.sqrt(3.0) * 30.0)
    corr2 = impl.amplitude_correction(0.15) ** 2
    for seed in (1, 2, 3, 4):
        lengths, times = impl.make_synthetic_pendulum_data(seed)
        assert times.shape == (lengths.size, 5)
        Ts = []
        for L, row in zip(lengths, times):
            T, uT = impl.period_from_timings(row, 10, u_single=u_frame)
            g, ug = impl.g_from_pendulum(L, T, 0.002, uT)
            assert abs(g * corr2 - 9.81) < 3 * ug
            Ts.append(T)
        A, p = impl.fit_power_law(lengths, np.array(Ts))
        assert p == pytest.approx(0.5, abs=0.01)
        assert 4 * np.pi ** 2 / A ** 2 * corr2 == pytest.approx(9.81, rel=4e-3)
    # without the amplitude correction the recovered g is biased low by about 2 theta0^2 / 16
    assert corr2 - 1 == pytest.approx(0.15 ** 2 / 8, rel=0.05)


def test_power_law_fit_recovers_exact_exponent():
    L = np.array([0.3, 0.5, 0.9, 1.4])
    A, p = impl.fit_power_law(L, 2 * np.pi * np.sqrt(L / 9.81))
    assert p == pytest.approx(0.5, abs=1e-10)
    assert A == pytest.approx(2 * np.pi / np.sqrt(9.81), rel=1e-10)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "oscillations_waves.png").exists()
