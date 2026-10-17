"""Public tests for P02.03 Forced Response & Resonance. They assert properties, with scipy as the oracle.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-02/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import solve_ivp

from conftest import load_impl

impl = load_impl(__file__)

P = impl.avionics_bracket()
M, C, K, F0 = P["m"], P["c"], P["k"], P["F0"]
WN = np.sqrt(K / M)
ZETA = C / (2 * np.sqrt(K * M))


def residual(x_of_t, force_of_t, t, m, c, k, dt=1e-6):
    """m x'' + c x' + k x - F(t) with derivatives by central differences."""
    x = x_of_t(t)
    xp, xm = x_of_t(t + dt), x_of_t(t - dt)
    v = (xp - xm) / (2 * dt)
    a = (xp - 2 * x + xm) / dt ** 2
    return m * a + c * v + k * x - force_of_t(t)


# ---------------------------------------------------------------- steady state

def test_steady_state_solves_the_forced_equation_at_several_frequencies():
    for f in (20.0, 54.0, 55.0, 70.0):
        w = 2 * np.pi * f
        t = np.linspace(0.011, 0.3, 7)
        res = residual(lambda tt: impl.steady_state_response(tt, w, F0, M, C, K), lambda tt: F0 * np.cos(w * tt), t, M, C, K)
        assert np.max(np.abs(res)) < 1e-4 * F0                                      # second difference with dt = 1e-6 on um-scale x


def test_undetermined_coefficients_agree_with_the_complex_gain():
    for w in (10.0, WN, 3 * WN):
        A, B = impl.undetermined_coefficients(w, F0, M, C, K)
        H = impl.complex_gain(w, M, C, K)
        assert A == pytest.approx(F0 * H.real, abs=1e-15)
        assert B == pytest.approx(-F0 * H.imag, abs=1e-15)                           # Re[H e^{iwt}] = Re H cos - Im H sin
        # A cos + B sin reproduces the steady-state function.
        t = np.linspace(0.0, 0.1, 50)
        npt.assert_allclose(A * np.cos(w * t) + B * np.sin(w * t), impl.steady_state_response(t, w, F0, M, C, K), atol=1e-15)


def test_gain_limits_at_low_and_high_frequency_and_at_resonance():
    assert float(impl.magnification(0.0, M, C, K)) == pytest.approx(1.0)              # static: X = F0 / k
    assert float(impl.phase_lag(0.0, M, C, K)) == pytest.approx(0.0)
    w_hi = 50 * WN
    assert float(impl.magnification(w_hi, M, C, K)) == pytest.approx(K / (M * w_hi ** 2), rel=1e-3)   # mass-controlled: 1/r^2
    assert float(impl.phase_lag(w_hi, M, C, K)) == pytest.approx(np.pi, abs=0.01)
    assert float(impl.phase_lag(WN, M, C, K)) == pytest.approx(np.pi / 2)                              # exactly 90 deg at omega_n
    assert float(impl.magnification(WN, M, C, K)) == pytest.approx(1 / (2 * ZETA))
    with pytest.raises(ValueError):
        impl.complex_gain(WN, M, 0.0, K)


def test_magnification_formula_in_r_and_zeta():
    r = np.linspace(0.1, 3.0, 30)
    w = r * WN
    expected = 1.0 / np.sqrt((1 - r ** 2) ** 2 + (2 * ZETA * r) ** 2)
    npt.assert_allclose(impl.magnification(w, M, C, K), expected, rtol=1e-12)
    # Phase lag is monotone increasing in frequency, from 0 to pi.
    phi = impl.phase_lag(w, M, C, K)
    assert np.all(np.diff(phi) > 0.0) and 0 < phi[0] < phi[-1] < np.pi


# ---------------------------------------------------------------- the peak

def test_peak_location_height_bandwidth_and_q():
    wp = impl.peak_frequency(M, C, K)
    assert wp == pytest.approx(WN * np.sqrt(1 - 2 * ZETA ** 2))
    # The analytic peak beats every sampled frequency.
    w = np.linspace(0.5 * WN, 1.5 * WN, 20001)
    Mw = impl.magnification(w, M, C, K)
    assert impl.peak_magnification(M, C, K) >= np.max(Mw) - 1e-9
    assert impl.peak_magnification(M, C, K) == pytest.approx(np.max(Mw), rel=1e-6)
    lo, hi = impl.half_power_bandwidth(M, C, K)
    assert lo < wp < hi
    for wb in (lo, hi):
        assert float(impl.magnification(wb, M, C, K)) ** 2 == pytest.approx(0.5 * impl.peak_magnification(M, C, K) ** 2, rel=1e-8)
    assert hi - lo == pytest.approx(2 * ZETA * WN, rel=0.01)                          # light-damping rule
    assert impl.quality_factor(M, C, K) == pytest.approx(1 / (2 * ZETA), rel=0.01)
    # Heavy damping: no peak.
    c_heavy = 2 * 0.8 * np.sqrt(K * M)
    assert impl.peak_frequency(M, c_heavy, K) == 0.0 and impl.peak_magnification(M, c_heavy, K) == 1.0
    with pytest.raises(ValueError):
        impl.half_power_bandwidth(M, c_heavy, K)


# ---------------------------------------------------------------- transients

def test_full_response_meets_initial_conditions_and_decays_to_steady_state():
    w = 2 * np.pi * 57.0
    x0, v0 = 2e-3, -0.1
    t = np.linspace(0.0, 1.5, 30001)
    x = impl.full_response(t, w, F0, M, C, K, x0, v0)
    assert x[0] == pytest.approx(x0, abs=1e-12)
    dt = 1e-7
    v_est = float((impl.full_response(dt, w, F0, M, C, K, x0, v0) - impl.full_response(-dt, w, F0, M, C, K, x0, v0)) / (2 * dt))
    assert v_est == pytest.approx(v0, rel=1e-5)
    T99 = impl.transient_decay_time(M, C, K, 0.01)
    assert T99 == pytest.approx(-np.log(0.01) / (ZETA * WN))
    late = t > 2.0 * T99  # transient is below 1e-4 of its initial size here
    npt.assert_allclose(x[late], impl.steady_state_response(t[late], w, F0, M, C, K), atol=1e-3 * F0 * abs(impl.complex_gain(w, M, C, K)))
    # The residual of the full response vanishes too.
    res = residual(lambda tt: impl.full_response(tt, w, F0, M, C, K, x0, v0), lambda tt: F0 * np.cos(w * tt), np.linspace(0.02, 0.2, 5), M, C, K)
    assert np.max(np.abs(res)) < 1e-3 * F0


def test_full_response_matches_rk4_and_solve_ivp():
    w = WN
    t_end, n = 0.5, 10000
    ts, ys = impl.rk4_forced(M, C, K, lambda tt: F0 * np.cos(w * tt), 0.0, 0.0, t_end, n)
    x_exact = impl.full_response(ts, w, F0, M, C, K, 0.0, 0.0)
    assert np.max(np.abs(ys[:, 0] - x_exact)) < 1e-9
    sol = solve_ivp(lambda tt, y: [y[1], (F0 * np.cos(w * tt) - C * y[1] - K * y[0]) / M], (0.0, t_end), [0.0, 0.0],
                    rtol=1e-10, atol=1e-13, t_eval=ts)
    assert np.max(np.abs(sol.y[0] - x_exact)) < 1e-10
    # RK4 is fourth order on this problem.
    errs = []
    for nn in (1000, 2000, 4000):
        tss, yss = impl.rk4_forced(M, C, K, lambda tt: F0 * np.cos(w * tt), 0.0, 0.0, 0.1, nn)
        errs.append(np.max(np.abs(yss[:, 0] - impl.full_response(tss, w, F0, M, C, K, 0.0, 0.0))))
    rates = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert np.all(np.abs(rates - 4.0) < 0.4), rates
    with pytest.raises(ValueError):
        impl.rk4_forced(M, C, K, lambda tt: 0.0, 0.0, 0.0, 1.0, 0)


def test_superposition_of_two_drives():
    w1, w2 = 2 * np.pi * 50.0, 2 * np.pi * 60.0
    t = np.linspace(0.0, 0.2, 2001)
    x12 = impl.steady_state_response(t, w1, F0, M, C, K) + impl.steady_state_response(t, w2, 0.5 * F0, M, C, K)
    ts, ys = impl.rk4_forced(M, C, K, lambda tt: F0 * np.cos(w1 * tt) + 0.5 * F0 * np.cos(w2 * tt),
                             float(x12[0]), float((x12[1] - x12[0]) / (t[1] - t[0])), 0.2, 2000)
    # Started on the combined steady state, the system stays on it (to the accuracy of the velocity estimate).
    assert np.max(np.abs(ys[:, 0] - x12)) < 2e-2 * np.max(np.abs(x12))


def test_undamped_resonance_grows_linearly_and_solves_the_equation():
    t = np.linspace(0.0, 0.5, 5001)
    x = impl.undamped_resonance(t, F0, M, K)
    res = residual(lambda tt: impl.undamped_resonance(tt, F0, M, K), lambda tt: F0 * np.cos(WN * tt), np.linspace(0.05, 0.4, 6), M, 0.0, K)
    assert np.max(np.abs(res)) < 1e-3 * F0
    env = F0 / (2 * M * WN) * t
    assert np.all(np.abs(x) <= env + 1e-15)                                           # inside the linear envelope
    # The envelope is reached: the peaks grow linearly with time.
    peaks = [np.max(np.abs(x[(t > a) & (t < a + 0.05)])) for a in (0.1, 0.2, 0.3, 0.4)]
    ratios = np.array(peaks[1:]) / np.array(peaks[:-1])
    expected = np.array([0.25, 0.35, 0.45]) / np.array([0.15, 0.25, 0.35])
    npt.assert_allclose(ratios, expected, rtol=0.1)
    # And it agrees with RK4 on the undamped system driven at omega_n.
    ts, ys = impl.rk4_forced(M, 0.0, K, lambda tt: F0 * np.cos(WN * tt), 0.0, 0.0, 0.5, 20000)
    assert np.max(np.abs(ys[:, 0] - impl.undamped_resonance(ts, F0, M, K))) < 1e-8


def test_beats_have_the_right_envelope_and_tend_to_resonance():
    w = 2 * np.pi * 52.0
    t = np.linspace(0.0, 1.0, 40001)
    x = impl.beats(t, w, F0, M, K)
    assert x[0] == pytest.approx(0.0, abs=1e-15)
    res = residual(lambda tt: impl.beats(tt, w, F0, M, K), lambda tt: F0 * np.cos(w * tt), np.linspace(0.05, 0.9, 6), M, 0.0, K)
    assert np.max(np.abs(res)) < 1e-3 * F0
    env = 2 * F0 / (M * abs(WN ** 2 - w ** 2)) * np.abs(np.sin((WN - w) * t / 2))
    assert np.all(np.abs(x) <= env + 1e-12)
    # Envelope zeros every 2 pi / |wn - w| = 1/3 s: the response is small there.
    Tb = 2 * np.pi / abs(WN - w)
    for n in (1, 2):
        i = np.argmin(np.abs(t - n * Tb))
        assert abs(x[i]) < 0.02 * np.max(np.abs(x))
    # Near omega_n, beats approach the t sin(wn t) solution over the first beat half-period.
    w_close = WN * (1 + 1e-4)
    t_short = np.linspace(0.0, 0.3, 3001)
    xr = impl.undamped_resonance(t_short, F0, M, K)
    npt.assert_allclose(impl.beats(t_short, w_close, F0, M, K), xr, atol=0.02 * np.max(np.abs(xr)))
    with pytest.raises(ValueError):
        impl.beats(t, WN, F0, M, K)


def test_step_response_final_value_and_overshoot():
    t = np.linspace(0.0, 3.0, 300001)
    x = impl.step_response(t, F0, M, C, K)
    assert x[0] == pytest.approx(0.0, abs=1e-15)
    assert x[-1] == pytest.approx(F0 / K, rel=1e-3)
    over = (np.max(x) - F0 / K) / (F0 / K)
    assert over == pytest.approx(impl.step_overshoot(ZETA), rel=1e-3)
    assert impl.step_overshoot(0.0) == 1.0
    assert impl.step_overshoot(0.5) == pytest.approx(np.exp(-np.pi * 0.5 / np.sqrt(0.75)))
    # Residual of the step response against a constant force.
    res = residual(lambda tt: impl.step_response(tt, F0, M, C, K), lambda tt: F0 * np.ones_like(np.asarray(tt)), np.linspace(0.01, 0.1, 5), M, C, K)
    assert np.max(np.abs(res)) < 1e-3 * F0
    with pytest.raises(ValueError):
        impl.step_overshoot(1.0)


# ---------------------------------------------------------------- base excitation

def test_transmissibility_crosses_one_at_root_two_for_every_damping():
    w_cross = np.sqrt(2.0) * WN
    for zeta in (0.02, 0.1, 0.3, 0.7):
        c = 2 * zeta * np.sqrt(K * M)
        assert float(impl.transmissibility(w_cross, M, c, K)) == pytest.approx(1.0, rel=1e-10)
        assert float(impl.transmissibility(0.0, M, c, K)) == pytest.approx(1.0)
    # Below the crossing more damping helps; above it more damping hurts.
    c_lo, c_hi = 2 * 0.05 * np.sqrt(K * M), 2 * 0.3 * np.sqrt(K * M)
    assert float(impl.transmissibility(WN, M, c_lo, K)) > float(impl.transmissibility(WN, M, c_hi, K))
    assert float(impl.transmissibility(3 * WN, M, c_lo, K)) < float(impl.transmissibility(3 * WN, M, c_hi, K))
    r = np.linspace(0.1, 4.0, 40)
    expected = np.sqrt((1 + (2 * ZETA * r) ** 2) / ((1 - r ** 2) ** 2 + (2 * ZETA * r) ** 2))
    npt.assert_allclose(impl.transmissibility(r * WN, M, C, K), expected, rtol=1e-12)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "resonance_explorer.png").stat().st_size > 1000
