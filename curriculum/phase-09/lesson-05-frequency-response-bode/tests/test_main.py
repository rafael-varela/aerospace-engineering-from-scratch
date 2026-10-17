"""Public tests for P09.05 Frequency Response & Bode Plots.

The independent checks: the Bode magnitude and phase computed from G(j omega) against
``scipy.signal.bode`` and ``scipy.signal.freqs`` (oracles only), and against the steady-state amplitude
and phase of a SIMULATED sine response (time-domain RK4, least-squares fit). Also: the factor-sum phase
equals the unwrapped phase; low- and high-frequency limits (DC gain, type, relative degree); the classic
corner corrections of the straight-line diagram (-3.01 dB, -20 log10 2 zeta); resonance and bandwidth
formulas; the all-pass (right-half-plane zero) case; the gain crossover of a servo loop.

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

CASES = {
    "servo closed loop": ([14.31], [0.05, 1.0, 14.31]),
    "two RC stages": ([1.0], [1.6e-3 ** 2, 2 * 1.6e-3, 1.0]),
    "pitch open loop": ([720.0, 73080.0, 108000.0], [1.0, 123.2, 2400.0, 15520.0, 42800.0, 0.0]),
    "lightly damped mode": ([25.0], [1.0, 0.5, 25.0]),
    "non-minimum phase": ([-2.0, 4.0], [1.0, 3.0, 2.0]),
}


# ---------------------------------------------------------------- against the scipy oracle

@pytest.mark.parametrize("name", list(CASES))
def test_bode_matches_scipy_signal_bode(name):
    num, den = CASES[name]
    G = impl.tf(num, den)
    w = np.logspace(-2, 4, 3001)
    m, ph = impl.bode(G, w)
    _, m_sp, ph_sp = signal.bode(signal.lti(num, den), w=w)
    npt.assert_allclose(m, m_sp, atol=1e-8)
    npt.assert_allclose(ph, ph_sp, atol=1e-7)
    _, H = signal.freqs(num, den, worN=w)
    npt.assert_allclose(impl.freq_response(G, w), H, rtol=1e-10)


@pytest.mark.parametrize("name", list(CASES))
def test_factor_sum_phase_equals_unwrapped_phase(name):
    G = impl.tf(*CASES[name])
    w = np.logspace(-2, 4, 3001)
    npt.assert_allclose(impl.phase_from_factors(G, w), impl.bode(G, w)[1], atol=1e-7)


def test_unwrap_removes_only_the_360_degree_jumps():
    true = np.linspace(0.0, -700.0, 400)
    wrapped = (true + 180.0) % 360.0 - 180.0
    out = impl.unwrap_deg(wrapped)
    npt.assert_allclose(out - out[0] + true[0], true, atol=1e-9)
    smooth = np.linspace(-10.0, -170.0, 50)
    npt.assert_allclose(impl.unwrap_deg(smooth), smooth)
    assert impl.mag_db(np.array([10.0, 0.1, 1j]))[0] == pytest.approx(20.0)


# ---------------------------------------------------------------- against the simulated sine test

@pytest.mark.parametrize("f_hz", [0.5, 2.0, 6.0])
def test_simulated_sine_response_matches_the_bode_servo(f_hz):
    G = impl.servo_closed_loop()
    w = 2 * np.pi * f_hz
    gain, ph = impl.sine_test(G, w)
    H = impl.freq_response(G, np.array([w]))[0]
    assert gain == pytest.approx(abs(H), rel=1e-5)
    assert ph == pytest.approx(np.degrees(np.angle(H)), abs=1e-3)


def test_simulated_sine_response_at_resonance_and_on_the_gyro_filter():
    G = impl.tf(*CASES["lightly damped mode"])           # zeta = 0.05, wn = 5: peak of 10 at resonance
    w = 5.0 * np.sqrt(1 - 2 * 0.05 ** 2)
    gain, ph = impl.sine_test(G, w)
    assert gain == pytest.approx(abs(impl.freq_response(G, np.array([w]))[0]), rel=1e-4)
    assert gain == pytest.approx(impl.resonance(0.05, 5.0)[1], rel=1e-4)
    H1 = impl.rc_filter(1.6e-3, 1)
    w1200 = 2 * np.pi * 1200.0
    g, p = impl.sine_test(H1, w1200)
    assert g == pytest.approx(1 / np.sqrt(1 + (w1200 * 1.6e-3) ** 2), rel=1e-5)
    assert p == pytest.approx(-np.degrees(np.arctan(w1200 * 1.6e-3)), abs=1e-3)


def test_fit_sinusoid_recovers_amplitude_and_phase_with_noise():
    t = np.linspace(0.0, 10.0, 5001)
    w = 2 * np.pi * 1.3
    A, ph = impl.fit_sinusoid(t, 2.5 * np.cos(w * t - 0.7) + 0.4, w)
    assert A == pytest.approx(2.5, rel=1e-12) and ph == pytest.approx(np.degrees(-0.7), abs=1e-10)
    rng = np.random.default_rng(5)
    y = 2.5 * np.cos(w * t - 0.7) + 0.1 * rng.standard_normal(t.size)
    A, ph = impl.fit_sinusoid(t, y, w)
    # least squares on N = 5001 points: sd(A) ~ sigma sqrt(2/N) = 0.002; allow 5 sd
    assert A == pytest.approx(2.5, abs=0.01)
    assert ph == pytest.approx(np.degrees(-0.7), abs=np.degrees(0.01 / 2.5))


# ---------------------------------------------------------------- limits and the straight-line diagram

def test_low_and_high_frequency_limits():
    G = impl.tf(*CASES["servo closed loop"])
    m, ph = impl.bode(G, np.array([1e-4, 1e5, 1e6]))
    assert m[0] == pytest.approx(0.0, abs=1e-6) and ph[0] == pytest.approx(0.0, abs=1e-3)
    assert m[2] - m[1] == pytest.approx(-40.0, abs=0.01)            # relative degree 2: -40 dB/decade
    assert ph[2] == pytest.approx(-180.0, abs=0.01)
    L = impl.tf(*CASES["pitch open loop"])
    mL, phL = impl.bode(L, np.array([1e-4, 1e-3]))
    assert mL[1] - mL[0] == pytest.approx(-20.0, abs=1e-3)           # type 1: -20 dB/decade, -90 deg
    assert phL[0] == pytest.approx(-90.0, abs=0.01)
    bf = impl.bode_form(L)
    assert bf["type"] == 1
    assert bf["K"] == pytest.approx(108000.0 / 42800.0, rel=1e-9)    # lim s L(s) as s -> 0


@pytest.mark.parametrize("zeta", [0.05, 0.2, 0.5])
def test_corner_corrections_of_the_straight_line_diagram(zeta):
    wn, wb = 10.0, 3.0
    P = impl.tf([wn ** 2], [1.0, 2 * zeta * wn, wn ** 2])
    err = impl.bode(P, np.array([wn]))[0][0] - impl.asymptotic_mag_db(P, np.array([wn]))[0]
    assert err == pytest.approx(-20 * np.log10(2 * zeta), abs=1e-9)
    R = impl.tf([1.0], [1 / wb, 1.0])
    w = np.array([wb / 10, wb, 10 * wb])
    e = impl.bode(R, w)[0] - impl.asymptotic_mag_db(R, w)
    npt.assert_allclose(e, [-0.0432, -3.0103, -0.0432], atol=1e-4)
    ww = np.logspace(-2, 3, 2001)
    assert np.max(np.abs(impl.asymptotic_phase_deg(R, ww) - impl.bode(R, ww)[1])) == pytest.approx(5.71, abs=0.02)


def test_straight_line_slopes_follow_the_bode_form():
    L = impl.tf(*CASES["pitch open loop"])
    w = np.array([0.01, 0.1, 1e4, 1e5])
    a = impl.asymptotic_mag_db(L, w)
    assert a[1] - a[0] == pytest.approx(-20.0)
    assert a[3] - a[2] == pytest.approx(-20.0 * 3)                   # 5 poles - 2 zeros
    assert a[0] == pytest.approx(20 * np.log10(108000.0 / 42800.0) + 40.0)
    npt.assert_allclose(impl.asymptotic_phase_deg(L, np.array([1e-3, 1e5])), [-90.0, -270.0])


def test_rhp_zero_has_the_magnitude_of_its_mirror_and_extra_lag():
    w = np.logspace(-2, 2, 401)
    ap = impl.tf([-1.0, 1.0], [1.0, 1.0])                            # all-pass (1 - s)/(1 + s)
    m, ph = impl.bode(ap, w)
    npt.assert_allclose(m, 0.0, atol=1e-12)
    npt.assert_allclose(ph, -2 * np.degrees(np.arctan(w)), atol=1e-9)
    lhp = impl.tf([1.0, 2.0], [1.0, 3.0, 2.0])
    rhp = impl.tf([-1.0, 2.0], [1.0, 3.0, 2.0])
    npt.assert_allclose(impl.bode(lhp, w)[0], impl.bode(rhp, w)[0], atol=1e-12)
    assert np.all(impl.bode(rhp, w)[1][1:] < impl.bode(lhp, w)[1][1:])


# ---------------------------------------------------------------- resonance, bandwidth, crossover

@pytest.mark.parametrize("zeta", [0.1, 0.3, 0.6])
def test_resonance_formula_matches_the_numerical_peak(zeta):
    wn = 4.0
    G = impl.tf([wn ** 2], [1.0, 2 * zeta * wn, wn ** 2])
    w = np.linspace(0.01, 10.0, 200001)
    mag = np.abs(impl.freq_response(G, w))
    wr, Mr = impl.resonance(zeta, wn)
    assert w[np.argmax(mag)] == pytest.approx(wr, abs=1e-3)
    assert mag.max() == pytest.approx(Mr, rel=1e-8)
    bw = wn * np.sqrt(1 - 2 * zeta ** 2 + np.sqrt(4 * zeta ** 4 - 4 * zeta ** 2 + 2))
    assert impl.bandwidth(G) == pytest.approx(bw, rel=1e-7)
    assert impl.resonance(0.8, wn) == (0.0, 1.0)


def test_bandwidth_of_first_order_and_crossover_of_the_servo_loop():
    tau = 0.05
    assert impl.bandwidth(impl.tf([1.0], [tau, 1.0])) == pytest.approx(1 / tau, rel=1e-8)
    K = 14.31
    open_loop = impl.tf([K], [tau, 1.0, 0.0])
    wc, ph = impl.gain_crossover(open_loop)
    # |K / (j w (j w tau + 1))| = 1  ->  w^2 (1 + w^2 tau^2) = K^2
    wc_exact = np.sqrt((-1 + np.sqrt(1 + 4 * K ** 2 * tau ** 2)) / (2 * tau ** 2))
    assert wc == pytest.approx(wc_exact, rel=1e-8)
    assert ph == pytest.approx(-90.0 - np.degrees(np.arctan(wc_exact * tau)), abs=1e-6)
    assert abs(impl.freq_response(open_loop, np.array([wc]))[0]) == pytest.approx(1.0, rel=1e-8)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "bode.png").stat().st_size > 1000
