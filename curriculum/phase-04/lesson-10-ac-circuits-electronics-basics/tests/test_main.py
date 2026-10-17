"""Public tests for P04.10 AC Circuits & Electronics Basics. They assert physics and circuit laws, not
magic numbers: impedances have the right phase and frequency dependence; the series RLC is purely
resistive at f0 = 1/(2 pi sqrt(LC)) with bandwidth f0/Q; RMS of a sine is A/sqrt(2) and the power
factor is cos(arg Z); the RC corner frequency is 1/(2 pi RC) at -3 dB and -45 degrees; the buffered
cascade matches ``scipy.signal.freqs`` and its step response matches ``scipy.signal.step`` and the
closed form; the unbuffered ladder matches nodal analysis; the closed-loop gain tends to the ideal
op-amp result; logic gates obey their truth tables; aliasing and quantisation behave as derived;
and the alias tone in SYNTHETIC sampled data is reduced by the filter's attenuation.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-04/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy import signal

from conftest import load_impl

impl = load_impl(__file__)


# ---------------------------------------------------------------- impedance and resonance

def test_impedances_have_the_right_phase_and_frequency_dependence():
    f = np.array([10.0, 100.0, 1000.0])
    ZR, ZL, ZC = impl.impedance_R(50.0, f), impl.impedance_L(1e-3, f), impl.impedance_C(1e-6, f)
    npt.assert_allclose(ZR, 50.0)
    npt.assert_allclose(np.angle(ZL), np.pi / 2)
    npt.assert_allclose(np.angle(ZC), -np.pi / 2)
    npt.assert_allclose(np.abs(ZL), 2 * np.pi * f * 1e-3)               # grows with f
    npt.assert_allclose(np.abs(ZC), 1 / (2 * np.pi * f * 1e-6))         # falls with f
    assert impl.series(10.0, 20.0) == pytest.approx(30.0)
    assert impl.parallel(10.0, 10.0) == pytest.approx(5.0)
    assert impl.parallel(10.0, 30.0, 15.0) == pytest.approx(1 / (1 / 10 + 1 / 30 + 1 / 15))


def test_series_rlc_is_resistive_at_resonance_with_bandwidth_f0_over_Q():
    R, L, C = 10.0, 10e-3, 1e-6
    f0 = impl.resonant_frequency(L, C)
    assert f0 == pytest.approx(1 / (2 * np.pi * np.sqrt(L * C)))
    I0 = impl.rlc_series_current(1.0, R, L, C, f0)
    assert I0 == pytest.approx(1.0 / R, rel=1e-12)                       # purely resistive, in phase
    Q = impl.quality_factor(R, L, C)
    assert Q == pytest.approx(2 * np.pi * f0 * L / R)
    f = np.linspace(0.5 * f0, 1.5 * f0, 200001)
    I = np.abs(impl.rlc_series_current(1.0, R, L, C, f))
    band = f[I >= I.max() / np.sqrt(2)]
    assert band[-1] - band[0] == pytest.approx(f0 / Q, rel=2e-3)
    # below resonance the circuit is capacitive (current leads), above inductive (current lags)
    assert np.angle(impl.rlc_series_current(1.0, R, L, C, 0.5 * f0)) > 0
    assert np.angle(impl.rlc_series_current(1.0, R, L, C, 2.0 * f0)) < 0


def test_rms_and_ac_power():
    t = np.linspace(0.0, 1.0, 40001)[:-1]
    assert impl.rms(3.0 * np.sin(2 * np.pi * 5 * t)) == pytest.approx(3.0 / np.sqrt(2), rel=1e-9)
    assert impl.rms(np.full(10, -2.0)) == pytest.approx(2.0)
    P, Q, pf = impl.ac_power(115.0, 10.0 + 0j)
    assert P == pytest.approx(115.0 ** 2 / 10.0) and Q == pytest.approx(0.0, abs=1e-9) and pf == pytest.approx(1.0)
    Z = complex(10.0, 5.0)
    P, Q, pf = impl.ac_power(115.0, Z)
    I = 115.0 / abs(Z)
    assert P == pytest.approx(I ** 2 * Z.real)                            # only the resistance dissipates
    assert Q == pytest.approx(I ** 2 * Z.imag)                            # inductive load absorbs vars
    assert pf == pytest.approx(np.cos(np.angle(Z)))
    _, Qc, _ = impl.ac_power(115.0, complex(10.0, -5.0))
    assert Qc < 0


# ---------------------------------------------------------------- the RC filter

def test_rc_corner_frequency_is_minus_3dB_and_minus_45_degrees():
    R, C = 3.3e3, 1e-6
    fc = impl.rc_corner_frequency(R, C)
    assert fc == pytest.approx(1 / (2 * np.pi * R * C))
    H = impl.rc_lowpass(fc, R, C)
    assert abs(H) == pytest.approx(1 / np.sqrt(2), rel=1e-12)
    assert np.degrees(np.angle(H)) == pytest.approx(-45.0)
    assert impl.gain_db(H) == pytest.approx(-3.0103, abs=1e-4)
    assert abs(impl.rc_lowpass(0.0, R, C)) == pytest.approx(1.0)
    # high-pass is the complement: H_lp + H_hp = 1 at every frequency
    f = np.logspace(-1, 4, 50)
    npt.assert_allclose(impl.rc_lowpass(f, R, C) + impl.rc_highpass(f, R, C), 1.0, atol=1e-12)


def test_rc_rolls_off_at_20_dB_per_decade_per_stage():
    R, C = 1e3, 1e-6
    fc = impl.rc_corner_frequency(R, C)
    for n in (1, 2, 3):
        g1 = impl.gain_db(impl.rc_lowpass(1e3 * fc, R, C, n))
        g2 = impl.gain_db(impl.rc_lowpass(1e4 * fc, R, C, n))
        assert g1 - g2 == pytest.approx(20.0 * n, rel=1e-5)


def test_buffered_cascade_matches_scipy_signal_freqs():
    R, C = 3.3e3, 1e-6
    f = np.logspace(0, 4, 60)
    for n in (1, 2, 3):
        lti = impl.rc_lti(R, C, n)
        _, H = signal.freqs(lti.num, lti.den, worN=2 * np.pi * f)
        npt.assert_allclose(impl.rc_lowpass(f, R, C, n), H, rtol=1e-10)
    # and an oracle built here from scratch: H(s) = 1/(RCs + 1)^2
    _, H2 = signal.freqs([1.0], np.polymul([R * C, 1.0], [R * C, 1.0]), worN=2 * np.pi * f)
    npt.assert_allclose(impl.rc_lowpass(f, R, C, 2), H2, rtol=1e-10)


def test_step_response_matches_closed_form_and_scipy_step():
    R, C = 10e3, 1e-6
    t = np.linspace(0.0, 5 * R * C, 501)
    v = impl.rc_step_response(t, R, C, V=3.3)
    assert v[0] == pytest.approx(0.0, abs=1e-15)
    assert impl.rc_step_response(R * C, R, C, 3.3) == pytest.approx(3.3 * (1 - np.exp(-1)))
    _, y = signal.step(([1.0], [R * C, 1.0]), T=t)
    npt.assert_allclose(v, 3.3 * y, atol=1e-9)
    # group delay at low frequency: -dphi/domega = n RC
    w = 2 * np.pi * np.array([0.01, 0.02])
    for n in (1, 2):
        ph = np.unwrap(np.angle(impl.rc_lowpass(w / (2 * np.pi), R, C, n)))
        assert -(ph[1] - ph[0]) / (w[1] - w[0]) == pytest.approx(impl.group_delay_dc(R, C, n), rel=1e-6)


def test_unbuffered_ladder_matches_nodal_analysis_and_differs_from_buffered():
    R, C = 3.3e3, 1e-6
    f = np.array([5.0, 48.0, 430.0])
    H = impl.rc2_unbuffered(f, R, C)
    for k, fk in enumerate(f):
        Y = 1j * 2 * np.pi * fk * C
        # unknowns v1, v2; source 1 V through R into node 1
        A = np.array([[1 / R + Y + 1 / R, -1 / R], [-1 / R, 1 / R + Y]])
        b = np.array([1 / R, 0.0])
        v = np.linalg.solve(A, b)
        assert H[k] == pytest.approx(v[1], rel=1e-12)
    assert abs(H[1]) < abs(impl.rc_lowpass(f[1], R, C, 2))                 # loading pulls the corner down
    # loaded single stage: DC gain is the divider R_L / (R + R_L)
    assert abs(impl.rc_loaded(0.0, 1e3, 1e-6, 9e3)) == pytest.approx(0.9)
    assert abs(impl.rc_loaded(100.0, 1e3, 1e-6, 1e15)) == pytest.approx(abs(impl.rc_lowpass(100.0, 1e3, 1e-6)), rel=1e-9)


def test_required_stages_meets_attenuation_and_droop():
    n, fc = impl.required_stages_and_corner(5.0, 430.0, 62.0, 0.5)
    assert n == 2
    H = impl.rc_lowpass(np.array([5.0, 430.0]), 1.0, 1.0 / (2 * np.pi * fc), n)
    assert abs(H[1]) == pytest.approx(1 / 62.0, rel=1e-9)                # attenuation met exactly at the largest fc
    assert -impl.gain_db(H[0]) <= 0.5
    n1, _ = impl.required_stages_and_corner(5.0, 430.0, 62.0, 3.0)        # a lax droop spec admits one stage
    assert n1 == 1
    with pytest.raises(ValueError):
        impl.required_stages_and_corner(40.0, 50.0, 1e4, 0.1)


# ---------------------------------------------------------------- op-amp, diode, transistor

def test_op_amp_gains_and_the_infinite_gain_limit():
    assert impl.noninverting_gain(27e3, 10e3) == pytest.approx(3.7)
    assert impl.noninverting_gain(0.0, 10e3) == pytest.approx(1.0)        # buffer
    assert impl.inverting_gain(47e3, 4.7e3) == pytest.approx(-10.0)
    beta = 10e3 / 37e3
    gains = [impl.closed_loop_gain(A, beta) for A in (1e2, 1e4, 1e6, 1e9)]
    errs = np.abs(np.array(gains) - 1 / beta)
    assert np.all(np.diff(errs) < 0) and errs[-1] < 1e-7
    # the error is about G^2/A: halving with doubling A
    A = 1e5
    assert 1 / beta - impl.closed_loop_gain(A, beta) == pytest.approx((1 / beta) ** 2 / A, rel=1e-3)


def test_diode_is_exponential_forward_and_saturates_reverse():
    VT = impl.thermal_voltage(300.0)
    assert VT == pytest.approx(0.02585, rel=1e-3)
    I1, I2 = impl.diode_current(0.6), impl.diode_current(0.6 + VT * np.log(10))
    assert I2 / I1 == pytest.approx(10.0, rel=1e-6)                        # 60 mV per decade at n = 1
    assert impl.diode_current(-1.0, I_s=1e-12) == pytest.approx(-1e-12, rel=1e-9)
    assert impl.diode_current(0.0) == 0.0
    Rb = impl.bjt_base_resistor(3.3, 0.12, 100.0)
    assert (3.3 - 0.7) / Rb == pytest.approx(5 * 0.12 / 100)


# ---------------------------------------------------------------- logic and the ADC

def test_logic_gates_truth_tables_and_half_adder():
    a = np.array([0, 0, 1, 1]); b = np.array([0, 1, 0, 1])
    npt.assert_array_equal(impl.gate("AND", a, b), [0, 0, 0, 1])
    npt.assert_array_equal(impl.gate("OR", a, b), [0, 1, 1, 1])
    npt.assert_array_equal(impl.gate("XOR", a, b), [0, 1, 1, 0])
    npt.assert_array_equal(impl.gate("NAND", a, b), [1, 1, 1, 0])
    npt.assert_array_equal(impl.gate("NOR", a, b), [1, 0, 0, 0])
    npt.assert_array_equal(impl.gate("NOT", a), [1, 1, 0, 0])
    # NAND is universal: NOT a = a NAND a; De Morgan
    npt.assert_array_equal(impl.gate("NAND", a, a), impl.gate("NOT", a))
    npt.assert_array_equal(impl.gate("NOR", a, b), impl.gate("AND", impl.gate("NOT", a), impl.gate("NOT", b)))
    for x in (0, 1):
        for y in (0, 1):
            s, c = impl.half_adder(x, y)
            assert 2 * c + s == x + y
    with pytest.raises(ValueError):
        impl.gate("XNAND", 1, 0)


def test_aliasing_folds_into_the_first_nyquist_zone():
    fs = 100.0
    assert impl.aliased_frequency(430.0, fs) == pytest.approx(30.0)
    assert impl.aliased_frequency(20.0, fs) == pytest.approx(20.0)
    assert impl.aliased_frequency(70.0, fs) == pytest.approx(30.0)
    assert impl.aliased_frequency(400.0, fs) == pytest.approx(0.0)
    f = np.linspace(0, 1000, 997)
    fa = impl.aliased_frequency(f, fs)
    assert np.all(fa >= 0) and np.all(fa <= fs / 2 + 1e-9)
    # the samples of a 430 Hz and a 30 Hz cosine are identical
    n = np.arange(50)
    npt.assert_allclose(np.cos(2 * np.pi * 430 * n / fs), np.cos(2 * np.pi * 30 * n / fs), atol=1e-9)


def test_adc_quantisation_lsb_clipping_and_error():
    assert impl.adc_lsb(3.3, 12) == pytest.approx(3.3 / 4096)
    npt.assert_array_equal(impl.adc_convert([-0.1, 0.0, 3.3, 5.0], 3.3, 12), [0, 0, 4095, 4095])
    rng = np.random.default_rng(1012)
    v = rng.uniform(0.01, 3.29, 20000)
    err = impl.counts_to_volts(impl.adc_convert(v, 3.3, 12), 3.3, 12) - v
    lsb = impl.adc_lsb(3.3, 12)
    assert np.max(np.abs(err)) <= 0.5 * lsb * (1 + 1e-9)
    # uniform error over one LSB: rms = LSB/sqrt(12); 20000 samples give about 0.5 % standard error
    assert np.sqrt(np.mean(err ** 2)) == pytest.approx(impl.quantisation_noise_rms(3.3, 12), rel=0.03)
    assert impl.adc_lsb(3.3, 16) == pytest.approx(impl.adc_lsb(3.3, 12) / 16)


def test_filter_reduces_the_alias_tone_in_synthetic_sampled_data():
    fs, f_vib, a_vib = 100.0, 430.0, 0.05
    R, C = 3.3e3, 1e-6
    t = np.arange(0.0, 3.0, 1.0 / 20000.0)
    v = impl.synthetic_sensor_signal(t, f_sig=2.0, a_sig=0.5, f_vib=f_vib, a_vib=a_vib)
    lsb = impl.adc_lsb(3.3, 12)
    amps = []
    for n in (0, 1, 2):
        ts, vs, counts = impl.filter_and_sample(t, v, R, C, n, fs)
        assert ts.size == 300 and counts.dtype.kind == "i"
        keep = ts >= 1.0
        amp = impl.tone_amplitude(counts[keep].astype(float), fs, 30.0)
        pred = a_vib * abs(impl.rc_lowpass(f_vib, R, C, n)) / lsb if n else a_vib / lsb
        assert amp == pytest.approx(pred, abs=0.6)                         # within about half an LSB of quantisation
        amps.append(amp)
    assert amps[0] > 50 and amps[2] < 1.0                                  # 62 LSB of alias without, under 1 LSB with two stages
    # the 2 Hz signal passes almost unchanged through the two stages
    ts, vs, _ = impl.filter_and_sample(t, v, R, C, 2, fs)
    keep = ts >= 1.0
    assert impl.tone_amplitude(vs[keep], fs, 2.0) == pytest.approx(0.5 * abs(impl.rc_lowpass(2.0, R, C, 2)), rel=1e-3)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "rc_filter.png").exists()
