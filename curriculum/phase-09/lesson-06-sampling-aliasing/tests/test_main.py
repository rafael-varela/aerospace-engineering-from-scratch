"""Public tests for P09.06 Sampling & Aliasing. They assert properties, not magic numbers: aliases fold
into [0, fs/2] and the samples of a tone and of its alias are identical; the FFT of a sampled tone peaks
at the predicted alias (independent check 1); the Butterworth prototype built from its poles equals
``scipy.signal.butter(analog=True)`` (independent check 2) and has |H|^2 = 1 / (1 + (w/wc)^2n); the order
rule meets both specifications and one order fewer fails; the group delay formula agrees with the slope
of the phase; the from-scratch analog simulator agrees with ``scipy.signal.lsim``, is linear and has the
right DC gain; sinc reconstruction interpolates the samples and recovers a band-limited signal; and in a
SYNTHETIC gyro record the alias measured after sampling equals a_vib |H(f_vib)|.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy import signal

from conftest import load_impl

impl = load_impl(__file__)


# ---------------------------------------------------------------- aliasing

def test_aliased_frequency_folds_into_the_first_nyquist_zone():
    fs = 400.0
    assert impl.aliased_frequency(390.0, fs) == pytest.approx(10.0)
    assert impl.aliased_frequency(1190.0, fs) == pytest.approx(10.0)
    assert impl.aliased_frequency(150.0, fs) == pytest.approx(150.0)        # below Nyquist: unchanged
    assert impl.aliased_frequency(800.0, fs) == pytest.approx(0.0, abs=1e-9)  # a multiple of fs looks like DC
    f = np.linspace(0.0, 5000.0, 2001)
    fa = impl.aliased_frequency(f, fs)
    assert fa.shape == f.shape
    assert np.all(fa >= 0.0) and np.all(fa <= fs / 2 + 1e-9)


@pytest.mark.parametrize("f", [130.0, 390.0, 610.0, 1190.0])
def test_samples_of_a_tone_and_of_its_alias_are_identical(f):
    fs = 400.0
    n = np.arange(64)
    fa = impl.aliased_frequency(f, fs)
    npt.assert_allclose(np.cos(2 * np.pi * f * n / fs), np.cos(2 * np.pi * fa * n / fs), atol=1e-9)


@pytest.mark.parametrize("f", [37.0, 390.0, 557.0, 1210.0])
def test_fft_of_a_sampled_tone_peaks_at_the_predicted_alias(f):
    # independent check: the spectrum of the samples, not the folding formula, says where the tone went
    fs, N = 400.0, 400
    n = np.arange(N)
    x = np.sin(2 * np.pi * f * n / fs + 0.4)
    assert impl.dominant_frequency(x, fs) == pytest.approx(impl.aliased_frequency(f, fs), abs=fs / N / 4)


# ---------------------------------------------------------------- Butterworth prototype

@pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6])
def test_butterworth_matches_scipy_butter_analog(n):
    wc = 2 * np.pi * 50.0
    b, a = impl.butterworth_analog(n, wc)
    bs, as_ = signal.butter(n, wc, analog=True)
    npt.assert_allclose(a / a[0], as_ / as_[0], rtol=1e-10)
    w = np.logspace(1, 4, 50)
    _, Hs = signal.freqs(bs, as_, worN=w)
    npt.assert_allclose(impl.freq_response(b, a, w), Hs, rtol=1e-9, atol=1e-14)


def test_butterworth_poles_lie_on_a_circle_in_the_left_half_plane():
    for n in (1, 2, 5, 8):
        p = impl.butterworth_poles(n, 3.0)
        assert p.size == n
        npt.assert_allclose(np.abs(p), 3.0, rtol=1e-12)
        assert np.all(p.real < 0)
        dist = np.abs(p[:, None] - np.conj(p)[None, :]).min(axis=1)
        assert np.all(dist < 1e-12)                                       # conjugate pairs: a(s) is real


def test_butterworth_magnitude_is_maximally_flat_and_minus_3db_at_the_corner():
    wc = 2 * np.pi * 40.0
    w = np.logspace(0, 4, 200)
    for n in (1, 2, 3, 4):
        b, a = impl.butterworth_analog(n, wc)
        H = impl.freq_response(b, a, w)
        npt.assert_allclose(np.abs(H) ** 2, 1.0 / (1.0 + (w / wc) ** (2 * n)), rtol=1e-9)
        assert 20 * np.log10(abs(impl.freq_response(b, a, wc))) == pytest.approx(-3.0103, abs=1e-3)
        assert abs(impl.freq_response(b, a, 0.0)) == pytest.approx(1.0)
        # roll-off of 20 n dB per decade far above the corner
        drop = 20 * np.log10(abs(impl.freq_response(b, a, 1e3 * wc)) / abs(impl.freq_response(b, a, 1e4 * wc)))
        assert drop == pytest.approx(20.0 * n, abs=1e-3)


def test_freq_response_matches_scipy_freqs_for_a_general_rational_function():
    b, a = [2.0, 3.0, 5.0], [1.0, 4.0, 9.0, 7.0]
    w = np.logspace(-2, 2, 40)
    _, Hs = signal.freqs(b, a, worN=w)
    npt.assert_allclose(impl.freq_response(np.array(b), np.array(a), w), Hs, rtol=1e-12)


@pytest.mark.parametrize("spec", [(30.0, 0.5, 370.0, 47.96), (20.0, 1.0, 180.0, 40.0), (5.0, 0.1, 25.0, 60.0)])
def test_order_rule_meets_both_specifications_and_one_order_fewer_cannot(spec):
    fp, droop, fst, att = spec
    n, lo, hi = impl.butterworth_order(fp, droop, fst, att)
    assert lo <= hi
    for fc in (lo, np.sqrt(lo * hi), hi):
        b, a = impl.butterworth_analog(n, 2 * np.pi * fc)
        assert -20 * np.log10(abs(impl.freq_response(b, a, 2 * np.pi * fp))) <= droop + 1e-9
        assert -20 * np.log10(abs(impl.freq_response(b, a, 2 * np.pi * fst))) >= att - 1e-9
    if n > 1:   # with n - 1 the passband needs a higher corner than the stopband allows
        ep, es = 10 ** (droop / 10) - 1, 10 ** (att / 10) - 1
        assert fp / ep ** (1 / (2 * (n - 1))) > fst / es ** (1 / (2 * (n - 1)))


def test_antialias_spec_puts_the_stop_edge_where_aliases_enter_the_band():
    d = impl.antialias_spec(30.0, 400.0, 5.0, 0.02, 0.5)
    assert d["f_stop"] == pytest.approx(370.0)
    assert d["atten_db"] == pytest.approx(20 * np.log10(250.0))
    assert d["fc_min"] <= d["fc"] <= d["fc_max"]
    # anything above f_stop is attenuated at least as much (Butterworth is monotonic), so the alias of
    # any line there stays below a_allow
    b, a = impl.butterworth_analog(d["n"], 2 * np.pi * d["fc"])
    f = np.linspace(370.0, 5000.0, 500)
    assert np.all(5.0 * np.abs(impl.freq_response(b, a, 2 * np.pi * f)) < 0.02)
    with pytest.raises(ValueError):
        impl.antialias_spec(30.0, 50.0, 5.0, 0.02, 0.5)        # fs < 2 f_band: no filter can help


def test_group_delay_at_dc_equals_the_slope_of_the_phase():
    for b, a in [impl.butterworth_analog(3, 2 * np.pi * 50.0), (np.array([1.0, 2.0]), np.array([1.0, 5.0, 6.0]))]:
        w = np.array([-1e-3, 1e-3])
        ph = np.unwrap(np.angle(impl.freq_response(b, a, w)))
        assert impl.group_delay_dc(b, a) == pytest.approx(-(ph[1] - ph[0]) / (w[1] - w[0]), rel=1e-5)
    # first-order low-pass 1/(tau s + 1): delay tau
    assert impl.group_delay_dc(np.array([1.0]), np.array([0.01, 1.0])) == pytest.approx(0.01)


# ---------------------------------------------------------------- time domain

def test_simulate_analog_matches_scipy_lsim_with_held_input():
    b, a = impl.butterworth_analog(3, 2 * np.pi * 50.0)
    t = np.arange(4000) / 20000.0
    rng = np.random.default_rng(3)
    u = np.sin(2 * np.pi * 37.0 * t) + 0.3 * rng.standard_normal(t.size)
    y = impl.simulate_analog(b, a, t, u)
    _, ys, _ = signal.lsim((b, a), u, t, interp=False)
    npt.assert_allclose(y, ys, atol=1e-9)


def test_simulate_analog_dc_gain_zero_input_and_superposition():
    b, a = impl.butterworth_analog(2, 2 * np.pi * 20.0)
    t = np.arange(3000) / 10000.0
    step = impl.simulate_analog(b, a, t, np.ones_like(t))
    assert step[-1] == pytest.approx(1.0, abs=1e-6)                         # unit DC gain
    npt.assert_allclose(impl.simulate_analog(b, a, t, np.zeros_like(t)), 0.0, atol=0.0)
    u1, u2 = np.sin(2 * np.pi * 5 * t), np.cos(2 * np.pi * 70 * t)
    npt.assert_allclose(impl.simulate_analog(b, a, t, 2 * u1 - 3 * u2),
                        2 * impl.simulate_analog(b, a, t, u1) - 3 * impl.simulate_analog(b, a, t, u2), atol=1e-12)


def test_sinc_reconstruction_interpolates_and_recovers_a_band_limited_signal():
    fs = 50.0
    n = np.arange(200)
    x = np.sin(2 * np.pi * 3.0 * n / fs) + 0.5 * np.cos(2 * np.pi * 11.0 * n / fs)
    npt.assert_allclose(impl.sinc_reconstruct(x, fs, n / fs), x, atol=1e-12)   # passes through every sample
    t = np.linspace(1.0, 3.0, 401)                                            # away from the record's ends
    xt = np.sin(2 * np.pi * 3.0 * t) + 0.5 * np.cos(2 * np.pi * 11.0 * t)
    assert np.max(np.abs(impl.sinc_reconstruct(x, fs, t) - xt)) < 0.01        # truncation error of a 4 s record


def test_reconstruction_of_an_undersampled_tone_is_its_alias():
    fs = 50.0
    n = np.arange(200)
    t = np.linspace(1.0, 3.0, 401)
    xr = impl.sinc_reconstruct(np.cos(2 * np.pi * 40.0 * n / fs), fs, t)
    assert np.max(np.abs(xr - np.cos(2 * np.pi * 10.0 * t))) < 0.01            # 40 Hz comes back as 10 Hz


# ---------------------------------------------------------------- the SYNTHETIC IMU experiment

def test_synthetic_gyro_alias_after_sampling_equals_the_predicted_attenuation():
    rows = impl.imu_experiment(orders=(0, 1, 3), duration=1.25)
    for r in rows:
        # tolerance: the seeded white noise puts about 0.005 deg/s into one DFT bin without a filter and
        # far less with one; 3 % covers the zero-order-hold amplitude error of the fine simulation grid
        assert r["measured"] == pytest.approx(r["predicted"], rel=0.03, abs=0.006)
    assert rows[0]["measured"] > 4.5                       # unprotected: the 5 deg/s line is all there, at 10 Hz
    assert rows[-1]["measured"] < 0.02                     # the designed order meets the 0.02 deg/s requirement
    assert rows[-1]["delay_ms"] > rows[1]["delay_ms"] > 0  # and pays for it in delay


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "sampling_aliasing.png").exists()
