"""Public tests for P09.07 Digital Filters. They assert properties, not magic numbers: the frequency
response from the coefficients equals ``scipy.signal.freqz``; the difference equation is linear and
time-invariant, returns the taps as the impulse response of an FIR and a^n for a first-order IIR, and
the DFT of its impulse response equals the frequency response (independent check); the moving average
has its zeros at multiples of fs/N; the windowed-sinc FIR is symmetric with a constant group delay of
(N - 1)/2 samples and the stopband its window promises; the bilinear transform equals
``scipy.signal.cont2discrete(method='bilinear')``, maps the left half plane inside the unit circle and,
prewarped, matches the analog response exactly at the chosen frequency; the digital Butterworth equals
``scipy.signal.butter(fs=...)``; and on a SYNTHETIC gyro record the filtered noise obeys Parseval.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy import signal

from conftest import load_impl

impl = load_impl(__file__)

FS = 1000.0


def some_filters():
    return [
        (np.array([0.2, 0.3, 0.5]), np.array([1.0])),                                  # FIR
        (np.array([0.1, 0.2]), np.array([1.0, -0.6])),                                 # first-order IIR
        signal.butter(3, 60.0, fs=FS),                                                 # oracle-designed IIR
        (np.array([1.0, -1.2, 0.5]), np.array([2.0, -1.0, 0.3])),                      # a[0] != 1
    ]


# ---------------------------------------------------------------- evaluation and the difference equation

@pytest.mark.parametrize("k", range(4))
def test_digital_freq_response_matches_scipy_freqz(k):
    b, a = some_filters()[k]
    f = np.linspace(0.0, FS / 2, 77)
    _, H = signal.freqz(b, a, worN=f, fs=FS)
    npt.assert_allclose(impl.digital_freq_response(b, a, f, FS), H, rtol=1e-10, atol=1e-13)


def test_fir_impulse_response_is_the_taps_and_first_order_iir_is_geometric():
    b = np.array([0.2, -0.5, 0.7, 0.1])
    npt.assert_allclose(impl.impulse_response(b, [1.0], 8), np.concatenate([b, np.zeros(4)]), atol=0.0)
    h = impl.impulse_response([1.0], [1.0, -0.8], 30)                     # y[n] = 0.8 y[n-1] + x[n]
    npt.assert_allclose(h, 0.8 ** np.arange(30), rtol=1e-12)
    h2 = impl.impulse_response([3.0], [2.0, -1.0], 10)                    # a[0] = 2: y = 1.5 (0.5)^n
    npt.assert_allclose(h2, 1.5 * 0.5 ** np.arange(10), rtol=1e-12)


@pytest.mark.parametrize("k", range(4))
def test_dft_of_the_impulse_response_equals_the_frequency_response(k):
    # independent check: the difference equation's own impulse response, transformed, gives freqz
    b, a = some_filters()[k]
    N = 4096                                                               # long enough for the IIR tails to die
    h = impl.impulse_response(b, a, N)
    H_dft = np.fft.rfft(h)
    f = np.fft.rfftfreq(N, 1 / FS)
    _, H = signal.freqz(b, a, worN=f, fs=FS)
    npt.assert_allclose(H_dft, H, atol=1e-9)


def test_difference_equation_is_linear_and_time_invariant():
    b, a = signal.butter(2, 40.0, fs=FS)
    rng = np.random.default_rng(1)
    x1, x2 = rng.standard_normal(300), rng.standard_normal(300)
    y1, y2 = impl.difference_equation(b, a, x1), impl.difference_equation(b, a, x2)
    npt.assert_allclose(impl.difference_equation(b, a, 2 * x1 - 0.5 * x2), 2 * y1 - 0.5 * y2, atol=1e-12)
    shifted = impl.difference_equation(b, a, np.concatenate([np.zeros(25), x1]))
    npt.assert_allclose(shifted[25:], y1, atol=1e-12)
    npt.assert_allclose(shifted[:25], 0.0, atol=0.0)                       # causal: nothing before the input


def test_group_delay_formula_equals_the_slope_of_the_unwrapped_phase():
    b, a = signal.butter(4, 40.0, fs=FS)
    f = np.linspace(1.0, 100.0, 50)
    df = 1e-4
    H = lambda ff: impl.digital_freq_response(b, a, ff, FS)
    dphase = np.angle(H(f + df) / H(f - df))                                # phase change over 2 df, no wrapping
    slope = dphase / (2 * df) / (2 * np.pi / FS)                            # d phase / d omega
    npt.assert_allclose(impl.group_delay(b, a, f, FS), -slope, rtol=1e-5)


def test_stability_from_the_poles():
    assert impl.is_stable_digital(np.array([1.0, -0.9]))
    assert not impl.is_stable_digital(np.array([1.0, -1.1]))
    assert not impl.is_stable_digital(np.array([1.0, -2.0, 1.0]))          # double pole on the unit circle
    assert impl.is_stable_digital(np.array([1.0]))                         # FIR: no poles


# ---------------------------------------------------------------- FIR design

def test_moving_average_has_unit_dc_gain_and_zeros_at_multiples_of_fs_over_n():
    b = impl.moving_average(50)
    npt.assert_allclose(b, 1 / 50)
    assert abs(impl.digital_freq_response(b, [1.0], 0.0, FS)) == pytest.approx(1.0)
    zeros = np.arange(1, 25) * FS / 50
    assert np.max(np.abs(impl.digital_freq_response(b, [1.0], zeros, FS))) < 1e-12
    npt.assert_allclose(impl.group_delay(b, [1.0], np.array([2.0, 7.0]), FS), 24.5, rtol=1e-9)


@pytest.mark.parametrize("name,ref", [("hann", np.hanning), ("hamming", np.hamming), ("blackman", np.blackman)])
def test_windows_match_their_definitions(name, ref):
    for n in (11, 64, 101):
        npt.assert_allclose(impl.window(name, n), ref(n), atol=1e-12)
    npt.assert_allclose(impl.window("rectangular", 9), 1.0)


@pytest.mark.parametrize("win,stop_db", [("hamming", -50.0), ("blackman", -70.0), ("hann", -40.0)])
def test_windowed_sinc_is_linear_phase_with_the_promised_stopband(win, stop_db):
    N, fc = 101, 40.0
    h = impl.windowed_sinc_lowpass(N, fc, FS, win)
    npt.assert_allclose(h, h[::-1], atol=1e-15)                              # symmetric taps
    assert h.sum() == pytest.approx(1.0)                                     # unit DC gain
    assert abs(impl.digital_freq_response(h, [1.0], fc, FS)) == pytest.approx(0.5, abs=0.02)   # -6 dB at fc
    f_pass = np.linspace(0.5, 15.0, 20)
    npt.assert_allclose(impl.group_delay(h, [1.0], f_pass, FS), (N - 1) / 2, rtol=1e-8)
    # beyond the transition band (about 6 fs/N for these windows) the gain stays below the window's floor
    f_stop = np.linspace(fc + 6 * FS / N, FS / 2, 2000)
    assert 20 * np.log10(np.max(np.abs(impl.digital_freq_response(h, [1.0], f_stop, FS)))) < stop_db


# ---------------------------------------------------------------- IIR design

@pytest.mark.parametrize("ba", [([1.0], [0.01, 1.0]), ([100.0, 0.0], [1.0, 20.0, 400.0]), ([2.0, 5.0], [1.0, 3.0, 7.0, 9.0])])
def test_bilinear_matches_scipy_cont2discrete(ba):
    b, a = (np.array(c, dtype=float) for c in ba)
    bz, az = impl.bilinear(b, a, FS)
    num, den, _ = signal.cont2discrete((b, a), 1 / FS, method="bilinear")
    num = np.ravel(num)
    npt.assert_allclose(az, den / den[0], rtol=1e-9, atol=1e-12)
    npt.assert_allclose(bz, num / den[0], rtol=1e-9, atol=1e-12)


def test_bilinear_maps_stable_poles_inside_and_unstable_outside_the_unit_circle():
    _, az = impl.bilinear(np.array([1.0]), np.array([1.0, 2.0, 50.0]), FS)
    assert impl.is_stable_digital(az)
    _, az_u = impl.bilinear(np.array([1.0]), np.array([1.0, -2.0, 50.0]), FS)
    assert not impl.is_stable_digital(az_u)


@pytest.mark.parametrize("fp", [40.0, 250.0])
def test_prewarping_matches_the_analog_response_exactly_at_the_chosen_frequency(fp):
    b, a = impl.butterworth_analog(2, 2 * np.pi * fp)
    bz, az = impl.bilinear(b, a, FS, f_prewarp=fp)
    assert impl.digital_freq_response(bz, az, fp, FS) == pytest.approx(complex(impl.analog_freq_response(b, a, fp)), abs=1e-12)
    bz0, az0 = impl.bilinear(b, a, FS)                                       # plain Tustin warps the corner down
    assert abs(impl.digital_freq_response(bz0, az0, fp, FS)) < 1 / np.sqrt(2)


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_digital_butterworth_equals_scipy_butter(n):
    bz, az = impl.butterworth_digital(n, 40.0, FS)
    bs, as_ = signal.butter(n, 40.0, fs=FS)
    npt.assert_allclose(az, as_, rtol=1e-9, atol=1e-12)
    npt.assert_allclose(bz, bs, rtol=1e-8, atol=1e-14)
    assert abs(impl.digital_freq_response(bz, az, FS / 2, FS)) < 1e-12        # z = -1 is an n-fold zero


# ---------------------------------------------------------------- the SYNTHETIC gyro comparison

def test_synthetic_gyro_comparison_obeys_parseval_and_the_delay_trade():
    rows = {r["name"]: r for r in impl.gyro_filter_comparison()}
    for r in rows.values():
        # filtered white noise: std = sigma sqrt(sum h^2); 3500 correlated samples, so allow 10 %
        assert r["noise_out"] == pytest.approx(r["noise_pred"], rel=0.10)
        assert r["noise_out"] < 0.3                                         # every filter reduces the noise
    assert rows["MA 50"]["residual_vib"] < 1e-6                            # 180 Hz sits on a zero of the MA
    assert rows["FIR 101 Hamming 40 Hz"]["delay_2hz_ms"] == pytest.approx(50.0)
    iir2, iir4 = rows["IIR Butterworth 2, 40 Hz"], rows["IIR Butterworth 4, 40 Hz"]
    assert iir2["delay_2hz_ms"] < iir4["delay_2hz_ms"] < rows["FIR 101 Hamming 40 Hz"]["delay_2hz_ms"]
    assert iir4["residual_vib"] < iir2["residual_vib"]
    # the residual vibration is the attenuation times the 3 deg/s line
    assert iir2["residual_vib"] == pytest.approx(3.0 * 10 ** (iir2["atten_180_db"] / 20), rel=1e-3)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "digital_filters.png").exists()
