"""Public tests for P02.10 Fourier Transform & the FFT. They assert properties, with scipy.fft as the oracle.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-02/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy import fft as sfft

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(210)


# ---------------------------------------------------------------- the DFT and its matrix

def test_dft_matrix_is_symmetric_and_scaled_unitary():
    for n in (1, 4, 7, 16):
        F = impl.dft_matrix(n)
        npt.assert_allclose(F, F.T, atol=1e-12)
        npt.assert_allclose(np.conj(F).T @ F, n * np.eye(n), atol=1e-10)    # F^H F = n I
        assert F[0, 0] == pytest.approx(1.0)
    with pytest.raises(ValueError):
        impl.dft_matrix(0)


@pytest.mark.parametrize("n", [1, 5, 16, 63, 128])
def test_dft_matches_scipy_for_any_length(n):
    x = RNG.standard_normal(n) + 1j * RNG.standard_normal(n)
    npt.assert_allclose(impl.dft(x), sfft.fft(x), rtol=1e-10, atol=1e-10)
    npt.assert_allclose(impl.idft(impl.dft(x)), x, atol=1e-10)             # round trip
    with pytest.raises(ValueError):
        impl.dft(np.ones((2, 2)))


def test_dft_is_linear_and_obeys_the_shift_theorem():
    n = 32
    x, y = RNG.standard_normal(n), RNG.standard_normal(n)
    a, b = 1.7, -0.4
    npt.assert_allclose(impl.dft(a * x + b * y), a * impl.dft(x) + b * impl.dft(y), atol=1e-10)
    # A circular shift by m multiplies bin k by exp(-2 pi i k m / N).
    m = 5
    k = np.arange(n)
    npt.assert_allclose(impl.dft(np.roll(x, m)), impl.dft(x) * np.exp(-2j * np.pi * k * m / n), atol=1e-10)
    # A real sequence has a conjugate-symmetric spectrum.
    X = impl.dft(x)
    npt.assert_allclose(X[1:], np.conj(X[1:][::-1]), atol=1e-10)


# ---------------------------------------------------------------- the FFT

@pytest.mark.parametrize("n", [1, 2, 8, 64, 1024])
def test_fft_radix2_matches_the_dft_and_scipy(n):
    x = RNG.standard_normal(n)
    X = impl.fft_radix2(x)
    npt.assert_allclose(X, impl.dft(x), rtol=1e-9, atol=1e-9)
    npt.assert_allclose(X, sfft.fft(x), rtol=1e-9, atol=1e-9)


def test_fft_radix2_rejects_lengths_that_are_not_powers_of_two():
    for n in (0, 3, 6, 100):
        with pytest.raises(ValueError):
            impl.fft_radix2(np.ones(max(n, 0)))
    # fft_any falls back to the DFT for those lengths.
    x = RNG.standard_normal(100)
    npt.assert_allclose(impl.fft_any(x), sfft.fft(x), atol=1e-9)
    npt.assert_allclose(impl.fft_any(x[:64]), sfft.fft(x[:64]), atol=1e-9)


def test_parseval_holds_for_dft_and_fft():
    x = RNG.standard_normal(256)
    e_t, e_f = impl.parseval(x, impl.fft_radix2(x))
    assert e_t == pytest.approx(e_f, rel=1e-12)
    e_t, e_f = impl.parseval(x[:100], impl.dft(x[:100]))
    assert e_t == pytest.approx(e_f, rel=1e-12)
    assert e_t == pytest.approx(np.sum(x[:100] ** 2))


def test_fft_frequencies_match_scipy_fftfreq():
    for n, dt in ((8, 1.0), (9, 0.01), (1024, 1 / 1024)):
        npt.assert_allclose(impl.fft_frequencies(n, dt), sfft.fftfreq(n, dt), atol=1e-12)
    f = impl.fft_frequencies(2048, 1 / 1024)
    assert f[1] == pytest.approx(0.5)                                      # resolution 1 / (N dt)
    assert np.max(np.abs(f)) == pytest.approx(512.0)                       # Nyquist f_s / 2
    with pytest.raises(ValueError):
        impl.fft_frequencies(8, 0.0)


# ---------------------------------------------------------------- amplitude spectrum, windows, leakage

def test_amplitude_spectrum_recovers_an_on_bin_sinusoid_exactly():
    fs, n = 1000.0, 1000
    t = np.arange(n) / fs
    for f0, a0, phase in ((10.0, 3.0, 0.0), (123.0, 0.7, 1.2), (499.0, 2.0, -0.4)):
        f, A = impl.amplitude_spectrum(1.5 + a0 * np.cos(2 * np.pi * f0 * t + phase), 1 / fs)
        k = np.argmin(np.abs(f - f0))
        assert f[k] == pytest.approx(f0)
        assert A[k] == pytest.approx(a0, rel=1e-10)                        # the sinusoid's amplitude, in its units
        assert A[0] == pytest.approx(1.5, rel=1e-10)                       # the mean
        others = np.delete(A, [0, k])
        assert np.max(others) < 1e-9
    assert f.size == n // 2 + 1 and f[-1] == pytest.approx(fs / 2)
    # An odd length has no Nyquist bin.
    f_odd, A_odd = impl.amplitude_spectrum(np.ones(101), 1.0)
    assert f_odd.size == 51


def test_hann_window_has_coherent_gain_one_half_and_limits_leakage():
    w = impl.hann_window(64)
    assert np.mean(w) == pytest.approx(0.5)
    assert w[0] == pytest.approx(0.0) and np.max(w) <= 1.0 + 1e-12
    with pytest.raises(ValueError):
        impl.hann_window(0)
    fs, n = 1024.0, 2048
    t = np.arange(n) / fs
    # A tone exactly half-way between bins: the rectangular window reads it 2/pi low, the Hann window about 15 % low.
    x = 2.0 * np.cos(2 * np.pi * 83.25 * t)
    _, A_rect = impl.amplitude_spectrum(x, 1 / fs)
    _, A_hann = impl.amplitude_spectrum(x, 1 / fs, window="hann")
    assert np.max(A_rect) == pytest.approx(2.0 * 2.0 / np.pi, rel=0.01)
    assert 0.80 * 2.0 < np.max(A_hann) < 0.90 * 2.0
    # On a bin both windows are exact.
    x = 2.0 * np.cos(2 * np.pi * 83.0 * t)
    assert np.max(impl.amplitude_spectrum(x, 1 / fs, window="hann")[1]) == pytest.approx(2.0, rel=1e-9)
    # Far from the tone, the Hann spectrum is lower (sidelobes fall faster).
    f, A_rect = impl.amplitude_spectrum(2.0 * np.cos(2 * np.pi * 83.25 * t), 1 / fs)
    _, A_hann = impl.amplitude_spectrum(2.0 * np.cos(2 * np.pi * 83.25 * t), 1 / fs, window="hann")
    far = np.abs(f - 83.25) > 10.0
    assert np.max(A_hann[far]) < 0.1 * np.max(A_rect[far])
    with pytest.raises(ValueError):
        impl.amplitude_spectrum(x, 1 / fs, window="blackman-harris-7")


def test_rms_from_spectrum_equals_time_domain_rms():
    fs, n = 1024.0, 2048
    t, a = impl.synthetic_vibration(fs, n / fs, seed=3)
    _, A = impl.amplitude_spectrum(a, 1 / fs)
    assert impl.rms_from_spectrum(A, a.size) == pytest.approx(np.sqrt(np.mean(a ** 2)), rel=1e-10)
    # Noise-free: RMS is sqrt(sum a_k^2 / 2) of the tones.
    _, a0 = impl.synthetic_vibration(fs, n / fs, noise_std=0.0)
    amps = np.array([amp for _, amp, _ in impl.vibration_components()])
    assert np.sqrt(np.mean(a0 ** 2)) == pytest.approx(np.sqrt(np.sum(amps ** 2) / 2), rel=2e-3)


# ---------------------------------------------------------------- the synthetic record and its peaks

def test_synthetic_record_is_reproducible_and_peaks_identify_the_tones():
    fs, dur = 1024.0, 2.0
    t1, a1 = impl.synthetic_vibration(fs, dur, seed=10)
    t2, a2 = impl.synthetic_vibration(fs, dur, seed=10)
    _, a3 = impl.synthetic_vibration(fs, dur, seed=11)
    npt.assert_array_equal(a1, a2)
    assert not np.allclose(a1, a3)
    assert t1.size == 2048 and t1[1] == pytest.approx(1 / fs)
    f, A = impl.amplitude_spectrum(a1, 1 / fs, window="hann")
    peaks = impl.dominant_peaks(f, A, 3, min_separation=5.0)
    found = sorted(fp for fp, _ in peaks)
    true = sorted(f0 for f0, _, _ in impl.vibration_components())
    npt.assert_allclose(found, true, atol=0.5)                              # to within one bin
    for fp, ap in peaks:
        a_true = [amp for f0, amp, _ in impl.vibration_components() if abs(f0 - fp) <= 0.5][0]
        assert 0.8 * a_true <= ap <= 1.05 * a_true                          # Hann: at most ~15 % low, noise 0.5 m/s^2
    with pytest.raises(ValueError):
        impl.dominant_peaks(f, A, 0)
    with pytest.raises(ValueError):
        impl.synthetic_vibration(0.0, 1.0)


def test_dominant_peaks_skips_sidelobes_within_the_separation():
    f = np.arange(0.0, 50.0, 0.5)
    A = np.zeros_like(f)
    A[20] = 5.0; A[19] = 1.0; A[21] = 1.2                                   # a tone with two leakage sidelobes
    A[80] = 2.0
    peaks = impl.dominant_peaks(f, A, 2, min_separation=2.0)
    assert [p[0] for p in peaks] == [10.0, 40.0]
    assert impl.dominant_peaks(f, A, 3, min_separation=0.0)[1] == (40.0, 2.0)


# ---------------------------------------------------------------- aliasing

def test_alias_frequency_folds_into_the_nyquist_band():
    fs = 120.0
    assert impl.alias_frequency(25.0, fs) == pytest.approx(25.0)
    assert impl.alias_frequency(83.3, fs) == pytest.approx(36.7)
    assert impl.alias_frequency(120.0, fs) == pytest.approx(0.0)
    assert impl.alias_frequency(60.0, fs) == pytest.approx(60.0)
    assert impl.alias_frequency(250.0, 100.0) == pytest.approx(50.0)
    assert impl.alias_frequency(130.0, 100.0) == pytest.approx(30.0)
    for f0 in RNG.uniform(0, 1000, 20):
        fa = impl.alias_frequency(f0, fs)
        assert 0.0 <= fa <= fs / 2
    with pytest.raises(ValueError):
        impl.alias_frequency(10.0, -1.0)
    # The spectrum of the undersampled record peaks where the formula says.
    t, a = impl.synthetic_vibration(fs, 2.0, noise_std=0.0)
    f, A = impl.amplitude_spectrum(a, 1 / fs)
    peaks = sorted(fp for fp, _ in impl.dominant_peaks(f, A, 3, min_separation=3.0))
    npt.assert_allclose(peaks, sorted(impl.alias_frequency(f0, fs) for f0, _, _ in impl.vibration_components()), atol=0.5)


# ---------------------------------------------------------------- the continuous transform

def test_sampled_transform_converges_to_the_continuous_one_at_first_and_second_order():
    a = 50.0
    errs_rect, errs_trap = [], []
    for dt in (1 / 1024, 1 / 2048, 1 / 4096):
        t = np.arange(0.0, 2.0, dt)
        x = np.exp(-a * t)
        for trap, store in ((False, errs_rect), (True, errs_trap)):
            f, X = impl.ft_from_samples(x, dt, trapezoid=trap)
            sel = (f >= 0) & (f <= 100.0)
            exact = impl.ft_exp_decay(f[sel], a)
            store.append(np.max(np.abs(X[sel] - exact) / np.abs(exact)))
    errs_rect, errs_trap = np.array(errs_rect), np.array(errs_trap)
    npt.assert_allclose(errs_rect[:-1] / errs_rect[1:], 2.0, rtol=0.05)   # rectangle sum: O(dt)
    npt.assert_allclose(errs_trap[:-1] / errs_trap[1:], 4.0, rtol=0.10)   # trapezoid: O(dt^2)
    assert errs_trap[-1] < 3e-3
    # At f = 0 the transform is the integral of x: 1/a.
    f, X = impl.ft_from_samples(np.exp(-a * np.arange(0.0, 2.0, 1 / 4096)), 1 / 4096, trapezoid=True)
    assert X[0].real == pytest.approx(1 / a, rel=1e-4)
    assert impl.ft_exp_decay(0.0, a) == pytest.approx(1 / a)
    assert abs(impl.ft_exp_decay(a / (2 * np.pi), a)) == pytest.approx(1 / (a * np.sqrt(2)))   # -3 dB at f = a / 2 pi
    with pytest.raises(ValueError):
        impl.ft_exp_decay(np.array([1.0]), -1.0)


def test_rectangular_pulse_transforms_to_a_sinc():
    T, dt = 0.1, 1 / 4096
    t = np.arange(0.0, 2.0, dt)
    x = (t < T).astype(float)
    f, X = impl.ft_from_samples(x, dt, trapezoid=True)
    sel = (f >= 0) & (f <= 100.0)
    exact = T * np.sinc(f[sel] * T) * np.exp(-1j * np.pi * f[sel] * T)
    npt.assert_allclose(X[sel], exact, atol=1e-4)
    assert abs(X[np.argmin(np.abs(f - 1 / T))]) < 1e-4                    # first zero of the sinc at f = 1 / T


# ---------------------------------------------------------------- convolution theorem

def test_convolution_via_fft_matches_numpy_convolve():
    for nx, nh in ((10, 3), (100, 37), (257, 64)):
        x, h = RNG.standard_normal(nx), RNG.standard_normal(nh)
        npt.assert_allclose(impl.convolve_fft(x, h), np.convolve(x, h), atol=1e-9)
    # Convolution with a unit impulse is the identity; with a delayed impulse, a delay.
    x = RNG.standard_normal(50)
    npt.assert_allclose(impl.convolve_fft(x, [1.0]), x, atol=1e-10)
    npt.assert_allclose(impl.convolve_fft(x, [0.0, 0.0, 1.0])[2:], x, atol=1e-10)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "vibration_fft.png").stat().st_size > 1000
