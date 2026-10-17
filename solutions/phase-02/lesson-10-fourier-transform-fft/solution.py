"""P02.10 Fourier Transform & the FFT: reference solution.

Spectrum of a vibration record from a (synthetic) accelerometer, computed with a discrete Fourier
transform written from scratch (the O(N^2) matrix form), a radix-2 Cooley-Tukey FFT written from
scratch, and ``scipy.fft`` as the oracle.

Conventions:
    * DFT:  X[k] = sum_{n=0}^{N-1} x[n] exp(-2 pi i k n / N),  k = 0 .. N-1  (no 1/N in front);
      inverse x[n] = (1/N) sum_k X[k] exp(+2 pi i k n / N). This is the ``scipy.fft`` convention.
    * Continuous transform:  X(f) = integral x(t) exp(-2 pi i f t) dt, frequency f in hertz.
      Sampling at spacing dt turns the integral into a Riemann sum, so  X(f_k) ~ dt * DFT[k].
    * Single-sided amplitude spectrum: A[k] = 2 |X[k]| / N for 0 < k < N/2 and A[0] = |X[0]| / N,
      so a sinusoid of amplitude a sitting exactly on bin k gives A[k] = a.
    * Every "measurement" in this lesson is synthetic (``synthetic_vibration``) and labelled as such.

Run:  python solution.py          (prints the spectrum table, writes vibration_fft.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-02/lesson-10-fourier-transform-fft
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


# ---------------------------------------------------------------- the discrete Fourier transform

def dft_matrix(n: int) -> np.ndarray:
    """Return the n x n DFT matrix F with entries F[k, m] = exp(-2 pi i k m / n).

    F is symmetric and F^H F = n I (it is sqrt(n) times a unitary matrix), so F^{-1} = conj(F) / n.
    Raise ``ValueError`` if n < 1.
    """
    if n < 1:
        raise ValueError("n must be a positive integer")
    k = np.arange(n)
    return np.exp(-2j * np.pi * np.outer(k, k) / n)


def dft(x: np.ndarray) -> np.ndarray:
    """Return the DFT of ``x`` by multiplying with the DFT matrix: O(N^2) operations, no tricks.

    X[k] = sum_n x[n] exp(-2 pi i k n / N). Works for any N; this is the definition the FFT must
    reproduce.
    """
    x = np.asarray(x, dtype=complex)
    if x.ndim != 1:
        raise ValueError("dft expects a one-dimensional sequence")
    return dft_matrix(x.size) @ x


def idft(X: np.ndarray) -> np.ndarray:
    """Return the inverse DFT, x[n] = (1/N) sum_k X[k] exp(+2 pi i k n / N), via conj(F) X / N."""
    X = np.asarray(X, dtype=complex)
    n = X.size
    return np.conj(dft_matrix(n)) @ X / n


def fft_radix2(x: np.ndarray) -> np.ndarray:
    """Return the DFT of ``x`` by the recursive radix-2 Cooley-Tukey algorithm.

    Split x into even and odd samples, transform each half (recursively), then combine with the
    twiddle factors w_k = exp(-2 pi i k / N):  X[k] = E[k] + w_k O[k],  X[k + N/2] = E[k] - w_k O[k].
    Cost O(N log N). Raise ``ValueError`` unless N is a power of two (N >= 1).
    """
    x = np.asarray(x, dtype=complex)
    n = x.size
    if n < 1 or (n & (n - 1)) != 0:
        raise ValueError("fft_radix2 needs a length that is a power of two")
    if n == 1:
        return x.copy()
    even = fft_radix2(x[0::2])
    odd = fft_radix2(x[1::2])
    twiddle = np.exp(-2j * np.pi * np.arange(n // 2) / n) * odd
    return np.concatenate([even + twiddle, even - twiddle])


def fft_any(x: np.ndarray) -> np.ndarray:
    """Return the DFT of ``x``: the radix-2 FFT when the length is a power of two, else the O(N^2) DFT."""
    x = np.asarray(x, dtype=complex)
    n = x.size
    if n >= 1 and (n & (n - 1)) == 0:
        return fft_radix2(x)
    return dft(x)


def fft_frequencies(n: int, dt: float) -> np.ndarray:
    """Return the frequency [Hz] of each DFT bin for N samples spaced dt apart, in scipy.fft.fftfreq order.

    f_k = k / (N dt) for k = 0 .. N/2 - 1, then the negative frequencies. The bin spacing 1 / (N dt)
    is the frequency resolution; f_s / 2 = 1 / (2 dt) is the Nyquist frequency.
    """
    if n < 1 or dt <= 0.0:
        raise ValueError("n must be positive and dt > 0")
    k = np.arange(n)
    k = np.where(k < (n + 1) // 2, k, k - n)
    return k / (n * dt)


# ---------------------------------------------------------------- windows and the amplitude spectrum

def hann_window(n: int) -> np.ndarray:
    """Return the periodic Hann window w[m] = 0.5 (1 - cos(2 pi m / N)), m = 0 .. N-1.

    Its mean (coherent gain) is exactly 0.5, so a windowed spectrum must be divided by 0.5 to recover
    sinusoid amplitudes.
    """
    if n < 1:
        raise ValueError("n must be positive")
    m = np.arange(n)
    return 0.5 * (1.0 - np.cos(2.0 * np.pi * m / n))


def amplitude_spectrum(x: np.ndarray, dt: float, window: str | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Return (f, A): single-sided amplitude spectrum of the real record ``x`` sampled at spacing dt.

    A[k] = 2 |X[k]| / (N g) for 0 < k < N/2 and A[0] = |X[0]| / (N g), where g is the coherent gain of
    the window (1 for none, 0.5 for 'hann'). A sinusoid of amplitude a exactly on a bin returns a.
    ``window`` is None (rectangular) or 'hann'. The DFT is computed with ``fft_any``.
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    if window is None:
        w, gain = np.ones(n), 1.0
    elif window == "hann":
        w, gain = hann_window(n), 0.5
    else:
        raise ValueError(f"unknown window {window!r}")
    X = fft_any(x * w)
    half = n // 2 + 1 if n % 2 == 0 else (n + 1) // 2
    f = np.abs(fft_frequencies(n, dt)[:half])
    if n % 2 == 0:
        f[-1] = 0.5 / dt                           # the Nyquist bin, which fftfreq lists as negative
    A = np.abs(X[:half]) / (n * gain)
    A[1:] *= 2.0
    if n % 2 == 0:
        A[-1] /= 2.0                               # Nyquist bin has no mirror image
    return f, A


def parseval(x: np.ndarray, X: np.ndarray) -> tuple[float, float]:
    """Return (sum |x[n]|^2, (1/N) sum |X[k]|^2): Parseval's theorem says they are equal."""
    x = np.asarray(x, dtype=complex)
    X = np.asarray(X, dtype=complex)
    return float(np.sum(np.abs(x) ** 2)), float(np.sum(np.abs(X) ** 2) / X.size)


# ---------------------------------------------------------------- the continuous transform, sampled

def ft_exp_decay(f: np.ndarray, a: float) -> np.ndarray:
    """Return the continuous Fourier transform of x(t) = exp(-a t) for t >= 0 (zero before): 1 / (a + 2 pi i f).

    The Laplace transform of the same function is 1 / (s + a); the Fourier transform is its value on
    the imaginary axis s = 2 pi i f. Raise ``ValueError`` unless a > 0 (otherwise the integral diverges).
    """
    if a <= 0.0:
        raise ValueError("a must be positive for the transform to exist")
    return 1.0 / (a + 2j * np.pi * np.asarray(f, dtype=float))


def ft_from_samples(x: np.ndarray, dt: float, trapezoid: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """Return (f, X_hat): the Riemann-sum estimate X_hat(f_k) = dt * DFT[k] of the continuous transform.

    Frequencies in fftfreq order. With ``trapezoid=True`` the first and last samples get half weight
    (the trapezoidal rule), which removes the first-order error a jump at t = 0 leaves in the plain
    rectangle sum. The estimate is exact for the sampled signal; for the underlying continuous signal
    it carries the sampling (aliasing) and truncation errors discussed in the lesson.
    """
    x = np.asarray(x, dtype=float).copy()
    if trapezoid:
        x[0] *= 0.5
        x[-1] *= 0.5
    return fft_frequencies(x.size, dt), dt * fft_any(x)


def alias_frequency(f: float, fs: float) -> float:
    """Return the apparent frequency [Hz] at which a sinusoid of frequency f appears when sampled at fs.

    The spectrum of a sampled signal repeats every fs, so f is indistinguishable from f + m fs for any
    integer m; the single-sided spectrum shows |f - m fs| folded into [0, fs / 2].
    """
    if fs <= 0.0:
        raise ValueError("sampling frequency must be positive")
    r = np.mod(f, fs)
    return float(min(r, fs - r))


# ---------------------------------------------------------------- the vibration record (SYNTHETIC)

# GIVEN
def vibration_components() -> list[tuple[float, float, float]]:
    """Return the (frequency [Hz], amplitude [m/s^2], phase [rad]) of each tone in the synthetic record.

    A propeller at 1500 rpm gives a once-per-revolution (1P) tone at 25 Hz and its 2P harmonic at 50 Hz;
    a lightly loaded skin panel responds at its first mode, 83.3 Hz, which does not sit on a 0.5 Hz
    bin. These numbers are chosen to be plausible for a light aircraft, not measured on one.
    """
    return [(25.0, 4.0, 0.3), (50.0, 1.5, -1.1), (83.3, 2.0, 2.0)]


def synthetic_vibration(fs: float, duration: float, seed: int = 10, noise_std: float = 0.5) -> tuple[np.ndarray, np.ndarray]:
    """Return (t, a): a SYNTHETIC accelerometer record [s, m/s^2] sampled at fs for ``duration`` seconds.

    a(t) = sum of the tones in ``vibration_components`` plus Gaussian white noise of standard deviation
    ``noise_std``, from ``numpy.random.default_rng(seed)``. No physical sensor was involved.
    """
    if fs <= 0.0 or duration <= 0.0:
        raise ValueError("fs and duration must be positive")
    n = int(round(fs * duration))
    t = np.arange(n) / fs
    a = np.zeros(n)
    for f0, amp, ph in vibration_components():
        a += amp * np.cos(2.0 * np.pi * f0 * t + ph)
    rng = np.random.default_rng(seed)
    a += noise_std * rng.standard_normal(n)
    return t, a


def dominant_peaks(f: np.ndarray, A: np.ndarray, n_peaks: int, min_separation: float = 0.0) -> list[tuple[float, float]]:
    """Return the ``n_peaks`` largest local maxima of A as (frequency, amplitude), strongest first.

    A local maximum is a bin whose amplitude exceeds both neighbours; peaks closer than
    ``min_separation`` [Hz] to an already chosen stronger peak are skipped (so leakage sidelobes are
    not reported as tones). Raise ``ValueError`` if ``n_peaks`` < 1.
    """
    if n_peaks < 1:
        raise ValueError("n_peaks must be at least 1")
    f = np.asarray(f, dtype=float)
    A = np.asarray(A, dtype=float)
    interior = np.arange(1, A.size - 1)
    locmax = interior[(A[interior] > A[interior - 1]) & (A[interior] >= A[interior + 1])]
    order = locmax[np.argsort(A[locmax])[::-1]]
    chosen: list[tuple[float, float]] = []
    for k in order:
        if all(abs(f[k] - fc) >= min_separation for fc, _ in chosen):
            chosen.append((float(f[k]), float(A[k])))
        if len(chosen) == n_peaks:
            break
    return chosen


def rms_from_spectrum(A: np.ndarray, n: int) -> float:
    """Return the RMS value implied by the single-sided amplitude spectrum ``A`` of an n-sample record.

    Each sinusoid of amplitude a contributes a^2 / 2 to the mean square; the DC term contributes A0^2,
    and for even n the last bin (the Nyquist bin, which has no mirror image) contributes A^2 as well.
    Parseval's theorem says the result equals the time-domain RMS of the record.
    """
    A = np.asarray(A, dtype=float)
    ms = A[0] ** 2 + 0.5 * np.sum(A[1:] ** 2)
    if n % 2 == 0:
        ms += 0.5 * A[-1] ** 2
    return float(np.sqrt(ms))


# ---------------------------------------------------------------- convolution

def convolve_fft(x: np.ndarray, h: np.ndarray) -> np.ndarray:
    """Return the linear convolution of x and h (length Nx + Nh - 1) through the convolution theorem.

    Zero-pad both to the next power of two at least Nx + Nh - 1, multiply the FFTs, inverse transform,
    take the real part and trim. Equals ``numpy.convolve(x, h)`` to rounding.
    """
    x = np.asarray(x, dtype=float)
    h = np.asarray(h, dtype=float)
    n_lin = x.size + h.size - 1
    n = 1
    while n < n_lin:
        n *= 2
    X = fft_radix2(np.concatenate([x, np.zeros(n - x.size)]))
    H = fft_radix2(np.concatenate([h, np.zeros(n - h.size)]))
    y = np.conj(fft_radix2(np.conj(X * H))) / n     # inverse FFT via the forward one: ifft(Y) = conj(fft(conj(Y)))/N
    return np.real(y[:n_lin])


# ---------------------------------------------------------------- main

# GIVEN
def main(out_dir: str | Path = ".") -> None:
    """Analyse the synthetic vibration record, demonstrate leakage, windows and aliasing, and plot."""
    import time

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy import fft as sfft

    out_dir = Path(out_dir)
    fs, duration = 1024.0, 2.0
    t, a = synthetic_vibration(fs, duration)
    n = a.size
    dt = 1.0 / fs
    print(f"SYNTHETIC record: fs = {fs:.0f} Hz, N = {n}, duration {duration} s, resolution {1 / (n * dt):.3f} Hz,"
          f" Nyquist {fs / 2:.0f} Hz")

    # From-scratch transforms against the library.
    t0 = time.perf_counter(); X_dft = dft(a); t_dft = time.perf_counter() - t0
    t0 = time.perf_counter(); X_fft = fft_radix2(a); t_fft = time.perf_counter() - t0
    sfft.fft(a)                                                   # warm up the library once
    t0 = time.perf_counter(); X_lib = sfft.fft(a); t_lib = time.perf_counter() - t0
    print(f"max |DFT - scipy| = {np.max(np.abs(X_dft - X_lib)):.2e}, max |FFT - scipy| = {np.max(np.abs(X_fft - X_lib)):.2e}")
    print(f"time: O(N^2) DFT {t_dft * 1e3:.1f} ms, Python radix-2 FFT {t_fft * 1e3:.1f} ms, scipy.fft {t_lib * 1e6:.0f} us")
    e_t, e_f = parseval(a, X_fft)
    print(f"Parseval: sum|x|^2 = {e_t:.4f}, (1/N) sum|X|^2 = {e_f:.4f}")

    f, A_rect = amplitude_spectrum(a, dt)
    _, A_hann = amplitude_spectrum(a, dt, window="hann")
    print("\ntrue tones (frequency Hz, amplitude m/s^2):", [(f0, amp) for f0, amp, _ in vibration_components()])
    print("rectangular window peaks:", [(round(fp, 1), round(ap, 3)) for fp, ap in dominant_peaks(f, A_rect, 3, 5.0)])
    print("Hann window peaks       :", [(round(fp, 1), round(ap, 3)) for fp, ap in dominant_peaks(f, A_hann, 3, 5.0)])
    rms_t = float(np.sqrt(np.mean(a ** 2)))
    print(f"RMS from time record {rms_t:.4f} m/s^2, from rectangular spectrum {rms_from_spectrum(A_rect, n):.4f} m/s^2")

    # Aliasing: the same signal sampled too slowly.
    fs_low = 120.0
    t_low, a_low = synthetic_vibration(fs_low, duration, noise_std=0.0)
    f_low, A_low = amplitude_spectrum(a_low, 1.0 / fs_low)
    print(f"\nsampled at {fs_low:.0f} Hz (Nyquist {fs_low / 2:.0f} Hz): peaks",
          [(round(fp, 1), round(ap, 3)) for fp, ap in dominant_peaks(f_low, A_low, 3, 3.0)])
    print("predicted alias of each tone:", [(f0, round(alias_frequency(f0, fs_low), 2)) for f0, _, _ in vibration_components()])

    # Continuous transform of a decaying exponential against dt * DFT (rectangle and trapezoid sums).
    a_decay = 50.0
    print()
    for dt_c in (1 / 1024, 1 / 2048, 1 / 4096):
        tt = np.arange(0.0, 2.0, dt_c)
        for trap in (False, True):
            fc, Xc = ft_from_samples(np.exp(-a_decay * tt), dt_c, trapezoid=trap)
            sel = (fc >= 0) & (fc <= 100.0)
            exact = ft_exp_decay(fc[sel], a_decay)
            err = np.max(np.abs(Xc[sel] - exact) / np.abs(exact))
            print(f"exp(-{a_decay:.0f} t): dt = 1/{round(1 / dt_c)} s, {'trapezoid' if trap else 'rectangle'} sum,"
                  f" max relative error vs 1/(a + 2 pi i f) up to 100 Hz = {err:.5f}")

    # Figure.
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    ax = axes[0, 0]
    ax.plot(t[: int(0.2 * fs)], a[: int(0.2 * fs)], lw=0.8)
    ax.set_xlabel("time [s]"); ax.set_ylabel("acceleration [m/s$^2$]")
    ax.set_title("synthetic accelerometer record (first 0.2 s)")
    ax = axes[0, 1]
    ax.plot(f, A_rect, lw=0.8, label="rectangular")
    ax.plot(f, A_hann, lw=0.8, label="Hann (gain corrected)")
    for f0, amp, _ in vibration_components():
        ax.plot(f0, amp, "kv", ms=5)
    ax.set_xlim(0, 150); ax.set_xlabel("frequency [Hz]"); ax.set_ylabel("amplitude [m/s$^2$]")
    ax.set_title("single-sided amplitude spectrum"); ax.legend()
    ax = axes[1, 0]
    ax.semilogy(f, A_rect + 1e-6, lw=0.8, label="rectangular")
    ax.semilogy(f, A_hann + 1e-6, lw=0.8, label="Hann")
    ax.set_xlim(60, 110); ax.set_ylim(1e-3, 5); ax.set_xlabel("frequency [Hz]"); ax.set_ylabel("amplitude [m/s$^2$]")
    ax.set_title("leakage around the off-bin 83.3 Hz tone"); ax.legend()
    ax = axes[1, 1]
    ax.plot(f_low, A_low, lw=1.0)
    for f0, amp, _ in vibration_components():
        ax.plot(alias_frequency(f0, fs_low), amp, "rx", ms=8)
    ax.set_xlabel("frequency [Hz]"); ax.set_ylabel("amplitude [m/s$^2$]")
    ax.set_title(f"same signal sampled at {fs_low:.0f} Hz: aliases (x = predicted)")
    fig.suptitle("P02.10  FFT of a synthetic vibration record")
    fig.tight_layout()
    fig.savefig(out_dir / "vibration_fft.png", dpi=150)
    print(f"\nsaved {out_dir / 'vibration_fft.png'}")


if __name__ == "__main__":
    main()
