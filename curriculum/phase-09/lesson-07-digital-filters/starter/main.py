"""P09.07 Digital Filters: starter code.

After the anti-aliasing filter and the sampler of P09.06, a gyro channel is a stream of numbers at
f_s = 1000 Hz. It still carries a 180 Hz vibration line (below Nyquist, so not an alias, but far above
the 0-10 Hz the attitude loops need) and broadband sensor noise. The flight software removes them with
a digital filter: a difference equation run once per sample. This module builds the common ones from
scratch and measures what each costs in delay:

* evaluation of a digital frequency response H(e^{j omega}) from the coefficients;
* the direct-form difference equation, and the impulse response obtained by running it on a unit
  impulse;
* FIR filters: the moving average (zeros at multiples of f_s/N) and the windowed-sinc low-pass
  (rectangular, Hann, Hamming, Blackman windows), both with exactly linear phase;
* IIR filters: the bilinear (Tustin) transform s = K (1 - z^-1)/(1 + z^-1), with K = 2 f_s or prewarped
  so that one chosen frequency maps exactly, and the digital Butterworth built from it;
* group delay from the coefficients, and pole-based stability of a recursive filter;
* a method comparison on a SYNTHETIC gyro record (seeded; not measured on any vehicle).

The tests check every response against ``scipy.signal.freqz`` and against the impulse response of the
difference equation (its DFT is the frequency response), the bilinear transform against
``scipy.signal.cont2discrete(method='bilinear')`` and the digital Butterworth against
``scipy.signal.butter(fs=...)``.

Conventions: coefficient arrays b, a are in powers of z^-1 (b[0] + b[1] z^-1 + ...), a[0] != 0;
analog coefficients are in powers of s, highest first. Frequencies f in Hz, f_s in Hz.

Run:  python main.py [out_dir]     (prints the lesson's tables, writes digital_filters.png)
Test: pytest curriculum/phase-09/lesson-07-digital-filters
Complete every function whose body is `raise NotImplementedError`; every other function
and main() are given and complete.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

FS_GYRO = 1000.0     # Hz, gyro sample rate (typical of a small-UAV flight controller)
F_VIB = 180.0        # Hz, vibration line (representative motor/propeller harmonic)


# ---------------------------------------------------------------- given helpers

def poly_from_roots(roots: np.ndarray) -> np.ndarray:
    """Monic polynomial (highest power first) with the given roots; real when the roots are conjugate pairs."""
    c = np.array([1.0 + 0j])
    for r in np.atleast_1d(roots):
        c = np.convolve(c, np.array([1.0, -r]))
    if np.max(np.abs(c.imag), initial=0.0) <= 1e-9 * np.max(np.abs(c)):
        return c.real.copy()
    return c


def butterworth_analog(n: int, wc: float) -> tuple[np.ndarray, np.ndarray]:
    """Analog Butterworth low-pass of order n and corner wc [rad/s] (the prototype of P09.06)."""
    k = np.arange(1, n + 1)
    poles = wc * np.exp(1j * np.pi * (2 * k + n - 1) / (2 * n))
    return np.array([wc ** n]), np.real(poly_from_roots(poles))


def analog_freq_response(b: np.ndarray, a: np.ndarray, f: np.ndarray | float) -> np.ndarray:
    """H(j 2 pi f) for analog coefficients in powers of s (highest first)."""
    s = 2j * np.pi * np.asarray(f, dtype=float)
    return np.polyval(b, s) / np.polyval(a, s)


def tone_amplitude(x: np.ndarray, fs: float, f: float) -> float:
    """Amplitude of the component of x at frequency f (single-bin DFT)."""
    x = np.asarray(x, dtype=float)
    n = np.arange(x.size)
    return float(2.0 / x.size * np.abs(np.sum(x * np.exp(-2j * np.pi * f * n / fs))))


def synthetic_gyro_record(seed: int = 7, duration: float = 4.0, fs: float = FS_GYRO,
                          a_vib: float = 3.0, noise_rms: float = 0.3) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (t, true_rate, vibration, noise) in s and deg/s. SYNTHETIC, not measured.

    true_rate = 5 sin(2 pi 0.5 t) + 2 sin(2 pi 2 t); vibration = a_vib sin(2 pi 180 t + 0.2);
    noise is white and Gaussian with standard deviation noise_rms (seeded). The sensed rate is the sum.
    """
    rng = np.random.default_rng(seed)
    t = np.arange(int(round(duration * fs))) / fs
    true = 5.0 * np.sin(2 * np.pi * 0.5 * t) + 2.0 * np.sin(2 * np.pi * 2.0 * t)
    vib = a_vib * np.sin(2 * np.pi * F_VIB * t + 0.2)
    noise = noise_rms * rng.standard_normal(t.size)
    return t, true, vib, noise


# ---------------------------------------------------------------- evaluation and the difference equation

def digital_freq_response(b: np.ndarray, a: np.ndarray, f: np.ndarray | float, fs: float) -> np.ndarray:
    """Return H(e^{j omega}) = B(e^{j omega}) / A(e^{j omega}), omega = 2 pi f / fs, coefficients in powers of z^-1."""
    raise NotImplementedError


def difference_equation(b: np.ndarray, a: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Run a[0] y[n] = sum_k b[k] x[n-k] - sum_{k>=1} a[k] y[n-k] from rest (direct form I)."""
    raise NotImplementedError


def impulse_response(b: np.ndarray, a: np.ndarray, n: int) -> np.ndarray:
    """First n samples of the impulse response, by running the difference equation on delta[n]."""
    raise NotImplementedError


def group_delay(b: np.ndarray, a: np.ndarray, f: np.ndarray | float, fs: float) -> np.ndarray:
    """Return the group delay -d(arg H)/d(omega) in samples at frequency f.

    For C(omega) = sum c_k e^{-j omega k}, -d arg C / d omega = Re[ sum k c_k e^{-j omega k} / C ];
    the delay of B/A is that of B minus that of A.
    """
    raise NotImplementedError


def is_stable_digital(a: np.ndarray) -> bool:
    """True if every root of a(z) (the poles of a recursive filter) lies strictly inside the unit circle."""
    raise NotImplementedError


# ---------------------------------------------------------------- FIR design

def moving_average(n_taps: int) -> np.ndarray:
    """Return the n_taps-point moving average b = [1/N, ..., 1/N] (a = [1])."""
    raise NotImplementedError


def window(name: str, n: int) -> np.ndarray:
    """Return the symmetric window of length n: 'rectangular', 'hann', 'hamming' or 'blackman'.

    w[k] for k = 0..n-1 with x = 2 pi k / (n - 1): hann 0.5 - 0.5 cos x; hamming 0.54 - 0.46 cos x;
    blackman 0.42 - 0.5 cos x + 0.08 cos 2x.
    """
    raise NotImplementedError


def windowed_sinc_lowpass(n_taps: int, fc: float, fs: float, win: str = "hamming") -> np.ndarray:
    """Return a linear-phase FIR low-pass of n_taps (odd) taps with cutoff fc [Hz].

    h[k] = 2 (fc/fs) sinc(2 (fc/fs) (k - M)), M = (n_taps - 1)/2, times the window, scaled to unit DC gain.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- IIR design: the bilinear transform

def bilinear(b: np.ndarray, a: np.ndarray, fs: float, f_prewarp: float | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Map the analog H(s) = b(s)/a(s) to a digital filter by s = K (1 - z^-1)/(1 + z^-1).

    K = 2 fs (Tustin) or, with f_prewarp, K = 2 pi f_p / tan(pi f_p / fs), so that the digital response
    at f_p equals the analog response at f_p exactly. Returns (bz, az) in powers of z^-1, az[0] = 1.
    """
    raise NotImplementedError


def butterworth_digital(n: int, fc: float, fs: float) -> tuple[np.ndarray, np.ndarray]:
    """Digital Butterworth low-pass of order n with its -3 dB point exactly at fc [Hz]:
    the analog prototype at 2 pi fc, mapped by the bilinear transform prewarped at fc."""
    raise NotImplementedError


# ---------------------------------------------------------------- the lesson's numbers

def lesson_filters(fs: float = FS_GYRO) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """The four filters compared in the lesson: (b, a) for each."""
    return {
        "MA 50": (moving_average(50), np.array([1.0])),
        "FIR 101 Hamming 40 Hz": (windowed_sinc_lowpass(101, 40.0, fs, "hamming"), np.array([1.0])),
        "IIR Butterworth 2, 40 Hz": butterworth_digital(2, 40.0, fs),
        "IIR Butterworth 4, 40 Hz": butterworth_digital(4, 40.0, fs),
    }


def gyro_filter_comparison(seed: int = 7) -> list[dict]:
    """Run each lesson filter on the SYNTHETIC gyro record; return per filter: name, atten_180_db,
    delay_2hz_ms, residual_vib (deg/s at 180 Hz in the output), noise_out (rms of the filtered noise),
    noise_pred (noise_rms * sqrt(sum h^2), Parseval), rms_err (output - true rate). The first 0.5 s
    (start-up transient) is discarded."""
    t, true, vib, noise = synthetic_gyro_record(seed)
    x = true + vib + noise
    keep = t >= 0.5
    rows = []
    for name, (b, a) in lesson_filters().items():
        y = difference_equation(b, a, x)
        yn = difference_equation(b, a, noise)
        h = impulse_response(b, a, 4000)
        rows.append({
            "name": name,
            "atten_180_db": float(20 * np.log10(abs(digital_freq_response(b, a, F_VIB, FS_GYRO)))),
            "delay_2hz_ms": float(1e3 * group_delay(b, a, 2.0, FS_GYRO) / FS_GYRO),
            "residual_vib": tone_amplitude(difference_equation(b, a, vib)[keep], FS_GYRO, F_VIB),
            "noise_out": float(np.std(yn[keep])),
            "noise_pred": float(0.3 * np.sqrt(np.sum(h ** 2))),
            "rms_err": float(np.sqrt(np.mean((y - true)[keep] ** 2))),
        })
    return rows


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's tables and save digital_filters.png in out_dir."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fs = FS_GYRO

    b_ma = moving_average(50)
    print(f"moving average, 50 taps at {fs:.0f} Hz: zeros at multiples of {fs / 50:.0f} Hz; "
          f"|H(180 Hz)| = {abs(digital_freq_response(b_ma, [1.0], 180.0, fs)):.2e}, |H(10 Hz)| = "
          f"{abs(digital_freq_response(b_ma, [1.0], 10.0, fs)):.4f}, delay {group_delay(b_ma, [1.0], 2.0, fs):.1f} samples")
    for win in ("rectangular", "hann", "hamming", "blackman"):
        h = windowed_sinc_lowpass(101, 40.0, fs, win)
        ff = np.linspace(80.0, 500.0, 4000)
        worst = 20 * np.log10(np.max(np.abs(digital_freq_response(h, [1.0], ff, fs))))
        print(f"windowed sinc, 101 taps, fc 40 Hz, {win:11s}: |H(40)| = {abs(digital_freq_response(h, [1.0], 40.0, fs)):.3f}, "
              f"worst gain above 80 Hz {worst:7.1f} dB, |H(180)| {20 * np.log10(abs(digital_freq_response(h, [1.0], 180.0, fs))):7.1f} dB")

    print("\nbilinear transform of the analog Butterworth 2, fc = 40 Hz and fc = 250 Hz, fs = 1000 Hz:")
    for fc in (40.0, 250.0):
        ba, aa = butterworth_analog(2, 2 * np.pi * fc)
        b0, a0 = bilinear(ba, aa, fs)
        b1, a1 = bilinear(ba, aa, fs, f_prewarp=fc)
        print(f"  fc = {fc:5.0f} Hz: |H(fc)| analog {abs(analog_freq_response(ba, aa, fc)):.4f}, Tustin "
              f"{abs(digital_freq_response(b0, a0, fc, fs)):.4f}, prewarped {abs(digital_freq_response(b1, a1, fc, fs)):.4f}; "
              f"Tustin maps the analog corner to {fs / np.pi * np.arctan(np.pi * fc / fs):.1f} Hz")
    bz, az = butterworth_digital(2, 40.0, fs)
    print("  Butterworth 2 at 40 Hz: b =", np.array2string(bz, precision=6), " a =", np.array2string(az, precision=6))
    print("  poles |z| =", np.array2string(np.abs(np.roots(az)), precision=5), " stable:", is_stable_digital(az))

    print("\nSYNTHETIC gyro record (not a measurement): 0.5 + 2 Hz manoeuvre, 3 deg/s at 180 Hz, 0.3 deg/s white noise, fs 1 kHz")
    rows = gyro_filter_comparison()
    for r in rows:
        att = "   zero  " if r["atten_180_db"] < -200 else f"{r['atten_180_db']:7.1f} dB"
        print(f"  {r['name']:26s} |H(180)| {att}, residual {r['residual_vib']:.4f} deg/s, "
              f"delay at 2 Hz {r['delay_2hz_ms']:5.1f} ms, noise {r['noise_out']:.4f} (Parseval {r['noise_pred']:.4f}) deg/s, "
              f"rms error vs true {r['rms_err']:.3f} deg/s")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    ff = np.linspace(0.1, fs / 2, 2000)
    for name, (b, a) in lesson_filters().items():
        axes[0].plot(ff, 20 * np.log10(np.abs(digital_freq_response(b, a, ff, fs)) + 1e-12), label=name, lw=1)
        axes[1].plot(ff[ff < 100], 1e3 * group_delay(b, a, ff[ff < 100], fs) / fs, label=name, lw=1)
    axes[0].axvline(F_VIB, color="r", ls=":", lw=1)
    axes[0].set_ylim(-120, 5); axes[0].set_xlabel("f [Hz]"); axes[0].set_ylabel("|H| [dB]"); axes[0].legend(fontsize=7)
    axes[0].grid(lw=0.3); axes[0].set_title("magnitude", fontsize=9)
    axes[1].set_xlabel("f [Hz]"); axes[1].set_ylabel("group delay [ms]"); axes[1].grid(lw=0.3); axes[1].set_title("delay", fontsize=9)
    t, true, vib, noise = synthetic_gyro_record()
    x = true + vib + noise
    sel = (t >= 1.0) & (t < 1.6)
    axes[2].plot(t[sel], x[sel], color="0.7", lw=0.5, label="sensed (SYNTHETIC)")
    axes[2].plot(t[sel], true[sel], "k--", lw=1, label="true rate")
    for name in ("FIR 101 Hamming 40 Hz", "IIR Butterworth 4, 40 Hz"):
        b, a = lesson_filters()[name]
        axes[2].plot(t[sel], difference_equation(b, a, x)[sel], lw=1, label=name)
    axes[2].set_xlabel("t [s]"); axes[2].set_ylabel("rate [deg/s]"); axes[2].legend(fontsize=7)
    axes[2].set_title("SYNTHETIC gyro record (not measured)", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "digital_filters.png", dpi=130)
    plt.close(fig)
    print(f"saved {out_dir / 'digital_filters.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
