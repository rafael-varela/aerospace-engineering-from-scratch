"""P09.06 Sampling & Aliasing: reference solution.

A MEMS rate gyro on a small UAV is sampled by the flight computer at f_s = 400 Hz. Its signal band
(the motion the attitude loops must see) is 0-30 Hz, but the airframe also shakes at the propeller
blade-pass frequency, here a 390 Hz line of 5 deg/s. Sampled without protection, that line lands at
10 Hz, in the middle of the control band, and nothing done after the sampler can remove it. This
module builds the tools to predict that and to prevent it:

* the aliased-frequency calculator f_a = |f - f_s round(f / f_s)|;
* the Butterworth analog low-pass prototype from its poles on a circle, checked against
  ``scipy.signal.butter(analog=True)`` in the tests;
* order and corner selection from a passband droop and a stopband attenuation, with the stopband
  edge at f_s - f_band (the lowest frequency that folds into the band);
* the group delay of a rational filter at DC (the price the control loop pays);
* an analog-filter simulator (companion-form state space, exact zero-order-hold steps on a fine
  grid) so the filter can be run on a SYNTHETIC gyro record and the alias measured after sampling;
* the FFT peak of a sampled tone and Whittaker-Shannon (sinc) reconstruction.

Every gyro record in this lesson is SYNTHETIC: generated from a seeded random-number generator and a
stated formula, not measured on any aircraft.

Units: hertz for frequencies (f), rad/s for angular frequencies (w), seconds, deg/s for gyro rates.

Run:  python solution.py [out_dir]     (prints the lesson's tables, writes sampling_aliasing.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-09/lesson-06-sampling-aliasing
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# The lesson's representative IMU channel (typical numbers for a small multirotor or UAV, not data
# from a specific vehicle).
FS_IMU = 400.0          # Hz, gyro sample rate seen by the flight computer
F_BAND = 30.0           # Hz, top of the signal band the attitude loops need
F_VIB = 390.0           # Hz, propeller blade-pass vibration line
A_VIB = 5.0             # deg/s, amplitude of that line at the gyro
A_ALLOW = 0.02          # deg/s, largest alias we accept (about the gyro's own noise in the band)
DROOP_DB = 0.5          # dB, largest loss allowed at F_BAND


# ---------------------------------------------------------------- given helpers

def poly_from_roots(roots: np.ndarray) -> np.ndarray:
    """Return the monic polynomial (highest power first) with the given roots.

    Built by repeated convolution with (s - r). When the roots come in conjugate pairs the result is
    real to rounding, and the imaginary part is dropped.
    """
    c = np.array([1.0 + 0j])
    for r in np.atleast_1d(roots):
        c = np.convolve(c, np.array([1.0, -r]))
    if np.max(np.abs(c.imag), initial=0.0) <= 1e-9 * np.max(np.abs(c)):
        return c.real.copy()
    return c


def tf_to_companion(b: np.ndarray, a: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Return (A, B, C, D) of the controllable companion form of b(s)/a(s) (proper, highest power first).

    x' = A x + B u, y = C x + D u, with A's first row -a[1:]/a[0] and ones on the subdiagonal.
    """
    b = np.atleast_1d(np.asarray(b, dtype=float))
    a = np.atleast_1d(np.asarray(a, dtype=float))
    n = len(a) - 1
    if len(b) > len(a):
        raise ValueError("improper transfer function")
    b = np.concatenate([np.zeros(len(a) - len(b)), b]) / a[0]
    a = a / a[0]
    D = float(b[0])
    A = np.zeros((n, n))
    if n:
        A[0, :] = -a[1:]
        A[1:, :-1] = np.eye(n - 1)
    B = np.zeros((n, 1))
    if n:
        B[0, 0] = 1.0
    C = (b[1:] - D * a[1:]).reshape(1, n)
    return A, B, C, D


def expm_series(M: np.ndarray, n_terms: int = 20) -> np.ndarray:
    """Matrix exponential by scaling and squaring of a truncated Taylor series (from P02.05)."""
    M = np.asarray(M, dtype=float)
    norm = np.linalg.norm(M, 1)
    s = int(np.ceil(np.log2(norm / 0.5))) if norm > 0.5 else 0
    X = M / 2.0 ** s
    term = np.eye(len(M))
    E = np.eye(len(M))
    for k in range(1, n_terms):
        term = term @ X / k
        E = E + term
    for _ in range(s):
        E = E @ E
    return E


def discretize_zoh(A: np.ndarray, B: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """Exact zero-order-hold step: (A_d, B_d) from exp([[A, B], [0, 0]] dt) (the augmented matrix of P02.05)."""
    n, m = B.shape
    M = np.zeros((n + m, n + m))
    M[:n, :n] = A
    M[:n, n:] = B
    E = expm_series(M * dt)
    return E[:n, :n], E[:n, n:]


def tone_amplitude(x: np.ndarray, fs: float, f: float) -> float:
    """Amplitude of the component of x at frequency f (single-bin DFT, 2/N |sum x e^{-j 2 pi f n / fs}|)."""
    x = np.asarray(x, dtype=float)
    n = np.arange(x.size)
    return float(2.0 / x.size * np.abs(np.sum(x * np.exp(-2j * np.pi * f * n / fs))))


def synthetic_imu_record(t: np.ndarray, seed: int = 9, a_vib: float = A_VIB, f_vib: float = F_VIB,
                         noise_rms: float = 0.05) -> tuple[np.ndarray, np.ndarray]:
    """Return (true_rate, sensed_rate) in deg/s on the fine time grid t. SYNTHETIC, not measured.

    true_rate = 2 sin(2 pi 1 t) + 1 sin(2 pi 3 t)   (a slow manoeuvre, periodic in any whole second)
    sensed    = true_rate + a_vib sin(2 pi f_vib t + 0.3) + 0.4 a_vib sin(2 pi 2 f_vib t + 1.1)
                + white noise of standard deviation noise_rms per fine-grid point (seeded).
    """
    rng = np.random.default_rng(seed)
    true = 2.0 * np.sin(2 * np.pi * 1.0 * t) + 1.0 * np.sin(2 * np.pi * 3.0 * t)
    vib = a_vib * np.sin(2 * np.pi * f_vib * t + 0.3) + 0.4 * a_vib * np.sin(2 * np.pi * 2 * f_vib * t + 1.1)
    return true, true + vib + noise_rms * rng.standard_normal(t.size)


# ---------------------------------------------------------------- aliasing

def aliased_frequency(f: np.ndarray | float, fs: float) -> np.ndarray | float:
    """Return the apparent frequency in [0, fs/2] of a sinusoid at f [Hz] sampled at fs [Hz].

    f_a = |f - fs round(f / fs)|: the sampled sinusoid is indistinguishable from one at f_a.
    """
    f = np.asarray(f, dtype=float)
    fa = np.abs(f - fs * np.round(f / fs))
    return float(fa) if fa.ndim == 0 else fa


def dominant_frequency(x: np.ndarray, fs: float, pad: int = 8) -> float:
    """Return the frequency [Hz] of the largest peak in the spectrum of x (mean removed, Hann window,
    zero-padded to pad * len(x) points so the peak is located to a fraction of fs / len(x))."""
    x = np.asarray(x, dtype=float)
    x = (x - x.mean()) * np.hanning(x.size)
    nfft = pad * x.size
    X = np.abs(np.fft.rfft(x, nfft))
    return float(np.argmax(X) * fs / nfft)


def sinc_reconstruct(samples: np.ndarray, fs: float, t: np.ndarray, t0: float = 0.0) -> np.ndarray:
    """Whittaker-Shannon reconstruction x(t) = sum_n x[n] sinc(fs (t - t0) - n), samples taken at t0 + n / fs."""
    samples = np.asarray(samples, dtype=float)
    n = np.arange(samples.size)
    arg = fs * (np.asarray(t, dtype=float)[:, None] - t0) - n[None, :]
    return np.sinc(arg) @ samples


# ---------------------------------------------------------------- Butterworth prototype

def butterworth_poles(n: int, wc: float) -> np.ndarray:
    """Return the n left-half-plane poles of an analog Butterworth low-pass of corner wc [rad/s].

    p_k = wc exp(j pi (2k + n - 1) / (2n)), k = 1..n: equally spaced on the circle |s| = wc.
    """
    if n < 1:
        raise ValueError("order must be at least 1")
    k = np.arange(1, n + 1)
    return wc * np.exp(1j * np.pi * (2 * k + n - 1) / (2 * n))


def butterworth_analog(n: int, wc: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (b, a), highest power first, of the unit-DC-gain analog Butterworth low-pass of order n, corner wc [rad/s]."""
    a = poly_from_roots(butterworth_poles(n, wc))
    return np.array([wc ** n]), np.real(a)


def freq_response(b: np.ndarray, a: np.ndarray, w: np.ndarray | float) -> np.ndarray:
    """Return H(jw) = b(jw) / a(jw) for analog polynomial coefficients (highest power first), w in rad/s."""
    s = 1j * np.asarray(w, dtype=float)
    return np.polyval(b, s) / np.polyval(a, s)


def butterworth_order(f_pass: float, droop_db: float, f_stop: float, atten_db: float) -> tuple[int, float, float]:
    """Smallest Butterworth order meeting |H(f_pass)| >= -droop_db dB and |H(f_stop)| <= -atten_db dB.

    Returns (n, fc_min, fc_max) [Hz]: any corner fc in [fc_min, fc_max] meets both specifications;
    fc_min meets the passband exactly, fc_max meets the stopband exactly.
    """
    if not (0 < f_pass < f_stop) or droop_db <= 0 or atten_db <= 0:
        raise ValueError("need 0 < f_pass < f_stop and positive dB specifications")
    ep = 10.0 ** (droop_db / 10.0) - 1.0
    es = 10.0 ** (atten_db / 10.0) - 1.0
    n = int(np.ceil(np.log10(es / ep) / (2.0 * np.log10(f_stop / f_pass)) - 1e-12))
    n = max(n, 1)
    fc_min = f_pass / ep ** (1.0 / (2 * n))
    fc_max = f_stop / es ** (1.0 / (2 * n))
    return n, float(fc_min), float(fc_max)


def antialias_spec(f_band: float, fs: float, a_vib: float, a_allow: float, droop_db: float) -> dict:
    """Design an analog Butterworth anti-aliasing filter for a channel sampled at fs.

    The stopband edge is fs - f_band, the lowest frequency that folds back into [0, f_band]. The
    attenuation needed there is 20 log10(a_vib / a_allow). Returns a dict with keys
    f_stop, atten_db, n, fc_min, fc_max, fc (the geometric mean of the allowed range).
    """
    f_stop = fs - f_band
    if f_stop <= f_band:
        raise ValueError("fs must exceed 2 f_band")
    atten = 20.0 * np.log10(a_vib / a_allow)
    n, lo, hi = butterworth_order(f_band, droop_db, f_stop, atten)
    return {"f_stop": f_stop, "atten_db": float(atten), "n": n, "fc_min": lo, "fc_max": hi, "fc": float(np.sqrt(lo * hi))}


def group_delay_dc(b: np.ndarray, a: np.ndarray) -> float:
    """Return the group delay -d(arg H)/dw at w = 0 [s] of H(s) = b(s)/a(s).

    With b(s) = ... + b1 s + b0 and a(s) = ... + a1 s + a0: tau(0) = a1/a0 - b1/b0.
    """
    b = np.atleast_1d(np.asarray(b, dtype=float))
    a = np.atleast_1d(np.asarray(a, dtype=float))
    b0, a0 = b[-1], a[-1]
    b1 = b[-2] if b.size > 1 else 0.0
    a1 = a[-2] if a.size > 1 else 0.0
    return float(a1 / a0 - b1 / b0)


# ---------------------------------------------------------------- time domain

def simulate_analog(b: np.ndarray, a: np.ndarray, t: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Return the output of the analog filter b(s)/a(s) driven by u on the uniform grid t, from rest.

    The input is held constant between grid points (zero-order hold); each step is exact for that
    held input: x_{k+1} = A_d x_k + B_d u_k, y_k = C x_k + D u_k.
    """
    t = np.asarray(t, dtype=float)
    u = np.asarray(u, dtype=float)
    A, B, C, D = tf_to_companion(b, a)
    n = A.shape[0]
    if n == 0:
        return D * u
    dt = t[1] - t[0]
    Ad, Bd = discretize_zoh(A, B, dt)
    bd = Bd[:, 0]
    X = np.empty((t.size, n))
    x = np.zeros(n)
    for k in range(t.size):
        X[k] = x
        x = Ad @ x + bd * u[k]
    return X @ C[0] + D * u


def sample_through_filter(t_fine: np.ndarray, v: np.ndarray, fs: float, b=None, a=None) -> tuple[np.ndarray, np.ndarray]:
    """Filter v (on the fine grid t_fine) with b/a if given, then sample at fs. Returns (t_samples, samples).

    The fine-grid step must divide 1/fs exactly (to rounding); every (1/fs)/dt-th point is kept.
    """
    dt = t_fine[1] - t_fine[0]
    step = int(round(1.0 / (fs * dt)))
    if abs(step * dt * fs - 1.0) > 1e-9:
        raise ValueError("fine grid must divide the sample period")
    y = v if b is None else simulate_analog(b, a, t_fine, v)
    return t_fine[::step], y[::step]


# ---------------------------------------------------------------- the lesson's numbers

def imu_experiment(orders=(0, 1, 2, 3, 4), fc: float | None = None, seed: int = 9, duration: float = 1.25,
                   f_sim: float = 20000.0, noise_rms: float = 0.05) -> list[dict]:
    """Run the SYNTHETIC gyro record through Butterworth filters of the given orders, sample at FS_IMU and
    measure the 10 Hz alias of the 390 Hz line (after a 0.25 s start-up transient is discarded).

    Returns one dict per order: n, fc, predicted (a_vib |H(f_vib)|), measured, delay_ms, phase_5hz_deg.
    """
    spec = antialias_spec(F_BAND, FS_IMU, A_VIB, A_ALLOW, DROOP_DB)
    fc = spec["fc"] if fc is None else fc
    t = np.arange(int(round(duration * f_sim))) / f_sim
    _, v = synthetic_imu_record(t, seed=seed, noise_rms=noise_rms)
    f_alias = aliased_frequency(F_VIB, FS_IMU)
    rows = []
    for n in orders:
        if n == 0:
            b = a = None
            H = 1.0 + 0j
            H5 = 1.0 + 0j
            tau = 0.0
        else:
            b, a = butterworth_analog(n, 2 * np.pi * fc)
            H = freq_response(b, a, 2 * np.pi * F_VIB)
            H5 = freq_response(b, a, 2 * np.pi * 5.0)
            tau = group_delay_dc(b, a)
        ts, ys = sample_through_filter(t, v, FS_IMU, b, a)
        keep = ts >= 0.25
        rows.append({"n": n, "fc": fc if n else np.inf, "predicted": A_VIB * abs(H),
                     "measured": tone_amplitude(ys[keep], FS_IMU, f_alias),
                     "delay_ms": 1e3 * tau, "phase_5hz_deg": float(np.degrees(np.angle(H5)))})
    return rows


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's tables and save sampling_aliasing.png in out_dir."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"aliases at fs = {FS_IMU:.0f} Hz (Nyquist {FS_IMU / 2:.0f} Hz):")
    for f in (10.0, 190.0, 210.0, 390.0, 410.0, 780.0, 1190.0):
        print(f"  {f:7.1f} Hz -> {aliased_frequency(f, FS_IMU):6.1f} Hz")
    N = 400
    n = np.arange(N)
    x = np.cos(2 * np.pi * F_VIB * n / FS_IMU)
    print(f"FFT peak of a {F_VIB:.0f} Hz tone sampled at {FS_IMU:.0f} Hz for 1 s: {dominant_frequency(x, FS_IMU):.2f} Hz")

    spec = antialias_spec(F_BAND, FS_IMU, A_VIB, A_ALLOW, DROOP_DB)
    print(f"\nspec: band 0-{F_BAND:.0f} Hz, droop <= {DROOP_DB} dB; stop edge fs - f_band = {spec['f_stop']:.0f} Hz, "
          f"attenuation >= 20 log10({A_VIB}/{A_ALLOW}) = {spec['atten_db']:.2f} dB")
    for nn in (1, 2, 3, 4):
        ep = 10 ** (DROOP_DB / 10) - 1
        es = 10 ** (spec['atten_db'] / 10) - 1
        lo, hi = F_BAND / ep ** (1 / (2 * nn)), spec["f_stop"] / es ** (1 / (2 * nn))
        print(f"  order {nn}: passband needs fc >= {lo:6.2f} Hz, stopband needs fc <= {hi:6.2f} Hz  -> {'OK' if lo <= hi else 'impossible'}")
    print(f"design: order {spec['n']}, fc in [{spec['fc_min']:.2f}, {spec['fc_max']:.2f}] Hz, chosen fc = {spec['fc']:.2f} Hz")
    bb, aa = butterworth_analog(spec["n"], 2 * np.pi * spec["fc"])
    print("  poles [rad/s]:", np.array2string(butterworth_poles(spec["n"], 2 * np.pi * spec["fc"]), precision=1))
    print("  a(s) =", np.array2string(aa, precision=6))
    for f in (F_BAND, spec["f_stop"], F_VIB):
        H = freq_response(bb, aa, 2 * np.pi * f)
        print(f"  |H({f:.0f} Hz)| = {abs(H):.5f} ({20 * np.log10(abs(H)):7.2f} dB), phase {np.degrees(np.angle(H)):7.2f} deg")
    print(f"  group delay at DC {1e3 * group_delay_dc(bb, aa):.2f} ms; phase at 5 Hz "
          f"{np.degrees(np.angle(freq_response(bb, aa, 2 * np.pi * 5.0))):.2f} deg")

    print("\nSYNTHETIC gyro record (not a measurement): 390 Hz line 5 deg/s, 780 Hz line 2 deg/s, seeded noise; "
          f"fc = {spec['fc']:.2f} Hz")
    rows = imu_experiment()
    for r in rows:
        print(f"  order {r['n']}: 10 Hz alias predicted {r['predicted']:.4f} deg/s, measured {r['measured']:.4f} deg/s, "
              f"DC delay {r['delay_ms']:.2f} ms, phase at 5 Hz {r['phase_5hz_deg']:.2f} deg")

    # sinc reconstruction of a band-limited signal from its samples
    fs = 50.0
    ns = np.arange(200)
    xs = np.sin(2 * np.pi * 3.0 * ns / fs) + 0.5 * np.cos(2 * np.pi * 11.0 * ns / fs)
    tt = np.linspace(1.0, 3.0, 801)
    xr = sinc_reconstruct(xs, fs, tt)
    xtrue = np.sin(2 * np.pi * 3.0 * tt) + 0.5 * np.cos(2 * np.pi * 11.0 * tt)
    print(f"\nsinc reconstruction (3 Hz + 11 Hz, fs = 50 Hz, 200 samples), max error on 1-3 s: {np.max(np.abs(xr - xtrue)):.4f}")
    xs2 = np.cos(2 * np.pi * 40.0 * ns / fs)
    xr2 = sinc_reconstruct(xs2, fs, tt)
    print(f"a 40 Hz cosine sampled at 50 Hz reconstructs as {dominant_frequency(xr2, 801 / 2.0):.2f} Hz "
          f"(alias {aliased_frequency(40.0, fs):.0f} Hz); max |x_r - cos(2 pi 10 t)| = "
          f"{np.max(np.abs(xr2 - np.cos(2 * np.pi * 10.0 * tt))):.4f}")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    ff = np.logspace(0, 3.3, 600)
    for nn in (1, 2, 3, 4):
        b_, a_ = butterworth_analog(nn, 2 * np.pi * spec["fc"])
        axes[0].semilogx(ff, 20 * np.log10(np.abs(freq_response(b_, a_, 2 * np.pi * ff))), label=f"order {nn}")
    axes[0].axvline(F_BAND, color="g", ls=":", lw=1)
    axes[0].axvline(spec["f_stop"], color="r", ls=":", lw=1)
    axes[0].axhline(-spec["atten_db"], color="r", ls="--", lw=0.8)
    axes[0].axvline(FS_IMU / 2, color="k", ls="-.", lw=0.8, label="Nyquist")
    axes[0].set_ylim(-100, 5); axes[0].set_xlabel("f [Hz]"); axes[0].set_ylabel("|H| [dB]")
    axes[0].set_title(f"Butterworth, fc = {spec['fc']:.1f} Hz", fontsize=9); axes[0].legend(fontsize=7); axes[0].grid(True, which="both", lw=0.3)
    t = np.arange(int(1.25 * 20000)) / 20000.0
    true, v = synthetic_imu_record(t)
    ts, y0 = sample_through_filter(t, v, FS_IMU)
    _, y3 = sample_through_filter(t, v, FS_IMU, bb, aa)
    axes[1].plot(ts, y0, ".", ms=2, label="sampled, no filter")
    axes[1].plot(ts, y3, "-", lw=1, label=f"sampled after order {spec['n']}")
    axes[1].plot(t[::50], true[::50], "k--", lw=0.8, label="true rate")
    axes[1].set_xlim(0.25, 1.25); axes[1].set_ylim(-12, 12)
    axes[1].set_xlabel("t [s]"); axes[1].set_ylabel("rate [deg/s]"); axes[1].legend(fontsize=7)
    axes[1].set_title("SYNTHETIC gyro record (not measured)", fontsize=9)
    keep = ts >= 0.25
    for y, lab in ((y0, "no filter"), (y3, f"order {spec['n']}")):
        Y = np.abs(np.fft.rfft(y[keep] - y[keep].mean())) * 2 / keep.sum()
        axes[2].semilogy(np.fft.rfftfreq(keep.sum(), 1 / FS_IMU), Y + 1e-6, label=lab)
    axes[2].axvline(aliased_frequency(F_VIB, FS_IMU), color="r", ls=":", lw=1)
    axes[2].set_xlabel("f [Hz] (after sampling)"); axes[2].set_ylabel("amplitude [deg/s]"); axes[2].legend(fontsize=7)
    axes[2].set_title("spectrum of the samples (SYNTHETIC)", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "sampling_aliasing.png", dpi=130)
    plt.close(fig)
    print(f"saved {out_dir / 'sampling_aliasing.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
