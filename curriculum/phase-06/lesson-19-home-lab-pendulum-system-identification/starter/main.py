"""P06.19 Home Lab: Pendulum System Identification: starter code.

Complete every function that raises NotImplementedError. The SYNTHETIC data tools (simulate_swing,
make_synthetic_swing, write_synthetic_csv, read_swing_csv), initial_guess and main() are given; the
generator runs once pendulum_rhs is written. When finished, running this file analyses the committed
SYNTHETIC record (or your own CSV) and writes the figure.

Identify the natural frequency and damping of a real string pendulum from a logged swing (a phone video
tracked frame by frame, or a gyroscope), with uncertainties, and decide what kind of damping it has.

The module covers four jobs.

* Prediction (before the swing): the small-angle frequency from L and g with its uncertainty, the
  large-amplitude correction from the exact (elliptic-integral) period, and an order-of-magnitude
  damping ratio from air drag on the bob.
* Model and SYNTHETIC data: the nonlinear pendulum with viscous, quadratic (air drag) and Coulomb
  (pivot friction) damping, integrated with ``solve_ivp``; a seeded SYNTHETIC swing record with sampling
  jitter and angle noise, for learners without the kit and for testing the analysis.
* Analysis: damped frequency from zero crossings and from an FFT peak; peaks with sub-sample
  refinement; damping by the logarithmic decrement over n cycles and by a straight-line fit of ln(peak)
  against cycle number; both together by nonlinear least squares of
  theta(t) = A exp(-zeta wn t) cos(wd t + phi) + c with ``scipy.optimize.curve_fit``, uncertainties from
  the covariance, checked by a seeded residual (block) bootstrap.
* Diagnosis: the decay envelope fitted with dA/dt = -kappa A^p (p = 0 Coulomb, 1 viscous, 2 quadratic)
  and the frequency measured in windows against amplitude, compared with the exact amplitude correction.

Units: SI; angles in radians in code (the CSV stores degrees, as a video tool would).

Run:  python main.py [out_dir] [csv]   (prints the numbers for the SYNTHETIC record or your CSV, writes
                                        pendulum_sysid.png)
Test: pytest curriculum/phase-06/lesson-19-home-lab-pendulum-system-identification
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
from scipy.special import ellipk

G_STANDARD = 9.80665
SYNTH_SEED = 619
SYNTH_DEFAULTS = dict(L=0.80, g=9.81, theta0_deg=20.0, b1=0.003, b2=0.004, b0=0.0, fs=30.0,
                      duration=120.0, jitter_s=0.002, noise_deg=0.3, offset_deg=0.5)


# ---------------------------------------------------------------- prediction

def amplitude_frequency_ratio(theta_amp: np.ndarray | float) -> np.ndarray | float:
    """omega(A) / omega_0 for an undamped pendulum of amplitude A [rad]: pi / (2 K(m)), m = sin^2(A/2)
    (scipy's ellipk takes the parameter m = k^2). Tends to 1 - A^2/16 for small A."""
    raise NotImplementedError


def predict_frequency(L: float, u_L: float, g: float, u_g: float, theta_amp: float = 0.0) -> dict[str, float]:
    """Predicted small-angle frequency f0 = sqrt(g/L)/(2 pi) [Hz] with its standard uncertainty
    (relative u = 0.5 sqrt((u_g/g)^2 + (u_L/L)^2)), and the frequency at amplitude ``theta_amp`` [rad]."""
    raise NotImplementedError


def quadratic_drag_coefficient(rho: float, Cd: float, d_bob: float, m: float, L: float) -> float:
    """b2 [1/rad] in theta'' = ... - b2 theta'|theta'| for air drag on a spherical bob of diameter d_bob:
    the drag 0.5 rho Cd (pi d^2/4) (L theta')^2 times L, over the inertia m L^2."""
    raise NotImplementedError


def envelope_rate(A: np.ndarray | float, omega: float, b1: float = 0.0, b2: float = 0.0,
                  b0: float = 0.0) -> np.ndarray | float:
    """Cycle-averaged dA/dt [rad/s] of a lightly damped pendulum (energy balance over one cycle):
    -(b1/2) A - (4 b2 omega / (3 pi)) A^2 - 2 b0 / (pi omega), for viscous (b1, 1/s), quadratic
    (b2, 1/rad) and Coulomb (b0, rad/s^2) damping."""
    raise NotImplementedError


def equivalent_zeta(A: np.ndarray | float, omega: float, b1: float = 0.0, b2: float = 0.0,
                    b0: float = 0.0) -> np.ndarray | float:
    """Equivalent viscous damping ratio at amplitude A: -envelope_rate / (omega A)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the model and SYNTHETIC data

def pendulum_rhs(t: float, y: np.ndarray, L: float, g: float, b1: float = 0.0, b2: float = 0.0,
                 b0: float = 0.0, eps: float = 1e-3) -> list[float]:
    """State y = (theta, theta'). theta'' = -(g/L) sin theta - b1 theta' - b2 theta'|theta'|
    - b0 tanh(theta'/eps): viscous, quadratic and (smoothed) Coulomb damping."""
    raise NotImplementedError


def simulate_swing(t: np.ndarray, L: float, g: float, theta0: float, b1: float = 0.0, b2: float = 0.0,
                   b0: float = 0.0, rtol: float = 1e-10, atol: float = 1e-12) -> np.ndarray:
    """theta(t) [rad] released from rest at theta0, by DOP853 with dense output, at the (sorted) times t."""
    sol = solve_ivp(pendulum_rhs, (0.0, float(np.max(t))), [theta0, 0.0], args=(L, g, b1, b2, b0),
                    method="DOP853", rtol=rtol, atol=atol, dense_output=True)
    return sol.sol(np.asarray(t, float))[0]


def make_synthetic_swing(seed: int, L: float = 0.80, g: float = 9.81, theta0_deg: float = 20.0,
                         b1: float = 0.003, b2: float = 0.004, b0: float = 0.0, fs: float = 30.0,
                         duration: float = 120.0, jitter_s: float = 0.002, noise_deg: float = 0.3,
                         offset_deg: float = 0.5) -> dict:
    """SYNTHETIC swing record, as a phone-video tracking tool would export it. The nonlinear pendulum is
    integrated with solve_ivp; samples are taken at frame instants k/fs plus normal jitter (sd jitter_s),
    with the timestamps rounded to 1 ms; the angle gets a constant offset (a tilted protractor) and normal
    noise (sd noise_deg). Returns dict(t [s], theta [rad], truth=dict(...)). NOT a measurement."""
    rng = np.random.default_rng(seed)
    n = int(round(duration * fs))
    t_true = np.arange(n) / fs + rng.normal(0.0, jitter_s, n)
    t_true = np.sort(np.clip(t_true, 0.0, None))
    th = simulate_swing(t_true, L, g, np.radians(theta0_deg), b1, b2, b0)
    theta = th + np.radians(offset_deg) + rng.normal(0.0, np.radians(noise_deg), n)
    truth = dict(seed=seed, L=L, g=g, theta0_deg=theta0_deg, b1=b1, b2=b2, b0=b0, fs=fs, duration=duration,
                 jitter_s=jitter_s, noise_deg=noise_deg, offset_deg=offset_deg, omega0=float(np.sqrt(g / L)))
    return {"t": np.round(t_true, 3), "theta": theta, "truth": truth}


def write_synthetic_csv(path: str | Path, seed: int = SYNTH_SEED, **kw) -> Path:
    """Write make_synthetic_swing(seed, **kw) as CSV (t_s, theta_deg) with a header that says SYNTHETIC."""
    d = make_synthetic_swing(seed, **kw)
    tr = d["truth"]
    path = Path(path)
    head = [f"# SYNTHETIC data from make_synthetic_swing(seed={seed}) in P06.19 solution.py: NOT a measurement",
            f"# model: nonlinear pendulum, L = {tr['L']} m, g = {tr['g']} m/s^2, released from rest at "
            f"{tr['theta0_deg']} deg; damping b1 = {tr['b1']} 1/s (viscous), b2 = {tr['b2']} 1/rad (quadratic), "
            f"b0 = {tr['b0']} rad/s^2 (Coulomb)",
            f"# sampling: {tr['fs']} frames/s for {tr['duration']} s, timing jitter sd {tr['jitter_s'] * 1e3:.1f} ms "
            f"(timestamps rounded to 1 ms); angle noise sd {tr['noise_deg']} deg; angle offset {tr['offset_deg']} deg",
            "t_s,theta_deg"]
    rows = [f"{t:.3f},{np.degrees(th):.3f}" for t, th in zip(d["t"], d["theta"])]
    path.write_text("\n".join(head + rows) + "\n")
    return path


def read_swing_csv(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Read a swing CSV with columns t_s, theta_deg (lines starting with # are comments).
    Returns (t [s], theta [rad])."""
    rows = [ln for ln in Path(path).read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    data = np.array([[float(v) for v in ln.split(",")[:2]] for ln in rows[1:]])   # rows[0] is the header
    return data[:, 0], np.radians(data[:, 1])


# ---------------------------------------------------------------- frequency

def zero_crossings(t: np.ndarray, y: np.ndarray, direction: int = 1) -> np.ndarray:
    """Times where y (mean removed) crosses zero upward (direction = 1) or downward (-1), by linear
    interpolation between the samples on either side."""
    raise NotImplementedError


def period_from_crossings(tc: np.ndarray) -> tuple[float, float]:
    """Damped period T_d and its standard uncertainty from a straight-line fit of crossing time against
    crossing number (the slope); the residual scatter sets the uncertainty."""
    raise NotImplementedError


def fft_peak_frequency(t: np.ndarray, y: np.ndarray, pad: int = 16) -> float:
    """Frequency [Hz] of the largest spectral peak: resample onto a uniform grid at the mean sample rate
    (np.interp), remove the mean, apply a Hann window, zero-pad by ``pad``, and refine the peak bin by a
    parabola through the log magnitudes of the three bins around it."""
    raise NotImplementedError


# ---------------------------------------------------------------- peaks and the logarithmic decrement

def find_peaks_interp(t: np.ndarray, y: np.ndarray, half_width: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """Positive peaks of y about its mean, one per cycle: for each lobe between an upward and the next
    downward zero crossing, take the largest sample and fit a parabola by least squares through it and
    ``half_width`` samples on each side (non-uniform times allowed). Returns (peak times, peak heights
    above the mean)."""
    raise NotImplementedError


def log_decrement(peaks: np.ndarray, n: int) -> tuple[float, float]:
    """Logarithmic decrement over n cycles, delta = ln(A_0 / A_n) / n, and the damping ratio
    zeta = delta / sqrt(4 pi^2 + delta^2)."""
    raise NotImplementedError


def log_decrement_fit(peaks: np.ndarray) -> dict[str, float]:
    """Straight-line least-squares fit of ln(A_k) against cycle number k: slope -delta. Returns delta,
    its standard uncertainty (from the residual scatter), zeta and u_zeta (zeta ~ delta/(2 pi))."""
    raise NotImplementedError


# ---------------------------------------------------------------- nonlinear least squares

def damped_cosine(t: np.ndarray, A: float, zeta: float, omega_n: float, phi: float, c: float) -> np.ndarray:
    """A exp(-zeta omega_n t) cos(omega_d t + phi) + c with omega_d = omega_n sqrt(1 - zeta^2)."""
    raise NotImplementedError


def initial_guess(t: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Starting values for the fit: omega from the FFT peak, zeta from the ln(peak) fit, A from the first
    peak, phi from the first peak's time, c from the mean."""
    w = 2 * np.pi * fft_peak_frequency(t, y)
    tp, ap = find_peaks_interp(t, y)
    z = max(log_decrement_fit(ap)["zeta"], 1e-5) if ap.size > 3 else 1e-3
    A0 = ap[0] * np.exp(z * w * tp[0])
    phi = -w * tp[0]
    return np.array([A0, z, w, (phi + np.pi) % (2 * np.pi) - np.pi, float(np.mean(y))])


def fit_damped_cosine(t: np.ndarray, y: np.ndarray, p0: np.ndarray | None = None) -> dict:
    """Fit damped_cosine by nonlinear least squares (scipy.optimize.curve_fit, unweighted). Returns
    p (A, zeta, omega_n, phi, c), the covariance (scaled by the residual variance, absolute_sigma=False),
    the standard errors perr, the residual RMS, and the damped frequency f_d = omega_d / (2 pi)."""
    raise NotImplementedError


def bootstrap_fit(t: np.ndarray, y: np.ndarray, fit: dict, n_boot: int = 200, seed: int = 0,
                  block: int = 1) -> np.ndarray:
    """Residual bootstrap of the fit: refit fitted + resampled residuals ``n_boot`` times. With
    block > 1, residuals are resampled in contiguous blocks of that many samples (a moving-block
    bootstrap), which keeps short-range correlation. Returns the (n_boot, 5) parameter samples."""
    raise NotImplementedError


# ---------------------------------------------------------------- diagnosis: kind of damping, amplitude effect

def envelope_power_law(t: np.ndarray, A0: float, kappa: float, p: float) -> np.ndarray:
    """Solution of dA/dt = -kappa A^p from A(0) = A0: A0 (1 - (1-p) kappa A0^(p-1) t)^(1/(1-p)), the
    exponential A0 exp(-kappa t) at p = 1; clipped at zero."""
    raise NotImplementedError


def fit_envelope(t_peaks: np.ndarray, A_peaks: np.ndarray) -> dict:
    """Fit A(t) = envelope_power_law(t; A0, kappa, p) to the peaks: p near 0 is Coulomb (constant loss
    per cycle), near 1 viscous (exponential), near 2 quadratic drag (1/A grows linearly). Returns p, its
    standard error, A0, kappa and the covariance."""
    raise NotImplementedError


def envelope_viscous_quadratic(t: np.ndarray, A0: float, alpha: float, beta: float) -> np.ndarray:
    """Solution of dA/dt = -alpha A - beta A^2 (viscous plus quadratic damping, cycle-averaged) from
    A(0) = A0: alpha A0 e^(-alpha t) / (alpha + beta A0 (1 - e^(-alpha t))). With envelope_rate,
    alpha = b1/2 and beta = 4 b2 omega / (3 pi)."""
    raise NotImplementedError


def fit_envelope_vq(t_peaks: np.ndarray, A_peaks: np.ndarray, omega: float) -> dict:
    """Fit envelope_viscous_quadratic to the peaks and convert to the damping coefficients
    b1 = 2 alpha [1/s] and b2 = 3 pi beta / (4 omega) [1/rad], with standard errors from the covariance."""
    raise NotImplementedError


def frequency_vs_amplitude(t: np.ndarray, y: np.ndarray, cycles: int = 6) -> tuple[np.ndarray, np.ndarray]:
    """Damped frequency [Hz] in consecutive windows of ``cycles`` cycles (from the upward zero crossings),
    with the mean peak amplitude [rad] in each window. Returns (amplitudes, frequencies)."""
    raise NotImplementedError


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".", csv_path: str | Path | None = None) -> None:
    """Analyse a swing record (the committed SYNTHETIC CSV by default, else regenerate it from the seed),
    print the lesson's numbers and save pendulum_sysid.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    if csv_path is not None and Path(csv_path).exists():
        t, th = read_swing_csv(csv_path)
        label = f"SYNTHETIC CSV {Path(csv_path).name}" if "SYNTHETIC" in Path(csv_path).name else f"record {Path(csv_path).name}"
    else:
        d = make_synthetic_swing(SYNTH_SEED, **SYNTH_DEFAULTS)
        t, th = d["t"], d["theta"]
        label = f"SYNTHETIC seed {SYNTH_SEED}"
    tr = SYNTH_DEFAULTS

    pred = predict_frequency(0.80, 0.003, 9.81, 0.01, np.radians(20.0))
    b2_pred = quadratic_drag_coefficient(1.2, 0.47, 0.05, 0.15, 0.80)
    w0 = pred["omega0"]
    print(f"prediction (L = 0.800 +/- 0.003 m, g = 9.81 +/- 0.01 m/s^2): f0 = {pred['f0']:.4f} +/- {pred['u_f0']:.4f} Hz; "
          f"at 20 deg {pred['f_amp']:.4f} Hz ({(pred['f_amp'] / pred['f0'] - 1) * 100:.2f} %)")
    print(f"  air drag on a 50 mm, 150 g bob (Cd 0.47): b2 = {b2_pred:.4f} /rad, zeta_eq at 20 deg = "
          f"{float(equivalent_zeta(np.radians(20.0), w0, b2=b2_pred)):.2e}, at 5 deg = {float(equivalent_zeta(np.radians(5.0), w0, b2=b2_pred)):.2e}")

    print(f"\ndata: {label}: {t.size} samples over {t[-1] - t[0]:.1f} s, mean interval {np.mean(np.diff(t)) * 1e3:.2f} ms "
          f"(sd {np.std(np.diff(t)) * 1e3:.2f} ms)")
    tc = zero_crossings(t, th)
    Td, uTd = period_from_crossings(tc)
    f_fft = fft_peak_frequency(t, th)
    print(f"damped frequency: zero crossings {1 / Td:.5f} +/- {uTd / Td ** 2:.5f} Hz ({tc.size} crossings); FFT peak {f_fft:.5f} Hz")

    tp, ap = find_peaks_interp(t, th)
    n_cyc = 50
    dl, zl = log_decrement(ap, n_cyc)
    lf = log_decrement_fit(ap)
    print(f"peaks: {ap.size}, first {np.degrees(ap[0]):.2f} deg, last {np.degrees(ap[-1]):.2f} deg")
    print(f"log decrement over {n_cyc} cycles: delta = {dl:.5f}, zeta = {zl:.3e}; "
          f"ln(peak) line over all peaks: zeta = {lf['zeta']:.3e} +/- {lf['u_zeta']:.1e}")

    fit = fit_damped_cosine(t, th)
    A, z, wn, ph, c = fit["p"]
    eA, ez, ew, eph, ec = fit["perr"]
    print(f"NLS fit: A = {np.degrees(A):.3f} +/- {np.degrees(eA):.3f} deg, zeta = {z:.4e} +/- {ez:.1e}, "
          f"omega_n = {wn:.5f} +/- {ew:.5f} rad/s (f_n = {wn / (2 * np.pi):.5f} Hz), c = {np.degrees(c):.3f} +/- {np.degrees(ec):.3f} deg; "
          f"residual RMS {np.degrees(fit['rms']):.3f} deg")
    boot = bootstrap_fit(t, th, fit, n_boot=100, seed=1, block=1)
    bootb = bootstrap_fit(t, th, fit, n_boot=100, seed=2, block=60)
    sd, sdb = boot.std(axis=0, ddof=1), bootb.std(axis=0, ddof=1)
    print(f"  bootstrap sd (100 resamples): zeta {sd[1]:.1e} (blocks of 2 s: {sdb[1]:.1e}), omega_n {sd[2]:.1e} "
          f"(blocks: {sdb[2]:.1e}) rad/s; covariance: {ez:.1e}, {ew:.1e}")

    env = fit_envelope(tp, ap)
    print(f"envelope dA/dt = -kappa A^p: p = {env['p']:.3f} +/- {env['u_p']:.3f} (0 Coulomb, 1 viscous, 2 quadratic)")
    amps, freqs = frequency_vs_amplitude(t, th, cycles=6)
    slope = np.polyfit(amps ** 2, freqs, 1)
    print(f"frequency against amplitude ({amps.size} windows of 6 cycles): f = {slope[1]:.5f} Hz x "
          f"(1 {slope[0] / slope[1]:+.4f} A^2); exact small-amplitude law: 1 - 0.0625 A^2")

    vq = fit_envelope_vq(tp, ap, wn)
    print(f"envelope dA/dt = -alpha A - beta A^2: b1 = {vq['b1']:.2e} +/- {vq['u_b1']:.1e} 1/s, "
          f"b2 = {vq['b2']:.2e} +/- {vq['u_b2']:.1e} 1/rad (correlation {vq['corr']:.2f})")
    ratio_mean = float(np.mean(amplitude_frequency_ratio(ap)))
    Aw = np.column_stack([amps ** 2, np.ones_like(amps)])
    cw, *_ = np.linalg.lstsq(Aw, freqs, rcond=None)
    rw = freqs - Aw @ cw
    u_int = float(np.sqrt((rw @ rw / (amps.size - 2)) * np.linalg.inv(Aw.T @ Aw)[1, 1]))
    print(f"small-angle frequency: NLS f_n / mean amplitude ratio = {wn / (2 * np.pi) / ratio_mean:.5f} Hz; "
          f"windowed intercept {cw[1]:.5f} +/- {u_int:.5f} Hz")

    zeta_true = equivalent_zeta(ap, np.sqrt(tr["g"] / tr["L"]), tr["b1"], tr["b2"], tr["b0"])
    print(f"\nSYNTHETIC truth (generator): omega_0 = {np.sqrt(tr['g'] / tr['L']):.5f} rad/s, f0 = "
          f"{np.sqrt(tr['g'] / tr['L']) / (2 * np.pi):.5f} Hz; zeta_eq from {zeta_true[0]:.3e} at the first peak "
          f"to {zeta_true[-1]:.3e} at the last (mean {np.mean(zeta_true):.3e})")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    ax[0].plot(t, np.degrees(th), ".", ms=1.5, label="data")
    ax[0].plot(t, np.degrees(damped_cosine(t, *fit["p"])), lw=0.6, label="NLS fit")
    ax[0].plot(tp, np.degrees(ap + np.mean(th)), "r.", ms=3, label="peaks")
    ax[0].set_xlim(0, 20)
    ax[0].set_xlabel("time (s)")
    ax[0].set_ylabel("angle (deg)")
    ax[0].set_title(f"{label}: first 20 s")
    ax[0].legend(fontsize=7)
    tt = tp - tp[0]
    ax[1].semilogy(tt, np.degrees(ap), "o", ms=3, label="peaks")
    ax[1].semilogy(tt, np.degrees(envelope_power_law(tt, env["A0"], env["kappa"], env["p"])), label=f"power law, p = {env['p']:.2f}")
    ax[1].semilogy(tt, np.degrees(ap[0] * np.exp(-lf["delta"] * np.arange(ap.size))), "--", label="exponential (log decrement)")
    ax[1].set_xlabel("time from first peak (s)")
    ax[1].set_ylabel("peak amplitude (deg)")
    ax[1].set_title(f"{label}: decay envelope")
    ax[1].legend(fontsize=7)
    Ag = np.linspace(0, amps.max() * 1.05, 50)
    ax[2].plot(np.degrees(amps), freqs, "o", ms=4, label="windowed zero crossings")
    ax[2].plot(np.degrees(Ag), slope[1] * amplitude_frequency_ratio(Ag), label="exact law, fitted f0")
    ax[2].set_xlabel("amplitude (deg)")
    ax[2].set_ylabel("damped frequency (Hz)")
    ax[2].set_title(f"{label}: amplitude-dependent frequency")
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out_dir / "pendulum_sysid.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'pendulum_sysid.png'}")


if __name__ == "__main__":
    # python main.py [out_dir] [csv]: analyse your own CSV (t_s, theta_deg) or, by default, the SYNTHETIC record
    args = sys.argv[1:]
    main(args[0] if args else ".", args[1] if len(args) > 1 else None)
