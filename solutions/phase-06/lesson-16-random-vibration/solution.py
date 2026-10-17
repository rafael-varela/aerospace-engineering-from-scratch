"""P06.16 Random Vibration: reference solution.

How hard does a launch shake an avionics box? A random-vibration test specification is a one-sided
acceleration power spectral density (PSD) G(f) in g^2/Hz, given as breakpoints joined by straight lines
on log-log axes (constant dB/octave slopes). This module

* evaluates such a spectrum and integrates it exactly, segment by segment, to the overall Grms;
* computes the response of a single-degree-of-freedom (SDOF) component on a base that moves with that
  acceleration: Miles' equation, the full integral of |H(f)|^2 G(f), the relative-displacement RMS and
  the zero-upcrossing rate from spectral moments;
* synthesizes a seeded Gaussian random acceleration record with the target PSD and drives an exactly
  discretized SDOF model with it, so the analytic Grms can be checked against a time-domain simulation;
* gives the Gaussian exceedance fraction for a k-sigma level and Davenport's expected peak factor,
  which says how far beyond "3 sigma" the largest peak of a long test is likely to go.

Conventions: f in Hz, G one-sided in g^2/Hz (so the variance is the integral of G over 0..inf),
accelerations in g unless a name says m/s^2, displacements in m. All random records are SYNTHETIC.

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes random_vibration.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-06/lesson-16-random-vibration
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.signal import lfilter
from scipy.special import erfc

G0 = 9.80665                                  # m/s^2 per g

# The component minimum workmanship random vibration level of NASA-STD-7001C, Table 2 (breakpoints in Hz,
# g^2/Hz): 0.01 at 20 Hz, +3 dB/oct to 80 Hz, 0.04 to 500 Hz, -3 dB/oct to 0.01 at 2000 Hz; 6.8 grms overall.
WORKMANSHIP_F = np.array([20.0, 80.0, 500.0, 2000.0])
WORKMANSHIP_G = np.array([0.01, 0.04, 0.04, 0.01])


# ---------------------------------------------------------------- the input spectrum

def segment_slopes_db_per_octave(f_bp: np.ndarray, G_bp: np.ndarray) -> np.ndarray:
    """Slope of each log-log segment in dB/octave: 10 log10(G2/G1) / log2(f2/f1)."""
    f_bp, G_bp = np.asarray(f_bp, float), np.asarray(G_bp, float)
    return 10.0 * np.log10(G_bp[1:] / G_bp[:-1]) / np.log2(f_bp[1:] / f_bp[:-1])


def psd_loglog(f: np.ndarray, f_bp: np.ndarray, G_bp: np.ndarray) -> np.ndarray:
    """Evaluate the breakpoint spectrum at frequencies ``f`` [Hz]: straight lines between breakpoints on
    log-log axes (G = G1 (f/f1)^n on each segment), and zero outside [f_bp[0], f_bp[-1]]."""
    f = np.asarray(f, float)
    f_bp, G_bp = np.asarray(f_bp, float), np.asarray(G_bp, float)
    out = np.zeros_like(f)
    inside = (f >= f_bp[0]) & (f <= f_bp[-1])
    out[inside] = np.exp(np.interp(np.log(f[inside]), np.log(f_bp), np.log(G_bp)))
    return out


def psd_area(f_bp: np.ndarray, G_bp: np.ndarray) -> float:
    """Exact area under the breakpoint spectrum (the mean square, g^2). On a segment G = G1 (f/f1)^n,
    the integral is G1 f1 [(f2/f1)^(n+1) - 1]/(n + 1), or G1 f1 ln(f2/f1) when n = -1 (-3 dB/oct)."""
    f_bp, G_bp = np.asarray(f_bp, float), np.asarray(G_bp, float)
    area = 0.0
    for f1, f2, G1, G2 in zip(f_bp[:-1], f_bp[1:], G_bp[:-1], G_bp[1:]):
        n = np.log(G2 / G1) / np.log(f2 / f1)
        if abs(n + 1.0) < 1e-12:
            area += G1 * f1 * np.log(f2 / f1)
        else:
            area += G1 * f1 * ((f2 / f1) ** (n + 1.0) - 1.0) / (n + 1.0)
    return float(area)


def grms(f_bp: np.ndarray, G_bp: np.ndarray) -> float:
    """Overall RMS acceleration [g] of the breakpoint spectrum: sqrt(psd_area)."""
    return float(np.sqrt(psd_area(f_bp, G_bp)))


# ---------------------------------------------------------------- SDOF response

def transmissibility_sq(f: np.ndarray, fn: float, zeta: float) -> np.ndarray:
    """|H(f)|^2 from base acceleration to the absolute acceleration of an SDOF mass on a spring and
    viscous damper: (1 + (2 zeta r)^2) / ((1 - r^2)^2 + (2 zeta r)^2), r = f/fn."""
    r = np.asarray(f, float) / fn
    return (1.0 + (2 * zeta * r) ** 2) / ((1.0 - r * r) ** 2 + (2 * zeta * r) ** 2)


def miles_grms(fn: float, zeta: float, G_at_fn: float) -> float:
    """Miles' equation: response Grms = sqrt((pi/2) fn Q G(fn)), Q = 1/(2 zeta). Assumes the input PSD is
    flat near fn and the response is dominated by the resonance (light damping)."""
    Q = 1.0 / (2.0 * zeta)
    return float(np.sqrt(0.5 * np.pi * fn * Q * G_at_fn))


def _log_grid(f_bp: np.ndarray, fn: float, zeta: float, n_per_decade: int) -> np.ndarray:
    """Integration grid: log-spaced over the spectrum, with extra points packed around the resonance."""
    lo, hi = float(f_bp[0]), float(f_bp[-1])
    base = np.geomspace(lo, hi, int(np.ceil(np.log10(hi / lo) * n_per_decade)) + 1)
    bw = max(zeta, 1e-4) * fn
    dense = fn + bw * np.sinh(np.linspace(-12.0, 12.0, 4001)) / np.sinh(12.0) * min(50.0, fn / bw * 0.99)
    pts = np.concatenate([base, dense[(dense > lo) & (dense < hi)], np.asarray(f_bp, float)])
    return np.unique(pts)


def response_grms(f_bp: np.ndarray, G_bp: np.ndarray, fn: float, zeta: float,
                  n_per_decade: int = 4000) -> float:
    """Response Grms = sqrt( integral |H(f)|^2 G(f) df ) over the spectrum's band, by the trapezoid rule
    on a log grid refined around fn (converged to about 1e-6 relative at the default grid)."""
    f = _log_grid(np.asarray(f_bp, float), fn, zeta, n_per_decade)
    integrand = transmissibility_sq(f, fn, zeta) * psd_loglog(f, f_bp, G_bp)
    return float(np.sqrt(np.trapezoid(integrand, f)))


def relative_displacement_rms(f_bp: np.ndarray, G_bp: np.ndarray, fn: float, zeta: float,
                              n_per_decade: int = 4000) -> float:
    """RMS of the relative displacement z = x - x_base [m] (spring deflection) for a base acceleration PSD
    G in g^2/Hz: |Z/A_base|^2 = 1 / ((omega_n^2 - omega^2)^2 + (2 zeta omega_n omega)^2)."""
    f = _log_grid(np.asarray(f_bp, float), fn, zeta, n_per_decade)
    wn, w = 2 * np.pi * fn, 2 * np.pi * f
    Hz2 = 1.0 / ((wn ** 2 - w ** 2) ** 2 + (2 * zeta * wn * w) ** 2)
    return float(G0 * np.sqrt(np.trapezoid(Hz2 * psd_loglog(f, f_bp, G_bp), f)))


def spectral_moments(f_bp: np.ndarray, G_bp: np.ndarray, fn: float, zeta: float,
                     n_per_decade: int = 4000) -> tuple[float, float]:
    """Moments m0 = integral S df and m2 = integral f^2 S df of the response acceleration PSD
    S = |H|^2 G (units g^2 and g^2 Hz^2)."""
    f = _log_grid(np.asarray(f_bp, float), fn, zeta, n_per_decade)
    S = transmissibility_sq(f, fn, zeta) * psd_loglog(f, f_bp, G_bp)
    return float(np.trapezoid(S, f)), float(np.trapezoid(f * f * S, f))


def upcrossing_rate(f_bp: np.ndarray, G_bp: np.ndarray, fn: float, zeta: float) -> float:
    """Rice's zero-upcrossing rate nu0+ = sqrt(m2/m0) [crossings per second] of the response."""
    m0, m2 = spectral_moments(f_bp, G_bp, fn, zeta)
    return float(np.sqrt(m2 / m0))


def gaussian_exceedance(k: float) -> float:
    """Fraction of time a zero-mean Gaussian signal spends beyond +/- k sigma: erfc(k / sqrt 2)."""
    return float(erfc(k / np.sqrt(2.0)))


def expected_peak_factor(nu0: float, duration: float) -> float:
    """Davenport's expected largest peak, in units of sigma, of a stationary Gaussian narrow-band signal
    with upcrossing rate nu0 over ``duration`` seconds: sqrt(2 ln N) + 0.5772 / sqrt(2 ln N), N = nu0 T."""
    N = nu0 * duration
    a = np.sqrt(2.0 * np.log(N))
    return float(a + 0.5772 / a)


# ---------------------------------------------------------------- time domain

def synthesize_random(f_bp: np.ndarray, G_bp: np.ndarray, fs: float, duration: float,
                      seed: int) -> tuple[np.ndarray, np.ndarray]:
    """SYNTHETIC stationary Gaussian record with one-sided PSD G (g^2/Hz): sum of cosines on the FFT grid
    df = 1/duration with independent normal in-phase and quadrature amplitudes of variance G(f_k) df,
    built with one inverse real FFT. Returns (t [s], a [g])."""
    rng = np.random.default_rng(seed)
    n = int(round(fs * duration))
    df = fs / n
    fk = np.arange(n // 2 + 1) * df
    sd = np.sqrt(psd_loglog(fk, f_bp, G_bp) * df)
    X = (rng.standard_normal(fk.size) + 1j * rng.standard_normal(fk.size)) * sd
    X[0] = 0.0
    if n % 2 == 0:
        X[-1] = 0.0
    # x[m] = sum_k Re(X_k e^{i 2 pi k m / n}) has variance sum_k G(f_k) df; irfft returns
    # (1/n) [X_0 + 2 sum_k Re(X_k e^{i 2 pi k m / n})], so multiply by n/2
    a = np.fft.irfft(X, n) * (n / 2)
    return np.arange(n) / fs, a


def sdof_discrete(fn: float, zeta: float, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """Exact discretization of the base-excited SDOF for a base acceleration that varies linearly between
    samples (first-order hold). State x = [z, z'] with z the relative displacement:
    z'' + 2 zeta wn z' + wn^2 z = -a_base. Returns IIR coefficients (b, a) for ``lfilter`` that map the
    base acceleration samples to the ABSOLUTE acceleration a_abs = -(2 zeta wn z' + wn^2 z)."""
    wn = 2 * np.pi * fn
    A = np.array([[0.0, 1.0], [-wn * wn, -2 * zeta * wn]])
    B = np.array([[0.0], [-1.0]])
    # first-order hold via the augmented exponential [[A, B, 0], [0, 0, I/dt], [0, 0, 0]]
    M = np.zeros((4, 4))
    M[:2, :2], M[:2, 2:3], M[2, 3] = A * dt, B * dt, 1.0
    E = expm(M)
    Phi, G1, G2 = E[:2, :2], E[:2, 2:3], E[:2, 3:4]
    # x[k+1] = Phi x[k] + (G1 - G2) u[k] + G2 u[k+1]; with u[k+1] the new input this is a state model
    # with direct feed-through. Transfer function from u to y = C x:
    Gm = G1 - G2
    C = np.array([[-wn * wn, -2 * zeta * wn]])
    # y(z) = C (zI - Phi)^-1 (Gm + G2 z) u(z)
    tr, det = np.trace(Phi), np.linalg.det(Phi)
    den = np.array([1.0, -tr, det])
    adj0 = np.array([[1.0, 0.0], [0.0, 1.0]])            # coefficient of z in adj(zI - Phi)
    adj1 = np.array([[-Phi[1, 1], Phi[0, 1]], [Phi[1, 0], -Phi[0, 0]]])   # constant term
    # numerator (in powers z^2, z, 1): C [adj0 z + adj1] [G2 z + Gm]
    n2 = (C @ adj0 @ G2).item()
    n1 = (C @ adj0 @ Gm + C @ adj1 @ G2).item()
    n0 = (C @ adj1 @ Gm).item()
    return np.array([n2, n1, n0]), den


def simulate_sdof_base(a_base: np.ndarray, dt: float, fn: float, zeta: float) -> np.ndarray:
    """Absolute acceleration of the SDOF mass (same units as ``a_base``) for a sampled base acceleration,
    starting from rest, with the exact first-order-hold discretization of ``sdof_discrete``."""
    b, a = sdof_discrete(fn, zeta, dt)
    return lfilter(b, a, np.asarray(a_base, float))


def variance_estimate_relative_sd(fn: float, zeta: float, duration: float) -> float:
    """Approximate relative standard deviation of the sample variance of a lightly damped SDOF response
    of a broadband input over ``duration``: 1/sqrt(zeta wn T) (from Var = (2/T) integral R^2 dtau with
    R(tau) ~ sigma^2 exp(-zeta wn |tau|) cos(wn tau))."""
    return float(1.0 / np.sqrt(zeta * 2 * np.pi * fn * duration))


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save random_vibration.png (SYNTHETIC time history)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.signal import welch

    out_dir = Path(out_dir)
    f_bp, G_bp = WORKMANSHIP_F, WORKMANSHIP_G
    print("input spectrum: NASA-STD-7001C Table 2 component minimum workmanship level:")
    print(f"  breakpoints {f_bp.tolist()} Hz; {G_bp.tolist()} g^2/Hz; slopes "
          f"{np.round(segment_slopes_db_per_octave(f_bp, G_bp), 2).tolist()} dB/oct")
    print(f"  overall {grms(f_bp, G_bp):.3f} Grms (area {psd_area(f_bp, G_bp):.3f} g^2)")

    fn, zeta, m_box = 150.0, 0.05, 2.0
    Gfn = float(psd_loglog(np.array([fn]), f_bp, G_bp)[0])
    g_miles = miles_grms(fn, zeta, Gfn)
    g_full = response_grms(f_bp, G_bp, fn, zeta)
    z_rms = relative_displacement_rms(f_bp, G_bp, fn, zeta)
    nu0 = upcrossing_rate(f_bp, G_bp, fn, zeta)
    print(f"\navionics box: fn = {fn} Hz, zeta = {zeta} (Q = {1 / (2 * zeta):.0f}), mass {m_box} kg")
    print(f"  Miles: {g_miles:.3f} Grms; full integral: {g_full:.3f} Grms (Miles/full = {g_miles / g_full:.4f})")
    print(f"  3-sigma acceleration {3 * g_full:.2f} g; 3-sigma mount force {3 * g_full * m_box * G0:.0f} N")
    print(f"  relative displacement {z_rms * 1e3:.4f} mm rms, 3-sigma {3 * z_rms * 1e3:.3f} mm "
          f"(Miles-type estimate g Grms/wn^2 = {G0 * g_miles / (2 * np.pi * fn) ** 2 * 1e3:.4f} mm)")
    T_test = 60.0
    pf = expected_peak_factor(nu0, T_test)
    print(f"  upcrossing rate {nu0:.1f} /s; expected largest peak in a {T_test:.0f} s test: {pf:.2f} sigma "
          f"= {pf * g_full:.1f} g; fraction of time beyond 3 sigma {gaussian_exceedance(3.0) * 100:.3f} %")

    # time-domain check (SYNTHETIC)
    fs, dur, seed = 8192.0, 32.0, 16
    t, a_in = synthesize_random(f_bp, G_bp, fs, dur, seed)
    a_out = simulate_sdof_base(a_in, 1 / fs, fn, zeta)
    skip = int(0.5 * fs)                                   # drop the start-up transient (about 0.5 s >> 1/(zeta wn))
    rms_in, rms_out = np.std(a_in), np.sqrt(np.mean(a_out[skip:] ** 2))
    u_rel = variance_estimate_relative_sd(fn, zeta, dur - 0.5) / 2
    print(f"\nSYNTHETIC time-domain check (seed {seed}, {dur:.0f} s at {fs:.0f} Hz):")
    print(f"  input {rms_in:.3f} Grms (target {grms(f_bp, G_bp):.3f}); response {rms_out:.3f} Grms "
          f"(integral {g_full:.3f}; expected scatter +/- {u_rel * 100:.1f} % at 1 sigma)")
    frac3 = np.mean(np.abs(a_out[skip:]) > 3 * rms_out)
    print(f"  fraction beyond 3 sigma {frac3 * 100:.3f} % (Gaussian {gaussian_exceedance(3) * 100:.3f} %); "
          f"largest |peak| {np.max(np.abs(a_out[skip:])) / rms_out:.2f} sigma "
          f"(Davenport for {dur - 0.5:.1f} s: {expected_peak_factor(nu0, dur - 0.5):.2f})")

    # Miles versus full integral over fn and zeta
    print("\nMiles / full integral:")
    print("   fn [Hz]   zeta=0.01   zeta=0.05   zeta=0.10")
    for fni in (30.0, 60.0, 150.0, 400.0, 1000.0, 1800.0):
        G_i = float(psd_loglog(np.array([fni]), f_bp, G_bp)[0])
        row = [miles_grms(fni, z, G_i) / response_grms(f_bp, G_bp, fni, z) for z in (0.01, 0.05, 0.10)]
        print(f"   {fni:7.0f}   " + "   ".join(f"{v:9.4f}" for v in row))

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    f = np.geomspace(10, 3000, 2000)
    ax[0].loglog(f, np.where(psd_loglog(f, f_bp, G_bp) > 0, psd_loglog(f, f_bp, G_bp), np.nan), label="input G(f)")
    ax[0].loglog(f, transmissibility_sq(f, fn, zeta) * psd_loglog(f, f_bp, G_bp), label="response |H|^2 G")
    fw, Pin = welch(a_in, fs=fs, nperseg=8192)
    fw2, Pout = welch(a_out[skip:], fs=fs, nperseg=8192)
    ax[0].loglog(fw[1:], Pin[1:], lw=0.5, alpha=0.6, label="Welch, SYNTHETIC input")
    ax[0].loglog(fw2[1:], Pout[1:], lw=0.5, alpha=0.6, label="Welch, SYNTHETIC response")
    ax[0].set_xlim(10, 3000)
    ax[0].set_ylim(1e-5, 10)
    ax[0].set_xlabel("frequency (Hz)")
    ax[0].set_ylabel("PSD (g^2/Hz)")
    ax[0].set_title("input and response PSD")
    ax[0].legend(fontsize=7)
    sl = slice(skip, skip + int(0.25 * fs))
    ax[1].plot(t[sl], a_in[sl], lw=0.5, label="base (input)")
    ax[1].plot(t[sl], a_out[sl], lw=0.8, label="box (response)")
    ax[1].axhline(3 * g_full, color="k", ls="--", lw=0.8)
    ax[1].axhline(-3 * g_full, color="k", ls="--", lw=0.8, label="+/- 3 sigma (analytic)")
    ax[1].set_xlabel("time (s)")
    ax[1].set_ylabel("acceleration (g)")
    ax[1].set_title(f"SYNTHETIC record, seed {seed}")
    ax[1].legend(fontsize=7)
    zs = np.geomspace(0.005, 0.2, 30)
    ax[2].semilogx(zs, [miles_grms(fn, z, Gfn) / response_grms(f_bp, G_bp, fn, z) for z in zs], "o-", ms=3)
    ax[2].set_xlabel("damping ratio zeta")
    ax[2].set_ylabel("Miles / full integral")
    ax[2].set_title(f"Miles' accuracy at fn = {fn:.0f} Hz")
    fig.tight_layout()
    fig.savefig(out_dir / "random_vibration.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'random_vibration.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
