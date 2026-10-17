"""P02.19 Random Processes & Power Spectra: starter code.

A Dryden gust generator: the one-sided power spectral density (PSD) and autocorrelation of the
longitudinal (u) and vertical (w) gust components, a shaping filter that turns seeded white noise
into a gust time history with that spectrum, from-scratch estimators (biased autocorrelation,
periodogram, Welch average) checked against ``scipy.signal`` as the oracle, and the response of an
aircraft heave mode (a first-order lag) to the gust: RMS load factor, level-crossing rate and the
record length needed to estimate a variance.

Conventions (all defined in the lesson, part 2):
    * omega is angular frequency [rad/s]; f = omega / (2 pi) is frequency [Hz]; V is airspeed [m/s]
      so that spatial frequency Omega = omega / V [rad/m].
    * The two-sided spectrum is S(omega) = integral R(tau) exp(-i omega tau) dtau.
    * The one-sided spectrum used throughout is Phi(omega) = S(omega) / pi for omega >= 0, so that
          sigma**2 = R(0) = integral_0^inf Phi(omega) d omega.
      In hertz, G(f) = 2 pi Phi(2 pi f) and sigma**2 = integral_0^inf G(f) df.
    * Dryden forms, with b = V / L the gust break frequency and x = omega / b:
          Phi_u(omega) = sigma_u**2 (2 L_u / (pi V)) / (1 + x**2)
          Phi_w(omega) = sigma_w**2 (L_w / (pi V)) (1 + 3 x**2) / (1 + x**2)**2
      (the lateral component v has the same form as w).
    * Sampled white noise of variance 1/dt stands for continuous white noise of two-sided PSD 1.

Complete every function marked NotImplementedError.

Run:  python main.py          (prints the gust statistics, writes gust_spectrum.png)
Test: pytest curriculum/phase-02/lesson-19-random-processes-power-spectra
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm
from scipy.special import erfc

G0 = 9.80665  # standard gravity [m/s^2]

# A light aircraft in low-altitude turbulence. Every number is a lesson parameter, not a measurement.
AIRCRAFT = {"mass": 1000.0, "area": 16.0, "cl_alpha": 5.0, "airspeed": 60.0, "rho": 1.19}
TURBULENCE = {"sigma_w": 1.5, "sigma_u": 1.5, "L_w": 300.0, "L_u": 300.0}


# ---------------------------------------------------------------- the Dryden model

def dryden_psd(omega: np.ndarray, sigma: float, L: float, V: float, component: str = "w") -> np.ndarray:
    """Return the one-sided Dryden PSD Phi(omega) [(m/s)^2 per rad/s] for omega >= 0.

    component "u": sigma^2 (2L/(pi V)) / (1 + (L omega/V)^2);
    component "w" (or "v"): sigma^2 (L/(pi V)) (1 + 3 (L omega/V)^2) / (1 + (L omega/V)^2)^2.
    Normalised so that the integral from 0 to infinity is sigma^2. Raise ``ValueError`` for a bad
    component or non-positive sigma, L, V.
    """
    raise NotImplementedError


def dryden_psd_hz(f: np.ndarray, sigma: float, L: float, V: float, component: str = "w") -> np.ndarray:
    """Return the one-sided PSD per hertz, G(f) = 2 pi Phi(2 pi f) [(m/s)^2/Hz], so that integral G df = sigma^2."""
    raise NotImplementedError


def dryden_autocorrelation(tau: np.ndarray, sigma: float, L: float, V: float, component: str = "w") -> np.ndarray:
    """Return the autocorrelation R(tau) [(m/s)^2] of the Dryden component (the inverse transform of Phi).

    "u": sigma^2 exp(-V|tau|/L);   "w"/"v": sigma^2 (1 - V|tau|/(2L)) exp(-V|tau|/L).
    """
    raise NotImplementedError


def psd_from_autocorrelation(omega: float, R: Callable[[float], float], tau_max: float) -> float:
    """Return Phi(omega) = (2/pi) integral_0^tau_max R(tau) cos(omega tau) dtau: the Wiener-Khinchin relation, numerically.

    ``tau_max`` must be several correlation times so that the truncated tail is negligible. Uses
    ``scipy.integrate.quad`` with the oscillatory cosine weight.
    """
    raise NotImplementedError


def integrated_variance(psd: Callable[[float], float]) -> float:
    """Return integral_0^inf psd(omega) d omega, the variance of the process with that one-sided PSD."""
    raise NotImplementedError


# ---------------------------------------------------------------- shaping filter and simulation

def shaping_filter(sigma: float, L: float, V: float, component: str = "w") -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return the continuous state-space (A, B, C, D) of the filter H(s) with |H(i omega)|^2 = pi Phi(omega).

    Driven by white noise of two-sided PSD 1, the output has one-sided PSD Phi and variance sigma^2.
    With b = V/L:
        "u":  H_u(s) = sigma sqrt(2L/V) b / (s + b)                       (one state)
        "w":  H_w(s) = sigma sqrt(L/V) b (sqrt(3) s + b) / (s + b)^2       (two states, controllable form)
    """
    raise NotImplementedError


def discretize_zoh(A: np.ndarray, B: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (Ad, Bd) for x[n+1] = Ad x[n] + Bd u[n] with u held constant over each step (zero-order hold).

    Both come from one matrix exponential of the block matrix [[A, B], [0, 0]] dt (P02.05).
    """
    raise NotImplementedError


def simulate_lti(Ad: np.ndarray, Bd: np.ndarray, C: np.ndarray, D: np.ndarray, u: np.ndarray,
                 x0: np.ndarray | None = None) -> np.ndarray:
    """Return the scalar output y[n] = C x[n] + D u[n] of the discrete system x[n+1] = Ad x[n] + Bd u[n].

    A plain loop written from scratch; ``scipy.signal.dlsim`` is the oracle in the tests.
    """
    raise NotImplementedError


def white_noise(n: int, dt: float, rng: np.random.Generator) -> np.ndarray:
    """Return n samples of Gaussian white noise with variance 1/dt: the sampled stand-in for continuous white noise of two-sided PSD 1."""
    raise NotImplementedError


def generate_gust(duration: float, dt: float, sigma: float, L: float, V: float,
                  component: str = "w", seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Return (t, gust) for a seeded Dryden gust record of the given duration [s] and step [s].

    White noise -> shaping filter (ZOH discretisation) -> gust velocity [m/s]. The first 10 L/V
    seconds are simulated and discarded so the record is stationary from t = 0.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- estimators

def autocorrelation_estimate(x: np.ndarray, max_lag: int) -> np.ndarray:
    """Return the biased autocorrelation estimate R[k] = (1/N) sum_{n} (x[n]-m)(x[n+k]-m) for k = 0..max_lag.

    The 1/N (not 1/(N-k)) normalisation keeps the estimate positive semi-definite, so its transform is
    a valid spectrum; the price is a bias of order k/N at lag k.
    """
    raise NotImplementedError


def periodogram(x: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (f [Hz], P [x^2/Hz]): the one-sided periodogram |X_k|^2 dt / N, doubled except at DC and Nyquist.

    Its mean is (nearly) the true PSD but its standard deviation equals its mean at every bin,
    whatever N: the periodogram is an inconsistent estimator, which is why Welch averages it.
    """
    raise NotImplementedError


def welch_psd(x: np.ndarray, dt: float, nperseg: int, overlap: float = 0.5) -> tuple[np.ndarray, np.ndarray]:
    """Return (f, P): Welch's averaged, Hann-windowed, mean-removed periodogram over segments with the given overlap.

    Each segment is detrended by its own mean, multiplied by the periodic Hann window
    w[n] = 0.5 - 0.5 cos(2 pi n / nperseg), transformed, and scaled by dt / sum(w^2); the one-sided
    spectrum is doubled except at DC and Nyquist, and the K segment spectra are averaged.
    Matches ``scipy.signal.welch(x, fs=1/dt, nperseg=nperseg)`` with its default window and detrend.
    """
    raise NotImplementedError


def band_average(f: np.ndarray, P: np.ndarray, f_lo: float, f_hi: float) -> float:
    """Return the mean of P over the bins with f_lo <= f < f_hi (a smoothed spectral level for comparisons)."""
    raise NotImplementedError


# ---------------------------------------------------------------- aircraft heave response

def heave_time_constant(mass: float, rho: float, V: float, S: float, cl_alpha: float) -> float:
    """Return tau = 2 m / (rho V S C_L_alpha) [s]: the time constant of the plunge-only (heave) mode.

    From m dv/dt = (rho V S C_L_alpha / 2)(w_g - v): the aircraft's vertical velocity v follows the
    gust w_g through a first-order lag, and its acceleration is a = (w_g - v)/tau.
    """
    raise NotImplementedError


def heave_response(omega: np.ndarray, tau: float) -> tuple[np.ndarray, np.ndarray]:
    """Return the complex frequency responses (H_v, H_a) from gust velocity to plunge velocity and to acceleration.

    H_v = 1 / (1 + i omega tau),  H_a = i omega / (1 + i omega tau) = i omega H_v.
    """
    raise NotImplementedError


def heave_response_variances(sigma_w: float, L: float, V: float, tau: float) -> tuple[float, float]:
    """Return (var_v, var_a): closed-form variances of plunge velocity and acceleration for Dryden w input.

    With b = V/L and c = 1/tau, integrating |H|^2 Phi_w over omega by partial fractions gives
        var_v = sigma_w^2 c (b + 2c) / (2 (b + c)^2),
        var_a = sigma_w^2 b c^2 (2b + 3c) / (2 (b + c)^2).
    Limits: c -> inf gives var_v -> sigma_w^2 (the aircraft follows the gust) and var_a -> inf (the
    Dryden gust is not differentiable); c -> 0 gives both -> 0.
    """
    raise NotImplementedError


def response_variance_numeric(sigma_w: float, L: float, V: float, tau: float, which: str = "a") -> float:
    """Return integral_0^inf |H(i omega)|^2 Phi_w(omega) d omega by quadrature, for which = "v" or "a"."""
    raise NotImplementedError


def simulate_heave(t: np.ndarray, w: np.ndarray, tau: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (v, a): plunge velocity and acceleration time histories for the gust record w(t).

    State v with dv/dt = (w - v)/tau, outputs v and a = (w - v)/tau. The gust is taken as *linear*
    between samples (first-order hold), for which the step is exact:
        v[n+1] = e^{-h} v[n] + (1 - e^{-h}) w[n] + (1 - (1 - e^{-h})/h) (w[n+1] - w[n]),   h = dt/tau.
    A zero-order hold would make w jump at every sample and the sampled a would carry the jumps:
    for a gust with a 1/omega^2 spectral tail those jumps scale as sqrt(dt) and bias sigma_a upward
    by several per cent at dt = 0.02 s (try it: replace the last term by zero).
    """
    raise NotImplementedError


# ---------------------------------------------------------------- what the statistics mean for the aircraft

def gaussian_exceedance_fraction(level: float, sigma: float) -> float:
    """Return P(|x| > level) for a zero-mean Gaussian of standard deviation sigma: the fraction of time spent beyond +-level."""
    raise NotImplementedError


def rice_upcrossing_rate(level: float, sigma_x: float, sigma_xdot: float) -> float:
    """Return Rice's expected rate [1/s] of upward crossings of ``level`` by a stationary Gaussian process.

    nu(level) = (sigma_xdot / (2 pi sigma_x)) exp(-level^2 / (2 sigma_x^2)); it needs a finite
    derivative variance, which the Dryden gust itself does not have.
    """
    raise NotImplementedError


def count_upcrossings(x: np.ndarray, level: float) -> int:
    """Return the number of samples n with x[n] < level <= x[n+1] (upward crossings of the level)."""
    raise NotImplementedError


def variance_estimate_relative_sd(T: float, L: float, V: float) -> float:
    """Return the relative standard deviation of the sample variance of a Dryden w record of length T [s].

    Var(sigma_hat^2) = (2/T) integral R^2 dtau = 5 sigma^4 L / (4 V T), so the relative SD is
    sqrt(5 L / (4 V T)); valid for T >> L/V.
    """
    raise NotImplementedError


def record_length_for_variance(rel_sd: float, L: float, V: float) -> float:
    """Return the record length T [s] that gives the requested relative SD of the variance estimate: T = 5L/(4 V rel_sd^2)."""
    raise NotImplementedError


def main(out_dir: str | Path = ".") -> None:
    """Generate a gust record, compare its statistics with the Dryden model, and report the aircraft's response."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    V, S, m, rho, cla = (AIRCRAFT[k] for k in ("airspeed", "area", "mass", "rho", "cl_alpha"))
    sw, L = TURBULENCE["sigma_w"], TURBULENCE["L_w"]
    tau = heave_time_constant(m, rho, V, S, cla)
    b = V / L
    print(f"aircraft: m = {m:.0f} kg, S = {S:.0f} m^2, V = {V:.0f} m/s, heave time constant tau = {tau:.3f} s")
    print(f"turbulence: sigma_w = {sw} m/s, L_w = {L:.0f} m, break frequency V/L = {b:.3f} rad/s = {b / (2 * np.pi):.4f} Hz")
    print(f"integral of Phi_w  = {integrated_variance(lambda om: dryden_psd(om, sw, L, V, 'w')):.6f}  (sigma_w^2 = {sw**2})")
    print(f"integral of Phi_u  = {integrated_variance(lambda om: dryden_psd(om, sw, L, V, 'u')):.6f}")
    for om in (0.0, 0.2, 2.0):
        wk = psd_from_autocorrelation(om, lambda s: float(dryden_autocorrelation(s, sw, L, V, 'w')), 40.0 * L / V)
        print(f"Wiener-Khinchin at omega = {om:4.1f} rad/s: from R(tau) {wk:.6f}, closed form {float(dryden_psd(om, sw, L, V, 'w')):.6f}")

    duration, dt, seed = 2000.0, 0.02, 2019
    t, w = generate_gust(duration, dt, sw, L, V, "w", seed=seed)
    print(f"\nrecord: {duration:.0f} s at {1 / dt:.0f} Hz, seed {seed}  (synthetic: generated, not measured)")
    print(f"sample mean {w.mean():+.4f} m/s, sample SD {w.std():.4f} m/s (model {sw});"
          f" predicted relative SD of the variance estimate {variance_estimate_relative_sd(duration, L, V):.3f}")
    print(f"record for 2 % variance accuracy: {record_length_for_variance(0.02, L, V) / 60:.0f} min")

    f_w, P_w = welch_psd(w, dt, nperseg=8192)
    for lo, hi in ((0.01, 0.1), (0.1, 1.0), (1.0, 10.0)):
        est = band_average(f_w, P_w, lo, hi)
        th = band_average(f_w, dryden_psd_hz(f_w, sw, L, V, "w"), lo, hi)
        print(f"Welch / Dryden in {lo:5.2f}-{hi:5.2f} Hz: {est / th:.3f}")
    f_per, P_per = periodogram(w, dt)
    band = (f_per > 0.1) & (f_per < 1.0)
    rel = P_per[band] / dryden_psd_hz(f_per[band], sw, L, V, "w")
    print(f"periodogram: SD/mean of (estimate / model) over the bins between 0.1 and 1 Hz = {np.std(rel) / np.mean(rel):.2f}"
          f"  (1.0 for an inconsistent estimator; Welch with 8192-point segments: "
          f"{np.std(P_w[(f_w > 0.1) & (f_w < 1.0)] / dryden_psd_hz(f_w[(f_w > 0.1) & (f_w < 1.0)], sw, L, V, 'w')) / np.mean(P_w[(f_w > 0.1) & (f_w < 1.0)] / dryden_psd_hz(f_w[(f_w > 0.1) & (f_w < 1.0)], sw, L, V, 'w')):.2f})")

    lags = np.arange(0, int(round(4 * L / V / dt)) + 1)
    R_est = autocorrelation_estimate(w, lags[-1])
    R_th = dryden_autocorrelation(lags * dt, sw, L, V, "w")
    print(f"autocorrelation: R(0) est {R_est[0]:.4f} vs {R_th[0]:.4f}; first zero of R_w at tau = 2L/V = {2 * L / V:.1f} s,"
          f" estimate crosses zero at {lags[np.argmax(R_est < 0)] * dt:.1f} s")

    var_v, var_a = heave_response_variances(sw, L, V, tau)
    v, a = simulate_heave(t, w, tau)
    sig_n = np.sqrt(var_a) / G0
    print(f"\nheave response: sigma_v = {np.sqrt(var_v):.4f} m/s (closed form) vs {v.std():.4f} simulated;"
          f" numeric quad {np.sqrt(response_variance_numeric(sw, L, V, tau, 'v')):.4f}")
    print(f"               sigma_a = {np.sqrt(var_a):.4f} m/s^2 (closed form) vs {a.std():.4f} simulated;"
          f" RMS load factor {sig_n:.4f} g; sharp-edged-gust value sigma_w/(g tau) = {sw / (G0 * tau):.4f} g")
    for lvl in (0.2, 0.3, 0.5):
        print(f"   fraction of time |dn| > {lvl} g: {gaussian_exceedance_fraction(lvl * G0, np.sqrt(var_a)):.2e}"
              f"  simulated {np.mean(np.abs(a) > lvl * G0):.2e}")
    nu0 = rice_upcrossing_rate(0.0, np.sqrt(var_v), np.sqrt(var_a))
    print(f"zero upcrossings of v: Rice {nu0 * duration:.0f} expected in {duration:.0f} s, counted {count_upcrossings(v, 0.0)}")
    nu1 = rice_upcrossing_rate(np.sqrt(var_v), np.sqrt(var_v), np.sqrt(var_a))
    print(f"upcrossings of v = +1 sigma: Rice {nu1 * duration:.0f}, counted {count_upcrossings(v, np.sqrt(var_v))}")

    # Computational experiment: scatter of the variance estimate against record length.
    print("\nvariance-estimate scatter over 24 seeds (synthetic records):")
    print("   T [s]   predicted rel SD   measured rel SD")
    for T in (100.0, 400.0, 1600.0):
        vars_ = [generate_gust(T, 0.1, sw, L, V, "w", seed=100 + k)[1].var() for k in range(24)]
        print(f"   {T:6.0f}   {variance_estimate_relative_sd(T, L, V):16.3f}   {np.std(vars_) / sw**2:15.3f}")

    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    axs[0, 0].plot(t[: int(120 / dt)], w[: int(120 / dt)], lw=0.6)
    axs[0, 0].set_xlabel("t [s]")
    axs[0, 0].set_ylabel("w_g [m/s]")
    axs[0, 0].set_title("synthetic Dryden gust record (first 120 s)")
    f_fit = f_w[1:]
    axs[0, 1].loglog(f_fit, P_w[1:], lw=0.8, label="Welch estimate (8192-point Hann)")
    axs[0, 1].loglog(f_fit, dryden_psd_hz(f_fit, sw, L, V, "w"), "k--", label="Dryden Phi_w")
    axs[0, 1].loglog(f_fit, dryden_psd_hz(f_fit, sw, L, V, "u"), "k:", label="Dryden Phi_u")
    axs[0, 1].set_xlabel("f [Hz]")
    axs[0, 1].set_ylabel("G [(m/s)^2/Hz]")
    axs[0, 1].set_ylim(1e-6, 1e2)
    axs[0, 1].legend(fontsize=8)
    axs[0, 1].set_title("gust PSD")
    axs[1, 0].plot(lags * dt, R_est, label="estimate")
    axs[1, 0].plot(lags * dt, R_th, "k--", label="Dryden R_w")
    axs[1, 0].axhline(0, color="gray", lw=0.5)
    axs[1, 0].set_xlabel("lag [s]")
    axs[1, 0].set_ylabel("R_w [(m/s)^2]")
    axs[1, 0].legend()
    axs[1, 0].set_title("autocorrelation")
    f_a, P_a = welch_psd(a / G0, dt, nperseg=8192)
    Hv, Ha = heave_response(2 * np.pi * f_fit, tau)
    axs[1, 1].loglog(f_fit, P_a[1:], lw=0.8, label="Welch estimate of load-factor PSD")
    axs[1, 1].loglog(f_fit, np.abs(Ha) ** 2 * dryden_psd_hz(f_fit, sw, L, V, "w") / G0**2, "k--", label="|H_a|^2 Phi_w / g^2")
    axs[1, 1].set_xlabel("f [Hz]")
    axs[1, 1].set_ylabel("G_n [g^2/Hz]")
    axs[1, 1].set_ylim(1e-8, 1e0)
    axs[1, 1].legend(fontsize=8)
    axs[1, 1].set_title("aircraft load-factor response")
    for ax in axs.ravel():
        ax.grid(True, which="both", lw=0.3)
    fig.suptitle("Dryden gust generator and heave response (synthetic data)")
    fig.tight_layout()
    fig.savefig(out_dir / "gust_spectrum.png", dpi=150)
    print(f"\nsaved {out_dir / 'gust_spectrum.png'}")


if __name__ == "__main__":
    main()
