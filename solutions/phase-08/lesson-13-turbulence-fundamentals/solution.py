"""P08.13 Turbulence Fundamentals: reference solution.

What are the turbulence intensity, integral scale, dissipation rate and Kolmogorov scale behind a
wind-tunnel turbulence grid, is a 1 mm hot wire small enough to see the smallest eddies, and how
long must the record be? This module analyses a velocity record u(t) from a single-point probe:

* Reynolds decomposition u = U + u' and the moments (variance, skewness, flatness);
* the autocorrelation rho(tau) and the integral time scale T = int rho dtau (to the first zero);
* Taylor's frozen-turbulence hypothesis, which turns times into lengths, x = U t, so L = U T;
* the power spectrum by Welch's method, its -5/3 inertial range, the dissipation rate epsilon from
  Kolmogorov's law E11(k) = C1 epsilon^(2/3) k^(-5/3) (C1 = 0.49), and the Kolmogorov scales;
* the record length needed for a given statistical uncertainty.

THE DATA ARE SYNTHETIC. ``synthetic_signal`` generates a seeded Gaussian random signal whose
spectrum is a stated model, the von Karman longitudinal spectrum

    E11(k) = (2 sigma^2 L / pi) / (1 + (1.339 k L)^2)^(5/6),   int_0^inf E11 dk = sigma^2,

with integral length scale exactly L, convected past the probe at U (Taylor's hypothesis built in).
Because the variance, the integral scale and the dissipation rate implied by the model are known,
the analysis can be checked: it must recover them within the statistical scatter of a finite
record. A synthetic Gaussian signal has no intermittency and no energy cascade; its skewness is 0
and its flatness 3, unlike real turbulence, and it has no dissipation range (the model spectrum
continues as k^(-5/3) to the Nyquist frequency). Simulation is never measurement.

Units: SI (m, s).

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes turbulence_synthetic.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-08/lesson-13-turbulence-fundamentals
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

C1_KOLMOGOROV = 0.49            # one-dimensional Kolmogorov constant for E11 (about 18/55 of C = 1.5)
NU_AIR_20C = 1.81e-5 / 1.204    # m^2/s, air at 20 C and 101 325 Pa (shared Phase 08 numbers)


# ---------------------------------------------------------------- given: the synthetic signal

def von_karman_e11(k: np.ndarray, sigma: float, L: float) -> np.ndarray:
    """One-sided von Karman longitudinal wavenumber spectrum (m^3/s^2) of a velocity component with
    variance sigma^2 and integral length scale L: (2 sigma^2 L / pi) / (1 + (1.339 k L)^2)^(5/6)."""
    k = np.asarray(k, dtype=float)
    return (2 * sigma ** 2 * L / np.pi) / (1 + (1.339 * k * L) ** 2) ** (5 / 6)


def model_epsilon(sigma: float, L: float) -> float:
    """Dissipation rate implied by the von Karman model: at large k, E11 -> (2 sigma^2 L / pi)
    (1.339 k L)^(-5/3), which equals C1 epsilon^(2/3) k^(-5/3) when
    epsilon = [(2 sigma^2 L / pi) (1.339 L)^(-5/3) / C1]^(3/2)."""
    return float(((2 * sigma ** 2 * L / np.pi) * (1.339 * L) ** (-5 / 3) / C1_KOLMOGOROV) ** 1.5)


def synthetic_signal(n: int, fs: float, U: float, sigma: float, L: float, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """SYNTHETIC velocity record (not a measurement): n samples at fs Hz of u(t) = U + u'(t), where
    u' is a zero-mean Gaussian random process with the von Karman spectrum convected at U, i.e. the
    one-sided frequency spectrum S(f) = (2 pi / U) E11(2 pi f / U). Built in the frequency domain:
    each Fourier coefficient f_k = k fs / n (0 < k < n/2) gets an independent complex Gaussian
    amplitude with E|X_k|^2 = n^2 S(f_k) df / 2, then an inverse real FFT. Seeded with
    numpy.random.default_rng(seed). Returns (t, u)."""
    rng = np.random.default_rng(seed)
    df = fs / n
    f = np.arange(n // 2 + 1) * df
    S = (2 * np.pi / U) * von_karman_e11(2 * np.pi * f / U, sigma, L)
    Z = (rng.standard_normal(f.size) + 1j * rng.standard_normal(f.size)) / np.sqrt(2)
    X = n * np.sqrt(S * df / 2) * Z
    X[0] = 0.0
    X[-1] = 0.0
    u = U + np.fft.irfft(X, n)
    return np.arange(n) / fs, u


def band_variance(n: int, fs: float, U: float, sigma: float, L: float) -> float:
    """Expected variance of ``synthetic_signal``: the model spectrum summed over the bins it
    synthesises (0 < f < fs/2). Slightly less than sigma^2: the energy above the Nyquist frequency
    and below the first bin is not in a finite, sampled record."""
    df = fs / n
    f = np.arange(1, n // 2) * df
    return float(np.sum((2 * np.pi / U) * von_karman_e11(2 * np.pi * f / U, sigma, L)) * df)


# ---------------------------------------------------------------- statistics

def reynolds_decomposition(u: np.ndarray) -> tuple[float, np.ndarray]:
    """Split a record into its time mean and fluctuation: u = U + u'. Returns (U, u')."""
    U = float(np.mean(u))
    return U, np.asarray(u, dtype=float) - U


def moments(u: np.ndarray) -> dict[str, float]:
    """Mean, variance <u'^2>, standard deviation, skewness <u'^3>/sigma^3 and flatness <u'^4>/sigma^4
    of a record (population moments, i.e. divided by n)."""
    U, up = reynolds_decomposition(u)
    var = float(np.mean(up ** 2))
    s = np.sqrt(var)
    return {"mean": U, "var": var, "std": float(s), "skew": float(np.mean(up ** 3) / s ** 3),
            "flat": float(np.mean(up ** 4) / var ** 2)}


def autocorrelation(up: np.ndarray, max_lag: int) -> np.ndarray:
    """Normalised autocorrelation rho(k) = R(k)/R(0), k = 0..max_lag, of a zero-mean record, with the
    biased estimator R(k) = (1/n) sum_{j=0}^{n-1-k} u'_j u'_{j+k}, computed by FFT with zero padding
    to at least 2n (so it is the linear, not circular, correlation)."""
    up = np.asarray(up, dtype=float)
    n = up.size
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    F = np.fft.rfft(up, nfft)
    R = np.fft.irfft(F * np.conj(F), nfft)[: max_lag + 1] / n
    return R / R[0]


def integral_time_scale(rho: np.ndarray, dt: float) -> float:
    """Integral time scale T = int_0^tau0 rho dtau by the trapezoid rule, up to the first lag tau0
    where rho falls to zero or below (the noisy tail beyond it is excluded)."""
    rho = np.asarray(rho, dtype=float)
    neg = np.nonzero(rho <= 0)[0]
    end = neg[0] + 1 if neg.size else rho.size
    return float(np.trapezoid(rho[:end], dx=dt))


def welch_psd(u: np.ndarray, fs: float, nperseg: int) -> tuple[np.ndarray, np.ndarray]:
    """One-sided power spectral density (units^2/Hz) by Welch's method: segments of nperseg samples
    with 50 % overlap, each with its mean removed and a Hann window applied, |FFT|^2 averaged and
    scaled so that the density integrates to the variance: S = 2 |X|^2 / (fs sum(w^2)) (not doubled
    at f = 0 and at the Nyquist frequency). Returns (f, S) with f = k fs / nperseg, k = 0..nperseg/2."""
    u = np.asarray(u, dtype=float)
    w = np.hanning(nperseg + 1)[:-1]            # periodic Hann window
    step = nperseg // 2
    starts = range(0, u.size - nperseg + 1, step)
    acc = np.zeros(nperseg // 2 + 1)
    for k, s in enumerate(starts, start=1):
        seg = u[s:s + nperseg]
        X = np.fft.rfft((seg - seg.mean()) * w)
        acc += np.abs(X) ** 2
    S = acc / k / (fs * np.sum(w ** 2))
    S[1:-1] *= 2
    return np.fft.rfftfreq(nperseg, 1 / fs), S


def integral_scale_from_spectrum(f: np.ndarray, S: np.ndarray, var: float, U: float, nbins: int = 3) -> float:
    """Integral length scale from the low-frequency plateau of the spectrum: for a one-sided
    frequency spectrum S(0) = 4 var T, so L = U T = U S(0) / (4 var), with S(0) estimated as the mean
    of the first ``nbins`` nonzero-frequency bins (the f = 0 bin of a mean-removed segment is empty)."""
    return float(U * np.mean(S[1:1 + nbins]) / (4 * var))


def fit_inertial_range(f: np.ndarray, S: np.ndarray, U: float, f_lo: float, f_hi: float) -> tuple[float, float]:
    """Fit the inertial range f_lo <= f <= f_hi of a frequency spectrum. Converts to wavenumber with
    Taylor's hypothesis (k = 2 pi f / U, E11 = U S / (2 pi)); returns the least-squares slope of
    log E11 against log k, and the dissipation rate from Kolmogorov's law,
    epsilon = mean over the band of (E11 k^(5/3) / C1)^(3/2). Returns (slope, epsilon)."""
    sel = (f >= f_lo) & (f <= f_hi)
    k = 2 * np.pi * f[sel] / U
    E = U * S[sel] / (2 * np.pi)
    slope = float(np.polyfit(np.log(k), np.log(E), 1)[0])
    eps = float(np.mean((E * k ** (5 / 3) / C1_KOLMOGOROV) ** 1.5))
    return slope, eps


def kolmogorov_scales(eps: float, nu: float) -> dict[str, float]:
    """Kolmogorov length eta = (nu^3/eps)^(1/4), time tau = (nu/eps)^(1/2) and velocity
    v = (nu eps)^(1/4): the scales at which the eddy Reynolds number eta v / nu = 1."""
    return {"eta": (nu ** 3 / eps) ** 0.25, "tau": (nu / eps) ** 0.5, "v": (nu * eps) ** 0.25}


def record_length_for(T_int: float, rel_err: float) -> float:
    """Record duration (s) for a one-standard-deviation relative error rel_err in the variance of a
    Gaussian signal with integral time scale T_int: var(sigma_hat^2)/sigma^4 ~ 2 T_int / T, so
    T = 2 T_int / rel_err^2 (valid for T >> T_int)."""
    return 2 * T_int / rel_err ** 2


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save turbulence_synthetic.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    U, sigma, L, fs, n, nu = 10.0, 0.3, 0.05, 20000.0, 2 ** 20, NU_AIR_20C
    t, u = synthetic_signal(n, fs, U, sigma, L, seed=2026)
    print(f"SYNTHETIC grid-turbulence record (von Karman spectrum; U = {U} m/s, sigma = {sigma} m/s, L = {L} m),")
    print(f"  {n} samples at {fs:.0f} Hz = {n / fs:.1f} s, seed 2026")
    m = moments(u)
    bv = band_variance(n, fs, U, sigma, L)
    print(f"mean {m['mean']:.4f} m/s; variance {m['var']:.5f} m^2/s^2 (model in band {bv:.5f}, total {sigma ** 2:.5f}); "
          f"Tu = {m['std'] / m['mean'] * 100:.2f} %")
    print(f"skewness {m['skew']:.3f}, flatness {m['flat']:.3f} (Gaussian: 0 and 3)")

    _, up = reynolds_decomposition(u)
    rho = autocorrelation(up, 4000)
    T = integral_time_scale(rho, 1 / fs)
    print(f"integral time scale {T * 1e3:.3f} ms -> L = U T = {U * T * 1e3:.2f} mm (model {L * 1e3:.1f} mm, T = {L / U * 1e3:.1f} ms)")
    f, S = welch_psd(u, fs, 2 ** 14)
    print(f"Welch PSD (nperseg 2^14): integral of S = {np.trapezoid(S, f):.5f} m^2/s^2; "
          f"L from S(0) = {integral_scale_from_spectrum(f, S, m['var'], U) * 1e3:.2f} mm")
    slope, eps = fit_inertial_range(f, S, U, 300.0, 3000.0)
    eps_m = model_epsilon(sigma, L)
    print(f"inertial range 300-3000 Hz: slope {slope:.3f} (Kolmogorov -5/3 = {-5 / 3:.3f}); "
          f"epsilon {eps:.4f} m^2/s^3 (model {eps_m:.4f})")
    ks = kolmogorov_scales(eps, nu)
    print(f"Kolmogorov scales (nu = {nu:.4e}): eta = {ks['eta'] * 1e3:.3f} mm, tau = {ks['tau'] * 1e3:.3f} ms, "
          f"v = {ks['v']:.4f} m/s; f_eta = U/(2 pi eta) = {U / (2 * np.pi * ks['eta']):.0f} Hz")
    print(f"ratio L/eta = {L / ks['eta']:.0f}; Re_L = sigma L / nu = {sigma * L / nu:.0f}; "
          f"a 1 mm wire spans {1e-3 / ks['eta']:.1f} eta")
    for e in (0.10, 0.05, 0.02):
        print(f"record length for {e * 100:.0f} % (1 sigma) on the variance: {record_length_for(L / U, e):.2f} s")

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.0))
    ax[0].plot(t[:4000] * 1e3, u[:4000], lw=0.6)
    ax[0].axhline(m["mean"], color="k", lw=0.8)
    ax[0].set_xlabel("t (ms)")
    ax[0].set_ylabel("u (m/s)")
    ax[0].set_title("SYNTHETIC signal (first 0.2 s)")
    tau = np.arange(rho.size) / fs
    ax[1].plot(tau * 1e3, rho, label="estimated ρ(τ), synthetic")
    ax[1].axvline(T * 1e3, color="k", ls=":", label=f"T = {T * 1e3:.2f} ms")
    ax[1].set_xlim(0, 60)
    ax[1].set_xlabel("τ (ms)")
    ax[1].set_ylabel("ρ")
    ax[1].set_title("autocorrelation (synthetic)")
    ax[1].legend(fontsize=8)
    ax[2].loglog(f[1:], S[1:], lw=0.7, label="Welch PSD, synthetic")
    ax[2].loglog(f[1:], (2 * np.pi / U) * von_karman_e11(2 * np.pi * f[1:] / U, sigma, L), "k--", lw=0.8, label="model (von Kármán)")
    ax[2].set_xlabel("f (Hz)")
    ax[2].set_ylabel("S (m²/s²/Hz)")
    ax[2].set_title("spectrum: −5/3 range (synthetic)")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "turbulence_synthetic.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'turbulence_synthetic.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
