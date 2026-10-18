"""P06.17 Vibration Isolation: starter code.

Complete every function that raises NotImplementedError. The P06.16 helpers (psd_loglog, psd_grms,
freq_grid, synthesize_random, simulate_sdof_base, miles_grms) and main() are given. When finished,
running this file prints the lesson's numbers and writes the figure.

How soft should the mounts under an inertial measurement unit (IMU) be? An isolator is a spring and a
damper between the vehicle structure and the sensitive box: above sqrt(2) times its natural frequency
it attenuates vibration, at its natural frequency it amplifies it, and the softer it is the further the
box sways. This module

* gives the base-excitation transmissibility of a viscous isolator (magnitude and phase) and of a
  hysteretic (elastomer-like) isolator, the exact resonant peak, and the corner frequency needed for a
  target attenuation at a disturbance frequency;
* checks the closed-form transmissibility against the steady state of a ``solve_ivp`` simulation;
* computes the random response of the isolated IMU (transmitted Grms and sway) by the lesson-16 integral
  of |H|^2 G, checked against Miles' equation and a seeded SYNTHETIC time-domain simulation;
* sweeps natural frequency and damping and picks the isolator that meets a set of requirements.

Conventions: f in Hz, r = f/fn, G one-sided in g^2/Hz, accelerations in g, displacements in m,
stiffness in N/m, damping in N s/m. Random records are SYNTHETIC.

Run:  python main.py [out_dir]         (prints the lesson's numbers, writes vibration_isolation.png)
Test: pytest curriculum/phase-06/lesson-17-vibration-isolation
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.optimize import brentq
from scipy.signal import lfilter

G0 = 9.80665

# NASA-STD-7001C Table 2 component minimum workmanship random vibration level (Hz, g^2/Hz), 6.78 Grms (P06.16).
SPEC_F = np.array([20.0, 80.0, 500.0, 2000.0])
SPEC_G = np.array([0.01, 0.04, 0.04, 0.01])


# ---------------------------------------------------------------- given helpers (from P06.16)

def psd_loglog(f: np.ndarray, f_bp: np.ndarray, G_bp: np.ndarray) -> np.ndarray:
    """Breakpoint spectrum evaluated at ``f``: straight lines on log-log axes, zero outside the band."""
    f = np.asarray(f, float)
    out = np.zeros_like(f)
    inside = (f >= f_bp[0]) & (f <= f_bp[-1])
    out[inside] = np.exp(np.interp(np.log(f[inside]), np.log(f_bp), np.log(G_bp)))
    return out


def psd_grms(f_bp: np.ndarray, G_bp: np.ndarray) -> float:
    """Exact overall Grms of a breakpoint spectrum (P06.16)."""
    area = 0.0
    for f1, f2, G1, G2 in zip(f_bp[:-1], f_bp[1:], G_bp[:-1], G_bp[1:]):
        n = np.log(G2 / G1) / np.log(f2 / f1)
        area += G1 * f1 * (np.log(f2 / f1) if abs(n + 1) < 1e-12 else ((f2 / f1) ** (n + 1) - 1) / (n + 1))
    return float(np.sqrt(area))


def freq_grid(f_bp: np.ndarray, fn: float, zeta: float, n_per_decade: int = 3000) -> np.ndarray:
    """Log grid over the spectrum band with extra points packed around fn (for integrals of |H|^2 G)."""
    lo, hi = float(f_bp[0]), float(f_bp[-1])
    base = np.geomspace(lo, hi, int(np.ceil(np.log10(hi / lo) * n_per_decade)) + 1)
    bw = max(zeta, 1e-3) * fn
    dense = fn + bw * np.linspace(-40.0, 40.0, 4001) ** 3 / 40.0 ** 2
    pts = np.concatenate([base, dense[(dense > lo) & (dense < hi)], np.asarray(f_bp, float)])
    return np.unique(pts)


def synthesize_random(f_bp: np.ndarray, G_bp: np.ndarray, fs: float, duration: float,
                      seed: int) -> tuple[np.ndarray, np.ndarray]:
    """SYNTHETIC stationary Gaussian base acceleration [g] with one-sided PSD G (P06.16)."""
    rng = np.random.default_rng(seed)
    n = int(round(fs * duration))
    fk = np.arange(n // 2 + 1) * fs / n
    X = (rng.standard_normal(fk.size) + 1j * rng.standard_normal(fk.size)) * np.sqrt(psd_loglog(fk, f_bp, G_bp) * fs / n)
    X[0] = 0.0
    if n % 2 == 0:
        X[-1] = 0.0
    return np.arange(n) / fs, np.fft.irfft(X, n) * (n / 2)


def simulate_sdof_base(a_base: np.ndarray, dt: float, fn: float, zeta: float) -> tuple[np.ndarray, np.ndarray]:
    """Exact first-order-hold simulation of z'' + 2 zeta wn z' + wn^2 z = -a_base (P06.16).
    Returns (absolute acceleration [units of a_base], relative displacement z [units of a_base * s^2])."""
    wn = 2 * np.pi * fn
    A = np.array([[0.0, 1.0], [-wn * wn, -2 * zeta * wn]])
    M = np.zeros((4, 4))
    M[:2, :2], M[:2, 2], M[2, 3] = A * dt, [0.0, -dt], 1.0
    E = expm(M)
    Phi, G1, G2 = E[:2, :2], E[:2, 2], E[:2, 3]
    Gm = G1 - G2
    den = np.array([1.0, -np.trace(Phi), np.linalg.det(Phi)])
    adj1 = np.array([[-Phi[1, 1], Phi[0, 1]], [Phi[1, 0], -Phi[0, 0]]])

    def num(C: np.ndarray) -> np.ndarray:
        return np.array([C @ G2, C @ Gm + C @ adj1 @ G2, C @ adj1 @ Gm])

    a_base = np.asarray(a_base, float)
    acc = lfilter(num(np.array([-wn * wn, -2 * zeta * wn])), den, a_base)
    z = lfilter(num(np.array([1.0, 0.0])), den, a_base)
    return acc, z


# ---------------------------------------------------------------- isolator basics

def natural_frequency(k: float, m: float) -> float:
    """fn = sqrt(k/m) / (2 pi) [Hz] for total mount stiffness k [N/m] and isolated mass m [kg]."""
    raise NotImplementedError


def stiffness_for_frequency(fn: float, m: float) -> float:
    """Total mount stiffness k = m (2 pi fn)^2 [N/m] that puts the isolation mode at fn."""
    raise NotImplementedError


def damping_coefficient(zeta: float, k: float, m: float) -> float:
    """Viscous damping coefficient c = 2 zeta sqrt(k m) [N s/m]."""
    raise NotImplementedError


def static_deflection(n_load: float, fn: float) -> float:
    """Mount deflection [m] under a steady acceleration of n_load g: n g0 / wn^2."""
    raise NotImplementedError


def transmissibility(r: np.ndarray, zeta: float) -> np.ndarray:
    """|X/Y| of a viscously damped isolator under base motion (equal for displacement and acceleration):
    sqrt( (1 + (2 zeta r)^2) / ((1 - r^2)^2 + (2 zeta r)^2) )."""
    raise NotImplementedError


def transmissibility_phase(r: np.ndarray, zeta: float) -> np.ndarray:
    """Phase lag [rad] of the isolated mass behind the base, in [0, pi): the argument of
    (1 + 2i zeta r) / (1 - r^2 + 2i zeta r), negated."""
    raise NotImplementedError


def transmissibility_hysteretic(r: np.ndarray, eta: float) -> np.ndarray:
    """|X/Y| for a hysteretic (structural) damper, complex stiffness k (1 + i eta), typical of elastomers:
    sqrt( (1 + eta^2) / ((1 - r^2)^2 + eta^2) ). Falls as 1/r^2 at high frequency."""
    raise NotImplementedError


def peak_transmissibility(zeta: float) -> tuple[float, float]:
    """Exact resonant peak of the viscous transmissibility: (r_peak, T_peak), with
    r_peak^2 = (sqrt(1 + 8 zeta^2) - 1) / (4 zeta^2) (tends to 1 as zeta -> 0)."""
    raise NotImplementedError


def zeta_for_peak(T_peak: float) -> float:
    """Damping ratio whose resonant peak transmissibility equals ``T_peak`` (> 1)."""
    raise NotImplementedError


def fn_for_attenuation(f_dist: float, T_target: float, zeta: float) -> float:
    """Highest isolator natural frequency [Hz] whose transmissibility at the disturbance frequency
    ``f_dist`` is ``T_target`` (< 1): solve T(f_dist/fn, zeta) = T_target for r > sqrt 2."""
    raise NotImplementedError


# ---------------------------------------------------------------- independent check: steady state by solve_ivp

def steady_state_ivp(r: float, zeta: float, fn: float = 1.0, n_cycles: int = 0,
                     rtol: float = 1e-10, atol: float = 1e-12) -> tuple[float, float]:
    """Integrate m x'' = -k (x - y) - c (x' - y') with y = sin(w t), w = r wn, from rest with
    ``solve_ivp`` (DOP853) until the transient has decayed (exp(-zeta wn t) < 1e-8, at least 30 cycles,
    or ``n_cycles`` if given), then fit x = a sin(w t) + b cos(w t) over the last 5 cycles.
    Returns (amplitude ratio |X/Y|, phase lag [rad])."""
    raise NotImplementedError


# ---------------------------------------------------------------- random response (P06.16 style)

def random_response(f_bp: np.ndarray, G_bp: np.ndarray, fn: float, zeta: float) -> tuple[float, float]:
    """(transmitted Grms [g], relative-displacement RMS [m]) of the isolated mass for a base acceleration
    PSD G (g^2/Hz): sqrt(integral |H|^2 G df) and g0 sqrt(integral G / ((wn^2 - w^2)^2 + (2 zeta wn w)^2) df)."""
    raise NotImplementedError


def miles_grms(fn: float, zeta: float, G_fn: float) -> float:
    """Miles' equation sqrt((pi/2) fn Q G(fn)), Q = 1/(2 zeta) (P06.16)."""
    return float(np.sqrt(0.5 * np.pi * fn * G_fn / (2 * zeta)))


# ---------------------------------------------------------------- design

def design_sweep(fns: np.ndarray, zetas: np.ndarray, f_bp: np.ndarray, G_bp: np.ndarray,
                 n_static: float = 0.0, f_hf: float = 500.0) -> dict[str, np.ndarray]:
    """Evaluate every (fn, zeta) pair. Returns arrays of shape (len(fns), len(zetas)):
    'grms' transmitted Grms, 'sway3' = static deflection under n_static g + 3 sigma relative
    displacement [m], 'tpeak' resonant peak transmissibility, 'thf' transmissibility at f_hf [Hz]."""
    raise NotImplementedError


def choose_isolator(fns: np.ndarray, zetas: np.ndarray, f_bp: np.ndarray, G_bp: np.ndarray, m: float,
                    fn_min: float, sway_max: float, tpeak_max: float, thf_max: float = np.inf,
                    n_static: float = 0.0, f_hf: float = 500.0) -> dict[str, float]:
    """Pick the (fn, zeta) on the grid with the lowest transmitted Grms that satisfies fn >= fn_min,
    static + 3 sigma sway <= sway_max, peak transmissibility <= tpeak_max and transmissibility at f_hf
    <= thf_max. Returns a dict with fn, zeta, k, c, grms, sway3, tpeak, thf, or raises ValueError if
    nothing on the grid is feasible."""
    raise NotImplementedError


# ---------------------------------------------------------------- main

M_IMU = 0.35          # kg, illustrative tactical-grade MEMS IMU with its bracket
FN_MIN = 60.0         # Hz, keep the mount mode well above the guidance bandwidth and its harmonics
SWAY_MAX = 0.80e-3    # m, static + 3 sigma, set by the connector and cable slack (illustrative)
TPEAK_MAX = 4.0       # resonant amplification allowed during a sine sweep (illustrative)
F_HF, THF_MAX = 500.0, 0.06   # transmissibility allowed at 500 Hz and above, where the IMU's sensors
                              # rectify vibration into bias (illustrative)
N_STATIC = 5.0        # g, quasi-static axial launch acceleration (typical)


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save vibration_isolation.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    f_bp, G_bp = SPEC_F, SPEC_G
    print(f"base input: NASA-STD-7001C Table 2 workmanship level, {psd_grms(f_bp, G_bp):.3f} Grms; IMU mass {M_IMU} kg")
    g_hard, _ = random_response(f_bp, G_bp, 3000.0, 0.05)
    print(f"rigidly mounted IMU: {psd_grms(f_bp, G_bp):.3f} Grms (the input itself); on a stiff bracket with its "
          f"mode at 3000 Hz, zeta = 0.05: {g_hard:.3f} Grms")

    print("\ntransmissibility: closed form vs solve_ivp steady state (zeta = 0.1)")
    print("     r     closed      ivp        |diff|     phase closed  phase ivp [deg]")
    r_check = np.array([0.5, 1.0, np.sqrt(2), 3.0, 10.0])
    ivp = [steady_state_ivp(float(r), 0.1) for r in r_check]
    for r, (Ti, pi_) in zip(r_check, ivp):
        Tc, pc = float(transmissibility(r, 0.1)), float(transmissibility_phase(r, 0.1))
        print(f"  {r:5.3f}   {Tc:.6f}   {Ti:.6f}   {abs(Tc - Ti):.1e}    {np.degrees(pc):8.3f}     {np.degrees(pi_):8.3f}")

    print("\npeak transmissibility and high-frequency isolation (r = 10):")
    for z in (0.02, 0.05, 0.1, 0.2, 0.4):
        rp, Tp = peak_transmissibility(z)
        print(f"  zeta {z:4.2f}: peak {Tp:6.2f} at r = {rp:.4f}; T(10) viscous {float(transmissibility(10, z)):.4f}; "
              f"hysteretic eta = 2 zeta: {float(transmissibility_hysteretic(10, 2 * z)):.4f}")
    print(f"  zeta giving a peak of {TPEAK_MAX}: {zeta_for_peak(TPEAK_MAX):.4f}")

    fns = np.arange(40.0, 201.0, 5.0)
    zetas = np.round(np.arange(0.05, 0.3001, 0.025), 3)
    best = choose_isolator(fns, zetas, f_bp, G_bp, M_IMU, FN_MIN, SWAY_MAX, TPEAK_MAX, THF_MAX, N_STATIC, F_HF)
    print(f"\nrequirements: fn >= {FN_MIN} Hz, static ({N_STATIC} g) + 3 sigma sway <= {SWAY_MAX * 1e3} mm, "
          f"peak T <= {TPEAK_MAX}, T({F_HF:.0f} Hz) <= {THF_MAX}")
    print(f"chosen isolator: fn = {best['fn']:.0f} Hz, zeta = {best['zeta']}, k = {best['k'] / 1e3:.1f} kN/m total "
          f"({best['k'] / 4e3:.2f} kN/m per mount of 4), c = {best['c']:.2f} N s/m")
    print(f"  transmitted {best['grms']:.3f} Grms ({best['grms'] / psd_grms(f_bp, G_bp) * 100:.0f} % of the input), "
          f"sway {best['sway3'] * 1e3:.3f} mm, peak T {best['tpeak']:.2f}, T({F_HF:.0f} Hz) {best['thf']:.4f}")
    print(f"  an elastomer with nearly the same resonant peak (eta = 2 zeta = {2 * best['zeta']}) has T({F_HF:.0f} Hz) = "
          f"{float(transmissibility_hysteretic(F_HF / best['fn'], 2 * best['zeta'])):.4f}")
    fn, z = best["fn"], best["zeta"]
    g_full, zr = random_response(f_bp, G_bp, fn, z)
    g_miles = miles_grms(fn, z, float(psd_loglog(np.array([fn]), f_bp, G_bp)[0]))
    print(f"  integral {g_full:.3f} Grms vs Miles {g_miles:.3f} Grms; sway: static {static_deflection(N_STATIC, fn) * 1e3:.3f} mm + "
          f"3 x {zr * 1e3:.4f} mm")

    fs, dur, seed = 8192.0, 32.0, 17
    _, a = synthesize_random(f_bp, G_bp, fs, dur + 0.5, seed)
    acc, zz = simulate_sdof_base(a, 1 / fs, fn, z)
    sk = int(0.5 * fs)
    rms_a, rms_z = np.sqrt(np.mean(acc[sk:] ** 2)), G0 * np.sqrt(np.mean(zz[sk:] ** 2))
    u = 0.5 / np.sqrt(z * 2 * np.pi * fn * dur)
    print(f"\nSYNTHETIC time-domain check (seed {seed}, {dur:.0f} s): transmitted {rms_a:.3f} Grms "
          f"({(rms_a / g_full - 1) * 100:+.1f} %), sway {rms_z * 1e3:.4f} mm rms ({(rms_z / zr - 1) * 100:+.1f} %); "
          f"expected scatter about +/- {u * 100:.1f} % (1 sigma)")

    print(f"\nsweep at zeta = {z}: fn [Hz], transmitted Grms, static + 3 sigma sway [mm], T({F_HF:.0f} Hz)")
    sw = design_sweep(fns, np.array([z]), f_bp, G_bp, N_STATIC, F_HF)
    for i in range(0, fns.size, 4):
        print(f"  {fns[i]:5.0f}   {sw['grms'][i, 0]:6.3f}   {sw['sway3'][i, 0] * 1e3:.3f}   {sw['thf'][i, 0]:.4f}")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    r = np.geomspace(0.1, 30, 1000)
    for zi in (0.02, 0.1, 0.3):
        ax[0].loglog(r, transmissibility(r, zi), label=f"viscous zeta = {zi}")
    ax[0].loglog(r, transmissibility_hysteretic(r, 0.2), "k--", label="hysteretic eta = 0.2")
    ax[0].loglog(r_check, [v[0] for v in ivp], "o", label="solve_ivp, zeta = 0.1")
    ax[0].axvline(np.sqrt(2), color="0.6", lw=0.8)
    ax[0].set_xlabel("frequency ratio r = f / fn")
    ax[0].set_ylabel("transmissibility |X/Y|")
    ax[0].set_title("isolator transmissibility")
    ax[0].legend(fontsize=7)
    sw_all = design_sweep(fns, zetas, f_bp, G_bp, N_STATIC, F_HF)
    for jj in (0, 2, 4, 6, 8, 10):
        ax[1].plot(sw_all["sway3"][:, jj] * 1e3, sw_all["grms"][:, jj], ".-", label=f"zeta = {zetas[jj]}")
    ax[1].axvline(SWAY_MAX * 1e3, color="k", ls="--", lw=0.8)
    ax[1].plot(best["sway3"] * 1e3, best["grms"], "k*", ms=12, label="chosen")
    ax[1].set_xlabel("static + 3 sigma sway (mm)")
    ax[1].set_ylabel("transmitted Grms")
    ax[1].set_title("trade: sway against transmitted vibration")
    ax[1].legend(fontsize=7)
    f = np.geomspace(15, 2500, 1500)
    G = psd_loglog(f, f_bp, G_bp)
    ax[2].loglog(f, np.where(G > 0, G, np.nan), label="base input")
    ax[2].loglog(f, np.where(G > 0, transmissibility(f / fn, z) ** 2 * G, np.nan), label="isolated IMU")
    ax[2].set_xlabel("frequency (Hz)")
    ax[2].set_ylabel("PSD (g^2/Hz)")
    ax[2].set_title(f"chosen isolator, fn = {fn:.0f} Hz, zeta = {z}")
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out_dir / "vibration_isolation.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'vibration_isolation.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
