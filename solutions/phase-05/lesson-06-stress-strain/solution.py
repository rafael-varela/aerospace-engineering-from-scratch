"""P05.06 Stress & Strain: reference solution.

A tensile test pulls a coupon of known cross-section and gauge length and records force and
elongation. Reduced to engineering stress and strain, that record gives the numbers every
structural calculation in this phase starts from: Young's modulus E, the 0.2 % offset yield
strength, the ultimate tensile strength, the elongation after fracture and the energy absorbed
per unit volume (toughness). This module builds the reduction from scratch:

* engineering and true stress and strain, and the constant-volume link between them;
* the isotropic Hooke's law in three dimensions (the constitutive law), with the shear and bulk
  moduli it implies;
* the Ramberg-Osgood curve, whose 0.2 % offset yield is exact by construction (the check the
  analyzer must pass), and Considere's criterion, which says where the engineering curve peaks;
* a seeded synthetic tensile-test generator (the stand-in dataset of the lesson, labelled
  synthetic wherever it appears) with load-cell noise, extensometer noise and machine compliance;
* the analyzer: modulus by least squares, 0.2 % offset yield, ultimate strength, elongation,
  toughness and modulus of resilience.

Units: stress in MPa (= N/mm^2), strain dimensionless, force in N, lengths in mm, energy per unit
volume in MJ/m^3 (= MPa x strain). Material values are typical handbook values, not design
allowables.

Run:  python solution.py [out_dir]     (prints the lesson's tables, writes stress_strain.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-05/lesson-06-stress-strain
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

OFFSET = 0.002          # the 0.2 % offset used to define the yield strength of metals

# Typical room-temperature properties (MPa, dimensionless). Typical values for orientation and
# for the synthetic data; they are not design allowables (those are statistically based minima
# from a materials handbook and depend on product form, thickness and grain direction).
MATERIALS: dict[str, dict[str, float]] = {
    "2024-T3":   {"E": 73.1e3, "nu": 0.33, "Fty": 345.0, "Ftu": 483.0, "elongation": 0.18},
    "7075-T6":   {"E": 71.7e3, "nu": 0.33, "Fty": 503.0, "Ftu": 572.0, "elongation": 0.11},
    "Ti-6Al-4V": {"E": 113.8e3, "nu": 0.342, "Fty": 880.0, "Ftu": 950.0, "elongation": 0.14},
    "4130":      {"E": 205.0e3, "nu": 0.29, "Fty": 435.0, "Ftu": 670.0, "elongation": 0.255},
}


# ---------------------------------------------------------------- stress and strain

def engineering_stress(force: np.ndarray | float, area0: float) -> np.ndarray | float:
    """Engineering (nominal) stress sigma = F / A0 [MPa], F in N, A0 the original area in mm^2."""
    return np.asarray(force, dtype=float) / area0


def engineering_strain(elongation: np.ndarray | float, gauge_length: float) -> np.ndarray | float:
    """Engineering strain eps = dL / L0 (dimensionless), elongation and gauge length in mm."""
    return np.asarray(elongation, dtype=float) / gauge_length


def true_stress_strain(stress: np.ndarray, strain: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """True stress and true (logarithmic) strain from engineering values, assuming constant volume
    and uniform deformation (valid up to the ultimate strength, before necking):

        sigma_t = sigma (1 + eps),   eps_t = ln(1 + eps).
    """
    stress = np.asarray(stress, dtype=float)
    strain = np.asarray(strain, dtype=float)
    return stress * (1.0 + strain), np.log1p(strain)


# ---------------------------------------------------------------- Hooke's law in 3-D

def shear_modulus(E: float, nu: float) -> float:
    """Shear modulus of an isotropic material, G = E / (2 (1 + nu))."""
    return E / (2.0 * (1.0 + nu))


def bulk_modulus(E: float, nu: float) -> float:
    """Bulk modulus of an isotropic material, K = E / (3 (1 - 2 nu)); infinite at nu = 0.5."""
    return E / (3.0 * (1.0 - 2.0 * nu)) if nu < 0.5 else np.inf


def hooke_strain(stress: np.ndarray, E: float, nu: float) -> np.ndarray:
    """Small strain tensor (3x3) from the stress tensor (3x3) for an isotropic linear-elastic
    material:  eps = ((1 + nu) sigma - nu tr(sigma) I) / E.  Shear entries are tensor strains
    (half the engineering shear strain gamma)."""
    s = np.asarray(stress, dtype=float)
    return ((1.0 + nu) * s - nu * np.trace(s) * np.eye(3)) / E


def hooke_stress(strain: np.ndarray, E: float, nu: float) -> np.ndarray:
    """Inverse of :func:`hooke_strain`: sigma = lambda tr(eps) I + 2 mu eps, with the Lame constants
    lambda = E nu / ((1 + nu)(1 - 2 nu)) and mu = G."""
    e = np.asarray(strain, dtype=float)
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = shear_modulus(E, nu)
    return lam * np.trace(e) * np.eye(3) + 2.0 * mu * e


# ---------------------------------------------------------------- material models

def ramberg_osgood_strain(sigma: np.ndarray | float, E: float, sigma_y: float, n: float,
                          offset: float = OFFSET) -> np.ndarray | float:
    """Ramberg-Osgood strain: eps = sigma / E + offset (sigma / sigma_y)^n. At sigma = sigma_y the
    plastic part is exactly the offset, so the 0.2 % offset yield of this curve is sigma_y."""
    s = np.asarray(sigma, dtype=float)
    return s / E + offset * np.sign(s) * (np.abs(s) / sigma_y) ** n


def considere_point(E: float, sigma_y: float, n: float) -> tuple[float, float]:
    """Ultimate point of a material whose TRUE stress-strain curve is Ramberg-Osgood.

    The engineering stress sigma_t exp(-eps_t) peaks where d sigma_t / d eps_t = sigma_t
    (Considere's criterion). With eps_t(sigma_t) Ramberg-Osgood that is
    sigma_t / E + offset n (sigma_t / sigma_y)^n = 1, solved here by bisection (brentq).
    Returns (eps_u, sigma_u): engineering strain and engineering stress at the maximum.
    """
    f = lambda s: s / E + OFFSET * n * (s / sigma_y) ** n - 1.0
    s_c = brentq(f, 1e-9 * sigma_y, 50.0 * sigma_y, xtol=1e-12, rtol=1e-14)
    eps_t = s_c / E + OFFSET * (s_c / sigma_y) ** n
    return float(np.expm1(eps_t)), float(s_c * np.exp(-eps_t))


def synthetic_tensile_test(E: float, sigma_y: float, n: float, eps_fracture: float, drop: float = 0.10,
                           n_points: int = 2000, noise_stress: float = 0.0, noise_strain: float = 0.0,
                           compliance: float = 0.0, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """SYNTHETIC tensile test record (engineering strain, engineering stress in MPa).

    Up to the ultimate point the material follows a Ramberg-Osgood TRUE stress-strain curve,
    converted to engineering values at constant volume; the maximum therefore sits at the
    Considere point. After it, necking is imitated by a smooth parabolic fall of engineering
    stress by the fraction ``drop`` of the ultimate, ending at fracture at ``eps_fracture``.
    The record is sampled at ``n_points`` equal strain steps (a constant crosshead rate).

    Measurement imperfections: Gaussian noise of standard deviation ``noise_stress`` (MPa, load
    cell) and ``noise_strain`` (extensometer), and a machine ``compliance`` (1/MPa) that adds
    compliance x sigma to the recorded strain, as happens when strain is taken from crosshead
    travel instead of an extensometer. Seeded with numpy.random.default_rng(seed).
    """
    eps_u, sig_u = considere_point(E, sigma_y, n)
    if eps_fracture <= eps_u:
        raise ValueError(f"eps_fracture {eps_fracture} must exceed the uniform strain {eps_u:.4f}")
    s_c = sig_u * (1.0 + eps_u)                                  # true stress at the maximum
    st = np.linspace(0.0, s_c, 40001)
    et = st / E + OFFSET * (st / sigma_y) ** n
    eps_pre, sig_pre = np.expm1(et), st * np.exp(-et)
    eps = np.linspace(0.0, eps_fracture, n_points)
    post = sig_u * (1.0 - drop * ((eps - eps_u) / (eps_fracture - eps_u)) ** 2)
    sig = np.where(eps <= eps_u, np.interp(eps, eps_pre, sig_pre), post)
    rng = np.random.default_rng(seed)
    sig_meas = sig + noise_stress * rng.standard_normal(n_points)
    eps_meas = eps + compliance * sig + noise_strain * rng.standard_normal(n_points)
    return eps_meas, sig_meas


# ---------------------------------------------------------------- the analyzer

def fit_modulus(strain: np.ndarray, stress: np.ndarray, lo: float = 0.10, hi: float = 0.40) -> tuple[float, float]:
    """Young's modulus by least squares on the points before the maximum whose stress lies between
    lo and hi times the maximum stress. Returns (E, eps0) where eps0 is the strain intercept at zero
    stress (the 'toe' of the curve: slack, seating or an extensometer zero offset)."""
    strain = np.asarray(strain, dtype=float)
    stress = np.asarray(stress, dtype=float)
    i_max = int(np.argmax(stress))
    smax = stress[i_max]
    sel = np.zeros(stress.size, bool)
    sel[:i_max + 1] = True
    sel &= (stress >= lo * smax) & (stress <= hi * smax)
    if sel.sum() < 3:
        raise ValueError("fewer than three points in the modulus window")
    slope, intercept = np.polyfit(strain[sel], stress[sel], 1)
    return float(slope), float(-intercept / slope)


def offset_yield(strain: np.ndarray, stress: np.ndarray, E: float, eps0: float = 0.0,
                 offset: float = OFFSET) -> tuple[float, float]:
    """0.2 % offset yield: where the curve first crosses the line of slope E through
    (eps0 + offset, 0). Searches before the maximum; the crossing is linearly interpolated.
    Returns (strain, stress) at yield."""
    strain = np.asarray(strain, dtype=float)
    stress = np.asarray(stress, dtype=float)
    i_max = int(np.argmax(stress))
    f = stress[:i_max + 1] - E * (strain[:i_max + 1] - eps0 - offset)
    idx = np.nonzero((f[:-1] > 0.0) & (f[1:] <= 0.0))[0]
    if idx.size == 0:
        raise ValueError("the offset line does not cross the curve before the maximum")
    i = int(idx[0])
    w = f[i] / (f[i] - f[i + 1])
    eps_y = strain[i] + w * (strain[i + 1] - strain[i])
    return float(eps_y), float(stress[i] + w * (stress[i + 1] - stress[i]))


def moving_average(x: np.ndarray, window: int) -> np.ndarray:
    """Centred moving average of odd length ``window`` (edges use the points available)."""
    x = np.asarray(x, dtype=float)
    if window <= 1:
        return x.copy()
    k = np.ones(window)
    return np.convolve(x, k, "same") / np.convolve(np.ones_like(x), k, "same")


def ultimate_strength(strain: np.ndarray, stress: np.ndarray, window: int = 1) -> tuple[float, float]:
    """Ultimate tensile strength: the maximum engineering stress, after an optional centred moving
    average of ``window`` points to keep load-cell noise from biasing the maximum upward.
    Returns (strain at the maximum, maximum stress)."""
    s = moving_average(stress, window)
    i = int(np.argmax(s))
    return float(np.asarray(strain)[i]), float(s[i])


def elongation_after_fracture(strain: np.ndarray, stress: np.ndarray, E: float, eps0: float = 0.0) -> float:
    """Permanent (plastic) strain after fracture: the last recorded strain minus the elastic
    recovery sigma_f / E and the toe offset. This is what fitting the broken halves back together
    measures."""
    return float(np.asarray(strain)[-1] - np.asarray(stress)[-1] / E - eps0)


def toughness(strain: np.ndarray, stress: np.ndarray) -> float:
    """Energy absorbed per unit volume up to fracture, the area under the engineering curve
    (trapezoid rule): MPa x strain = MJ/m^3."""
    return float(np.trapezoid(stress, strain))


def resilience(sigma_y: float, E: float) -> float:
    """Modulus of resilience, the elastic energy per unit volume at yield, sigma_y^2 / (2 E) [MJ/m^3]."""
    return sigma_y ** 2 / (2.0 * E)


def remove_compliance(strain: np.ndarray, stress: np.ndarray, compliance: float) -> np.ndarray:
    """Correct strain taken from crosshead travel: subtract the machine's compliance x stress."""
    return np.asarray(strain, dtype=float) - compliance * np.asarray(stress, dtype=float)


@dataclass
class TensileResult:
    """Reduced tensile-test properties (MPa, strain, MJ/m^3)."""
    E: float
    eps0: float
    eps_y: float
    sigma_y: float
    eps_u: float
    sigma_u: float
    elongation: float
    toughness: float
    resilience: float


def analyze_tensile(strain: np.ndarray, stress: np.ndarray, lo: float = 0.10, hi: float = 0.40,
                    window: int = 1) -> TensileResult:
    """Reduce one engineering stress-strain record to its properties."""
    E, eps0 = fit_modulus(strain, stress, lo, hi)
    eps_y, sig_y = offset_yield(strain, stress, E, eps0)
    eps_u, sig_u = ultimate_strength(strain, stress, window)
    return TensileResult(E, eps0, eps_y, sig_y, eps_u, sig_u,
                         elongation_after_fracture(strain, stress, E, eps0),
                         toughness(strain, stress), resilience(sig_y, E))


# ---------------------------------------------------------------- main

COUPON = {"E": 73.1e3, "sigma_y": 345.0, "n": 9.0, "eps_fracture": 0.18, "drop": 0.10}   # synthetic 2024-T3


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save stress_strain.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)

    print("typical properties (not design allowables)")
    print("material      E (GPa)   nu     G (GPa)  Fty (MPa)  Ftu (MPa)  elong")
    for name, m in MATERIALS.items():
        print(f"{name:12s} {m['E'] / 1e3:7.1f}  {m['nu']:.3f}  {shear_modulus(m['E'], m['nu']) / 1e3:7.1f}"
              f"  {m['Fty']:9.0f}  {m['Ftu']:9.0f}  {m['elongation']:.3f}")

    # ---- a 2024-T3 strap in tension, and a pressurised fuselage skin (Hooke's law in 2-D)
    al = MATERIALS["2024-T3"]
    F, b, t, L0 = 20e3, 25.0, 3.0, 200.0
    sig = engineering_stress(F, b * t)
    eps = sig / al["E"]
    print(f"\nstrap {b:.0f} x {t:.0f} mm, {F / 1e3:.0f} kN: sigma = {sig:.1f} MPa, eps = {eps:.3e}, "
          f"elongation over {L0:.0f} mm = {eps * L0:.3f} mm, width change = {-al['nu'] * eps * b * 1e3:.1f} um, "
          f"stress / Fty = {sig / al['Fty']:.3f}")
    p, R, ts = 0.0552, 1975.0, 1.6                       # MPa, mm, mm
    s_hoop, s_axial = p * R / ts, p * R / (2 * ts)
    e = hooke_strain(np.diag([s_hoop, s_axial, 0.0]), al["E"], al["nu"])
    print(f"fuselage skin R = {R / 1e3:.3f} m, t = {ts} mm, p = {p * 1e3:.1f} kPa: hoop {s_hoop:.1f} MPa, axial "
          f"{s_axial:.1f} MPa; strains hoop {e[0, 0]:.3e}, axial {e[1, 1]:.3e}, thickness {e[2, 2]:.3e}; "
          f"radius grows {e[0, 0] * R:.2f} mm (uniaxial hoop estimate {s_hoop / al['E'] * R:.2f} mm)")

    # ---- Ramberg-Osgood: the offset yield is exact
    sg = np.linspace(0, 500.0, 5001)
    ey, sy = offset_yield(ramberg_osgood_strain(sg, COUPON["E"], COUPON["sigma_y"], COUPON["n"]), sg, COUPON["E"])
    print(f"\nRamberg-Osgood check: offset yield of the exact curve {sy:.4f} MPa (input {COUPON['sigma_y']:.1f})")
    eu, su = considere_point(COUPON["E"], COUPON["sigma_y"], COUPON["n"])
    print(f"Considere point of the synthetic 2024-T3 material: eps_u = {eu:.4f}, sigma_u = {su:.1f} MPa")

    # ---- the synthetic coupon, clean and with noise
    clean = synthetic_tensile_test(**COUPON)
    noisy = synthetic_tensile_test(**COUPON, noise_stress=1.0, noise_strain=2e-5, seed=2024)
    print("\nSYNTHETIC 2024-T3 coupon     E (GPa)  Fty (MPa)  Ftu (MPa)  eps_u   elong   toughness  resilience")
    for label, (x, y), w in (("clean", clean, 1), ("noisy, seed 2024", noisy, 21)):
        r = analyze_tensile(x, y, window=w)
        print(f"  {label:26s} {r.E / 1e3:7.2f}  {r.sigma_y:9.1f}  {r.sigma_u:9.1f}  {r.eps_u:.4f}  {r.elongation:.4f}"
              f"  {r.toughness:6.1f}     {r.resilience:.3f}")

    print("Ftu of the noisy record against the smoothing window: " + ", ".join(
        f"{w} pt {ultimate_strength(*noisy, window=w)[1]:.1f}" for w in (1, 5, 21, 51)) + " MPa")

    # ---- experiment 1: modulus window
    print("\nmodulus window (noisy coupon)   E (GPa)   Fty (MPa)")
    for lo, hi in ((0.05, 0.20), (0.10, 0.40), (0.20, 0.60), (0.30, 0.80)):
        r = analyze_tensile(*noisy, lo=lo, hi=hi, window=21)
        print(f"  {lo:.2f}-{hi:.2f} of max            {r.E / 1e3:7.2f}   {r.sigma_y:7.1f}")

    # ---- experiment 2: noise Monte Carlo
    print("\nnoise study, 200 seeds          mean E   sd E (GPa)  mean Fty  sd Fty (MPa)")
    for ns, ne in ((0.5, 1e-5), (1.0, 2e-5), (3.0, 5e-5)):
        Es, Ys = [], []
        for s in range(200):
            r = analyze_tensile(*synthetic_tensile_test(**COUPON, noise_stress=ns, noise_strain=ne, seed=s), window=21)
            Es.append(r.E)
            Ys.append(r.sigma_y)
        print(f"  {ns:.1f} MPa, {ne:.0e}           {np.mean(Es) / 1e3:7.2f}  {np.std(Es) / 1e3:6.2f}"
              f"      {np.mean(Ys):7.1f}  {np.std(Ys):6.2f}")

    # ---- experiment 3: crosshead compliance
    print("\ncrosshead compliance   apparent E  Fty with      Fty with      elong   corrected")
    print("  (1/MPa)               (GPa)       fitted E      handbook E            E (GPa)")
    for c in (0.0, 2e-6, 5e-6, 1e-5):
        x, y = synthetic_tensile_test(**COUPON, compliance=c)
        r = analyze_tensile(x, y)
        _, sy_hb = offset_yield(x, y, MATERIALS["2024-T3"]["E"], r.eps0)
        rc = analyze_tensile(remove_compliance(x, y, c), y)
        print(f"  {c:.0e}               {r.E / 1e3:7.2f}     {r.sigma_y:7.1f}       {sy_hb:7.1f}"
              f"     {r.elongation:.4f}  {rc.E / 1e3:7.2f}")

    # ---- figure
    r = analyze_tensile(*noisy, window=21)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].plot(noisy[0] * 100, noisy[1], lw=0.8, label="synthetic 2024-T3 record (seed 2024)")
    st, et = true_stress_strain(clean[1], clean[0])
    m = clean[0] <= r.eps_u
    ax[0].plot(et[m] * 100, st[m], "--", lw=1, label="true stress (to the maximum)")
    xx = np.linspace(r.eps0 + OFFSET, r.eps_y + 0.002, 50)
    ax[0].plot(xx * 100, r.E * (xx - r.eps0 - OFFSET), "k:", lw=1, label="0.2 % offset line")
    ax[0].plot(r.eps_y * 100, r.sigma_y, "o", label=f"Fty {r.sigma_y:.0f} MPa")
    ax[0].plot(r.eps_u * 100, r.sigma_u, "s", label=f"Ftu {r.sigma_u:.0f} MPa")
    ax[0].set_xlabel("strain (%)"); ax[0].set_ylabel("stress (MPa)"); ax[0].legend(fontsize=8)
    ax[0].set_title("SYNTHETIC tensile record and its reduction")
    ax[1].plot(noisy[0] * 100, noisy[1], ".", ms=2)
    ax[1].plot(xx * 100, r.E * (xx - r.eps0 - OFFSET), "k:", lw=1)
    ax[1].set_xlim(0, 1.2); ax[1].set_ylim(0, 420)
    ax[1].set_xlabel("strain (%)"); ax[1].set_title("elastic region and offset construction")
    fig.tight_layout()
    fig.savefig(out_dir / "stress_strain.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'stress_strain.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
