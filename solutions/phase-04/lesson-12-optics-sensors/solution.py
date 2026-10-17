"""P04.12 Optics & Sensors: reference solution.

Part one: the optics budget of a star tracker, the camera that tells a spacecraft where it is pointing
by photographing stars.

* reflection, refraction (Snell's law), total internal reflection and normal-incidence Fresnel loss;
* the thin lens: image distance, magnification, f-number;
* pixel scale and field of view of a detector behind a lens;
* the diffraction-limited spot (Airy radius 1.22 lambda N) against the pixel size;
* how many stars fall in the field (an empirical whole-sky count, labelled as an assumption) and how
  many photons a star delivers;
* centroiding: a Gaussian spot integrated over pixels, photon and read noise (seeded Monte Carlo,
  SYNTHETIC images), the 1/sqrt(N) law, pixel-locking bias for a too-sharp spot, and the attitude
  accuracy that results.

Part two: four sensor archetypes and the ADC behind them.

* the pitot-static probe (incompressible and subsonic compressible airspeed);
* the strain gauge in a Wheatstone bridge (quarter and full bridge, exact and linearised);
* the NTC thermistor (beta model) in a divider;
* the spring-mass accelerometer;
* ADC resolution in the physical unit, and the bits needed for a resolution.

Units: SI; angles in radians unless the name says deg or arcsec; magnitudes are astronomical (a
difference of 5 magnitudes is a factor 100 in flux).

Run:  python solution.py [out_dir]     (prints the verification tables, writes star_tracker.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-04/lesson-12-optics-sensors
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.special import erf

ARCSEC = np.pi / (180.0 * 3600.0)      # rad per arcsecond
G0 = 9.80665
GAMMA_AIR = 1.4
P0_SL = 101325.0
RHO0_SL = 1.225

# Assumed photometric zero point for a broad silicon band (about 400-900 nm): photons per second per
# square metre from a magnitude-0 star. An order-of-magnitude design value, not a calibrated standard.
PHOTON_FLUX_MAG0 = 1.0e10
# Assumed empirical whole-sky star count brighter than magnitude m: N(m) = A exp(B m). A rough fit used
# for star-tracker sizing; good to tens of per cent for 3 < m < 7 and not a catalogue.
STAR_COUNT_A, STAR_COUNT_B = 6.57, 1.08


# ---------------------------------------------------------------- reflection and refraction

def refraction_angle(n1: float, n2: float, theta1: np.ndarray | float) -> np.ndarray | float:
    """Snell's law n1 sin(theta1) = n2 sin(theta2): the refracted angle, or nan where there is total
    internal reflection (sin theta2 > 1)."""
    s = n1 * np.sin(np.asarray(theta1, dtype=float)) / n2
    with np.errstate(invalid="ignore"):
        return np.where(np.abs(s) <= 1.0, np.arcsin(np.clip(s, -1.0, 1.0)), np.nan)


def critical_angle(n1: float, n2: float) -> float:
    """Angle of incidence beyond which light going from n1 into n2 < n1 is totally reflected:
    arcsin(n2 / n1). Raises ValueError if n2 >= n1 (no total internal reflection)."""
    if n2 >= n1:
        raise ValueError("total internal reflection needs n2 < n1")
    return float(np.arcsin(n2 / n1))


def fresnel_reflectance_normal(n1: float, n2: float) -> float:
    """Fraction of power reflected at normal incidence: ((n1 - n2) / (n1 + n2))^2 (4 % for glass in air)."""
    return ((n1 - n2) / (n1 + n2)) ** 2


# ---------------------------------------------------------------- the thin lens and the camera

def thin_lens_image_distance(f: float, d_o: np.ndarray | float) -> np.ndarray | float:
    """Image distance from 1/d_o + 1/d_i = 1/f (real-is-positive convention). d_o = inf gives f;
    an object inside the focal length gives a negative (virtual) image distance."""
    d_o = np.asarray(d_o, dtype=float)
    with np.errstate(divide="ignore"):
        return 1.0 / (1.0 / f - 1.0 / d_o)


def magnification(d_o: float, d_i: float) -> float:
    """Lateral magnification m = -d_i / d_o (negative: inverted real image)."""
    return -d_i / d_o


def f_number(f: float, aperture_diameter: float) -> float:
    """N = f / D. Image irradiance scales as 1/N^2; the diffraction spot scales as N."""
    return f / aperture_diameter


def pixel_scale(pixel_pitch: float, f: float) -> float:
    """Angle subtended by one pixel at the centre of the field, pitch / f [rad] (small-angle)."""
    return pixel_pitch / f


def field_of_view(n_pixels: int, pixel_pitch: float, f: float) -> float:
    """Full field of view along one side of the detector: 2 arctan(n p / (2 f)) [rad]."""
    return 2.0 * np.arctan(n_pixels * pixel_pitch / (2.0 * f))


def airy_radius(wavelength: float, N: float) -> float:
    """Radius of the first dark ring of the diffraction pattern of a circular aperture in the focal
    plane: 1.22 lambda N. A perfect lens cannot focus a star to a smaller spot."""
    return 1.22 * wavelength * N


def airy_gaussian_sigma(wavelength: float, N: float) -> float:
    """Standard deviation of the Gaussian that best matches the Airy core: about 0.42 lambda N."""
    return 0.42 * wavelength * N


# ---------------------------------------------------------------- stars and photons

def stars_brighter_than(m: float) -> float:
    """Whole-sky count of stars brighter than visual magnitude m from the assumed fit A exp(B m)."""
    return STAR_COUNT_A * np.exp(STAR_COUNT_B * m)


def fov_solid_angle(fov_x: float, fov_y: float) -> float:
    """Solid angle of a rectangular field of view of full angles fov_x, fov_y [sr]:
    4 arcsin(sin(fov_x/2) sin(fov_y/2))."""
    return 4.0 * np.arcsin(np.sin(fov_x / 2.0) * np.sin(fov_y / 2.0))


def mean_stars_in_fov(m_limit: float, solid_angle: float) -> float:
    """Expected number of stars brighter than m_limit in a field of the given solid angle (stars taken
    as uniformly spread over the sky, which they are not: the galactic plane has several times more)."""
    return stars_brighter_than(m_limit) * solid_angle / (4.0 * np.pi)


def star_electrons(m: float, aperture_diameter: float, t_exp: float, qe: float = 0.5,
                   transmission: float = 0.8) -> float:
    """Photo-electrons collected from a star of magnitude m: flux0 10^(-0.4 m) x pi D^2/4 x t x QE x optics."""
    area = np.pi * aperture_diameter ** 2 / 4.0
    return PHOTON_FLUX_MAG0 * 10.0 ** (-0.4 * m) * area * t_exp * qe * transmission


# ---------------------------------------------------------------- centroiding (SYNTHETIC images)

def gaussian_spot(x0: float, y0: float, sigma: float, flux: float, size: int = 9) -> np.ndarray:
    """Expected electrons per pixel of a Gaussian spot of total ``flux`` centred at (x0, y0) pixels on a
    size x size window whose pixel (i, j) covers [j - 1/2, j + 1/2] x [i - 1/2, i + 1/2] in (x, y),
    with the window centre at pixel (size // 2, size // 2). Integrated exactly over each pixel with erf."""
    edges = np.arange(size + 1) - 0.5 - size // 2
    s2 = np.sqrt(2.0) * sigma
    fx = 0.5 * np.diff(erf((edges - x0) / s2))
    fy = 0.5 * np.diff(erf((edges - y0) / s2))
    return flux * np.outer(fy, fx)


def centroid(img: np.ndarray, background: float = 0.0) -> tuple[float, float]:
    """Centre of mass (x, y) of a window, in pixels relative to its central pixel, after subtracting a
    constant background."""
    w = np.asarray(img, dtype=float) - background
    size = w.shape[0]
    c = np.arange(size) - size // 2
    tot = w.sum()
    return float((w.sum(axis=0) @ c) / tot), float((w.sum(axis=1) @ c) / tot)


def simulate_centroid_errors(flux: float, sigma: float, read_noise: float = 0.0, n_trials: int = 500,
                             size: int = 9, seed: int = 0, subpixel: bool = True) -> np.ndarray:
    """SYNTHETIC Monte Carlo: n_trials star images with Poisson photon noise and Gaussian read noise,
    the true centre drawn uniformly within the central pixel (``subpixel``) or at its centre; returns the
    x-centroid errors in pixels."""
    rng = np.random.default_rng(seed)
    errs = np.empty(n_trials)
    for k in range(n_trials):
        x0, y0 = (rng.uniform(-0.5, 0.5, 2) if subpixel else (0.0, 0.0))
        lam = gaussian_spot(x0, y0, sigma, flux, size)
        img = rng.poisson(lam).astype(float) + rng.normal(0.0, read_noise, lam.shape)
        errs[k] = centroid(img)[0] - x0
    return errs


def centroid_sigma_photon_limit(sigma_psf: float, n_electrons: float) -> float:
    """Photon-noise-limited centroid error of a Gaussian spot: sigma_psf / sqrt(N) (pixels if sigma is)."""
    return sigma_psf / np.sqrt(n_electrons)


def pixel_locking_bias(sigma: float, n_pos: int = 101, size: int = 9) -> tuple[np.ndarray, np.ndarray]:
    """Noise-free centroid error against true sub-pixel position x0 in [-0.5, 0.5]: the systematic
    S-curve that appears when the spot is not much wider than a pixel. Returns (x0, bias)."""
    x0 = np.linspace(-0.5, 0.5, n_pos)
    bias = np.array([centroid(gaussian_spot(x, 0.0, sigma, 1.0, size))[0] - x for x in x0])
    return x0, bias


def attitude_error_cross_boresight(centroid_px: float, scale_rad_per_px: float, n_stars: float) -> float:
    """Cross-boresight attitude error from n_stars independent centroids: sigma_c x scale / sqrt(n) [rad]."""
    return centroid_px * scale_rad_per_px / np.sqrt(n_stars)


@dataclass
class StarTracker:
    """A star-tracker camera."""
    f: float = 25e-3                 # m focal length
    aperture: float = 25e-3 / 1.2    # m (f/1.2)
    n_pixels: int = 1024
    pitch: float = 5.5e-6            # m
    wavelength: float = 600e-9       # m, effective
    m_limit: float = 6.0
    t_exp: float = 0.1               # s

    @property
    def N(self) -> float:
        return f_number(self.f, self.aperture)

    @property
    def scale(self) -> float:
        return pixel_scale(self.pitch, self.f)

    @property
    def fov(self) -> float:
        return field_of_view(self.n_pixels, self.pitch, self.f)

    @property
    def n_stars(self) -> float:
        return mean_stars_in_fov(self.m_limit, fov_solid_angle(self.fov, self.fov))


# ---------------------------------------------------------------- sensor archetypes

def pitot_airspeed(delta_p: np.ndarray | float, rho: float = RHO0_SL) -> np.ndarray | float:
    """Incompressible Bernoulli: V = sqrt(2 (p_total - p_static) / rho). Valid below about Mach 0.3."""
    return np.sqrt(2.0 * np.asarray(delta_p, dtype=float) / rho)


def pitot_delta_p_compressible(mach: np.ndarray | float, p_static: float = P0_SL, gamma: float = GAMMA_AIR) -> np.ndarray | float:
    """Subsonic isentropic impact pressure: q_c = p ((1 + (gamma - 1)/2 M^2)^(gamma/(gamma - 1)) - 1)."""
    M = np.asarray(mach, dtype=float)
    return p_static * ((1.0 + 0.5 * (gamma - 1.0) * M * M) ** (gamma / (gamma - 1.0)) - 1.0)


def calibrated_airspeed(q_c: np.ndarray | float, gamma: float = GAMMA_AIR) -> np.ndarray | float:
    """Calibrated airspeed from impact pressure with sea-level standard constants (subsonic):
    V_c = a0 sqrt(2/(gamma - 1) ((q_c / p0 + 1)^((gamma - 1)/gamma) - 1)), a0 = sqrt(gamma p0 / rho0)."""
    a0 = np.sqrt(gamma * P0_SL / RHO0_SL)
    q = np.asarray(q_c, dtype=float)
    return a0 * np.sqrt(2.0 / (gamma - 1.0) * ((q / P0_SL + 1.0) ** ((gamma - 1.0) / gamma) - 1.0))


def bridge_output(V_ex: float, R1: float, R2: float, R3: float, R4: float) -> float:
    """Wheatstone bridge output V_out = V_ex (R3 / (R3 + R4) - R2 / (R1 + R2)): the difference between
    the mid-points of two dividers (R1 top, R2 bottom) and (R4 top, R3 bottom). Zero when R1 R3 = R2 R4."""
    return V_ex * (R3 / (R3 + R4) - R2 / (R1 + R2))


def quarter_bridge(V_ex: float, strain: np.ndarray | float, gauge_factor: float = 2.0, R: float = 350.0) -> np.ndarray | float:
    """One active gauge (R3 = R (1 + GF eps)) and three fixed resistors R: exact output
    V_ex (x / (2 (2 + x))), x = GF eps; for small x about V_ex GF eps / 4."""
    x = gauge_factor * np.asarray(strain, dtype=float)
    return V_ex * x / (2.0 * (2.0 + x))


def full_bridge(V_ex: float, strain: np.ndarray | float, gauge_factor: float = 2.0) -> np.ndarray | float:
    """Four active gauges, two in tension and two in compression on a bending beam: V_ex GF eps,
    exactly linear and four times the quarter bridge's small-signal sensitivity."""
    return V_ex * gauge_factor * np.asarray(strain, dtype=float)


def thermistor_resistance(T: np.ndarray | float, R0: float = 10e3, T0: float = 298.15, beta: float = 3950.0) -> np.ndarray | float:
    """NTC thermistor beta model R(T) = R0 exp(beta (1/T - 1/T0)), T in kelvin."""
    return R0 * np.exp(beta * (1.0 / np.asarray(T, dtype=float) - 1.0 / T0))


def thermistor_temperature(R: np.ndarray | float, R0: float = 10e3, T0: float = 298.15, beta: float = 3950.0) -> np.ndarray | float:
    """Inverse of the beta model: 1/T = 1/T0 + ln(R / R0) / beta."""
    return 1.0 / (1.0 / T0 + np.log(np.asarray(R, dtype=float) / R0) / beta)


def divider_voltage(V_s: float, R_top: float, R_bottom: np.ndarray | float) -> np.ndarray | float:
    """Output of a voltage divider measured across R_bottom."""
    Rb = np.asarray(R_bottom, dtype=float)
    return V_s * Rb / (R_top + Rb)


def accelerometer_deflection(accel: np.ndarray | float, f_n: float) -> np.ndarray | float:
    """Static deflection of a spring-mass accelerometer of natural frequency f_n: x = a / omega_n^2."""
    return np.asarray(accel, dtype=float) / (2.0 * np.pi * f_n) ** 2


def adc_resolution(full_scale_range: float, bits: int) -> float:
    """Smallest step of the measured quantity resolved by an ADC spanning full_scale_range: range / 2^bits."""
    return full_scale_range / 2 ** bits


def bits_needed(full_scale_range: float, resolution: float) -> int:
    """Smallest number of bits with range / 2^bits <= resolution."""
    return int(np.ceil(np.log2(full_scale_range / resolution)))


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Run the lesson's numbers, print them and save star_tracker.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)

    # ---- refraction
    print(f"glass n = 1.5 to air: critical angle {np.degrees(critical_angle(1.5, 1.0)):.2f} deg; normal reflectance "
          f"{fresnel_reflectance_normal(1.0, 1.5) * 100:.2f} %; 30 deg in air -> {np.degrees(refraction_angle(1.0, 1.5, np.radians(30))):.2f} deg in glass")
    print(f"fibre core 1.48 / cladding 1.46: critical angle {np.degrees(critical_angle(1.48, 1.46)):.2f} deg")
    print(f"thin lens f = 25 mm: object at 1 m -> image at {thin_lens_image_distance(25e-3, 1.0) * 1e3:.3f} mm; at infinity -> "
          f"{thin_lens_image_distance(25e-3, np.inf) * 1e3:.3f} mm; at 50 mm -> {thin_lens_image_distance(25e-3, 50e-3) * 1e3:.1f} mm, m = "
          f"{magnification(50e-3, thin_lens_image_distance(25e-3, 50e-3)):.2f}")

    # ---- the star tracker
    st = StarTracker()
    print(f"\nstar tracker: f = {st.f * 1e3:.0f} mm, f/{st.N:.1f} (D = {st.aperture * 1e3:.1f} mm), {st.n_pixels} x {st.n_pixels} px of {st.pitch * 1e6:.1f} um")
    print(f"  pixel scale {st.scale * 1e6:.0f} urad = {st.scale / ARCSEC:.1f} arcsec/px; FOV {np.degrees(st.fov):.2f} deg square, "
          f"{fov_solid_angle(st.fov, st.fov):.4f} sr")
    r_airy = airy_radius(st.wavelength, st.N)
    print(f"  Airy radius at {st.wavelength * 1e9:.0f} nm: {r_airy * 1e6:.2f} um = {r_airy / st.pitch:.2f} px "
          f"(Gaussian sigma {airy_gaussian_sigma(st.wavelength, st.N) / st.pitch:.3f} px): the focused star is smaller than a pixel")
    print(f"  stars brighter than m = {st.m_limit:.0f} over the sky (assumed fit): {stars_brighter_than(st.m_limit):.0f}; expected in FOV: {st.n_stars:.1f}")
    for m in (3.0, 5.0, 6.0):
        print(f"  m = {m:.0f} star, {st.t_exp * 1e3:.0f} ms: {star_electrons(m, st.aperture, st.t_exp):.0f} electrons")
    Ne = star_electrons(st.m_limit, st.aperture, st.t_exp)
    sig = 1.0
    errs0 = simulate_centroid_errors(Ne, sig, 0.0, 2000, seed=12)
    print(f"  SYNTHETIC centroiding, m = 6, defocused sigma = 1 px, 9 x 9 window, no read noise: rms {np.std(errs0):.4f} px "
          f"(sigma/sqrt(N) = {centroid_sigma_photon_limit(sig, Ne):.4f}; with the pixel's own 1/12 px^2: "
          f"{centroid_sigma_photon_limit(np.sqrt(sig ** 2 + 1 / 12), Ne):.4f})")
    print("  SYNTHETIC centroid rms (px) against window size and read noise (2000 images each, seed 12):")
    print("    read noise   3x3     5x5     7x7     9x9")
    for rn in (0.0, 3.0, 10.0):
        row = [np.sqrt(np.mean(simulate_centroid_errors(Ne, sig, rn, 2000, size=sz, seed=12) ** 2)) for sz in (3, 5, 7, 9)]
        print(f"    {rn:4.0f} e-    " + "  ".join(f"{v:.4f}" for v in row))
    errs = simulate_centroid_errors(Ne, sig, 3.0, 2000, size=5, seed=12)
    for s_ in (0.2, 0.5, 1.0):
        x0, b = pixel_locking_bias(s_)
        print(f"  pixel-locking bias, sigma = {s_} px: max |bias| {np.max(np.abs(b)):.4f} px")
    c_px = float(np.sqrt(np.mean(errs ** 2)))
    att = attitude_error_cross_boresight(c_px, st.scale, st.n_stars)
    print(f"  baseline 5 x 5 window, 3 e- read noise: {c_px:.4f} px = {c_px * st.scale / ARCSEC:.2f} arcsec per star; "
          f"cross-boresight attitude error with {st.n_stars:.0f} stars: {att / ARCSEC:.2f} arcsec")

    # ---- sensors
    print(f"\npitot: 0.6 kPa -> {pitot_airspeed(600.0):.2f} m/s; Mach 0.6 at sea level: q_c = {pitot_delta_p_compressible(0.6) / 1e3:.2f} kPa, "
          f"V_c = {calibrated_airspeed(pitot_delta_p_compressible(0.6)):.2f} m/s (true {0.6 * np.sqrt(1.4 * P0_SL / RHO0_SL):.2f}), "
          f"incompressible formula gives {pitot_airspeed(pitot_delta_p_compressible(0.6)):.2f} m/s")
    eps = 1000e-6
    print(f"strain gauge GF 2, 1000 microstrain, 3.3 V excitation: quarter bridge {quarter_bridge(3.3, eps) * 1e3:.4f} mV "
          f"(linear {3.3 * 2 * eps / 4 * 1e3:.4f} mV); full bridge {full_bridge(3.3, eps) * 1e3:.3f} mV")
    for TC in (-20.0, 25.0, 60.0):
        R = thermistor_resistance(TC + 273.15)
        print(f"thermistor 10 k / B 3950 at {TC:+.0f} C: {R:.0f} ohm, divider (10 k top, 3.3 V) {divider_voltage(3.3, 10e3, R):.3f} V")
    print(f"accelerometer f_n = 5 kHz: 1 g deflects {accelerometer_deflection(G0, 5e3) * 1e9:.2f} nm")
    print(f"12-bit ADC over +-16 g: {adc_resolution(32.0, 12) * 1000:.2f} mg per count; bits for 1 mg over 32 g: {bits_needed(32.0, 1e-3)}")

    # ---- figure
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))
    ax = axs[0]
    img = np.random.default_rng(3).poisson(gaussian_spot(0.23, -0.31, 1.0, Ne, 9)) + 0.0
    ax.imshow(img, origin="lower", extent=(-4.5, 4.5, -4.5, 4.5), cmap="gray")
    cx, cy = centroid(img)
    ax.plot(0.23, -0.31, "r+", ms=12, label="true"); ax.plot(cx, cy, "cx", ms=10, label="centroid")
    ax.set_title(f"SYNTHETIC m = 6 star, {Ne:.0f} e-"); ax.legend(fontsize=8)
    ax = axs[1]
    for s_ in (0.2, 0.35, 0.5, 1.0):
        x0, b = pixel_locking_bias(s_)
        ax.plot(x0, b, label=f"sigma = {s_} px")
    ax.set_xlabel("true sub-pixel position (px)"); ax.set_ylabel("centroid bias (px)"); ax.set_title("pixel locking (noise-free)"); ax.legend()
    ax = axs[2]
    fluxes = np.array([100, 300, 1000, 3000, 10000])
    rms = [np.std(simulate_centroid_errors(F, 1.0, 0.0, 400, seed=7)) for F in fluxes]
    ax.loglog(fluxes, rms, "o-", label="SYNTHETIC Monte Carlo")
    ax.loglog(fluxes, 1.0 / np.sqrt(fluxes), "--", label="sigma / sqrt(N)")
    ax.set_xlabel("electrons"); ax.set_ylabel("centroid rms (px)"); ax.set_title("photon-limited centroiding"); ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "star_tracker.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
