"""P04.11 Electromagnetic Waves & Antennas: starter code.

From Maxwell's equations to a CubeSat downlink:

* the speed of light from mu0 and eps0; a plane wave's E/B ratio, its Poynting intensity and the
  inverse-square law for a point source;
* radiation pressure, and the thrust of a solar sail;
* the far-field patterns of the short (Hertzian) dipole, the half-wave dipole and a dipole of any
  length, the last also computed numerically by summing Hertzian elements with a sinusoidal current;
* directivity by integrating the pattern over the sphere (``scipy.integrate.quad`` as the oracle),
  half-power beamwidth, radiation resistance, gain = efficiency x directivity, effective aperture;
* polarisation mismatch between linear and circular antennas;
* the Friis transmission equation and a link budget for a 437 MHz CubeSat downlink at 500 km,
  with slant range from the elevation angle, thermal noise kTB and the Eb/N0 margin.

Units: SI; powers also in dBW/dBm, gains in dBi (decibels relative to an isotropic antenna).

Run:  python main.py [out_dir]     (prints the verification tables, writes dipole_pattern.png)
Test: pytest curriculum/phase-04/lesson-11-electromagnetic-waves-antennas
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import quad

MU0 = 4.0e-7 * np.pi                  # T m / A
EPS0 = 8.8541878128e-12               # F / m
C_LIGHT = 299_792_458.0               # m / s (exact by definition of the metre)
ETA0 = MU0 * C_LIGHT                  # ohm, impedance of free space, about 376.73
K_BOLTZMANN = 1.380649e-23            # J / K
SOLAR_CONSTANT = 1361.0               # W / m^2 at 1 au
AU = 1.495978707e11                   # m
R_EARTH = 6.371e6                     # m


# ---------------------------------------------------------------- plane waves

def wave_speed(mu: float = MU0, eps: float = EPS0) -> float:
    """Speed of an electromagnetic wave from Maxwell's equations: 1 / sqrt(mu eps)."""
    raise NotImplementedError


def b_from_e(E: float, c: float = C_LIGHT) -> float:
    """Magnetic amplitude of a plane wave with electric amplitude E: B = E / c."""
    raise NotImplementedError


def plane_wave_intensity(E0: float) -> float:
    """Time-averaged Poynting intensity of a plane wave of peak field E0: I = E0^2 / (2 eta0) = c eps0 E0^2 / 2."""
    raise NotImplementedError


def e_from_intensity(I: float) -> float:
    """Peak electric field of a plane wave of intensity I: sqrt(2 eta0 I)."""
    raise NotImplementedError


def isotropic_intensity(P: float, r: np.ndarray | float) -> np.ndarray | float:
    """Intensity at distance r of power P radiated equally in all directions: P / (4 pi r^2)."""
    raise NotImplementedError


def radiation_pressure(I: float, reflectivity: float = 0.0) -> float:
    """Pressure of light of intensity I on a flat surface facing it: (1 + reflectivity) I / c.
    0 for a perfect absorber gives I/c; 1 for a perfect mirror gives 2 I/c."""
    raise NotImplementedError


def sail_acceleration(area: float, mass: float, r_au: float = 1.0, reflectivity: float = 1.0,
                      cone_deg: float = 0.0) -> float:
    """Acceleration of a flat sail facing the Sun at angle ``cone_deg`` (between the normal and the
    sunline): a = (1 + rho) I(r) A cos^2(alpha) / (m c), I(r) = SOLAR_CONSTANT / r_au^2.
    (For a perfect mirror the force is along the normal; one cos for the projected area, one for the
    normal component of the momentum change. For rho < 1 this is the normal component only.)"""
    raise NotImplementedError


# ---------------------------------------------------------------- dipole patterns

def hertzian_pattern(theta: np.ndarray | float) -> np.ndarray | float:
    """Normalised radiation intensity of an electrically short dipole along z: sin^2(theta)."""
    raise NotImplementedError


def dipole_pattern(theta: np.ndarray | float, length_wl: float) -> np.ndarray | float:
    """Unnormalised radiation intensity of a centre-fed thin dipole of total length L = length_wl
    wavelengths with a sinusoidal current: F(theta) = [(cos(kL/2 cos t) - cos(kL/2)) / sin t]^2.
    For L = lambda/2 this is [cos(pi/2 cos t) / sin t]^2. Zero on the axis (limit taken)."""
    raise NotImplementedError


def halfwave_pattern(theta: np.ndarray | float) -> np.ndarray | float:
    """Normalised half-wave dipole intensity [cos(pi/2 cos theta) / sin theta]^2 (max 1 at broadside)."""
    raise NotImplementedError


def dipole_pattern_numeric(theta: np.ndarray, length_wl: float, n_seg: int = 200) -> np.ndarray:
    """The same far-field intensity computed by superposition: the dipole is cut into n_seg Hertzian
    elements at z_k with current I(z) = sin(k (L/2 - |z|)); each contributes I_k dz exp(j k z_k cos t),
    and the sum is multiplied by the element factor sin t. Returns |sum|^2 scaled so that it equals
    ``dipole_pattern`` (which differs from the integral by the constant factor (2/k)^2 ... see docstring of
    the test). Proves the closed form is the array sum of short dipoles."""
    raise NotImplementedError


def directivity(pattern, n_theta: int = 2001) -> float:
    """Maximum directivity D = 4 pi U_max / P_rad of an axisymmetric pattern U(theta), with
    P_rad = 2 pi int_0^pi U sin(theta) d theta, by the trapezoid rule on n_theta points."""
    raise NotImplementedError


def directivity_quad(pattern, theta_max: float = np.pi / 2) -> float:
    """The same with ``scipy.integrate.quad`` (the oracle): U_max taken at theta_max."""
    raise NotImplementedError


def half_power_beamwidth(pattern, n_theta: int = 180001) -> float:
    """Full width in degrees between the two directions where an axisymmetric pattern U(theta), maximal
    at broadside, falls to half its maximum."""
    raise NotImplementedError


def radiation_resistance_short(length_wl: float) -> float:
    """Radiation resistance of a short dipole with uniform current: R = (2 pi / 3) eta0 (L / lambda)^2,
    about 80 pi^2 (L / lambda)^2 ohm (with the triangular current of a real centre-fed short dipole it is a
    quarter of this)."""
    raise NotImplementedError


def radiation_resistance_dipole(length_wl: float) -> float:
    """Radiation resistance referred to the current maximum, R = (eta0 / (2 pi)) int_0^pi F(theta) sin(theta) d theta
    with F from ``dipole_pattern`` (73.1 ohm for the half-wave dipole), by ``quad``. From P = I0^2 R / 2 with
    the far field |E| = eta0 I0 sqrt(F) / (2 pi r)."""
    raise NotImplementedError


def to_db(x: np.ndarray | float) -> np.ndarray | float:
    """10 log10 x (power ratios)."""
    raise NotImplementedError


def from_db(x_db: np.ndarray | float) -> np.ndarray | float:
    """10^(x / 10)."""
    raise NotImplementedError


def gain(directivity_value: float, efficiency: float) -> float:
    """Gain = radiation efficiency x directivity (an antenna that dissipates power in its own
    resistance has less gain than its pattern alone suggests)."""
    raise NotImplementedError


def effective_aperture(G: float, wavelength: float) -> float:
    """Effective area of a receiving antenna of gain G: A_e = G lambda^2 / (4 pi). An isotropic
    antenna has A_e = lambda^2 / (4 pi) whatever its size: small at short wavelength."""
    raise NotImplementedError


def dish_gain(diameter: float, wavelength: float, efficiency: float = 0.55) -> float:
    """Gain of a circular aperture (parabolic dish): eta (pi D / lambda)^2, from A_e = eta pi D^2 / 4."""
    raise NotImplementedError


def polarization_loss(tx: str, rx: str, angle_deg: float = 0.0) -> float:
    """Polarisation loss factor (power fraction received) for ideal antennas:
    linear-linear cos^2(angle between them); linear-circular 1/2 at any angle; circular-circular of the
    same hand 1, opposite hand 0. ``tx``/``rx`` are 'linear', 'rhcp' or 'lhcp'."""
    raise NotImplementedError


# ---------------------------------------------------------------- the link budget

def free_space_path_loss(distance: np.ndarray | float, wavelength: float) -> np.ndarray | float:
    """(4 pi d / lambda)^2, the Friis spreading factor as a loss (> 1)."""
    raise NotImplementedError


def friis_received_power(P_t: float, G_t: float, G_r: float, distance: np.ndarray | float,
                         wavelength: float) -> np.ndarray | float:
    """Friis transmission equation, P_r = P_t G_t G_r (lambda / (4 pi d))^2, for antennas in each other's
    far field, matched in polarisation and impedance."""
    raise NotImplementedError


def slant_range(altitude: float, elevation_deg: np.ndarray | float, R: float = R_EARTH) -> np.ndarray | float:
    """Distance from a ground station to a satellite at ``altitude`` seen at elevation angle e over a
    spherical Earth: d = sqrt((R + h)^2 - R^2 cos^2 e) - R sin e. Equals h at zenith."""
    raise NotImplementedError


@dataclass
class Link:
    """A downlink budget; every gain and loss as a linear power ratio except where named _dB."""
    frequency: float = 437e6          # Hz (UHF amateur-satellite band)
    P_t: float = 1.0                  # W transmitter output
    G_t_dBi: float = 2.15             # half-wave dipole
    G_r_dBi: float = 14.0             # ground Yagi
    altitude: float = 500e3           # m
    elevation_deg: float = 10.0
    other_losses_dB: float = 3.0 + 1.0 + 2.0   # polarisation (linear to circular) + cable + pointing/atmosphere
    T_sys: float = 500.0              # K, system noise temperature at UHF (sky, ground, receiver)
    bit_rate: float = 9600.0          # bit/s
    required_EbN0_dB: float = 10.0    # for the chosen modulation and coding (assumption)

    @property
    def wavelength(self) -> float:
        return C_LIGHT / self.frequency

    @property
    def distance(self) -> float:
        return float(slant_range(self.altitude, self.elevation_deg))

    def budget(self) -> dict:
        """Received power, noise density, Eb/N0 and margin, in dB units."""
        Pt_dBW = to_db(self.P_t)
        fspl_dB = to_db(free_space_path_loss(self.distance, self.wavelength))
        Pr_dBW = Pt_dBW + self.G_t_dBi + self.G_r_dBi - fspl_dB - self.other_losses_dB
        N0_dBW_Hz = to_db(K_BOLTZMANN * self.T_sys)
        EbN0_dB = Pr_dBW - N0_dBW_Hz - to_db(self.bit_rate)
        return dict(Pt_dBW=Pt_dBW, fspl_dB=float(fspl_dB), Pr_dBW=float(Pr_dBW), Pr_dBm=float(Pr_dBW + 30.0),
                    N0_dBW_Hz=float(N0_dBW_Hz), EbN0_dB=float(EbN0_dB),
                    margin_dB=float(EbN0_dB - self.required_EbN0_dB), distance=self.distance)


# ---------------------------------------------------------------- the Hertzian dipole field and the inverse square

def hertzian_far_field(I0: float, dl: float, wavelength: float, r: float, theta: np.ndarray) -> np.ndarray:
    """Peak far-field E_theta [V/m] of a short dipole of length dl carrying peak current I0:
    |E| = eta0 k I0 dl sin(theta) / (4 pi r)."""
    raise NotImplementedError


def hertzian_radiated_power(I0: float, dl: float, wavelength: float) -> float:
    """Total power of a short dipole with uniform peak current I0: P = (1/2) I0^2 R_rad,
    R_rad = (2 pi / 3) eta0 (dl/lambda)^2 (equivalently P = eta0 k^2 I0^2 dl^2 / (12 pi))."""
    raise NotImplementedError


def power_through_sphere(I0: float, dl: float, wavelength: float, r: float, n: int = 2001) -> float:
    """Integrate the Poynting intensity |E|^2 / (2 eta0) of ``hertzian_far_field`` over a sphere of
    radius r: the same at every r (energy conservation is the inverse-square law)."""
    raise NotImplementedError


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Run the lesson's numbers, print them and save dipole_pattern.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)

    # ---- plane waves and the Sun
    print(f"1/sqrt(mu0 eps0) = {wave_speed():.6e} m/s (c = {C_LIGHT:.6e}); eta0 = {ETA0:.2f} ohm")
    E_sun = e_from_intensity(SOLAR_CONSTANT)
    print(f"sunlight at 1 au, {SOLAR_CONSTANT:.0f} W/m^2: E0 = {E_sun:.0f} V/m, B0 = {b_from_e(E_sun) * 1e6:.2f} uT; "
          f"pressure absorbed {radiation_pressure(SOLAR_CONSTANT, 0) * 1e6:.2f} uPa, reflected {radiation_pressure(SOLAR_CONSTANT, 1) * 1e6:.2f} uPa")
    A_sail, m_sail = 32.0, 5.0
    a1 = sail_acceleration(A_sail, m_sail)
    print(f"sail {A_sail:.0f} m^2, {m_sail:.0f} kg, perfect mirror at 1 au: force {a1 * m_sail * 1e6:.0f} uN, "
          f"a = {a1:.2e} m/s^2, dv per day {a1 * 86400:.2f} m/s; at 35 deg cone {sail_acceleration(A_sail, m_sail, cone_deg=35):.2e} m/s^2; "
          f"at 1.52 au {sail_acceleration(A_sail, m_sail, 1.52):.2e} m/s^2")
    print(f"1 W isotropic at 500 km: {isotropic_intensity(1.0, 500e3):.2e} W/m^2, at 1000 km: {isotropic_intensity(1.0, 1000e3):.2e} W/m^2")

    # ---- dipoles
    D_h, D_hq = directivity(hertzian_pattern), directivity_quad(hertzian_pattern)
    D_hw, D_hwq = directivity(halfwave_pattern), directivity_quad(halfwave_pattern)
    print(f"\nshort dipole: D = {D_h:.5f} (quad {D_hq:.5f}, exact 1.5) = {to_db(D_h):.2f} dBi, HPBW {half_power_beamwidth(hertzian_pattern):.2f} deg")
    print(f"half-wave dipole: D = {D_hw:.5f} (quad {D_hwq:.5f}) = {to_db(D_hw):.2f} dBi, HPBW {half_power_beamwidth(halfwave_pattern):.2f} deg, "
          f"R_rad = {radiation_resistance_dipole(0.5):.2f} ohm")
    t = np.linspace(0.01, np.pi - 0.01, 7)
    err = np.max(np.abs(dipole_pattern_numeric(t, 0.5) - dipole_pattern(t, 0.5)))
    print(f"half-wave pattern by summing 200 Hertzian elements: max |difference| {err:.1e}")
    for Lw in (0.1, 0.5, 1.0, 1.25, 1.5):
        f = lambda th, L=Lw: dipole_pattern(th, L)
        th = np.linspace(1e-4, np.pi - 1e-4, 20001)
        print(f"  L = {Lw:4.2f} lambda: D = {directivity(f):.3f} ({to_db(directivity(f)):.2f} dBi), max at {np.degrees(th[np.argmax(f(th))]):.1f} deg")
    print(f"short dipole R_rad (uniform current) at L = 0.05 lambda: {radiation_resistance_short(0.05):.2f} ohm")

    # ---- inverse square from the Hertzian field
    lam = C_LIGHT / 437e6
    for r in (10.0, 1e3, 1e5):
        print(f"power through a sphere of r = {r:.0e} m: {power_through_sphere(1.0, 0.05 * lam, lam, r):.6f} W "
              f"(radiated {hertzian_radiated_power(1.0, 0.05 * lam, lam):.6f} W)")

    # ---- apertures
    for name, G, lamx in (("isotropic, 437 MHz", 1.0, lam), ("half-wave dipole, 437 MHz", D_hw, lam),
                          ("isotropic, 2.25 GHz", 1.0, C_LIGHT / 2.25e9), ("3 m dish, 2.25 GHz", dish_gain(3.0, C_LIGHT / 2.25e9), C_LIGHT / 2.25e9)):
        print(f"  {name}: G = {to_db(G):.1f} dBi, A_e = {effective_aperture(G, lamx):.4f} m^2")

    # ---- the link budget
    link = Link()
    b = link.budget()
    print(f"\nCubeSat downlink {link.frequency / 1e6:.0f} MHz, {link.P_t} W, {link.G_t_dBi} dBi -> {link.G_r_dBi} dBi, "
          f"{link.altitude / 1e3:.0f} km at {link.elevation_deg:.0f} deg elevation: slant range {b['distance'] / 1e3:.0f} km")
    print(f"  FSPL {b['fspl_dB']:.1f} dB, other losses {link.other_losses_dB:.1f} dB, P_r = {b['Pr_dBm']:.1f} dBm, "
          f"N0 = {b['N0_dBW_Hz'] + 30:.1f} dBm/Hz, Eb/N0 = {b['EbN0_dB']:.1f} dB at {link.bit_rate:.0f} bit/s, margin {b['margin_dB']:.1f} dB")
    for e in (5.0, 10.0, 30.0, 90.0):
        bb = Link(elevation_deg=e).budget()
        print(f"  elevation {e:4.0f} deg: range {bb['distance'] / 1e3:6.0f} km, FSPL {bb['fspl_dB']:.1f} dB, margin {bb['margin_dB']:5.1f} dB")
    print(f"  polarisation: linear-linear at 90 deg {polarization_loss('linear', 'linear', 90):.2f}, linear-circular "
          f"{to_db(polarization_loss('linear', 'rhcp')):.2f} dB, RHCP-LHCP {polarization_loss('rhcp', 'lhcp'):.0f}")
    print(f"  highest bit rate with margin >= 3 dB at 10 deg: {link.bit_rate * from_db(b['margin_dB'] - 3.0):.0f} bit/s")

    # ---- figure
    fig = plt.figure(figsize=(14, 4.6))
    th = np.linspace(0, 2 * np.pi, 721)
    ax = fig.add_subplot(1, 3, 1, projection="polar")
    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    ax.plot(th, hertzian_pattern(th), label="short dipole")
    ax.plot(th, halfwave_pattern(th), label="half-wave")
    ax.plot(th, dipole_pattern(th, 1.25) / dipole_pattern(np.pi / 2, 1.25) * 1.0, label="1.25 lambda (norm.)", lw=0.8)
    ax.set_title("power patterns (dipole along 0 deg)"); ax.legend(loc="lower left", fontsize=7)
    ax = fig.add_subplot(1, 3, 2)
    th2 = np.linspace(0.01, np.pi - 0.01, 181)
    ax.plot(np.degrees(th2), to_db(halfwave_pattern(th2) + 1e-12), label="closed form")
    ax.plot(np.degrees(th2[::10]), to_db(dipole_pattern_numeric(th2[::10], 0.5) + 1e-12), "o", ms=4, label="sum of 200 elements")
    ax.set_ylim(-30, 1); ax.set_xlabel("theta (deg)"); ax.set_ylabel("relative intensity (dB)"); ax.set_title("half-wave dipole"); ax.legend()
    ax = fig.add_subplot(1, 3, 3)
    el = np.linspace(2, 90, 200)
    ax.plot(el, [Link(elevation_deg=e).budget()["margin_dB"] for e in el])
    ax.axhline(3.0, color="r", lw=0.8)
    ax.set_xlabel("elevation (deg)"); ax.set_ylabel("Eb/N0 margin (dB)"); ax.set_title("437 MHz CubeSat downlink, 500 km")
    fig.tight_layout()
    fig.savefig(out_dir / "dipole_pattern.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
