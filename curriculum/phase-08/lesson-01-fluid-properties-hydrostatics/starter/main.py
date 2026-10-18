"""P08.01 Fluid Properties & Hydrostatics: starter code.

What are the temperature, pressure, density and viscosity of the air the light aeroplane flies in
at 3000 m, and of the air a sounding rocket climbs through up to 47 km? This module builds the
International Standard Atmosphere (ISA) from one law of physics and one assumption:

* hydrostatic equilibrium, dp/dh = -rho g0 (a fluid at rest carries no shear stress, so the only
  force balancing the weight of a slab of air is the pressure difference across it);
* the ideal-gas law p = rho R T (`P07.01`), with the ISA's prescribed temperature profile T(h):
  piecewise linear in geopotential altitude h.

It computes the pressure two independent ways: the closed forms obtained by integrating the
hydrostatic equation by hand layer by layer (power law in a gradient layer, exponential in an
isothermal one), and a fixed-step RK4 integration of dp/dh = -g0 p / (R T(h)) that knows nothing
about the closed forms. Both are checked against the published ISA anchor values.

Also: Sutherland's law for the viscosity of air, the U-tube and inclined manometers, the pressure
at the bottom of an accelerating propellant tank, and Archimedes' principle checked by integrating
the pressure over the surface of a submerged sphere (liquid) and of a balloon (compressible ISA).

Units: SI throughout (m, s, kg, K, Pa). Altitudes are geopotential unless a name says geometric.

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes isa_atmosphere.png)
Test: pytest curriculum/phase-08/lesson-01-fluid-properties-hydrostatics
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable

import numpy as np

# ---------------------------------------------------------------- constants (ISA, published)

G0 = 9.80665                 # m/s^2, standard gravity
R_AIR = 287.05287            # J/(kg K), specific gas constant of dry air used by the ISA
GAMMA = 1.4                  # ratio of specific heats of air
T0 = 288.15                  # K, sea-level temperature
P0 = 101325.0                # Pa, sea-level pressure
R_EARTH = 6356766.0          # m, the Earth radius used to define geopotential altitude
SUTH_C1 = 1.458e-6           # kg/(m s K^0.5), Sutherland constant for air
SUTH_S = 110.4               # K, Sutherland temperature for air

# Layer base geopotential altitudes (m) and lapse rates dT/dh (K/m), sea level to 47 km.
LAYER_BASES = np.array([0.0, 11000.0, 20000.0, 32000.0, 47000.0])
LAYER_LAPSE = np.array([-6.5e-3, 0.0, 1.0e-3, 2.8e-3])
H_TOP = 47000.0

# Water and air at 20 C (lesson-wide shared numbers).
RHO_WATER = 998.2            # kg/m^3
RHO_AIR_20C = 1.204          # kg/m^3


def _ret(a: np.ndarray):
    """Return a Python-friendly scalar for 0-d input, the array otherwise."""
    a = np.asarray(a)
    return a[()] if a.ndim == 0 else a


def _layer_index(h: np.ndarray) -> np.ndarray:
    """Index of the ISA layer containing each geopotential altitude h (0 <= h <= 47 km)."""
    h = np.asarray(h, dtype=float)
    if np.any(h < -1e-9) or np.any(h > H_TOP + 1e-9):
        raise ValueError("this model covers geopotential altitudes 0 to 47 000 m")
    return np.clip(np.searchsorted(LAYER_BASES, h, side="right") - 1, 0, len(LAYER_LAPSE) - 1)


def _base_temperatures() -> np.ndarray:
    """Temperature at the base of each layer, built by stepping the lapse rates up from T0."""
    Tb = [T0]
    for k in range(len(LAYER_LAPSE)):
        Tb.append(Tb[-1] + LAYER_LAPSE[k] * (LAYER_BASES[k + 1] - LAYER_BASES[k]))
    return np.array(Tb)


def layer_base_pressures() -> np.ndarray:
    """Pressure (Pa) at the base of each layer (0, 11, 20, 32, 47 km) from the closed forms of
    ``isa_pressure``, stepping up from P0 one layer at a time."""
    raise NotImplementedError


# ---------------------------------------------------------------- the standard atmosphere

def geopotential_altitude(z: np.ndarray) -> np.ndarray:
    """Geopotential altitude h = r0 z / (r0 + z) from geometric altitude z (m), r0 = 6 356 766 m.
    h is the height that has the same potential energy per unit mass at constant g0 as z has in the
    real inverse-square field, so dp/dh = -rho g0 holds exactly in h."""
    z = np.asarray(z, dtype=float)
    return _ret(R_EARTH * z / (R_EARTH + z))


def isa_temperature(h: np.ndarray) -> np.ndarray:
    """ISA temperature (K) at geopotential altitude h (m), 0-47 km: T = T_b + a (h - h_b) in the
    layer whose base is h_b, with lapse rates -6.5, 0, +1.0, +2.8 K/km."""
    raise NotImplementedError


def isa_pressure(h: np.ndarray) -> np.ndarray:
    """ISA pressure (Pa) at geopotential altitude h (m) from the closed-form integrals of
    dp/dh = -g0 p / (R T):

        gradient layer (a != 0):  p = p_b (T / T_b)^(-g0 / (a R))
        isothermal layer (a = 0): p = p_b exp(-g0 (h - h_b) / (R T_b))
    """
    raise NotImplementedError


def isa_density(h: np.ndarray) -> np.ndarray:
    """ISA density (kg/m^3) from the ideal-gas law rho = p / (R T)."""
    raise NotImplementedError


def integrate_hydrostatic(h_max: float, dh: float, p_start: float = P0,
                          T_of_h: Callable[[np.ndarray], np.ndarray] | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Integrate dp/dh = -g0 p / (R T(h)) from h = 0 (p = p_start) to h_max with classical RK4 at a
    fixed step dh. ``T_of_h`` defaults to ``isa_temperature``. This is the independent route: it
    uses the ODE only, not the closed forms. Choose dh so that the layer boundaries fall on steps
    (T(h) has kinks there; RK4 keeps fourth order only if no step straddles one).
    Returns (h, p) on the n + 1 nodes."""
    raise NotImplementedError


def incompressible_pressure(h: np.ndarray, rho: float = 1.225, p_ref: float = P0) -> np.ndarray:
    """Hydrostatic pressure of a constant-density fluid: p = p_ref - rho g0 h (h up). Right for a
    liquid; for air it is the first term of the Taylor series of the ISA and fails with height."""
    raise NotImplementedError


def scale_height(T: float) -> float:
    """Pressure scale height H = R T / g0 (m) of an isothermal atmosphere: p falls by e every H."""
    raise NotImplementedError


# ---------------------------------------------------------------- viscosity

def sutherland_viscosity(T: np.ndarray) -> np.ndarray:
    """Dynamic viscosity of air (Pa s) by Sutherland's law, mu = C1 T^1.5 / (T + S), with
    C1 = 1.458e-6 kg/(m s K^0.5) and S = 110.4 K. Gas viscosity rises with temperature (molecules
    carry momentum between layers faster), unlike a liquid's."""
    raise NotImplementedError


def kinematic_viscosity(h: np.ndarray) -> np.ndarray:
    """Kinematic viscosity nu = mu / rho (m^2/s) of ISA air at geopotential altitude h."""
    raise NotImplementedError


# ---------------------------------------------------------------- manometers, tanks, buoyancy

def manometer_delta_p(dh: np.ndarray, rho_m: float = RHO_WATER, rho_f: float = RHO_AIR_20C,
                      angle_deg: float = 90.0, g: float = G0) -> np.ndarray:
    """Pressure difference (Pa) read by a U-tube manometer whose liquid (density rho_m) stands
    with a reading dh between the legs, measured along the tube; the legs are filled above the
    liquid by a fluid of density rho_f. For a tube inclined at angle_deg to the horizontal the
    vertical height is dh sin(angle): dp = (rho_m - rho_f) g dh sin(angle)."""
    raise NotImplementedError


def tank_bottom_pressure(p_ullage: float, rho: float, fill_height: float, accel: float = 0.0,
                         g: float = G0) -> float:
    """Pressure (Pa) at the outlet of a propellant tank whose liquid column of height fill_height
    sits under an ullage pressure p_ullage while the vehicle accelerates along its axis at accel
    (m/s^2, positive up). In the tank's frame the liquid feels an effective gravity g + accel, so
    the column adds rho (g + accel) fill_height. In free fall (accel = -g) it adds nothing."""
    raise NotImplementedError


def pressure_force_on_sphere(p_of_z: Callable[[np.ndarray], np.ndarray], z_c: float, radius: float,
                             n: int = 200) -> float:
    """Net vertical pressure force (N, positive up) on a sphere of the given radius centred at
    height z_c in a fluid at rest whose pressure is p_of_z(z), by integrating -p n dA over the
    surface with the midpoint rule on n polar bands:

        F_z = -sum p(z_c + R cos th) cos th * 2 pi R^2 sin th d th.

    This uses only the surface pressures; Archimedes' principle (F = weight of displaced fluid)
    is what it should reproduce."""
    raise NotImplementedError


def displaced_weight(rho_of_z: Callable[[np.ndarray], np.ndarray], z_c: float, radius: float,
                     n: int = 2000, g: float = G0) -> float:
    """Weight (N) of the fluid a sphere displaces, integral of rho(z) g over its volume, by the
    midpoint rule on n horizontal discs of area pi (R^2 - (z - z_c)^2)."""
    raise NotImplementedError


def observed_order(errors: np.ndarray, ratio: float = 2.0) -> np.ndarray:
    """Observed convergence orders log(e_k / e_{k+1}) / log(ratio) for grids refined by ``ratio``."""
    e = np.asarray(errors, dtype=float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


# ---------------------------------------------------------------- the lesson's run

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save isa_atmosphere.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    published = {11000.0: 22632.06, 20000.0: 5474.889, 32000.0: 868.019, 47000.0: 110.906}
    print("ISA anchors (published data) against the closed forms and RK4 (dh = 100 m):")
    h_rk, p_rk = integrate_hydrostatic(H_TOP, 100.0)
    print("   h (m)    T (K)    published p (Pa)  closed form (Pa)   RK4 (Pa)       rel. diff closed/published")
    for hh, pp in published.items():
        pc = float(isa_pressure(hh))
        pr = float(p_rk[int(round(hh / 100.0))])
        print(f"  {hh:6.0f}  {float(isa_temperature(hh)):7.2f}  {pp:14.3f}  {pc:14.3f}  {pr:14.3f}   {pc / pp - 1:+.1e}")
    print(f"rho(0) = {float(isa_density(0.0)):.5f} kg/m^3, rho(11 km) = {float(isa_density(11000.0)):.5f} kg/m^3 "
          f"(published 1.225, 0.36392)")
    print(f"mu(T0) = {float(sutherland_viscosity(T0)):.5e} Pa s (published 1.7894e-5)")

    h_cr = 3000.0
    T, p, rho = float(isa_temperature(h_cr)), float(isa_pressure(h_cr)), float(isa_density(h_cr))
    mu = float(sutherland_viscosity(T))
    print(f"\nlight aeroplane at {h_cr:.0f} m: T = {T:.2f} K, p = {p:.1f} Pa, rho = {rho:.4f} kg/m^3, "
          f"sigma = rho/rho0 = {rho / float(isa_density(0.0)):.4f}, mu = {mu:.4e} Pa s, nu = {mu / rho:.4e} m^2/s")
    print(f"  incompressible estimate p0 - rho0 g h = {float(incompressible_pressure(h_cr)):.1f} Pa "
          f"({float(incompressible_pressure(h_cr)) / p - 1:+.2%})")

    print("\nRK4 step study, error at 47 km against the closed form:")
    errs = []
    steps = [1000.0, 500.0, 250.0, 125.0]
    for dh in steps:
        _, pp = integrate_hydrostatic(H_TOP, dh)
        errs.append(abs(pp[-1] - float(isa_pressure(H_TOP))) / float(isa_pressure(H_TOP)))
        print(f"  dh = {dh:6.1f} m: relative error {errs[-1]:.3e}")
    print(f"  observed orders {np.array2string(observed_order(errs), precision=2)}")

    print(f"\nscale height at T0: {scale_height(T0):.0f} m; at 216.65 K: {scale_height(216.65):.0f} m")

    # manometer: 25 mm of water
    print(f"\nU-tube, 25.0 mm of water under air: {float(manometer_delta_p(0.025)):.2f} Pa; inclined at 10 deg, "
          f"the same pressure gives {0.025 / np.sin(np.radians(10.0)) * 1e3:.1f} mm along the tube")

    # propellant tank (LOX, typical density 1141 kg/m^3) under 3 g of thrust acceleration
    p_t = tank_bottom_pressure(3.0e5, 1141.0, 2.0, accel=2 * G0)
    print(f"LOX tank (1141 kg/m^3 typical), 2 m column, 3.0 bar ullage, 3 g total: outlet {p_t / 1e5:.3f} bar "
          f"(on the pad {tank_bottom_pressure(3.0e5, 1141.0, 2.0) / 1e5:.3f} bar)")

    # buoyancy: a 2 m radius helium balloon at 3000 m
    Rb = 2.0
    F = pressure_force_on_sphere(lambda z: isa_pressure(z), h_cr, Rb, n=400)
    W = displaced_weight(lambda z: isa_density(z), h_cr, Rb)
    print(f"balloon r = {Rb} m at {h_cr:.0f} m: surface-pressure force {F:.3f} N, displaced weight {W:.3f} N, "
          f"rho g V = {rho * G0 * 4 / 3 * np.pi * Rb ** 3:.3f} N")

    h = np.linspace(0.0, H_TOP, 471)
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.5), sharey=True)
    ax[0].plot(isa_temperature(h), h / 1e3)
    ax[0].set_xlabel("T (K)")
    ax[0].set_ylabel("geopotential altitude h (km)")
    ax[0].set_title("ISA temperature")
    ax[1].semilogx(isa_pressure(h), h / 1e3, label="closed form")
    ax[1].semilogx(p_rk[::10], h_rk[::10] / 1e3, ".", ms=3, label="RK4, dh = 100 m")
    ax[1].semilogx(list(published.values()), np.array(list(published.keys())) / 1e3, "ks", mfc="none",
                   label="ISA anchors (published)")
    ax[1].set_xlabel("p (Pa)")
    ax[1].set_title("pressure: two routes")
    ax[1].legend(fontsize=8)
    ax[2].plot(sutherland_viscosity(isa_temperature(h)) * 1e5, h / 1e3, label="mu (1e-5 Pa s)")
    ax[2].semilogx(kinematic_viscosity(h) * 1e5, h / 1e3, label="nu (1e-5 m^2/s)")
    ax[2].set_xscale("log")
    ax[2].set_title("viscosity (Sutherland)")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "isa_atmosphere.png", dpi=110)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'isa_atmosphere.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
