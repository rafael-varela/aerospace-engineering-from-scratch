"""P04.07 Kinetic Theory of Gases: starter code.

Complete every function marked NotImplementedError. When finished, running this file produces the
lesson's figure and prints the verification tables.

Pressure and temperature from molecular motion, applied to the air an aircraft flies through:

* pressure from momentum transfer to a wall, p = n m <v^2> / 3, and the ideal-gas law p = n k T;
* the Maxwell-Boltzmann speed distribution: its density, the most probable, mean and rms speeds, the
  fraction above a given speed (atmospheric escape of hydrogen against nitrogen);
* equipartition: molar heat capacities and gamma = (f + 2) / f from the number of degrees of freedom;
  the speed of sound a = sqrt(gamma R T / M) and its ratio to the rms speed, sqrt(gamma / 3);
* mean free path and the Knudsen number: where the continuum model of the air stops;
* two SIMULATIONS (labelled as simulations wherever they appear): free molecules in a box, whose wall
  momentum transfer reproduces n k T; and a random-pair elastic-collision model that relaxes any speed
  distribution to Maxwell-Boltzmann while conserving energy and momentum exactly;
* the International Standard Atmosphere to 86 km (seven layers) as the application: T, p, rho and the
  speed of sound against altitude, checked against the hydrostatic equation.

Units: SI (K, Pa, kg/m^3, m/s, J). Molar masses in kg/mol.

Run:  python main.py          (prints the verification tables, writes kinetic_theory.png)
Test: pytest curriculum/phase-04/lesson-07-kinetic-theory-of-gases
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.integrate import quad

K_B = 1.380649e-23                   # J/K
N_A = 6.02214076e23                  # 1/mol
R_GAS = K_B * N_A                    # J/(mol K)
M_AIR = 28.9647e-3                   # kg/mol
M_N2 = 28.0134e-3
M_O2 = 31.998e-3
M_H2 = 2.01588e-3
M_HE = 4.0026e-3
G0 = 9.80665                         # m/s^2
R_AIR = R_GAS / M_AIR                # J/(kg K) = 287.05
GAMMA_AIR = 1.4
D_N2 = 3.7e-10                       # m, effective hard-sphere diameter of N2 for the mean free path

# ISA layers: base geopotential altitude (m), lapse rate (K/m); base temperature and pressure computed from them
ISA_BASES = np.array([0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0, 86000.0])
ISA_LAPSE = np.array([-0.0065, 0.0, 0.001, 0.0028, 0.0, -0.0028, -0.002])
ISA_T0, ISA_P0 = 288.15, 101325.0


# ---------------------------------------------------------------- pressure from molecular motion

def pressure_from_momentum_transfer(n_density: float, m: float, v_rms: float) -> float:
    """p = n m v_rms^2 / 3: each wall collision transfers 2 m v_x, the flux of molecules at the wall is
    n <|v_x|> / 2, and <v_x^2> = <v^2> / 3 by isotropy. n_density in 1/m^3, m in kg, v_rms in m/s."""
    raise NotImplementedError


def ideal_gas_pressure(n_density: float, T: float) -> float:
    """p = n k T."""
    raise NotImplementedError


def number_density(p: float, T: float) -> float:
    """n = p / (k T) [1/m^3]."""
    raise NotImplementedError


def rms_speed(T: float, M: float) -> float:
    """v_rms = sqrt(3 R T / M) for molar mass M [kg/mol]."""
    raise NotImplementedError


def mean_speed(T: float, M: float) -> float:
    """<v> = sqrt(8 R T / (pi M))."""
    raise NotImplementedError


def most_probable_speed(T: float, M: float) -> float:
    """v_p = sqrt(2 R T / M): the peak of the Maxwell-Boltzmann density."""
    raise NotImplementedError


def maxwell_boltzmann_pdf(v: np.ndarray | float, T: float, M: float) -> np.ndarray | float:
    """Density of molecular speeds, f(v) = 4 pi (m / 2 pi k T)^{3/2} v^2 exp(-m v^2 / 2 k T), with m = M / N_A."""
    raise NotImplementedError


def fraction_above(v_min: float, T: float, M: float) -> float:
    """Fraction of molecules with speed above v_min, by numerical integration of the density."""
    raise NotImplementedError


def sample_speeds(n: int, T: float, M: float, rng: np.random.Generator) -> np.ndarray:
    """Speeds of n molecules drawn from the Maxwell-Boltzmann distribution: each velocity component is
    normal with variance k T / m, and the speed is the norm."""
    raise NotImplementedError


# ---------------------------------------------------------------- equipartition, gamma, speed of sound

def gamma_from_dof(f: int | float) -> float:
    """gamma = Cp / Cv = (f + 2) / f for f quadratic degrees of freedom per molecule (3 monatomic, 5 diatomic
    at room temperature, 7 with vibration excited)."""
    raise NotImplementedError


def molar_heat_capacities(f: int | float) -> tuple[float, float]:
    """(Cv, Cp) per mole from equipartition: Cv = f R / 2, Cp = Cv + R."""
    raise NotImplementedError


def speed_of_sound(T: float, gamma: float = GAMMA_AIR, M: float = M_AIR) -> float:
    """a = sqrt(gamma R T / M)."""
    raise NotImplementedError


def sound_to_rms_ratio(gamma: float) -> float:
    """a / v_rms = sqrt(gamma / 3): sound travels at a fixed fraction of the molecular rms speed."""
    raise NotImplementedError


# ---------------------------------------------------------------- mean free path and Knudsen number

def mean_free_path(T: float, p: float, d: float = D_N2) -> float:
    """lambda = k T / (sqrt(2) pi d^2 p) for hard spheres of diameter d."""
    raise NotImplementedError


def collision_frequency(T: float, p: float, M: float = M_AIR, d: float = D_N2) -> float:
    """Collisions per second per molecule, <v> / lambda."""
    raise NotImplementedError


def knudsen_number(T: float, p: float, length: float, d: float = D_N2) -> float:
    """Kn = lambda / length. Below about 0.01 the air is a continuum; above about 10 it is free molecules."""
    raise NotImplementedError


# ---------------------------------------------------------------- the International Standard Atmosphere

def isa_temperature(h: np.ndarray | float) -> np.ndarray | float:
    """ISA temperature [K] at geopotential altitude h [m], 0 <= h <= 86 km (seven linear layers)."""
    raise NotImplementedError


def isa_pressure(h: np.ndarray | float) -> np.ndarray | float:
    """ISA pressure [Pa] from the hydrostatic equation dp/dh = -rho g0 and the ideal gas law in each layer:
    p = pb (Tb / T)^{g0 / (L R)} for a lapse-rate layer and pb exp(-g0 (h - hb) / (R Tb)) for an isothermal one."""
    raise NotImplementedError


def isa_density(h: np.ndarray | float) -> np.ndarray | float:
    """rho = p / (R_air T) [kg/m^3]."""
    raise NotImplementedError


def isa_speed_of_sound(h: np.ndarray | float) -> np.ndarray | float:
    """a(h) = sqrt(gamma R_air T(h))."""
    raise NotImplementedError


def hydrostatic_residual(h: np.ndarray, dh: float = 1.0) -> np.ndarray:
    """Relative residual of dp/dh + rho g0 = 0 by central differences, a check that the layers are consistent."""
    raise NotImplementedError


# ---------------------------------------------------------------- simulations (labelled as such)

def simulate_box_pressure(n_particles: int, T: float, M: float, L: float, t_total: float,
                          rng: np.random.Generator) -> tuple[float, float]:
    """SIMULATION: n non-interacting molecules in a cubic box of side L, velocities drawn from the
    Maxwell-Boltzmann distribution, bouncing elastically off the walls. The pressure on the wall x = L is
    the total momentum it receives, 2 m |v_x| per hit, divided by the time and the wall area.

    Hits are counted analytically from free flight: a molecule with |v_x| hits the x = L wall every
    2 L / |v_x| seconds, so over t_total it hits floor((|v_x| t_total + phase) / 2L) times with a random
    starting phase. Returns (p_simulated, p_ideal = n k T) with n = n_particles / L^3.
    """
    raise NotImplementedError


def relax_to_maxwellian(velocities: np.ndarray, n_collisions: int, rng: np.random.Generator) -> np.ndarray:
    """SIMULATION of elastic collisions: repeatedly pick two molecules of equal mass at random and give
    their relative velocity a random new direction, keeping its magnitude and the centre-of-mass velocity.
    Each collision conserves momentum and kinetic energy exactly; whatever the initial velocities, the
    speeds relax to the Maxwell-Boltzmann distribution at the temperature set by the total energy.

    ``velocities`` has shape (n, 3); a modified copy is returned."""
    raise NotImplementedError


def temperature_from_velocities(velocities: np.ndarray, M: float) -> float:
    """Kinetic temperature from equipartition: (3/2) k T = (1/2) m <v^2>."""
    raise NotImplementedError


# ---------------------------------------------------------------- driver

def main(out_dir: str | Path = ".") -> None:
    """Print the verification tables and write kinetic_theory.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    T, p = 288.15, 101325.0
    n = number_density(p, T)
    m = M_AIR / N_A
    print(f"sea-level air: n = {n:.3e} /m^3, m = {m:.3e} kg, v_rms = {rms_speed(T, M_AIR):.1f} m/s, "
          f"<v> = {mean_speed(T, M_AIR):.1f} m/s, v_p = {most_probable_speed(T, M_AIR):.1f} m/s")
    print(f"pressure from momentum transfer n m v_rms^2 / 3 = {pressure_from_momentum_transfer(n, m, rms_speed(T, M_AIR)):.1f} Pa "
          f"= n k T = {ideal_gas_pressure(n, T):.1f} Pa")
    print(f"wall collision rate n <v> / 4 = {n * mean_speed(T, M_AIR) / 4:.3e} per m^2 per s")
    norm, _ = quad(lambda v: maxwell_boltzmann_pdf(v, T, M_AIR), 0.0, np.inf)
    print(f"Maxwell-Boltzmann density integrates to {norm:.6f}; fraction above 1000 m/s: {fraction_above(1000.0, T, M_AIR):.3e}")
    cv, cp = molar_heat_capacities(5)
    print(f"diatomic (f = 5): Cv = {cv:.2f}, Cp = {cp:.2f} J/(mol K), gamma = {gamma_from_dof(5):.3f}; monatomic gamma = {gamma_from_dof(3):.4f}")
    a = speed_of_sound(T)
    print(f"speed of sound at 288.15 K: {a:.1f} m/s = {sound_to_rms_ratio(GAMMA_AIR):.4f} x v_rms ({sound_to_rms_ratio(GAMMA_AIR) * rms_speed(T, M_AIR):.1f} m/s)")
    print(f"mean free path at sea level {mean_free_path(T, p)*1e9:.1f} nm, collision frequency {collision_frequency(T, p):.2e} /s")

    print("\natmospheric escape at the exobase (T = 1000 K, escape speed at 500 km altitude = 10.8 km/s):")
    for name, M in (("H2", M_H2), ("He", M_HE), ("N2", M_N2), ("O2", M_O2)):
        print(f"  {name:>3}: v_rms = {rms_speed(1000.0, M):>6.0f} m/s, fraction above 10 800 m/s = {fraction_above(10800.0, 1000.0, M):.2e}")

    print("\nISA (geopotential altitude), with the hydrostatic residual |dp/dh + rho g| / (rho g):")
    print(f"{'h (km)':>7} {'T (K)':>8} {'p (Pa)':>10} {'rho (kg/m^3)':>13} {'a (m/s)':>8} {'lambda (m)':>11} {'Kn (1 m)':>10} {'hydro':>9}")
    for h in (0.0, 5e3, 11e3, 20e3, 32e3, 47e3, 71e3, 86e3):
        Th, ph, rh = isa_temperature(h), isa_pressure(h), isa_density(h)
        res = abs(float(hydrostatic_residual(np.array([min(max(h, 1.0), 85999.0)]))[0]))
        print(f"{h/1e3:>7.0f} {Th:>8.2f} {ph:>10.2f} {rh:>13.4e} {isa_speed_of_sound(h):>8.1f} {mean_free_path(Th, ph):>11.2e} "
              f"{knudsen_number(Th, ph, 1.0):>10.2e} {res:>9.1e}")
    hs = np.linspace(0.0, 86e3, 8601)
    kn = knudsen_number(isa_temperature(hs), isa_pressure(hs), 1.0)
    print(f"Kn = 0.01 for a 1 m body at {hs[np.argmax(kn > 0.01)]/1e3:.1f} km (ISA stops at 86 km; the continuum model is already marginal there)")

    rng = np.random.default_rng(407)
    print("\nSIMULATION 1: free molecules in a box (seed 407)")
    for N in (1000, 10000, 100000):
        ps, pi = simulate_box_pressure(N, 300.0, M_N2, 1e-6, 1e-7, rng)
        print(f"  N = {N:>6d}: p_sim / (n k T) = {ps / pi:.4f}")
    print("SIMULATION 2: elastic collisions relax a single-speed gas to Maxwell-Boltzmann (seed 407)")
    N = 20000
    dirs = rng.normal(size=(N, 3))
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    v0 = rms_speed(300.0, M_N2) * dirs                                    # every molecule at the same speed
    v1 = relax_to_maxwellian(v0, 10 * N, rng)
    s1 = np.linalg.norm(v1, axis=1)
    Tk = temperature_from_velocities(v1, M_N2)
    print(f"  energy ratio after/before {np.sum(v1**2)/np.sum(v0**2):.12f}; momentum change {np.linalg.norm(v1.sum(0) - v0.sum(0)):.2e} m/s")
    print(f"  kinetic temperature {Tk:.2f} K; mean speed {s1.mean():.1f} vs {mean_speed(Tk, M_N2):.1f}; "
          f"most-probable bin {np.histogram(s1, 60)[1][np.argmax(np.histogram(s1, 60)[0])]:.0f} vs v_p {most_probable_speed(Tk, M_N2):.1f}")
    print(f"  <v^4>/<v^2>^2 = {np.mean(s1**4)/np.mean(s1**2)**2:.4f} (Maxwell-Boltzmann: 5/3 = 1.6667; single speed: 1)")

    fig, axs = plt.subplots(1, 3, figsize=(14, 4.3))
    ax = axs[0]
    vv = np.linspace(0.0, 2000.0, 400)
    for Tp in (200.0, 300.0, 600.0):
        ax.plot(vv, maxwell_boltzmann_pdf(vv, Tp, M_N2) * 1e3, label=f"N2 at {Tp:.0f} K")
    ax.plot(vv, maxwell_boltzmann_pdf(vv, 300.0, M_H2) * 1e3, "--", label="H2 at 300 K")
    ax.set_xlabel("speed (m/s)")
    ax.set_ylabel("density (per km/s)")
    ax.set_title("Maxwell-Boltzmann")
    ax.legend(fontsize=8)
    ax = axs[1]
    ax.hist(s1, bins=60, density=True, alpha=0.6, label="SIMULATION after 10 N collisions")
    ax.plot(vv, maxwell_boltzmann_pdf(vv, Tk, M_N2), "k-", label=f"Maxwell-Boltzmann at {Tk:.0f} K")
    ax.axvline(rms_speed(300.0, M_N2), color="r", ls=":", label="initial single speed")
    ax.set_xlabel("speed (m/s)")
    ax.set_title("relaxation to equilibrium")
    ax.legend(fontsize=8)
    ax = axs[2]
    ax.plot(isa_temperature(hs), hs / 1e3, label="T (K)")
    ax.plot(isa_speed_of_sound(hs), hs / 1e3, label="a (m/s)")
    ax.set_xlabel("T (K), a (m/s)")
    ax.set_ylabel("geopotential altitude (km)")
    ax.set_title("ISA to 86 km")
    ax2 = ax.twiny()
    ax2.semilogx(isa_pressure(hs), hs / 1e3, "g--", label="p (Pa)")
    ax2.set_xlabel("p (Pa)")
    ax.legend(fontsize=8, loc="center right")
    fig.tight_layout()
    fig.savefig(out_dir / "kinetic_theory.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
