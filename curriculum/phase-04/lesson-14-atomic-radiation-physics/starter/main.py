"""P04.14 Atomic & Radiation Physics: starter code.

Complete every function marked NotImplementedError. When finished, running this file produces the lesson's output and figure.

From the hydrogen atom to the radiation budget of a spacecraft:

* Bohr energy levels with the reduced-mass correction, transition wavelengths (Lyman, Balmer), photon
  energies and photon counts of a lidar pulse, lidar range from round-trip time;
* nuclear binding energy from measured atomic masses, the semi-empirical mass formula as a model, and
  the Q-value of the Pu-238 alpha decay;
* the half-life law, activity, specific activity and decay heat of Pu-238;
* a two-member decay chain (Bateman solution) checked against ``scipy.integrate.solve_ivp``;
* absorbed dose (gray), equivalent dose (sievert), dose from a particle fluence and stopping power,
  total ionizing dose with a design margin;
* single-event upset rate and the Poisson probability of at least one upset;
* RTG electrical power over a mission (Pu-238 decay plus an assumed thermocouple degradation);
* a first-order exponential-attenuation shielding MODEL: dose behind an areal density, the attenuation
  length fitted from two dose-depth points, and the thickness for a dose target.

Nothing here is measured: every input is a stated constant, a published atomic mass, or an explicitly
labelled illustrative number. The shielding law is a model, not a transport calculation.

Run:  python main.py          (prints the tables, writes radiation_dose.png)
Test: pytest curriculum/phase-04/lesson-14-atomic-radiation-physics
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------- constants (SI unless stated)

H_PLANCK = 6.62607015e-34       # J s
C = 299_792_458.0               # m/s
E_CHARGE = 1.602176634e-19      # C (J per eV)
N_A = 6.02214076e23             # 1/mol
RYDBERG_EV = 13.605693122994    # eV, Rydberg energy (infinite nuclear mass)
M_E_U = 5.48579909065e-4        # electron mass in u
M_P_U = 1.007276466621          # proton mass in u
M_H_U = 1.00782503207           # hydrogen-1 atomic mass in u
M_N_U = 1.00866491595           # neutron mass in u
U_MEV = 931.49410242            # MeV per u
MEV_J = 1.602176634e-13         # J per MeV
YEAR_S = 365.25 * 86400.0       # s (Julian year)

# published atomic masses (u), AME 2020 values rounded to 6 decimals
ATOMIC_MASS_U = {
    "H-2": 2.014102, "He-4": 4.002603, "C-12": 12.000000, "O-16": 15.994915,
    "Fe-56": 55.934936, "Ni-62": 61.928345, "U-234": 234.040952, "U-235": 235.043930,
    "U-238": 238.050788, "Pu-238": 238.049560,
}
PU238_HALF_LIFE_Y = 87.7         # years
PU238_MOLAR_MASS = 238.049560    # g/mol

RADIATION_WEIGHTING = {"photon": 1.0, "electron": 1.0, "proton": 2.0, "alpha": 20.0}


# ---------------------------------------------------------------- the hydrogen atom and photons

def energy_level(n: int | np.ndarray, Z: int = 1, reduced_mass: bool = True) -> np.ndarray | float:
    """Bohr energy of level n [eV]: -Z^2 R_y / n^2, with R_y scaled by mu/m_e for a hydrogen nucleus
    when ``reduced_mass`` is True (13.5983 eV instead of 13.6057 eV for the ground state)."""
    raise NotImplementedError


def photon_energy_ev(wavelength_m: np.ndarray | float) -> np.ndarray | float:
    """Photon energy h c / lambda in eV."""
    raise NotImplementedError


def transition_wavelength(n_upper: int, n_lower: int, Z: int = 1, reduced_mass: bool = True) -> float:
    """Vacuum wavelength [m] of the photon emitted in the n_upper -> n_lower transition."""
    raise NotImplementedError


def photons_per_pulse(pulse_energy_j: float, wavelength_m: float) -> float:
    """Number of photons in a laser pulse of the given energy and wavelength."""
    raise NotImplementedError


def lidar_range(round_trip_s: np.ndarray | float) -> np.ndarray | float:
    """Range c t / 2 [m] from the round-trip time of a lidar pulse."""
    raise NotImplementedError


# ---------------------------------------------------------------- nuclei

def binding_energy_mev(Z: int, A: int, atomic_mass_u: float) -> float:
    """Nuclear binding energy [MeV] from the atomic mass: (Z m_H + N m_n - M) c^2.

    Using the hydrogen-atom mass (not the proton mass) cancels the Z electron masses."""
    raise NotImplementedError


def semf_binding_energy_mev(Z: int, A: int) -> float:
    """Semi-empirical (liquid-drop) mass formula [MeV], a MODEL:
    a_V A - a_S A^(2/3) - a_C Z(Z-1)/A^(1/3) - a_A (A-2Z)^2/A + pairing,
    with a_V = 15.75, a_S = 17.8, a_C = 0.711, a_A = 23.7, a_P = 11.18 MeV; pairing +a_P/sqrt(A) for
    even-even, -a_P/sqrt(A) for odd-odd, 0 for odd A. Good to about 1 % for A > 30, poor for light nuclei."""
    raise NotImplementedError


def q_value_mev(parent_u: float, *products_u: float) -> float:
    """Q-value [MeV] of a decay from atomic masses: (M_parent - sum M_products) c^2. Positive = allowed."""
    raise NotImplementedError


# ---------------------------------------------------------------- decay

def decay_constant(half_life: float) -> float:
    """lambda = ln 2 / T_half (same time unit as the half-life, inverted)."""
    raise NotImplementedError


def remaining_fraction(t: np.ndarray | float, half_life: float) -> np.ndarray | float:
    """N(t)/N0 = exp(-lambda t) = 2^(-t/T_half)."""
    raise NotImplementedError


def activity_bq(n_atoms: float, half_life_s: float) -> float:
    """Activity A = lambda N [Bq] for N atoms and a half-life in seconds."""
    raise NotImplementedError


def specific_activity_bq_per_g(half_life_s: float, molar_mass_g: float) -> float:
    """Activity of one gram of a pure nuclide [Bq/g]: lambda N_A / M."""
    raise NotImplementedError


def decay_heat_w_per_g(half_life_s: float, molar_mass_g: float, q_mev: float) -> float:
    """Thermal power of one gram [W/g], assuming all the decay energy Q is deposited locally
    (true for alpha decay in a solid fuel pellet; the recoil is included in Q)."""
    raise NotImplementedError


def bateman_two(t: np.ndarray | float, N1_0: float, lam1: float, lam2: float) -> tuple[np.ndarray, np.ndarray]:
    """Parent and daughter populations of a chain 1 -> 2 -> (stable), daughter initially absent:
    N1 = N1_0 e^{-l1 t},  N2 = N1_0 l1/(l2 - l1) (e^{-l1 t} - e^{-l2 t}). Requires l1 != l2."""
    raise NotImplementedError


def chain_numeric(t: np.ndarray, N1_0: float, lam1: float, lam2: float) -> tuple[np.ndarray, np.ndarray]:
    """The same chain integrated with solve_ivp (Radau, tight tolerances): the oracle for ``bateman_two``."""
    raise NotImplementedError


# ---------------------------------------------------------------- dose

def absorbed_dose_gy(energy_j: float, mass_kg: float) -> float:
    """Absorbed dose [Gy = J/kg]."""
    raise NotImplementedError


def equivalent_dose_sv(doses_gy: dict[str, float]) -> float:
    """Equivalent dose [Sv] = sum over radiation types of w_R D_R, with w_R from RADIATION_WEIGHTING
    (photons and electrons 1, protons 2, alphas 20)."""
    raise NotImplementedError


def dose_from_fluence_gy(fluence_per_cm2: float, stopping_power_mev_cm2_per_g: float) -> float:
    """Dose [Gy] deposited by a fluence of particles crossing a thin layer: Phi S, with
    1 MeV/g = 1.602e-10 Gy. Thin layer: the particle loses little energy in it."""
    raise NotImplementedError


def gy_to_krad(d_gy: float) -> float:
    """1 rad = 0.01 Gy, so 1 krad = 10 Gy."""
    raise NotImplementedError


def total_ionizing_dose(dose_rate_per_year: float, years: float, margin: float = 2.0) -> float:
    """Design TID = dose rate x mission duration x radiation design margin (same dose unit as the rate)."""
    raise NotImplementedError


def seu_rate_per_day(cross_section_cm2_per_bit: float, flux_per_cm2_s: float, n_bits: float) -> float:
    """Expected single-event upsets per day: sigma x flux x bits x 86 400 s (a constant-flux, constant-sigma model)."""
    raise NotImplementedError


def prob_at_least_one(rate: float, t: float) -> float:
    """Poisson probability of at least one event in time t at mean rate ``rate``: 1 - exp(-rate t)."""
    raise NotImplementedError


# ---------------------------------------------------------------- RTG

def rtg_electrical_power(t_years: np.ndarray | float, p_thermal_0: float, efficiency: float,
                         half_life_y: float = PU238_HALF_LIFE_Y, degradation_per_year: float = 0.0) -> np.ndarray | float:
    """Electrical power [W] after t years: P_th0 2^(-t/T) x efficiency x (1 - d)^t.

    The 2^(-t/T) factor is the Pu-238 decay; ``degradation_per_year`` is an ASSUMED fractional loss of the
    thermocouples' conversion efficiency per year (0 by default)."""
    raise NotImplementedError


def rtg_fuel_mass_g(p_electric_eol: float, mission_years: float, efficiency: float, q_mev: float,
                    degradation_per_year: float = 0.0) -> float:
    """Mass of pure Pu-238 [g] needed so that the electrical power at end of mission is ``p_electric_eol``."""
    raise NotImplementedError


# ---------------------------------------------------------------- shielding (exponential MODEL)

def areal_density_g_cm2(thickness_mm: np.ndarray | float, density_g_cm3: float = 2.70) -> np.ndarray | float:
    """Areal density rho x [g/cm^2] of a slab (aluminium by default; 1 mm Al = 0.27 g/cm^2)."""
    raise NotImplementedError


def shielded_dose(d0: float, areal_density: np.ndarray | float, attenuation_length: float) -> np.ndarray | float:
    """MODEL: dose behind areal density x is D0 exp(-x / L). D0 is the dose at x = 0 of the model
    (not the unshielded dose, which a real spectrum makes much larger), L the effective attenuation length
    in g/cm^2. Exact for narrow-beam monoenergetic photons (L = 1/(mu/rho)); a fit, valid only over the
    fitted range, for the electron and proton spectra of space."""
    raise NotImplementedError


def fit_attenuation(x1: float, d1: float, x2: float, d2: float) -> tuple[float, float]:
    """Fit D0 and L of the exponential model through two (areal density, dose) points."""
    raise NotImplementedError


def areal_density_for_dose(d_target: float, d0: float, attenuation_length: float) -> float:
    """Areal density [g/cm^2] at which the model dose falls to d_target (0 if no shielding is needed)."""
    raise NotImplementedError


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson tables and write radiation_dose.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    print(f"hydrogen ground state: {energy_level(1):.4f} eV (reduced mass), {energy_level(1, reduced_mass=False):.4f} eV (infinite mass)")
    for (u, l, name) in ((2, 1, "Lyman alpha"), (3, 2, "H-alpha"), (4, 2, "H-beta")):
        print(f"  {name:12s} {transition_wavelength(u, l) * 1e9:8.2f} nm")
    print(f"Nd:YAG lidar 1064 nm: photon {photon_energy_ev(1064e-9):.3f} eV, 10 mJ pulse = {photons_per_pulse(10e-3, 1064e-9):.3e} photons; "
          f"6.67 us round trip = {lidar_range(6.67e-6):.0f} m")
    for name, Z, A in (("He-4", 2, 4), ("Fe-56", 26, 56), ("U-238", 92, 238), ("Pu-238", 94, 238)):
        B = binding_energy_mev(Z, A, ATOMIC_MASS_U[name])
        print(f"  {name:7s} B = {B:8.2f} MeV, B/A = {B / A:.3f} MeV; SEMF B/A = {semf_binding_energy_mev(Z, A) / A:.3f}")
    q = q_value_mev(ATOMIC_MASS_U["Pu-238"], ATOMIC_MASS_U["U-234"], ATOMIC_MASS_U["He-4"])
    T = PU238_HALF_LIFE_Y * YEAR_S
    sa = specific_activity_bq_per_g(T, PU238_MOLAR_MASS)
    print(f"Pu-238 alpha decay Q = {q:.3f} MeV; specific activity {sa:.3e} Bq/g = {sa / 3.7e10:.2f} Ci/g; "
          f"decay heat {decay_heat_w_per_g(T, PU238_MOLAR_MASS, q):.4f} W/g")
    for yrs in (0, 14, 17, 50):
        print(f"  after {yrs:2d} y: {remaining_fraction(yrs, PU238_HALF_LIFE_Y) * 100:6.2f} % of the heat")
    m = rtg_fuel_mass_g(110.0, 14.0, 0.06, q)
    print(f"Pu-238 for 110 We at end of 14 y, 6 % efficiency: {m / 1e3:.2f} kg; with 1.5 %/y degradation: "
          f"{rtg_fuel_mass_g(110.0, 14.0, 0.06, q, 0.015) / 1e3:.2f} kg; thermal power at launch "
          f"{m * decay_heat_w_per_g(T, PU238_MOLAR_MASS, q):.0f} W, electrical {rtg_electrical_power(0.0, m * decay_heat_w_per_g(T, PU238_MOLAR_MASS, q), 0.06):.1f} W")
    d0, L = fit_attenuation(areal_density_g_cm2(1.0), 20.0, areal_density_g_cm2(3.0), 4.0)
    x_need = areal_density_for_dose(30.0 / (10 * 2.0), d0, L)
    print(f"shielding MODEL (illustrative dose-depth points 20 and 4 krad/yr behind 1 and 3 mm Al): "
          f"D0 = {d0:.1f} krad/yr, L = {L:.4f} g/cm^2; 30 krad part, 10 y, margin 2 -> "
          f"{x_need:.3f} g/cm^2 = {x_need / 0.27:.2f} mm Al")
    r = seu_rate_per_day(1e-14, 1.0, 8 * 2 ** 30)
    print(f"SEU (illustrative): sigma 1e-14 cm^2/bit, 1 /cm^2/s, 1 GiB -> {r:.2f} upsets/day; "
          f"P(>=1 in 1 h) = {prob_at_least_one(r / 24, 1.0):.3f}")

    t = np.linspace(0, 30, 301)
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    p_th0 = m * decay_heat_w_per_g(T, PU238_MOLAR_MASS, q)
    ax[0].plot(t, rtg_electrical_power(t, p_th0, 0.06), label="decay only")
    ax[0].plot(t, rtg_electrical_power(t, p_th0, 0.06, degradation_per_year=0.015), label="+1.5 %/y (assumed)")
    ax[0].axhline(110.0, color="gray", ls="--", lw=0.8)
    ax[0].set_xlabel("mission time (years)"); ax[0].set_ylabel("electrical power (W)"); ax[0].legend()
    ax[0].set_title(f"RTG, {p_th0:.0f} W thermal at launch")
    x = np.linspace(0, 2.0, 200)
    ax[1].semilogy(x, shielded_dose(d0, x, L))
    ax[1].semilogy([areal_density_g_cm2(1.0), areal_density_g_cm2(3.0)], [20.0, 4.0], "o")
    ax[1].set_xlabel("aluminium areal density (g/cm^2)"); ax[1].set_ylabel("dose (krad/yr), MODEL")
    ax[1].set_title("exponential shielding model (illustrative)")
    tt = np.linspace(0, 5, 200)
    lam1, lam2 = decay_constant(2.0), decay_constant(0.3)
    N1, N2 = bateman_two(tt, 1.0, lam1, lam2)
    ax[2].plot(tt, lam1 * N1, label="parent activity")
    ax[2].plot(tt, lam2 * N2, label="daughter activity")
    ax[2].set_xlabel("time (units of the half-lives)"); ax[2].set_ylabel("activity (arb.)"); ax[2].legend()
    ax[2].set_title("two-member chain, T1 = 2, T2 = 0.3")
    fig.tight_layout()
    fig.savefig(out_dir / "radiation_dose.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
