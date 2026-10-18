"""P07.01 Systems, Properties & State: starter code.

Complete every function marked NotImplementedError. When finished, running this file prints the
lesson's tables and writes its figure.

An ideal-gas property calculator built on the NASA Glenn 9-coefficient polynomials:

* cp, h and s at 1 bar for single species from the three NASA polynomial forms, with the interval
  (200-1000 K or 1000-6000 K) picked by temperature;
* dry air as a fixed mixture (N2 0.7808, O2 0.2095, Ar 0.0093, CO2 0.0004 by mole), on a mass basis:
  molar mass, gas constant, cp, h, u, s(T, p) including the partial-pressure (mixing) term, gamma;
* the ideal-gas equation of state p v = R T, checked against the kinetic-theory pressure rho v_rms^2 / 3
  of P04.07;
* the state postulate for a simple compressible substance: any two independent properties of
  (p, T, v, h, u, s) fix the state; (T, h), (T, u) and (h, u) are NOT independent for an ideal gas;
* the comparison that motivates the lesson: constant cp = 1005 J/(kg K) against the variable-cp model
  between compressor delivery and turbine inlet temperatures.

Units: SI. Temperatures in K, pressures in Pa, molar quantities per mol, mixture quantities per kg.
R = 8.31446261815324 J/(mol K) (CODATA 2018); p deg = 1 bar = 1e5 Pa for the NASA data.

Run:  python main.py              (prints the verification tables, writes ideal_gas_properties.png)
Test: pytest curriculum/phase-07/lesson-01-systems-properties-state
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.optimize import brentq

R_U = 8.31446261815324          # J/(mol K), CODATA 2018 (the NASA database itself was fitted with 8.314510)
K_B = 1.380649e-23              # J/K
N_A = 6.02214076e23             # 1/mol
P_REF = 1.0e5                   # Pa, standard-state pressure of the NASA data (1 bar)
P_ATM = 101325.0                # Pa, 1 standard atmosphere
CP_CONST = 1005.0               # J/(kg K), the "textbook constant" cp of air the lesson tests

# NASA Glenn 9-coefficient polynomials, copied unrounded from data/thermo.inp of NASA CEA
# (https://github.com/nasa/cea, Apache-2.0), the database documented in McBride, Zehe & Gordon,
# NASA/TP-2002-211556 (2002). Each interval is (T_low, T_high, [a1..a7], [b1, b2]); M in g/mol; Hf298 in J/mol.
NASA9 = {
    "N2": {"M": 28.0134, "Hf298": 0.0, "intervals": [
        (200.0, 1000.0, [22103.71497, -381.846182, 6.08273836, -0.00853091441, 1.384646189e-05, -9.62579362e-09, 2.519705809e-12], [710.846086, -10.76003744]),
        (1000.0, 6000.0, [587712.406, -2239.249073, 6.06694922, -0.00061396855, 1.491806679e-07, -1.923105485e-11, 1.061954386e-15], [12832.10415, -15.86640027])]},
    "O2": {"M": 31.9988, "Hf298": 0.0, "intervals": [
        (200.0, 1000.0, [-34255.6342, 484.700097, 1.119010961, 0.00429388924, -6.83630052e-07, -2.0233727e-09, 1.039040018e-12], [-3391.45487, 18.4969947]),
        (1000.0, 6000.0, [-1037939.022, 2344.830282, 1.819732036, 0.001267847582, -2.188067988e-07, 2.053719572e-11, -8.19346705e-16], [-16890.10929, 17.38716506])]},
    "Ar": {"M": 39.948, "Hf298": 0.0, "intervals": [
        (200.0, 1000.0, [0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0], [-745.375, 4.37967491]),
        (1000.0, 6000.0, [20.10538475, -0.0599266107, 2.500069401, -3.99214116e-08, 1.20527214e-11, -1.819015576e-15, 1.078576636e-19], [-744.993961, 4.37918011])]},
    "CO2": {"M": 44.0095, "Hf298": -393510.0, "intervals": [
        (200.0, 1000.0, [49436.5054, -626.411601, 5.30172524, 0.002503813816, -2.127308728e-07, -7.68998878e-10, 2.849677801e-13], [-45281.9846, -7.04827944]),
        (1000.0, 6000.0, [117696.2419, -1788.791477, 8.29152319, -9.22315678e-05, 4.86367688e-09, -1.891053312e-12, 6.33003659e-16], [-39083.5059, -26.52669281])]},
    "H2O": {"M": 18.01528, "Hf298": -241826.0, "intervals": [
        (200.0, 1000.0, [-39479.6083, 575.573102, 0.931782653, 0.00722271286, -7.34255737e-06, 4.95504349e-09, -1.336933246e-12], [-33039.7431, 17.24205775]),
        (1000.0, 6000.0, [1034972.096, -2412.698562, 4.64611078, 0.002291998307, -6.83683048e-07, 9.42646893e-11, -4.82238053e-15], [-13842.86509, -7.97814851])]},
    "He": {"M": 4.002602, "Hf298": 0.0, "intervals": [
        (300.0, 1000.0, [0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0], [-745.375, 0.928723974]),
        (1000.0, 6000.0, [0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0], [-745.375, 0.928723974])]},
    "H2": {"M": 2.01588, "Hf298": 0.0, "intervals": [
        (200.0, 1000.0, [40783.2321, -800.918604, 8.21470201, -0.01269714457, 1.753605076e-05, -1.20286027e-08, 3.36809349e-12], [2682.484665, -30.43788844]),
        (1000.0, 6000.0, [560812.801, -837.150474, 2.975364532, 0.001252249124, -3.74071619e-07, 5.9366252e-11, -3.6069941e-15], [5339.82441, -2.202774769])]},
}

# Dry air, fixed for the whole of Phase 07 (mole fractions).
AIR = {"N2": 0.7808, "O2": 0.2095, "Ar": 0.0093, "CO2": 0.0004}

# Lienhard & Lienhard, A Heat Transfer Textbook, 6th ed. (v6.00, 2024), Appendix A, Table A.6 (air at
# 101 325 Pa): cp in J/(kg K). Their stated uncertainty for gas cp is 0.2-2 %.
LIENHARD_AIR_T = np.array([200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0,
                           800.0, 850.0, 900.0, 950.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0])
LIENHARD_AIR_CP = np.array([1007.0, 1006.0, 1006.0, 1009.0, 1014.0, 1021.0, 1030.0, 1040.0, 1051.0, 1063.0, 1075.0,
                            1087.0, 1099.0, 1110.0, 1121.0, 1131.0, 1141.0, 1159.0, 1174.0, 1188.0, 1200.0, 1210.0])


# ---------------------------------------------------------------- given helpers

def coefficients(species: str, T: float, interval: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Return (a[0:7], b[0:2]) of the NASA interval that contains T (the lower interval at exactly
    1000 K). ``interval`` = 0 or 1 forces one interval (used to test continuity at the break)."""
    ivs = NASA9[species]["intervals"]
    if interval is not None:
        lo, hi, a, b = ivs[interval]
        return np.array(a), np.array(b)
    for lo, hi, a, b in ivs:
        if lo <= T <= hi:
            return np.array(a), np.array(b)
    raise ValueError(f"T = {T} K is outside the NASA range of {species} ({ivs[0][0]}-{ivs[-1][1]} K)")


def _per_T(fn, species: str, T, interval: int | None):
    """Evaluate fn(a, b, T) element by element, picking the interval for each temperature."""
    Ta = np.asarray(T, dtype=float)
    out = np.array([fn(*coefficients(species, t, interval), t) for t in Ta.ravel()]).reshape(Ta.shape)
    return float(out) if out.ndim == 0 else out


def molar_mass(species: str) -> float:
    """Molar mass in kg/mol."""
    return NASA9[species]["M"] * 1e-3


def rms_speed(T: float, M: float) -> float:
    """v_rms = sqrt(3 R T / M) from kinetic theory (P04.07); M in kg/mol."""
    return float(np.sqrt(3.0 * R_U * T / M))


# ---------------------------------------------------------------- the three NASA forms

def cp_over_R(a: np.ndarray, T: float) -> float:
    """cp/R = a1 T^-2 + a2 T^-1 + a3 + a4 T + a5 T^2 + a6 T^3 + a7 T^4."""
    raise NotImplementedError


def h_over_RT(a: np.ndarray, b1: float, T: float) -> float:
    """h/(RT) = -a1 T^-2 + a2 ln(T)/T + a3 + a4 T/2 + a5 T^2/3 + a6 T^3/4 + a7 T^4/5 + b1/T."""
    raise NotImplementedError


def s0_over_R(a: np.ndarray, b2: float, T: float) -> float:
    """s deg/R at 1 bar = -a1 T^-2/2 - a2 T^-1 + a3 ln T + a4 T + a5 T^2/2 + a6 T^3/3 + a7 T^4/4 + b2."""
    raise NotImplementedError


def cp_molar(species: str, T, interval: int | None = None):
    """Molar heat capacity at constant pressure [J/(mol K)]; T may be an array."""
    raise NotImplementedError


def h_molar(species: str, T, interval: int | None = None):
    """Molar enthalpy [J/mol], including the enthalpy of formation (h(298.15 K) = Hf298)."""
    raise NotImplementedError


def s0_molar(species: str, T, interval: int | None = None):
    """Molar entropy at p deg = 1 bar [J/(mol K)]."""
    raise NotImplementedError


def fit_jump(species: str, T_break: float = 1000.0) -> tuple[float, float, float]:
    """Mismatch of the two NASA intervals at the break temperature, as (dcp/R, dh/(R T), ds/R):
    upper-interval value minus lower-interval value, made dimensionless."""
    raise NotImplementedError


# ---------------------------------------------------------------- mixtures (mass basis)

def mixture_molar_mass(x: dict[str, float]) -> float:
    """M_mix = sum x_i M_i [kg/mol] for mole fractions x."""
    raise NotImplementedError


def mixture_cp(x: dict[str, float], T):
    """cp of an ideal-gas mixture per unit mass [J/(kg K)]: sum x_i cp_i / M_mix."""
    raise NotImplementedError


def mixture_h(x: dict[str, float], T):
    """h of an ideal-gas mixture per unit mass [J/kg], on the NASA (formation-enthalpy) datum."""
    raise NotImplementedError


def mixture_s(x: dict[str, float], T, p):
    """s of an ideal-gas mixture per unit mass [J/(kg K)]: sum x_i (s_i deg(T) - R ln(x_i p / p deg)) / M_mix.
    The -R ln x_i part is the entropy of mixing; the -R ln(p/p deg) part is the pressure dependence."""
    raise NotImplementedError


M_AIR = sum(xi * NASA9[sp]["M"] for sp, xi in AIR.items()) * 1e-3     # kg/mol (given, so the module imports)
R_AIR = R_U / M_AIR                                                  # J/(kg K)


def air_cp(T):
    """cp of dry air [J/(kg K)]."""
    raise NotImplementedError


def air_h(T):
    """h of dry air [J/kg] (NASA datum: elements in their reference states have h = 0 at 298.15 K)."""
    raise NotImplementedError


def air_u(T):
    """u = h - R T of dry air [J/kg]."""
    raise NotImplementedError


def air_s(T, p):
    """s of dry air [J/(kg K)] at temperature T and pressure p."""
    raise NotImplementedError


def air_gamma(T):
    """gamma = cp / cv = cp / (cp - R) of dry air."""
    raise NotImplementedError


# ---------------------------------------------------------------- equation of state and kinetic theory

def ideal_gas_pressure(rho: float, T: float, M: float) -> float:
    """p = rho R T / M [Pa]."""
    raise NotImplementedError


def specific_volume(p: float, T: float, M: float) -> float:
    """v = R T / (M p) [m^3/kg]."""
    raise NotImplementedError


def kinetic_pressure(rho: float, T: float, M: float) -> float:
    """Pressure from molecular momentum transfer, p = rho v_rms^2 / 3 (P04.07)."""
    raise NotImplementedError


# ---------------------------------------------------------------- inversion and the state postulate

def temperature_from_h(h: float, x: dict[str, float] = AIR) -> float:
    """Temperature [K] at which the mixture's enthalpy equals h [J/kg], by bracketed root finding
    (h is strictly increasing in T because cp > 0)."""
    raise NotImplementedError


def temperature_from_s(s: float, p: float, x: dict[str, float] = AIR) -> float:
    """Temperature [K] at which the mixture's entropy at pressure p equals s [J/(kg K)]."""
    raise NotImplementedError


def air_state(**given: float) -> dict[str, float]:
    """Complete the state of dry air from exactly two of p [Pa], T [K], v [m^3/kg], h [J/kg],
    u [J/kg], s [J/(kg K)]. Returns all six. Raises ValueError for a pair that is not independent for
    an ideal gas ((T, h), (T, u), (h, u): each is a function of T alone) or for the wrong count."""
    raise NotImplementedError


# ---------------------------------------------------------------- constant cp against variable cp

def delta_h_variable(T1: float, T2: float) -> float:
    """h(T2) - h(T1) of dry air from the NASA polynomials [J/kg]."""
    raise NotImplementedError


def delta_h_constant(T1: float, T2: float, cp: float = CP_CONST) -> float:
    """cp (T2 - T1) with a constant cp [J/kg]."""
    raise NotImplementedError


def exit_temperature(T1: float, q: float, constant_cp: bool = False, cp: float = CP_CONST) -> float:
    """Temperature after adding heat q [J/kg] to air at constant pressure (q = h2 - h1), with the
    variable-cp model (default) or a constant cp."""
    raise NotImplementedError


# ---------------------------------------------------------------- driver

def main(out_dir: str | Path = ".") -> None:
    """Print the verification tables and the engineering answer; write ideal_gas_properties.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    print(f"dry air: M = {M_AIR*1e3:.4f} g/mol, R = {R_AIR:.3f} J/(kg K)")
    print("\nh(298.15 K) against the tabulated enthalpy of formation (J/mol):")
    for sp in NASA9:
        if NASA9[sp]["intervals"][0][0] <= 298.15:
            print(f"  {sp:>4}: h = {h_molar(sp, 298.15):12.2f}, Hf298 = {NASA9[sp]['Hf298']:12.2f}, "
                  f"difference {h_molar(sp, 298.15) - NASA9[sp]['Hf298']:+.2f}")
    print("\nfit continuity at 1000 K (upper minus lower interval): dcp/R, dh/(RT), ds/R")
    for sp in ("N2", "O2", "Ar", "CO2", "H2O", "H2"):
        dcp, dh, ds = fit_jump(sp)
        print(f"  {sp:>4}: {dcp:+.2e} {dh:+.2e} {ds:+.2e}")
    print("\nequipartition limits: cp/R of He at 500 K = "
          f"{cp_molar('He', 500.0)/R_U:.4f} (5/2); N2 at 300 K = {cp_molar('N2', 300.0)/R_U:.4f} (7/2); "
          f"N2 at 3000 K = {cp_molar('N2', 3000.0)/R_U:.4f} (below 9/2: vibration, theta_v = 3390 K, not yet fully excited)")

    print("\nair cp: NASA mixture against Lienhard Table A.6 (J/(kg K))")
    cp_n = air_cp(LIENHARD_AIR_T)
    dev = (cp_n - LIENHARD_AIR_CP) / LIENHARD_AIR_CP
    for T, a, b, d in zip(LIENHARD_AIR_T, cp_n, LIENHARD_AIR_CP, dev):
        if T in (200, 300, 500, 800, 1000, 1200, 1500):
            print(f"  {T:6.0f} K: NASA {a:7.1f}, Lienhard {b:6.0f}, deviation {100*d:+.2f} %")
    print(f"  largest |deviation| over all {len(dev)} rows: {100*np.max(np.abs(dev)):.2f} %")

    rho, T = 1.225, 288.15
    print(f"\nkinetic theory: rho v_rms^2 / 3 = {kinetic_pressure(rho, T, M_AIR):.2f} Pa; rho R T = {ideal_gas_pressure(rho, T, M_AIR):.2f} Pa")

    print("\nstate postulate: (p, T) -> state -> back from other pairs")
    st = air_state(p=P_ATM, T=T)
    print("  " + ", ".join(f"{k} = {v:.6g}" for k, v in st.items()))
    for pair in (("p", "v"), ("p", "h"), ("v", "s"), ("u", "s"), ("p", "s")):
        st2 = air_state(**{k: st[k] for k in pair})
        print(f"  from {pair}: T = {st2['T']:.6f} K, p = {st2['p']:.4f} Pa")
    try:
        air_state(T=T, h=st["h"])
    except ValueError as e:
        print(f"  from ('T', 'h'): ValueError: {e}")

    print("\nconstant cp = 1005 J/(kg K) against variable cp")
    print(f"  cp(T) of air: 300 K {air_cp(300.0):.1f}, 850 K {air_cp(850.0):.1f}, 1700 K {air_cp(1700.0):.1f} J/(kg K); "
          f"gamma: {air_gamma(300.0):.4f}, {air_gamma(850.0):.4f}, {air_gamma(1700.0):.4f}")
    for T1, T2, label in ((288.15, 850.0, "compression 288.15 -> 850 K"),
                          (850.0, 1700.0, "heat addition 850 -> 1700 K"),
                          (1700.0, 1000.0, "turbine expansion 1700 -> 1000 K")):
        dv, dc = delta_h_variable(T1, T2), delta_h_constant(T1, T2)
        print(f"  {label}: variable {dv/1e3:8.2f} kJ/kg, constant {dc/1e3:8.2f} kJ/kg, "
              f"constant-cp error {100*(dc-dv)/dv:+.2f} %")
    q = delta_h_variable(850.0, 1700.0)
    print(f"  the same heat {q/1e3:.2f} kJ/kg added at 850 K with constant cp gives "
          f"{exit_temperature(850.0, q, constant_cp=True):.1f} K instead of {exit_temperature(850.0, q):.1f} K")
    Tm = (850.0 + 1700.0) / 2
    print(f"  mean cp between 850 and 1700 K = {q/850.0:.1f} J/(kg K); cp at the mean temperature {Tm:.0f} K = {air_cp(Tm):.1f}")

    Ts = np.linspace(200.0, 2500.0, 461)
    cps = air_cp(Ts)
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
    axs[0].plot(Ts, cps, label="NASA 9-coefficient mixture")
    axs[0].plot(LIENHARD_AIR_T, LIENHARD_AIR_CP, "o", ms=4, label="Lienhard Table A.6")
    axs[0].axhline(CP_CONST, color="k", ls=":", label="constant 1005")
    axs[0].set_xlabel("T (K)")
    axs[0].set_ylabel("cp (J/(kg K))")
    axs[0].set_title("dry air cp")
    axs[0].legend(fontsize=8)
    axs[1].plot(Ts, air_gamma(Ts))
    axs[1].axhline(1.4, color="k", ls=":")
    axs[1].set_xlabel("T (K)")
    axs[1].set_ylabel("gamma")
    axs[1].set_title("gamma of dry air")
    T2s = np.linspace(400.0, 2000.0, 81)
    err = [100 * (delta_h_constant(298.15, t) - delta_h_variable(298.15, t)) / delta_h_variable(298.15, t) for t in T2s]
    axs[2].plot(T2s, err)
    axs[2].set_xlabel("T2 (K)")
    axs[2].set_ylabel("error of 1005 (T2 - 298.15) (%)")
    axs[2].set_title("constant-cp enthalpy error from 298.15 K")
    fig.tight_layout()
    fig.savefig(out_dir / "ideal_gas_properties.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
