"""P07.08 Exergy: reference solution.

Where does a turbojet core lose the work its fuel could have done? This module audits the exergy of a
single-spool turbojet core, station by station, with real gas properties and real combustion:

* properties of air, Jet-A vapour (C12H23) and the products of its complete combustion from the NASA
  Glenn 9-coefficient polynomials, as ideal-gas mixtures (molar enthalpy including formation, molar
  entropy with the mixing term -R ln(y_i p / p°));
* the engine: inlet (pressure recovery), compressor and turbine (isentropic efficiencies, the turbine
  driving the compressor), combustor (pressure loss, fuel-air ratio from the energy balance to the
  turbine inlet temperature), and a fully expanded nozzle (isentropic efficiency);
* per component: the first law, the entropy generated, the exergy destroyed T0 S_gen, and the same
  destruction from an exergy balance (the Gouy-Stodola theorem says they are equal);
* the fuel exergy (the reversible work of burning the fuel with the engine's air at the dead state),
  the exergy leaving in the hot exhaust, the useful product (the jet's kinetic-energy gain), and the
  audit: exergy in = product + exhaust loss + sum of destructions, closing to round-off.

Assumptions (stated in the lesson): steady flow, adiabatic components, complete combustion to CO2 and
H2O (no dissociation, no CO, no NOx), fuel supplied as Jet-A vapour at T0 and the combustor pressure
(the pump work for liquid fuel and its heat of vaporisation are neglected), ideal-gas mixtures with no
condensation even at the dead state, and a "restricted" dead state: the products' exergy is measured
relative to the products themselves at (T0, p0), not after they have diffused into the atmosphere.

Property data: NASA Glenn 9-coefficient polynomials, copied unrounded from data/thermo.inp of NASA CEA
(https://github.com/nasa/cea, Apache-2.0), documented in McBride, Zehe & Gordon, NASA/TP-2002-211556.
Dry air: N2 0.7808, O2 0.2095, Ar 0.0093, CO2 0.0004 (mole fractions). R = 8.31446261815324 J/(mol K)
(CODATA 2018); s° at p° = 1 bar.

Units: SI; amounts per kilogram of air entering the engine (mol/kg, J/kg, J/(kg K)), so a power in W
is the value times the air mass flow in kg/s.

Run:  python solution.py [out_dir]      (prints the lesson's numbers, writes exergy_audit.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-07/lesson-08-exergy
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

R_UNIVERSAL = 8.31446261815324          # J/(mol K), CODATA 2018 exact
P_STD = 1.0e5                           # Pa, reference pressure of the NASA s° data (1 bar)

# ---------------------------------------------------------------- property data (given)
# NASA Glenn 9-coefficient polynomials (NASA CEA data/thermo.inp, NASA/TP-2002-211556), unrounded.
# Each interval: (T_low, T_high, [a1..a7], [b1, b2]); h includes the enthalpy of formation.
NASA9 = {
    "N2": {"M": 28.0134e-3, "intervals": [
        (200.0, 1000.0, [22103.71497, -381.846182, 6.08273836, -0.00853091441, 1.384646189e-05, -9.62579362e-09, 2.519705809e-12], [710.846086, -10.76003744]),
        (1000.0, 6000.0, [587712.406, -2239.249073, 6.06694922, -0.00061396855, 1.491806679e-07, -1.923105485e-11, 1.061954386e-15], [12832.10415, -15.86640027])]},
    "O2": {"M": 31.9988e-3, "intervals": [
        (200.0, 1000.0, [-34255.6342, 484.700097, 1.119010961, 0.00429388924, -6.83630052e-07, -2.0233727e-09, 1.039040018e-12], [-3391.45487, 18.4969947]),
        (1000.0, 6000.0, [-1037939.022, 2344.830282, 1.819732036, 0.001267847582, -2.188067988e-07, 2.053719572e-11, -8.19346705e-16], [-16890.10929, 17.38716506])]},
    "Ar": {"M": 39.948e-3, "intervals": [
        (200.0, 1000.0, [0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0], [-745.375, 4.37967491]),
        (1000.0, 6000.0, [20.10538475, -0.0599266107, 2.500069401, -3.99214116e-08, 1.20527214e-11, -1.819015576e-15, 1.078576636e-19], [-744.993961, 4.37918011])]},
    "CO2": {"M": 44.0095e-3, "intervals": [
        (200.0, 1000.0, [49436.5054, -626.411601, 5.30172524, 0.002503813816, -2.127308728e-07, -7.68998878e-10, 2.849677801e-13], [-45281.9846, -7.04827944]),
        (1000.0, 6000.0, [117696.2419, -1788.791477, 8.29152319, -9.22315678e-05, 4.86367688e-09, -1.891053312e-12, 6.33003659e-16], [-39083.5059, -26.52669281])]},
    "H2O": {"M": 18.01528e-3, "intervals": [
        (200.0, 1000.0, [-39479.6083, 575.573102, 0.931782653, 0.00722271286, -7.34255737e-06, 4.95504349e-09, -1.336933246e-12], [-33039.7431, 17.24205775]),
        (1000.0, 6000.0, [1034972.096, -2412.698562, 4.64611078, 0.002291998307, -6.83683048e-07, 9.42646893e-11, -4.82238053e-15], [-13842.86509, -7.97814851])]},
    "Jet-A(g)": {"M": 167.31102e-3, "intervals": [
        (273.15, 1000.0, [-606869.559, 8328.25959, -43.1232127, 0.2572390455, -0.000262931604, 1.64498894e-07, -4.64533514e-11], [-76069.6276, 279.4305937]),
        (1000.0, 6000.0, [18583561.02, -76772.1989, 141.9826133, -0.00743752453, 5.85620255e-07, 1.223955647e-11, -3.149201922e-15], [422198.952, -898.606104])]},
}
SPECIES = ("N2", "O2", "Ar", "CO2", "H2O")                       # the gas-mixture species
AIR_X = np.array([0.7808, 0.2095, 0.0093, 0.0004, 0.0])            # Phase 07 dry-air model (mole fractions)
M_SPECIES = np.array([NASA9[k]["M"] for k in SPECIES])
M_AIR = float(AIR_X @ M_SPECIES)                                   # kg/mol, 28.965 g/mol
N_AIR = AIR_X / M_AIR                                              # mol of each species per kg of air
FUEL = "Jet-A(g)"
M_FUEL = NASA9[FUEL]["M"]                                          # kg/mol, C12H23
NU_FUEL = np.array([0.0, -(12.0 + 23.0 / 4.0), 0.0, 12.0, 23.0 / 2.0])   # C12H23 + 17.75 O2 -> 12 CO2 + 11.5 H2O


def _nasa(species: str, T: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """cp/R, h/(RT) and s°/R of one species, choosing the interval by T."""
    T = np.asarray(T, dtype=float)
    out = []
    for (tl, th, a, b) in NASA9[species]["intervals"]:
        a1, a2, a3, a4, a5, a6, a7 = a
        b1, b2 = b
        cp = a1 / T**2 + a2 / T + a3 + a4 * T + a5 * T**2 + a6 * T**3 + a7 * T**4
        h = -a1 / T**2 + a2 * np.log(T) / T + a3 + a4 * T / 2 + a5 * T**2 / 3 + a6 * T**3 / 4 + a7 * T**4 / 5 + b1 / T
        s = -a1 / T**2 / 2 - a2 / T + a3 * np.log(T) + a4 * T + a5 * T**2 / 2 + a6 * T**3 / 3 + a7 * T**4 / 4 + b2
        out.append((cp, h, s))
    lo = T <= 1000.0
    return tuple(np.where(lo, out[0][i], out[1][i]) for i in range(3))


def species_props(species: str, T: float) -> tuple[float, float, float]:
    """Molar cp (J/(mol K)), h (J/mol, including the enthalpy of formation) and s° (J/(mol K), at p°)."""
    cp, h, s = _nasa(species, T)
    return float(cp) * R_UNIVERSAL, float(h) * R_UNIVERSAL * T, float(s) * R_UNIVERSAL


def mixture_enthalpy(n: np.ndarray, T: float) -> float:
    """Enthalpy of a mixture of n_i moles of SPECIES at T, J (per kg of air when n is per kg of air)."""
    return float(sum(ni * species_props(k, T)[1] for k, ni in zip(SPECIES, n) if ni != 0.0))


def mixture_cp(n: np.ndarray, T: float) -> float:
    """Heat capacity of the mixture, J/K (per kg of air)."""
    return float(sum(ni * species_props(k, T)[0] for k, ni in zip(SPECIES, n) if ni != 0.0))


def mixture_s0(n: np.ndarray, T: float) -> float:
    """Sum of n_i s°_i(T), J/K: the mixture entropy without its pressure and mixing terms."""
    return float(sum(ni * species_props(k, T)[2] for k, ni in zip(SPECIES, n) if ni != 0.0))


def mixture_entropy(n: np.ndarray, T: float, p: float) -> float:
    """Entropy of an ideal-gas mixture, S = sum n_i [s°_i(T) - R ln(y_i p / p°)], J/K (per kg of air)."""
    n = np.asarray(n, dtype=float)
    y = n / n.sum()
    mix = sum(ni * np.log(yi * p / P_STD) for ni, yi in zip(n, y) if ni > 0.0)
    return mixture_s0(n, T) - R_UNIVERSAL * float(mix)


def temperature_from_enthalpy(n: np.ndarray, H: float, T_guess: float = 800.0) -> float:
    """Temperature at which mixture_enthalpy(n, T) = H, by Newton's method (dH/dT = mixture_cp > 0)."""
    T = T_guess
    for _ in range(100):
        step = (mixture_enthalpy(n, T) - H) / mixture_cp(n, T)
        T = min(max(T - step, 200.0), 6000.0)
        if abs(step) < 1e-10:
            return T
    raise RuntimeError("temperature inversion did not converge")


def temperature_from_s0(n: np.ndarray, S0: float, T_guess: float = 800.0) -> float:
    """Temperature at which mixture_s0(n, T) = S0, by Newton's method (dS0/dT = mixture_cp / T > 0)."""
    T = T_guess
    for _ in range(100):
        step = (mixture_s0(n, T) - S0) / (mixture_cp(n, T) / T)
        T = min(max(T - step, 200.0), 6000.0)
        if abs(step) < 1e-10:
            return T
    raise RuntimeError("temperature inversion did not converge")


def fuel_enthalpy(T: float) -> float:
    """Molar enthalpy of Jet-A vapour, J/mol (NASA 'Jet-A(g)', formation enthalpy -249 657 J/mol at 298.15 K)."""
    return species_props(FUEL, T)[1]


def fuel_entropy(T: float, p: float) -> float:
    """Molar entropy of pure Jet-A vapour at (T, p), J/(mol K)."""
    return species_props(FUEL, T)[2] - R_UNIVERSAL * np.log(p / P_STD)


# ---------------------------------------------------------------- fuel properties

def lower_heating_value(T: float = 298.15) -> float:
    """Lower heating value of Jet-A vapour, J/kg of fuel: the enthalpy released when it burns completely
    with O2 to CO2 and H2O vapour, reactants and products at the same T."""
    h = np.array([species_props(k, T)[1] for k in SPECIES])
    return float((fuel_enthalpy(T) - NU_FUEL @ h) / M_FUEL)


def standard_reaction_exergy(T: float = 298.15) -> float:
    """-Delta G° of the combustion reaction per kg of fuel, every species pure at p° and T: the
    reversible work of the reaction itself, J/kg."""
    g = np.array([species_props(k, T)[1] - T * species_props(k, T)[2] for k in SPECIES])
    g_fuel = fuel_enthalpy(T) - T * species_props(FUEL, T)[2]
    return float((g_fuel - NU_FUEL @ g) / M_FUEL)


def products(n_fuel: float) -> np.ndarray:
    """Moles of each species per kg of air after complete combustion of n_fuel mol of C12H23 with 1 kg
    of air: CO2 + 12 n_f, H2O + 11.5 n_f, O2 - 17.75 n_f (element conservation)."""
    n = N_AIR + n_fuel * NU_FUEL
    if n[1] < 0.0:
        raise ValueError("not enough oxygen: fuel-air ratio above stoichiometric")
    return n


def element_counts(n: np.ndarray, n_fuel: float = 0.0) -> dict[str, float]:
    """Moles of C, H, O, N and Ar atoms in a mixture n plus n_fuel mol of C12H23."""
    N2, O2, Ar, CO2, H2O = n
    return {"C": CO2 + 12.0 * n_fuel, "H": 2.0 * H2O + 23.0 * n_fuel, "O": 2.0 * O2 + 2.0 * CO2 + H2O,
            "N": 2.0 * N2, "Ar": Ar}


# ---------------------------------------------------------------- the engine

@dataclass
class EngineSpec:
    """Turbojet core inputs (typical, illustrative values as defaults)."""
    T0: float = 288.15        # ambient (dead-state) temperature, K
    p0: float = 101325.0      # ambient (dead-state) pressure, Pa
    V0: float = 0.0           # flight speed, m/s (0: sea-level static)
    pi_d: float = 0.98        # inlet total-pressure recovery
    pi_c: float = 12.0        # compressor pressure ratio
    eta_c: float = 0.85       # compressor isentropic efficiency
    T_tit: float = 1400.0     # turbine inlet temperature, K
    dp_cc: float = 0.04       # combustor fractional pressure loss
    eta_t: float = 0.88       # turbine isentropic efficiency
    eta_n: float = 0.97       # nozzle isentropic efficiency (fully expanded to p0)


def engine_states(spec: EngineSpec) -> dict:
    """Station states per kg of air. Station 0 ambient (static, moving at V0 in the engine frame);
    2 compressor face (total; kinetic energy taken as negligible); 3 compressor exit; 4 turbine inlet;
    5 turbine exit (total); 9 nozzle exit (static, at p0, velocity V9). Returns a dict with:
    'T', 'p' (dicts by station), 'n_air' (air moles), 'n_prod' (product moles), 'n_fuel' (mol of fuel),
    'f' (fuel-air mass ratio), 'w_c' (J/kg air), 'V9' (m/s).
    Inlet: h2 = h0 + V0^2/2 (adiabatic); p2 = pi_d times the isentropic ram pressure, which satisfies
    S0(T2) - S0(T0) = R_mix ln(p2s/p0). Compressor: isentropic efficiency. Combustor: n_fuel from the energy balance
    H_air(T3) + n_fuel h_fuel(T0) = H_prod(n_fuel, T_tit), linear in n_fuel because the product
    moles are. Turbine: (power) H_prod(T4) - H_prod(T5) = w_c; p5 from the isentropic exit T5s.
    Nozzle: T9 from h9 = h5 - eta_n (h5 - h9s), V9 = sqrt(2 (h5 - h9) / mass of products)."""
    s = spec
    na = N_AIR
    T, p = {0: s.T0}, {0: s.p0}
    # inlet: adiabatic, total enthalpy conserved; ideal ram pressure from the isentrope, then recovery
    H0 = mixture_enthalpy(na, s.T0)
    T[2] = temperature_from_enthalpy(na, H0 + 0.5 * s.V0 ** 2, s.T0)
    R_a = R_UNIVERSAL * na.sum()
    p[2] = s.pi_d * s.p0 * np.exp((mixture_s0(na, T[2]) - mixture_s0(na, s.T0)) / R_a)
    # compressor
    p[3] = s.pi_c * p[2]
    T3s = temperature_from_s0(na, mixture_s0(na, T[2]) + R_a * np.log(s.pi_c), T[2])
    H2 = mixture_enthalpy(na, T[2])
    w_c = (mixture_enthalpy(na, T3s) - H2) / s.eta_c
    T[3] = temperature_from_enthalpy(na, H2 + w_c, T3s)
    # combustor: H_air(T3) + n_f h_f(T0) = sum (N_AIR + n_f nu) h(T4)  ->  linear in n_f
    p[4] = p[3] * (1.0 - s.dp_cc)
    T[4] = s.T_tit
    h4 = np.array([species_props(k, s.T_tit)[1] for k in SPECIES])
    H3 = mixture_enthalpy(na, T[3])
    n_fuel = (na @ h4 - H3) / (fuel_enthalpy(s.T0) - NU_FUEL @ h4)
    npr = products(n_fuel)
    R_p = R_UNIVERSAL * npr.sum()
    # turbine drives the compressor
    H4 = mixture_enthalpy(npr, T[4])
    T[5] = temperature_from_enthalpy(npr, H4 - w_c, T[4])
    T5s = temperature_from_enthalpy(npr, H4 - w_c / s.eta_t, T[5])
    p[5] = p[4] * np.exp((mixture_s0(npr, T5s) - mixture_s0(npr, T[4])) / R_p)
    if p[5] <= s.p0:
        raise ValueError("turbine exit pressure below ambient: no nozzle pressure ratio left")
    # nozzle, fully expanded to p0
    p[9] = s.p0
    H5 = mixture_enthalpy(npr, T[5])
    T9s = temperature_from_s0(npr, mixture_s0(npr, T[5]) + R_p * np.log(s.p0 / p[5]), T[5])
    H9 = H5 - s.eta_n * (H5 - mixture_enthalpy(npr, T9s))
    T[9] = temperature_from_enthalpy(npr, H9, T9s)
    m_p = 1.0 + n_fuel * M_FUEL
    V9 = float(np.sqrt(2.0 * (H5 - H9) / m_p))
    return {"T": T, "p": p, "n_air": na, "n_prod": npr, "n_fuel": float(n_fuel),
            "f": float(n_fuel * M_FUEL), "w_c": float(w_c), "V9": V9}


def first_law_residuals(spec: EngineSpec) -> dict[str, float]:
    """Energy in minus energy out (J/kg air) of each component, from engine_states: inlet
    (H0 + V0^2/2 - H2), compressor (H2 + w_c - H3), combustor (H3 + n_f h_f - H4), turbine
    (H4 - w_c - H5), nozzle (H5 - H9 - m_p V9^2/2). All zero to round-off."""
    st = engine_states(spec)
    T, na, npr, nf = st["T"], st["n_air"], st["n_prod"], st["n_fuel"]
    m_p = 1.0 + st["f"]
    return {
        "inlet": mixture_enthalpy(na, T[0]) + 0.5 * spec.V0 ** 2 - mixture_enthalpy(na, T[2]),
        "compressor": mixture_enthalpy(na, T[2]) + st["w_c"] - mixture_enthalpy(na, T[3]),
        "combustor": mixture_enthalpy(na, T[3]) + nf * fuel_enthalpy(spec.T0) - mixture_enthalpy(npr, T[4]),
        "turbine": mixture_enthalpy(npr, T[4]) - st["w_c"] - mixture_enthalpy(npr, T[5]),
        "nozzle": mixture_enthalpy(npr, T[5]) - mixture_enthalpy(npr, T[9]) - 0.5 * m_p * st["V9"] ** 2,
    }


def flow_exergy(n: np.ndarray, T: float, p: float, T0: float, p0: float, V: float = 0.0, m: float = 1.0) -> float:
    """Thermomechanical flow exergy of a stream of fixed composition n (mass m), relative to the same
    mixture at the dead state (T0, p0): (H - H0) - T0 (S - S0) + m V^2/2, J (per kg of air)."""
    return (mixture_enthalpy(n, T) - mixture_enthalpy(n, T0)
            - T0 * (mixture_entropy(n, T, p) - mixture_entropy(n, T0, p0)) + 0.5 * m * V ** 2)


def fuel_exergy(spec: EngineSpec, st: dict | None = None) -> float:
    """Exergy supplied by the fuel, J/kg air: the reversible work of turning the air and fuel at the dead
    state (fuel at T0 and the combustor pressure) into the products at the dead state,
    B_air(T0, p0) + n_f b_fuel(T0, p3) - B_prod(T0, p0), with B = H - T0 S."""
    st = engine_states(spec) if st is None else st
    T0, p0 = spec.T0, spec.p0
    B_air = mixture_enthalpy(st["n_air"], T0) - T0 * mixture_entropy(st["n_air"], T0, p0)
    b_f = fuel_enthalpy(T0) - T0 * fuel_entropy(T0, st["p"][3])
    B_prod = mixture_enthalpy(st["n_prod"], T0) - T0 * mixture_entropy(st["n_prod"], T0, p0)
    return float(B_air + st["n_fuel"] * b_f - B_prod)


def exergy_audit(spec: EngineSpec) -> dict:
    """The exergy audit, J/kg of air. Returns a dict with:
    'S_gen' and 'D' (dicts by component: inlet, compressor, combustor, turbine, nozzle) with
    D = T0 S_gen; 'D_balance' (the same destructions from exergy balances, exergy in - exergy out - work
    out, which Gouy-Stodola says equal T0 S_gen); 'Ex_fuel'; 'Ex_ram' = V0^2/2 (the kinetic exergy of
    the incoming air in the engine frame); 'product' = m_p V9^2/2 - V0^2/2 (the jet's kinetic-energy
    gain); 'exhaust' (thermomechanical exergy left in the jet, which the atmosphere destroys);
    'residual' = Ex_fuel + Ex_ram - product - exhaust - sum(D), which is zero to round-off;
    'eta_II' = product / (Ex_fuel + Ex_ram); 'eta_th' = product / (m_f LHV); 'KE9' = m_p V9^2/2; and
    component second-law efficiencies 'eta_II_comp': compressor (exergy gained / work in), turbine
    (work out / exergy given up), nozzle (kinetic energy out / exergy given up); 'states' is engine_states."""
    st = engine_states(spec)
    T, p, na, npr, nf = st["T"], st["p"], st["n_air"], st["n_prod"], st["n_fuel"]
    T0, p0 = spec.T0, spec.p0
    m_p = 1.0 + st["f"]
    S = {k: mixture_entropy(na, T[k], p[k]) for k in (0, 2, 3)}
    S.update({k: mixture_entropy(npr, T[k], p[k]) for k in (4, 5, 9)})
    s_fuel = nf * fuel_entropy(T0, p[3])
    S_gen = {"inlet": S[2] - S[0], "compressor": S[3] - S[2], "combustor": S[4] - S[3] - s_fuel,
             "turbine": S[5] - S[4], "nozzle": S[9] - S[5]}
    D = {k: T0 * v for k, v in S_gen.items()}
    # exergy-balance route: psi relative to each stream's own dead state, fuel exergy across the combustor
    psi_a = {k: flow_exergy(na, T[k], p[k], T0, p0) for k in (0, 2, 3)}
    psi_p = {k: flow_exergy(npr, T[k], p[k], T0, p0) for k in (4, 5, 9)}
    Ex_f = fuel_exergy(spec, st)
    w_c = st["w_c"]
    KE9 = 0.5 * m_p * st["V9"] ** 2
    Ex_ram = 0.5 * spec.V0 ** 2
    D_bal = {"inlet": psi_a[0] + Ex_ram - psi_a[2],
             "compressor": psi_a[2] + w_c - psi_a[3],
             "combustor": psi_a[3] + Ex_f - psi_p[4],
             "turbine": psi_p[4] - w_c - psi_p[5],
             "nozzle": psi_p[5] - psi_p[9] - KE9}
    product = KE9 - Ex_ram
    exhaust = psi_p[9]
    residual = Ex_f + Ex_ram - KE9 - exhaust - sum(D.values())
    lhv = lower_heating_value(298.15)
    eta_comp = {"compressor": (psi_a[3] - psi_a[2]) / w_c,
                "turbine": w_c / (psi_p[4] - psi_p[5]),
                "nozzle": KE9 / (psi_p[5] - psi_p[9])}
    return {"S_gen": S_gen, "D": D, "D_balance": D_bal, "Ex_fuel": Ex_f, "Ex_ram": Ex_ram,
            "product": product, "KE9": KE9, "exhaust": exhaust, "residual": float(residual),
            "eta_II": product / (Ex_f + Ex_ram), "eta_th": product / (st["f"] * lhv),
            "eta_II_comp": eta_comp, "states": st}


# ---------------------------------------------------------------- the lesson's numbers

def _print_audit(spec: EngineSpec, label: str, mdot: float) -> dict:
    a = exergy_audit(spec)
    st = a["states"]
    tot = a["Ex_fuel"] + a["Ex_ram"]
    print(f"\n{label}")
    print(f"  T0 = {spec.T0} K, p0 = {spec.p0:.0f} Pa, V0 = {spec.V0:.1f} m/s; pi_d = {spec.pi_d}, pi_c = {spec.pi_c}, eta_c = {spec.eta_c}, "
          f"TIT = {spec.T_tit} K, dp_cc = {spec.dp_cc}, eta_t = {spec.eta_t}, eta_n = {spec.eta_n}")
    for k in (0, 2, 3, 4, 5, 9):
        print(f"    station {k}: T = {st['T'][k]:8.2f} K, p = {st['p'][k] / 1e3:8.2f} kPa")
    print(f"  f = {st['f']:.5f} (fuel-air ratio), w_c = {st['w_c'] / 1e3:.2f} kJ/kg, V9 = {st['V9']:.1f} m/s, "
          f"specific thrust (m_p V9 - V0) = {(1 + st['f']) * st['V9'] - spec.V0:.1f} N s/kg")
    print(f"  fuel exergy {a['Ex_fuel'] / 1e3:.2f} kJ/kg air = {a['Ex_fuel'] / st['f'] / 1e6:.3f} MJ/kg fuel; ram exergy {a['Ex_ram'] / 1e3:.2f} kJ/kg")
    print("  where the exergy goes (kJ/kg air, % of exergy in), destruction by T0 S_gen and by exergy balance:")
    print(f"    {'product: jet KE gain':24s} {a['product'] / 1e3:8.2f}  {a['product'] / tot * 100:5.1f} %")
    print(f"    {'exhaust (lost)':24s} {a['exhaust'] / 1e3:8.2f}  {a['exhaust'] / tot * 100:5.1f} %")
    for k in a["D"]:
        print(f"    {'destroyed: ' + k:24s} {a['D'][k] / 1e3:8.2f}  {a['D'][k] / tot * 100:5.1f} %   (balance {a['D_balance'][k] / 1e3:8.2f}, "
              f"S_gen {a['S_gen'][k]:.2f} J/(kg K))")
    print(f"  Gouy-Stodola closure: in - product - exhaust - sum D = {a['residual']:.2e} J/kg "
          f"(relative {a['residual'] / tot:.1e}); max |T0 S_gen - balance| = {max(abs(a['D'][k] - a['D_balance'][k]) for k in a['D']):.1e} J/kg")
    print(f"  eta_II = {a['eta_II']:.4f}, eta_th = {a['eta_th']:.4f}; component eta_II: " +
          ", ".join(f"{k} {v:.4f}" for k, v in a["eta_II_comp"].items()))
    print(f"  at {mdot} kg/s of air: fuel {mdot * st['f']:.4f} kg/s, jet power {mdot * a['product'] / 1e6:.3f} MW, "
          f"combustor destruction {mdot * a['D']['combustor'] / 1e6:.3f} MW")
    return a


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save exergy_audit.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    lhv = lower_heating_value(298.15)
    exs = standard_reaction_exergy(298.15)
    print(f"Jet-A(g) C12H23 (NASA): LHV(298.15 K) = {lhv / 1e6:.3f} MJ/kg, -dG°(298.15 K) = {exs / 1e6:.3f} MJ/kg, ratio {exs / lhv:.4f}")
    print(f"stoichiometric fuel-air ratio = {N_AIR[1] / (12 + 23 / 4) * M_FUEL:.5f}")

    mdot = 20.0                                                     # kg/s, illustrative turbojet core
    base = EngineSpec()
    a = _print_audit(base, "turbojet core, sea-level static (typical, illustrative)", mdot)
    print(f"  fuel exergy per kg fuel / LHV = {a['Ex_fuel'] / a['states']['f'] / lhv:.4f}")

    rev = EngineSpec(pi_d=1.0, eta_c=1.0, eta_t=1.0, eta_n=1.0, dp_cc=0.0)
    r = _print_audit(rev, "reversible turbomachinery limit (pi_d = 1, eta = 1, no combustor pressure loss)", mdot)

    cruise = EngineSpec(T0=216.65, p0=22632.0, V0=0.8 * np.sqrt(1.4 * R_UNIVERSAL / M_AIR * 216.65))
    _print_audit(cruise, "cruise, 11 km ISA, Mach 0.8 (typical, illustrative)", mdot)

    print("\nsweep at sea-level static: TIT (K), pi_c, eta_II, eta_th, combustor D / exergy in, exhaust / exergy in")
    rows = []
    for tit in (1200.0, 1400.0, 1600.0):
        for pc in (6.0, 12.0, 20.0, 30.0):
            aa = exergy_audit(EngineSpec(T_tit=tit, pi_c=pc))
            tot = aa["Ex_fuel"] + aa["Ex_ram"]
            rows.append((tit, pc, aa))
            print(f"  {tit:6.0f}  {pc:4.0f}   {aa['eta_II']:.4f}   {aa['eta_th']:.4f}   {aa['D']['combustor'] / tot:.4f}   {aa['exhaust'] / tot:.4f}")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    labels = ["product (jet KE)", "exhaust loss", "D inlet", "D compressor", "D combustor", "D turbine", "D nozzle"]
    for i, (lab, aa) in enumerate((("sea-level static", a), ("reversible turbomachinery", r))):
        vals = [aa["product"], aa["exhaust"]] + [aa["D"][k] for k in ("inlet", "compressor", "combustor", "turbine", "nozzle")]
        tot = aa["Ex_fuel"] + aa["Ex_ram"]
        left = 0.0
        for j, v in enumerate(vals):
            ax[0].barh(i, v / tot * 100, left=left, color=f"C{j}", label=labels[j] if i == 0 else None)
            left += v / tot * 100
    ax[0].set_yticks([0, 1])
    ax[0].set_yticklabels(["sea-level static", "reversible\nturbomachinery"])
    ax[0].set_xlabel("% of exergy in")
    ax[0].set_title("turbojet core exergy audit (illustrative)")
    ax[0].legend(fontsize=7, loc="lower right")
    for tit in (1200.0, 1400.0, 1600.0):
        pcs = [rr[1] for rr in rows if rr[0] == tit]
        ax[1].plot(pcs, [rr[2]["eta_II"] for rr in rows if rr[0] == tit], "o-", label=f"eta_II, TIT {tit:.0f} K")
        ax[1].plot(pcs, [rr[2]["D"]["combustor"] / (rr[2]["Ex_fuel"] + rr[2]["Ex_ram"]) for rr in rows if rr[0] == tit], "s--",
                   label=f"combustor D share, TIT {tit:.0f} K")
    ax[1].set_xlabel("compressor pressure ratio")
    ax[1].set_ylabel("fraction of exergy in")
    ax[1].set_title("second-law efficiency and combustor destruction")
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out_dir / "exergy_audit.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'exergy_audit.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
