"""P07.09 Gas Mixtures & Psychrometrics: reference solution.

What humidity and dew point does an airliner cabin settle at in cruise, when every kilogram of
fresh air comes from 11 km where the air holds almost no water? This module builds the answer from
the ideal-gas mixture rules and the two conservation laws of moist air:

* ideal-gas mixtures: molar mass, mole and mass fractions, Dalton's partial pressures and Amagat's
  partial volumes (identical for ideal gases), mixture enthalpy and entropy from NASA polynomials,
  and the ideal entropy of mixing -R sum x ln x;
* water: the saturation pressure from the steam table of P07.03 (NIST SRD 69 values via Yan,
  Appendix A), checked against an independent estimate from the NASA polynomials, the Gibbs
  equality g_vapour(T, p_sat) = g_liquid(T);
* moist air: humidity ratio, relative humidity, dew point, enthalpy per kilogram of dry air,
  adiabatic mixing of two streams, cooling with condensation;
* the environmental control system (ECS): outside air at cruise, bleed, pack, mixing manifold with
  recirculated air, cabin, outflow; dry-air, water and energy balances on every box.

Units: SI (K, Pa, kg, J, W). Mass flows are of DRY air (kg/s) unless named otherwise; humidity
ratios are kg of water vapour per kg of dry air; enthalpies of moist air are J per kg of dry air.
Pressures: the NASA data are at p_ref = 1 bar; the ISA sea-level pressure is 101 325 Pa.

Run:  python solution.py [out_dir]   (prints the lesson's numbers, writes psychrometrics.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-07/lesson-09-gas-mixtures-psychrometrics
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

R_U = 8.31446261815324          # J/(mol K), CODATA 2018 (exact)
P_REF = 1.0e5                   # Pa, standard-state pressure of the NASA data
P_ATM = 101325.0                # Pa
G0 = 9.80665                    # m/s^2

# Dry air model fixed for Phase 07 (mole fractions; M = 28.965 g/mol)
AIR_X = {"N2": 0.7808, "O2": 0.2095, "Ar": 0.0093, "CO2": 0.0004}

# NASA Glenn 9-coefficient polynomials, copied unrounded from data/thermo.inp of NASA CEA
# (https://github.com/nasa/cea, Apache-2.0; McBride, Zehe & Gordon, NASA/TP-2002-211556).
# Only the 200-1000 K intervals are needed here (H2O(L): 273.15-373.15 K, and 373.15-600 K).
# cp/R = a1 T^-2 + a2 T^-1 + a3 + a4 T + a5 T^2 + a6 T^3 + a7 T^4; h includes Hf298; s at p = 1 bar.
NASA9: dict[str, dict] = {
    'N2': {
        "M": 28.0134,  # g/mol
        "Hf298": 0.0,  # J/mol
        "elements": {'N': 2},
        "intervals": [
            (200.0, 1000.0, [22103.71497, -381.846182, 6.08273836, -0.00853091441, 1.384646189e-05, -9.62579362e-09, 2.519705809e-12], [710.846086, -10.76003744]),
        ],
    },
    'O2': {
        "M": 31.9988,  # g/mol
        "Hf298": 0.0,  # J/mol
        "elements": {'O': 2},
        "intervals": [
            (200.0, 1000.0, [-34255.6342, 484.700097, 1.119010961, 0.00429388924, -6.83630052e-07, -2.0233727e-09, 1.039040018e-12], [-3391.45487, 18.4969947]),
        ],
    },
    'Ar': {
        "M": 39.948,  # g/mol
        "Hf298": 0.0,  # J/mol
        "elements": {'Ar': 1},
        "intervals": [
            (200.0, 1000.0, [0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0], [-745.375, 4.37967491]),
        ],
    },
    'CO2': {
        "M": 44.0095,  # g/mol
        "Hf298": -393510.0,  # J/mol
        "elements": {'C': 1, 'O': 2},
        "intervals": [
            (200.0, 1000.0, [49436.5054, -626.411601, 5.30172524, 0.002503813816, -2.127308728e-07, -7.68998878e-10, 2.849677801e-13], [-45281.9846, -7.04827944]),
        ],
    },
    'H2O': {
        "M": 18.01528,  # g/mol
        "Hf298": -241826.0,  # J/mol
        "elements": {'H': 2, 'O': 1},
        "intervals": [
            (200.0, 1000.0, [-39479.6083, 575.573102, 0.931782653, 0.00722271286, -7.34255737e-06, 4.95504349e-09, -1.336933246e-12], [-33039.7431, 17.24205775]),
        ],
    },
    'H2O(L)': {
        "M": 18.01528,  # g/mol
        "Hf298": -285830.0,  # J/mol
        "elements": {'H': 2, 'O': 1},
        "intervals": [
            (273.15, 373.15, [1326371304.0, -24482953.88, 187942.8776, -767.899505, 1.761556813, -0.002151167128, 1.092570813e-06], [110176047.6, -977970.097]),
            (373.15, 600.0, [1263631001.0, -16803802.49, 92782.3479, -272.237395, 0.447924376, -0.000391939743, 1.425743266e-07], [81131768.8, -513441.808]),
        ],
    },
}

# Saturated water, temperature table: T (deg C), p_sat (MPa), h_f (kJ/kg), h_g (kJ/kg).
# Yan, Introduction to Engineering Thermodynamics (LibreTexts), Appendix A, Table A1; values generated
# from the NIST Chemistry WebBook (SRD 69, Lemmon, McLinden & Friend). Reference state: u_f = s_f = 0
# for liquid at the triple point. The printed pressures have a resolution of 10 Pa.
SAT_WATER = np.array([
    (0.01, 0.00061, 0.00, 2500.92),
    (5.00, 0.00087, 21.02, 2510.06),
    (10.00, 0.00123, 42.02, 2519.21),
    (15.00, 0.00171, 62.98, 2528.33),
    (20.00, 0.00234, 83.91, 2537.43),
    (25.00, 0.00317, 104.83, 2546.51),
    (30.00, 0.00425, 125.73, 2555.55),
    (35.00, 0.00563, 146.63, 2564.55),
    (40.00, 0.00738, 167.53, 2573.51),
    (45.00, 0.00960, 188.44, 2582.43),
    (50.00, 0.01235, 209.34, 2591.29),
    (55.00, 0.01576, 230.26, 2600.09),
    (60.00, 0.01995, 251.18, 2608.83),
])
T_TRIPLE = 273.16               # K
M_W = NASA9["H2O"]["M"] / 1000.0   # kg/mol


# ---------------------------------------------------------------- given: NASA polynomials and the ISA

def _coeffs(species: str, T: float) -> tuple[list[float], list[float]]:
    """The (a1..a7, b1 b2) coefficients of the interval holding T. Raises ValueError outside the fits."""
    for lo, hi, a, b in NASA9[species]["intervals"]:
        if lo <= T <= hi:
            return a, b
    raise ValueError(f"T = {T} K outside the fitted range for {species}")


def cp_species(species: str, T: float) -> float:
    """Molar heat capacity cp(T) [J/(mol K)] from the NASA 9-coefficient polynomial."""
    a, _ = _coeffs(species, T)
    return R_U * (a[0] / T ** 2 + a[1] / T + a[2] + a[3] * T + a[4] * T ** 2 + a[5] * T ** 3 + a[6] * T ** 4)


def h_species(species: str, T: float) -> float:
    """Molar enthalpy h(T) [J/mol] including the enthalpy of formation (h(298.15 K) = Hf298)."""
    a, b = _coeffs(species, T)
    return R_U * T * (-a[0] / T ** 2 + a[1] * np.log(T) / T + a[2] + a[3] * T / 2 + a[4] * T ** 2 / 3
                      + a[5] * T ** 3 / 4 + a[6] * T ** 4 / 5 + b[0] / T)


def s_species(species: str, T: float) -> float:
    """Molar standard-state entropy s°(T) [J/(mol K)] at p_ref = 1 bar."""
    a, b = _coeffs(species, T)
    return R_U * (-a[0] / T ** 2 / 2 - a[1] / T + a[2] * np.log(T) + a[3] * T + a[4] * T ** 2 / 2
                  + a[5] * T ** 3 / 3 + a[6] * T ** 4 / 4 + b[1])


def g_species(species: str, T: float) -> float:
    """Molar standard-state Gibbs energy g°(T) = h - T s° [J/mol] at 1 bar."""
    return h_species(species, T) - T * s_species(species, T)


def isa_state(z: float) -> tuple[float, float]:
    """ISA temperature [K] and pressure [Pa] at geopotential altitude z [m], 0-20 km (troposphere
    with lapse 6.5 K/km, isothermal tropopause above 11 km), with the ICAO R_air = 287.05287 J/(kg K)."""
    R_air, L, T0 = 287.05287, 0.0065, 288.15
    if z <= 11000.0:
        T = T0 - L * z
        return T, P_ATM * (T / T0) ** (G0 / (R_air * L))
    T11, p11 = isa_state(11000.0)
    return T11, p11 * np.exp(-G0 * (z - 11000.0) / (R_air * T11))


# ---------------------------------------------------------------- ideal-gas mixtures

def molar_mass(x: dict[str, float]) -> float:
    """Mixture molar mass M = sum x_i M_i [kg/mol] for mole fractions x (normalised here)."""
    tot = sum(x.values())
    return sum(xi * NASA9[s]["M"] for s, xi in x.items()) / tot / 1000.0


def mass_fractions(x: dict[str, float]) -> dict[str, float]:
    """Mass fractions y_i = x_i M_i / M from mole fractions x."""
    M = molar_mass(x) * 1000.0
    tot = sum(x.values())
    return {s: xi / tot * NASA9[s]["M"] / M for s, xi in x.items()}


def mole_fractions(y: dict[str, float]) -> dict[str, float]:
    """Mole fractions x_i = (y_i / M_i) / sum(y_j / M_j) from mass fractions y."""
    n = {s: yi / NASA9[s]["M"] for s, yi in y.items()}
    tot = sum(n.values())
    return {s: ni / tot for s, ni in n.items()}


def partial_pressures(n: dict[str, float], T: float, V: float) -> dict[str, float]:
    """Dalton: the pressure each component would exert ALONE in the whole volume V at T,
    p_i = n_i R T / V [Pa], for amounts n [mol]. For ideal gases they add to the mixture pressure."""
    return {s: ni * R_U * T / V for s, ni in n.items()}


def partial_volumes(n: dict[str, float], T: float, p: float) -> dict[str, float]:
    """Amagat: the volume each component would occupy ALONE at the mixture's T and p,
    V_i = n_i R T / p [m^3]. For ideal gases they add to the mixture volume."""
    return {s: ni * R_U * T / p for s, ni in n.items()}


def mixture_cp(x: dict[str, float], T: float) -> float:
    """Molar cp of an ideal mixture, sum x_i cp_i(T) [J/(mol K)]."""
    tot = sum(x.values())
    return sum(xi * cp_species(s, T) for s, xi in x.items()) / tot


def mixture_h(x: dict[str, float], T: float) -> float:
    """Molar enthalpy of an ideal mixture, sum x_i h_i(T) [J/mol]: no enthalpy of mixing."""
    tot = sum(x.values())
    return sum(xi * h_species(s, T) for s, xi in x.items()) / tot


def mixture_s(x: dict[str, float], T: float, p: float) -> float:
    """Molar entropy of an ideal mixture at T, p: sum x_i [s°_i(T) - R ln(x_i p / p_ref)] [J/(mol K)].
    Each component is evaluated at its own partial pressure x_i p."""
    tot = sum(x.values())
    out = 0.0
    for s, xi in x.items():
        if xi > 0.0:
            xf = xi / tot
            out += xf * (s_species(s, T) - R_U * np.log(xf * p / P_REF))
    return out


def entropy_of_mixing(x: dict[str, float]) -> float:
    """Ideal entropy of mixing per mole of mixture, -R sum x_i ln x_i [J/(mol K)], when pure components
    at the same T and p are mixed to the same T and p. Zero for a pure gas; R ln N for N equal parts."""
    tot = sum(x.values())
    return -R_U * sum((xi / tot) * np.log(xi / tot) for xi in x.values() if xi > 0.0)


# ---------------------------------------------------------------- water: saturation pressure

def p_sat(T: float) -> float:
    """Saturation pressure of water [Pa] from the steam table (P07.03), 273.16-333.15 K, by linear
    interpolation of ln p against 1/T between table rows (exact at the rows; Clausius-Clapeyron says
    ln p is nearly linear in 1/T). Raises ValueError outside the table."""
    Tk = SAT_WATER[:, 0] + 273.15
    Tk[0] = T_TRIPLE
    if not (Tk[0] - 1e-9 <= T <= Tk[-1] + 1e-9):
        raise ValueError(f"T = {T} K outside the saturation table (273.16-333.15 K)")
    x = 1.0 / Tk[::-1]
    y = np.log(SAT_WATER[::-1, 1] * 1e6)
    return float(np.exp(np.interp(1.0 / T, x, y)))


def p_sat_gibbs(T: float) -> float:
    """Independent estimate of the saturation pressure [Pa] from the NASA polynomials alone: at
    saturation g_vapour(T, p_sat) = g_liquid(T), and for an ideal-gas vapour and an incompressible
    liquid (whose g barely depends on p) that gives p_sat = p_ref exp(-(g°_vap - g°_liq) / (R T)).
    Valid where both fits are, 273.15-600 K."""
    dg = g_species("H2O", T) - g_species("H2O(L)", T)
    return float(P_REF * np.exp(-dg / (R_U * T)))


def latent_heat_triple() -> float:
    """h_fg at the triple point from the table [J/kg]."""
    return float((SAT_WATER[0, 3] - SAT_WATER[0, 2]) * 1e3)


def p_sat_extrapolated(T: float) -> float:
    """ESTIMATE of the saturation pressure over (supercooled) liquid water below the triple point [Pa],
    by Clausius-Clapeyron with constant h_fg from the triple point:
    ln(p / p_tp) = -(h_fg M_w / R)(1/T - 1/T_tp). Not checked against data (no verified sub-zero
    table is in this course's sources); the frost point over ice is a little higher."""
    p_tp = SAT_WATER[0, 1] * 1e6
    return float(p_tp * np.exp(-latent_heat_triple() * M_W / R_U * (1.0 / T - 1.0 / T_TRIPLE)))


def p_sat_any(T: float) -> float:
    """Saturation pressure from the table at and above the triple point, the Clausius-Clapeyron
    estimate below it (continuous at 273.16 K)."""
    return p_sat(T) if T >= T_TRIPLE else p_sat_extrapolated(T)


# ---------------------------------------------------------------- moist air

M_AIR = sum(xi * NASA9[k]["M"] for k, xi in AIR_X.items()) / 1000.0   # kg/mol, 28.9657 g/mol
EPS = M_W / M_AIR                 # 0.62195: ratio of the molar masses of water and dry air


def humidity_ratio(p_v: float, p: float) -> float:
    """omega = m_v / m_a = eps p_v / (p - p_v) [kg water / kg dry air], from Dalton (p = p_a + p_v)."""
    return EPS * p_v / (p - p_v)


def vapour_pressure(omega: float, p: float) -> float:
    """Partial pressure of water vapour p_v = omega p / (eps + omega) [Pa] (inverse of humidity_ratio)."""
    return omega * p / (EPS + omega)


def relative_humidity(T: float, omega: float, p: float) -> float:
    """phi = p_v / p_sat(T) (dimensionless; 1 = saturated). Uses the table at and above 273.16 K and
    the Clausius-Clapeyron estimate below."""
    return vapour_pressure(omega, p) / p_sat_any(T)


def omega_saturated(T: float, p: float) -> float:
    """Humidity ratio of saturated moist air at T and total pressure p."""
    return humidity_ratio(p_sat_any(T), p)


def dew_point(omega: float, p: float) -> tuple[float, bool]:
    """Dew point [K]: the temperature at which air of humidity ratio omega at total pressure p becomes
    saturated when cooled at constant p and omega, p_sat(T_dp) = p_v. Returns (T_dp, is_estimate):
    is_estimate is True when T_dp lies below the triple point and comes from p_sat_extrapolated."""
    pv = vapour_pressure(omega, p)
    p_tp = SAT_WATER[0, 1] * 1e6
    if pv >= p_tp:
        return float(brentq(lambda T: p_sat(T) - pv, T_TRIPLE, SAT_WATER[-1, 0] + 273.15, xtol=1e-12)), False
    # invert the Clausius-Clapeyron estimate in closed form
    inv = 1.0 / T_TRIPLE - R_U / (latent_heat_triple() * M_W) * np.log(pv / p_tp)
    return float(1.0 / inv), True


def dry_air_enthalpy(T: float) -> float:
    """Specific enthalpy of dry air [J/kg] relative to 0 deg C (273.15 K), from the NASA mixture."""
    return (mixture_h(AIR_X, T) - mixture_h(AIR_X, 273.15)) / M_AIR


def vapour_enthalpy(T: float) -> float:
    """Specific enthalpy of water vapour as an ideal gas [J/kg] on the steam-table reference (liquid
    at the triple point = 0): h_g(0.01 deg C) from the table plus the NASA ideal-gas change from 273.16 K."""
    return SAT_WATER[0, 3] * 1e3 + (h_species("H2O", T) - h_species("H2O", T_TRIPLE)) / M_W


def liquid_enthalpy(T: float) -> float:
    """Specific enthalpy of liquid water [J/kg] from the table's h_f (linear interpolation), 0.01-60 deg C."""
    return float(np.interp(T - 273.15, SAT_WATER[:, 0], SAT_WATER[:, 2]) * 1e3)


def moist_air_enthalpy(T: float, omega: float) -> float:
    """Enthalpy of moist air per kilogram of DRY air [J/kg_a]: h = h_a(T) + omega h_v(T)."""
    return dry_air_enthalpy(T) + omega * vapour_enthalpy(T)


def temperature_from_enthalpy(h: float, omega: float) -> float:
    """Invert moist_air_enthalpy for T [K] (200-1000 K) at fixed omega."""
    return float(brentq(lambda T: moist_air_enthalpy(T, omega) - h, 200.0, 1000.0, xtol=1e-12))


def mix_streams(m1: float, T1: float, w1: float, m2: float, T2: float, w2: float) -> tuple[float, float, float]:
    """Adiabatic, steady mixing of two moist-air streams (dry-air mass flows m1, m2 [kg/s]) with no work:
    dry air  m3 = m1 + m2;  water  m3 w3 = m1 w1 + m2 w2;  energy  m3 h3 = m1 h1 + m2 h2.
    Returns (m3, T3, w3). Both w3 and h3 are mass-weighted averages, so the mixed state lies on the
    straight line joining the two inlet states in the (omega, h) plane."""
    m3 = m1 + m2
    w3 = (m1 * w1 + m2 * w2) / m3
    h3 = (m1 * moist_air_enthalpy(T1, w1) + m2 * moist_air_enthalpy(T2, w2)) / m3
    return m3, temperature_from_enthalpy(h3, w3), w3


def cool_dehumidify(m_a: float, T_in: float, w_in: float, T_out: float, p: float) -> dict[str, float]:
    """Steady cooling of moist air at constant p to T_out (a pack heat exchanger with a water separator).
    If T_out is above the dew point the humidity ratio is unchanged; otherwise the air leaves saturated,
    w_out = w_sat(T_out, p), and the excess water leaves as liquid at T_out.
    water: m_a w_in = m_a w_out + m_cond;  energy: Q = m_a (h_in - h_out) - m_cond h_f(T_out).
    Returns {'w_out', 'm_cond' [kg/s], 'Q' [W, heat removed, positive], 'phi_out'}."""
    w_sat = omega_saturated(T_out, p)
    w_out = min(w_in, w_sat)
    m_cond = m_a * (w_in - w_out)
    h_cond = liquid_enthalpy(T_out) if m_cond > 0.0 else 0.0
    Q = m_a * (moist_air_enthalpy(T_in, w_in) - moist_air_enthalpy(T_out, w_out)) - m_cond * h_cond
    return {"w_out": w_out, "m_cond": m_cond, "Q": Q, "phi_out": relative_humidity(T_out, w_out, p)}


def cabin_humidity(m_fresh: float, m_recirc: float, m_water: float, w_out: float) -> tuple[float, float]:
    """Steady water balance of the cabin and the mixing manifold. Fresh air (dry-air flow m_fresh,
    humidity w_out) mixes with recirculated cabin air (m_recirc, w_cabin); the occupants add m_water
    [kg/s]; the outflow valve exhausts m_fresh of cabin air. Returns (w_cabin, w_supply).
    Cabin + manifold:  m_fresh w_out + m_water = m_fresh w_cabin   =>  w_cabin = w_out + m_water / m_fresh,
    independent of m_recirc. Manifold:  w_supply = (m_fresh w_out + m_recirc w_cabin) / (m_fresh + m_recirc)."""
    w_c = w_out + m_water / m_fresh
    w_s = (m_fresh * w_out + m_recirc * w_c) / (m_fresh + m_recirc)
    return w_c, w_s


def ecs_cruise(n_occ: int = 186, m_fresh_per_occ: float = 0.25 / 60, recirc_fraction: float = 0.5,
               water_per_occ: float = 0.060 / 3600, Q_sens: float = 15.0e3, T_cabin: float = 295.15,
               z_cabin: float = 2438.4, z_flight: float = 11000.0, w_out: float = 5.0e-5,
               T_bleed: float = 473.15) -> dict[str, float]:
    """The cruise environmental control system, steady state, with every balance written out.

    Outside air at the ISA state at z_flight (humidity ratio w_out) is taken as engine bleed at T_bleed,
    cooled by the pack (no water to condense) to T_pack, mixed in the manifold with recirculated cabin
    air (recirc_fraction of the supply), and delivered to the cabin at T_supply; the cabin, at T_cabin
    and the ISA pressure of the cabin altitude z_cabin, gains the sensible load Q_sens [W] and the
    occupants' water m_water as vapour at T_cabin; the outflow valve exhausts the fresh-air flow.

    Unknowns from the balances: w_cabin, w_supply (water), h_supply from the cabin energy balance,
    T_pack from the manifold energy balance, the pack heat duty Q_pack = m_f (h(T_bleed) - h(T_pack)).
    Returns a dict of states, flows and residuals (dry-air, water and energy; all should be ~0)."""
    T_o, p_o = isa_state(z_flight)
    _, p_c = isa_state(z_cabin)
    m_f = n_occ * m_fresh_per_occ
    m_s = m_f / (1.0 - recirc_fraction)
    m_r = m_s - m_f
    m_w = n_occ * water_per_occ
    w_c, w_s = cabin_humidity(m_f, m_r, m_w, w_out)
    h_c = moist_air_enthalpy(T_cabin, w_c)
    # cabin: m_s h_s + Q_sens + m_w h_v(T_c) = m_s h_c
    h_s = h_c - (Q_sens + m_w * vapour_enthalpy(T_cabin)) / m_s
    T_s = temperature_from_enthalpy(h_s, w_s)
    # manifold: m_f h(T_pack, w_out) + m_r h_c = m_s h_s
    h_p = (m_s * h_s - m_r * h_c) / m_f
    T_p = temperature_from_enthalpy(h_p, w_out)
    Q_pack = m_f * (moist_air_enthalpy(T_bleed, w_out) - h_p)
    T_dp, est = dew_point(w_c, p_c)
    res_air = (m_f + m_r) - m_s
    res_water = m_f * w_out + m_w - m_f * w_c
    res_energy = m_f * h_p + Q_sens + m_w * vapour_enthalpy(T_cabin) - m_f * h_c
    return {"T_out": T_o, "p_out": p_o, "p_cabin": p_c, "m_fresh": m_f, "m_recirc": m_r, "m_supply": m_s,
            "m_water": m_w, "w_out": w_out, "w_cabin": w_c, "w_supply": w_s,
            "phi_cabin": relative_humidity(T_cabin, w_c, p_c),
            "phi_fresh_at_cabin": relative_humidity(T_cabin, w_out, p_c),
            "pv_out": vapour_pressure(w_out, p_o),
            "T_dew_cabin": T_dp, "dew_is_estimate": est, "T_supply": T_s, "T_pack": T_p, "Q_pack": Q_pack,
            "res_air": res_air, "res_water": res_water, "res_energy": res_energy}


# ---------------------------------------------------------------- the lesson's numbers

def ground_case(T_amb: float = 308.15, phi_amb: float = 0.60, T_pack: float = 275.15,
                m_fresh: float = 186 * 0.25 / 60) -> dict[str, float]:
    """On the ground on a hot, humid day: outside air at T_amb, phi_amb and 101 325 Pa through a pack
    that cools it to T_pack; returns cool_dehumidify's result plus the inlet humidity and dew point."""
    w_in = humidity_ratio(phi_amb * p_sat(T_amb), P_ATM)
    out = cool_dehumidify(m_fresh, T_amb, w_in, T_pack, P_ATM)
    out.update({"w_in": w_in, "T_dew_in": dew_point(w_in, P_ATM)[0], "m_fresh": m_fresh})
    return out


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save psychrometrics.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    M_air = molar_mass(AIR_X)
    y = mass_fractions(AIR_X)
    print(f"dry air: M = {M_air * 1e3:.4f} g/mol, eps = M_w/M_a = {EPS:.5f}; mass fractions "
          + ", ".join(f"{k} {v:.4f}" for k, v in y.items()))
    n = {k: v for k, v in AIR_X.items()}          # 1 mol of air
    T, V = 300.0, R_U * 300.0 / P_ATM
    pp = partial_pressures(n, T, V)
    pv = partial_volumes(n, T, P_ATM)
    print(f"1 mol of air at 300 K, 1 atm: sum of Dalton partial pressures {sum(pp.values()):.3f} Pa, "
          f"sum of Amagat partial volumes {sum(pv.values()) * 1e3:.5f} L (V = {V * 1e3:.5f} L)")
    print(f"entropy of mixing of air from pure N2, O2, Ar, CO2: {entropy_of_mixing(AIR_X):.4f} J/(mol K) "
          f"= {entropy_of_mixing(AIR_X) / M_air:.2f} J/(kg K); R ln 2 = {R_U * np.log(2):.4f}")
    print("\nsaturation pressure: table (NIST SRD 69 via Yan) against the NASA Gibbs-equality estimate")
    print("  T (C)   table (Pa)   NASA (Pa)   diff (%)")
    worst = 0.0
    for row in SAT_WATER:
        Tk = max(row[0] + 273.15, T_TRIPLE)
        a, b = p_sat(Tk), p_sat_gibbs(Tk)
        worst = max(worst, abs(b / a - 1))
        if row[0] in (0.01, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0):
            print(f"  {row[0]:5.2f}   {a:9.1f}   {b:9.1f}   {100 * (b / a - 1):+6.2f}")
    print(f"  largest difference 0.01-60 C: {100 * worst:.2f} %")
    print(f"  h_fg at 0.01 C: table {latent_heat_triple() / 1e3:.2f} kJ/kg, NASA "
          f"{(h_species('H2O', T_TRIPLE) - h_species('H2O(L)', T_TRIPLE)) / M_W / 1e3:.2f} kJ/kg")

    e = ecs_cruise()
    print(f"\nECS at cruise (typical, illustrative inputs): outside {e['T_out']:.2f} K, {e['p_out']:.0f} Pa, "
          f"w = {e['w_out']:.1e} (p_v = {e['pv_out']:.2f} Pa); cabin 22.0 C at {e['p_cabin']:.0f} Pa (8000 ft)")
    print(f"  fresh {e['m_fresh']:.3f} kg/s, recirculated {e['m_recirc']:.3f} kg/s, supply {e['m_supply']:.3f} kg/s, "
          f"occupant water {e['m_water'] * 3600:.2f} kg/h")
    print(f"  outside air brought to the cabin unchanged: RH {100 * e['phi_fresh_at_cabin']:.2f} %")
    print(f"  cabin w = {e['w_cabin'] * 1e3:.3f} g/kg, RH {100 * e['phi_cabin']:.1f} %, dew point "
          f"{e['T_dew_cabin'] - 273.15:.1f} C ({'Clausius-Clapeyron estimate below 0.01 C' if e['dew_is_estimate'] else 'table'})")
    print(f"  supply w = {e['w_supply'] * 1e3:.3f} g/kg at {e['T_supply'] - 273.15:.2f} C; pack outlet "
          f"{e['T_pack'] - 273.15:.2f} C; pack duty {e['Q_pack'] / 1e3:.1f} kW from {473.15 - 273.15:.0f} C bleed")
    print(f"  residuals: dry air {e['res_air']:.1e} kg/s, water {e['res_water']:.1e} kg/s, energy {e['res_energy']:.1e} W")

    print("\nrecirculation sweep (same fresh air): fraction, cabin RH (%), supply T (C), pack outlet T (C)")
    for r in (0.0, 0.3, 0.5, 0.6):
        er = ecs_cruise(recirc_fraction=r)
        print(f"  {r:.1f}   {100 * er['phi_cabin']:.2f}   {er['T_supply'] - 273.15:6.2f}   {er['T_pack'] - 273.15:6.2f}")
    print("fresh-air sweep (recirc 0.5): kg/min per occupant, cabin RH (%), dew point (C)")
    for mf in (0.15, 0.25, 0.35, 0.5):
        er = ecs_cruise(m_fresh_per_occ=mf / 60)
        print(f"  {mf:.2f}   {100 * er['phi_cabin']:.2f}   {er['T_dew_cabin'] - 273.15:6.2f}")
    print("outside-humidity sweep: w_out, cabin RH (%)")
    for wo in (0.0, 1e-5, 5e-5, 1e-4, 3e-4):
        print(f"  {wo:.0e}   {100 * ecs_cruise(w_out=wo)['phi_cabin']:.2f}")

    g = ground_case()
    print(f"\nground, 35 C and 60 % RH: w = {g['w_in'] * 1e3:.2f} g/kg, dew point {g['T_dew_in'] - 273.15:.1f} C; "
          f"pack to 2 C: w_out = {g['w_out'] * 1e3:.3f} g/kg, condensate {g['m_cond'] * 3600:.1f} kg/h, "
          f"heat removed {g['Q'] / 1e3:.1f} kW")
    s1, s2 = (g['m_fresh'], 275.15, g['w_out']), (g['m_fresh'], 297.15, 8e-3)
    m3, T3, w3 = mix_streams(*s1, *s2)
    print(f"mixing equal flows (2 C, {g['w_out'] * 1e3:.2f} g/kg) and (24 C, 8.00 g/kg): {T3 - 273.15:.3f} C, "
          f"{w3 * 1e3:.3f} g/kg (arithmetic mean of T would be 13.000 C)")

    # psychrometric chart at the cabin pressure, with the cruise states
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
    p_c = e["p_cabin"]
    Ts = np.linspace(258.15, 313.15, 200)
    for phi in (0.1, 0.2, 0.4, 0.6, 0.8, 1.0):
        w = [humidity_ratio(phi * p_sat_any(t), p_c) * 1e3 for t in Ts]
        ax[0].plot(Ts - 273.15, w, "k-" if phi == 1.0 else "k:", lw=1.2 if phi == 1.0 else 0.7)
        k = int(np.searchsorted(w, 24.0)) - 1 if w[-1] > 24.0 else len(w) - 1
        ax[0].text(Ts[k] - 273.15, w[k], f" {int(phi * 100)} %", fontsize=7, va="bottom")
    ax[0].axvspan(Ts[0] - 273.15, 0.01, color="0.9", label="below 0.01 C: CC estimate")
    pts = {"cabin": (295.15, e["w_cabin"]), "supply": (e["T_supply"], e["w_supply"]), "pack out": (e["T_pack"], e["w_out"])}
    for k, (t, w) in pts.items():
        ax[0].plot(t - 273.15, w * 1e3, "o")
        ax[0].annotate(k, (t - 273.15, w * 1e3), textcoords="offset points", xytext=(4, 4), fontsize=8)
    ax[0].plot([e["T_pack"] - 273.15, 22.0], [e["w_out"] * 1e3, e["w_cabin"] * 1e3], "C0--", lw=0.8)
    ax[0].set_xlabel("dry-bulb temperature (C)")
    ax[0].set_ylabel("humidity ratio (g water / kg dry air)")
    ax[0].set_ylim(0, 25)
    ax[0].set_title(f"cruise ECS states, p = {p_c / 1e3:.1f} kPa")
    ax[0].legend(fontsize=7, loc="upper left")
    Tt = np.linspace(T_TRIPLE, 333.15, 120)
    ax[1].plot(Tt - 273.15, [100 * (p_sat_gibbs(t) / p_sat(t) - 1) for t in Tt], label="NASA Gibbs estimate vs table")
    ax[1].axhline(0.0, color="k", lw=0.6)
    ax[1].set_xlabel("T (C)")
    ax[1].set_ylabel("difference in p_sat (%)")
    ax[1].set_title("independent check of the saturation pressure")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "psychrometrics.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'psychrometrics.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
