"""P07.11 Chemical Equilibrium: starter code.

At what temperature and composition do combustion gases actually settle? Lesson P07.10 burned every
fuel completely to CO2 and H2O. Above about 2000 K that is wrong: the products dissociate into CO,
H2, OH, H, O and NO, which absorbs energy and lowers the flame temperature. This module finds the
true (equilibrium) composition by minimising the Gibbs energy at fixed T and p, subject to the
conservation of atoms:

* the element-potential method of NASA RP-1311 (Gordon & McBride; the algorithm behind NASA CEA):
  Newton's method on the Lagrange conditions mu_j / RT = sum_i a_ij pi_i, in the logarithms of the
  amounts, reduced to a small linear system in the element potentials pi_i;
* an independent check by a general-purpose constrained minimiser (scipy.optimize.minimize, SLSQP),
  which must reach the same Gibbs energy;
* the law of mass action: for every independent reaction, the equilibrium constant from the standard
  Gibbs energy change, K_p = exp(-Delta G° / RT), must match the solver's composition;
* the equilibrium adiabatic flame temperature (energy balance with the equilibrium products at each
  trial T), compared with the complete-combustion value of P07.10;
* applications: H2/O2 and C/H/O/N (Jet-A/air, CH4/air) flames and an illustrative rocket chamber.

COMPUTATIONAL ONLY: there is no combustion, flame, burning, fuel handling or fuel of any kind in this
lesson. Every result is computed.

Ideal-gas mixture; gas-phase species only (no soot, no liquid); NASA polynomials to 6000 K.
Units: SI (K, Pa, mol, J). Amounts are in mol; element amounts b are mol of atoms.
Standard-state pressure p° = 1 bar (the NASA data).

Run:  python main.py [out_dir]   (prints the lesson's numbers, writes equilibrium.png)
Test: pytest curriculum/phase-07/lesson-11-chemical-equilibrium
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize

R_U = 8.31446261815324          # J/(mol K), CODATA 2018 (exact)
P_REF = 1.0e5                   # Pa, standard-state pressure of the NASA data
T_REF = 298.15                  # K

# Dry air model fixed for Phase 07 (mole fractions; M = 28.965 g/mol)
AIR_X = {"N2": 0.7808, "O2": 0.2095, "Ar": 0.0093, "CO2": 0.0004}

# NASA Glenn 9-coefficient polynomials, copied unrounded from data/thermo.inp of NASA CEA
# (https://github.com/nasa/cea, Apache-2.0; McBride, Zehe & Gordon, NASA/TP-2002-211556).
# Intervals 200-1000 K and 1000-6000 K. Jet-A(g) is the CEA surrogate C12H23 as a vapour.
# cp/R = a1 T^-2 + a2 T^-1 + a3 + a4 T + a5 T^2 + a6 T^3 + a7 T^4; h includes Hf298; s at p = 1 bar.
NASA9: dict[str, dict] = {
    'N2': {
        "M": 28.0134,  # g/mol
        "Hf298": 0.0,  # J/mol
        "elements": {'N': 2},
        "intervals": [
            (200.0, 1000.0, [22103.71497, -381.846182, 6.08273836, -0.00853091441, 1.384646189e-05, -9.62579362e-09, 2.519705809e-12], [710.846086, -10.76003744]),
            (1000.0, 6000.0, [587712.406, -2239.249073, 6.06694922, -0.00061396855, 1.491806679e-07, -1.923105485e-11, 1.061954386e-15], [12832.10415, -15.86640027]),
        ],
    },
    'O2': {
        "M": 31.9988,  # g/mol
        "Hf298": 0.0,  # J/mol
        "elements": {'O': 2},
        "intervals": [
            (200.0, 1000.0, [-34255.6342, 484.700097, 1.119010961, 0.00429388924, -6.83630052e-07, -2.0233727e-09, 1.039040018e-12], [-3391.45487, 18.4969947]),
            (1000.0, 6000.0, [-1037939.022, 2344.830282, 1.819732036, 0.001267847582, -2.188067988e-07, 2.053719572e-11, -8.19346705e-16], [-16890.10929, 17.38716506]),
        ],
    },
    'Ar': {
        "M": 39.948,  # g/mol
        "Hf298": 0.0,  # J/mol
        "elements": {'Ar': 1},
        "intervals": [
            (200.0, 1000.0, [0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0], [-745.375, 4.37967491]),
            (1000.0, 6000.0, [20.10538475, -0.0599266107, 2.500069401, -3.99214116e-08, 1.20527214e-11, -1.819015576e-15, 1.078576636e-19], [-744.993961, 4.37918011]),
        ],
    },
    'CO2': {
        "M": 44.0095,  # g/mol
        "Hf298": -393510.0,  # J/mol
        "elements": {'C': 1, 'O': 2},
        "intervals": [
            (200.0, 1000.0, [49436.5054, -626.411601, 5.30172524, 0.002503813816, -2.127308728e-07, -7.68998878e-10, 2.849677801e-13], [-45281.9846, -7.04827944]),
            (1000.0, 6000.0, [117696.2419, -1788.791477, 8.29152319, -9.22315678e-05, 4.86367688e-09, -1.891053312e-12, 6.33003659e-16], [-39083.5059, -26.52669281]),
        ],
    },
    'H2O': {
        "M": 18.01528,  # g/mol
        "Hf298": -241826.0,  # J/mol
        "elements": {'H': 2, 'O': 1},
        "intervals": [
            (200.0, 1000.0, [-39479.6083, 575.573102, 0.931782653, 0.00722271286, -7.34255737e-06, 4.95504349e-09, -1.336933246e-12], [-33039.7431, 17.24205775]),
            (1000.0, 6000.0, [1034972.096, -2412.698562, 4.64611078, 0.002291998307, -6.83683048e-07, 9.42646893e-11, -4.82238053e-15], [-13842.86509, -7.97814851]),
        ],
    },
    'CO': {
        "M": 28.0101,  # g/mol
        "Hf298": -110535.196,  # J/mol
        "elements": {'C': 1, 'O': 1},
        "intervals": [
            (200.0, 1000.0, [14890.45326, -292.2285939, 5.72452717, -0.00817623503, 1.456903469e-05, -1.087746302e-08, 3.027941827e-12], [-13031.31878, -7.85924135]),
            (1000.0, 6000.0, [461919.725, -1944.704863, 5.91671418, -0.000566428283, 1.39881454e-07, -1.787680361e-11, 9.62093557e-16], [-2466.261084, -13.87413108]),
        ],
    },
    'H2': {
        "M": 2.01588,  # g/mol
        "Hf298": 0.0,  # J/mol
        "elements": {'H': 2},
        "intervals": [
            (200.0, 1000.0, [40783.2321, -800.918604, 8.21470201, -0.01269714457, 1.753605076e-05, -1.20286027e-08, 3.36809349e-12], [2682.484665, -30.43788844]),
            (1000.0, 6000.0, [560812.801, -837.150474, 2.975364532, 0.001252249124, -3.74071619e-07, 5.9366252e-11, -3.6069941e-15], [5339.82441, -2.202774769]),
        ],
    },
    'OH': {
        "M": 17.00734,  # g/mol
        "Hf298": 37278.206,  # J/mol
        "elements": {'O': 1, 'H': 1},
        "intervals": [
            (200.0, 1000.0, [-1998.85899, 93.0013616, 3.050854229, 0.001529529288, -3.157890998e-06, 3.31544618e-09, -1.138762683e-12], [2991.214235, 4.67411079]),
            (1000.0, 6000.0, [1017393.379, -2509.957276, 5.11654786, 0.000130529993, -8.28432226e-08, 2.006475941e-11, -1.556993656e-15], [20196.40206, -11.01282337]),
        ],
    },
    'H': {
        "M": 1.00794,  # g/mol
        "Hf298": 217998.828,  # J/mol
        "elements": {'H': 1},
        "intervals": [
            (200.0, 1000.0, [0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0], [25473.70801, -0.446682853]),
            (1000.0, 6000.0, [60.7877425, -0.1819354417, 2.500211817, -1.226512864e-07, 3.73287633e-11, -5.68774456e-15, 3.410210197e-19], [25474.86398, -0.448191777]),
        ],
    },
    'O': {
        "M": 15.9994,  # g/mol
        "Hf298": 249175.003,  # J/mol
        "elements": {'O': 1},
        "intervals": [
            (200.0, 1000.0, [-7953.6113, 160.7177787, 1.966226438, 0.00101367031, -1.110415423e-06, 6.5175075e-10, -1.584779251e-13], [28403.62437, 8.40424182]),
            (1000.0, 6000.0, [261902.0262, -729.872203, 3.31717727, -0.000428133436, 1.036104594e-07, -9.43830433e-12, 2.725038297e-16], [33924.2806, -0.667958535]),
        ],
    },
    'NO': {
        "M": 30.0061,  # g/mol
        "Hf298": 91271.31,  # J/mol
        "elements": {'N': 1, 'O': 1},
        "intervals": [
            (200.0, 1000.0, [-11439.16503, 153.6467592, 3.43146873, -0.002668592368, 8.48139912e-06, -7.68511105e-09, 2.386797655e-12], [9098.21441, 6.72872549]),
            (1000.0, 6000.0, [223901.8716, -1289.651623, 5.43393603, -0.00036560349, 9.88096645e-08, -1.416076856e-11, 9.38018462e-16], [17503.17656, -8.50166909]),
        ],
    },
    'N': {
        "M": 14.0067,  # g/mol
        "Hf298": 472680.0,  # J/mol
        "elements": {'N': 1},
        "intervals": [
            (200.0, 1000.0, [0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0], [56104.6378, 4.193905036]),
            (1000.0, 6000.0, [88765.0138, -107.12315, 2.362188287, 0.0002916720081, -1.7295151e-07, 4.01265788e-11, -2.677227571e-15], [56973.5133, 4.865231506]),
        ],
    },
    'CH4': {
        "M": 16.04246,  # g/mol
        "Hf298": -74600.0,  # J/mol
        "elements": {'C': 1, 'H': 4},
        "intervals": [
            (200.0, 1000.0, [-176685.0998, 2786.18102, -12.0257785, 0.0391761929, -3.61905443e-05, 2.026853043e-08, -4.97670549e-12], [-23313.1436, 89.0432275]),
            (1000.0, 6000.0, [3730042.76, -13835.01485, 20.49107091, -0.001961974759, 4.72731304e-07, -3.72881469e-11, 1.623737207e-15], [75320.6691, -121.9124889]),
        ],
    },
    'Jet-A(g)': {
        "M": 167.31102,  # g/mol
        "Hf298": -249657.0,  # J/mol
        "elements": {'C': 12, 'H': 23},
        "intervals": [
            (273.15, 1000.0, [-606869.559, 8328.25959, -43.1232127, 0.2572390455, -0.000262931604, 1.64498894e-07, -4.64533514e-11], [-76069.6276, 279.4305937]),
            (1000.0, 6000.0, [18583561.02, -76772.1989, 141.9826133, -0.00743752453, 5.85620255e-07, 1.223955647e-11, -3.149201922e-15], [422198.952, -898.606104]),
        ],
    },
}

ELEMENTS = ("C", "H", "O", "N", "Ar")
PRODUCT_CANDIDATES = ("CO2", "H2O", "CO", "H2", "O2", "N2", "Ar", "OH", "H", "O", "NO", "N", "CH4")


# ---------------------------------------------------------------- given: NASA polynomials

def _coeffs(species: str, T: float) -> tuple[list[float], list[float]]:
    """The (a1..a7, b1 b2) coefficients of the interval holding T. Raises ValueError outside the fits."""
    for lo, hi, a, b in NASA9[species]["intervals"]:
        if lo <= T <= hi:
            return a, b
    raise ValueError(f"T = {T} K outside the fitted range for {species}")


def h_species(species: str, T: float) -> float:
    """Molar enthalpy h(T) [J/mol] measured from the elements: h(298.15 K) = Hf298."""
    a, b = _coeffs(species, T)
    return R_U * T * (-a[0] / T ** 2 + a[1] * np.log(T) / T + a[2] + a[3] * T / 2 + a[4] * T ** 2 / 3
                      + a[5] * T ** 3 / 4 + a[6] * T ** 4 / 5 + b[0] / T)


def s_species(species: str, T: float) -> float:
    """Molar standard-state entropy s°(T) [J/(mol K)] at p° = 1 bar."""
    a, b = _coeffs(species, T)
    return R_U * (-a[0] / T ** 2 / 2 - a[1] / T + a[2] * np.log(T) + a[3] * T + a[4] * T ** 2 / 2
                  + a[5] * T ** 3 / 3 + a[6] * T ** 4 / 4 + b[1])


def g_species(species: str, T: float) -> float:
    """Molar standard-state Gibbs energy g°(T) = h - T s° [J/mol] at 1 bar."""
    return h_species(species, T) - T * s_species(species, T)


def molar_mass(species: str) -> float:
    """Molar mass [kg/mol]."""
    return NASA9[species]["M"] / 1000.0


def element_totals(n: dict[str, float]) -> dict[str, float]:
    """Moles of each element in the amounts n [mol]."""
    out: dict[str, float] = {}
    for s, nj in n.items():
        for el, a in NASA9[s]["elements"].items():
            out[el] = out.get(el, 0.0) + a * nj
    return out


def element_matrix(species: list[str], elements: list[str]) -> np.ndarray:
    """A[i, j] = atoms of element i in one molecule of species j."""
    return np.array([[NASA9[s]["elements"].get(e, 0) for s in species] for e in elements], dtype=float)


def species_for(b0: dict[str, float], candidates: tuple[str, ...] = PRODUCT_CANDIDATES) -> tuple[list[str], list[str]]:
    """The elements present in b0 (in ELEMENTS order) and the candidate species made only of them."""
    els = [e for e in ELEMENTS if b0.get(e, 0.0) > 0.0]
    sp = [s for s in candidates if set(NASA9[s]["elements"]) <= set(els)]
    return sp, els


# ---------------------------------------------------------------- given: complete combustion (from P07.10)

def o2_stoich(fuel: str) -> float:
    """Moles of O2 to burn one mole of fuel CxHyOz to CO2 and H2O: x + y/4 - z/2."""
    el = NASA9[fuel]["elements"]
    return el.get("C", 0) + el.get("H", 0) / 4.0 - el.get("O", 0) / 2.0


def reactants(fuel: str, phi: float, oxidizer: str = "air") -> dict[str, float]:
    """Reactant amounts per mole of fuel at equivalence ratio phi, with dry air or pure O2."""
    n_o2 = o2_stoich(fuel) / phi
    out = {fuel: 1.0}
    if oxidizer == "O2":
        out["O2"] = n_o2
    else:
        for k, x in AIR_X.items():
            out[k] = out.get(k, 0.0) + n_o2 * x / AIR_X["O2"]
    return out


def reactant_enthalpy(fuel: str, phi: float, oxidizer: str = "air", T_fuel: float = T_REF,
                      T_ox: float = T_REF) -> float:
    """Enthalpy [J per mole of fuel] of the reactants, fuel at T_fuel and oxidizer at T_ox."""
    n = reactants(fuel, phi, oxidizer)
    return h_species(fuel, T_fuel) + sum(v * h_species(k, T_ox) for k, v in n.items() if k != fuel)


def complete_products(fuel: str, phi: float, oxidizer: str = "air") -> dict[str, float]:
    """Products of complete combustion (P07.10): C to CO first, then H to H2O, then CO to CO2, spare O2."""
    el = element_totals(reactants(fuel, phi, oxidizer))
    C, H, O = el.get("C", 0.0), el.get("H", 0.0), el.get("O", 0.0)
    o_left = O - C
    h2o = min(H / 2.0, o_left)
    o_left -= h2o
    co2 = min(C, o_left)
    out = {"CO2": co2, "H2O": h2o, "CO": C - co2, "H2": H / 2.0 - h2o, "O2": max((o_left - co2) / 2.0, 0.0),
           "N2": el.get("N", 0.0) / 2.0, "Ar": el.get("Ar", 0.0)}
    return {k: v for k, v in out.items() if v > 0.0}


def complete_flame_temperature(fuel: str, phi: float, oxidizer: str = "air", T_fuel: float = T_REF,
                               T_ox: float = T_REF) -> float:
    """Complete-combustion adiabatic flame temperature (P07.10) [K]."""
    H_r = reactant_enthalpy(fuel, phi, oxidizer, T_fuel, T_ox)
    prod = complete_products(fuel, phi, oxidizer)
    return float(brentq(lambda T: sum(v * h_species(k, T) for k, v in prod.items()) - H_r, 200.0, 6000.0, xtol=1e-10))


# ---------------------------------------------------------------- Gibbs energy

def gibbs_rt(n: np.ndarray, species: list[str], T: float, p: float) -> float:
    """Dimensionless Gibbs energy G/(R T) of an ideal-gas mixture of amounts n [mol] at T, p:
    G/RT = sum_j n_j [g°_j(T)/RT + ln(n_j / n_tot) + ln(p / p°)]  (terms with n_j = 0 contribute 0)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the element-potential method

def equilibrium_tp(b0: dict[str, float], T: float, p: float, species: list[str] | None = None,
                   ln_n0: np.ndarray | None = None, tol: float = 1e-12, max_iter: int = 500) -> dict:
    """Equilibrium composition at fixed T [K] and p [Pa] for element amounts b0 {element: mol of atoms},
    by minimising G subject to A n = b0 with the element-potential method of NASA RP-1311.

    Unknowns: Delta ln n_j (each species), pi_i (element potentials, the Lagrange multipliers / RT) and
    Delta ln n (the total amount n is carried as its own variable). Newton's method on
        Delta ln n_j - sum_i a_ij pi_i - Delta ln n = -mu_j/RT                    (one per species)
        sum_j a_kj n_j Delta ln n_j = b0_k - b_k                                 (one per element)
        sum_j n_j Delta ln n_j - n Delta ln n = n - sum_j n_j                     (total amount)
    Substituting the first into the others leaves a (n_elements + 1) square system:
        sum_i [sum_j a_kj a_ij n_j] pi_i + b_k Delta ln n = b0_k - b_k + sum_j a_kj n_j mu_j/RT
        sum_i b_i pi_i + (sum_j n_j - n) Delta ln n     = n - sum_j n_j + sum_j n_j mu_j/RT
    then Delta ln n_j = -mu_j/RT + sum_i a_ij pi_i + Delta ln n. Steps are damped as in RP-1311:
    lambda = min(1, 2 / max(5 |Delta ln n|, |Delta ln n_j| of major species), and a limit that stops a
    trace species (n_j/n < 1e-8) from jumping above 1e-4 in one step). Converged, as in RP-1311 but
    tighter, when every (n_j/n)|Delta ln n_j| and |Delta ln n| is below tol and every |Delta ln n_j|
    (trace species included) is below 1e-6. (The matrix becomes ill-conditioned when one species holds
    nearly all of two elements, e.g. H2O at low T, so Delta ln n_j of trace species stalls near 1e-8.)

    Returns {'species', 'elements', 'n' (mol), 'ln_n', 'x' (mole fractions), 'pi', 'iterations',
    'G_RT', 'element_residual' (max |A n - b0| / max b0)}."""
    raise NotImplementedError


def equilibrium_scipy(b0: dict[str, float], T: float, p: float, species: list[str] | None = None) -> dict:
    """The independent check: minimise G/RT directly with scipy.optimize.minimize (SLSQP), subject to
    A n = b0 and n_j >= 1e-30, with the analytic gradient dG/dn_j / RT = mu_j/RT. Returns
    {'species', 'n', 'x', 'G_RT', 'element_residual', 'success'}."""
    raise NotImplementedError


# ---------------------------------------------------------------- law of mass action

def equilibrium_constant(reaction: dict[str, float], T: float) -> float:
    """K_p = exp(-Delta G° / RT) for a reaction {species: stoichiometric coefficient}, products
    positive and reactants negative; Delta G° = sum nu_j g°_j(T) at p° = 1 bar."""
    raise NotImplementedError


def formation_reactions(species: list[str], base: list[str]) -> list[dict[str, float]]:
    """Independent reactions: one per non-base species j, forming j from the base species,
    sum_b nu_b B -> j, with nu solved from the element balance A_base nu = a_j. Returned as
    {j: +1, B: -nu_B, ...}. The base must contain one species per element with an invertible A_base."""
    raise NotImplementedError


def mass_action_residuals(result: dict, T: float, p: float, base: list[str]) -> list[tuple[dict, float, float]]:
    """For each formation reaction, (reaction, ln K_p from Delta G°, ln Q_p from the composition), with
    ln Q_p = sum_j nu_j ln(x_j p / p°). At equilibrium ln Q_p = ln K_p for every reaction."""
    raise NotImplementedError


# ---------------------------------------------------------------- equilibrium flame temperature

def equilibrium_flame_temperature(fuel: str, phi: float, oxidizer: str = "air", p: float = 101325.0,
                                  T_fuel: float = T_REF, T_ox: float = T_REF) -> tuple[float, dict]:
    """Adiabatic flame temperature with equilibrium products [K]: the T at which the enthalpy of the
    equilibrium mixture at (T, p) equals the reactants' enthalpy. Each trial T needs an equilibrium
    solve (warm-started from the last one); brentq on T between 1000 K and the complete-combustion
    value + 300 K. (For lean and stoichiometric mixtures the complete-combustion value is an upper
    bound, because dissociation only absorbs energy; for rich ones the fixed product order of P07.10 is
    not the equilibrium split between CO, CO2, H2 and H2O, and equilibrium can come out slightly hotter.)
    Returns (T_ad, equilibrium result at T_ad)."""
    raise NotImplementedError


def oxidizer_fuel_to_phi(fuel: str, of_mass: float, oxidizer: str = "O2") -> float:
    """Equivalence ratio from a rocket's mixture ratio O/F (mass of oxidizer per mass of fuel)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the lesson's numbers

def _fmt(r: dict, cut: float = 1e-5) -> str:
    return ", ".join(f"{s} {x:.4f}" if x >= 1e-3 else f"{s} {x:.2e}" for s, x in zip(r["species"], r["x"]) if x >= cut)


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save equilibrium.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    print("COMPUTATIONAL ONLY: no combustion, flame or fuel of any kind; every number below is computed.\n")

    b_ho = element_totals(reactants("H2", 1.0, "O2"))
    r = equilibrium_tp(b_ho, 3000.0, 1e5)
    s = equilibrium_scipy(b_ho, 3000.0, 1e5)
    print(f"H2 + O2/2 at 3000 K, 1 bar: {r['iterations']} Newton iterations, element residual {r['element_residual']:.1e}")
    print(f"  mole fractions: {_fmt(r)}")
    print(f"  element potentials pi_H = {r['pi'][0]:.4f}, pi_O = {r['pi'][1]:.4f}")
    print(f"  G/RT: element potentials {r['G_RT']:.12f}, SLSQP {s['G_RT']:.12f} (difference {s['G_RT'] - r['G_RT']:.1e})")
    print("  law of mass action (base species H2, O2):")
    for rx, lnK, lnQ in mass_action_residuals(r, 3000.0, 1e5, ["H2", "O2"]):
        eq = " + ".join(f"{-v:g} {k}" for k, v in rx.items() if v < 0) + " -> " + " + ".join(f"{v:g} {k}" for k, v in rx.items() if v > 0)
        print(f"    {eq:28s} K_p = {np.exp(lnK):.4e}, ln K = {lnK:+.6f}, ln Q = {lnQ:+.6f}")
    print(f"  complete combustion would be pure H2O; equilibrium x_H2O = {r['x'][r['species'].index('H2O')]:.4f}")

    print("\nH2/O2 (phi = 1) mole fraction of H2O against T and p:")
    for T in (1000.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0):
        row = [equilibrium_tp(b_ho, T, pp) for pp in (1e4, 1e5, 1e6, 1e7)]
        print(f"  {T:6.0f} K: " + "  ".join(f"{q['x'][q['species'].index('H2O')]:.6f}" for q in row) + "   (0.1, 1, 10, 100 bar)")
    r1 = equilibrium_tp(b_ho, 1000.0, 1e5)
    print(f"  low-T limit, 1000 K and 1 bar: 1 - x_H2O = {1 - r1['x'][r1['species'].index('H2O')]:.2e}")
    b_air = element_totals(reactants("CH4", 0.8, "air"))
    r2 = equilibrium_tp(b_air, 1000.0, 101325.0)
    comp = complete_products("CH4", 0.8, "air")
    ntot = sum(comp.values())
    dev = max(abs(r2["x"][r2["species"].index(k)] - v / ntot) for k, v in comp.items())
    print(f"  CH4/air phi = 0.8 at 1000 K: largest |x_eq - x_complete| = {dev:.2e} (NO {r2['x'][r2['species'].index('NO')]:.2e})")

    print("\nadiabatic flame temperatures: complete combustion (P07.10) against equilibrium")
    cases = [("CH4", "air", 1.0, 101325.0, 298.15), ("Jet-A(g)", "air", 1.0, 101325.0, 298.15),
             ("Jet-A(g)", "air", 1.0, 30e5, 800.0), ("Jet-A(g)", "air", 0.396, 30e5, 800.0),
             ("H2", "O2", 1.0, 1e5, 298.15), ("CH4", "O2", 1.0, 1e5, 298.15)]
    for fu, ox, ph, pp, Tox in cases:
        Tc = complete_flame_temperature(fu, ph, ox, T_REF, Tox)
        Te, re = equilibrium_flame_temperature(fu, ph, ox, pp, T_REF, Tox)
        print(f"  {fu}/{ox} phi {ph:.3f}, {pp / 1e5:5.2f} bar, oxidizer at {Tox:.0f} K: complete {Tc:7.1f} K, "
              f"equilibrium {Te:7.1f} K (lower by {Tc - Te:6.1f} K)")
    print("  Jet-A/air, air at 800 K, 30 bar (typical, illustrative combustor): equilibrium T_ad against phi")
    phis = np.round(np.arange(0.8, 1.41, 0.05), 3)
    Te_list = [equilibrium_flame_temperature("Jet-A(g)", ph, "air", 30e5, T_REF, 800.0)[0] for ph in phis]
    for ph, Te in zip(phis, Te_list):
        print(f"    phi {ph:.2f}: {Te:.1f} K (complete {complete_flame_temperature('Jet-A(g)', ph, 'air', T_REF, 800.0):.1f} K)")
    k = int(np.argmax(Te_list))
    print(f"    equilibrium peak near phi = {phis[k]:.2f}")

    print("\nillustrative rocket chamber: H2/O2 at O/F = 6 (mass), gaseous reactants at 298.15 K")
    phi_r = oxidizer_fuel_to_phi("H2", 6.0)
    for pc in (1e5, 20e6):
        Tc, rc = equilibrium_flame_temperature("H2", phi_r, "O2", pc)
        Mbar = float(np.sum(rc["x"] * np.array([molar_mass(q) for q in rc["species"]])))
        print(f"  p_c = {pc / 1e6:5.2f} MPa: phi = {phi_r:.4f}, T_c = {Tc:.1f} K (complete {complete_flame_temperature('H2', phi_r, 'O2'):.1f} K), "
              f"M = {Mbar * 1e3:.3f} g/mol")
        print(f"    {_fmt(rc)}")
    s2 = equilibrium_scipy(element_totals(reactants("H2", phi_r, "O2")), Tc, 20e6)
    print(f"  SLSQP at the chamber state: G/RT {s2['G_RT']:.10f} vs {rc['G_RT']:.10f}")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    Ts = np.linspace(1500.0, 4500.0, 61)
    sols, ln0 = [], None
    for T in Ts:
        q = equilibrium_tp(b_ho, T, 1e5, None, ln0)
        ln0 = q["ln_n"]
        sols.append(q["x"])
    sols = np.array(sols)
    for j, sp in enumerate(r["species"]):
        ax[0].semilogy(Ts, sols[:, j], label=sp)
    ax[0].set_ylim(1e-4, 1.2)
    ax[0].set_xlabel("T (K)")
    ax[0].set_ylabel("mole fraction")
    ax[0].set_title("H2 + O2/2 at equilibrium, 1 bar")
    ax[0].legend(fontsize=8)
    pp = np.linspace(0.6, 1.4, 33)
    ax[1].plot(pp, [complete_flame_temperature("Jet-A(g)", q, "air", T_REF, 800.0) for q in pp], "--", label="complete combustion")
    for pr, st in ((1e5, ":"), (30e5, "-")):
        ax[1].plot(pp, [equilibrium_flame_temperature("Jet-A(g)", q, "air", pr, T_REF, 800.0)[0] for q in pp], st,
                   label=f"equilibrium, {pr / 1e5:.0f} bar")
    ax[1].set_xlabel("equivalence ratio phi")
    ax[1].set_ylabel("adiabatic flame temperature (K)")
    ax[1].set_title("Jet-A/air, air at 800 K")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "equilibrium.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'equilibrium.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
