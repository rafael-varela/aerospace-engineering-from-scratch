"""P07.10 Chemical Reactions & Combustion: starter code.

How hot can a combustor make the air it burns in, and how much fuel does it take to reach a turbine
inlet temperature? This module answers with the two conservation laws of a reacting system:

* atoms are conserved (stoichiometry): the oxygen a fuel CxHyOz needs, the fuel-air ratio, the
  equivalence ratio phi, and the products of complete combustion (lean: CO2, H2O and spare O2;
  rich: the oxygen goes to H2O first, then CO, then CO2, leaving H2 and CO);
* energy is conserved (first law for a steady, adiabatic, work-free combustor): the enthalpy of the
  reactants equals the enthalpy of the products, with every enthalpy measured from the elements
  (NASA polynomials include the enthalpy of formation), so the heat of reaction needs no separate
  bookkeeping;
* heating values: LHV and HHV at 298.15 K from enthalpies of formation, with water as vapour (LHV)
  or liquid (HHV); their difference checked against h_fg from the steam table of P07.03;
* the complete-combustion adiabatic flame temperature for Jet-A/air at combustor-inlet conditions,
  CH4/O2 and H2/O2, and the fuel-air ratio that gives a chosen turbine inlet temperature.

COMPUTATIONAL ONLY: there is no combustion, flame, burning, fuel handling or fuel of any kind in this
lesson. Every result is computed.

Complete combustion ignores dissociation (CO2 -> CO + O2/2, H2O -> OH, H, O ...); above about 2000 K
it overestimates the flame temperature, by more than 1000 K for the oxygen flames. P07.11 computes
the equilibrium composition and the lower, physical flame temperature.

Units: SI (K, Pa, mol, kg, J). Amounts are per mole of fuel unless stated. Enthalpies of formation
at 298.15 K and 1 bar (the NASA reference state).

Run:  python main.py [out_dir]   (prints the lesson's numbers, writes flame_temperature.png)
Test: pytest curriculum/phase-07/lesson-10-combustion
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

R_U = 8.31446261815324          # J/(mol K), CODATA 2018 (exact)
T_REF = 298.15                  # K, reference temperature of the enthalpies of formation

# Dry air model fixed for Phase 07 (mole fractions; M = 28.965 g/mol)
AIR_X = {"N2": 0.7808, "O2": 0.2095, "Ar": 0.0093, "CO2": 0.0004}

# Steam table row at 25 C (Yan, Introduction to Engineering Thermodynamics (LibreTexts), Appendix A,
# Table A1, values from NIST SRD 69): h_f and h_g in kJ/kg. Used only for the independent HHV - LHV check.
H_F_25C, H_G_25C = 104.83, 2546.51

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

FUELS = ("CH4", "H2", "Jet-A(g)")
PRODUCT_SPECIES = ("CO2", "H2O", "CO", "H2", "O2", "N2", "Ar")


# ---------------------------------------------------------------- given: NASA polynomials

def _coeffs(species: str, T: float) -> tuple[list[float], list[float]]:
    """The (a1..a7, b1 b2) coefficients of the interval holding T. Raises ValueError outside the fits."""
    for lo, hi, a, b in NASA9[species]["intervals"]:
        if lo <= T <= hi:
            return a, b
    raise ValueError(f"T = {T} K outside the fitted range for {species}")


def cp_species(species: str, T: float) -> float:
    """Molar heat capacity cp(T) [J/(mol K)]."""
    a, _ = _coeffs(species, T)
    return R_U * (a[0] / T ** 2 + a[1] / T + a[2] + a[3] * T + a[4] * T ** 2 + a[5] * T ** 3 + a[6] * T ** 4)


def h_species(species: str, T: float) -> float:
    """Molar enthalpy h(T) [J/mol] measured from the elements: h(298.15 K) = Hf298."""
    a, b = _coeffs(species, T)
    return R_U * T * (-a[0] / T ** 2 + a[1] * np.log(T) / T + a[2] + a[3] * T / 2 + a[4] * T ** 2 / 3
                      + a[5] * T ** 3 / 4 + a[6] * T ** 4 / 5 + b[0] / T)


def molar_mass(species: str) -> float:
    """Molar mass [kg/mol]."""
    return NASA9[species]["M"] / 1000.0


def mixture_mass(n: dict[str, float]) -> float:
    """Mass [kg] of the amounts n [mol]."""
    return sum(ni * molar_mass(s) for s, ni in n.items())


def element_totals(n: dict[str, float]) -> dict[str, float]:
    """Moles of each element in the amounts n [mol]: sum_j a_ij n_j."""
    out: dict[str, float] = {}
    for s, nj in n.items():
        for el, a in NASA9[s]["elements"].items():
            out[el] = out.get(el, 0.0) + a * nj
    return out


# ---------------------------------------------------------------- stoichiometry

def o2_stoich(fuel: str) -> float:
    """Moles of O2 needed to burn one mole of fuel CxHyOz completely to CO2 and H2O: x + y/4 - z/2."""
    raise NotImplementedError


def reactants(fuel: str, phi: float, oxidizer: str = "air") -> dict[str, float]:
    """Reactant amounts per mole of fuel at equivalence ratio phi (phi < 1 lean, > 1 rich):
    O2 = o2_stoich / phi; with oxidizer 'air' the N2, Ar and CO2 of the dry-air model come along in
    proportion (AIR_X[k] / AIR_X['O2'] per mole of O2); with 'O2' the oxidizer is pure oxygen."""
    raise NotImplementedError


def fuel_oxidizer_ratio(fuel: str, phi: float, oxidizer: str = "air") -> float:
    """Mass of fuel per mass of oxidizer (the fuel-air ratio f for air): f = phi f_stoich."""
    raise NotImplementedError


def equivalence_ratio(fuel: str, f: float, oxidizer: str = "air") -> float:
    """phi = f / f_stoich: the inverse of fuel_oxidizer_ratio."""
    raise NotImplementedError


def complete_products(fuel: str, phi: float, oxidizer: str = "air") -> dict[str, float]:
    """Products of complete combustion per mole of fuel, from atom conservation alone.
    Lean or stoichiometric: all C to CO2, all H to H2O, spare O2. Rich, a fixed order with no
    equilibrium: every carbon atom first takes one O (CO), the hydrogen then takes what it can (H2O,
    the rest left as H2), and any oxygen still left turns CO into CO2. Raises ValueError if there are
    fewer O atoms than C atoms (solid carbon would form; outside this model)."""
    raise NotImplementedError


# ---------------------------------------------------------------- energy

def mixture_enthalpy(n: dict[str, float], T: float) -> float:
    """Enthalpy [J] of the amounts n [mol] at T, sum n_j h_j(T), measured from the elements."""
    raise NotImplementedError


def reactant_enthalpy(fuel: str, phi: float, oxidizer: str = "air", T_fuel: float = T_REF,
                      T_ox: float = T_REF) -> float:
    """Enthalpy [J per mole of fuel] of the reactants, fuel at T_fuel and oxidizer at T_ox."""
    raise NotImplementedError


def heating_values(fuel: str) -> tuple[float, float]:
    """Lower and higher heating values [J/kg of fuel] at 298.15 K: the enthalpy released when one
    mole of fuel and its stoichiometric O2 at 298.15 K become CO2 and H2O at 298.15 K,
    LHV = Hf(fuel) - sum nu Hf(products, H2O vapour), HHV with liquid water H2O(L)."""
    raise NotImplementedError


def adiabatic_flame_temperature(fuel: str, phi: float, oxidizer: str = "air", T_fuel: float = T_REF,
                                T_ox: float = T_REF) -> float:
    """Complete-combustion adiabatic flame temperature [K]: the T at which the products' enthalpy
    equals the reactants' (steady flow, no heat loss, no work, no kinetic energy change; for ideal
    gases h does not depend on pressure, so neither does this T). Solved by brentq on
    H_products(T) - H_reactants = 0 between 200 K and 6000 K."""
    raise NotImplementedError


def phi_for_temperature(fuel: str, T_target: float, oxidizer: str = "air", T_fuel: float = T_REF,
                        T_ox: float = T_REF) -> float:
    """Lean equivalence ratio at which the complete-combustion products reach T_target (a turbine
    inlet temperature), by brentq on phi in (1e-4, 1)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the lesson's numbers

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save flame_temperature.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    print("COMPUTATIONAL ONLY: no combustion, flame or fuel of any kind; every number below is computed.\n")
    print("fuel       O2/mol fuel   f_stoich (air)   AFR_stoich   LHV (MJ/kg)   HHV (MJ/kg)   HHV-LHV vs table h_fg")
    for fu in FUELS:
        lhv, hhv = heating_values(fu)
        el = NASA9[fu]["elements"]
        w = el.get("H", 0) / 2.0 * molar_mass("H2O") / molar_mass(fu)          # kg water per kg fuel
        check = w * (H_G_25C - H_F_25C) * 1e3
        fst = fuel_oxidizer_ratio(fu, 1.0)
        print(f"{fu:9s}  {o2_stoich(fu):8.2f}      {fst:.5f}          {1 / fst:6.2f}      {lhv / 1e6:7.3f}      "
              f"{hhv / 1e6:7.3f}       {(hhv - lhv) / 1e6:.4f} vs {check / 1e6:.4f} ({100 * ((hhv - lhv) / check - 1):+.3f} %)")

    print("\ncomplete-combustion adiabatic flame temperatures, reactants at 298.15 K, phi = 1:")
    for fu, ox in (("CH4", "air"), ("H2", "air"), ("Jet-A(g)", "air"), ("CH4", "O2"), ("H2", "O2")):
        print(f"  {fu}/{ox}: {adiabatic_flame_temperature(fu, 1.0, ox):.1f} K")

    T3, T4 = 800.0, 1700.0
    print(f"\nJet-A(g)/air, air at T3 = {T3:.0f} K, fuel vapour at 298.15 K (typical, illustrative combustor inlet):")
    print("  phi     f         T_ad (K)   products (mol per mol fuel)")
    phis = np.round(np.arange(0.3, 1.51, 0.1), 2)
    Tads = []
    for ph in phis:
        Tad = adiabatic_flame_temperature("Jet-A(g)", ph, "air", T_REF, T3)
        Tads.append(Tad)
        prod = complete_products("Jet-A(g)", ph)
        print(f"  {ph:.1f}   {fuel_oxidizer_ratio('Jet-A(g)', ph):.5f}   {Tad:8.1f}   "
              + ", ".join(f"{k} {v:.2f}" for k, v in prod.items() if k in ("CO2", "H2O", "CO", "H2", "O2")))
    fine = np.linspace(0.9, 1.2, 301)
    Tf = [adiabatic_flame_temperature("Jet-A(g)", p, "air", T_REF, T3) for p in fine]
    print(f"  peak {max(Tf):.1f} K at phi = {fine[int(np.argmax(Tf))]:.3f}")
    ph4 = phi_for_temperature("Jet-A(g)", T4, "air", T_REF, T3)
    f4 = fuel_oxidizer_ratio("Jet-A(g)", ph4)
    print(f"  turbine inlet {T4:.0f} K needs phi = {ph4:.4f}, f = {f4:.5f} (1 kg fuel per {1 / f4:.1f} kg air); "
          f"energy check f LHV = {f4 * heating_values('Jet-A(g)')[0] / 1e3:.1f} kJ/kg air")
    H_r = reactant_enthalpy("Jet-A(g)", ph4, "air", T_REF, T3)
    H_p = mixture_enthalpy(complete_products("Jet-A(g)", ph4), T4)
    print(f"  energy balance residual at that point: {H_p - H_r:.2e} J per mol fuel (|H_r| = {abs(H_r):.4e})")
    print("  sensitivity to the inlet temperature at phi = 1 and at the turbine-inlet phi:")
    for Tin in (298.15, 600.0, 800.0, 900.0):
        print(f"    T3 = {Tin:6.1f} K: T_ad(phi=1) = {adiabatic_flame_temperature('Jet-A(g)', 1.0, 'air', T_REF, Tin):.1f} K, "
              f"T_ad(phi={ph4:.3f}) = {adiabatic_flame_temperature('Jet-A(g)', ph4, 'air', T_REF, Tin):.1f} K")

    print("\nH2/O2 and CH4/O2 against phi (reactants at 298.15 K):")
    for ph in (0.5, 0.8, 1.0, 1.2, 1.5, 2.0):
        a = adiabatic_flame_temperature("H2", ph, "O2")
        b = adiabatic_flame_temperature("CH4", ph, "O2")
        print(f"  phi {ph:.1f}: H2/O2 {a:.0f} K, CH4/O2 {b:.0f} K")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    pp = np.linspace(0.3, 1.5, 121)
    ax[0].plot(pp, [adiabatic_flame_temperature("Jet-A(g)", p, "air", T_REF, T3) for p in pp], label="Jet-A/air, air at 800 K")
    ax[0].plot(pp, [adiabatic_flame_temperature("Jet-A(g)", p, "air") for p in pp], label="Jet-A/air, 298 K")
    ax[0].plot(pp, [adiabatic_flame_temperature("CH4", p, "air") for p in pp], "--", label="CH4/air, 298 K")
    ax[0].plot(pp, [adiabatic_flame_temperature("H2", p, "air") for p in pp], ":", label="H2/air, 298 K")
    ax[0].axhline(T4, color="k", lw=0.6)
    ax[0].plot([ph4], [T4], "ko")
    ax[0].annotate(f"turbine inlet {T4:.0f} K at phi = {ph4:.3f}", (ph4, T4), textcoords="offset points", xytext=(8, -14), fontsize=8)
    ax[0].set_xlabel("equivalence ratio phi")
    ax[0].set_ylabel("complete-combustion T_ad (K)")
    ax[0].set_title("adiabatic flame temperature (no dissociation)")
    ax[0].legend(fontsize=8)
    for sp in ("CO2", "H2O", "O2", "CO", "H2"):
        ax[1].plot(pp, [complete_products("Jet-A(g)", p).get(sp, 0.0) for p in pp], label=sp)
    ax[1].set_xlabel("equivalence ratio phi")
    ax[1].set_ylabel("mol per mol Jet-A")
    ax[1].set_title("complete-combustion products, Jet-A/air")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "flame_temperature.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'flame_temperature.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
