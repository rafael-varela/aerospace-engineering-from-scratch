"""Public tests for P07.10 Chemical Reactions & Combustion. They assert conservation, not magic numbers:
atoms balance exactly between reactants and products for lean, stoichiometric and rich mixtures;
mass is conserved; the fuel-air ratio and the equivalence ratio invert each other; the heating values
from enthalpies of formation match the standard enthalpies of combustion of methane and hydrogen
(OpenStax Chemistry 2e, Table 5.2) and obey HHV - LHV = (water per kg fuel) x h_fg from the steam table
(an independent data set); the energy balance H_reactants = H_products(T_ad) closes to round-off; T_ad
peaks at stoichiometric, rises with the inlet temperature and returns to the inlet temperature as
phi -> 0; and the fuel-air ratio solved for a turbine inlet temperature reproduces it.

COMPUTATIONAL ONLY: nothing here involves real combustion or fuel.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import pytest

from conftest import load_impl

impl = load_impl(__file__)

CASES = [("CH4", "air"), ("H2", "air"), ("Jet-A(g)", "air"), ("CH4", "O2"), ("H2", "O2")]


# ---------------------------------------------------------------- stoichiometry

def test_stoichiometric_oxygen():
    assert impl.o2_stoich("CH4") == 2.0
    assert impl.o2_stoich("H2") == 0.5
    assert impl.o2_stoich("Jet-A(g)") == pytest.approx(12 + 23 / 4)


def test_reactants_carry_air_in_the_fixed_proportions():
    r = impl.reactants("CH4", 0.8, "air")
    assert r["CH4"] == 1.0 and r["O2"] == pytest.approx(2.0 / 0.8)
    assert r["N2"] / r["O2"] == pytest.approx(0.7808 / 0.2095, rel=1e-14)
    assert r["Ar"] / r["O2"] == pytest.approx(0.0093 / 0.2095, rel=1e-14)
    r = impl.reactants("H2", 1.0, "O2")
    assert set(r) == {"H2", "O2"} and r["O2"] == 0.5


@pytest.mark.parametrize("fuel,ox", CASES)
def test_atoms_and_mass_balance_exactly(fuel, ox):
    for phi in (0.3, 0.7, 1.0, 1.2, 1.5):
        r = impl.reactants(fuel, phi, ox)
        p = impl.complete_products(fuel, phi, ox)
        er, ep = impl.element_totals(r), impl.element_totals(p)
        assert set(er) == set(ep)
        for el in er:
            assert ep[el] == pytest.approx(er[el], rel=1e-13, abs=1e-15)
        assert impl.mixture_mass(p) == pytest.approx(impl.mixture_mass(r), rel=1e-13)
        assert all(v >= 0.0 for v in p.values())


def test_product_regimes_lean_stoichiometric_rich():
    lean = impl.complete_products("CH4", 0.5, "O2")
    assert lean["CO2"] == pytest.approx(1.0) and lean["H2O"] == pytest.approx(2.0) and lean["O2"] == pytest.approx(2.0)
    st = impl.complete_products("CH4", 1.0, "O2")
    assert set(st) == {"CO2", "H2O"}                                          # nothing left over
    rich = impl.complete_products("H2", 2.0, "O2")
    assert rich["H2O"] == pytest.approx(0.5) and rich["H2"] == pytest.approx(0.5) and "O2" not in rich
    rich = impl.complete_products("CH4", 1.2, "O2")
    assert "O2" not in rich and rich["CO"] > 0 and rich["H2O"] == pytest.approx(2.0)
    with pytest.raises(ValueError):
        impl.complete_products("CH4", 5.0, "O2")                              # fewer O than C atoms: soot


def test_fuel_air_ratio_and_equivalence_ratio_invert():
    fst = impl.fuel_oxidizer_ratio("CH4", 1.0)
    # by hand: 16.04246 / (2 / 0.2095 * 28.9657) with the phase's air model
    assert fst == pytest.approx(16.04246 / (2 / 0.2095 * 28.96573), rel=1e-5)
    for fuel in impl.FUELS:
        for phi in (0.25, 0.8, 1.0, 1.3):
            f = impl.fuel_oxidizer_ratio(fuel, phi)
            assert f == pytest.approx(phi * impl.fuel_oxidizer_ratio(fuel, 1.0), rel=1e-13)
            assert impl.equivalence_ratio(fuel, f) == pytest.approx(phi, rel=1e-13)


# ---------------------------------------------------------------- heating values

def test_heating_values_from_enthalpies_of_formation():
    lhv, hhv = impl.heating_values("CH4")
    # Hess's law by hand from the tabulated Hf (J/mol): CH4 -74600, CO2 -393510, H2O -241826, H2O(L) -285830.
    # The polynomials reproduce Hf to about 2 J/mol (they were fitted with R = 8.314510).
    assert lhv * 16.04246e-3 == pytest.approx(-74600 + 393510 + 2 * 241826, abs=10.0)
    assert hhv * 16.04246e-3 == pytest.approx(-74600 + 393510 + 2 * 285830, abs=10.0)
    lhv2, hhv2 = impl.heating_values("H2")
    assert lhv2 * 2.01588e-3 == pytest.approx(241826, abs=5.0)
    assert hhv2 * 2.01588e-3 == pytest.approx(285830, abs=5.0)
    for fuel in impl.FUELS:
        l, h = impl.heating_values(fuel)
        assert h > l > 0


def test_higher_heating_values_against_standard_enthalpies_of_combustion():
    # OpenStax Chemistry 2e, section 5.3, Table 5.2 "Standard Molar Enthalpies of Combustion" (25 C, liquid
    # water): methane -890.8 kJ/mol, hydrogen -285.8 kJ/mol. An independent compilation, quoted to 0.1 kJ/mol.
    _, hhv = impl.heating_values("CH4")
    assert hhv * impl.molar_mass("CH4") == pytest.approx(890.8e3, rel=1e-3)
    _, hhv = impl.heating_values("H2")
    assert hhv * impl.molar_mass("H2") == pytest.approx(285.8e3, rel=1e-3)


@pytest.mark.parametrize("fuel", ["CH4", "H2", "Jet-A(g)"])
def test_hhv_minus_lhv_matches_steam_table_latent_heat(fuel):
    # Independent data: h_fg at 25 C from NIST SRD 69 (via Yan) times the water formed per kg of fuel.
    # The NASA Hf difference (44 004 J/mol) and the table agree to 0.04 %.
    l, h = impl.heating_values(fuel)
    w = impl.NASA9[fuel]["elements"]["H"] / 2 * impl.molar_mass("H2O") / impl.molar_mass(fuel)
    assert (h - l) == pytest.approx(w * (impl.H_G_25C - impl.H_F_25C) * 1e3, rel=1e-3)


# ---------------------------------------------------------------- adiabatic flame temperature

@pytest.mark.parametrize("fuel,ox", CASES)
def test_energy_balance_closes_at_the_flame_temperature(fuel, ox):
    for phi, Tox in ((0.6, 298.15), (1.0, 700.0), (1.3, 298.15)):
        Tad = impl.adiabatic_flame_temperature(fuel, phi, ox, 298.15, Tox)
        H_r = impl.reactant_enthalpy(fuel, phi, ox, 298.15, Tox)
        H_p = impl.mixture_enthalpy(impl.complete_products(fuel, phi, ox), Tad)
        assert abs(H_p - H_r) < 1e-9 * max(abs(H_r), 1e5)
        assert Tad > Tox


def test_flame_temperature_peaks_at_stoichiometric():
    for fuel, ox, Tin in (("Jet-A(g)", "air", 800.0), ("CH4", "air", 298.15), ("H2", "O2", 298.15)):
        phis = np.linspace(0.6, 1.4, 81)
        T = [impl.adiabatic_flame_temperature(fuel, p, ox, 298.15, Tin) for p in phis]
        assert 0.95 <= phis[int(np.argmax(T))] <= 1.1


def test_flame_temperature_limits_and_monotonicity():
    # phi -> 0: almost pure air, T_ad -> T_air
    assert impl.adiabatic_flame_temperature("Jet-A(g)", 1e-6, "air", 298.15, 800.0) == pytest.approx(800.0, abs=0.05)
    # hotter inlet air, hotter flame, but by less than the inlet rise at phi = 1 (products' cp is larger)
    a = impl.adiabatic_flame_temperature("Jet-A(g)", 1.0, "air", 298.15, 600.0)
    b = impl.adiabatic_flame_temperature("Jet-A(g)", 1.0, "air", 298.15, 800.0)
    assert 0 < b - a < 200.0
    # oxygen instead of air: no nitrogen to heat, a much hotter flame
    assert impl.adiabatic_flame_temperature("CH4", 1.0, "O2") > impl.adiabatic_flame_temperature("CH4", 1.0, "air") + 2000


def test_lean_flame_rise_matches_heating_value_estimate():
    # The first law split into two steps (burn at 298.15 K releasing the LHV, then heat the products):
    # T_ad - T0 = m_fuel LHV / (m_products cp_bar), with cp_bar the products' mean cp over [T0, T_ad].
    fuel, phi, T0 = "CH4", 0.05, 298.15
    Tad = impl.adiabatic_flame_temperature(fuel, phi)
    f = impl.fuel_oxidizer_ratio(fuel, phi)
    lhv, _ = impl.heating_values(fuel)
    p = impl.complete_products(fuel, phi)
    cp_bar = (impl.mixture_enthalpy(p, Tad) - impl.mixture_enthalpy(p, T0)) / (impl.mixture_mass(p) * (Tad - T0))
    m_air = impl.mixture_mass(p) - impl.molar_mass(fuel)
    assert (Tad - T0) == pytest.approx(f * lhv * m_air / (impl.mixture_mass(p) * cp_bar), rel=1e-9)


def test_phi_for_turbine_inlet_temperature():
    phi = impl.phi_for_temperature("Jet-A(g)", 1700.0, "air", 298.15, 800.0)
    assert 0.2 < phi < 0.6                                                    # gas turbines run lean overall
    assert impl.adiabatic_flame_temperature("Jet-A(g)", phi, "air", 298.15, 800.0) == pytest.approx(1700.0, abs=1e-6)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "flame_temperature.png").stat().st_size > 10_000
