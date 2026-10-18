"""Public tests for P07.11 Chemical Equilibrium. They assert thermodynamics, not magic numbers: atoms
are conserved to round-off; the element-potential solution satisfies the Lagrange conditions
mu_j/RT = sum_i a_ij pi_i; an independent minimiser (scipy SLSQP) reaches the same Gibbs energy and
no feasible perturbation lowers it; the law of mass action holds for every independent reaction with
K_p = exp(-Delta G°/RT); the answer does not depend on the starting guess; Le Chatelier's principle
holds in T and p; at low temperature the composition is complete combustion; and the equilibrium
adiabatic flame temperature closes its energy balance and lies below the complete-combustion value
of P07.10 for lean and stoichiometric mixtures.

COMPUTATIONAL ONLY: nothing here involves real combustion or fuel.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import pytest

from conftest import load_impl

impl = load_impl(__file__)

B_HO = {"H": 2.0, "O": 1.0}                                     # H2 + O2/2


def b_chon(phi=1.0, fuel="CH4"):
    return impl.element_totals(impl.reactants(fuel, phi, "air"))


# ---------------------------------------------------------------- Gibbs energy and reactions

def test_gibbs_energy_of_pure_gas_and_of_a_mixture():
    T, p = 2000.0, 3e5
    g = impl.g_species("H2O", T) / (8.31446261815324 * T)
    assert impl.gibbs_rt(np.array([2.0]), ["H2O"], T, p) == pytest.approx(2 * (g + np.log(3.0)), rel=1e-13)
    # mixing lowers G: two pure gases at T, p versus their mixture at T, p
    n = np.array([1.0, 1.0])
    pure = impl.gibbs_rt(n[:1], ["H2"], T, p) + impl.gibbs_rt(n[1:], ["O2"], T, p)
    assert impl.gibbs_rt(n, ["H2", "O2"], T, p) == pytest.approx(pure - 2 * np.log(2), rel=1e-12)


def test_equilibrium_constants_compose_like_reactions():
    T = 2500.0
    K1 = impl.equilibrium_constant({"H2O": 1, "H2": -1, "O2": -0.5}, T)
    K2 = impl.equilibrium_constant({"H2": 1, "H2O": -1, "O2": 0.5}, T)
    assert K1 * K2 == pytest.approx(1.0, rel=1e-12)                         # reverse reaction
    Ka = impl.equilibrium_constant({"OH": 1, "H2": -0.5, "O2": -0.5}, T)
    Kb = impl.equilibrium_constant({"H": 1, "H2": -0.5}, T)
    Kc = impl.equilibrium_constant({"OH": 1, "H": -1, "O2": -0.5}, T)
    assert Kc == pytest.approx(Ka / Kb, rel=1e-12)                           # Hess's law for K
    # dissociation favoured by temperature
    assert impl.equilibrium_constant({"H": 2, "H2": -1}, 3500.0) > impl.equilibrium_constant({"H": 2, "H2": -1}, 2500.0)


def test_formation_reactions_balance_atoms():
    sp, els = impl.species_for(b_chon())
    base = ["CO2", "H2O", "O2", "N2", "Ar"]
    rx = impl.formation_reactions(sp, base)
    assert len(rx) == len(sp) - len(base)                                     # one per non-base species
    for r in rx:
        tot = {}
        for s, nu in r.items():
            for e, a in impl.NASA9[s]["elements"].items():
                tot[e] = tot.get(e, 0.0) + nu * a
        assert all(abs(v) < 1e-12 for v in tot.values())


# ---------------------------------------------------------------- the solver

@pytest.mark.parametrize("b0,T,p", [(B_HO, 3000.0, 1e5), (B_HO, 1500.0, 2e7), (b_chon(), 2200.0, 101325.0),
                                    (b_chon(1.3, "Jet-A(g)"), 2500.0, 30e5)])
def test_elements_conserved_and_lagrange_conditions_hold(b0, T, p):
    r = impl.equilibrium_tp(b0, T, p)
    A = impl.element_matrix(r["species"], r["elements"])
    b = np.array([b0[e] for e in r["elements"]])
    assert np.max(np.abs(A @ r["n"] - b)) < 1e-12 * np.max(b)
    assert r["element_residual"] < 1e-12
    assert np.all(r["n"] > 0) and r["x"].sum() == pytest.approx(1.0, abs=1e-14)
    # Lagrange: every species' chemical potential is the sum of its atoms' element potentials
    g0 = np.array([impl.g_species(s, T) for s in r["species"]]) / (8.31446261815324 * T)
    mu = g0 + np.log(r["x"]) + np.log(p / impl.P_REF)
    npt_ok = np.abs(mu - A.T @ r["pi"])
    assert np.max(npt_ok) < 1e-6


@pytest.mark.parametrize("b0,T,p", [(B_HO, 3000.0, 1e5), (b_chon(), 2200.0, 101325.0), (b_chon(0.7, "Jet-A(g)"), 2400.0, 30e5)])
def test_second_minimiser_reaches_the_same_gibbs_energy(b0, T, p):
    r = impl.equilibrium_tp(b0, T, p)
    s = impl.equilibrium_scipy(b0, T, p)
    assert s["element_residual"] < 1e-8
    assert abs(s["G_RT"] - r["G_RT"]) < 1e-9 * abs(r["G_RT"])
    assert s["G_RT"] >= r["G_RT"] - 1e-9 * abs(r["G_RT"])                     # nothing beats the true minimum
    major = r["x"] > 1e-3
    np.testing.assert_allclose(s["x"][major], r["x"][major], rtol=1e-4)


def test_no_feasible_perturbation_lowers_gibbs_energy():
    T, p = 2800.0, 1e5
    r = impl.equilibrium_tp(B_HO, T, p)
    A = impl.element_matrix(r["species"], r["elements"])
    _, _, Vt = np.linalg.svd(A)
    null = Vt[A.shape[0]:]                                                   # directions that keep A n = b
    rng = np.random.default_rng(7)
    for _ in range(20):
        d = rng.standard_normal(null.shape[0]) @ null
        for eps in (1e-3, 1e-2):
            n = r["n"] + eps * d * r["n"].min()
            if np.all(n > 0):
                assert impl.gibbs_rt(n, r["species"], T, p) > r["G_RT"]


@pytest.mark.parametrize("b0,T,p,base", [(B_HO, 3000.0, 1e5, ["H2", "O2"]), (B_HO, 2000.0, 1e7, ["H2O", "O2"]),
                                         (b_chon(), 2200.0, 101325.0, ["CO2", "H2O", "O2", "N2", "Ar"]),
                                         (b_chon(1.2, "Jet-A(g)"), 2400.0, 30e5, ["CO", "H2", "H2O", "N2", "Ar"])])
def test_law_of_mass_action_for_every_independent_reaction(b0, T, p, base):
    r = impl.equilibrium_tp(b0, T, p)
    res = impl.mass_action_residuals(r, T, p, base)
    assert len(res) == len(r["species"]) - len(base)
    for rx, lnK, lnQ in res:
        assert lnQ == pytest.approx(lnK, abs=1e-6)


def test_answer_does_not_depend_on_the_starting_guess():
    b0, T, p = b_chon(0.9), 2300.0, 5e5
    a = impl.equilibrium_tp(b0, T, p)
    rng = np.random.default_rng(3)
    b = impl.equilibrium_tp(b0, T, p, a["species"], rng.uniform(-12, 1, len(a["species"])))
    np.testing.assert_allclose(b["x"], a["x"], rtol=1e-6, atol=1e-14)


def test_le_chatelier_in_pressure_and_temperature():
    # H2O <-> H2 + O2/2 makes more moles: higher p favours H2O, higher T favours dissociation
    xw = lambda T, p: (lambda r: r["x"][r["species"].index("H2O")])(impl.equilibrium_tp(B_HO, T, p))
    assert xw(3000.0, 1e4) < xw(3000.0, 1e5) < xw(3000.0, 1e6) < xw(3000.0, 1e7)
    assert xw(2000.0, 1e5) > xw(2500.0, 1e5) > xw(3000.0, 1e5) > xw(3500.0, 1e5)


def test_low_temperature_limit_is_complete_combustion():
    r = impl.equilibrium_tp(B_HO, 1000.0, 1e5)
    assert r["x"][r["species"].index("H2O")] > 1 - 1e-6
    r = impl.equilibrium_tp(b_chon(0.8), 1000.0, 101325.0)
    comp = impl.complete_products("CH4", 0.8, "air")
    ntot = sum(comp.values())
    for k, v in comp.items():
        assert r["x"][r["species"].index(k)] == pytest.approx(v / ntot, abs=1e-4)


# ---------------------------------------------------------------- the flame

@pytest.mark.parametrize("fuel,ox,phi,p,Tox", [("CH4", "air", 1.0, 101325.0, 298.15), ("H2", "O2", 1.0, 1e5, 298.15),
                                              ("Jet-A(g)", "air", 0.8, 30e5, 800.0)])
def test_equilibrium_flame_energy_balance_and_below_complete_combustion(fuel, ox, phi, p, Tox):
    T, r = impl.equilibrium_flame_temperature(fuel, phi, ox, p, 298.15, Tox)
    H_r = impl.reactant_enthalpy(fuel, phi, ox, 298.15, Tox)
    H_p = float(np.sum(r["n"] * np.array([impl.h_species(s, T) for s in r["species"]])))
    assert abs(H_p - H_r) < 1e-6 * max(abs(H_r), 1e5)
    assert T < impl.complete_flame_temperature(fuel, phi, ox, 298.15, Tox)


def test_flame_temperature_rises_with_pressure_and_gap_is_largest_near_stoichiometric():
    lo = impl.equilibrium_flame_temperature("H2", 1.0, "O2", 1e5)[0]
    hi = impl.equilibrium_flame_temperature("H2", 1.0, "O2", 1e7)[0]
    assert hi > lo + 100.0                                                    # less dissociation at high p
    gaps = []
    for phi in (0.4, 0.7, 1.0):
        Tc = impl.complete_flame_temperature("Jet-A(g)", phi, "air", 298.15, 800.0)
        Te = impl.equilibrium_flame_temperature("Jet-A(g)", phi, "air", 30e5, 298.15, 800.0)[0]
        gaps.append(Tc - Te)
    assert 0 < gaps[0] < gaps[1] < gaps[2]


def test_rocket_chamber_mixture_ratio():
    phi = impl.oxidizer_fuel_to_phi("H2", 6.0)
    assert phi == pytest.approx(0.5 * 31.9988 / 2.01588 / 6.0, rel=1e-12)       # stoichiometric O/F = 7.94
    Tc, r = impl.equilibrium_flame_temperature("H2", phi, "O2", 20e6)
    assert 3000.0 < Tc < impl.complete_flame_temperature("H2", phi, "O2")
    assert r["x"][r["species"].index("H2")] > 0.2                              # fuel-rich: much free H2


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "equilibrium.png").stat().st_size > 10_000
