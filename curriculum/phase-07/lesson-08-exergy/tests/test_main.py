"""Public tests for P07.08 Exergy. They assert thermodynamics, not magic numbers: atoms are conserved by
the combustion; the heating value follows from the NASA enthalpies of formation; the first law closes
for every component; the exergy destroyed in each component by T0 S_gen equals the destruction from an
exergy balance (the Gouy-Stodola theorem) and the whole audit closes to round-off; every destruction is
non-negative for any reasonable engine; the flow exergy vanishes at the dead state and equals R T0 ln(p/p0)
at T0; second-law efficiencies are at most 1 and equal 1 in the reversible limit; and every loss
moves exergy from the product into destruction.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import pytest

from conftest import load_impl

impl = load_impl(__file__)

BASE = impl.EngineSpec()
CRUISE = impl.EngineSpec(T0=216.65, p0=22632.0, V0=236.1)
REV = impl.EngineSpec(pi_d=1.0, eta_c=1.0, eta_t=1.0, eta_n=1.0, dp_cc=0.0)
COMPONENTS = ("inlet", "compressor", "combustor", "turbine", "nozzle")


# ---------------------------------------------------------------- mixtures and fuel

def test_combustion_conserves_atoms():
    for n_f in (0.0, 0.05, 0.2):
        before = impl.element_counts(impl.N_AIR, n_f)
        after = impl.element_counts(impl.products(n_f))
        for el in ("C", "H", "O", "N", "Ar"):
            assert after[el] == pytest.approx(before[el], rel=1e-14, abs=1e-14), el
    with pytest.raises(ValueError):
        impl.products(1.0)                                     # far richer than stoichiometric


def test_heating_value_from_formation_enthalpies():
    # LHV = [Hf(fuel) - 12 Hf(CO2) - 11.5 Hf(H2O)] / M_fuel with the database's Hf298 values (O2: 0).
    # NASA h(298.15 K) reproduces each Hf298 to about 2 J/mol (R = 8.314510 in the database, CODATA
    # 2018 here), so allow 40 J/mol of fuel: a relative error of 1e-5.
    lhv = (-249657.0 + 12 * 393510.0 + 11.5 * 241826.0) / impl.M_FUEL
    assert impl.lower_heating_value(298.15) == pytest.approx(lhv, rel=1e-5)
    # the reaction makes more moles of gas than it consumes, so Delta S > 0 and -Delta G > LHV
    assert impl.standard_reaction_exergy(298.15) > impl.lower_heating_value(298.15)
    assert impl.standard_reaction_exergy(298.15) < 1.05 * impl.lower_heating_value(298.15)


def test_mixture_entropy_includes_a_positive_mixing_term():
    T, p = 500.0, 2.0e5
    S_mix = impl.mixture_entropy(impl.N_AIR, T, p)
    S_pure = sum(impl.mixture_entropy(np.where(np.arange(5) == i, ni, 0.0), T, p) for i, ni in enumerate(impl.N_AIR) if ni > 0)
    assert S_mix > S_pure
    y = impl.N_AIR / impl.N_AIR.sum()
    mixing = -impl.R_UNIVERSAL * sum(ni * np.log(yi) for ni, yi in zip(impl.N_AIR, y) if ni > 0)
    assert S_mix - S_pure == pytest.approx(mixing, rel=1e-12)


def test_flow_exergy_vanishes_at_the_dead_state_and_has_the_right_sign_elsewhere():
    n, T0, p0 = impl.N_AIR, 288.15, 101325.0
    assert impl.flow_exergy(n, T0, p0, T0, p0) == pytest.approx(0.0, abs=1e-9)
    # at p0, any temperature other than T0 has positive exergy (hot or cold)
    for T in (200.0, 250.0, 500.0, 1400.0):
        assert impl.flow_exergy(n, T, p0, T0, p0) > 0.0
    # at T0 the flow exergy is R T0 ln(p/p0) per mole: positive above p0, negative below
    # (a stream below ambient pressure needs work to be pushed out into the atmosphere)
    R_mix = impl.R_UNIVERSAL * n.sum()
    for p in (0.3 * p0, 5 * p0):
        assert impl.flow_exergy(n, T0, p, T0, p0) == pytest.approx(R_mix * T0 * np.log(p / p0), rel=1e-10)
    # kinetic energy adds one for one
    assert impl.flow_exergy(n, T0, p0, T0, p0, V=100.0) == pytest.approx(5000.0, rel=1e-12)


# ---------------------------------------------------------------- the engine: first law

@pytest.mark.parametrize("spec", [BASE, CRUISE, REV], ids=["static", "cruise", "reversible"])
def test_first_law_closes_for_every_component(spec):
    res = impl.first_law_residuals(spec)
    assert set(res) == set(COMPONENTS)
    for k, v in res.items():
        assert abs(v) < 1e-8, k                         # J/kg of air, against energies of order 1e6
    st = impl.engine_states(spec)
    assert 0 < st["f"] < impl.N_AIR[1] / 17.75 * impl.M_FUEL       # lean of stoichiometric
    assert st["p"][5] > spec.p0 and st["V9"] > spec.V0
    assert st["T"][4] == spec.T_tit


def test_fuel_air_ratio_rises_with_turbine_inlet_temperature():
    f = [impl.engine_states(impl.EngineSpec(T_tit=t))["f"] for t in (1200.0, 1400.0, 1600.0)]
    assert np.all(np.diff(f) > 0)


# ---------------------------------------------------------------- Gouy-Stodola

@pytest.mark.parametrize("spec", [BASE, CRUISE, REV], ids=["static", "cruise", "reversible"])
def test_gouy_stodola_closes_to_round_off(spec):
    a = impl.exergy_audit(spec)
    total_in = a["Ex_fuel"] + a["Ex_ram"]
    assert abs(a["residual"]) < 1e-12 * total_in
    assert total_in == pytest.approx(a["product"] + a["Ex_ram"] + a["exhaust"] + sum(a["D"].values()), rel=1e-12)
    for k in COMPONENTS:
        assert a["D"][k] == pytest.approx(spec.T0 * a["S_gen"][k], rel=1e-14, abs=1e-12)
        assert a["D"][k] == pytest.approx(a["D_balance"][k], abs=1e-12 * total_in), k   # two routes, one answer


def test_every_destruction_is_non_negative_for_random_engines():
    rng = np.random.default_rng(7)
    for _ in range(12):
        spec = impl.EngineSpec(T0=rng.uniform(216.65, 310.0), V0=rng.uniform(0.0, 260.0), pi_d=rng.uniform(0.9, 1.0),
                               pi_c=rng.uniform(4.0, 25.0), eta_c=rng.uniform(0.75, 1.0), T_tit=rng.uniform(1100.0, 1700.0),
                               dp_cc=rng.uniform(0.0, 0.06), eta_t=rng.uniform(0.8, 1.0), eta_n=rng.uniform(0.9, 1.0))
        a = impl.exergy_audit(spec)
        for k in COMPONENTS:
            assert a["D"][k] >= -1e-9, (k, spec)
        assert a["exhaust"] > 0 and a["product"] > 0
        assert abs(a["residual"]) < 1e-12 * a["Ex_fuel"]


# ---------------------------------------------------------------- second-law efficiencies

def test_second_law_efficiencies_are_at_most_one_and_one_when_reversible():
    a = impl.exergy_audit(BASE)
    assert 0 < a["eta_II"] < 1
    assert a["eta_II"] < a["eta_th"]                    # fuel exergy exceeds the heating value here
    for k, v in a["eta_II_comp"].items():
        assert 0 < v < 1, k
    r = impl.exergy_audit(REV)
    for k, v in r["eta_II_comp"].items():
        assert v == pytest.approx(1.0, abs=1e-12), k
    for k in ("inlet", "compressor", "turbine", "nozzle"):
        assert abs(r["D"][k]) < 1e-9, k
    # combustion and the hot exhaust stay irreversible: the core still cannot reach eta_II = 1
    assert r["D"]["combustor"] > 0 and r["exhaust"] > 0 and r["eta_II"] < 1
    assert r["eta_II"] == pytest.approx(1 - (r["D"]["combustor"] + r["exhaust"]) / r["Ex_fuel"], rel=1e-12)


def test_fuel_exergy_is_close_to_the_reaction_gibbs_energy():
    for spec in (BASE, CRUISE):
        a = impl.exergy_audit(spec)
        per_kg_fuel = a["Ex_fuel"] / a["states"]["f"]
        assert impl.lower_heating_value() < per_kg_fuel < 1.06 * impl.lower_heating_value()
        assert per_kg_fuel == pytest.approx(impl.standard_reaction_exergy(298.15), rel=0.03)


def test_each_loss_raises_its_own_destruction_and_lowers_eta_II():
    a0 = impl.exergy_audit(BASE)
    for change, comp in (({"eta_c": 0.80}, "compressor"), ({"eta_t": 0.83}, "turbine"),
                         ({"eta_n": 0.92}, "nozzle"), ({"pi_d": 0.95}, "inlet")):
        a = impl.exergy_audit(impl.EngineSpec(**{**BASE.__dict__, **change}))
        assert a["D"][comp] > a0["D"][comp], comp
        assert a["eta_II"] < a0["eta_II"], comp


def test_combustor_dominates_and_shrinks_as_TIT_rises():
    a = impl.exergy_audit(BASE)
    assert a["D"]["combustor"] == max(a["D"].values())
    share = []
    for t in (1200.0, 1400.0, 1600.0):
        aa = impl.exergy_audit(impl.EngineSpec(T_tit=t))
        share.append(aa["D"]["combustor"] / aa["Ex_fuel"])
    assert np.all(np.diff(share) < 0)                   # hotter combustion, less destruction per unit fuel exergy


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "exergy_audit.png").stat().st_size > 10_000
