"""Public tests for P07.09 Gas Mixtures & Psychrometrics. They assert thermodynamics, not magic numbers:
Dalton's and Amagat's rules agree for ideal gases; mole and mass fractions convert both ways; the
entropy of mixing is non-negative with the right limits; the steam-table saturation pressure agrees
with an independent estimate from the NASA polynomials (the Gibbs equality); humidity ratio, relative
humidity and dew point are consistent inverses; adiabatic mixing conserves dry air, water and energy
and lands on the straight line between the inlet states in (omega, h); cooling below the dew point
removes exactly the excess water; and the cabin's water and energy balances close, with the cabin
humidity independent of recirculation.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import pytest

from conftest import load_impl

impl = load_impl(__file__)

R = 8.31446261815324


# ---------------------------------------------------------------- ideal-gas mixtures

def test_dry_air_molar_mass_and_fraction_round_trip():
    M = impl.molar_mass(impl.AIR_X)
    assert M == pytest.approx(28.965e-3, abs=2e-6)                   # the phase's dry-air model
    y = impl.mass_fractions(impl.AIR_X)
    assert sum(y.values()) == pytest.approx(1.0, abs=1e-14)
    assert y["O2"] > impl.AIR_X["O2"] and y["N2"] < impl.AIR_X["N2"]   # O2 heavier than the mean, N2 lighter
    x = impl.mole_fractions(y)
    for k in impl.AIR_X:
        assert x[k] == pytest.approx(impl.AIR_X[k], rel=1e-12)


def test_dalton_and_amagat_agree_for_ideal_gases():
    n = {"N2": 2.0, "O2": 0.5, "CO2": 0.1, "H2O": 0.4}
    T, V = 350.0, 0.08
    p = sum(n.values()) * R * T / V
    pp = impl.partial_pressures(n, T, V)
    pv = impl.partial_volumes(n, T, p)
    assert sum(pp.values()) == pytest.approx(p, rel=1e-13)              # Dalton: pressures add
    assert sum(pv.values()) == pytest.approx(V, rel=1e-13)              # Amagat: volumes add
    ntot = sum(n.values())
    for k in n:
        assert pp[k] / p == pytest.approx(n[k] / ntot, rel=1e-13)      # both fractions are mole fractions
        assert pv[k] / V == pytest.approx(n[k] / ntot, rel=1e-13)


def test_entropy_of_mixing_limits():
    assert impl.entropy_of_mixing({"N2": 1.0}) == 0.0                   # a pure gas does not "mix"
    assert impl.entropy_of_mixing({"N2": 0.5, "O2": 0.5}) == pytest.approx(R * np.log(2), rel=1e-13)
    assert impl.entropy_of_mixing({"N2": 1, "O2": 1, "Ar": 1}) == pytest.approx(R * np.log(3), rel=1e-13)
    vals = [impl.entropy_of_mixing({"N2": 1 - f, "O2": f}) for f in np.linspace(0.01, 0.99, 99)]
    assert min(vals) > 0 and np.argmax(vals) == 49                     # positive, largest at 50/50


def test_mixture_properties_are_mole_weighted_and_entropy_includes_mixing():
    T, p = 500.0, 2.0e5
    x = impl.AIR_X
    h = sum(xi * impl.h_species(k, T) for k, xi in x.items())
    assert impl.mixture_h(x, T) == pytest.approx(h, rel=1e-13)
    # mixture entropy minus the pure components at the same T, p equals the entropy of mixing
    s_pure = sum(xi * (impl.s_species(k, T) - R * np.log(p / impl.P_REF)) for k, xi in x.items())
    assert impl.mixture_s(x, T, p) - s_pure == pytest.approx(impl.entropy_of_mixing(x), rel=1e-10)
    # cp is dh/dT of the mixture
    dT = 1e-3
    assert impl.mixture_cp(x, T) == pytest.approx((impl.mixture_h(x, T + dT) - impl.mixture_h(x, T - dT)) / (2 * dT), rel=1e-7)


# ---------------------------------------------------------------- saturation pressure: independent check

def test_table_interpolation_reproduces_the_rows():
    for row in impl.SAT_WATER:
        T = max(row[0] + 273.15, impl.T_TRIPLE)
        assert impl.p_sat(T) == pytest.approx(row[1] * 1e6, rel=1e-12)
    with pytest.raises(ValueError):
        impl.p_sat(250.0)


def test_saturation_pressure_table_agrees_with_nasa_gibbs_equality():
    # Two independent routes: NIST SRD 69 table values (printed to 10 Pa) and g_vap = g_liq from the
    # NASA ideal-gas and liquid fits. Their difference is real-gas non-ideality of the vapour plus table
    # rounding; 1 % covers 0.01-60 C (the largest is 0.52 % at 60 C).
    for T in np.linspace(impl.T_TRIPLE, 333.15, 25):
        assert impl.p_sat_gibbs(T) == pytest.approx(impl.p_sat(T), rel=1e-2)
    # Clausius-Clapeyron: d ln p / d(1/T) = -h_fg M / R; the NASA estimate's slope gives h_fg near 2500 kJ/kg
    T1, T2 = 273.16, 278.15
    hfg = -R / impl.M_W * (np.log(impl.p_sat_gibbs(T2)) - np.log(impl.p_sat_gibbs(T1))) / (1 / T2 - 1 / T1)
    assert hfg == pytest.approx(impl.latent_heat_triple(), rel=5e-3)


def test_extrapolation_below_the_triple_point_is_continuous_and_monotone():
    assert impl.p_sat_extrapolated(impl.T_TRIPLE) == pytest.approx(impl.p_sat(impl.T_TRIPLE), rel=1e-12)
    T = np.linspace(250.0, 273.16, 30)
    p = [impl.p_sat_any(t) for t in T]
    assert np.all(np.diff(p) > 0)


# ---------------------------------------------------------------- moist-air state

def test_humidity_ratio_vapour_pressure_and_relative_humidity_are_consistent():
    p = 75_000.0
    for pv in (10.0, 500.0, 3000.0):
        w = impl.humidity_ratio(pv, p)
        assert impl.vapour_pressure(w, p) == pytest.approx(pv, rel=1e-13)
    T = 300.0
    assert impl.relative_humidity(T, impl.omega_saturated(T, p), p) == pytest.approx(1.0, rel=1e-13)
    assert impl.relative_humidity(T, 0.0, p) == 0.0
    assert impl.EPS == pytest.approx(18.01528 / 28.965, rel=1e-4)


def test_dew_point_inverts_saturation_and_flags_the_estimate():
    p = impl.P_ATM
    for T in (280.0, 300.0, 330.0):
        Td, est = impl.dew_point(impl.omega_saturated(T, p), p)          # saturated air: dew point = T
        assert Td == pytest.approx(T, abs=1e-6) and not est
    w = impl.humidity_ratio(0.5 * impl.p_sat(300.0), p)
    Td, _ = impl.dew_point(w, p)
    assert Td < 300.0 and impl.p_sat(Td) == pytest.approx(0.5 * impl.p_sat(300.0), rel=1e-9)
    Td, est = impl.dew_point(impl.humidity_ratio(300.0, p), p)          # 300 Pa: below the triple point
    assert est and Td < impl.T_TRIPLE
    assert impl.p_sat_extrapolated(Td) == pytest.approx(300.0, rel=1e-10)


def test_moist_enthalpy_references_and_vapour_against_table():
    assert impl.dry_air_enthalpy(273.15) == pytest.approx(0.0, abs=1e-9)
    assert impl.vapour_enthalpy(impl.T_TRIPLE) == pytest.approx(impl.SAT_WATER[0, 3] * 1e3, rel=1e-13)
    # ideal-gas vapour enthalpy against the table's saturated-vapour h_g (low pressure, so within 0.3 %)
    for row in impl.SAT_WATER[1:]:
        assert impl.vapour_enthalpy(row[0] + 273.15) == pytest.approx(row[3] * 1e3, rel=3e-3)
    # dry-air cp near 1005 J/(kg K) around room temperature (NASA mixture)
    cp = (impl.dry_air_enthalpy(305.0) - impl.dry_air_enthalpy(285.0)) / 20.0
    assert cp == pytest.approx(1005.0, rel=3e-3)
    h = impl.moist_air_enthalpy(310.0, 0.012)
    assert impl.temperature_from_enthalpy(h, 0.012) == pytest.approx(310.0, abs=1e-8)


# ---------------------------------------------------------------- processes and balances

def test_adiabatic_mixing_closes_balances_and_lies_on_the_straight_line():
    m1, T1, w1 = 0.7, 277.0, 0.004
    m2, T2, w2 = 0.4, 300.0, 0.009
    m3, T3, w3 = impl.mix_streams(m1, T1, w1, m2, T2, w2)
    assert m3 == pytest.approx(m1 + m2, rel=1e-15)                                    # dry air
    assert m3 * w3 == pytest.approx(m1 * w1 + m2 * w2, rel=1e-13)                     # water
    h1, h2, h3 = (impl.moist_air_enthalpy(T, w) for T, w in ((T1, w1), (T2, w2), (T3, w3)))
    assert m3 * h3 == pytest.approx(m1 * h1 + m2 * h2, rel=1e-11)                     # energy
    # collinear in (omega, h), with the lever rule: distances in inverse proportion to the flows
    cross = (w2 - w1) * (h3 - h1) - (w3 - w1) * (h2 - h1)
    assert abs(cross) < 1e-9 * abs((w2 - w1) * (h2 - h1))
    assert (w3 - w1) / (w2 - w3) == pytest.approx(m2 / m1, rel=1e-12)
    assert min(T1, T2) < T3 < max(T1, T2)
    # mixing a stream with itself changes nothing
    _, T4, w4 = impl.mix_streams(0.3, T1, w1, 0.5, T1, w1)
    assert T4 == pytest.approx(T1, abs=1e-8) and w4 == pytest.approx(w1, rel=1e-14)


def test_cooling_with_and_without_condensation_closes_water_and_energy():
    p, m = impl.P_ATM, 0.8
    w_in = impl.humidity_ratio(0.6 * impl.p_sat(308.15), p)
    Td, _ = impl.dew_point(w_in, p)
    dry = impl.cool_dehumidify(m, 308.15, w_in, Td + 2.0, p)                   # stays above the dew point
    assert dry["m_cond"] == 0.0 and dry["w_out"] == w_in and dry["phi_out"] < 1.0
    assert dry["Q"] == pytest.approx(m * (impl.moist_air_enthalpy(308.15, w_in) - impl.moist_air_enthalpy(Td + 2.0, w_in)), rel=1e-12)
    wet = impl.cool_dehumidify(m, 308.15, w_in, 275.15, p)                       # well below it
    assert wet["phi_out"] == pytest.approx(1.0, rel=1e-12)
    assert m * w_in == pytest.approx(m * wet["w_out"] + wet["m_cond"], rel=1e-13)   # water balance
    e_in = m * impl.moist_air_enthalpy(308.15, w_in)
    e_out = m * impl.moist_air_enthalpy(275.15, wet["w_out"]) + wet["m_cond"] * impl.liquid_enthalpy(275.15) + wet["Q"]
    assert e_out == pytest.approx(e_in, rel=1e-12)                                   # energy balance
    assert wet["Q"] > dry["Q"] > 0.0                                                 # latent heat makes it costlier


def test_cabin_humidity_set_by_fresh_air_not_recirculation():
    m_f, m_w, w_o = 0.775, 3.1e-3, 5e-5
    ws = [impl.cabin_humidity(m_f, m_r, m_w, w_o) for m_r in (0.0, 0.3, 0.775, 1.5)]
    for w_c, w_s in ws:
        assert w_c == pytest.approx(ws[0][0], rel=1e-14)
        assert m_f * w_o + m_w == pytest.approx(m_f * w_c, rel=1e-13)               # cabin + manifold water
        assert w_o * (1 - 1e-12) <= w_s <= w_c * (1 + 1e-12)
    assert ws[0][1] == pytest.approx(w_o, rel=1e-14)                                # no recirculation: supply = fresh
    w_a, _ = impl.cabin_humidity(m_f, 0.5, m_w, w_o)
    w_b, _ = impl.cabin_humidity(2 * m_f, 0.5, m_w, w_o)
    assert (w_b - w_o) == pytest.approx((w_a - w_o) / 2, rel=1e-12)                 # double the fresh air, half the rise


def test_ecs_cruise_balances_close_and_states_are_physical():
    e = impl.ecs_cruise()
    assert abs(e["res_air"]) < 1e-12 and abs(e["res_water"]) < 1e-15 and abs(e["res_energy"]) < 1e-6
    assert e["T_out"] == pytest.approx(216.65, abs=1e-9) and e["p_out"] == pytest.approx(22632.0, rel=1e-4)
    assert 0.0 < e["phi_fresh_at_cabin"] < 0.01                                     # outside air is bone dry
    assert 0.0 < e["phi_cabin"] < 1.0
    assert e["T_pack"] < e["T_supply"] < 295.15                                       # supply between pack and cabin
    assert e["T_dew_cabin"] < e["T_pack"]                                             # no condensation in the ducts
    assert e["Q_pack"] > 0.0
    # overall control volume, bleed port to outflow valve (pack, manifold and cabin together): the pack
    # removes the bleed enthalpy above the exhausted cabin air plus the cabin's sensible and vapour loads
    m_f, m_w, Tc = 186 * 0.25 / 60, 186 * 0.060 / 3600, 295.15
    q_overall = (m_f * (impl.moist_air_enthalpy(473.15, 5.0e-5) - impl.moist_air_enthalpy(Tc, e["w_cabin"]))
                 + 15.0e3 + m_w * impl.vapour_enthalpy(Tc))
    assert e["Q_pack"] == pytest.approx(q_overall, rel=1e-9)


def test_ecs_limits_recirculation_and_moisture():
    a = impl.ecs_cruise(recirc_fraction=0.0)
    b = impl.ecs_cruise(recirc_fraction=0.5)
    assert a["T_supply"] == pytest.approx(a["T_pack"], abs=1e-8)                      # no recirculation: supply is pack air
    assert a["phi_cabin"] == pytest.approx(b["phi_cabin"], rel=1e-12)                 # humidity unaffected
    assert a["T_pack"] == pytest.approx(b["T_pack"], abs=1e-6)                        # pack outlet set by the cabin load
    dry = impl.ecs_cruise(water_per_occ=0.0)
    assert dry["w_cabin"] == pytest.approx(dry["w_out"], rel=1e-14)                   # no people, no humidity rise
    more = impl.ecs_cruise(Q_sens=30e3)
    assert more["T_pack"] < b["T_pack"]                                               # more load needs colder pack air


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "psychrometrics.png").stat().st_size > 10_000
