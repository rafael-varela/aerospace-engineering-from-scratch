"""Public tests for P07.04 The Second Law & Entropy. They assert the second law and the property nature of
entropy, not magic numbers: S_gen >= 0 for every process over seeded random inputs; S_gen -> 0 as the
temperature difference -> 0 (linearly for heat across dT, quadratically for two bodies equilibrating);
entropy change is path-independent (three different paths and the end-state formula agree, the trapezoid
rule converging at second order); energy balances of the exchanger and the mixer close; the ledger refuses
negative entries; and every cycle the tracker accepts obeys the Carnot bound, while every cycle that beats
Carnot is rejected.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. The only file written goes to tmp_path.
"""
import numpy as np
import pytest

from conftest import load_impl

impl = load_impl(__file__)

R = impl.R_AIR


# ---------------------------------------------------------------- entropy is a property

def test_end_state_formula_matches_constant_cp_limits():
    # isothermal: -R ln(p2/p1) exactly; small isobaric step: cp ln(T2/T1)
    assert impl.ds_ideal_gas(400.0, 1e5, 400.0, 3e5) == pytest.approx(-R * np.log(3.0), rel=1e-12)
    T1, T2 = 500.0, 501.0
    assert impl.ds_ideal_gas(T1, 2e5, T2, 2e5) == pytest.approx(impl.air_cp(500.5) * np.log(T2 / T1), rel=1e-6)
    # antisymmetric
    assert impl.ds_ideal_gas(300.0, 1e5, 700.0, 8e5) == pytest.approx(-impl.ds_ideal_gas(700.0, 8e5, 300.0, 1e5))


def test_entropy_change_is_path_independent():
    T1, p1, T2, p2 = 320.0, 0.9e5, 760.0, 1.2e6
    exact = impl.ds_ideal_gas(T1, p1, T2, p2)
    paths = [
        (lambda l: np.where(l < 0.5, T1 + (T2 - T1) * 2 * l, T2), lambda l: np.where(l < 0.5, p1, p1 + (p2 - p1) * (2 * l - 1))),
        (lambda l: np.where(l < 0.5, T1, T1 + (T2 - T1) * (2 * l - 1)), lambda l: np.where(l < 0.5, p1 + (p2 - p1) * 2 * l, p2)),
        (lambda l: T1 + (T2 - T1) * l + 120 * np.sin(3 * np.pi * l), lambda l: p1 * (p2 / p1) ** l * (1 + 0.4 * np.sin(2 * np.pi * l) ** 2)),
    ]
    for Tp, pp in paths:
        assert impl.ds_path_integral(Tp, pp, 20001) == pytest.approx(exact, rel=1e-6)
    # a closed loop returns to the start: zero entropy change for the gas
    loop = impl.ds_path_integral(lambda l: 500 + 200 * np.sin(2 * np.pi * l), lambda l: 3e5 * (1 + 0.5 * np.sin(4 * np.pi * l)), 20001)
    assert abs(loop) < 1e-5


def test_path_integral_converges_at_second_order():
    T1, p1, T2, p2 = 300.0, 1e5, 800.0, 1e6
    exact = impl.ds_ideal_gas(T1, p1, T2, p2)
    Tp = lambda l: T1 + (T2 - T1) * l + 150 * np.sin(3 * np.pi * l)
    pp = lambda l: p1 * (p2 / p1) ** l * (1 + 0.5 * np.sin(2 * np.pi * l) ** 2)
    errs = np.array([abs(impl.ds_path_integral(Tp, pp, n) - exact) for n in (101, 201, 401, 801)])
    orders = np.log2(errs[:-1] / errs[1:])
    assert np.all((orders > 1.8) & (orders < 2.2))


def test_irreversible_and_reversible_routes_give_the_same_property_change():
    # adiabatic throttling (irreversible, q = 0) and a reversible isothermal expansion between the same
    # states give the same Delta s; in the throttle all of it is generated, in the reversible path all of it is q/T
    p1, p2, T = 1.0e6, 0.25e6, 600.0
    ds = impl.ds_ideal_gas(T, p1, T, p2)
    assert impl.throttle_sgen(p1, p2) == pytest.approx(ds, rel=1e-12)
    q_rev = R * T * np.log(p1 / p2)                                          # heat into a reversible isothermal expansion
    assert q_rev / T == pytest.approx(ds, rel=1e-12)


# ---------------------------------------------------------------- S_gen >= 0 and -> 0

def test_every_process_generates_non_negative_entropy_over_random_inputs():
    rng = np.random.default_rng(404)
    for _ in range(200):
        Ta, Tb = rng.uniform(210.0, 990.0, 2)
        assert impl.heat_transfer_sgen(rng.uniform(1, 1e4), max(Ta, Tb), min(Ta, Tb)) >= 0
        assert impl.bodies_equilibrate(rng.uniform(0.1, 5), rng.uniform(100, 4000), Ta, rng.uniform(0.1, 5), rng.uniform(100, 4000), Tb)[1] >= 0
        pa, pb = rng.uniform(1e4, 3e6, 2)
        assert impl.throttle_sgen(max(pa, pb), min(pa, pb)) >= 0
        m = rng.uniform(0.1, 3, 3)
        Ts = rng.uniform(220.0, 900.0, 3)
        ps = rng.uniform(1e5, 5e5, 3)
        assert impl.mix_streams(m, Ts, ps, ps.min())[1] >= -1e-9
        assert impl.liquid_line_friction(rng.uniform(0.1, 10), rng.uniform(0, 1e6), rng.uniform(500, 1500), rng.uniform(1000, 4000), Ta)[1] >= 0
        n = rng.uniform(0.0, 2.0, 3)
        assert impl.mix_gases_sgen({"a": n[0], "b": n[1], "c": n[2]}) >= 0


def test_sgen_vanishes_as_the_temperature_difference_vanishes():
    dT = np.array([10.0, 1.0, 0.1, 0.01])
    s_q = np.array([impl.heat_transfer_sgen(1000.0, 300.0 + d, 300.0) for d in dT])
    s_b = np.array([impl.bodies_equilibrate(1.0, 900.0, 300 + d / 2, 1.0, 900.0, 300 - d / 2)[1] for d in dT])
    assert np.all(np.diff(s_q) < 0) and np.all(np.diff(s_b) < 0)
    # slopes: S_gen ~ dT for fixed Q, ~ dT^2 for two bodies (C dT^2 / (4 T^2) for equal bodies)
    assert np.allclose(np.log10(s_q[:-1] / s_q[1:]), 1.0, atol=0.05)
    assert np.allclose(np.log10(s_b[:-1] / s_b[1:]), 2.0, atol=0.01)
    assert s_b[-1] == pytest.approx(900.0 * 0.01 ** 2 / (4 * 300.0 ** 2), rel=1e-3)
    assert impl.heat_transfer_sgen(500.0, 300.0, 300.0) == 0.0
    assert impl.mix_gases_sgen({"N2": 1.0}) == 0.0
    assert impl.liquid_line_friction(1.0, 0.0, 800.0, 2000.0, 300.0)[1] == 0.0


def test_throttle_refuses_a_pressure_rise():
    with pytest.raises(ValueError):
        impl.throttle_sgen(1e5, 2e5)


# ---------------------------------------------------------------- open-system balances

def test_heat_exchanger_balances_and_refuses_a_temperature_cross():
    r = impl.heat_exchanger(1.0, 700.0, 470.0, 1.2, 320.0)
    assert r["Q"] == pytest.approx(1.0 * (impl.air_h(700.0) - impl.air_h(470.0)))
    assert 1.2 * (impl.air_h(r["T_c_out"]) - impl.air_h(320.0)) == pytest.approx(r["Q"], rel=1e-10)
    assert r["S_gen"] > 0
    # a warmer cold stream (smaller temperature difference) lowers the generation; more cold flow at the
    # same inlet keeps the cold side colder on average and raises it
    assert impl.heat_exchanger(1.0, 700.0, 470.0, 1.2, 400.0)["S_gen"] < r["S_gen"]
    assert impl.heat_exchanger(1.0, 700.0, 470.0, 5.0, 320.0)["S_gen"] > r["S_gen"]
    with pytest.raises(ValueError):
        impl.heat_exchanger(1.0, 700.0, 300.0, 1.0, 320.0)                   # hot outlet below cold inlet


def test_mixing_energy_balance_and_limits():
    T_out, S = impl.mix_streams([0.5, 0.5], [270.0, 297.0], [1.05e5, 1.05e5], 1.05e5)
    assert 0.5 * impl.air_h(270.0) + 0.5 * impl.air_h(297.0) == pytest.approx(impl.air_h(T_out), rel=1e-12)
    assert 270.0 < T_out < 297.0 and S > 0
    T_eq, S_eq = impl.mix_streams([0.3, 0.7], [400.0, 400.0], [2e5, 2e5], 2e5)
    assert T_eq == pytest.approx(400.0) and abs(S_eq) < 1e-9
    # mixing of different gases: -R_U sum n ln x; two equal halves give R_U ln 2 per mole
    assert impl.mix_gases_sgen({"a": 0.5, "b": 0.5}) == pytest.approx(impl.R_U * np.log(2.0))


def test_liquid_line_friction_matches_the_small_dp_formula():
    mdot, dp, rho, c, T = 2.0, 1.5e5, 800.0, 2000.0, 300.0
    T_out, S = impl.liquid_line_friction(mdot, dp, rho, c, T)
    assert c * (T_out - T) == pytest.approx(dp / rho)                         # h = c T + p / rho constant
    assert S == pytest.approx(mdot * dp / (rho * T), rel=1e-3)
    assert T * S < mdot * dp / rho                                            # lost work at T below the pumping power


# ---------------------------------------------------------------- the ledger

def test_ledger_ranks_totals_and_refuses_negative_entries():
    L = impl.EntropyLedger(T0=300.0)
    L.add("a", 2.0)
    L.add("b", 5.0)
    L.add("c", 0.0)
    assert L.total() == pytest.approx(7.0)
    assert L.lost_work() == pytest.approx(2100.0)
    r = L.ranked()
    assert [e[0] for e in r] == ["b", "a", "c"] and sum(e[2] for e in r) == pytest.approx(1.0)
    with pytest.raises(ValueError):
        L.add("perpetual motion", -1.0)
    assert L.total() == pytest.approx(7.0)                                    # the refused entry is not recorded


def test_bleed_air_ledger_is_complete_and_positive():
    L = impl.bleed_air_ledger()
    assert len(L.entries) == 4 and all(s > 0 for _, s in L.entries)
    c = impl.BLEED
    # the throttle entry equals mdot R ln(p1/p2); the throttle dominates in this illustrative case
    assert dict(L.entries)["pressure-regulating valve (throttling)"] == pytest.approx(c["mdot"] * R * np.log(c["p_bleed"] / c["p_prv"]))
    assert L.ranked()[0][0].startswith("pressure-regulating valve")


# ---------------------------------------------------------------- cycles and the Carnot bound

def test_clausius_integral_of_a_carnot_cycle_is_zero():
    TH, TL, QH = 900.0, 300.0, 1000.0
    QL = -QH * TL / TH
    assert impl.clausius_integral([QH, QL], [TH, TL]) == pytest.approx(0.0, abs=1e-12)
    r = impl.check_cycle([QH, QL], [TH, TL])
    assert r["eta"] == pytest.approx(impl.carnot_efficiency(TH, TL)) and abs(r["S_gen"]) < 1e-12


def test_every_accepted_cycle_obeys_the_carnot_bound():
    rng = np.random.default_rng(7)
    accepted = rejected = 0
    for _ in range(2000):
        TH = rng.uniform(400.0, 2000.0)
        TL = rng.uniform(150.0, TH - 10.0)
        QH = rng.uniform(1.0, 1000.0)
        eta = rng.uniform(0.0, 0.99)
        Q, T = [QH, -(1 - eta) * QH], [TH, TL]
        if eta > impl.carnot_efficiency(TH, TL) * (1 + 1e-9):
            with pytest.raises(ValueError):
                impl.check_cycle(Q, T)
            rejected += 1
        elif eta < impl.carnot_efficiency(TH, TL) * (1 - 1e-9):
            r = impl.check_cycle(Q, T)
            assert r["eta"] <= r["eta_carnot"] and r["S_gen"] >= 0
            accepted += 1
    assert accepted > 100 and rejected > 100
    # a cycle with several reservoirs is still bounded by Carnot between its extremes
    r = impl.check_cycle([400.0, 300.0, -350.0, -200.0], [1200.0, 900.0, 400.0, 300.0])
    assert r["eta"] <= r["eta_carnot"]


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "entropy_generation.png").exists()
