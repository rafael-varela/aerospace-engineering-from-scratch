"""Public tests for P07.07 Real Brayton Cycles. They assert thermodynamics, not magic numbers: with no
recuperator, no intercooler and ideal components the model is the P07.06 cycle; with an ideal
recuperator it is the closed form 1 - (T1/T3) pi^((gamma-1)/gamma), and with an intercooler too, the
intercooled-recuperated closed form; the first law closes for every component and for the cycle; the
entropy generated in every component is non-negative and sums to the Clausius integral; efficiency
stays below Carnot; component losses only ever lower efficiency; the recuperator helps only while the
turbine exit is hotter than the compressor exit; and sqrt(pi) is the best intercooler split.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import pytest

from conftest import load_impl

impl = load_impl(__file__)

AIR = impl.AIR
COLD = impl.constant_cp_gas(1.4)
T1, T3 = 288.15, 1400.0


def spec(**kw):
    return impl.CycleSpec(**kw)


REAL = dict(eta_c=0.85, eta_t=0.88, dp_cc=0.04)
ICR = dict(intercool=True, dp_ic=0.02, eps=0.85, dp_rc=0.03, dp_rh=0.04)


# ---------------------------------------------------------------- ideal limits

@pytest.mark.parametrize("pi", [4.0, 12.0, 30.0])
def test_no_recuperator_ideal_components_recover_lesson_06(pi):
    pf = impl.cycle_performance(COLD, spec(pi=pi, T_tit=T3))
    assert pf["eta"] == pytest.approx(impl.simple_cycle_efficiency(pi), rel=1e-11)
    x = pi ** (0.4 / 1.4)
    cp = 1.4 * impl.R_AIR / 0.4
    assert pf["w_net"] == pytest.approx(cp * T1 * ((T3 / T1) * (1 - 1 / x) - (x - 1)), rel=1e-10)


@pytest.mark.parametrize("pi", [3.0, 8.0, 12.0])
def test_ideal_recuperator_recovers_the_closed_form(pi):
    pf = impl.cycle_performance(COLD, spec(pi=pi, T_tit=T3, eps=1.0))
    assert pf["eta"] == pytest.approx(impl.ideal_recuperated_efficiency(pi, T1, T3), rel=1e-10)
    st = impl.cycle_states(COLD, spec(pi=pi, T_tit=T3, eps=1.0))
    assert st["T"][4] == pytest.approx(st["T"][6], rel=1e-12)     # cold side leaves at the turbine exit T
    assert st["T"][7] == pytest.approx(st["T"][3], rel=1e-12)     # exhaust leaves at the compressor exit T
    # with an intercooler as well
    pf2 = impl.cycle_performance(COLD, spec(pi=pi, T_tit=T3, eps=1.0, intercool=True))
    assert pf2["eta"] == pytest.approx(impl.ideal_icr_efficiency(pi, T1, T3), rel=1e-10)
    assert pf2["eta"] > pf["eta"]


def test_recuperation_helps_only_below_the_crossover_pressure_ratio():
    # ideal: the recuperated and simple efficiencies cross where T2 = T4, pi = (T3/T1)^(gamma/(2(gamma-1)))
    pi_x = (T3 / T1) ** 1.75
    assert impl.ideal_recuperated_efficiency(pi_x, T1, T3) == pytest.approx(impl.simple_cycle_efficiency(pi_x), rel=1e-12)
    for pi, better in ((8.0, True), (25.0, False)):
        e_rec = impl.cycle_performance(COLD, spec(pi=pi, T_tit=T3, eps=0.8))["eta"]
        e_sim = impl.cycle_performance(COLD, spec(pi=pi, T_tit=T3))["eta"]
        assert bool(e_rec > e_sim) is better
        st = impl.cycle_states(COLD, spec(pi=pi, T_tit=T3, eps=0.8))
        assert bool(st["T"][6] > st["T"][3]) is better                # turbine exit hotter than compressor exit


# ---------------------------------------------------------------- conservation

@pytest.mark.parametrize("gas", [AIR, COLD], ids=["nasa", "constant-cp"])
def test_first_law_closes_per_component_and_overall(gas):
    sp = spec(pi=12.0, T_tit=T3, **REAL, **ICR)
    bal = impl.component_balances(gas, sp)
    assert set(bal) >= {"compressors", "intercooler", "recuperator", "combustor", "turbine", "exhaust", "cycle"}
    for k, v in bal.items():
        assert abs(v) < 1e-9 * 1e6, k                                 # J/kg against energies of order 1e6
    pf = impl.cycle_performance(gas, sp)
    assert pf["q_in"] - pf["q_ic"] - pf["q_ex"] == pytest.approx(pf["w_net"], rel=1e-12)
    # the recuperator moves heat from the exhaust to the compressed air; it adds none
    st = impl.cycle_states(gas, sp)
    assert pf["q_rec"] == pytest.approx(st["h"][6] - st["h"][7], rel=1e-12) and pf["q_rec"] > 0


@pytest.mark.parametrize("gas", [AIR, COLD], ids=["nasa", "constant-cp"])
def test_entropy_generation_non_negative_and_equals_clausius(gas):
    sp = spec(pi=12.0, T_tit=T3, **REAL, **ICR)
    sg = impl.entropy_generation(gas, sp)
    for k, v in sg.items():
        assert v >= -1e-9, k
    pf = impl.cycle_performance(gas, sp)
    assert sg["total"] == pytest.approx((pf["q_ic"] + pf["q_ex"]) / T1 - pf["q_in"] / T3, rel=1e-10)
    # an ideal simple cycle generates entropy only where heat crosses a finite temperature difference
    sg0 = impl.entropy_generation(gas, spec(pi=12.0, T_tit=T3))
    assert abs(sg0["LP compressor"]) < 1e-9 and abs(sg0["turbine"]) < 1e-9
    assert sg0["combustor"] > 0 and sg0["exhaust"] > 0


def test_efficiency_below_carnot_for_every_configuration():
    for gas in (AIR, COLD):
        for extra in ({}, dict(intercool=True), dict(eps=0.9), dict(intercool=True, eps=1.0)):
            for pi in (3.0, 12.0, 30.0):
                pf = impl.cycle_performance(gas, spec(pi=pi, T_tit=T3, **extra))
                assert 0 < pf["eta"] < pf["eta_carnot"]
                assert pf["eta_carnot"] == pytest.approx(1 - T1 / T3)


# ---------------------------------------------------------------- losses and design choices

def test_every_loss_lowers_efficiency():
    base = dict(pi=12.0, T_tit=T3, **REAL, **ICR)
    e0 = impl.cycle_performance(AIR, spec(**base))["eta"]
    for change in ({"eta_c": 0.80}, {"eta_t": 0.83}, {"eps": 0.70}, {"dp_cc": 0.08}, {"dp_rc": 0.06},
                   {"dp_rh": 0.08}, {"dp_ic": 0.05}, {"dp_ex": 0.03}, {"T_ic": 320.0}):
        e = impl.cycle_performance(AIR, spec(**{**base, **change}))["eta"]
        assert e < e0, change
    # real components: lower efficiency and work than the ideal cycle
    ideal = impl.cycle_performance(AIR, spec(pi=12.0, T_tit=T3))
    real = impl.cycle_performance(AIR, spec(pi=12.0, T_tit=T3, **REAL))
    assert real["eta"] < ideal["eta"] and real["w_net"] < ideal["w_net"]


def test_intercooling_raises_work_and_needs_a_recuperator_to_raise_efficiency():
    a = impl.cycle_performance(AIR, spec(pi=12.0, T_tit=T3, **REAL))
    b = impl.cycle_performance(AIR, spec(pi=12.0, T_tit=T3, **REAL, intercool=True, dp_ic=0.02))
    c = impl.cycle_performance(AIR, spec(pi=12.0, T_tit=T3, **REAL, **ICR))
    assert b["w_net"] > a["w_net"] and b["w_c"] < a["w_c"]
    assert b["eta"] < a["eta"]            # the extra heat to raise colder air to TIT costs more than the work saved
    assert c["eta"] > a["eta"]            # with a recuperator the exhaust supplies that heat


def test_sqrt_pi_is_the_best_intercooler_split_for_constant_cp():
    for pi in (9.0, 16.0):
        sp = spec(pi=pi, T_tit=T3, eta_c=0.85, intercool=True)
        assert impl.best_split(COLD, sp) == pytest.approx(np.sqrt(pi), rel=1e-5)
    # with pressure loss in the intercooler the best split moves above sqrt(pi)
    assert impl.best_split(COLD, spec(pi=16.0, T_tit=T3, eta_c=0.85, intercool=True, dp_ic=0.05)) > 4.0


def test_station_pressures_follow_the_losses():
    sp = spec(pi=12.0, T_tit=T3, **REAL, **ICR, dp_ex=0.01)
    p = impl.cycle_states(AIR, sp)["p"]
    assert p[3] == pytest.approx(12.0 * sp.p1)
    assert p[2] == pytest.approx(p[1] * 0.98) and p[4] == pytest.approx(p[3] * 0.97)
    assert p[5] == pytest.approx(p[4] * 0.96) and p[7] == pytest.approx(sp.p1 / 0.99)
    assert p[7] == pytest.approx(p[6] * 0.96)
    assert p[1] == pytest.approx(np.sqrt(12.0) * sp.p1)


def test_sweep_shapes():
    pis = np.array([3.0, 5.0, 8.0, 12.0, 20.0, 30.0])
    e_s, w_s = impl.sweep(AIR, spec(T_tit=T3, **REAL), pis)
    e_r, _ = impl.sweep(AIR, spec(T_tit=T3, **REAL, eps=0.85, dp_rc=0.03, dp_rh=0.04), pis)
    assert np.argmax(e_r) < np.argmax(e_s)          # recuperated cycle peaks at a lower pressure ratio
    assert 0 < np.argmax(w_s) < len(pis) - 1        # work has an interior maximum


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "real_brayton.png").stat().st_size > 10_000
