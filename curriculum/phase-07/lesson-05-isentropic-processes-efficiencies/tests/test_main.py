"""Public tests for P07.05 Isentropic Processes & Efficiencies. They assert thermodynamics, not magic
numbers: the variable-cp isentropic relation reduces to T1 pi^((gamma-1)/gamma) for a constant-cp gas;
the property inversions are exact inverses; eta = 1 generates no entropy and eta < 1 always generates
some; isentropic and polytropic efficiencies convert into each other, by the closed form for constant
cp and by the solver for NASA air; many small stages at a stage efficiency reproduce the polytropic
closed form at first order; and the first law closes three ways (h2 - h1, the integral of cp dT, and
the integral of v dp + T ds along the path).

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

AIR = impl.AIR
COLD = impl.constant_cp_gas(1.4)
T1 = 288.15


# ---------------------------------------------------------------- property inversions

def test_temperature_inversions_are_exact_inverses():
    T = np.array([250.0, 700.0, 999.0, 1001.0, 1800.0])
    npt.assert_allclose(impl.temperature_from_s0(AIR, AIR.s0(T)), T, rtol=1e-11)
    npt.assert_allclose(impl.temperature_from_h(AIR, AIR.h(T)), T, rtol=1e-11)
    assert isinstance(impl.temperature_from_h(AIR, float(AIR.h(500.0))), float)


def test_variable_cp_isentrope_reduces_to_constant_cp_closed_form():
    # the same s0-based solver, given a constant-cp gas, must give T1 pi^((gamma-1)/gamma)
    for gamma in (1.4, 5.0 / 3.0, 1.3):
        gas = impl.constant_cp_gas(gamma)
        for Ta, pi in ((T1, 1.2), (T1, 4.0), (T1, 40.0), (1400.0, 1.0 / 12.0)):
            assert impl.isentropic_temperature(gas, Ta, pi) == pytest.approx(
                impl.isentropic_temperature_closed_form(Ta, pi, gamma), rel=1e-11)


def test_nasa_air_isentrope_is_isentropic_and_runs_cooler_than_gamma_1_4():
    for pi in (2.0, 10.0, 40.0):
        T2s = impl.isentropic_temperature(AIR, T1, pi)
        assert impl.entropy_generation(AIR, T1, T2s, pi) == pytest.approx(0.0, abs=1e-9)
        # cp rises with T, so gamma falls below 1.4 and the real isentrope is cooler
        assert T2s < impl.isentropic_temperature_closed_form(T1, pi, 1.4)
    # at a small pressure ratio near room temperature, gamma = 1.4 is almost exact
    assert impl.isentropic_temperature(AIR, T1, 1.2) == pytest.approx(
        impl.isentropic_temperature_closed_form(T1, 1.2, 1.4), abs=0.02)


# ---------------------------------------------------------------- efficiencies and the second law

@pytest.mark.parametrize("gas", [AIR, COLD], ids=["nasa", "constant-cp"])
def test_eta_one_generates_no_entropy_and_eta_below_one_always_does(gas):
    pi = 12.0
    assert impl.compressor_exit_temperature(gas, T1, pi, 1.0) == pytest.approx(impl.isentropic_temperature(gas, T1, pi), rel=1e-12)
    assert impl.entropy_generation(gas, T1, impl.compressor_exit_temperature(gas, T1, pi, 1.0), pi) == pytest.approx(0.0, abs=1e-9)
    assert impl.entropy_generation(gas, 1400.0, impl.turbine_exit_temperature(gas, 1400.0, pi, 1.0), 1 / pi) == pytest.approx(0.0, abs=1e-9)
    sgen_c = [impl.entropy_generation(gas, T1, impl.compressor_exit_temperature(gas, T1, pi, e), pi) for e in (0.95, 0.9, 0.8)]
    sgen_t = [impl.entropy_generation(gas, 1400.0, impl.turbine_exit_temperature(gas, 1400.0, pi, e), 1 / pi) for e in (0.95, 0.9, 0.8)]
    assert np.all(np.array(sgen_c) > 0) and np.all(np.diff(sgen_c) > 0)
    assert np.all(np.array(sgen_t) > 0) and np.all(np.diff(sgen_t) > 0)


def test_efficiency_definitions_round_trip():
    pi = 25.0
    T2 = impl.compressor_exit_temperature(AIR, T1, pi, 0.86)
    assert impl.isentropic_efficiency(AIR, T1, T2, pi) == pytest.approx(0.86, rel=1e-10)
    T2p = impl.polytropic_exit_temperature(AIR, T1, pi, 0.91)
    assert impl.polytropic_efficiency(AIR, T1, T2p, pi) == pytest.approx(0.91, rel=1e-10)
    Tt = impl.turbine_exit_temperature(AIR, 1500.0, 8.0, 0.9)
    assert impl.isentropic_efficiency(AIR, 1500.0, Tt, 1 / 8.0) == pytest.approx(0.9, rel=1e-10)
    Ttp = impl.polytropic_exit_temperature(AIR, 1500.0, 1 / 8.0, 0.88)
    assert impl.polytropic_efficiency(AIR, 1500.0, Ttp, 1 / 8.0) == pytest.approx(0.88, rel=1e-10)


def test_polytropic_closed_form_and_the_compressor_turbine_asymmetry():
    for pi in (2.0, 10.0, 40.0):
        for eta_p in (0.85, 0.9):
            T2 = impl.polytropic_exit_temperature(COLD, T1, pi, eta_p)
            assert T2 == pytest.approx(T1 * pi ** (0.4 / 1.4 / eta_p), rel=1e-11)
            eta_s = impl.isentropic_efficiency(COLD, T1, T2, pi)
            assert eta_s == pytest.approx(impl.isentropic_from_polytropic(pi, eta_p), rel=1e-10)
            assert eta_s < eta_p                                                  # compressor: eta_s below eta_p
            Tt = impl.polytropic_exit_temperature(COLD, 1400.0, 1 / pi, eta_p)
            eta_t = impl.isentropic_efficiency(COLD, 1400.0, Tt, 1 / pi)
            assert eta_t == pytest.approx(impl.isentropic_from_polytropic_turbine(pi, eta_p), rel=1e-10)
            assert eta_t > eta_p                                                  # turbine: eta_s above eta_p
    # the two efficiencies coincide as the pressure ratio goes to 1
    assert impl.isentropic_from_polytropic(1.0001, 0.9) == pytest.approx(0.9, abs=1e-4)
    assert impl.isentropic_from_polytropic_turbine(1.0001, 0.9) == pytest.approx(0.9, abs=1e-4)
    # and eta_s falls as the pressure ratio rises at fixed eta_p
    assert np.all(np.diff([impl.isentropic_from_polytropic(p, 0.9) for p in (2, 5, 10, 20, 40)]) < 0)


def test_polytropic_entropy_generation_is_independent_of_cp():
    # along a polytropic path s_gen = R ln(pi) (1/eta_p - 1) for any ideal gas
    pi, eta_p = 40.0, 0.9
    expected = impl.R_AIR * np.log(pi) * (1 / eta_p - 1)
    for gas in (AIR, COLD):
        T2 = impl.polytropic_exit_temperature(gas, T1, pi, eta_p)
        assert impl.entropy_generation(gas, T1, T2, pi) == pytest.approx(expected, rel=1e-9)


# ---------------------------------------------------------------- stages

@pytest.mark.parametrize("gas", [AIR, COLD], ids=["nasa", "constant-cp"])
def test_many_small_stages_reproduce_the_polytropic_closed_form(gas):
    pi, eta = 40.0, 0.9
    T_lim = impl.polytropic_exit_temperature(gas, T1, pi, eta)
    errs = []
    for n in (8, 16, 32, 64):
        T = impl.multistage_compressor(gas, T1, pi, n, eta)
        assert len(T) == n + 1 and T[0] == T1 and np.all(np.diff(T) > 0)
        errs.append(abs(T[-1] - T_lim))
    orders = impl.observed_order(errs)
    assert np.all(orders > 0.85) and np.all(orders < 1.15)                    # first order in 1/n
    T = impl.multistage_compressor(gas, T1, pi, 256, eta)
    assert impl.polytropic_efficiency(gas, T1, T[-1], pi) == pytest.approx(eta, abs=5e-4)
    # one stage is the whole compressor at its isentropic efficiency
    assert impl.multistage_compressor(gas, T1, pi, 1, eta)[-1] == pytest.approx(
        impl.compressor_exit_temperature(gas, T1, pi, eta), rel=1e-12)


def test_multistage_overall_efficiency_is_below_stage_efficiency():
    # the preheat effect: later stages work on hotter air, so the same stage eta_s costs more work
    pi = 30.0
    T = impl.multistage_compressor(AIR, T1, pi, 10, 0.9)
    assert impl.isentropic_efficiency(AIR, T1, T[-1], pi) < 0.9
    pi_st = pi ** 0.1
    w_stage = [impl.specific_work(AIR, T[i], T[i + 1]) for i in range(10)]
    assert np.all(np.diff(w_stage) > 0)
    # stage works add to the overall work (first law, stage by stage)
    assert sum(w_stage) == pytest.approx(impl.specific_work(AIR, T[0], T[-1]), rel=1e-12)
    for i in range(10):
        assert impl.isentropic_efficiency(AIR, T[i], T[i + 1], pi_st) == pytest.approx(0.9, rel=1e-9)


# ---------------------------------------------------------------- first law

@pytest.mark.parametrize("gas", [AIR, COLD], ids=["nasa", "constant-cp"])
def test_first_law_closes_three_ways(gas):
    dh, cpdT, vdp, Tds = impl.first_law_integrals(gas, T1, 40.0, 0.9)
    assert cpdT == pytest.approx(dh, rel=1e-9)
    assert vdp + Tds == pytest.approx(dh, rel=1e-5)            # trapezoid on 2001 points along the path
    assert vdp > 0 and Tds > 0                                   # friction heats the gas: T ds > 0
    # the friction heat is the fraction (1 - eta_p) of the work, exactly, along a polytropic path
    assert Tds / dh == pytest.approx(0.1, rel=1e-5)
    # a reversible path: T ds = 0 and all the work is v dp
    dh1, _, vdp1, Tds1 = impl.first_law_integrals(gas, T1, 40.0, 1.0)
    assert abs(Tds1) < 1e-6 * dh1 and vdp1 == pytest.approx(dh1, rel=1e-5)


def test_work_and_dimensional_sanity():
    w = impl.specific_work(AIR, T1, impl.polytropic_exit_temperature(AIR, T1, 40.0, 0.9))
    assert 5e5 < w < 8e5                                          # J/kg: several hundred kJ/kg for pi = 40
    assert impl.specific_work(AIR, 1400.0, impl.turbine_exit_temperature(AIR, 1400.0, 12.0, 0.9)) < 0
    assert 1.39 < float(AIR.cp(300.0)) / (float(AIR.cp(300.0)) - impl.R_AIR) < 1.41


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "compressor_efficiency.png").stat().st_size > 10_000
