"""Public tests for P07.06 Power Cycles: Brayton. They assert thermodynamics, not magic numbers: the
state table is isentropic and isobaric where it should be; the first law closes for the cycle; the
constant-cp efficiency and work match their closed forms; the net work equals the enclosed area in
both the p-v and the T-s plane by a numerical contour integral, converging at second order; the
efficiency is below Carnot between the same extreme temperatures; and the pressure ratio of maximum
work matches (T3/T1)^(gamma/(2(gamma-1))), with T2 = T4 there for any ideal gas.

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
T1, P1, PI, T3 = 288.15, 101325.0, 12.0, 1400.0
GASES = pytest.mark.parametrize("gas", [AIR, COLD], ids=["nasa", "constant-cp"])


# ---------------------------------------------------------------- state table

@GASES
def test_state_table_has_the_four_processes(gas):
    st = impl.brayton_states(gas, T1, P1, PI, T3)
    T, p, v, s = st["T"], st["p"], st["v"], st["s"]
    assert T[0] == T1 and T[2] == T3
    npt.assert_allclose(p, [P1, PI * P1, PI * P1, P1], rtol=1e-14)
    npt.assert_allclose(p * v, gas.R * T, rtol=1e-14)                    # ideal-gas equation of state
    assert s[1] == pytest.approx(s[0], abs=1e-8) and s[3] == pytest.approx(s[2], abs=1e-8)   # isentropic legs
    assert T[0] < T[1] < T[2] and T[3] > T[0] and T[3] < T[2]
    npt.assert_allclose(st["h"], gas.h(T), rtol=1e-14)


@GASES
def test_first_law_closes_for_the_cycle(gas):
    pf = impl.brayton_performance(gas, T1, P1, PI, T3)
    assert pf["q_in"] - pf["q_out"] == pytest.approx(pf["w_net"], rel=1e-12)
    assert pf["w_t"] - pf["w_c"] == pytest.approx(pf["w_net"], rel=1e-12)
    assert pf["eta"] == pytest.approx(pf["w_net"] / pf["q_in"], rel=1e-14)
    assert pf["bwr"] == pytest.approx(pf["w_c"] / pf["w_t"], rel=1e-14)
    assert 0 < pf["bwr"] < 1 and pf["q_out"] > 0


# ---------------------------------------------------------------- closed forms

def test_constant_cp_efficiency_and_work_match_closed_forms():
    for pi in (2.0, 5.0, 12.0, 30.0):
        for T3s in (1100.0, 1400.0, 1700.0):
            pf = impl.brayton_performance(COLD, T1, P1, pi, T3s)
            assert pf["eta"] == pytest.approx(float(impl.efficiency_closed_form(pi, 1.4)), rel=1e-11)
            assert pf["w_net"] == pytest.approx(float(impl.net_work_closed_form(T1, T3s, pi, 1.4)), rel=1e-10)
    # eta depends only on the pressure ratio: same for every T3; 1 - T1/T2
    st = impl.brayton_states(COLD, T1, P1, PI, T3)
    assert float(impl.efficiency_closed_form(PI)) == pytest.approx(1 - st["T"][0] / st["T"][1], rel=1e-12)
    # monatomic gas: larger gamma, higher efficiency at the same pi
    assert impl.efficiency_closed_form(PI, 5 / 3) > impl.efficiency_closed_form(PI, 1.4)


def test_work_limits_zero_at_pi_one_and_at_pi_max():
    assert float(impl.net_work_closed_form(T1, T3, 1.0)) == pytest.approx(0.0, abs=1e-9)
    pi_max = (T3 / T1) ** 3.5                     # T2s = T3: no room to add heat
    assert float(impl.net_work_closed_form(T1, T3, pi_max)) == pytest.approx(0.0, abs=1e-6)
    pf = impl.brayton_performance(COLD, T1, P1, pi_max * (1 - 1e-9), T3)
    assert pf["eta"] == pytest.approx(impl.carnot_efficiency(T1, T3), rel=1e-6)   # Carnot reached with no work
    assert abs(pf["w_net"]) < 1e-3


# ---------------------------------------------------------------- enclosed area

@GASES
def test_net_work_equals_enclosed_area_in_p_v_and_T_s(gas):
    wn = impl.brayton_performance(gas, T1, P1, PI, T3)["w_net"]
    path = impl.cycle_path(gas, T1, P1, PI, T3, 400)
    assert len(path["p"]) == 4 * 400 + 1
    assert path["p"][0] == path["p"][-1] and path["T"][0] == path["T"][-1]
    assert impl.enclosed_area(path["v"], path["p"]) == pytest.approx(wn, rel=1e-5)
    assert impl.enclosed_area(path["s"], path["T"]) == pytest.approx(wn, rel=1e-5)


def test_enclosed_area_is_a_signed_contour_integral():
    # unit square traversed clockwise in (x, y) has area +1; counter-clockwise -1
    x = np.array([0.0, 0.0, 1.0, 1.0, 0.0])
    y = np.array([0.0, 1.0, 1.0, 0.0, 0.0])
    assert impl.enclosed_area(x, y) == pytest.approx(1.0)
    assert impl.enclosed_area(x[::-1], y[::-1]) == pytest.approx(-1.0)
    # a circle of radius 2 by many points, clockwise
    th = np.linspace(0, 2 * np.pi, 4001)
    assert impl.enclosed_area(2 * np.sin(th), 2 * np.cos(th)) == pytest.approx(4 * np.pi, rel=1e-5)


def test_contour_integrals_converge_at_second_order():
    wn = impl.brayton_performance(AIR, T1, P1, PI, T3)["w_net"]
    e_pv, e_ts = [], []
    for n in (25, 50, 100, 200):
        path = impl.cycle_path(AIR, T1, P1, PI, T3, n)
        e_pv.append(abs(impl.enclosed_area(path["v"], path["p"]) - wn))
        e_ts.append(abs(impl.enclosed_area(path["s"], path["T"]) - wn))
    for e in (e_pv, e_ts):
        o = impl.observed_order(e)
        assert np.all(o > 1.8) and np.all(o < 2.2)


def test_isentropic_legs_contribute_nothing_to_T_ds():
    path = impl.cycle_path(AIR, T1, P1, PI, T3, 100)
    s = path["s"]
    npt.assert_allclose(s[:101], s[0], atol=1e-8)                 # 1-2 isentropic
    npt.assert_allclose(s[200:301], s[200], atol=1e-8)            # 3-4 isentropic
    npt.assert_allclose(path["p"][100:201], PI * P1, rtol=1e-14)  # 2-3 isobaric


# ---------------------------------------------------------------- Carnot and the optimum

@GASES
def test_efficiency_below_carnot(gas):
    for pi in (2.0, 8.0, 20.0, 50.0):
        for T3s in (1100.0, 1700.0):
            pf = impl.brayton_performance(gas, T1, P1, pi, T3s)
            assert 0 < pf["eta"] < pf["eta_carnot"]
            assert pf["eta_carnot"] == pytest.approx(impl.carnot_efficiency(T1, T3s))
    assert impl.carnot_efficiency(300.0, 1200.0) == pytest.approx(0.75)


def test_optimum_pressure_ratio_matches_closed_form():
    for T3s in (1100.0, 1400.0, 1700.0):
        pio = impl.optimum_pressure_ratio(COLD, T1, P1, T3s)
        assert pio == pytest.approx(impl.optimum_pressure_ratio_closed_form(T1, T3s, 1.4), rel=1e-5)
        st = impl.brayton_states(COLD, T1, P1, pio, T3s)
        assert st["T"][1] == pytest.approx(np.sqrt(T1 * T3s), rel=1e-5)          # T2 = sqrt(T1 T3)


@GASES
def test_at_maximum_work_T2_equals_T4_for_any_ideal_gas(gas):
    # dw_net/d(ln pi) = R (T4 - T2) for an ideal gas, so the maximum is where T2 = T4
    pio = impl.optimum_pressure_ratio(gas, T1, P1, T3)
    st = impl.brayton_states(gas, T1, P1, pio, T3)
    assert st["T"][1] == pytest.approx(st["T"][3], rel=1e-5)
    w = [impl.brayton_performance(gas, T1, P1, pio * f, T3)["w_net"] for f in (0.9, 1.0, 1.1)]
    assert w[1] > w[0] and w[1] > w[2]


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "brayton_cycle.png").stat().st_size > 10_000
