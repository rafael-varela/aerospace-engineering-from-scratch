"""Public tests for P07.02 The First Law. They assert conservation and closed forms, not magic numbers:
the orifice flow is continuous at the choking ratio and scales with p0/sqrt(T0); the control-volume
energy balance U + int(mdot h) - int(Qdot) = U0 closes to round-off in every heat mode; the adiabatic tank state stays on p ~ rho^gamma and its internal-energy drop equals the closed-system
boundary work; the ODE matches the closed-form choked adiabatic blowdown and the isothermal exponential;
a finite wall conductance lies between the two limits; and the blowdown ends at the back pressure.

This is a COMPUTATIONAL model only: no pressure vessel is built, filled or vented at home in this course.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. The only file written goes to tmp_path.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

TANK = impl.Tank(V=0.02, p0=20e6, T0=290.0, gas="N2")
ORF = impl.Orifice(Cd=0.75, A=np.pi * 0.0015 ** 2 / 4)
PB = 0.2e6


# ---------------------------------------------------------------- orifice

@pytest.mark.parametrize("gamma, expected", [(1.4, 0.528282), (5 / 3, 0.487139)])
def test_critical_pressure_ratio(gamma, expected):
    assert impl.critical_pressure_ratio(gamma) == pytest.approx(expected, abs=1e-6)


def test_mass_flow_is_continuous_at_choking_and_has_the_right_scaling():
    R, g, cv, cp = impl.gas_properties("N2")
    p0, T0, Cd, A = 1e6, 300.0, 0.8, 1e-6
    pc = p0 * impl.critical_pressure_ratio(g)
    just_choked = impl.mass_flow(p0, T0, pc * (1 - 1e-9), Cd, A, R, g)
    just_sub = impl.mass_flow(p0, T0, pc * (1 + 1e-9), Cd, A, R, g)
    assert just_sub == pytest.approx(just_choked, rel=1e-6)
    # choked flow is independent of the back pressure and proportional to p0 / sqrt(T0)
    assert impl.mass_flow(p0, T0, 0.1 * pc, Cd, A, R, g) == pytest.approx(just_choked, rel=1e-12)
    assert impl.mass_flow(2 * p0, 4 * T0, 0.0, Cd, A, R, g) == pytest.approx(impl.mass_flow(p0, T0, 0.0, Cd, A, R, g))
    # subsonic flow falls to zero as the back pressure rises to p0, and there is no reverse flow
    sub = [impl.mass_flow(p0, T0, r * p0, Cd, A, R, g) for r in (0.6, 0.8, 0.95, 0.999)]
    assert np.all(np.diff(sub) < 0) and sub[-1] < 0.1 * just_choked
    assert impl.mass_flow(p0, T0, p0, Cd, A, R, g) == 0.0
    assert impl.mass_flow(p0, T0, 2 * p0, Cd, A, R, g) == 0.0


# ---------------------------------------------------------------- conservation

@pytest.mark.parametrize("mode, hA", [("adiabatic", 0.0), ("isothermal", 0.0), ("wall", 20.0)])
def test_energy_and_mass_balances_close_to_round_off(mode, hA):
    res = impl.simulate_blowdown(TANK, ORF, PB, mode, t_end=400.0, hA=hA)
    assert impl.energy_residual(res) < 1e-12
    m0 = impl.initial_mass(TANK)
    # mass only leaves: m falls monotonically and stays positive
    assert np.all(np.diff(res["m"]) <= 1e-12 * m0) and res["m"][-1] > 0
    # the energy that left as enthalpy is positive and bounded by U0 + heat in
    assert 0 < res["H_out"][-1] <= res["U"][0] + res["Q_in"][-1] + 1e-9 * res["U"][0]


def test_isothermal_mode_holds_t_and_its_heat_is_r_t_times_the_mass_out():
    res = impl.simulate_blowdown(TANK, ORF, PB, "isothermal", t_end=300.0)
    npt.assert_allclose(res["T"], TANK.T0, rtol=1e-9)
    R = impl.gas_properties("N2")[0]
    m0 = impl.initial_mass(TANK)
    assert res["Q_in"][-1] == pytest.approx(R * TANK.T0 * (m0 - res["m"][-1]), rel=1e-8)


# ---------------------------------------------------------------- the adiabatic tank and closed-system work

def test_adiabatic_tank_state_follows_p_proportional_to_rho_gamma():
    res = impl.simulate_blowdown(TANK, ORF, PB, "adiabatic", t_end=800.0)
    g = impl.gas_properties("N2")[1]
    rho = res["m"] / TANK.V
    k = res["p"] / rho ** g
    assert np.ptp(k) / k[0] < 1e-8
    assert res["T"][-1] < 0.5 * TANK.T0                                   # strong cooling


def test_internal_energy_drop_of_the_remaining_gas_equals_its_boundary_work():
    res = impl.simulate_blowdown(TANK, ORF, PB, "adiabatic", t_end=150.0)
    R, g, cv, cp = impl.gas_properties("N2")
    m0 = impl.initial_mass(TANK)
    w = impl.polytropic_work(TANK.p0, TANK.V / m0, res["p"][-1], TANK.V / res["m"][-1], g)
    assert w == pytest.approx(cv * (TANK.T0 - res["T"][-1]), rel=1e-8)    # q = 0: u1 - u2 = w
    # the isothermal formula is the n -> 1 limit of the polytropic one
    p1, v1, v2 = 1e6, 0.1, 0.3
    w1 = impl.polytropic_work(p1, v1, p1 * v1 / v2, v2, 1.0)
    assert w1 == pytest.approx(p1 * v1 * np.log(3.0))
    n = 1.0 + 1e-6
    assert impl.polytropic_work(p1, v1, p1 * (v1 / v2) ** n, v2, n) == pytest.approx(w1, rel=1e-5)


# ---------------------------------------------------------------- closed forms

@pytest.mark.parametrize("gas", ["N2", "He"])
def test_ode_matches_the_closed_form_choked_adiabatic_blowdown(gas):
    tank = impl.Tank(V=TANK.V, p0=TANK.p0, T0=TANK.T0, gas=gas)
    g = impl.gas_properties(gas)[1]
    p_unchoke = PB / impl.critical_pressure_ratio(g)
    t_unchoke = impl.adiabatic_time_to_pressure(1.05 * p_unchoke, tank, ORF)
    res = impl.simulate_blowdown(tank, ORF, PB, "adiabatic", t_end=t_unchoke)
    p, T, m = impl.adiabatic_choked_closed_form(res["t"], tank, ORF)
    npt.assert_allclose(res["p"], p, rtol=1e-7)
    npt.assert_allclose(res["T"], T, rtol=1e-7)
    npt.assert_allclose(res["m"], m, rtol=1e-7)
    assert res["p"][-1] == pytest.approx(1.05 * p_unchoke, rel=1e-6)


def test_ode_matches_the_isothermal_exponential_and_the_time_constant():
    tau = impl.blowdown_time_constant(TANK, ORF)
    p_unchoke = PB / impl.critical_pressure_ratio(1.4)
    t_end = tau * np.log(TANK.p0 / (1.05 * p_unchoke))
    res = impl.simulate_blowdown(TANK, ORF, PB, "isothermal", t_end=t_end)
    p, m = impl.isothermal_choked_closed_form(res["t"], TANK, ORF)
    npt.assert_allclose(res["p"], p, rtol=1e-7)
    # tau is the initial mass over the initial mass flow
    R, g, cv, cp = impl.gas_properties("N2")
    mdot0 = impl.mass_flow(TANK.p0, TANK.T0, PB, ORF.Cd, ORF.A, R, g)
    assert tau == pytest.approx(impl.initial_mass(TANK) / mdot0, rel=1e-12)


def test_p_stop_event_and_the_time_to_pressure_inverse():
    p_stop = 4e6
    res = impl.simulate_blowdown(TANK, ORF, PB, "adiabatic", t_end=1e4, p_stop=p_stop)
    assert res["p"][-1] == pytest.approx(p_stop, rel=1e-8)
    assert res["t_stop"] == pytest.approx(impl.adiabatic_time_to_pressure(p_stop, TANK, ORF), rel=1e-7)
    assert impl.adiabatic_time_to_pressure(TANK.p0, TANK, ORF) == pytest.approx(0.0, abs=1e-12)


# ---------------------------------------------------------------- limits and the engineering consequences

def test_wall_heat_transfer_lies_between_the_limits():
    p_stop = 4e6
    T_end = []
    for hA in (0.0, 5.0, 50.0, 1e5):
        r = impl.simulate_blowdown(TANK, ORF, PB, "wall", t_end=1e4, p_stop=p_stop, hA=hA)
        T_end.append(r["T"][-1])
    ad = impl.simulate_blowdown(TANK, ORF, PB, "adiabatic", t_end=1e4, p_stop=p_stop)
    assert T_end[0] == pytest.approx(ad["T"][-1], rel=1e-9)               # hA = 0 is the adiabatic case
    assert np.all(np.diff(T_end) > 0)                                       # more conductance, warmer gas
    assert T_end[-1] == pytest.approx(TANK.T0, rel=2e-3)                   # very large hA tends to isothermal


def test_remaining_mass_and_recovery_pressure():
    g = 1.4
    p_ratio = 0.15
    res = impl.simulate_blowdown(TANK, ORF, PB, "adiabatic", t_end=1e4, p_stop=p_ratio * TANK.p0)
    m0 = impl.initial_mass(TANK)
    assert res["m"][-1] / m0 == pytest.approx(impl.remaining_mass_fraction(p_ratio, g, "adiabatic"), rel=1e-8)
    assert impl.remaining_mass_fraction(p_ratio, g, "adiabatic") > impl.remaining_mass_fraction(p_ratio, g, "isothermal")
    # warming the stranded gas back to T0 at constant volume: p = m R T0 / V
    R = impl.gas_properties("N2")[0]
    assert impl.recovery_pressure(res["p"][-1], res["T"][-1], TANK.T0) == pytest.approx(res["m"][-1] * R * TANK.T0 / TANK.V, rel=1e-8)


def test_long_blowdown_ends_at_the_back_pressure_and_helium_is_faster():
    res = impl.simulate_blowdown(TANK, ORF, PB, "isothermal", t_end=3000.0)
    assert res["p"][-1] == pytest.approx(PB, rel=1e-3)
    assert np.all(res["p"] >= PB * (1 - 1e-9))
    he = impl.Tank(V=TANK.V, p0=TANK.p0, T0=TANK.T0, gas="He")
    assert impl.blowdown_time_constant(he, ORF) < 0.5 * impl.blowdown_time_constant(TANK, ORF)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "tank_blowdown.png").exists()
