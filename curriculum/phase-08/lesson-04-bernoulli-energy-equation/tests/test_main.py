"""Public tests for P08.04 Bernoulli & Energy Equation. They assert physics, not magic numbers:
every airspeed conversion inverts to round-off (TAS -> CAS -> TAS, CAS <-> q_c, venturi dp <-> Q);
the compressible pitot law tends to Bernoulli as M -> 0 with the M^2/4 correction; at sea level
CAS = EAS = TAS; a static-port error moves the altimeter by dp/(rho g), consistent with the
hydrostatic gradient of the atmosphere; the Euler equation integrated through a contraction
converges to Bernoulli at second order; and the energy equation balances with losses and shaft
work in the limits that must hold.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)


# ---------------------------------------------------------------- atmosphere and Bernoulli

def test_troposphere_matches_published_isa_values():
    T, p, rho = impl.isa_troposphere(np.array([0.0, 11000.0]))
    assert rho[0] == pytest.approx(1.225, rel=1e-4)
    assert p[1] == pytest.approx(22632.06, rel=1e-5)                     # published ISA anchor
    assert T[1] == pytest.approx(216.65)
    with pytest.raises(ValueError):
        impl.isa_troposphere(12000.0)


def test_bernoulli_airspeed_inverts_dynamic_pressure_and_scales():
    V = np.array([0.0, 10.0, 60.0])
    npt.assert_allclose(impl.bernoulli_airspeed(impl.dynamic_pressure(0.9091, V), 0.9091), V, rtol=1e-14)
    assert float(impl.bernoulli_airspeed(4 * 500.0, 1.2)) == pytest.approx(2 * float(impl.bernoulli_airspeed(500.0, 1.2)))
    with pytest.raises(ValueError):
        impl.bernoulli_airspeed(-1.0, 1.2)


def test_compressible_pitot_tends_to_bernoulli_with_the_m2_over_4_term():
    T, p = 268.65, 70108.5
    a = np.sqrt(1.4 * impl.R_AIR * T)
    rho = p / (impl.R_AIR * T)
    for M in (0.01, 0.05, 0.1):
        V = M * a
        ratio = float(impl.impact_pressure(V, p, T)) / float(impl.dynamic_pressure(rho, V))
        assert ratio - 1.0 == pytest.approx(M ** 2 / 4, rel=0.02)
    assert float(impl.impact_pressure(0.0, p, T)) == 0.0
    with pytest.raises(ValueError):
        impl.impact_pressure(1.01 * a, p, T)


# ---------------------------------------------------------------- the airspeed chain

def test_at_sea_level_cas_equals_eas_equals_tas():
    r = impl.airspeeds_from_tas(75.0, 0.0)
    assert r["CAS"] == pytest.approx(75.0, rel=1e-12)
    assert r["EAS"] == pytest.approx(75.0, rel=1e-4)
    assert r["IAS"] == r["CAS"]


def test_cas_and_impact_pressure_are_inverse():
    cas = np.array([20.0, 51.75, 120.0, 250.0])
    npt.assert_allclose(impl.cas_from_impact_pressure(impl.impact_pressure_from_cas(cas)), cas, rtol=1e-13)


@pytest.mark.parametrize("h", [0.0, 1500.0, 3000.0, 9000.0])
def test_tas_cas_round_trip(h):
    for tas in (30.0, 60.0, 150.0):
        assert impl.tas_from_cas(impl.airspeeds_from_tas(tas, h)["CAS"], h) == pytest.approx(tas, rel=1e-12)


def test_eas_keeps_the_dynamic_pressure_and_cas_exceeds_it_slightly():
    for h in (1000.0, 3000.0, 6000.0):
        r = impl.airspeeds_from_tas(60.0, h)
        assert impl.RHO0 * r["EAS"] ** 2 == pytest.approx(r["rho"] * 60.0 ** 2, rel=1e-12)
        assert r["EAS"] < r["CAS"] < 60.0
        assert (r["CAS"] - r["EAS"]) / r["EAS"] < 3e-3                    # a light aeroplane: CAS ~ EAS
    gaps = [impl.airspeeds_from_tas(150.0, h)["CAS"] - impl.airspeeds_from_tas(150.0, h)["EAS"] for h in (0.0, 3000.0, 6000.0)]
    assert np.all(np.diff(gaps) > 0)                                       # grows with altitude


def test_static_port_error_signs_and_small_error_linearisation():
    base = impl.airspeeds_from_tas(60.0, 3000.0)
    low = impl.airspeeds_from_tas(60.0, 3000.0, cp_static=-0.02)            # port in fast air: low pressure
    assert low["IAS"] > base["IAS"] and low["p_static_sensed"] < base["p"]
    e = impl.static_error_effects(-0.02, base["q"], base["rho"])
    assert e["altitude_error"] > 0                                         # altimeter reads high
    # dCAS ~ -cp q / (rho0 CAS) for a small error at low Mach number
    assert low["IAS"] - base["IAS"] == pytest.approx(0.02 * base["q"] / (impl.RHO0 * base["CAS"]), rel=0.02)
    assert impl.airspeeds_from_tas(60.0, 3000.0, instrument_error=0.5)["IAS"] == pytest.approx(base["IAS"] + 0.5)


def test_altimeter_error_is_consistent_with_the_hydrostatic_gradient():
    r = impl.airspeeds_from_tas(60.0, 3000.0)
    e = impl.static_error_effects(-0.03, r["q"], r["rho"])
    _, p_up, _ = impl.isa_troposphere(3000.0 + e["altitude_error"])
    assert float(p_up) - r["p"] == pytest.approx(e["dp"], rel=2e-3)          # the sensed pressure, found higher up


# ---------------------------------------------------------------- Euler equation, contraction, venturi

def test_euler_integration_converges_to_bernoulli_at_second_order():
    A = lambda x: impl.contraction_area(x, 4.0, 0.5, 2.0)
    errs = []
    for n in (25, 50, 100):
        x = np.linspace(0.0, 2.0, n + 1)
        errs.append(np.max(np.abs(impl.euler_pressure_along_duct(A, x, 15.0, 1.204, 0.0)
                                  - impl.bernoulli_pressure_along_duct(A, x, 15.0, 1.204, 0.0))))
    orders = impl.observed_order(errs)
    assert np.all(orders > 1.8) and np.all(orders < 2.2)


def test_bernoulli_keeps_total_pressure_and_a_straight_duct_keeps_static():
    A = lambda x: impl.contraction_area(x, 4.0, 0.5, 2.0)
    x = np.linspace(0.0, 2.0, 101)
    p = impl.bernoulli_pressure_along_duct(A, x, 15.0, 1.204, 1000.0)
    V = 15.0 / A(x)
    npt.assert_allclose(p + 0.5 * 1.204 * V ** 2, 1000.0 + 0.5 * 1.204 * V[0] ** 2, rtol=1e-14)
    flat = lambda xx: np.full_like(xx, 2.0)
    npt.assert_allclose(impl.euler_pressure_along_duct(flat, x, 15.0, 1.204, 1000.0), 1000.0, atol=1e-10)


def test_contraction_area_and_venturi_inversion():
    x = np.array([0.0, 1.0, 2.0])
    npt.assert_allclose(impl.contraction_area(x, 4.0, 0.5, 2.0), [4.0, 2.25, 0.5])
    A1, A2, rho, Q = 4.0, 0.5, 1.204, 15.0
    dp = 0.5 * rho * ((Q / A2) ** 2 - (Q / A1) ** 2)
    assert impl.venturi_flow_rate(dp, A1, A2, rho) == pytest.approx(Q, rel=1e-14)
    assert impl.venturi_flow_rate(dp, A1, A2, rho, Cd=0.98) == pytest.approx(0.98 * Q)


# ---------------------------------------------------------------- the energy equation

def test_energy_equation_reduces_to_bernoulli_without_losses_or_work():
    rho, V1, V2, p1, z1, z2 = 998.2, 2.0, 6.0, 2e5, 0.0, 3.0
    p2 = p1 + 0.5 * rho * (V1 ** 2 - V2 ** 2) + rho * impl.G0 * (z1 - z2)   # Bernoulli
    assert impl.shaft_head(p1, V1, z1, p2, V2, z2, rho) == pytest.approx(0.0, abs=1e-12)
    assert impl.shaft_head(p1, V1, z1, p2, V2, z2, rho, h_loss=1.7) == pytest.approx(1.7)
    assert impl.shaft_head(p2, V2, z2, p1, V1, z1, rho, h_loss=0.0) == pytest.approx(0.0, abs=1e-12)
    # a turbine: the outlet pressure is lower than Bernoulli allows, so the head is negative
    assert impl.shaft_head(p1, V1, z1, p2 - 5e4, V2, z2, rho) < 0


def test_kinetic_energy_correction_and_power():
    h1 = impl.shaft_head(0.0, 3.0, 0.0, 0.0, 3.0, 0.0, 998.2, alpha1=1.0, alpha2=2.0)
    assert h1 == pytest.approx(9.0 / (2 * impl.G0))                       # laminar outlet alpha = 2
    P = impl.shaft_power(10.0, 998.2, 0.01)
    assert P == pytest.approx(998.2 * impl.G0 * 0.01 * 10.0)
    assert impl.shaft_power(10.0, 998.2, 0.01, efficiency=0.5) == pytest.approx(2 * P)
    with pytest.raises(ValueError):
        impl.shaft_power(10.0, 998.2, 0.01, efficiency=0.0)


def test_open_tunnel_power_scales_with_cube_of_speed_and_falls_with_diffuser_area():
    P1 = impl.open_tunnel_fan_power(20.0, 0.5, 1.204, 0.6, 1.0)
    P2 = impl.open_tunnel_fan_power(40.0, 0.5, 1.204, 0.6, 1.0)
    assert P2 / P1 == pytest.approx(8.0, rel=1e-12)
    assert impl.open_tunnel_fan_power(20.0, 0.5, 1.204, 0.6, 2.0) < P1
    # no losses and exit = test section: the fan supplies exactly the kinetic-energy flux rho Q V^2 / 2
    assert impl.open_tunnel_fan_power(20.0, 0.5, 1.204, 0.0, 0.5) == pytest.approx(0.5 * 1.204 * 0.5 * 20.0 ** 3)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "pitot_static.png").stat().st_size > 10_000
