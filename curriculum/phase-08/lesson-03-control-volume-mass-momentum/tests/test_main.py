"""Public tests for P08.03 Control Volume: Mass & Momentum. They assert conservation, not magic
numbers: a uniform pressure exerts no net force on a closed control surface; the box around a
turbojet balances mass to round-off and returns the textbook thrust whatever its size; the thrust
is unchanged in a frame moving along the axis (Galilean invariance) because mass balances, and is
not when a face is left out; limiting cases (no flow, a static engine, a perfectly expanded
nozzle); momentum-flux correction factors from quadrature against their closed forms and never
below 1; the rocket equation by RK4 against Tsiolkovsky, at fourth order; and a moving vane whose
power peaks at a third of the jet speed.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import dataclasses

import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

CRUISE = dict(mdot_a=10.0, f=0.02, V0=236.0, Ve=650.0, pe=40000.0, Te=900.0, pa=22632.0, rho_a=0.36392)


def box(A_box, **kw):
    c = dict(CRUISE, **kw)
    return impl.turbojet_control_volume(c["mdot_a"], c["f"], c["V0"], c["Ve"], c["pe"], c["Te"], c["pa"], c["rho_a"], A_box)


def formula(**kw):
    c = dict(CRUISE, **kw)
    Ae = impl.exit_area((1 + c["f"]) * c["mdot_a"], c["Ve"], c["pe"], c["Te"])
    return impl.turbojet_thrust(c["mdot_a"], c["f"], c["V0"], c["Ve"], c["pe"], c["pa"], Ae)


# ---------------------------------------------------------------- face sums

def test_uniform_pressure_exerts_no_net_force_on_a_closed_surface():
    p = 7.3e4
    faces = [impl.Face(2.0, -1.0, 1.2, 0.0, 0.0, p), impl.Face(2.0, 1.0, 1.2, 0.0, 0.0, p),
             impl.Face(5.0, 0.0, 1.2, 0.0, 0.0, p)]
    assert impl.pressure_force_x(faces) == 0.0
    assert impl.mass_outflow(faces) == 0.0 and impl.solid_force_x(faces) == 0.0


def test_face_sums_for_a_straight_duct():
    # uniform flow in one end and out of the other: no net mass, no net momentum, no force
    rho, V, A = 1.2, 30.0, 0.5
    faces = [impl.Face(A, -1.0, rho, -V, V, 1e5), impl.Face(A, 1.0, rho, V, V, 1e5)]
    assert impl.mass_outflow(faces) == pytest.approx(0.0, abs=1e-12)
    assert impl.momentum_outflow_x(faces) == pytest.approx(0.0, abs=1e-9)
    assert impl.solid_force_x(faces) == pytest.approx(0.0, abs=1e-9)


# ---------------------------------------------------------------- the turbojet box

@pytest.mark.parametrize("A_box", [0.3, 1.0, 10.0, 1000.0])
def test_box_mass_balance_closes_and_thrust_is_independent_of_box_size(A_box):
    faces = box(A_box)
    m_scale = CRUISE["rho_a"] * CRUISE["V0"] * A_box
    assert abs(impl.mass_outflow(faces)) < 1e-13 * m_scale
    F = formula()["F"]
    assert impl.solid_force_x(faces) == pytest.approx(F, rel=1e-12, abs=1e-14 * m_scale * CRUISE["V0"])


def test_exit_area_satisfies_continuity():
    mdot, Ve, pe, Te = 10.2, 650.0, 4e4, 900.0
    Ae = impl.exit_area(mdot, Ve, pe, Te)
    assert pe / (impl.R_AIR * Te) * Ve * Ae == pytest.approx(mdot, rel=1e-14)


def test_thrust_is_galilean_invariant_because_mass_balances():
    faces = box(5.0)
    W = 50.0                                                  # observer moving at -W along x
    shifted = [dataclasses.replace(f, ux=f.ux + W) for f in faces]
    assert impl.solid_force_x(shifted) == pytest.approx(impl.solid_force_x(faces), rel=1e-12)
    # leave out the side face and the balance no longer closes: the "force" then depends on the frame
    broken = faces[:3] + faces[4:]
    broken_shifted = shifted[:3] + shifted[4:]
    assert abs(impl.mass_outflow(broken)) > 1.0
    assert abs(impl.solid_force_x(broken_shifted) - impl.solid_force_x(broken)) > 10.0


def test_breakdown_adds_up_and_limits():
    r = formula()
    assert r["F"] == pytest.approx(r["gross"] - r["ram_drag"] + r["pressure"])
    assert r["pressure"] > 0                                                  # underexpanded nozzle
    s = formula(V0=0.0, pe=22632.0)
    assert s["ram_drag"] == 0.0 and s["pressure"] == 0.0
    assert s["F"] == pytest.approx(1.02 * 10.0 * 650.0)
    z = impl.turbojet_thrust(10.0, 0.0, 300.0, 300.0, 5e4, 5e4, 0.1)         # nothing changes: no thrust
    assert z["F"] == 0.0


def test_box_reproduces_pressure_thrust_slope():
    # with A_e held fixed, dF/dp_e = A_e
    Ae = 0.1
    F1 = impl.turbojet_thrust(10.0, 0.02, 236.0, 650.0, 4.0e4, 22632.0, Ae)["F"]
    F2 = impl.turbojet_thrust(10.0, 0.02, 236.0, 650.0, 4.1e4, 22632.0, Ae)["F"]
    assert (F2 - F1) / 1e3 == pytest.approx(Ae)


def test_box_rejects_a_box_smaller_than_the_jet():
    with pytest.raises(ValueError):
        box(0.01)


# ---------------------------------------------------------------- non-uniform profiles

def test_profile_fluxes_uniform_and_laminar():
    u = impl.profile_fluxes(lambda r: np.full_like(r, 3.0), 0.2, 1.1)
    assert u["beta"] == pytest.approx(1.0, rel=1e-12)
    assert u["mdot"] == pytest.approx(1.1 * 3.0 * np.pi * 0.04, rel=1e-12)
    lam = impl.profile_fluxes(lambda r: 2.0 * (1 - (r / 0.2) ** 2), 0.2, 1.1)
    assert lam["V_mean"] == pytest.approx(1.0, rel=1e-12)                     # mean = u_max / 2
    assert lam["beta"] == pytest.approx(impl.beta_laminar(), rel=1e-12)


@pytest.mark.parametrize("n", [6.0, 7.0, 9.0])
def test_power_law_beta_matches_closed_form(n):
    b = impl.profile_fluxes(lambda r: (1 - r) ** (1 / n), 1.0, 1.0, n=200)["beta"]
    assert b == pytest.approx(impl.beta_power_law(n), rel=1e-5)
    assert impl.beta_power_law(1e6) == pytest.approx(1.0, rel=1e-5)          # flat profile limit


def test_beta_is_never_below_one():
    rng = np.random.default_rng(8)                                           # seeded random positive profiles
    for _ in range(20):
        c = rng.uniform(0.1, 1.0, 4)
        b = impl.profile_fluxes(lambda r: c[0] + c[1] * r + c[2] * r ** 2 + c[3] * np.cos(5 * r), 1.0, 1.0)["beta"]
        assert b >= 1.0 - 1e-12                                              # Cauchy-Schwarz


# ---------------------------------------------------------------- rockets and moving control volumes

def test_rocket_thrust_vacuum_minus_sea_level_is_pa_Ae():
    Fsl = impl.rocket_thrust(5.0, 2500.0, 8e4, 101325.0, 0.05)
    Fvac = impl.rocket_thrust(5.0, 2500.0, 8e4, 0.0, 0.05)
    assert Fvac - Fsl == pytest.approx(101325.0 * 0.05)


def test_rocket_burn_matches_tsiolkovsky_at_fourth_order():
    m0, mp, tb, c = 100.0, 60.0, 8.0, 2000.0
    md = mp / tb
    exact = impl.tsiolkovsky_delta_v(c, m0, m0 - mp, tb)
    errs = [abs(impl.rocket_burn(m0, md, md * c, tb, n)[1][-1] - exact) for n in (10, 20, 40)]
    orders = impl.observed_order(errs)
    assert np.all(orders > 3.8) and np.all(orders < 4.2)
    t, v, h = impl.rocket_burn(m0, md, md * c, tb, 400)
    assert v[0] == 0.0 and h[0] == 0.0
    npt.assert_allclose(np.diff(h) > 0, True)                                # climbing all the way


def test_without_gravity_delta_v_is_c_log_mass_ratio():
    t, v, h = impl.rocket_burn(50.0, 2.0, 2.0 * 1800.0, 10.0, 200, g=0.0)
    assert v[-1] == pytest.approx(1800.0 * np.log(50.0 / 30.0), rel=1e-9)
    assert impl.tsiolkovsky_delta_v(1800.0, 50.0, 30.0) == pytest.approx(1800.0 * np.log(5 / 3))


def test_vane_forces_and_the_one_third_rule():
    rho, A, V = 998.2, 1e-4, 20.0
    Fx, Fy = impl.vane_force(rho, A, V, 0.0, 180.0)
    assert Fx == pytest.approx(2 * rho * A * V ** 2) and abs(Fy) < 1e-9
    assert impl.vane_force(rho, A, V, 0.0, 0.0) == (0.0, 0.0)
    Fx90, Fy90 = impl.vane_force(rho, A, V, 0.0, 90.0)
    assert Fx90 == pytest.approx(rho * A * V ** 2) and Fy90 == pytest.approx(-rho * A * V ** 2)
    us = np.linspace(0.0, V, 3001)
    P = [impl.vane_force(rho, A, V, u, 180.0)[0] * u for u in us]
    assert us[int(np.argmax(P))] == pytest.approx(V / 3, abs=V / 3000)


def test_thrust_reverser_fraction():
    assert impl.thrust_reverser_fraction(90.0) == pytest.approx(0.0, abs=1e-15)
    assert impl.thrust_reverser_fraction(180.0) == pytest.approx(1.0)
    assert impl.thrust_reverser_fraction(135.0) == pytest.approx(np.sqrt(0.5))


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "jet_thrust.png").stat().st_size > 10_000
