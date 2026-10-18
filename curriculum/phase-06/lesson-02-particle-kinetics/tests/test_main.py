"""Public tests for P06.02 Particle Kinetics. They assert mechanics, not magic numbers: the pull-up load
factor n = 1 + V^2/(gR) and its top-of-loop and straight-line limits; the level-turn bank, radius and rate
satisfy both force balances; seat and wing loads scale with n; the closed-form dive-recovery altitude
loss and time solve their ODEs and agree with `solve_ivp`; the load factor recovered from integrated
positions alone equals the commanded one; a loop flown with n = cos(gamma) + V^2/(gR) is a circle;
with thrust equal to drag the specific energy is conserved; a slower g-onset always costs altitude.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import quad

from conftest import load_impl

impl = load_impl(__file__)

G = 9.80665


# ---------------------------------------------------------------- free-body results

def test_bottom_of_pull_up_load_factor_and_its_inverse():
    V, R = 250.0, 1000.0
    assert float(impl.loop_load_factor(V, R)) == pytest.approx(1 + V ** 2 / (G * R))
    for n in (2.0, 5.0, 9.0):
        R_n = float(impl.radius_for_load_factor(V, n))
        assert float(impl.loop_load_factor(V, R_n)) == pytest.approx(n, rel=1e-12)
    # dimensionless: doubling V and quadrupling R leaves n unchanged
    assert float(impl.loop_load_factor(2 * V, 4 * R)) == pytest.approx(float(impl.loop_load_factor(V, R)))


def test_loop_limits_top_bottom_and_straight_line():
    V, R = 200.0, 1500.0
    assert float(impl.loop_load_factor(V, R, np.pi)) == pytest.approx(V ** 2 / (G * R) - 1)
    # bottom minus top is always 2 g at constant speed
    assert float(impl.loop_load_factor(V, R) - impl.loop_load_factor(V, R, np.pi)) == pytest.approx(2.0)
    # an infinite radius is straight flight: n = cos(gamma)
    assert float(impl.loop_load_factor(V, np.inf, 0.3)) == pytest.approx(np.cos(0.3))
    assert np.isinf(float(impl.radius_for_load_factor(V, 1.0)))


def test_level_turn_satisfies_both_force_balances():
    V, n = 180.0, 4.0
    phi, R, om = impl.level_turn(V, n)
    W = 1.0
    L = n * W
    assert L * np.cos(phi) == pytest.approx(W)                           # vertical balance
    assert L * np.sin(phi) == pytest.approx((W / G) * V ** 2 / R)       # centripetal
    assert om == pytest.approx(V / R)
    assert impl.level_turn(V, 1.0)[1] == np.inf
    with pytest.raises(ValueError):
        impl.level_turn(V, 0.8)


def test_seat_and_wing_loads_scale_with_load_factor():
    assert float(impl.apparent_weight(100.0, 1.0)) == pytest.approx(100.0 * G)
    npt.assert_allclose(impl.apparent_weight(80.0, np.array([1.0, 4.5, 9.0])), 80.0 * G * np.array([1.0, 4.5, 9.0]))
    assert float(impl.apparent_weight(100.0, 0.0)) == 0.0                # weightless at n = 0


# ---------------------------------------------------------------- closed-form recovery

def test_closed_form_recovery_matches_quadrature():
    V, n, g0 = 220.0, 6.0, np.radians(-45.0)
    dh_q, _ = quad(lambda g: V ** 2 * np.sin(g) / (G * (n - np.cos(g))), g0, 0.0)
    t_q, _ = quad(lambda g: V / (G * (n - np.cos(g))), g0, 0.0)
    assert impl.pullout_altitude_loss(V, n, g0) == pytest.approx(-dh_q, rel=1e-10)
    assert impl.pullout_time(V, n, g0) == pytest.approx(t_q, rel=1e-10)


def test_closed_form_recovery_limits():
    V = 250.0
    assert impl.pullout_altitude_loss(V, 5.0, 0.0) == pytest.approx(0.0, abs=1e-12)    # already level
    # small dive angle: delta_h ~ (V^2/g) gamma0^2 / (2 (n - 1)), a circle of radius V^2/(g(n-1))
    g0 = -0.01
    assert impl.pullout_altitude_loss(V, 5.0, g0) == pytest.approx(V ** 2 / G * g0 ** 2 / (2 * 4.0), rel=1e-3)
    # more g, less altitude; scales as V^2
    assert impl.pullout_altitude_loss(V, 8.0, -1.0) < impl.pullout_altitude_loss(V, 4.0, -1.0)
    assert impl.pullout_altitude_loss(2 * V, 6.0, -1.0) == pytest.approx(4 * impl.pullout_altitude_loss(V, 6.0, -1.0))
    with pytest.raises(ValueError):
        impl.pullout_altitude_loss(V, 1.0, -0.5)


# ---------------------------------------------------------------- the point-mass model

def test_rhs_straight_and_level_is_an_equilibrium():
    n1 = impl.ramp_load_factor(1.0, 1.0, np.inf)
    d = impl.flight_path_rhs(0.0, np.array([0.0, 1000.0, 200.0, 0.0]), n1)
    npt.assert_allclose(d, [200.0, 0.0, 0.0, 0.0], atol=1e-14)
    # unloaded (n = 0) the path curves down at g/V; without speed hold a 30 deg climb loses g/2
    d = impl.flight_path_rhs(0.0, np.array([0.0, 0.0, 100.0, np.radians(30)]), lambda t, y: 0.0, speed_hold=False)
    assert d[3] == pytest.approx(-G * np.cos(np.radians(30)) / 100.0)
    assert d[2] == pytest.approx(-G / 2)


def test_ramp_load_factor():
    f = impl.ramp_load_factor(1.0, 7.0, 3.0)
    assert f(0.0, None) == 1.0 and f(1.0, None) == pytest.approx(4.0) and f(5.0, None) == 7.0
    assert impl.ramp_load_factor(1.0, 7.0, np.inf)(0.0, None) == 7.0


def test_simulated_recovery_matches_closed_form():
    V, n, g0 = 250.0, 7.0, np.radians(-60.0)
    sol = impl.simulate_pullout(V, g0, impl.ramp_load_factor(n, n, np.inf), h0=0.0)
    assert sol.status == 1                                               # stopped by the level-off event
    assert -sol.y_events[0][0][1] == pytest.approx(impl.pullout_altitude_loss(V, n, g0), rel=1e-8)
    assert sol.t_events[0][0] == pytest.approx(impl.pullout_time(V, n, g0), rel=1e-8)
    assert sol.y_events[0][0][3] == pytest.approx(0.0, abs=1e-9)


def test_load_factor_recovered_from_positions_alone():
    V, g0 = 230.0, np.radians(-40.0)
    cmd = impl.ramp_load_factor(1.0, 6.0, 5.0)
    sol = impl.simulate_pullout(V, g0, cmd, h0=2000.0)
    t = np.linspace(0.0, sol.t_events[0][0], 3001)
    y = sol.sol(t)
    n_kin = impl.kinematic_load_factor(t, y[0], y[1])
    n_cmd = np.array([cmd(tt, None) for tt in t])
    inner = slice(20, -20)                       # np.gradient is one-sided at the ends; the ramp has a kink at 1 s
    kink = np.abs(t - 1.0) < 0.02
    ok = np.zeros_like(t, dtype=bool)
    ok[inner] = True
    ok &= ~kink
    npt.assert_allclose(n_kin[ok], n_cmd[ok], atol=2e-3)


def test_constant_speed_loop_is_a_circle():
    V, R = 200.0, 1200.0
    law = lambda t, y: float(impl.loop_load_factor(V, R, y[3]))
    sol = impl.simulate_pullout(V, 0.0, law, h0=500.0, t_max=2 * np.pi * R / V, stop_level=False)
    rad = np.hypot(sol.y[0], sol.y[1] - (500.0 + R))
    npt.assert_allclose(rad, R, atol=1e-5)
    assert sol.y[3, -1] == pytest.approx(2 * np.pi, rel=1e-8)              # one full turn of the velocity


def test_energy_conserved_when_thrust_equals_drag():
    sol = impl.simulate_pullout(250.0, np.radians(-50.0), impl.ramp_load_factor(1.0, 6.0, 4.0),
                                speed_hold=False, h0=0.0)
    E = 0.5 * sol.y[2] ** 2 + G * sol.y[1]
    assert np.max(np.abs(E - E[0])) < 1e-8 * E[0]
    assert sol.y[2, -1] > 250.0                                          # the aircraft speeds up in the dive
    # excess thrust adds energy at the rate g V (T - D)/W
    s2 = impl.simulate_pullout(200.0, 0.0, impl.ramp_load_factor(1.0, 1.0, np.inf), speed_hold=False,
                               excess_thrust=0.1, h0=0.0, t_max=5.0, stop_level=False)
    E2 = 0.5 * s2.y[2] ** 2 + G * s2.y[1]
    assert E2[-1] - E2[0] == pytest.approx(G * 0.1 * (s2.y[0, -1] - s2.y[0, 0]), rel=1e-8)


def test_slow_onset_costs_altitude_monotonically():
    V, g0, n = 250.0, np.radians(-60.0), 7.0
    losses = [impl.recovery_altitude_loss(V, g0, n, r) for r in (1.0, 3.0, 10.0, np.inf)]
    assert np.all(np.diff(losses) < 0)
    assert losses[-1] == pytest.approx(impl.pullout_altitude_loss(V, n, g0), rel=1e-7)
    # without speed hold the aircraft accelerates in the dive and needs more height
    assert impl.recovery_altitude_loss(V, g0, n, 6.0, speed_hold=False) > impl.recovery_altitude_loss(V, g0, n, 6.0)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "pullout.png").stat().st_size > 10_000
