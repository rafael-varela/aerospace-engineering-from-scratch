"""Public tests for P04.02 Work, Energy & Power. They assert physics, not magic numbers: the work of
a conservative force is path-independent and that of drag is not; the energy height is conserved
exactly in flight without thrust and drag and changes at the rate P_s = (T - D)V/W with them; the
power and drag minima sit where the closed forms put them; the phugoid period matches Lanchester's
formula at constant density and the density-gradient correction in the ISA; the ceilings are
consistent with P_s; the energy-method climb is never slower than the alternatives.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-04/parameterized-cases/.
"""
import math
from dataclasses import replace

import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import solve_ivp

from conftest import load_impl

impl = load_impl(__file__)

G = impl.G0
AC = impl.Aircraft()


# ---------------------------------------------------------------- the atmosphere

def test_isa_anchor_values_continuity_and_range():
    T0, p0, r0 = impl.isa(0.0)
    assert (T0, p0) == (pytest.approx(288.15), pytest.approx(101325.0))
    assert r0 == pytest.approx(1.225, rel=1e-3)
    Ta, pa, ra = impl.isa(11_000.0 - 1e-6)
    Tb, pb, rb = impl.isa(11_000.0 + 1e-6)
    assert Ta == pytest.approx(216.65) and Tb == pytest.approx(216.65)
    assert pa == pytest.approx(pb, rel=1e-9)                                         # continuous at the tropopause
    assert impl.isa(11_000.0)[1] == pytest.approx(22632.0, rel=1e-3)
    assert impl.speed_of_sound(0.0) == pytest.approx(340.29, rel=1e-3)
    with pytest.raises(ValueError):
        impl.isa(25_000.0)


# ---------------------------------------------------------------- work and energy

def test_work_of_gravity_is_path_independent_and_drag_work_is_not():
    m = 3.0
    grav = lambda r, dr: np.array([0.0, -m * G])
    s = np.linspace(0.0, 1.0, 401)
    straight = np.column_stack([10 * s, 4 * s])
    wavy = np.column_stack([10 * s, 4 * s + 3 * np.sin(5 * np.pi * s)])
    w1, w2 = impl.work_along_path(grav, straight), impl.work_along_path(grav, wavy)
    assert w1 == pytest.approx(-m * G * 4.0, rel=1e-12)
    assert w2 == pytest.approx(w1, rel=1e-10)
    # closed loop: zero work for the conservative force, negative for drag
    loop = np.column_stack([np.cos(2 * np.pi * s), np.sin(2 * np.pi * s)])
    assert impl.work_along_path(grav, loop) == pytest.approx(0.0, abs=1e-12)
    D = 2.0
    drag = lambda r, dr: -D * dr / np.linalg.norm(dr)
    assert impl.work_along_path(drag, loop) == pytest.approx(-D * 2 * np.pi, rel=1e-4)    # -D x arc length
    assert impl.work_along_path(drag, wavy) < impl.work_along_path(drag, straight) < 0      # longer path, more loss


def test_work_energy_theorem_on_a_spring():
    # F = -c x along a straight line: W from x0 to x1 equals the change of -(c/2) x^2 (trapezoid is exact for linear F)
    c = 50.0
    path = np.linspace(0.3, -0.1, 7)[:, None]
    W = impl.work_along_path(lambda r, dr: -c * r, path)
    assert W == pytest.approx(-(c / 2) * (0.1 ** 2 - 0.3 ** 2), rel=1e-12)


def test_kinetic_energy_and_energy_height():
    assert impl.kinetic_energy(5000.0, 150.0) == pytest.approx(56.25e6)
    assert impl.energy_height(1000.0, 0.0) == pytest.approx(1000.0)
    assert impl.energy_height(0.0, math.sqrt(2 * G * 500.0)) == pytest.approx(500.0)
    npt.assert_allclose(impl.energy_height(np.array([0.0, 100.0]), np.array([10.0, 0.0])), [100.0 / (2 * G), 100.0])
    # kinetic energy per unit weight is the speed part of h_E
    assert impl.kinetic_energy(AC.mass, 120.0) / AC.W == pytest.approx(impl.energy_height(0.0, 120.0))
    assert impl.zoom_altitude(0.0, 150.0, 80.0) == pytest.approx((150.0 ** 2 - 80.0 ** 2) / (2 * G))
    assert impl.zoom_altitude(2000.0, 100.0, 100.0) == pytest.approx(2000.0)


# ---------------------------------------------------------------- forces and power

def test_drag_minimum_and_power_minimum_are_at_the_closed_form_speeds():
    h = 3000.0
    rho = impl.isa(h)[2]
    V = np.linspace(50.0, 200.0, 30001)
    D = impl.drag(AC, rho, V)
    Pr, Pa = impl.power_required_available(AC, h, V)
    assert V[np.argmin(D)] == pytest.approx(impl.min_drag_speed(AC, h), abs=0.02)
    assert V[np.argmin(Pr)] == pytest.approx(impl.min_power_speed(AC, h), abs=0.02)
    assert impl.min_power_speed(AC, h) == pytest.approx(impl.min_drag_speed(AC, h) / 3 ** 0.25)
    assert D.min() == pytest.approx(AC.W / impl.max_lift_to_drag(AC), rel=1e-6)          # D_min = W / (L/D)max
    npt.assert_allclose(Pa, impl.thrust(AC, rho) * V)
    # at the minimum-drag speed the induced and zero-lift drags are equal
    Vmd = impl.min_drag_speed(AC, h)
    q = 0.5 * rho * Vmd ** 2
    CL = AC.W / (q * AC.S)
    assert AC.CD0 == pytest.approx(AC.k * CL ** 2, rel=1e-12)


def test_drag_grows_with_load_factor_and_thrust_lapses_with_density():
    rho = impl.isa(0.0)[2]
    D1 = float(impl.drag(AC, rho, 100.0, n=1.0))
    D2 = float(impl.drag(AC, rho, 100.0, n=2.0))
    q = 0.5 * rho * 100.0 ** 2
    CL1 = AC.W / (q * AC.S)
    assert D2 - D1 == pytest.approx(q * AC.S * AC.k * 3 * CL1 ** 2, rel=1e-12)       # induced drag x n^2
    assert impl.thrust(AC, impl.RHO0) == pytest.approx(AC.T0)
    assert impl.thrust(AC, impl.RHO0 / 2) == pytest.approx(AC.T0 * 0.5 ** AC.thrust_exp)
    assert impl.max_lift_to_drag(AC) == pytest.approx(1 / (2 * math.sqrt(AC.k * AC.CD0)))


def test_specific_excess_power_definition_and_vectorisation():
    h, V = 2000.0, 140.0
    rho = impl.isa(h)[2]
    ps = float(impl.specific_excess_power(AC, h, V))
    assert ps == pytest.approx((impl.thrust(AC, rho) - float(impl.drag(AC, rho, V))) * V / AC.W, rel=1e-12)
    grid = impl.specific_excess_power(AC, np.array([[0.0], [5000.0]]), np.array([100.0, 150.0, 200.0]))
    assert grid.shape == (2, 3)
    assert np.all(grid[1] < grid[0])                                                # less thrust, less excess power at altitude
    # a turn (n = 2) costs excess power
    assert float(impl.specific_excess_power(AC, h, V, n=2.0)) < ps


def test_best_climb_and_ceilings_are_consistent():
    v0, p0 = impl.best_climb(AC, 0.0)
    v9, p9 = impl.best_climb(AC, 9000.0)
    assert p0 > p9 > 0
    assert v9 > v0                                                                    # best-climb TAS rises with altitude
    # the search found a maximum
    for dv in (-3.0, 3.0):
        assert float(impl.specific_excess_power(AC, 0.0, v0 + dv)) <= p0 + 1e-9
    h_abs, h_srv = impl.ceiling(AC), impl.ceiling(AC, 0.508)
    assert h_srv < h_abs
    assert impl.best_climb(AC, h_abs)[1] == pytest.approx(0.0, abs=1e-3)
    assert impl.best_climb(AC, h_srv)[1] == pytest.approx(0.508, abs=1e-3)
    # heavier: lower ceiling; more thrust: higher
    assert impl.ceiling(replace(AC, mass=1.1 * AC.mass), 0.508) < h_srv < impl.ceiling(replace(AC, T0=1.1 * AC.T0), 0.508)
    with pytest.raises(ValueError):
        impl.ceiling(replace(AC, T0=1000.0))


# ---------------------------------------------------------------- equations of motion

def test_energy_height_is_conserved_without_thrust_and_drag():
    h0, V0 = 3000.0, 150.0
    CL = AC.W / (0.5 * impl.isa(h0)[2] * V0 ** 2 * AC.S)
    fl = impl.simulate_flight(AC, h0, 1.08 * V0, 0.05, CL, 150.0, thrust_on=False, drag_on=False)
    hE = impl.energy_height(fl.h, fl.V)
    assert np.ptp(hE) < 1e-6                                                         # metres, over 150 s
    assert np.ptp(fl.h) > 50.0                                                       # while the altitude swings


def test_energy_rate_equals_specific_excess_power_along_a_trajectory():
    h0, V0 = 2000.0, 130.0
    CL = AC.W / (0.5 * impl.isa(h0)[2] * V0 ** 2 * AC.S)
    fl = impl.simulate_flight(AC, h0, V0, 0.0, CL, 60.0, n_samples=6001)
    hE = impl.energy_height(fl.h, fl.V)
    dhE = np.gradient(hE, fl.t)
    q = 0.5 * np.array([impl.isa(h)[2] for h in fl.h]) * fl.V ** 2
    T = np.array([impl.thrust(AC, impl.isa(h)[2]) for h in fl.h])
    D = q * AC.S * (AC.CD0 + AC.k * CL ** 2)
    npt.assert_allclose(dhE[2:-2], ((T - D) * fl.V / AC.W)[2:-2], atol=1e-4)
    # the right-hand side against an independent solve_ivp call
    ref = solve_ivp(impl.point_mass_rhs, (0, 60.0), [0.0, h0, V0, 0.0], args=(AC, CL), method="RK45", rtol=1e-11,
                    atol=1e-9, t_eval=fl.t[::600])
    npt.assert_allclose(ref.y[1], fl.h[::600], atol=1e-4)


def test_phugoid_matches_lanchester_at_constant_density_and_the_gradient_formula_in_the_isa():
    V0, h0 = 150.0, 3000.0
    rho = impl.isa(h0)[2]
    CL = AC.W / (0.5 * rho * V0 ** 2 * AC.S)
    flc = impl.simulate_flight(AC, h0, 1.01 * V0, 0.0, CL, 300.0, thrust_on=False, drag_on=False, const_rho=rho)
    assert impl.oscillation_period(flc.t, flc.h) == pytest.approx(impl.phugoid_period_lanchester(V0), rel=2e-3)
    assert impl.phugoid_period_lanchester(V0) == pytest.approx(math.sqrt(2) * math.pi * V0 / G)
    fl = impl.simulate_flight(AC, h0, 1.01 * V0, 0.0, CL, 300.0, thrust_on=False, drag_on=False)
    beta = impl.density_scale_gradient(h0)
    assert 0.9e-4 < beta < 1.2e-4                                                    # 1/m: a density scale height near 10 km
    T_isa = impl.oscillation_period(fl.t, fl.h)
    assert T_isa == pytest.approx(impl.phugoid_period_with_gradient(V0, beta), rel=5e-3)
    assert T_isa < impl.phugoid_period_lanchester(V0)                                 # the density gradient stiffens the phugoid
    assert impl.phugoid_period_with_gradient(V0, 0.0) == pytest.approx(impl.phugoid_period_lanchester(V0))


def test_oscillation_period_estimator():
    t = np.linspace(0, 10, 2001)
    assert impl.oscillation_period(t, 3 + np.sin(2 * np.pi * t / 2.5)) == pytest.approx(2.5, rel=1e-4)
    with pytest.raises(ValueError):
        impl.oscillation_period(t, t)


def test_glide_at_best_lift_to_drag():
    # a steady glide at constant density: L = W cos(gamma), D = W sin(-gamma), so x / height lost = L/D exactly
    h0 = 3000.0
    rho = impl.isa(h0)[2]
    CL = math.sqrt(AC.CD0 / AC.k)
    gam = -math.atan(1 / impl.max_lift_to_drag(AC))
    V = math.sqrt(2 * AC.W * math.cos(gam) / (rho * AC.S * CL))
    fl = impl.simulate_flight(AC, h0, V, gam, CL, 60.0, thrust_on=False, const_rho=rho)
    assert fl.x[-1] / (h0 - fl.h[-1]) == pytest.approx(impl.max_lift_to_drag(AC), rel=1e-8)
    npt.assert_allclose(fl.V, V, rtol=1e-8)
    assert np.all(np.diff(impl.energy_height(fl.h, fl.V)) < 0)                     # a glider only loses energy
    # any other C_L glides worse
    for CL2 in (0.7 * CL, 1.4 * CL):
        assert CL2 / (AC.CD0 + AC.k * CL2 ** 2) < impl.max_lift_to_drag(AC)


# ---------------------------------------------------------------- the energy method

def test_rutowski_path_lies_on_the_energy_curves_and_maximises_ps():
    grid = np.linspace(2000.0, 9000.0, 8)
    V, h, ps = impl.rutowski_path(AC, grid)
    npt.assert_allclose(impl.energy_height(h, V), grid, rtol=1e-10)
    assert np.all(h >= 0) and np.all(ps > 0)
    assert np.all(np.diff(ps) < 0)                                                   # P_s falls as the energy rises
    for e, v, p in zip(grid, V, ps):
        for dv in (-5.0, 5.0):
            hh = e - (v + dv) ** 2 / (2 * G)
            if hh >= 0:
                assert float(impl.specific_excess_power(AC, hh, v + dv)) <= p + 1e-6


def test_time_to_energy_closed_form_and_energy_method_is_fastest():
    grid = np.linspace(0.0, 1000.0, 11)
    assert impl.time_to_energy(grid, np.full(11, 10.0)) == pytest.approx(100.0)
    with pytest.raises(ValueError):
        impl.time_to_energy(grid, np.linspace(5.0, -1.0, 11))
    hE0, hE1 = float(impl.energy_height(0.0, 150.0)), 8000.0
    g2 = np.linspace(hE0, hE1, 81)
    t_opt = impl.time_to_energy(g2, impl.rutowski_path(AC, g2)[2])
    for Vc in (120.0, 150.0, 170.0):
        assert impl.climb_time_accelerate_then_constant(AC, 150.0, Vc, hE1, n=81) >= 0.99 * t_opt
    t_slow = impl.time_to_energy(g2, impl.constant_speed_climb_ps(AC, 150.0, g2))
    assert t_slow > t_opt


def test_acceleration_factor():
    assert impl.climb_rate_with_acceleration(20.0, 150.0, 0.0) == 20.0
    assert impl.climb_rate_with_acceleration(20.0, 150.0, 5e-3) < 20.0                 # accelerating while climbing costs climb rate
    assert impl.climb_rate_with_acceleration(20.0, 150.0, -5e-3) > 20.0                # decelerating feeds the climb


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "energy_height.png").exists()
