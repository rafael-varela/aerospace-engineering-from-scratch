"""Public tests for P04.05 Gravitation & Kepler's Laws. They assert the inverse-square law and its
consequences (g at altitude, circular and escape speeds, vis-viva), conservation of energy and angular
momentum by the leapfrog, its second-order convergence against a solve_ivp oracle, Kepler's three laws
on the simulated orbit, centre-of-mass motion in the full two-body problem, and the Hohmann transfer.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. The only file written goes to tmp_path.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import solve_ivp

from conftest import load_impl

impl = load_impl(__file__)

MU = impl.MU_EARTH
RE = impl.R_EARTH
A, E = 7500e3, 0.1


# ---------------------------------------------------------------- the law

def test_acceleration_is_inverse_square_and_points_at_the_centre():
    r = np.array([7000e3, 0.0])
    a1 = impl.gravity_acceleration(r)
    a2 = impl.gravity_acceleration(2 * r)
    assert a1[1] == 0.0 and a1[0] < 0.0                               # toward the origin
    assert np.linalg.norm(a1) == pytest.approx(MU / 7000e3 ** 2, rel=1e-12)
    assert np.linalg.norm(a1) / np.linalg.norm(a2) == pytest.approx(4.0, rel=1e-12)
    stack = impl.gravity_acceleration(np.array([r, 2 * r]))            # stacks of states
    assert stack.shape == (2, 2)
    npt.assert_allclose(stack[1], a2)


def test_newtons_third_law_for_the_pair_force():
    r1, r2 = np.array([0.0, 0.0, 0.0]), np.array([3.0, 4.0, 0.0])
    f21 = impl.gravitational_force(1e5, 2e5, r1, r2)
    f12 = impl.gravitational_force(2e5, 1e5, r2, r1)
    npt.assert_allclose(f21, -f12, rtol=1e-14)
    assert np.linalg.norm(f21) == pytest.approx(impl.G_NEWTON * 1e5 * 2e5 / 25.0, rel=1e-12)
    assert np.dot(f21, r2 - r1) < 0.0                                 # attractive


def test_g_at_altitude_limits():
    assert impl.g_at_altitude(0.0) == pytest.approx(MU / RE ** 2)
    assert impl.g_at_altitude(0.0) == pytest.approx(9.82, abs=0.01)
    assert impl.g_at_altitude(400e3) < impl.g_at_altitude(0.0)
    assert impl.g_at_altitude(RE) == pytest.approx(impl.g_at_altitude(0.0) / 4)
    assert impl.g_at_altitude(np.array([0.0, 400e3])).shape == (2,)


def test_circular_period_and_escape_speed():
    r = RE + 400e3
    v = impl.circular_speed(r)
    assert v == pytest.approx(np.sqrt(MU / r))
    assert 7600 < v < 7700                                            # ISS-like orbit
    assert impl.orbital_period(r) == pytest.approx(2 * np.pi * r / v)
    assert impl.escape_speed(r) == pytest.approx(np.sqrt(2.0) * v)
    # a circular orbit is a centripetal balance: v^2 / r = g(h)
    assert v ** 2 / r == pytest.approx(impl.g_at_altitude(400e3))


def test_vis_viva_reduces_to_circular_and_perigee_speeds():
    assert impl.vis_viva(1e7, 1e7) == pytest.approx(impl.circular_speed(1e7))
    r0, v0 = impl.perigee_state(A, E)
    assert np.linalg.norm(v0) == pytest.approx(impl.vis_viva(A * (1 - E), A))
    assert np.linalg.norm(v0) == pytest.approx(np.sqrt(MU * (1 + E) / (A * (1 - E))))
    assert impl.vis_viva(A * (1 + E), A) < impl.vis_viva(A * (1 - E), A)   # slower at apogee
    assert impl.vis_viva(2 * A, np.inf) == pytest.approx(impl.escape_speed(2 * A))   # parabolic: a = inf


def test_energy_and_elements_round_trip():
    r0, v0 = impl.perigee_state(A, E)
    assert impl.specific_energy(r0, v0) == pytest.approx(-MU / (2 * A), rel=1e-12)
    assert impl.specific_angular_momentum(r0, v0) == pytest.approx(np.sqrt(MU * A * (1 - E ** 2)), rel=1e-12)
    a, e, rp, ra = impl.orbital_elements(r0, v0)
    assert a == pytest.approx(A, rel=1e-10)
    assert e == pytest.approx(E, abs=1e-10)
    assert rp == pytest.approx(A * (1 - E), rel=1e-10)
    assert ra == pytest.approx(A * (1 + E), rel=1e-10)
    # the same orbit from a point away from perigee gives the same elements
    t, R, V = impl.propagate_leapfrog(A, E, 0.3, 4000)
    a2, e2, _, _ = impl.orbital_elements(R[-1], V[-1])
    assert a2 == pytest.approx(A, rel=1e-5)
    assert e2 == pytest.approx(E, abs=1e-4)
    # an escaping state has non-negative energy and infinite apoapsis
    _, e3, _, ra3 = impl.orbital_elements(r0, 1.5 * v0)
    assert e3 > 1.0 and ra3 == np.inf


# ---------------------------------------------------------------- the leapfrog

def _oscillator_accel(x):
    return -4.0 * x                                                   # omega = 2


def test_leapfrog_is_second_order_on_the_harmonic_oscillator():
    errs, hs = [], []
    for n in (100, 200, 400, 800):
        h = 3.0 / n                                                   # t_end = 3 is not a whole period
        t, X, V = impl.leapfrog(_oscillator_accel, np.array([1.0]), np.array([0.0]), h, n)
        assert t.shape == (n + 1,) and X.shape == (n + 1, 1)
        assert t[-1] == pytest.approx(3.0, rel=1e-12)
        errs.append(abs(X[-1, 0] - np.cos(2 * t[-1])))
        hs.append(h)
    assert 1.9 < impl.convergence_order(np.array(hs), np.array(errs)) < 2.1


def test_leapfrog_is_time_reversible():
    r0, v0 = impl.perigee_state(A, E)
    h = impl.orbital_period(A) / 500
    accel = lambda r: impl.gravity_acceleration(r)
    _, R, V = impl.leapfrog(accel, r0, v0, h, 137)
    _, Rb, Vb = impl.leapfrog(accel, R[-1], -V[-1], h, 137)          # reverse the velocity and go back
    npt.assert_allclose(Rb[-1], r0, rtol=0, atol=1e-5)                 # metres on a 7000 km orbit: rounding
    npt.assert_allclose(-Vb[-1], v0, rtol=0, atol=1e-8)


def test_leapfrog_conserves_angular_momentum_to_rounding_and_bounds_the_energy_error():
    spo = 500
    t, R, V = impl.propagate_leapfrog(A, E, 10.0, spo)
    h = impl.specific_angular_momentum(R, V)
    assert np.max(np.abs(h - h[0]) / abs(h[0])) < 1e-12                # exact for a central force
    dE = impl.energy_drift(R, V)
    assert np.max(np.abs(dE)) < 1e-4
    # bounded: the error over ten orbits is no larger than over the first one (no secular drift)
    assert np.max(np.abs(dE)) <= 1.05 * np.max(np.abs(dE[:spo + 1]))
    # and it shrinks four-fold when the step halves
    _, R2, V2 = impl.propagate_leapfrog(A, E, 10.0, 2 * spo)
    assert np.max(np.abs(impl.energy_drift(R2, V2))) == pytest.approx(np.max(np.abs(dE)) / 4, rel=0.1)


def test_leapfrog_converges_to_the_solve_ivp_oracle_at_order_two():
    t_ref, R_ref, _ = impl.propagate_solve_ivp(A, E, 1.0, n_out=501)
    hs, errs = [], []
    for spo in (250, 500, 1000):
        t, R, V = impl.propagate_leapfrog(A, E, 1.0, spo)
        Ri = np.column_stack([np.interp(t_ref, t, R[:, 0]), np.interp(t_ref, t, R[:, 1])])
        errs.append(impl.position_error(Ri, R_ref))
        hs.append(impl.orbital_period(A) / spo)
    assert errs[-1] < 2e3                                             # under 2 km after one orbit
    assert 1.8 < impl.convergence_order(np.array(hs), np.array(errs)) < 2.2


def test_two_body_rhs_matches_solve_ivp_closed_orbit():
    r0, v0 = impl.perigee_state(A, E)
    T = impl.orbital_period(A)
    sol = solve_ivp(impl.two_body_rhs, (0.0, T), np.concatenate([r0, v0]), method="DOP853", rtol=1e-12, atol=1e-6)
    npt.assert_allclose(sol.y[:2, -1], r0, atol=1.0)                   # back at perigee within 1 m
    npt.assert_allclose(sol.y[2:, -1], v0, atol=1e-3)


# ---------------------------------------------------------------- Kepler's laws

def test_keplers_first_law_the_orbit_is_an_ellipse_with_the_primary_at_a_focus():
    t, R, V = impl.propagate_leapfrog(A, E, 1.0, 2000)
    rn = np.linalg.norm(R, axis=1)
    assert rn.min() == pytest.approx(A * (1 - E), rel=1e-6)
    assert rn.max() == pytest.approx(A * (1 + E), rel=1e-5)        # apogee shifts with the bounded energy error
    # sum of distances to the two foci is 2a everywhere on an ellipse (second focus at (-2ae, 0))
    d2 = np.linalg.norm(R - np.array([-2 * A * E, 0.0]), axis=1)
    npt.assert_allclose(rn + d2, 2 * A, rtol=1e-4)


def test_keplers_second_law_equal_areas_in_equal_times():
    spo = 1000
    t, R, V = impl.propagate_leapfrog(A, 0.5, 2.0, spo)                 # a very eccentric orbit
    dA = impl.areal_velocity(R, t[1] - t[0])
    assert dA.shape == (R.shape[0] - 1,)
    assert (dA.max() - dA.min()) / dA.mean() < 1e-10
    assert dA.mean() == pytest.approx(0.5 * abs(impl.specific_angular_momentum(R[0], V[0])), rel=1e-6)


def test_keplers_third_law_from_simulated_periods():
    ratios = []
    for a in (7000e3, 12000e3, 30000e3):
        t, R, V = impl.propagate_leapfrog(a, 0.3, 3.0, 800)
        T = impl.period_from_trajectory(t, R)
        assert T == pytest.approx(impl.orbital_period(a), rel=1e-4)
        ratios.append(impl.keplers_third_law_ratio(a, T))
    ratios = np.array(ratios)
    npt.assert_allclose(ratios, 4 * np.pi ** 2 / MU, rtol=2e-4)
    with pytest.raises(ValueError):
        impl.period_from_trajectory(*impl.propagate_leapfrog(A, E, 0.5, 200)[:2])


# ---------------------------------------------------------------- the full two-body problem

def test_full_two_body_conserves_momentum_and_reduces_to_the_one_body_problem():
    m1, m2 = impl.M_EARTH, impl.M_MOON
    mu = impl.G_NEWTON * (m1 + m2)
    d = impl.D_EARTH_MOON
    v = impl.circular_speed(d, mu)
    T = impl.orbital_period(d, mu)
    sim = impl.simulate_two_body(m1, m2, np.array([d, 0.0]), np.array([0.0, v]), T, n_out=201)
    p_scale = m2 * v
    assert np.max(np.linalg.norm(sim["momentum"], axis=1)) / p_scale < 1e-9
    assert np.max(np.linalg.norm(sim["com"] - sim["com"][0], axis=1)) < 1.0
    # the relative motion is the one-body problem with mu = G (m1 + m2): a circle of radius d
    rn = np.linalg.norm(sim["r_rel"], axis=1)
    npt.assert_allclose(rn, d, rtol=1e-8)
    npt.assert_allclose(sim["r_rel"][-1], sim["r_rel"][0], atol=d * 1e-7)   # one period later: back
    # the heavy body barely moves: its orbit radius is m2/(m1+m2) of the separation
    assert np.max(np.linalg.norm(sim["r1"], axis=1)) == pytest.approx(m2 / (m1 + m2) * d, rel=1e-6)
    # the rhs is antisymmetric in the forces (Newton's third law)
    y = np.array([0.0, 0.0, d, 0.0, 0.0, 0.0, 0.0, v])
    dy = impl.two_body_full_rhs(0.0, y, m1, m2)
    npt.assert_allclose(m1 * dy[4:6] + m2 * dy[6:8], 0.0, atol=1e-12 * m2 * np.linalg.norm(dy[6:8]))


# ---------------------------------------------------------------- the transfer

def test_hohmann_transfer_known_case_and_limits():
    r1 = RE + 400e3
    r_geo = (MU * (impl.DAY_SIDEREAL / (2 * np.pi)) ** 2) ** (1 / 3)
    dv1, dv2, tof = impl.hohmann_transfer(r1, r_geo)
    assert dv1 + dv2 == pytest.approx(3857.0, rel=0.01)                 # about 3.9 km/s, LEO to GEO
    assert dv1 > dv2 > 0.0
    assert tof == pytest.approx(np.pi * np.sqrt((0.5 * (r1 + r_geo)) ** 3 / MU))
    assert 5.0 * 3600 < tof < 5.5 * 3600
    # the transfer orbit's energy sits between the two circles
    eps_t = -MU / (r1 + r_geo)
    assert -MU / (2 * r1) < eps_t < -MU / (2 * r_geo)
    # no transfer, no burn; and the reverse transfer costs the same
    assert impl.hohmann_transfer(r1, r1) == pytest.approx((0.0, 0.0, 0.5 * impl.orbital_period(r1)))
    back = impl.hohmann_transfer(r_geo, r1)
    assert back[0] + back[1] == pytest.approx(dv1 + dv2)
    # the first burn alone raises the apogee to r2 exactly
    v_after = impl.circular_speed(r1) + dv1
    _, _, _, ra = impl.orbital_elements(np.array([r1, 0.0]), np.array([0.0, v_after]))
    assert ra == pytest.approx(r_geo, rel=1e-9)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "two_body.png").exists()
