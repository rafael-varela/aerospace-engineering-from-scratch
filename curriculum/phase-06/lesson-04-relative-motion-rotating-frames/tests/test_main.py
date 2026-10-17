"""Public tests for P06.04 Relative Motion & Rotating Frames. They assert mechanics, not magic numbers:
the five-term acceleration equals the second derivative of a point's inertial position computed by
brute force; its limits (no rotation, a point fixed in a turning frame); Coriolis does no work; frame
conversions round-trip and the local ENU basis is orthonormal; the independent check: the same shot
integrated in the rotating Earth-fixed frame and in the inertial frame lands in the same place, and the
rotating run conserves the Jacobi integral; the flat-Earth Coriolis integration agrees with the
first-order closed form; deflection is to the right in the north, to the left in the south, and a
northward shot at the equator lands west by Omega T^2 v_up / 3; the Foucault pendulum's swing plane
turns at exactly -Omega sin(latitude), matching the closed form from the non-rotating view.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

OM = 7.2921159e-5
G = 9.80665


def _Rz(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])


# ---------------------------------------------------------------- the five terms

def test_five_term_acceleration_equals_brute_force_second_derivative():
    # a frame whose origin accelerates and whose angle about z is th(t) = 0.3 + 0.8 t + 0.25 t^2
    RO = lambda t: np.array([1.0 + 2.0 * t + 0.5 * 3.0 * t ** 2, -t, 0.2 * t ** 2])
    th = lambda t: 0.3 + 0.8 * t + 0.25 * t ** 2
    rb = lambda t: np.array([2.0 + np.sin(t), 0.5 * t ** 2, 1.0 - t])       # position in the moving frame
    rI = lambda t: RO(t) + _Rz(th(t)) @ rb(t)
    t, h = 0.7, 1e-4
    a_brute = (rI(t + h) - 2 * rI(t) + rI(t - h)) / h ** 2
    Q = _Rz(th(t))
    w = np.array([0, 0, 0.8 + 0.5 * t])
    al = np.array([0, 0, 0.5])
    vrel = Q @ ((rb(t + h) - rb(t - h)) / (2 * h))
    arel = Q @ ((rb(t + h) - 2 * rb(t) + rb(t - h)) / h ** 2)
    aO = np.array([3.0, 0.0, 0.4])
    a, terms = impl.five_term_acceleration(aO, al, w, Q @ rb(t), vrel, arel)
    npt.assert_allclose(a, a_brute, rtol=1e-5, atol=1e-5)
    assert set(terms) == {"origin", "euler", "centripetal", "coriolis", "relative"}
    npt.assert_allclose(sum(terms.values()), a)


def test_five_term_limits():
    z = np.zeros(3)
    a, _ = impl.five_term_acceleration([1.0, 2, 3], z, z, [5.0, 0, 0], [0, 1.0, 0], [0, 0, 4.0])
    npt.assert_allclose(a, [1.0, 2.0, 7.0])                                  # no rotation: a_O + a_rel
    w = np.array([0, 0, 2.0])
    a, t = impl.five_term_acceleration(z, z, w, [3.0, 0, 0], z, z)          # point fixed in a turning frame
    npt.assert_allclose(a, [-12.0, 0, 0])                                    # omega^2 r toward the axis
    npt.assert_allclose(t["coriolis"], 0.0)


def test_coriolis_does_no_work_and_is_perpendicular_to_spin():
    rng = np.random.default_rng(2)
    for _ in range(5):
        r, v = rng.normal(size=3), rng.normal(size=3)
        f = np.zeros(3)
        a = impl.rotating_frame_acceleration(r, v, f, omega=np.array([0, 0, 1.3]))
        cor = a + np.cross([0, 0, 1.3], np.cross([0, 0, 1.3], r))            # remove the centrifugal part
        assert abs(cor @ v) < 1e-12 and abs(cor[2]) < 1e-12


def test_frame_conversions_round_trip():
    r, v, t = np.array([6.4e6, 1e5, 2e6]), np.array([100.0, -300.0, 50.0]), 1234.0
    rR, vR = impl.inertial_to_rotating(r, v, t)
    rI, vI = impl.rotating_to_inertial(rR, vR, t)
    npt.assert_allclose(rI, r, rtol=1e-14)
    npt.assert_allclose(vI, v, rtol=1e-12, atol=1e-9)
    r0, v0 = impl.inertial_to_rotating(r, v, 0.0)
    npt.assert_allclose(r0, r)
    npt.assert_allclose(v0, v - np.cross([0, 0, OM], r))                      # the ground's velocity removed


def test_enu_basis_and_launch_state():
    for lat, lon in ((0.3, 1.0), (-0.9, -2.0), (1.5, 0.0)):
        B = impl.enu_basis(lat, lon)
        npt.assert_allclose(B @ B.T, np.eye(3), atol=1e-14)
        npt.assert_allclose(np.cross(B[0], B[1]), B[2], atol=1e-14)          # E x N = U
        assert B[1][2] > 0                                                    # north climbs toward the pole
    r0, v0 = impl.launch_state(0.5, 0.2, 800.0, np.radians(30), np.radians(90))
    E, N, U = impl.enu_basis(0.5, 0.2)
    assert np.linalg.norm(r0) == pytest.approx(impl.R_E)
    assert np.linalg.norm(v0) == pytest.approx(800.0)
    assert v0 @ U == pytest.approx(400.0) and v0 @ N == pytest.approx(0.0, abs=1e-9)
    assert v0 @ E > 0                                                         # azimuth 90 deg is east


# ---------------------------------------------------------------- independent check: two frames

@pytest.mark.parametrize("az", [np.radians(30.0), np.radians(250.0)])
def test_rotating_and_inertial_integrations_land_together(az):
    lat, lon = np.radians(40.0), np.radians(-75.0)
    r0, v0 = impl.launch_state(lat, lon, 700.0, np.radians(50), az)
    tr, rr, vr = impl.simulate_rotating(r0, v0)
    ti, ri, vi = impl.simulate_inertial(r0, v0)
    assert tr == pytest.approx(ti, rel=1e-9)
    assert np.linalg.norm(rr - ri) < 1e-3                                     # m, after ~100 s of flight
    npt.assert_allclose(vr, vi, atol=1e-6)
    # Jacobi integral of the rotating frame: v^2/2 - mu/r - |omega x r|^2/2 is conserved
    J = lambda r, v: 0.5 * v @ v - impl.MU_E / np.linalg.norm(r) - 0.5 * np.sum(np.cross([0, 0, OM], r) ** 2)
    assert J(rr, vr) == pytest.approx(J(r0, v0), rel=1e-10)


def test_non_rotating_earth_has_no_cross_range():
    lat, lon = np.radians(49.0), 0.0
    for az in (0.0, np.radians(77.0)):
        r0, v0 = impl.launch_state(lat, lon, 600.0, np.radians(45), az)
        _, r, _ = impl.simulate_rotating(r0, v0, Omega=0.0)
        down, cross = impl.range_and_cross(r, lat, lon, az)
        assert abs(cross) < 1e-6 and down > 0


# ---------------------------------------------------------------- flat Earth, first order

@pytest.mark.parametrize("az", [0.0, np.pi / 2, np.pi])
def test_flat_earth_integration_matches_first_order_closed_form(az):
    lat, V, el = np.radians(49.0), 1000.0, np.radians(45.0)
    v = np.array([np.sin(az) * np.cos(el), np.cos(az) * np.cos(el), np.sin(el)]) * V
    _, r = impl.flat_earth_rotating(v, lat)
    _, r0 = impl.flat_earth_rotating(v, lat, Omega=0.0)
    shift = r[:2] - r0[:2]
    pred = impl.first_order_impact_shift(v, lat)
    # neglected terms are O(Omega T) ~ 1 % of the first-order shift
    assert np.linalg.norm(shift - pred) < 0.02 * np.linalg.norm(pred) + 2.0
    assert np.linalg.norm(pred) > 300.0


def test_deflection_right_in_the_north_left_in_the_south_and_linear_in_omega():
    v = np.array([0.0, 300.0, 0.0])                                          # horizontal, northward
    for lat, sign in ((np.radians(45), 1), (np.radians(-45), -1)):
        d = impl.first_order_deflection(v, lat, 10.0, g=0.0)
        assert sign * d[0] > 0                                               # east is to the right of north
        assert d[0] == pytest.approx(sign * OM * np.sin(abs(lat)) * 300.0 * 100.0)
    d1 = impl.first_order_deflection(v, 0.7, 10.0, Omega=OM)
    d2 = impl.first_order_deflection(v, 0.7, 10.0, Omega=2 * OM)
    npt.assert_allclose(d2, 2 * d1)


def test_northward_shot_at_the_equator_lands_west():
    vu = 600.0
    v = np.array([0.0, 400.0, vu])
    T = 2 * vu / G
    sh = impl.first_order_impact_shift(v, 0.0)
    assert sh[0] == pytest.approx(-OM * T ** 2 * vu / 3, rel=1e-12)
    assert abs(sh[1]) < 1e-9


def test_spherical_coriolis_alone_close_to_flat_earth():
    lat, lon = np.radians(49.0), 0.0
    r0, v0 = impl.launch_state(lat, lon, 1000.0, np.radians(45), 0.0)
    _, r, _ = impl.simulate_rotating(r0, v0)
    _, rn, _ = impl.simulate_rotating(r0, v0, coriolis=False)
    cross = impl.range_and_cross(r, lat, lon, 0.0)[1] - impl.range_and_cross(rn, lat, lon, 0.0)[1]
    v = np.array([0.0, 1000 * np.cos(np.radians(45)), 1000 * np.sin(np.radians(45))])
    pred = impl.first_order_impact_shift(v, lat)[0]
    # spherical flight is ~1.4 % longer (gravity weakens with height, the ground curves away): ~3 % more
    assert pred < cross < 1.06 * pred


# ---------------------------------------------------------------- Foucault pendulum

def test_foucault_closed_form_solves_the_equations():
    L, lat, amp = 67.0, np.radians(48.85), 3.0
    w0, Oz = np.sqrt(G / L), OM * np.sin(lat)
    t, h = np.array([0.0, 5.0, 123.4]), 1e-3
    x, y = impl.foucault_closed_form(t, L, lat, amp)
    xp, yp = impl.foucault_closed_form(t + h, L, lat, amp)
    xm, ym = impl.foucault_closed_form(t - h, L, lat, amp)
    vx, vy = (xp - xm) / (2 * h), (yp - ym) / (2 * h)
    ax, ay = (xp - 2 * x + xm) / h ** 2, (yp - 2 * y + ym) / h ** 2
    npt.assert_allclose(ax, -w0 ** 2 * x + 2 * Oz * vy, atol=1e-6)
    npt.assert_allclose(ay, -w0 ** 2 * y - 2 * Oz * vx, atol=1e-6)
    assert x[0] == pytest.approx(amp) and y[0] == 0.0
    assert abs(vx[0]) < 1e-8 and abs(vy[0]) < 1e-8                           # released from rest on the floor


def test_foucault_integration_matches_closed_form_and_precesses_at_omega_sin_lat():
    L, lat, amp = 67.0, np.radians(48.85), 3.0
    t, x, y = impl.simulate_foucault(L, lat, amp, 600.0, n_out=12001)
    xc, yc = impl.foucault_closed_form(t, L, lat, amp)
    assert np.max(np.hypot(x - xc, y - yc)) < 1e-6
    rate = impl.swing_plane_rate(t, x, y)
    assert rate == pytest.approx(-OM * np.sin(lat), rel=1e-4)                 # clockwise in the north


def test_foucault_limits_equator_pole_and_no_rotation():
    L, amp = 30.0, 1.0
    t, x, y = impl.simulate_foucault(L, 0.0, amp, 300.0, n_out=6001)
    assert abs(impl.swing_plane_rate(t, x, y)) < 1e-9                         # equator: no precession
    t, x, y = impl.simulate_foucault(L, np.radians(-90.0), amp, 300.0, n_out=6001)
    assert impl.swing_plane_rate(t, x, y) == pytest.approx(OM, rel=1e-4)      # south pole: counterclockwise
    t, x, y = impl.simulate_foucault(L, 0.8, amp, 300.0, Omega=0.0, n_out=6001)
    assert np.max(np.abs(y)) < 1e-9


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "rotating_earth.png").stat().st_size > 10_000
