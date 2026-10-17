"""Public tests for P04.01 Kinematics & Newton's Laws. They assert physics, not magic numbers:
the drag-free limit of the simulation against the constant-acceleration closed forms, the purely
vertical and purely horizontal limits against their closed forms, energy and symmetry properties,
convergence with the solver tolerance, the wind triangle as vector addition, and the recovery of
a known terminal speed from SYNTHETIC frame-quantised drop data.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-04/parameterized-cases/.
"""
import math

import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import solve_ivp

from conftest import load_impl

impl = load_impl(__file__)

G = impl.G0
# the ping-pong ball's drag parameter, computed here (not through impl) so that collection never depends on learner code
K_PP = impl.RHO_AIR * impl.PING_PONG["cd"] * math.pi * impl.PING_PONG["diameter"] ** 2 / 4 / (2 * impl.PING_PONG["mass"])


# ---------------------------------------------------------------- Newton's second law and the drag parameter

def test_drag_parameter_scales_like_area_over_mass():
    k = impl.drag_parameter(0.0027, 0.040, cd=0.5, rho=1.2)
    assert k == pytest.approx(1.2 * 0.5 * math.pi * 0.02 ** 2 / (2 * 0.0027))
    assert impl.drag_parameter(0.0054, 0.040, cd=0.5, rho=1.2) == pytest.approx(k / 2)      # double the mass
    assert impl.drag_parameter(0.0027, 0.080, cd=0.5, rho=1.2) == pytest.approx(4 * k)      # double the diameter
    assert impl.drag_parameter(0.0027, 0.040, cd=1.0, rho=1.2) == pytest.approx(2 * k)


def test_terminal_speed_balances_drag_and_weight():
    v_t = impl.terminal_speed(K_PP)
    assert K_PP * v_t ** 2 == pytest.approx(G, rel=1e-12)          # drag acceleration equals g
    assert 7.0 < v_t < 10.0                                          # a ping-pong ball: about 8 m/s
    assert math.isinf(impl.terminal_speed(0.0))
    # the rhs is zero acceleration for a body falling at v_t
    rhs = impl.projectile_rhs(0.0, np.array([0.0, 10.0, 0.0, -v_t]), K_PP)
    npt.assert_allclose(rhs, [0.0, -v_t, 0.0, 0.0], atol=1e-12)


def test_rhs_without_drag_is_free_fall_and_drag_opposes_velocity():
    y = np.array([1.0, 2.0, 3.0, -4.0])
    npt.assert_allclose(impl.projectile_rhs(0.0, y, 0.0), [3.0, -4.0, 0.0, -G])
    a = impl.projectile_rhs(0.0, y, 0.2)[2:]
    drag = a - np.array([0.0, -G])
    v = y[2:]
    assert np.dot(drag, v) < 0                                       # opposes motion
    assert abs(drag[0] * v[1] - drag[1] * v[0]) < 1e-12              # and is parallel to it
    assert np.linalg.norm(drag) == pytest.approx(0.2 * np.linalg.norm(v) ** 2)


# ---------------------------------------------------------------- closed forms

def test_vacuum_closed_forms_are_consistent():
    v0, ang = 12.0, 37.0
    T = impl.vacuum_time_of_flight(v0, ang)
    R = impl.vacuum_range(v0, ang)
    y = impl.vacuum_state(T, v0, ang)
    assert y.shape == (4,)
    assert y[1] == pytest.approx(0.0, abs=1e-12)
    assert y[0] == pytest.approx(R)
    assert y[3] == pytest.approx(-v0 * math.sin(math.radians(ang)))       # lands at the launch speed, mirrored
    Y = impl.vacuum_state(np.linspace(0, T, 7), v0, ang)
    assert Y.shape == (7, 4)
    npt.assert_allclose(Y[:, 2], v0 * math.cos(math.radians(ang)))       # vx constant
    assert impl.vacuum_range(v0, 45.0) == pytest.approx(v0 ** 2 / G)
    assert impl.vacuum_range(v0, 30.0) == pytest.approx(impl.vacuum_range(v0, 60.0))   # complementary angles


def test_vertical_drop_closed_form_satisfies_the_ode():
    v_t = 8.0
    t = np.linspace(0.0, 3.0, 61)
    s, v = impl.vertical_drop_state(t, v_t)
    # ds/dt = v and dv/dt = g (1 - v^2/v_t^2), checked by central differences
    ds = np.gradient(s, t)
    npt.assert_allclose(ds[1:-1], v[1:-1], rtol=2e-3)
    dv = np.gradient(v, t)
    npt.assert_allclose(dv[1:-1], G * (1 - v[1:-1] ** 2 / v_t ** 2), rtol=5e-3, atol=1e-3)
    assert v[-1] < v_t and v[-1] > 0.99 * v_t                             # approaching terminal speed
    s_inf, v_inf = impl.vertical_drop_state(t, math.inf)
    npt.assert_allclose(s_inf, 0.5 * G * t ** 2)
    npt.assert_allclose(v_inf, G * t)


def test_vertical_drop_time_inverts_the_distance_law_and_has_the_vacuum_limit():
    for h in (0.3, 1.0, 2.0, 50.0):
        t = impl.vertical_drop_time(h, 8.0)
        s, _ = impl.vertical_drop_state(t, 8.0)
        assert float(s) == pytest.approx(h, rel=1e-10)
        assert t > math.sqrt(2 * h / G)                                    # drag lengthens every fall
    assert impl.vertical_drop_time(1.5, math.inf) == pytest.approx(math.sqrt(3.0 / G))
    assert impl.vertical_drop_time(1.5, 1e6) == pytest.approx(math.sqrt(3.0 / G), rel=1e-9)      # no cancellation
    # first-order expansion t ~ t0 (1 + g h / (6 v_t^2)) for weak drag
    h, v_t = 1.0, 40.0
    t0 = math.sqrt(2 * h / G)
    assert impl.vertical_drop_time(h, v_t) == pytest.approx(t0 * (1 + G * h / (6 * v_t ** 2)), rel=1e-5)
    # a very long fall: distance ~ v_t t, no overflow
    assert impl.vertical_drop_time(1e6, 8.0) == pytest.approx(1e6 / 8.0, rel=1e-3)


def test_horizontal_decay_closed_form_satisfies_the_ode():
    v0, k = 6.0, 0.14
    t = np.linspace(0.0, 2.0, 41)
    x, v = impl.horizontal_decay_state(t, v0, k)
    npt.assert_allclose(np.gradient(x, t)[1:-1], v[1:-1], rtol=2e-3)
    npt.assert_allclose(np.gradient(v, t)[1:-1], -k * v[1:-1] ** 2, rtol=5e-3)
    assert x[-1] == pytest.approx(math.log(1 + k * v0 * 2.0) / k)
    x0, v00 = impl.horizontal_decay_state(t, v0, 0.0)
    npt.assert_allclose(x0, v0 * t)
    npt.assert_allclose(v00, v0)


# ---------------------------------------------------------------- the simulation against its limits

def test_simulation_without_drag_reproduces_the_vacuum_formulas():
    v0, ang = 9.0, 50.0
    tr = impl.simulate_projectile(v0, ang, 0.0, rtol=1e-11)
    assert tr.time_of_flight == pytest.approx(impl.vacuum_time_of_flight(v0, ang), rel=1e-8)
    assert tr.range == pytest.approx(impl.vacuum_range(v0, ang), rel=1e-8)
    assert tr.impact_speed == pytest.approx(v0, rel=1e-8)
    assert tr.apex_height == pytest.approx((v0 * math.sin(math.radians(ang))) ** 2 / (2 * G), rel=1e-4)
    Y = impl.vacuum_state(tr.t, v0, ang)
    npt.assert_allclose(tr.x, Y[:, 0], atol=1e-7)
    npt.assert_allclose(tr.z, Y[:, 1], atol=1e-7)


def test_simulation_of_a_drop_matches_the_vertical_closed_form():
    v_t = impl.terminal_speed(K_PP)
    tr = impl.simulate_projectile(0.0, 90.0, K_PP, z0=2.0, rtol=1e-12)
    assert tr.time_of_flight == pytest.approx(impl.vertical_drop_time(2.0, v_t), rel=1e-9)
    assert tr.range == pytest.approx(0.0, abs=1e-12)
    s, v = impl.vertical_drop_state(tr.t, v_t)
    npt.assert_allclose(2.0 - tr.z, s, atol=1e-8)
    npt.assert_allclose(-tr.vz, v, atol=1e-8)
    assert tr.impact_speed < v_t


def test_simulation_horizontal_limit_and_solve_ivp_oracle():
    # launch horizontally from a great height: for a short time x(t) follows the pure-drag decay
    v0, k = 6.0, K_PP
    tr = impl.simulate_projectile(v0, 0.0, k, z0=1000.0, rtol=1e-11, max_time=200.0)
    t_short = 0.05
    i = np.searchsorted(tr.t, t_short)
    x_cf, _ = impl.horizontal_decay_state(tr.t[:i], v0, k)
    npt.assert_allclose(tr.x[:i], x_cf, rtol=2e-3)
    # and the whole trajectory agrees with an independent solve_ivp (DOP853) call
    def rhs(t, y):
        s = math.hypot(y[2], y[3])
        return [y[2], y[3], -k * s * y[2], -G - k * s * y[3]]
    ref = solve_ivp(rhs, (0, tr.time_of_flight), [0, 1000.0, v0, 0.0], method="DOP853", rtol=1e-12, atol=1e-12,
                    t_eval=tr.t)
    npt.assert_allclose(tr.x, ref.y[0], atol=1e-6)
    npt.assert_allclose(tr.z, ref.y[1], atol=1e-6)


def test_drag_shortens_range_and_flight_and_removes_the_symmetry():
    v0, ang = 6.0, 45.0
    vac = impl.simulate_projectile(v0, ang, 0.0)
    dr = impl.simulate_projectile(v0, ang, K_PP)
    assert dr.range < vac.range
    assert dr.time_of_flight < vac.time_of_flight          # a tossed ball with drag lands sooner (lower apex)
    assert dr.apex_height < vac.apex_height
    assert dr.impact_speed < v0                             # drag dissipates kinetic energy
    assert vac.impact_speed == pytest.approx(v0, rel=1e-6)
    # with drag the descent is steeper than the ascent
    impact_angle = math.degrees(math.atan2(-dr.vz[-1], dr.vx[-1]))
    assert impact_angle > ang
    # the apex comes before the midpoint of the flight
    assert dr.t[np.argmax(dr.z)] < 0.5 * dr.time_of_flight
    # range is monotone in k
    ranges = [impl.drag_range(v0, ang, kk) for kk in (0.0, 0.05, 0.1, 0.2, 0.4)]
    assert all(a > b for a, b in zip(ranges, ranges[1:]))


def test_range_converges_with_tolerance():
    r = [impl.simulate_projectile(6.0, 45.0, K_PP, rtol=rt).range for rt in (1e-4, 1e-7, 1e-10)]
    assert abs(r[1] - r[2]) < abs(r[0] - r[2])
    assert abs(r[1] - r[2]) < 1e-6


def test_optimal_angle_is_45_without_drag_and_lower_with_drag():
    a0, r0 = impl.optimal_angle(6.0, 0.0)
    assert a0 == pytest.approx(45.0, abs=0.05)
    assert r0 == pytest.approx(36.0 / G, rel=1e-5)
    a1, r1 = impl.optimal_angle(6.0, K_PP)
    assert 35.0 < a1 < 45.0
    assert r1 >= impl.drag_range(6.0, 45.0, K_PP)
    assert r1 >= impl.drag_range(6.0, a1 - 2.0, K_PP) and r1 >= impl.drag_range(6.0, a1 + 2.0, K_PP)
    a2, _ = impl.optimal_angle(6.0, 4 * K_PP)
    assert a2 < a1                                          # more drag, flatter optimum


# ---------------------------------------------------------------- relative motion

def test_wind_triangle_is_vector_addition():
    assert impl.wind_triangle(60.0, 90.0, 0.0, 0.0) == pytest.approx((60.0, 90.0))
    gs, trk = impl.wind_triangle(60.0, 90.0, 15.0, 90.0)        # headwind from the east
    assert (gs, trk) == pytest.approx((45.0, 90.0))
    gs, trk = impl.wind_triangle(60.0, 90.0, 15.0, 270.0)       # tailwind
    assert (gs, trk) == pytest.approx((75.0, 90.0))
    gs, trk = impl.wind_triangle(60.0, 0.0, 15.0, 270.0)        # pure crosswind from the west, heading north
    assert gs == pytest.approx(math.hypot(60.0, 15.0))
    assert trk == pytest.approx(math.degrees(math.atan2(15.0, 60.0)))   # drifts east
    # symmetry: mirror heading and wind about north, track mirrors
    gs1, t1 = impl.wind_triangle(50.0, 30.0, 10.0, 100.0)
    gs2, t2 = impl.wind_triangle(50.0, 330.0, 10.0, 260.0)
    assert gs1 == pytest.approx(gs2)
    assert t1 == pytest.approx((360.0 - t2) % 360.0)


def test_heading_for_track_inverts_the_wind_triangle():
    for tas, trk, w, wf in [(60.0, 90.0, 15.0, 30.0), (80.0, 200.0, 25.0, 315.0), (50.0, 0.0, 20.0, 90.0)]:
        hd, gs = impl.heading_for_track(tas, trk, w, wf)
        gs2, trk2 = impl.wind_triangle(tas, hd, w, wf)
        assert gs2 == pytest.approx(gs, rel=1e-10)
        assert trk2 == pytest.approx(trk, abs=1e-8)
    hd, gs = impl.heading_for_track(60.0, 90.0, 15.0, 90.0)       # straight headwind: no correction
    assert (hd, gs) == pytest.approx((90.0, 45.0))
    with pytest.raises(ValueError):
        impl.heading_for_track(10.0, 90.0, 15.0, 0.0)              # crosswind exceeds airspeed


# ---------------------------------------------------------------- the experiment pipeline on SYNTHETIC data

def test_predicted_fall_times_order_and_limits():
    hs = np.array([0.6, 1.2, 1.8])
    t_vac, t_drag = impl.predict_fall_times(hs, K_PP)
    npt.assert_allclose(t_vac, np.sqrt(2 * hs / G))
    assert np.all(t_drag > t_vac)
    assert np.all(np.diff(t_drag - t_vac) > 0)                        # the excess grows with height
    t_vac0, t_drag0 = impl.predict_fall_times(hs, 0.0)
    npt.assert_allclose(t_drag0, t_vac0)


def test_fit_recovers_the_terminal_speed_from_synthetic_drops():
    # SYNTHETIC frame-quantised data at 120 fps, consistent with the model by construction; the fitted
    # terminal speed must agree with the one used to generate it within 2.5 standard uncertainties,
    # and the recovered C_D with the generating value within 25 % (quantisation at 120 fps is ~3 ms
    # against a drag excess of 5-26 ms, so the fit is well determined but not exact).
    hs = np.array([0.6, 0.9, 1.2, 1.5, 1.8])
    v_true = impl.terminal_speed(K_PP)
    u_frame = impl.frame_timing(1, 120.0)[1]
    for seed in (1, 2, 3):
        synth = impl.make_synthetic_drop_data(K_PP, hs, 120.0, n_drops=8, seed=seed)
        means, us = map(np.array, zip(*[impl.reduce_timings(s, u_frame) for s in synth]))
        v_fit, u_v = impl.fit_terminal_speed(hs, means, us)
        assert abs(v_fit - v_true) < 2.5 * u_v + 0.05
        assert u_v < 0.6
        cd = impl.drag_coefficient_from_terminal_speed(v_fit, impl.PING_PONG["mass"], impl.PING_PONG["diameter"])
        assert cd == pytest.approx(impl.PING_PONG["cd"], rel=0.25)
    # the inversion is exact on the nominal value
    assert impl.drag_coefficient_from_terminal_speed(v_true, impl.PING_PONG["mass"], impl.PING_PONG["diameter"]) == pytest.approx(impl.PING_PONG["cd"], rel=1e-12)


def test_fit_on_exact_times_is_exact_and_a_dense_object_gives_a_large_terminal_speed():
    hs = np.array([0.5, 1.0, 1.5, 2.0])
    _, t_exact = impl.predict_fall_times(hs, K_PP)
    v_fit, u_v = impl.fit_terminal_speed(hs, t_exact)
    assert v_fit == pytest.approx(impl.terminal_speed(K_PP), rel=1e-6)
    assert u_v < 1e-3
    k_bb = impl.drag_parameter(**impl.BEAN_BAG)
    _, t_bb = impl.predict_fall_times(hs, k_bb)
    v_bb, _ = impl.fit_terminal_speed(hs, t_bb)
    assert v_bb > 2 * impl.terminal_speed(K_PP)


def test_synthetic_data_is_consistent_with_the_drag_model_and_not_with_vacuum_at_1p8m():
    hs = np.array([0.6, 1.2, 1.8])
    u_frame = impl.frame_timing(1, 60.0)[1]
    t_vac, t_drag = impl.predict_fall_times(hs, K_PP)
    n_bad = 0
    for seed in range(4):
        synth = impl.make_synthetic_drop_data(K_PP, hs, 60.0, n_drops=6, seed=seed)
        for s, tv, td in zip(synth, t_vac, t_drag):
            m, u = impl.reduce_timings(s, u_frame)
            assert impl.normalized_error(m, u, td, 0.0) < 2.5           # consistent with the generating model
        m, u = impl.reduce_timings(synth[-1], u_frame)
        n_bad += impl.normalized_error(m, u, t_vac[-1], 0.0) > 2.0
    assert n_bad >= 3                                                    # vacuum is rejected at 1.8 m in most seeds
    assert impl.normalized_error(1.0, 0.1, 1.0, 0.1) == 0.0
    assert impl.normalized_error(1.0, 0.0, 1.2, 0.0) == math.inf


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "projectile_drag.png").exists()
