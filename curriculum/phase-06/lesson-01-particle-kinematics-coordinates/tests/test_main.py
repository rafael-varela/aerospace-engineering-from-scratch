"""Public tests for P06.01 Particle Kinematics in Many Coordinates. They assert kinematics, not magic
numbers: the polar and cylindrical formulas against direct differentiation of x = r cos(theta),
y = r sin(theta); straight-line flight has zero acceleration in polar coordinates even though r and
theta both curve; Cartesian <-> cylindrical conversions round-trip; the path description satisfies
a_t^2 + a_n^2 = |a|^2, has orthonormal e_t, e_n, e_b and the helix and circle curvatures; and the
independent check: finite differences of sampled radar data converge to the closed-form trajectory at
second order by both the cylindrical and the Cartesian route, while synthetic noise makes the error
grow again at small steps, with the optimum near the predicted step.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)


def _polar_numeric(rf, thf, t, h=1e-4):
    """Cartesian velocity and acceleration of x = r cos th, y = r sin th by central differences."""
    xy = lambda s: np.array([rf(s) * np.cos(thf(s)), rf(s) * np.sin(thf(s))])
    v = (xy(t + h) - xy(t - h)) / (2 * h)
    a = (xy(t + h) - 2 * xy(t) + xy(t - h)) / h ** 2
    return v, a


# ---------------------------------------------------------------- polar and cylindrical

def test_polar_components_match_direct_differentiation():
    rf = lambda s: 3.0 + 0.5 * s ** 2
    thf = lambda s: 0.4 * s + 0.1 * s ** 3
    t = 1.3
    r, rd, rdd = rf(t), t * 1.0, 1.0
    th, thd, thdd = thf(t), 0.4 + 0.3 * t ** 2, 0.6 * t
    vr, vth, ar, ath = impl.polar_kinematics(r, rd, rdd, th, thd, thdd)
    v_num, a_num = _polar_numeric(rf, thf, t)
    er = np.array([np.cos(th), np.sin(th)])
    et = np.array([-np.sin(th), np.cos(th)])
    npt.assert_allclose([vr, vth], [v_num @ er, v_num @ et], rtol=1e-7)
    npt.assert_allclose([ar, ath], [a_num @ er, a_num @ et], rtol=1e-5)


def test_uniform_circular_motion_is_purely_centripetal():
    R, w = 1500.0, 0.12
    vr, vth, ar, ath = impl.polar_kinematics(R, 0.0, 0.0, 0.7, w, 0.0)
    assert vr == 0.0 and vth == pytest.approx(R * w)
    assert ar == pytest.approx(-R * w ** 2) and ath == pytest.approx(0.0)


def test_straight_line_flight_has_zero_acceleration_in_polar_coordinates():
    # an aircraft passing a radar at constant velocity: r and theta both vary nonlinearly
    t = np.linspace(-60.0, 60.0, 7)
    p0, v = np.array([-4000.0, 3000.0, 2000.0]), np.array([220.0, 15.0, 0.0])
    pos = p0 + t[:, None] * v
    r, rd, rdd, th, thd, thdd, z, zd, zdd = impl.cartesian_to_cylindrical(pos, np.tile(v, (7, 1)), np.zeros((7, 3)))
    assert np.ptp(rd) > 100.0 and np.ptp(thd) > 1e-3            # the rates are far from constant
    vr, vth, ar, ath = impl.polar_kinematics(r, rd, rdd, th, thd, thdd)
    npt.assert_allclose(ar, 0.0, atol=1e-10)                    # r'' = r theta'^2 exactly cancels
    npt.assert_allclose(ath, 0.0, atol=1e-10)                   # r theta'' = -2 r' theta'
    npt.assert_allclose(np.hypot(vr, vth), np.hypot(220.0, 15.0), rtol=1e-12)


def test_cylindrical_round_trip_on_a_general_motion():
    rng = np.random.default_rng(3)
    pos = rng.uniform(-5000, 5000, (20, 3))
    vel = rng.uniform(-300, 300, (20, 3))
    acc = rng.uniform(-50, 50, (20, 3))
    r, rd, rdd, th, thd, thdd, z, zd, zdd = impl.cartesian_to_cylindrical(pos, vel, acc)
    v_cyl, a_cyl = impl.cylindrical_kinematics(r, rd, rdd, th, thd, thdd, z, zd, zdd)
    assert v_cyl.shape == (20, 3) and a_cyl.shape == (20, 3)
    npt.assert_allclose(impl.cylindrical_to_cartesian(v_cyl, th), vel, rtol=1e-10, atol=1e-9)
    npt.assert_allclose(impl.cylindrical_to_cartesian(a_cyl, th), acc, rtol=1e-10, atol=1e-9)
    npt.assert_allclose(r * np.cos(th), pos[:, 0], rtol=1e-12)
    # the speed is the same in every basis
    npt.assert_allclose(np.linalg.norm(v_cyl, axis=1), np.linalg.norm(vel, axis=1), rtol=1e-12)


def test_cylindrical_to_cartesian_is_a_rotation():
    th = 0.9
    comp = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    E = impl.cylindrical_to_cartesian(comp, th)
    npt.assert_allclose(E @ E.T, np.eye(3), atol=1e-15)
    assert np.linalg.det(E) == pytest.approx(1.0)
    npt.assert_allclose(E[0], [np.cos(th), np.sin(th), 0.0])


# ---------------------------------------------------------------- path coordinates

def test_path_decomposition_and_orthonormal_triad():
    rng = np.random.default_rng(7)
    vel = rng.normal(0, 100, (50, 3))
    acc = rng.normal(0, 10, (50, 3))
    pk = impl.path_kinematics(vel, acc)
    npt.assert_allclose(pk["a_t"] ** 2 + pk["a_n"] ** 2, np.sum(acc ** 2, axis=1), rtol=1e-10)
    assert np.all(pk["a_n"] >= 0)
    for a, b in (("e_t", "e_t"), ("e_n", "e_n"), ("e_b", "e_b")):
        npt.assert_allclose(np.sum(pk[a] * pk[b], axis=1), 1.0, rtol=1e-12)
    for a, b in (("e_t", "e_n"), ("e_t", "e_b"), ("e_n", "e_b")):
        npt.assert_allclose(np.sum(pk[a] * pk[b], axis=1), 0.0, atol=1e-12)
    # a is in the osculating plane: no binormal component; and a = a_t e_t + a_n e_n
    npt.assert_allclose(np.sum(acc * pk["e_b"], axis=1), 0.0, atol=1e-10)
    npt.assert_allclose(pk["a_t"][:, None] * pk["e_t"] + pk["a_n"][:, None] * pk["e_n"], acc, atol=1e-10)
    npt.assert_allclose(pk["kappa"] * pk["speed"] ** 2, pk["a_n"], rtol=1e-12)


def test_path_kinematics_limits_line_circle_and_speed_change():
    pk = impl.path_kinematics(np.array([[100.0, 0, 0]]), np.array([[3.0, 0, 0]]))
    assert pk["a_t"][0] == pytest.approx(3.0) and pk["a_n"][0] == pytest.approx(0.0)
    assert np.isinf(pk["rho"][0])
    V, R = 250.0, 3000.0
    pk = impl.path_kinematics(np.array([[0.0, V, 0]]), np.array([[-V ** 2 / R, 0, 0]]))
    assert pk["rho"][0] == pytest.approx(R) and pk["a_t"][0] == pytest.approx(0.0, abs=1e-12)
    npt.assert_allclose(pk["e_n"][0], [-1.0, 0.0, 0.0])
    npt.assert_allclose(pk["e_b"][0], [0.0, 0.0, 1.0])          # counterclockwise turn: binormal up


def test_helix_closed_form_is_self_consistent():
    t = np.linspace(0, 60, 601)
    pos, vel, acc = impl.helix_turn(t)
    tm, v_fd, a_fd = impl.central_differences(t, pos, 1)
    # truncation: (h^2/6) V omega^2 = 3.3e-3 m/s on v, (h^2/12) R omega^4 = 1.7e-4 m/s^2 on a
    npt.assert_allclose(v_fd, vel[1:-1], atol=5e-3)
    npt.assert_allclose(a_fd, acc[1:-1], atol=3e-4)
    kap = impl.helix_curvature(impl.R_TURN, impl.OMEGA_TURN, impl.CLIMB)
    pk = impl.path_kinematics(vel, acc)
    npt.assert_allclose(pk["kappa"], kap, rtol=1e-12)
    npt.assert_allclose(pk["a_t"], 0.0, atol=1e-10)
    npt.assert_allclose(pk["a_n"], impl.R_TURN * impl.OMEGA_TURN ** 2, rtol=1e-12)
    # zero climb: the helix is a circle of radius R
    assert impl.helix_curvature(2000.0, 0.1, 0.0) == pytest.approx(1 / 2000.0)


# ---------------------------------------------------------------- independent check: sampled data

def test_central_differences_exact_for_quadratics_and_shapes():
    t = np.linspace(0.0, 2.0, 21)
    f = np.stack([3 * t ** 2 - t + 2, -t ** 2], axis=-1)
    tm, fd, fdd = impl.central_differences(t, f, 2)
    assert tm.shape == (17,) and fd.shape == (17, 2)
    npt.assert_allclose(fd[:, 0], 6 * tm - 1, atol=1e-12)
    npt.assert_allclose(fdd[:, 0], 6.0, atol=1e-10)
    npt.assert_allclose(fdd[:, 1], -2.0, atol=1e-10)
    with pytest.raises(ValueError):
        impl.central_differences(t, f, 11)


@pytest.mark.parametrize("route", ["track_from_radar", "track_cartesian"])
def test_sampled_radar_converges_at_second_order(route):
    dt = 0.1
    t = np.arange(0.0, 80.0 + dt / 2, dt)
    pos, vel, acc = impl.helix_turn(t)
    r, th, z = impl.radar_samples(pos)
    ev, ea = [], []
    for k in (40, 20, 10, 5):
        tm, v, a = getattr(impl, route)(t, r, th, z, k)
        assert len(tm) == len(t) - 2 * k
        ev.append(np.max(np.abs(v - vel[k:-k])))
        ea.append(np.max(np.abs(a - acc[k:-k])))
    assert np.all(np.abs(impl.observed_order(ev) - 2) < 0.1)
    assert np.all(np.abs(impl.observed_order(ea) - 2) < 0.1)
    assert ea[-1] < 1e-3 * impl.R_TURN * impl.OMEGA_TURN ** 2


def test_two_routes_agree_and_recover_the_turn():
    dt = 0.1
    t = np.arange(0.0, 60.0 + dt / 2, dt)
    pos, vel, acc = impl.helix_turn(t, phase=1.0)
    r, th, z = impl.radar_samples(pos)
    _, v1, a1 = impl.track_from_radar(t, r, th, z, 2)
    _, v2, a2 = impl.track_cartesian(t, r, th, z, 2)
    npt.assert_allclose(a1, a2, atol=2e-3)                      # both O(h^2) from the truth, h = 0.2 s
    pk = impl.path_kinematics(v1, a1)
    kap = impl.helix_curvature(impl.R_TURN, impl.OMEGA_TURN, impl.CLIMB)
    npt.assert_allclose(pk["rho"], 1 / kap, rtol=1e-3)
    npt.assert_allclose(pk["a_n"], impl.R_TURN * impl.OMEGA_TURN ** 2, rtol=1e-3)


def test_azimuth_is_unwrapped_through_pi():
    # a track crossing the negative x axis (theta jumps from +pi to -pi without unwrapping)
    t = np.arange(0.0, 40.0, 0.1)
    pos, vel, acc = impl.helix_turn(t, centre=(-6000.0, 0.0), phase=0.0)
    r, th, z = impl.radar_samples(pos)
    assert np.max(np.abs(np.diff(th))) < 0.1
    _, v, a = impl.track_from_radar(t, r, th, z, 2)
    npt.assert_allclose(a, acc[2:-2], atol=5e-3)


def test_noise_makes_small_steps_worse_and_the_optimum_is_near_prediction():
    dt = 0.1
    t = np.arange(0.0, 120.0 + dt / 2, dt)
    pos, vel, acc = impl.helix_turn(t)
    strides = np.array([5, 10, 20, 40, 60, 80, 140])
    err = impl.accel_error_vs_stride(t, pos, acc, strides, 5.0, 1e-3, 5.0, seed=11)
    assert err[0] > 20 * err.min()                              # noise / h^2 dominates at h = 0.5 s
    assert err[-1] > 1.5 * err.min()                             # truncation h^2 dominates at h = 14 s
    h_opt, e_opt = impl.optimal_step(8.0, impl.R_TURN * impl.OMEGA_TURN ** 4)
    h_best = strides[np.argmin(err)] * dt
    # one seeded record: the measured optimum is within a factor 2 of the model (the error curve is flat there)
    assert h_opt / 2 < h_best < 2 * h_opt
    assert e_opt / 3 < err.min() < 3 * e_opt


def test_noise_is_seeded_and_reproducible():
    pos, _, _ = impl.helix_turn(np.arange(0, 10, 0.1))
    a = impl.radar_samples(pos, 5.0, 1e-3, 5.0, np.random.default_rng(4))
    b = impl.radar_samples(pos, 5.0, 1e-3, 5.0, np.random.default_rng(4))
    npt.assert_array_equal(a[0], b[0])
    clean = impl.radar_samples(pos)
    assert 3.0 < np.std(a[0] - clean[0]) < 7.0                  # 100 samples of sigma = 5 m


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "radar_track.png").stat().st_size > 10_000
