"""Public tests for P04.04 Rotation & Angular Momentum. They assert physics, not magic numbers: the
inertia formulas and the parallel-axis theorem, conservation of the inertial angular momentum and the
kinetic energy in torque-free motion, stable spin about the extreme axes and exponential growth about
the intermediate axis at the predicted rate, the body-frame nutation of a symmetric spinner, the
reaction-wheel slew against its closed form with total angular momentum zero, and the heavy top's
steady precession, fast-top limit, three invariants and nutation turning point.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-04/parameterized-cases/.
"""
import math

import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

G = 9.80665
Q_ID = np.array([1.0, 0.0, 0.0, 0.0])
# the lesson's top: 100 g disk, r 30 mm, 4 mm thick, centre 40 mm from the pivot (computed here, not by impl)
M_TOP, R_TOP, L_TOP, TH_TOP = 0.10, 0.030, 0.040, 0.004
I1_TOP = M_TOP * (3 * R_TOP ** 2 + TH_TOP ** 2) / 12 + M_TOP * L_TOP ** 2
I3_TOP = 0.5 * M_TOP * R_TOP ** 2
I_TOP = np.diag([I1_TOP, I1_TOP, I3_TOP])
I_SAT = np.diag([260.0, 208.0, 135.0])


def _gravity(t, q, w):
    return impl.gravity_torque_body(q, M_TOP, L_TOP)


# ---------------------------------------------------------------- inertia

def test_inertia_formulas_and_their_limits():
    m, r, L = 2.0, 0.1, 0.4
    Ic = impl.inertia_solid_cylinder(m, r, L)
    npt.assert_allclose(np.diag(Ic), [m * (3 * r * r + L * L) / 12] * 2 + [m * r * r / 2], rtol=1e-14)
    npt.assert_allclose(Ic, np.diag(np.diag(Ic)), atol=0)                                # principal axes
    rod = impl.inertia_solid_cylinder(m, 0.0, L)
    assert rod[0, 0] == pytest.approx(m * L * L / 12) and rod[2, 2] == 0.0             # thin rod
    disk = impl.inertia_solid_cylinder(m, r, 0.0)
    assert disk[0, 0] + disk[1, 1] == pytest.approx(disk[2, 2])                        # perpendicular-axis theorem
    Ib = impl.inertia_box(6.0, 1.0, 2.0, 3.0)
    npt.assert_allclose(np.diag(Ib), [6.5, 5.0, 2.5], rtol=1e-14)
    cube = impl.inertia_box(12.0, 1.0, 1.0, 1.0)
    npt.assert_allclose(cube, 2.0 * np.eye(3), rtol=1e-14)


def test_parallel_axis_theorem_in_tensor_form():
    Ic = impl.inertia_box(3.0, 0.2, 0.5, 0.9)
    d = np.array([0.3, -0.2, 0.7])
    Ip = impl.parallel_axis(Ic, 3.0, d)
    npt.assert_allclose(Ip, Ip.T, atol=1e-15)                                          # symmetric
    n = np.array([1.0, 0.0, 0.0])
    # moment about a line through the new point = moment about the parallel line through the centre + m r_perp^2
    r_perp2 = d @ d - (d @ n) ** 2
    assert n @ Ip @ n == pytest.approx(n @ Ic @ n + 3.0 * r_perp2, rel=1e-13)
    assert Ip[0, 1] == pytest.approx(-3.0 * d[0] * d[1])                               # products of inertia appear
    npt.assert_allclose(impl.parallel_axis(Ic, 3.0, np.zeros(3)), Ic, atol=0)
    assert np.all(np.linalg.eigvalsh(Ip) > 0)


# ---------------------------------------------------------------- energy, momentum, kinematics

def test_energy_and_momentum_of_a_known_state():
    w = np.array([0.1, -0.2, 0.3])
    assert impl.rotational_energy(I_SAT, w) == pytest.approx(0.5 * (260 * 0.01 + 208 * 0.04 + 135 * 0.09))
    stack = np.vstack([w, 2 * w])
    npt.assert_allclose(impl.rotational_energy(I_SAT, stack), np.array([1.0, 4.0]) * impl.rotational_energy(I_SAT, w), rtol=1e-14)
    npt.assert_allclose(impl.angular_momentum_inertial(Q_ID, w, I_SAT), I_SAT @ w, rtol=1e-14)
    q90 = impl.quat_from_axis_angle([0, 0, 1], math.pi / 2)                            # body x points along inertial y
    npt.assert_allclose(impl.angular_momentum_inertial(q90, np.array([1.0, 0, 0]), I_SAT), [0, 260, 0], atol=1e-12)
    h = np.array([0.0, 0.0, 5.0])
    npt.assert_allclose(impl.angular_momentum_inertial(Q_ID, w, I_SAT, h), I_SAT @ w + h, rtol=1e-14)


def test_euler_rates_reduce_to_known_cases():
    w = np.array([0.3, -0.1, 0.2])
    npt.assert_allclose(impl.euler_rates(w, 50.0 * np.eye(3), np.zeros(3)), 0.0, atol=1e-15)   # sphere: no gyroscopic term
    tau = np.array([1.0, 2.0, 3.0])
    npt.assert_allclose(impl.euler_rates(np.zeros(3), I_SAT, tau), tau / np.diag(I_SAT), rtol=1e-14)
    wd = impl.euler_rates(w, I_SAT, np.zeros(3))
    I1, I2, I3 = np.diag(I_SAT)
    npt.assert_allclose(wd, [(I2 - I3) * w[1] * w[2] / I1, (I3 - I1) * w[2] * w[0] / I2, (I1 - I2) * w[0] * w[1] / I3], rtol=1e-12)
    # the wheel torque appears on the body with the opposite sign
    hd = np.array([0.0, 0.0, 0.5])
    npt.assert_allclose(impl.euler_rates(np.zeros(3), I_SAT, np.zeros(3), np.zeros(3), hd), -hd / np.diag(I_SAT), rtol=1e-14)


def test_quaternion_rate_preserves_the_norm_and_rotates_about_the_axis():
    rng = np.random.default_rng(4)
    for _ in range(5):
        q = rng.normal(size=4); q /= np.linalg.norm(q)
        w = rng.normal(size=3)
        assert q @ impl.quaternion_rate(q, w) == pytest.approx(0.0, abs=1e-14)       # d|q|^2/dt = 0
    # constant rate about body z from identity: q(t) = [cos(wt/2), 0, 0, sin(wt/2)], so dq/dt at 0 is [0, 0, 0, w/2]
    npt.assert_allclose(impl.quaternion_rate(Q_ID, np.array([0, 0, 0.8])), [0, 0, 0, 0.4], atol=1e-15)


# ---------------------------------------------------------------- torque-free motion

def test_torque_free_motion_conserves_inertial_momentum_and_energy():
    w0 = np.array([0.2, -0.4, 0.3])
    q0 = impl.quat_from_axis_angle([1, 2, 3], 0.7)
    rot = impl.simulate_rigid_body(I_SAT, w0, q0, 60.0, n_samples=301)
    L = impl.angular_momentum_inertial(rot.q, rot.w, I_SAT)
    npt.assert_allclose(L, np.broadcast_to(L[0], L.shape), atol=1e-8 * np.linalg.norm(L[0]))
    T = impl.rotational_energy(I_SAT, rot.w)
    assert np.ptp(T) < 1e-9 * T[0]
    npt.assert_allclose(np.linalg.norm(rot.q, axis=1), 1.0, atol=1e-12)
    assert np.ptp(np.linalg.norm(rot.w, axis=1)) > 1e-3                               # yet the body rates change


def test_spin_about_the_extreme_axes_is_stable_and_wobbles_at_the_predicted_frequency():
    eps = 1e-4
    for ax in (0, 2):
        w0 = np.full(3, eps); w0[ax] = 0.5
        rot = impl.simulate_rigid_body(I_SAT, w0, Q_ID, 100.0, n_samples=2001, rtol=1e-10)
        others = [i for i in range(3) if i != ax]
        assert np.max(np.abs(rot.w[:, others])) < 10 * eps
        f = impl.principal_axis_wobble_frequency(I_SAT, ax, 0.5)
        # dominant frequency of the small rate oscillation
        sig = rot.w[:, others[0]] - rot.w[:, others[0]].mean()
        freqs = np.fft.rfftfreq(sig.size, rot.t[1] - rot.t[0]) * 2 * math.pi
        f_meas = freqs[np.argmax(np.abs(np.fft.rfft(sig)))]
        assert f_meas == pytest.approx(f, abs=2 * math.pi / 100.0)                    # one FFT bin
    with pytest.raises(ValueError):
        impl.principal_axis_wobble_frequency(I_SAT, 1, 0.5)


def test_intermediate_axis_perturbation_grows_at_the_predicted_rate():
    w2 = 0.5
    sigma = impl.intermediate_axis_growth_rate(I_SAT, w2)
    I1, I2, I3 = np.diag(I_SAT)
    assert sigma == pytest.approx(w2 * math.sqrt((I1 - I2) * (I2 - I3) / (I1 * I3)))   # order of the axes does not matter
    rot = impl.simulate_rigid_body(I_SAT, np.array([1e-7, w2, 1e-7]), Q_ID, 60.0, n_samples=601)
    pert = np.hypot(rot.w[:, 0], rot.w[:, 2])
    win = (rot.t > 15.0) & (pert < 1e-3)
    rate = np.polyfit(rot.t[win], np.log(pert[win]), 1)[0]
    assert rate == pytest.approx(sigma, rel=0.02)


def test_symmetric_body_nutates_at_the_body_frame_rate():
    I = np.diag([800.0, 800.0, 1200.0])
    rot = impl.simulate_rigid_body(I, np.array([0.01, 0.0, 1.0]), Q_ID, 20.0, n_samples=801)
    lam = impl.symmetric_nutation_rate(800.0, 1200.0, 1.0)
    assert lam == pytest.approx(0.5)
    ph = np.unwrap(np.arctan2(rot.w[:, 1], rot.w[:, 0]))
    assert np.polyfit(rot.t, ph, 1)[0] == pytest.approx(lam, rel=1e-6)
    npt.assert_allclose(rot.w[:, 2], 1.0, atol=1e-12)                                  # w3 constant for I1 = I2
    npt.assert_allclose(np.hypot(rot.w[:, 0], rot.w[:, 1]), 0.01, rtol=1e-8)
    prolate = impl.symmetric_nutation_rate(1200.0, 800.0, 1.0)
    assert prolate < 0                                                                  # opposite sense for a prolate body


# ---------------------------------------------------------------- reaction wheel

def test_bang_bang_slew_matches_the_simulation_and_total_momentum_stays_zero():
    Iz, ang, T = I_SAT[2, 2], math.radians(30.0), 60.0
    tau, wmax, hmax = impl.bang_bang_slew(Iz, ang, T)
    assert tau * (T / 2) ** 2 / Iz == pytest.approx(ang)                               # two halves of (1/2) a t^2
    assert hmax == pytest.approx(Iz * wmax)
    wheel = lambda t, q, w: np.array([0.0, 0.0, -tau if t < T / 2 else tau])
    rot = impl.simulate_rigid_body(I_SAT, np.zeros(3), Q_ID, T, wheel_torque=wheel, n_samples=241, max_step=0.5)
    yaw = 2 * np.arctan2(rot.q[:, 3], rot.q[:, 0])
    assert yaw[-1] == pytest.approx(ang, rel=1e-6)
    assert np.max(rot.w[:, 2]) == pytest.approx(wmax, rel=1e-3)
    assert abs(rot.w[-1, 2]) < 1e-9
    body = rot.w @ I_SAT.T
    npt.assert_allclose(body + rot.h, 0.0, atol=1e-12)                                 # body and wheel cancel
    assert np.max(np.abs(rot.h[:, 2])) == pytest.approx(hmax, rel=1e-3)


# ---------------------------------------------------------------- the heavy top

def test_gravity_torque_is_horizontal_and_vanishes_upright():
    npt.assert_allclose(impl.gravity_torque_body(Q_ID, M_TOP, L_TOP), 0.0, atol=1e-15)
    th = 0.4
    q = impl.quat_from_axis_angle([1, 0, 0], th)
    tb = impl.gravity_torque_body(q, M_TOP, L_TOP)
    assert np.linalg.norm(tb) == pytest.approx(M_TOP * G * L_TOP * math.sin(th))
    assert tb[2] == pytest.approx(0.0, abs=1e-15)                                      # no torque about the symmetry axis
    ti = impl.quat_to_matrix(q) @ tb
    assert ti[2] == pytest.approx(0.0, abs=1e-15)                                      # nor about the vertical


def test_steady_precession_root_limits_and_validation():
    th, w3 = math.radians(30.0), 2 * math.pi * 50
    Om = impl.steady_precession_rate(I1_TOP, I3_TOP, M_TOP, L_TOP, th, w3)
    assert I1_TOP * math.cos(th) * Om ** 2 - I3_TOP * w3 * Om + M_TOP * G * L_TOP == pytest.approx(0.0, abs=1e-12)
    fast = impl.fast_top_precession(I3_TOP, M_TOP, L_TOP, w3)
    assert fast == pytest.approx(M_TOP * G * L_TOP / (I3_TOP * w3))
    assert Om > fast                                                                     # below 90 deg tilt the steady root lies above the fast-top value
    big = 2 * math.pi * 5000
    assert impl.steady_precession_rate(I1_TOP, I3_TOP, M_TOP, L_TOP, th, big) == pytest.approx(
        impl.fast_top_precession(I3_TOP, M_TOP, L_TOP, big), rel=1e-5)                 # fast-top limit
    assert impl.steady_precession_rate(I1_TOP, I3_TOP, M_TOP, L_TOP, math.pi / 2, w3) == pytest.approx(fast, rel=1e-12)
    with pytest.raises(ValueError):
        impl.steady_precession_rate(I1_TOP, I3_TOP, M_TOP, L_TOP, th, 2 * math.pi * 10)


def test_simulated_top_precesses_steadily_at_the_root_and_keeps_its_invariants():
    th, w3 = math.radians(30.0), 2 * math.pi * 50
    Om = impl.steady_precession_rate(I1_TOP, I3_TOP, M_TOP, L_TOP, th, w3)
    q0, wb = impl.top_initial_state(th, w3, Om)
    npt.assert_allclose(wb, [0.0, Om * math.sin(th), w3])
    rot = impl.simulate_rigid_body(I_TOP, wb, q0, 0.6, torque=_gravity, n_samples=301, rtol=1e-9)
    tilt, az = impl.axis_tilt_and_azimuth(rot)
    assert np.ptp(tilt) < 1e-6
    assert impl.measured_precession_rate(rot) == pytest.approx(Om, rel=1e-5)
    E, Lz, L3 = impl.top_invariants(rot, I_TOP, M_TOP, L_TOP)
    for x in (E, Lz, L3):
        assert np.ptp(x) < 1e-7 * abs(x[0])


def test_top_released_without_precession_nods_down_to_the_turning_point():
    th, w3 = math.radians(30.0), 2 * math.pi * 50
    q0, wb = impl.top_initial_state(th, w3, 0.0)
    rot = impl.simulate_rigid_body(I_TOP, wb, q0, 0.3, torque=_gravity, n_samples=3001, rtol=1e-9)
    tilt, _ = impl.axis_tilt_and_azimuth(rot)
    th_turn = impl.nutation_turning_tilt(I1_TOP, I3_TOP * w3, M_TOP, L_TOP, th)
    assert tilt.min() == pytest.approx(th, abs=1e-7)                                   # it starts at the top turning point
    assert tilt.max() == pytest.approx(th_turn, rel=1e-4)
    small = 2 * I1_TOP * M_TOP * G * L_TOP * math.sin(th) / (I3_TOP * w3) ** 2
    assert th_turn - th == pytest.approx(small, rel=0.15)                              # fast-top estimate is close
    assert th_turn - th > small                                                         # and slightly low here
    E, Lz, L3 = impl.top_invariants(rot, I_TOP, M_TOP, L_TOP)
    assert np.ptp(E) < 1e-7 * abs(E[0]) and np.ptp(Lz) < 1e-7 * abs(Lz[0])


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "spinning_top.png").stat().st_size > 10_000
