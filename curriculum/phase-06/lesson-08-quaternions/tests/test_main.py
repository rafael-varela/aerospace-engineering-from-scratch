"""Public tests for P06.08 Quaternions. They assert properties, not magic numbers: the Hamilton product's
rules (i j = k, non-commutative, associative, norm multiplicative), conjugate and inverse; unit
quaternions of axis-angle turns and back, with q and -q the same attitude; the quaternion DCM equal to
the P06.07 DCM route (Rodrigues, 3-2-1) and to an independent oracle (scipy.spatial.transform.Rotation,
used only here); Shepperd's DCM-to-quaternion conversion in every branch; composition order; Euler angles
both ways, including gimbal lock; the kinematic equation keeping |q| = 1 to first order and agreeing with
dC/dt = -[w x] C; integration against the DCM route and against closed forms, fourth-order RK4, the exact
exponential step, renormalisation; and SLERP's endpoints, constant rate and shortest arc.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import math

import numpy as np
import numpy.testing as npt
import pytest
from scipy.spatial.transform import Rotation

from conftest import load_impl

impl = load_impl(__file__)

DEG = math.pi / 180.0
I4 = np.array([1.0, 0.0, 0.0, 0.0])


def _rand_q(rng):
    q = rng.normal(size=4)
    return q / np.linalg.norm(q)


# ---------------------------------------------------------------- algebra

def test_hamilton_rules():
    i, j, k = np.eye(4)[1], np.eye(4)[2], np.eye(4)[3]
    npt.assert_allclose(impl.quat_multiply(i, j), k, atol=0)
    npt.assert_allclose(impl.quat_multiply(j, k), i, atol=0)
    npt.assert_allclose(impl.quat_multiply(k, i), j, atol=0)
    npt.assert_allclose(impl.quat_multiply(j, i), -k, atol=0)                  # not commutative
    for u in (i, j, k):
        npt.assert_allclose(impl.quat_multiply(u, u), -I4, atol=0)             # i^2 = j^2 = k^2 = -1
    npt.assert_allclose(impl.quat_multiply(impl.quat_multiply(i, j), k), -I4, atol=0)   # ijk = -1


def test_product_is_associative_and_norm_multiplicative():
    rng = np.random.default_rng(1)
    for _ in range(20):
        p, q, r = rng.normal(size=(3, 4))
        npt.assert_allclose(impl.quat_multiply(impl.quat_multiply(p, q), r), impl.quat_multiply(p, impl.quat_multiply(q, r)), atol=1e-13)
        assert impl.quat_norm(impl.quat_multiply(p, q)) == pytest.approx(impl.quat_norm(p) * impl.quat_norm(q), rel=1e-13)
        npt.assert_allclose(impl.quat_multiply(p, I4), p, atol=0)


def test_conjugate_and_inverse():
    rng = np.random.default_rng(2)
    p, q = rng.normal(size=(2, 4))
    npt.assert_allclose(impl.quat_conjugate(impl.quat_multiply(p, q)),
                        impl.quat_multiply(impl.quat_conjugate(q), impl.quat_conjugate(p)), atol=1e-14)
    npt.assert_allclose(impl.quat_multiply(q, impl.quat_inverse(q)), I4, atol=1e-14)
    npt.assert_allclose(impl.quat_multiply(impl.quat_inverse(q), q), I4, atol=1e-14)
    assert impl.quat_norm(impl.quat_normalize(q)) == pytest.approx(1.0, abs=1e-15)


def test_axis_angle_round_trip_and_double_cover():
    rng = np.random.default_rng(3)
    for _ in range(20):
        n = rng.normal(size=3)
        n /= np.linalg.norm(n)
        a = rng.uniform(0.01, math.pi - 0.01)
        q = impl.quat_from_axis_angle(n, a)
        assert impl.quat_norm(q) == pytest.approx(1.0, abs=1e-15)
        for qq in (q, -q):                                                      # q and -q: same attitude
            n2, a2 = impl.axis_angle_from_quat(qq)
            npt.assert_allclose(n2, n, atol=1e-12)
            assert a2 == pytest.approx(a, abs=1e-12)
    n0, a0 = impl.axis_angle_from_quat(I4)
    assert a0 == 0.0 and np.linalg.norm(n0) == pytest.approx(1.0)


# ---------------------------------------------------------------- quaternion <-> DCM <-> Euler

def test_quaternion_dcm_equals_rodrigues_and_scipy():
    rng = np.random.default_rng(4)
    for _ in range(20):
        n = rng.normal(size=3)
        a = rng.uniform(-math.pi, math.pi)
        q = impl.quat_from_axis_angle(n, a)
        C = impl.dcm_from_quat(q)
        npt.assert_allclose(C, impl.dcm_from_axis_angle(n, a), atol=1e-14)    # the P06.07 route
        npt.assert_allclose(C.T @ C, np.eye(3), atol=1e-14)
        assert np.linalg.det(C) == pytest.approx(1.0, abs=1e-14)
        # SciPy: scalar-last and active (body to reference), i.e. the transpose of C_b/n
        npt.assert_allclose(Rotation.from_quat([q[1], q[2], q[3], q[0]]).as_matrix().T, C, atol=1e-14)
        npt.assert_allclose(impl.dcm_from_quat(-q), C, atol=0)


def test_unnormalised_quaternion_scales_the_dcm():
    q = impl.quat_from_axis_angle([1.0, -2.0, 0.5], 0.8)
    C = impl.dcm_from_quat(1.01 * q)
    npt.assert_allclose(C, 1.01 ** 2 * impl.dcm_from_quat(q), atol=1e-14)


def test_quat_from_dcm_round_trip_in_every_branch():
    rng = np.random.default_rng(5)
    cases = [_rand_q(rng) for _ in range(30)]
    # near-180-degree turns about each axis exercise the q1, q2, q3 branches of Shepperd's method
    cases += [impl.quat_from_axis_angle(e, math.pi - 1e-3) for e in np.eye(3)]
    cases += [impl.quat_from_axis_angle(e, math.pi) for e in np.eye(3)]
    for q in cases:
        qb = impl.quat_from_dcm(impl.dcm_from_quat(q))
        assert qb[0] >= 0 and impl.quat_norm(qb) == pytest.approx(1.0, abs=1e-14)
        assert min(np.abs(qb - q).max(), np.abs(qb + q).max()) < 1e-13


def test_transform_matches_dcm_and_preserves_length():
    rng = np.random.default_rng(6)
    for _ in range(10):
        q, v = _rand_q(rng), rng.normal(size=3)
        vb = impl.transform(q, v)
        npt.assert_allclose(vb, impl.dcm_from_quat(q) @ v, atol=1e-14)
        assert np.linalg.norm(vb) == pytest.approx(np.linalg.norm(v), rel=1e-14)
    # a frame yawed 30 deg sees a vector 30 deg east of north straight ahead (P01.15's sign check)
    q30 = impl.quat_from_axis_angle([0, 0, 1], 30 * DEG)
    npt.assert_allclose(impl.transform(q30, np.array([math.cos(30 * DEG), math.sin(30 * DEG), 0.0])), [1, 0, 0], atol=1e-15)


def test_composition_order_mirrors_the_dcm():
    rng = np.random.default_rng(7)
    q_bn, q_cb = _rand_q(rng), _rand_q(rng)
    npt.assert_allclose(impl.dcm_from_quat(impl.quat_multiply(q_bn, q_cb)),
                        impl.dcm_from_quat(q_cb) @ impl.dcm_from_quat(q_bn), atol=1e-14)


def test_euler_321_both_ways_against_the_dcm_route():
    rng = np.random.default_rng(8)
    for _ in range(30):
        psi, th, ph = rng.uniform(-math.pi, math.pi), rng.uniform(-1.5, 1.5), rng.uniform(-math.pi, math.pi)
        q = impl.quat_from_euler_321(psi, th, ph)
        npt.assert_allclose(impl.dcm_from_quat(q), impl.dcm_from_euler_321(psi, th, ph), atol=1e-14)
        npt.assert_allclose(impl.euler_321_from_quat(q), [psi, th, ph], atol=1e-10)
    for th in (math.pi / 2, -math.pi / 2):                                     # gimbal lock
        q = impl.quat_from_euler_321(0.4, th, -0.3)
        e = impl.euler_321_from_quat(q)
        npt.assert_allclose(impl.dcm_from_euler_321(*e), impl.dcm_from_quat(q), atol=1e-12)


def test_error_angle_properties():
    rng = np.random.default_rng(9)
    qa, qb = _rand_q(rng), _rand_q(rng)
    e = impl.quat_error_angle(qa, qb)
    assert e == pytest.approx(impl.quat_error_angle(qb, qa), rel=1e-12)
    assert e == pytest.approx(impl.quat_error_angle(-qa, qb), rel=1e-12)       # sign does not matter
    assert e == pytest.approx(impl.misalignment_angle(impl.dcm_from_quat(qa), impl.dcm_from_quat(qb)), abs=1e-12)
    assert 0.0 <= e <= math.pi
    assert impl.quat_error_angle(qa, qa) == pytest.approx(0.0, abs=1e-7)


# ---------------------------------------------------------------- kinematics and integration

def test_quat_rate_is_orthogonal_to_q_and_matches_dcm_kinematics():
    rng = np.random.default_rng(10)
    h = 1e-6
    for _ in range(10):
        q, w = _rand_q(rng), rng.normal(size=3)
        qd = impl.quat_rate(q, w)
        assert q @ qd == pytest.approx(0.0, abs=1e-15)                         # d|q|^2/dt = 0
        Cd_fd = (impl.dcm_from_quat(q + h * qd) - impl.dcm_from_quat(q - h * qd)) / (2 * h)
        npt.assert_allclose(Cd_fd, -impl.skew(w) @ impl.dcm_from_quat(q), atol=1e-8)


def test_quaternion_route_matches_dcm_route_for_coning():
    pd, sd, beta = 0.4, 2.0, 20 * DEG
    w = impl.coning_rate(pd, sd, beta)
    C0 = impl.coning_dcm(0.0, pd, sd, beta)
    tq, Q = impl.integrate_quat(impl.quat_from_dcm(C0), w, 20.0, n_samples=81)
    tc, Cs = impl.integrate_dcm(C0, w, 20.0, n_samples=81)
    for m in range(tq.size):
        Cq = impl.dcm_from_quat(impl.quat_normalize(Q[m]))
        assert impl.misalignment_angle(Cq, Cs[m]) < 1e-9                       # P06.07 route
        assert impl.misalignment_angle(Cq, impl.coning_dcm(tq[m], pd, sd, beta)) < 1e-9   # closed form
    assert np.abs(np.linalg.norm(Q, axis=1) - 1).max() < 1e-9


def test_exact_step_reproduces_fixed_axis_rotation():
    w = np.array([1.0, 2.0, 2.0]) / 3.0 * (30 * DEG)
    q = I4.copy()
    for _ in range(500):
        q = impl.quat_exp_step(q, w, 0.04)
    exact = impl.quat_from_axis_angle(w, np.linalg.norm(w) * 20.0)
    assert impl.quat_error_angle(q, exact) < 1e-13
    assert abs(impl.quat_norm(q) - 1) < 1e-13
    npt.assert_allclose(impl.quat_exp_step(I4, np.zeros(3), 0.1), I4, atol=0)


def test_rk4_order_and_norm_drift():
    pd, sd, beta = 0.4, 2.0, 20 * DEG
    w = impl.coning_rate(pd, sd, beta)
    q0 = impl.quat_from_dcm(impl.coning_dcm(0.0, pd, sd, beta))
    errs, drift = [], []
    for dt in [0.1, 0.05, 0.025]:
        t, Q = impl.integrate_quat_rk4(q0, w, 10.0, dt)
        errs.append(impl.quat_error_angle(Q[-1], impl.quat_from_dcm(impl.coning_dcm(t[-1], pd, sd, beta))))
        drift.append(abs(impl.quat_norm(Q[-1]) - 1))
    o_e = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    o_d = np.log2(np.array(drift[:-1]) / np.array(drift[1:]))
    assert np.all((o_e > 3.7) & (o_e < 4.3)), o_e
    assert np.all((o_d > 4.5) & (o_d < 5.5)), o_d                              # norm defect: one order higher
    t, Qn = impl.integrate_quat_rk4(q0, w, 10.0, 0.1, normalize=True)
    assert np.abs(np.linalg.norm(Qn, axis=1) - 1).max() < 1e-14
    e_n = impl.quat_error_angle(Qn[-1], impl.quat_from_dcm(impl.coning_dcm(t[-1], pd, sd, beta)))
    assert e_n == pytest.approx(errs[0], rel=0.05)                             # renormalising does not fix attitude


# ---------------------------------------------------------------- interpolation

def test_slerp_endpoints_unit_norm_and_constant_rate():
    q0 = impl.quat_from_euler_321(0.2, 0.1, -0.3)
    q1 = impl.quat_from_euler_321(1.4, -0.5, 0.9)
    npt.assert_allclose(impl.slerp(q0, q1, 0.0), q0, atol=1e-14)
    assert min(np.abs(impl.slerp(q0, q1, 1.0) - q1).max(), np.abs(impl.slerp(q0, q1, 1.0) + q1).max()) < 1e-14
    us = np.linspace(0, 1, 11)
    qs = [impl.slerp(q0, q1, u) for u in us]
    for q in qs:
        assert impl.quat_norm(q) == pytest.approx(1.0, abs=1e-14)
    steps = [impl.quat_error_angle(qs[m + 1], qs[m]) for m in range(10)]
    npt.assert_allclose(steps, steps[0], rtol=1e-10)                           # equal angle per equal step
    assert sum(steps) == pytest.approx(impl.quat_error_angle(q1, q0), rel=1e-10)   # a single great arc
    # nlerp follows the same arc but not at constant rate
    qn = impl.nlerp(q0, q1, 0.25)
    assert impl.quat_error_angle(qn, q0) != pytest.approx(0.25 * impl.quat_error_angle(q1, q0), rel=1e-3)
    assert impl.quat_error_angle(impl.nlerp(q0, q1, 0.5), impl.slerp(q0, q1, 0.5)) < 1e-12


def test_slerp_takes_the_shorter_arc():
    qa = I4
    qb = impl.quat_from_axis_angle([0, 0, 1], 170 * DEG)
    for qq in (qb, -qb):
        _, a = impl.axis_angle_from_quat(impl.slerp(qa, qq, 0.5))
        assert a == pytest.approx(85 * DEG, rel=1e-12)
    _, a_long = impl.axis_angle_from_quat(impl.nlerp(qa, -qb, 0.5, shortest=False))
    assert a_long == pytest.approx(95 * DEG, rel=1e-12)                        # the long way round
    npt.assert_allclose(impl.slerp(qa, qa, 0.3), qa, atol=1e-14)


def test_star_tracker_slerp_interpolation_converges_at_second_order():
    pd, sd, beta = 0.4, 2.0, 20 * DEG
    t_fine = np.linspace(0.0, 6.0, 241)
    truth = [impl.quat_from_dcm(impl.coning_dcm(t, pd, sd, beta)) for t in t_fine]
    errs = []
    for dts in [0.4, 0.2, 0.1]:
        ts = np.arange(0.0, 6.0 + 1e-9, dts)
        qs = np.array([impl.quat_from_dcm(impl.coning_dcm(t, pd, sd, beta)) for t in ts])
        qi = impl.interpolate_samples(ts, qs, t_fine, "slerp")
        errs.append(max(impl.quat_error_angle(a, b) for a, b in zip(qi, truth)))
    orders = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert np.all((orders > 1.8) & (orders < 2.3)), orders


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "quaternions.png").stat().st_size > 10_000
