"""Public tests for P06.07 3D Rotations: Euler Angles & DCMs. They assert properties, not magic numbers:
every DCM is a rotation (C^T C = E, det C = +1); the twelve Euler sequences agree with an independent
oracle (scipy.spatial.transform.Rotation, used only here) and with the written-out 3-2-1 matrix; angles
round-trip through the DCM, and at gimbal lock the rebuilt matrix is still right; the Euler-rate matrix
matches the textbook 3-2-1 kinematics, has |det B| = |cos a2| (|sin a2| for symmetric sequences) and agrees
with a finite difference of the DCM; the DCM kinematic equation integrated by solve_ivp and by RK4
reproduces two closed-form attitude histories (a fixed-axis rotation and steady coning), RK4 converges at
fourth order, SVD re-orthonormalisation restores C^T C = E without changing the attitude, and Euler-angle
integration loses accuracy near gimbal lock where the DCM does not.

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
SEQS = ["121", "123", "131", "132", "212", "213", "231", "232", "312", "313", "321", "323"]
AXIS_LETTER = {"1": "X", "2": "Y", "3": "Z"}


def _random_angles(rng, seq):
    a1, a3 = rng.uniform(-math.pi, math.pi, 2)
    a2 = rng.uniform(0.2, math.pi - 0.2) if seq[0] == seq[2] else rng.uniform(-1.3, 1.3)
    return np.array([a1, a2, a3])


# ---------------------------------------------------------------- DCMs

def test_every_sequence_gives_a_rotation_matrix():
    rng = np.random.default_rng(7)
    for seq in SEQS:
        for _ in range(5):
            C = impl.dcm_from_euler(_random_angles(rng, seq), seq)
            npt.assert_allclose(C.T @ C, np.eye(3), atol=1e-14)
            assert np.linalg.det(C) == pytest.approx(1.0, abs=1e-14)
            assert impl.orthogonality_error(C) < 1e-14


def test_321_matches_the_written_out_matrix():
    psi, th, ph = 30 * DEG, 10 * DEG, -5 * DEG
    c, s = math.cos, math.sin
    C_ref = np.array([
        [c(th) * c(psi), c(th) * s(psi), -s(th)],
        [s(ph) * s(th) * c(psi) - c(ph) * s(psi), s(ph) * s(th) * s(psi) + c(ph) * c(psi), s(ph) * c(th)],
        [c(ph) * s(th) * c(psi) + s(ph) * s(psi), c(ph) * s(th) * s(psi) - s(ph) * c(psi), c(ph) * c(th)]])
    npt.assert_allclose(impl.dcm_from_euler([psi, th, ph], "321"), C_ref, atol=1e-15)


def test_all_twelve_sequences_against_scipy_oracle():
    rng = np.random.default_rng(12)
    for seq in SEQS:
        a = _random_angles(rng, seq)
        letters = "".join(AXIS_LETTER[ch] for ch in seq)            # upper case: intrinsic (body-axis) turns
        C_scipy = Rotation.from_euler(letters, a).as_matrix().T      # SciPy is active; the lesson is passive
        npt.assert_allclose(impl.dcm_from_euler(a, seq), C_scipy, atol=1e-14)


def test_bad_sequences_are_refused():
    for bad in ["112", "331", "12", "1234", "124", "abc"]:
        with pytest.raises(ValueError):
            impl.dcm_from_euler([0.1, 0.2, 0.3], bad)
    assert impl.is_symmetric("313") and not impl.is_symmetric("321")


def test_euler_angles_round_trip_for_all_sequences():
    rng = np.random.default_rng(3)
    for seq in SEQS:
        for _ in range(10):
            a = _random_angles(rng, seq)
            npt.assert_allclose(impl.euler_from_dcm(impl.dcm_from_euler(a, seq), seq), a, atol=1e-11)


def test_extraction_of_an_arbitrary_attitude_rebuilds_it():
    rng = np.random.default_rng(4)
    for _ in range(20):
        C = Rotation.random(random_state=rng).as_matrix().T
        for seq in SEQS:
            a = impl.euler_from_dcm(C, seq)
            npt.assert_allclose(impl.dcm_from_euler(a, seq), C, atol=1e-12)
            if seq[0] == seq[2]:
                assert 0.0 <= a[1] <= math.pi
            else:
                assert -math.pi / 2 <= a[1] <= math.pi / 2


def test_gimbal_lock_rebuilds_the_matrix():
    for seq, a2 in [("321", math.pi / 2), ("321", -math.pi / 2), ("123", math.pi / 2), ("313", 0.0), ("313", math.pi), ("212", 0.0)]:
        a = np.array([0.7, a2, -0.4])
        C = impl.dcm_from_euler(a, seq)
        b = impl.euler_from_dcm(C, seq)
        npt.assert_allclose(impl.dcm_from_euler(b, seq), C, atol=1e-12)      # angles differ, attitude does not
        assert b[2] == 0.0                                                    # the convention: a3 set to zero


# ---------------------------------------------------------------- Euler-angle kinematics

def test_rate_matrix_matches_textbook_321_kinematics():
    psi, th, ph = 0.4, -0.3, 1.1
    dpsi, dth, dph = 0.05, -0.2, 0.3
    p = dph - dpsi * math.sin(th)
    q = dth * math.cos(ph) + dpsi * math.sin(ph) * math.cos(th)
    r = -dth * math.sin(ph) + dpsi * math.cos(ph) * math.cos(th)
    npt.assert_allclose(impl.euler_rate_matrix([psi, th, ph], "321") @ [dpsi, dth, dph], [p, q, r], atol=1e-15)
    npt.assert_allclose(impl.euler_angle_rates([psi, th, ph], np.array([p, q, r]), "321"), [dpsi, dth, dph], atol=1e-14)


def test_rate_matrix_determinant_shows_the_singularity():
    rng = np.random.default_rng(9)
    for seq in SEQS:
        a = _random_angles(rng, seq)
        d = abs(np.linalg.det(impl.euler_rate_matrix(a, seq)))
        assert d == pytest.approx(abs(math.sin(a[1])) if seq[0] == seq[2] else abs(math.cos(a[1])), rel=1e-12)
    with pytest.raises(ValueError):
        impl.euler_angle_rates([0.3, math.pi / 2, 0.2], np.array([0.0, 0.0, 0.1]), "321")
    with pytest.raises(ValueError):
        impl.euler_angle_rates([0.3, 0.0, 0.2], np.array([0.1, 0.0, 0.0]), "313")


def test_rate_matrix_agrees_with_finite_difference_of_the_dcm():
    """d/dt C(a(t)) by central differences must equal -[w x] C with w = B(a) da/dt (the kinematic equation)."""
    rng = np.random.default_rng(10)
    h = 1e-6
    for seq in SEQS:
        a, ad = _random_angles(rng, seq), rng.normal(size=3)
        Cdot_fd = (impl.dcm_from_euler(a + h * ad, seq) - impl.dcm_from_euler(a - h * ad, seq)) / (2 * h)
        w = impl.euler_rate_matrix(a, seq) @ ad
        npt.assert_allclose(impl.dcm_rate(impl.dcm_from_euler(a, seq), w), Cdot_fd, atol=1e-8)


def test_dcm_rate_preserves_orthogonality_to_first_order():
    C = impl.dcm_from_euler([0.3, -0.5, 1.2], "321")
    Cd = impl.dcm_rate(C, np.array([0.2, -0.7, 0.4]))
    npt.assert_allclose(Cd.T @ C + C.T @ Cd, np.zeros((3, 3)), atol=1e-15)   # d/dt (C^T C) = 0


# ---------------------------------------------------------------- integration against closed forms

def test_solve_ivp_dcm_matches_single_axis_closed_form():
    axis = np.array([1.0, 2.0, 2.0]) / 3.0
    rate = 30 * DEG
    C0 = impl.dcm_from_euler([0.2, 0.1, -0.3], "321")
    t, Cs = impl.integrate_dcm(C0, lambda t: rate * axis, 20.0, n_samples=101)
    for m in range(0, t.size, 10):
        assert impl.misalignment_angle(Cs[m], impl.single_axis_dcm(C0, axis, rate, t[m])) < 1e-9
        assert impl.orthogonality_error(Cs[m]) < 1e-8
    # the rotation axis has fixed body components, so its reference-frame components C^T axis never change
    npt.assert_allclose(Cs[-1].T @ axis, C0.T @ axis, atol=1e-9)


def test_solve_ivp_dcm_matches_steady_coning():
    pd, sd, beta = 0.4, 2.0, 20 * DEG
    C0 = impl.coning_dcm(0.0, pd, sd, beta)
    t, Cs = impl.integrate_dcm(C0, lambda t: impl.coning_body_rate(t, pd, sd, beta), 20.0, n_samples=81)
    for m in range(t.size):
        assert impl.misalignment_angle(Cs[m], impl.coning_dcm(t[m], pd, sd, beta)) < 1e-8
    # the inertial 3 axis seen from the body stays at the half-cone angle beta from body axis 3
    tilt = [math.acos(np.clip(C[2, 2], -1, 1)) for C in Cs]
    npt.assert_allclose(tilt, beta, atol=1e-8)


def test_coning_body_rate_is_consistent_with_the_rate_matrix():
    pd, sd, beta, t = 0.4, 2.0, 0.35, 1.7
    a = np.array([pd * t, beta, sd * t])
    npt.assert_allclose(impl.coning_body_rate(t, pd, sd, beta), impl.euler_rate_matrix(a, "313") @ [pd, 0.0, sd], atol=1e-14)


def test_rk4_converges_at_fourth_order():
    pd, sd, beta = 0.4, 2.0, 20 * DEG
    C0 = impl.coning_dcm(0.0, pd, sd, beta)
    w = lambda t: impl.coning_body_rate(t, pd, sd, beta)
    errs = []
    for dt in [0.1, 0.05, 0.025]:
        t, Cs = impl.integrate_dcm_rk4(C0, w, 10.0, dt)
        errs.append(impl.misalignment_angle(Cs[-1], impl.coning_dcm(t[-1], pd, sd, beta)))
    orders = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert np.all((orders > 3.7) & (orders < 4.3)), orders


def test_svd_reorthonormalisation_fixes_drift_without_moving_the_attitude():
    pd, sd, beta = 0.4, 2.0, 20 * DEG
    C0 = impl.coning_dcm(0.0, pd, sd, beta)
    w = lambda t: impl.coning_body_rate(t, pd, sd, beta)
    t, raw = impl.integrate_dcm_rk4(C0, w, 30.0, 0.1)
    _, ren = impl.integrate_dcm_rk4(C0, w, 30.0, 0.1, renormalize_every=1)
    exact = impl.coning_dcm(t[-1], pd, sd, beta)
    assert impl.orthogonality_error(raw[-1]) > 1e-4                         # RK4 alone drifts
    assert impl.orthogonality_error(ren[-1]) < 1e-13                        # renormalised every step
    e_raw = impl.misalignment_angle(impl.orthonormalize_svd(raw[-1]), exact)
    e_ren = impl.misalignment_angle(ren[-1], exact)
    assert e_ren == pytest.approx(e_raw, rel=0.05)                          # the attitude error is not cured


def test_svd_gives_the_nearest_rotation():
    rng = np.random.default_rng(5)
    R0 = Rotation.random(random_state=rng).as_matrix()
    A = R0 + 1e-3 * rng.normal(size=(3, 3))
    R = impl.orthonormalize_svd(A)
    npt.assert_allclose(R.T @ R, np.eye(3), atol=1e-14)
    assert np.linalg.det(R) == pytest.approx(1.0, abs=1e-14)
    npt.assert_allclose(impl.orthonormalize_svd(R0), R0, atol=1e-14)        # a rotation is left alone
    d = np.linalg.norm(A - R)
    for _ in range(50):                                                      # no nearby rotation is closer
        Q = Rotation.from_rotvec(1e-3 * rng.normal(size=3)).as_matrix() @ R
        assert np.linalg.norm(A - Q) >= d - 1e-15


def test_misalignment_angle_properties():
    C1 = impl.dcm_from_euler([0.3, 0.2, -0.1], "321")
    assert impl.misalignment_angle(C1, C1) == pytest.approx(0.0, abs=1e-15)
    C2 = impl.dcm_from_axis_angle(np.array([0.0, 0.6, 0.8]), 0.4) @ C1
    assert impl.misalignment_angle(C2, C1) == pytest.approx(0.4, rel=1e-12)
    assert impl.misalignment_angle(C1, C2) == pytest.approx(0.4, rel=1e-12)
    C3 = impl.dcm_from_axis_angle(np.array([1.0, 0.0, 0.0]), math.pi) @ C1
    assert impl.misalignment_angle(C3, C1) == pytest.approx(math.pi, rel=1e-12)


def test_euler_angle_integration_agrees_away_from_gimbal_lock():
    pd, sd, beta = 0.4, 2.0, 20 * DEG
    w = lambda t: impl.coning_body_rate(t, pd, sd, beta)
    res = impl.integrate_euler_angles([0.0, beta, 0.0], w, 10.0, seq="313", n_samples=51)
    assert res["success"]
    for m in range(res["t"].size):
        C = impl.dcm_from_euler(res["angles"][m], "313")
        assert impl.misalignment_angle(C, impl.coning_dcm(res["t"][m], pd, sd, beta)) < 1e-8


def test_hammerhead_euler_angles_fail_near_vertical_while_dcm_does_not():
    errs_e, errs_c = [], []
    for delta in [10 * DEG, 0.1 * DEG]:
        segs, C_end = impl.hammerhead(0.5, 0.5, math.pi / 2 - delta)
        C = impl.integrate_segments_rk4(np.eye(3), segs, 0.01, "dcm")
        a = impl.integrate_segments_rk4(np.zeros(3), segs, 0.01, "euler")
        errs_c.append(impl.misalignment_angle(C, C_end))
        errs_e.append(impl.misalignment_angle(impl.dcm_from_euler(a, "321"), C_end))
    assert max(errs_c) < 1e-9                                  # the DCM does not care how close to vertical
    assert errs_e[0] < 1e-8                                    # 10 deg from vertical: angles are fine
    assert errs_e[1] > 1e3 * errs_e[0]                         # 0.1 deg: the angle rates blow up
    # the closed form: heading reversed, pitch back to theta_top below the horizon... seen as yaw 180 deg
    psi, th, ph = impl.euler_from_dcm(C_end, "321")
    assert abs(abs(psi) - math.pi) < 1e-9 and abs(ph) < 1e-9


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "attitude_kinematics.png").stat().st_size > 10_000
