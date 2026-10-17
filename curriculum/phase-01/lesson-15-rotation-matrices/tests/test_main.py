"""Public tests for P01.15 Matrix Transformations & Rotation Matrices. They assert properties, not magic numbers.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-01/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.spatial.transform import Rotation

from conftest import load_impl

impl = load_impl(__file__)

DEG = np.pi / 180.0
RNG = np.random.default_rng(1515)
ATTITUDES = [(30 * DEG, 10 * DEG, -5 * DEG), (-120 * DEG, -40 * DEG, 170 * DEG), (0.0, 0.0, 0.0),
             (179 * DEG, 89 * DEG, -179 * DEG)]


# ---------------------------------------------------------------- elementary rotations

def test_elementary_rotations_are_proper_orthogonal_for_any_angle():
    for a in np.linspace(-2 * np.pi, 2 * np.pi, 17):
        for R in (impl.rot_x(a), impl.rot_y(a), impl.rot_z(a)):
            npt.assert_allclose(R.T @ R, np.eye(3), atol=1e-14)
            assert np.linalg.det(R) == pytest.approx(1.0, abs=1e-14)
            assert impl.is_rotation_matrix(R)


def test_elementary_rotations_reduce_to_identity_and_invert_by_negative_angle():
    for R in (impl.rot_x, impl.rot_y, impl.rot_z):
        npt.assert_allclose(R(0.0), np.eye(3), atol=1e-15)
        npt.assert_allclose(R(0.7) @ R(-0.7), np.eye(3), atol=1e-15)
        npt.assert_allclose(R(0.3) @ R(0.4), R(0.7), atol=1e-15)          # angles add about one axis


def test_rot_z_is_passive_a_frame_turned_toward_a_vector_sees_it_on_its_x_axis():
    """A vector 30 deg east of north has components (cos30, sin30, 0) in NED. A frame yawed by 30 deg
    looks straight at it, so in that frame the vector lies along +x. That is the passive convention."""
    v_ned = np.array([np.cos(30 * DEG), np.sin(30 * DEG), 0.0])
    npt.assert_allclose(impl.rot_z(30 * DEG) @ v_ned, [1.0, 0.0, 0.0], atol=1e-15)
    # Right-hand rule about x: a frame rolled by +90 deg has its new y where old z was and its new z where
    # old -y was, so the old +y axis has components (0, 0, -1) in the rolled frame.
    npt.assert_allclose(impl.rot_x(90 * DEG) @ np.array([0.0, 1.0, 0.0]), [0.0, 0.0, -1.0], atol=1e-15)


# ---------------------------------------------------------------- the 3-2-1 DCM

@pytest.mark.parametrize("psi,theta,phi", ATTITUDES)
def test_dcm_matches_scipy_rotation_oracle(psi, theta, phi):
    C = impl.dcm_from_euler_321(psi, theta, phi)
    # scipy builds the *active* body-to-NED matrix for the intrinsic Z-Y-X sequence; our C_bn is its transpose.
    R_nb = Rotation.from_euler("ZYX", [psi, theta, phi]).as_matrix()
    npt.assert_allclose(C, R_nb.T, atol=1e-14)
    assert impl.is_rotation_matrix(C)
    npt.assert_allclose(np.linalg.inv(C), C.T, atol=1e-13)                 # inverse is the transpose


def test_dcm_closed_form_entries():
    psi, theta, phi = 0.4, -0.3, 1.1
    C = impl.dcm_from_euler_321(psi, theta, phi)
    assert C[0, 2] == pytest.approx(-np.sin(theta))
    assert C[0, 0] == pytest.approx(np.cos(theta) * np.cos(psi))
    assert C[0, 1] == pytest.approx(np.cos(theta) * np.sin(psi))
    assert C[1, 2] == pytest.approx(np.sin(phi) * np.cos(theta))
    assert C[2, 2] == pytest.approx(np.cos(phi) * np.cos(theta))


def test_yaw_only_pitch_only_roll_only_reduce_to_elementary_rotations():
    npt.assert_allclose(impl.dcm_from_euler_321(0.5, 0.0, 0.0), impl.rot_z(0.5), atol=1e-15)
    npt.assert_allclose(impl.dcm_from_euler_321(0.0, 0.5, 0.0), impl.rot_y(0.5), atol=1e-15)
    npt.assert_allclose(impl.dcm_from_euler_321(0.0, 0.0, 0.5), impl.rot_x(0.5), atol=1e-15)


def test_is_rotation_matrix_rejects_reflections_scalings_and_wrong_shapes():
    reflection = np.diag([1.0, 1.0, -1.0])                 # orthogonal, det = -1
    assert np.allclose(reflection.T @ reflection, np.eye(3)) and not impl.is_rotation_matrix(reflection)
    assert not impl.is_rotation_matrix(1.01 * np.eye(3))
    assert not impl.is_rotation_matrix(np.eye(2))
    assert not impl.is_rotation_matrix(np.eye(3) + 1e-6, tol=1e-9)
    assert impl.is_rotation_matrix(np.eye(3) + 1e-12, tol=1e-9)


def test_transform_preserves_length_and_angles_and_handles_batches():
    C = impl.dcm_from_euler_321(*ATTITUDES[1])
    v = RNG.normal(size=(50, 3))
    w = impl.transform(C, v)
    assert w.shape == v.shape
    npt.assert_allclose(np.linalg.norm(w, axis=1), np.linalg.norm(v, axis=1), rtol=1e-13)
    npt.assert_allclose(w @ w.T, v @ v.T, atol=1e-11)                      # all dot products preserved
    npt.assert_allclose(impl.transform(C, v[3]), w[3], atol=1e-14)
    npt.assert_allclose(impl.transform(C.T, w), v, atol=1e-13)              # transpose undoes it
    with pytest.raises(ValueError):
        impl.transform(C, np.ones((3, 4)))


def test_compose_order_and_inverse():
    C_ba = impl.dcm_from_euler_321(0.2, 0.1, -0.3)
    C_cb = impl.dcm_from_euler_321(-1.0, 0.4, 0.9)
    C_ca = impl.compose(C_cb, C_ba)
    v_a = np.array([1.0, -2.0, 0.5])
    npt.assert_allclose(C_ca @ v_a, C_cb @ (C_ba @ v_a), atol=1e-14)
    assert impl.is_rotation_matrix(C_ca)
    assert not np.allclose(C_ca, impl.compose(C_ba, C_cb))                 # rotations do not commute
    npt.assert_allclose(impl.compose(C_ca.T, C_ca), np.eye(3), atol=1e-14)


def test_dcm_is_not_the_sum_of_its_angles_but_a_product():
    """A 3-2-1 DCM equals the elementary rotations multiplied right-to-left; left-to-right is a different attitude."""
    psi, theta, phi = ATTITUDES[0]
    npt.assert_allclose(impl.dcm_from_euler_321(psi, theta, phi), impl.rot_x(phi) @ impl.rot_y(theta) @ impl.rot_z(psi))
    assert not np.allclose(impl.dcm_from_euler_321(psi, theta, phi), impl.rot_z(psi) @ impl.rot_y(theta) @ impl.rot_x(phi))


# ---------------------------------------------------------------- Euler-angle extraction

@pytest.mark.parametrize("psi,theta,phi", ATTITUDES)
def test_euler_extraction_round_trips(psi, theta, phi):
    out = impl.euler_321_from_dcm(impl.dcm_from_euler_321(psi, theta, phi))
    npt.assert_allclose(out, [psi, theta, phi], atol=1e-12)


def test_euler_extraction_random_round_trip_through_scipy():
    for _ in range(20):
        q = Rotation.random(random_state=int(RNG.integers(1 << 31)))
        C = q.as_matrix().T
        psi, theta, phi = impl.euler_321_from_dcm(C)
        assert -np.pi / 2 <= theta <= np.pi / 2
        npt.assert_allclose(impl.dcm_from_euler_321(psi, theta, phi), C, atol=1e-12)


@pytest.mark.parametrize("theta_sign", [1.0, -1.0])
def test_euler_extraction_at_gimbal_lock_reproduces_the_dcm(theta_sign):
    psi0, phi0 = 0.8, -0.5
    C = impl.dcm_from_euler_321(psi0, theta_sign * np.pi / 2, phi0)
    psi, theta, phi = impl.euler_321_from_dcm(C)
    assert phi == 0.0 and theta == pytest.approx(theta_sign * np.pi / 2)
    npt.assert_allclose(impl.dcm_from_euler_321(psi, theta, phi), C, atol=1e-12)   # same attitude, different angles
    assert psi == pytest.approx(psi0 - theta_sign * phi0)                          # only psi -/+ phi survives


# ---------------------------------------------------------------- axis-angle and Rodrigues

def test_skew_matrix_implements_the_cross_product():
    for _ in range(5):
        a, b = RNG.normal(size=3), RNG.normal(size=3)
        npt.assert_allclose(impl.skew(a) @ b, np.cross(a, b), atol=1e-14)
        assert np.allclose(impl.skew(a).T, -impl.skew(a))
        assert impl.skew(a) @ a == pytest.approx(0.0, abs=1e-14)


def test_rodrigues_reproduces_elementary_rotations_and_euler_dcms():
    npt.assert_allclose(impl.dcm_from_axis_angle([0, 0, 1], 0.6), impl.rot_z(0.6), atol=1e-15)
    npt.assert_allclose(impl.dcm_from_axis_angle([0, 5, 0], 0.6), impl.rot_y(0.6), atol=1e-15)   # axis normalised
    npt.assert_allclose(impl.dcm_from_axis_angle([1, 0, 0], 0.6), impl.rot_x(0.6), atol=1e-15)
    for att in ATTITUDES:
        C = impl.dcm_from_euler_321(*att)
        n, a = impl.axis_angle_from_dcm(C)
        assert 0.0 <= a <= np.pi and np.linalg.norm(n) == pytest.approx(1.0)
        npt.assert_allclose(impl.dcm_from_axis_angle(n, a), C, atol=1e-12)
        npt.assert_allclose(C @ n, n, atol=1e-12)                                   # the axis is unchanged
        assert np.cos(a) == pytest.approx((np.trace(C) - 1) / 2)
    with pytest.raises(ValueError):
        impl.dcm_from_axis_angle([0.0, 0.0, 0.0], 1.0)


def test_axis_angle_handles_identity_and_half_turn():
    n, a = impl.axis_angle_from_dcm(np.eye(3))
    assert a == 0.0 and np.linalg.norm(n) == pytest.approx(1.0)
    for axis in ([0, 0, 1.0], [1, 1, 0.0], [1, -2, 3.0]):
        C = impl.dcm_from_axis_angle(axis, np.pi)
        n, a = impl.axis_angle_from_dcm(C)
        assert a == pytest.approx(np.pi, abs=1e-7)
        unit = np.asarray(axis) / np.linalg.norm(axis)
        assert abs(n @ unit) == pytest.approx(1.0, abs=1e-6)                        # sign is arbitrary at pi
        npt.assert_allclose(impl.dcm_from_axis_angle(n, a), C, atol=1e-6)


def test_small_angle_dcm_is_first_order_accurate_and_second_order_non_orthogonal():
    axis = np.array([0.6, 0.0, 0.8])
    errs, orths = [], []
    for d in (1e-1, 5e-2, 2.5e-2, 1.25e-2):
        approx = impl.small_angle_dcm(axis * d)
        errs.append(np.max(np.abs(approx - impl.dcm_from_axis_angle(axis, d))))
        orths.append(impl.orthogonality_error(approx))
    rates = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert np.all(rates > 1.9), f"small-angle error should fall as dtheta^2, got rates {rates}"
    npt.assert_allclose(orths, [d * d for d in (1e-1, 5e-2, 2.5e-2, 1.25e-2)], rtol=1e-12)


def test_orthonormalize_restores_a_drifted_dcm_quadratically():
    C_true = impl.dcm_from_euler_321(*ATTITUDES[1])
    C = C_true + 1e-3 * RNG.normal(size=(3, 3))
    e0 = impl.orthogonality_error(C)
    e1 = impl.orthogonality_error(impl.orthonormalize(C))
    e2 = impl.orthogonality_error(impl.orthonormalize(C, iterations=2))
    assert e1 < 1e-2 * e0 and e2 < 1e-2 * e1
    assert impl.is_rotation_matrix(impl.orthonormalize(C, iterations=4), tol=1e-12)
    npt.assert_allclose(impl.orthonormalize(C_true), C_true, atol=1e-15)           # a rotation is a fixed point


# ---------------------------------------------------------------- applications

def test_gravity_in_body_closed_form_and_magnitude():
    g = 9.80665
    for theta, phi in [(0.0, 0.0), (10 * DEG, -5 * DEG), (-30 * DEG, 60 * DEG), (89 * DEG, 0.0)]:
        gb = impl.gravity_in_body(theta, phi, g)
        npt.assert_allclose(gb, g * np.array([-np.sin(theta), np.sin(phi) * np.cos(theta), np.cos(phi) * np.cos(theta)]), atol=1e-13)
        assert np.linalg.norm(gb) == pytest.approx(g)
    npt.assert_allclose(impl.gravity_in_body(0.0, 0.0), [0.0, 0.0, 9.80665])
    # Nose up: gravity has a component along -x (pointing backwards along the fuselage).
    assert impl.gravity_in_body(20 * DEG, 0.0)[0] < 0


def test_body_velocity_to_ned_and_flight_path_angle():
    V, alpha = 70.0, 4 * DEG
    v_b = V * np.array([np.cos(alpha), 0.0, np.sin(alpha)])
    # Wings level: gamma = theta - alpha exactly.
    v_n = impl.body_velocity_to_ned(30 * DEG, 10 * DEG, 0.0, v_b)
    assert np.linalg.norm(v_n) == pytest.approx(V)
    assert impl.flight_path_angle(v_n) == pytest.approx(10 * DEG - alpha)
    assert np.arctan2(v_n[1], v_n[0]) == pytest.approx(30 * DEG)                  # track = heading with no sideslip
    # With bank at fixed theta and alpha, sin(gamma) = cos(a) sin(th) - sin(a) cos(phi) cos(th): part of the
    # angle of attack now tilts sideways instead of pushing the path below the nose, so gamma rises.
    v_n_bank = impl.body_velocity_to_ned(30 * DEG, 10 * DEG, 40 * DEG, v_b)
    gamma_bank = impl.flight_path_angle(v_n_bank)
    assert gamma_bank > impl.flight_path_angle(v_n)
    assert np.sin(gamma_bank) == pytest.approx(np.cos(alpha) * np.sin(10 * DEG) - np.sin(alpha) * np.cos(40 * DEG) * np.cos(10 * DEG))
    with pytest.raises(ValueError):
        impl.flight_path_angle(np.zeros(3))


def test_integration_converges_at_second_order_and_renormalization_removes_orthogonality_drift():
    axis, rate, T = np.array([1.0, 2.0, 2.0]) / 3.0, 30 * DEG, 20.0
    raw = np.array([impl.integration_errors(axis, rate, T, n) for n in (200, 400, 800)])
    ren = np.array([impl.integration_errors(axis, rate, T, n, renormalize=True) for n in (200, 400, 800)])
    rates = np.log2(raw[:-1, 1] / raw[1:, 1])
    assert np.all(rates > 1.9), f"misalignment should fall as dt^2, got {rates}"
    assert np.all(raw[:, 0] > 1e-2) and np.all(ren[:, 0] < 1e-5)                 # drift removed
    npt.assert_allclose(ren[:, 1], raw[:, 1], rtol=1e-6)                           # but misalignment unchanged
    # The misalignment agrees with the analytic estimate rate*T*(rate*dt)^2/3 to a few percent.
    for n, (_, mis) in zip((200, 400, 800), raw):
        assert mis == pytest.approx(rate * T * (rate * T / n) ** 2 / 3, rel=0.05)
    with pytest.raises(ValueError):
        impl.integrate_constant_rate(axis, rate, T, 0)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "dcm_integration_error.png").stat().st_size > 1000
