"""Public tests for P02.14 Tensor Notation & Index Gymnastics. They assert properties: the epsilon-delta
identity, index-notation cross product and determinant against NumPy, the transformation laws against
matrix forms and against the 2-D Mohr formulas, invariance of the stress invariants and of von Mises stress
under random rotations, principal stresses against ``numpy.linalg.eigh``, isotropic stiffness invariant
under every rotation (and an anisotropic one not), and Hooke's law round trips.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-02/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(214)
CASE = impl.fuselage_case()


def random_symmetric(scale: float = 100e6) -> np.ndarray:
    A = RNG.normal(size=(3, 3)) * scale
    return 0.5 * (A + A.T)


def random_rotation() -> np.ndarray:
    return impl.rotation_matrix(RNG.normal(size=3), RNG.uniform(-np.pi, np.pi))


# ---------------------------------------------------------------- the two symbols

def test_levi_civita_is_totally_antisymmetric_and_the_epsilon_delta_identity_holds():
    eps = impl.levi_civita()
    assert eps.shape == (3, 3, 3)
    assert eps[0, 1, 2] == 1.0 and eps[2, 1, 0] == -1.0 and eps[0, 0, 1] == 0.0
    npt.assert_array_equal(eps, -np.transpose(eps, (1, 0, 2)))
    npt.assert_array_equal(eps, -np.transpose(eps, (0, 2, 1)))
    npt.assert_array_equal(eps, np.transpose(eps, (1, 2, 0)))  # cyclic
    assert impl.epsilon_delta_identity_residual() == 0.0
    npt.assert_array_equal(impl.kronecker_delta(), np.eye(3))


def test_cross_product_and_determinant_from_epsilon_agree_with_numpy():
    for _ in range(10):
        a, b = RNG.normal(size=3), RNG.normal(size=3)
        npt.assert_allclose(impl.cross(a, b), np.cross(a, b), atol=1e-12)
        assert abs(impl.dot(a, impl.cross(a, b))) < 1e-12
        A = RNG.normal(size=(3, 3))
        npt.assert_allclose(impl.determinant(A), np.linalg.det(A), rtol=1e-12)
        npt.assert_allclose(impl.matrix_product(A, A.T), A @ A.T, rtol=1e-12)
        npt.assert_allclose(impl.trace(A), np.trace(A), rtol=1e-12)
        npt.assert_allclose(impl.double_dot(A, A), np.sum(A * A), rtol=1e-12)


def test_symmetric_and_antisymmetric_parts_and_the_axial_vector():
    A = RNG.normal(size=(3, 3))
    S, W = impl.symmetric_part(A), impl.antisymmetric_part(A)
    npt.assert_allclose(S + W, A, atol=1e-14)
    npt.assert_allclose(S, S.T, atol=1e-14)
    npt.assert_allclose(W, -W.T, atol=1e-14)
    assert abs(impl.double_dot(S, W)) < 1e-12  # symmetric : antisymmetric = 0
    w = impl.axial_vector(W)
    v = RNG.normal(size=3)
    npt.assert_allclose(W @ v, np.cross(w, v), atol=1e-12)


# ---------------------------------------------------------------- rotations

def test_rodrigues_matrix_is_a_proper_rotation_and_rotates_axes_passively():
    for _ in range(10):
        Q = random_rotation()
        assert impl.is_orthogonal(Q)
        npt.assert_allclose(Q @ Q.T, np.eye(3), atol=1e-12)
        npt.assert_allclose(np.linalg.det(Q), 1.0, rtol=1e-12)
    # Axes rotated by +90 deg about x3: the old e1 is the new -e2' direction, i.e. v = e1 has components (0, -1, 0).
    Q = impl.rotation_matrix([0, 0, 1], np.pi / 2)
    npt.assert_allclose(impl.transform_vector(np.array([1.0, 0, 0]), Q), [0, -1, 0], atol=1e-12)
    # The axis itself is unchanged and rotating by 2 pi is the identity.
    n = RNG.normal(size=3)
    n /= np.linalg.norm(n)
    npt.assert_allclose(impl.transform_vector(n, impl.rotation_matrix(n, 1.3)), n, atol=1e-12)
    npt.assert_allclose(impl.rotation_matrix(n, 2 * np.pi), np.eye(3), atol=1e-12)
    with pytest.raises(ValueError):
        impl.rotation_matrix([0, 0, 0], 1.0)


def test_transformation_laws_match_matrix_forms_and_compose():
    S = random_symmetric()
    v = RNG.normal(size=3)
    Q1, Q2 = random_rotation(), random_rotation()
    npt.assert_allclose(impl.transform_vector(v, Q1), Q1 @ v, rtol=1e-12)
    npt.assert_allclose(impl.transform_tensor2(S, Q1), Q1 @ S @ Q1.T, rtol=1e-12)
    # Two successive rotations equal one rotation by the product.
    npt.assert_allclose(impl.transform_tensor2(impl.transform_tensor2(S, Q1), Q2),
                        impl.transform_tensor2(S, Q2 @ Q1), rtol=1e-12)
    # Lengths and inner products are preserved.
    npt.assert_allclose(impl.dot(impl.transform_vector(v, Q1), impl.transform_vector(v, Q1)), impl.dot(v, v), rtol=1e-12)
    # The identity leaves everything alone and the inverse rotation undoes it.
    npt.assert_allclose(impl.transform_tensor2(impl.transform_tensor2(S, Q1), Q1.T), S, rtol=1e-12)
    C = impl.isotropic_stiffness(70e9, 0.3)
    npt.assert_allclose(impl.transform_tensor4(C, np.eye(3)), C, rtol=1e-12)


def test_traction_transforms_as_a_vector_so_cauchy_formula_is_frame_independent():
    S = random_symmetric()
    n = RNG.normal(size=3)
    n /= np.linalg.norm(n)
    Q = random_rotation()
    t_old = impl.traction(S, n)
    t_new = impl.traction(impl.transform_tensor2(S, Q), impl.transform_vector(n, Q))
    npt.assert_allclose(t_new, impl.transform_vector(t_old, Q), rtol=1e-12)
    sn, tau = impl.normal_and_shear_on_plane(S, n)
    npt.assert_allclose(sn ** 2 + tau ** 2, impl.dot(t_old, t_old), rtol=1e-12)
    with pytest.raises(ValueError):
        impl.traction(RNG.normal(size=(3, 3)) * 1e6, n)


# ---------------------------------------------------------------- invariants, principal stresses

def test_invariants_are_unchanged_by_rotation_and_agree_with_the_characteristic_polynomial():
    S = random_symmetric()
    I1, I2, I3 = impl.stress_invariants(S)
    for _ in range(5):
        J1, J2, J3 = impl.stress_invariants(impl.transform_tensor2(S, random_rotation()))
        npt.assert_allclose([J1, J2, J3], [I1, I2, I3], rtol=1e-9)
    vals = np.linalg.eigvalsh(S)
    npt.assert_allclose(I1, vals.sum(), rtol=1e-12)
    npt.assert_allclose(I2, vals[0] * vals[1] + vals[1] * vals[2] + vals[2] * vals[0], rtol=1e-9)
    npt.assert_allclose(I3, np.prod(vals), rtol=1e-9)
    # Every principal stress is a root of s^3 - I1 s^2 + I2 s - I3.
    for s in vals:
        assert abs(s ** 3 - I1 * s ** 2 + I2 * s - I3) < 1e-9 * max(abs(vals)) ** 3


def test_principal_stresses_match_eigh_and_their_planes_carry_no_shear():
    S = random_symmetric()
    vals, vecs = impl.principal_stresses(S)
    npt.assert_allclose(vals, np.sort(np.linalg.eigvalsh(S))[::-1], rtol=1e-12)
    npt.assert_allclose(vecs.T @ vecs, np.eye(3), atol=1e-12)
    npt.assert_allclose(np.linalg.det(vecs), 1.0, rtol=1e-12)
    for k in range(3):
        sn, tau = impl.normal_and_shear_on_plane(S, vecs[:, k])
        npt.assert_allclose(sn, vals[k], rtol=1e-10)
        assert tau < 1e-6 * np.max(np.abs(vals))
    # In the principal frame (rows of Q = principal directions) the tensor is diagonal.
    npt.assert_allclose(impl.transform_tensor2(S, vecs.T), np.diag(vals), atol=1e-6 * np.max(np.abs(vals)))


def test_deviator_von_mises_and_max_shear_properties():
    S = random_symmetric()
    p, s = impl.hydrostatic_and_deviatoric(S)
    assert abs(impl.trace(s)) < 1e-6 * np.max(np.abs(S))
    npt.assert_allclose(p * np.eye(3) + s, S, rtol=1e-12)
    vm = impl.von_mises(S)
    for _ in range(5):
        npt.assert_allclose(impl.von_mises(impl.transform_tensor2(S, random_rotation())), vm, rtol=1e-9)
    # Hydrostatic pressure adds nothing to von Mises; uniaxial gives sigma; pure shear gives sqrt(3) tau.
    npt.assert_allclose(impl.von_mises(S + 50e6 * np.eye(3)), vm, rtol=1e-9)
    npt.assert_allclose(impl.von_mises(np.diag([200e6, 0, 0])), 200e6, rtol=1e-12)
    tau = 30e6
    npt.assert_allclose(impl.von_mises(np.array([[0, tau, 0], [tau, 0, 0], [0, 0, 0]])), np.sqrt(3) * tau, rtol=1e-12)
    npt.assert_allclose(impl.max_shear(np.array([[0, tau, 0], [tau, 0, 0], [0, 0, 0]])), tau, rtol=1e-12)
    # max shear is the largest shear found on any of many random planes.
    tmax = impl.max_shear(S)
    shears = []
    for _ in range(400):
        n = RNG.normal(size=3)
        shears.append(impl.normal_and_shear_on_plane(S, n / np.linalg.norm(n))[1])
    assert max(shears) <= tmax * (1 + 1e-9) and max(shears) > 0.9 * tmax


def test_mohr_formulas_are_the_tensor_law_restricted_to_the_plane():
    for _ in range(10):
        sx, sy, txy = RNG.normal(size=3) * 80e6
        theta = RNG.uniform(-np.pi, np.pi)
        S = np.array([[sx, txy, 0.0], [txy, sy, 0.0], [0.0, 0.0, 0.0]])
        Sp = impl.transform_tensor2(S, impl.rotation_matrix([0, 0, 1], theta))
        sxp, syp, txyp = impl.mohr_2d(sx, sy, txy, theta)
        npt.assert_allclose([sxp, syp, txyp], [Sp[0, 0], Sp[1, 1], Sp[0, 1]], rtol=1e-10, atol=1e-3)
        # Mohr's circle: the rotated point is at the same distance from the centre.
        c = 0.5 * (sx + sy)
        R = np.hypot(0.5 * (sx - sy), txy)
        npt.assert_allclose(np.hypot(sxp - c, txyp), R, rtol=1e-10)


# ---------------------------------------------------------------- strain and Hooke

def test_strain_from_displacement_gradient_and_rigid_rotation_gives_no_strain():
    H = RNG.normal(size=(3, 3)) * 1e-3
    eps, omega = impl.strain_from_displacement_gradient(H)
    npt.assert_allclose(eps + omega, H, atol=1e-16)
    npt.assert_allclose(eps, eps.T, atol=1e-16)
    # A small rigid rotation u = w x r has gradient H_ij = -eps_ijk w_k: antisymmetric, zero strain.
    w = RNG.normal(size=3) * 1e-3
    Hrot = -np.einsum("ijk,k->ij", impl.levi_civita(), w)
    e2, o2 = impl.strain_from_displacement_gradient(Hrot)
    assert np.max(np.abs(e2)) < 1e-16
    npt.assert_allclose(impl.axial_vector(o2), w, atol=1e-16)


def test_isotropic_stiffness_is_invariant_under_every_rotation_but_an_anisotropic_one_is_not():
    C = impl.isotropic_stiffness(73e9, 0.33)
    assert C.shape == (3, 3, 3, 3)
    # minor and major symmetries
    npt.assert_allclose(C, np.transpose(C, (1, 0, 2, 3)), atol=1e-3)
    npt.assert_allclose(C, np.transpose(C, (0, 1, 3, 2)), atol=1e-3)
    npt.assert_allclose(C, np.transpose(C, (2, 3, 0, 1)), atol=1e-3)
    for _ in range(5):
        npt.assert_allclose(impl.transform_tensor4(C, random_rotation()), C, rtol=1e-9, atol=1e-3)
    Ca = C.copy()
    Ca[0, 0, 0, 0] *= 1.5
    assert np.max(np.abs(impl.transform_tensor4(Ca, random_rotation()) - Ca)) > 1e-3 * np.max(np.abs(Ca))
    with pytest.raises(ValueError):
        impl.isotropic_stiffness(73e9, 0.5)
    with pytest.raises(ValueError):
        impl.isotropic_stiffness(-1.0, 0.3)


def test_hooke_law_round_trips_and_reproduces_uniaxial_and_shear_moduli():
    E, nu = 73e9, 0.33
    C = impl.isotropic_stiffness(E, nu)
    S = random_symmetric()
    eps = impl.isotropic_strain(S, E, nu)
    npt.assert_allclose(impl.hooke_stress(C, eps), S, rtol=1e-9)
    # Uniaxial stress: eps_11 = sigma/E, eps_22 = eps_33 = -nu sigma/E.
    e = impl.isotropic_strain(np.diag([100e6, 0, 0]), E, nu)
    npt.assert_allclose([e[0, 0], e[1, 1], e[2, 2]], [100e6 / E, -nu * 100e6 / E, -nu * 100e6 / E], rtol=1e-12)
    # Pure shear: sigma_12 = 2 mu eps_12 = G gamma.
    G = E / (2 * (1 + nu))
    eps12 = np.array([[0, 1e-3, 0], [1e-3, 0, 0], [0, 0, 0]])
    npt.assert_allclose(impl.hooke_stress(C, eps12)[0, 1], 2 * G * 1e-3, rtol=1e-12)
    # Hooke's law is frame independent: rotating eps and C gives the rotated sigma.
    Q = random_rotation()
    npt.assert_allclose(impl.hooke_stress(impl.transform_tensor4(C, Q), impl.transform_tensor2(eps, Q)),
                        impl.transform_tensor2(impl.hooke_stress(C, eps), Q), rtol=1e-9)


# ---------------------------------------------------------------- the fuselage skin

def test_fuselage_skin_stress_scales_correctly_and_rejects_thick_walls():
    p, r, t, T = CASE["p"], CASE["r"], CASE["t"], CASE["torque"]
    S = impl.fuselage_skin_stress(p, r, t, T)
    npt.assert_allclose(S, S.T)
    npt.assert_allclose(S[1, 1], 2 * S[0, 0], rtol=1e-12)        # hoop = twice longitudinal
    npt.assert_allclose(S[1, 1], p * r / t, rtol=1e-12)
    npt.assert_allclose(S[0, 1], T / (2 * np.pi * r ** 2 * t), rtol=1e-12)
    assert S[2, 2] == 0.0 and S[0, 2] == 0.0 and S[1, 2] == 0.0
    # Doubling the pressure doubles the normal stresses and leaves the shear; doubling t halves everything.
    S2 = impl.fuselage_skin_stress(2 * p, r, t, T)
    npt.assert_allclose(np.diag(S2), 2 * np.diag(S), rtol=1e-12)
    npt.assert_allclose(S2[0, 1], S[0, 1], rtol=1e-12)
    npt.assert_allclose(impl.fuselage_skin_stress(p, r, 2 * t, T), 0.5 * S, rtol=1e-12)
    # No torque: axes 1-2 are principal and von Mises is sqrt(3)/2 of the hoop stress.
    S0 = impl.fuselage_skin_stress(p, r, t, 0.0)
    npt.assert_allclose(impl.von_mises(S0), np.sqrt(3) / 2 * S0[1, 1], rtol=1e-12)
    with pytest.raises(ValueError):
        impl.fuselage_skin_stress(p, r, 0.3 * r, T)
    with pytest.raises(ValueError):
        impl.fuselage_skin_stress(p, -r, t, T)


def test_seam_stresses_match_the_tensor_law_and_are_periodic_in_the_seam_angle():
    S = impl.fuselage_skin_stress(CASE["p"], CASE["r"], CASE["t"], CASE["torque"])
    a = np.radians(CASE["seam_angle_deg"])
    sn, tau = impl.seam_stresses(S, a)
    # A seam along axis 1 (a = 0) is opened by the hoop stress; one around the hoop (a = 90 deg) by the longitudinal stress.
    npt.assert_allclose(impl.seam_stresses(S, 0.0)[0], S[1, 1], rtol=1e-12)
    npt.assert_allclose(impl.seam_stresses(S, np.pi / 2)[0], S[0, 0], rtol=1e-12)
    npt.assert_allclose(impl.seam_stresses(S, 0.0)[1], abs(S[0, 1]), rtol=1e-12)
    # Same from rotating the tensor to seam axes.
    Sp = impl.transform_tensor2(S, impl.rotation_matrix([0, 0, 1], a))
    npt.assert_allclose(sn, Sp[1, 1], rtol=1e-10)
    npt.assert_allclose(tau, abs(Sp[0, 1]), rtol=1e-10)
    # Periodic in 180 deg; the normal stress is bounded by the principal stresses.
    npt.assert_allclose(impl.seam_stresses(S, a + np.pi), (sn, tau), rtol=1e-10)
    vals, _ = impl.principal_stresses(S)
    for ang in np.linspace(0, np.pi, 37):
        s_ang, t_ang = impl.seam_stresses(S, ang)
        assert vals[2] - 1e-6 <= s_ang <= vals[0] + 1e-6
        assert t_ang <= impl.max_shear(S) * (1 + 1e-9)


def test_yield_margin_is_frame_independent_and_signs_are_right():
    S = impl.fuselage_skin_stress(CASE["p"], CASE["r"], CASE["t"], CASE["torque"])
    m = impl.yield_margin(S, CASE["sigma_yield"])
    assert m > 0.0
    npt.assert_allclose(impl.yield_margin(impl.transform_tensor2(S, random_rotation()), CASE["sigma_yield"]), m, rtol=1e-9)
    npt.assert_allclose(impl.yield_margin(np.diag([325e6, 0, 0]), 325e6), 0.0, atol=1e-12)
    assert impl.yield_margin(np.diag([400e6, 0, 0]), 325e6) < 0.0


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "stress_rotation.png").stat().st_size > 1000
