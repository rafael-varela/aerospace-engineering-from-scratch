"""Public tests for P01.16 Eigenvalues & Eigenvectors. They assert properties, with numpy.linalg as the oracle.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-01/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.linalg import expm

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(1616)
MASSES, POSITIONS = impl.satellite_point_masses()



def sat_tensor() -> np.ndarray:
    """The satellite inertia tensor about its centre of mass (built inside each test, not at import time)."""
    return impl.inertia_tensor(MASSES, POSITIONS)


def random_symmetric(n: int, spread: float = 10.0) -> np.ndarray:
    B = RNG.normal(size=(n, n))
    return B + B.T + spread * np.eye(n)


# ---------------------------------------------------------------- the inertia tensor

def test_inertia_tensor_is_symmetric_positive_with_correct_diagonal():
    I_SAT = sat_tensor()
    I = I_SAT
    npt.assert_allclose(I, I.T, atol=1e-12)
    cm = impl.center_of_mass(MASSES, POSITIONS)
    d = POSITIONS - cm
    # Diagonal entries are sum m (y^2 + z^2) etc.; off-diagonals are -sum m x y etc.
    assert I[0, 0] == pytest.approx(np.sum(MASSES * (d[:, 1] ** 2 + d[:, 2] ** 2)))
    assert I[0, 1] == pytest.approx(-np.sum(MASSES * d[:, 0] * d[:, 1]))
    assert I[1, 2] == pytest.approx(-np.sum(MASSES * d[:, 1] * d[:, 2]))
    assert np.all(np.diag(I) > 0.0) and np.any(np.abs(I - np.diag(np.diag(I))) > 1e-3)   # products of inertia present
    assert np.all(np.linalg.eigvalsh(I) > 0.0)                                            # positive definite


def test_parallel_axis_theorem_and_input_validation():
    cm = impl.center_of_mass(MASSES, POSITIONS)
    shift = np.array([0.5, -0.2, 0.3])
    I_cm = impl.inertia_tensor(MASSES, POSITIONS)
    I_o = impl.inertia_tensor(MASSES, POSITIONS, about=cm + shift)
    M = MASSES.sum()
    npt.assert_allclose(I_o, I_cm + M * (shift @ shift * np.eye(3) - np.outer(shift, shift)), rtol=1e-12)
    assert np.all(np.diag(I_o) >= np.diag(I_cm))                                         # CM gives the minimum
    with pytest.raises(ValueError):
        impl.inertia_tensor(np.array([1.0, -1.0]), np.zeros((2, 3)))
    with pytest.raises(ValueError):
        impl.inertia_tensor(np.array([1.0, 1.0]), np.zeros((3, 3)))


def test_rayleigh_quotient_is_the_moment_about_the_axis():
    I_SAT = sat_tensor()
    for _ in range(5):
        n = RNG.normal(size=3)
        assert impl.rayleigh_quotient(I_SAT, n) == pytest.approx(impl.moment_about_axis(MASSES, POSITIONS, n), rel=1e-12)
    lam = np.linalg.eigvalsh(I_SAT)
    for _ in range(20):
        q = impl.rayleigh_quotient(I_SAT, RNG.normal(size=3))
        assert lam[0] - 1e-12 <= q <= lam[-1] + 1e-12                                     # bounded by the extreme eigenvalues


# ---------------------------------------------------------------- characteristic polynomial, closed form

def test_characteristic_coefficients_match_numpy_poly_and_the_invariants():
    I_SAT = sat_tensor()
    for A in (I_SAT, random_symmetric(3), RNG.normal(size=(3, 3))):
        c = impl.characteristic_coefficients(A)
        npt.assert_allclose(c, np.poly(A), rtol=1e-12, atol=1e-12)
        lam = np.linalg.eigvals(A)
        assert -c[1] == pytest.approx(np.trace(A))
        assert -c[3] == pytest.approx(np.linalg.det(A))
        assert c[2] == pytest.approx(np.real(lam[0] * lam[1] + lam[0] * lam[2] + lam[1] * lam[2]))
        # Every eigenvalue is a root of the polynomial.
        for l in lam:
            assert abs(np.polyval(c, l)) < 1e-9 * max(1.0, abs(l) ** 3)
    with pytest.raises(ValueError):
        impl.characteristic_coefficients(np.eye(2))


def test_closed_form_symmetric_eigenvalues_match_eigvalsh():
    I_SAT = sat_tensor()
    for A in (I_SAT, random_symmetric(3), random_symmetric(3, spread=0.0), np.diag([3.0, 1.0, 2.0])):
        npt.assert_allclose(impl.symmetric_eigenvalues_3x3(A), np.linalg.eigvalsh(A), rtol=1e-10, atol=1e-10)
    npt.assert_allclose(impl.symmetric_eigenvalues_3x3(4.0 * np.eye(3)), [4.0, 4.0, 4.0])
    with pytest.raises(ValueError):
        impl.symmetric_eigenvalues_3x3(np.array([[1.0, 2.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]))


# ---------------------------------------------------------------- power and inverse iteration

def test_power_iteration_finds_the_largest_eigenpair():
    I_SAT = sat_tensor()
    lam, v, n = impl.power_iteration(I_SAT)
    assert lam == pytest.approx(np.max(np.linalg.eigvalsh(I_SAT)), rel=1e-10)
    assert np.linalg.norm(v) == pytest.approx(1.0)
    assert np.linalg.norm(I_SAT @ v - lam * v) < 1e-8 * lam                             # residual of the eigen-equation
    assert 1 < n < 10000
    # Convergence is linear with ratio lambda_2 / lambda_1: a tighter tolerance costs log(tol)/log(ratio) more steps.
    lam_s = np.sort(np.linalg.eigvalsh(I_SAT))
    ratio = lam_s[1] / lam_s[2]
    _, _, n_loose = impl.power_iteration(I_SAT, tol=1e-6)
    _, _, n_tight = impl.power_iteration(I_SAT, tol=1e-12)
    assert n_tight - n_loose == pytest.approx(np.log(1e-6) / np.log(ratio), rel=0.15)


def test_power_iteration_handles_negative_dominant_eigenvalue_and_raises_when_stuck():
    A = np.diag([-5.0, 2.0, 1.0])
    lam, v, _ = impl.power_iteration(A, x0=np.array([1.0, 1.0, 1.0]))
    assert lam == pytest.approx(-5.0)
    assert abs(v[0]) == pytest.approx(1.0)
    # Two eigenvalues of equal magnitude and opposite sign: the iterate never settles.
    with pytest.raises(RuntimeError):
        impl.power_iteration(np.diag([3.0, -3.0, 1.0]), x0=np.array([1.0, 1.0, 1.0]), max_iter=500)


def test_shifted_inverse_iteration_targets_any_eigenvalue():
    I_SAT = sat_tensor()
    lam_s = np.sort(np.linalg.eigvalsh(I_SAT))
    for target in lam_s:
        lam, v, n = impl.shifted_inverse_iteration(I_SAT, target + 0.3)
        assert lam == pytest.approx(target, rel=1e-10)
        assert np.linalg.norm(I_SAT @ v - lam * v) < 1e-8 * lam
    lam0, _, n0 = impl.shifted_inverse_iteration(I_SAT, 0.0)
    assert lam0 == pytest.approx(lam_s[0], rel=1e-10)                                     # shift 0 gives the smallest
    with pytest.raises(RuntimeError):
        impl.shifted_inverse_iteration(np.diag([1.0, 2.0, 3.0]), 2.0)                     # exact eigenvalue: singular


# ---------------------------------------------------------------- Jacobi rotations

@pytest.mark.parametrize("n", [3, 6, 12])
def test_jacobi_matches_eigh_and_gives_orthonormal_eigenvectors(n):
    A = random_symmetric(n)
    lam, V, sweeps = impl.jacobi_eigen(A)
    npt.assert_allclose(lam, np.linalg.eigvalsh(A), rtol=1e-10, atol=1e-10)
    npt.assert_allclose(V.T @ V, np.eye(n), atol=1e-12)
    npt.assert_allclose(A @ V, V @ np.diag(lam), atol=1e-9)
    npt.assert_allclose(V @ np.diag(lam) @ V.T, A, atol=1e-9)                             # spectral decomposition
    assert np.all(np.diff(lam) >= 0.0)
    # Eigenvectors agree with numpy up to sign.
    _, V_np = np.linalg.eigh(A)
    npt.assert_allclose(np.abs(np.sum(V * V_np, axis=0)), 1.0, atol=1e-8)
    assert sweeps <= 12


def test_jacobi_converges_quadratically_and_preserves_the_invariants():
    A = random_symmetric(8)
    hist = []
    lam, V, _ = impl.jacobi_eigen(A, tol=1e-13, history=hist)
    hist = np.array(hist)
    assert np.sum(lam) == pytest.approx(np.trace(A))                                      # trace invariant
    assert np.prod(lam) == pytest.approx(np.linalg.det(A), rel=1e-9)                      # determinant invariant
    assert np.sum(lam ** 2) == pytest.approx(np.sum(A ** 2))                              # Frobenius invariant
    # The off-diagonal norm falls monotonically, and the reduction factor per sweep itself shrinks
    # (quadratic convergence), so the final sweep gains far more than the first.
    assert np.all(np.diff(hist) < 0.0)
    ratios = hist[1:] / hist[:-1]
    assert len(ratios) >= 2 and ratios[-1] < 0.1 * ratios[0]
    with pytest.raises(ValueError):
        impl.jacobi_eigen(np.array([[1.0, 2.0], [0.0, 1.0]]))
    npt.assert_allclose(impl.off_diagonal_norm(np.array([[1.0, 2.0], [2.0, 1.0]])), np.sqrt(8.0))


def test_jacobi_on_a_diagonal_matrix_takes_one_sweep_and_sorts():
    lam, V, sweeps = impl.jacobi_eigen(np.diag([3.0, 1.0, 2.0]))
    npt.assert_allclose(lam, [1.0, 2.0, 3.0])
    assert sweeps == 1
    npt.assert_allclose(np.abs(V), np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]]), atol=1e-15)


# ---------------------------------------------------------------- principal axes

def test_principal_axes_diagonalize_the_tensor_with_a_proper_rotation():
    I_SAT = sat_tensor()
    moments, C = impl.principal_axes(I_SAT)
    npt.assert_allclose(C.T @ C, np.eye(3), atol=1e-12)
    assert np.linalg.det(C) == pytest.approx(1.0, abs=1e-12)
    I_p = impl.rotate_tensor(I_SAT, C)
    npt.assert_allclose(I_p, np.diag(moments), atol=1e-9)
    npt.assert_allclose(moments, np.linalg.eigvalsh(I_SAT), rtol=1e-10)
    # Each row is a principal axis: the moment about it is its eigenvalue (geometry check, no linear algebra).
    for k in range(3):
        assert impl.moment_about_axis(MASSES, POSITIONS, C[k]) == pytest.approx(moments[k], rel=1e-10)
    # Rotating back recovers the body tensor.
    npt.assert_allclose(impl.rotate_tensor(I_p, C.T), I_SAT, atol=1e-9)


def test_principal_moments_are_rotation_invariant():
    I_SAT = sat_tensor()
    """Expressing the same body in another frame changes every entry but not the principal moments."""
    angle = 0.7
    c, s = np.cos(angle), np.sin(angle)
    R = np.array([[c, 0, -s], [0, 1, 0], [s, 0, c]]) @ np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])
    I_rot = impl.rotate_tensor(I_SAT, R)
    assert not np.allclose(I_rot, I_SAT)
    npt.assert_allclose(impl.principal_axes(I_rot)[0], impl.principal_axes(I_SAT)[0], rtol=1e-10)
    npt.assert_allclose(impl.characteristic_coefficients(I_rot), impl.characteristic_coefficients(I_SAT), rtol=1e-10)


def test_triangle_inequalities_and_spin_stability_rules():
    I_SAT = sat_tensor()
    moments, _ = impl.principal_axes(I_SAT)
    assert impl.triangle_inequalities_hold(moments)
    assert impl.triangle_inequalities_hold([1.0, 1.0, 2.0])                               # flat plate: equality
    assert not impl.triangle_inequalities_hold([1.0, 1.0, 2.5])                           # no body can do this
    assert not impl.triangle_inequalities_hold([-1.0, 2.0, 2.0])
    # Point masses always give a realizable tensor, whatever the geometry.
    for _ in range(5):
        I_rand = impl.inertia_tensor(RNG.uniform(0.5, 5.0, size=6), RNG.normal(size=(6, 3)))
        assert impl.triangle_inequalities_hold(np.linalg.eigvalsh(I_rand))
    assert impl.spin_axis_stability(moments, 2) == {"rigid": True, "with_dissipation": True}
    assert impl.spin_axis_stability(moments, 0) == {"rigid": True, "with_dissipation": False}
    assert impl.spin_axis_stability(moments, 1) == {"rigid": False, "with_dissipation": False}
    with pytest.raises(ValueError):
        impl.spin_axis_stability([1.0, 1.0, 2.0], 0)
    with pytest.raises(ValueError):
        impl.spin_axis_stability(moments, 3)


# ---------------------------------------------------------------- functions of a matrix, stability

def test_matrix_function_reproduces_powers_square_roots_and_the_exponential():
    A = random_symmetric(5)
    npt.assert_allclose(impl.matrix_function_symmetric(A, lambda l: l ** 2), A @ A, rtol=1e-10)
    S = impl.matrix_function_symmetric(A, np.sqrt)
    npt.assert_allclose(S @ S, A, rtol=1e-10)
    npt.assert_allclose(impl.matrix_function_symmetric(A, np.exp), expm(A), rtol=1e-9)
    npt.assert_allclose(impl.matrix_function_symmetric(A, lambda l: 1.0 / l), np.linalg.inv(A), rtol=1e-10)


def test_decay_rates_of_a_damped_oscillator_match_the_analytic_eigenvalues():
    w, zeta = 2.0, 0.1
    A = np.array([[0.0, 1.0], [-w * w, -2.0 * zeta * w]])
    rates = impl.linear_system_decay_rates(A)
    npt.assert_allclose(rates, [-zeta * w, -zeta * w])
    # The actual solution decays at that rate: |u(t)| ~ exp(-zeta w t).
    t = 5.0
    u = expm(A * t) @ np.array([1.0, 0.0])
    assert np.linalg.norm(u) < 1.5 * np.exp(-zeta * w * t)
    # Negative damping flips the sign: unstable.
    assert np.all(impl.linear_system_decay_rates(np.array([[0.0, 1.0], [-w * w, +2.0 * zeta * w]])) > 0.0)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "eigen_convergence.png").stat().st_size > 1000
