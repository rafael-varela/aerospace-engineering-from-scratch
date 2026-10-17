"""Public tests for P01.13 Systems of Linear Equations. The oracle is numpy.linalg; the assertions are properties.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-01/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(113)


def random_system(n: int) -> tuple[np.ndarray, np.ndarray]:
    return RNG.standard_normal((n, n)), RNG.standard_normal(n)


# ---------------------------------------------------------------- triangular solves

def test_back_and_forward_substitution_match_numpy_on_triangular_systems():
    for n in (1, 4, 9):
        A, b = random_system(n)
        U = np.triu(A) + n * np.eye(n)
        L = np.tril(A) + n * np.eye(n)
        npt.assert_allclose(impl.back_substitution(U, b), np.linalg.solve(U, b), rtol=1e-12, atol=1e-12)
        npt.assert_allclose(impl.forward_substitution(L, b), np.linalg.solve(L, b), rtol=1e-12, atol=1e-12)
    with pytest.raises(ValueError):
        impl.back_substitution(np.array([[1.0, 2.0], [0.0, 0.0]]), np.array([1.0, 1.0]))
    with pytest.raises(ValueError):
        impl.forward_substitution(np.array([[0.0, 0.0], [1.0, 2.0]]), np.array([1.0, 1.0]))


def test_substitution_does_not_modify_its_inputs():
    U = np.array([[2.0, 1.0], [0.0, 3.0]])
    c = np.array([1.0, 6.0])
    U0, c0 = U.copy(), c.copy()
    impl.back_substitution(U, c)
    npt.assert_array_equal(U, U0)
    npt.assert_array_equal(c, c0)


# ---------------------------------------------------------------- elimination

def test_forward_elimination_produces_an_upper_triangular_system_with_the_same_solution():
    A, b = random_system(6)
    U, c, perm = impl.forward_elimination(A, b)
    assert np.allclose(np.tril(U, -1), 0.0)
    assert sorted(perm.tolist()) == list(range(6))
    npt.assert_allclose(impl.back_substitution(U, c), np.linalg.solve(A, b), rtol=1e-10)
    # Elimination must not change A or b in place.
    A2, b2 = random_system(4)
    A3, b3 = A2.copy(), b2.copy()
    impl.forward_elimination(A2, b2)
    npt.assert_array_equal(A2, A3)
    npt.assert_array_equal(b2, b3)


def test_gaussian_elimination_matches_numpy_with_tiny_residual():
    for n in (2, 5, 20, 60):
        A, b = random_system(n)
        x = impl.gaussian_elimination(A, b)
        npt.assert_allclose(x, np.linalg.solve(A, b), rtol=1e-9, atol=1e-12)
        assert impl.residual_norm(A, x, b) < 1e-12


def test_identity_linearity_and_superposition():
    A, b = random_system(5)
    npt.assert_allclose(impl.gaussian_elimination(np.eye(5), b), b)
    x = impl.gaussian_elimination(A, b)
    npt.assert_allclose(impl.gaussian_elimination(A, 3.0 * b), 3.0 * x, rtol=1e-12)
    b2 = RNG.standard_normal(5)
    x2 = impl.gaussian_elimination(A, b2)
    npt.assert_allclose(impl.gaussian_elimination(A, b + b2), x + x2, rtol=1e-10)


def test_pivoting_handles_a_zero_pivot_and_a_tiny_one():
    A = np.array([[0.0, 1.0], [1.0, 0.0]])
    npt.assert_allclose(impl.gaussian_elimination(A, np.array([2.0, 3.0])), [3.0, 2.0])
    with pytest.raises(ValueError):
        impl.gaussian_elimination(A, np.array([2.0, 3.0]), pivot=False)
    T = np.array([[1e-13, 1.0], [1.0, 1.0]])
    t = np.array([1.0, 2.0])
    exact = np.linalg.solve(T, t)
    assert np.max(np.abs(impl.gaussian_elimination(T, t) - exact)) < 1e-12
    assert np.max(np.abs(impl.gaussian_elimination(T, t, pivot=False) - exact)) > 1e-6   # catastrophic without pivoting


def test_singular_and_malformed_systems_raise():
    with pytest.raises(ValueError):
        impl.gaussian_elimination(np.array([[1.0, 2.0], [2.0, 4.0]]), np.array([1.0, 2.0]))
    u = RNG.standard_normal(5)
    v = RNG.standard_normal(5)
    with pytest.raises(ValueError):
        impl.gaussian_elimination(np.outer(u, v) + np.outer(v, u), RNG.standard_normal(5))   # rank 2
    with pytest.raises(ValueError):
        impl.gaussian_elimination(np.ones((2, 3)), np.ones(2))
    with pytest.raises(ValueError):
        impl.gaussian_elimination(np.eye(3), np.ones(2))


# ---------------------------------------------------------------- LU, determinant, inverse

def test_lu_factors_reconstruct_pa_and_have_the_right_structure():
    for n in (3, 8, 25):
        A, _ = random_system(n)
        P, L, U = impl.lu_decompose(A)
        npt.assert_allclose(P @ A, L @ U, atol=1e-12)
        npt.assert_allclose(np.diag(L), 1.0)
        assert np.allclose(np.triu(L, 1), 0.0) and np.allclose(np.tril(U, -1), 0.0)
        assert np.allclose(P @ P.T, np.eye(n)) and np.allclose(P.sum(axis=0), 1.0) and set(np.unique(P)) <= {0.0, 1.0}
        assert np.all(np.abs(L) <= 1.0 + 1e-14)                 # partial pivoting bounds the multipliers


def test_lu_solve_agrees_with_elimination_for_several_right_hand_sides():
    A, _ = random_system(7)
    P, L, U = impl.lu_decompose(A)
    for _ in range(3):
        b = RNG.standard_normal(7)
        npt.assert_allclose(impl.lu_solve(P, L, U, b), impl.gaussian_elimination(A, b), rtol=1e-10)


def test_determinant_matches_numpy_and_obeys_its_rules():
    for n in (2, 5, 12):
        A, _ = random_system(n)
        assert impl.determinant(A) == pytest.approx(np.linalg.det(A), rel=1e-10)
    A, _ = random_system(4)
    B, _ = random_system(4)
    assert impl.determinant(A @ B) == pytest.approx(impl.determinant(A) * impl.determinant(B), rel=1e-10)
    assert impl.determinant(2.0 * A) == pytest.approx(16.0 * impl.determinant(A), rel=1e-10)
    assert impl.determinant(np.eye(6)) == pytest.approx(1.0)
    assert impl.determinant(np.array([[1.0, 2.0], [2.0, 4.0]])) == 0.0
    swapped = A[[1, 0, 2, 3]]
    assert impl.determinant(swapped) == pytest.approx(-impl.determinant(A), rel=1e-10)


def test_gauss_jordan_inverse_matches_numpy_and_inverts():
    for n in (1, 3, 10):
        A, _ = random_system(n)
        Ainv = impl.inverse_gauss_jordan(A)
        npt.assert_allclose(Ainv, np.linalg.inv(A), rtol=1e-9, atol=1e-12)
        npt.assert_allclose(A @ Ainv, np.eye(n), atol=1e-11)
        npt.assert_allclose(Ainv @ A, np.eye(n), atol=1e-11)
    with pytest.raises(ValueError):
        impl.inverse_gauss_jordan(np.array([[1.0, 2.0], [3.0, 6.0]]))


# ---------------------------------------------------------------- refinement, conditioning, cost

def test_iterative_refinement_repairs_an_unpivoted_solve_but_not_an_ill_conditioned_problem():
    T = np.array([[1e-13, 1.0], [1.0, 1.0]])
    t = np.array([1.0, 2.0])
    exact = np.linalg.solve(T, t)
    assert np.max(np.abs(impl.iterative_refinement(T, t, steps=2, pivot=False) - exact)) < 1e-12
    H = impl.hilbert(10)
    ones = np.ones(10)
    x = impl.gaussian_elimination(H, H @ ones)
    assert impl.residual_norm(H, x, H @ ones) < 1e-12          # the residual is tiny ...
    assert np.max(np.abs(x - ones)) > 1e-7                     # ... and the answer is still wrong
    x_ref = impl.iterative_refinement(H, H @ ones, steps=3)
    assert impl.residual_norm(H, x_ref, H @ ones) < 1e-12
    assert np.max(np.abs(x_ref - ones)) > 1e-9                 # refinement cannot fix conditioning


def test_hilbert_matrix_is_symmetric_with_growing_condition_number():
    for n in (3, 6, 9):
        H = impl.hilbert(n)
        assert H.shape == (n, n)
        npt.assert_array_equal(H, H.T)
        assert H[0, 0] == 1.0 and H[n - 1, n - 1] == pytest.approx(1.0 / (2 * n - 1))
    assert np.linalg.cond(impl.hilbert(9)) > 1e3 * np.linalg.cond(impl.hilbert(6))


def test_residual_norm_is_zero_for_the_exact_solution_and_scale_free():
    A, b = random_system(4)
    x = np.linalg.solve(A, b)
    assert impl.residual_norm(A, x, b) < 1e-13
    assert impl.residual_norm(A, np.zeros(4), b) == pytest.approx(1.0)
    assert impl.residual_norm(10 * A, x, 10 * b) == pytest.approx(impl.residual_norm(A, x, b), abs=1e-14)


def test_elimination_flops_scale_as_two_thirds_n_cubed():
    assert impl.elimination_flops(1) == pytest.approx(1.0)
    ratio = impl.elimination_flops(800) / impl.elimination_flops(400)
    assert 7.5 < ratio < 8.5
    assert impl.elimination_flops(1000) == pytest.approx(2.0 / 3.0 * 1000**3, rel=0.01)
    with pytest.raises(ValueError):
        impl.elimination_flops(0)


# ---------------------------------------------------------------- the truss

def test_truss_equilibrium_matrix_has_the_right_shape_and_newtons_third_law_in_every_member_column():
    nodes, members, supports = impl.two_bay_truss()
    A = impl.truss_equilibrium_matrix(nodes, members, supports)
    assert A.shape == (12, 12)
    m = len(members)
    for col in range(m):
        fx, fy = A[0::2, col].sum(), A[1::2, col].sum()        # internal force: no net force on the structure
        assert abs(fx) < 1e-14 and abs(fy) < 1e-14
        assert np.linalg.norm(A[:, col]) == pytest.approx(np.sqrt(2.0))   # +e at one end, -e at the other
    for col, (node, dof) in enumerate(supports):
        assert A[2 * node + dof, m + col] == 1.0 and np.count_nonzero(A[:, m + col]) == 1
    with pytest.raises(ValueError):
        impl.truss_equilibrium_matrix(nodes, [(0, 0)], supports)


def test_truss_solution_satisfies_global_equilibrium_and_joint_equilibrium():
    nodes, members, supports = impl.two_bay_truss()
    loads = {1: (0.0, -15000.0)}
    forces, reactions = impl.solve_truss(nodes, members, supports, loads)
    assert forces.shape == (9,) and reactions.shape == (3,)
    # Global equilibrium: reactions carry the load, no horizontal reaction, moment balance about node 0.
    assert reactions[0] == pytest.approx(0.0, abs=1e-9)
    assert reactions[1] + reactions[2] == pytest.approx(15000.0)
    assert reactions[2] * 4.0 == pytest.approx(15000.0 * 2.0)            # roller at x = 4 m, load at x = 2 m
    # Joint equilibrium: residual of the full system.
    A = impl.truss_equilibrium_matrix(nodes, members, supports)
    b = impl.truss_load_vector(len(nodes), loads)
    assert impl.residual_norm(A, np.concatenate([forces, reactions]), b) < 1e-12
    npt.assert_allclose(np.concatenate([forces, reactions]), np.linalg.solve(A, b), rtol=1e-10, atol=1e-8)
    # Hand check at node 0: reaction 7.5 kN up is balanced by the diagonal 0-4 (e = (0.8, 0.6)) in compression.
    assert forces[7] == pytest.approx(-7500.0 / 0.6)
    assert forces[0] == pytest.approx(-0.8 * forces[7])


def test_truss_forces_are_linear_in_the_load_and_zero_force_members_are_detected():
    nodes, members, supports = impl.two_bay_truss()
    f1, r1 = impl.solve_truss(nodes, members, supports, {1: (0.0, -15000.0)})
    f2, r2 = impl.solve_truss(nodes, members, supports, {1: (0.0, -30000.0)})
    npt.assert_allclose(f2, 2 * f1, atol=1e-8)
    npt.assert_allclose(r2, 2 * r1, atol=1e-8)
    # Node 3 joins only two non-collinear members and carries no load: both must be zero-force.
    assert abs(f1[2]) < 1e-9 and abs(f1[4]) < 1e-9
    assert impl.truss_load_vector(6, {1: (0.0, -15000.0)})[3] == 15000.0


def test_truss_with_a_missing_member_is_a_mechanism_and_raises():
    nodes, members, supports = impl.two_bay_truss()
    with pytest.raises(ValueError):
        impl.solve_truss(nodes, members[:-1], supports, {1: (0.0, -15000.0)})            # 11 unknowns, 12 equations
    with pytest.raises(ValueError):
        impl.solve_truss(nodes, members[:-1] + [(3, 5)], supports, {1: (0.0, -15000.0)})  # square count, but singular
