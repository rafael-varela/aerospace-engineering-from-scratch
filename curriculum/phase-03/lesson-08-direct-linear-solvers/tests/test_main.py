"""Public tests for P03.08 Direct Linear Solvers. They assert factorisation identities, residuals,
limiting cases, conditioning laws and a convergence order, with ``scipy.linalg`` as the oracle.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-03/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy import linalg as sla

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(308)


def spd_matrix(n):
    B = RNG.standard_normal((n, n))
    return B @ B.T + n * np.eye(n)


# ---------------------------------------------------------------- triangular solves

def test_triangular_substitutions_match_scipy():
    n = 7
    L = np.tril(RNG.standard_normal((n, n))) + 3 * np.eye(n)
    U = np.triu(RNG.standard_normal((n, n))) + 3 * np.eye(n)
    b = RNG.standard_normal(n)
    npt.assert_allclose(impl.forward_substitution(L, b), sla.solve_triangular(L, b, lower=True), rtol=1e-12)
    npt.assert_allclose(impl.back_substitution(U, b), sla.solve_triangular(U, b, lower=False), rtol=1e-12)
    L1 = np.tril(RNG.standard_normal((n, n)), -1) + np.eye(n)
    y = impl.forward_substitution(L1, b, unit_diagonal=True)
    npt.assert_allclose(L1 @ y, b, atol=1e-12)
    with pytest.raises(ValueError):
        impl.back_substitution(np.array([[1.0, 2.0], [0.0, 0.0]]), np.ones(2))


# ---------------------------------------------------------------- LU

@pytest.mark.parametrize("n", [1, 3, 6, 12])
def test_lu_reconstructs_the_permuted_matrix_and_matches_scipy(n):
    A = RNG.standard_normal((n, n))
    L, U, perm = impl.lu_decompose(A)
    npt.assert_allclose(L @ U, A[perm], atol=1e-13)
    npt.assert_allclose(np.diag(L), 1.0)
    assert np.all(np.abs(np.tril(L, -1)) <= 1.0 + 1e-15)        # partial pivoting: multipliers <= 1
    assert np.allclose(np.tril(U, -1), 0.0) and np.allclose(np.triu(L, 1), 0.0)
    P, Ls, Us = sla.lu(A)
    npt.assert_allclose(L, Ls, atol=1e-12)                       # same pivots, same factors
    npt.assert_allclose(U, Us, atol=1e-12)
    b = RNG.standard_normal(n)
    x = impl.lu_solve(L, U, perm, b)
    npt.assert_allclose(x, sla.solve(A, b), rtol=1e-10, atol=1e-12)
    assert impl.backward_error(A, x, b) < 1e-15


def test_pivoting_is_needed_and_singular_matrices_are_refused():
    A = np.array([[0.0, 1.0], [1.0, 0.0]])
    with pytest.raises(ValueError):
        impl.lu_decompose(A, pivot=False)
    L, U, perm = impl.lu_decompose(A)
    npt.assert_allclose(impl.lu_solve(L, U, perm, np.array([2.0, 3.0])), [3.0, 2.0])
    # A tiny pivot without pivoting loses everything; with pivoting it is harmless.
    eps = 1e-17
    A = np.array([[eps, 1.0], [1.0, 1.0]])
    b = np.array([1.0, 2.0])
    L, U, perm = impl.lu_decompose(A, pivot=False)
    x_bad = impl.lu_solve(L, U, perm, b)
    L, U, perm = impl.lu_decompose(A)
    x_good = impl.lu_solve(L, U, perm, b)
    npt.assert_allclose(x_good, [1.0, 1.0], rtol=1e-12)
    assert np.max(np.abs(x_bad - [1.0, 1.0])) > 0.1
    with pytest.raises(ValueError):
        impl.lu_decompose(np.ones((3, 3)))
    with pytest.raises(ValueError):
        impl.lu_decompose(np.ones((2, 3)))
    # Several right-hand sides at once.
    A = RNG.standard_normal((5, 5)); B = RNG.standard_normal((5, 3))
    npt.assert_allclose(impl.solve(A, B), sla.solve(A, B), rtol=1e-10)


def test_determinant_inverse_and_condition_number():
    A = RNG.standard_normal((6, 6))
    L, U, perm = impl.lu_decompose(A)
    assert impl.determinant_from_lu(U, perm) == pytest.approx(np.linalg.det(A), rel=1e-10)
    A2 = A[[1, 0, 2, 3, 4, 5]]                                    # one row swap flips the sign
    L2, U2, perm2 = impl.lu_decompose(A2)
    assert impl.determinant_from_lu(U2, perm2) == pytest.approx(-np.linalg.det(A), rel=1e-10)
    npt.assert_allclose(impl.inverse_from_lu(L, U, perm), np.linalg.inv(A), rtol=1e-9, atol=1e-12)
    assert impl.condition_number_1(A) == pytest.approx(np.linalg.cond(A, 1), rel=1e-9)
    assert impl.condition_number_1(np.eye(4)) == pytest.approx(1.0)
    assert impl.norm1(np.array([[1.0, -2.0], [3.0, 4.0]])) == 6.0


def test_growth_factor_is_two_to_the_n_minus_one_for_wilkinson_and_small_for_random():
    for n in (4, 6, 10):
        assert impl.growth_factor(impl.wilkinson_matrix(n)) == pytest.approx(2.0 ** (n - 1))
    W = impl.wilkinson_matrix(5)
    assert np.all(np.diag(W) == 1.0) and np.all(W[:, -1] == 1.0) and np.all(np.tril(W, -1)[:, :-1] <= 0.0)
    assert impl.growth_factor(RNG.standard_normal((50, 50))) < 10.0


# ---------------------------------------------------------------- Cholesky

@pytest.mark.parametrize("n", [1, 4, 9])
def test_cholesky_reconstructs_and_matches_scipy(n):
    A = spd_matrix(n)
    L = impl.cholesky(A)
    npt.assert_allclose(L @ L.T, A, rtol=1e-12, atol=1e-12)
    assert np.allclose(np.triu(L, 1), 0.0) and np.all(np.diag(L) > 0)
    npt.assert_allclose(L, sla.cholesky(A, lower=True), rtol=1e-10, atol=1e-12)
    b = RNG.standard_normal(n)
    npt.assert_allclose(impl.cholesky_solve(L, b), sla.solve(A, b, assume_a="pos"), rtol=1e-10, atol=1e-12)


def test_cholesky_refuses_unsymmetric_and_indefinite_matrices():
    with pytest.raises(ValueError):
        impl.cholesky(np.array([[2.0, 1.0], [0.0, 2.0]]))
    with pytest.raises(ValueError):
        impl.cholesky(np.array([[1.0, 2.0], [2.0, 1.0]]))         # eigenvalues 3 and -1
    with pytest.raises(ValueError):
        impl.cholesky(np.array([[1.0, 1.0], [1.0, 1.0]]))         # singular (semi-definite)
    # The pivot test *is* the definiteness test: the smallest pivot vanishes with the smallest eigenvalue.
    A = np.array([[4.0, 2.0], [2.0, 1.0 + 1e-8]])
    L = impl.cholesky(A)
    assert L[1, 1] ** 2 == pytest.approx(1e-8, rel=1e-4)


# ---------------------------------------------------------------- tridiagonal

def test_thomas_matches_dense_solve_and_costs_little():
    n = 50
    lower = RNG.uniform(-1, 1, n); upper = RNG.uniform(-1, 1, n); diag = 3.0 + RNG.uniform(0, 1, n)
    d = RNG.standard_normal(n)
    A = np.diag(diag) + np.diag(lower[1:], -1) + np.diag(upper[:-1], 1)
    npt.assert_allclose(impl.thomas_solve(lower, diag, upper, d), sla.solve(A, d), rtol=1e-12)
    with pytest.raises(ValueError):
        impl.thomas_solve(np.array([0.0, 1.0]), np.array([0.0, 1.0]), np.array([1.0, 0.0]), np.ones(2))


# ---------------------------------------------------------------- conditioning and refinement

def test_backward_error_is_tiny_and_forward_error_follows_the_condition_number():
    conds, fwd = [], []
    for n in (4, 6, 8, 10):
        H = impl.hilbert(n)
        npt.assert_allclose(H, sla.hilbert(n))
        x_true = np.ones(n)
        b = H @ x_true
        L, U, perm = impl.lu_decompose(H)
        x = impl.lu_solve(L, U, perm, b)
        kappa = impl.condition_number_1(H)
        err = impl.norm1(x - x_true) / impl.norm1(x_true)
        assert impl.backward_error(H, x, b) < 1e-15                # stable: the residual is at round-off
        assert err < kappa * np.finfo(float).eps                   # the forward-error bound holds
        conds.append(kappa); fwd.append(err)
    # The forward error grows (nearly) linearly with kappa: an observed "order" of 1 in kappa.
    slope = np.polyfit(np.log(conds), np.log(fwd), 1)[0]
    assert 0.75 < slope < 1.15, (conds, fwd)
    assert fwd[-1] > 1e-6 and fwd[0] < 1e-11


def test_iterative_refinement_reduces_the_error_on_an_ill_conditioned_system():
    H = impl.hilbert(8)
    x_true = np.ones(8)
    b = H @ x_true
    L, U, perm = impl.lu_decompose(H)
    x0 = impl.lu_solve(L, U, perm, b)
    x, hist = impl.iterative_refinement(H, L, U, perm, b, 3)
    assert hist.size == 3
    assert impl.norm1(x - x_true) <= impl.norm1(x0 - x_true)
    assert impl.norm1(H @ x - b) <= 10 * impl.norm1(H @ x0 - b) + 1e-16


# ---------------------------------------------------------------- the truss

def test_truss_stiffness_is_symmetric_with_three_rigid_body_modes_and_spd_when_supported():
    nodes, elements = impl.spar_truss()
    assert nodes.shape == (10, 2) and elements.shape == (17, 2)
    K = impl.truss_stiffness(nodes, elements, 70e9, 2e-4)
    npt.assert_allclose(K, K.T, atol=1e-6)
    lam = np.linalg.eigvalsh(K)
    assert np.sum(np.abs(lam) < 1e-6 * lam[-1]) == 3                 # two translations and a rotation
    Kf, free = impl.apply_supports(K, [0, 1, 10, 11])
    assert Kf.shape == (16, 16) and free.size == 16
    assert np.linalg.eigvalsh(Kf)[0] > 0
    L = impl.cholesky(Kf)
    npt.assert_allclose(L @ L.T, Kf, rtol=1e-12, atol=1e-12 * np.max(np.abs(Kf)))
    with pytest.raises(ValueError):
        impl.cholesky(K)
    # A single bar along x: stiffness EA/L on the x degrees of freedom only.
    K1 = impl.truss_stiffness(np.array([[0.0, 0.0], [2.0, 0.0]]), np.array([[0, 1]]), 70e9, 2e-4)
    npt.assert_allclose(K1, 7e6 * np.array([[1, 0, -1, 0], [0, 0, 0, 0], [-1, 0, 1, 0], [0, 0, 0, 0]]))


def test_truss_solution_is_in_equilibrium_and_matches_statics_and_scipy():
    E, A = 70e9, 2e-4
    nodes, elements = impl.spar_truss()
    fixed = [0, 1, 10, 11]
    f = np.zeros(20); f[9] = -4e3                                    # 4 kN down at the lower tip node
    u, R = impl.solve_truss(nodes, elements, E, A, fixed, f)
    assert np.allclose(u[fixed], 0.0)
    assert np.sum(R[1::2]) == pytest.approx(4e3, rel=1e-10)          # vertical reactions balance the load
    assert np.sum(R[0::2]) == pytest.approx(0.0, abs=1e-6)
    assert np.sum(R[[1, 11]]) == pytest.approx(4e3) and np.allclose(np.delete(R, [0, 1, 10, 11]), 0.0, atol=1e-6)
    # Moment about the root: chord forces = M / depth; shear through the diagonal.
    N = impl.member_forces(nodes, elements, E, A, u)
    assert N[1] == pytest.approx(4e3 * 2.0 / 0.3, rel=1e-9)           # upper chord at the root, tension
    assert N[0] == pytest.approx(-4e3 * 1.5 / 0.3, rel=1e-9)          # lower chord at the root, compression
    assert N[2] == pytest.approx(-4e3 * np.hypot(0.5, 0.3) / 0.3, rel=1e-9)   # root diagonal
    assert u[9] < 0 and abs(u[9]) < 0.05                              # tip goes down, by centimetres
    # Against scipy, and linear superposition of two load cases with one factorisation.
    K = impl.truss_stiffness(nodes, elements, E, A)
    Kf, free = impl.apply_supports(K, fixed)
    npt.assert_allclose(u[free], sla.solve(Kf, f[free], assume_a="pos"), rtol=1e-10)
    g = np.zeros(20); g[13] = -1e3; g[15] = -1e3
    U, _ = impl.solve_truss(nodes, elements, E, A, fixed, np.column_stack([f, g, f + g]))
    npt.assert_allclose(U[:, 2], U[:, 0] + U[:, 1], rtol=1e-10, atol=1e-15)
    npt.assert_allclose(U[:, 0], u, rtol=1e-12)


# ---------------------------------------------------------------- the heated rod: observed order

def test_heated_rod_converges_at_second_order_and_matches_the_exact_solution():
    ns = [10, 20, 40, 80]
    errs = impl.heated_rod_errors(ns)
    hs = 0.5 / (np.array(ns) + 1)
    assert np.all(np.diff(errs) < 0)
    assert abs(impl.convergence_order(hs, errs) - 2.0) < 0.1, errs
    x, lo, di, up, rhs = impl.heated_rod_system(40, 0.5, 15.0, 2e5, 300.0)
    T = impl.thomas_solve(lo, di, up, rhs)
    T_ex = impl.heated_rod_exact(x, 0.5, 15.0, 2e5, 300.0)
    assert np.max(np.abs(T - T_ex)) < 0.2
    assert T_ex[0] > 300.0 and np.max(T_ex) == pytest.approx(300.0 + 2e5 * 0.25 / (15.0 * np.pi ** 2), rel=1e-3)
    # The tridiagonal matrix is SPD (k/h^2 times [2, -1]): Cholesky accepts it.
    Ad = np.diag(di) + np.diag(lo[1:], -1) + np.diag(up[:-1], 1)
    npt.assert_allclose(impl.cholesky_solve(impl.cholesky(Ad), rhs), T, rtol=1e-10)
    with pytest.raises(ValueError):
        impl.convergence_order(hs, np.zeros(4))


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "direct_solvers.png").stat().st_size > 1000
