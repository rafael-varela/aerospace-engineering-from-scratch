"""Public tests for P03.09 Iterative & Sparse Solvers. They assert matrix identities, convergence to
the direct solution, the theoretical contraction rates, the CG error bound, preconditioner
behaviour, the scaling of iteration counts with resolution, and a discretisation order, with
``scipy.sparse.linalg`` as the oracle.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-03/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
import scipy.sparse as sp
from scipy import linalg as sla
from scipy.sparse.linalg import cg as scipy_cg, spsolve, spsolve_triangular

from conftest import load_impl

impl = load_impl(__file__)

L = 0.2


def manufactured_system(n, k=1.0):
    X, Y, h = impl.grid(n, L)
    T_ex, f = impl.manufactured_solution(X, Y, L, k)
    A = impl.poisson_matrix(n, h, k)
    b = impl.poisson_rhs(n, h, f, 0.0, k)
    return A, b, T_ex, h


# ---------------------------------------------------------------- the matrix

def test_grid_spacing_and_layout():
    X, Y, h = impl.grid(4, L)
    assert h == pytest.approx(L / 5)
    assert X.shape == (4, 4)
    npt.assert_allclose(X[0], h * np.arange(1, 5))       # X[j, i] = x_i
    npt.assert_allclose(Y[:, 0], h * np.arange(1, 5))


def test_poisson_matrix_is_the_kronecker_laplacian():
    n, h, k = 6, 0.1, 3.0
    A = impl.poisson_matrix(n, h, k)
    assert sp.issparse(A)
    T = sp.diags([-np.ones(n - 1), 2 * np.ones(n), -np.ones(n - 1)], [-1, 0, 1])
    ref = (k / h ** 2) * (sp.kron(sp.eye(n), T) + sp.kron(T, sp.eye(n)))
    assert abs(A - ref).max() < 1e-12
    assert A.nnz == 5 * n * n - 4 * n
    assert abs(A - A.T).max() == 0.0


def test_variable_conductivity_matrix_is_symmetric_and_uses_harmonic_means():
    n, h = 5, 0.05
    kf = np.full((n + 2, n + 2), 2.0)
    kf[:, 4:] = 8.0                                        # a jump in k across a vertical line
    A = impl.poisson_matrix(n, h, kf)
    assert abs(A - A.T).max() < 1e-12
    # unknown (i=2, j=2) sits at grid column 3 (k = 2) and its right neighbour at column 4 (k = 8)
    p = 2 + n * 2
    k_face = 2 * 2.0 * 8.0 / (2.0 + 8.0)                   # harmonic mean 3.2
    assert A[p, p + 1] == pytest.approx(-k_face / h ** 2)
    assert A[p, p - 1] == pytest.approx(-2.0 / h ** 2)
    assert np.all(sla.eigvalsh(A.toarray()) > 0)
    uniform = impl.poisson_matrix(n, h, np.full((n + 2, n + 2), 2.0))
    assert abs(uniform - impl.poisson_matrix(n, h, 2.0)).max() < 1e-12


def test_wall_temperature_with_no_source_is_the_solution():
    n, h = 7, 0.02
    A = impl.poisson_matrix(n, h, 5.0)
    b = impl.poisson_rhs(n, h, np.zeros((n, n)), 350.0, 5.0)
    T = spsolve(sp.csc_matrix(A), b)
    npt.assert_allclose(T, 350.0, rtol=1e-12)


def test_eigenvalue_formula_and_jacobi_radius():
    n, h = 6, 0.03
    A = impl.poisson_matrix(n, h, 1.0).toarray()
    ev = sla.eigvalsh(A)
    lo, hi = impl.poisson_eigen_bounds(n, h)
    assert lo == pytest.approx(ev[0], rel=1e-10)
    assert hi == pytest.approx(ev[-1], rel=1e-10)
    G = np.eye(n * n) - A / np.diag(A)[:, None]
    rho = np.max(np.abs(np.linalg.eigvals(G)))
    assert impl.jacobi_spectral_radius(n) == pytest.approx(rho, rel=1e-10)
    assert impl.jacobi_spectral_radius(n) == pytest.approx(np.cos(np.pi / 7))


# ---------------------------------------------------------------- stationary methods

def test_jacobi_converges_at_the_predicted_rate():
    A, b, _, _ = manufactured_system(10)
    x_ref = spsolve(sp.csc_matrix(A), b)
    x, hist = impl.jacobi(A, b, tol=1e-8)
    assert hist[-1] <= 1e-8 * np.linalg.norm(b)
    assert np.max(np.abs(x - x_ref)) < 1e-6
    ratio = np.mean(hist[-40:] / hist[-41:-1])
    assert ratio == pytest.approx(impl.jacobi_spectral_radius(10), abs=2e-4)
    predicted = np.log(1e-8 * np.linalg.norm(b) / hist[0]) / np.log(impl.jacobi_spectral_radius(10))
    assert hist.size - 1 == pytest.approx(predicted, rel=0.05)


def test_gauss_seidel_takes_half_the_jacobi_iterations():
    A, b, _, _ = manufactured_system(8)
    x_ref = spsolve(sp.csc_matrix(A), b)
    x_j, h_j = impl.jacobi(A, b, tol=1e-7)
    x_gs, h_gs = impl.gauss_seidel(A, b, tol=1e-7)
    assert np.max(np.abs(x_gs - x_ref)) < 1e-6
    assert h_gs[-1] <= 1e-7 * np.linalg.norm(b)
    assert 0.42 < (h_gs.size - 1) / (h_j.size - 1) < 0.58


def test_initial_guess_is_honoured():
    A, b, _, _ = manufactured_system(6)
    x_ref = spsolve(sp.csc_matrix(A), b)
    for method in (impl.jacobi, impl.gauss_seidel, impl.steepest_descent, impl.conjugate_gradient):
        x, hist = method(A, b, x0=x_ref, tol=1e-10)
        assert hist.size == 1                              # already converged: no iteration taken
        npt.assert_allclose(x, x_ref)


# ---------------------------------------------------------------- Krylov methods

def test_steepest_descent_converges_and_residual_decreases():
    A, b, _, _ = manufactured_system(8)
    x_ref = spsolve(sp.csc_matrix(A), b)
    x, hist = impl.steepest_descent(A, b, tol=1e-8)
    assert np.max(np.abs(x - x_ref)) < 1e-6
    assert hist[-1] <= 1e-8 * np.linalg.norm(b)
    # SD minimises phi(x) = (1/2) ||x - x*||_A^2 along each step, so the A-norm error must fall monotonically
    # (the residual norm itself need not: it zigzags)
    errs = [impl.a_norm(A, impl.steepest_descent(A, b, tol=0.0, maxiter=k)[0] - x_ref) for k in range(0, 20)]
    assert np.all(np.diff(errs) < 0)


def test_conjugate_gradient_matches_direct_solve_and_scipy():
    A, b, _, _ = manufactured_system(12)
    x_ref = spsolve(sp.csc_matrix(A), b)
    x, hist = impl.conjugate_gradient(A, b, tol=1e-10)
    assert np.max(np.abs(x - x_ref)) < 1e-9 * np.max(np.abs(x_ref))
    assert hist[0] == pytest.approx(np.linalg.norm(b))
    x_sp, info = scipy_cg(A, b, rtol=1e-10, atol=0.0)
    assert info == 0
    npt.assert_allclose(x, x_sp, rtol=1e-8, atol=1e-12)
    assert hist.size - 1 < 40                              # far fewer than the N = 144 finite-termination bound


def test_conjugate_gradient_terminates_in_at_most_n_steps():
    rng = np.random.default_rng(309)
    B = rng.standard_normal((9, 9))
    A = sp.csr_matrix(B @ B.T + 9 * np.eye(9))
    b = rng.standard_normal(9)
    x, hist = impl.conjugate_gradient(A, b, tol=1e-14, maxiter=9)
    assert hist[-1] / hist[0] < 1e-10
    npt.assert_allclose(x, np.linalg.solve(A.toarray(), b), rtol=1e-8)


def test_cg_a_norm_error_decreases_monotonically_and_obeys_the_bound():
    A, b, _, _ = manufactured_system(12)
    x_ref = spsolve(sp.csc_matrix(A), b)
    lo, hi = impl.poisson_eigen_bounds(12, L / 13)
    kappa = hi / lo
    e0 = impl.a_norm(A, x_ref)
    errs = []
    for k in range(0, 10):
        xk, _ = impl.conjugate_gradient(A, b, tol=0.0, maxiter=k)
        errs.append(impl.a_norm(A, xk - x_ref) / e0)
    errs = np.array(errs)
    assert np.all(np.diff(errs) < 0)
    assert np.all(errs[1:] <= impl.cg_error_bound(kappa, np.arange(1, 10)) * (1 + 1e-12))
    assert impl.cg_error_bound(1.0, 5) == 0.0
    assert impl.cg_error_bound(100.0, 0) == pytest.approx(2.0)


def test_steepest_descent_is_much_slower_than_cg():
    A, b, _, _ = manufactured_system(12)
    _, h_sd = impl.steepest_descent(A, b, tol=1e-6)
    _, h_cg = impl.conjugate_gradient(A, b, tol=1e-6)
    assert (h_sd.size - 1) > 10 * (h_cg.size - 1)


# ---------------------------------------------------------------- preconditioners

def test_jacobi_preconditioner_is_the_identity_up_to_scaling_for_constant_k():
    A, b, _, _ = manufactured_system(8)
    M = impl.jacobi_preconditioner(A)
    npt.assert_allclose(M(b), b / A.diagonal())
    _, h0 = impl.conjugate_gradient(A, b, tol=1e-9)
    _, h1 = impl.conjugate_gradient(A, b, tol=1e-9, M=M)
    assert h0.size == h1.size
    npt.assert_allclose(h0[:-1], h1[:-1], rtol=1e-8)       # identical iterates (the last entry is round-off)
    assert h1[-1] <= 1e-9 * np.linalg.norm(b)


def test_incomplete_cholesky_is_exact_when_there_is_no_fill():
    n = 9
    T = sp.csr_matrix(sp.diags([-np.ones(n - 1), 4 * np.ones(n), -np.ones(n - 1)], [-1, 0, 1]))
    Lic = impl.incomplete_cholesky(T)
    assert sp.issparse(Lic)
    L_exact = sla.cholesky(T.toarray(), lower=True)
    npt.assert_allclose(Lic.toarray(), L_exact, atol=1e-13)
    A = impl.poisson_matrix(6, 0.1, 1.0)
    Lic = impl.incomplete_cholesky(A)
    assert abs(sp.triu(Lic, 1)).max() == 0.0               # lower triangular
    pattern_A = set(zip(*sp.tril(A).nonzero()))
    assert set(zip(*Lic.nonzero())) <= pattern_A            # no fill-in kept
    with pytest.raises(ValueError):
        impl.incomplete_cholesky(sp.csr_matrix(np.array([[1.0, 2.0], [2.0, 1.0]])))


def test_sparse_triangular_solves_match_scipy():
    rng = np.random.default_rng(310)
    n = 12
    Ld = np.tril(rng.standard_normal((n, n))) * (rng.random((n, n)) < 0.4) + 3 * np.eye(n)
    Ls = sp.csr_matrix(Ld)
    b = rng.standard_normal(n)
    npt.assert_allclose(impl.sparse_triangular_solve(Ls, b, lower=True), spsolve_triangular(Ls, b, lower=True), rtol=1e-12)
    Us = sp.csr_matrix(Ld.T)
    npt.assert_allclose(impl.sparse_triangular_solve(Us, b, lower=False), spsolve_triangular(Us, b, lower=False), rtol=1e-12)


def test_ic_preconditioner_cuts_the_iteration_count_and_keeps_the_answer():
    kfield = impl.plate_conductivity(12, L, 1.0, 0.01, (0.06, 0.14, 0.06, 0.14))
    h = L / 13
    A = impl.poisson_matrix(12, h, kfield)
    b = impl.poisson_rhs(12, h, np.full((12, 12), 2000.0), 300.0, kfield)
    x_ref = spsolve(sp.csc_matrix(A), b)
    x0, h0 = impl.conjugate_gradient(A, b, tol=1e-9)
    x1, h1 = impl.conjugate_gradient(A, b, tol=1e-9, M=impl.jacobi_preconditioner(A))
    x2, h2 = impl.conjugate_gradient(A, b, tol=1e-9, M=impl.ic_preconditioner(A))
    for x in (x0, x1, x2):
        assert np.max(np.abs(x - x_ref)) < 1e-5 * np.max(np.abs(x_ref - 300.0))
    assert h1.size < 0.8 * h0.size                           # a 100:1 jump in k: diagonal scaling pays
    assert h2.size < 0.85 * h1.size                          # IC(0) pays more
    assert h2.size < 0.6 * h0.size


# ---------------------------------------------------------------- the plate and the studies

def test_plate_centre_rise_matches_the_series_solution_and_is_symmetric():
    kt, q = 170.0 * 2e-3, 2000.0
    X, Y, T, its = impl.solve_plate(32, L, kt, q, 300.0)
    assert T.shape == (32, 32)
    rise = T.max() - 300.0
    assert rise == pytest.approx(0.07367 * q * L ** 2 / kt, rel=3e-3)
    npt.assert_allclose(T, T.T, rtol=1e-8)                 # symmetric about the diagonal
    npt.assert_allclose(T, T[::-1, :], rtol=1e-8)          # and about the mid-plane
    assert 40 < its < 120


def test_plate_conductivity_field_marks_the_insert():
    k = impl.plate_conductivity(8, L, 1.0, 0.1, (0.05, 0.1, 0.0, 0.2))
    assert k.shape == (10, 10)
    assert k[0, 0] == 1.0
    x = L / 9 * np.arange(10)
    inside = (x >= 0.05) & (x <= 0.1)
    npt.assert_array_equal(k[:, inside], 0.1)
    npt.assert_array_equal(k[:, ~inside], 1.0)


def test_discretisation_order_is_two():
    ns = [6, 12, 24]
    errs = impl.discretisation_errors(ns, L=L, k=1.0)
    hs = L / (np.array(ns) + 1)
    assert np.all(np.diff(errs) < 0)
    assert 1.85 < impl.convergence_order(hs, errs) < 2.15
    with pytest.raises(ValueError):
        impl.convergence_order(hs, np.array([1.0, 0.0, -1.0]))


def test_cg_iterations_grow_linearly_and_jacobi_quadratically_with_n():
    ns = [8, 16, 32]
    its_cg = impl.iterations_vs_n(impl.conjugate_gradient, ns, tol=1e-6, L=L)
    its_j = impl.iterations_vs_n(impl.jacobi, ns, tol=1e-6, L=L)
    assert 0.8 < impl.scaling_exponent(ns, its_cg) < 1.25
    assert 1.7 < impl.scaling_exponent(ns, its_j) < 2.2
    assert its_j[-1] > 50 * its_cg[-1]


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "iterative_solvers.png").exists()
