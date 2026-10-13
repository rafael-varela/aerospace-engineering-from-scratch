"""Public tests for P01.18 SVD & Matrix Norms. They assert properties, with numpy.linalg as the oracle.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-01/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(1818)


def bar_matrix(n: int = 10, damaged: int | None = None, ratio: float = 1.0) -> np.ndarray:
    k = impl.element_stiffness(impl.E_ALUMINIUM, impl.BAR_AREA, impl.BAR_LENGTH / n)
    ks = np.full(n, k)
    if damaged is not None:
        ks[damaged] *= ratio
    return impl.chain_stiffness_matrix(ks)


# ---------------------------------------------------------------- norms

def test_vector_norms_match_numpy_and_obey_the_norm_axioms():
    x, y = RNG.normal(size=7), RNG.normal(size=7)
    for p in (1, 2, np.inf):
        assert impl.vector_norm(x, p) == pytest.approx(np.linalg.norm(x, p))
        assert impl.vector_norm(x + y, p) <= impl.vector_norm(x, p) + impl.vector_norm(y, p) + 1e-14   # triangle
        assert impl.vector_norm(3.0 * x, p) == pytest.approx(3.0 * impl.vector_norm(x, p))           # homogeneity
        assert impl.vector_norm(np.zeros(3), p) == 0.0
    # Ordering of the three norms: |x|_inf <= |x|_2 <= |x|_1 <= sqrt(n) |x|_2.
    assert impl.vector_norm(x, np.inf) <= impl.vector_norm(x, 2) <= impl.vector_norm(x, 1) <= np.sqrt(7) * impl.vector_norm(x, 2)
    with pytest.raises(ValueError):
        impl.vector_norm(x, 3)


def test_matrix_norms_match_numpy_for_random_and_structural_matrices():
    for A in (RNG.normal(size=(5, 5)), RNG.normal(size=(8, 3)), bar_matrix(6)):
        assert impl.matrix_norm_1(A) == pytest.approx(np.linalg.norm(A, 1))
        assert impl.matrix_norm_inf(A) == pytest.approx(np.linalg.norm(A, np.inf))
        assert impl.frobenius_norm(A) == pytest.approx(np.linalg.norm(A, "fro"))
        assert impl.matrix_norm_2_power(A) == pytest.approx(np.linalg.norm(A, 2), rel=1e-9)
        # Transposing swaps the 1- and inf-norms and leaves the others alone.
        assert impl.matrix_norm_1(A.T) == pytest.approx(impl.matrix_norm_inf(A))
        assert impl.frobenius_norm(A.T) == pytest.approx(impl.frobenius_norm(A))


def test_induced_norms_are_consistent_and_submultiplicative():
    A, B = RNG.normal(size=(6, 6)), RNG.normal(size=(6, 6))
    for norm_m, p in ((impl.matrix_norm_1, 1), (impl.matrix_norm_inf, np.inf), (impl.matrix_norm_2_power, 2)):
        for _ in range(10):
            x = RNG.normal(size=6)
            assert impl.vector_norm(A @ x, p) <= norm_m(A) * impl.vector_norm(x, p) * (1 + 1e-12)   # |Ax| <= |A||x|
        assert norm_m(A @ B) <= norm_m(A) * norm_m(B) * (1 + 1e-12)                                   # |AB| <= |A||B|
    # The induced norms are attained: the 1-norm by a unit coordinate vector, the inf-norm by a sign vector.
    j = int(np.argmax(np.sum(np.abs(A), axis=0)))
    assert impl.vector_norm(A @ np.eye(6)[j], 1) == pytest.approx(impl.matrix_norm_1(A))
    i = int(np.argmax(np.sum(np.abs(A), axis=1)))
    assert impl.vector_norm(A @ np.sign(A[i]), np.inf) == pytest.approx(impl.matrix_norm_inf(A))


def test_two_norm_equals_largest_singular_value_and_frobenius_equals_root_sum_of_squares():
    A = RNG.normal(size=(7, 4))
    _, s, _ = impl.svd_jacobi(A)
    assert impl.matrix_norm_2_power(A) == pytest.approx(s[0], rel=1e-9)
    assert impl.frobenius_norm(A) == pytest.approx(np.sqrt(np.sum(s ** 2)))
    assert s[0] <= impl.frobenius_norm(A) <= np.sqrt(4) * s[0]                                      # |A|_2 <= |A|_F <= sqrt(n)|A|_2
    Q, _ = np.linalg.qr(RNG.normal(size=(5, 5)))
    assert impl.matrix_norm_2_power(Q) == pytest.approx(1.0)                                        # orthogonal: no stretch
    assert impl.condition_number(Q) == pytest.approx(1.0, rel=1e-10)


# ---------------------------------------------------------------- SVD

@pytest.mark.parametrize("m,n", [(4, 4), (9, 5), (12, 12)])
def test_svd_jacobi_reconstructs_and_matches_numpy(m, n):
    A = RNG.normal(size=(m, n))
    U, s, Vt = impl.svd_jacobi(A)
    assert U.shape == (m, n) and s.shape == (n,) and Vt.shape == (n, n)
    npt.assert_allclose((U * s) @ Vt, A, atol=1e-12)
    npt.assert_allclose(U.T @ U, np.eye(n), atol=1e-12)
    npt.assert_allclose(Vt @ Vt.T, np.eye(n), atol=1e-12)
    assert np.all(np.diff(s) <= 0.0) and np.all(s >= 0.0)
    npt.assert_allclose(s, np.linalg.svd(A, compute_uv=False), rtol=1e-10)
    # Singular vectors agree with numpy up to sign.
    Un, _, Vtn = np.linalg.svd(A, full_matrices=False)
    npt.assert_allclose(np.abs(np.sum(U * Un, axis=0)), 1.0, atol=1e-8)
    npt.assert_allclose(np.abs(np.sum(Vt * Vtn, axis=1)), 1.0, atol=1e-8)
    # A^T A has eigenvalues s^2 with eigenvectors the rows of Vt.
    npt.assert_allclose(A.T @ A @ Vt.T, Vt.T * s ** 2, atol=1e-9)


def test_svd_of_symmetric_positive_definite_matrix_is_its_eigendecomposition():
    K = bar_matrix(10)
    U, s, Vt = impl.svd_jacobi(K)
    npt.assert_allclose(s, np.sort(np.linalg.eigvalsh(K))[::-1], rtol=1e-10)
    npt.assert_allclose(U, Vt.T, atol=1e-8)                                                         # U = V for SPD
    with pytest.raises(ValueError):
        impl.svd_jacobi(np.ones((2, 3)))


def test_svd_handles_rank_deficient_and_diagonal_matrices():
    D = np.diag([3.0, 1.0, 2.0])
    U, s, Vt = impl.svd_jacobi(D)
    npt.assert_allclose(s, [3.0, 2.0, 1.0])
    npt.assert_allclose((U * s) @ Vt, D, atol=1e-14)
    A = np.outer(RNG.normal(size=6), RNG.normal(size=4)) + np.outer(RNG.normal(size=6), RNG.normal(size=4))   # rank 2
    U, s, Vt = impl.svd_jacobi(A)
    assert np.sum(s > 1e-10 * s[0]) == 2
    npt.assert_allclose(s[2:], 0.0, atol=1e-12)
    npt.assert_allclose((U * s) @ Vt, A, atol=1e-12)
    npt.assert_allclose(U.T @ U, np.eye(4), atol=1e-10)                                             # U completed to orthonormal
    assert impl.condition_number(A) == np.inf


# ---------------------------------------------------------------- condition number and the structural model

def test_chain_stiffness_assembly_and_physics():
    ks = np.array([1.0, 2.0, 3.0])
    K = impl.chain_stiffness_matrix(ks)
    npt.assert_allclose(K, [[3.0, -2.0, 0.0], [-2.0, 5.0, -3.0], [0.0, -3.0, 3.0]])
    npt.assert_allclose(K, K.T)
    assert np.all(np.linalg.eigvalsh(K) > 0.0)                                                      # clamped: no mechanism
    npt.assert_allclose(K @ np.ones(3), [1.0, 0.0, 0.0])                                             # rigid shift loads only the clamp
    # Series springs: a tip load gives tip deflection F * sum(1/k_i).
    f = np.array([0.0, 0.0, 10.0])
    assert np.linalg.solve(K, f)[-1] == pytest.approx(10.0 * np.sum(1.0 / ks))
    with pytest.raises(ValueError):
        impl.chain_stiffness_matrix(np.array([1.0, -1.0]))
    with pytest.raises(ValueError):
        impl.element_stiffness(70e9, 0.0, 0.1)


def test_uniform_bar_condition_number_matches_analytic_and_grows_as_n_squared():
    for n in (3, 10, 25):
        K = bar_matrix(n)
        c = impl.condition_number(K)
        assert c == pytest.approx(np.linalg.cond(K), rel=1e-9)
        assert c == pytest.approx(impl.uniform_chain_condition_analytic(n), rel=1e-9)
    c16, c32, c64 = (impl.condition_number(bar_matrix(n)) for n in (16, 32, 64))
    rates = np.log2([c32 / c16, c64 / c32])
    assert np.all(rates > 1.9) and np.all(rates < 2.05)                                            # second order in n
    assert impl.uniform_chain_condition_analytic(64) / (16 * 64 ** 2 / np.pi ** 2) == pytest.approx(1.0, rel=0.02)
    # Scaling every stiffness scales K but not its condition number.
    assert impl.condition_number(1e-9 * bar_matrix(10)) == pytest.approx(impl.condition_number(bar_matrix(10)), rel=1e-9)
    with pytest.raises(ValueError):
        impl.uniform_chain_condition_analytic(0)


def test_tip_deflection_against_closed_form():
    K = bar_matrix(10)
    f = np.zeros(10)
    f[-1] = 1000.0
    u = np.linalg.solve(K, f)
    assert u[-1] == pytest.approx(impl.tip_deflection_analytic(1000.0, impl.E_ALUMINIUM, impl.BAR_AREA, impl.BAR_LENGTH), rel=1e-12)
    npt.assert_allclose(np.diff(u), u[0], rtol=1e-10)                                               # uniform strain


def test_error_amplification_bound_holds_and_is_attained_in_the_worst_case():
    K = bar_matrix(10)
    f = np.zeros(10)
    f[-1] = 1000.0
    for _ in range(10):
        df = 10.0 * RNG.normal(size=10)
        actual, bound = impl.error_amplification(K, f, df)
        assert actual <= bound * (1 + 1e-12)
    f_w, df_w = impl.worst_case_loads(K)
    assert impl.vector_norm(f_w) == pytest.approx(1.0) and impl.vector_norm(df_w) == pytest.approx(1.0)
    actual, bound = impl.error_amplification(K, f_w, 0.01 * df_w)
    assert actual == pytest.approx(bound, rel=1e-8)
    assert bound == pytest.approx(0.01 * impl.condition_number(K), rel=1e-8)
    # The reverse pairing gives the smallest possible amplification, 1/cond.
    actual_r, _ = impl.error_amplification(K, df_w, 0.01 * f_w)
    assert actual_r == pytest.approx(0.01 / impl.condition_number(K), rel=1e-8)


def test_damaged_element_raises_the_condition_number_as_one_over_the_ratio():
    conds = np.array([impl.condition_number(bar_matrix(10, damaged=4, ratio=r)) for r in (1e-2, 1e-3, 1e-4)])
    assert np.all(np.diff(conds) > 0.0)
    rates = np.log10(conds[1:] / conds[:-1])
    assert np.all(rates > 0.95) and np.all(rates < 1.05)                                           # cond ~ 1/ratio
    modes = impl.mechanism_modes(bar_matrix(10, damaged=4, ratio=0.0))
    assert modes.shape == (10, 1)
    npt.assert_allclose(modes[:4, 0], 0.0, atol=1e-10)                                              # nodes before the break fixed
    npt.assert_allclose(np.abs(modes[4:, 0]), 1.0 / np.sqrt(6.0), rtol=1e-8)                         # nodes after move together
    assert impl.mechanism_modes(bar_matrix(10)).shape == (10, 0)


# ---------------------------------------------------------------- pseudoinverse and low rank

def test_pseudoinverse_matches_numpy_and_the_four_cases():
    A_tall = RNG.normal(size=(8, 3))
    A_wide = RNG.normal(size=(3, 8))
    S = bar_matrix(5)
    for A in (A_tall, A_wide, S):
        npt.assert_allclose(impl.pseudoinverse(A), np.linalg.pinv(A), atol=1e-10)
    npt.assert_allclose(impl.pseudoinverse(A_tall), np.linalg.solve(A_tall.T @ A_tall, A_tall.T), atol=1e-10)   # least squares
    npt.assert_allclose(impl.pseudoinverse(A_wide), A_wide.T @ np.linalg.inv(A_wide @ A_wide.T), atol=1e-10)   # min norm
    npt.assert_allclose(impl.pseudoinverse(S), np.linalg.inv(S), rtol=1e-8)                                   # invertible
    b = RNG.normal(size=3)
    x_min = impl.pseudoinverse(A_wide) @ b
    npt.assert_allclose(A_wide @ x_min, b, atol=1e-10)                                              # solves exactly
    z = impl.null_space(A_wide) @ RNG.normal(size=5)
    assert impl.vector_norm(x_min + z) > impl.vector_norm(x_min)                                    # any other solution is longer
    # Singular square matrix: A A^+ A = A still holds (Penrose), and the null-space part is dropped.
    K0 = bar_matrix(6, damaged=2, ratio=0.0)
    P = impl.pseudoinverse(K0)
    npt.assert_allclose(K0 @ P @ K0, K0, atol=1e-6 * impl.matrix_norm_2_power(K0))
    npt.assert_allclose(P, np.linalg.pinv(K0), rtol=1e-6, atol=1e-18)


def test_low_rank_approximation_satisfies_eckart_young():
    A = RNG.normal(size=(10, 6))
    _, s, _ = impl.svd_jacobi(A)
    for k in (1, 2, 4):
        Ak = impl.low_rank_approximation(A, k)
        assert np.linalg.matrix_rank(Ak) == k
        assert impl.matrix_norm_2_power(A - Ak) == pytest.approx(s[k], rel=1e-8)                    # 2-norm error = s_{k+1}
        assert impl.frobenius_norm(A - Ak) == pytest.approx(np.sqrt(np.sum(s[k:] ** 2)))             # Frobenius error
        # Any other rank-k matrix (here: a random projection of the same rank) does worse.
        B = impl.low_rank_approximation(A + 0.3 * RNG.normal(size=A.shape), k)
        assert impl.frobenius_norm(A - B) > impl.frobenius_norm(A - Ak)
    npt.assert_allclose(impl.low_rank_approximation(A, 6), A, atol=1e-12)
    with pytest.raises(ValueError):
        impl.low_rank_approximation(A, 0)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "stiffness_conditioning.png").stat().st_size > 1000
