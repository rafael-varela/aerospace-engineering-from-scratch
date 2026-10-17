"""Public tests for P01.17 Orthogonality, Least Squares & QR. They assert properties, with numpy.linalg as the oracle.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-01/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(1717)


def random_tall(m: int, n: int) -> np.ndarray:
    return RNG.normal(size=(m, n))


def same_up_to_column_signs(Q1: np.ndarray, Q2: np.ndarray, atol: float = 1e-10) -> bool:
    signs = np.sign(np.sum(Q1 * Q2, axis=0))
    return np.allclose(Q1 * signs, Q2, atol=atol)


# ---------------------------------------------------------------- Gram-Schmidt QR

@pytest.mark.parametrize("factor", ["classical_gram_schmidt", "modified_gram_schmidt"])
@pytest.mark.parametrize("m,n", [(5, 2), (17, 3), (40, 6)])
def test_gram_schmidt_gives_a_valid_reduced_qr(factor, m, n):
    A = random_tall(m, n)
    Q, R = getattr(impl, factor)(A)
    assert Q.shape == (m, n) and R.shape == (n, n)
    npt.assert_allclose(Q.T @ Q, np.eye(n), atol=1e-12)                       # orthonormal columns
    npt.assert_allclose(Q @ R, A, atol=1e-12)                                 # exact factorisation
    npt.assert_allclose(R, np.triu(R), atol=0.0)                              # upper triangular
    assert np.all(np.diag(R) > 0.0)                                           # positive diagonal: unique QR
    # numpy's Householder QR may differ by column signs only.
    Qn, Rn = np.linalg.qr(A)
    assert same_up_to_column_signs(Q, Qn)
    npt.assert_allclose(np.abs(R), np.abs(Rn), atol=1e-10)


def test_gram_schmidt_first_column_and_a_simple_hand_case():
    A = np.array([[1.0, 1.0], [1.0, 0.0], [0.0, 1.0]])
    Q, R = impl.modified_gram_schmidt(A)
    npt.assert_allclose(Q[:, 0], A[:, 0] / np.sqrt(2.0))                      # q_0 is a_0 normalised
    npt.assert_allclose(R[0, 0], np.sqrt(2.0))
    npt.assert_allclose(R[0, 1], 1.0 / np.sqrt(2.0))                          # projection of a_1 on q_0
    npt.assert_allclose(Q[:, 1], np.array([1.0, -1.0, 2.0]) / np.sqrt(6.0))
    npt.assert_allclose(R[1, 1], np.sqrt(1.5))


def test_gram_schmidt_rejects_dependent_columns_and_wide_matrices():
    A = np.column_stack([np.ones(5), 2.0 * np.ones(5)])
    for fn in (impl.classical_gram_schmidt, impl.modified_gram_schmidt):
        with pytest.raises(ValueError):
            fn(A)
        with pytest.raises(ValueError):
            fn(np.ones((2, 3)))


def test_modified_beats_classical_on_an_ill_conditioned_matrix():
    A, _ = impl.laeuchli_matrix(1e-8)
    Qc, _ = impl.classical_gram_schmidt(A)
    Qm, _ = impl.modified_gram_schmidt(A)
    loss_c, loss_m = impl.orthogonality_loss(Qc), impl.orthogonality_loss(Qm)
    assert loss_c > 1e-2                                                      # classical: Q is not orthogonal at all
    assert loss_m < 1e-6 and loss_m < 1e-4 * loss_c                            # modified: orthogonal to ~eps*cond
    assert impl.orthogonality_loss(np.linalg.qr(A)[0]) < 1e-14                 # Householder: orthogonal to rounding
    # Both still reproduce A = QR: loss of orthogonality is not loss of the factorisation.
    for fn in (impl.classical_gram_schmidt, impl.modified_gram_schmidt):
        Q, R = fn(A)
        npt.assert_allclose(Q @ R, A, atol=1e-15)


def test_back_substitution_matches_solve_and_rejects_singular():
    for n in (2, 5, 9):
        R = np.triu(random_tall(n, n)) + 3.0 * np.eye(n)
        y = RNG.normal(size=n)
        npt.assert_allclose(impl.back_substitution(R, y), np.linalg.solve(R, y), rtol=1e-12)
    with pytest.raises(ValueError):
        impl.back_substitution(np.array([[1.0, 2.0], [0.0, 0.0]]), np.ones(2))


# ---------------------------------------------------------------- projections

def test_projection_matrix_is_an_orthogonal_projector_onto_the_column_space():
    A = random_tall(12, 3)
    P = impl.projection_matrix(A)
    npt.assert_allclose(P @ P, P, atol=1e-12)                                 # idempotent
    npt.assert_allclose(P.T, P, atol=1e-12)                                   # symmetric
    npt.assert_allclose(P @ A, A, atol=1e-12)                                 # leaves the columns alone
    assert np.trace(P) == pytest.approx(3.0)                                  # rank = number of columns
    npt.assert_allclose(P, A @ np.linalg.solve(A.T @ A, A.T), atol=1e-11)     # the textbook formula
    b = RNG.normal(size=12)
    e = b - P @ b
    npt.assert_allclose(A.T @ e, 0.0, atol=1e-11)                             # error is orthogonal to every column
    assert np.linalg.norm(P @ b) ** 2 + np.linalg.norm(e) ** 2 == pytest.approx(np.linalg.norm(b) ** 2)   # Pythagoras


# ---------------------------------------------------------------- least squares

def test_least_squares_solvers_agree_with_lstsq_and_satisfy_the_normal_equations():
    A = random_tall(30, 4)
    b = RNG.normal(size=30)
    ref = np.linalg.lstsq(A, b, rcond=None)[0]
    x_n = impl.least_squares_normal(A, b)
    x_q = impl.least_squares_qr(A, b)
    npt.assert_allclose(x_n, ref, rtol=1e-10)
    npt.assert_allclose(x_q, ref, rtol=1e-10)
    assert impl.residual_orthogonality(A, b, x_q) < 1e-11
    assert impl.residual_orthogonality(A, b, x_q + 1e-3) > 1e-4               # a wrong x leaves a non-orthogonal residual
    # The solution minimises the residual: any perturbation makes it larger.
    r0 = np.linalg.norm(b - A @ x_q)
    for _ in range(5):
        assert np.linalg.norm(b - A @ (x_q + 1e-2 * RNG.normal(size=4))) > r0


def test_consistent_system_is_solved_exactly_and_square_system_reduces_to_solve():
    A = random_tall(20, 3)
    x_true = np.array([0.3, -1.2, 2.0])
    npt.assert_allclose(impl.least_squares_qr(A, A @ x_true), x_true, rtol=1e-12)
    S = random_tall(4, 4) + 4.0 * np.eye(4)
    b = RNG.normal(size=4)
    npt.assert_allclose(impl.least_squares_qr(S, b), np.linalg.solve(S, b), rtol=1e-10)


def test_normal_equations_square_the_condition_number_and_fail_first():
    for eps in (1e-3, 1e-5):
        r = impl.laeuchli_comparison(eps)
        assert r["cond_AtA"] == pytest.approx(r["cond_A"] ** 2, rel=1e-6)
        assert r["err_qr"] < 1e-13 and r["err_lstsq"] < 1e-13
    # Below eps = 1e-8, eps**2 is under machine epsilon: A^T A rounds to the singular all-ones matrix.
    broken = impl.laeuchli_comparison(1e-8)
    assert broken["err_normal"] > 1e-2 or np.isinf(broken["err_normal"])
    assert broken["err_qr"] < 1e-13 and broken["err_lstsq"] < 1e-13
    # Forming Q^T b explicitly inherits the loss of orthogonality of Q, divided by the tiny R diagonal.
    assert impl.laeuchli_comparison(1e-7)["err_qr_explicit"] > 1e-4
    c_a, c_ata = impl.condition_squaring(impl.design_matrix(np.linspace(0.0, 1.5, 17)))
    assert c_ata == pytest.approx(c_a ** 2, rel=1e-8)
    with pytest.raises(ValueError):
        impl.laeuchli_matrix(0.0)


def test_explicit_and_augmented_qr_agree_for_a_well_conditioned_problem():
    A = random_tall(25, 3)
    b = RNG.normal(size=25)
    npt.assert_allclose(impl.least_squares_qr(A, b, explicit_qtb=True), impl.least_squares_qr(A, b), rtol=1e-12)


# ---------------------------------------------------------------- the drag polar

def test_design_matrix_columns_and_model_names():
    cl = np.array([-0.2, 0.0, 0.5, 1.0])
    A = impl.design_matrix(cl)
    npt.assert_allclose(A, np.column_stack([np.ones(4), cl ** 2]))
    npt.assert_allclose(impl.design_matrix(cl, "general"), np.column_stack([np.ones(4), cl, cl ** 2]))
    with pytest.raises(ValueError):
        impl.design_matrix(cl, "cubic")


def test_fit_recovers_the_exact_polar_without_noise():
    _, cl, cd = impl.wind_tunnel_data(noise_cl=0.0, noise_cd=0.0)
    fit = impl.fit_drag_polar(cl, cd)
    assert fit["cd0"] == pytest.approx(impl.TRUE_CD0, rel=1e-12)
    assert fit["k"] == pytest.approx(impl.TRUE_K, rel=1e-12)
    assert fit["rms"] < 1e-15
    assert impl.oswald_efficiency(fit["k"], impl.ASPECT_RATIO) == pytest.approx(impl.TRUE_E, rel=1e-12)


def test_fit_with_noise_recovers_the_polar_within_its_error_bars():
    for seed in (1, 2, 3, 17):
        _, cl, cd = impl.wind_tunnel_data(seed=seed)
        fit = impl.fit_drag_polar(cl, cd)
        assert abs(fit["cd0"] - impl.TRUE_CD0) < 4.0 * fit["cd0_se"]
        assert abs(fit["k"] - impl.TRUE_K) < 4.0 * fit["k_se"]
        assert 0.0 < fit["cd0_se"] < 0.1 * fit["cd0"] and 0.0 < fit["k_se"] < 0.1 * fit["k"]
        assert fit["residuals"].shape == cl.shape
        assert np.abs(np.mean(fit["residuals"])) < 1e-12                       # a constant column forces zero-mean residuals
        assert abs(np.sum(fit["residuals"] * cl ** 2)) < 1e-12                 # and orthogonality to the C_L^2 column
    with pytest.raises(ValueError):
        impl.fit_drag_polar(np.array([0.1, 0.2]), np.array([0.02, 0.03]))


def test_standard_errors_scale_like_one_over_sqrt_n_and_match_monte_carlo():
    # With the same residual level, four times the points over the same alpha range halves the standard error
    # (the small extra factor is the m/(m-2) degrees-of-freedom correction).
    _, cl, _ = impl.wind_tunnel_data(n_points=17, noise_cl=0.0)
    _, cl4, _ = impl.wind_tunnel_data(n_points=68, noise_cl=0.0)
    se17 = impl.parameter_standard_errors(impl.design_matrix(cl), 1e-3 * np.ones(17))
    se68 = impl.parameter_standard_errors(impl.design_matrix(cl4), 1e-3 * np.ones(68))
    expected = 2.0 * np.sqrt((17.0 / 15.0) / (68.0 / 66.0))
    npt.assert_allclose(se17 / se68, expected, rtol=0.05)
    sd, se = impl.monte_carlo_k_scatter(n_trials=150, seed=5, noise_cl=0.0)
    assert sd == pytest.approx(se, rel=0.2)                                   # formula holds when the regressor is exact
    sd_cl, se_cl = impl.monte_carlo_k_scatter(n_trials=150, seed=5, noise_cl=0.01)
    assert sd_cl > sd                                                         # regressor noise adds scatter
    with pytest.raises(ValueError):
        impl.parameter_standard_errors(np.ones((2, 2)), np.zeros(2))


def test_max_lift_to_drag_formula_and_oswald_validation():
    ld, cl_opt = impl.max_lift_to_drag(0.025, 0.05)
    assert cl_opt == pytest.approx(np.sqrt(0.5))
    assert ld == pytest.approx(1.0 / (2.0 * np.sqrt(0.025 * 0.05)))
    # Brute force: L/D = cl / (cd0 + k cl^2) peaks at the same place.
    cls = np.linspace(0.05, 2.0, 20001)
    ld_curve = cls / (0.025 + 0.05 * cls ** 2)
    assert ld_curve.max() == pytest.approx(ld, rel=1e-6)
    assert cls[np.argmax(ld_curve)] == pytest.approx(cl_opt, abs=2e-4)
    for bad in ((0.0, 0.05), (0.025, -1.0)):
        with pytest.raises(ValueError):
            impl.max_lift_to_drag(*bad)
    with pytest.raises(ValueError):
        impl.oswald_efficiency(-0.01, 8.0)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "drag_polar_fit.png").stat().st_size > 1000
