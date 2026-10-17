"""Public tests for P01.14 Vector Spaces, Basis & Rank. The oracle is numpy.linalg (matrix_rank, SVD); assertions are properties.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-01/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(114)


def random_rank(m: int, n: int, r: int) -> np.ndarray:
    """An (m, n) matrix of exact rank r (r <= min(m, n)), or zeros for r = 0."""
    if r == 0:
        return np.zeros((m, n))
    return RNG.standard_normal((m, r)) @ RNG.standard_normal((r, n))


def projector(basis: np.ndarray) -> np.ndarray:
    """Orthogonal projector onto the column span of ``basis`` (zero matrix for an empty basis)."""
    if basis.shape[1] == 0:
        return np.zeros((basis.shape[0], basis.shape[0]))
    return basis @ np.linalg.pinv(basis)


# ---------------------------------------------------------------- rref

def test_rref_of_a_small_matrix_by_hand():
    A = np.array([[1.0, 2.0, 2.0, 2.0], [2.0, 4.0, 6.0, 8.0], [3.0, 6.0, 8.0, 10.0]])   # Strang's example
    R, pivots = impl.rref(A)
    assert pivots == [0, 2]
    npt.assert_allclose(R, [[1.0, 2.0, 0.0, -2.0], [0.0, 0.0, 1.0, 2.0], [0.0, 0.0, 0.0, 0.0]], atol=1e-12)
    npt.assert_array_equal(A, [[1.0, 2.0, 2.0, 2.0], [2.0, 4.0, 6.0, 8.0], [3.0, 6.0, 8.0, 10.0]])   # input untouched


def test_rref_structure_identity_zero_and_idempotence():
    R, p = impl.rref(np.eye(4))
    npt.assert_array_equal(R, np.eye(4))
    assert p == [0, 1, 2, 3]
    R0, p0 = impl.rref(np.zeros((3, 5)))
    assert p0 == [] and np.all(R0 == 0.0)
    B = random_rank(6, 8, 3)
    R, piv = impl.rref(B)
    assert len(piv) == 3 and piv == sorted(piv)
    for i, c in enumerate(piv):                     # pivot columns are unit vectors e_i
        npt.assert_allclose(R[:, c], np.eye(6)[:, i], atol=1e-12)
    assert np.allclose(R[3:], 0.0)                  # rows beyond the rank are zero
    R2, piv2 = impl.rref(R)
    npt.assert_allclose(R2, R, atol=1e-12)
    assert piv2 == piv
    with pytest.raises(ValueError):
        impl.rref(np.ones(3))


# ---------------------------------------------------------------- rank

def test_rank_matches_numpy_on_matrices_of_known_rank():
    for m, n, r in ((5, 5, 5), (5, 5, 2), (8, 3, 3), (3, 8, 3), (7, 9, 0), (6, 6, 1), (10, 4, 2)):
        B = random_rank(m, n, r)
        assert impl.rank(B) == r == np.linalg.matrix_rank(B)


def test_rank_is_invariant_under_transpose_scaling_and_gram_products():
    B = random_rank(7, 5, 3)
    assert impl.rank(B.T) == impl.rank(B) == 3
    assert impl.rank(1e-6 * B) == 3 and impl.rank(1e6 * B) == 3        # relative tolerance
    assert impl.rank(B.T @ B) == 3 and impl.rank(B @ B.T) == 3
    assert impl.rank(np.hstack([B, B])) == 3                             # repeating columns adds nothing
    assert impl.rank(np.hstack([B, RNG.standard_normal((7, 1))])) == 4   # a generic new column does


def test_rank_tolerance_decides_nearly_dependent_columns():
    v = RNG.standard_normal(6)
    A = np.column_stack([v, 2 * v + 1e-13 * RNG.standard_normal(6)])
    assert impl.rank(A) == 1                     # default tolerance 1e-10: dependent
    assert impl.rank(A, tol=1e-16) == 2          # tighter than the perturbation: independent


# ---------------------------------------------------------------- null spaces and bases

def test_nullspace_is_annihilated_has_the_right_dimension_and_matches_the_svd_span():
    for m, n, r in ((5, 8, 3), (6, 6, 4), (9, 4, 2), (4, 4, 4)):
        B = random_rank(m, n, r)
        N = impl.nullspace(B)
        assert N.shape == (n, n - r)
        npt.assert_allclose(B @ N, 0.0, atol=1e-10)
        assert np.linalg.matrix_rank(N) == n - r if n > r else N.size == 0
        _, sv, vt = np.linalg.svd(B)
        N_svd = vt[r:].T
        npt.assert_allclose(projector(N), N_svd @ N_svd.T, atol=1e-10)
    # Special solutions carry a 1 in their own free variable and 0 in the others.
    B = random_rank(5, 8, 3)
    _, pivots = impl.rref(B)
    free = [c for c in range(8) if c not in pivots]
    N = impl.nullspace(B)
    npt.assert_allclose(N[free], np.eye(5), atol=1e-12)


def test_left_nullspace_annihilates_from_the_left():
    B = random_rank(7, 5, 3)
    Y = impl.left_nullspace(B)
    assert Y.shape == (7, 4)
    npt.assert_allclose(Y.T @ B, 0.0, atol=1e-10)
    assert np.linalg.matrix_rank(Y) == 4


def test_column_and_row_space_bases_span_and_are_independent():
    B = random_rank(6, 7, 3)
    C = impl.column_space_basis(B)
    assert C.shape == (6, 3) and impl.is_independent(C)
    assert np.linalg.matrix_rank(np.hstack([C, B])) == 3            # every column of B is in span(C)
    for col in C.T:                                                  # the basis is made of actual columns of B
        assert any(np.allclose(col, b) for b in B.T)
    Rw = impl.row_space_basis(B)
    assert Rw.shape == (7, 3) and impl.is_independent(Rw)
    assert np.linalg.matrix_rank(np.hstack([Rw, B.T])) == 3
    npt.assert_allclose(Rw.T @ impl.nullspace(B), 0.0, atol=1e-10)   # row space is orthogonal to the null space


def test_is_independent_on_hand_cases():
    assert impl.is_independent(np.eye(3))
    assert not impl.is_independent(np.array([[1.0, 2.0], [1.0, 2.0]]))
    assert not impl.is_independent(RNG.standard_normal((2, 3)))      # three vectors in R^2
    assert impl.is_independent(RNG.standard_normal((5, 3)))
    with pytest.raises(ValueError):
        impl.is_independent(np.ones(4))


def test_four_subspace_dimensions_obey_rank_nullity_both_ways():
    for m, n, r in ((5, 8, 3), (8, 5, 5), (6, 6, 0)):
        d = impl.four_subspace_dimensions(random_rank(m, n, r))
        assert d["rank"] == r
        assert d["rank"] + d["nullity"] == n
        assert d["rank"] + d["left_nullity"] == m


# ---------------------------------------------------------------- Ax = b

def test_solve_general_returns_a_particular_solution_plus_the_null_space_or_raises():
    B = random_rank(5, 8, 3)
    x_true = RNG.standard_normal(8)
    b = B @ x_true
    assert impl.is_consistent(B, b)
    xp, N = impl.solve_general(B, b)
    npt.assert_allclose(B @ xp, b, atol=1e-9)
    assert N.shape == (8, 5)
    for _ in range(3):
        npt.assert_allclose(B @ (xp + N @ RNG.standard_normal(5)), b, atol=1e-9)
    b_bad = b + RNG.standard_normal(5)            # generic perturbation leaves the 3-D column space
    assert not impl.is_consistent(B, b_bad)
    with pytest.raises(ValueError):
        impl.solve_general(B, b_bad)
    with pytest.raises(ValueError):
        impl.solve_general(B, b[:4])
    A = RNG.standard_normal((4, 4))
    xp, N = impl.solve_general(A, A @ x_true[:4])
    assert N.shape == (4, 0)
    npt.assert_allclose(xp, x_true[:4], atol=1e-9)


# ---------------------------------------------------------------- the truss

VARIANTS = {"determinate": (9, 0, 12, 0, 0), "missing-diagonal": (8, -1, 11, 1, 0),
            "extra-diagonal": (10, 1, 12, 0, 1), "maxwell-trap": (9, 0, 11, 1, 1)}


def test_classification_of_the_four_trusses():
    for name, (m, maxwell, rk, k, s) in VARIANTS.items():
        nodes, members, supports = impl.truss_variant(name)
        info = impl.classify_truss(nodes, members, supports)
        A = impl.truss_equilibrium_matrix(nodes, members, supports)
        assert (info["m"], info["maxwell"], info["rank"], info["k"], info["s"]) == (m, maxwell, rk, k, s), name
        assert info["rank"] == np.linalg.matrix_rank(A)
        assert info["maxwell"] == info["s"] - info["k"] == impl.maxwell_count(m, 3, 6)   # rank-nullity in disguise
    assert impl.classify_truss(*impl.truss_variant("determinate"))["status"] == "determinate"
    assert impl.classify_truss(*impl.truss_variant("missing-diagonal"))["status"] == "mechanism"
    assert impl.classify_truss(*impl.truss_variant("extra-diagonal"))["status"] == "indeterminate"
    assert impl.classify_truss(*impl.truss_variant("maxwell-trap"))["status"] == "mechanism and indeterminate"
    with pytest.raises(ValueError):
        impl.truss_variant("warren")


def test_mechanism_mode_stretches_no_member_and_moves_no_support():
    for name in ("missing-diagonal", "maxwell-trap"):
        nodes, members, supports = impl.truss_variant(name)
        A = impl.truss_equilibrium_matrix(nodes, members, supports)
        D = impl.truss_mechanisms(A)
        assert D.shape == (12, 1)
        d = D[:, 0]
        assert np.max(np.abs(d)) > 0.0
        npt.assert_allclose(A.T @ d, 0.0, atol=1e-10)
        npt.assert_allclose(impl.member_elongations(nodes, members, d), 0.0, atol=1e-10)
        for node, dof in supports:
            assert abs(d[2 * node + dof]) < 1e-12
        # The motion is in the bay without a diagonal: the loaded joint 1 moves while 0 and 2 do not.
        assert abs(d[2 * 1 + 1]) > 1e-3 * np.max(np.abs(d))
    nodes, members, supports = impl.truss_variant("determinate")
    assert impl.truss_mechanisms(impl.truss_equilibrium_matrix(nodes, members, supports)).shape == (12, 0)


def test_self_stress_state_is_in_equilibrium_without_load_and_lives_in_the_doubled_bay():
    for name in ("extra-diagonal", "maxwell-trap"):
        nodes, members, supports = impl.truss_variant(name)
        A = impl.truss_equilibrium_matrix(nodes, members, supports)
        S = impl.truss_self_stress(A)
        assert S.shape == (len(members) + 3, 1)
        x = S[:, 0] / np.max(np.abs(S[:, 0]))
        npt.assert_allclose(A @ x, 0.0, atol=1e-10)
        npt.assert_allclose(x[len(members):], 0.0, atol=1e-10)                 # no reactions are needed
        left_bay = {(0, 1), (3, 4), (0, 3), (1, 4), (0, 4), (1, 3)}
        for (i, k), f in zip(members, x):
            if (i, k) in left_bay:
                assert abs(f) > 0.1
            else:
                assert abs(f) < 1e-10
        # Both diagonals in tension, chords and verticals in compression, with the cos/sin ratio of the bay.
        diag = [f for (i, k), f in zip(members, x) if (i, k) in {(0, 4), (1, 3)}]
        chord = [f for (i, k), f in zip(members, x) if (i, k) in {(0, 1), (3, 4)}]
        assert np.sign(diag[0]) == np.sign(diag[1]) and np.sign(chord[0]) == -np.sign(diag[0])
        assert abs(chord[0] / diag[0]) == pytest.approx(0.8)


def test_member_elongations_are_linear_and_vanish_for_rigid_translation():
    nodes, members, supports = impl.truss_variant("determinate")
    d = RNG.standard_normal(12)
    e1 = impl.member_elongations(nodes, members, d)
    npt.assert_allclose(impl.member_elongations(nodes, members, 2 * d), 2 * e1)
    rigid = np.tile([0.3, -0.7], 6)
    npt.assert_allclose(impl.member_elongations(nodes, members, rigid), 0.0, atol=1e-12)
    # Consistency with the equilibrium matrix: elongation = -A_members^T d.
    A = impl.truss_equilibrium_matrix(nodes, members, supports)
    npt.assert_allclose(e1, -(A[:, :len(members)].T @ d), atol=1e-12)
