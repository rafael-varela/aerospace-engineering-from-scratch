"""Public tests for P08.05 Dimensional Analysis & Similarity. They assert algebra and physics, not
magic numbers: every group found is dimensionless (D e = 0); the number of groups is n - rank D;
the SymPy (exact) and SVD (floating-point) nullspaces span the same space, also for seeded random
dimension matrices; the method of repeating variables returns the Reynolds, Mach and force-
coefficient groups and refuses dependent repeating variables; Reynolds and Mach numbers do not
change when the inputs are expressed in other units; a 1:k model at the Reynolds-matching speed
runs at k times the full-scale Mach number; and the pressurised and cryogenic tunnel solutions,
substituted back, match both numbers.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

WING = ["F", "rho", "V", "c", "mu", "a"]
PIPE = ["dp", "rho", "V", "d", "mu", "eps", "L"]


# ---------------------------------------------------------------- the dimension matrix

def test_dimension_matrix_columns():
    D = impl.dimension_matrix(WING)
    assert D.shape == (4, 6)
    npt.assert_array_equal(D[:, 0], [1, 1, -2, 0])                  # force: M L T^-2
    npt.assert_array_equal(D[:, 4], [1, -1, -1, 0])                 # viscosity: M L^-1 T^-1


@pytest.mark.parametrize("names,count", [(WING, 3), (PIPE, 4), (["P", "rho", "omega", "d", "V"], 2),
                                         (["V", "g", "L"], 1), (["k", "mu", "cp"], 1)])
def test_pi_count_is_n_minus_rank(names, count):
    D = impl.dimension_matrix(names)
    assert impl.pi_count(D) == count == len(names) - np.linalg.matrix_rank(D)


def test_sympy_groups_are_integer_and_dimensionless():
    for names in (WING, PIPE):
        D = impl.dimension_matrix(names)
        groups = impl.pi_groups_sympy(D)
        assert len(groups) == impl.pi_count(D)
        for e in groups:
            assert all(float(x).is_integer() for x in e)
            assert impl.is_dimensionless(D, np.array(e, dtype=float))


def test_svd_basis_is_orthonormal_and_dimensionless():
    D = impl.dimension_matrix(PIPE)
    B = impl.pi_groups_svd(D)
    assert B.shape == (7, 4)
    npt.assert_allclose(B.T @ B, np.eye(4), atol=1e-12)
    assert np.max(np.abs(D @ B)) < 1e-12


@pytest.mark.parametrize("names", [WING, PIPE, ["P", "rho", "omega", "d", "V", "mu"]])
def test_sympy_and_svd_span_the_same_space(names):
    D = impl.dimension_matrix(names)
    B_sym = np.array(impl.pi_groups_sympy(D), dtype=float).T
    assert impl.same_span(B_sym, impl.pi_groups_svd(D))


def test_random_dimension_matrices_agree():
    rng = np.random.default_rng(5)                                  # seeded integer dimension matrices
    for _ in range(10):
        n = int(rng.integers(4, 8))
        D = rng.integers(-3, 4, size=(4, n)).astype(float)
        D[3] = D[0] + D[1] if rng.random() < 0.5 else D[3]          # sometimes rank deficient
        B_svd = impl.pi_groups_svd(D)
        B_sym = np.array(impl.pi_groups_sympy(D), dtype=float).T
        assert B_svd.shape[1] == len(impl.pi_groups_sympy(D)) == impl.pi_count(D)
        if B_svd.shape[1]:
            assert impl.same_span(B_sym, B_svd)


def test_same_span_detects_different_spaces():
    a = np.array([[1.0], [0.0], [0.0]])
    b = np.array([[0.0], [1.0], [0.0]])
    assert not impl.same_span(a, b)
    assert impl.same_span(a, 3.0 * a)


# ---------------------------------------------------------------- repeating variables

def test_repeating_variables_give_force_coefficient_reynolds_and_mach():
    D = impl.dimension_matrix(WING)
    g = impl.pi_groups_repeating(D, WING, ["rho", "V", "c"])
    assert set(g) == {"F", "mu", "a"}
    npt.assert_allclose(g["F"], [1, -1, -2, -2, 0, 0], atol=1e-12)       # F / (rho V^2 c^2)
    npt.assert_allclose(g["mu"], [0, -1, -1, -1, 1, 0], atol=1e-12)      # 1 / Re
    npt.assert_allclose(g["a"], [0, 0, -1, 0, 0, 1], atol=1e-12)         # 1 / M
    B_rep = np.array(list(g.values())).T
    assert impl.same_span(B_rep, impl.pi_groups_svd(D))


def test_repeating_variables_must_be_independent():
    D = impl.dimension_matrix(WING)
    with pytest.raises(ValueError):
        impl.pi_groups_repeating(D, WING, ["V", "a", "c"])            # V and a have the same dimensions


def test_is_dimensionless():
    D = impl.dimension_matrix(WING)
    assert impl.is_dimensionless(D, np.array([0, 1, 1, 1, -1, 0.0]))   # Re
    assert not impl.is_dimensionless(D, np.array([0, 1, 1, 0, 0, 0.0]))


def test_format_group():
    assert impl.format_group(["rho", "V", "c", "mu"], np.array([1, 1, 1, -1.0])) == "rho^1 V^1 c^1 mu^-1"


# ---------------------------------------------------------------- similarity numbers

def test_reynolds_and_mach_are_unit_invariant():
    rho, V, c, mu, a = 0.9091, 60.0, 1.5, 1.6937e-5, 328.6
    Re_si = impl.reynolds(rho, V, c, mu)
    # the same flow in g, cm, s: rho in g/cm^3, V in cm/s, c in cm, mu in poise (g/(cm s))
    Re_cgs = impl.reynolds(rho * 1e-3, V * 100.0, c * 100.0, mu * 10.0)
    assert Re_cgs == pytest.approx(Re_si, rel=1e-12)
    # knots and feet for speed and length, slugs per cubic foot and lbf s/ft^2 for rho and mu
    kt, ft, slug_ft3, lbfs_ft2 = 1852.0 / 3600.0, 0.3048, 515.3788, 47.88026
    Re_imp = impl.reynolds(rho / slug_ft3, V / ft, c / ft, mu / lbfs_ft2)
    assert Re_imp == pytest.approx(Re_si, rel=1e-6)
    assert impl.mach(V / kt, a / kt) == pytest.approx(impl.mach(V, a), rel=1e-14)


def test_full_scale_numbers_use_the_shared_atmosphere():
    fs = impl.full_scale_numbers(3000.0, 60.0, 1.5)
    assert fs["rho"] == pytest.approx(0.9091, rel=1e-4)
    assert fs["mu"] == pytest.approx(1.6937e-5, rel=1e-4)
    assert fs["Re"] == pytest.approx(fs["rho"] * 60.0 * 1.5 / fs["mu"], rel=1e-14)
    assert 0.18 < fs["M"] < 0.19


def test_reynolds_matching_speed_scales_inversely_with_model_size_so_mach_cannot_match():
    fs = impl.full_scale_numbers(3000.0, 60.0, 1.5)
    T, p, rho = impl.isa_troposphere(3000.0)
    for k in (2.0, 5.0):
        V = impl.speed_for_reynolds(fs["Re"], 1.5 / k, rho, fs["mu"])
        assert impl.reynolds(rho, V, 1.5 / k, fs["mu"]) == pytest.approx(fs["Re"], rel=1e-12)
        # same air as full scale: the Re-matched model runs at k times the Mach number
        assert impl.mach(V, impl.speed_of_sound(T)) == pytest.approx(k * fs["M"], rel=1e-12)


def test_pressurised_tunnel_matches_both_and_scales_with_model_size():
    fs = impl.full_scale_numbers(3000.0, 60.0, 1.5)
    out = [impl.pressurised_tunnel(fs["Re"], fs["M"], 1.5 / k, 288.15) for k in (2.0, 4.0)]
    for k, o in zip((2.0, 4.0), out):
        assert impl.reynolds(o["rho"], o["V"], 1.5 / k, o["mu"]) == pytest.approx(fs["Re"], rel=1e-12)
        assert impl.mach(o["V"], impl.speed_of_sound(288.15)) == pytest.approx(fs["M"], rel=1e-12)
    assert out[1]["p"] / out[0]["p"] == pytest.approx(2.0, rel=1e-12)


def test_cryogenic_tunnel_matches_both_at_one_atmosphere():
    fs = impl.full_scale_numbers(3000.0, 60.0, 1.5)
    o = impl.cryogenic_tunnel(fs["Re"], fs["M"], 0.3, 101325.0)
    assert o["T"] < 288.15
    assert impl.reynolds(o["rho"], o["V"], 0.3, o["mu"]) == pytest.approx(fs["Re"], rel=1e-9)
    assert impl.mach(o["V"], impl.speed_of_sound(o["T"])) == pytest.approx(fs["M"], rel=1e-12)
    with pytest.raises(ValueError):
        impl.cryogenic_tunnel(fs["Re"] * 1e4, fs["M"], 0.003, 101325.0)


def test_stokes_drag_collapses_to_cd_re_equals_24():
    # any fluid, any small sphere: dimensional analysis says C_D depends on Re alone
    for rho, mu, V, d in ((998.2, 1.002e-3, 1e-4, 1e-4), (1.204, 1.81e-5, 1e-3, 1e-5), (1260.0, 1.41, 0.01, 0.003)):
        cd = impl.drag_coefficient(impl.stokes_drag(mu, V, d), rho, V, np.pi * d ** 2 / 4)
        assert cd * impl.reynolds(rho, V, d, mu) == pytest.approx(24.0, rel=1e-12)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "similarity.png").stat().st_size > 10_000
