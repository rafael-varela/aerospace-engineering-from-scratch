"""Public tests for P05.16 Energy Methods. They assert properties, not magic numbers: the three energy
formulas; Castigliano's dU/dP equals the closed-form deflection; the flexibility matrix is symmetric
(Maxwell) and positive definite and its quadratic form is the energy (Clapeyron); the unit-load tip
deflection of a spar under any lift distribution equals direct double integration of M/EI; the L-bracket
reproduces its bending-plus-torsion closed form and the rigid-torsion limit; the force method reproduces
the direct stiffness method on determinate and indeterminate trusses, with joint equilibrium to
round-off; the braced spar satisfies compatibility and global equilibrium, recovers 3qL/8 in the rigid
limit, and its strut force grows monotonically with strut stiffness from zero.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

EI = 73.1e9 * 4e-6
L = 5.0


# ---------------------------------------------------------------- energy and Castigliano

def test_energy_formulas():
    assert impl.strain_energy_axial(1e3, 2.0, 1e7) == pytest.approx(1e6 * 2 / 2e7)
    assert impl.strain_energy_torsion(50.0, 0.5, 2e3) == pytest.approx(2500 * 0.5 / 4e3)
    x = np.linspace(0, L, 201)
    M = 1e3 * (L - x)                                        # tip load 1 kN
    assert impl.strain_energy_bending(x, M, EI) == pytest.approx(1e6 * L ** 3 / (6 * EI), rel=1e-12)


def test_castigliano_on_a_cantilever_gives_PL3_over_3EI():
    x = np.linspace(0, L, 401)

    def U(P):
        return impl.strain_energy_bending(x, P[0] * (L - x), EI)

    for P in (100.0, 2e3):
        d = impl.castigliano_deflection(U, np.array([P]), 0)
        assert d == pytest.approx(P * L ** 3 / (3 * EI), rel=1e-9)


def test_unit_load_tip_deflection_equals_direct_integration_for_any_lift():
    x = np.linspace(0, L, 4001)
    rng = np.random.default_rng(1601)
    for q in (impl.elliptic_lift(x, L, 18.6e3), np.full_like(x, 3.7e3), 1e3 * (1 + rng.uniform(0, 1) * x / L)):
        M = impl.cantilever_moment(x, q)
        d1 = impl.unit_load_deflection(x, M, L - x, EI)
        w = impl.deflection_by_integration(x, M, EI)
        assert d1 == pytest.approx(w[-1], rel=1e-5)          # trapezoid O(h^2) at 4001 points
        # equilibrium: M(0) is the moment of the lift about the root; M(L) = 0
        assert M[0] == pytest.approx(np.trapezoid(q * x, x), rel=1e-6)
        assert abs(M[-1]) < 1e-9 * abs(M[0])


def test_uniform_and_elliptic_closed_forms():
    x = np.linspace(0, L, 4001)
    q0 = 3.7e3
    M = impl.cantilever_moment(x, np.full_like(x, q0))
    npt.assert_allclose(M, q0 * (L - x) ** 2 / 2, rtol=1e-6, atol=1e-6 * q0 * L ** 2)
    w = impl.deflection_by_integration(x, M, EI)
    assert w[-1] == pytest.approx(q0 * L ** 4 / (8 * EI), rel=1e-6)
    q = impl.elliptic_lift(x, L, 18.6e3)
    assert np.trapezoid(q, x) == pytest.approx(18.6e3, rel=1e-3)   # sqrt singularity at the tip: slow
    Me = impl.cantilever_moment(x, q)
    assert Me[0] == pytest.approx(4 * 18.6e3 * L / (3 * np.pi), rel=1e-3)


# ---------------------------------------------------------------- flexibility matrix

def bracket():
    I = np.pi / 64 * (0.0254 ** 4 - (0.0254 - 2 * 0.00124) ** 4)
    return 0.3, 0.5, 205e9 * I, 80e9 * 2 * I


def test_l_bracket_matches_closed_form_and_castigliano():
    a, b, EIb, GJb = bracket()
    mem = impl.l_bracket_members(a, b, EIb, GJb)
    F = impl.flexibility_matrix(mem)
    P = 500.0
    exact = impl.l_bracket_tip_deflection(P, a, b, EIb, GJb)
    assert F[0, 0] * P == pytest.approx(exact, rel=1e-10)
    d = impl.castigliano_deflection(lambda Pv: impl.strain_energy_members(mem, Pv), np.array([P, 0.0]), 0)
    assert d == pytest.approx(exact, rel=1e-8)
    # torsionally rigid limit: only the two bending terms remain
    F_rigid = impl.flexibility_matrix(impl.l_bracket_members(a, b, EIb, np.inf))
    assert F_rigid[0, 0] == pytest.approx((a ** 3 + b ** 3) / (3 * EIb), rel=1e-10)


def test_flexibility_is_symmetric_positive_definite_and_gives_the_energy():
    a, b, EIb, GJb = bracket()
    mem = impl.l_bracket_members(a, b, EIb, GJb)
    F = impl.flexibility_matrix(mem)
    assert F[0, 1] == pytest.approx(F[1, 0], rel=1e-12)              # Maxwell's reciprocal theorem
    assert F[0, 1] == pytest.approx(b ** 3 / (3 * EIb), rel=1e-10)
    assert np.all(np.linalg.eigvalsh(F) > 0)
    rng = np.random.default_rng(1602)
    for _ in range(3):
        P = rng.normal(scale=500, size=2)
        U = impl.strain_energy_members(mem, P)
        assert U == pytest.approx(0.5 * P @ F @ P, rel=1e-10)          # Clapeyron: U = (1/2) sum P delta


def test_superposition_of_deflections():
    a, b, EIb, GJb = bracket()
    mem = impl.l_bracket_members(a, b, EIb, GJb)
    U = lambda Pv: impl.strain_energy_members(mem, Pv)
    d1 = impl.castigliano_deflection(U, np.array([300.0, 0.0]), 0)
    d2 = impl.castigliano_deflection(U, np.array([0.0, 200.0]), 0)
    d12 = impl.castigliano_deflection(U, np.array([300.0, 200.0]), 0)
    assert d12 == pytest.approx(d1 + d2, rel=1e-8)


# ---------------------------------------------------------------- trusses

NODES = np.array([[0.0, 0.0], [0.0, 0.30], [0.0, 0.60], [0.50, 0.25]])
FREE = np.array([6, 7])
LOAD = np.array([2.0e3, -5.9e3])


@pytest.mark.parametrize("bars", [[(0, 3), (2, 3)], [(0, 3), (1, 3), (2, 3)]])
def test_force_method_equals_stiffness_method_and_balances_the_joint(bars):
    EA = np.array([10e6, 7e6, 12e6])[:len(bars)]
    N = impl.truss_force_method(NODES, bars, EA, LOAD, FREE)
    u, Ns = impl.truss_stiffness_method(NODES, bars, EA, LOAD, FREE)
    npt.assert_allclose(N, Ns, rtol=1e-9, atol=1e-6)
    C = impl.truss_equilibrium_matrix(NODES, bars, FREE)
    npt.assert_allclose(C @ N, LOAD, atol=1e-8)                        # joint equilibrium
    for dof in (0, 1):
        d = impl.truss_unit_load_displacement(NODES, bars, EA, N, FREE, dof)
        assert d == pytest.approx(u[dof], rel=1e-9)


def test_redundant_bar_takes_load_in_proportion_to_its_stiffness():
    bars = [(0, 3), (1, 3), (2, 3)]
    N_soft = impl.truss_force_method(NODES, bars, np.array([10e6, 1e-3, 10e6]), LOAD, FREE)
    N_det = impl.truss_force_method(NODES, [(0, 3), (2, 3)], np.array([10e6, 10e6]), LOAD, FREE)
    assert abs(N_soft[1]) < 1e-6 * np.max(np.abs(N_det))               # a vanishing bar carries nothing
    npt.assert_allclose(N_soft[[0, 2]], N_det, rtol=1e-6)
    # the truss with the extra bar is stiffer
    u3, _ = impl.truss_stiffness_method(NODES, bars, np.full(3, 10e6), LOAD, FREE)
    u2, _ = impl.truss_stiffness_method(NODES, [(0, 3), (2, 3)], np.full(2, 10e6), LOAD, FREE)
    assert np.linalg.norm(u3) < np.linalg.norm(u2)


# ---------------------------------------------------------------- the braced spar

def lift(x):
    return impl.elliptic_lift(x, L, 18.6e3)


def test_braced_spar_compatibility_and_equilibrium():
    r = impl.braced_spar(L, 3.0, 1.2, EI, 73.1e9 * 1.2e-3, 73.1e9 * 3e-4, lift)
    assert abs(r["compatibility"]) < 1e-4
    x, M = r["x"], r["M"]
    # moment equilibrium about the root: lift moment minus the strut's vertical pull times a
    assert M[0] == pytest.approx(np.trapezoid(r["q"] * x, x) - r["R"] * 3.0, rel=1e-6)
    assert r["root_shear"] == pytest.approx(r["lift"] - r["R"])
    assert abs(M[-1]) < 1e-6 * abs(M).max()
    assert r["S"] > 0                                                     # lift puts the strut in tension


def test_braced_spar_rigid_limit_is_the_propped_cantilever():
    for a in (5.0, 3.0):
        r = impl.braced_spar(L, a, 1.0, EI, 1e30, 1e30, lambda x: np.full_like(x, 1e3))
        # no deflection at a rigid strut; S comes from Simpson and w from the trapezoid rule, so the two
        # agree to the trapezoid's O(h^2) error, about 1e-6 of the peak deflection at 2001 stations
        assert abs(r["w_strut"]) < 1e-5 * abs(r["w"]).max()
    r = impl.braced_spar(L, L, 1.0, EI, 1e30, 1e30, lambda x: np.full_like(x, 1e3))
    assert r["R"] == pytest.approx(3 * 1e3 * L / 8, rel=1e-6)


def test_strut_force_grows_monotonically_with_strut_stiffness():
    S = [impl.braced_spar(L, 3.0, 1.2, EI, 87.7e6, k, lift)["S"] for k in (1e-6, 1e5, 1e6, 1e7, 1e8, 1e12)]
    assert S[0] == pytest.approx(0.0, abs=1e-3)
    assert np.all(np.diff(S) > 0)
    tips = [impl.braced_spar(L, 3.0, 1.2, EI, 87.7e6, k, lift)["tip"] for k in (1e5, 1e7)]
    assert tips[1] < tips[0]


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "energy_methods.png").stat().st_size > 10_000
