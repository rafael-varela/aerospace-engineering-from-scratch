"""Public tests for P03.13 Finite Element Basics. They assert properties, with SciPy as the oracle.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-03/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import quad
from scipy.linalg import eigh, solve_banded

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(313)
E, RHO, L, A = 70.0e9, 2700.0, 5.0, 6.0e-4
OM = 400.0 * 2 * np.pi / 60


# ---------------------------------------------------------------- exact solutions

def test_exact_solutions_satisfy_their_differential_equations_and_boundary_conditions():
    x = np.linspace(0.0, L, 2001)
    h = x[1] - x[0]
    P = 17.0e3
    u, s = impl.exact_rotating_bar(x, E, RHO, OM, L, A, P)
    assert u[0] == 0.0 and s[-1] == pytest.approx(P / A, rel=1e-12)              # fixed root, tip load
    npt.assert_allclose(E * np.gradient(u, h)[5:-5], s[5:-5], rtol=1e-5)          # sigma = E u'
    npt.assert_allclose(A * np.gradient(s, h)[5:-5], -RHO * A * OM ** 2 * x[5:-5], rtol=1e-4)   # (A sigma)' = -b
    u0, s0 = impl.exact_rotating_bar(x, E, RHO, OM, L, A, 0.0)
    assert s0[0] == pytest.approx(0.5 * RHO * OM ** 2 * L ** 2)                   # root stress = total CF load / A
    x2 = np.linspace(0.0, 1.5, 2001)
    h2 = x2[1] - x2[0]
    u2, s2 = impl.exact_tapered_bar(x2, E, 8e-4, 0.6, 1.5, 5e4)
    assert u2[0] == 0.0
    npt.assert_allclose(s2 * 8e-4 * (1 - 0.6 * x2 / 1.5), np.full_like(x2, 5e4), rtol=1e-12)   # constant axial force
    npt.assert_allclose(E * np.gradient(u2, h2)[5:-5], s2[5:-5], rtol=1e-4)
    with pytest.raises(ValueError):
        impl.exact_tapered_bar(x, E, 8e-4, 1.0, 1.5, 5e4)


# ---------------------------------------------------------------- element level

def test_shape_functions_partition_unity_and_element_matrices_have_the_textbook_form():
    xi = np.linspace(-1, 1, 7)
    N = impl.shape_functions(xi)
    assert N.shape == (7, 2)
    npt.assert_allclose(N.sum(axis=1), 1.0)
    npt.assert_allclose(N[0], [1.0, 0.0])
    npt.assert_allclose(N[-1], [0.0, 1.0])
    k = impl.element_stiffness(E, A, 0.5)
    npt.assert_allclose(k, E * A / 0.5 * np.array([[1, -1], [-1, 1]]))
    npt.assert_allclose(k @ np.ones(2), 0.0, atol=1e-6)                         # rigid translation: no force
    m = impl.element_mass(RHO, A, 0.5)
    assert m.sum() == pytest.approx(RHO * A * 0.5)                               # total mass of the element
    npt.assert_allclose(m, m.T)
    with pytest.raises(ValueError):
        impl.element_stiffness(E, A, 0.0)
    with pytest.raises(ValueError):
        impl.element_load(lambda x: x, 1.0, 1.0)


def test_element_load_matches_scipy_quad_and_the_lumped_limit():
    x1, x2 = 1.2, 1.7
    b = lambda x: 3.0 * x ** 2 - x + 2.0                                        # noqa: E731
    f = impl.element_load(b, x1, x2, n_gauss=3)
    N1 = lambda x: (x2 - x) / (x2 - x1)                                         # noqa: E731
    N2 = lambda x: (x - x1) / (x2 - x1)                                         # noqa: E731
    assert f[0] == pytest.approx(quad(lambda x: N1(x) * b(x), x1, x2)[0], rel=1e-12)
    assert f[1] == pytest.approx(quad(lambda x: N2(x) * b(x), x1, x2)[0], rel=1e-12)
    assert f.sum() == pytest.approx(quad(b, x1, x2)[0], rel=1e-12)             # total load preserved
    npt.assert_allclose(impl.element_load(lambda x: np.full_like(x, 4.0), 0.0, 2.0), [4.0, 4.0])   # uniform load splits evenly
    npt.assert_allclose(impl.element_load(lambda x: x, 0.0, 1.0), [1 / 6, 1 / 3])                  # linear load: 1/3 - 2/3 split


# ---------------------------------------------------------------- assembly and solve

def test_assembled_stiffness_is_symmetric_singular_before_and_spd_after_constraint():
    nodes = np.sort(RNG.uniform(0.0, L, 9))
    nodes[0], nodes[-1] = 0.0, L
    K, f, M = impl.assemble(nodes, E, A, lambda x: np.ones_like(x), RHO)
    npt.assert_allclose(K, K.T) and npt.assert_allclose(M, M.T)
    npt.assert_allclose(K @ np.ones(nodes.size), 0.0, atol=1e-3)                # one rigid-body mode
    assert np.linalg.matrix_rank(K) == nodes.size - 1
    assert np.all(np.linalg.eigvalsh(K[1:, 1:]) > 0.0)                          # SPD once the root is fixed
    assert np.all(np.linalg.eigvalsh(M) > 0.0)
    assert f.sum() == pytest.approx(L)                                          # unit body force: total L
    assert M.sum() == pytest.approx(RHO * A * L)
    # Tridiagonal: nothing couples non-adjacent nodes.
    assert np.max(np.abs(np.triu(K, 2))) == 0.0
    with pytest.raises(ValueError):
        impl.assemble(np.array([0.0, 1.0, 0.5]), E, A)
    with pytest.raises(ValueError):
        impl.solve_bar(K, f, {})


def test_solution_equals_scipy_banded_solve_and_reaction_balances_the_load():
    nodes = np.linspace(0.0, L, 13)
    b = lambda x: RHO * A * OM ** 2 * x                                          # noqa: E731
    K, f, _ = impl.assemble(nodes, E, A, b)
    P = 17.0e3
    f = impl.apply_point_load(f, 12, P)
    u, r = impl.solve_bar(K, f, {0: 0.0})
    assert u[0] == 0.0 and np.all(u[1:] > 0.0) and np.all(np.diff(u) > 0.0)
    npt.assert_allclose(r[1:], 0.0, atol=1e-6)                                  # reactions only at the fixed dof
    assert r[0] == pytest.approx(-(0.5 * RHO * A * OM ** 2 * L ** 2 + P), rel=1e-10)   # equilibrium
    Kf = K[1:, 1:]
    ab = np.zeros((3, Kf.shape[0]))
    ab[0, 1:] = np.diag(Kf, 1)
    ab[1] = np.diag(Kf)
    ab[2, :-1] = np.diag(Kf, -1)
    npt.assert_allclose(u[1:], solve_banded((1, 1), ab, f[1:]), rtol=1e-10)
    # A prescribed non-zero end displacement with no load gives a uniform strain.
    u2, r2 = impl.solve_bar(K, np.zeros(13), {0: 0.0, 12: 1e-3})
    npt.assert_allclose(impl.element_stress(u2, nodes, E), E * 1e-3 / L, rtol=1e-10)
    assert r2[12] == pytest.approx(E * A * 1e-3 / L, rel=1e-10)


def test_patch_test_linear_field_is_reproduced_exactly_on_a_nonuniform_mesh():
    nodes = np.sort(np.concatenate([[0.0, L], RNG.uniform(0.0, L, 10)]))
    K, f, _ = impl.assemble(nodes, E, A)
    eps = 2.5e-4
    u, r = impl.solve_bar(K, impl.apply_point_load(f, nodes.size - 1, E * A * eps), {0: 0.0})
    npt.assert_allclose(u, eps * nodes, rtol=1e-12, atol=1e-18)                 # exact for a linear field
    npt.assert_allclose(impl.element_stress(u, nodes, E), E * eps, rtol=1e-12)
    npt.assert_allclose(impl.interpolate(u, nodes, [0.37 * L, 0.9 * L]), eps * np.array([0.37 * L, 0.9 * L]), rtol=1e-12)


def test_uniform_bar_is_exact_at_the_nodes_and_minimises_potential_energy():
    c = impl.rotor_blade()
    P = c["m_tip"] * c["omega"] ** 2 * c["L"]
    for n in (3, 10, 25):
        nodes, u, sig, K, f = impl.solve_rotor_blade(c, n)
        u_ex, s_ex = impl.exact_rotating_bar(nodes, c["E"], c["rho"], c["omega"], c["L"], c["A"], P)
        npt.assert_allclose(u, u_ex, rtol=1e-10, atol=1e-18)                    # nodal exactness
        # The element stress is the exact stress averaged over the element (nodal values are exact), which
        # differs from the exact stress at the element centre by h^2 sigma'' / 24.
        npt.assert_allclose(sig, c["E"] * np.diff(u_ex) / np.diff(nodes), rtol=1e-10)
        mid = 0.5 * (nodes[:-1] + nodes[1:])
        s_mid = impl.exact_rotating_bar(mid, c["E"], c["rho"], c["omega"], c["L"], c["A"], P)[1]
        h = c["L"] / n
        npt.assert_allclose(sig - s_mid, -c["rho"] * c["omega"] ** 2 * h ** 2 / 24, rtol=1e-8)
        # Galerkin: the FE solution minimises Pi over admissible fields; any perturbation raises it.
        pi0 = impl.potential_energy(K, f, u)
        for _ in range(5):
            d = RNG.standard_normal(u.size) * 1e-4 * np.max(u)
            d[0] = 0.0
            assert impl.potential_energy(K, f, u + d) > pi0
        # Clapeyron: at equilibrium, strain energy = half the external work, Pi = -U.
        assert impl.strain_energy(K, u) == pytest.approx(0.5 * f @ u, rel=1e-12)
        assert pi0 == pytest.approx(-impl.strain_energy(K, u), rel=1e-12)


def test_tapered_bar_matches_scipy_bvp_and_converges():
    A0, c, L2, P = 8.0e-4, 0.6, 1.5, 5.0e4
    errs = []
    for n in (4, 8, 16, 32):
        nodes = np.linspace(0.0, L2, n + 1)
        K, f, _ = impl.assemble(nodes, E, lambda x: A0 * (1 - c * x / L2))
        u, _ = impl.solve_bar(K, impl.apply_point_load(f, n, P), {0: 0.0})
        u_ex = impl.exact_tapered_bar(nodes, E, A0, c, L2, P)[0]
        npt.assert_allclose(impl.tapered_bar_bvp_oracle(E, A0, c, L2, P, nodes), u_ex, rtol=1e-7)
        errs.append(np.max(np.abs(u - u_ex)))
    npt.assert_allclose(impl.observed_order(np.array(errs)), 2.0, atol=0.1)      # mid-point area: O(h^2) at the nodes
    assert errs[-1] < 1e-4 * u_ex[-1]


def test_observed_orders_are_two_in_l2_and_one_in_energy():
    e_l2, e_en = impl.convergence_study(impl.rotor_blade(), [4, 8, 16, 32])
    npt.assert_allclose(impl.observed_order(e_l2), 2.0, atol=0.05)
    npt.assert_allclose(impl.observed_order(e_en), 1.0, atol=0.05)
    assert np.all(np.diff(e_l2) < 0) and np.all(np.diff(e_en) < 0)
    # Error norms are zero for an exactly representable field.
    nodes = np.linspace(0.0, 1.0, 6)
    assert impl.l2_error(2.0 * nodes, nodes, lambda x: 2.0 * x) == pytest.approx(0.0, abs=1e-14)
    assert impl.energy_error(2.0 * nodes, nodes, E, lambda x: np.full_like(x, 2.0 * E)) == pytest.approx(0.0, abs=1e-6)


# ---------------------------------------------------------------- vibration

def test_axial_frequencies_match_scipy_eigh_and_converge_to_the_exact_ones_from_above():
    w_ex = impl.exact_axial_frequencies(E, RHO, L, 3)
    assert w_ex[0] == pytest.approx(np.pi / (2 * L) * np.sqrt(E / RHO))
    prev = None
    for n in (4, 8, 16, 32):
        nodes = np.linspace(0.0, L, n + 1)
        K, _, M = impl.assemble(nodes, E, A, None, RHO)
        w = impl.axial_frequencies(K, M, [0])
        lam = eigh(K[1:, 1:], M[1:, 1:], eigvals_only=True)
        npt.assert_allclose(w, np.sqrt(lam), rtol=1e-9)
        assert np.all(w[:3] >= w_ex - 1e-9)                                     # consistent mass: upper bounds
        if prev is not None:
            assert np.all(w[:3] <= prev[:3] + 1e-9)                              # monotone with refinement
        prev = w
    err = []
    for n in (8, 16, 32):
        nodes = np.linspace(0.0, L, n + 1)
        K, _, M = impl.assemble(nodes, E, A, None, RHO)
        err.append(impl.axial_frequencies(K, M, [0])[0] - w_ex[0])
    npt.assert_allclose(impl.observed_order(np.array(err)), 2.0, atol=0.1)       # first mode converges as h^2


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "bar_fem.png").stat().st_size > 1000
