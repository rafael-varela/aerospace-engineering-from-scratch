"""Public tests for P05.03 Trusses: Joints & Sections. They assert statics, not magic numbers: the joint
matrix has one column per member and reaction and its member columns sum to zero force; counting and rank
classify determinate, indeterminate and mechanism trusses; the method of joints leaves zero force at every
joint and its reactions balance the loads; the method of sections, with reactions from whole-body
equilibrium, reproduces the joint forces in every panel; so does the beam analogy M/h and V/sin; the
direct stiffness method gives the same forces whatever EA is for a determinate truss and EA-dependent
forces for an indeterminate one; symmetry, superposition, zero-force members, a vanishing load and the
1/depth scaling of chord forces hold; the 3-D engine mount is determinate and in equilibrium.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)


def _triangle(P=1000.0):
    """Three-member truss: pin at 0, roller at 1, load at the apex 2."""
    nodes = [[0.0, 0.0], [4.0, 0.0], [1.0, 3.0]]
    members = [[0, 1], [0, 2], [1, 2]]
    fixed = [[True, True], [False, True], [False, False]]
    loads = np.zeros((3, 2)); loads[2] = [0.0, -P]
    return impl.Truss(nodes, members, fixed, loads)


# ---------------------------------------------------------------- the joint matrix and counting

def test_joint_matrix_shape_and_member_columns_are_self_equilibrated():
    t = impl.pratt_truss()
    B = impl.joint_matrix(t)
    m, r = len(t.members), int(t.fixed.sum())
    assert B.shape == (t.nodes.size, m + r)
    d = t.dim
    for k in range(m):                            # a member pulls its two joints equally and oppositely
        col = B[:, k].reshape(-1, d)
        npt.assert_allclose(col.sum(axis=0), 0.0, atol=1e-15)
        assert np.count_nonzero(np.abs(col).sum(axis=1)) == 2


def test_determinacy_classification():
    assert impl.determinacy(impl.pratt_truss())["status"] == "determinate"
    ic = impl.determinacy(impl.pratt_truss(n_panels=4, cross_diagonals=True))
    assert ic["status"] == "indeterminate" and ic["indeterminate_degree"] == 4
    # remove a diagonal: m + r = 2j - 1, a mechanism (the panel shears)
    t = impl.pratt_truss(n_panels=4)
    g = impl.pratt_members(4)
    keep = [k for k in range(len(t.members)) if k != g["diagonals"][1]]
    mech = impl.Truss(t.nodes, t.members[keep], t.fixed, t.loads)
    im = impl.determinacy(mech)
    assert im["status"] == "mechanism" and im["mechanism_degree"] == 1
    with pytest.raises(impl.TrussError):
        impl.solve_joints(mech)
    with pytest.raises(impl.TrussError):
        impl.stiffness_solve(mech, 1e6)
    # right count, wrong geometry: three collinear joints, the middle one free, is a mechanism
    col = impl.Truss([[0, 0], [1, 0], [2, 0]], [[0, 1], [1, 2]], [[True, True], [False, False], [True, True]])
    assert impl.determinacy(col)["status"] == "mechanism"


# ---------------------------------------------------------------- the method of joints

def test_triangle_matches_hand_solution():
    P = 1000.0
    t = _triangle(P)
    f, R = impl.solve_joints(t)
    # reactions by moments about joint 0: R1 * 4 = P * 1
    assert R[1, 1] == pytest.approx(P / 4) and R[0, 1] == pytest.approx(3 * P / 4) and R[0, 0] == pytest.approx(0.0, abs=1e-12)
    # joint 1: vertical equilibrium of member 1-2 gives its force, then the base member
    u12 = np.array([-3.0, 3.0]) / np.sqrt(18)
    f12 = -(P / 4) / u12[1]
    assert f[2] == pytest.approx(f12)
    assert f[0] == pytest.approx(f12 * u12[0])
    assert f[1] < 0 and f[2] < 0 and f[0] > 0                 # rafters in compression, tie in tension


def test_joint_equilibrium_and_global_balance():
    for t in (impl.pratt_truss(), impl.pratt_truss(8, 6.0, 0.9, 1500.0), impl.engine_mount_truss()):
        f, R = impl.solve_joints(t)
        scale = np.abs(t.loads).max()
        npt.assert_allclose(impl.joint_residual(t, f, R), 0.0, atol=1e-10 * scale)
        npt.assert_allclose(t.loads.sum(axis=0) + R.sum(axis=0), 0.0, atol=1e-9 * scale)
        assert np.all(R[~t.fixed] == 0.0)


def test_symmetric_load_gives_symmetric_forces_and_superposition_holds():
    n = 6
    t = impl.pratt_truss(n)
    f, _ = impl.solve_joints(t)
    g = impl.pratt_members(n)
    for name in ("bottom", "top", "diagonals", "verticals"):
        npt.assert_allclose(f[g[name]], f[g[name]][::-1], atol=1e-9)
    t2 = impl.pratt_truss(n, P=2 * 2000.0)
    npt.assert_allclose(impl.solve_joints(t2)[0], 2 * f, atol=1e-9)
    t0 = impl.pratt_truss(n, P=0.0)
    npt.assert_allclose(impl.solve_joints(t0)[0], 0.0, atol=1e-12)
    # signs of a Pratt under gravity: bottom chord tension, top chord compression, diagonals tension
    assert np.all(f[g["top"]] < 0) and np.all(f[g["diagonals"]] > 0) and np.all(f[g["bottom"]] >= -1e-9)


def test_zero_force_members_at_an_unloaded_two_member_joint():
    # joint 3 joins two collinear members and an unloaded third: the third carries nothing
    nodes = [[0, 0], [2, 0], [1, 1], [1, 0]]
    members = [[0, 3], [3, 1], [3, 2], [0, 2], [1, 2]]
    fixed = [[True, True], [False, True], [False, False], [False, False]]
    loads = np.zeros((4, 2)); loads[2] = [0.0, -500.0]
    t = impl.Truss(nodes, members, fixed, loads)
    f, _ = impl.solve_joints(t)
    assert list(impl.zero_force_members(f, 500.0)) == [2]
    zf = impl.zero_force_members(impl.solve_joints(impl.pratt_truss(6))[0], 2000.0)
    assert len(zf) == 3


# ---------------------------------------------------------------- the method of sections

def test_external_reactions_from_whole_body_equal_joint_reactions():
    t = impl.pratt_truss(8, 6.4, 0.8, 1200.0)
    t.loads[3] += [400.0, -900.0]                             # an unsymmetric extra load
    _, R = impl.solve_joints(t)
    npt.assert_allclose(impl.external_reactions_2d(t), R, atol=1e-9)


def test_method_of_sections_reproduces_joint_forces_in_every_panel():
    for n, span, h, P in ((6, 4.8, 0.6, 2000.0), (8, 6.4, 0.5, 1000.0)):
        t = impl.pratt_truss(n, span, h, P)
        f, _ = impl.solve_joints(t)
        Rext = impl.external_reactions_2d(t)
        for panel in range(n):
            cut, side = impl.pratt_section(n, panel)
            npt.assert_allclose(impl.solve_section(t, cut, side, Rext), f[cut], atol=1e-9 * P)
            # the free body on the other side gives the same forces
            other = [j for j in range(len(t.nodes)) if j not in side]
            npt.assert_allclose(impl.solve_section(t, cut, other, Rext), f[cut], atol=1e-9 * P)


def test_beam_analogy_chord_forces_and_one_over_depth_scaling():
    n, span, P = 6, 4.8, 2000.0
    for h in (0.4, 0.6, 1.2):
        t = impl.pratt_truss(n, span, h, P)
        f, _ = impl.solve_joints(t)
        for panel in range(n):
            cut, _ = impl.pratt_section(n, panel)
            npt.assert_allclose(impl.pratt_beam_analogy(t, n, panel), f[cut], atol=1e-9 * P)
    depths = np.array([0.3, 0.6, 1.2])
    fmax, _ = impl.depth_sweep(span, n, P, depths)
    npt.assert_allclose(fmax * depths, fmax[0] * depths[0], rtol=1e-12)   # chord force = M_max / h


def test_section_refuses_a_member_that_does_not_cross_the_cut():
    t = impl.pratt_truss(6)
    cut, side = impl.pratt_section(6, 2)
    with pytest.raises(impl.TrussError):
        impl.solve_section(t, [cut[0], 0], side, impl.external_reactions_2d(t))  # member 0 lies inside


# ---------------------------------------------------------------- the stiffness method

def test_stiffness_method_matches_joints_whatever_EA_for_a_determinate_truss():
    rng = np.random.default_rng(5031)
    for t in (impl.pratt_truss(), impl.engine_mount_truss()):
        f, R = impl.solve_joints(t)
        scale = np.abs(f).max()
        for EA in (1e6, rng.uniform(1e5, 1e8, len(t.members))):
            _, fk, Rk = impl.stiffness_solve(t, EA)
            npt.assert_allclose(fk, f, atol=1e-8 * scale)
            npt.assert_allclose(Rk, R, atol=1e-8 * scale)


def test_stiffness_method_on_indeterminate_truss_balances_and_depends_on_EA():
    t = impl.pratt_truss(4, cross_diagonals=True)
    A = impl.tube_area(0.01905, 0.000889)
    _, f1, R1 = impl.stiffness_solve(t, 200e9 * A)
    npt.assert_allclose(impl.joint_residual(t, f1, R1), 0.0, atol=1e-8 * 2000)
    EA = np.full(len(t.members), 200e9 * A); EA[-1] *= 10
    _, f2, R2 = impl.stiffness_solve(t, EA)
    npt.assert_allclose(impl.joint_residual(t, f2, R2), 0.0, atol=1e-8 * 2000)
    assert abs(f2[-1] - f1[-1]) > 100.0                      # the stiffer diagonal attracts load
    npt.assert_allclose(R1, R2, atol=1e-8 * 2000)            # but the external reactions are determinate
    # uniform scaling of every EA scales displacements and leaves forces alone
    Q1, f3, _ = impl.stiffness_solve(t, 2 * 200e9 * A)
    Q0, _, _ = impl.stiffness_solve(t, 200e9 * A)
    npt.assert_allclose(f3, f1, atol=1e-8 * 2000)
    npt.assert_allclose(Q1, Q0 / 2, atol=1e-15)


# ---------------------------------------------------------------- 3-D engine mount and helpers

def test_engine_mount_is_determinate_and_loads_are_statically_equivalent():
    t = impl.engine_mount_truss()
    info = impl.determinacy(t)
    assert info["status"] == "determinate" and info["m"] == 12 and info["r"] == 12
    ring = t.nodes[4:]
    cg = np.array([0.55 + 0.30, 0.0, 0.0])
    F, M = np.array([2000.0, 0.0, -4000.0]), np.array([-300.0, 50.0, 0.0])
    f = impl.rigid_body_nodal_loads(ring, cg, F, M)
    npt.assert_allclose(f.sum(axis=0), F, atol=1e-9)
    npt.assert_allclose(np.cross(ring - cg, f).sum(axis=0), M, atol=1e-9)


def test_tube_formulas():
    od, wall = 0.02, 0.001
    assert impl.tube_area(od, wall) == pytest.approx(np.pi * (od - wall) * wall)
    assert impl.tube_area(od, od / 2) == pytest.approx(np.pi * od ** 2 / 4)          # solid bar limit
    assert impl.tube_second_moment(od, od / 2) == pytest.approx(np.pi * od ** 4 / 64)
    assert impl.euler_load_pinned(200e9, 1e-9, 2.0) == pytest.approx(impl.euler_load_pinned(200e9, 1e-9, 1.0) / 4)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "truss.png").exists()
