"""Public tests for P05.02 Equilibrium of Rigid Bodies. They assert statics, not magic numbers: the
equilibrium matrix has the right shape and its columns move correctly with the reference point; the
classifier counts determinate, indeterminate, partially and improperly constrained support sets by rank;
solved reactions leave zero force and moment about any point; the strut-braced wing matches the
one-equation moment check about its root pin; the aileron's pushrod force matches the moment about the
hinge line; the engine mount passes a virtual-work check whose link-length changes are measured from the
moved geometry; superposition, load reversal and the near-improper blow-up hold.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)


def _beam(supports):
    """A 2-D beam from x = 0 to 4 m with the given support list."""
    return [u for s in supports for u in s]


# ---------------------------------------------------------------- matrix and loads

def test_equilibrium_matrix_columns_for_force_and_couple():
    U = [impl.Unknown("f", (1.0, 2.0, 3.0), (0.0, 0.0, 2.0)), impl.Unknown("c", (5.0, 5.0, 5.0), (1.0, 0.0, 0.0), "couple")]
    A = impl.equilibrium_matrix(U)
    assert A.shape == (6, 2)
    npt.assert_allclose(A[:, 0], [0, 0, 1, 2, -1, 0])          # [u; r x u] with u normalised
    npt.assert_allclose(A[:, 1], [0, 0, 0, 1, 0, 0])           # a couple has no force and no lever
    A2 = impl.equilibrium_matrix([impl.Unknown("f", (1.0, 2.0), (0.0, 1.0))], dim=2)
    npt.assert_allclose(A2[:, 0], [0, 1, 1])                   # rows Fx, Fy, Mz


def test_changing_the_reference_point_changes_moments_by_O_cross_R_and_not_the_solution():
    rng = np.random.default_rng(5021)
    U, _, _ = impl.hexapod_mount()
    pts, F = rng.normal(size=(3, 3)), rng.normal(size=(3, 3)) * 100
    O = rng.normal(size=3)
    b0, bO = impl.load_vector(pts, F), impl.load_vector(pts, F, ref=O)
    npt.assert_allclose(bO[:3], b0[:3])
    npt.assert_allclose(bO[3:], b0[3:] - np.cross(O, b0[:3]), atol=1e-10)
    x0 = impl.solve_reactions(U, pts, F)
    xO = impl.solve_reactions(U, pts, F, ref=O)
    for k in x0:
        assert xO[k] == pytest.approx(x0[k], rel=1e-9, abs=1e-8)


# ---------------------------------------------------------------- classification

def test_classification_of_planar_beams():
    cases = {
        "determinate": [impl.pin("A", (0, 0)), impl.roller("B", (4, 0), (0, 1))],
        "indeterminate": [impl.pin("A", (0, 0)), impl.pin("B", (4, 0))],
        "partially constrained": [impl.roller("A", (0, 0), (0, 1)), impl.roller("B", (4, 0), (0, 1))],
        "improperly constrained": [impl.roller("A", (0, 0), (0, 1)), impl.roller("B", (2, 0), (0, 1)),
                                   impl.roller("C", (4, 0), (0, 1))],
    }
    for status, sup in cases.items():
        info = impl.classify_supports(impl.equilibrium_matrix(_beam(sup), dim=2))
        assert info["status"] == status, status
        assert info["n_unknowns"] - info["rank"] == info["indeterminate_degree"]
        assert info["n_eq"] - info["rank"] == info["freedom_degree"]
    # three concurrent links: enough unknowns, rank 2: a rotation about the meeting point is free
    conc = _beam([impl.link("a", (1, 0), (0, 0)), impl.link("b", (0, 1), (0, 0)), impl.link("c", (1, 1), (0, 0))])
    info = impl.classify_supports(impl.equilibrium_matrix(conc, dim=2))
    assert info["status"] == "improperly constrained" and info["freedom_degree"] == 1
    fixed = impl.classify_supports(impl.equilibrium_matrix(impl.fixed_2d("A", (0, 0)), dim=2))
    assert fixed["status"] == "determinate"                    # a cantilever


def test_classification_of_3d_supports_and_mount_variants():
    assert impl.classify_supports(impl.equilibrium_matrix(impl.fixed_3d("A", (0, 0, 0))))["status"] == "determinate"
    a = impl.Aileron()
    assert impl.classify_supports(impl.equilibrium_matrix(a.unknowns()))["status"] == "determinate"
    two_balls = impl.ball_socket("A", (0, 0, 0)) + impl.ball_socket("B", (0, 1, 0))
    info = impl.classify_supports(impl.equilibrium_matrix(two_balls))
    assert info["status"] == "improperly constrained" and info["freedom_degree"] == 1   # spins about A-B
    assert info["indeterminate_degree"] == 1                                            # and fights itself along A-B
    expect = {(): ("determinate", 0, 0), ("extra_link",): ("indeterminate", 1, 0),
              ("drop_link",): ("partially constrained", 0, 1), ("spread0",): ("improperly constrained", 1, 1)}
    kw = {(): {}, ("extra_link",): dict(extra_link=True), ("drop_link",): dict(drop_link=2),
          ("spread0",): dict(spread_deg=0.0)}
    for key, (status, ind, free) in expect.items():
        U, _, _ = impl.hexapod_mount(**kw[key])
        info = impl.classify_supports(impl.equilibrium_matrix(U))
        assert (info["status"], info["indeterminate_degree"], info["freedom_degree"]) == (status, ind, free), key


def test_hinge_unknowns_are_normal_to_the_axis_unless_axial():
    h = impl.hinge("h", (0, 0, 0), (0, 1, 0), axial=False)
    assert len(h) == 2 and all(abs(u.direction[1]) < 1e-15 for u in h)
    h = impl.hinge("h", (0, 0, 0), (1, 1, 0), axial=True, moments=True)
    a = np.array([1, 1, 0]) / np.sqrt(2)
    assert len(h) == 5
    forces = [u for u in h if u.kind == "force"]
    assert sum(abs(np.dot(u.direction, a)) > 0.99 for u in forces) == 1


def test_solver_refuses_non_determinate_supports():
    pts, F = np.array([[2.0, 0.0]]), np.array([[0.0, -1000.0]])
    for sup in ([impl.pin("A", (0, 0)), impl.pin("B", (4, 0))],
                [impl.roller("A", (0, 0), (0, 1)), impl.roller("B", (4, 0), (0, 1))]):
        with pytest.raises(impl.SupportError):
            impl.solve_reactions(_beam(sup), pts, F, dim=2)
    assert issubclass(impl.SupportError, ValueError)


# ---------------------------------------------------------------- equilibrium of solved bodies

def test_simply_supported_beam_matches_the_lever_rule():
    U = _beam([impl.pin("A", (0, 0)), impl.roller("B", (4, 0), (0, 1))])
    P, a = 1000.0, 1.0
    r = impl.solve_reactions(U, [[a, 0.0]], [[0.0, -P]], dim=2)
    assert r["B"] == pytest.approx(P * a / 4.0)
    assert r["A_y"] == pytest.approx(P * (4.0 - a) / 4.0)
    assert r["A_x"] == pytest.approx(0.0, abs=1e-12)
    # a pure couple on the beam is reacted by an equal and opposite vertical couple
    r = impl.solve_reactions(U, [[a, 0.0]], [[0.0, 0.0]], couples=[[500.0]], dim=2)
    assert r["B"] == pytest.approx(-125.0) and r["A_y"] == pytest.approx(125.0)


def test_strut_braced_wing_matches_moments_about_the_root_and_is_in_equilibrium():
    w = impl.StrutWing(lift_half=20000.0, wing_weight=600.0, n=3.8)
    r = impl.strut_wing_reactions(w)
    assert r["strut"] == pytest.approx(impl.strut_force_closed_form(w), rel=1e-12)
    pts, F = w.loads()
    for ref in ((0.0, 0.0), (3.0, 1.0), (w.y_strut, 0.0)):
        npt.assert_allclose(impl.residual(w.unknowns(), r, pts, F, dim=2, ref=ref), 0.0, atol=1e-9 * 2e4)
    assert r["strut"] > 0.0                                     # positive lift: strut in tension
    rneg = impl.strut_wing_reactions(impl.StrutWing(lift_half=-8000.0, wing_weight=600.0, n=-1.5))
    assert rneg["strut"] < 0.0                                  # negative g: strut in compression


def test_strut_wing_superposition_and_balanced_lift():
    base = dict(semispan=5.0, y_strut=2.6)
    a = impl.strut_wing_reactions(impl.StrutWing(lift_half=10000.0, **base))
    b = impl.strut_wing_reactions(impl.StrutWing(wing_weight=500.0, n=2.0, **base))
    c = impl.strut_wing_reactions(impl.StrutWing(lift_half=10000.0, wing_weight=500.0, n=2.0, **base))
    for k in c:
        assert c[k] == pytest.approx(a[k] + b[k], abs=1e-8 * 1e4)
    # lift moment about the root balancing the weight moment: the strut idles
    W = 1000.0
    L = W * 2.5 / (4 * 5.0 / (3 * np.pi))
    r = impl.strut_wing_reactions(impl.StrutWing(lift_half=L, wing_weight=W, n=1.0, **base))
    assert r["strut"] == pytest.approx(0.0, abs=1e-9 * L)


def test_aileron_pushrod_equals_moment_about_hinge_line_and_body_is_in_equilibrium():
    rng = np.random.default_rng(5022)
    for _ in range(4):
        a = impl.Aileron(N=rng.uniform(-1500, 1500), x_cp=rng.uniform(-0.15, -0.02), y_cp=rng.uniform(0.2, 1.3),
                         Fx=rng.uniform(-100, 0))
        r = impl.aileron_reactions(a)
        assert r["pushrod"] == pytest.approx(impl.pushrod_force_by_axis_moment(a), rel=1e-10, abs=1e-9)
        assert r["pushrod"] == pytest.approx(impl.hinge_moment(a) / a.horn, rel=1e-10, abs=1e-9)
        pts, F = a.loads()
        for ref in (None, (0.3, -1.0, 2.0)):
            npt.assert_allclose(impl.residual(a.unknowns(), r, pts, F, ref=ref), 0.0, atol=1e-9 * 1500)
        assert r["hinge1_y"] == pytest.approx(0.0, abs=1e-9)    # no spanwise load: no axial reaction


def test_engine_mount_passes_the_virtual_work_check():
    U, eng, frm = impl.hexapod_mount()
    pts, F, C = impl.engine_loads(mass=120.0, nz=4.0, thrust=3000.0, torque=-450.0)
    r = impl.solve_reactions(U, pts, F, C)
    x = np.array([r[u.name] for u in U])
    scale = np.abs(x).max() * 0.5                               # force x a typical lever
    rng = np.random.default_rng(5023)
    for _ in range(6):
        t, th = rng.normal(size=3), rng.normal(size=3)
        assert abs(impl.virtual_work_residual(eng, frm, x, pts, F, C, t, th)) < 1e-8 * scale
    # a wrong answer fails it: flip one link force
    x_bad = x.copy(); x_bad[0] *= -1
    worst = max(abs(impl.virtual_work_residual(eng, frm, x_bad, pts, F, C, *rng.normal(size=(2, 3)))) for _ in range(6))
    assert worst > 1e-3 * scale


def test_link_lengths_are_unchanged_by_zero_motion_and_rotation_is_rigid():
    U, eng, frm = impl.hexapod_mount()
    l0 = np.linalg.norm(frm - eng, axis=1)
    npt.assert_allclose(impl.link_lengths(eng, frm, np.zeros(3), np.zeros(3)), l0)
    # rotating the whole engine about the thrust axis by 2 pi brings every link back
    npt.assert_allclose(impl.link_lengths(eng, frm, np.zeros(3), np.array([2 * np.pi, 0, 0])), l0, atol=1e-12)


def test_mount_forces_reverse_with_the_load_and_vanish_without_it():
    U, _, _ = impl.hexapod_mount()
    pts, F, C = impl.engine_loads(120.0, 4.0, 3000.0, -450.0)
    r1 = impl.solve_reactions(U, pts, F, C)
    r2 = impl.solve_reactions(U, pts, -F, -C)
    r0 = impl.solve_reactions(U, pts, 0 * F, 0 * C)
    for k in r1:
        assert r2[k] == pytest.approx(-r1[k])
        assert r0[k] == pytest.approx(0.0, abs=1e-12)


def test_link_forces_blow_up_as_the_mount_approaches_an_improper_constraint():
    loads = impl.engine_loads(120.0, 4.0, 3000.0, -450.0)
    spreads = np.array([30.0, 10.0, 3.0, 1.0])
    fmax, cond = impl.spread_sweep(spreads, loads)
    assert np.all(np.diff(fmax) > 0) and np.all(np.diff(cond) > 0)
    assert cond[-1] > 20 * cond[0]
    # with no torque the radial-plane mount is still loaded fine: only the torque has no path
    U0, _, _ = impl.hexapod_mount(spread_deg=0.0)
    A0 = impl.equilibrium_matrix(U0)
    b_no_torque = impl.load_vector(loads[0], loads[1])
    b_torque = impl.load_vector(loads[0], loads[1], loads[2])
    res = lambda b: np.linalg.norm(A0 @ np.linalg.lstsq(A0, -b, rcond=None)[0] + b)
    assert res(b_no_torque) < 1e-9 * np.linalg.norm(b_no_torque)
    assert res(b_torque) == pytest.approx(450.0, rel=1e-9)     # the unreachable torque about the thrust axis


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "equilibrium.png").exists()
