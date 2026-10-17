"""Public tests for P05.01 Forces, Moments & Free-Body Diagrams. They assert statics, not magic numbers:
a moment is perpendicular to the lever and the force and vanishes for a force through the point;
Varignon's theorem; the moment about an axis does not depend on the point chosen on the axis; moving a
force system changes its moment by (A - B) x R and a couple's moment not at all; the central axis
carries a moment parallel to the resultant; distributed loads reproduce the uniform, triangular and
elliptic closed forms; the free-body solver leaves zero residual about any point and refuses concurrent
unknowns; the landing-gear brace force from the 3 x 3 solve equals the one-equation moment check about
the trunnion; superposition and the D = 0 limit hold.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)


# ---------------------------------------------------------------- moments

def test_moment_is_perpendicular_to_lever_and_force_with_magnitude_r_f_sin():
    rng = np.random.default_rng(5011)
    for _ in range(10):
        r, F, P = rng.normal(size=3), rng.normal(size=3) * 1e3, rng.normal(size=3)
        M = impl.moment_about_point(r, F, P)
        d = r - P
        assert abs(np.dot(M, d)) < 1e-9 * np.linalg.norm(M) * np.linalg.norm(d) + 1e-12
        assert abs(np.dot(M, F)) < 1e-9 * np.linalg.norm(M) * np.linalg.norm(F) + 1e-12
        sin = np.linalg.norm(np.cross(d, F)) / (np.linalg.norm(d) * np.linalg.norm(F))
        assert np.linalg.norm(M) == pytest.approx(np.linalg.norm(d) * np.linalg.norm(F) * sin)
    # a force whose line passes through P has no moment about P
    P = np.array([1.0, 2.0, 3.0])
    F = np.array([2.0, -1.0, 4.0])
    npt.assert_allclose(impl.moment_about_point(P + 2.5 * F, F, P), 0.0, atol=1e-12)


def test_moment_2d_matches_z_component_of_3d_and_sign_convention():
    rng = np.random.default_rng(5012)
    r, F, P = rng.normal(size=2), rng.normal(size=2), rng.normal(size=2)
    M3 = impl.moment_about_point(np.append(r, 0), np.append(F, 0), np.append(P, 0))
    assert impl.moment_2d(r, F, P) == pytest.approx(M3[2])
    assert M3[0] == 0.0 and M3[1] == 0.0
    assert impl.moment_2d([1.0, 0.0], [0.0, 10.0]) == pytest.approx(10.0)       # anticlockwise positive


def test_varignon_moment_of_the_sum_equals_sum_of_moments():
    rng = np.random.default_rng(5013)
    r, P = rng.normal(size=3), rng.normal(size=3)
    Fs = rng.normal(size=(5, 3))
    lhs = impl.moment_about_point(r, Fs.sum(axis=0), P)
    rhs = sum(impl.moment_about_point(r, f, P) for f in Fs)
    npt.assert_allclose(lhs, rhs, atol=1e-12)


def test_moment_about_axis_is_independent_of_the_point_on_the_axis_and_ignores_parallel_forces():
    rng = np.random.default_rng(5014)
    r, F, P, u = rng.normal(size=3), rng.normal(size=3), rng.normal(size=3), rng.normal(size=3)
    m1 = impl.moment_about_axis(r, F, P, u)
    m2 = impl.moment_about_axis(r, F, P + 3.7 * u, 5.0 * u)                    # other point, unscaled u
    assert m1 == pytest.approx(m2)
    assert m1 == pytest.approx(np.dot(u / np.linalg.norm(u), np.cross(r - P, F)))
    assert impl.moment_about_axis(r, 4.0 * u, P, u) == pytest.approx(0.0, abs=1e-12)   # F parallel to axis
    assert impl.moment_about_axis(P + 2.0 * u, F, P, u) == pytest.approx(0.0, abs=1e-12)  # F cuts the axis


def test_transport_moment_and_couples():
    rng = np.random.default_rng(5015)
    pts, Fs = rng.normal(size=(4, 3)), rng.normal(size=(4, 3))
    A, B = rng.normal(size=3), rng.normal(size=3)
    R_A, M_A = impl.resultant(pts, Fs, None, A)
    R_B, M_B = impl.resultant(pts, Fs, None, B)
    npt.assert_allclose(R_A, R_B)
    npt.assert_allclose(impl.transport_moment(M_A, R_A, A, B), M_B, atol=1e-12)
    # a couple (equal and opposite forces) has the same moment about every point
    F = rng.normal(size=3)
    p1, p2 = rng.normal(size=3), rng.normal(size=3)
    R1, M1 = impl.resultant([p1, p2], [F, -F], None, A)
    R2, M2 = impl.resultant([p1, p2], [F, -F], None, B)
    npt.assert_allclose(R1, 0.0, atol=1e-14)
    npt.assert_allclose(M1, M2, atol=1e-12)
    npt.assert_allclose(M1, np.cross(p1 - p2, F), atol=1e-12)
    # a free couple adds to the moment and not to the force
    R3, M3 = impl.resultant([p1], [F], [[0.0, 0.0, 7.0]], A)
    npt.assert_allclose(M3 - impl.moment_about_point(p1, F, A), [0, 0, 7.0], atol=1e-12)


def test_central_axis_carries_a_moment_parallel_to_the_resultant():
    rng = np.random.default_rng(5016)
    pts, Fs = rng.normal(size=(3, 3)), rng.normal(size=(3, 3))
    R, M_O = impl.resultant(pts, Fs)
    p, pitch = impl.central_axis(R, M_O)
    M_p = impl.transport_moment(M_O, R, np.zeros(3), p)
    npt.assert_allclose(np.cross(M_p, R), 0.0, atol=1e-10)
    npt.assert_allclose(M_p, pitch * R, atol=1e-10)
    # a planar system has zero pitch: it reduces to a single force on the central axis
    pts2, Fs2 = rng.normal(size=(3, 2)), rng.normal(size=(3, 2))
    R2, M2 = impl.resultant(pts2, Fs2)
    p2, pitch2 = impl.central_axis(R2, M2)
    assert pitch2 == pytest.approx(0.0, abs=1e-12)
    npt.assert_allclose(impl.moment_about_point(p2, R2), M2, atol=1e-10)    # Varignon: R at p2 reproduces M
    with pytest.raises(ValueError):
        impl.central_axis(np.zeros(3), np.array([0.0, 0.0, 1.0]))


# ---------------------------------------------------------------- distributed loads

def test_uniform_and_triangular_distributed_loads():
    L, w0 = 2.0, 300.0
    x = np.linspace(0.0, L, 11)
    W, xb = impl.distributed_load_resultant(x, np.full_like(x, w0))
    assert W == pytest.approx(w0 * L) and xb == pytest.approx(L / 2)
    W, xb = impl.distributed_load_resultant(x, w0 * x / L)                     # linear: trapezoid exact
    assert W == pytest.approx(w0 * L / 2)
    x = np.linspace(0.0, L, 2001)
    W, xb = impl.distributed_load_resultant(x, w0 * x / L)
    assert xb == pytest.approx(2 * L / 3, rel=1e-6)                             # x w is quadratic: O(h^2)


def test_elliptic_lift_total_and_centroid_converge():
    s, Lh = 5.0, 5000.0
    errs = []
    for n in (41, 161, 641, 2561):
        y = np.linspace(0.0, s, n)
        W, yb = impl.distributed_load_resultant(y, impl.elliptic_lift(y, s, Lh))
        errs.append(abs(yb - 4 * s / (3 * np.pi)))
        assert W == pytest.approx(Lh, rel=2e-2)
    assert np.all(np.diff(errs) < 0)
    assert errs[-1] < 1e-4 * s
    # the square-root tip slows the trapezoid rule below O(h^2) but it still converges
    order = np.log(errs[0] / errs[-1]) / np.log(4 ** 3)
    assert order > 1.2


# ---------------------------------------------------------------- the free-body solver

def test_solve_fbd_2d_leaves_zero_residual_about_any_point():
    rng = np.random.default_rng(5017)
    kp, kf = rng.normal(size=(3, 2)), rng.normal(size=(3, 2)) * 1e3
    up = np.array([[0.0, 0.0], [0.0, 0.0], [2.0, 1.0]])
    ud = np.array([[1.0, 0.0], [0.0, 1.0], [0.3, 1.0]])
    C = 250.0
    x = impl.solve_fbd_2d(kp, kf, up, ud, known_couple=C)
    U = ud / np.linalg.norm(ud, axis=1, keepdims=True)
    pts = np.vstack([kp, up])
    Fs = np.vstack([kf, x[:, None] * U])
    for ref in (np.zeros(2), np.array([5.0, -3.0]), up[2]):
        npt.assert_allclose(impl.equilibrium_residual_2d(pts, Fs, C, ref), 0.0, atol=1e-8)
    x2 = impl.solve_fbd_2d(kp, kf, up, ud, known_couple=C, ref=np.array([-4.0, 9.0]))
    npt.assert_allclose(x, x2, rtol=1e-10, atol=1e-8)                           # reference point is free


def test_solve_fbd_2d_refuses_concurrent_or_parallel_unknowns():
    kp, kf = np.array([[1.0, 1.0]]), np.array([[0.0, -100.0]])
    concurrent = (np.zeros((3, 2)), np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]))
    parallel = (np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]]), np.array([[0.0, 1.0]] * 3))
    for up, ud in (concurrent, parallel):
        with pytest.raises(ValueError):
            impl.solve_fbd_2d(kp, kf, up, ud)


# ---------------------------------------------------------------- the landing gear

def test_gear_brace_force_by_system_equals_single_moment_equation():
    geom = impl.GearGeometry()
    rng = np.random.default_rng(5018)
    for _ in range(5):
        V, D = rng.uniform(5e3, 3e4), rng.uniform(-5e3, 1e4)
        r = impl.gear_reactions(geom, V, D)
        assert r["brace"] == pytest.approx(impl.brace_force_by_moments(geom, V, D), rel=1e-12)
    # the hand number: lever of the drag load about A is the leg length; brace lever from geometry
    V, D = 14563.0, 3641.0
    lever_D = geom.trunnion[1] - geom.contact[1]
    assert impl.brace_force_by_moments(geom, V, D) == pytest.approx(D * lever_D / impl.brace_lever_arm(geom))


def test_gear_free_body_is_in_equilibrium_and_ay_balances_vertical_loads():
    geom = impl.GearGeometry()
    V, D = impl.landing_load_case(1100.0)
    pts, Fs = impl.gear_fbd_forces(geom, V, D)
    scale = V * 1.2
    for ref in (geom.trunnion, geom.contact, geom.brace_frame, np.array([7.0, -3.0])):
        npt.assert_allclose(impl.equilibrium_residual_2d(pts, Fs, 0.0, ref), 0.0, atol=1e-10 * scale)
    r = impl.gear_reactions(geom, V, D)
    u = geom.brace_dir()
    assert r["Ay"] == pytest.approx(-V - r["brace"] * u[1])
    assert r["brace"] > 0.0                                                     # aft drag puts the brace in tension


def test_gear_superposition_and_zero_drag_limit():
    geom = impl.GearGeometry()
    V, D = 12000.0, 3000.0
    a, b, c = impl.gear_reactions(geom, V, 0.0), impl.gear_reactions(geom, 0.0, D), impl.gear_reactions(geom, V, D)
    for k in ("Ax", "Ay", "brace"):
        assert c[k] == pytest.approx(a[k] + b[k], abs=1e-8 * V)
        assert impl.gear_reactions(geom, 2 * V, 2 * D)[k] == pytest.approx(2 * c[k], abs=1e-8 * V)
    # contact directly below the trunnion and no drag: the ground force passes through A, the brace idles
    assert a["brace"] == pytest.approx(0.0, abs=1e-9 * V)
    assert a["Ay"] == pytest.approx(-V)
    assert impl.gear_reactions(geom, 0.0, 0.0)["brace"] == pytest.approx(0.0, abs=1e-12)


def test_brace_force_grows_without_bound_as_the_brace_line_approaches_the_trunnion():
    geom = impl.GearGeometry()
    V, D = 14563.0, 3641.0
    h = np.array([0.5, 0.8, 1.0, 1.1, 1.15, 1.19])
    T = impl.brace_sweep(geom, V, D, h)
    assert np.all(np.diff(T) > 0)
    assert T[-1] > 10 * T[0]
    lever = np.array([impl.brace_lever_arm(impl.GearGeometry(geom.contact, geom.trunnion, np.array([0.0, hh]),
                                                             geom.brace_frame)) for hh in h])
    npt.assert_allclose(T * lever, D * 1.2, rtol=1e-10)                        # T times its lever is constant
    bad = impl.GearGeometry(geom.contact, geom.trunnion, np.array([0.0, 0.5]), np.array([0.0, 1.5]))
    with pytest.raises(ValueError):
        impl.brace_force_by_moments(bad, V, D)                                  # brace along the leg: no lever


def test_landing_load_case_scales_with_mass_and_factor():
    V, D = impl.landing_load_case(1000.0, main_share=0.9, n_gear=2.0, drag_ratio=0.3)
    assert V == pytest.approx(2.0 * 0.45 * 1000.0 * 9.80665)
    assert D == pytest.approx(0.3 * V)
    V2, _ = impl.landing_load_case(2000.0, main_share=0.9, n_gear=2.0, drag_ratio=0.3)
    assert V2 == pytest.approx(2 * V)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "landing_gear_fbd.png").exists()
