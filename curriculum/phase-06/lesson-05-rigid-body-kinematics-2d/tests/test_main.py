"""Public tests for P06.05 Rigid Body Kinematics in 2D. They assert kinematics, not magic numbers: the
Grashof classes; Newton's method closes the loop, converges quadratically and agrees with the
circle-intersection closed form on both branches; link lengths are preserved; the Jacobian is singular
at toggle; the velocity solution satisfies v_B = v_A + w3 x r_BA and agrees with instant centres
(Kennedy's theorem: I13 and I24 on their lines); a parallelogram has w3 = 0 and w4 = w2; power is
conserved through the gearing; and the independent check: finite differences of the position solution
converge at second order to the analytic angular velocities and accelerations, including a coupler point.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

LINKS = (0.100, 0.040, 0.120, 0.080)          # crank-rocker: 40 + 120 < 100 + 80


def _cr(w, r):
    return w * np.array([-r[1], r[0]])


# ---------------------------------------------------------------- classification

@pytest.mark.parametrize("links,kind", [
    ((0.100, 0.040, 0.120, 0.080), "crank-rocker"),
    ((0.040, 0.100, 0.120, 0.080), "double-crank"),
    ((0.100, 0.080, 0.040, 0.120), "double-rocker"),
    ((0.100, 0.080, 0.120, 0.040), "rocker-crank"),
    ((0.100, 0.050, 0.100, 0.050), "change-point"),
    ((0.100, 0.090, 0.120, 0.095), "triple-rocker"),
])
def test_grashof_classes(links, kind):
    assert impl.grashof(*links) == kind


# ---------------------------------------------------------------- position

def test_newton_closes_the_loop_and_preserves_link_lengths():
    for th2 in np.radians([0.0, 37.0, 140.0, 260.0]):
        g = impl.position_circle_intersection(th2, LINKS, branch=1)
        th3, th4, it = impl.solve_position_newton(th2, LINKS, (g[0] + 0.2, g[1] - 0.2))
        assert np.linalg.norm(impl.loop_residual(th2, th3, th4, LINKS)) < 1e-13
        P = impl.joint_positions(th2, th3, th4, LINKS)
        assert np.linalg.norm(P["A"] - P["O2"]) == pytest.approx(LINKS[1])
        assert np.linalg.norm(P["B"] - P["A"]) == pytest.approx(LINKS[2])
        assert np.linalg.norm(P["B"] - P["O4"]) == pytest.approx(LINKS[3])
        assert it <= 8


def test_newton_agrees_with_circle_intersection_on_both_branches():
    th2 = np.radians(75.0)
    for branch in (1, -1):
        g3, g4 = impl.position_circle_intersection(th2, LINKS, branch)
        th3, th4, _ = impl.solve_position_newton(th2, LINKS, (g3 + 0.1, g4 + 0.1))
        assert np.cos(th3 - g3) == pytest.approx(1.0, abs=1e-12)
        assert np.cos(th4 - g4) == pytest.approx(1.0, abs=1e-12)
    # the two branches are mirror images about the line A -> O4: different rocker angles
    assert abs(impl.position_circle_intersection(th2, LINKS, 1)[1] - impl.position_circle_intersection(th2, LINKS, -1)[1]) > 0.5


def test_newton_converges_quadratically():
    th2 = np.radians(50.0)
    ex = np.array(impl.position_circle_intersection(th2, LINKS, 1))
    x = ex + np.array([0.15, -0.1])
    errs = []
    for _ in range(4):
        x = x - np.linalg.solve(impl.loop_jacobian(x[0], x[1], LINKS), impl.loop_residual(th2, x[0], x[1], LINKS))
        errs.append(np.linalg.norm(x - ex))
    # e_{k+1} ~ C e_k^2: the ratio stays bounded while the errors fall by orders of magnitude
    assert errs[2] < 1e-6 and errs[2] / errs[1] ** 2 < 10.0 and errs[1] / errs[0] ** 2 < 10.0


def test_jacobian_against_finite_differences_and_singular_at_toggle():
    th2, th3, th4 = 0.7, 0.3, 1.4
    J = impl.loop_jacobian(th3, th4, LINKS)
    h = 1e-7
    Jn = np.column_stack([(impl.loop_residual(th2, th3 + h, th4, LINKS) - impl.loop_residual(th2, th3 - h, th4, LINKS)) / (2 * h),
                          (impl.loop_residual(th2, th3, th4 + h, LINKS) - impl.loop_residual(th2, th3, th4 - h, LINKS)) / (2 * h)])
    npt.assert_allclose(J, Jn, atol=1e-8)
    assert abs(np.linalg.det(impl.loop_jacobian(0.9, 0.9, LINKS))) < 1e-15      # coupler and rocker in line


def test_assembly_failure_is_reported():
    # a triple-rocker cannot reach every crank angle
    L = (0.100, 0.090, 0.120, 0.095)
    with pytest.raises(ValueError):
        impl.position_circle_intersection(np.pi, (0.100, 0.090, 0.050, 0.030), 1)
    with pytest.raises(ValueError):                    # at th2 = 0, |r3 - r4| > |O4 - A|: no triangle
        impl.position_circle_intersection(0.0, L, 1)
    th3, th4 = impl.sweep_positions(np.radians(np.linspace(60, 120, 7)), L,
                                    impl.position_circle_intersection(np.radians(60), L, 1))
    assert np.all(np.isfinite(th4))


# ---------------------------------------------------------------- velocity and acceleration

def test_velocity_satisfies_relative_velocity_equation_and_instant_centres():
    th2, w2 = np.radians(120.0), 5.0
    th3, th4 = impl.position_circle_intersection(th2, LINKS, 1)
    w3, w4 = impl.velocity_analysis(th2, th3, th4, w2, LINKS)
    P = impl.joint_positions(th2, th3, th4, LINKS)
    vA = _cr(w2, P["A"])
    vB_coupler = vA + _cr(w3, P["B"] - P["A"])
    vB_rocker = _cr(w4, P["B"] - P["O4"])
    npt.assert_allclose(vB_coupler, vB_rocker, atol=1e-12)
    r32, r42 = impl.velocity_ratio_ic(th2, th3, th4, LINKS)
    assert r32 == pytest.approx(w3 / w2, rel=1e-10) and r42 == pytest.approx(w4 / w2, rel=1e-10)
    ic = impl.instant_centres(th2, th3, th4, LINKS)
    # Kennedy: I13 on line O2-A and on line O4-B; I24 on the ground line and on line A-B
    c = lambda p, q, r: (q - p)[0] * (r - p)[1] - (q - p)[1] * (r - p)[0]
    for p, q in ((ic["I12"], ic["I23"]), (ic["I14"], ic["I34"])):
        assert abs(c(p, q, ic["I13"])) < 1e-12
    assert abs(ic["I24"][1]) < 1e-14 and abs(c(ic["I23"], ic["I34"], ic["I24"])) < 1e-12
    # every coupler point moves as if pinned at I13
    rP, vP, _ = impl.coupler_point(th2, th3, w2, 0.0, w3, 0.0, LINKS, 0.05, 0.6)
    npt.assert_allclose(vP, _cr(w3, rP - ic["I13"]), atol=1e-12)


def test_parallelogram_coupler_translates():
    L = (0.100, 0.040, 0.100, 0.040)
    for th2 in (0.3, 1.2, 2.5):
        th3, th4 = impl.position_circle_intersection(th2, L, 1)
        npt.assert_allclose([np.cos(th3), np.sin(th3)], [1.0, 0.0], atol=1e-12)   # coupler stays parallel to ground
        assert np.cos(th4 - th2) == pytest.approx(1.0)
        w3, w4 = impl.velocity_analysis(th2, th3, th4, 3.0, L)
        assert abs(w3) < 1e-12 and w4 == pytest.approx(3.0)
        a3, a4 = impl.acceleration_analysis(th2, th3, th4, 3.0, 2.0, w3, w4, L)
        assert abs(a3) < 1e-10 and a4 == pytest.approx(2.0)


def test_velocities_are_linear_in_crank_rate_and_power_is_conserved():
    th2 = 0.9
    th3, th4 = impl.position_circle_intersection(th2, LINKS, 1)
    w3a, w4a = impl.velocity_analysis(th2, th3, th4, 1.0, LINKS)
    w3b, w4b = impl.velocity_analysis(th2, th3, th4, -2.5, LINKS)
    assert w3b == pytest.approx(-2.5 * w3a) and w4b == pytest.approx(-2.5 * w4a)
    assert impl.velocity_analysis(th2, th3, th4, 0.0, LINKS) == (0.0, 0.0)
    # ideal massless linkage: T2 w2 = M4 w4, so the torque ratio is w4/w2 (checked by virtual work
    # on the loop: a small crank rotation moves the rocker by (w4/w2) times as much)
    d = 1e-6
    t3b, t4b = impl.solve_position_newton(th2 + d, LINKS, (th3, th4))[:2]
    assert (t4b - th4) / d == pytest.approx(w4a, rel=1e-5)


def test_transmission_angle():
    assert impl.transmission_angle(0.2, 0.2 + np.pi / 2) == pytest.approx(np.pi / 2)
    assert impl.transmission_angle(0.2, 0.2 + np.radians(130)) == pytest.approx(np.radians(50))
    assert impl.transmission_angle(1.0, 1.0) == pytest.approx(0.0)
    # law of cosines across the diagonal B-O2 ... via the diagonal A-O4
    th2 = 1.1
    th3, th4 = impl.position_circle_intersection(th2, LINKS, 1)
    r1, r2, r3, r4 = LINKS
    cg = (r3 ** 2 + r4 ** 2 - r1 ** 2 - r2 ** 2 + 2 * r1 * r2 * np.cos(th2)) / (2 * r3 * r4)
    g = np.arccos(cg)
    assert impl.transmission_angle(th3, th4) == pytest.approx(min(g, np.pi - g), rel=1e-10)


# ---------------------------------------------------------------- independent check: finite differences

def _history(t, w2, al2, th20=0.4):
    th2 = th20 + w2 * t + 0.5 * al2 * t ** 2
    g = impl.position_circle_intersection(float(th2[0]), LINKS, 1)
    th3, th4 = impl.sweep_positions(th2, LINKS, g)
    return th2, th3, th4


def test_finite_differences_converge_to_analytic_rates_at_second_order():
    w2, al2, t0 = 4.0, 6.0, 0.1
    th2c, th3c, th4c = _history(np.array([t0]), w2, al2)
    W3, W4 = impl.velocity_analysis(th2c[0], th3c[0], th4c[0], w2 + al2 * t0, LINKS)
    A3, A4 = impl.acceleration_analysis(th2c[0], th3c[0], th4c[0], w2 + al2 * t0, al2, W3, W4, LINKS)
    ev, ea = [], []
    for h in (0.008, 0.004, 0.002, 0.001):
        t = np.array([t0 - h, t0, t0 + h])
        _, th3, th4 = _history(t, w2, al2)
        _, d1, d2 = impl.central_rates(t, th4)
        _, e1, e2 = impl.central_rates(t, th3)
        ev.append(abs(d1[0] - W4) + abs(e1[0] - W3))
        ea.append(abs(d2[0] - A4) + abs(e2[0] - A3))
    assert np.all(np.abs(impl.observed_order(ev) - 2) < 0.1)
    assert np.all(np.abs(impl.observed_order(ea) - 2) < 0.1)
    assert ev[-1] < 1e-4 * abs(W4) and ea[-1] < 1e-3 * abs(A4)          # h = 1 ms


def test_coupler_point_acceleration_matches_finite_differences():
    w2, al2, t0, h = 3.0, -2.0, 0.2, 1e-3
    t = np.array([t0 - h, t0, t0 + h])
    th2, th3, th4 = _history(t, w2, al2)
    p, dl = 0.07, -0.5
    pos = np.array([impl.coupler_point(th2[i], th3[i], 0, 0, 0, 0, LINKS, p, dl)[0] for i in range(3)])
    w2c = w2 + al2 * t0
    w3, w4 = impl.velocity_analysis(th2[1], th3[1], th4[1], w2c, LINKS)
    a3, a4 = impl.acceleration_analysis(th2[1], th3[1], th4[1], w2c, al2, w3, w4, LINKS)
    _, vP, aP = impl.coupler_point(th2[1], th3[1], w2c, al2, w3, a3, LINKS, p, dl)
    npt.assert_allclose((pos[2] - pos[0]) / (2 * h), vP, rtol=1e-5, atol=1e-7)
    npt.assert_allclose((pos[2] - 2 * pos[1] + pos[0]) / h ** 2, aP, rtol=1e-4, atol=1e-5)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "four_bar.png").stat().st_size > 10_000
