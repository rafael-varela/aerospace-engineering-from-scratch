"""Public tests for P08.09 Vorticity & Circulation. They assert the theorems, not magic numbers: the
Lamb-Oseen vortex has the vorticity that is the curl of its velocity, holds total circulation Gamma,
tends to a point vortex outside its core, and satisfies the vorticity diffusion equation; Stokes's
theorem holds (line integral = area integral of vorticity); circulation round an irrotational
region is zero and does not depend on the loop's shape; RK4 reproduces the exact particle paths at
fourth order; Kelvin's theorem holds for an advected material loop in inviscid flow, to a
discretisation error that falls as the square of the segment length, and fails as Lamb-Oseen
predicts with viscosity; Biot-Savart reproduces the infinite and semi-infinite line vortices, closes
Stokes's theorem in 3-D, and gives the horseshoe downwash formulas; and a starting vortex makes the
total circulation zero.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import quad

from conftest import load_impl

impl = load_impl(__file__)

G, RC2 = 2.5, 0.49


def vel(rc2=RC2, Gamma=G):
    return lambda p: impl.vortex_velocity(p, Gamma, rc2)


# ---------------------------------------------------------------- the Lamb-Oseen vortex

def test_lamb_oseen_vorticity_is_the_curl_of_its_velocity_and_integrates_to_gamma():
    h = 1e-5
    pts = np.array([[0.3, 0.2], [-0.7, 0.5], [1.1, -0.9], [1e-3, 0.0]])
    dvdx = (impl.vortex_velocity(pts + [h, 0], G, RC2)[:, 1] - impl.vortex_velocity(pts - [h, 0], G, RC2)[:, 1]) / (2 * h)
    dudy = (impl.vortex_velocity(pts + [0, h], G, RC2)[:, 0] - impl.vortex_velocity(pts - [0, h], G, RC2)[:, 0]) / (2 * h)
    npt.assert_allclose(dvdx - dudy, impl.lamb_oseen_vorticity(np.hypot(pts[:, 0], pts[:, 1]), G, RC2), rtol=1e-6)
    # divergence-free
    dudx = (impl.vortex_velocity(pts + [h, 0], G, RC2)[:, 0] - impl.vortex_velocity(pts - [h, 0], G, RC2)[:, 0]) / (2 * h)
    dvdy = (impl.vortex_velocity(pts + [0, h], G, RC2)[:, 1] - impl.vortex_velocity(pts - [0, h], G, RC2)[:, 1]) / (2 * h)
    assert np.max(np.abs(dudx + dvdy)) < 1e-6
    total = quad(lambda r: 2 * np.pi * r * float(impl.lamb_oseen_vorticity(r, G, RC2)), 0, np.inf)[0]
    assert total == pytest.approx(G, rel=1e-10)


def test_lamb_oseen_limits_axis_and_point_vortex():
    ut = impl.lamb_oseen_utheta(np.array([0.0, 1e-9, 20.0]), G, RC2)
    assert ut[0] == 0.0 and abs(ut[1]) < 1e-8                          # solid-body rotation on the axis
    assert ut[2] == pytest.approx(G / (2 * np.pi * 20.0), rel=1e-12)    # point vortex far away
    r_max = impl.R_MAX_FACTOR * np.sqrt(RC2)                            # peak swirl
    e = 1e-4
    u = impl.lamb_oseen_utheta(np.array([r_max - e, r_max, r_max + e]), G, RC2)
    assert u[1] > u[0] and u[1] > u[2]
    assert impl.core_radius2(0.3, 1.5e-5, 100.0) == pytest.approx(0.09 + 4 * 1.5e-5 * 100.0)


def test_lamb_oseen_satisfies_the_vorticity_diffusion_equation():
    # omega_t = nu (omega_rr + omega_r / r) with rc2 = rc0^2 + 4 nu t
    nu, rc0, t, dt, dr = 0.02, 0.5, 3.0, 1e-4, 1e-4
    w = lambda r, tt: impl.lamb_oseen_vorticity(r, G, impl.core_radius2(rc0, nu, tt))
    r = np.array([0.2, 0.6, 1.3])
    w_t = (w(r, t + dt) - w(r, t - dt)) / (2 * dt)
    w_rr = (w(r + dr, t) - 2 * w(r, t) + w(r - dr, t)) / dr ** 2
    w_r = (w(r + dr, t) - w(r - dr, t)) / (2 * dr)
    npt.assert_allclose(w_t, nu * (w_rr + w_r / r), rtol=1e-5)


# ---------------------------------------------------------------- Stokes's theorem and circulation

def test_stokes_theorem_line_integral_equals_area_integral():
    for c, R in (((0.0, 0.0), 0.8), ((0.9, 0.3), 0.6), ((-1.5, 2.0), 0.4)):
        line = impl.circulation_polygon(impl.circle_loop(c, R, 2000), vel())
        area = impl.disk_area_integral(lambda p: impl.lamb_oseen_vorticity(np.hypot(p[:, 0], p[:, 1]), G, RC2), c, R)
        assert line == pytest.approx(area, rel=1e-5)      # midpoint rule on chords: O((2 pi / n)^2)
    # a concentric circle: the closed form
    assert impl.circulation_polygon(impl.circle_loop((0, 0), 0.8, 2000), vel()) == pytest.approx(
        float(impl.lamb_oseen_circulation(0.8, G, RC2)), rel=1e-5)


def test_disk_quadrature_is_exact_for_polynomials():
    # x^2 + y^2 over the unit disk = pi / 2; a linear function gives area * value at the centre
    assert impl.disk_area_integral(lambda p: p[:, 0] ** 2 + p[:, 1] ** 2, (0, 0), 1.0) == pytest.approx(np.pi / 2, rel=1e-13)
    assert impl.disk_area_integral(lambda p: 3 + p[:, 0], (2.0, 1.0), 0.5) == pytest.approx(5 * np.pi * 0.25, rel=1e-13)


def test_circulation_in_irrotational_flow_depends_only_on_what_is_enclosed():
    vort = [(0.0, 0.0, 1.7)]
    v = lambda p: impl.point_vortices_velocity(p, vort)
    sq = np.array([[x, -2.0] for x in np.linspace(-2, 2, 400, endpoint=False)]
                  + [[2.0, y] for y in np.linspace(-2, 2, 400, endpoint=False)]
                  + [[x, 2.0] for x in np.linspace(2, -2, 400, endpoint=False)]
                  + [[-2.0, y] for y in np.linspace(2, -2, 400, endpoint=False)])
    for loop in (impl.circle_loop((0.2, -0.1), 1.0, 2000), impl.circle_loop((0.0, 0.0), 5.0, 2000), sq):
        assert impl.circulation_polygon(loop, v) == pytest.approx(1.7, rel=1e-5)
    assert abs(impl.circulation_polygon(impl.circle_loop((3.0, 0.0), 1.0, 2000), v)) < 1e-10


# ---------------------------------------------------------------- material loops: Kelvin

def test_turning_angle_matches_quadrature_and_rk4_converges_at_fourth_order():
    nu, rc0, t = 0.01, 1.0, 40.0
    r = np.array([0.3, 0.9, 1.5])
    ref = [quad(lambda s: float(impl.lamb_oseen_utheta(np.array([ri]), 1.0, impl.core_radius2(rc0, nu, s))[0]) / ri, 0, t)[0]
           for ri in r]
    npt.assert_allclose(impl.turning_angle(r, t, 1.0, rc0, nu), ref, rtol=1e-10)
    pts = impl.circle_loop((0.9, 0.0), 0.6, 12)
    exact = impl.material_loop_exact(pts, t, 1.0, rc0, nu)
    errs = [np.max(np.linalg.norm(impl.advect_rk4(pts, 0.0, t, n, 1.0, rc0, nu) - exact, axis=1)) for n in (10, 20, 40)]
    assert np.all(np.abs(impl.observed_order(errs) - 4.0) < 0.3)
    # particles keep their distance from the axis (purely azimuthal flow)
    npt.assert_allclose(np.hypot(*exact.T), np.hypot(*pts.T), rtol=1e-13)


def test_kelvin_circulation_of_a_material_loop_is_constant_without_viscosity():
    c0, R0, T = (0.9, 0.0), 0.6, 120.0
    changes = []
    for n in (1000, 2000):
        loop0 = impl.circle_loop(c0, R0, n)
        G0 = impl.circulation_polygon(loop0, vel(1.0, 1.0))
        lp = impl.material_loop_exact(loop0, T, 1.0, 1.0, 0.0)
        G1 = impl.circulation_polygon(lp, vel(1.0, 1.0))
        changes.append(abs(G1 - G0))
        length0 = np.sum(np.linalg.norm(np.roll(loop0, -1, 0) - loop0, axis=1))
        length1 = np.sum(np.linalg.norm(np.roll(lp, -1, 0) - lp, axis=1))
        assert length1 > 4 * length0                                    # the loop has been stretched hard
    assert changes[-1] < 2e-5 * G0
    assert changes[0] / changes[1] == pytest.approx(4.0, rel=0.1)        # the change is discretisation, ~ ds^2
    # the same with the loop advected by RK4 instead of the exact map
    loop0 = impl.circle_loop(c0, R0, 2000)
    lp = impl.advect_rk4(loop0, 0.0, T, 200, 1.0, 1.0, 0.0)
    assert impl.circulation_polygon(lp, vel(1.0, 1.0)) == pytest.approx(G0, rel=2e-5)


def test_viscosity_changes_the_circulation_as_lamb_oseen_predicts():
    nu, t = 0.01, 60.0
    rc2 = impl.core_radius2(1.0, nu, t)
    # concentric material loop: stays a circle, circulation Gamma (1 - exp(-R^2 / rc2))
    lp = impl.material_loop_exact(impl.circle_loop((0, 0), 0.8, 2000), t, 1.0, 1.0, nu)
    G_line = impl.circulation_polygon(lp, vel(rc2, 1.0))
    assert G_line == pytest.approx(float(impl.lamb_oseen_circulation(0.8, 1.0, rc2)), rel=1e-5)
    assert G_line < float(impl.lamb_oseen_circulation(0.8, 1.0, 1.0))
    # off-centre loop: line integral round the material loop = area integral of the current vorticity
    c0, R0 = (0.9, 0.0), 0.6
    lp = impl.material_loop_exact(impl.circle_loop(c0, R0, 4000), t, 1.0, 1.0, nu)
    G_line = impl.circulation_polygon(lp, vel(rc2, 1.0))
    G_area = impl.material_circulation_by_area(c0, R0, t, 1.0, 1.0, nu)
    assert G_line == pytest.approx(G_area, rel=1e-5)
    assert G_area < 0.8 * impl.material_circulation_by_area(c0, R0, 0.0, 1.0, 1.0, nu)


# ---------------------------------------------------------------- Biot-Savart

def test_biot_savart_limits_and_symmetry():
    P = np.array([[0.0, 0.0, 0.7]])
    u = impl.biot_savart_segment(P, np.array([-1e5, 0, 0]), np.array([1e5, 0, 0]), 3.0)[0]
    # infinite line along +x: speed Gamma/(2 pi d), direction x_hat x z_hat = -y_hat
    assert u[1] == pytest.approx(-3.0 / (2 * np.pi * 0.7), rel=1e-8) and abs(u[0]) < 1e-15 and abs(u[2]) < 1e-15
    u_half = impl.biot_savart_segment(P, np.array([0.0, 0, 0]), np.array([1e5, 0, 0]), 3.0)[0]
    assert u_half[1] == pytest.approx(u[1] / 2, rel=1e-6)                 # semi-infinite: half
    u_rev = impl.biot_savart_segment(P, np.array([1e5, 0, 0]), np.array([-1e5, 0, 0]), 3.0)[0]
    npt.assert_allclose(u_rev, -u, rtol=1e-12)                            # reversing the filament
    # a point on the line of the filament gets nothing (no self-induction)
    assert np.all(impl.biot_savart_segment(np.array([[2.0, 0, 0]]), np.array([0.0, 0, 0]), np.array([1.0, 0, 0]), 3.0) == 0)


def test_biot_savart_satisfies_stokes_round_the_filament():
    A, B, Gam = np.array([0.0, 0.0, -50.0]), np.array([0.0, 0.0, 50.0]), 1.3
    n = 2000
    th = 2 * np.pi * np.arange(n) / n
    loop = np.column_stack([0.5 * np.cos(th) + 0.1, 0.5 * np.sin(th), np.full(n, 0.2)])
    nxt = np.roll(loop, -1, 0)
    u = impl.biot_savart_segment(0.5 * (loop + nxt), A, B, Gam)
    assert np.sum(u * (nxt - loop)) == pytest.approx(Gam, rel=1e-3)        # finite length 100: ~1e-4 short
    far = np.column_stack([0.5 * np.cos(th) + 3.0, 0.5 * np.sin(th), np.full(n, 0.2)])
    nf = np.roll(far, -1, 0)
    assert abs(np.sum(impl.biot_savart_segment(0.5 * (far + nf), A, B, Gam) * (nf - far))) < 1e-4 * Gam


def test_horseshoe_downwash_formulas():
    bp, Gam = 7.854, 17.0
    w0 = impl.horseshoe_velocity(np.array([[0.0, 0.0, 0.0]]), bp, Gam)[0]
    assert w0[2] == pytest.approx(-Gam / (np.pi * bp), rel=1e-6)           # two half-infinite legs
    wf = impl.horseshoe_velocity(np.array([[5e3, 0.0, 0.0]]), bp, Gam, length=1e7)[0]
    assert wf[2] == pytest.approx(-2 * Gam / (np.pi * bp), rel=1e-4)       # two infinite legs
    up = impl.horseshoe_velocity(np.array([[-3.0, 0.0, 0.0]]), bp, Gam)[0]
    assert up[2] > 0                                                       # upwash ahead of the wing
    outside = impl.horseshoe_velocity(np.array([[2.0, 1.2 * bp, 0.0]]), bp, Gam)[0]
    assert outside[2] > 0                                                  # upwash outboard of the tips


def test_starting_vortex_keeps_the_total_circulation_zero():
    Gam = 17.0
    vort = [(0.0, 0.0, -Gam), (50.0, 0.0, Gam)]
    v = lambda p: impl.point_vortices_velocity(p, vort)
    assert abs(impl.circulation_polygon(impl.circle_loop((25.0, 0.0), 40.0, 4000), v)) < 1e-9
    assert impl.circulation_polygon(impl.circle_loop((0.0, 0.0), 2.0, 2000), v) == pytest.approx(-Gam, rel=1e-5)
    # the starting vortex induces a downwash at the wing that fades as it is left behind
    w = [impl.point_vortices_velocity(np.array([[0.0, 0.0]]), [(d, 0.0, Gam)])[0, 1] for d in (5.0, 50.0)]
    assert w[0] == pytest.approx(-Gam / (2 * np.pi * 5.0)) and abs(w[1]) < abs(w[0])


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "vorticity.png").stat().st_size > 10_000
