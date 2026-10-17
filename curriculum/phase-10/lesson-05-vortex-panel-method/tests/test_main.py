"""Public tests for P10.05 Vortex Panel Method. They assert properties of the method, not one magic
number: the panel influence equals a direct quadrature of point vortices along the panel and looks
like a point vortex from far away; the solved system satisfies flow tangency and the Kutta
condition; a symmetric section at zero incidence has zero lift; the solution is linear in alpha;
cl from the circulation and from the surface pressure agree and the pressure drag tends to zero;
against the exact Karman-Trefftz flow cl converges at an observed order near 2 and Cp away from the
trailing edge converges; and as the thickness goes to zero the panel cl and cm tend to thin-airfoil
theory.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy import integrate

from conftest import load_impl

impl = load_impl(__file__)

MU = complex(-0.08, 0.06)
TAU = np.radians(10.0)


def polygon_circle(n):
    th = np.linspace(0, 2 * np.pi, n + 1)
    return np.cos(th), np.sin(th)


# ---------------------------------------------------------------- geometry and influence

def test_panel_geometry_normals_are_unit_outward_and_perpendicular():
    x, y = polygon_circle(24)
    g = impl.panel_geometry(x, y)
    npt.assert_allclose(np.hypot(g["nx"], g["ny"]), 1.0)
    npt.assert_allclose(g["nx"] * g["tx"] + g["ny"] * g["ty"], 0.0, atol=1e-15)
    assert np.all(g["nx"] * g["xc"] + g["ny"] * g["yc"] > 0)                # outward from the centre
    assert g["L"].sum() == pytest.approx(2 * np.pi, rel=0.02)


def test_linear_vortex_influence_matches_direct_quadrature():
    x, y = np.array([0.2, 1.1]), np.array([0.1, 0.5])                    # one panel
    g = impl.panel_geometry(x, y)
    pts = np.array([[0.9, -0.4], [0.3, 0.9], [-1.0, 0.2]])
    ua, va, ub, vb = impl.linear_vortex_influence(pts[:, 0], pts[:, 1], g)
    L = g["L"][0]
    for k, (px, py) in enumerate(pts):
        def vel(s, which):
            gs = (1 - s / L) if which == "a" else s / L
            xs, ys = x[0] + s * g["tx"][0], y[0] + s * g["ty"][0]
            rx, ry = px - xs, py - ys
            r2 = rx ** 2 + ry ** 2
            return gs * ry / (2 * np.pi * r2), -gs * rx / (2 * np.pi * r2)   # clockwise vortex
        for which, (u, v) in (("a", (ua, va)), ("b", (ub, vb))):
            uq = integrate.quad(lambda s: vel(s, which)[0], 0, L)[0]
            vq = integrate.quad(lambda s: vel(s, which)[1], 0, L)[0]
            assert u[k, 0] == pytest.approx(uq, abs=1e-10)
            assert v[k, 0] == pytest.approx(vq, abs=1e-10)


def test_panel_looks_like_a_point_vortex_from_far_away():
    x, y = np.array([0.0, 0.4]), np.array([0.0, 0.3])
    g = impl.panel_geometry(x, y)
    P = np.array([300.0]), np.array([-200.0])
    ua, va, ub, vb = impl.linear_vortex_influence(*P, g)
    G = g["L"][0]                                                          # uniform unit strength
    rx, ry = P[0] - g["xc"][0], P[1] - g["yc"][0]
    r2 = rx ** 2 + ry ** 2
    assert (ua + ub)[0, 0] == pytest.approx(G * ry / (2 * np.pi * r2), rel=1e-5)
    assert (va + vb)[0, 0] == pytest.approx(-G * rx / (2 * np.pi * r2), rel=1e-5)


# ---------------------------------------------------------------- the solution

@pytest.fixture(scope="module")
def naca2412():
    x, y = impl.naca4("2412", 80)
    return x, y, impl.vortex_panel(x, y, np.radians(4.0))


def test_boundary_condition_and_kutta_condition_hold(naca2412):
    x, y, s = naca2412
    assert s["gamma"][0] + s["gamma"][-1] == pytest.approx(0.0, abs=1e-12)
    g = s["geo"]
    ua, va, ub, vb = impl.linear_vortex_influence(g["xc"], g["yc"], g, outside_self=True)
    U, V = impl.node_influence(ua, va, ub, vb)
    a = s["alpha"]
    vn = (U @ s["gamma"] + np.cos(a)) * g["nx"] + (V @ s["gamma"] + np.sin(a)) * g["ny"]
    npt.assert_allclose(vn, 0.0, atol=1e-10)
    assert s["residual"] < 1e-10


def test_symmetric_section_at_zero_incidence_has_no_lift():
    x, y = impl.naca4("0012", 60)
    s = impl.vortex_panel(x, y, 0.0)
    assert abs(s["cl_gamma"]) < 1e-12 and abs(s["cl_p"]) < 1e-10 and abs(s["cm"]) < 1e-10
    npt.assert_allclose(s["cp"], s["cp"][::-1], atol=1e-10)


def test_solution_superposes_two_free_streams():
    # the system is linear in the free stream (cos alpha, sin alpha): the circulation at any alpha is
    # the alpha = 0 solution times cos(alpha) plus the alpha = 90 deg solution times sin(alpha)
    x, y = impl.naca4("4412", 50)
    G0, G90 = (impl.vortex_panel(x, y, np.radians(a))["Gamma"] for a in (0.0, 90.0))
    for a in (-5.0, 3.0, 11.0):
        G = impl.vortex_panel(x, y, np.radians(a))["Gamma"]
        assert G == pytest.approx(G0 * np.cos(np.radians(a)) + G90 * np.sin(np.radians(a)), rel=1e-10)
    assert G0 > 0                                                         # camber: lift at zero alpha


def test_lift_from_gamma_and_pressure_agree_and_pressure_drag_vanishes():
    prev = None
    for n in (40, 80, 160):
        x, y = impl.naca4("2412", n)
        s = impl.vortex_panel(x, y, np.radians(6.0))
        assert s["cl_p"] == pytest.approx(s["cl_gamma"], rel=5e-3)
        if prev is not None:
            assert abs(s["cd_p"]) < abs(prev)                             # d'Alembert, approached
        prev = s["cd_p"]
    assert abs(prev) < 2e-3


# ---------------------------------------------------------------- against the exact solution

@pytest.fixture(scope="module")
def kt_study():
    return impl.kt_convergence(1.0, MU, TAU, np.radians(4.0), [20, 40, 80, 160, 320])


def test_cl_converges_to_karman_trefftz_at_second_order(kt_study):
    p = impl.observed_order(kt_study["e_clg"])
    assert 1.6 < np.mean(p) < 2.4                                         # second order on average
    assert kt_study["e_clg"][-1] < 1e-4                                   # N = 320
    assert 1.6 < np.mean(impl.observed_order(kt_study["e_clp"])) < 2.4


def test_cp_converges_away_from_the_trailing_edge(kt_study):
    e = kt_study["e_cp"]
    assert np.all(np.diff(e) < 0)                                         # monotone
    assert e[-1] < 1e-3
    assert np.mean(impl.observed_order(e)[-2:]) > 1.6                     # second order once resolved


def test_exact_reference_is_self_consistent():
    al = np.radians(3.0)
    a, beta = impl.circle_geometry(1.0, MU)
    th = -beta + 2 * np.pi * (np.arange(2048) + 0.5) / 2048
    ex = impl.kt_exact(1.0, MU, TAU, al, th)
    assert ex["cl"] == pytest.approx(8 * np.pi * a * np.sin(al + beta) / ex["chord"])
    assert ex["cp"].max() <= 1.0 + 1e-12


# ---------------------------------------------------------------- thin-airfoil limit

def test_thin_airfoil_limit_of_cambered_sections():
    al = np.radians(4.0)
    cl_ta, cm_ta, aL0 = impl.thin_airfoil_naca4(0.02, 0.4, al)
    dcl, dcm = [], []
    for t in (0.12, 0.06, 0.03, 0.015):
        x, y = impl.naca4_mpt(0.02, 0.4, t, 120)
        s = impl.vortex_panel(x, y, al)
        dcl.append(s["cl_gamma"] - cl_ta)
        dcm.append(s["cm"] - cm_ta)
    dcl, dcm = np.array(dcl), np.array(dcm)
    assert np.all(dcl > 0) and np.all(np.diff(dcl) < 0)                   # thickness adds lift, less and less
    assert dcl[-1] < 0.015 * cl_ta                                        # within 1.5 % at t/c = 1.5 %
    assert np.all(np.diff(np.abs(dcm)) < 0) and abs(dcm[-1]) < 0.03 * abs(cm_ta)


def test_thin_symmetric_section_tends_to_two_pi_alpha():
    al = np.radians(3.0)
    x, y = impl.naca4_mpt(0.0, 0.0, 0.01, 120)
    s = impl.vortex_panel(x, y, al)
    assert s["cl_gamma"] == pytest.approx(2 * np.pi * al, rel=0.01)
    assert abs(s["cm"]) < 2e-3


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path, n_max=160)
    assert (tmp_path / "vortex_panel.png").stat().st_size > 1000
