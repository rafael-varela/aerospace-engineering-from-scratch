"""Public tests for P10.04 Thin Airfoil Theory. They assert the theory's structure and its checks: a
flat plate gives cl = 2 pi alpha with the centre of pressure at the quarter chord; a parabolic camber
line reproduces its closed forms; cl is linear in alpha with slope 2 pi and cm_c/4 does not depend on
alpha; the zero-lift angle zeroes cl; quadrature split at a kink is exact to round-off; the plain
flap reproduces tau = 1 - (theta_h - sin theta_h)/pi and its moment closed form, with the right
limits; camber and flap superpose; the 230 camber line has its design ideal cl; an independent
discrete-vortex model converges to the Glauert result at second order; and the vortex panel method
tends to thin-airfoil theory as the thickness goes to zero.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy import integrate

from conftest import load_impl

impl = load_impl(__file__)

FLAT = lambda x: np.zeros_like(np.asarray(x, float))


def test_flat_plate():
    for a in (-0.1, 0.05, 0.2):
        r = impl.thin_airfoil(FLAT, a)
        assert r["cl"] == pytest.approx(2 * np.pi * a)
        assert abs(r["cm_c4"]) < 1e-15 and r["alpha_L0"] == 0.0
        assert r["x_cp"] == pytest.approx(0.25)
        npt.assert_allclose(r["A"][1:], 0.0, atol=1e-15)


def test_parabolic_camber_closed_forms():
    h = 0.03                                     # z = 4 h x (1 - x), dz/dx = 4 h cos(theta)
    f = lambda x: 4 * h * (1 - 2 * np.asarray(x))
    r = impl.thin_airfoil(f, np.radians(3.0))
    assert r["A"][1] == pytest.approx(4 * h)
    assert abs(r["A"][2]) < 1e-14 and abs(r["A"][3]) < 1e-14
    assert r["alpha_L0"] == pytest.approx(-2 * h)
    assert r["cm_c4"] == pytest.approx(-np.pi * h)
    assert r["cl"] == pytest.approx(2 * np.pi * (np.radians(3.0) + 2 * h))


def test_lift_is_linear_with_2pi_slope_and_cm_is_constant():
    f, br = impl.naca4_camber_slope(0.04, 0.3)
    rs = [impl.thin_airfoil(f, np.radians(a), br) for a in (-2.0, 1.0, 7.0)]
    cl = np.array([r["cl"] for r in rs])
    npt.assert_allclose(np.diff(cl) / np.radians([3.0, 6.0]), 2 * np.pi, rtol=1e-12)
    npt.assert_allclose([r["cm_c4"] for r in rs], rs[0]["cm_c4"], atol=1e-15)   # c/4 is the a.c.
    assert impl.thin_airfoil(f, rs[0]["alpha_L0"], br)["cl"] == pytest.approx(0.0, abs=1e-13)


def test_centre_of_pressure_and_leading_edge_moment():
    f, br = impl.naca4_camber_slope(0.02, 0.4)
    prev = None
    for a in (1.0, 3.0, 6.0, 10.0):
        r = impl.thin_airfoil(f, np.radians(a), br)
        assert r["cm_le"] == pytest.approx(-r["x_cp"] * r["cl"])
        assert r["cm_le"] == pytest.approx(r["cm_c4"] - 0.25 * r["cl"])          # moment transfer
        if prev is not None:
            assert 0.25 < r["x_cp"] < prev                                       # moves forward with cl
        prev = r["x_cp"]


def test_quadrature_matches_adaptive_integration_and_split_is_exact():
    f, br = impl.naca4_camber_slope(0.02, 0.4)
    tp = np.arccos(1 - 2 * 0.4)
    g = lambda th, n: f(0.5 * (1 - np.cos(th))) * np.cos(n * th)
    for n in range(4):
        ref = integrate.quad(g, 0, tp, args=(n,), epsabs=1e-14)[0] + integrate.quad(g, tp, np.pi, args=(n,), epsabs=1e-14)[0]
        assert impl.glauert_integrals(f, br, 3, 16)[n] == pytest.approx(ref, abs=1e-12)  # smooth pieces
    unsplit = impl.glauert_integrals(f, [], 3, 16)
    split = impl.glauert_integrals(f, br, 3, 16)
    assert np.max(np.abs(unsplit - split)) > 1e-7                                # the kink costs accuracy


def test_naca2412_against_its_closed_form():
    # alpha_L0 of the 4-digit line in closed form (integrate dz/dx (cos th - 1) by hand on each piece)
    m, p = 0.02, 0.4
    tp = np.arccos(1 - 2 * p)
    F = lambda th, k: (2 * m / k) * ((p - 0.5) * (np.sin(th) - th) + 0.25 * (np.sin(th) * np.cos(th) + th) - 0.5 * np.sin(th) * 1)
    a1 = (F(tp, p ** 2) - F(0, p ** 2)) + (F(np.pi, (1 - p) ** 2) - F(tp, (1 - p) ** 2))
    f, br = impl.naca4_camber_slope(m, p)
    assert impl.thin_airfoil(f, 0.0, br)["alpha_L0"] == pytest.approx(-a1 / np.pi, rel=1e-10)
    assert np.degrees(impl.thin_airfoil(f, 0.0, br)["alpha_L0"]) == pytest.approx(-2.077, abs=0.001)


@pytest.mark.parametrize("E", [0.1, 0.25, 0.4, 0.7])
def test_flap_effectiveness_and_moment_closed_forms(E):
    d = np.radians(5.0)
    f, br = impl.flap_camber_slope(E, d)
    r0 = impl.thin_airfoil(f, 0.0, br)
    assert -r0["alpha_L0"] / d == pytest.approx(impl.flap_effectiveness(E), rel=1e-12)
    assert r0["cm_c4"] / d == pytest.approx(impl.flap_cm_derivative(E), rel=1e-12)
    assert r0["cl"] == pytest.approx(2 * np.pi * impl.flap_effectiveness(E) * d, rel=1e-12)


def test_flap_limits():
    assert impl.flap_effectiveness(1.0) == pytest.approx(1.0)                    # whole plate rotates
    assert impl.flap_effectiveness(1e-4) == pytest.approx(4 / np.pi * np.sqrt(1e-4), rel=1e-2)  # tiny tab
    assert impl.flap_cm_derivative(0.5) == pytest.approx(-0.5)
    Es = np.linspace(0.05, 0.95, 19)
    assert np.all(np.diff([impl.flap_effectiveness(e) for e in Es]) > 0)


def test_camber_and_flap_superpose():
    fc, bc = impl.naca4_camber_slope(0.02, 0.4)
    ff, bf = impl.flap_camber_slope(0.25, np.radians(10.0))
    fs, bs = impl.combine_slopes((fc, bc), (ff, bf))
    a = np.radians(2.0)
    r, rc, rf = impl.thin_airfoil(fs, a, bs), impl.thin_airfoil(fc, a, bc), impl.thin_airfoil(ff, 0.0, bf)
    assert r["cl"] == pytest.approx(rc["cl"] + rf["cl"], rel=1e-12)
    assert r["cm_c4"] == pytest.approx(rc["cm_c4"] + rf["cm_c4"], rel=1e-12)


def test_naca230_camber_has_its_design_ideal_cl():
    f, br = impl.naca5_camber_slope(0.2025, 15.957)
    r = impl.thin_airfoil(f, 0.0, br)
    assert np.pi * r["A"][1] == pytest.approx(0.30, abs=5e-4)                   # ideal cl = pi A1
    a_ideal = impl.glauert_integrals(f, br, 0)[0] / np.pi                         # A0 = 0
    assert impl.fourier_coefficients(f, a_ideal, br)[0] == pytest.approx(0.0, abs=1e-15)


def test_lumped_vortex_is_exact_for_plate_and_parabola_and_converges_at_order_two():
    a = np.radians(4.0)
    assert impl.lumped_vortex(FLAT, a, 5)["cl"] == pytest.approx(2 * np.pi * a, rel=1e-12)
    par = lambda x: 0.12 * (1 - 2 * np.asarray(x))
    assert impl.lumped_vortex(par, a, 6)["cl"] == pytest.approx(impl.thin_airfoil(par, a)["cl"], rel=1e-12)
    f, br = impl.naca4_camber_slope(0.02, 0.4)
    ex = impl.thin_airfoil(f, a, br)
    e_cl = [abs(impl.lumped_vortex(f, a, n)["cl"] - ex["cl"]) for n in (20, 40, 80, 160)]
    e_cm = [abs(impl.lumped_vortex(f, a, n)["cm_c4"] - ex["cm_c4"]) for n in (20, 40, 80, 160)]
    assert np.all(np.abs(impl.observed_order(e_cl) - 2) < 0.1)
    assert np.all(np.abs(impl.observed_order(e_cm) - 2) < 0.1)


def test_panel_method_tends_to_thin_airfoil_as_thickness_vanishes():
    a = np.radians(4.0)
    f, br = impl.naca4_camber_slope(0.02, 0.4)
    ex = impl.thin_airfoil(f, a, br)
    d = []
    for t in (0.12, 0.06, 0.03, 0.015):
        cl, cm = impl.panel_cl_cm(*impl.naca4_coords(0.02, 0.4, t, 100), a)
        d.append((cl - ex["cl"]) / ex["cl"])
    d = np.array(d)
    assert np.all(d > 0) and np.all(np.diff(d) < 0) and d[-1] < 0.015
    npt.assert_allclose(d[:-1] / d[1:], 2.0, rtol=0.15)                         # difference ~ proportional to t


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "thin_airfoil.png").stat().st_size > 1000
