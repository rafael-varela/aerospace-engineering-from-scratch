"""Public tests for P08.11 Integral Boundary Layer Methods. They assert physics, not magic numbers:
the closure fits are continuous and pass through l = 0 at lambda = -0.09; Thwaites's quadrature
reproduces its own closed forms (flat plate, wedge flows, the Hiemenz stagnation point, Howarth's
decelerating flow) and Blasius (published constants, labelled) to the method's known accuracy; the
momentum integral equation is satisfied by Thwaites's output; the exact Falkner-Skan solutions
from solve_bvp reduce to Blasius and Hiemenz and lose their wall shear near beta = -0.1988; an
independent marching solution of the boundary-layer PDEs agrees on separation; and separation on a
cylinder does not depend on speed, size or viscosity.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

BLASIUS_FPP0 = 0.332057            # published (Blasius flat plate; derived in P08.10)
BLASIUS_THETA = 0.664115
HOWARTH_XSEP = 0.1199              # published separation point of Ue = U0 (1 - x/L)
HIEMENZ_FPP0 = 1.23259             # published stagnation-point wall curvature (beta = 1)

NU, U = 1.863e-5, 60.0


# ---------------------------------------------------------------- closure

def test_closure_fits_are_continuous_and_separate_at_minus_0_09():
    l0m, H0m = impl.thwaites_l_H(-1e-12)
    l0p, H0p = impl.thwaites_l_H(0.0)
    assert l0m == pytest.approx(l0p, abs=1e-6) and l0p == pytest.approx(0.22)
    assert abs(H0m - H0p) < 0.01                               # the two fits meet (to their own accuracy)
    l_sep, H_sep = impl.thwaites_l_H(-0.09)
    assert abs(l_sep) < 2e-3 and 3.4 < H_sep < 3.8
    lam = np.linspace(-0.09, 0.1, 50)
    l, H = impl.thwaites_l_H(lam)
    assert np.all(np.diff(l) > 0) and np.all(np.diff(H) < 0)   # more favourable: more shear, fuller profile


# ---------------------------------------------------------------- Thwaites's quadrature

def test_flat_plate_against_blasius():
    x = np.linspace(0.0, 1.5, 1501)
    r = impl.thwaites(x, np.full_like(x, U), np.zeros_like(x), NU)
    Re = U * x[1:] / NU
    npt.assert_allclose(r["theta"][1:] * np.sqrt(Re) / x[1:], np.sqrt(0.45), rtol=1e-12)    # its own closed form
    assert np.sqrt(0.45) == pytest.approx(BLASIUS_THETA, rel=0.011)                       # 1 % from Blasius
    npt.assert_allclose(r["cf"][1:] * np.sqrt(Re), 2 * BLASIUS_FPP0, rtol=0.015)
    npt.assert_allclose(r["H"], 2.61) and np.all(r["lam"] == 0.0)


def test_output_satisfies_thwaites_ode_and_the_momentum_integral():
    # Thwaites's ODE, Ue d(theta^2)/dx + 6 theta^2 Ue' = 0.45 nu, holds to discretisation error;
    # von Karman's momentum integral d(theta)/dx + (2 + H) theta Ue'/Ue = c_f/2 holds to the accuracy
    # of the linear fit F(lambda) = 0.45 - 6 lambda, a few per cent away from separation
    x = np.linspace(0.0, 0.1, 2001)
    Ue, dU = U * (1 - x), np.full_like(x, -U)
    r = impl.thwaites(x, Ue, dU, NU)
    th2 = r["theta"] ** 2
    ode = Ue * np.gradient(th2, x) + 6 * th2 * dU
    npt.assert_allclose(ode[20:-1], 0.45 * NU, rtol=1e-3)
    dth = np.gradient(r["theta"], x)
    lhs = dth + (2 + r["H"]) * r["theta"] * dU / Ue
    ok = (x > 0.005) & (r["lam"] > -0.05)
    assert ok.sum() > 500
    npt.assert_allclose(lhs[ok], r["cf"][ok] / 2, rtol=0.06)


def test_wedge_flows_have_constant_lambda_equal_to_the_closed_form():
    for m in (0.1, 0.5, 1.0):
        x = np.linspace(0.0, 0.5, 2001)[1:]
        Ue, dU = 10.0 * x ** m, 10.0 * m * x ** (m - 1)
        theta0 = np.sqrt(0.45 * NU * x[0] / ((5 * m + 1) * Ue[0]))           # start on the similar solution
        r = impl.thwaites(x, Ue, dU, NU, theta0=theta0)
        w = impl.wedge_thwaites(m)
        npt.assert_allclose(r["lam"][200:], w["lam"], rtol=2e-3)
        assert w["lam"] == pytest.approx(0.45 * m / (5 * m + 1))


def test_stagnation_point_limit_and_hiemenz():
    a = 400.0
    x = np.linspace(0.0, 0.01, 201)
    r = impl.thwaites(x, a * x, np.full_like(x, a), NU)
    assert r["theta"][0] == pytest.approx(np.sqrt(0.075 * NU / a), rel=1e-12)   # the limit at x = 0
    # theta is constant at a stagnation point; the trapezoid rule for int x^5 dx is poor in the first
    # few cells (it gives h^6/2 instead of h^6/6 on the first), so check from the 20th station on
    npt.assert_allclose(r["theta"][20:], np.sqrt(0.075 * NU / a), rtol=5e-3)
    npt.assert_allclose(r["lam"][20:], 0.075, rtol=1e-2)
    hz = impl.falkner_skan(1.0)
    assert hz["fpp0"] == pytest.approx(HIEMENZ_FPP0, abs=1e-4)
    # Thwaites under-predicts Hiemenz's momentum thickness by about 6 %
    assert np.sqrt(0.075) / hz["theta"] == pytest.approx(0.937, abs=0.01)


def test_howarth_separation_closed_form_and_published():
    # Thwaites for Ue = U (1 - x/L): lambda = -0.45 (1 - (1 - x/L)^6) / (6 (1 - x/L)^6)
    x = np.linspace(0.0, 0.2, 2001)
    r = impl.thwaites(x, U * (1 - x), np.full_like(x, -U), NU)
    xi = 1 - x
    npt.assert_allclose(r["lam"], -0.45 * (1 - xi ** 6) / (6 * xi ** 6), rtol=1e-5, atol=1e-12)
    xs = impl.separation_point(x, r["lam"])
    assert -0.45 * (1 - (1 - xs) ** 6) / (6 * (1 - xs) ** 6) == pytest.approx(-0.09, abs=1e-6)
    assert xs == pytest.approx(HOWARTH_XSEP, rel=0.04)                     # within 4 % of the exact value


def test_separation_point_interpolates_and_reports_none():
    x = np.array([0.0, 1.0, 2.0, 3.0])
    assert impl.separation_point(x, np.array([0.0, -0.05, -0.13, -0.2])) == pytest.approx(1.5)
    assert np.isnan(impl.separation_point(x, np.zeros(4)))


def test_cylinder_separation_is_independent_of_speed_size_and_viscosity():
    angles = []
    for Uc, R, nu in ((20.0, 0.05, 1.5e-5), (3.0, 0.004, 1.0e-6), (60.0, 1.0, NU)):
        x = np.linspace(0.0, 0.75 * np.pi * R, 3001)
        Ue, dU = impl.cylinder_edge_velocity(x, Uc, R)
        r = impl.thwaites(x, Ue, dU, nu)
        angles.append(np.degrees(impl.separation_point(x, r["lam"]) / R))
    npt.assert_allclose(angles, angles[0], rtol=1e-6)
    assert 90.0 < angles[0] < 110.0                                        # behind the shoulder, not at it


def test_no_pressure_gradient_means_no_separation():
    x = np.linspace(0.0, 3.0, 301)
    r = impl.thwaites(x, np.full_like(x, U), np.zeros_like(x), NU)
    assert np.isnan(impl.separation_point(x, r["lam"]))


def test_thwaites_is_unit_invariant():
    x = np.linspace(0.0, 0.1, 401)
    r1 = impl.thwaites(x, U * (1 - x), np.full_like(x, -U), NU)
    r2 = impl.thwaites(x * 1e3, U * 1e3 * (1 - x), np.full_like(x, -U), NU * 1e6)   # millimetres
    npt.assert_allclose(r2["lam"], r1["lam"], rtol=1e-9, atol=1e-15)
    npt.assert_allclose(r2["theta"][1:], r1["theta"][1:] * 1e3, rtol=1e-9)
    npt.assert_allclose(r2["cf"][1:], r1["cf"][1:], rtol=1e-9)


# ---------------------------------------------------------------- exact similarity solutions

def test_falkner_skan_reduces_to_blasius():
    fs = impl.falkner_skan(0.0)
    # with eta = y sqrt(Ue/(2 nu x)) the Blasius wall curvature is 0.332057 sqrt(2)
    assert fs["fpp0"] == pytest.approx(BLASIUS_FPP0 * np.sqrt(2), abs=2e-5)
    assert fs["H"] == pytest.approx(2.591, abs=2e-3)
    assert fs["lam"] == 0.0 and fs["l"] == pytest.approx(BLASIUS_FPP0 * BLASIUS_THETA, rel=1e-3)


def test_falkner_skan_wall_shear_vanishes_near_beta_minus_0_1988():
    a, b = impl.falkner_skan(-0.15), impl.falkner_skan(-0.198)
    assert a["fpp0"] > b["fpp0"] > 0 and b["fpp0"] < 0.03
    assert b["H"] > 3.5                                                   # separating profiles are thick and empty
    assert -0.075 < b["lam"] < -0.06                                      # exact separation lambda ~ -0.068, not -0.09


def test_thwaites_closure_tracks_exact_wedge_flows():
    for beta in (-0.1, 0.0, 0.3, 1.0):
        fs = impl.falkner_skan(beta)
        l, H = impl.thwaites_l_H(fs["lam"])
        assert l == pytest.approx(fs["l"], rel=0.08) and H == pytest.approx(fs["H"], rel=0.07)
        w = impl.wedge_thwaites(fs["m"])
        assert w["theta_eta"] == pytest.approx(fs["theta"], rel=0.07)


# ---------------------------------------------------------------- independent check: marching the PDEs

def test_marched_pde_matches_blasius_shear_on_a_flat_plate():
    x = np.geomspace(1e-4, 1.0, 801)
    y_max = 15 * np.sqrt(NU * 1.0 / U)
    xd, y, tau = impl.march_boundary_layer(x, np.full_like(x, U), np.zeros_like(x), NU, y_max, ny=200)
    assert xd.size == x.size                                              # no separation on a flat plate
    cf = 2 * tau[-1] / U ** 2
    assert cf * np.sqrt(U * xd[-1] / NU) == pytest.approx(2 * BLASIUS_FPP0, rel=0.01)


def test_marched_pde_finds_howarth_separation():
    x = np.linspace(1e-4, 0.2, 2001)
    y_max = 20 * np.sqrt(NU * 0.12 / U)
    xd, y, tau = impl.march_boundary_layer(x, U * (1 - x), np.full_like(x, -U), NU, y_max, ny=200)
    assert xd[-1] < 0.2 and np.all(tau > 0) and np.all(np.diff(tau[len(tau) // 2:]) < 0)
    xs = impl.march_separation(xd, tau)
    assert xs == pytest.approx(HOWARTH_XSEP, rel=0.005)                   # published, to 0.5 %
    # and Thwaites (an integral closure) lands within 3 % of the marched PDE
    r = impl.thwaites(x, U * (1 - x), np.full_like(x, -U), NU)
    assert impl.separation_point(x, r["lam"]) == pytest.approx(xs, rel=0.03)


def test_march_separation_extrapolates_a_square_root():
    xs = 0.37
    x = np.linspace(0.2, 0.36, 40)
    assert impl.march_separation(x, 3.0 * np.sqrt(xs - x)) == pytest.approx(xs, rel=1e-10)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "thwaites.png").stat().st_size > 10_000
