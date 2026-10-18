"""Public tests for P10.10 Prandtl Lifting-Line Theory. They assert the theory, not magic numbers:
the elliptic wing is reproduced exactly (A_n = 0 for n > 1, e = 1, the closed-form CL); the solution
is linear in alpha and zero at the zero-lift angle; a symmetric wing has no even coefficients and a
symmetric load that vanishes at the tips; CL and CDi from the coefficients equal the span integrals
of the load and of load times induced angle; collocation converges at the expected observed order
(4 for the rectangular wing, 2 when the chord or twist has a kink at the root); the induced-drag
factor has its minimum at a moderate taper; CL_alpha tends to a0 as AR grows; twist superposes;
an anti-symmetric twist rolls the wing the right way without changing CL.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

B, S, C = 10.0, 12.0, 1.2


def test_collocation_points():
    th = impl.collocation_points(8, True)
    assert len(th) == 8 and th[-1] == pytest.approx(np.pi / 2) and th[0] > 0
    th = impl.collocation_points(9, False)
    assert len(th) == 9 and 0 < th[0] and th[-1] < np.pi
    npt.assert_allclose(th + th[::-1], np.pi)                                 # symmetric about the root


@pytest.mark.parametrize("AR,a0,alpha_deg", [(8.3333, 2 * np.pi, 5.0), (4.0, 5.7, 3.0), (20.0, 2 * np.pi, -2.0)])
def test_elliptic_wing_matches_the_closed_form_exactly(AR, a0, alpha_deg):
    w = impl.elliptic_wing(B, B ** 2 / AR, a0=a0)
    sol = impl.solve_monoplane(w, np.radians(alpha_deg), 12)
    r = impl.wing_coefficients(w, sol)
    assert np.max(np.abs(sol["A"][1:])) < 1e-12
    assert r["e"] == pytest.approx(1.0, abs=1e-12)
    cf = impl.elliptic_closed_form(a0, AR, np.radians(alpha_deg))
    assert r["CL"] == pytest.approx(cf["CL"], rel=1e-6)
    assert r["CDi"] == pytest.approx(cf["CDi"], rel=1e-6)


def test_linear_in_alpha_and_zero_at_zero_lift_angle():
    aL0 = np.radians(-2.0)
    w = impl.rectangular_wing(B, C, alpha_L0=lambda y: np.full_like(np.asarray(y, dtype=float), aL0))
    cl = [impl.wing_coefficients(w, impl.solve_monoplane(w, a))["CL"] for a in (aL0, 0.0, 0.05, 0.1)]
    assert abs(cl[0]) < 1e-14
    slope = (cl[3] - cl[1]) / 0.1
    assert cl[2] == pytest.approx(cl[1] + slope * 0.05, rel=1e-12)
    assert cl[1] == pytest.approx(slope * (0.0 - aL0), rel=1e-12)


def test_symmetric_wing_has_no_even_terms_and_a_symmetric_load():
    w = impl.tapered_wing(B, S, 0.45, twist=impl.linear_twist(B, np.radians(-3.0)))
    full = impl.solve_monoplane(w, 0.08, 80, symmetric=False)
    npt.assert_allclose(full["A"][1::2], 0.0, atol=1e-12)
    half = impl.solve_monoplane(w, 0.08, 40, symmetric=True)
    assert impl.wing_coefficients(w, full)["CL"] == pytest.approx(impl.wing_coefficients(w, half)["CL"], rel=1e-3)
    y = np.linspace(-4.9, 4.9, 15)
    g = impl.span_loading(w, half, y)["gamma"]
    npt.assert_allclose(g, g[::-1], rtol=1e-12)
    assert np.all(np.abs(impl.span_loading(w, half, np.array([-B / 2, B / 2]))["gamma"]) < 1e-12)


def test_coefficients_equal_span_integrals_of_the_load():
    w = impl.tapered_wing(B, S, 0.3, twist=impl.linear_twist(B, np.radians(-2.0)))
    V = 40.0
    sol = impl.solve_monoplane(w, 0.07, 40)
    r = impl.wing_coefficients(w, sol)
    th = (np.arange(4000) + 0.5) * np.pi / 4000
    y = -0.5 * B * np.cos(th)
    dy = 0.5 * B * np.sin(th) * np.pi / 4000
    sl = impl.span_loading(w, sol, y, V)
    Sref = w.area()
    assert 2 / (V * Sref) * np.sum(sl["gamma"] * dy) == pytest.approx(r["CL"], rel=1e-10)          # Kutta-Joukowski
    assert 2 / (V * Sref) * np.sum(sl["gamma"] * sl["alpha_i"] * dy) == pytest.approx(r["CDi"], rel=1e-8)   # near-field drag


def test_collocation_converges_at_the_observed_orders():
    cs = impl.convergence_study(impl.rectangular_wing(B, C), np.radians(5.0), [10, 20, 40, 80], N_ref=640)
    assert np.all(cs["order_CL"] > 3.7) and np.all(cs["order_CL"] < 4.3)       # smooth chord: fast
    cs = impl.convergence_study(impl.tapered_wing(B, S, 0.4), np.radians(5.0), [10, 20, 40, 80], N_ref=640)
    assert np.all(cs["order_CL"] > 1.8) and np.all(cs["order_CL"] < 2.3)       # root kink in c(y): order 2
    assert cs["err_CL"][-1] < 1e-4


def test_induced_drag_factor_has_a_minimum_at_moderate_taper():
    tapers = np.linspace(0.0, 1.0, 41)
    d = [impl.wing_coefficients(impl.tapered_wing(B, S, t), impl.solve_monoplane(impl.tapered_wing(B, S, t), 0.1, 60))["delta"]
         for t in tapers]
    assert np.all(np.array(d) > 0)
    assert 0.3 <= tapers[int(np.argmin(d))] <= 0.45
    assert d[-1] > 4 * min(d) and d[0] > d[-1]                                # rectangle and triangle both worse
    # at fixed (rectangular) planform a larger AR departs further from elliptic
    dr = [impl.wing_coefficients(impl.rectangular_wing(B, B / AR), impl.solve_monoplane(impl.rectangular_wing(B, B / AR), 0.1, 60))["delta"]
          for AR in (4, 8, 16)]
    assert dr[0] < dr[1] < dr[2]


def test_lift_slope_limits_and_tau_definition():
    a0 = 2 * np.pi
    mid = impl.lift_slope_and_tau(impl.rectangular_wing(B, B / 200.0), 60)
    big = impl.lift_slope_and_tau(impl.rectangular_wing(B, B / 2000.0), 60)
    assert mid["CL_alpha"] < big["CL_alpha"] < a0 and big["CL_alpha"] / a0 > 0.995   # 2-D limit
    r = impl.lift_slope_and_tau(impl.rectangular_wing(B, C), 60)
    AR = B / C
    assert r["CL_alpha"] < a0 / (1 + a0 / (np.pi * AR))                       # below the elliptic slope
    assert r["CL_alpha"] == pytest.approx(a0 / (1 + a0 / (np.pi * AR) * (1 + r["tau"])), rel=1e-8)
    assert impl.lift_slope_and_tau(impl.elliptic_wing(B, S), 20)["tau"] == pytest.approx(0.0, abs=1e-10)


def test_twist_superposes_and_washout_unloads_the_tips():
    flat = impl.tapered_wing(B, S, 0.6)
    tw = impl.tapered_wing(B, S, 0.6, twist=impl.linear_twist(B, np.radians(-4.0)))
    a = 0.08
    s_flat, s_tw = impl.solve_monoplane(flat, a), impl.solve_monoplane(tw, a)
    s_tw0 = impl.solve_monoplane(tw, 0.0)
    s_flat0 = impl.solve_monoplane(flat, 0.0)
    npt.assert_allclose(s_tw["A"], s_flat["A"] + (s_tw0["A"] - s_flat0["A"]), atol=1e-14)   # linear superposition
    CL = 0.5
    a_f, a_t = impl.alpha_for_cl(flat, CL), impl.alpha_for_cl(tw, CL)
    assert a_t > a_f
    y = np.array([0.0, 4.5])
    clf = impl.span_loading(flat, impl.solve_monoplane(flat, a_f), y)["cl"]
    clt = impl.span_loading(tw, impl.solve_monoplane(tw, a_t), y)["cl"]
    assert clt[1] / clt[0] < clf[1] / clf[0]                                  # washout lowers tip cl / root cl
    assert impl.wing_coefficients(tw, impl.solve_monoplane(tw, a_t))["CL"] == pytest.approx(CL, rel=1e-10)


def test_antisymmetric_twist_rolls_without_changing_lift():
    flat = impl.rectangular_wing(B, C)
    ail = impl.rectangular_wing(B, C, twist=lambda y: np.radians(2.0) * np.sign(y) * (np.abs(y) > 3.0))
    s0 = impl.solve_monoplane(flat, 0.05, 40, symmetric=False)
    s1 = impl.solve_monoplane(ail, 0.05, 40, symmetric=False)
    assert impl.wing_coefficients(ail, s1)["CL"] == pytest.approx(impl.wing_coefficients(flat, s0)["CL"], rel=1e-10)
    assert impl.rolling_moment_coefficient(ail, s1) < 0                         # more lift on the right: right wing up
    assert impl.rolling_moment_coefficient(flat, s0) == pytest.approx(0.0, abs=1e-14)
    s2 = impl.solve_monoplane(impl.rectangular_wing(B, C, twist=lambda y: np.radians(4.0) * np.sign(y) * (np.abs(y) > 3.0)), 0.05, 40, symmetric=False)
    assert impl.rolling_moment_coefficient(ail, s2) == pytest.approx(2 * impl.rolling_moment_coefficient(ail, s1), rel=1e-10)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "lifting_line.png").stat().st_size > 10_000
