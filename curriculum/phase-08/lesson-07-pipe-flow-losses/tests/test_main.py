"""Public tests for P08.07 Pipe Flow & Losses. They assert fluid mechanics and numerics, not magic
numbers: the laminar friction factor follows from the exact Hagen-Poiseuille profile; Newton's
method solves Colebrook to round-off, converges quadratically and agrees with the explicit Haaland
formula to its published accuracy; Colebrook reduces to the smooth-pipe law and the fully rough law
in its limits; the Darcy-Weisbach drop reproduces Hagen-Poiseuille in laminar flow; the pressure
budget adds up, conserves mass through segments of different size, and is unchanged by a change of
units; and the line-sizing root-find meets its target.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import pytest
from scipy.integrate import quad

from conftest import load_impl

impl = load_impl(__file__)

RES = np.logspace(np.log10(4e3), 8, 25)
RELS = [0.0, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2]


# ---------------------------------------------------------------- laminar

def test_laminar_friction_factor_from_the_poiseuille_profile():
    # f = 8 tau_w / (rho V^2) with tau_w and V from the exact profile w = G (R^2 - r^2) / (4 mu)
    rho, mu, R, G = 900.0, 0.05, 4e-3, 2e3
    w = lambda r: G * (R ** 2 - r ** 2) / (4 * mu)
    V = quad(lambda r: 2 * np.pi * r * w(r), 0, R)[0] / (np.pi * R ** 2)
    tau_w = G * R / 2                                  # mu |dw/dr| at r = R, = force balance G pi R^2 = tau 2 pi R
    f = 8 * tau_w / (rho * V ** 2)
    Re = impl.reynolds(rho, V, 2 * R, mu)
    assert f == pytest.approx(impl.friction_laminar(Re), rel=1e-10)
    assert f * Re == pytest.approx(64.0, rel=1e-10)


def test_darcy_weisbach_reproduces_hagen_poiseuille_in_laminar_flow():
    rho, mu, D, L, Q = 1000.0, 1e-3, 6e-3, 2.0, 2e-6
    V = Q / (np.pi * D ** 2 / 4)
    Re = impl.reynolds(rho, V, D, mu)
    assert Re < 2300
    dp = impl.darcy_weisbach_dp(impl.friction_factor(Re, 1e-3), L, D, rho, V)
    assert dp == pytest.approx(128 * mu * L * Q / (np.pi * D ** 4), rel=1e-12)   # roughness does not matter


# ---------------------------------------------------------------- Colebrook

def test_newton_solves_colebrook_to_round_off_across_the_moody_chart():
    for e in RELS:
        for Re in RES:
            f, hist = impl.friction_colebrook(Re, e)
            assert abs(impl.colebrook_residual(f, Re, e)) < 1e-12
            assert 0.005 < f < 0.08
            assert hist.size - 1 <= 6                    # from the Haaland start


def test_newton_converges_quadratically_from_a_poor_start():
    Re, e = 2e5, 1e-4
    f, hist = impl.friction_colebrook(Re, e, x0=1.0)
    x_star = 1 / np.sqrt(f)
    err = np.abs(hist - x_star)
    err = err[err > 1e-12]
    # e_{k+1} / e_k^2 stays bounded (quadratic), so the number of correct digits roughly doubles
    ratios = err[1:] / err[:-1] ** 2
    assert np.all(ratios < 1.0)
    assert err[-1] < 1e-6 and hist.size - 1 <= 10


def test_newton_iterates_after_the_first_lie_left_of_the_root_and_increase():
    # g(x) = x + 2 log10(a + b x) is increasing and concave, so each tangent zero has g <= 0:
    # from the first step on, the iterates sit at or left of the root and climb to it
    for Re, e, x0 in [(2e5, 1e-4, None), (2e5, 1e-4, 1.0), (2e5, 1e-4, 50.0), (1e7, 0.0, 30.0),
                      (4e3, 0.05, 20.0)]:
        f, hist = impl.friction_colebrook(Re, e, x0=x0)
        a, b = e / 3.7, 2.51 / Re
        x_star = 1 / np.sqrt(f)
        assert abs(x_star + 2 * np.log10(a + b * x_star)) < 1e-12
        later = hist[1:]
        assert np.all(later <= x_star + 1e-12 * x_star)
        assert np.all(np.diff(later) >= -1e-12 * x_star)
        if x0 is not None and x0 > x_star:
            assert later[0] < x_star                     # thrown to the left by the first step


def test_haaland_agrees_with_colebrook_to_its_published_accuracy():
    worst = max(abs(impl.friction_haaland(Re, e) / impl.friction_colebrook(Re, e)[0] - 1)
                for Re in RES for e in RELS)
    assert worst < 0.02                                  # Haaland 1983: within about 1.5 %


def test_colebrook_limits_smooth_and_fully_rough():
    # smooth pipe: Colebrook with eps = 0 is Prandtl's law 1/sqrt f = 2 log10(Re sqrt f) - 0.8 (2 log10 2.51)
    for Re in (1e4, 1e6, 1e8):
        f = impl.friction_colebrook(Re, 0.0)[0]
        assert 1 / np.sqrt(f) == pytest.approx(2 * np.log10(Re * np.sqrt(f)) - 2 * np.log10(2.51), rel=1e-12)
    # fully rough: as Re -> infinity, 1/sqrt f -> -2 log10(eps / 3.7), independent of Re
    for e in (1e-3, 1e-2):
        f_inf = (2 * np.log10(3.7 / e)) ** -2
        assert impl.friction_colebrook(1e12, e)[0] == pytest.approx(f_inf, rel=1e-4)


def test_blasius_and_monotonic_trends():
    for Re in (5e3, 2e4, 1e5):
        assert impl.friction_blasius(Re) == pytest.approx(impl.friction_colebrook(Re, 0.0)[0], rel=0.03)
    f_re = [impl.friction_colebrook(Re, 0.0)[0] for Re in RES]
    assert np.all(np.diff(f_re) < 0)                     # smooth: f falls with Re
    f_e = [impl.friction_colebrook(1e6, e)[0] for e in RELS]
    assert np.all(np.diff(f_e) > 0)                      # rougher: larger f


def test_friction_factor_regimes_are_continuous():
    e = 1e-4
    assert impl.friction_factor(1000.0, e) == pytest.approx(0.064)
    for Re0 in (impl.RE_LAMINAR, impl.RE_TURBULENT):
        assert impl.friction_factor(Re0 * (1 - 1e-9), e) == pytest.approx(impl.friction_factor(Re0 * (1 + 1e-9), e), rel=1e-6)
    assert impl.friction_factor(1e5, e) == pytest.approx(impl.friction_colebrook(1e5, e)[0])
    with pytest.raises(ValueError):
        impl.friction_factor(0.0, e)


# ---------------------------------------------------------------- losses and the budget

def test_loss_scaling_with_length_and_flow():
    rho, mu, D, eps = impl.RHO_RP1, impl.MU_RP1, 0.02, impl.ROUGHNESS_DRAWN
    a = impl.segment_loss(1.0, rho, mu, D, 2.0, eps, K=3.0)
    b = impl.segment_loss(1.0, rho, mu, D, 4.0, eps, K=3.0)
    assert b["dp_friction"] == pytest.approx(2 * a["dp_friction"])          # proportional to L
    assert b["dp_minor"] == pytest.approx(a["dp_minor"])
    c = impl.segment_loss(2.0, rho, mu, D, 2.0, eps, K=3.0)
    assert c["dp_minor"] == pytest.approx(4 * a["dp_minor"])                 # K rho V^2 / 2
    # turbulent friction grows a little slower than V^2 (f falls with Re), a lot faster than V
    ratio = c["dp_friction"] / a["dp_friction"]
    assert 3.0 < ratio < 4.0


def test_injector_orifice_inverts():
    A = impl.injector_area_for_dp(1.6, 810.0, 0.7, 0.4e6)
    assert impl.injector_dp(1.6, 810.0, 0.7, A) == pytest.approx(0.4e6, rel=1e-12)


def test_budget_adds_up_and_conserves_mass_through_segments():
    rho, mu, eps = impl.RHO_RP1, impl.MU_RP1, impl.ROUGHNESS_DRAWN
    segs = [{"D": 0.025, "L": 1.0, "roughness": eps, "K": 0.5},
            {"D": 0.016, "L": 1.5, "roughness": eps, "K": 1.3}]
    A_inj = impl.injector_area_for_dp(1.6, rho, 0.7, 0.4e6)
    b = impl.feed_line_budget(1.6, rho, mu, segs, 2.0e6, 0.7, A_inj, head=1.0, accel_g=3.0)
    total = b["p_chamber"] + b["dp_injector"] + b["dp_friction"] + b["dp_minor"] - b["dp_head"]
    assert b["p_tank"] == pytest.approx(total, rel=1e-14)
    s1, s2 = b["segments"]
    assert rho * s1["V"] * np.pi * 0.025 ** 2 / 4 == pytest.approx(1.6)      # mass flow in each segment
    assert rho * s2["V"] * np.pi * 0.016 ** 2 / 4 == pytest.approx(1.6)
    assert s1["V"] * 0.025 ** 2 == pytest.approx(s2["V"] * 0.016 ** 2)
    assert b["dp_head"] == pytest.approx(rho * 3.0 * impl.G0 * 1.0)


def test_budget_at_zero_flow_is_hydrostatic():
    rho, mu = impl.RHO_RP1, impl.MU_RP1
    b = impl.feed_line_budget(1e-12, rho, mu, impl.fuel_line(), 2.0e6, 0.7, 1e-4, head=1.2, accel_g=2.0)
    assert b["p_tank"] == pytest.approx(2.0e6 - rho * 2.0 * impl.G0 * 1.2, rel=1e-9)


def test_results_do_not_depend_on_the_units():
    # the same line in mm, g and ms: Re and f are unchanged, and dp in g/(mm ms^2) is dp in Pa / 1e6
    rho, mu, D, L, eps, mdot = 810.0, 2.0e-3, 0.02, 2.5, 1.5e-6, 1.6
    si = impl.segment_loss(mdot, rho, mu, D, L, eps, K=6.0)
    mm = impl.segment_loss(mdot * 1e3 / 1e3, rho * 1e3 / 1e9, mu * 1e3 / (1e3 * 1e3), D * 1e3, L * 1e3, eps * 1e3, K=6.0)
    assert mm["Re"] == pytest.approx(si["Re"], rel=1e-12)
    assert mm["f"] == pytest.approx(si["f"], rel=1e-12)
    assert mm["dp_friction"] * 1e6 == pytest.approx(si["dp_friction"], rel=1e-12)
    assert mm["V"] == pytest.approx(si["V"], rel=1e-12)                    # mm/ms = m/s


def test_line_diameter_meets_the_loss_target():
    rho, mu, eps = impl.RHO_RP1, impl.MU_RP1, impl.ROUGHNESS_DRAWN
    K = impl.fuel_line()[0]["K"]
    D = impl.line_diameter_for_loss(1.6, rho, mu, 2.5, eps, K, 0.2e6)
    s = impl.segment_loss(1.6, rho, mu, D, 2.5, eps, K)
    assert s["dp_friction"] + s["dp_minor"] == pytest.approx(0.2e6, rel=1e-8)
    s_small = impl.segment_loss(1.6, rho, mu, 0.95 * D, 2.5, eps, K)
    assert s_small["dp_friction"] + s_small["dp_minor"] > 0.2e6
    # loss scales roughly as D^-5 at fixed mdot (V^2 ~ D^-4, 1/D, a slowly varying f)
    big = impl.segment_loss(1.6, rho, mu, 2 * D, 2.5, eps, 0.0)
    ref = impl.segment_loss(1.6, rho, mu, D, 2.5, eps, 0.0)
    assert 4.5 < np.log2(ref["dp_friction"] / big["dp_friction"]) < 5.0


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "feed_line.png").stat().st_size > 10_000
