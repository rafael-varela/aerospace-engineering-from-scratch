"""Public tests for P10.14 Low Reynolds Number Aerodynamics. They assert physics, not magic numbers:
Thwaites' method reproduces its closed forms for Falkner-Skan (power-law) flows, the flat plate
within 1 % of Blasius, the stagnation-point value lambda = 0.075 and Howarth's separation point within
3 %; the laminar separation point does not depend on Re while the boundary layer thins as 1/sqrt(Re);
Michel's flat-plate transition Reynolds number does not depend on the speed; the bubble length scales
as 1/Re and its closure switches at Re_min; friction laws meet at the transition point; the authentic
NACA 0012 data load completely and show the effect of the grit strip.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import pytest

from conftest import load_impl

impl = load_impl(__file__)

CL = 0.445


# ---------------------------------------------------------------- Thwaites

@pytest.mark.parametrize("m", [0.0, 0.1, 0.3, 1.0])
def test_thwaites_power_law_closed_form(m):
    # Ue = s^m: theta^2 = 0.45 nu s^(1 - m) / (5 m + 1), lambda = 0.45 m / (5 m + 1)
    nu = 1e-5
    s = np.linspace(0, 1, 20001)
    ue = s ** m if m > 0 else np.ones_like(s)
    th, lam = impl.thwaites(s, ue, nu)
    exact = np.sqrt(0.45 * nu * s[1000:] ** (1 - m) / (5 * m + 1))
    np.testing.assert_allclose(th[1000:], exact, rtol=2e-3)
    np.testing.assert_allclose(lam[1000:-1], 0.45 * m / (5 * m + 1), rtol=1e-2, atol=1e-6)


def test_flat_plate_close_to_blasius_and_scales_with_nu():
    s = np.linspace(0, 2, 4001)
    th1, _ = impl.thwaites(s, np.ones_like(s), 1e-6)
    th2, _ = impl.thwaites(s, np.ones_like(s), 4e-6)
    assert th1[-1] / np.sqrt(1e-6 * 2) == pytest.approx(0.664, rel=0.012)     # Thwaites 0.671 vs Blasius 0.664
    np.testing.assert_allclose(th2[1:], 2 * th1[1:], rtol=1e-12)              # theta ~ sqrt(nu)


def test_stagnation_point_limit():
    s = np.linspace(0, 0.1, 1001)
    th, lam = impl.thwaites(s, 3.0 * s, 1e-5)
    np.testing.assert_allclose(lam, 0.075, rtol=1e-3)
    np.testing.assert_allclose(th, np.sqrt(0.075 * 1e-5 / 3.0), rtol=1e-3)    # constant thickness


def test_howarth_separation_and_reynolds_independence():
    x = np.linspace(0, 1, 4001)
    seps = []
    for nu in (1e-7, 1e-5, 1e-3):
        _, lam = impl.thwaites(x, 1 - x, nu)
        seps.append(impl.laminar_separation(x, lam))
    assert seps[0] == pytest.approx(0.120, rel=0.03)
    np.testing.assert_allclose(seps, seps[0], rtol=1e-9)
    _, lam = impl.thwaites(x, 1 + x, 1e-5)
    assert np.isnan(impl.laminar_separation(x, lam))                           # favourable: no separation


def test_michel_flat_plate_transition_reynolds_number_is_speed_independent():
    s = np.linspace(0, 5, 20001)
    nu = 1.5e-5
    Rex = []
    for U in (10.0, 20.0, 40.0):
        ue = np.full_like(s, U)
        th, _ = impl.thwaites(s, ue, nu)
        s_tr = impl.michel_transition(s, ue, th, nu)
        Rex.append(U * s_tr / nu)
    np.testing.assert_allclose(Rex, Rex[0], rtol=2e-3)
    assert 5e5 < Rex[0] < 5e6
    ue = np.full_like(s, 1.0)                                                  # very slow: never transitions
    th, _ = impl.thwaites(s, ue, nu)
    assert np.isnan(impl.michel_transition(s, ue, th, nu))


# ---------------------------------------------------------------- bubbles on airfoils

def test_bubble_length_scaling():
    assert impl.bubble_length(1.2, 2e-6) == pytest.approx(2 * impl.bubble_length(1.2, 1e-6))
    assert impl.bubble_length(1.2, 1e-6) == pytest.approx(impl.BUBBLE_RE * 1e-6 / 1.2)


def test_separation_point_is_reynolds_independent_and_bubble_shrinks():
    a = impl.bubble_analysis("0012", CL, 1e5)
    b = impl.bubble_analysis("0012", CL, 4e5)
    assert a["x_sep"] == pytest.approx(b["x_sep"], rel=1e-9)
    assert a["length"] == pytest.approx(4 * b["length"], rel=1e-9)
    assert a["theta_sep_over_c"] == pytest.approx(2 * b["theta_sep_over_c"], rel=1e-9)
    assert 0 < a["x_sep"] < 1


def test_closure_switches_at_re_min():
    r = impl.bubble_analysis("0012", CL, 1.64e5)
    lo = impl.bubble_analysis("0012", CL, 0.95 * r["Re_min"])
    hi = impl.bubble_analysis("0012", CL, 1.05 * r["Re_min"])
    assert not lo["closes"] and hi["closes"]


def test_thin_cambered_plate_leaves_less_recovery_after_separation():
    t = impl.bubble_analysis("0012", CL, 5e4)
    p = impl.bubble_analysis("4402", CL, 5e4)
    assert p["x_sep"] > t["x_sep"]                                             # separates later
    assert p["dCp_rec"] < 0.5 * t["dCp_rec"]                                   # and little recovery is left
    assert p["ue_peak"] < t["ue_peak"]                                         # gentler suction peak


def test_natural_transition_reynolds_number_separates_regimes():
    Re_nat = impl.natural_transition_re("0012", CL)
    assert 1e5 < Re_nat < 1e7
    assert impl.bubble_analysis("0012", CL, 0.8 * Re_nat)["bubble"]
    assert not impl.bubble_analysis("0012", CL, 1.25 * Re_nat)["bubble"]


# ---------------------------------------------------------------- friction drag vs Re

def test_friction_laws_and_mixed_flow():
    assert float(impl.cf_laminar(1e6)) == pytest.approx(1.328e-3)
    assert float(impl.cf_turbulent(1e7)) == pytest.approx(0.455 / 7 ** 2.58)
    Re_tr = 5e5
    # continuous at the transition point
    assert float(impl.cf_mixed(Re_tr * (1 + 1e-9), Re_tr)) == pytest.approx(float(impl.cf_laminar(Re_tr)), rel=1e-6)
    Re = np.array([1e6, 3e6, 1e7, 1e8])
    deficit = Re * (impl.cf_turbulent(Re) - impl.cf_mixed(Re, Re_tr))
    np.testing.assert_allclose(deficit, deficit[0], rtol=1e-10)                # a constant (about 1700 here)
    assert 1500 < deficit[0] < 1900
    assert float(impl.section_cd_friction(1e6, t_c=0.0, Re_tr=None)) == pytest.approx(2 * float(impl.cf_turbulent(1e6)))
    assert impl.reynolds_number(12.0, 0.2) == pytest.approx(12.0 * 0.2 / 1.461e-5)


def test_laminar_drag_rises_steeply_at_low_re():
    cd = impl.section_cd_friction(np.array([2e4, 8e4]), Re_tr=1e12)
    assert cd[0] / cd[1] == pytest.approx(2.0, rel=1e-10)                      # 1/sqrt(Re)


# ---------------------------------------------------------------- authentic data

def test_ladson_data_loads_and_grit_effects():
    rows = impl.read_ladson()
    assert len(rows) == 161
    blocks = impl.ladson_blocks(rows)
    assert len([k for k in blocks if k[0] == "I"]) == 6 and len([k for k in blocks if k[0] == "VII"]) == 5
    ReI, cdI = impl.cd_at_zero_lift(blocks, "I")
    ReV, cdV = impl.cd_at_zero_lift(blocks, "VII")
    assert np.all(np.diff(ReI) > 0) and np.all(cdI > 0.005) and np.all(cdI < 0.008)
    assert np.all(cdV[:4] > cdI[:4])                                           # grit costs drag at every Re
    _, clI = impl.clmax_vs_re(blocks, "I")
    _, clV = impl.clmax_vs_re(blocks, "VII")
    assert np.all(clV[:4] < clI[:4])                                           # and maximum lift


def test_panel_helper_sanity():
    X, Y = impl.naca4("0012")
    assert abs(impl.hess_smith(X, Y, 0.0)["cl"]) < 1e-10
    a = impl.alpha_for_cl("0012", 0.5)
    assert np.degrees(a) == pytest.approx(0.5 / (2 * np.pi * 1.09) * 180 / np.pi, rel=0.1)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "low_re.png").stat().st_size > 10_000
