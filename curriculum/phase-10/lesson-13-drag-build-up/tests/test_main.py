"""Public tests for P10.13 Drag Build-Up. They assert physics and bookkeeping, not magic numbers:
the friction laws and their mixed form reduce to the pure laminar and turbulent limits; earlier
transition costs drag; form factors tend to 1 for thin sections and slender bodies; wetted areas
match closed forms; each component's drag scales with its wetted area; C_D0 adds up; lifting line
gives e within 2 % of 1 for a taper ratio of 0.4 and below 1 for a rectangle; the parabolic polar's best L/D
matches a brute-force maximum; level flight has L = W; and on the authentic NACA 0012 data the
fully turbulent estimate brackets the measurements from the side the physics predicts.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import pytest

from conftest import load_impl

impl = load_impl(__file__)


# ---------------------------------------------------------------- skin friction

def test_friction_laws():
    assert float(impl.cf_laminar(1e6)) == pytest.approx(1.328e-3)
    assert float(impl.cf_turbulent(1e6)) == pytest.approx(0.455 / 6 ** 2.58)
    Re = np.logspace(5.5, 8, 20)
    assert np.all(impl.cf_turbulent(Re) > impl.cf_laminar(Re))
    assert np.all(np.diff(impl.cf_turbulent(Re)) < 0)


def test_mixed_friction_limits_and_monotonic_in_transition():
    Re = 5e6
    assert float(impl.cf_mixed(Re, 0.0)) == pytest.approx(float(impl.cf_turbulent(Re)))
    assert float(impl.cf_mixed(Re, 1.0)) == pytest.approx(float(impl.cf_laminar(Re)))
    xs = np.linspace(0.02, 0.98, 25)
    cf = np.array([float(impl.cf_mixed(Re, x)) for x in xs])
    assert np.all(np.diff(cf) < 0)                                     # later transition, less friction
    # the laminar-run correction Re (Cf_t - Cf_mixed) depends only on the transition Reynolds number
    a = 5e6 * (impl.cf_turbulent(5e6) - impl.cf_mixed(5e6, 0.1))
    b = 2.5e6 * (impl.cf_turbulent(2.5e6) - impl.cf_mixed(2.5e6, 0.2))
    assert float(a) == pytest.approx(float(b), rel=1e-12)


# ---------------------------------------------------------------- form factors and areas

def test_form_factors():
    assert impl.form_factor_section(0.0) == 1.0
    assert impl.form_factor_section(0.12) == pytest.approx(1 + 0.24 + 60 * 0.12 ** 4)
    assert impl.form_factor_section(0.20) > impl.form_factor_section(0.12)
    assert impl.form_factor_body(4.0) > impl.form_factor_body(8.0)
    assert impl.form_factor_body(8.0) == pytest.approx(1 + 60 / 512 + 8 / 400)
    assert float(impl.section_cd(5e6, 0.0, 0.0)) == pytest.approx(2 * float(impl.cf_turbulent(5e6)))


def test_perimeter_and_body_wetted_area():
    assert impl.perimeter_ratio("0001") == pytest.approx(1.0, abs=2e-3)
    assert 1.0 < impl.perimeter_ratio("0012") < impl.perimeter_ratio("0018")
    x = np.linspace(0, 3.0, 301)
    assert impl.body_wetted_area(x, np.full_like(x, 0.4)) == pytest.approx(2 * np.pi * 0.4 * 3.0)   # cylinder
    t = np.linspace(0, np.pi, 4001)
    xs = 0.5 * (1 - np.cos(t))
    assert impl.body_wetted_area(xs, 0.5 * np.sin(t)) == pytest.approx(np.pi, rel=1e-5)             # sphere 4 pi R^2


def test_component_scaling():
    a = impl.component("x", 1.0, 10.0, 1.2, 1.0)
    b = impl.component("x", 1.0, 20.0, 1.2, 1.1)
    assert b["D_over_q"] == pytest.approx(2.2 * a["D_over_q"])
    assert a["CD"] == pytest.approx(a["D_over_q"] / 12.0)
    assert a["Re"] == pytest.approx(55.0 / 1.461e-5)
    fast = impl.component("x", 1.0, 10.0, 1.2, 1.0, V=110.0)
    assert fast["Cf"] < a["Cf"]                                             # higher Re, lower friction coefficient


def test_cd0_bookkeeping():
    comps = impl.trainer_components()
    assert [c["name"] for c in comps][:2] == ["wing", "fuselage"]
    r = impl.cd0_build_up(comps, gear_dq=0.0, cooling_dq=0.0, misc=0.0)
    assert r["CD0"] == pytest.approx(sum(c["CD"] for c in comps))
    r2 = impl.cd0_build_up(comps, gear_dq=0.12, cooling_dq=0.0, misc=0.0)
    assert r2["CD0"] - r["CD0"] == pytest.approx(0.12 / 12.0)
    r3 = impl.cd0_build_up(comps)
    assert r3["CD0"] == pytest.approx(r3["components"] + r3["gear"] + r3["cooling"] + r3["misc"])
    assert 0.015 < r3["CD0"] < 0.05                                        # a light aeroplane's range
    wing = comps[0]
    assert 2.0 * 10.68 < wing["S_wet"] < 2.1 * 10.68                        # about twice the exposed area


# ---------------------------------------------------------------- induced drag and polar

def test_lifting_line_span_efficiency():
    e_rect = impl.lifting_line_e(8.33, 1.0)
    e_tap = impl.lifting_line_e(8.33, 0.4)
    assert 0.9 < e_rect < 0.99
    assert e_tap > 0.98 and e_tap > e_rect and e_tap > impl.lifting_line_e(8.33, 0.1)
    # an elliptic wing is the e = 1 limit; no planform beats it
    assert e_tap < 1.0
    assert impl.lifting_line_e(4.0, 1.0) > e_rect > impl.lifting_line_e(16.0, 1.0)   # rectangle: worse at high AR


def test_oswald_estimate_limits():
    e_inv = impl.oswald_e(0.03, k_visc=0.0)
    assert e_inv == pytest.approx(impl.lifting_line_e(100 / 12) * (1 - 2 * 0.11 ** 2))
    assert impl.oswald_e(0.0) == pytest.approx(e_inv)
    assert impl.oswald_e(0.04) < impl.oswald_e(0.02) < e_inv


def test_best_ld_matches_brute_force():
    CD0, e = 0.03, 0.75
    b = impl.best_ld(CD0, e)
    CL = np.linspace(0.05, 2.0, 200001)
    LD = CL / impl.polar(CL, CD0, e)
    assert b["LD"] == pytest.approx(LD.max(), rel=1e-8)
    assert b["CL"] == pytest.approx(CL[np.argmax(LD)], abs=1e-4)
    # at best L/D the induced drag equals the zero-lift drag
    assert impl.polar(b["CL"], CD0, e) == pytest.approx(2 * CD0)
    assert 0.5 * 1.225 * b["V"] ** 2 * 12.0 * b["CL"] == pytest.approx(750 * 9.81)


def test_cruise_drag_lift_equals_weight():
    r = impl.cruise_drag(55.0, 0.03, 0.75)
    assert 0.5 * 1.225 * 55.0 ** 2 * 12.0 * r["CL"] == pytest.approx(750 * 9.81)
    assert r["P"] == pytest.approx(r["D"] * 55.0)
    assert r["CDi"] == pytest.approx(r["CL"] ** 2 / (np.pi * 0.75 * 100 / 12))


# ---------------------------------------------------------------- data-based investigation

def test_ladson_blocks_near_zero_lift():
    rows = impl.read_ladson()
    Re, al, cd = impl.ladson_cd_near_zero_lift(rows, "I")
    assert len(Re) == 6 and np.all(np.abs(al) < 0.2) and np.all(np.diff(Re) > 0)
    ReG, alG, cdG = impl.ladson_cd_near_zero_lift(rows, "VII")
    assert len(ReG) == 5 and np.all(np.abs(alG) < 0.2)


def test_fully_turbulent_estimate_brackets_the_data():
    rows = impl.read_ladson()
    Re, _, cd = impl.ladson_cd_near_zero_lift(rows, "I")
    est = impl.section_cd(Re, 0.12, 0.0)
    assert np.all(est > cd)                                 # free transition: laminar runs lower the drag
    gap = (est - cd) / cd
    assert gap[0] > gap[-1]                                 # the gap closes as Re rises
    ReG, _, cdG = impl.ladson_cd_near_zero_lift(rows, "VII")
    assert np.all(impl.section_cd(ReG, 0.12, 0.05) < cdG)   # grit adds drag beyond a smooth turbulent layer


def test_effective_transition_reproduces_measurement():
    rows = impl.read_ladson()
    Re, _, cd = impl.ladson_cd_near_zero_lift(rows, "I")
    xt = np.array([impl.effective_transition(r, c) for r, c in zip(Re, cd)])
    for r, c, x in zip(Re, cd, xt):
        assert float(impl.section_cd(r, 0.12, x)) == pytest.approx(c, rel=1e-8)
    assert np.all((xt > 0) & (xt < 1))
    assert np.all(np.diff(xt) < 0)                          # transition moves forward with Re
    assert np.isnan(impl.effective_transition(5e6, 0.02))  # no transition location gives that much drag


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "drag_build_up.png").stat().st_size > 10_000
