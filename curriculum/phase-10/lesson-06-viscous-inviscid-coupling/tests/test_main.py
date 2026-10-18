"""Public tests for P10.06 Viscous-Inviscid Coupling. They assert the boundary-layer physics and the
checks the lesson relies on: the Falkner-Skan solutions reproduce Blasius and Hiemenz; Thwaites gives
theta sqrt(Re_x)/x = sqrt(0.45) on a flat plate (within 1.5 % of Blasius) and a constant theta in
stagnation flow (within 7 % of Hiemenz); Michel's criterion puts flat-plate transition at a few
million; Head's method reproduces the 1/7-power-law flat plate and separates in a strong adverse
gradient; Squire-Young reduces to 2 theta when the trailing-edge velocity equals the free stream; the
stagnation split is symmetric for a symmetric section; drag falls with Re and rises with a trip; the
displacement thickness decambers the section; and, when the authentic Ladson data file is present,
the NACA 0012 prediction is within 25 % of the measured free-transition drag.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

N = 60      # panels per surface: enough for these properties, fast


# ---------------------------------------------------------------- similarity solutions and Thwaites

@pytest.fixture(scope="module")
def blasius():
    return impl.falkner_skan(0.0)


def test_falkner_skan_reproduces_blasius_and_hiemenz(blasius):
    assert blasius["fpp0"] == pytest.approx(0.4696, abs=2e-4)
    assert blasius["I_theta"] * np.sqrt(2) == pytest.approx(0.664, abs=1e-3)     # theta sqrt(Re_x)/x
    assert blasius["H"] == pytest.approx(2.59, abs=0.01)
    hz = impl.falkner_skan(1.0)
    assert hz["fpp0"] == pytest.approx(1.2326, abs=2e-4)
    assert hz["H"] == pytest.approx(2.216, abs=0.005)


def test_thwaites_flat_plate_against_blasius(blasius):
    s = np.linspace(0, 2.0, 1001)
    nu = 2e-6
    r = impl.thwaites(s, np.ones_like(s), nu)
    k = r["theta"][1:] * np.sqrt(1 / (nu * s[1:]))                            # theta sqrt(Re_x)/x
    npt.assert_allclose(k, np.sqrt(0.45), rtol=1e-10)
    assert k[-1] == pytest.approx(blasius["I_theta"] * np.sqrt(2), rel=0.015)
    npt.assert_allclose(r["lam"], 0.0, atol=1e-12)
    cf_blasius = 2 * blasius["fpp0"] / np.sqrt(2)                              # cf sqrt(Re_x)
    assert r["cf"][-1] * np.sqrt(s[-1] / nu) == pytest.approx(cf_blasius, rel=0.015)


def test_thwaites_stagnation_flow_against_hiemenz():
    s = np.linspace(0, 0.5, 801)
    C, nu = 4.0, 1e-6
    r = impl.thwaites(s, C * s, nu)
    npt.assert_allclose(r["theta"], np.sqrt(0.075 * nu / C), rtol=1e-6)        # constant thickness
    npt.assert_allclose(r["lam"], 0.075, rtol=1e-6)
    hz = impl.falkner_skan(1.0)
    assert r["theta"][400] == pytest.approx(hz["I_theta"] * np.sqrt(nu / C), rel=0.07)


def test_thwaites_correlations_are_continuous_and_signal_separation():
    H, l = impl.thwaites_H_l(np.array([-1e-9, 0.0, 1e-9]))
    assert np.ptp(H) < 1e-3 and np.ptp(l) < 1e-6
    lam = np.linspace(-0.09, 0.1, 50)
    H, l = impl.thwaites_H_l(lam)
    assert np.all(np.diff(H) < 0) and np.all(np.diff(l) > 0)                  # adverse gradient: fuller, less shear
    s = np.linspace(0, 1, 801)
    r = impl.thwaites(s, 1.0 - 0.6 * s, 1e-6)                                  # strongly decelerating
    assert r["i_sep"] is not None and r["lam"][r["i_sep"]] < -0.09


def test_michel_transition_on_a_flat_plate(blasius):
    Re_x = np.logspace(4, 8, 400)
    Re_th = blasius["I_theta"] * np.sqrt(2) * np.sqrt(Re_x)
    tr = impl.michel_transition(Re_x, Re_th)
    assert not tr[0] and tr[-1]
    i = int(np.argmax(tr))
    assert np.all(tr[i:])                                                      # once crossed, stays crossed
    assert 2e6 < Re_x[i] < 1e7


# ---------------------------------------------------------------- Head and Squire-Young

def test_head_shape_factor_inverse():
    H = np.linspace(1.2, 2.8, 40)
    npt.assert_allclose(impl.head_H_of_H1(impl.head_H1(H)), H, rtol=1e-10)


def test_head_flat_plate_against_power_law():
    Re_L = 1e7
    x = np.linspace(0, 1, 2001)
    lam = impl.thwaites(x, np.ones_like(x), 1 / Re_L)
    i0 = 20
    t = impl.head_turbulent(x, np.ones_like(x), 1 / Re_L, x[i0], lam["theta"][i0])
    assert t["s_sep"] is None
    assert t["theta"][-1] == pytest.approx(0.036 * Re_L ** -0.2, rel=0.05)    # 1/7-power law
    assert t["cf"][-1] == pytest.approx(0.0592 * Re_L ** -0.2, rel=0.05)
    assert 1.25 < t["H"][-1] < 1.45 and np.all(np.diff(t["theta"]) > 0)


def test_head_separates_in_a_strong_adverse_gradient():
    s = np.linspace(0, 1, 801)
    ue = 1.0 - 0.5 * s
    t = impl.head_turbulent(s, ue, 1e-6, 0.1, 3e-4)
    assert t["s_sep"] is not None and 0.1 < t["s_sep"] < 1.0
    assert t["H"][-1] == pytest.approx(2.4)


def test_squire_young_limits():
    assert impl.squire_young(0.002, 1.4, 1.0) == pytest.approx(0.004)
    assert impl.squire_young(0.002, 2.0, 1.0) == pytest.approx(0.004)         # H drops out at ue = V
    assert impl.squire_young(0.002, 1.4, 0.9) < impl.squire_young(0.002, 1.4, 1.0)
    assert impl.squire_young(0.002, 1.4, 0.9, chord=2.0) == pytest.approx(0.5 * impl.squire_young(0.002, 1.4, 0.9))


# ---------------------------------------------------------------- coupling

def test_stagnation_split_of_a_symmetric_section():
    x, y = impl.naca4("0012", N)
    ps = impl.panel_solve(x, y, 0.0)
    up, lo = impl.split_at_stagnation(x, y, ps["Vt"])
    for br in (up, lo):
        assert br["s"][0] == 0.0 and br["ue"][0] == 0.0
        assert np.all(np.diff(br["s"]) > 0)
    assert up["x"][0] == pytest.approx(0.0, abs=1e-3)                          # stagnation at the nose
    npt.assert_allclose(up["ue"], lo["ue"], atol=1e-9)
    up4, lo4 = impl.split_at_stagnation(x, y, impl.panel_solve(x, y, np.radians(4.0))["Vt"])
    assert up4["s"][-1] > lo4["s"][-1]                                         # stagnation moves onto the lower surface


def test_symmetric_section_drag_splits_equally_and_falls_with_Re():
    cds = []
    for Re in (2e6, 6e6, 1.8e7):
        r = impl.section_drag("0012", 0.0, Re, N)
        assert r["upper"]["cd"] == pytest.approx(r["lower"]["cd"], rel=1e-5)        # ODE tolerance
        assert abs(r["cl"]) < 1e-10
        cds.append(r["cd"])
    assert cds[0] > cds[1] > cds[2]


def test_trip_moves_transition_forward_and_raises_drag():
    free = impl.section_drag("0012", 0.0, 6e6, N)
    trip = impl.section_drag("0012", 0.0, 6e6, N, trip_x=0.05)
    assert trip["upper"]["mechanism"] == "forced"
    assert trip["upper"]["x_tr"] == pytest.approx(0.05, abs=0.02)
    assert free["upper"]["x_tr"] > 0.15
    assert trip["cd"] > 1.15 * free["cd"]


def test_displacement_thickness_decambers_the_section():
    r = impl.displacement_iteration("2412", np.radians(2.0), 4.5e6, N)
    assert np.all(r["dstar"] >= 0)
    assert 0.0 < (r["cl0"] - r["cl1"]) / r["cl0"] < 0.06


DATA = impl.find_data_file()


@pytest.mark.skipif(DATA is None, reason="authentic Ladson data file not present")
def test_naca0012_drag_against_ladson_free_transition():
    d = impl.read_ladson(DATA)
    Re, cd = impl.ladson_zero_lift_drag(d, "I")
    assert Re.size >= 4 and np.all(cd > 0.005) and np.all(cd < 0.008)
    for R, c in zip(Re, cd):
        if R <= 1.25e7:                                                        # beyond, transition modelling dominates
            pred = impl.section_drag("0012", 0.0, R, N)["cd"]
            assert pred == pytest.approx(c, rel=0.25)
    Ref, cdf = impl.ladson_zero_lift_drag(d, "VII")
    assert np.all(cdf > cd.max())                                              # the wraparound grit adds drag


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path, n=40)
    assert (tmp_path / "viscous_drag.png").stat().st_size > 1000
