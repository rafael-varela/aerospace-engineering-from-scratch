"""Public tests for P10.08 High-Lift Devices. They assert flap and slat physics, not magic numbers:
the hinge sits where cos(theta_h) = 2E - 1; the flap effectiveness runs from 0 to 1 and agrees with
quadrature of the kinked camber line by the given thin-airfoil solver; flap increments are linear in
the deflection and superpose on camber; the largest nose-down moment per degree belongs to a 25 %
flap; the hinge moment is restoring, linear and has the right whole-plate limit; the panel method
with a deflected flap approaches thin-airfoil theory as the thickness goes to zero (the independent
route); the multi-element panel model shows the slat effect (lower main-element suction peak) and
the circulation effect (more main-element lift ahead of a slotted flap); the suction-peak cl_max
criterion is self-consistent on the clean section; stall speed scales as it must.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)


# ---------------------------------------------------------------- thin-airfoil flap theory

def test_hinge_angle_places_the_hinge():
    assert impl.hinge_angle(0.5) == pytest.approx(np.pi / 2)
    assert impl.hinge_angle(1.0) == pytest.approx(0.0, abs=1e-12)
    assert impl.hinge_angle(0.0) == pytest.approx(np.pi)
    for E in (0.1, 0.3, 0.7):
        assert 0.5 * (1 - np.cos(impl.hinge_angle(E))) == pytest.approx(1 - E)


def test_flap_effectiveness_limits_and_monotonic():
    assert impl.flap_effectiveness(1e-6) == pytest.approx(4 / np.pi * 1e-3, rel=1e-3)   # small flap: (4/pi) sqrt(E)
    assert impl.flap_effectiveness(1.0) == pytest.approx(1.0)                 # all-moving: delta acts as alpha
    assert impl.flap_effectiveness(0.25) == pytest.approx(1 - (2 * np.pi / 3 - np.sqrt(3) / 2) / np.pi)
    t = [impl.flap_effectiveness(E) for E in np.linspace(0.01, 0.99, 50)]
    assert np.all(np.diff(t) > 0)
    assert impl.flap_effectiveness(0.2) > 0.2                                 # small flaps punch above their chord


def test_flap_increments_match_quadrature_of_the_kinked_camber_line():
    for E in (0.15, 0.25, 0.4):
        d = np.radians(10.0)
        inc = impl.flap_increments(E, d)
        num = impl.thin_airfoil(impl.flap_slope(E, d), 0.0, n_terms=5, m=40000)
        assert num["cl"] == pytest.approx(inc["dcl"], rel=1e-4)
        assert num["cm_c4"] == pytest.approx(inc["dcm_c4"], rel=1e-3)
        assert num["alpha_L0"] == pytest.approx(inc["dalpha_L0"], rel=1e-4)


def test_flap_is_linear_and_superposes_on_camber():
    m, p = 0.02, 0.4
    base = lambda x: np.where(x < p, 2 * m / p ** 2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))
    a = np.radians(3.0)
    c0 = impl.thin_airfoil(base, a)["cl"]
    for deg in (5.0, 10.0, 20.0):
        d = np.radians(deg)
        c = impl.thin_airfoil(impl.flap_slope(0.25, d, base), a, m=40000)["cl"]
        assert c - c0 == pytest.approx(impl.flap_increments(0.25, d)["dcl"], rel=1e-4)
    assert impl.flap_increments(0.25, 0.2)["dcl"] == pytest.approx(2 * impl.flap_increments(0.25, 0.1)["dcl"])


def test_largest_nose_down_moment_per_degree_is_a_quarter_chord_flap():
    Es = np.linspace(0.05, 0.95, 181)
    dcm = np.array([impl.flap_increments(E, 1.0)["dcm_c4"] for E in Es])
    assert np.all(dcm < 0)                                                    # TE-down flap pitches nose down
    assert Es[np.argmin(dcm)] == pytest.approx(0.25, abs=0.006)
    assert dcm.min() == pytest.approx(-3 * np.sqrt(3) / 8, rel=1e-4)


def test_hinge_moment_is_restoring_linear_and_has_the_whole_plate_limit():
    E = 0.25
    cha = (impl.hinge_moment(E, 0.02, 0.0) - impl.hinge_moment(E, -0.02, 0.0)) / 0.04
    chd = (impl.hinge_moment(E, 0.0, 0.02) - impl.hinge_moment(E, 0.0, -0.02)) / 0.04
    assert cha < 0 and chd < 0
    assert impl.hinge_moment(E, 0.0, 0.0) == pytest.approx(0.0, abs=1e-10)
    combo = impl.hinge_moment(E, 0.03, 0.05)
    assert combo == pytest.approx(0.03 * cha + 0.05 * chd, rel=1e-3)
    # a "flap" of the whole chord is the plate itself turning about its leading edge: cm_LE = -cl/4
    a, d = 0.02, 0.03
    assert impl.hinge_moment(0.999, a, d) == pytest.approx(-np.pi / 2 * (a + d), rel=0.02)
    # quadrature converged
    assert impl.hinge_moment(E, 0.0, 0.1, n_terms=300) == pytest.approx(
        impl.hinge_moment(E, 0.0, 0.1, n_terms=600, m=16000), rel=2e-3)


# ---------------------------------------------------------------- panel model of a plain flap

def test_flapped_section_geometry():
    X0, Y0 = impl.naca4("2412", 60)
    X, Y = impl.flapped_section("2412", 0.25, 0.0, 60)
    npt.assert_allclose(X, X0, atol=1e-15); npt.assert_allclose(Y, Y0, atol=1e-15)
    X, Y = impl.flapped_section("2412", 0.25, 20.0, 60)
    assert X[0] == pytest.approx(X[-1]) and Y[0] == pytest.approx(Y[-1])     # still closed
    assert Y[0] == pytest.approx(-0.25 * np.sin(np.radians(20.0)), abs=0.01)  # TE drops by about cf sin(delta)
    i_le = int(np.argmin(X))
    assert np.all(np.diff(X[: i_le + 1]) < 0) and np.all(np.diff(X[i_le:]) > 0)   # no fold in the contour


def test_panel_flap_derivative_approaches_thin_airfoil_as_thickness_vanishes():
    tau = impl.flap_effectiveness(0.25)
    d = [impl.panel_flap_derivative(t, 0.25, n=100) for t in (0.02, 0.06, 0.12)]
    assert d[0] == pytest.approx(2 * np.pi * tau, rel=0.01)                   # independent route
    assert d[0] < d[1] < d[2]                                                  # thickness adds effectiveness
    assert d[2] / d[0] < 1.15


def test_multi_element_solver_far_apart_elements_are_independent():
    A = impl.naca4("0012", 60)
    B = (A[0] + 1000.0, A[1] + 600.0)                                        # mutual upwash ~ Gamma/(2 pi d)
    alone = impl.panel_solve([A], 5.0)["cl"]
    r = impl.panel_solve([A, B], 5.0)
    assert r["parts"][0]["cl"] == pytest.approx(alone, rel=2e-3)
    assert r["parts"][1]["cl"] == pytest.approx(alone, rel=2e-3)
    assert r["cl"] == pytest.approx(r["cl_kj"], rel=0.01)                    # pressure vs Kutta-Joukowski


def test_element_peaks_bookkeeping():
    b = impl.panel_basis(impl.slat_configuration())
    ep = impl.element_peaks(b, 8.0)
    assert sum(ep["cls"]) == pytest.approx(ep["cl"])
    assert len(ep["peaks"]) == 2
    for pk, part in zip(ep["peaks"], ep["result"]["parts"]):
        assert pk >= part["cp"].min() - 1e-12 and pk < 1.0


def test_slat_effect_lowers_the_main_element_peak():
    clean = impl.panel_basis([impl.naca4("2412", 100)])
    slat = impl.panel_basis(impl.slat_configuration())
    for a in (8.0, 12.0):
        c, s = impl.element_peaks(clean, a), impl.element_peaks(slat, a)
        assert s["peaks"][1] > 0.5 * c["peaks"][0]                            # main peak more than halved
        assert s["cls"][1] < c["cl"]                                          # the slat unloads the main element
        assert s["cl"] > c["cl"]                                              # but the pair lifts more


def test_circulation_effect_of_a_slotted_flap():
    clean = impl.element_peaks(impl.panel_basis([impl.naca4("2412", 100)]), 2.0)
    fl = impl.element_peaks(impl.panel_basis(impl.slotted_flap_configuration(25.0)), 2.0)
    assert fl["cls"][0] > clean["cl"] + 0.5                                   # flap raises the main element's lift
    te_clean = clean["result"]["parts"][0]["cp"][[0, -1]].mean()
    te_flap = fl["result"]["parts"][0]["cp"][[0, -1]].mean()
    assert te_flap < te_clean                                                 # dumping: main TE in faster flow


# ---------------------------------------------------------------- cl_max and the landing payoff

def test_peak_criterion_is_self_consistent_on_the_clean_section():
    clean = impl.naca4("2412", 100)
    a1, cp1 = impl.critical_peak(clean, 1.4)
    a2, cp2 = impl.critical_peak(clean, 1.6)
    assert a2 > a1 and cp2 < cp1 < 0
    r = impl.clmax_by_peak([clean], cp2)
    assert r["cl"] == pytest.approx(1.6, abs=1e-4) and r["alpha"] == pytest.approx(a2, abs=1e-3)


def test_devices_raise_clmax_plain_flap_less_than_its_lift_shift():
    clean = impl.naca4("2412", 100)
    a_c, cpc = impl.critical_peak(clean, 1.6)
    pf = impl.clmax_by_peak([impl.flapped_section("2412", 0.25, 10.0, 100)], cpc)
    shift = impl.flap_increments(0.25, np.radians(10.0))["dcl"]
    assert 0 < pf["cl"] - 1.6 < shift                                         # flap lowers the stall angle
    assert pf["alpha"] < a_c
    sl = impl.clmax_by_peak(impl.slat_configuration(), cpc)
    assert sl["cl"] > 1.6 and sl["alpha"] > a_c                               # slat delays the stall angle


def test_stall_speed_and_wing_clmax():
    v = impl.stall_speed(750.0, 12.0, 1.5)
    assert v == pytest.approx(np.sqrt(2 * 750 * 9.81 / (1.225 * 12.0 * 1.5)))
    assert impl.stall_speed(4 * 750.0, 12.0, 1.5) == pytest.approx(2 * v)    # V_s ~ sqrt(W)
    assert impl.stall_speed(750.0, 12.0, 6.0) == pytest.approx(v / 2)        # V_s ~ 1/sqrt(CL_max)
    assert impl.wing_clmax(1.6, 0.0, 0.9, 0.6) == pytest.approx(1.44)
    assert impl.wing_clmax(1.6, 1.0, 0.9, 0.6) == pytest.approx(0.9 * (1.6 + 0.6))


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "high_lift.png").stat().st_size > 10_000
