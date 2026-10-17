"""Public tests for P08.14 Turbulent Boundary Layers. They assert physics, not magic numbers:
Spalding's composite law tends to u+ = y+ at the wall and to the log law far from it, and its
Newton inversion is exact to round-off; a Clauser fit recovers the friction velocity used to build a
SEEDED SYNTHETIC profile (labelled synthetic); the 1/7th-law plate satisfies the momentum integral
and integrates to 5/4 of its local c_f; the correlations (labelled correlations) agree with each
other to a few per cent, an integrated local correlation agrees with the integrated ones, and the
Karman-Schoenherr root satisfies its own equation; the mixed plate lies between laminar and
turbulent and reduces to them in the limits; and a friction law built from the log-wake profile and
the momentum integral alone, with no correlation, agrees with the correlations.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

NU, U = 1.863e-5, 60.0           # the light aeroplane at 3000 m


# ---------------------------------------------------------------- law of the wall

def test_log_law_values_and_slope():
    assert impl.log_law(1.0) == pytest.approx(5.0)
    y = np.array([100.0, 100.0 * np.e])
    assert np.diff(impl.log_law(y))[0] == pytest.approx(1 / 0.41)    # one e-fold of y+ adds 1/kappa
    assert impl.log_law(np.e, kappa=0.384, B=4.17) == pytest.approx(1 / 0.384 + 4.17)


def test_spalding_limits_sublayer_and_log_law():
    up = np.array([0.01, 0.1, 1.0])
    npt.assert_allclose(impl.spalding_yplus(up), up, rtol=2e-3)       # u+ = y+ at the wall
    yp = np.array([1e3, 1e4, 1e5])
    npt.assert_allclose(impl.spalding_uplus(yp), impl.log_law(yp), atol=0.05)   # log law far from it
    yb = np.array([5.0, 10.0, 30.0])
    assert np.all(impl.spalding_uplus(yb) < np.minimum(yb, impl.log_law(yb)))   # buffer layer below both


def test_spalding_inversion_is_exact_and_monotone():
    yp = np.concatenate([[0.0], np.geomspace(1e-3, 1e6, 200)])
    up = impl.spalding_uplus(yp)
    assert up[0] == 0.0
    npt.assert_allclose(impl.spalding_yplus(up), yp, rtol=1e-10, atol=1e-12)
    assert np.all(np.diff(up) > 0)


def test_clauser_fit_recovers_synthetic_u_tau():
    # SYNTHETIC profile (Spalding + Coles wake, 1 % noise, seeded). Over 40 seeds the fit is biased
    # by -0.6 % (buffer layer and wake inside the window) with a scatter of 0.23 %; 2 % is > 5 sigma
    ut, delta = 2.0, 0.025
    y = np.geomspace(2 * NU / ut, delta, 40)
    Uy = impl.synthetic_profile(y, ut, NU, delta, noise=0.01, seed=3)
    ut_fit, npts = impl.clauser_fit(y, Uy, NU, delta)
    assert ut_fit == pytest.approx(ut, rel=0.02) and npts >= 5
    assert (y * ut_fit / NU >= 30).sum() >= npts                      # every point used is in the log region
    # a clean log-law profile is fitted exactly
    yl = np.geomspace(1e-3, 5e-3, 12)
    ul = ut * impl.log_law(yl * ut / NU)
    assert impl.clauser_fit(yl, ul, NU, delta=1.0)[0] == pytest.approx(ut, rel=1e-8)


def test_clauser_fit_refuses_a_profile_with_no_log_region():
    y = np.geomspace(1e-6, 2e-6, 10)                                  # all inside the sublayer
    with pytest.raises(ValueError):
        impl.clauser_fit(y, 1.0 * y / NU * 1.0, NU, delta=0.02)


# ---------------------------------------------------------------- the 1/7th-law plate

def test_power_law_plate_satisfies_the_momentum_integral():
    x = np.linspace(0.05, 1.5, 3001)
    r = impl.power_law_bl(x, U, NU)
    dth = np.gradient(r["theta"], x)
    npt.assert_allclose(dth[5:-5], r["cf"][5:-5] / 2, rtol=1e-4)
    npt.assert_allclose(r["delta"] / x * r["Re_x"] ** 0.2, 0.3707, atol=1e-4)
    npt.assert_allclose(r["H"], 9 / 7) and npt.assert_allclose(r["theta"] / r["delta"], 7 / 72)


def test_integrated_cf_is_five_quarters_of_local_for_a_power_law():
    ReL = U * 1.5 / NU
    # C_F = (1/L) int c_f dx with c_f ~ x^(-1/5) gives 5/4 c_f(L); 0.074/0.0592 = 1.25 exactly
    assert float(impl.CF_plate(ReL, "power")) == pytest.approx(1.25 * float(impl.cf_local(ReL, "power")), rel=1e-12)
    assert impl.integrate_cf(ReL, "power") == pytest.approx(float(impl.CF_plate(ReL, "power")), rel=3e-3)


# ---------------------------------------------------------------- correlations

def test_correlations_agree_with_each_other():
    for R in (1e6, 5e6, 1e7, 1e8):
        kinds = ("power", "prandtl-schlichting", "karman-schoenherr", "ittc") if R <= 1e7 else \
            ("prandtl-schlichting", "karman-schoenherr", "ittc")       # the 1/7th-law fit is for Re < 1e7
        vals = [float(impl.CF_plate(R, k)) for k in kinds]
        assert max(vals) / min(vals) < 1.08                           # within 8 %
        loc = [float(impl.cf_local(R, k)) for k in ("power", "schultz-grunow", "white")[(R > 1e7):]]
        assert max(loc) / min(loc) < 1.10                             # local laws scatter a little more
    assert float(impl.cf_local(1e8, "power")) < 0.9 * float(impl.cf_local(1e8, "white"))   # power law fails high


def test_integrated_local_correlation_matches_integrated_correlations():
    # integrating a local correlation reproduces the integrated ones to the scatter of the
    # correlations themselves (8 %; the largest gap is White's local law against Karman-Schoenherr at 1e6)
    for R in (1e6, 1e7, 1e8):
        for loc in ("schultz-grunow", "white"):
            for tot in ("prandtl-schlichting", "karman-schoenherr"):
                assert impl.integrate_cf(R, loc) == pytest.approx(float(impl.CF_plate(R, tot)), rel=0.08)


def test_karman_schoenherr_root_satisfies_its_equation():
    R = np.array([1e5, 4.83e6, 1e9])
    C = impl.CF_plate(R, "karman-schoenherr")
    npt.assert_allclose(0.242 / np.sqrt(C), np.log10(R * C), rtol=1e-12)
    assert np.all(np.diff(C) < 0)
    with pytest.raises(ValueError):
        impl.CF_plate(1e6, "no-such-correlation")


def test_turbulent_friction_far_exceeds_laminar_and_falls_more_slowly():
    R = np.array([1e6, 1e7])
    turb = impl.CF_plate(R, "prandtl-schlichting")
    lam = impl.laminar_CF(R)
    assert np.all(turb > 3 * lam)
    assert turb[0] / turb[1] < lam[0] / lam[1]                       # Re^-0.2-ish against Re^-0.5


# ---------------------------------------------------------------- the mixed plate

def test_mixed_plate_limits_and_bounds():
    ReL = U * 1.5 / NU
    for Rtr in (5e5, 2e6):
        cp = impl.mixed_plate_prandtl(ReL, Rtr)
        cm, rx0 = impl.mixed_plate_momentum(ReL, Rtr)
        lam, turb = float(impl.laminar_CF(ReL)), float(impl.CF_plate(ReL))
        assert lam < cp < turb and lam < cm < turb * 1.02
        assert 0 < rx0 < Rtr                                         # virtual origin between leading edge and transition
    # transition at the trailing edge: all laminar; at the leading edge: all turbulent
    assert impl.mixed_plate_prandtl(ReL, ReL) == pytest.approx(float(impl.laminar_CF(ReL)), rel=1e-12)
    assert impl.mixed_plate_prandtl(ReL, 1e3) == pytest.approx(float(impl.CF_plate(ReL)), rel=1e-3)
    cm0, _ = impl.mixed_plate_momentum(ReL, 1.0)
    th = float(impl.power_law_bl(1.5, U, NU)["theta"])
    assert cm0 == pytest.approx(2 * th / 1.5, rel=1e-4)               # = the 1/7th-law plate's 2 theta/L


def test_later_transition_means_less_drag():
    ReL = 1e7
    cps = [impl.mixed_plate_prandtl(ReL, R) for R in (3e5, 1e6, 3e6)]
    cms = [impl.mixed_plate_momentum(ReL, R)[0] for R in (3e5, 1e6, 3e6)]
    assert cps[0] > cps[1] > cps[2] and cms[0] > cms[1] > cms[2]


# ---------------------------------------------------------------- independent friction law

def test_log_wake_state_is_consistent():
    st = impl.log_wake_state(1e4)
    assert st["cf"] == pytest.approx(2 / st["Ue_plus"] ** 2)
    assert 1.2 < st["H"] < 1.45                                       # turbulent: much fuller than laminar 2.59
    assert st["Ue_plus"] == pytest.approx(float(impl.spalding_uplus(1e4)) + 2 * 0.55 / 0.41, rel=1e-12)


def test_log_wake_friction_law_agrees_with_correlations():
    R = np.array([2e6, 1e7, 1e8])
    c = impl.log_wake_plate(R)
    for kind in ("white", "schultz-grunow"):
        npt.assert_allclose(c, impl.cf_local(R, kind), rtol=0.08)    # no correlation went into c
    assert np.all(np.diff(c) < 0)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "turbulent_bl.png").stat().st_size > 10_000
