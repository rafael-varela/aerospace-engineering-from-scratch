"""Public tests for P07.13 Heat Transfer: Convection. They assert properties, not magic numbers:
the interpolated air and water properties reproduce every node of Lienhard's tables and are
internally consistent (nu = mu/rho, Pr = mu cp/k); the laminar flat plate satisfies the
Reynolds-Colburn analogy St Pr^(2/3) = c_f/2 exactly and its average Nusselt number is twice the
local value at the trailing edge; the turbulent power laws obey the analogy to the factor
Pr^(0.6 - 1/3) that their Prandtl exponents imply; the mixed laminar-turbulent average equals the
numerical average of the local coefficient and is continuous at transition; correlations have the
right limits (Churchill-Bernstein -> 0.3, Churchill-Chu -> 0.68 and 0.36 as Re, Ra -> 0) and
scalings (Re^(1/2), Ra^(1/4), Ra^(1/3) with h independent of L); Gnielinski and Dittus-Boelter
agree within 20 % for gases; and the three applications close their energy balances.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)


# ---------------------------------------------------------------- properties

def test_air_properties_reproduce_table_nodes_and_are_consistent():
    tab = impl.AIR_TABLE_A6
    p = impl.air_properties(tab[:, 0])
    npt.assert_allclose(p["rho"], tab[:, 1], rtol=1e-12)
    npt.assert_allclose(p["k"], tab[:, 5], rtol=1e-12)
    npt.assert_allclose(p["nu"], tab[:, 4] * 1e-5, rtol=1e-12)
    npt.assert_allclose(p["Pr"], tab[:, 7], rtol=1e-12)
    # the table's own columns agree with each other to its rounding (4 significant figures)
    npt.assert_allclose(p["mu"] / p["rho"], p["nu"], rtol=3e-3)
    npt.assert_allclose(p["mu"] * p["cp"] / p["k"], p["Pr"], rtol=1e-2)
    npt.assert_allclose(p["beta"], 1.0 / tab[:, 0])
    # between nodes, linear interpolation lies between the neighbours
    q = impl.air_properties(305.0)
    assert 0.0264 < float(q["k"]) < 0.0271 and float(q["k"]) == pytest.approx(0.02675)
    with pytest.raises(ValueError):
        impl.air_properties(150.0)


def test_water_properties_reproduce_table_nodes():
    tab = impl.WATER_TABLE_A3
    w = impl.water_properties(tab[:, 0])
    npt.assert_allclose(w["k"], tab[:, 3], rtol=1e-12)
    npt.assert_allclose(w["Pr"], tab[:, 6], rtol=1e-12)
    npt.assert_allclose(w["nu"] * w["rho"] * w["cp"] / w["k"], w["Pr"], rtol=1.5e-2)
    with pytest.raises(ValueError):
        impl.water_properties(400.0)


# ---------------------------------------------------------------- flat plate: the independent checks

@pytest.mark.parametrize("Pr", [0.7, 1.0, 5.85])
def test_laminar_reynolds_colburn_analogy_holds_exactly(Pr):
    Re = np.logspace(3, 5.6, 9)
    St = impl.stanton(impl.plate_laminar_local_nu(Re, Pr), Re, Pr)
    npt.assert_allclose(St * Pr ** (2 / 3), impl.plate_laminar_cf_local(Re) / 2, rtol=1e-12)


def test_turbulent_analogy_differs_by_the_prandtl_factor_only():
    Re = np.logspace(5.8, 7.5, 7)
    for Pr in (0.7, 1.0):
        St = impl.stanton(impl.plate_turbulent_local_nu(Re, Pr), Re, Pr)
        ratio = St * Pr ** (2 / 3) / (impl.plate_turbulent_cf_local(Re) / 2)
        npt.assert_allclose(ratio, Pr ** (0.6 - 1 / 3), rtol=1e-12)            # exact at Pr = 1
    assert abs(0.7 ** (0.6 - 1 / 3) - 1) < 0.1                                  # within 10 % for air


def test_laminar_average_is_twice_the_trailing_edge_local_value():
    for Re in (1e3, 4e4, 3e5):
        for Pr in (0.71, 7.0):
            assert impl.plate_laminar_average_nu(Re, Pr) == pytest.approx(2 * impl.plate_laminar_local_nu(Re, Pr), rel=1e-12)


def test_mixed_average_matches_numerical_integration_and_is_continuous():
    Pr = 0.707
    for ReL in (1e4, 3e5, 8e5, 3e6, 1e7):
        assert impl.plate_average_nu(ReL, Pr) == pytest.approx(impl.plate_average_nu_numerical(ReL, Pr), rel=1e-7)
    Rc = impl.RE_CRIT_PLATE
    lo, hi = impl.plate_average_nu(Rc * (1 - 1e-9), Pr), impl.plate_average_nu(Rc * (1 + 1e-9), Pr)
    assert hi == pytest.approx(lo, rel=1e-7)                                    # no jump at transition
    # beyond transition the average rises faster than the laminar law would predict
    assert impl.plate_average_nu(4 * Rc, Pr) > impl.plate_laminar_average_nu(4 * Rc, Pr)
    # the local value does jump at a sharp transition (turbulent above laminar)
    assert impl.plate_local_nu(Rc * 1.0001, Pr) > 1.5 * impl.plate_local_nu(Rc * 0.9999, Pr)


def test_laminar_plate_scales_as_root_Re():
    Pr = 0.71
    assert impl.plate_laminar_average_nu(4e4, Pr) / impl.plate_laminar_average_nu(1e4, Pr) == pytest.approx(2.0)
    # h ~ U^(1/2) L^(-1/2): doubling L at fixed U lowers h by 1/sqrt(2)
    k, nu, U = 0.0264, 1.575e-5, 3.0
    h = lambda L: impl.h_from_nu(impl.plate_laminar_average_nu(impl.reynolds(U, L, nu), Pr), k, L)
    assert h(0.4) / h(0.2) == pytest.approx(1 / np.sqrt(2))


# ---------------------------------------------------------------- cylinder and pipe

def test_churchill_bernstein_limits():
    Pr = 0.71
    assert impl.cylinder_crossflow_nu(1e-10, Pr) == pytest.approx(0.3, abs=1e-4)
    # at low Re the high-Re bracket [1 + (Re/282000)^(5/8)]^(4/5) is ~1, so the simplified form (7.66)
    # agrees: within 2.5 % up to Re = 1000 and 6 % at Re = 4000, where Lienhard draws the line
    for Re, tol in ((10.0, 0.025), (100.0, 0.025), (1000.0, 0.025), (4000.0, 0.06)):
        simple = 0.3 + 0.62 * np.sqrt(Re) * Pr ** (1 / 3) / (1 + (0.4 / Pr) ** (2 / 3)) ** 0.25
        assert impl.cylinder_crossflow_nu(Re, Pr) == pytest.approx(simple, rel=tol)
    Re = np.logspace(0, 7, 50)
    assert np.all(np.diff(impl.cylinder_crossflow_nu(Re, Pr)) > 0)


def test_pipe_correlations_and_their_range():
    assert impl.pipe_nu(1000, 0.7, wall="q") == pytest.approx(4.364)
    assert impl.pipe_nu(1000, 0.7, wall="T") == pytest.approx(3.657)
    with pytest.raises(ValueError):
        impl.gnielinski_nu(1500.0, 0.7)
    Re = np.array([1e4, 3e4, 1e5, 1e6])
    # for gases Dittus-Boelter and Gnielinski agree within 20 % in fully turbulent flow
    r = impl.dittus_boelter_nu(Re, 0.71) / impl.gnielinski_nu(Re, 0.71)
    assert np.all(np.abs(r - 1) < 0.2)
    # Filonenko's friction factor near Blasius' 0.316 Re^-1/4 at Re = 1e4 (within 6 %)
    assert impl.filonenko_f(1e4) == pytest.approx(0.316 * 1e4 ** -0.25, rel=0.06)
    assert np.all(np.diff(impl.gnielinski_nu(np.logspace(3.5, 6, 20), 3.0)) > 0)
    # the correlations written out (Filonenko f, Gnielinski 7.41, Dittus-Boelter) pin every coefficient
    for Re_, Pr_ in ((1e4, 0.71), (1e5, 7.0)):
        f = (1.82 * np.log10(Re_) - 1.64) ** -2
        assert impl.filonenko_f(Re_) == pytest.approx(f, rel=1e-12)
        nu_g = (f / 8) * (Re_ - 1000) * Pr_ / (1 + 12.7 * np.sqrt(f / 8) * (Pr_ ** (2 / 3) - 1))
        assert impl.gnielinski_nu(Re_, Pr_) == pytest.approx(nu_g, rel=1e-10)
        assert impl.dittus_boelter_nu(Re_, Pr_) == pytest.approx(0.023 * Re_ ** 0.8 * Pr_ ** 0.4, rel=1e-12)
        assert impl.dittus_boelter_nu(Re_, Pr_, heating=False) == pytest.approx(0.023 * Re_ ** 0.8 * Pr_ ** 0.3, rel=1e-12)


# ---------------------------------------------------------------- natural convection

def test_natural_convection_limits_and_scalings():
    Pr = 0.71
    assert impl.vertical_plate_nu(1e-12, Pr) == pytest.approx(0.68, abs=1e-3)
    assert impl.horizontal_cylinder_nu(1e-12, Pr) == pytest.approx(0.36, abs=1e-3)
    # laminar boundary layers: Nu - constant scales as Ra^(1/4)
    for f, c in ((impl.vertical_plate_nu, 0.68), (impl.horizontal_cylinder_nu, 0.36)):
        assert (f(1.6e7, Pr) - c) / (f(1e6, Pr) - c) == pytest.approx(2.0, rel=1e-10)
    # turbulent upper face: Nu ~ Ra^(1/3), so h does not depend on the plate size
    assert impl.horizontal_plate_up_nu(8e8, Pr) / impl.horizontal_plate_up_nu(1e8, Pr) == pytest.approx(2.0, rel=1e-10)
    # the all-Ra form of Churchill-Chu stays within 25 % of the laminar form for 1e4 <= Ra <= 1e8
    Ra = np.logspace(4, 8, 9)
    assert np.all(np.abs(impl.vertical_plate_nu_allra(Ra, Pr) / impl.vertical_plate_nu(Ra, Pr) - 1) < 0.25)
    # thick-boundary-layer correction raises Nu when it applies, and is continuous where Nu = 10
    assert impl.horizontal_plate_up_nu(10.0, Pr) > 0.560 * 10 ** 0.25 / (1 + (0.492 / Pr) ** (9 / 16)) ** (4 / 9)


def test_h_natural_is_symmetric_for_cold_and_hot_vertical_plates_and_zero_at_zero_dT():
    Tinf = 300.0
    assert impl.h_natural("vertical", 300.0, Tinf, 0.2) == 0.0
    h_hot = impl.h_natural("vertical", 320.0, Tinf, 0.2)
    h_cold = impl.h_natural("vertical", 280.0, Tinf, 0.2)
    assert h_cold == pytest.approx(h_hot, rel=0.15)        # same |dT|, properties at slightly different film T
    # a hot plate facing up loses more than one facing down (unstable vs stable)
    assert impl.h_natural("up", 330.0, Tinf, 0.1) > impl.h_natural("down", 330.0, Tinf, 0.1)
    # a cold plate facing up behaves like a hot plate facing down
    assert impl.h_natural("up", 270.0, Tinf, 0.1) == pytest.approx(impl.h_natural("down", 330.0, Tinf, 0.1), rel=0.15)


def test_radiation_linearisation_is_exact():
    eps, Tw, Ts = 0.9, 330.0, 295.0
    assert impl.h_radiation_linear(eps, Tw, Ts) * (Tw - Ts) == pytest.approx(eps * impl.SIGMA * (Tw ** 4 - Ts ** 4), rel=1e-12)
    assert impl.h_radiation_linear(eps, 300.0, 300.0) == pytest.approx(4 * eps * impl.SIGMA * 300.0 ** 3)


# ---------------------------------------------------------------- applications

def test_box_surface_temperature_closes_the_energy_balance_and_is_monotone():
    L, W, H, Tinf = 0.30, 0.12, 0.20, 298.15
    for eps in (0.0, 0.9):
        Tw = impl.box_surface_temperature(60.0, Tinf, L, W, H, eps)
        assert impl.box_heat_loss_natural(Tw, Tinf, L, W, H, eps)["Q"] == pytest.approx(60.0, rel=1e-9)
    Ts = [impl.box_surface_temperature(q, Tinf, L, W, H, 0.9) for q in (10, 30, 60, 90)]
    assert np.all(np.diff(Ts) > 0)
    # radiation helps: the painted box runs cooler than the bare one
    assert impl.box_surface_temperature(60.0, Tinf, L, W, H, 0.9) < impl.box_surface_temperature(60.0, Tinf, L, W, H, 0.0)
    # natural convection coefficients in still air are a few W/m^2 K
    r = impl.box_heat_loss_natural(330.0, Tinf, L, W, H, 0.0)
    assert 2.0 < r["h_conv"] < 10.0


def test_required_speed_meets_the_load_and_vanishes_when_not_needed():
    L, W, H, Tinf = 0.30, 0.12, 0.20, 298.15
    U = impl.box_required_speed(60.0, 323.15, Tinf, L, W, H, 0.9)
    assert U > 0
    assert impl.box_heat_loss_forced(U, 323.15, Tinf, L, W, H, 0.9)["Q"] == pytest.approx(60.0, rel=1e-8)
    assert impl.box_required_speed(5.0, 323.15, Tinf, L, W, H, 0.9) == 0.0


def test_cold_plate_energy_balance_and_regime_jump():
    Q, D, Lc, Tin = 400.0, 0.004, 1.0, 293.15
    for mdot in (0.005, 0.02):
        c = impl.cold_plate(Q, mdot, D, Lc, Tin)
        Tb = 0.5 * (Tin + c["T_out"])
        cp = float(impl.water_properties(Tb)["cp"])
        assert mdot * cp * (c["T_out"] - Tin) == pytest.approx(Q, rel=1e-9)     # first law on the coolant
        assert c["q_flux"] * np.pi * D * Lc == pytest.approx(Q)
        assert c["T_wall_out"] > c["T_out"] > Tin
    lam, turb = impl.cold_plate(Q, 0.005, D, Lc, Tin), impl.cold_plate(Q, 0.02, D, Lc, Tin)
    assert lam["Re"] < 2300 < turb["Re"] and lam["Nu"] == pytest.approx(4.364)
    assert turb["h"] > 5 * lam["h"]


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "convection.png").stat().st_size > 10_000
