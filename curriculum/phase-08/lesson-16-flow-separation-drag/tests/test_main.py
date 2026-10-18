"""Public tests for P08.16 Flow Separation & Drag.

They assert physics, not one magic number: attached potential flow round a cylinder has zero drag and
zero side force (d'Alembert's paradox); the pressure-drag integral is exact for a Fourier mode and
reproduces the closed form of the separated-flow model; the sphere correlations tend to Stokes's 24/Re
in the creeping limit, sit in the Newton range at moderate Re and (Morrison only) show the drag crisis;
the Strouhal number is absent below the shedding threshold and near 0.21 in the subcritical range;
terminal velocities satisfy the force balance, reduce to Stokes's law for tiny particles, match a time
integration of the equation of motion, and obey Reynolds similarity through the Archimedes number;
a fairing cuts the drag of a round strut by roughly an order of magnitude.

The drag coefficients here are CORRELATIONS (curve fits to experiments), never measurements.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

TH = 2 * np.pi * np.arange(256) / 256


# ---------------------------------------------------------------- pressure drag

def test_dalembert_attached_potential_flow_has_no_drag_and_no_side_force():
    cp = impl.cp_potential_cylinder(TH)
    assert cp[0] == pytest.approx(1.0) and cp[64] == pytest.approx(-3.0)     # stagnation; shoulder
    assert abs(impl.pressure_drag_coefficient(TH, cp)) < 1e-13
    assert abs(impl.pressure_lift_coefficient(TH, cp)) < 1e-13


def test_pressure_drag_integral_is_exact_for_fourier_modes():
    th = 2 * np.pi * np.arange(16) / 16
    for a, b in ((0.3, 0.7), (-1.0, 2.0)):
        cp = a + b * np.cos(th) + 0.4 * np.cos(3 * th)
        assert impl.pressure_drag_coefficient(th, cp) == pytest.approx(np.pi * b / 2, rel=1e-12)


def test_separated_model_matches_its_closed_form_and_base_suction_adds_drag():
    th = 2 * np.pi * np.arange(400_000) / 400_000
    for deg, cpb in ((82.0, -1.2), (100.0, -0.6), (120.0, -0.2)):
        ts = np.radians(deg)
        cd = impl.pressure_drag_coefficient(th, impl.cp_separated_model(th, ts, cpb))
        s = np.sin(ts)
        # front: (1/2) integral of (1 - 4 sin^2) cos over |theta| < theta_s; wake: -cp_b sin(theta_s)
        assert cd == pytest.approx(s - 4 * s ** 3 / 3 - cpb * s, abs=1e-4)
    cds = [impl.pressure_drag_coefficient(TH, impl.cp_separated_model(TH, np.radians(82), c)) for c in (-0.6, -1.0, -1.4)]
    assert cds[0] < cds[1] < cds[2]
    # a wake that is symmetric about the stagnation line makes no side force
    assert abs(impl.pressure_lift_coefficient(TH, impl.cp_separated_model(TH, np.radians(82), -1.2))) < 1e-12


# ---------------------------------------------------------------- correlations

def test_sphere_correlations_tend_to_stokes_in_the_creeping_limit():
    re = np.array([1e-4, 1e-3, 1e-2])
    npt.assert_allclose(impl.cd_sphere_morrison(re), 24 / re, rtol=0.01)
    npt.assert_allclose(impl.cd_sphere_schiller_naumann(re), 24 / re, rtol=0.01)
    npt.assert_allclose(impl.cd_sphere_stokes(re), 24 / re, rtol=1e-15)


def test_newton_regime_and_drag_crisis():
    re = np.logspace(3.3, 5, 20)
    cdm = impl.cd_sphere_morrison(re)
    assert np.all((cdm > 0.35) & (cdm < 0.5))                     # roughly constant, the Newton regime
    sn = impl.cd_sphere_schiller_naumann(np.array([999.0, 1001.0, 1e6]))
    assert abs(sn[0] - sn[1]) < 0.01 and sn[2] == 0.44            # nearly continuous; no crisis
    re_c, cd_min = impl.drag_crisis()
    assert 2e5 < re_c < 6e5 and cd_min < 0.15                     # Morrison: the drag crisis
    assert float(impl.cd_sphere_morrison(1e5)) > 3 * cd_min
    cyl = impl.cd_cylinder_white(np.array([1e2, 1e3, 1e4, 1e5]))
    assert np.all(np.diff(cyl) < 0) and 1.0 < cyl[-1] < 1.3 and cyl[0] > cyl[-1]


def test_strouhal_threshold_continuity_and_subcritical_value():
    assert np.isnan(impl.strouhal_cylinder(30.0))                  # steady wake: no shedding
    re = np.array([60.0, 100.0, 149.999, 150.001, 299.999, 300.001, 2000.0, 1e5])
    st = impl.strouhal_cylinder(re)
    assert st[2] == pytest.approx(st[3], abs=1e-5) and st[4] == pytest.approx(st[5], abs=1e-5)
    assert np.all(np.diff(st) >= -1e-9)
    assert 0.20 < st[-1] < 0.215
    # f = St V / d: at fixed Re (V d fixed) doubling V and halving d quadruples the frequency
    f1 = impl.shedding_frequency(10.0, 0.02, impl.RHO_AIR, impl.MU_AIR)
    f2 = impl.shedding_frequency(20.0, 0.01, impl.RHO_AIR, impl.MU_AIR)
    assert f2 == pytest.approx(4 * f1, rel=1e-12)


# ---------------------------------------------------------------- terminal velocity

def test_terminal_velocity_satisfies_the_force_balance():
    d, rp = 0.040, 2.7e-3 / (np.pi * 0.040 ** 3 / 6)               # a ping-pong ball
    V, re, cd = impl.terminal_velocity(d, rp, impl.RHO_AIR, impl.MU_AIR)
    W = (rp - impl.RHO_AIR) * impl.G0 * np.pi * d ** 3 / 6
    D = 0.5 * impl.RHO_AIR * V ** 2 * cd * np.pi * d ** 2 / 4
    assert D == pytest.approx(W, rel=1e-10)
    assert re == pytest.approx(impl.RHO_AIR * V * d / impl.MU_AIR, rel=1e-12)
    assert 8.0 < V < 11.0


def test_stokes_limit_of_terminal_velocity():
    args = (1e-5, 2500.0, impl.RHO_WATER, impl.MU_WATER)            # a 10 um glass bead in water
    Vs = impl.stokes_terminal_velocity(*args)
    V_exact, re, _ = impl.terminal_velocity(*args, cd=impl.cd_sphere_stokes)
    assert V_exact == pytest.approx(Vs, rel=1e-10)                   # 24/Re is Stokes's law
    V, re, _ = impl.terminal_velocity(*args)
    assert re < 1e-2 and V == pytest.approx(Vs, rel=0.01)


def test_time_integration_approaches_the_terminal_velocity():
    d, rp = 0.010, 7850.0                                            # a steel ball in water
    V, _, _ = impl.terminal_velocity(d, rp, impl.RHO_WATER, impl.MU_WATER)
    t, Vt, z = impl.simulate_fall(d, rp, impl.RHO_WATER, impl.MU_WATER, 3.0)
    assert np.all(np.diff(Vt) >= -1e-8)                              # accelerates monotonically (to solver tolerance)
    assert Vt[-1] == pytest.approx(V, rel=1e-4)
    # initial acceleration is g (1 - rho_f/rho_p): drag is zero at rest
    a0 = (Vt[1] - Vt[0]) / (t[1] - t[0])
    assert a0 == pytest.approx(impl.G0 * (1 - impl.RHO_WATER / rp), rel=0.02)
    assert z[-1] == pytest.approx(np.trapezoid(Vt, t), rel=1e-4)


def test_reynolds_similarity_through_the_archimedes_number():
    # c_d depends only on Re, so the terminal Re depends only on Ar = (rho_p - rho_f) rho_f g d^3 / mu^2
    d1, rp, rf, mu1 = 0.004, 2500.0, impl.RHO_WATER, impl.MU_WATER
    d2 = 0.008
    mu2 = mu1 * (d2 / d1) ** 1.5                                     # same Ar
    _, re1, cd1 = impl.terminal_velocity(d1, rp, rf, mu1)
    _, re2, cd2 = impl.terminal_velocity(d2, rp, rf, mu2)
    assert re2 == pytest.approx(re1, rel=1e-9) and cd2 == pytest.approx(cd1, rel=1e-9)


# ---------------------------------------------------------------- streamlining

def test_fairing_saves_an_order_of_magnitude():
    T, p, rho, mu = impl.isa_troposphere(3000.0)
    s = impl.strut_drag(0.025, 1.0, 60.0, rho, mu)
    assert s["re_d"] == pytest.approx(rho * 60.0 * 0.025 / mu)
    assert 5 < s["ratio"] < 25
    assert s["power_round"] == pytest.approx(s["drag_round"] * 60.0)
    assert s["drag_faired"] == pytest.approx(0.5 * rho * 60.0 ** 2 * s["cd_faired"] * 0.025, rel=1e-12)
    assert s["equivalent_wire_d"] < 0.25 * 0.025
    # the round strut's drag goes almost as V^2 (c_d varies only weakly with Re)
    s2 = impl.strut_drag(0.025, 1.0, 120.0, rho, mu)
    assert 3.8 < s2["drag_round"] / s["drag_round"] < 4.0
    # a thinner fairing has less drag per unit thickness... until friction on its longer chord dominates
    cds = [impl.fairing_cd(tc, 1e6) for tc in (0.1, 0.25, 0.5)]
    assert cds[1] < cds[0] and cds[1] < cds[2]


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "drag.png").stat().st_size > 10_000
