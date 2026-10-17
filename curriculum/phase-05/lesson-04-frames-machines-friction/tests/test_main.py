"""Public tests for P05.04 Frames, Machines & Friction. They assert statics, not magic numbers: circle
intersections lie on both circles; the kinematics keeps every pushrod at its length and returns to
neutral; the small-angle gearing equals the product of the arm ratios; each free body (stick, bellcrank,
elevator) is in force and moment equilibrium; the frictionless stick force from the member-by-member
chain equals the virtual-work force computed from the kinematics alone; with pin friction the driving
force exceeds the ideal one, the back-driven force falls below it, the band vanishes as mu -> 0, and the
power balance closes; the capstan law matches a numerical integration; the pulley ratio reduces to the
friction-circle closed form at 180 degrees and to 1 + 2 rho sin(theta/2) for small friction.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

LK = impl.Linkage()


def _cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


# ---------------------------------------------------------------- kinematics

def test_circle_intersections_lie_on_both_circles_and_refuse_disjoint_circles():
    rng = np.random.default_rng(5041)
    for _ in range(5):
        c0, c1 = rng.normal(size=2), rng.normal(size=2)
        d = np.linalg.norm(c1 - c0)
        r0 = rng.uniform(0.6, 1.0) * d
        r1 = rng.uniform(0.6, 1.0) * d
        P = impl.circle_intersections(c0, r0, c1, r1)
        npt.assert_allclose(np.linalg.norm(P - c0, axis=1), r0, rtol=1e-12)
        npt.assert_allclose(np.linalg.norm(P - c1, axis=1), r1, rtol=1e-12)
    with pytest.raises(ValueError):
        impl.circle_intersections([0, 0], 1.0, [3, 0], 1.0)


def test_kinematics_keeps_pushrod_lengths_and_neutral_is_zero():
    assert impl.solve_kinematics(LK, 0.0) == pytest.approx((0.0, 0.0), abs=1e-12)
    for phi in np.radians([-20.0, -7.0, 3.0, 15.0, 20.0]):
        psi, delta = impl.solve_kinematics(LK, phi)
        _, Ps = LK.stick_points(phi)
        B1, B2 = LK.bellcrank_points(psi)
        H = LK.horn_point(delta)
        assert np.linalg.norm(B1 - Ps) == pytest.approx(LK.L1, rel=1e-12)
        assert np.linalg.norm(H - B2) == pytest.approx(LK.L2, rel=1e-12)
        assert np.linalg.norm(H - LK.E0) == pytest.approx(LK.h, rel=1e-12)


def test_small_angle_gearing_is_the_product_of_arm_ratios_and_reverses_with_the_horn():
    dpsi, ddelta = impl.gearing(LK, 0.0)
    assert dpsi == pytest.approx(LK.a_s / LK.b1, rel=1e-7)
    assert ddelta == pytest.approx(-(LK.a_s / LK.b1) * (LK.b2 / LK.h), rel=1e-7)
    lk2 = impl.Linkage(b2=0.05)                    # halve the output arm: half the elevator travel
    assert impl.gearing(lk2, 0.0)[1] == pytest.approx(0.5 * ddelta, rel=1e-7)
    assert impl.mechanical_advantage(lk2, 0.0) == pytest.approx(2 * impl.mechanical_advantage(LK, 0.0), rel=1e-7)


# ---------------------------------------------------------------- statics

def test_every_free_body_is_in_equilibrium():
    for phi in np.radians([-18.0, 5.0, 20.0]):
        for mu, dirn in ((0.0, 1), (0.3, 1), (0.3, -1)):
            s = impl.linkage_forces(LK, phi, mu, direction=dirn)
            p = s["_parts"]
            Mf = s["M_friction"]
            scale = abs(s["P1"]) + 1.0
            # stick: grip force, pushrod 1, pivot reaction
            Fs = s["F"] * p["tg"] + s["P1"] * p["u1"] + s["R_stick"]
            Ms = _cross(p["grip"] - LK.S0, s["F"] * p["tg"]) + _cross(p["Ps"] - LK.S0, s["P1"] * p["u1"]) + Mf[0]
            # bellcrank: both pushrods and the pivot reaction
            Fb = -s["P1"] * p["u1"] + s["P2"] * p["u2"] + s["R_bell"]
            Mb = _cross(p["B1"] - LK.C0, -s["P1"] * p["u1"]) + _cross(p["B2"] - LK.C0, s["P2"] * p["u2"]) + Mf[1]
            # elevator: aerodynamic resultant, pushrod 2, hinge reaction
            Fe = p["F_aero"] - s["P2"] * p["u2"] + s["R_hinge"]
            Me = s["hinge_moment"] + _cross(p["H"] - LK.E0, -s["P2"] * p["u2"]) + Mf[2]
            npt.assert_allclose(np.concatenate([Fs, Fb, Fe]), 0.0, atol=1e-9 * scale)
            npt.assert_allclose([Ms, Mb, Me], 0.0, atol=1e-9 * scale)


def test_frictionless_chain_equals_virtual_work_over_the_travel():
    for phi in np.radians(np.linspace(-20.0, 20.0, 9)):
        F_chain = impl.linkage_forces(LK, phi)["F"]
        F_vw = impl.stick_force_virtual_work(LK, phi)
        assert F_chain == pytest.approx(F_vw, rel=1e-7, abs=1e-9)


def test_aerodynamic_resultant_reproduces_the_hinge_moment():
    s = impl.linkage_forces(LK, np.radians(12.0))
    p = s["_parts"]
    delta = s["delta"]
    cp = LK.E0 + LK.e_cp * np.array([np.cos(delta), np.sin(delta)])
    assert _cross(cp - LK.E0, p["F_aero"]) == pytest.approx(s["hinge_moment"], rel=1e-12)


def test_zero_hinge_moment_needs_no_force_and_force_is_linear_in_k():
    lk0 = impl.Linkage(k_hinge=0.0)
    s = impl.linkage_forces(lk0, np.radians(15.0), mu=0.3)
    assert s["F"] == pytest.approx(0.0, abs=1e-12) and s["P1"] == pytest.approx(0.0, abs=1e-12)
    F1 = impl.linkage_forces(impl.Linkage(k_hinge=50.0), np.radians(15.0))["F"]
    F2 = impl.linkage_forces(impl.Linkage(k_hinge=100.0), np.radians(15.0))["F"]
    assert F2 == pytest.approx(2 * F1, rel=1e-12)


# ---------------------------------------------------------------- pin friction

def test_friction_band_brackets_the_ideal_force_and_vanishes_with_mu():
    phi = np.radians(20.0)
    F0 = impl.linkage_forces(LK, phi)["F"]
    bands = []
    for mu in (0.4, 0.2, 0.1, 0.05):
        Fd = impl.linkage_forces(LK, phi, mu, direction=+1)["F"]
        Fb = impl.linkage_forces(LK, phi, mu, direction=-1)["F"]
        assert Fb < F0 < Fd
        bands.append(Fd - Fb)
    assert np.all(np.diff(bands) < 0)
    assert bands[-1] / bands[-2] == pytest.approx(0.5, rel=0.02)     # band ~ friction radius ~ mu for small mu
    assert impl.efficiency(LK, phi, 1e-9) == pytest.approx(1.0, abs=1e-8)
    assert impl.efficiency(LK, phi, 0.4) < impl.efficiency(LK, phi, 0.2) < 1.0


def test_friction_solution_closes_the_power_balance():
    for phi in np.radians([-15.0, 10.0, 20.0]):
        for dirn in (1, -1):
            s = impl.linkage_forces(LK, phi, 0.3, direction=dirn)
            assert abs(impl.friction_power_balance(LK, phi, s, 0.3, direction=dirn)) < 1e-6 * abs(s["F"]) * LK.Lg


def test_friction_radius_limits():
    assert impl.friction_radius(0.005, 0.0) == 0.0
    assert impl.friction_radius(0.005, 1e-4) == pytest.approx(0.005 * 1e-4, rel=1e-7)
    assert impl.friction_radius(0.005, 1e6) == pytest.approx(0.005, rel=1e-9)        # never beyond the pin


# ---------------------------------------------------------------- cables

def test_capstan_closed_form_matches_numerical_integration():
    for mu, wrap in ((0.1, np.pi / 6), (0.2, np.pi), (0.3, 3 * np.pi)):
        assert impl.capstan_numeric(mu, wrap) == pytest.approx(impl.capstan_ratio(mu, wrap), rel=1e-8)
    assert impl.capstan_ratio(0.0, 5.0) == 1.0
    assert impl.capstan_ratio(0.2, 2 * np.pi) == pytest.approx(impl.capstan_ratio(0.2, np.pi) ** 2)


def test_pulley_ratio_limits():
    mu_b, r_b, R = 0.15, 0.005, 0.0381
    rho = impl.friction_radius(r_b, mu_b) / R
    assert impl.pulley_ratio(mu_b, r_b, R, np.pi) == pytest.approx((1 + rho) / (1 - rho), rel=1e-12)
    for th in (np.pi / 6, np.pi / 2):
        k = impl.pulley_ratio(1e-4, r_b, R, th)
        rho_s = impl.friction_radius(r_b, 1e-4) / R
        assert k - 1 == pytest.approx(2 * rho_s * np.sin(th / 2), rel=1e-3)
    assert impl.pulley_ratio(0.0, r_b, R, 1.0) == pytest.approx(1.0)
    k = impl.pulley_ratio(mu_b, r_b, R, np.pi / 2)                # the root satisfies the moment balance
    assert (k - 1) * R == pytest.approx(impl.friction_radius(r_b, mu_b) * np.sqrt(k * k + 1 - 2 * k * np.cos(np.pi / 2)))
    with pytest.raises(ValueError):
        impl.pulley_ratio(10.0, 0.05, 0.04, 1.0)


def test_cable_run_is_the_product_of_its_elements():
    a = ("pulley", 0.15, 0.005, 0.0381, np.pi / 2)
    b = ("fairlead", 0.2, np.radians(3.0))
    assert impl.cable_run_ratio([a, b, a]) == pytest.approx(impl.pulley_ratio(*a[1:]) ** 2 * impl.capstan_ratio(0.2, np.radians(3.0)))
    assert impl.cable_run_ratio([]) == 1.0
    with pytest.raises(ValueError):
        impl.cable_run_ratio([("spring", 1.0)])


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "control_linkage.png").exists()
