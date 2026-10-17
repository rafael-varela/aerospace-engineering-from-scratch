"""Public tests for P10.17 Propeller & Rotor Theory. They assert physics, not magic numbers:
the actuator disk conserves momentum and energy and its efficiency formula follows from them;
the static limit gives P = T^(3/2)/sqrt(2 rho A); the BEM residual is zero at the solution; BEM's
ideal limit (no drag, no swirl, no tip loss) reproduces momentum theory exactly for a blade designed
for uniform inflow, and never beats it otherwise; each loss lowers the efficiency; the
coefficients collapse at fixed J (thrust ~ n^2); the relaxed fixed-point iteration and the
bracketing solver agree; the element study converges at about second order.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

RHO = 1.225
R = 0.127
A = np.pi * R ** 2


@pytest.fixture(scope="module")
def blade():
    return impl.reference_propeller()


# ---------------------------------------------------------------- momentum theory

@pytest.mark.parametrize("T,V", [(1.2, 12.0), (5.0, 3.0), (20.0, 30.0), (3.0, 0.0)])
def test_actuator_disk_momentum_and_energy(T, V):
    ad = impl.actuator_disk(T, V, A)
    v = ad["v"]
    mdot = RHO * A * (V + v)
    # thrust = mass flow x velocity change; far wake at V + 2v
    assert mdot * (ad["v_wake"] - V) == pytest.approx(T, rel=1e-10)
    assert ad["v_wake"] == pytest.approx(V + 2 * v)
    # power = kinetic energy added to the stream per second
    assert ad["P"] == pytest.approx(0.5 * mdot * (ad["v_wake"] ** 2 - V ** 2), rel=1e-10)


def test_static_limit_and_power_law():
    T = 4.0
    ad = impl.actuator_disk(T, 0.0, A)
    assert ad["P"] == pytest.approx(impl.static_ideal_power(T, A), rel=1e-12)
    assert ad["v"] == pytest.approx(np.sqrt(T / (2 * RHO * A)), rel=1e-12)
    # doubling thrust costs 2^(3/2) the power; doubling the disk area divides it by sqrt(2)
    assert impl.static_ideal_power(2 * T, A) / impl.static_ideal_power(T, A) == pytest.approx(2 ** 1.5)
    assert impl.static_ideal_power(T, 2 * A) / impl.static_ideal_power(T, A) == pytest.approx(2 ** -0.5)


def test_ideal_efficiency_formula_matches_disk_and_limits():
    D, n = 2 * R, 100.0
    for T, V in [(1.2, 12.0), (3.0, 8.0), (0.2, 20.0)]:
        CT = T / (RHO * n ** 2 * D ** 4)
        J = V / (n * D)
        assert float(impl.ideal_efficiency(CT, J)) == pytest.approx(impl.actuator_disk(T, V, A)["eta"], rel=1e-12)
    assert float(impl.ideal_efficiency(1e-9, 0.5)) == pytest.approx(1.0, abs=1e-6)      # light loading -> 1
    etas = impl.ideal_efficiency(np.array([0.01, 0.05, 0.1]), 0.4)
    assert np.all(np.diff(etas) < 0)                                                     # heavier loading costs


# ---------------------------------------------------------------- BEM building blocks

def test_tip_loss_factor_limits():
    r = np.linspace(0.2 * R, R, 50)
    F = impl.prandtl_tip_loss(r, R, 2, np.full_like(r, 0.2))
    assert np.all((F >= 0) & (F <= 1 + 1e-12))
    assert F[-1] == pytest.approx(0.0, abs=1e-9)                                         # zero at the tip
    assert np.all(np.diff(F) <= 1e-12)                                                   # falls towards the tip
    F_many = impl.prandtl_tip_loss(r[:-1], R, 200, np.full(49, 0.2))
    assert np.all(F_many > 0.99)                                                         # B -> infinity: no loss


def test_bem_residual_vanishes_and_velocity_triangle_holds(blade):
    res = impl.bem(blade, 12.0, 100.0)
    Om = 2 * np.pi * 100.0
    Rr, k, kp, F = impl.element_residual(res["phi"], res["r"], blade.chord(res["r"]), blade.beta(res["r"]),
                                         blade, 12.0, Om)
    assert np.max(np.abs(Rr)) < 1e-9
    npt.assert_allclose(np.tan(res["phi"]), (12.0 + res["ua"]) / (Om * res["r"] - res["ut"]), rtol=1e-9)
    # annulus momentum balance: dT/dr = 4 pi r rho F (V + u_a) u_a
    npt.assert_allclose(res["dT"], 4 * np.pi * res["r"] * RHO * res["F"] * (12.0 + res["ua"]) * res["ua"],
                        rtol=1e-7, atol=1e-9)


def test_static_case_converges_and_coefficients_are_scale_free(blade):
    a = impl.bem(blade, 0.0, 80.0)
    b = impl.bem(blade, 0.0, 120.0)
    assert a["T"] > 0 and a["P"] > 0
    assert b["T"] / a["T"] == pytest.approx((120 / 80) ** 2, rel=1e-8)                  # T ~ n^2 at J = 0
    assert b["CT"] == pytest.approx(a["CT"], rel=1e-8) and b["CP"] == pytest.approx(a["CP"], rel=1e-8)
    # figure of merit below 1: the real blade needs more than the ideal power
    assert impl.static_ideal_power(a["T"], A) < a["P"]


def test_coefficients_collapse_at_fixed_advance_ratio(blade):
    J = 0.35
    D = 2 * R
    r1 = impl.bem(blade, J * 90.0 * D, 90.0)
    r2 = impl.bem(blade, J * 130.0 * D, 130.0)
    assert r1["CT"] == pytest.approx(r2["CT"], rel=1e-8)
    assert r1["CQ"] == pytest.approx(r2["CQ"], rel=1e-8)
    assert r1["CP"] == pytest.approx(2 * np.pi * r1["CQ"], rel=1e-10)
    assert r1["eta"] == pytest.approx(J * r1["CT"] / r1["CP"], rel=1e-10)


# ---------------------------------------------------------------- independent check: momentum theory

def test_uniform_inflow_blade_reproduces_momentum_theory():
    T, V, n = 1.2, 12.0, 100.0
    ub = impl.design_uniform_inflow_blade(T, V, n, R)
    res = impl.bem(ub, V, n, drag=False, swirl=False, tip_loss=False)
    A_ann = np.pi * (R ** 2 - ub.r_hub ** 2)
    ad = impl.actuator_disk(T, V, A_ann)
    assert res["T"] == pytest.approx(T, rel=1e-6)
    assert np.ptp(res["ua"]) < 1e-6 * ad["v"]                                            # uniform inflow
    npt.assert_allclose(res["ua"], ad["v"], rtol=1e-6)
    assert res["eta"] == pytest.approx(ad["eta"], rel=1e-7)                              # ideal efficiency
    assert res["P"] == pytest.approx(ad["P"], rel=1e-6)                                  # energy


def test_ideal_bem_energy_per_annulus():
    """With cd = 0 and no swirl, each annulus's power is exactly dT (V + u_a)."""
    blade = impl.reference_propeller()
    res = impl.bem(blade, 10.0, 100.0, drag=False, swirl=False, tip_loss=False)
    Om = 2 * np.pi * 100.0
    npt.assert_allclose(res["dQ"] * Om, res["dT"] * (10.0 + res["ua"]), rtol=1e-8, atol=1e-10)


def test_real_blade_never_beats_momentum_theory(blade):
    D = 2 * R
    for J in (0.15, 0.3, 0.45):
        for kw in ({}, dict(drag=False, swirl=False, tip_loss=False)):
            res = impl.bem(blade, J * 100.0 * D, 100.0, **kw)
            assert res["T"] > 0
            assert res["eta"] <= float(impl.ideal_efficiency(res["CT"], J)) + 1e-9


def test_each_loss_lowers_efficiency(blade):
    V, n = 12.0, 103.0
    steps = [dict(drag=False, swirl=False, tip_loss=False), dict(drag=False, swirl=True, tip_loss=False),
             dict(drag=False, swirl=True, tip_loss=True), dict(drag=True, swirl=True, tip_loss=True)]
    etas = [impl.bem(blade, V, n, **kw)["eta"] for kw in steps]
    assert np.all(np.diff(etas) < 0)
    Ts = [impl.bem(blade, V, n, **kw)["T"] for kw in steps[:3]]
    assert Ts[2] < Ts[1] < Ts[0]                                                         # swirl, tip loss cost thrust


def test_thrust_falls_with_advance_ratio_and_crosses_zero(blade):
    sw = impl.sweep_advance_ratio(blade, 100.0, np.linspace(0.05, 0.75, 15))
    assert np.all(np.diff(sw["CT"]) < 0)
    assert sw["CT"][0] > 0 and sw["CT"][-1] < 0                                          # windmilling at high J
    i = int(np.argmax(sw["eta"]))
    assert 0 < i < len(sw["eta"]) - 1                                                    # interior efficiency peak


def test_rpm_for_thrust_hits_the_target(blade):
    T_req = impl.uav_cruise_thrust()
    assert 0.5 < T_req < 3.0                                                             # a 1.2 kg UAV with L/D ~ 10
    n = impl.rpm_for_thrust(blade, T_req, 12.0)
    assert impl.bem(blade, 12.0, n)["T"] == pytest.approx(T_req, rel=1e-6)


# ---------------------------------------------------------------- numerics

def test_fixed_point_agrees_with_bracketing(blade):
    for kw in ({}, dict(swirl=False), dict(V=0.0)):
        V = kw.pop("V", 12.0)
        a = impl.bem(blade, V, 100.0, n_el=30, **kw)
        b = impl.bem_fixed_point(blade, V, 100.0, n_el=30, relax=0.2, **kw)
        assert b["residual"] < 1e-9
        assert b["T"] == pytest.approx(a["T"], rel=1e-7)
        assert b["P"] == pytest.approx(a["P"], rel=1e-7)


def test_element_convergence_is_about_second_order(blade):
    ref = impl.bem(blade, 12.0, 100.0, n_el=640)["T"]
    errs = [abs(impl.bem(blade, 12.0, 100.0, n_el=ne)["T"] - ref) for ne in (20, 40, 80)]
    orders = impl.observed_order(errs)
    assert np.all(orders > 1.5) and np.all(orders < 2.5)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "propeller_bem.png").stat().st_size > 10_000
