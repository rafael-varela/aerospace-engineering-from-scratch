"""Public tests for P06.06 Rigid Body Kinetics in 2D. They assert physics, not magic numbers: plane
kinetics of a body on an incline (rolling against sliding, the friction needed to roll), the touchdown
spin-up of a landing-gear wheel by impulse-momentum (final rolling state, spin-up time, slip distance)
and its friction energy computed two ways (integral of friction x slip velocity, and the loss of kinetic
energy), the independence of that energy from the friction coefficient and the load history, the limits
of a very heavy aircraft and of pre-spin, and a hybrid sliding/rolling/locked time integration that must
reproduce the closed forms and close its energy audit at every sample.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import math

import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

G = 9.80665
WHEEL = impl.Wheel()                     # the illustrative main-gear wheel (given)
I_W, R_W, M_W = WHEEL.I, WHEEL.R, WHEEL.M


# ---------------------------------------------------------------- warm-up: incline

def test_incline_rolling_acceleration_and_ordering_of_shapes():
    th, m, R = math.radians(12.0), 4.0, 0.2
    accs = []
    for k2 in [1.0, 0.5, 0.4]:                                  # hoop, disk, sphere: I = k2 m R^2
        r = impl.incline_motion(th, m, k2 * m * R * R, R, mu_s=0.9, mu_k=0.7)
        assert r["rolls"]
        assert r["a"] == pytest.approx(G * math.sin(th) / (1.0 + k2), rel=1e-13)
        assert r["alpha"] == pytest.approx(r["a"] / R, rel=1e-13)               # rolling constraint
        # Newton for the centre and for rotation about it, with the returned friction
        assert m * r["a"] == pytest.approx(m * G * math.sin(th) - r["friction"], rel=1e-12)
        assert k2 * m * R * R * r["alpha"] == pytest.approx(r["friction"] * R, rel=1e-12)
        accs.append(r["a"])
    assert accs[0] < accs[1] < accs[2]                          # the hoop is slowest, the sphere fastest
    assert impl.incline_motion(th, m, 0.5 * m * R * R, R, 0.9, 0.7)["normal"] == pytest.approx(m * G * math.cos(th))


def test_incline_slides_below_the_friction_needed_to_roll():
    th, m, R = math.radians(35.0), 2.0, 0.1
    I = 0.5 * m * R * R
    mu_min = impl.min_friction_to_roll(th, m, I, R)
    assert mu_min == pytest.approx(math.tan(th) / 3.0, rel=1e-13)              # disk: tan(theta)/(1 + 2)
    assert impl.incline_motion(th, m, I, R, mu_s=1.001 * mu_min, mu_k=0.1)["rolls"]
    r = impl.incline_motion(th, m, I, R, mu_s=0.999 * mu_min, mu_k=0.1)
    assert not r["rolls"]
    assert r["a"] == pytest.approx(G * (math.sin(th) - 0.1 * math.cos(th)), rel=1e-13)
    assert r["alpha"] == pytest.approx(0.1 * m * G * math.cos(th) * R / I, rel=1e-13)
    assert r["a"] > r["alpha"] * R                               # the contact point slips forward


def test_frictionless_incline_does_not_spin():
    r = impl.incline_motion(math.radians(20.0), 1.0, 0.01, 0.2, mu_s=0.0, mu_k=0.0)
    assert not r["rolls"] and r["alpha"] == 0.0
    assert r["a"] == pytest.approx(G * math.sin(math.radians(20.0)), rel=1e-14)


# ---------------------------------------------------------------- closed-form spin-up

def test_effective_mass_limits():
    assert impl.effective_mass(1e12, 3.0, 0.5) == pytest.approx(3.0 / 0.25, rel=1e-9)   # heavy aircraft: I/R^2
    assert impl.effective_mass(7.0, 1e12, 0.5) == pytest.approx(7.0, rel=1e-9)          # very heavy wheel: M
    assert impl.effective_mass(M_W, I_W, R_W) < min(M_W, I_W / R_W ** 2)


def test_spin_up_ends_rolling_with_consistent_impulse():
    cf = impl.spin_up_closed_form(70.0, R_W, I_W, M_W, 0.55, 120e3)
    assert cf["V_final"] == pytest.approx(cf["omega_final"] * R_W, rel=1e-13)        # rolling at the end
    assert M_W * (70.0 - cf["V_final"]) == pytest.approx(cf["impulse"], rel=1e-10)   # linear impulse on the axle
    assert I_W * cf["omega_final"] == pytest.approx(cf["impulse"] * R_W, rel=1e-12)  # angular impulse on the wheel
    assert cf["t_spin"] == pytest.approx(cf["impulse"] / (0.55 * 120e3), rel=1e-13)
    assert cf["slip_distance"] == pytest.approx(0.5 * 70.0 * cf["t_spin"], rel=1e-12)  # slip falls linearly
    assert cf["peak_drag"] == pytest.approx(0.55 * 120e3)


def test_friction_energy_two_routes_agree():
    rng = np.random.default_rng(6)
    for _ in range(20):
        V0, R, I, M = rng.uniform(20, 90), rng.uniform(0.2, 0.7), rng.uniform(1, 40), rng.uniform(500, 3e4)
        w0 = rng.uniform(0, 0.8) * V0 / R
        cf = impl.spin_up_closed_form(V0, R, I, M, rng.uniform(0.2, 0.9), rng.uniform(1e4, 2e5), rng.uniform(0, 0.2), w0)
        assert cf["energy_friction"] == pytest.approx(impl.kinetic_energy_loss(V0, R, I, M, w0), rel=1e-9)


def test_friction_energy_independent_of_mu_and_load_history():
    E = [impl.spin_up_closed_form(70.0, R_W, I_W, M_W, mu, N, tr)["energy_friction"]
         for mu in (0.3, 0.55, 0.8) for N in (60e3, 120e3) for tr in (0.0, 0.05, 0.3)]
    npt.assert_allclose(E, E[0], rtol=1e-13)
    t1 = impl.spin_up_closed_form(70.0, R_W, I_W, M_W, 0.4, 100e3)["t_spin"]
    t2 = impl.spin_up_closed_form(70.0, R_W, I_W, M_W, 0.4, 200e3)["t_spin"]
    assert t2 == pytest.approx(t1 / 2, rel=1e-13)                # twice the load, half the spin-up time


def test_heavy_aircraft_limit_energy_equals_wheel_kinetic_energy():
    cf = impl.spin_up_closed_form(70.0, R_W, I_W, 1e12, 0.55, 120e3)
    assert cf["V_final"] == pytest.approx(70.0, rel=1e-9)
    assert cf["energy_friction"] == pytest.approx(0.5 * I_W * (70.0 / R_W) ** 2, rel=1e-9)


def test_pre_spin_reduces_energy_quadratically():
    E0 = impl.spin_up_closed_form(70.0, R_W, I_W, M_W, 0.55, 120e3, 0.1)["energy_friction"]
    for f in (0.25, 0.5, 0.9):
        E = impl.spin_up_closed_form(70.0, R_W, I_W, M_W, 0.55, 120e3, 0.1, omega0=f * 70.0 / R_W)["energy_friction"]
        assert E == pytest.approx(E0 * (1 - f) ** 2, rel=1e-12)
    with pytest.raises(ValueError):
        impl.spin_up_closed_form(70.0, R_W, I_W, M_W, 0.55, 120e3, 0.0, omega0=1.1 * 70.0 / R_W)


def test_ramp_branches_meet_and_tend_to_the_constant_load():
    mu, N = 0.55, 120e3
    J = impl.effective_mass(M_W, I_W, R_W) * 70.0
    tr_star = 2 * J / (mu * N)                                    # spin-up ends exactly at the end of the ramp
    a = impl.spin_up_closed_form(70.0, R_W, I_W, M_W, mu, N, tr_star * (1 - 1e-9))
    b = impl.spin_up_closed_form(70.0, R_W, I_W, M_W, mu, N, tr_star * (1 + 1e-9))
    assert a["t_spin"] == pytest.approx(b["t_spin"], rel=1e-7)
    assert a["slip_distance"] == pytest.approx(b["slip_distance"], rel=1e-7)
    assert a["t_spin"] == pytest.approx(tr_star, rel=1e-7)
    c0 = impl.spin_up_closed_form(70.0, R_W, I_W, M_W, mu, N, 0.0)
    c1 = impl.spin_up_closed_form(70.0, R_W, I_W, M_W, mu, N, 1e-7)
    assert c1["t_spin"] == pytest.approx(c0["t_spin"], rel=1e-5)
    assert c1["slip_distance"] == pytest.approx(c0["slip_distance"], rel=1e-5)


# ---------------------------------------------------------------- time integration

@pytest.mark.parametrize("t_ramp", [0.0, 0.03, 0.15])            # constant load, spin-up after and during the ramp
def test_simulation_reproduces_closed_form(t_ramp):
    cf = impl.spin_up_closed_form(70.0, R_W, I_W, M_W, 0.55, 120e3, t_ramp)
    res = impl.simulate_touchdown(WHEEL, 70.0, 0.55, 0.8, impl.ramp_load(120e3, t_ramp), 0.25, n_samples=501)
    assert res.t_spin == pytest.approx(cf["t_spin"], rel=1e-8)
    assert res.V[-1] == pytest.approx(cf["V_final"], rel=1e-10)
    assert res.omega[-1] == pytest.approx(cf["omega_final"], rel=1e-9)
    assert res.E_friction[-1] == pytest.approx(cf["energy_friction"], rel=1e-8)
    assert res.slip_distance[-1] == pytest.approx(cf["slip_distance"], rel=1e-8)


def test_energy_audit_closes_and_friction_heat_equals_kinetic_energy_loss():
    res = impl.simulate_touchdown(WHEEL, 70.0, 0.55, 0.8, impl.ramp_load(120e3, 0.1), 0.25, n_samples=501)
    au = impl.energy_audit(res, WHEEL)
    assert au["max_rel_residual"] < 1e-9
    assert np.all(np.diff(res.E_friction) >= -1e-9)              # dissipation never runs backwards
    assert res.E_friction[-1] == pytest.approx(au["KE"][0] - au["KE"][-1], rel=1e-9)
    assert np.all(res.E_brake == 0.0)


def test_rolling_after_spin_up_without_brake():
    res = impl.simulate_touchdown(WHEEL, 70.0, 0.55, 0.8, impl.ramp_load(120e3, 0.0), 0.2, n_samples=401)
    roll = res.rolling
    assert roll[-1] and not roll[1]
    assert res.transitions[0][1] == "sliding" and res.transitions[1][1] == "rolling"
    npt.assert_allclose(res.slip[roll], 0.0, atol=1e-8)
    npt.assert_allclose(res.friction[roll], 0.0, atol=1e-12)                 # free rolling needs no friction
    assert np.ptp(res.V[roll]) < 1e-9
    assert np.all(res.slip[~roll][1:] > 0)                                    # forward slip until spin-up


def test_brake_below_limit_keeps_rolling_with_predicted_deceleration():
    N = 120e3
    Tmax = impl.max_brake_torque(0.8, N, I_W, M_W, R_W)
    assert impl.rolling_friction_for_brake(Tmax, I_W, M_W, R_W) == pytest.approx(0.8 * N, rel=1e-13)
    Tb = 0.9 * Tmax
    res = impl.simulate_touchdown(WHEEL, 70.0, 0.55, 0.8, impl.ramp_load(N, 0.0), 0.3,
                                  brake_torque=lambda t: Tb if t >= 0.1 else 0.0, n_samples=601)
    assert [p for _, p in res.transitions] == ["sliding", "rolling"]
    late = res.t > 0.12
    F = impl.rolling_friction_for_brake(Tb, I_W, M_W, R_W)
    dec = -np.diff(res.V[late]) / np.diff(res.t[late])
    npt.assert_allclose(dec, F / M_W, rtol=1e-6)
    npt.assert_allclose(res.slip[late], 0.0, atol=1e-8)
    assert impl.energy_audit(res, WHEEL)["max_rel_residual"] < 1e-9
    assert res.E_brake[-1] > 0


def test_brake_above_limit_skids_then_locks():
    N = 120e3
    Tb = 1.3 * impl.max_brake_torque(0.8, N, I_W, M_W, R_W)
    res = impl.simulate_touchdown(WHEEL, 70.0, 0.55, 0.8, impl.ramp_load(N, 0.0), 0.4,
                                  brake_torque=lambda t: Tb if t >= 0.1 else 0.0, n_samples=801)
    assert [p for _, p in res.transitions] == ["sliding", "rolling", "sliding", "locked"]
    t_lock = res.transitions[3][0]
    # after the skid starts, omega falls at (mu_k N R - T_b)/I until it reaches zero
    assert t_lock == pytest.approx(0.1 + res.omega[np.searchsorted(res.t, 0.1)] * I_W / (Tb - 0.55 * N * R_W), rel=1e-3)
    locked = res.t > t_lock + 1e-3
    npt.assert_allclose(res.omega[locked], 0.0, atol=1e-9)
    dec = -np.diff(res.V[locked]) / np.diff(res.t[locked])
    npt.assert_allclose(dec, 0.55 * N / M_W, rtol=1e-6)                     # a locked tyre drags at mu_k N
    assert 0.55 * N < impl.rolling_friction_for_brake(0.9 * Tb / 1.3, I_W, M_W, R_W)   # less than rolling braking
    assert impl.energy_audit(res, WHEEL)["max_rel_residual"] < 1e-9


def test_overspun_wheel_spins_down_and_pushes_the_axle_forward():
    V0, w0 = 50.0, 1.2 * 50.0 / R_W
    res = impl.simulate_touchdown(WHEEL, V0, 0.55, 0.8, impl.ramp_load(80e3, 0.0), 0.2, omega0=w0, n_samples=401)
    Me = impl.effective_mass(M_W, I_W, R_W)
    J = Me * (V0 - w0 * R_W)                                     # negative: friction acts forward on the axle
    assert res.V[-1] == pytest.approx(V0 - J / M_W, rel=1e-10)
    assert res.V[-1] > V0
    assert np.all(res.friction[1:][~res.rolling[1:]] < 0)
    assert res.E_friction[-1] == pytest.approx(0.5 * Me * (V0 - w0 * R_W) ** 2, rel=1e-8)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "wheel_spinup.png").stat().st_size > 10_000
