"""Public tests for P06.10 Gyroscopes & Spin Stability. They assert physics, not magic numbers: the
nutation of an axisymmetric spinner (body rate (I3 - It) w3 / It, inertial coning rate H / It, constant
half-cone angle) against simulation; gyroscopic precession tau / (H sin theta) and the tilt of H under an
inertially fixed torque; spin averaging of a body-fixed torque (wander 2 tau / (I3 w3^2), a quarter when
the spin doubles); the energy at fixed |H| and the major/minor-axis rule; the energy-sink sign; and a Kane
damper whose closed-form linear nutation rate matches the nonlinear simulation, which conserves H, never
gains energy and ends in a flat spin with T = H^2 / (2 (It + J)).

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import math

import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

I_OB = np.diag([800.0, 800.0, 1200.0])          # oblate spinner
I_PR = np.diag([1500.0, 1500.0, 600.0])         # prolate spinner


def _angle_to(v, h):
    """Angle between v (one vector or rows) and h."""
    v = np.asarray(v, dtype=float)
    v = v / np.linalg.norm(v, axis=-1, keepdims=True)
    return np.arccos(np.clip(v @ (h / np.linalg.norm(h)), -1, 1))


# ---------------------------------------------------------------- nutation

def test_nutation_closed_forms():
    nu = impl.nutation(800.0, 1200.0, 1.0, 0.05)
    assert nu["H"] == pytest.approx(math.hypot(1200.0, 40.0), rel=1e-15)
    assert nu["body_rate"] == pytest.approx(0.5, rel=1e-15)
    assert nu["inertial_rate"] == pytest.approx(nu["H"] / 800.0, rel=1e-15)
    assert math.tan(nu["theta"]) == pytest.approx(40.0 / 1200.0, rel=1e-13)
    assert impl.nutation(1500.0, 600.0, 3.0, 0.1)["body_rate"] < 0          # prolate: against the spin
    assert impl.nutation(800.0, 1200.0, 1.0, 0.0)["theta"] == 0.0


@pytest.mark.parametrize("I, w3, wt", [(I_OB, 1.0, 0.05), (I_PR, 3.0, 0.2)])
def test_simulated_nutation_matches_closed_form(I, w3, wt):
    It, I3 = I[0, 0], I[2, 2]
    nu = impl.nutation(It, I3, w3, wt)
    res = impl.simulate_rigid(I, np.array([wt, 0.0, w3]), 30.0, n_samples=1501)
    Hi = impl.inertial_momentum(I, res)
    npt.assert_allclose(Hi - Hi[0], 0.0, atol=1e-9 * nu["H"])
    cone = _angle_to(impl.spin_axis_inertial(res), Hi[0])
    npt.assert_allclose(cone, nu["theta"], atol=1e-8)                        # constant half-cone angle
    assert impl.axis_coning_rate(res, Hi[0]) == pytest.approx(nu["inertial_rate"], rel=1e-7)
    phase = np.unwrap(np.arctan2(res.w[:, 1], res.w[:, 0]))                 # transverse rate in the body
    assert np.polyfit(res.t, phase, 1)[0] == pytest.approx(nu["body_rate"], rel=1e-7)


# ---------------------------------------------------------------- precession under torque

def test_gyroscope_precesses_at_tau_over_H_sin_theta():
    th, tau0 = math.radians(30.0), 2.0
    q0 = np.array([math.cos(th / 2), math.sin(th / 2), 0.0, 0.0])          # spin axis tilted 30 deg from z
    res = impl.simulate_rigid(I_OB, np.array([0.0, 0.0, 1.0]), 60.0, torque=impl.gyroscope_torque(I_OB, tau0),
                              q0=q0, n_samples=601)
    H = impl.inertial_momentum(I_OB, res)
    Hn = np.linalg.norm(H, axis=1)
    npt.assert_allclose(Hn, 1200.0, rtol=1e-10)                             # the torque is perpendicular to H
    npt.assert_allclose(np.arccos(H[:, 2] / Hn), th, atol=1e-9)             # and to z: the tilt is fixed
    az = np.unwrap(np.arctan2(H[:, 1], H[:, 0]))
    assert np.polyfit(res.t, az, 1)[0] == pytest.approx(impl.precession_rate(tau0, 1200.0, th), rel=1e-7)
    assert impl.precession_rate(tau0, 1200.0, math.pi / 2) == pytest.approx(tau0 / 1200.0)   # OpenStax 11.4


def test_fixed_torque_tilts_H_towards_itself():
    tau = np.array([0.0, 5.0, 0.0])
    res = impl.simulate_rigid(I_OB, np.array([0.0, 0.0, 1.0]), 40.0, torque=lambda t, q, w: impl.dcm_from_quat(q) @ tau,
                              n_samples=401)
    H = impl.inertial_momentum(I_OB, res)
    tilt = np.arctan2(H[:, 1], H[:, 2])
    npt.assert_allclose(tilt, impl.tilt_under_fixed_torque(5.0, 1200.0, res.t), atol=1e-9)
    npt.assert_allclose(H[:, 0], 0.0, atol=1e-8)                            # it moves along the torque, not across


def test_spin_averaging_bounds_the_wander_and_scales_as_one_over_spin_squared():
    out = []
    for spin in (1.0, 2.0):
        res = impl.simulate_rigid(I_OB, np.array([0.0, 0.0, spin]), 30.0, torque=lambda t, q, w: np.array([1.0, 0.0, 0.0]),
                                  n_samples=3001)
        H = impl.inertial_momentum(I_OB, res)
        out.append(_angle_to(H, H[0]).max())
        assert out[-1] == pytest.approx(impl.spin_averaged_wander(1.0, 1200.0, spin), rel=0.02)
    assert out[1] / out[0] == pytest.approx(0.25, rel=0.03)


# ---------------------------------------------------------------- energy, dissipation, the major-axis rule

def test_energy_at_fixed_momentum_is_lowest_about_the_major_axis():
    H = 10.0
    assert impl.energy_at_fixed_momentum(H, 4.0, 1.0, 0.0) == pytest.approx(H * H / 2.0)
    assert impl.energy_at_fixed_momentum(H, 4.0, 1.0, math.pi / 2) == pytest.approx(H * H / 8.0)
    th = np.linspace(0, math.pi / 2, 50)
    T = impl.energy_at_fixed_momentum(H, 4.0, 1.0, th)
    assert np.all(np.diff(T) < 0)                                            # prolate: T falls as theta grows
    T_ob = impl.energy_at_fixed_momentum(H, 1.0, 4.0, th)
    assert np.all(np.diff(T_ob) > 0)                                         # oblate: T is least at theta = 0
    assert impl.energy_sink_theta_rate(1.0, H, 4.0, 1.0, 0.3) > 0
    assert impl.energy_sink_theta_rate(1.0, H, 1.0, 4.0, 0.3) < 0
    fs = impl.flat_spin(H, 4.0)
    assert fs["rate"] == pytest.approx(2.5) and fs["T"] == pytest.approx(12.5)


def test_spin_stability_rule():
    I = np.diag([100.0, 200.0, 300.0])
    assert [impl.spin_stability(I, k, False) for k in range(3)] == ["stable", "unstable", "stable"]
    assert [impl.spin_stability(I, k, True) for k in range(3)] == ["unstable", "unstable", "stable"]
    E = impl.explorer_like()
    assert impl.spin_stability(E, 2, False) == "stable" and impl.spin_stability(E, 2, True) == "unstable"
    assert impl.spin_stability(I_OB, 2, True) == "stable"


# ---------------------------------------------------------------- the Kane damper

def test_kane_linear_matrix_limits():
    It, I3, J, Om = 1500.0, 600.0, 60.0, 3.0
    ev = np.sort_complex(np.linalg.eigvals(impl.kane_linear_matrix(It, I3, J, 0.0, Om)))
    npt.assert_allclose(ev, np.sort_complex(np.array([1j * (I3 - It) * Om / It, -1j * Om])), atol=1e-14)
    lam0 = impl.kane_nutation_mode(It, I3, J, 1e-9, Om)
    assert lam0.imag == pytest.approx((I3 - It) * Om / It, rel=1e-6)
    r1 = impl.kane_nutation_mode(It, I3, J, 1e-3, Om).real
    r2 = impl.kane_nutation_mode(It, I3, J, 2e-3, Om).real
    assert r2 == pytest.approx(2 * r1, rel=1e-3)                             # weak damping: rate proportional to c


def test_damper_destabilises_prolate_and_stabilises_oblate():
    for c in (1.0, 10.0, 100.0, 1000.0):
        assert impl.kane_nutation_mode(1500.0, 600.0, 60.0, c, 3.0).real > 0
        assert impl.kane_nutation_mode(800.0, 1200.0, 60.0, c, 1.0).real < 0


def test_damper_simulation_conserves_momentum_and_dissipates_energy():
    J, c = 60.0, 100.0
    res = impl.simulate_with_damper(I_PR, J, c, np.array([0.05, 0.0, 3.0]), np.array([0.0, 0.0, 3.0]), 100.0, n_samples=1001)
    H = impl.inertial_momentum(I_PR, res, J)
    npt.assert_allclose(H - H[0], 0.0, atol=1e-8 * np.linalg.norm(H[0]))
    E = impl.damper_energy(I_PR, J, res)
    assert np.all(np.diff(E) <= 1e-10 * E[0]) and E[-1] < E[0]
    res0 = impl.simulate_with_damper(I_PR, J, 0.0, np.array([0.05, 0.0, 3.0]), np.array([0.0, 0.0, 3.0]), 20.0, n_samples=201)
    E0 = impl.damper_energy(I_PR, J, res0)
    npt.assert_allclose(E0, E0[0], rtol=1e-10)                               # no viscosity, no loss


def test_damped_prolate_spinner_grows_at_the_linear_rate_and_ends_in_a_flat_spin():
    J, c, spin = 60.0, 100.0, 3.0
    lam = impl.kane_nutation_mode(1500.0, 600.0, J, c, spin)
    res = impl.simulate_with_damper(I_PR, J, c, np.array([1e-3, 0.0, spin]), np.array([0.0, 0.0, spin]), 600.0, n_samples=3001)
    amp = np.hypot(res.w[:, 0], res.w[:, 1])
    assert impl.fit_growth_rate(res.t, amp, 1e-2, 0.2) == pytest.approx(lam.real, rel=0.01)
    H = impl.inertial_momentum(I_PR, res, J)
    Hn = np.linalg.norm(H[0])
    E = impl.damper_energy(I_PR, J, res)
    assert E[-1] == pytest.approx(Hn ** 2 / (2 * (1500.0 + J)), rel=1e-6)    # the flat spin
    assert _angle_to(impl.spin_axis_inertial(res)[-1], H[0]) == pytest.approx(math.pi / 2, abs=1e-3)
    assert np.linalg.norm(res.w[-1]) == pytest.approx(Hn / (1500.0 + J), rel=1e-6)


def test_damped_oblate_spinner_loses_its_nutation():
    J, c = 60.0, 100.0
    lam = impl.kane_nutation_mode(800.0, 1200.0, J, c, 1.0)
    res = impl.simulate_with_damper(I_OB, J, c, np.array([0.05, 0.0, 1.0]), np.array([0.0, 0.0, 1.0]), 300.0, n_samples=1501)
    amp = np.hypot(res.w[:, 0], res.w[:, 1])
    assert amp[-1] < 0.01 * amp[0]
    assert impl.fit_growth_rate(res.t, amp, 1e-4, 0.03) == pytest.approx(lam.real, rel=0.03)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "spin_stability.png").stat().st_size > 10_000
