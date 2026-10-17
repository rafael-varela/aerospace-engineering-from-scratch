"""Public tests for P06.09 3D Rigid Body Dynamics: Euler's Equations. They assert physics, not magic
numbers: principal axes as a proper eigenbasis of the inertia tensor; Euler's equations in principal and
non-principal axes and the power balance; conservation of kinetic energy, of |H| and of the inertial
angular-momentum vector in torque-free motion, and Poinsot's invariable plane; dH/dt = torque in the
inertial frame; the polhode classification; the exact Jacobi-elliptic solution satisfying Euler's equations
and matching the simulation on both sides of the separatrix, with any axis order and signs; its period and
flip interval against the simulated flips; and its limits (the small-wobble frequency of P04.04 and the
logarithmic growth of the flip interval near the separatrix).

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import math

import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

I_ABC = np.diag([1.0, 2.0, 3.0])
I_T = np.diag([2.2e-6, 6.0e-6, 7.4e-6])          # an illustrative T-handle


def _rand_rot(rng):
    q = rng.normal(size=4)
    return impl.dcm_from_quat(q / np.linalg.norm(q))


# ---------------------------------------------------------------- inertia and Euler's equations

def test_principal_axes_form_a_proper_eigenbasis():
    rng = np.random.default_rng(1)
    R0 = _rand_rot(rng)
    I = R0 @ np.diag([135.4, 208.3, 260.4]) @ R0.T
    mom, R = impl.principal_axes(I)
    npt.assert_allclose(mom, [135.4, 208.3, 260.4], rtol=1e-12)
    npt.assert_allclose(R.T @ I @ R, np.diag(mom), atol=1e-10)
    npt.assert_allclose(R.T @ R, np.eye(3), atol=1e-13)
    assert np.linalg.det(R) == pytest.approx(1.0, abs=1e-13)


def test_euler_rhs_component_form_and_sphere():
    w, tau = np.array([0.3, -1.1, 0.7]), np.array([0.05, 0.02, -0.04])
    d = np.diag(I_ABC)
    expected = [((d[1] - d[2]) * w[1] * w[2] + tau[0]) / d[0], ((d[2] - d[0]) * w[2] * w[0] + tau[1]) / d[1],
                ((d[0] - d[1]) * w[0] * w[1] + tau[2]) / d[2]]
    npt.assert_allclose(impl.euler_rhs(I_ABC, w, tau), expected, rtol=1e-14)
    npt.assert_allclose(impl.euler_rhs(2.0 * np.eye(3), w), np.zeros(3), atol=1e-16)   # no gyroscopic term


def test_euler_rhs_is_frame_covariant_and_balances_power():
    rng = np.random.default_rng(2)
    R = _rand_rot(rng)
    I = R @ I_ABC @ R.T                                   # the same body described in rotated axes
    w, tau = rng.normal(size=3), rng.normal(size=3)
    npt.assert_allclose(impl.euler_rhs(I, R @ w, R @ tau), R @ impl.euler_rhs(I_ABC, w, tau), atol=1e-13)
    wd = impl.euler_rhs(I, w, tau)
    assert w @ (I @ wd) == pytest.approx(w @ tau, rel=1e-12)             # dT/dt = w . tau


def test_energy_and_momentum_of_known_states():
    w = np.array([0.1, -0.2, 0.3])
    assert impl.kinetic_energy(I_ABC, w) == pytest.approx(0.5 * (0.01 + 2 * 0.04 + 3 * 0.09))
    npt.assert_allclose(impl.angular_momentum_body(I_ABC, w), [0.1, -0.4, 0.9], rtol=1e-15)
    npt.assert_allclose(impl.kinetic_energy(I_ABC, np.vstack([w, 2 * w])), [0.18, 0.72], rtol=1e-14)
    q90 = np.array([math.cos(math.pi / 4), 0.0, 0.0, math.sin(math.pi / 4)])  # body = reference yawed 90 deg
    npt.assert_allclose(impl.angular_momentum_inertial(I_ABC, np.array([1.0, 0, 0]), q90), [0, 1, 0], atol=1e-15)


# ---------------------------------------------------------------- torque-free motion

def test_torque_free_conserves_T_H_and_inertial_H():
    w0 = np.array([0.3, 1.0, 0.8])
    res = impl.simulate(I_ABC, w0, 20.0, n_samples=401)
    dr = impl.invariant_drift(I_ABC, res)
    assert dr["T"] < 1e-9 and dr["H"] < 1e-9 and dr["H_inertial"] < 1e-9
    assert np.ptp(res.w[:, 0]) > 0.1                                       # while the body rates change a lot
    npt.assert_allclose(np.linalg.norm(res.q, axis=1), 1.0, atol=1e-14)


def test_invariable_plane():
    res = impl.simulate(I_ABC, np.array([1.0, 0.5, 0.2]), 15.0, n_samples=301)
    assert impl.invariable_plane_residual(I_ABC, res) < 1e-9


def test_torque_changes_inertial_momentum_at_the_applied_rate():
    tau_i = np.array([0.0, 0.02, -0.01])                                   # a torque fixed in inertial space
    torque = lambda t, q, w: impl.dcm_from_quat(q) @ tau_i
    w0 = np.array([0.4, -0.3, 0.9])
    res = impl.simulate(I_ABC, w0, 10.0, torque=torque, n_samples=201)
    Hi = impl.angular_momentum_inertial(I_ABC, res.w, res.q)
    npt.assert_allclose(Hi, Hi[0] + np.outer(res.t, tau_i), atol=1e-9)


def test_non_principal_axes_give_the_same_motion():
    rng = np.random.default_rng(3)
    R = _rand_rot(rng)
    I_rot = R @ I_ABC @ R.T
    w0 = np.array([0.3, 1.0, 0.8])
    a = impl.simulate(I_ABC, w0, 8.0, n_samples=81)
    b = impl.simulate(I_rot, R @ w0, 8.0, n_samples=81)
    npt.assert_allclose(b.w, a.w @ R.T, atol=1e-9)
    assert impl.kinetic_energy(I_rot, b.w[-1]) == pytest.approx(impl.kinetic_energy(I_ABC, w0), rel=1e-10)


# ---------------------------------------------------------------- polhodes and the elliptic solution

def test_polhode_classification():
    assert impl.polhode_type(I_ABC, np.array([0.0, 0.1, 1.0])) == "major"
    assert impl.polhode_type(I_ABC, np.array([1.0, 0.1, 0.0])) == "minor"
    assert impl.polhode_type(I_ABC, np.array([0.0, 1.0, 0.0])) == "separatrix"
    # H^2 = 2 T I2 exactly: I1 w1^2 (I1 - I2) + I3 w3^2 (I3 - I2) = 0
    w1 = 1.0
    w3 = math.sqrt(1.0 * (2.0 - 1.0) / (3.0 * (3.0 - 2.0))) * w1
    assert impl.polhode_type(I_ABC, np.array([w1, 0.7, w3])) == "separatrix"
    with pytest.raises(ValueError):
        impl.elliptic_parameters(I_ABC, np.array([0.0, 1.0, 0.0]))


@pytest.mark.parametrize("w0", [[0.3, 1.0, 0.8], [1.0, 0.5, 0.2], [-0.3, 1.0, -0.8], [0.6, -0.9, 0.1], [-1.0, -0.4, 0.3]])
def test_elliptic_solution_matches_simulation(w0):
    w0 = np.array(w0)
    res = impl.simulate(I_ABC, w0, 12.0, n_samples=241)
    npt.assert_allclose(impl.elliptic_solution(I_ABC, w0, res.t), res.w, atol=1e-8)


def test_elliptic_solution_handles_any_axis_order():
    perm = [2, 0, 1]
    I_p = np.diag(np.diag(I_ABC)[perm])
    w0 = np.array([0.3, 1.0, 0.8])
    t = np.linspace(0, 10, 101)
    w_ref = impl.elliptic_solution(I_ABC, w0, t)
    npt.assert_allclose(impl.elliptic_solution(I_p, w0[perm], t), w_ref[:, perm], atol=1e-12)


def test_elliptic_solution_satisfies_euler_and_conservation():
    w0 = np.array([0.5, -0.8, 0.6])
    t = np.linspace(0, 6, 61)
    h = 1e-6
    w = impl.elliptic_solution(I_ABC, w0, t)
    wd = (impl.elliptic_solution(I_ABC, w0, t + h) - impl.elliptic_solution(I_ABC, w0, t - h)) / (2 * h)
    for k in range(t.size):
        npt.assert_allclose(wd[k], impl.euler_rhs(I_ABC, w[k]), atol=1e-7)
    npt.assert_allclose(impl.kinetic_energy(I_ABC, w), impl.kinetic_energy(I_ABC, w0), rtol=1e-12)
    npt.assert_allclose(np.linalg.norm(impl.angular_momentum_body(I_ABC, w), axis=1),
                        np.linalg.norm(I_ABC @ w0), rtol=1e-12)


def test_period_and_flip_interval_match_the_simulated_flips():
    spin = 2 * math.pi * 4.0
    w0 = np.array([0.01 * spin, spin, 0.0])                     # T-handle near its intermediate axis
    p = impl.elliptic_parameters(I_T, w0)
    assert p["flip_interval"] == pytest.approx(0.5 * p["period"], rel=1e-14)
    res = impl.simulate(I_T, w0, 2.2 * p["period"], n_samples=4001)
    fl = impl.flip_times(res.t, res.w[:, 1])
    assert fl.size >= 4
    npt.assert_allclose(np.diff(fl), p["flip_interval"], rtol=1e-4)
    assert impl.invariant_drift(I_T, res)["H_inertial"] < 1e-8


def test_small_amplitude_period_tends_to_the_wobble_frequency():
    spin = 3.0
    for eps in [1e-2, 1e-3]:
        p = impl.elliptic_parameters(I_ABC, np.array([0.0, eps * spin, spin]))
        f = impl.small_wobble_frequency(I_ABC, 2, spin)
        assert p["period"] == pytest.approx(2 * math.pi / f, rel=5 * eps)   # amplitude correction O(eps)
    assert p["m"] < 1e-5
    with pytest.raises(ValueError):
        impl.small_wobble_frequency(I_ABC, 1, spin)


def test_flip_interval_grows_logarithmically_near_the_separatrix():
    spin = 2 * math.pi * 4.0
    sig = impl.intermediate_growth_rate(I_T, spin)
    prev = None
    for eps in [1e-3, 1e-4, 1e-5, 1e-6]:
        p = impl.elliptic_parameters(I_T, np.array([eps * spin, spin, 0.0]))
        if prev is not None:
            # 2K/lambda ~ (2/sigma) ln(4/k') with k' proportional to eps: 2 ln(10)/sigma per decade
            assert p["flip_interval"] - prev == pytest.approx(2 * math.log(10) / sig, rel=1e-3)
        prev = p["flip_interval"]
        assert p["lam"] == pytest.approx(sig, rel=1e-3)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "tumbling_body.png").stat().st_size > 10_000
