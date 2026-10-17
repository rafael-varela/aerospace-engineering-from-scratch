"""Public tests for P03.07 Symplectic Integrators. They assert conservation, symmetry, symplecticity
and convergence order, with the closed-form Kepler solution and ``scipy.integrate.solve_ivp`` as oracles.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-03/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import solve_ivp

from conftest import load_impl

impl = load_impl(__file__)

MU = impl.MU_EARTH
A, E = 7000e3, 0.1


def kepler_rhs(t, y):
    return np.concatenate([y[2:], impl.kepler_acceleration(y[:2], MU)])


# ---------------------------------------------------------------- orbit bookkeeping and the oracle

def test_perigee_state_energy_momentum_and_period():
    q0, p0 = impl.perigee_state(A, E, MU)
    assert impl.energy(q0, p0, MU) == pytest.approx(-MU / (2 * A), rel=1e-12)
    assert impl.angular_momentum(q0, p0) == pytest.approx(np.sqrt(MU * A * (1 - E ** 2)), rel=1e-12)
    assert impl.orbital_period(q0, p0, MU) == pytest.approx(2 * np.pi * np.sqrt(A ** 3 / MU), rel=1e-12)
    assert np.linalg.norm(q0) == pytest.approx(A * (1 - E))
    a = impl.kepler_acceleration(q0, MU)
    assert np.linalg.norm(a) == pytest.approx(MU / np.linalg.norm(q0) ** 2, rel=1e-12)
    assert a @ q0 < 0                                            # points toward the centre
    with pytest.raises(ValueError):
        impl.orbital_period(q0, 2 * p0, MU)                      # escape speed: unbound
    with pytest.raises(ValueError):
        impl.perigee_state(A, 1.0, MU)
    # Array form of the invariants.
    Es = impl.energy(np.vstack([q0, q0]), np.vstack([p0, p0]), MU)
    assert Es.shape == (2,) and Es[0] == pytest.approx(-MU / (2 * A))


def test_solve_kepler_satisfies_keplers_equation():
    for e in (0.0, 0.1, 0.5, 0.9):
        for M in (0.0, 0.3, 2.0, 3.1, 5.9):
            Ecc = impl.solve_kepler(M, e)
            assert abs(Ecc - e * np.sin(Ecc) - M) < 1e-12


def test_kepler_exact_matches_solve_ivp_and_is_periodic():
    q0, p0 = impl.perigee_state(A, E, MU)
    P = impl.orbital_period(q0, p0, MU)
    ts = np.array([0.25, 0.5, 1.0, 1.5, 2.3]) * P
    q_ex, p_ex = impl.kepler_exact(ts, q0, p0, MU)
    sol = solve_ivp(kepler_rhs, (0.0, ts[-1]), np.concatenate([q0, p0]), method="DOP853", rtol=1e-13, atol=1e-6, t_eval=ts)
    npt.assert_allclose(q_ex, sol.y[:2].T, atol=1e-3)            # m, on a 7000 km orbit
    npt.assert_allclose(p_ex, sol.y[2:].T, atol=1e-6)            # m/s
    npt.assert_allclose(q_ex[2], q0, atol=1e-6)                  # back at perigee after one period
    q_half = impl.kepler_exact(0.5 * P, q0, p0, MU)[0][0]
    assert np.linalg.norm(q_half) == pytest.approx(A * (1 + E), rel=1e-10)   # apogee
    # A circular orbit keeps its radius; an unbound state is refused.
    qc, pc = impl.perigee_state(A, 0.0, MU)
    qs, _ = impl.kepler_exact(np.linspace(0, P, 7), qc, pc, MU)
    npt.assert_allclose(np.linalg.norm(qs, axis=1), A, rtol=1e-10)
    with pytest.raises(ValueError):
        impl.kepler_exact(100.0, q0, 2 * p0, MU)


# ---------------------------------------------------------------- one step of each method on the oscillator

def test_one_step_formulas_on_the_harmonic_oscillator():
    w, h = 2.0, 0.3
    acc = lambda q: -w ** 2 * q
    q, p = np.array([1.0, 0.0]), np.array([0.5, 0.0])
    q1, p1 = impl.euler_step(acc, q, p, h)
    npt.assert_allclose([q1[0], p1[0]], [q[0] + h * p[0], p[0] - h * w ** 2 * q[0]], rtol=1e-14)
    q1, p1 = impl.symplectic_euler_step(acc, q, p, h)
    p_new = p[0] - h * w ** 2 * q[0]
    npt.assert_allclose([q1[0], p1[0]], [q[0] + h * p_new, p_new], rtol=1e-14)
    q1, p1 = impl.stormer_verlet_step(acc, q, p, h)
    q_new = q[0] + h * p[0] - 0.5 * h ** 2 * w ** 2 * q[0]
    npt.assert_allclose([q1[0], p1[0]], [q_new, p[0] - 0.5 * h * w ** 2 * (q[0] + q_new)], rtol=1e-14)
    # RK4 on the oscillator is the degree-4 Taylor polynomial of the rotation.
    z = -1j * w * h                                            # y = q + i p/w obeys y' = -i w y
    R = 1 + z + z ** 2 / 2 + z ** 3 / 6 + z ** 4 / 24
    q1, p1 = impl.rk4_step(acc, q, p, h)
    assert q1[0] + 1j * p1[0] / w == pytest.approx((q[0] + 1j * p[0] / w) * R, rel=1e-12)
    # Yoshida's weights: w0 + 2 w1 = 1 and w0^3 + 2 w1^3 = 0 (the fourth-order conditions).
    assert impl.YOSHIDA_W0 + 2 * impl.YOSHIDA_W1 == pytest.approx(1.0)
    assert impl.YOSHIDA_W0 ** 3 + 2 * impl.YOSHIDA_W1 ** 3 == pytest.approx(0.0, abs=1e-14)
    # Second components stay zero: the force is linear and per-component.
    assert q1[1] == 0.0 and p1[1] == 0.0


# ---------------------------------------------------------------- observed convergence orders

@pytest.mark.parametrize("name, order, tol", [("symplectic_euler", 1.0, 0.2), ("verlet", 2.0, 0.15),
                                               ("rk4", 4.0, 0.3), ("yoshida4", 4.0, 0.3)])
def test_observed_order_of_each_method_over_one_and_a_half_orbits(name, order, tol):
    q0, p0 = impl.perigee_state(A, E, MU)
    P = impl.orbital_period(q0, p0, MU)
    ns = [300, 600, 1200]
    errs = impl.position_errors(impl.STEPPERS[name], q0, p0, 1.5 * P, ns, MU)
    assert np.all(np.diff(errs) < 0)
    assert abs(impl.convergence_order(1.5 * P / np.array(ns), errs) - order) < tol, errs


def test_symplectic_euler_looks_second_order_over_exactly_one_period():
    # Symplectic Euler is Verlet conjugated by a half-step shift; after a whole period the O(h)
    # shift cancels and only Verlet's O(h^2) error remains.
    q0, p0 = impl.perigee_state(A, E, MU)
    P = impl.orbital_period(q0, p0, MU)
    ns = [200, 400, 800]
    errs = impl.position_errors(impl.symplectic_euler_step, q0, p0, P, ns, MU)
    assert abs(impl.convergence_order(P / np.array(ns), errs) - 2.0) < 0.15, errs


def test_convergence_order_helper_and_propagate_bookkeeping():
    hs = np.array([0.4, 0.2, 0.1])
    assert impl.convergence_order(hs, 5.0 * hs ** 3) == pytest.approx(3.0, abs=1e-10)
    with pytest.raises(ValueError):
        impl.convergence_order(hs, np.array([1.0, 0.0, 1.0]))
    q0, p0 = impl.perigee_state(A, E, MU)
    t, qs, ps = impl.propagate(impl.stormer_verlet_step, lambda q: impl.kepler_acceleration(q, MU), q0, p0, 10.0, 5)
    assert qs.shape == (6, 2) and ps.shape == (6, 2) and t[-1] == pytest.approx(50.0)
    npt.assert_array_equal(qs[0], q0)
    with pytest.raises(ValueError):
        impl.propagate(impl.stormer_verlet_step, lambda q: q, q0, p0, 1.0, 0)


# ---------------------------------------------------------------- symplecticity, reversibility, invariants

def test_symplectic_methods_preserve_phase_area_and_rk4_does_not():
    w = 2 * np.pi
    acc = lambda q: -w ** 2 * q
    q, p = np.array([1.0, 0.0]), np.array([0.0, 0.0])
    h = 1.0 / w                                                   # h omega = 1: a coarse step
    for name in ("symplectic_euler", "verlet", "yoshida4"):
        M = impl.step_jacobian(impl.STEPPERS[name], acc, q, p, h)
        assert abs(np.linalg.det(M) - 1.0) < 1e-8, name          # the map is linear here: FD is exact
        assert impl.symplectic_defect(M) < 1e-7, name
    M = impl.step_jacobian(impl.rk4_step, acc, q, p, h)
    assert np.linalg.det(M) < 0.99                                # 0.9758: RK4 shrinks phase area
    assert impl.symplectic_defect(M) > 1e-2
    M = impl.step_jacobian(impl.euler_step, acc, q, p, h)
    assert np.linalg.det(M) > 1.01                                # explicit Euler inflates it
    # Nonlinear case: the Kepler problem in units mu = 1.
    acc1 = lambda q: impl.kepler_acceleration(q, 1.0)
    q1, p1 = impl.perigee_state(1.0 / (1 - 0.1), 0.1, mu=1.0)
    for name in ("symplectic_euler", "verlet", "yoshida4"):
        assert impl.symplectic_defect(impl.step_jacobian(impl.STEPPERS[name], acc1, q1, p1, 0.2)) < 1e-7, name
    assert impl.symplectic_defect(impl.step_jacobian(impl.rk4_step, acc1, q1, p1, 0.2)) > 1e-7
    assert impl.symplectic_defect(np.eye(4)) == 0.0


def test_verlet_and_yoshida_are_time_reversible_and_the_others_are_not():
    acc1 = lambda q: impl.kepler_acceleration(q, 1.0)
    q1, p1 = impl.perigee_state(1.0 / (1 - 0.1), 0.1, mu=1.0)
    for name in ("verlet", "yoshida4"):
        assert impl.reversibility_defect(impl.STEPPERS[name], acc1, q1, p1, 0.05) < 1e-12, name
    for name in ("euler", "symplectic_euler", "rk4"):
        assert impl.reversibility_defect(impl.STEPPERS[name], acc1, q1, p1, 0.05) > 1e-11, name


def test_verlet_conserves_the_modified_energy_exactly_on_the_oscillator():
    w = 2 * np.pi
    acc = lambda q: -w ** 2 * q
    h = 0.5 / w                                                   # h omega = 0.5
    _, qs, ps = impl.propagate(impl.stormer_verlet_step, acc, np.array([1.0, 0.0]), np.array([0.0, 0.0]), h, 400)
    H = 0.5 * ps[:, 0] ** 2 + 0.5 * w ** 2 * qs[:, 0] ** 2
    Hh = impl.modified_energy_oscillator(qs[:, 0], ps[:, 0], w, h)
    assert np.ptp(Hh) / Hh[0] < 1e-12                             # exactly conserved (to round-off)
    assert np.ptp(H) / H[0] == pytest.approx((h * w) ** 2 / 4, rel=1e-3)   # the true energy swings by (h w)^2 / 4
    assert np.max(np.abs(qs[:, 0])) < 1.0 + 1e-9                  # and the motion stays bounded
    # Beyond h omega = 2 the Verlet oscillator is unstable.
    _, qs_bad, _ = impl.propagate(impl.stormer_verlet_step, acc, np.array([1.0, 0.0]), np.array([0.0, 0.0]), 2.1 / w, 200)
    assert np.max(np.abs(qs_bad[:, 0])) > 1e3


def test_angular_momentum_is_exact_for_kick_drift_methods_only():
    q0, p0 = impl.perigee_state(A, E, MU)
    P = impl.orbital_period(q0, p0, MU)
    acc = lambda q: impl.kepler_acceleration(q, MU)
    L0 = impl.angular_momentum(q0, p0)
    for name in ("symplectic_euler", "verlet", "yoshida4"):
        _, qs, ps = impl.propagate(impl.STEPPERS[name], acc, q0, p0, P / 200, 10 * 200)
        assert np.max(np.abs(impl.angular_momentum(qs, ps) - L0)) / L0 < 1e-12, name
    _, qs, ps = impl.propagate(impl.rk4_step, acc, q0, p0, P / 200, 10 * 200)
    assert np.max(np.abs(impl.angular_momentum(qs, ps) - L0)) / L0 > 1e-9


def test_energy_error_is_bounded_for_verlet_and_secular_for_rk4_and_euler():
    q0, p0 = impl.perigee_state(A, E, MU)
    n_orb = 40
    t, dE = impl.energy_drift(impl.stormer_verlet_step, q0, p0, n_orb, 200, MU)
    assert np.max(np.abs(dE)) < 2e-4
    first, last = dE[: 200], dE[-200:]
    assert np.ptp(last) == pytest.approx(np.ptp(first), rel=0.05)   # same oscillation amplitude after 40 orbits
    assert abs(impl.drift_slope(t, dE)) * n_orb < 0.02 * np.max(np.abs(dE))   # no trend
    t, dE = impl.energy_drift(impl.rk4_step, q0, p0, n_orb, 200, MU)
    slope = impl.drift_slope(t, dE)
    assert slope < 0                                              # RK4 loses energy here
    assert dE[-1] == pytest.approx(slope * n_orb, rel=0.1)        # linear growth: the drift is secular
    assert abs(dE[-1]) > 10 * np.max(np.abs(dE[: 200]))          # far beyond the first orbit's swing
    t, dE = impl.energy_drift(impl.euler_step, q0, p0, 3, 200, MU)
    assert dE[-1] > 0.1 and np.all(np.diff(dE) >= -1e-12)         # explicit Euler pumps energy monotonically


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "symplectic_orbit.png").stat().st_size > 1000
