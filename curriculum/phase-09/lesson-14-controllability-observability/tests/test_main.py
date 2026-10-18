"""Public tests for P09.14 Controllability & Observability.

They assert structural properties, not magic numbers: the Kalman rank tests and the PBH eigenvalue test
agree on random controllable and constructed uncontrollable systems (the lesson's independent check);
observability is controllability of the dual; ranks survive a change of coordinates; the satellite loses
an axis when an orthogonal wheel fails but not when one pyramid wheel fails; the unobservable direction of
a two-axis star tracker is the boresight angle; the Lyapunov solver agrees with
``scipy.linalg.solve_continuous_lyapunov`` and the Gramians agree with quadrature of their defining
integrals (the other independent check); minimum-energy steering reaches its target with the energy the
Gramian predicts and beats a bang-bang manoeuvre.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import simpson, solve_ivp
from scipy.linalg import expm, solve_continuous_lyapunov

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(91401)
J = impl.satellite_inertia()


def random_system(n, m, controllable=True, seed=0):
    """A random (A, B); if not controllable, built as T blockdiag(A11, A22) T^-1 with B = T [B1; 0]."""
    rng = np.random.default_rng(seed)
    if controllable:
        return rng.standard_normal((n, n)), rng.standard_normal((n, m))
    k = n - 2
    A = np.zeros((n, n))
    A[:k, :k] = rng.standard_normal((k, k))
    A[:k, k:] = rng.standard_normal((k, n - k))
    A[k:, k:] = rng.standard_normal((n - k, n - k))
    B = np.zeros((n, m))
    B[:k] = rng.standard_normal((k, m))
    T = rng.standard_normal((n, n)) + 3.0 * np.eye(n)
    return T @ A @ np.linalg.inv(T), T @ B


def quad_gramian(A, B, T, n_pts):
    """Simpson quadrature (error ~ h^4) of int_0^T e^{At} B B^T e^{A^T t} dt with scipy's expm (independent of impl)."""
    t = np.linspace(0.0, T, n_pts)
    vals = np.array([expm(A * ti) @ B @ B.T @ expm(A.T * ti) for ti in t])
    return simpson(vals, x=t, axis=0)


# ---------------------------------------------------------------- rank tests

def test_ctrb_and_obsv_have_the_kalman_structure_and_are_dual():
    A, B = random_system(4, 2, seed=1)
    C = RNG.standard_normal((1, 4))
    Wc = impl.ctrb(A, B)
    assert Wc.shape == (4, 8)
    npt.assert_allclose(Wc[:, 4:6], A @ A @ B, rtol=1e-12)
    npt.assert_allclose(Wc[:, 6:8], A @ A @ A @ B, rtol=1e-12)
    O = impl.obsv(A, C)
    assert O.shape == (4, 4)
    npt.assert_allclose(O, impl.ctrb(A.T, C.T).T, rtol=1e-12)
    assert impl.is_observable(A, C) == impl.is_controllable(A.T, C.T)


@pytest.mark.parametrize("seed,controllable,m", [(2, True, 1), (3, True, 2), (4, True, 1),
                                                  (5, False, 1), (6, False, 2), (7, False, 1)])
def test_pbh_agrees_with_the_rank_test(seed, controllable, m):
    A, B = random_system(5, m, controllable, seed)
    rank_says = impl.is_controllable(A, B)
    pbh_says = len(impl.pbh_uncontrollable_modes(A, B)) == 0
    assert rank_says == pbh_says == controllable


def test_pbh_names_the_uncontrollable_eigenvalue():
    # x2 is driven by nothing: its eigenvalue -3 is the uncontrollable one
    A = np.array([[-1.0, 1.0], [0.0, -3.0]])
    B = np.array([[1.0], [0.0]])
    bad = impl.pbh_uncontrollable_modes(A, B)
    assert len(bad) == 1 and bad[0] == pytest.approx(-3.0)
    assert not impl.is_controllable(A, B)
    # dually, measuring only x1 of the transposed system leaves -3 unobservable
    assert impl.pbh_unobservable_modes(A.T, B.T)[0] == pytest.approx(-3.0)


def test_controllability_survives_a_change_of_coordinates():
    A, B = random_system(4, 1, seed=8)
    T = RNG.standard_normal((4, 4)) + 2.0 * np.eye(4)
    At, Bt = np.linalg.solve(T, A @ T), np.linalg.solve(T, B)
    assert impl.is_controllable(At, Bt) and impl.is_controllable(A, B)
    npt.assert_allclose(impl.ctrb(At, Bt), np.linalg.solve(T, impl.ctrb(A, B)), atol=1e-8)
    A2, B2 = random_system(4, 1, False, 9)
    assert not impl.is_controllable(np.linalg.solve(T, A2 @ T), np.linalg.solve(T, B2))


# ---------------------------------------------------------------- the satellite

def test_orthogonal_wheels_lose_an_axis_when_one_fails():
    A, B = impl.satellite_model(J, impl.wheel_axes("orthogonal"))
    assert impl.is_controllable(A, B)
    A, B = impl.satellite_model(J, impl.wheel_axes("orthogonal", failed=(2,)))
    assert not impl.is_controllable(A, B)
    assert impl.numerical_rank(impl.ctrb(A, B)) == 4
    bad = impl.pbh_uncontrollable_modes(A, B)
    assert len(bad) == 1 and abs(bad[0]) < 1e-12                          # the rigid-body mode at s = 0


def test_pyramid_survives_any_single_wheel_failure_but_not_two():
    for failed in ((), (0,), (1,), (2,), (3,)):
        A, B = impl.satellite_model(J, impl.wheel_axes("pyramid", failed))
        assert impl.is_controllable(A, B), failed
        assert impl.pbh_uncontrollable_modes(A, B) == []
    for failed in ((0, 1), (0, 2), (1, 3)):
        A, B = impl.satellite_model(J, impl.wheel_axes("pyramid", failed))
        assert not impl.is_controllable(A, B), failed


def test_two_axis_star_tracker_cannot_see_its_boresight_angle():
    A, _ = impl.satellite_model(J, impl.wheel_axes("orthogonal"))
    I6 = np.eye(6)
    assert impl.is_observable(A, I6[:3])                                 # 3-axis tracker alone
    assert impl.is_observable(A, I6[[0, 1, 2, 4, 5]])                    # x gyro failed: still observable
    C = I6[[0, 1, 3, 4, 5]]                                              # tracker x, y + all gyros
    assert not impl.is_observable(A, C)
    N = impl.unobservable_subspace(A, C)
    assert N.shape == (6, 1)
    npt.assert_allclose(np.abs(N[:, 0]), I6[2], atol=1e-12)              # the boresight (z) angle
    npt.assert_allclose(impl.obsv(A, C) @ N, 0.0, atol=1e-12)
    assert impl.unobservable_subspace(A, I6[3:]).shape == (6, 3)        # gyros only: all three angles


# ---------------------------------------------------------------- Lyapunov and Gramians

@pytest.mark.parametrize("n", [2, 4, 6])
def test_lyapunov_solver_matches_scipy_and_has_zero_residual(n):
    A = RNG.standard_normal((n, n))
    A = A - (np.max(np.linalg.eigvals(A).real) + 0.3) * np.eye(n)
    Bm = RNG.standard_normal((n, 2))
    Q = Bm @ Bm.T
    X = impl.lyap(A, Q)
    npt.assert_allclose(A @ X + X @ A.T + Q, 0.0, atol=1e-9 * np.abs(Q).max())
    npt.assert_allclose(X, solve_continuous_lyapunov(A, -Q), rtol=1e-8, atol=1e-10)


def test_infinite_gramian_equals_quadrature_of_its_integral():
    A, B, C = impl.short_period_model()
    W = impl.ctrb_gramian(A, B)
    # slowest decay of e^{At} B B^T e^{A^T t} is e^{-2 * 1.85 t}: by 12 s the tail is ~e^-44
    Wq = quad_gramian(A, B, 12.0, 4001)
    npt.assert_allclose(W, Wq, rtol=1e-7)
    assert np.all(np.linalg.eigvalsh(W) > 0)                             # controllable <=> W > 0
    M = impl.obsv_gramian(A, C)
    npt.assert_allclose(A.T @ M + M @ A + C.T @ C, 0.0, atol=1e-12)
    with pytest.raises(ValueError):
        impl.ctrb_gramian(*impl.satellite_model(J, impl.wheel_axes("orthogonal")))   # not stable


def test_finite_gramian_equals_quadrature_and_the_lyapunov_identity():
    A, B = impl.satellite_model(J, impl.wheel_axes("pyramid"))
    T = 8.0
    W = impl.ctrb_gramian_finite(A, B, T)
    npt.assert_allclose(W, quad_gramian(A, B, T, 801), rtol=1e-5, atol=1e-9)
    # d/dT W(T) = e^{AT} B B^T e^{A^T T}, and A W + W A^T + B B^T = e^{AT} B B^T e^{A^T T}
    E = expm(A * T)
    npt.assert_allclose(A @ W + W @ A.T + B @ B.T, E @ B @ B.T @ E.T, atol=1e-10)


def test_double_integrator_gramian_in_closed_form():
    Jz, T = 3.0, 5.0
    A = np.array([[0.0, 1.0], [0.0, 0.0]])
    B = np.array([[0.0], [1.0 / Jz]])
    W = impl.ctrb_gramian_finite(A, B, T)
    npt.assert_allclose(W, np.array([[T ** 3 / 3, T ** 2 / 2], [T ** 2 / 2, T]]) / Jz ** 2, rtol=1e-12)


def test_gramian_is_singular_exactly_when_a_wheel_axis_is_lost():
    A, B = impl.satellite_model(J, impl.wheel_axes("orthogonal", failed=(1,)))
    W = impl.ctrb_gramian_finite(A, B, 10.0)
    ev = np.linalg.eigvalsh(W)
    assert ev.min() < 1e-12 * ev.max()
    with pytest.raises(ValueError):
        impl.min_energy_input(A, B, np.zeros(6), np.ones(6), 10.0)


# ---------------------------------------------------------------- minimum-energy steering

def test_minimum_energy_input_reaches_the_target_with_the_predicted_energy():
    A, B = impl.satellite_model(J, impl.wheel_axes("pyramid", failed=(2,)))
    x0 = np.zeros(6)
    xf = np.array([np.radians(5.0), np.radians(-3.0), np.radians(8.0), 0.0, 0.0, 0.0])
    T = 30.0
    u, E = impl.min_energy_input(A, B, x0, xf, T)
    sol = solve_ivp(lambda t, x: A @ x + B @ u(t), (0.0, T), x0, rtol=1e-10, atol=1e-13)
    npt.assert_allclose(sol.y[:, -1], xf, atol=1e-8)
    t = np.linspace(0.0, T, 3001)
    Eq = np.trapezoid(np.sum(u(t) ** 2, axis=1), t)
    assert Eq == pytest.approx(E, rel=1e-6)


def test_rest_to_rest_slew_matches_the_closed_form_and_beats_bang_bang():
    Jz, dth = 3.0, np.radians(10.0)
    A = np.array([[0.0, 1.0], [0.0, 0.0]])
    B = np.array([[0.0], [1.0 / Jz]])
    energies = []
    for T in (10.0, 20.0):
        u, E = impl.min_energy_input(A, B, np.zeros(2), np.array([dth, 0.0]), T)
        t = np.linspace(0.0, T, 11)
        npt.assert_allclose(u(t)[:, 0], 6 * Jz * dth / T ** 2 * (1 - 2 * t / T), rtol=1e-9, atol=1e-14)
        assert E == pytest.approx(12 * Jz ** 2 * dth ** 2 / T ** 3, rel=1e-10)
        # bang-bang (+tau for T/2, -tau for T/2) does the same slew with tau = 4 J dth / T^2 and energy tau^2 T
        tau = 4 * Jz * dth / T ** 2
        assert tau ** 2 * T == pytest.approx(4.0 / 3.0 * E, rel=1e-10)
        energies.append(E)
    assert energies[0] / energies[1] == pytest.approx(8.0, rel=1e-10)    # energy ~ 1/T^3


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "controllability.png").stat().st_size > 1000
