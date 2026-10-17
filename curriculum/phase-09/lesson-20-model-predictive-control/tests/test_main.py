"""Tests for P09.20 Model Predictive Control.

Run against the starter (default) or the reference solution (AEFS_IMPL=solution).
The independent checks: the unconstrained MPC gain equals the finite-horizon LQR gain from the Riccati
recursion; with the DARE terminal cost it equals the infinite-horizon LQR gain (DARE checked against
``scipy.linalg.solve_discrete_are``); the constraints are never violated in closed loop.
"""
from __future__ import annotations

import numpy as np
import numpy.testing as npt
import pytest
from conftest import load_impl
from scipy.linalg import expm, solve_discrete_are

impl = load_impl(__file__)

RNG = np.random.default_rng(20250920)
_, _, AD, BD = impl.altitude_model()
Q, R = impl.altitude_weights()
LIM = impl.altitude_limits()
X0 = np.array([-20.0, 0.0, 0.0, 0.0])


def random_system(n: int, m: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    A = RNG.standard_normal((n, n)) * 0.6
    B = RNG.standard_normal((n, m))
    C = RNG.standard_normal((n, n))
    return A, B, C @ C.T + 0.1 * np.eye(n), np.eye(m) * RNG.uniform(0.5, 2.0)


# ---------------------------------------------------------------- model and Riccati

def test_altitude_model_discretisation_and_physics():
    A, B, Ad, Bd = impl.altitude_model()
    npt.assert_allclose(Ad, expm(A * 0.1), rtol=1e-10, atol=1e-12)
    lam = np.linalg.eigvals(A)
    assert np.sum(np.abs(lam) < 1e-9) == 2                              # altitude and trimmed flight-path integrators
    # a steady climb at gamma = theta needs no elevator: (0, g, g, 0) is an equilibrium of the gamma/theta/q rows
    g = 0.05
    npt.assert_allclose((A @ np.array([0.0, g, g, 0.0]))[1:], 0.0, atol=1e-15)


@pytest.mark.parametrize("n,m", [(2, 1), (4, 1), (4, 2)])
def test_dare_matches_scipy(n, m):
    A, B, Qr, Rr = random_system(n, m)
    P = impl.dare(A, B, Qr, Rr)
    ref = solve_discrete_are(A, B, Qr, Rr)
    npt.assert_allclose(P, ref, rtol=1e-8, atol=1e-10 * np.abs(ref).max())
    res = Qr + A.T @ P @ A - A.T @ P @ B @ np.linalg.solve(Rr + B.T @ P @ B, B.T @ P @ A) - P
    assert np.abs(res).max() < 1e-8 * np.abs(P).max()
    K = impl.lqr_gain(A, B, Qr, Rr, P)
    assert np.max(np.abs(np.linalg.eigvals(A - B @ K))) < 1.0


def test_riccati_recursion_converges_to_the_dare():
    Ks, P0 = impl.finite_horizon_lqr(AD, BD, Q, R, np.zeros((4, 4)), 400)
    P = impl.dare(AD, BD, Q, R)
    npt.assert_allclose(P0, P, rtol=1e-8)
    npt.assert_allclose(Ks[0], impl.lqr_gain(AD, BD, Q, R, P), rtol=1e-8)
    # the gains near the end of the horizon (short time-to-go) differ from the stationary one
    assert np.abs(Ks[-1] - Ks[0]).max() > 1e-2


# ---------------------------------------------------------------- condensed MPC

def test_prediction_matrices_reproduce_a_direct_simulation():
    N = 12
    Sx, Su = impl.prediction_matrices(AD, BD, N)
    assert Sx.shape == (48, 4) and Su.shape == (48, 12)
    U = RNG.standard_normal(N) * 0.1
    x = X0.copy()
    for k in range(N):
        x = AD @ x + BD[:, 0] * U[k]
        npt.assert_allclose((Sx @ X0 + Su @ U)[4 * k:4 * k + 4], x, rtol=1e-12, atol=1e-12)
    npt.assert_allclose(Su[:4, 1:], 0.0)                               # causality: x_1 does not depend on u_1...


def test_condensed_cost_equals_the_stagewise_sum():
    N = 8
    P = impl.dare(AD, BD, Q, R)
    H, F, Y = impl.condensed_cost(AD, BD, Q, R, P, N)
    for _ in range(5):
        U = RNG.standard_normal(N) * 0.05
        x0 = RNG.standard_normal(4) * [5, 0.05, 0.05, 0.1]
        J, x = 0.0, x0.copy()
        for k in range(N):
            J += x @ Q @ x + U[k] * R[0, 0] * U[k]
            x = AD @ x + BD[:, 0] * U[k]
        J += x @ P @ x
        npt.assert_allclose(0.5 * U @ H @ U + x0 @ F.T @ U + x0 @ Y @ x0, J, rtol=1e-10)
    assert np.all(np.linalg.eigvalsh(H) > 0)


@pytest.mark.parametrize("N", [1, 3, 10, 40])
def test_unconstrained_mpc_equals_finite_horizon_lqr(N):
    Km = impl.unconstrained_mpc_gain(AD, BD, Q, R, Q, N)
    Ks, _ = impl.finite_horizon_lqr(AD, BD, Q, R, Q, N)
    npt.assert_allclose(Km, Ks[0], rtol=1e-9, atol=1e-12)


@pytest.mark.parametrize("N", [1, 5, 25])
def test_dare_terminal_cost_makes_mpc_equal_infinite_horizon_lqr(N):
    P = impl.dare(AD, BD, Q, R)
    K_inf = impl.lqr_gain(AD, BD, Q, R, P)
    npt.assert_allclose(impl.unconstrained_mpc_gain(AD, BD, Q, R, P, N), K_inf, rtol=1e-8, atol=1e-11)
    npt.assert_allclose(impl.lqr_gain(AD, BD, Q, R, solve_discrete_are(AD, BD, Q, R)), K_inf, rtol=1e-8)


# ---------------------------------------------------------------- the QP solver

def test_qp_small_hand_case():
    """min (z1 - 1)^2 + (z2 - 2)^2 s.t. z1 + z2 <= 1, z >= 0: solution (0, 1), multiplier of the sum 2."""
    H = 2.0 * np.eye(2)
    f = np.array([-2.0, -4.0])
    G = np.array([[1.0, 1.0], [-1.0, 0.0], [0.0, -1.0]])
    w = np.array([1.0, 0.0, 0.0])
    z, lam, _ = impl.solve_qp(H, f, G, w)
    npt.assert_allclose(z, [0.0, 1.0], atol=1e-12)
    npt.assert_allclose(lam, [2.0, 0.0, 0.0], atol=1e-12)             # 2(z - c) + 2 (1, 1) = 0 at z = (0, 1)
    k = impl.kkt_residuals(H, f, G, w, z, lam)
    assert max(k.values()) < 1e-12


@pytest.mark.parametrize("seed", [1, 2, 3, 4, 5])
def test_qp_satisfies_kkt_on_random_problems(seed):
    rng = np.random.default_rng(seed)
    n, mc = 12, 30
    M = rng.standard_normal((n, n))
    H = M @ M.T + 0.5 * np.eye(n)
    f = rng.standard_normal(n) * 5
    G = rng.standard_normal((mc, n))
    w = rng.uniform(0.1, 1.0, mc)                                     # z = 0 is strictly feasible
    z, lam, _ = impl.solve_qp(H, f, G, w)
    k = impl.kkt_residuals(H, f, G, w, z, lam)
    assert k["stationarity"] < 1e-9 and k["primal"] < 1e-9 and k["dual"] == 0.0 and k["complementarity"] < 1e-9
    # optimality against feasible perturbations (a convex QP has no better feasible point)
    J = lambda zz: 0.5 * zz @ H @ zz + f @ zz
    for _ in range(50):
        zt = z + rng.standard_normal(n) * 0.05
        if np.all(G @ zt <= w):
            assert J(zt) >= J(z) - 1e-12


def test_qp_without_active_constraints_is_the_unconstrained_minimum_and_infeasible_raises():
    H = np.diag([1.0, 2.0])
    f = np.array([-0.1, 0.2])
    z, lam, it = impl.solve_qp(H, f, np.array([[1.0, 0.0]]), np.array([5.0]))
    npt.assert_allclose(z, [0.1, -0.1], atol=1e-15)
    assert lam[0] == 0.0 and it == 0
    with pytest.raises(ValueError):
        impl.solve_qp(H, f, np.array([[1.0, 0.0], [-1.0, 0.0]]), np.array([-1.0, -1.0]))     # z1 <= -1 and z1 >= 1


# ---------------------------------------------------------------- receding horizon on the altitude hold

@pytest.fixture(scope="module")
def mpc_run():
    return impl.simulate_mpc(X0, 150, N=20)


def test_closed_loop_never_violates_a_constraint(mpc_run):
    u, x = mpc_run["u"], mpc_run["x"]
    tol = 1e-9
    assert np.abs(u).max() <= LIM["u_max"] + tol
    du = np.diff(np.concatenate([[0.0], u]))
    assert np.abs(du).max() <= LIM["du_max"] + tol
    assert np.abs(LIM["V"] * x[:, 1]).max() <= LIM["climb_max"] + tol
    # and the limits are actually used (the problem is genuinely constrained)
    assert np.abs(LIM["V"] * x[:, 1]).max() > LIM["climb_max"] - 1e-6
    assert np.abs(du).max() > LIM["du_max"] - 1e-9


def test_mpc_reaches_the_commanded_altitude(mpc_run):
    x = mpc_run["x"]
    assert abs(x[-1, 0]) < 0.05
    assert x[:, 0].max() < 0.05                                         # essentially no overshoot


def test_unconstrained_lqr_would_break_the_limits():
    lq = impl.simulate_lqr(X0, 100)
    assert np.abs(lq["u"]).max() > 3 * LIM["u_max"]
    assert np.abs(LIM["V"] * lq["x"][:, 1]).max() > 3 * LIM["climb_max"]
    clipped = impl.simulate_lqr(X0, 100, clip=LIM)
    assert np.abs(LIM["V"] * clipped["x"][:, 1]).max() > 2 * LIM["climb_max"]   # clipping the actuator ignores the state limit


def test_with_loose_limits_mpc_is_the_lqr_loop():
    loose = {"u_max": 10.0, "du_max": 10.0, "climb_max": 100.0, "V": LIM["V"]}
    x0 = np.array([-2.0, 0.0, 0.0, 0.0])
    r = impl.simulate_mpc(x0, 40, N=10, limits=loose)
    lq = impl.simulate_lqr(x0, 40)
    npt.assert_allclose(r["x"], lq["x"], atol=1e-9)
    npt.assert_allclose(r["u"], lq["u"], atol=1e-10)


def test_too_short_a_horizon_loses_feasibility_and_says_so():
    """With N = 3 (0.3 s) the controller cannot see that the climb must be stopped within the elevator rate
    limit: a few samples in, no input sequence meets the climb-rate limit and the QP must report it."""
    with pytest.raises(ValueError):
        impl.simulate_mpc(X0, 30, N=3)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "mpc_altitude.png").stat().st_size > 1000
