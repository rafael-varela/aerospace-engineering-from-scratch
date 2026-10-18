"""Tests for P09.16 LQR.

Run against the starter (default) or the reference solution (AEFS_IMPL=solution).
The independent checks: the Riccati solution against ``scipy.linalg.solve_continuous_are``; the residual
of the ARE; the Hamiltonian and Newton-Kleinman solutions agree with each other.
"""
from __future__ import annotations

import numpy as np
import numpy.testing as npt
import pytest
from conftest import load_impl
from scipy.linalg import expm, solve_continuous_are, solve_continuous_lyapunov

impl = load_impl(__file__)

RNG = np.random.default_rng(20250916)


def random_system(n: int, m: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """A random (generically controllable) system, possibly open-loop unstable, with Q > 0 and R > 0."""
    A = RNG.standard_normal((n, n))
    B = RNG.standard_normal((n, m))
    C = RNG.standard_normal((n, n))
    Q = C @ C.T + 0.1 * np.eye(n)
    D = RNG.standard_normal((m, m))
    R = D @ D.T + 0.5 * np.eye(m)
    return A, B, Q, R


SAT_A, SAT_B = impl.satellite_model(10.0)
SAT_Q = np.diag([1.0 / np.radians(5.0) ** 2, 1.0 / np.radians(0.5) ** 2])   # Bryson: 5 deg, 0.5 deg/s
SAT_R = np.array([[1.0 / 0.02 ** 2]])                                     # Bryson: 0.02 N m


# ---------------------------------------------------------------- Lyapunov

@pytest.mark.parametrize("n", [2, 3, 5])
def test_lyapunov_solve_residual_and_scipy(n):
    A = RNG.standard_normal((n, n)) - 3.0 * np.eye(n)        # shifted to be Hurwitz
    C = RNG.standard_normal((n, n))
    Q = C @ C.T
    X = impl.lyapunov_solve(A, Q)
    npt.assert_allclose(A.T @ X + X @ A + Q, 0.0, atol=1e-10 * np.abs(Q).max())
    # scipy solves a X + X a^H = q, so pass a = A^T, q = -Q
    npt.assert_allclose(X, solve_continuous_lyapunov(A.T, -Q), rtol=1e-9, atol=1e-12)
    npt.assert_allclose(X, X.T, atol=1e-14 * np.abs(X).max())
    assert np.all(np.linalg.eigvalsh(X) > -1e-10)             # Hurwitz A and Q >= 0 give X >= 0


def test_lyapunov_singular_case_raises():
    with pytest.raises(ValueError):
        impl.lyapunov_solve(np.array([[1.0, 0.0], [0.0, -1.0]]), np.eye(2))   # eigenvalues sum to 0


# ---------------------------------------------------------------- Riccati: the independent checks

@pytest.mark.parametrize("n,m", [(2, 1), (3, 1), (4, 2), (6, 3)])
def test_hamiltonian_care_matches_scipy_and_has_small_residual(n, m):
    A, B, Q, R = random_system(n, m)
    P = impl.care_hamiltonian(A, B, Q, R)
    ref = solve_continuous_are(A, B, Q, R)
    npt.assert_allclose(P, ref, rtol=1e-8, atol=1e-9 * np.abs(ref).max())
    assert impl.are_residual(A, B, Q, R, P) < 1e-10
    assert np.all(np.linalg.eigvalsh(P) > 0)                    # Q > 0 makes P positive definite
    K = np.linalg.solve(R, B.T @ P)
    assert np.all(np.linalg.eigvals(A - B @ K).real < 0)          # the stabilising solution


@pytest.mark.parametrize("n,m", [(2, 1), (3, 2), (5, 1)])
def test_newton_kleinman_agrees_with_hamiltonian_and_scipy(n, m):
    A, B, Q, R = random_system(n, m)
    P_n, iters = impl.care_newton_kleinman(A, B, Q, R)
    P_h = impl.care_hamiltonian(A, B, Q, R)
    npt.assert_allclose(P_n, P_h, rtol=1e-8, atol=1e-9 * np.abs(P_h).max())
    npt.assert_allclose(P_n, solve_continuous_are(A, B, Q, R), rtol=1e-8, atol=1e-9 * np.abs(P_h).max())
    assert impl.are_residual(A, B, Q, R, P_n) < 1e-10
    assert iters < 60


def test_newton_kleinman_is_monotone_and_quadratic():
    """Kleinman's theorem: P_1 >= P_2 >= ... >= P, and the error roughly squares each step near the end."""
    A, B = SAT_A, SAT_B
    P_star = solve_continuous_are(A, B, SAT_Q, SAT_R)
    K = impl.bass_stabilizing_gain(A, B)
    assert np.all(np.linalg.eigvals(A - B @ K).real < 0)
    errs, prev = [], None
    for _ in range(14):
        P = impl.lyapunov_solve(A - B @ K, SAT_Q + K.T @ SAT_R @ K)
        if prev is not None:
            assert np.all(np.linalg.eigvalsh(prev - P) > -1e-8 * np.abs(P).max())   # decreasing
        assert np.all(np.linalg.eigvalsh(P - P_star) > -1e-8 * np.abs(P_star).max())  # bounded below by P
        errs.append(np.linalg.norm(P - P_star) / np.linalg.norm(P_star))
        K = np.linalg.solve(SAT_R, B.T @ P)
        prev = P
    # quadratic tail: once the error is below 0.1, e_{k+1} <= C e_k^2 (C ~ 1.5 here; 3 allows rounding)
    pairs = [(a, b) for a, b in zip(errs[:-1], errs[1:]) if a < 0.1 and b > 1e-13]
    assert len(pairs) >= 3
    assert all(b <= 3.0 * a ** 2 for a, b in pairs)
    assert errs[-1] < 1e-13


def test_newton_kleinman_rejects_destabilising_start():
    with pytest.raises(ValueError):
        impl.care_newton_kleinman(SAT_A, SAT_B, SAT_Q, SAT_R, K0=np.zeros((1, 2)))   # double integrator, K = 0


def test_hamiltonian_spectrum_is_symmetric_about_the_imaginary_axis():
    A, B, Q, R = random_system(4, 2)
    lam = np.linalg.eigvals(impl.hamiltonian_matrix(A, B, Q, R))
    npt.assert_allclose(np.sort_complex(lam), np.sort_complex(-lam), atol=1e-9)
    # the stable half are the LQR closed-loop poles
    K, P, cl = impl.lqr(A, B, Q, R)
    npt.assert_allclose(np.sort_complex(cl), np.sort_complex(lam[lam.real < 0]), atol=1e-8)


def test_unstabilisable_system_has_no_stabilising_solution():
    A = np.array([[1.0, 0.0], [0.0, -1.0]])
    B = np.array([[0.0], [1.0]])                                  # the unstable mode is uncontrollable
    with pytest.raises(ValueError):
        impl.care_hamiltonian(A, B, np.eye(2), np.eye(1))


# ---------------------------------------------------------------- the satellite regulator

def test_double_integrator_closed_form():
    """For J theta'' = u with Q = diag(q1, q2), R = r: k1 = sqrt(q1/r), k2 = sqrt(q2/r + 2 J k1)."""
    J = 10.0
    q1, q2, r = SAT_Q[0, 0], SAT_Q[1, 1], SAT_R[0, 0]
    for method in ("hamiltonian", "newton"):
        K, P, lam = impl.lqr(SAT_A, SAT_B, SAT_Q, SAT_R, method=method)
        k1 = np.sqrt(q1 / r)
        npt.assert_allclose(K[0], [k1, np.sqrt(q2 / r + 2 * J * k1)], rtol=1e-9)
        npt.assert_allclose(np.sort(lam.real), np.sort(np.linalg.eigvals(SAT_A - SAT_B @ K).real), rtol=1e-12)


def test_bryson_weights():
    Q, R = impl.bryson_weights([0.1, 2.0], [0.5])
    npt.assert_allclose(Q, np.diag([100.0, 0.25]))
    npt.assert_allclose(R, [[4.0]])


def test_optimal_cost_is_x0_P_x0_and_beats_perturbed_gains():
    K, P, _ = impl.lqr(SAT_A, SAT_B, SAT_Q, SAT_R)
    x0 = np.array([np.radians(5.0), np.radians(0.1)])
    J_opt = impl.cost_of_gain(SAT_A, SAT_B, SAT_Q, SAT_R, K, x0)
    npt.assert_allclose(J_opt, x0 @ P @ x0, rtol=1e-9)
    # simulated cost by quadrature equals the Lyapunov value (a third route)
    t, x, u = impl.simulate(SAT_A, SAT_B, K, x0, 300.0, 0.02)
    npt.assert_allclose(impl.quadratic_cost(t, x, u, SAT_Q, SAT_R), J_opt, rtol=1e-4)
    for dK in RNG.standard_normal((20, 2)) * 0.1 * np.abs(K):
        assert impl.cost_of_gain(SAT_A, SAT_B, SAT_Q, SAT_R, K + dK, x0) >= J_opt * (1 - 1e-12)
    assert impl.cost_of_gain(SAT_A, SAT_B, SAT_Q, SAT_R, -K, x0) == np.inf


def test_unsaturated_simulation_matches_the_matrix_exponential():
    K, _, _ = impl.lqr(SAT_A, SAT_B, SAT_Q, SAT_R)
    x0 = np.array([np.radians(5.0), 0.0])
    t, x, u = impl.simulate(SAT_A, SAT_B, K, x0, 60.0, 0.05)
    Acl = SAT_A - SAT_B @ K
    ref = np.array([expm(Acl * ti) @ x0 for ti in t[::50]])
    npt.assert_allclose(x[::50], ref, atol=1e-10)
    npt.assert_allclose(u[:, 0], -(x @ K.T)[:, 0], atol=1e-15)
    npt.assert_allclose(impl.expm(Acl * 7.0), expm(Acl * 7.0), rtol=1e-9, atol=1e-12)


def test_saturated_simulation_respects_the_torque_limit_and_still_converges():
    K, _, _ = impl.lqr(SAT_A, SAT_B, SAT_Q, SAT_R)
    x0 = np.array([np.radians(20.0), 0.0])
    t, x, u = impl.simulate(SAT_A, SAT_B, K, x0, 300.0, 0.05, u_max=0.02)
    assert np.abs(u).max() <= 0.02 + 1e-15
    assert abs(x[-1, 0]) < 1e-3 * x0[0]
    # zero initial state: nothing happens (no spurious input)
    t0, xz, uz = impl.simulate(SAT_A, SAT_B, K, np.zeros(2), 10.0, 0.05, u_max=0.02)
    assert np.abs(xz).max() == 0.0 and np.abs(uz).max() == 0.0


# ---------------------------------------------------------------- guaranteed margins

@pytest.mark.parametrize("seed", [1, 2, 3, 4])
def test_kalman_return_difference_and_guaranteed_margins(seed):
    """Single-input LQR with R > 0: |1 + L(jw)| >= 1 at every w, hence PM >= 60 deg and stability for k in (1/2, inf)."""
    rng = np.random.default_rng(seed)
    n = 3
    A = rng.standard_normal((n, n))
    B = rng.standard_normal((n, 1))
    C = rng.standard_normal((n, n))
    Q = C @ C.T + 0.01 * np.eye(n)
    R = np.array([[rng.uniform(0.1, 10.0)]])
    K, _, _ = impl.lqr(A, B, Q, R)
    w = np.logspace(-3, 3, 3000)
    L = impl.loop_transfer(A, B, K, w)
    assert np.abs(1 + L).min() >= 1 - 1e-9
    pm, wc = impl.phase_margin(w, L)
    assert pm >= 60.0 - 1e-6 or np.isinf(pm)
    lo, hi = impl.gain_margin_interval(A, B, K, np.concatenate([np.linspace(0.51, 1.0, 50), np.logspace(0, 3, 60)]))
    assert lo <= 0.51 + 1e-12 and hi >= 1e3 - 1e-9


def test_loop_transfer_limits():
    K, _, _ = impl.lqr(SAT_A, SAT_B, SAT_Q, SAT_R)
    w = np.array([1e-6, 1e4])
    L = impl.loop_transfer(SAT_A, SAT_B, K, w)
    # low frequency: L ~ k1 / (J s^2) -> huge, phase -180; high frequency: L ~ k2 / (J s), phase -90
    J = 10.0
    npt.assert_allclose(L[0], K[0, 0] / (J * (1j * w[0]) ** 2), rtol=1e-3)
    npt.assert_allclose(L[1], K[0, 1] / (J * 1j * w[1]), rtol=1e-3)


def test_heavier_control_weight_gives_lower_torque_and_slower_response():
    rows = impl.weight_sweep(np.array([0.3, 1.0, 3.0]))
    torques = [r["peak_torque"] for r in rows]
    settles = [r["settle"] for r in rows]
    assert torques[0] > torques[1] > torques[2]
    assert settles[0] < settles[2]
    assert all(r["pm"] >= 60.0 for r in rows)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "lqr_satellite.png").stat().st_size > 1000
