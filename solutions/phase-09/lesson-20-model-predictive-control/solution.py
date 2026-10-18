"""P09.20 Model Predictive Control: reference solution.

Linear MPC from scratch for an altitude-hold autopilot with an elevator travel limit, an elevator rate limit
and a climb-rate limit. At every sample the controller predicts the next N states as a linear function of
the N future inputs, X = Sx x0 + Su U, writes the quadratic cost as a condensed quadratic program in U,
solves it under linear inequality constraints with a dual active-set method (Goldfarb-Idnani), applies the
first input, and repeats (receding horizon).

Independent checks (in the tests):
    * with no active constraint the first MPC move is u0 = -K_mpc x0 with K_mpc equal to the first gain of the
      finite-horizon LQR from the backward Riccati recursion;
    * with the terminal cost P from the discrete algebraic Riccati equation (DARE, solved here by iterating the
      recursion and checked against ``scipy.linalg.solve_discrete_are``) K_mpc equals the infinite-horizon LQR
      gain for every horizon N;
    * the QP solution satisfies the KKT conditions, and the closed loop never violates a constraint.

Conventions:
    * State x = (dh, gamma, theta, q): altitude error h - h_ref [m], flight-path angle [rad], pitch angle [rad],
      pitch rate [rad/s]. Input u = elevator [rad], positive nose-up. Speed constant (autothrottle).
    * Cost J = sum_{k=0}^{N-1} (x_k^T Q x_k + u_k^T R u_k) + x_N^T P x_N; condensed: J = 1/2 U^T H U + x0^T F^T U + const.
    * QP: minimise 1/2 z^T H z + f^T z subject to G z <= w; multipliers lam >= 0 with H z + f + G^T lam = 0.

Run:  python solution.py          (prints the altitude-hold study, writes mpc_altitude.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-09/lesson-20-model-predictive-control
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


# ---------------------------------------------------------------- given: model and discretisation

def expm(M: np.ndarray, n_terms: int = 20) -> np.ndarray:
    """Return e^M by scaling and squaring of a truncated Taylor series (as built in P02.05)."""
    M = np.asarray(M, dtype=float)
    norm = np.linalg.norm(M, 1)
    s = max(0, int(np.ceil(np.log2(norm / 0.5)))) if norm > 0.5 else 0
    Ms = M / 2.0 ** s
    term = np.eye(M.shape[0])
    E = np.eye(M.shape[0])
    for k in range(1, n_terms):
        term = term @ Ms / k
        E = E + term
    for _ in range(s):
        E = E @ E
    return E


def c2d_zoh(A: np.ndarray, B: np.ndarray, T: float) -> tuple[np.ndarray, np.ndarray]:
    """Return the exact zero-order-hold (Ad, Bd) from the augmented matrix exponential (P09.17)."""
    n, m = B.shape
    M = np.zeros((n + m, n + m))
    M[:n, :n] = A
    M[:n, n:] = B
    E = expm(M * T)
    return E[:n, :n], E[:n, n:]


def altitude_model(V: float = 25.0, dt: float = 0.1) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (A, B, Ad, Bd) of a representative small-UAV longitudinal model at constant speed V [m/s].

    dh' = V gamma;  gamma' = L_alpha (theta - gamma);  theta' = q;
    q' = M_alpha (theta - gamma) + M_q q + M_delta delta_e,
    with L_alpha = 1.0 1/s, M_alpha = -9 1/s^2, M_q = -2.5 1/s, M_delta = 12 1/s^2 (alpha = theta - gamma).
    Discretised by ZOH at dt (default 0.1 s, a 10 Hz outer loop). Representative numbers, not one aircraft.
    """
    La, Ma, Mq, Md = 1.0, -9.0, -2.5, 12.0
    A = np.array([[0.0, V, 0.0, 0.0],
                  [0.0, -La, La, 0.0],
                  [0.0, 0.0, 0.0, 1.0],
                  [0.0, -Ma, Ma, Mq]])
    B = np.array([[0.0], [0.0], [0.0], [Md]])
    Ad, Bd = c2d_zoh(A, B, dt)
    return A, B, Ad, Bd


def altitude_limits(V: float = 25.0, dt: float = 0.1) -> dict[str, float]:
    """Return the constraint set: elevator +-10 deg, elevator rate +-30 deg/s (per-sample change +-3 deg at
    10 Hz), climb rate |V gamma| <= 2.5 m/s."""
    return {"u_max": np.radians(10.0), "du_max": np.radians(30.0) * dt, "climb_max": 2.5, "V": V}


def altitude_weights() -> tuple[np.ndarray, np.ndarray]:
    """Return (Q, R) from Bryson-like scales: 2 m, 0.1 rad, 0.2 rad, 0.5 rad/s; elevator 0.17 rad."""
    return np.diag([1 / 2.0 ** 2, 1 / 0.1 ** 2, 1 / 0.2 ** 2, 1 / 0.5 ** 2]), np.array([[1 / 0.17 ** 2]])


# ---------------------------------------------------------------- Riccati

def finite_horizon_lqr(Ad: np.ndarray, Bd: np.ndarray, Q: np.ndarray, R: np.ndarray, P_N: np.ndarray,
                       N: int) -> tuple[list[np.ndarray], np.ndarray]:
    """Return ([K_0, ..., K_{N-1}], P_0) of the finite-horizon discrete LQR by the backward Riccati recursion
    P_N given;  K_k = (R + B^T P_{k+1} B)^{-1} B^T P_{k+1} A;  P_k = Q + A^T P_{k+1} (A - B K_k)."""
    P = np.asarray(P_N, dtype=float)
    Ks = [None] * N
    for k in range(N - 1, -1, -1):
        K = np.linalg.solve(R + Bd.T @ P @ Bd, Bd.T @ P @ Ad)
        P = Q + Ad.T @ P @ (Ad - Bd @ K)
        P = 0.5 * (P + P.T)
        Ks[k] = K
    return Ks, P


def dare(Ad: np.ndarray, Bd: np.ndarray, Q: np.ndarray, R: np.ndarray, tol: float = 1e-12,
         max_iter: int = 100000) -> np.ndarray:
    """Return the stabilising solution of the DARE P = Q + A^T P A - A^T P B (R + B^T P B)^{-1} B^T P A by
    iterating the Riccati recursion from P = Q until ||P_{k+1} - P_k||_F <= tol ||P_k||_F.
    Raise ``RuntimeError`` if it has not converged in max_iter steps."""
    P = np.asarray(Q, dtype=float).copy()
    for _ in range(max_iter):
        K = np.linalg.solve(R + Bd.T @ P @ Bd, Bd.T @ P @ Ad)
        Pn = Q + Ad.T @ P @ (Ad - Bd @ K)
        Pn = 0.5 * (Pn + Pn.T)
        if np.linalg.norm(Pn - P) <= tol * np.linalg.norm(Pn):
            return Pn
        P = Pn
    raise RuntimeError("Riccati iteration did not converge")


def lqr_gain(Ad: np.ndarray, Bd: np.ndarray, Q: np.ndarray, R: np.ndarray, P: np.ndarray) -> np.ndarray:
    """Return K = (R + B^T P B)^{-1} B^T P A."""
    return np.linalg.solve(R + Bd.T @ P @ Bd, Bd.T @ P @ Ad)


# ---------------------------------------------------------------- condensed MPC

def prediction_matrices(Ad: np.ndarray, Bd: np.ndarray, N: int) -> tuple[np.ndarray, np.ndarray]:
    """Return (Sx, Su) with X = (x_1, ..., x_N) stacked = Sx x0 + Su U, U = (u_0, ..., u_{N-1}) stacked.

    Block row k (k = 1..N): Sx_k = A^k, Su_{k,j} = A^{k-1-j} B for j < k, zero otherwise.
    """
    n, m = Bd.shape
    Sx = np.zeros((N * n, n))
    Su = np.zeros((N * n, N * m))
    Ak = np.eye(n)
    powers = [np.eye(n)]
    for k in range(1, N + 1):
        Ak = Ad @ Ak
        powers.append(Ak)
        Sx[(k - 1) * n:k * n, :] = Ak
    for k in range(1, N + 1):
        for j in range(k):
            Su[(k - 1) * n:k * n, j * m:(j + 1) * m] = powers[k - 1 - j] @ Bd
    return Sx, Su


def condensed_cost(Ad: np.ndarray, Bd: np.ndarray, Q: np.ndarray, R: np.ndarray, P: np.ndarray,
                   N: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (H, F, Y) with J(U; x0) = 1/2 U^T H U + x0^T F^T U + x0^T Y x0.

    With Qbar = blkdiag(Q, ..., Q, P) (N blocks, the last P) and Rbar = blkdiag(R, ..., R):
    H = 2 (Su^T Qbar Su + Rbar),  F = 2 Su^T Qbar Sx,  Y = Q + Sx^T Qbar Sx.
    """
    n, m = Bd.shape
    Sx, Su = prediction_matrices(Ad, Bd, N)
    Qbar = np.zeros((N * n, N * n))
    for k in range(N - 1):
        Qbar[k * n:(k + 1) * n, k * n:(k + 1) * n] = Q
    Qbar[(N - 1) * n:, (N - 1) * n:] = P
    Rbar = np.kron(np.eye(N), R)
    H = 2.0 * (Su.T @ Qbar @ Su + Rbar)
    F = 2.0 * Su.T @ Qbar @ Sx
    Y = Q + Sx.T @ Qbar @ Sx
    return 0.5 * (H + H.T), F, Y


def unconstrained_mpc_gain(Ad: np.ndarray, Bd: np.ndarray, Q: np.ndarray, R: np.ndarray, P: np.ndarray,
                           N: int) -> np.ndarray:
    """Return K_mpc (m x n): the first m rows of H^{-1} F, so the unconstrained MPC applies u0 = -K_mpc x0."""
    H, F, _ = condensed_cost(Ad, Bd, Q, R, P, N)
    m = Bd.shape[1]
    return np.linalg.solve(H, F)[:m, :]


def mpc_constraints(Ad: np.ndarray, Bd: np.ndarray, N: int, x0: np.ndarray, u_prev: float, u_max: float,
                    du_max: float, c_row: np.ndarray | None = None, c_max: float | None = None,
                    pred: tuple[np.ndarray, np.ndarray] | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Return (G, w) with G U <= w encoding, for a single input:
        |u_k| <= u_max (k = 0..N-1);  |u_k - u_{k-1}| <= du_max with u_{-1} = u_prev;
        and, if c_row is given, |c_row x_k| <= c_max for k = 1..N (x_k = Sx_k x0 + Su_k U).
    Rows in that order: +u, -u, +du, -du, +c, -c. ``pred`` = (Sx, Su) may be passed to avoid recomputing them.
    """
    I = np.eye(N)
    D = I - np.eye(N, k=-1)                          # (D U)_k = u_k - u_{k-1}, first row u_0
    e = np.zeros(N)
    e[0] = u_prev
    Gs = [I, -I, D, -D]
    ws = [np.full(N, u_max), np.full(N, u_max), du_max + e, du_max - e]
    if c_row is not None:
        Sx, Su = prediction_matrices(Ad, Bd, N) if pred is None else pred
        Cbar = np.kron(np.eye(N), np.atleast_2d(c_row))   # picks c_row x_k from the stacked X
        CS = Cbar @ Su
        cx = Cbar @ Sx @ np.asarray(x0, dtype=float)
        Gs += [CS, -CS]
        ws += [c_max - cx, c_max + cx]
    return np.vstack(Gs), np.concatenate(ws)


# ---------------------------------------------------------------- QP: dual active set (Goldfarb-Idnani)

def solve_qp(H: np.ndarray, f: np.ndarray, G: np.ndarray, w: np.ndarray, tol: float = 1e-10,
             max_iter: int = 2000) -> tuple[np.ndarray, np.ndarray, int]:
    """Return (z, lam, iterations) minimising 1/2 z^T H z + f^T z subject to G z <= w (H positive definite).

    Dual active-set method of Goldfarb and Idnani: start at the unconstrained minimum z = -H^{-1} f (dual
    feasible, lam = 0); repeatedly pick the most violated constraint p and move along
        z_dir = -H^{-1} (g_p + N^T r),   r = -(N H^{-1} N^T)^{-1} N H^{-1} g_p
    (N = rows of G in the active set), which keeps the active constraints at equality while lam_p grows by t
    and the active multipliers change by t r. The step t is the smaller of the full step that makes
    constraint p active and the partial step at which an active multiplier reaches zero (that constraint is
    then dropped). Every iterate is dual feasible (lam >= 0) and the objective rises monotonically; the method
    ends when no constraint is violated by more than ``tol`` (scaled). Raise ``ValueError`` if infeasible.
    """
    H = np.asarray(H, dtype=float)
    f = np.asarray(f, dtype=float)
    G = np.atleast_2d(np.asarray(G, dtype=float))
    w = np.asarray(w, dtype=float)
    Hinv = np.linalg.inv(H)
    z = -Hinv @ f
    m_c = G.shape[0]
    lam = np.zeros(m_c)
    active: list[int] = []
    scale = 1.0 + np.abs(w).max() if m_c else 1.0
    it = 0
    while True:
        viol = G @ z - w
        p = int(np.argmax(viol)) if m_c else -1
        if m_c == 0 or viol[p] <= tol * scale:
            return z, lam, it
        while True:
            it += 1
            if it > max_iter:
                raise RuntimeError("QP iteration limit")
            gp = G[p]
            if active:
                Nm = G[active]
                NHi = Nm @ Hinv
                r = -np.linalg.solve(NHi @ Nm.T, NHi @ gp)
                zd = -Hinv @ (gp + Nm.T @ r)
            else:
                r = np.zeros(0)
                zd = -Hinv @ gp
            # partial (dual) step: an active multiplier hits zero
            t1, k_drop = np.inf, -1
            for j, rj in enumerate(r):
                if rj < -1e-14:
                    tj = lam[active[j]] / (-rj)
                    if tj < t1:
                        t1, k_drop = tj, j
            # full (primal) step: constraint p becomes active
            gz = gp @ zd
            t2 = (gp @ z - w[p]) / (-gz) if gz < -1e-14 * (1.0 + np.linalg.norm(gp) ** 2 * np.abs(Hinv).max()) else np.inf
            t = min(t1, t2)
            if not np.isfinite(t):
                raise ValueError("QP is infeasible")
            if np.isfinite(t2):
                z = z + t * zd
            for j, idx in enumerate(active):
                lam[idx] += t * r[j]
            lam[p] += t
            if t2 <= t1:
                active.append(p)
                break
            lam[active[k_drop]] = 0.0
            active.pop(k_drop)


def kkt_residuals(H: np.ndarray, f: np.ndarray, G: np.ndarray, w: np.ndarray, z: np.ndarray,
                  lam: np.ndarray) -> dict[str, float]:
    """Return the four KKT residuals: stationarity ||H z + f + G^T lam||inf, primal max(G z - w, 0),
    dual max(-lam, 0) and complementarity max |lam_i (G z - w)_i|."""
    s = G @ z - w
    return {"stationarity": float(np.abs(H @ z + f + G.T @ lam).max()),
            "primal": float(max(0.0, s.max())), "dual": float(max(0.0, (-lam).max())),
            "complementarity": float(np.abs(lam * s).max())}


# ---------------------------------------------------------------- receding horizon

def mpc_step(Ad: np.ndarray, Bd: np.ndarray, Q: np.ndarray, R: np.ndarray, P: np.ndarray, N: int,
             x0: np.ndarray, u_prev: float, limits: dict[str, float], climb_constraint: bool = True,
             cache: dict | None = None) -> tuple[float, np.ndarray, int]:
    """Return (u0, U, QP iterations): solve the constrained MPC problem at state x0 and return its first move.

    ``cache`` (a dict, filled on first use) keeps H, F, Sx, Su between calls with the same model and horizon.
    """
    if cache is None:
        cache = {}
    if "H" not in cache:
        cache["H"], cache["F"], _ = condensed_cost(Ad, Bd, Q, R, P, N)
        cache["pred"] = prediction_matrices(Ad, Bd, N)
    H, F = cache["H"], cache["F"]
    c_row = np.array([0.0, limits["V"], 0.0, 0.0]) if climb_constraint else None
    G, w = mpc_constraints(Ad, Bd, N, x0, u_prev, limits["u_max"], limits["du_max"], c_row,
                           limits["climb_max"] if climb_constraint else None, pred=cache["pred"])
    U, lam, it = solve_qp(H, F @ x0, G, w)
    return float(U[0]), U, it


def simulate_mpc(x0: np.ndarray, n_steps: int, N: int = 30, terminal: str = "dare",
                 limits: dict[str, float] | None = None, climb_constraint: bool = True) -> dict[str, np.ndarray]:
    """Run the receding-horizon loop on the altitude model (the plant equals the prediction model).

    terminal "dare" uses P from ``dare``; "Q" uses P = Q. Returns arrays t, x (n_steps + 1, 4), u (n_steps,),
    iters (n_steps,). Raise ``ValueError`` naming the step if a QP is infeasible (no input sequence can meet
    every constraint over the horizon from the current state).
    """
    _, _, Ad, Bd = altitude_model()
    Q, R = altitude_weights()
    lim = altitude_limits() if limits is None else limits
    P = dare(Ad, Bd, Q, R) if terminal == "dare" else Q
    x = np.asarray(x0, dtype=float).copy()
    xs, us, its = [x.copy()], [], []
    u_prev = 0.0
    cache: dict = {}
    for k in range(n_steps):
        try:
            u, _, it = mpc_step(Ad, Bd, Q, R, P, N, x, u_prev, lim, climb_constraint, cache)
        except ValueError as exc:
            raise ValueError(f"MPC QP infeasible at step {k} (t = {0.1 * k:.1f} s)") from exc
        u = float(np.clip(u, -lim["u_max"], lim["u_max"]))     # guards rounding only (|change| < 1e-12)
        x = Ad @ x + Bd[:, 0] * u
        xs.append(x.copy())
        us.append(u)
        its.append(it)
        u_prev = u
    return {"t": np.arange(n_steps + 1) * 0.1, "x": np.array(xs), "u": np.array(us), "iters": np.array(its)}


def simulate_lqr(x0: np.ndarray, n_steps: int, clip: dict[str, float] | None = None) -> dict[str, np.ndarray]:
    """Run u = -K_inf x on the altitude model; if ``clip`` is given, saturate u and its per-sample change to the
    limits (a naive way to respect the actuator, which ignores the climb-rate limit)."""
    _, _, Ad, Bd = altitude_model()
    Q, R = altitude_weights()
    K = lqr_gain(Ad, Bd, Q, R, dare(Ad, Bd, Q, R))
    x = np.asarray(x0, dtype=float).copy()
    xs, us = [x.copy()], []
    u_prev = 0.0
    for _ in range(n_steps):
        u = float((-K @ x)[0])
        if clip is not None:
            u = float(np.clip(u, u_prev - clip["du_max"], u_prev + clip["du_max"]))
            u = float(np.clip(u, -clip["u_max"], clip["u_max"]))
        x = Ad @ x + Bd[:, 0] * u
        xs.append(x.copy())
        us.append(u)
        u_prev = u
    return {"t": np.arange(n_steps + 1) * 0.1, "x": np.array(xs), "u": np.array(us)}


def settling_time(t: np.ndarray, dh: np.ndarray, band: float = 0.4) -> float:
    """Return the first time after which |dh| stays within ``band`` metres (inf if it is outside at the end)."""
    out = np.where(np.abs(dh) > band)[0]
    if len(out) == 0:
        return 0.0
    return float(t[out[-1] + 1]) if out[-1] + 1 < len(t) else float("inf")


def main(out_dir: str | Path = ".") -> None:
    """Print the MPC checks and the altitude-hold study, and save mpc_altitude.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    A, B, Ad, Bd = altitude_model()
    Q, R = altitude_weights()
    lim = altitude_limits()
    print("altitude model eigenvalues (continuous):", np.round(np.linalg.eigvals(A), 3))
    P = dare(Ad, Bd, Q, R)
    K = lqr_gain(Ad, Bd, Q, R, P)
    print(f"infinite-horizon LQR K = {np.round(K[0], 4).tolist()}; closed-loop |z| = {np.round(np.abs(np.linalg.eigvals(Ad - Bd @ K)), 3)}")
    for N in (5, 20):
        Km = unconstrained_mpc_gain(Ad, Bd, Q, R, Q, N)
        Ks, _ = finite_horizon_lqr(Ad, Bd, Q, R, Q, N)
        Kp = unconstrained_mpc_gain(Ad, Bd, Q, R, P, N)
        print(f"N = {N:2d}: |K_mpc(P=Q) - K_0 Riccati| = {np.abs(Km - Ks[0]).max():.1e}; |K_mpc(P=DARE) - K_inf| = {np.abs(Kp - K).max():.1e};"
              f" K_mpc(P=Q) = {np.round(Km[0], 4).tolist()}")

    x0 = np.array([-20.0, 0.0, 0.0, 0.0])          # 20 m below the commanded altitude
    n = 250
    lq = simulate_lqr(x0, n)
    lq_dh = np.diff(np.concatenate([[0.0], lq["u"]]))
    print(f"\n20 m climb, unconstrained LQR: peak elevator {np.degrees(np.abs(lq['u']).max()):.1f} deg, peak rate {np.degrees(np.abs(lq_dh).max()) / 0.1:.0f} deg/s,"
          f" peak climb {np.abs(lim['V'] * lq['x'][:, 1]).max():.2f} m/s, settling (0.4 m) {settling_time(lq['t'], lq['x'][:, 0]):.1f} s")
    lc = simulate_lqr(x0, n, clip=lim)
    print(f"   LQR clipped to the actuator: peak climb {np.abs(lim['V'] * lc['x'][:, 1]).max():.2f} m/s (limit 2.5),"
          f" overshoot {max(0.0, lc['x'][:, 0].max()):.2f} m, settling {settling_time(lc['t'], lc['x'][:, 0]):.1f} s")
    for N in (10, 20, 30, 50):
        r = simulate_mpc(x0, n, N=N)
        du = np.diff(np.concatenate([[0.0], r["u"]]))
        print(f"   MPC N = {N:2d} (P = DARE): peak elevator {np.degrees(np.abs(r['u']).max()):.2f} deg, peak rate {np.degrees(np.abs(du).max()) / 0.1:.2f} deg/s,"
              f" peak climb {np.abs(lim['V'] * r['x'][:, 1]).max():.4f} m/s, overshoot {max(0.0, r['x'][:, 0].max()):.3f} m,"
              f" settling {settling_time(r['t'], r['x'][:, 0]):.1f} s, QP iterations max {r['iters'].max()}")
    r30 = simulate_mpc(x0, n, N=30)
    print("   horizon study, terminal cost P = Q (no tail) and P = DARE:")
    for N in (3, 5, 10, 20):
        for term in ("Q", "dare"):
            try:
                r = simulate_mpc(x0, n, N=N, terminal=term)
                print(f"      N = {N:2d}, P = {term:4s}: overshoot {max(0.0, r['x'][:, 0].max()):.3f} m, settling {settling_time(r['t'], r['x'][:, 0]):.1f} s,"
                      f" error at 25 s {r['x'][-1, 0]:+.3f} m")
            except ValueError as exc:
                print(f"      N = {N:2d}, P = {term:4s}: {exc}")
    for N in (3, 5, 10, 20):
        Km = unconstrained_mpc_gain(Ad, Bd, Q, R, Q, N)
        print(f"      unconstrained N = {N:2d}, P = Q: closed-loop spectral radius {np.max(np.abs(np.linalg.eigvals(Ad - Bd @ Km))):.4f}")

    # KKT certificate for the first QP
    H, F, _ = condensed_cost(Ad, Bd, Q, R, P, 30)
    G, w = mpc_constraints(Ad, Bd, 30, x0, 0.0, lim["u_max"], lim["du_max"], np.array([0.0, lim["V"], 0.0, 0.0]), lim["climb_max"])
    U, lam, it = solve_qp(H, F @ x0, G, w)
    kk = kkt_residuals(H, F @ x0, G, w, U, lam)
    print(f"first QP: {G.shape[0]} constraints, {int(np.sum(lam > 1e-9))} active, {it} iterations;"
          f" KKT residuals {', '.join(f'{k} {v:.1e}' for k, v in kk.items())}")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.0))
    for res, lab, st in ((lq, "LQR, unconstrained", "--"), (lc, "LQR, clipped", ":"), (r30, "MPC, N = 30", "-")):
        axes[0].plot(res["t"], 20.0 + res["x"][:, 0], st, label=lab)
        axes[1].plot(res["t"], lim["V"] * res["x"][:, 1], st, label=lab)
        axes[2].step(res["t"][:-1], np.degrees(res["u"]), st, where="post", label=lab)
    axes[0].set_ylabel("altitude gained [m]")
    axes[1].axhline(lim["climb_max"], color="k", lw=0.6); axes[1].axhline(-lim["climb_max"], color="k", lw=0.6)
    axes[1].set_ylabel("climb rate [m/s]"); axes[1].set_ylim(-4, 11)
    axes[2].axhline(10.0, color="k", lw=0.6); axes[2].axhline(-10.0, color="k", lw=0.6)
    axes[2].set_ylabel("elevator [deg]"); axes[2].set_ylim(-25, 25)
    for ax in axes:
        ax.set_xlabel("t [s]"); ax.grid(lw=0.4); ax.set_xlim(0, 15)
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "mpc_altitude.png", dpi=150)
    plt.close(fig)
    print(f"saved {out_dir / 'mpc_altitude.png'}")


if __name__ == "__main__":
    main()
