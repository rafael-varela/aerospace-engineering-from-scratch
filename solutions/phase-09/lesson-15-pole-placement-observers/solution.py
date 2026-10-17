"""P09.15 Pole Placement & Observers: reference solution.

Full-state feedback u = N r - K x with K from Ackermann's formula, a Luenberger observer whose gain L is
Ackermann's formula applied to the dual system (A^T, C^T), the separation principle (the closed loop with
the observer has exactly the regulator poles and the observer poles), a reference gain for zero
steady-state error, and the observer error decaying at the rate of the slowest placed observer pole.

The aerospace case is single-axis attitude control of a small satellite with a reaction wheel,
J theta'' = tau, after its gyro has failed: only the star-tracker angle is measured, and the observer
supplies the rate. A second, four-state example places the poles of the pitch-attitude model of P09.11
(short period, integrator, elevator servo). Star-tracker noise in the experiment is SYNTHETIC and seeded.

Conventions:
    * Single input for Ackermann; B is (n, 1), C is (1, n) for the observer, K is (1, n), L is (n, 1).
    * Satellite state x = (theta, omega) [rad, rad/s]; input wheel torque tau [N m].

Run:  python solution.py          (prints the design, writes pole_placement.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-09/lesson-15-pole-placement-observers
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


# ---------------------------------------------------------------- given

def expm(M: np.ndarray, n_terms: int = 20) -> np.ndarray:
    """Return e^M by scaling and squaring a truncated Taylor series (the method of P02.05)."""
    M = np.asarray(M, dtype=float)
    norm = np.linalg.norm(M, 1)
    s = max(0, int(np.ceil(np.log2(norm / 0.5)))) if norm > 0.5 else 0
    X = M / 2 ** s
    term = np.eye(len(M))
    E = np.eye(len(M))
    for k in range(1, n_terms):
        term = term @ X / k
        E = E + term
    for _ in range(s):
        E = E @ E
    return E


def discretize_zoh(A: np.ndarray, B: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (A_d, B_d) of the exact zero-order-hold discretisation from e^{[[A, B], [0, 0]] dt}."""
    n, m = B.shape
    M = np.zeros((n + m, n + m))
    M[:n, :n] = A
    M[:n, n:] = B
    E = expm(M * dt)
    return E[:n, :n], E[:n, n:]


def satellite_axis(J: float = 3.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (A, B, C) of one satellite axis, J theta'' = tau, measured by the star tracker only (gyro failed)."""
    return np.array([[0.0, 1.0], [0.0, 0.0]]), np.array([[0.0], [1.0 / J]]), np.array([[1.0, 0.0]])


def pitch_model() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (A, B, C) of the representative pitch model of P09.11: x = (alpha, q, theta, delta), input delta_cmd, output theta."""
    A = np.array([[-1.2, 1.0, 0.0, 0.0],
                  [-6.0, -2.5, 0.0, 12.0],
                  [0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 0.0, -20.0]])
    return A, np.array([[0.0], [0.0], [0.0], [20.0]]), np.array([[0.0, 0.0, 1.0, 0.0]])


def second_order_poles(zeta: float, wn: float) -> np.ndarray:
    """Return the pair -zeta wn +- j wn sqrt(1 - zeta^2) (a repeated real pole -wn when zeta = 1)."""
    if zeta >= 1.0:
        r = np.sqrt(zeta ** 2 - 1.0)
        return np.array([-zeta * wn + wn * r, -zeta * wn - wn * r], dtype=complex)
    wd = wn * np.sqrt(1.0 - zeta ** 2)
    return np.array([-zeta * wn + 1j * wd, -zeta * wn - 1j * wd])


# ---------------------------------------------------------------- design

def char_poly(poles: np.ndarray) -> np.ndarray:
    """Return the real monic polynomial with the given roots (complex poles must come in conjugate pairs)."""
    c = np.poly(np.asarray(poles))
    if np.max(np.abs(np.imag(c))) > 1e-9 * max(1.0, np.max(np.abs(c))):
        raise ValueError("complex poles must come in conjugate pairs")
    return np.real(c)


def ackermann(A: np.ndarray, B: np.ndarray, poles: np.ndarray) -> np.ndarray:
    """Return K (1 x n) with eig(A - B K) = poles, by Ackermann's formula K = e_n^T W_c^{-1} phi_d(A).

    W_c = [B, AB, ..., A^{n-1} B] and phi_d(A) = A^n + a_1 A^{n-1} + ... + a_n I is the desired characteristic
    polynomial evaluated at A. Raise ``ValueError`` if (A, B) is not controllable (rank W_c < n): no gain can
    move an uncontrollable pole.
    """
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float).reshape(-1, 1)
    n = A.shape[0]
    Wc = np.hstack([np.linalg.matrix_power(A, k) @ B for k in range(n)])
    sv = np.linalg.svd(Wc, compute_uv=False)
    if sv[-1] <= 1e-12 * sv[0]:
        raise ValueError("(A, B) is not controllable: rank W_c < n")
    a = char_poly(poles)
    phi = np.zeros((n, n))
    for c in a:                                     # Horner: phi = (...((A + a1 I) A + a2 I) A ...) + a_n I
        phi = phi @ A + c * np.eye(n)
    en = np.zeros(n)
    en[-1] = 1.0
    row = np.linalg.solve(Wc.T, en)                 # e_n^T W_c^{-1}
    return (row @ phi)[None, :]


def observer_gain(A: np.ndarray, C: np.ndarray, poles: np.ndarray) -> np.ndarray:
    """Return L (n x 1) with eig(A - L C) = poles, by duality: L = ackermann(A^T, C^T, poles)^T.

    Raise ``ValueError`` (from ``ackermann``) if (A, C) is not observable.
    """
    try:
        return ackermann(A.T, C.T, poles).T
    except ValueError as exc:
        raise ValueError("(A, C) is not observable: rank of the observability matrix < n") from exc


def reference_gain(A: np.ndarray, B: np.ndarray, C: np.ndarray, K: np.ndarray) -> float:
    """Return N so that u = N r - K x gives unit DC gain from r to y: N = -1 / (C (A - B K)^{-1} B)."""
    g = (C @ np.linalg.solve(A - B @ K, B))[0, 0]
    return float(-1.0 / g)


def observer_closed_loop(A: np.ndarray, B: np.ndarray, C: np.ndarray, K: np.ndarray, L: np.ndarray,
                         N: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (A_cl, B_cl) of the plant plus observer with u = N r - K x_hat, state (x, x_hat), inputs (r, v).

    x' = A x + B (N r - K x_hat);  x_hat' = A x_hat + B (N r - K x_hat) + L (C x + v - C x_hat), v = sensor noise.
    In the coordinates (x, e = x - x_hat) the matrix is block triangular [[A - BK, BK], [0, A - LC]]: the
    separation principle.
    """
    n = A.shape[0]
    Acl = np.block([[A, -B @ K], [L @ C, A - B @ K - L @ C]])
    Bcl = np.zeros((2 * n, 2))
    Bcl[:n, 0] = (B * N)[:, 0]
    Bcl[n:, 0] = (B * N)[:, 0]
    Bcl[n:, 1] = L[:, 0]
    return Acl, Bcl


def simulate_observer_loop(A: np.ndarray, B: np.ndarray, C: np.ndarray, K: np.ndarray, L: np.ndarray, N: float,
                           x0: np.ndarray, xhat0: np.ndarray, r: np.ndarray, v: np.ndarray, dt: float) -> dict[str, np.ndarray]:
    """Simulate the observer-based loop with r and v held over each dt (exact ZOH); return x, x_hat, u, e = x - x_hat.

    r and v are arrays of length N_steps (the reference and the measurement noise samples); outputs have
    N_steps rows, the state at the start of each hold.
    """
    n = A.shape[0]
    Acl, Bcl = observer_closed_loop(A, B, C, K, L, N)
    Ad, Bd = discretize_zoh(Acl, Bcl, dt)
    z = np.concatenate([np.asarray(x0, dtype=float), np.asarray(xhat0, dtype=float)])
    Z = np.empty((len(r), 2 * n))
    w = np.column_stack([r, v])
    for k in range(len(r)):
        Z[k] = z
        z = Ad @ z + Bd @ w[k]
    x, xh = Z[:, :n], Z[:, n:]
    u = N * np.asarray(r) - xh @ K[0]
    return {"x": x, "x_hat": xh, "u": u, "e": x - xh}


def fit_decay_rate(t: np.ndarray, e: np.ndarray, t_start: float, t_end: float = np.inf) -> float:
    """Return the decay rate sigma [1/s] of a least-squares fit log|e| = c - sigma t over t_start <= t <= t_end.

    e may be (N,) or (N, n) (the Euclidean norm of each row is used). Choose the window after the faster modes
    have died out and before |e| reaches the rounding floor of the states it is the difference of; samples below
    1e-12 x the first sample in the window are dropped as a safeguard.
    """
    e = np.asarray(e, dtype=float)
    mag = np.linalg.norm(e, axis=1) if e.ndim == 2 else np.abs(e)
    sel = (t >= t_start) & (t <= t_end)
    tt, mm = t[sel], mag[sel]
    keep = mm > 1e-12 * mm[0]
    slope = np.polyfit(tt[keep], np.log(mm[keep]), 1)[0]
    return float(-slope)


# ---------------------------------------------------------------- the example

def main(out_dir: str | Path = ".") -> None:
    """Design the satellite regulator and observer, check separation and decay, run the noise study, save a figure."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    J = 3.0
    A, B, C = satellite_axis(J)
    print(f"one satellite axis, J = {J} kg m^2 (representative), star tracker only (gyro failed); wheel limit 10 mN m")
    reg = second_order_poles(0.707, 0.1)
    obs = np.array([-0.3, -0.6])
    K = ackermann(A, B, reg)
    L = observer_gain(A, C, obs)
    N = reference_gain(A, B, C, K)
    print(f"regulator poles {np.round(reg, 4)} -> K = {np.array2string(K[0], precision=5)} (= J (wn^2, 2 zeta wn))")
    print(f"observer poles {obs} -> L = {np.array2string(L[:, 0], precision=4)}; reference gain N = {N:.5f} (= K[0])")
    print("eig(A - BK) =", np.round(np.linalg.eigvals(A - B @ K), 6), " eig(A - LC) =", np.round(np.linalg.eigvals(A - L @ C), 6))
    Acl, _ = observer_closed_loop(A, B, C, K, L, N)
    print("eig of the 4x4 loop with the observer:", np.round(np.sort_complex(np.linalg.eigvals(Acl)), 6))

    dt = 0.1
    t = np.arange(0.0, 120.0, dt)
    ref = np.radians(10.0) * np.ones_like(t)
    x0 = np.array([0.0, np.radians(0.05)])            # unknown tip-off rate 0.05 deg/s
    run = simulate_observer_loop(A, B, C, K, L, N, x0, np.zeros(2), ref, np.zeros_like(t), dt)
    sf = simulate_observer_loop(A, B, C, K, L, N, x0, x0, ref, np.zeros_like(t), dt)
    theta = run["x"][:, 0]
    os_pct = 100 * (theta.max() - ref[0]) / ref[0]
    print(f"10 deg slew from a 0.05 deg/s tip-off: overshoot {os_pct:.2f} %, peak torque {1e3 * np.abs(run['u']).max():.3f} mN m,"
          f" final error {np.degrees(ref[0] - theta[-1]):.2e} deg")
    print(f"same slew with a perfect initial estimate: peak torque {1e3 * np.abs(sf['u']).max():.3f} mN m;"
          f" max |theta_observer - theta_state_feedback| = {np.degrees(np.abs(theta - sf['x'][:, 0]).max()):.4f} deg")
    rate = fit_decay_rate(t, run["e"], 20.0, 60.0)
    print(f"observer error: fitted decay rate {rate:.5f} 1/s over 20-60 s (slowest placed observer pole {-obs.max():.2f});"
          f" |e| falls from {np.linalg.norm(run['e'][0]):.2e} to {np.linalg.norm(run['e'][600]):.2e} at 60 s")

    print("\nregulator speed vs wheel torque (10 deg step, perfect estimate):")
    for wn in (0.05, 0.1, 0.2, 0.4):
        Ki = ackermann(A, B, second_order_poles(0.707, wn))
        Ni = reference_gain(A, B, C, Ki)
        ri = simulate_observer_loop(A, B, C, Ki, L, Ni, np.zeros(2), np.zeros(2), ref, np.zeros_like(t), dt)
        print(f"  wn = {wn:.2f} rad/s: peak torque {1e3 * np.abs(ri['u']).max():7.3f} mN m (J wn^2 x 10 deg = {1e3 * J * wn ** 2 * ref[0]:.3f})")

    sigma = np.radians(10.0 / 3600.0)                 # 10 arcsec star-tracker noise, SYNTHETIC
    rng = np.random.default_rng(1509)
    v = sigma * rng.standard_normal(len(t))
    print(f"\nSYNTHETIC star-tracker noise, sigma = 10 arcsec at {1 / dt:.0f} Hz, seed 1509 (holding 0 deg):")
    noise_rows = []
    for wo in (0.15, 0.3, 0.6, 1.2):
        Li = observer_gain(A, C, np.array([-wo, -2 * wo]))
        ri = simulate_observer_loop(A, B, C, K, Li, N, np.zeros(2), np.zeros(2), np.zeros_like(t), v, dt)
        late = t >= 40.0
        rate_err = np.std(ri["e"][late, 1])
        tq = np.std(ri["u"][late])
        noise_rows.append((wo, rate_err, tq))
        print(f"  observer poles (-{wo}, -{2 * wo}): rms rate-estimate error {np.degrees(rate_err) * 3600:.3f} arcsec/s, rms torque {1e6 * tq:.3f} uN m")

    Ap, Bp, Cp = pitch_model()
    pitch_poles = np.array([-4 + 4j, -4 - 4j, -6.0, -20.0])
    Kp = ackermann(Ap, Bp, pitch_poles)
    print(f"\npitch model (4 states): K = {np.array2string(Kp[0], precision=4)}; eig(A - BK) =",
          np.round(np.sort_complex(np.linalg.eigvals(Ap - Bp @ Kp)), 6))
    try:
        observer_gain(A, np.array([[0.0, 1.0]]), obs)
    except ValueError as exc:
        print("observer from the gyro alone (star tracker failed instead):", exc)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    axes[0].plot(t, np.degrees(theta), label="with observer, x_hat(0) = 0")
    axes[0].plot(t, np.degrees(sf["x"][:, 0]), "--", label="perfect estimate")
    axes[0].set_xlabel("t [s]"); axes[0].set_ylabel("theta [deg]"); axes[0].legend(fontsize=8); axes[0].grid(lw=0.4)
    axes[0].set_title("10 deg slew (simulation)", fontsize=9)
    en = np.linalg.norm(run["e"], axis=1)
    axes[1].semilogy(t, en, label="|x - x_hat|")
    axes[1].semilogy(t, en[200] * np.exp(-rate * (t - t[200])), "k:", label=f"fit 20-60 s, rate {rate:.3f} 1/s")
    axes[1].set_xlabel("t [s]"); axes[1].legend(fontsize=8); axes[1].grid(True, which="both", lw=0.3)
    axes[1].set_title("observer error", fontsize=9)
    wo_, re_, tq_ = zip(*noise_rows)
    axes[2].loglog(wo_, np.degrees(re_) * 3600, "o-", label="rate error [arcsec/s]")
    axes[2].loglog(wo_, 1e6 * np.array(tq_), "s-", label="torque [uN m]")
    axes[2].set_xlabel("slowest observer pole [1/s]"); axes[2].legend(fontsize=8); axes[2].grid(True, which="both", lw=0.3)
    axes[2].set_title("SYNTHETIC tracker noise: rms", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "pole_placement.png", dpi=150)
    plt.close(fig)
    print(f"saved {out_dir / 'pole_placement.png'}")


if __name__ == "__main__":
    main()
