"""P09.16 LQR: starter code.

The linear-quadratic regulator for a single-axis satellite attitude loop driven by a reaction wheel.
The continuous algebraic Riccati equation (CARE)

    A^T P + P A - P B R^{-1} B^T P + Q = 0

is solved from scratch two independent ways: from the stable invariant subspace of the Hamiltonian
matrix, and by Newton-Kleinman iteration (a sequence of Lyapunov equations, each solved by a
Kronecker-product linear solve). Both are checked against ``scipy.linalg.solve_continuous_are`` in the
tests. Then: the optimal gain K = R^{-1} B^T P, the optimal cost x0^T P x0, Bryson's rule for the
weights, the Kalman return-difference inequality and the guaranteed margins of single-input LQR
(gain margin [1/2, inf), phase margin >= 60 deg), and a simulation with the wheel's torque limit.

Conventions:
    * Satellite state x = (theta, omega): attitude error [rad] and body rate [rad/s] about one axis.
    * Input u = wheel torque on the body [N m]; the wheel absorbs the opposite torque, so the wheel
      momentum change is h_w(t) = -J (omega(t) - omega(0)) (total angular momentum is conserved).
    * u = -K x is the control law; A - B K is the closed-loop matrix.

Complete every function marked NotImplementedError.

Run:  python main.py          (prints the satellite design, writes lqr_satellite.png)
Test: pytest curriculum/phase-09/lesson-16-lqr
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


# ---------------------------------------------------------------- the plant (given)

def satellite_model(J: float = 10.0) -> tuple[np.ndarray, np.ndarray]:
    """Return (A, B) of the rigid single-axis attitude model J theta'' = u, state (theta, omega).

    J = 10 kg m^2 is representative of a 50-100 kg small satellite about one axis.
    """
    A = np.array([[0.0, 1.0], [0.0, 0.0]])
    B = np.array([[0.0], [1.0 / J]])
    return A, B


def wheel_limits() -> dict[str, float]:
    """Return representative small-satellite reaction-wheel limits: torque 0.02 N m, momentum 0.2 N m s."""
    return {"torque": 0.02, "momentum": 0.2}


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


# ---------------------------------------------------------------- Riccati and Lyapunov solvers

def lyapunov_solve(A: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Return X solving A^T X + X A + Q = 0, by the Kronecker form (I kron A^T + A^T kron I) vec X = -vec Q.

    Unique when no two eigenvalues of A sum to zero (in particular when A is Hurwitz). Raise ``ValueError``
    if some lambda_i + lambda_j is zero to within 1e-10 max(1, max|lambda|): the Kronecker matrix, whose
    eigenvalues are exactly these sums, is then singular. The result is symmetrised when Q is symmetric.
    """
    raise NotImplementedError


def hamiltonian_matrix(A: np.ndarray, B: np.ndarray, Q: np.ndarray, R: np.ndarray) -> np.ndarray:
    """Return the 2n x 2n Hamiltonian matrix H = [[A, -B R^{-1} B^T], [-Q, -A^T]] of the LQR problem."""
    raise NotImplementedError


def care_hamiltonian(A: np.ndarray, B: np.ndarray, Q: np.ndarray, R: np.ndarray) -> np.ndarray:
    """Return the stabilising solution P of the CARE from the stable invariant subspace of H.

    Take the n eigenvectors of H whose eigenvalues have negative real part, stack them as [X1; X2], and
    return P = X2 X1^{-1} (real part, symmetrised). Raise ``ValueError`` if H has eigenvalues on the
    imaginary axis (no stabilising solution: (A, B) not stabilisable or (A, Q) has an unobservable mode on
    the axis) or X1 is singular.
    """
    raise NotImplementedError


def bass_stabilizing_gain(A: np.ndarray, B: np.ndarray, beta: float | None = None) -> np.ndarray:
    """Return a stabilising gain K0 = B^T W^{-1} by Bass's method.

    W solves (A + beta I) W + W (A + beta I)^T = 2 B B^T; for controllable (A, B) and beta larger than the
    largest |Re lambda(A)|, every eigenvalue of A - B K0 has real part at most -beta. Default
    beta = 1 + ||A||_1. Used to start Newton-Kleinman when A itself is not Hurwitz.
    """
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    if beta is None:
        beta = 1.0 + np.linalg.norm(A, 1)
    Ab = A + beta * np.eye(A.shape[0])
    # (A + bI) W + W (A + bI)^T = 2 B B^T  <=>  M^T W + W M + Qw = 0 with M = -(A + bI)^T, Qw = 2 B B^T
    W = lyapunov_solve(-Ab.T, 2.0 * B @ B.T)
    return B.T @ np.linalg.inv(W)


def care_newton_kleinman(A: np.ndarray, B: np.ndarray, Q: np.ndarray, R: np.ndarray,
                         K0: np.ndarray | None = None, tol: float = 1e-12,
                         max_iter: int = 100) -> tuple[np.ndarray, int]:
    """Return (P, iterations): the stabilising CARE solution by Newton-Kleinman iteration.

    From a stabilising K_k, solve the Lyapunov equation
        (A - B K_k)^T P_k + P_k (A - B K_k) + Q + K_k^T R K_k = 0,
    then set K_{k+1} = R^{-1} B^T P_k. Stop when ||P_k - P_{k-1}||_F <= tol ||P_k||_F. K0 defaults to
    zero if A is Hurwitz and to ``bass_stabilizing_gain`` otherwise. Raise ``ValueError`` if K0, or any later
    iterate, is not stabilising (in floating point a huge K0 can lead the iteration to a non-stabilising root
    of the ARE whose residual is just as small). Convergence is monotone (P_k decreases) and quadratic near the solution. The iteration also
    stops when the update has reached rounding level (below 1e-9 ||P_k||_F) and has stopped shrinking.
    """
    raise NotImplementedError


def are_residual(A: np.ndarray, B: np.ndarray, Q: np.ndarray, R: np.ndarray, P: np.ndarray) -> float:
    """Return ||A^T P + P A - P B R^{-1} B^T P + Q||_F / max(||Q||_F, ||A^T P + P A||_F): a relative CARE residual."""
    raise NotImplementedError


# ---------------------------------------------------------------- the regulator

def lqr(A: np.ndarray, B: np.ndarray, Q: np.ndarray, R: np.ndarray,
        method: str = "hamiltonian") -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (K, P, closed-loop eigenvalues) of the infinite-horizon LQR, K = R^{-1} B^T P.

    ``method`` is "hamiltonian" or "newton".
    """
    raise NotImplementedError


def bryson_weights(x_max: np.ndarray | list[float], u_max: np.ndarray | list[float]) -> tuple[np.ndarray, np.ndarray]:
    """Return (Q, R) by Bryson's rule: Q = diag(1/x_max_i^2), R = diag(1/u_max_j^2)."""
    raise NotImplementedError


def cost_of_gain(A: np.ndarray, B: np.ndarray, Q: np.ndarray, R: np.ndarray, K: np.ndarray,
                 x0: np.ndarray) -> float:
    """Return J(K) = int_0^inf (x^T Q x + u^T R u) dt for u = -K x from x0, as x0^T X x0.

    X solves (A - B K)^T X + X (A - B K) + Q + K^T R K = 0. Return ``inf`` if A - B K is not Hurwitz.
    """
    raise NotImplementedError


def loop_transfer(A: np.ndarray, B: np.ndarray, K: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Return the loop transfer at the plant input L(jw) = K (jw I - A)^{-1} B for a single input, shape (len(w),)."""
    raise NotImplementedError


def phase_margin(w: np.ndarray, L: np.ndarray) -> tuple[float, float]:
    """Return (phase margin [deg], crossover frequency [rad/s]) from samples of L(jw) on an increasing grid.

    The crossover is where |L| falls through 1 (linear interpolation in log|L| against log w); the margin is
    180 deg + angle L there. Return (inf, nan) if |L| never crosses 1.
    """
    mag = np.abs(L)
    idx = np.where((mag[:-1] >= 1.0) & (mag[1:] < 1.0))[0]
    if len(idx) == 0:
        return float("inf"), float("nan")
    i = idx[0]
    lw0, lw1 = np.log(w[i]), np.log(w[i + 1])
    lm0, lm1 = np.log(mag[i]), np.log(mag[i + 1])
    f = lm0 / (lm0 - lm1)
    wc = float(np.exp(lw0 + f * (lw1 - lw0)))
    ph = np.unwrap(np.angle(L))
    phc = ph[i] + f * (ph[i + 1] - ph[i])
    pm = np.degrees(phc) + 180.0
    pm = (pm + 180.0) % 360.0 - 180.0
    return float(pm), wc


def gain_margin_interval(A: np.ndarray, B: np.ndarray, K: np.ndarray, k_grid: np.ndarray) -> tuple[float, float]:
    """Return (k_low, k_high): the smallest and largest k in ``k_grid`` (sorted, containing 1) for which every
    gain from that k up to 1, or from 1 up to it, keeps A - k B K Hurwitz. Scans outward from k = 1."""
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    K = np.atleast_2d(K)
    k_grid = np.sort(np.asarray(k_grid, dtype=float))

    def stable(k: float) -> bool:
        return bool(np.all(np.linalg.eigvals(A - k * B @ K).real < 0))

    i1 = int(np.argmin(np.abs(k_grid - 1.0)))
    lo = k_grid[i1]
    for k in k_grid[i1::-1]:
        if not stable(k):
            break
        lo = k
    hi = k_grid[i1]
    for k in k_grid[i1:]:
        if not stable(k):
            break
        hi = k
    return float(lo), float(hi)


def simulate(A: np.ndarray, B: np.ndarray, K: np.ndarray, x0: np.ndarray, t_end: float, dt: float,
             u_max: float | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (t, x, u) for u = -K x (clipped to +-u_max if given) by RK4 with fixed step dt.

    The control is evaluated inside every RK4 stage, so without a limit this integrates the linear closed
    loop dx/dt = (A - B K) x to fourth order; with a limit it is the saturated (nonlinear) loop.
    """
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    K = np.atleast_2d(K)
    n_steps = int(round(t_end / dt))
    t = np.arange(n_steps + 1) * dt
    x = np.zeros((n_steps + 1, A.shape[0]))
    x[0] = x0

    def ctrl(xx: np.ndarray) -> np.ndarray:
        u = -K @ xx
        if u_max is not None:
            u = np.clip(u, -u_max, u_max)
        return u

    def f(xx: np.ndarray) -> np.ndarray:
        return A @ xx + B @ ctrl(xx)

    for k in range(n_steps):
        xk = x[k]
        k1 = f(xk)
        k2 = f(xk + 0.5 * dt * k1)
        k3 = f(xk + 0.5 * dt * k2)
        k4 = f(xk + dt * k3)
        x[k + 1] = xk + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
    u = np.array([ctrl(xx) for xx in x])
    return t, x, u


def quadratic_cost(t: np.ndarray, x: np.ndarray, u: np.ndarray, Q: np.ndarray, R: np.ndarray) -> float:
    """Return the trapezoidal integral of x^T Q x + u^T R u over the sampled trajectory."""
    integrand = np.einsum("ij,jk,ik->i", x, Q, x) + np.einsum("ij,jk,ik->i", u, np.atleast_2d(R), u)
    return float(np.sum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(t)))


# ---------------------------------------------------------------- the example

def weight_sweep(rhos: np.ndarray, J: float = 10.0, theta0_deg: float = 5.0) -> list[dict[str, float]]:
    """For R = rho R_Bryson, return per rho: closed-loop poles' slowest real part, 2 % settling time of theta
    (from a simulation), the peak torque |u(0)| for a theta0 error, and the phase margin."""
    A, B = satellite_model(J)
    lim = wheel_limits()
    Q, R0 = bryson_weights([np.radians(5.0), np.radians(0.5)], [lim["torque"]])
    x0 = np.array([np.radians(theta0_deg), 0.0])
    w = np.logspace(-4, 2, 3000)
    rows = []
    for rho in rhos:
        K, P, lam = lqr(A, B, Q, rho * R0)
        t, x, u = simulate(A, B, K, x0, 400.0, 0.05)
        out = np.where(np.abs(x[:, 0]) > 0.02 * abs(x0[0]))[0]
        ts = float(t[min(out[-1] + 1, len(t) - 1)]) if len(out) else 0.0
        pm, wc = phase_margin(w, loop_transfer(A, B, K, w))
        rows.append({"rho": float(rho), "slowest": float(lam.real.max()), "settle": ts,
                     "peak_torque": float(np.abs(u).max()), "pm": pm, "wc": wc})
    return rows


def main(out_dir: str | Path = ".") -> None:
    """Print the satellite LQR design and its checks, and save lqr_satellite.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    J = 10.0
    A, B = satellite_model(J)
    lim = wheel_limits()
    th_max, om_max = np.radians(5.0), np.radians(0.5)
    Q, R = bryson_weights([th_max, om_max], [lim["torque"]])
    print(f"satellite: J = {J} kg m^2, wheel torque limit {lim['torque']} N m, momentum {lim['momentum']} N m s")
    print(f"Bryson: theta_max = 5 deg, omega_max = 0.5 deg/s, u_max = {lim['torque']} N m ->"
          f" Q = diag({Q[0, 0]:.2f}, {Q[1, 1]:.1f}), R = {R[0, 0]:.0f}")

    P_h = care_hamiltonian(A, B, Q, R)
    P_n, iters = care_newton_kleinman(A, B, Q, R)
    print(f"P (Hamiltonian) = [[{P_h[0, 0]:.5f}, {P_h[0, 1]:.5f}], [{P_h[1, 0]:.5f}, {P_h[1, 1]:.4f}]]")
    print(f"Newton-Kleinman: {iters} iterations from the Bass gain, max |P_n - P_h| / max|P_h| = {np.abs(P_n - P_h).max() / np.abs(P_h).max():.1e}")
    print(f"ARE residuals: Hamiltonian {are_residual(A, B, Q, R, P_h):.1e}, Newton {are_residual(A, B, Q, R, P_n):.1e}")
    K, P, lam = lqr(A, B, Q, R)
    zeta = -lam.real / np.abs(lam)
    print(f"K = [{K[0, 0]:.4f} N m/rad, {K[0, 1]:.4f} N m s/rad]; closed-loop poles {np.round(lam, 4)} rad/s")
    # hand formula for the double integrator: K = [sqrt(q1/r), sqrt(J^2... )] (see lesson 2.4)
    q1, q2, r = Q[0, 0], Q[1, 1], R[0, 0]
    k1 = np.sqrt(q1 / r)
    k2 = np.sqrt(q2 / r + 2.0 * J * k1)
    print(f"hand formula: k1 = sqrt(q1/r) = {k1:.4f}, k2 = sqrt(q2/r + 2 J k1) = {k2:.4f}; zeta = {zeta.min():.3f}")

    x0 = np.array([np.radians(5.0), 0.0])
    t, x, u = simulate(A, B, K, x0, 200.0, 0.05)
    J_sim = quadratic_cost(t, x, u, Q, R)
    print(f"5 deg error: peak torque {np.abs(u).max():.4f} N m, cost x0'Px0 = {x0 @ P @ x0:.4f}, simulated cost {J_sim:.4f}")
    out = np.where(np.abs(x[:, 0]) > 0.02 * x0[0])[0]
    print(f"   2 % settling {t[out[-1] + 1]:.1f} s, peak body rate {np.degrees(np.abs(x[:, 1]).max()):.3f} deg/s,"
          f" peak wheel momentum change {J * np.abs(x[:, 1]).max():.4f} N m s")
    for scale in (0.5, 2.0):
        print(f"   cost with gain {scale} K: {cost_of_gain(A, B, Q, R, scale * K, x0):.4f}")

    w = np.logspace(-4, 2, 4000)
    L = loop_transfer(A, B, K, w)
    pm, wc = phase_margin(w, L)
    print(f"loop at the plant input: phase margin {pm:.1f} deg at {wc:.4f} rad/s; min |1 + L| = {np.abs(1 + L).min():.4f};"
          f" delay margin {np.radians(pm) / wc:.1f} s")
    lo, hi = gain_margin_interval(A, B, K, np.concatenate([np.linspace(0.01, 1.0, 100), np.logspace(0, 4, 200)]))
    print(f"gain scan: stable for k in [{lo:.2f}, {hi:.0f}] (grid limits 0.01 and 1e4; theory [1/2, inf), here the double integrator keeps any k > 0)")

    # saturation: a 20 deg error with the 0.02 N m limit
    x20 = np.array([np.radians(20.0), 0.0])
    t2, xs, us = simulate(A, B, K, x20, 300.0, 0.05, u_max=lim["torque"])
    t2l, xl, ul = simulate(A, B, K, x20, 300.0, 0.05)
    def settle(tt: np.ndarray, th: np.ndarray) -> float:
        o = np.where(np.abs(th) > 0.02 * abs(th[0]))[0]
        return float(tt[o[-1] + 1])

    print(f"20 deg error: linear loop asks {np.abs(ul).max():.4f} N m; 2 % settling torque-limited {settle(t2, xs[:, 0]):.1f} s"
          f" vs unlimited {settle(t2l, xl[:, 0]):.1f} s; time at the limit {np.sum(np.abs(us[:, 0]) >= lim['torque'] - 1e-12) * 0.05:.1f} s;"
          f" peak rate {np.degrees(np.abs(xs[:, 1]).max()):.3f} vs {np.degrees(np.abs(xl[:, 1]).max()):.3f} deg/s,"
          f" peak wheel momentum change {J * np.abs(xs[:, 1]).max():.4f} vs {J * np.abs(xl[:, 1]).max():.4f} N m s (capacity {lim['momentum']})")

    rows = weight_sweep(np.array([0.1, 0.3, 1.0, 3.0, 10.0]))
    print("weight sweep R = rho R_Bryson: rho, slowest pole, 2% settling [s], peak torque [N m], PM [deg], wc [rad/s]")
    for row in rows:
        print(f"   {row['rho']:5.1f} {row['slowest']:8.4f} {row['settle']:7.1f} {row['peak_torque']:.4f} {row['pm']:6.1f} {row['wc']:.4f}")

    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    axes[0].plot(t, np.degrees(x[:, 0]), label="5 deg, linear")
    axes[0].plot(t2, np.degrees(xs[:, 0]), label="20 deg, torque-limited")
    axes[0].plot(t2l, np.degrees(xl[:, 0]), "--", label="20 deg, unlimited")
    axes[0].set_xlabel("t [s]"); axes[0].set_ylabel("theta [deg]"); axes[0].legend(fontsize=8); axes[0].grid(lw=0.4)
    axes[1].plot(t2, us[:, 0], label="torque-limited")
    axes[1].plot(t2l, ul[:, 0], "--", label="unlimited")
    axes[1].axhline(lim["torque"], color="k", lw=0.6); axes[1].axhline(-lim["torque"], color="k", lw=0.6)
    axes[1].set_xlabel("t [s]"); axes[1].set_ylabel("wheel torque [N m]"); axes[1].legend(fontsize=8); axes[1].grid(lw=0.4)
    axes[2].plot((1 + L).real, (1 + L).imag, label="1 + L(jw)")
    th = np.linspace(0, 2 * np.pi, 200)
    axes[2].plot(np.cos(th), np.sin(th), "k:", lw=0.8, label="unit circle about 0")
    axes[2].set_aspect("equal"); axes[2].set_xlim(-0.5, 4); axes[2].set_ylim(-3, 1)
    axes[2].set_title("Kalman: |1 + L| >= 1", fontsize=9); axes[2].legend(fontsize=8); axes[2].grid(lw=0.4)
    fig.tight_layout()
    fig.savefig(out_dir / "lqr_satellite.png", dpi=150)
    plt.close(fig)
    print(f"saved {out_dir / 'lqr_satellite.png'}")


if __name__ == "__main__":
    main()
