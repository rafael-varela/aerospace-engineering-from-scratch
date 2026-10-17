"""P02.05 Matrix Exponential & Linear Systems: starter code.

The state transition matrix Phi(t) = e^{At} of a linear system dx/dt = A x + B u, built from scratch
(Taylor series, scaling and squaring, eigen-decomposition) and checked against ``scipy.linalg.expm``;
its group, inverse, derivative and determinant properties; the forced solution by variation of
parameters and the augmented-matrix trick for a constant input; zero-order-hold discretisation; and
stability from the eigenvalues. The example is the Clohessy-Wiltshire (Hill) relative-motion model of a
chaser spacecraft near a target in circular orbit, whose closed-form transition matrix is the oracle,
and a two-impulse rendezvous solved as a linear system.

Conventions:
    * CW state x = (x, y, xdot, ydot): x radial (outward), y along-track (direction of motion), in the
      target's rotating frame; n = sqrt(mu / a^3) is the orbital mean motion [rad/s].
    * Phi(t) maps x(0) to x(t): x(t) = Phi(t) x(0).

Complete every function marked NotImplementedError.

Run:  python main.py          (prints the rendezvous example, writes state_transition.png)
Test: pytest curriculum/phase-02/lesson-05-matrix-exponential
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

MU_EARTH = 3.986004418e14  # m^3/s^2


def iss_orbit() -> dict[str, float]:
    """Return the circular-orbit parameters of a station-class target: a = 6778 km (400 km altitude), mu of Earth, mean motion n."""
    a = 6778e3
    return {"a": a, "mu": MU_EARTH, "n": float(np.sqrt(MU_EARTH / a ** 3))}


# ---------------------------------------------------------------- the matrix exponential, from scratch

def expm_taylor(A: np.ndarray, n_terms: int = 30) -> np.ndarray:
    """Return sum_{k=0}^{n_terms-1} A^k / k!, the truncated Taylor series of e^A (accurate only when ||A|| is small)."""
    raise NotImplementedError


def expm_scaling_squaring(A: np.ndarray, n_terms: int = 20) -> np.ndarray:
    """Return e^A by scaling and squaring: pick s with ||A||/2^s <= 0.5, take the Taylor series of e^{A/2^s}, square s times.

    Uses the identity e^A = (e^{A/2^s})^(2^s); the Taylor series converges fast for the scaled matrix and the
    squarings undo the scaling exactly. This is the core idea behind every library ``expm``.
    """
    raise NotImplementedError


def expm_eig(A: np.ndarray, cond_limit: float = 1e8) -> np.ndarray:
    """Return e^A = V diag(e^lambda) V^{-1} from the eigen-decomposition (real part), for a diagonalisable A.

    Raise ``ValueError`` if the eigenvector matrix has condition number above ``cond_limit``: a defective or
    nearly defective A (repeated eigenvalue with too few eigenvectors, like Clohessy-Wiltshire) cannot be
    exponentiated this way. The lesson explains why the other two methods do not care.
    """
    raise NotImplementedError


def state_transition(A: np.ndarray, t: float) -> np.ndarray:
    """Return Phi(t) = e^{At} by ``expm_scaling_squaring``."""
    raise NotImplementedError


def propagate(A: np.ndarray, x0: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Return x(t_i) = Phi(t_i) x0 for each t_i, shape (len(t), n), for the homogeneous system dx/dt = A x."""
    raise NotImplementedError


def is_stable(A: np.ndarray) -> bool:
    """Return True if every eigenvalue of A has strictly negative real part (asymptotic stability of dx/dt = A x)."""
    raise NotImplementedError


# ---------------------------------------------------------------- forced systems

def forced_response_zoh(A: np.ndarray, B: np.ndarray, x0: np.ndarray, u: np.ndarray, dt: float) -> np.ndarray:
    """Return the states x_0, x_1, ..., x_N of dx/dt = A x + B u with u held constant at u[k] over [k dt, (k+1) dt].

    Uses the augmented matrix M = [[A, B], [0, 0]]: e^{M dt} = [[A_d, B_d], [0, I]] with A_d = e^{A dt} and
    B_d = int_0^dt e^{A s} ds B, so x_{k+1} = A_d x_k + B_d u_k exactly (zero-order hold). Returns shape (N+1, n).
    """
    raise NotImplementedError


def discretize(A: np.ndarray, B: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (A_d, B_d) of the exact zero-order-hold discretisation, from the exponential of the augmented matrix."""
    raise NotImplementedError


def modal_decomposition(A: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (eigenvalues, V, V^{-1}) so that x(t) = sum_k (V^{-1} x0)_k e^{lambda_k t} v_k: the normal modes of dx/dt = A x."""
    raise NotImplementedError


# ---------------------------------------------------------------- Clohessy-Wiltshire relative motion

def cw_matrix(n: float) -> np.ndarray:
    """Return A of the in-plane Clohessy-Wiltshire equations for x = (x, y, xdot, ydot):
    xddot = 3 n^2 x + 2 n ydot,  yddot = -2 n xdot.

    Linearised relative motion of a chaser about a target in a circular orbit of mean motion n (radial x,
    along-track y). A has eigenvalues 0, 0, +-i n and is defective: the double zero carries the secular drift.
    """
    return np.array([[0.0, 0.0, 1.0, 0.0],
                     [0.0, 0.0, 0.0, 1.0],
                     [3.0 * n * n, 0.0, 0.0, 2.0 * n],
                     [0.0, 0.0, -2.0 * n, 0.0]])


def cw_stm_closed_form(n: float, t: float) -> np.ndarray:
    """Return the closed-form Clohessy-Wiltshire state transition matrix Phi(t) for x = (x, y, xdot, ydot).

    x(t)    = (4 - 3 cos nt) x0 + (sin nt / n) xdot0 + (2 (1 - cos nt) / n) ydot0
    y(t)    = 6 (sin nt - nt) x0 + y0 - (2 (1 - cos nt) / n) xdot0 + ((4 sin nt - 3 nt) / n) ydot0
    xdot(t) = 3 n sin nt x0 + cos nt xdot0 + 2 sin nt ydot0
    ydot(t) = 6 n (cos nt - 1) x0 - 2 sin nt xdot0 + (4 cos nt - 3) ydot0
    """
    raise NotImplementedError


def cw_rhs(t: float, x: np.ndarray, n: float, accel: np.ndarray | None = None) -> np.ndarray:
    """Return dx/dt for the CW equations with an optional constant thrust acceleration (ax, ay) [m/s^2]."""
    a = np.zeros(2) if accel is None else np.asarray(accel, dtype=float)
    return np.array([x[2], x[3], 3.0 * n * n * x[0] + 2.0 * n * x[3] + a[0], -2.0 * n * x[2] + a[1]])


def rendezvous_delta_v(n: float, r0: np.ndarray, v0: np.ndarray, T: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (dv1, dv2): the two impulsive velocity changes [m/s] that take the chaser from (r0, v0) to the origin at rest in time T.

    Partition Phi(T) into blocks: r(T) = Phi_rr r0 + Phi_rv v0' with v0' = v0 + dv1. Setting r(T) = 0 gives the
    linear system Phi_rv dv1 = -(Phi_rr r0 + Phi_rv v0), solved with ``numpy.linalg.solve``; dv2 then cancels
    the arrival velocity. Raise ``ValueError`` if Phi_rv is singular (T a multiple of the orbital period).
    """
    raise NotImplementedError


# ---------------------------------------------------------------- the example

def main(out_dir: str | Path = ".") -> None:
    """Print the state-transition checks and the rendezvous, and save state_transition.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.integrate import solve_ivp
    from scipy.linalg import expm

    out_dir = Path(out_dir)
    orb = iss_orbit()
    n = orb["n"]
    P = 2 * np.pi / n
    print(f"target orbit: a = {orb['a'] / 1e3:.0f} km, n = {n:.6e} rad/s, period {P / 60:.2f} min")
    A = cw_matrix(n)
    print("CW eigenvalues:", np.round(np.linalg.eigvals(A), 8), " cond(V) =", f"{np.linalg.cond(np.linalg.eig(A)[1]):.2e}")

    # Three ways to the state transition matrix for half an orbit.
    T = P / 2
    Phi_ss = state_transition(A, T)
    Phi_cf = cw_stm_closed_form(n, T)
    Phi_sp = expm(A * T)
    print(f"Phi(P/2): max |scaling-squaring - closed form| = {np.max(np.abs(Phi_ss - Phi_cf)):.2e},"
          f" max |scipy expm - closed form| = {np.max(np.abs(Phi_sp - Phi_cf)):.2e}")
    print(f"plain Taylor with 30 terms at ||A T|| = {np.linalg.norm(A * T, 1):.2f}: error {np.max(np.abs(expm_taylor(A * T) - Phi_cf)):.2e};"
          f" with 60 terms: {np.max(np.abs(expm_taylor(A * T, 60) - Phi_cf)):.2e}")
    try:
        expm_eig(A * T)
    except ValueError as exc:
        print("expm_eig on CW:", exc)
    print(f"det Phi(P/2) = {np.linalg.det(Phi_cf):.12f} (= e^{{tr(A) T}} = 1), Phi(P/2) Phi(P/2) - Phi(P) max = {np.max(np.abs(Phi_cf @ Phi_cf - cw_stm_closed_form(n, P))):.2e}")

    # A chaser 1 km behind, 50 m below, at rest in the rotating frame: where is it after half an orbit?
    x0 = np.array([-50.0, -1000.0, 0.0, 0.0])
    t = np.linspace(0.0, P, 721)
    xs = propagate(A, x0, t)
    xT = Phi_cf @ x0
    print(f"\nchaser from (x, y) = (-50, -1000) m at rest: after P/2 at ({xT[0]:.1f}, {xT[1]:.1f}) m moving ({xT[2]:.3f}, {xT[3]:.3f}) m/s;"
          f" after one orbit at ({xs[-1, 0]:.1f}, {xs[-1, 1]:.1f}) m (drift {xs[-1, 1] - x0[1]:.1f} m per orbit)")
    sol = solve_ivp(lambda tt, x: cw_rhs(tt, x, n), (0.0, P), x0, rtol=1e-10, atol=1e-9, t_eval=t)
    print(f"max |propagate - solve_ivp| over one orbit: {np.max(np.abs(xs - sol.y.T)):.2e} m")

    # Two-impulse rendezvous in a quarter orbit and in half an orbit.
    for frac in (0.25, 0.5, 0.75):
        Tr = frac * P
        dv1, dv2 = rendezvous_delta_v(n, x0[:2], x0[2:], Tr)
        x_after = np.concatenate([x0[:2], x0[2:] + dv1])
        arrive = cw_stm_closed_form(n, Tr) @ x_after
        print(f"rendezvous in {frac:.2f} orbit ({Tr / 60:.1f} min): dv1 = ({dv1[0]:.4f}, {dv1[1]:.4f}) m/s, dv2 = ({dv2[0]:.4f}, {dv2[1]:.4f}) m/s,"
              f" total {np.linalg.norm(dv1) + np.linalg.norm(dv2):.4f} m/s; arrival position error {np.linalg.norm(arrive[:2]):.2e} m")

    # Constant thrust: ZOH vs solve_ivp, and the discretisation.
    accel = np.array([0.0, 2e-4])        # 0.2 mm/s^2 along-track for 10 minutes
    dt = 10.0
    N = 60
    B = np.array([[0.0, 0.0], [0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    u = np.tile(accel[:, None], (1, N))
    x_zoh = forced_response_zoh(A, B, x0, u, dt)
    sol_f = solve_ivp(lambda tt, x: cw_rhs(tt, x, n, accel), (0.0, N * dt), x0, rtol=1e-11, atol=1e-10, t_eval=np.arange(N + 1) * dt)
    print(f"\n10 min of 0.2 mm/s^2 along-track thrust: end state ZOH ({x_zoh[-1, 0]:.3f}, {x_zoh[-1, 1]:.3f}) m, solve_ivp ({sol_f.y[0, -1]:.3f}, {sol_f.y[1, -1]:.3f}) m,"
          f" max difference {np.max(np.abs(x_zoh - sol_f.y.T)):.2e} m")
    Ad, Bd = discretize(A, B, dt)
    print(f"A_d(10 s) = Phi(10 s): {np.max(np.abs(Ad - cw_stm_closed_form(n, dt))):.2e}; B_d ~ dt B for small n dt: max |B_d - dt B| = {np.max(np.abs(Bd - dt * B)):.2e}")
    print("stable?", is_stable(A), " (CW: marginal, eigenvalues on the imaginary axis; a damped short period:",
          is_stable(np.array([[-0.8, 1.0], [-4.0, -1.2]])), ")")

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    axes[0].plot(xs[:, 1] / 1e3, xs[:, 0] / 1e3, label="free drift, one orbit")
    for frac, col in ((0.25, "C1"), (0.5, "C2")):
        Tr = frac * P
        dv1, _ = rendezvous_delta_v(n, x0[:2], x0[2:], Tr)
        traj = propagate(A, np.concatenate([x0[:2], x0[2:] + dv1]), np.linspace(0.0, Tr, 200))
        axes[0].plot(traj[:, 1] / 1e3, traj[:, 0] / 1e3, col, label=f"rendezvous in {frac} orbit")
    axes[0].plot(0, 0, "k*", ms=10, label="target")
    axes[0].set_xlabel("along-track y [km]"); axes[0].set_ylabel("radial x [km]"); axes[0].legend(fontsize=8); axes[0].grid(lw=0.4)
    axes[0].set_title("Clohessy-Wiltshire relative motion", fontsize=9)
    ks = np.arange(2, 41)
    errs = [np.max(np.abs(expm_taylor(A * T, k) - Phi_cf)) for k in ks]
    axes[1].semilogy(ks, errs, "o-", ms=3, label=f"plain Taylor, ||A T|| = {np.linalg.norm(A * T, 1):.1f}")
    errs_ss = [np.max(np.abs(expm_scaling_squaring(A * T, k) - Phi_cf)) for k in ks]
    axes[1].semilogy(ks, errs_ss, "s-", ms=3, label="scaling and squaring")
    axes[1].set_xlabel("Taylor terms"); axes[1].set_ylabel("max |error| in Phi(P/2)"); axes[1].legend(fontsize=8)
    axes[1].grid(True, which="both", lw=0.4); axes[1].set_ylim(1e-17, 1e3)
    fig.tight_layout()
    fig.savefig(out_dir / "state_transition.png", dpi=150)
    print(f"saved {out_dir / 'state_transition.png'}")


if __name__ == "__main__":
    main()
