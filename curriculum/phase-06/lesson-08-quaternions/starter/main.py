"""P06.08 Quaternions: starter code.

A quaternion library from scratch: the Hamilton product, conjugate, norm and inverse; the unit
quaternion of an axis and angle; the transformation of vector components; the conversions
quaternion <-> DCM <-> 3-2-1 Euler angles; the kinematic equation dq/dt = q * (0, w) / 2 and its
integration (solve_ivp, fixed-step RK4 with or without renormalisation, and the exact exponential step);
the attitude error between two quaternions; and SLERP, with the normalised linear interpolation it
improves on. Everything is checked against the P06.07 DCM route, re-implemented here in the helpers
marked "carried over from P06.07" so that this lesson stands alone.

Conventions (the passive ones of P01.15 and P06.07): q = [q0, q1, q2, q3] = [cos(a/2), sin(a/2) n] is
the attitude of the body frame b relative to the reference frame n: the frame b is n turned by the angle a
about the unit axis n (the same in both frames). Its DCM is

    C_b/n(q) = (q0^2 - v.v) E + 2 v v^T - 2 q0 [v x],    v = (q1, q2, q3),

equal to P01.15's passive Rodrigues matrix. Components transform as v_b = q^* (0, v_n) q. Successive turns
compose as q_c/n = q_b/n * q_c/b (the turn about the newer frame on the right), mirroring
C_c/n = C_c/b C_b/n. w is the body rate of b relative to n in body components, and dq/dt = q * (0, w) / 2.
q and -q are the same attitude.

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes quaternions.png)
Test: pytest curriculum/phase-06/lesson-08-quaternions
"""
from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

RateFn = Callable[[float], np.ndarray]
Q_ID = np.array([1.0, 0.0, 0.0, 0.0])


# ---------------------------------------------------------------- carried over from P06.07 (given)

def rot(axis: int, angle: float) -> np.ndarray:
    """Passive elementary rotation about axis 1, 2 or 3 (P01.15, P06.07)."""
    c, s = math.cos(angle), math.sin(angle)
    if axis == 1:
        return np.array([[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]])
    if axis == 2:
        return np.array([[c, 0.0, -s], [0.0, 1.0, 0.0], [s, 0.0, c]])
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])


def skew(v: np.ndarray) -> np.ndarray:
    """The cross-product matrix [v x]."""
    return np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])


def dcm_from_euler_321(psi: float, theta: float, phi: float) -> np.ndarray:
    """C_b/n = R1(phi) R2(theta) R3(psi) (P06.07)."""
    return rot(1, phi) @ rot(2, theta) @ rot(3, psi)


def dcm_from_axis_angle(axis: np.ndarray, angle: float) -> np.ndarray:
    """Passive Rodrigues: cos(a) E + (1 - cos(a)) n n^T - sin(a) [n x] (P01.15)."""
    n = np.asarray(axis, dtype=float)
    n = n / np.linalg.norm(n)
    c, s = math.cos(angle), math.sin(angle)
    return c * np.eye(3) + (1.0 - c) * np.outer(n, n) - s * skew(n)


def misalignment_angle(Ca: np.ndarray, Cb: np.ndarray) -> float:
    """Rotation angle of Ca Cb^T (P06.07)."""
    E = Ca @ Cb.T
    w = np.array([E[2, 1] - E[1, 2], E[0, 2] - E[2, 0], E[1, 0] - E[0, 1]])
    return math.atan2(float(np.linalg.norm(w)), float(np.trace(E) - 1.0))


def integrate_dcm(C0: np.ndarray, w_fn: RateFn, t_end: float, n_samples: int = 201,
                  rtol: float = 1e-10, atol: float = 1e-12) -> tuple[np.ndarray, np.ndarray]:
    """dC/dt = -[w x] C by solve_ivp (DOP853), the P06.07 route. Returns (t, C (n, 3, 3))."""
    def rhs(t, y):
        return (-skew(w_fn(t)) @ y.reshape(3, 3)).ravel()
    t_eval = np.linspace(0.0, t_end, n_samples)
    sol = solve_ivp(rhs, (0.0, t_end), np.asarray(C0, dtype=float).ravel(), method="DOP853",
                    t_eval=t_eval, rtol=rtol, atol=atol)
    return sol.t, sol.y.T.reshape(-1, 3, 3)


def coning_rate(psi_dot: float, phi_dot: float, beta: float) -> RateFn:
    """Body rates of steady 3-1-3 coning (P06.07, eq. 6)."""
    def w(t: float) -> np.ndarray:
        ph = phi_dot * t
        return np.array([psi_dot * math.sin(beta) * math.sin(ph), psi_dot * math.sin(beta) * math.cos(ph),
                         phi_dot + psi_dot * math.cos(beta)])
    return w


def coning_dcm(t: float, psi_dot: float, phi_dot: float, beta: float) -> np.ndarray:
    """Closed-form coning attitude R3(phi_dot t) R1(beta) R3(psi_dot t) (P06.07)."""
    return rot(3, phi_dot * t) @ rot(1, beta) @ rot(3, psi_dot * t)


# ---------------------------------------------------------------- quaternion algebra

def quat_multiply(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Hamilton product p * q = (p0 q0 - pv.qv, p0 qv + q0 pv + pv x qv)."""
    raise NotImplementedError


def quat_conjugate(q: np.ndarray) -> np.ndarray:
    """q^* = (q0, -qv)."""
    raise NotImplementedError


def quat_norm(q: np.ndarray) -> float:
    """|q| = sqrt(q0^2 + q1^2 + q2^2 + q3^2)."""
    raise NotImplementedError


def quat_normalize(q: np.ndarray) -> np.ndarray:
    """q / |q|."""
    raise NotImplementedError


def quat_inverse(q: np.ndarray) -> np.ndarray:
    """q^-1 = q^* / |q|^2, so that q * q^-1 = (1, 0, 0, 0)."""
    raise NotImplementedError


def quat_from_axis_angle(axis: np.ndarray, angle: float) -> np.ndarray:
    """Unit quaternion (cos(a/2), sin(a/2) n) of a turn by `angle` about `axis` (normalised here)."""
    raise NotImplementedError


def axis_angle_from_quat(q: np.ndarray) -> tuple[np.ndarray, float]:
    """Axis and angle in [0, pi] of a unit quaternion; q and -q give the same answer (the sign of q is
    chosen so that q0 >= 0). For the identity the axis is returned as (1, 0, 0)."""
    raise NotImplementedError


def transform(q: np.ndarray, v_n: np.ndarray) -> np.ndarray:
    """Body components of a vector given in the reference frame: (0, v_b) = q^* (0, v_n) q."""
    raise NotImplementedError


# ---------------------------------------------------------------- quaternion <-> DCM <-> Euler

def dcm_from_quat(q: np.ndarray) -> np.ndarray:
    """C_b/n = (q0^2 - v.v) E + 2 v v^T - 2 q0 [v x]. No normalisation: a non-unit q gives |q|^2 times a
    rotation, which is how a drifting norm shows up in the DCM."""
    raise NotImplementedError


def quat_from_dcm(C: np.ndarray) -> np.ndarray:
    """Unit quaternion (q0 >= 0) of a rotation matrix, by Shepperd's method: compute the largest of
    4 q0^2 = 1 + tr C, 4 q1^2 = 1 + 2 C11 - tr C, ... and divide the off-diagonal sums and differences by it,
    so that no square root of a small number and no division by a small number occurs."""
    raise NotImplementedError


def quat_from_euler_321(psi: float, theta: float, phi: float) -> np.ndarray:
    """q_b/n of the 3-2-1 sequence: the three turns composed in order, q3(psi) * q2(theta) * q1(phi)."""
    raise NotImplementedError


def euler_321_from_quat(q: np.ndarray) -> np.ndarray:
    """(psi, theta, phi) of a unit quaternion, from the DCM entries written in q:
    psi = atan2(2(q0 q3 + q1 q2), 1 - 2(q2^2 + q3^2)), theta = asin(2(q0 q2 - q1 q3)),
    phi = atan2(2(q0 q1 + q2 q3), 1 - 2(q1^2 + q2^2)). At gimbal lock (|theta| = 90 deg) phi is set to 0
    and the whole rotation goes into psi, as in P06.07."""
    raise NotImplementedError


# ---------------------------------------------------------------- kinematics

def quat_rate(q: np.ndarray, w: np.ndarray) -> np.ndarray:
    """dq/dt = q * (0, w) / 2 with w the body rate in body axes."""
    raise NotImplementedError


def integrate_quat(q0: np.ndarray, w_fn: RateFn, t_end: float, n_samples: int = 201,
                   rtol: float = 1e-10, atol: float = 1e-12) -> tuple[np.ndarray, np.ndarray]:
    """Integrate dq/dt = q * (0, w(t)) / 2 with solve_ivp (DOP853). The output is NOT renormalised, so
    |q| - 1 shows the integrator's drift. Returns (t, q (n, 4))."""
    raise NotImplementedError


def integrate_quat_rk4(q0: np.ndarray, w_fn: RateFn, t_end: float, dt: float,
                       normalize: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """Fixed-step RK4 on the quaternion kinematic equation, optionally renormalising after every step.
    Returns (t, q) at every step."""
    raise NotImplementedError


def quat_exp_step(q: np.ndarray, w: np.ndarray, dt: float) -> np.ndarray:
    """Exact update for a body rate w held constant over dt: q * (cos(|w| dt / 2), sin(|w| dt / 2) w/|w|)."""
    raise NotImplementedError


def quat_error_angle(qa: np.ndarray, qb: np.ndarray) -> float:
    """The rotation angle (rad, in [0, pi]) between two attitudes: the angle of qb^-1 qa, insensitive to the
    signs of qa and qb: 2 atan2(|vector part|, |scalar part|) of qb^* qa (both normalised)."""
    raise NotImplementedError


# ---------------------------------------------------------------- interpolation

def slerp(q0: np.ndarray, q1: np.ndarray, u: float) -> np.ndarray:
    """Spherical linear interpolation from q0 (u = 0) to q1 (u = 1) along the shorter arc: if q0.q1 < 0,
    q1 is replaced by -q1 (the same attitude). With cos(Om) = q0.q1:
    q(u) = (sin((1 - u) Om) q0 + sin(u Om) q1) / sin(Om); a constant rotation rate about a fixed axis."""
    raise NotImplementedError


def nlerp(q0: np.ndarray, q1: np.ndarray, u: float, shortest: bool = True) -> np.ndarray:
    """Normalised linear interpolation ((1 - u) q0 + u q1) / |...|, on the shorter arc unless
    shortest=False. Same path as SLERP, but not at constant rate."""
    raise NotImplementedError


def interpolate_samples(t_s: np.ndarray, q_s: np.ndarray, t: np.ndarray, method: str = "slerp") -> np.ndarray:
    """Interpolate sampled attitudes (a star tracker's output at times t_s) at the times t, sample interval
    by sample interval, with "slerp" or "nlerp". Returns q (len(t), 4)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the lesson's output

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save quaternions.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    deg = math.pi / 180.0
    psi, th, ph = 30 * deg, 10 * deg, -5 * deg
    q = quat_from_euler_321(psi, th, ph)
    ax_, ang = axis_angle_from_quat(q)
    print(f"attitude yaw 30, pitch 10, roll -5 deg: q = {np.round(q, 6)}, |q| = {quat_norm(q):.15f}")
    print(f"  axis {np.round(ax_, 4)}, angle {ang / deg:.4f} deg; |C(q) - C_321| = {np.abs(dcm_from_quat(q) - dcm_from_euler_321(psi, th, ph)).max():.1e}")
    print(f"  back from the DCM: {np.round(quat_from_dcm(dcm_from_euler_321(psi, th, ph)), 6)}; "
          f"Euler angles back {np.round(euler_321_from_quat(q) / deg, 6)} deg")
    g_b = transform(q, np.array([0.0, 0.0, 9.80665]))
    print(f"  gravity in body axes by q* g q: {np.round(g_b, 4)} m/s^2 (by the DCM: {np.round(dcm_euler := dcm_from_euler_321(psi, th, ph) @ [0, 0, 9.80665], 4)})")

    # quaternion route against the DCM route: coning spacecraft
    pd, sd, beta = 0.4, 2.0, 20 * deg
    w = coning_rate(pd, sd, beta)
    C0 = coning_dcm(0.0, pd, sd, beta)
    q0 = quat_from_dcm(C0)
    tq, Q = integrate_quat(q0, w, 30.0, n_samples=601)
    tc, Cs = integrate_dcm(C0, w, 30.0, n_samples=601)
    mis_qc = np.array([misalignment_angle(dcm_from_quat(quat_normalize(Q[m])), Cs[m]) for m in range(tq.size)])
    mis_qx = np.array([misalignment_angle(dcm_from_quat(quat_normalize(Q[m])), coning_dcm(tq[m], pd, sd, beta)) for m in range(tq.size)])
    print(f"coning 30 s, solve_ivp rtol 1e-10: quaternion vs DCM route max {mis_qc.max():.2e} rad; quaternion vs closed form "
          f"{mis_qx.max():.2e} rad; max | |q| - 1 | = {np.abs(np.linalg.norm(Q, axis=1) - 1).max():.2e}")
    print("RK4 on the quaternion, coning 30 s:")
    rows = []
    for dt in [0.1, 0.05, 0.025, 0.0125]:
        t1, Q1 = integrate_quat_rk4(q0, w, 30.0, dt)
        t2, Q2 = integrate_quat_rk4(q0, w, 30.0, dt, normalize=True)
        ex = coning_dcm(t1[-1], pd, sd, beta)
        e1 = misalignment_angle(dcm_from_quat(quat_normalize(Q1[-1])), ex)
        e2 = misalignment_angle(dcm_from_quat(Q2[-1]), ex)
        n1 = abs(quat_norm(Q1[-1]) - 1)
        o1 = float(np.linalg.norm(dcm_from_quat(Q1[-1]).T @ dcm_from_quat(Q1[-1]) - np.eye(3)))
        rows.append((dt, e1, n1, e2, o1))
        print(f"  dt {dt:7.4f}: error {e1:.3e} rad, | |q| - 1 | {n1:.3e} (DCM of the raw q off orthogonal by {o1:.2e}); "
              f"renormalised each step: error {e2:.3e} rad")
    print(f"  observed order of the error: {np.round([math.log2(rows[k][1] / rows[k + 1][1]) for k in range(3)], 2)}; "
          f"of the norm drift: {np.round([math.log2(rows[k][2] / rows[k + 1][2]) for k in range(3)], 2)}")
    # exact exponential step for a constant rate
    wc = np.array([0.1, 0.2, 0.2]) * 3 * 30 * deg     # 30 deg/s about (1, 2, 2)/3
    qe = Q_ID.copy()
    for _ in range(2000):
        qe = quat_exp_step(qe, wc, 0.01)
    print(f"exact step, constant 30 deg/s for 20 s in 2000 steps: error vs closed form "
          f"{quat_error_angle(qe, quat_from_axis_angle(wc, np.linalg.norm(wc) * 20.0)):.1e} rad, | |q| - 1 | = {abs(quat_norm(qe) - 1):.1e}")

    # star tracker samples and interpolation
    print("star-tracker samples of the coning truth, interpolated (max error over 30 s, rad):")
    t_fine = np.linspace(0.0, 30.0, 3001)
    truth = [quat_from_dcm(coning_dcm(t, pd, sd, beta)) for t in t_fine]
    interp_rows = []
    for dts in [1.0, 0.5, 0.25, 0.125]:
        ts = np.arange(0.0, 30.0 + 1e-9, dts)
        qs = np.array([quat_from_dcm(coning_dcm(t, pd, sd, beta)) for t in ts])
        es = max(quat_error_angle(a, b) for a, b in zip(interpolate_samples(ts, qs, t_fine, "slerp"), truth))
        en = max(quat_error_angle(a, b) for a, b in zip(interpolate_samples(ts, qs, t_fine, "nlerp"), truth))
        interp_rows.append((dts, es, en))
        print(f"  sample interval {dts:5.3f} s: SLERP {es:.3e}, NLERP {en:.3e}")
    print(f"  SLERP order {np.round([math.log2(interp_rows[k][1] / interp_rows[k + 1][1]) for k in range(3)], 2)}")
    qa = Q_ID
    qb = quat_from_axis_angle([0, 0, 1], 170 * deg)
    qb_neg = -qb
    print(f"double cover: q = {np.round(qb, 4)} and -q are one attitude; halfway by SLERP: "
          f"{axis_angle_from_quat(slerp(qa, qb_neg, 0.5))[1] / deg:.1f} deg; by NLERP without the shortest-arc check "
          f"{axis_angle_from_quat(nlerp(qa, qb_neg, 0.5, shortest=False))[1] / deg:.1f} deg the other way round")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    ax[0].semilogy(tq, np.maximum(mis_qc, 1e-17), label="quaternion vs DCM route")
    ax[0].semilogy(tq, np.maximum(np.abs(np.linalg.norm(Q, axis=1) - 1), 1e-17), label="| |q| - 1 |")
    ax[0].set_xlabel("t (s)"); ax[0].legend(); ax[0].set_title("coning, solve_ivp")
    dts = [r[0] for r in rows]
    ax[1].loglog(dts, [r[1] for r in rows], "o-", label="attitude error")
    ax[1].loglog(dts, [r[2] for r in rows], "s-", label="norm drift")
    ax[1].set_xlabel("RK4 step (s)"); ax[1].legend(); ax[1].set_title("RK4 on q, 30 s")
    ax[2].loglog([r[0] for r in interp_rows], [r[1] for r in interp_rows], "o-", label="SLERP")
    ax[2].loglog([r[0] for r in interp_rows], [r[2] for r in interp_rows], "s-", label="NLERP")
    ax[2].set_xlabel("star-tracker sample interval (s)"); ax[2].set_ylabel("max error (rad)"); ax[2].legend()
    ax[2].set_title("interpolating sampled attitude")
    fig.tight_layout()
    fig.savefig(out_dir / "quaternions.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out_dir / 'quaternions.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
