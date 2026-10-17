"""P06.07 3D Rotations: Euler Angles & DCMs: reference solution.

An attitude kinematics library: the direction cosine matrix (DCM) of any of the twelve Euler-angle
sequences, the angles back from a DCM (with the gimbal-lock case handled), the matrix that turns
Euler-angle rates into body rates and its inverse, the DCM kinematic equation dC/dt = -[w x] C, its
integration with `scipy.integrate.solve_ivp` and with a fixed-step RK4 that re-orthonormalises, and
two closed-form attitude histories to check the integration against: a constant rate about a fixed
axis, and steady coning (a 3-1-3 motion with constant rates, the free motion of a symmetric spinner).

Conventions (those of P01.15): every matrix is passive. C = C_b/n maps the components of a vector in the
reference frame n (north-east-down for an aircraft, inertial for a spacecraft) to its components in the
body frame b: v_b = C v_n. The elementary rotations are R1, R2, R3 about axes 1, 2, 3 (x, y, z). An Euler
sequence is a string of three axis digits in the order the turns are made, e.g. "321" = yaw psi about 3,
then pitch theta about the new 2, then roll phi about the newest 1, and

    C = R_{s3}(a3) R_{s2}(a2) R_{s1}(a1)          with angles (a1, a2, a3) in the order of the turns.

w is the angular velocity of b relative to n, in body components (what a strapdown gyro triad measures).
Angles are in radians.

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes attitude_kinematics.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-06/lesson-07-euler-angles-dcms
"""
from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

SEQUENCES = ["121", "123", "131", "132", "212", "213", "231", "232", "312", "313", "321", "323"]
RateFn = Callable[[float], np.ndarray]


# ---------------------------------------------------------------- elementary rotations (given, from P01.15)

def rot(axis: int, angle: float) -> np.ndarray:
    """Passive elementary rotation about body axis 1, 2 or 3 by `angle` (rad): the components of a fixed
    vector in a frame turned by `angle` about that axis. R3(psi) = [[c, s, 0], [-s, c, 0], [0, 0, 1]]."""
    c, s = math.cos(angle), math.sin(angle)
    if axis == 1:
        return np.array([[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]])
    if axis == 2:
        return np.array([[c, 0.0, -s], [0.0, 1.0, 0.0], [s, 0.0, c]])
    if axis == 3:
        return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])
    raise ValueError("axis must be 1, 2 or 3")


def skew(v: np.ndarray) -> np.ndarray:
    """The cross-product matrix [v x]: skew(v) @ u = v x u."""
    return np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])


def dcm_from_axis_angle(axis: np.ndarray, angle: float) -> np.ndarray:
    """Passive Rodrigues formula: C = cos(a) E + (1 - cos(a)) n n^T - sin(a) [n x] (P01.15)."""
    n = np.asarray(axis, dtype=float)
    n = n / np.linalg.norm(n)
    c, s = math.cos(angle), math.sin(angle)
    return c * np.eye(3) + (1.0 - c) * np.outer(n, n) - s * skew(n)


def _parse(seq: str) -> tuple[int, int, int]:
    """Validate an Euler sequence string and return its three axis digits."""
    if len(seq) != 3 or any(ch not in "123" for ch in seq) or seq[0] == seq[1] or seq[1] == seq[2]:
        raise ValueError(f"not an Euler sequence: {seq!r} (three axes from 1, 2, 3, no two neighbours equal)")
    return int(seq[0]), int(seq[1]), int(seq[2])


def is_symmetric(seq: str) -> bool:
    """True for the six 'proper' (symmetric) sequences such as 3-1-3, whose first and last axes coincide."""
    a, _, c = _parse(seq)
    return a == c


# ---------------------------------------------------------------- Euler angles <-> DCM

def dcm_from_euler(angles: np.ndarray, seq: str = "321") -> np.ndarray:
    """C_b/n = R_{s3}(a3) R_{s2}(a2) R_{s1}(a1) for the sequence `seq` and angles (a1, a2, a3)."""
    i, j, k = _parse(seq)
    a1, a2, a3 = angles
    return rot(k, a3) @ rot(j, a2) @ rot(i, a1)


def euler_from_dcm(C: np.ndarray, seq: str = "321", tol: float = 1e-10) -> np.ndarray:
    """Angles (a1, a2, a3) of `seq` that reproduce C.

    Ranges: asymmetric sequences a2 in [-pi/2, pi/2]; symmetric sequences a2 in [0, pi]; a1, a3 in (-pi, pi].
    Method: M = C^T = A_i(a1) A_j(a2) A_k(a3) with A the active rotations, eps = +1 if (i, j, k) is a cyclic
    order of (1, 2, 3) and -1 otherwise (for a symmetric sequence k is the axis not used).
      asymmetric: sin a2 = eps M[i,k], a1 = atan2(-eps M[j,k], M[k,k]), a3 = atan2(-eps M[i,j], M[i,i]);
      symmetric:  cos a2 = M[i,i],     a1 = atan2(M[j,i], -eps M[k,i]), a3 = atan2(M[i,j], eps M[i,k]).
    At gimbal lock (cos a2 = 0, or sin a2 = 0 for a symmetric sequence) only a combination of a1 and a3
    is defined; the function then sets a3 = 0 and puts the whole rotation into a1.
    """
    s1, s2, s3 = _parse(seq)
    i, j = s1 - 1, s2 - 1
    sym = s1 == s3
    k = 3 - i - j if sym else s3 - 1
    eps = 1.0 if (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] else -1.0
    M = np.asarray(C, dtype=float).T
    if sym:
        a2 = math.acos(max(-1.0, min(1.0, M[i, i])))
        locked = abs(math.sin(a2)) < tol
        if not locked:
            a1 = math.atan2(M[j, i], -eps * M[k, i])
            a3 = math.atan2(M[i, j], eps * M[i, k])
    else:
        a2 = math.asin(max(-1.0, min(1.0, eps * M[i, k])))
        locked = abs(math.cos(a2)) < tol
        if not locked:
            a1 = math.atan2(-eps * M[j, k], M[k, k])
            a3 = math.atan2(-eps * M[i, j], M[i, i])
    if locked:
        a3 = 0.0
        Ni = M @ rot(s2, a2)                        # = A_i(a1), since A_j(a2)^T = R_j(a2)
        p, q = (i + 1) % 3, (i + 2) % 3
        a1 = math.atan2(Ni[q, p], Ni[p, p])
    return np.array([a1, a2, a3])


# ---------------------------------------------------------------- Euler-angle kinematics

def euler_rate_matrix(angles: np.ndarray, seq: str = "321") -> np.ndarray:
    """The 3 x 3 matrix B with w_b = B @ (da1/dt, da2/dt, da3/dt).

    Each angle rate is a rotation rate about an intermediate axis; carried into body axes:
    w_b = R_{s3}(a3) R_{s2}(a2) e_{s1} da1 + R_{s3}(a3) e_{s2} da2 + e_{s3} da3.
    For 3-2-1: p = dphi - dpsi sin(theta), q = dtheta cos(phi) + dpsi sin(phi) cos(theta),
    r = -dtheta sin(phi) + dpsi cos(phi) cos(theta)."""
    i, j, k = _parse(seq)
    a1, a2, a3 = angles
    E = np.eye(3)
    R3 = rot(k, a3)
    return np.column_stack([R3 @ rot(j, a2) @ E[i - 1], R3 @ E[j - 1], E[k - 1]])


def euler_angle_rates(angles: np.ndarray, w: np.ndarray, seq: str = "321", tol: float = 1e-9) -> np.ndarray:
    """Euler-angle rates from body rates: solve B(angles) x = w. Raises ValueError at gimbal lock
    (|det B| < tol; |det B| = |cos a2| for an asymmetric sequence and |sin a2| for a symmetric one)."""
    B = euler_rate_matrix(angles, seq)
    if abs(np.linalg.det(B)) < tol:
        raise ValueError("gimbal lock: the Euler-angle rates are undefined at this attitude")
    return np.linalg.solve(B, np.asarray(w, dtype=float))


# ---------------------------------------------------------------- the DCM kinematic equation

def dcm_rate(C: np.ndarray, w: np.ndarray) -> np.ndarray:
    """dC_b/n/dt = -[w x] C_b/n with w the body rate of b relative to n in body axes."""
    return -skew(w) @ C


def orthogonality_error(C: np.ndarray) -> float:
    """Frobenius norm of C^T C - E: zero for a rotation matrix."""
    return float(np.linalg.norm(C.T @ C - np.eye(3)))


def orthonormalize_svd(C: np.ndarray) -> np.ndarray:
    """The rotation matrix nearest to C in the Frobenius norm: C = U S V^T  ->  U diag(1, 1, det(U V^T)) V^T.
    This is the polar factor of C (the exact version of P01.15's first-order step C (3E - C^T C)/2)."""
    U, _, Vt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(U @ Vt))
    return U @ np.diag([1.0, 1.0, d]) @ Vt


def misalignment_angle(Ca: np.ndarray, Cb: np.ndarray) -> float:
    """The single rotation angle (rad) that takes attitude Cb to attitude Ca: the angle of E = Ca Cb^T,
    atan2(|w|, tr E - 1) with w = (E32 - E23, E13 - E31, E21 - E12), accurate near 0 and near pi."""
    E = Ca @ Cb.T
    w = np.array([E[2, 1] - E[1, 2], E[0, 2] - E[2, 0], E[1, 0] - E[0, 1]])
    return math.atan2(float(np.linalg.norm(w)), float(np.trace(E) - 1.0))


def integrate_dcm(C0: np.ndarray, w_fn: RateFn, t_end: float, n_samples: int = 201,
                  rtol: float = 1e-10, atol: float = 1e-12) -> tuple[np.ndarray, np.ndarray]:
    """Integrate dC/dt = -[w(t) x] C as nine coupled ODEs with solve_ivp (DOP853).
    Returns (t, C) with C of shape (n_samples, 3, 3). No re-orthonormalisation: the drift of C^T C from E
    measures the integrator's error."""
    def rhs(t, y):
        return dcm_rate(y.reshape(3, 3), w_fn(t)).ravel()
    t_eval = np.linspace(0.0, t_end, n_samples)
    sol = solve_ivp(rhs, (0.0, t_end), np.asarray(C0, dtype=float).ravel(), method="DOP853", t_eval=t_eval,
                    rtol=rtol, atol=atol)
    return sol.t, sol.y.T.reshape(-1, 3, 3)


def integrate_dcm_rk4(C0: np.ndarray, w_fn: RateFn, t_end: float, dt: float,
                      renormalize_every: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Fixed-step classical RK4 on dC/dt = -[w x] C (what a flight computer does at its gyro rate).
    If renormalize_every = k > 0, the matrix is replaced by `orthonormalize_svd` every k steps.
    Returns (t, C) at every step."""
    n = int(round(t_end / dt))
    C = np.asarray(C0, dtype=float).copy()
    ts = np.linspace(0.0, n * dt, n + 1)
    out = np.empty((n + 1, 3, 3))
    out[0] = C
    for m in range(n):
        t = ts[m]
        k1 = dcm_rate(C, w_fn(t))
        k2 = dcm_rate(C + 0.5 * dt * k1, w_fn(t + 0.5 * dt))
        k3 = dcm_rate(C + 0.5 * dt * k2, w_fn(t + 0.5 * dt))
        k4 = dcm_rate(C + dt * k3, w_fn(t + dt))
        C = C + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        if renormalize_every and (m + 1) % renormalize_every == 0:
            C = orthonormalize_svd(C)
        out[m + 1] = C
    return ts, out


def integrate_euler_angles(angles0: np.ndarray, w_fn: RateFn, t_end: float, seq: str = "321",
                           n_samples: int = 201, rtol: float = 1e-10, atol: float = 1e-12) -> dict:
    """Integrate the Euler-angle rates B(angles)^-1 w(t) with solve_ivp (DOP853), the alternative a DCM
    replaces. Near gimbal lock the rates grow without bound; the integration may stop early or lose accuracy.
    Returns {"t", "angles" (n, 3), "success", "nfev", "message"}."""
    def rhs(t, y):
        B = euler_rate_matrix(y, seq)
        return np.linalg.solve(B, w_fn(t))
    t_eval = np.linspace(0.0, t_end, n_samples)
    sol = solve_ivp(rhs, (0.0, t_end), np.asarray(angles0, dtype=float), method="DOP853", t_eval=t_eval,
                    rtol=rtol, atol=atol)
    return {"t": sol.t, "angles": sol.y.T, "success": bool(sol.success), "nfev": int(sol.nfev), "message": sol.message}


def integrate_euler_angles_rk4(angles0: np.ndarray, w_fn: RateFn, t_end: float, dt: float,
                               seq: str = "321") -> tuple[np.ndarray, np.ndarray]:
    """Fixed-step classical RK4 on the Euler-angle rates, the angle-based counterpart of
    `integrate_dcm_rk4`. Returns (t, angles) at every step; raises ValueError if a stage lands exactly on
    gimbal lock."""
    n = int(round(t_end / dt))
    y = np.asarray(angles0, dtype=float).copy()
    ts = np.linspace(0.0, n * dt, n + 1)
    out = np.empty((n + 1, 3))
    out[0] = y
    f = lambda t, y: euler_angle_rates(y, w_fn(t), seq, tol=0.0)
    for m in range(n):
        t = ts[m]
        k1 = f(t, y)
        k2 = f(t + 0.5 * dt, y + 0.5 * dt * k1)
        k3 = f(t + 0.5 * dt, y + 0.5 * dt * k2)
        k4 = f(t + dt, y + dt * k3)
        y = y + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        out[m + 1] = y
    return ts, out


# ---------------------------------------------------------------- closed-form attitude histories

def single_axis_dcm(C0: np.ndarray, axis: np.ndarray, rate: float, t: float) -> np.ndarray:
    """Exact attitude after a constant body rate `rate` about a fixed unit `axis` for time t, from C0.
    The axis is fixed in both frames, so C(t) = C(axis, rate t) C0 (passive Rodrigues)."""
    return dcm_from_axis_angle(axis, rate * t) @ C0


def coning_dcm(t: float, psi_dot: float, phi_dot: float, beta: float, psi0: float = 0.0, phi0: float = 0.0) -> np.ndarray:
    """Exact attitude of steady coning: a 3-1-3 motion with constant precession rate psi_dot, constant
    nutation (half-cone) angle beta and constant spin rate phi_dot:
    C(t) = R3(phi0 + phi_dot t) R1(beta) R3(psi0 + psi_dot t)."""
    return dcm_from_euler([psi0 + psi_dot * t, beta, phi0 + phi_dot * t], "313")


def coning_body_rate(t: float, psi_dot: float, phi_dot: float, beta: float, phi0: float = 0.0) -> np.ndarray:
    """Body rates of steady coning, from B_313 (psi_dot, 0, phi_dot):
    w = (psi_dot sin(beta) sin(phi), psi_dot sin(beta) cos(phi), phi_dot + psi_dot cos(beta)), phi = phi0 + phi_dot t."""
    ph = phi0 + phi_dot * t
    return np.array([psi_dot * math.sin(beta) * math.sin(ph), psi_dot * math.sin(beta) * math.cos(ph),
                     phi_dot + psi_dot * math.cos(beta)])


def hammerhead(q: float, r: float, theta_top: float) -> tuple[list[tuple[np.ndarray, float]], np.ndarray]:
    """The lesson's stall turn (hammerhead): from level flight, pitch up at body rate q until the pitch is
    theta_top, then yaw at body rate r through 180 deg (the nose swings across the vertical), then pitch at q
    for the same time again (the nose comes down to theta_top below the horizon, heading reversed).
    Returns (segments, C_end): a list of (constant body rate, duration) and the exact final attitude
    from three single-axis rotations."""
    t1, t2 = theta_top / q, math.pi / abs(r)
    segs = [(np.array([0.0, q, 0.0]), t1), (np.array([0.0, 0.0, r]), t2), (np.array([0.0, q, 0.0]), t1)]
    C = np.eye(3)
    for w, T in segs:
        C = single_axis_dcm(C, w / np.linalg.norm(w), float(np.linalg.norm(w)), T)
    return segs, C


def piecewise_rate(segments: list[tuple[np.ndarray, float]]) -> tuple[RateFn, float]:
    """A body-rate function w(t) from a list of (constant rate, duration), and the total duration."""
    ends = np.cumsum([T for _, T in segments])

    def w_fn(t: float) -> np.ndarray:
        m = min(int(np.searchsorted(ends, t, side="right")), len(segments) - 1)
        return segments[m][0]
    return w_fn, float(ends[-1])


def integrate_segments_rk4(x0: np.ndarray, segments: list[tuple[np.ndarray, float]], dt: float,
                           kind: str = "dcm", seq: str = "321") -> np.ndarray:
    """RK4 through a piecewise-constant rate history, restarting at every switch so that no step straddles
    one (each segment uses the step closest to dt that divides it). kind = "dcm" carries C, kind = "euler"
    carries the Euler angles of `seq`. Returns the final C (3 x 3) or the final angles (3,)."""
    x = np.asarray(x0, dtype=float).copy()
    for w, T in segments:
        n = max(1, int(math.ceil(T / dt - 1e-12)))
        h = T / n
        wf = lambda t, w=w: w
        if kind == "dcm":
            _, out = integrate_dcm_rk4(x, wf, n * h, h)
        else:
            _, out = integrate_euler_angles_rk4(x, wf, n * h, h, seq)
        x = out[-1]
    return x


# ---------------------------------------------------------------- the lesson's output

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save attitude_kinematics.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    deg = math.pi / 180.0
    ang = np.array([30.0, 10.0, -5.0]) * deg
    C = dcm_from_euler(ang, "321")
    print("3-2-1 DCM for yaw 30, pitch 10, roll -5 deg:")
    print(np.array2string(C, precision=5, suppress_small=True))
    print(f"  orthogonality error {orthogonality_error(C):.1e}, det {np.linalg.det(C):.15f}")
    for seq in ["313", "123", "213"]:
        a = euler_from_dcm(C, seq)
        print(f"  same attitude as {seq}: {np.round(a / deg, 4)} deg; rebuilt error {np.abs(dcm_from_euler(a, seq) - C).max():.1e}")
    B = euler_rate_matrix(ang, "321")
    print(f"3-2-1 rate matrix B at that attitude, det B = {np.linalg.det(B):.6f} (cos theta = {math.cos(ang[1]):.6f})")
    w = np.array([0.1, 0.05, -0.02])
    print(f"  body rate {w} rad/s -> Euler rates (psi, theta, phi) dot = {np.round(euler_angle_rates(ang, w), 6)} rad/s")
    for th in [80.0, 89.0, 89.9]:
        r = euler_angle_rates([0.0, th * deg, 30.0 * deg], np.array([0.0, 0.0, 0.1]))
        print(f"  pitch {th:5.1f} deg, roll 30 deg, yaw rate r = 0.1 rad/s: psi dot = {r[0]:9.4f}, phi dot = {r[2]:9.4f} rad/s")

    # DCM integration checks
    axis = np.array([1.0, 2.0, 2.0]) / 3.0
    rate = 30.0 * deg
    t, Cs = integrate_dcm(np.eye(3), lambda t: rate * axis, 20.0, n_samples=401)
    err = max(misalignment_angle(Cs[m], single_axis_dcm(np.eye(3), axis, rate, t[m])) for m in range(t.size))
    print(f"single axis 30 deg/s for 20 s (solve_ivp, rtol 1e-10): max misalignment {err:.2e} rad, "
          f"max orthogonality error {max(orthogonality_error(c) for c in Cs):.2e}")
    pd, sd, beta = 0.4, 2.0, 20.0 * deg
    wc = lambda t: coning_body_rate(t, pd, sd, beta)
    C0 = coning_dcm(0.0, pd, sd, beta)
    t, Cs = integrate_dcm(C0, wc, 30.0, n_samples=601)
    err = max(misalignment_angle(Cs[m], coning_dcm(t[m], pd, sd, beta)) for m in range(t.size))
    print(f"coning psi_dot {pd}, phi_dot {sd} rad/s, beta 20 deg, 30 s (solve_ivp): max misalignment {err:.2e} rad, "
          f"orthogonality {max(orthogonality_error(c) for c in Cs):.2e}")
    print("RK4 at fixed step on the coning case, 30 s:")
    rows = []
    for dt in [0.1, 0.05, 0.025, 0.0125]:
        tr, Cr = integrate_dcm_rk4(C0, wc, 30.0, dt)
        tr2, Cr2 = integrate_dcm_rk4(C0, wc, 30.0, dt, renormalize_every=1)
        e1 = misalignment_angle(Cr[-1], coning_dcm(tr[-1], pd, sd, beta))
        e2 = misalignment_angle(orthonormalize_svd(Cr2[-1]), coning_dcm(tr2[-1], pd, sd, beta))
        rows.append((dt, e1, orthogonality_error(Cr[-1]), e2, orthogonality_error(Cr2[-1])))
        print(f"  dt {dt:7.4f}: misalignment {e1:.3e} rad, orthogonality {rows[-1][2]:.3e}; "
              f"renormalised each step: misalignment {e2:.3e}, orthogonality {rows[-1][4]:.1e}")
    orders = [math.log2(rows[m][1] / rows[m + 1][1]) for m in range(len(rows) - 1)]
    print(f"  observed order of the misalignment: {np.round(orders, 2)}")

    # the hammerhead: Euler angles against the DCM
    q, r = 0.5, 0.5
    print("hammerhead, q = r = 0.5 rad/s, pitch at the top 90 deg - delta:")
    for dl in [10.0, 1.0, 0.1, 0.01]:
        segs, Cend = hammerhead(q, r, (90.0 - dl) * deg)
        wf, T = piecewise_rate(segs)
        th, Ch = integrate_dcm(np.eye(3), wf, T, n_samples=2001)
        ei = integrate_euler_angles([0.0, 0.0, 0.0], wf, T, n_samples=2001)
        e_dcm = misalignment_angle(Ch[-1], Cend)
        e_eul = misalignment_angle(dcm_from_euler(ei["angles"][-1], "321"), Cend) if ei["success"] else float("nan")
        psid = max(abs(euler_angle_rates(euler_from_dcm(c), wf(t), tol=0.0)[0]) for t, c in zip(th, Ch))
        e_dcm4 = misalignment_angle(integrate_segments_rk4(np.eye(3), segs, 0.01, "dcm"), Cend)
        e_eul4 = misalignment_angle(dcm_from_euler(integrate_segments_rk4(np.zeros(3), segs, 0.01, "euler"), "321"), Cend)
        print(f"  delta {dl:6.2f} deg: largest |psi dot| {psid:7.1f} rad/s; adaptive (rtol 1e-10): DCM error {e_dcm:.1e} rad, "
              f"Euler-angle error {e_eul:.1e} rad with {ei['nfev']} rhs evaluations; RK4 at 100 Hz: DCM error {e_dcm4:.1e} rad, "
              f"Euler-angle error {e_eul4:.1e} rad")
        if dl == 1.0:
            tl, Cl, eul_int = th, Ch, ei
    eul = np.array([euler_from_dcm(c, "321") for c in Cl])
    mis = np.array([misalignment_angle(dcm_from_euler(eul_int["angles"][m], "321"), Cl[m]) for m in range(eul_int["t"].size)])
    uy = np.unwrap(eul[:, 0])
    j = int(np.argmax(np.abs(np.diff(uy))))
    print(f"  delta 1 deg: max pitch {np.degrees(eul[:, 1].max()):.4f} deg; yaw from the DCM changes by "
          f"{np.degrees(uy[j + 1] - uy[j]):.2f} deg between t = {tl[j]:.4f} and {tl[j + 1]:.4f} s "
          f"while the attitude moves {np.degrees(misalignment_angle(Cl[j + 1], Cl[j])):.3f} deg; "
          f"yaw at start {np.degrees(eul[0, 0]):.1f}, at end {np.degrees(eul[-1, 0]):.1f} deg")
    ok = eul_int["t"].size

    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    ax[0].plot(tl, np.degrees(eul[:, 0]), label="yaw psi")
    ax[0].plot(tl, np.degrees(eul[:, 1]), label="pitch theta")
    ax[0].plot(tl, np.degrees(eul[:, 2]), label="roll phi")
    ax[0].set_xlabel("t (s)"); ax[0].set_ylabel("deg"); ax[0].legend(); ax[0].set_title("3-2-1 angles, hammerhead 1 deg from vertical")
    ax[1].semilogy(tl[:ok], np.maximum(mis, 1e-17))
    ax[1].set_xlabel("t (s)"); ax[1].set_ylabel("rad"); ax[1].set_title("Euler-angle integration vs DCM")
    dts = [r[0] for r in rows]
    ax[2].loglog(dts, [r[1] for r in rows], "o-", label="misalignment")
    ax[2].loglog(dts, [r[2] for r in rows], "s-", label="orthogonality error")
    ax[2].loglog(dts, [r[4] for r in rows], "^-", label="orthogonality, renormalised")
    ax[2].set_xlabel("RK4 step (s)"); ax[2].legend(); ax[2].set_title("coning, 30 s")
    fig.tight_layout()
    fig.savefig(out_dir / "attitude_kinematics.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out_dir / 'attitude_kinematics.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
