"""P01.15 Matrix Transformations & Rotation Matrices: reference solution.

Direction cosine matrix (DCM) tools for aircraft and spacecraft attitude, built from scratch
with NumPy and checked against ``numpy.linalg`` and ``scipy.spatial.transform.Rotation``.

Convention (aerospace standard, used throughout the curriculum):
    * Every matrix here is a *passive* rotation, also called a coordinate transformation.
      ``C_bn = dcm_from_euler_321(psi, theta, phi)`` maps the components of one fixed vector
      from the navigation (NED) frame to the body frame: ``v_b = C_bn @ v_n``.
    * The 3-2-1 sequence rotates by yaw ``psi`` about z, then pitch ``theta`` about the new y,
      then roll ``phi`` about the newest x:  C_bn = R1(phi) @ R2(theta) @ R3(psi).
    * Angles are in radians. Body axes: x forward, y right wing, z down. NED: x north,
      y east, z down.

Run:  python solution.py          (prints the attitude example, writes dcm_integration_error.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-01/lesson-15-rotation-matrices
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

# Default tolerance used when deciding whether a matrix "is" a rotation matrix.
ROTATION_TOL = 1e-9


# ---------------------------------------------------------------- elementary rotations

def rot_x(angle: float) -> np.ndarray:
    """Return R1(angle): the passive rotation about the x axis.

    Transforms vector components from a frame into one rotated by ``angle`` (radians, right-hand rule)
    about the shared x axis. Rows are the new basis vectors written in the old frame.
    """
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[1.0, 0.0, 0.0],
                     [0.0, c, s],
                     [0.0, -s, c]])


def rot_y(angle: float) -> np.ndarray:
    """Return R2(angle): the passive rotation about the y axis (same convention as ``rot_x``)."""
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, 0.0, -s],
                     [0.0, 1.0, 0.0],
                     [s, 0.0, c]])


def rot_z(angle: float) -> np.ndarray:
    """Return R3(angle): the passive rotation about the z axis (same convention as ``rot_x``)."""
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, s, 0.0],
                     [-s, c, 0.0],
                     [0.0, 0.0, 1.0]])


# ---------------------------------------------------------------- the DCM and its properties

def dcm_from_euler_321(psi: float, theta: float, phi: float) -> np.ndarray:
    """Return C_bn = R1(phi) @ R2(theta) @ R3(psi), the NED-to-body DCM for yaw-pitch-roll angles.

    ``v_body = C_bn @ v_ned``. The transpose maps body components back to NED.
    """
    return rot_x(phi) @ rot_y(theta) @ rot_z(psi)


def is_rotation_matrix(C: np.ndarray, tol: float = ROTATION_TOL) -> bool:
    """Return True if ``C`` is 3x3, C^T C = I to within ``tol`` (max-abs entry) and det C = +1 to within ``tol``.

    An orthogonal matrix with det = -1 is a reflection, not a rotation, and must return False.
    """
    C = np.asarray(C, dtype=float)
    if C.shape != (3, 3):
        return False
    orthogonal = np.max(np.abs(C.T @ C - np.eye(3))) < tol
    proper = abs(np.linalg.det(C) - 1.0) < tol
    return bool(orthogonal and proper)


def transform(C: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Apply the DCM to one vector of shape (3,) or to a batch of shape (N, 3), row by row.

    Returns an array of the same shape as ``v``. Raise ``ValueError`` for any other shape.
    """
    C = np.asarray(C, dtype=float)
    v = np.asarray(v, dtype=float)
    if v.shape == (3,):
        return C @ v
    if v.ndim == 2 and v.shape[1] == 3:
        return v @ C.T
    raise ValueError(f"v must have shape (3,) or (N, 3), got {v.shape}")


def compose(C_cb: np.ndarray, C_ba: np.ndarray) -> np.ndarray:
    """Return C_ca = C_cb @ C_ba: the DCM from frame a to frame c given a->b and b->c.

    Order matters: the first rotation applied is the right-hand factor.
    """
    return np.asarray(C_cb, dtype=float) @ np.asarray(C_ba, dtype=float)


def euler_321_from_dcm(C: np.ndarray) -> tuple[float, float, float]:
    """Return (psi, theta, phi) in radians such that ``dcm_from_euler_321(psi, theta, phi)`` reproduces ``C``.

    theta is taken in [-pi/2, pi/2]; psi and phi in (-pi, pi]. At gimbal lock (|theta| = 90 deg) only
    the sum or difference of psi and phi is defined: return phi = 0 and put the whole rotation in psi.
    """
    C = np.asarray(C, dtype=float)
    s_theta = -C[0, 2]
    if abs(s_theta) >= 1.0 - 1e-12:            # gimbal lock
        theta = np.copysign(np.pi / 2, s_theta)
        if s_theta > 0:                           # C[1,0] = sin(phi - psi), C[1,1] = cos(phi - psi)
            psi = -np.arctan2(C[1, 0], C[1, 1])
        else:                                     # C[1,0] = -sin(phi + psi), C[1,1] = cos(phi + psi)
            psi = np.arctan2(-C[1, 0], C[1, 1])
        return float(psi), float(theta), 0.0
    theta = np.arcsin(np.clip(s_theta, -1.0, 1.0))
    psi = np.arctan2(C[0, 1], C[0, 0])
    phi = np.arctan2(C[1, 2], C[2, 2])
    return float(psi), float(theta), float(phi)


# ---------------------------------------------------------------- axis-angle form

def skew(v: np.ndarray) -> np.ndarray:
    """Return the 3x3 skew-symmetric matrix [v x] such that skew(v) @ w == np.cross(v, w)."""
    v = np.asarray(v, dtype=float)
    return np.array([[0.0, -v[2], v[1]],
                     [v[2], 0.0, -v[0]],
                     [-v[1], v[0], 0.0]])


def dcm_from_axis_angle(axis: np.ndarray, angle: float) -> np.ndarray:
    """Return the passive DCM for a frame rotated by ``angle`` about the unit vector ``axis`` (Rodrigues formula).

    C = cos(a) I + (1 - cos(a)) n n^T - sin(a) [n x]. The axis is normalised internally; raise
    ``ValueError`` if it has zero length.
    """
    n = np.asarray(axis, dtype=float)
    norm = np.linalg.norm(n)
    if norm == 0.0:
        raise ValueError("axis must have non-zero length")
    n = n / norm
    c, s = np.cos(angle), np.sin(angle)
    return c * np.eye(3) + (1.0 - c) * np.outer(n, n) - s * skew(n)


def axis_angle_from_dcm(C: np.ndarray) -> tuple[np.ndarray, float]:
    """Return (unit axis, angle in [0, pi]) of the rotation ``C`` (inverse of ``dcm_from_axis_angle``).

    The skew part gives w = [C12 - C21, C20 - C02, C01 - C10] = 2 sin(a) n and the trace gives
    trace C - 1 = 2 cos(a), so angle = atan2(|w|, trace C - 1), which keeps full precision near 0
    and near pi (arccos of the trace alone would not). For 0 < angle < pi the axis is w / |w|;
    at angle = 0 return the x axis; near angle = pi use the symmetric part, where the sign of the
    axis is arbitrary (both signs give the same C).
    """
    C = np.asarray(C, dtype=float)
    w = np.array([C[1, 2] - C[2, 1], C[2, 0] - C[0, 2], C[0, 1] - C[1, 0]])
    two_sin = np.linalg.norm(w)
    two_cos = np.trace(C) - 1.0
    angle = float(np.arctan2(two_sin, two_cos))
    if angle < 1e-12:
        return np.array([1.0, 0.0, 0.0]), 0.0
    if two_sin > 2e-6:
        return w / two_sin, angle
    # angle ~ pi: n n^T = (sym(C) - cos(a) I) / (1 - cos(a)); take the column with the largest diagonal
    cos_a = 0.5 * two_cos
    nnT = (0.5 * (C + C.T) - cos_a * np.eye(3)) / (1.0 - cos_a)
    k = int(np.argmax(np.diag(nnT)))
    n = nnT[:, k] / np.sqrt(max(nnT[k, k], 1e-300))
    return n / np.linalg.norm(n), angle


def small_angle_dcm(dtheta: np.ndarray) -> np.ndarray:
    """Return the first-order DCM I - [dtheta x] for a small rotation vector ``dtheta`` (radians).

    This is what a strapdown navigator uses for one gyro sample. It is not exactly orthogonal:
    C^T C = I + [dtheta x]^T [dtheta x], an error of order |dtheta|**2.
    """
    return np.eye(3) - skew(dtheta)


def orthonormalize(C: np.ndarray, iterations: int = 1) -> np.ndarray:
    """Return C pulled back onto the rotation group by ``iterations`` steps of C <- C (3I - C^T C) / 2.

    One step removes the first-order orthogonality error; each further step squares the remaining error.
    """
    C = np.asarray(C, dtype=float).copy()
    for _ in range(iterations):
        C = 0.5 * C @ (3.0 * np.eye(3) - C.T @ C)
    return C


def orthogonality_error(C: np.ndarray) -> float:
    """Return max |C^T C - I|, the amount by which ``C`` fails to be orthogonal."""
    C = np.asarray(C, dtype=float)
    return float(np.max(np.abs(C.T @ C - np.eye(3))))


# ---------------------------------------------------------------- applications

def gravity_in_body(theta: float, phi: float, g: float = 9.80665) -> np.ndarray:
    """Return the gravity vector (NED: [0, 0, g]) expressed in body axes for pitch ``theta`` and roll ``phi``.

    Closed form: g * [-sin(theta), sin(phi) cos(theta), cos(phi) cos(theta)]. Yaw does not enter
    because gravity is along the yaw axis.
    """
    return dcm_from_euler_321(0.0, theta, phi) @ np.array([0.0, 0.0, g])


def body_velocity_to_ned(psi: float, theta: float, phi: float, v_body: np.ndarray) -> np.ndarray:
    """Return the velocity in NED components given its body components and the 3-2-1 attitude."""
    return dcm_from_euler_321(psi, theta, phi).T @ np.asarray(v_body, dtype=float)


def flight_path_angle(v_ned: np.ndarray) -> float:
    """Return the flight-path angle gamma = asin(-v_down / |v|) in radians (positive = climbing)."""
    v = np.asarray(v_ned, dtype=float)
    speed = np.linalg.norm(v)
    if speed == 0.0:
        raise ValueError("flight-path angle is undefined at zero speed")
    return float(np.arcsin(np.clip(-v[2] / speed, -1.0, 1.0)))


def integrate_constant_rate(axis: np.ndarray, rate: float, duration: float, n_steps: int,
                            renormalize: bool = False) -> np.ndarray:
    """Compose ``n_steps`` small-angle DCMs for a frame turning at constant ``rate`` (rad/s) about ``axis``.

    Each step uses ``small_angle_dcm(axis * rate * dt)``; with ``renormalize`` the running product
    is passed through one ``orthonormalize`` step after every multiplication. The exact answer is
    ``dcm_from_axis_angle(axis, rate * duration)``.
    """
    n = np.asarray(axis, dtype=float)
    n = n / np.linalg.norm(n)
    if n_steps < 1:
        raise ValueError("n_steps must be at least 1")
    dt = duration / n_steps
    step = small_angle_dcm(n * rate * dt)
    C = np.eye(3)
    for _ in range(n_steps):
        C = step @ C
        if renormalize:
            C = orthonormalize(C)
    return C


def misalignment_angle(C_a: np.ndarray, C_b: np.ndarray) -> float:
    """Return the angle (rad) of the rotation C_a @ C_b^T: how far the frame described by ``C_a`` is from ``C_b``.

    Zero when the two DCMs agree; insensitive to how many full turns either rotation contains.
    """
    _, ang = axis_angle_from_dcm(np.asarray(C_a, dtype=float) @ np.asarray(C_b, dtype=float).T)
    return ang


def integration_errors(axis: np.ndarray, rate: float, duration: float, n_steps: int,
                       renormalize: bool = False) -> tuple[float, float]:
    """Return (orthogonality error, misalignment angle in rad) of ``integrate_constant_rate`` against the exact DCM.

    The misalignment is ``misalignment_angle`` between the numerical DCM (pulled back onto the rotation
    group by five ``orthonormalize`` steps) and ``dcm_from_axis_angle(axis, rate * duration)``.
    """
    C_num = integrate_constant_rate(axis, rate, duration, n_steps, renormalize)
    orth = orthogonality_error(C_num)
    C_exact = dcm_from_axis_angle(axis, rate * duration)
    return orth, misalignment_angle(orthonormalize(C_num, iterations=5), C_exact)


def main(out_dir: str | Path = ".") -> None:
    """Print the attitude example and plot the small-angle integration errors against step size."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    deg = np.pi / 180.0
    psi, theta, phi = 30.0 * deg, 10.0 * deg, -5.0 * deg
    C = dcm_from_euler_321(psi, theta, phi)
    np.set_printoptions(precision=5, suppress=True)
    print("attitude: yaw 30 deg, pitch 10 deg, roll -5 deg")
    print("C_bn =\n", C)
    print(f"orthogonality error |C^T C - I| = {orthogonality_error(C):.1e}, det = {np.linalg.det(C):.12f}")
    print("recovered angles [deg]:", np.array(euler_321_from_dcm(C)) / deg)
    n, a = axis_angle_from_dcm(C)
    print(f"single-axis form: axis = {n}, angle = {a / deg:.4f} deg")

    print("\ngravity in body axes [m/s^2]:", gravity_in_body(theta, phi))
    alpha = 4.0 * deg
    V = 70.0
    v_b = V * np.array([np.cos(alpha), 0.0, np.sin(alpha)])
    v_n = body_velocity_to_ned(psi, theta, phi, v_b)
    print(f"airspeed {V} m/s at alpha = 4 deg -> v_body = {v_b}")
    print(f"v_ned = {v_n}  (north, east, down)")
    print(f"flight-path angle = {flight_path_angle(v_n) / deg:.4f} deg  (theta - alpha = {(theta - alpha) / deg:.1f} deg)")

    # Computational experiment: strapdown integration at 30 deg/s for 20 s (600 deg) about a skewed axis.
    axis = np.array([1.0, 2.0, 2.0]) / 3.0
    rate, duration = 30.0 * deg, 20.0
    steps = np.array([200, 400, 800, 1600, 3200, 6400])
    rows = []
    for n_steps in steps:
        orth_raw, ang_raw = integration_errors(axis, rate, duration, int(n_steps))
        orth_ren, ang_ren = integration_errors(axis, rate, duration, int(n_steps), renormalize=True)
        rows.append((duration / n_steps, orth_raw, ang_raw, orth_ren, ang_ren))
    print("\n   dt [s]   orth err (raw)   misalign (raw) [rad]   orth err (renorm)   misalign (renorm) [rad]")
    for dt, o1, a1, o2, a2 in rows:
        print(f"  {dt:7.4f}   {o1:12.3e}   {a1:18.3e}   {o2:15.3e}   {a2:20.3e}")
    rows = np.array(rows)
    print("observed convergence order of the misalignment:", np.round(np.log2(rows[:-1, 2] / rows[1:, 2]), 2))
    print("observed convergence order of the raw orthogonality error:", np.round(np.log2(rows[:-1, 1] / rows[1:, 1]), 2))
    print(f"predicted misalignment at dt = {rows[0, 0]} s: rate*T*(rate*dt)^2/3 = "
          f"{rate * duration * (rate * rows[0, 0]) ** 2 / 3:.3e} rad")

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.loglog(rows[:, 0], rows[:, 1], "o-", label="orthogonality error, raw")
    ax.loglog(rows[:, 0], rows[:, 2], "s-", label="misalignment angle, raw")
    ax.loglog(rows[:, 0], rows[:, 3], "o--", label="orthogonality error, renormalized")
    ax.loglog(rows[:, 0], rows[:, 4], "s--", label="misalignment angle, renormalized")
    ax.set_xlabel("step size dt [s]")
    ax.set_ylabel("error after 20 s")
    ax.set_title("Small-angle DCM integration, 30 deg/s")
    ax.grid(True, which="both", lw=0.4)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "dcm_integration_error.png", dpi=150)
    print(f"saved {out_dir / 'dcm_integration_error.png'}")


if __name__ == "__main__":
    main()
