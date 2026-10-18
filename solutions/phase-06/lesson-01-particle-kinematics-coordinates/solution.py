"""P06.01 Particle Kinematics in Many Coordinates: reference solution.

A ground radar tracks an aircraft flying a climbing turn. The radar reports where the aircraft is in
cylindrical coordinates about the radar (horizontal range r, azimuth theta, height z). From those
samples, what are the aircraft's velocity and acceleration, its speed, its normal acceleration and the
radius of the turn it is flying? This module describes the same motion three ways and converts between
them:

* polar and cylindrical coordinates (r, theta, z) with the unit vectors e_r, e_theta, e_z that turn
  with the particle: v = r' e_r + r theta' e_theta + z' e_z and
  a = (r'' - r theta'^2) e_r + (r theta'' + 2 r' theta') e_theta + z'' e_z;
* Cartesian coordinates (x, y, z) in a fixed frame;
* path (normal-tangential) coordinates: the speed, the tangential acceleration dv/dt along the
  velocity, the normal acceleration v^2 / rho toward the centre of curvature, and the curvature
  kappa = |v x a| / |v|^3.

The independent check is the analytic trajectory itself: a helix whose position, velocity and
acceleration are known in closed form. Finite differences of the sampled radar data must converge to it
at second order; the cylindrical and Cartesian routes must agree; and the path description must satisfy
a_t^2 + a_n^2 = |a|^2 in any frame.

Frames: x east, y north, z up, origin at the radar; theta measured from x toward y (counterclockwise
seen from above). Units: SI (m, s, rad).

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes radar_track.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-06/lesson-01-particle-kinematics-coordinates
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

G0 = 9.80665

# the lesson's climbing turn (illustrative numbers for a fighter-sized aircraft)
CENTRE = (8000.0, 6000.0)   # m, turn centre east and north of the radar
R_TURN = 2000.0             # m, horizontal turn radius
OMEGA_TURN = 0.1            # rad/s, turn rate (counterclockwise seen from above): 200 m/s horizontal
CLIMB = 10.0                # m/s, constant climb rate
Z0 = 3000.0                 # m, height at t = 0


# ---------------------------------------------------------------- polar and cylindrical

def polar_kinematics(r, rd, rdd, th, thd, thdd):
    """Velocity and acceleration components in the turning polar basis (e_r, e_theta):

        v_r = r',  v_theta = r theta',
        a_r = r'' - r theta'^2,  a_theta = r theta'' + 2 r' theta'.

    ``th`` is not needed for the components (only for converting them); it is accepted so that the
    signature lists the full state. Returns (v_r, v_theta, a_r, a_theta), each shaped like ``r``."""
    r, rd, rdd, thd, thdd = (np.asarray(a, dtype=float) for a in (r, rd, rdd, thd, thdd))
    return rd, r * thd, rdd - r * thd ** 2, r * thdd + 2.0 * rd * thd


def cylindrical_kinematics(r, rd, rdd, th, thd, thdd, z, zd, zdd):
    """Velocity and acceleration in cylindrical components (e_r, e_theta, e_z): the polar result of
    ``polar_kinematics`` plus z' and z'' along the fixed axis e_z. Returns (v, a), each of shape
    (..., 3) with the last axis ordered (r, theta, z)."""
    vr, vth, ar, ath = polar_kinematics(r, rd, rdd, th, thd, thdd)
    v = np.stack([vr, vth, np.asarray(zd, dtype=float) * np.ones_like(vr)], axis=-1)
    a = np.stack([ar, ath, np.asarray(zdd, dtype=float) * np.ones_like(ar)], axis=-1)
    return v, a


def cylindrical_to_cartesian(comp, th):
    """Convert a vector given by cylindrical components (c_r, c_theta, c_z), shape (..., 3), at azimuth
    ``th`` to Cartesian components: e_r = (cos th, sin th, 0), e_theta = (-sin th, cos th, 0)."""
    comp = np.asarray(comp, dtype=float)
    c, s = np.cos(th), np.sin(th)
    return np.stack([comp[..., 0] * c - comp[..., 1] * s,
                     comp[..., 0] * s + comp[..., 1] * c,
                     comp[..., 2]], axis=-1)


def cartesian_to_cylindrical(pos, vel, acc):
    """The inverse conversion: from Cartesian position, velocity and acceleration (each (..., 3)) to the
    cylindrical coordinates and their first two time derivatives,
    (r, r', r'', theta, theta', theta'', z, z', z'').

    r = sqrt(x^2 + y^2), theta = atan2(y, x); the rates follow from projecting v and a on e_r and
    e_theta and inverting ``polar_kinematics``: r' = v . e_r, theta' = (v . e_theta) / r,
    r'' = a . e_r + r theta'^2, theta'' = (a . e_theta - 2 r' theta') / r. Undefined on the axis r = 0."""
    pos, vel, acc = (np.asarray(a, dtype=float) for a in (pos, vel, acc))
    x, y = pos[..., 0], pos[..., 1]
    r = np.hypot(x, y)
    th = np.arctan2(y, x)
    er = np.stack([np.cos(th), np.sin(th)], axis=-1)
    et = np.stack([-np.sin(th), np.cos(th)], axis=-1)
    v_r = np.sum(vel[..., :2] * er, axis=-1)
    v_t = np.sum(vel[..., :2] * et, axis=-1)
    a_r = np.sum(acc[..., :2] * er, axis=-1)
    a_t = np.sum(acc[..., :2] * et, axis=-1)
    rd, thd = v_r, v_t / r
    rdd = a_r + r * thd ** 2
    thdd = (a_t - 2.0 * rd * thd) / r
    return r, rd, rdd, th, thd, thdd, pos[..., 2], vel[..., 2], acc[..., 2]


# ---------------------------------------------------------------- path coordinates

def path_kinematics(vel, acc):
    """Path (normal-tangential) description of the motion from Cartesian v and a, each (..., 3).

    Returns a dict with
      speed  |v|;
      a_t    tangential acceleration v . a / |v| (the rate of change of speed);
      a_n    normal acceleration |v x a| / |v| (>= 0, toward the centre of curvature);
      kappa  curvature |v x a| / |v|^3 (1/m) and rho = 1 / kappa (m, inf on a straight line);
      e_t, e_n, e_b  unit tangent, principal normal and binormal (e_b = e_t x e_n). e_n is undefined on
                     a straight line; it is returned as NaN there."""
    vel, acc = np.asarray(vel, dtype=float), np.asarray(acc, dtype=float)
    speed = np.linalg.norm(vel, axis=-1)
    e_t = vel / speed[..., None]
    a_t = np.sum(acc * e_t, axis=-1)
    a_perp = acc - a_t[..., None] * e_t
    a_n = np.linalg.norm(np.cross(vel, acc), axis=-1) / speed
    with np.errstate(divide="ignore", invalid="ignore"):
        kappa = a_n / speed ** 2
        rho = np.where(kappa > 0, 1.0 / kappa, np.inf)
        e_n = a_perp / np.linalg.norm(a_perp, axis=-1)[..., None]
    e_b = np.cross(e_t, e_n)
    return {"speed": speed, "a_t": a_t, "a_n": a_n, "kappa": kappa, "rho": rho,
            "e_t": e_t, "e_n": e_n, "e_b": e_b}


# ---------------------------------------------------------------- the analytic trajectory

def helix_turn(t, centre=CENTRE, R=R_TURN, omega=OMEGA_TURN, climb=CLIMB, z0=Z0, phase=0.0):
    """Closed-form climbing turn: a helix about the vertical axis through ``centre``,

        x = cx + R cos(omega t + phase),  y = cy + R sin(omega t + phase),  z = z0 + climb t.

    Returns (pos, vel, acc), each of shape (len(t), 3), with the derivatives differentiated by hand:
    the horizontal acceleration is R omega^2 toward the turn axis and there is no vertical acceleration."""
    t = np.atleast_1d(np.asarray(t, dtype=float))
    ph = omega * t + phase
    c, s = np.cos(ph), np.sin(ph)
    one = np.ones_like(t)
    pos = np.stack([centre[0] + R * c, centre[1] + R * s, z0 + climb * t], axis=-1)
    vel = np.stack([-R * omega * s, R * omega * c, climb * one], axis=-1)
    acc = np.stack([-R * omega ** 2 * c, -R * omega ** 2 * s, 0.0 * one], axis=-1)
    return pos, vel, acc


def helix_curvature(R: float, omega: float, climb: float) -> float:
    """Curvature of a helix of radius R and pitch 2 pi b, b = climb / omega: kappa = R / (R^2 + b^2)."""
    b = climb / omega
    return R / (R ** 2 + b ** 2)


# ---------------------------------------------------------------- samples and differences

def central_differences(t, f, k: int = 1):
    """First and second derivatives of uniformly sampled data by central differences over a stride of
    k samples (step h = k dt):

        f'(t_i)  ~ (f_{i+k} - f_{i-k}) / (2h),     f''(t_i) ~ (f_{i+k} - 2 f_i + f_{i-k}) / h^2,

    both second-order accurate. ``f`` may have trailing dimensions. Only the points that have both
    neighbours are returned: (t[k:-k], f'[k:-k], f''[k:-k])."""
    t = np.asarray(t, dtype=float)
    f = np.asarray(f, dtype=float)
    if k < 1 or 2 * k >= len(t):
        raise ValueError("stride must satisfy 1 <= k < len(t) / 2")
    h = k * (t[1] - t[0])
    fp, f0, fm = f[2 * k:], f[k:-k], f[:-2 * k]
    return t[k:-k], (fp - fm) / (2.0 * h), (fp - 2.0 * f0 + fm) / h ** 2


def radar_samples(pos, sigma_r: float = 0.0, sigma_th: float = 0.0, sigma_z: float = 0.0, rng=None):
    """What the radar reports: cylindrical coordinates (r, theta, z) of each Cartesian position (..., 3),
    with the azimuth unwrapped so that it is continuous. Optional SYNTHETIC Gaussian noise of standard
    deviations sigma_r (m), sigma_th (rad) and sigma_z (m) is added with ``rng``
    (a numpy Generator; required if any sigma > 0). The noise is synthetic: no real radar is modelled."""
    pos = np.asarray(pos, dtype=float)
    r = np.hypot(pos[:, 0], pos[:, 1])
    th = np.unwrap(np.arctan2(pos[:, 1], pos[:, 0]))
    z = pos[:, 2].copy()
    if sigma_r or sigma_th or sigma_z:
        r = r + sigma_r * rng.standard_normal(r.shape)
        th = th + sigma_th * rng.standard_normal(th.shape)
        z = z + sigma_z * rng.standard_normal(z.shape)
    return r, th, z


def track_from_radar(t, r, th, z, k: int = 1):
    """Velocity and acceleration from radar samples by the cylindrical route: difference r, theta and z
    (``central_differences`` with stride k), put the rates into ``cylindrical_kinematics``, and convert
    the components to Cartesian. Returns (t_mid, vel, acc) with vel and acc of shape (n, 3)."""
    tm, rd, rdd = central_differences(t, r, k)
    _, thd, thdd = central_differences(t, th, k)
    _, zd, zdd = central_differences(t, z, k)
    rm, thm = r[k:-k], th[k:-k]
    v_cyl, a_cyl = cylindrical_kinematics(rm, rd, rdd, thm, thd, thdd, z[k:-k], zd, zdd)
    return tm, cylindrical_to_cartesian(v_cyl, thm), cylindrical_to_cartesian(a_cyl, thm)


def track_cartesian(t, r, th, z, k: int = 1):
    """The same estimate by the Cartesian route: convert the samples to x, y, z first and difference
    those. Returns (t_mid, vel, acc). Agrees with ``track_from_radar`` to O(h^2), not exactly, because
    differencing does not commute with the nonlinear map (r, theta) -> (x, y)."""
    pos = np.stack([r * np.cos(th), r * np.sin(th), z], axis=-1)
    tm, vel, acc = central_differences(t, pos, k)
    return tm, vel, acc


def observed_order(errors, ratio: float = 2.0) -> np.ndarray:
    """Observed convergence orders log(e_k / e_{k+1}) / log(ratio) of errors on steps refined by ``ratio``."""
    e = np.asarray(errors, dtype=float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


def optimal_step(sigma: float, f4: float) -> tuple[float, float]:
    """Step that minimises the error of a second derivative by central differences from noisy samples.

    Truncation error ~ (h^2 / 12) |f''''|; noise error ~ sqrt(6) sigma / h^2 (the standard deviation
    of (f_+ - 2 f_0 + f_-) / h^2 with independent noise sigma on each sample). The sum is smallest at
    h* = (12 sqrt(6) sigma / |f''''|)^(1/4), where the two parts are equal and the error is
    2 sqrt(sqrt(6) sigma |f''''| / 12). Returns (h*, error at h*)."""
    A = abs(f4) / 12.0
    B = np.sqrt(6.0) * sigma
    h = (B / A) ** 0.25
    return float(h), float(2.0 * np.sqrt(A * B))


def accel_error_vs_stride(t, pos_true, acc_true, strides, sigma_r=0.0, sigma_th=0.0, sigma_z=0.0, seed=0):
    """RMS error of the acceleration estimated by ``track_from_radar`` for each stride k, from one
    SYNTHETIC noisy radar record (seeded). Returns an array of RMS |a_est - a_true| (m/s^2)."""
    rng = np.random.default_rng(seed)
    r, th, z = radar_samples(pos_true, sigma_r, sigma_th, sigma_z, rng)
    out = []
    for k in strides:
        tm, _, acc = track_from_radar(t, r, th, z, int(k))
        out.append(np.sqrt(np.mean(np.sum((acc - acc_true[k:-k]) ** 2, axis=-1))))
    return np.array(out)


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save radar_track.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    V_h = R_TURN * OMEGA_TURN
    print(f"climbing turn: centre {CENTRE} m, R = {R_TURN:.0f} m, omega = {OMEGA_TURN} rad/s, "
          f"climb {CLIMB} m/s; horizontal speed {V_h:.1f} m/s")
    kap = helix_curvature(R_TURN, OMEGA_TURN, CLIMB)
    V = np.hypot(V_h, CLIMB)
    print(f"closed form: speed {V:.3f} m/s, curvature {kap:.6e} 1/m, radius of curvature {1 / kap:.2f} m, "
          f"a_n = V^2 kappa = {V ** 2 * kap:.4f} m/s^2 ({V ** 2 * kap / G0:.4f} g), a_t = 0")

    # one instant, described three ways
    t1 = np.array([12.0])
    p, v, a = helix_turn(t1)
    r, rd, rdd, th, thd, thdd, z, zd, zdd = cartesian_to_cylindrical(p, v, a)
    vc, ac = cylindrical_kinematics(r, rd, rdd, th, thd, thdd, z, zd, zdd)
    pk = path_kinematics(v, a)
    print(f"\nat t = {t1[0]:.0f} s:")
    print(f"  Cartesian    pos {np.round(p[0], 1)} m, v {np.round(v[0], 3)} m/s, a {np.round(a[0], 4)} m/s^2")
    print(f"  cylindrical  r = {r[0]:.1f} m, theta = {np.degrees(th[0]):.3f} deg, r' = {rd[0]:.3f} m/s, "
          f"theta' = {thd[0]:.6f} rad/s, r'' = {rdd[0]:.4f} m/s^2, theta'' = {thdd[0]:.4e} rad/s^2")
    print(f"               v = ({vc[0, 0]:.3f}, {vc[0, 1]:.3f}, {vc[0, 2]:.3f}), "
          f"a = ({ac[0, 0]:.4f}, {ac[0, 1]:.4f}, {ac[0, 2]:.4f}) in (e_r, e_theta, e_z)")
    print(f"               r theta'^2 = {r[0] * thd[0] ** 2:.4f} m/s^2, 2 r' theta' = {2 * rd[0] * thd[0]:.4f} m/s^2")
    print(f"  path         speed {pk['speed'][0]:.3f} m/s, a_t = {pk['a_t'][0]:.2e}, a_n = {pk['a_n'][0]:.4f} m/s^2, "
          f"rho = {pk['rho'][0]:.2f} m")
    back = cylindrical_to_cartesian(ac, th)
    print(f"  round trip   |a(cyl -> Cartesian) - a| = {np.linalg.norm(back - a):.1e} m/s^2")

    # sampled radar track, noise-free: convergence
    dt = 0.1
    t = np.arange(0.0, 120.0 + dt / 2, dt)
    pos, vel, acc = helix_turn(t)
    r, th, z = radar_samples(pos)
    print("\nnoise-free radar at 10 Hz, acceleration error against the closed form (RMS, m/s^2)")
    print("   stride   h (s)   cylindrical route   Cartesian route")
    ks = [40, 20, 10, 5]               # step halved each time
    e_cyl, e_car = [], []
    for k in ks:
        tm, vv, aa = track_from_radar(t, r, th, z, k)
        _, vv2, aa2 = track_cartesian(t, r, th, z, k)
        e_cyl.append(np.sqrt(np.mean(np.sum((aa - acc[k:-k]) ** 2, axis=-1))))
        e_car.append(np.sqrt(np.mean(np.sum((aa2 - acc[k:-k]) ** 2, axis=-1))))
        print(f"   {k:5d}   {k * dt:5.1f}   {e_cyl[-1]:.4e}          {e_car[-1]:.4e}")
    print(f"   observed orders: cylindrical {np.array2string(observed_order(e_cyl), precision=3)}, "
          f"Cartesian {np.array2string(observed_order(e_car), precision=3)}")
    tm, vv, aa = track_from_radar(t, r, th, z, 10)
    pkm = path_kinematics(vv, aa)
    print(f"   stride 10: mean a_n {np.mean(pkm['a_n']):.4f} m/s^2, mean rho {np.mean(pkm['rho']):.2f} m, "
          f"max |a_t| {np.max(np.abs(pkm['a_t'])):.2e} m/s^2")

    # noisy radar (SYNTHETIC noise): the step trade-off
    sig_r, sig_th, sig_z = 5.0, 1.0e-3, 5.0
    strides = np.array([5, 10, 20, 30, 40, 50, 60, 80, 100, 140])
    err = accel_error_vs_stride(t, pos, acc, strides, sig_r, sig_th, sig_z, seed=1)
    r_mean = float(np.mean(r))
    sig_cross = sig_th * r_mean
    f4 = R_TURN * OMEGA_TURN ** 4
    h_opt, e_opt = optimal_step(np.hypot(sig_r, sig_cross) / np.sqrt(2), f4)
    print(f"\nSYNTHETIC noise sigma_r = {sig_r} m, sigma_theta = {sig_th * 1e3:.1f} mrad "
          f"(= {sig_cross:.1f} m cross-range at the mean range {r_mean:.0f} m), sigma_z = {sig_z} m")
    print("   stride   h (s)   RMS acceleration error (m/s^2)")
    for k, e in zip(strides, err):
        print(f"   {k:5d}   {k * dt:5.1f}   {e:9.3f}")
    i = int(np.argmin(err))
    print(f"   best stride {strides[i]} (h = {strides[i] * dt:.1f} s), error {err[i]:.3f} m/s^2 "
          f"= {err[i] / (R_TURN * OMEGA_TURN ** 2) * 100:.1f} % of the 20 m/s^2 turn acceleration")
    print(f"   error model: f'''' = R omega^4 = {f4:.3g} m/s^4, per-axis sigma {np.hypot(sig_r, sig_cross) / np.sqrt(2):.2f} m: "
          f"h* = {h_opt:.2f} s, error per axis ~ {e_opt:.3f} m/s^2")

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    ax[0].plot(pos[:, 0] / 1e3, pos[:, 1] / 1e3, label="aircraft (helix, plan view)")
    ax[0].plot(0, 0, "k^", ms=9, label="radar")
    ax[0].plot(CENTRE[0] / 1e3, CENTRE[1] / 1e3, "k+", ms=9, label="turn axis")
    ax[0].set_aspect("equal")
    ax[0].set_xlabel("east (km)")
    ax[0].set_ylabel("north (km)")
    ax[0].set_title("radar and track")
    ax[0].legend(fontsize=8)
    ax[1].plot(t, r / 1e3, label="r (km)")
    ax1b = ax[1].twinx()
    ax1b.plot(t, np.degrees(th), "C1", label="theta (deg)")
    ax[1].set_xlabel("t (s)")
    ax[1].set_ylabel("horizontal range r (km)")
    ax1b.set_ylabel("azimuth theta (deg)")
    ax[1].set_title("what the radar sees (cylindrical)")
    ax[2].loglog(np.array(ks) * dt, e_cyl, "o-", label="noise-free, cylindrical route")
    ax[2].loglog(np.array(ks) * dt, e_car, "x--", label="noise-free, Cartesian route")
    ax[2].loglog(strides * dt, err, "s-", label="SYNTHETIC noise")
    ax[2].loglog(np.array(ks) * dt, e_cyl[0] * (np.array(ks) / ks[0]) ** 2, "k:", lw=0.8, label="slope 2")
    ax[2].set_xlabel("difference step h (s)")
    ax[2].set_ylabel("RMS acceleration error (m/s^2)")
    ax[2].set_title("truncation against noise")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "radar_track.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'radar_track.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
