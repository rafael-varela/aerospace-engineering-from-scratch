"""P06.04 Relative Motion & Rotating Frames: reference solution.

A shell is fired 100 km across a rotating Earth. Where does it land compared with the same shot on an
Earth that does not turn, and why? And why does the swing plane of a long pendulum turn once every
31.8 hours in Paris? Both questions are about describing motion in a frame that rotates.

For a point P seen from a frame with origin O, angular velocity omega and angular acceleration alpha,
the inertial acceleration of P is the five-term sum

    a_P = a_O + alpha x r + omega x (omega x r) + 2 omega x v_rel + a_rel,

(origin, Euler, centripetal, Coriolis, relative). On the Earth (a_O = 0 for an Earth-centred frame,
alpha = 0), Newton's law in the rotating frame becomes
a_rel = f - 2 omega x v_rel - omega x (omega x r) with f the gravitational force per unit mass.

The independent check is to integrate the same motion twice: in the rotating Earth-fixed frame with
the Coriolis and centrifugal terms, and in a non-rotating frame with no fictitious terms at all, then
rotate the second answer into the Earth-fixed frame. They must agree to the solver tolerance. A
first-order closed form on a flat Earth and the Foucault pendulum's exact precession rate
Omega sin(latitude) are further checks.

Frames: ECEF-like Earth-fixed axes (z along the spin axis, x through latitude 0, longitude 0) on a
spherical, non-rotating-gravity Earth; local ENU (east, north, up) at the launch point. The Earth is a
sphere of the mean radius with point-mass gravity: no oblateness, no air. Units: SI, radians.

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes rotating_earth.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-06/lesson-04-relative-motion-rotating-frames
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

OMEGA_E = 7.2921159e-5     # rad/s, Earth's sidereal rotation rate
R_E = 6.371e6              # m, mean Earth radius (spherical model)
MU_E = 3.986004418e14      # m^3/s^2, Earth's gravitational parameter
G_SURF = 9.80665           # m/s^2, standard gravity (flat-Earth model)
OMEGA_VEC = np.array([0.0, 0.0, OMEGA_E])


# ---------------------------------------------------------------- the five-term acceleration

def five_term_acceleration(a_O, alpha, omega, r, v_rel, a_rel) -> tuple[np.ndarray, dict]:
    """Inertial acceleration of a point P described in a moving frame (all vectors in one basis):

        a_P = a_O + alpha x r + omega x (omega x r) + 2 omega x v_rel + a_rel.

    Returns (a_P, terms) where terms names each contribution: 'origin', 'euler', 'centripetal',
    'coriolis', 'relative'."""
    a_O, alpha, omega, r, v_rel, a_rel = (np.asarray(a, dtype=float) for a in (a_O, alpha, omega, r, v_rel, a_rel))
    terms = {"origin": a_O,
             "euler": np.cross(alpha, r),
             "centripetal": np.cross(omega, np.cross(omega, r)),
             "coriolis": 2.0 * np.cross(omega, v_rel),
             "relative": a_rel}
    return sum(terms.values()), terms


def rotating_frame_acceleration(r, v, f, omega=OMEGA_VEC) -> np.ndarray:
    """Acceleration seen in a frame turning at constant omega about a fixed origin, from the five-term
    law solved for a_rel: a_rel = f - 2 omega x v - omega x (omega x r), with f the real force per unit
    mass. The two subtracted terms are the Coriolis and centrifugal 'fictitious' accelerations."""
    r, v, f, omega = (np.asarray(a, dtype=float) for a in (r, v, f, omega))
    return f - 2.0 * np.cross(omega, v) - np.cross(omega, np.cross(omega, r))


# ---------------------------------------------------------------- frames and states

def rot_z(angle: float) -> np.ndarray:
    """Rotation matrix by ``angle`` about z (active, counterclockwise)."""
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def inertial_to_rotating(r_I, v_I, t: float, Omega: float = OMEGA_E) -> tuple[np.ndarray, np.ndarray]:
    """Position and velocity in the Earth-fixed frame from the non-rotating frame (the two coincide at
    t = 0): r_R = Rz(-Omega t) r_I and v_R = Rz(-Omega t) v_I - omega x r_R (the transport theorem)."""
    Q = rot_z(-Omega * t)
    r_R = Q @ np.asarray(r_I, dtype=float)
    v_R = Q @ np.asarray(v_I, dtype=float) - np.cross([0.0, 0.0, Omega], r_R)
    return r_R, v_R


def rotating_to_inertial(r_R, v_R, t: float, Omega: float = OMEGA_E) -> tuple[np.ndarray, np.ndarray]:
    """The inverse of ``inertial_to_rotating``: v_I = Rz(Omega t) (v_R + omega x r_R)."""
    r_R, v_R = np.asarray(r_R, dtype=float), np.asarray(v_R, dtype=float)
    Q = rot_z(Omega * t)
    return Q @ r_R, Q @ (v_R + np.cross([0.0, 0.0, Omega], r_R))


def enu_basis(lat: float, lon: float) -> np.ndarray:
    """Rows: the east, north and up unit vectors at (lat, lon) on a sphere, in Earth-fixed components."""
    sl, cl, so, co = np.sin(lat), np.cos(lat), np.sin(lon), np.cos(lon)
    return np.array([[-so, co, 0.0], [-sl * co, -sl * so, cl], [cl * co, cl * so, sl]])


def launch_state(lat: float, lon: float, speed: float, elevation: float, azimuth: float,
                 R: float = R_E) -> tuple[np.ndarray, np.ndarray]:
    """Earth-fixed position on the sphere and Earth-relative velocity of a launch at ``speed``,
    ``elevation`` above the horizon and ``azimuth`` clockwise from north. Returns (r0, v0)."""
    E, N, U = enu_basis(lat, lon)
    d = np.cos(elevation) * (np.sin(azimuth) * E + np.cos(azimuth) * N) + np.sin(elevation) * U
    return R * U, speed * d


def gravity(r, mu: float = MU_E) -> np.ndarray:
    """Point-mass gravity -mu r / |r|^3."""
    r = np.asarray(r, dtype=float)
    return -mu * r / np.linalg.norm(r) ** 3


# ---------------------------------------------------------------- the two integrations

def _ground(R):
    def hit(t, y):
        return np.linalg.norm(y[:3]) - R
    hit.terminal = True
    hit.direction = -1.0
    return hit


def simulate_rotating(r0, v0, Omega: float = OMEGA_E, mu: float = MU_E, R: float = R_E, t_max: float = 3600.0,
                      rtol: float = 1e-11, atol: float = 1e-6, coriolis: bool = True, centrifugal: bool = True):
    """Integrate in the Earth-fixed (rotating) frame: r'' = g(r) - 2 omega x r' - omega x (omega x r),
    from the Earth-relative state (r0, v0), until the projectile returns to |r| = R. Returns
    (t_impact, r_impact, v_impact) in the Earth-fixed frame.

    ``coriolis`` and ``centrifugal`` switch the two fictitious terms off one at a time, for experiments
    only: the physical answer needs both. (Switching off only the Coriolis term isolates the deflection
    that the flat-Earth formula predicts; the centrifugal term is a static field that on the real Earth
    is part of what a plumb line calls gravity.)"""
    w = np.array([0.0, 0.0, Omega])

    def rhs(t, y):
        r, v = y[:3], y[3:]
        a = gravity(r, mu)
        if coriolis:
            a = a - 2.0 * np.cross(w, v)
        if centrifugal:
            a = a - np.cross(w, np.cross(w, r))
        return np.concatenate([v, a])
    sol = solve_ivp(rhs, (0.0, t_max), np.concatenate([r0, v0]), rtol=rtol, atol=atol,
                    events=_ground(R), method="DOP853")
    if sol.status != 1:
        raise RuntimeError("no impact within t_max")
    y = sol.y_events[0][0]
    return float(sol.t_events[0][0]), y[:3], y[3:]


def simulate_inertial(r0, v0, Omega: float = OMEGA_E, mu: float = MU_E, R: float = R_E, t_max: float = 3600.0,
                      rtol: float = 1e-11, atol: float = 1e-6):
    """The same shot integrated in the non-rotating frame, where r'' = g(r) and nothing else: start with
    the inertial velocity v_I = v0 + omega x r0 (the gun moves with the ground), integrate to impact and
    rotate the impact point back into the Earth-fixed frame, which has turned by Omega t meanwhile.
    Returns (t_impact, r_impact, v_impact) in the Earth-fixed frame."""
    rI0, vI0 = rotating_to_inertial(r0, v0, 0.0, Omega)

    def rhs(t, y):
        return np.concatenate([y[3:], gravity(y[:3], mu)])
    sol = solve_ivp(rhs, (0.0, t_max), np.concatenate([rI0, vI0]), rtol=rtol, atol=atol,
                    events=_ground(R), method="DOP853")
    if sol.status != 1:
        raise RuntimeError("no impact within t_max")
    t = float(sol.t_events[0][0])
    rR, vR = inertial_to_rotating(sol.y_events[0][0][:3], sol.y_events[0][0][3:], t, Omega)
    return t, rR, vR


def range_and_cross(r_imp, lat: float, lon: float, azimuth: float, R: float = R_E) -> tuple[float, float]:
    """Down-range (along the great circle of the launch azimuth) and cross-range (positive to the right
    of the line of fire) of an impact point, as arc lengths on the sphere (m)."""
    E, N, U = enu_basis(lat, lon)
    h = np.sin(azimuth) * E + np.cos(azimuth) * N          # horizontal direction of fire
    right = np.cross(h, U)                                  # normal of the launch plane, to the right
    p = np.asarray(r_imp, dtype=float) / np.linalg.norm(r_imp)
    cross = R * np.arcsin(np.clip(p @ right, -1, 1))
    down = R * np.arctan2(p @ h, p @ U)
    return float(down), float(cross)


# ---------------------------------------------------------------- flat Earth, first order

def omega_enu(lat: float, Omega: float = OMEGA_E) -> np.ndarray:
    """Earth's angular velocity in local ENU components: Omega (0, cos lat, sin lat)."""
    return Omega * np.array([0.0, np.cos(lat), np.sin(lat)])


def flat_earth_rotating(v0_enu, lat: float, g: float = G_SURF, Omega: float = OMEGA_E, rtol: float = 1e-11):
    """Flat Earth, uniform g (the centrifugal term absorbed into g), Coriolis only, in local ENU:
    r'' = (0, 0, -g) - 2 omega x r', from the origin until the projectile is back at z = 0.
    Returns (t_impact, r_impact)."""
    w = omega_enu(lat, Omega)
    gv = np.array([0.0, 0.0, -g])

    def rhs(t, y):
        return np.concatenate([y[3:], gv - 2.0 * np.cross(w, y[3:])])

    def ground(t, y):
        return y[2]
    ground.terminal = True
    ground.direction = -1.0
    T0 = 2.0 * v0_enu[2] / g
    sol = solve_ivp(rhs, (0.0, 2 * T0), np.concatenate([[0, 0, 0], v0_enu]), rtol=rtol, atol=1e-9,
                    events=ground, method="DOP853", first_step=1e-3)
    return float(sol.t_events[0][0]), sol.y_events[0][0][:3]


def first_order_deflection(v0_enu, lat: float, t: float, g: float = G_SURF, Omega: float = OMEGA_E) -> np.ndarray:
    """First-order (in Omega) displacement caused by the Coriolis term on a flat Earth at time t:
    r1'' = -2 omega x (v0 + g t) with g = (0, 0, -g) integrates to r1 = -omega x (v0 t^2 + g t^3 / 3)."""
    w = omega_enu(lat, Omega)
    gv = np.array([0.0, 0.0, -g])
    return -np.cross(w, np.asarray(v0_enu, dtype=float) * t ** 2 + gv * t ** 3 / 3.0)


def first_order_impact_shift(v0_enu, lat: float, g: float = G_SURF, Omega: float = OMEGA_E) -> np.ndarray:
    """Horizontal shift (east, north) of the impact point to first order in Omega. At the vacuum flight
    time T = 2 v_up / g the Coriolis displacement is r1(T); its vertical part also changes the flight
    time by dT = -r1_z / v_z(T) = r1_z / v_up, which moves the impact along the horizontal velocity
    (the Eotvos effect on range). Returns r1_h(T) + v_h dT."""
    v0 = np.asarray(v0_enu, dtype=float)
    T = 2.0 * v0[2] / g
    r1 = first_order_deflection(v0, lat, T, g, Omega)
    dT = r1[2] / v0[2]
    return r1[:2] + v0[:2] * dT


# ---------------------------------------------------------------- Foucault pendulum

def foucault_rhs(t, y, w0: float, Om_z: float):
    """Small-amplitude pendulum in the local horizontal plane of a rotating Earth, y = [x, y, x', y']
    (x east, y north): r'' = -w0^2 r - 2 Om_z e_z x r', with Om_z = Omega sin(lat) the vertical component
    of the Earth's rotation (the horizontal component only changes the tension, to this order) and the
    centrifugal term absorbed into g."""
    x, yy, vx, vy = y
    return [vx, vy, -w0 ** 2 * x + 2.0 * Om_z * vy, -w0 ** 2 * yy - 2.0 * Om_z * vx]


def simulate_foucault(L: float, lat: float, amp: float, t_end: float, g: float = G_SURF,
                      Omega: float = OMEGA_E, n_out: int = 20001, rtol: float = 1e-10, atol: float = 1e-10):
    """Release the pendulum from rest (Earth-relative) at (amp, 0) and integrate in the rotating frame.
    Returns (t, x, y)."""
    w0 = np.sqrt(g / L)
    t = np.linspace(0.0, t_end, n_out)
    sol = solve_ivp(foucault_rhs, (0.0, t_end), [amp, 0.0, 0.0, 0.0], args=(w0, Omega * np.sin(lat)),
                    t_eval=t, rtol=rtol, atol=atol, method="DOP853")
    return sol.t, sol.y[0], sol.y[1]


def foucault_closed_form(t, L: float, lat: float, amp: float, g: float = G_SURF, Omega: float = OMEGA_E):
    """Exact solution of ``foucault_rhs`` from the non-rotating view. In axes that do not turn about the
    vertical, s'' = -(w0^2 + Om_z^2) s: an isotropic oscillator at W = sqrt(w0^2 + Om_z^2) started at
    (amp, 0) with velocity (0, Om_z amp) (the floor carries the bob), so
    s = (amp cos W t, (Om_z amp / W) sin W t), a very thin ellipse. Seen from the Earth, r = Rz(-Om_z t) s:
    the swing plane turns clockwise (in the north) at exactly Om_z. Returns (x, y)."""
    t = np.asarray(t, dtype=float)
    Oz = Omega * np.sin(lat)
    W = np.sqrt(g / L + Oz ** 2)
    sx, sy = amp * np.cos(W * t), Oz * amp / W * np.sin(W * t)
    c, s = np.cos(-Oz * t), np.sin(-Oz * t)
    return c * sx - s * sy, s * sx + c * sy


def swing_plane_rate(t, x, y) -> float:
    """Precession rate of the swing plane (rad/s, counterclockwise positive), from the azimuth of the
    bob at each turning point (local maxima of x^2 + y^2), taken mod pi, unwrapped, and fitted by a
    straight line."""
    r2 = x * x + y * y
    i = np.where((r2[1:-1] > r2[:-2]) & (r2[1:-1] >= r2[2:]))[0] + 1
    ang = np.unwrap(2.0 * np.arctan2(y[i], x[i])) / 2.0
    return float(np.polyfit(t[i], ang, 1)[0])


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save rotating_earth.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    lat, lon = np.radians(49.0), np.radians(2.5)
    speed, elev = 1000.0, np.radians(45.0)
    print(f"Earth rotation {OMEGA_E:.7e} rad/s; launch at {np.degrees(lat):.1f} deg N, "
          f"{speed:.0f} m/s at {np.degrees(elev):.0f} deg elevation, vacuum, spherical Earth")

    # the five terms for an airliner and for the shell
    E, N, U = enu_basis(lat, lon)
    r0, v_air = R_E * U, 250.0 * N
    a, terms = five_term_acceleration(np.zeros(3), np.zeros(3), OMEGA_VEC, r0, v_air, np.zeros(3))
    print(f"\nairliner flying north at 250 m/s: |Coriolis| = {np.linalg.norm(terms['coriolis']):.4f} m/s^2 "
          f"({np.linalg.norm(terms['coriolis']) / G_SURF * 1e3:.2f} mg), "
          f"|centripetal| = {np.linalg.norm(terms['centripetal']):.4f} m/s^2; "
          f"2 Omega V sin(lat) = {2 * OMEGA_E * 250 * np.sin(lat):.4f}")

    print("\nlong-range shot: rotating-frame vs inertial-frame integration (the independent check)")
    print(" azimuth  flight (s)  range (km)  |rotating - inertial| (m)")
    shots = (("north", 0.0), ("east", np.pi / 2), ("south", np.pi), ("west", 3 * np.pi / 2))
    for name, az in shots:
        r0, v0 = launch_state(lat, lon, speed, elev, az)
        tr, rr, _ = simulate_rotating(r0, v0)
        ti, ri, _ = simulate_inertial(r0, v0)
        print(f" {name:7s}  {tr:9.3f}  {range_and_cross(rr, lat, lon, az)[0] / 1e3:9.3f}  {np.linalg.norm(rr - ri):19.2e}"
              f"   (flight times differ by {abs(tr - ti):.1e} s)")
    r0, v0 = launch_state(lat, lon, speed, elev, 0.0)
    t0, rn, _ = simulate_rotating(r0, v0, Omega=0.0)
    print(f" non-rotating Earth: range {range_and_cross(rn, lat, lon, 0.0)[0] / 1e3:.3f} km in {t0:.3f} s; "
          f"flat Earth, uniform g: {speed ** 2 * np.sin(2 * elev) / G_SURF / 1e3:.3f} km in {2 * speed * np.sin(elev) / G_SURF:.3f} s")

    print("\nimpact shifts (m): Coriolis alone (full minus no-Coriolis), centrifugal alone, and the flat-Earth Coriolis models")
    print(" azimuth   Coriolis: down  cross(+right)   centrifugal: down  cross   flat 1st order: down  cross"
          "   flat sim: down  cross   no rotation at all: down  cross")
    for name, az in shots:
        r0, v0 = launch_state(lat, lon, speed, elev, az)
        _, rr, _ = simulate_rotating(r0, v0)
        _, rnc, _ = simulate_rotating(r0, v0, coriolis=False)
        _, rncf, _ = simulate_rotating(r0, v0, centrifugal=False)
        _, rn, _ = simulate_rotating(r0, v0, Omega=0.0)
        d, c = range_and_cross(rr, lat, lon, az)
        dc, cc = range_and_cross(rnc, lat, lon, az)
        df, cf = range_and_cross(rncf, lat, lon, az)
        dn, cn = range_and_cross(rn, lat, lon, az)
        v_enu = np.array([np.sin(az), np.cos(az), 0.0]) * speed * np.cos(elev) + np.array([0, 0, speed * np.sin(elev)])
        sh = first_order_impact_shift(v_enu, lat)
        h = np.array([np.sin(az), np.cos(az)])
        rgt = np.array([np.cos(az), -np.sin(az)])
        _, rf = flat_earth_rotating(v_enu, lat)
        _, rf0 = flat_earth_rotating(v_enu, lat, Omega=0.0)
        dfl = rf[:2] - rf0[:2]
        print(f" {name:7s}  {d - dc:15.1f}  {c - cc:13.1f}   {d - df:17.1f}  {c - cf:5.1f}   {sh @ h:20.1f}  {sh @ rgt:5.1f}"
              f"   {dfl @ h:14.1f}  {dfl @ rgt:5.1f}   {d - dn:24.1f}  {c - cn:5.1f}")
    print(f" centrifugal field at {np.degrees(lat):.0f} deg: horizontal (toward the equator) "
          f"{OMEGA_E ** 2 * R_E * np.cos(lat) * np.sin(lat):.4f} m/s^2, vertical {OMEGA_E ** 2 * R_E * np.cos(lat) ** 2:.4f} m/s^2; "
          f"(1/2) a_h t^2 over 146 s = {0.5 * OMEGA_E ** 2 * R_E * np.cos(lat) * np.sin(lat) * 146.17 ** 2:.0f} m")

    # southern hemisphere and the equator
    print("\nCoriolis cross-range of a northward shot against latitude (spherical: full minus no-Coriolis)")
    print("   lat (deg)   cross (m, +right)   flat 1st order (m)")
    lats = np.radians([-60, -30, 0, 30, 49, 60, 80])
    cr = []
    for la in lats:
        r0, v0 = launch_state(la, lon, speed, elev, 0.0)
        _, rr, _ = simulate_rotating(r0, v0)
        _, rn, _ = simulate_rotating(r0, v0, coriolis=False)
        c = range_and_cross(rr, la, lon, 0.0)[1] - range_and_cross(rn, la, lon, 0.0)[1]
        v_enu = np.array([0.0, speed * np.cos(elev), speed * np.sin(elev)])
        sh = first_order_impact_shift(v_enu, la)
        cr.append(c)
        print(f"   {np.degrees(la):9.0f}   {c:17.1f}   {sh[0]:18.1f}")

    # Foucault pendulum
    Lp, latp, amp = 67.0, np.radians(48.85), 3.0
    Oz = OMEGA_E * np.sin(latp)
    print(f"\nFoucault pendulum L = {Lp} m at {np.degrees(latp):.2f} deg N, amplitude {amp} m: "
          f"period {2 * np.pi * np.sqrt(Lp / G_SURF):.3f} s")
    print(f"  predicted precession -Omega sin(lat) = {-Oz:.6e} rad/s = {-np.degrees(Oz) * 3600:.3f} deg/h; "
          f"full turn in {2 * np.pi / Oz / 3600:.2f} h")
    t, x, y = simulate_foucault(Lp, latp, amp, 3600.0, n_out=72001)
    rate = swing_plane_rate(t, x, y)
    xc, yc = foucault_closed_form(t, Lp, latp, amp)
    print(f"  rotating-frame integration over 1 h: rate {rate:.6e} rad/s ({rate / -Oz - 1:+.1e} relative), "
          f"max |r - closed form| = {np.max(np.hypot(x - xc, y - yc)):.1e} m")
    print(f"  swing plane turned {np.degrees(rate * 3600):.3f} deg in 1 h; minor axis of the ellipse "
          f"Om_z amp / W = {Oz * amp / np.sqrt(G_SURF / Lp + Oz ** 2) * 1e3:.3f} mm")

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.4))
    ax[0].plot(np.degrees(lats), cr, "o-", label="spherical Earth, Coriolis term")
    ax[0].plot(np.degrees(lats), [first_order_impact_shift(np.array([0, speed * np.cos(elev), speed * np.sin(elev)]), la)[0]
                                   for la in lats], "k--", label="flat Earth, first order")
    ax[0].axhline(0, color="0.6", lw=0.6)
    ax[0].set_xlabel("latitude (deg)")
    ax[0].set_ylabel("cross-range of a northward shot (m, + right)")
    ax[0].set_title("Coriolis deflection, 1000 m/s at 45 deg")
    ax[0].legend(fontsize=8)
    ax[1].plot(x[:3000], y[:3000] * 1e3, lw=0.5)
    ax[1].set_xlabel("x east (m)")
    ax[1].set_ylabel("y north (mm)")
    ax[1].set_title("Foucault bob, first 150 s (rotating frame)")
    k = slice(0, None, 50)
    ax[2].plot(x[k], y[k], ",", alpha=0.6)
    ax[2].set_aspect("equal")
    ax[2].set_xlabel("x east (m)")
    ax[2].set_ylabel("y north (m)")
    ax[2].set_title("swing plane over 1 h (turns clockwise)")
    fig.tight_layout()
    fig.savefig(out_dir / "rotating_earth.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'rotating_earth.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
