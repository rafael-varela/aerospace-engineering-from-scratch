"""P04.01 Kinematics & Newton's Laws: reference solution.

A projectile with quadratic air drag, built from Newton's second law and checked against the
closed forms that exist in three limits (no drag, purely vertical motion, purely horizontal motion)
and against `scipy.integrate.solve_ivp` with tight tolerances. The same code predicts the drop and
toss experiment of `projects/p04-drag-drop/` (ping-pong ball against a dense object) and fits the
terminal speed, hence the drag coefficient, from measured fall times. The wind triangle (relative
motion) is the second aerospace example.

State vector for the projectile: y = [x, z, vx, vz], x horizontal, z vertical upward, SI units.
Drag parameter k = rho C_D A / (2 m) [1/m]: the acceleration from drag is -k |v| v.

Run:  python solution.py [out_dir]      (prints the verification tables, writes projectile_drag.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-04/lesson-01-kinematics-newtons-laws
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares, minimize_scalar

G0 = 9.80665          # m/s^2, standard gravity
RHO_AIR = 1.20        # kg/m^3, room air at about 20 C and sea level

# The two objects of the experiment (nominal values; measure your own and replace them).
PING_PONG = dict(mass=0.0027, diameter=0.040, cd=0.50)     # kg, m, dimensionless
BEAN_BAG = dict(mass=0.100, diameter=0.080, cd=1.0)        # kg, m, a soft lumpy object: cd is a guess


# ---------------------------------------------------------------- Newton's second law with quadratic drag

# LEARNER
def drag_parameter(mass: float, diameter: float, cd: float = 0.5, rho: float = RHO_AIR) -> float:
    """k = rho C_D A / (2 m) in 1/m, with A = pi d^2 / 4 the frontal area of a sphere.

    The drag acceleration is k |v| v, so k has units of 1/length: 1/k is the distance over which
    drag alone (no gravity) halves the speed of a body launched at any speed (see
    `horizontal_decay_state`), the natural length scale of the problem.
    """
    area = math.pi * diameter ** 2 / 4.0
    return rho * cd * area / (2.0 * mass)


# LEARNER
def terminal_speed(k: float, g: float = G0) -> float:
    """Speed at which drag balances weight for a falling body: v_t = sqrt(g / k).

    Returns inf when k == 0 (no drag: nothing stops the fall from accelerating).
    """
    if k <= 0.0:
        return math.inf
    return math.sqrt(g / k)


# LEARNER
def projectile_rhs(t: float, y: np.ndarray, k: float, g: float = G0) -> np.ndarray:
    """Right-hand side of the equations of motion for y = [x, z, vx, vz]:

        dx/dt = vx,  dz/dt = vz,  dvx/dt = -k |v| vx,  dvz/dt = -g - k |v| vz,

    with |v| = sqrt(vx^2 + vz^2). This is Newton's second law divided by the mass: the only forces
    are the weight (down) and the drag (opposite to the velocity, magnitude k m |v|^2).
    """
    vx, vz = y[2], y[3]
    speed = math.hypot(vx, vz)
    return np.array([vx, vz, -k * speed * vx, -g - k * speed * vz])


# ---------------------------------------------------------------- closed forms in three limits

# LEARNER
def vacuum_state(t: np.ndarray | float, v0: float, angle_deg: float, g: float = G0) -> np.ndarray:
    """Drag-free projectile launched from the origin: constant acceleration in z, none in x.

    Returns an array of shape (4,) for scalar t or (n, 4) for an array, columns [x, z, vx, vz].
    """
    t = np.asarray(t, dtype=float)
    a = math.radians(angle_deg)
    vx0, vz0 = v0 * math.cos(a), v0 * math.sin(a)
    x = vx0 * t
    z = vz0 * t - 0.5 * g * t ** 2
    vx = np.full_like(t, vx0)
    vz = vz0 - g * t
    return np.stack([x, z, vx, vz], axis=-1)


# LEARNER
def vacuum_time_of_flight(v0: float, angle_deg: float, g: float = G0) -> float:
    """Time to return to z = 0 without drag: 2 v0 sin(angle) / g."""
    return 2.0 * v0 * math.sin(math.radians(angle_deg)) / g


# LEARNER
def vacuum_range(v0: float, angle_deg: float, g: float = G0) -> float:
    """Horizontal distance to z = 0 without drag: v0^2 sin(2 angle) / g."""
    return v0 ** 2 * math.sin(math.radians(2.0 * angle_deg)) / g


# LEARNER
def vertical_drop_state(t: np.ndarray | float, v_t: float, g: float = G0) -> tuple[np.ndarray, np.ndarray]:
    """Fall from rest with quadratic drag, distance fallen and speed (both positive downward):

        v(t) = v_t tanh(g t / v_t),   s(t) = (v_t^2 / g) ln cosh(g t / v_t).

    With v_t = inf (no drag) returns g t^2 / 2 and g t. Works on arrays.
    """
    t = np.asarray(t, dtype=float)
    if math.isinf(v_t):
        return 0.5 * g * t ** 2, g * t
    u = g * t / v_t
    # ln cosh u without overflow: for large u, ln cosh u = u - ln 2 + ln(1 + e^{-2u})
    lncosh = np.where(u > 20.0, u - math.log(2.0) + np.log1p(np.exp(-2.0 * np.minimum(u, 700.0))), np.log(np.cosh(np.minimum(u, 20.0))))
    return v_t ** 2 / g * lncosh, v_t * np.tanh(u)


# LEARNER
def vertical_drop_time(h: float, v_t: float, g: float = G0) -> float:
    """Time to fall a height h from rest with quadratic drag (terminal speed v_t):

        t = (v_t / g) arccosh(exp(g h / v_t^2)),

    evaluated as (v_t / g) log1p(x + sqrt(x (x + 2))) with x = expm1(g h / v_t^2) so that the weak-drag
    limit keeps its digits; returns sqrt(2 h / g) for v_t = inf.
    """
    if math.isinf(v_t):
        return math.sqrt(2.0 * h / g)
    a = g * h / v_t ** 2
    if a > 300.0:                      # exp(a) would overflow: arccosh(e^a) = a + ln 2 to double precision
        return v_t / g * (a + math.log(2.0))
    x = math.expm1(a)
    return v_t / g * math.log1p(x + math.sqrt(x * (x + 2.0)))


# LEARNER
def horizontal_decay_state(t: np.ndarray | float, v0: float, k: float) -> tuple[np.ndarray, np.ndarray]:
    """Horizontal motion with quadratic drag and no gravity (a puck on ice in air, or the x-motion
    of a projectile in the instant before gravity matters):

        v(t) = v0 / (1 + k v0 t),   x(t) = ln(1 + k v0 t) / k.

    Returns (x, v); for k = 0 returns (v0 t, v0).
    """
    t = np.asarray(t, dtype=float)
    if k == 0.0:
        return v0 * t, np.full_like(t, v0)
    return np.log1p(k * v0 * t) / k, v0 / (1.0 + k * v0 * t)


# ---------------------------------------------------------------- the simulation

@dataclass
class Trajectory:
    """Result of `simulate_projectile`: sampled states plus the impact quantities."""
    t: np.ndarray
    x: np.ndarray
    z: np.ndarray
    vx: np.ndarray
    vz: np.ndarray
    time_of_flight: float
    range: float
    impact_speed: float
    apex_height: float


# LEARNER
def simulate_projectile(v0: float, angle_deg: float, k: float, g: float = G0, z0: float = 0.0,
                        rtol: float = 1e-10, n_samples: int = 400, max_time: float = 60.0) -> Trajectory:
    """Integrate `projectile_rhs` from (0, z0) with speed v0 at `angle_deg` above the horizontal until
    the projectile reaches z = 0 going down, with `scipy.integrate.solve_ivp` (RK45, rtol given,
    atol = 1e-12) and a terminal event on z = 0. Samples the dense output at `n_samples` times.

    The event is only armed once the body is moving downward (direction = -1), so a launch from
    z0 = 0 does not trigger it at t = 0. Raises ValueError if the ground is not reached within
    `max_time` seconds.
    """
    a = math.radians(angle_deg)
    y0 = np.array([0.0, z0, v0 * math.cos(a), v0 * math.sin(a)])

    def hit_ground(t, y, *args):
        return y[1]

    hit_ground.terminal = True
    hit_ground.direction = -1
    sol = solve_ivp(projectile_rhs, (0.0, max_time), y0, args=(k, g), method="RK45", rtol=rtol, atol=1e-12,
                    events=hit_ground, dense_output=True)
    if not sol.t_events[0].size:
        raise ValueError("the projectile did not reach the ground within max_time")
    t_end = float(sol.t_events[0][0])
    t = np.linspace(0.0, t_end, n_samples)
    Y = sol.sol(t)
    y_end = sol.y_events[0][0]
    return Trajectory(t=t, x=Y[0], z=Y[1], vx=Y[2], vz=Y[3], time_of_flight=t_end, range=float(y_end[0]),
                      impact_speed=float(math.hypot(y_end[2], y_end[3])), apex_height=float(Y[1].max()))


# LEARNER
def drag_range(v0: float, angle_deg: float, k: float, g: float = G0) -> float:
    """Horizontal range to z = 0 with drag, from `simulate_projectile` (rtol 1e-8 is enough here)."""
    return simulate_projectile(v0, angle_deg, k, g, rtol=1e-8, n_samples=2).range


# LEARNER
def optimal_angle(v0: float, k: float, g: float = G0) -> tuple[float, float]:
    """Launch angle (degrees) that maximises the range, and that range, by a bounded scalar search
    on [5, 85] degrees (`scipy.optimize.minimize_scalar`, xatol 1e-3 degrees). Without drag the
    answer is 45 degrees; with drag it is lower.
    """
    res = minimize_scalar(lambda a: -drag_range(v0, a, k, g), bounds=(5.0, 85.0), method="bounded",
                          options={"xatol": 1e-3})
    return float(res.x), float(-res.fun)


# ---------------------------------------------------------------- relative motion: the wind triangle

# LEARNER
def wind_triangle(tas: float, heading_deg: float, wind_speed: float, wind_from_deg: float) -> tuple[float, float]:
    """Ground speed and track (degrees clockwise from north) of an aircraft flying at true airspeed
    `tas` on `heading_deg`, in a wind of `wind_speed` blowing FROM `wind_from_deg` (meteorological
    convention). Velocity addition: v_ground = v_air + v_wind, in east/north components.
    """
    h, w = math.radians(heading_deg), math.radians(wind_from_deg)
    air = np.array([tas * math.sin(h), tas * math.cos(h)])
    wind = np.array([-wind_speed * math.sin(w), -wind_speed * math.cos(w)])
    ground = air + wind
    gs = float(np.hypot(*ground))
    track = math.degrees(math.atan2(ground[0], ground[1])) % 360.0
    return gs, track


# LEARNER
def heading_for_track(tas: float, track_deg: float, wind_speed: float, wind_from_deg: float) -> tuple[float, float]:
    """Heading to fly (degrees) so that the track over the ground is `track_deg`, and the resulting
    ground speed. Closed form: the wind correction angle satisfies
    sin(WCA) = (W / TAS) sin(wind_from - track), heading = track + WCA,
    ground speed = TAS cos(WCA) - W cos(wind_from - track).
    Raises ValueError when the crosswind exceeds the airspeed (no heading can hold the track).
    """
    d = math.radians(wind_from_deg - track_deg)
    s = wind_speed * math.sin(d) / tas
    if abs(s) > 1.0:
        raise ValueError("crosswind component exceeds the true airspeed: the track cannot be held")
    wca = math.asin(s)
    gs = tas * math.cos(wca) - wind_speed * math.cos(d)
    if gs <= 0.0:
        raise ValueError("the wind is stronger than the airspeed along the track")
    return (track_deg + math.degrees(wca)) % 360.0, gs


# ---------------------------------------------------------------- the experiment: predict, reduce, fit

def frame_timing(n_frames: int, fps: float) -> tuple[float, float]:
    """Time from a frame count and its Type B standard uncertainty 0.41/fps (from P00.15)."""
    return n_frames / fps, math.sqrt(2.0) * 0.5 / (math.sqrt(3.0) * fps)


def reduce_timings(times: np.ndarray, u_single: float = 0.0) -> tuple[float, float]:
    """Mean and combined standard uncertainty of repeated timings: Type A s/sqrt(n) in quadrature
    with the Type B term `u_single` that does not average down (from P00.15)."""
    x = np.asarray(times, dtype=float)
    u_a = x.std(ddof=1) / math.sqrt(x.size) if x.size > 1 else 0.0
    return float(x.mean()), math.hypot(u_a, u_single)


def normalized_error(a: float, u_a: float, b: float, u_b: float) -> float:
    """E_n = |a - b| / sqrt(u_a^2 + u_b^2); inf when both uncertainties are zero and a != b."""
    u = math.hypot(u_a, u_b)
    if u == 0.0:
        return 0.0 if a == b else math.inf
    return abs(a - b) / u


# LEARNER
def predict_fall_times(heights: np.ndarray, k: float, g: float = G0) -> tuple[np.ndarray, np.ndarray]:
    """Vacuum and drag fall times for each height: (sqrt(2h/g), vertical_drop_time(h, v_t))."""
    heights = np.asarray(heights, dtype=float)
    v_t = terminal_speed(k, g)
    t_vac = np.sqrt(2.0 * heights / g)
    t_drag = np.array([vertical_drop_time(h, v_t, g) for h in heights])
    return t_vac, t_drag


# LEARNER
def fit_terminal_speed(heights: np.ndarray, times: np.ndarray, u_times: np.ndarray | None = None,
                       g: float = G0) -> tuple[float, float]:
    """Least-squares fit of the terminal speed to measured mean fall times at several heights,
    minimising sum ((t_i - vertical_drop_time(h_i, v_t)) / u_i)^2 with `scipy.optimize.least_squares`
    (u_i = 1 when not given). Returns (v_t, u_vt) with the standard uncertainty from the Jacobian
    at the solution: u_vt^2 = s^2 (J^T J)^-1, s^2 the reduced chi-square (set to 1 when weights are given
    and the fit has no degrees of freedom to estimate it from). The fit is done in w = 1/v_t^2 so that
    the no-drag limit (w = 0) is a finite parameter value.
    """
    heights = np.asarray(heights, dtype=float)
    times = np.asarray(times, dtype=float)
    u = np.ones_like(times) if u_times is None else np.asarray(u_times, dtype=float)

    def model(w):
        v_t = math.inf if w[0] <= 0.0 else 1.0 / math.sqrt(w[0])
        return np.array([vertical_drop_time(h, v_t, g) for h in heights])

    def resid(w):
        return (times - model(w)) / u

    w0 = np.array([1.0 / 100.0])         # start at v_t = 10 m/s
    res = least_squares(resid, w0, bounds=(0.0, np.inf), xtol=1e-14, ftol=1e-14, gtol=1e-14)
    w = float(res.x[0])
    J = res.jac                          # d resid / d w, shape (n, 1)
    dof = max(times.size - 1, 1)
    s2 = float(res.fun @ res.fun) / dof if u_times is None else max(float(res.fun @ res.fun) / dof, 1.0)
    cov_w = s2 / float((J.T @ J)[0, 0])
    v_t = 1.0 / math.sqrt(w) if w > 0.0 else math.inf
    # v_t = w^-1/2: dv/dw = -v_t / (2 w)
    u_vt = math.sqrt(cov_w) * v_t / (2.0 * w) if w > 0.0 else math.inf
    return v_t, u_vt


# LEARNER
def drag_coefficient_from_terminal_speed(v_t: float, mass: float, diameter: float, rho: float = RHO_AIR,
                                         g: float = G0) -> float:
    """Invert v_t = sqrt(2 m g / (rho C_D A)): C_D = 2 m g / (rho A v_t^2)."""
    area = math.pi * diameter ** 2 / 4.0
    return 2.0 * mass * g / (rho * area * v_t ** 2)


def make_synthetic_drop_data(k: float, heights: np.ndarray, fps: float, n_drops: int = 5, seed: int = 0,
                             release_bias_s: float = 0.0, g: float = G0) -> list[np.ndarray]:
    """SYNTHETIC frame-counted fall times, for testing the analysis pipeline only. Not a measurement.

    For each height the true drag fall time (plus an optional constant release bias) is quantised to
    whole frames at a random phase, so the data are consistent with the drag model by construction.
    Nothing produced here may appear in a report as data.
    """
    rng = np.random.default_rng(seed)
    v_t = terminal_speed(k, g)
    out = []
    for h in np.asarray(heights, dtype=float):
        t_true = vertical_drop_time(h, v_t, g) + release_bias_s
        phase = rng.uniform(0.0, 1.0, n_drops)             # where inside a frame the release happened
        frames = np.floor((t_true * fps + phase)) - np.floor(phase)
        out.append(frames / fps)
    return out


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Run the lesson's numbers, print them and save projectile_drag.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    k_pp = drag_parameter(**PING_PONG)
    k_bb = drag_parameter(**BEAN_BAG)
    vt_pp, vt_bb = terminal_speed(k_pp), terminal_speed(k_bb)
    print(f"ping-pong ball: k = {k_pp:.4f} 1/m, terminal speed {vt_pp:.2f} m/s; bean bag: k = {k_bb:.4f} 1/m, v_t = {vt_bb:.1f} m/s")

    # ---- the drop: closed form against solve_ivp, and the experiment's prediction table
    heights = np.array([0.6, 0.9, 1.2, 1.5, 1.8, 2.0])
    t_vac, t_pp = predict_fall_times(heights, k_pp)
    _, t_bb = predict_fall_times(heights, k_bb)
    print("\nfall from rest: height, vacuum time, ping-pong time (drag), bean-bag time (drag), ping-pong excess (ms, %)")
    for h, tv, tp, tb in zip(heights, t_vac, t_pp, t_bb):
        print(f"  {h:4.1f} m   {tv:.4f} s   {tp:.4f} s   {tb:.4f} s   {1e3 * (tp - tv):5.1f} ms  {100 * (tp / tv - 1):4.1f} %")
    tr = simulate_projectile(0.0, 90.0, k_pp, z0=2.0, rtol=1e-12)
    print(f"solve_ivp drop of the ping-pong ball from 2.0 m: {tr.time_of_flight:.6f} s; closed form {vertical_drop_time(2.0, vt_pp):.6f} s; "
          f"difference {abs(tr.time_of_flight - vertical_drop_time(2.0, vt_pp)):.1e} s; impact speed {tr.impact_speed:.3f} m/s "
          f"(vacuum {math.sqrt(2 * G0 * 2.0):.3f} m/s, {100 * tr.impact_speed / vt_pp:.0f} % of v_t)")
    u30, u60 = frame_timing(1, 30.0)[1], frame_timing(1, 60.0)[1]
    excess = t_pp[-1] - t_vac[-1]
    print(f"frame uncertainty 0.41/fps: {1e3 * u30:.1f} ms at 30 fps, {1e3 * u60:.1f} ms at 60 fps; the 2.0 m drag excess {1e3 * excess:.1f} ms "
          f"is {excess / u30:.1f} and {excess / u60:.1f} of those")

    # ---- the toss
    v0, ang = 6.0, 45.0
    tv = simulate_projectile(v0, ang, 0.0)
    td = simulate_projectile(v0, ang, k_pp)
    print(f"\ntoss at {v0} m/s, {ang} deg: vacuum range {vacuum_range(v0, ang):.3f} m (code {tv.range:.3f}), time {vacuum_time_of_flight(v0, ang):.3f} s, "
          f"apex {v0 ** 2 * math.sin(math.radians(ang)) ** 2 / (2 * G0):.3f} m")
    print(f"  with drag (ping-pong): range {td.range:.3f} m ({100 * td.range / tv.range:.0f} % of vacuum), time {td.time_of_flight:.3f} s, "
          f"apex {td.apex_height:.3f} m, impact speed {td.impact_speed:.2f} m/s (launch {v0} m/s), impact angle "
          f"{math.degrees(math.atan2(-td.vz[-1], td.vx[-1])):.1f} deg below horizontal (launch {ang})")
    ia, ir = optimal_angle(v0, k_pp)
    print(f"  best angle with drag {ia:.1f} deg (range {ir:.3f} m); without drag 45 deg; at 45 deg the drag range is {td.range:.3f} m")
    print("  range against angle (vacuum / drag), m:")
    for a in (15, 25, 30, 35, 40, 45, 55, 65):
        print(f"    {a:2d} deg  {vacuum_range(v0, a):.3f}  {drag_range(v0, a, k_pp):.3f}")
    tb = simulate_projectile(v0, ang, k_bb)
    print(f"  bean bag at the same launch: range {tb.range:.3f} m ({100 * tb.range / tv.range:.0f} % of vacuum)")
    # convergence in the tolerance
    rs = [simulate_projectile(v0, ang, k_pp, rtol=r).range for r in (1e-4, 1e-6, 1e-8, 1e-10)]
    print(f"  range at rtol 1e-4, 1e-6, 1e-8, 1e-10: " + ", ".join(f"{r:.8f}" for r in rs) + f" m; spread {max(rs) - min(rs):.1e} m")

    # ---- the wind triangle
    gs, trk = wind_triangle(60.0, 90.0, 15.0, 30.0)
    hd, gs2 = heading_for_track(60.0, 90.0, 15.0, 30.0)
    print(f"\nwind triangle: TAS 60 m/s heading 090, wind 15 m/s from 030: ground speed {gs:.2f} m/s, track {trk:.1f} deg "
          f"(drift {trk - 90:.1f} deg); to make good 090 fly heading {hd:.1f} deg at {gs2:.2f} m/s ground speed")
    gs_h, _ = wind_triangle(60.0, 90.0, 15.0, 90.0)
    gs_t, _ = wind_triangle(60.0, 90.0, 15.0, 270.0)
    print(f"  pure headwind: {gs_h:.0f} m/s; pure tailwind: {gs_t:.0f} m/s; round trip 100 km each way in that wind takes "
          f"{(1e5 / gs_h + 1e5 / gs_t) / 60:.1f} min against {2e5 / 60 / 60:.1f} min in still air")

    # ---- SYNTHETIC data through the analysis pipeline (labelled synthetic; not a measurement)
    fps = 60.0
    hs = np.array([0.6, 0.9, 1.2, 1.5, 1.8])
    synth = make_synthetic_drop_data(k_pp, hs, fps, n_drops=5, seed=41)
    u_frame = frame_timing(1, fps)[1]
    means, us = zip(*[reduce_timings(s, u_frame) for s in synth])
    means, us = np.array(means), np.array(us)
    vt_fit, u_vt = fit_terminal_speed(hs, means, us)
    cd_fit = drag_coefficient_from_terminal_speed(vt_fit, PING_PONG["mass"], PING_PONG["diameter"])
    print(f"\nSYNTHETIC ping-pong drops at {fps:.0f} fps, 5 per height (seed 41), through the pipeline:")
    tvac, tdr = predict_fall_times(hs, k_pp)
    for h, m, u, a, b in zip(hs, means, us, tvac, tdr):
        print(f"  {h:3.1f} m: mean {m:.4f} +- {u:.4f} s; E_n against vacuum {normalized_error(m, u, a, 0.0):.2f}, against drag {normalized_error(m, u, b, 0.0):.2f}")
    print(f"  fitted v_t = {vt_fit:.2f} +- {u_vt:.2f} m/s (true {vt_pp:.2f}); C_D = {cd_fit:.2f} (true {PING_PONG['cd']})")
    synth30 = make_synthetic_drop_data(k_pp, hs, 30.0, n_drops=5, seed=41)
    m30, u30s = zip(*[reduce_timings(s, frame_timing(1, 30.0)[1]) for s in synth30])
    vt30, uvt30 = fit_terminal_speed(hs, np.array(m30), np.array(u30s))
    print(f"  the same at 30 fps: v_t = {vt30:.1f} +- {uvt30:.1f} m/s")

    # ---- figure
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    ax = axes[0]
    ax.plot(tv.x, tv.z, "k--", label="no drag")
    ax.plot(td.x, td.z, "C0", label="ping-pong ball")
    ax.plot(tb.x, tb.z, "C1", label="bean bag")
    ax.set_xlabel("x (m)"); ax.set_ylabel("z (m)"); ax.set_title(f"toss at {v0} m/s, {ang:.0f} deg"); ax.legend(); ax.set_aspect("equal")
    ax = axes[1]
    angs = np.arange(5.0, 86.0, 2.5)
    ax.plot(angs, [vacuum_range(v0, a) for a in angs], "k--", label="no drag")
    ax.plot(angs, [drag_range(v0, a, k_pp) for a in angs], "C0", label="ping-pong ball")
    ax.axvline(ia, color="C0", ls=":"); ax.axvline(45, color="k", ls=":")
    ax.set_xlabel("launch angle (deg)"); ax.set_ylabel("range (m)"); ax.set_title("range against angle"); ax.legend()
    ax = axes[2]
    hh = np.linspace(0.3, 2.2, 100)
    a_vac, a_pp = predict_fall_times(hh, k_pp)
    _, a_bb = predict_fall_times(hh, k_bb)
    ax.plot(hh, a_vac, "k--", label="vacuum")
    ax.plot(hh, a_pp, "C0", label="ping-pong (model)")
    ax.plot(hh, a_bb, "C1", label="bean bag (model)")
    ax.errorbar(hs, means, yerr=2 * us, fmt="o", color="C0", label="SYNTHETIC drops, k = 2")
    ax.set_xlabel("height (m)"); ax.set_ylabel("fall time (s)"); ax.set_title("drop: model and SYNTHETIC data"); ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "projectile_drag.png", dpi=120)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'projectile_drag.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
