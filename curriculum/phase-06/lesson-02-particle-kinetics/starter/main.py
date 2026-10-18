"""P06.02 Particle Kinetics: starter code.

How hard does a fighter pull at the bottom of a pull-up, what do the pilot and the wing feel, and how
much altitude does a dive recovery need when the g comes on at a finite rate? This module treats the
aircraft as a particle of mass m in the vertical plane and writes Newton's second law in path
(normal-tangential) coordinates, the frame that moves with the velocity:

    tangential:  m dV/dt      = T - D - W sin(gamma)
    normal:      m V dgamma/dt = L - W cos(gamma)

with V the speed, gamma the flight-path angle (positive climbing), W = m g. The load factor n = L / W
is what an accelerometer along the aircraft's normal axis reads in g. At the bottom of a circular
pull-up (gamma = 0, constant V, radius R) the normal equation gives n = 1 + V^2 / (g R).

The independent checks: the closed-form altitude loss and time of a constant-n, constant-speed dive
recovery against `solve_ivp`; the load factor recovered from the integrated *positions alone*
(differentiated twice, then projected on the path normal) against the commanded load factor; a
constant-speed loop commanded with n = cos(gamma) + V^2/(g R) staying on a circle of radius R; and,
with thrust equal to drag, conservation of the specific energy V^2/2 + g h.

You write every function whose body is `raise NotImplementedError`. `ramp_load_factor`,
`recovery_altitude_loss` and `main` are given.

Frame: x horizontal (direction of flight at gamma = 0), h up; the normal unit vector
e_n = (-sin gamma, cos gamma) points to the centre of curvature of a pull-up. Units: SI.

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes pullout.png)
Test: pytest curriculum/phase-06/lesson-02-particle-kinetics
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

G0 = 9.80665
M_AIRCRAFT = 12000.0      # kg, illustrative fighter mass in combat configuration
M_PILOT = 100.0           # kg, illustrative pilot plus helmet, suit and seat-pack share
N_LIMIT = 9.0             # illustrative structural and human limit load factor


# ---------------------------------------------------------------- free-body results in path coordinates

def loop_load_factor(V, R, gamma=0.0):
    """Load factor on a circular vertical path of radius R flown at constant speed V:
    the normal equation L - W cos(gamma) = m V^2 / R gives n = cos(gamma) + V^2 / (g R).
    At the bottom (gamma = 0) n = 1 + V^2/(g R); at the top (gamma = pi) n = V^2/(g R) - 1."""
    raise NotImplementedError


def radius_for_load_factor(V, n, gamma=0.0):
    """Radius of curvature of the flight path at load factor n, speed V and flight-path angle gamma:
    R = V^2 / (g (n - cos gamma)); at the bottom of a pull-up R = V^2 / (g (n - 1)). Infinite when the
    lift only balances the weight component (a straight path)."""
    raise NotImplementedError


def level_turn(V: float, n: float) -> tuple[float, float, float]:
    """Sustained level turn at load factor n and speed V. Vertical balance L cos(phi) = W gives the bank
    angle phi = arccos(1/n); the horizontal component L sin(phi) = m V^2 / R gives
    R = V^2 / (g sqrt(n^2 - 1)) and the turn rate omega = V / R. Returns (phi, R, omega)."""
    raise NotImplementedError


def apparent_weight(m: float, n) -> np.ndarray:
    """Force between a body of mass m and its support (the seat) when the aircraft pulls n g along its
    normal axis: n m g (N). It is also the lift the wing of an aircraft of mass m must produce."""
    raise NotImplementedError


# ---------------------------------------------------------------- the dive recovery in closed form

def pullout_altitude_loss(V: float, n: float, gamma0: float) -> float:
    """Altitude lost recovering from a dive at flight-path angle gamma0 < 0 to level flight at constant
    speed V and constant load factor n > 1, from dh/dgamma = V^2 sin(gamma) / (g (n - cos gamma)):

        delta_h = (V^2 / g) ln((n - cos gamma0) / (n - 1))     (positive, m)."""
    raise NotImplementedError


def pullout_time(V: float, n: float, gamma0: float) -> float:
    """Time for the same recovery, from dt = V dgamma / (g (n - cos gamma)). With
    F(gamma) = 2 / sqrt(n^2 - 1) atan(sqrt((n + 1)/(n - 1)) tan(gamma / 2)),
    t = (V / g) (F(0) - F(gamma0)), valid for |gamma0| < pi."""
    raise NotImplementedError


# ---------------------------------------------------------------- the point-mass flight path

def ramp_load_factor(n0: float, n_max: float, ndot: float) -> Callable[[float, np.ndarray], float]:
    """Commanded load factor rising from n0 at the g-onset rate ndot (g/s) and holding at n_max:
    n(t) = min(n0 + ndot t, n_max). ndot = inf applies n_max at once. Returns n_cmd(t, y)."""
    if np.isinf(ndot):
        return lambda t, y: n_max
    return lambda t, y: min(n0 + ndot * t, n_max)


def flight_path_rhs(t: float, y: np.ndarray, n_cmd: Callable, speed_hold: bool = True,
                    excess_thrust: float = 0.0) -> list[float]:
    """Point-mass equations of motion in path coordinates, state y = [x, h, V, gamma]:

        x' = V cos(gamma),  h' = V sin(gamma),
        V' = g ((T - D)/W - sin(gamma)),   gamma' = g (n - cos(gamma)) / V.

    With ``speed_hold`` the throttle and speed brake hold V constant (V' = 0); otherwise
    (T - D)/W = ``excess_thrust`` (0: thrust equals drag, so V^2/2 + g h is conserved).
    ``n_cmd(t, y)`` is the commanded load factor."""
    raise NotImplementedError


def simulate_pullout(V0: float, gamma0: float, n_cmd: Callable, speed_hold: bool = True,
                     excess_thrust: float = 0.0, h0: float = 3000.0, t_max: float = 60.0,
                     stop_level: bool = True, rtol: float = 1e-10, atol: float = 1e-9, max_step: float = np.inf):
    """Integrate ``flight_path_rhs`` from (0, h0, V0, gamma0) with `solve_ivp` (RK45, dense output).
    With ``stop_level`` the run ends when gamma rises through 0 (level flight). Returns the
    `solve_ivp` result; ``sol.t_events[0]`` holds the level-off time."""
    raise NotImplementedError


def kinematic_load_factor(t: np.ndarray, x: np.ndarray, h: np.ndarray) -> np.ndarray:
    """Load factor recovered from positions alone, the way an observer with a theodolite would do it:
    differentiate x(t), h(t) twice (``np.gradient``, second order), form the specific force
    f = a + g e_h (what an accelerometer measures) and project it on the path normal
    e_n = (-sin gamma, cos gamma) with gamma = atan2(h', x'):  n = f . e_n / g."""
    raise NotImplementedError


def recovery_altitude_loss(V0: float, gamma0: float, n_max: float, ndot: float, speed_hold: bool = True) -> float:
    """Altitude lost from the start of the pull (n = 1) to level flight with a finite g-onset rate,
    by simulation (m, positive)."""
    sol = simulate_pullout(V0, gamma0, ramp_load_factor(1.0, n_max, ndot), speed_hold=speed_hold, h0=0.0)
    if sol.status != 1:
        raise RuntimeError("the aircraft did not level off within t_max")
    return float(-sol.y_events[0][0][1])


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save pullout.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    V = 250.0
    W = M_AIRCRAFT * G0
    print(f"aircraft {M_AIRCRAFT:.0f} kg (W = {W / 1e3:.2f} kN), pilot with gear {M_PILOT:.0f} kg, V = {V:.0f} m/s")
    R9 = float(radius_for_load_factor(V, N_LIMIT))
    print(f"\npull-up: radius for n = {N_LIMIT:.0f} at the bottom R = V^2/(g(n-1)) = {R9:.1f} m")
    for R in (800.0, 1000.0, 2000.0, 4000.0):
        print(f"  R = {R:6.0f} m: n bottom {float(loop_load_factor(V, R)):.3f}, top (same V) {float(loop_load_factor(V, R, np.pi)):.3f}")
    print(f"  at n = {N_LIMIT:.0f}: seat force {float(apparent_weight(M_PILOT, N_LIMIT)) / 1e3:.2f} kN, "
          f"lift {float(apparent_weight(M_AIRCRAFT, N_LIMIT)) / 1e3:.1f} kN, each wing {float(apparent_weight(M_AIRCRAFT, N_LIMIT)) / 2e3:.1f} kN")
    print(f"  at n = 7 (the recovery below): seat force {float(apparent_weight(M_PILOT, 7.0)) / 1e3:.2f} kN, "
          f"lift {float(apparent_weight(M_AIRCRAFT, 7.0)) / 1e3:.1f} kN")
    phi, Rt, om = level_turn(V, N_LIMIT)
    print(f"  level turn at the same n: bank {np.degrees(phi):.2f} deg, radius {Rt:.1f} m, rate {np.degrees(om):.2f} deg/s")

    # dive recovery: closed form against solve_ivp, instant g
    g0 = np.radians(-60.0)
    n = 7.0
    dh = pullout_altitude_loss(V, n, g0)
    tp = pullout_time(V, n, g0)
    sol = simulate_pullout(V, g0, ramp_load_factor(n, n, np.inf), h0=0.0)
    dh_sim, tp_sim = -sol.y_events[0][0][1], sol.t_events[0][0]
    print(f"\ndive recovery from {np.degrees(g0):.0f} deg at {V:.0f} m/s and n = {n:.0f}, constant speed, g applied at once:")
    print(f"  closed form: altitude loss {dh:.3f} m in {tp:.4f} s")
    print(f"  solve_ivp:   altitude loss {dh_sim:.3f} m in {tp_sim:.4f} s   (differences {dh_sim - dh:.1e} m, {tp_sim - tp:.1e} s)")

    # load factor from positions alone
    t = np.linspace(0.0, tp_sim, 2001)
    yy = sol.sol(t)
    n_kin = kinematic_load_factor(t, yy[0], yy[1])
    print(f"  load factor from positions alone: mean {np.mean(n_kin):.6f}, max |n_kin - n| = {np.max(np.abs(n_kin - n)):.2e}")

    # finite g-onset
    print(f"\ng-onset rate: altitude lost from {np.degrees(g0):.0f} deg at {V:.0f} m/s pulling to {n:.0f} g (constant speed)")
    print("   onset (g/s)   altitude loss (m)   extra (m)")
    rates = [1.0, 2.0, 4.0, 6.0, 10.0, np.inf]
    losses = []
    for r in rates:
        losses.append(recovery_altitude_loss(V, g0, n, r))
        print(f"   {r:9.0f}     {losses[-1]:9.1f}        {losses[-1] - dh:7.1f}")
    print(f"   estimate for 6 g/s: the ramp takes (n-1)/ndot = {(n - 1) / 6:.2f} s, roughly a delay of half that, "
          f"V sin|gamma0| x {(n - 1) / 12:.2f} s = {V * np.sin(-g0) * (n - 1) / 12:.1f} m")
    dh_e = recovery_altitude_loss(V, g0, n, 6.0, speed_hold=False)
    print(f"   thrust = drag instead of constant speed (6 g/s): altitude loss {dh_e:.1f} m")
    sol_e = simulate_pullout(V, g0, ramp_load_factor(1.0, n, 6.0), speed_hold=False, h0=0.0)
    E = 0.5 * sol_e.y[2] ** 2 + G0 * sol_e.y[1]
    print(f"   specific energy drift V^2/2 + g h: {np.max(np.abs(E - E[0])) / E[0]:.1e} relative; "
          f"speed at level-off {sol_e.y_events[0][0][2]:.1f} m/s")

    # a constant-speed loop with n = cos(gamma) + V^2/(gR)
    R = 1500.0
    loop = lambda tt, y: float(loop_load_factor(V, R, y[3]))
    sl = simulate_pullout(V, 0.0, loop, h0=1000.0, t_max=2 * np.pi * R / V, stop_level=False)
    rad = np.hypot(sl.y[0], sl.y[1] - (1000.0 + R))
    print(f"\nconstant-speed loop, R = {R:.0f} m: n from {float(loop_load_factor(V, R, np.pi)):.3f} (top) "
          f"to {float(loop_load_factor(V, R)):.3f} (bottom); max |radius - R| over the loop {np.max(np.abs(rad - R)):.1e} m")

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    for r, c in zip([1.0, 4.0, np.inf], ["C0", "C1", "C2"]):
        s = simulate_pullout(V, g0, ramp_load_factor(1.0 if np.isfinite(r) else n, n, r), h0=0.0)
        tt = np.linspace(0, s.t[-1], 400)
        y = s.sol(tt)
        lab = "instant" if np.isinf(r) else f"{r:.0f} g/s"
        ax[0].plot(y[0], y[1], c, label=lab)
        ax[1].plot(tt, [ramp_load_factor(1.0 if np.isfinite(r) else n, n, r)(q, None) for q in tt], c, label=lab)
    ax[1].plot(t, n_kin, "k:", label="from positions (instant)")
    ax[0].set_xlabel("x (m)")
    ax[0].set_ylabel("h (m)")
    ax[0].set_aspect("equal")
    ax[0].set_title("60 deg dive recovery at 7 g, 250 m/s")
    ax[0].legend()
    ax[1].set_xlabel("t (s)")
    ax[1].set_ylabel("load factor n")
    ax[1].set_title("commanded and recovered n")
    ax[1].legend()
    fr = [r for r in rates if np.isfinite(r)]
    ax[2].plot(fr, losses[:-1], "o-", label="simulation")
    ax[2].axhline(dh, color="k", ls="--", lw=0.8, label="closed form, instant g")
    ax[2].set_xscale("log")
    ax[2].set_xlabel("g-onset rate (g/s)")
    ax[2].set_ylabel("altitude loss (m)")
    ax[2].set_title("cost of a slow g-onset")
    ax[2].legend()
    fig.tight_layout()
    fig.savefig(out_dir / "pullout.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'pullout.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
