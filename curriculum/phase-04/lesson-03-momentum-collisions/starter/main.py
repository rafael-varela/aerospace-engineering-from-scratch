"""P04.03 Momentum & Collisions: starter code.

Linear momentum, impulse and collisions in one dimension, then the variable-mass system: a rocket
treated first as a cart that throws pellets backwards one at a time (pure momentum conservation, no
calculus), then in the limit of infinitely many infinitesimal pellets, which is the Tsiolkovsky
rocket equation. The continuous rocket is integrated with `scipy.integrate.solve_ivp` (vertical
ascent with gravity and optional drag), checked against the closed-form burnout velocity and
altitude, and checked for conservation of the momentum of rocket plus exhaust. Staging follows from
the same equation.

Conventions: one dimension, positive up (or forward). Exhaust speed u > 0 is measured relative to the
rocket *after* each ejection. Mass flow mdot > 0 is the rate at which the rocket loses mass.

Run:  python main.py [out_dir]     (prints the tables, writes rocket_momentum.png)
Test: pytest curriculum/phase-04/lesson-03-momentum-collisions
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

G0 = 9.80665     # m/s^2, standard gravity (also defines specific impulse: u = Isp g0)


# ---------------------------------------------------------------- momentum, impulse, collisions

def impulse(t: np.ndarray, force: np.ndarray) -> float:
    """J = integral of F dt [N s] from a sampled force history, by the trapezoid rule."""
    raise NotImplementedError


def centre_of_mass(masses: np.ndarray, positions: np.ndarray) -> np.ndarray:
    """x_cm = sum(m_i x_i) / sum(m_i). `positions` has shape (n,) or (n, d); returns a scalar array or (d,)."""
    raise NotImplementedError


def collide_1d(m1: float, v1: float, m2: float, v2: float, e: float = 1.0) -> tuple[float, float]:
    """Velocities after a head-on collision with coefficient of restitution e (0 <= e <= 1):
    momentum is conserved, m1 v1 + m2 v2 = m1 v1' + m2 v2', and the separation speed is e times the
    approach speed, v2' - v1' = -e (v2 - v1). e = 1: elastic; e = 0: perfectly inelastic (the bodies stick).
    """
    raise NotImplementedError


def collision_energy_loss(m1: float, v1: float, m2: float, v2: float, e: float) -> float:
    """Kinetic energy dissipated in the collision: (1/2) mu (1 - e^2) (v1 - v2)^2 with the reduced mass
    mu = m1 m2 / (m1 + m2). Zero for e = 1; for e = 0 it is all the kinetic energy of the relative motion."""
    raise NotImplementedError


# ---------------------------------------------------------------- the pellet rocket and Tsiolkovsky

def pellet_rocket_velocities(m0: float, mf: float, u: float, n: int) -> np.ndarray:
    """Velocity after each of n equal ejections that take the mass from m0 to mf, starting at rest.

    Ejection i throws a pellet of mass dm = (m0 - mf)/n backward at speed u relative to the rocket
    after the throw. Momentum before = momentum after in the frame where the rocket moved at v:
        m v = (m - dm)(v + dv) + dm (v + dv - u)   =>   dv = u dm / m,
    with m the mass *before* the throw. Returns an array of n + 1 velocities (v_0 = 0).
    """
    raise NotImplementedError


def tsiolkovsky_dv(u: float, m0: float, mf: float) -> float:
    """Ideal velocity change of a rocket burning from m0 to mf with exhaust speed u: u ln(m0/mf)."""
    raise NotImplementedError


def mass_ratio_for_dv(dv: float, u: float) -> float:
    """m0/mf needed for a velocity change dv: exp(dv/u)."""
    raise NotImplementedError


def pellet_convergence(m0: float, mf: float, u: float, ns: list[int]) -> np.ndarray:
    """Error of the n-pellet rocket's final velocity against Tsiolkovsky, u ln(m0/mf) - v_n, for each n.
    Positive (the pellet rocket is always slower) and O(1/n)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the continuous rocket

@dataclass
class RocketFlight:
    t: np.ndarray
    z: np.ndarray
    v: np.ndarray
    m: np.ndarray
    p_exhaust: np.ndarray
    burn_time: float


def rocket_rhs(t: float, y: np.ndarray, mdot: float, u: float, g: float, k_drag: float,
               scale_height: float = 8500.0) -> np.ndarray:
    """Vertical rocket with state y = [z, v, m, P_ex]:

        dz/dt = v,
        dv/dt = (mdot u - m g - k_drag exp(-z/H) v |v|) / m,
        dm/dt = -mdot,
        dP_ex/dt = mdot (v - u)       (momentum carried away by the exhaust, in the ground frame).

    The thrust mdot u is the momentum flux of the exhaust relative to the rocket. The drag is
    k_drag exp(-z/H) v |v| with k_drag = rho0 C_D A / 2 [kg/m] and an exponential atmosphere of scale
    height H (a crude model, for the size of the drag loss only; H = inf gives constant density).
    """
    raise NotImplementedError


def simulate_rocket(m0: float, mf: float, mdot: float, u: float, g: float = G0, k_drag: float = 0.0,
                    v0: float = 0.0, rtol: float = 1e-11, n_samples: int = 501, scale_height: float = 8500.0) -> RocketFlight:
    """Integrate `rocket_rhs` from (z, v, m, P_ex) = (0, v0, m0, 0) to burnout at t_b = (m0 - mf)/mdot with
    `solve_ivp` (DOP853, atol 1e-9), sampled at n_samples times. Raises ValueError if the initial
    thrust-to-weight ratio mdot u/(m0 g) is not above 1 while starting from rest under gravity (the rocket
    would sit on the pad)."""
    raise NotImplementedError


def burnout_closed_form(m0: float, mf: float, mdot: float, u: float, g: float = G0) -> tuple[float, float]:
    """Burnout velocity and altitude of a vertical, drag-free rocket from rest with constant mdot and g:

        v_b = u ln(m0/mf) - g t_b,
        z_b = u t_b - u (mf/mdot) ln(m0/mf) - g t_b^2 / 2,   t_b = (m0 - mf)/mdot.
    """
    raise NotImplementedError


def velocity_losses(flight: RocketFlight, u: float, g: float = G0) -> tuple[float, float, float]:
    """Split the ideal velocity change into what was achieved and what was lost:
    returns (ideal dv = u ln(m0/mf), gravity loss = g t_b, drag loss = ideal - gravity loss - (v_b - v_0))."""
    raise NotImplementedError


# ---------------------------------------------------------------- staging

def stack_dv(payload: float, stages: list[tuple[float, float, float]]) -> tuple[float, list[float]]:
    """Ideal velocity change of a stack of stages fired in order (first in the list fires first).
    Each stage is (propellant mass, structural mass, exhaust speed u). Stage i burns with everything above
    it (later stages and the payload) as dead weight, and drops its structure at burnout.
    Returns (total dv, list of per-stage dv)."""
    raise NotImplementedError


def structure_mass(m_prop: float, eps: float) -> float:
    """Structural mass of a stage with structural coefficient eps = m_s/(m_s + m_p): m_s = eps m_p/(1 - eps)."""
    raise NotImplementedError


def best_two_stage_split(payload: float, m_prop_total: float, eps: float, u: float) -> tuple[float, float]:
    """Fraction of the total propellant to put in the first stage that maximises the two-stage dv (both
    stages with structural coefficient eps and exhaust speed u), by a bounded scalar search on [0.01, 0.99].
    Returns (fraction, dv)."""
    raise NotImplementedError


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save rocket_momentum.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)

    # ---- docking: a perfectly inelastic collision
    m_ship, m_station, v_rel = 12_000.0, 420_000.0, 0.10
    v1, v2 = collide_1d(m_ship, v_rel, m_station, 0.0, e=0.0)
    print(f"docking: {m_ship / 1e3:.0f} t vehicle at {v_rel} m/s into a {m_station / 1e3:.0f} t station: common velocity {1e3 * v1:.3f} mm/s; "
          f"energy absorbed by the docking mechanism {collision_energy_loss(m_ship, v_rel, m_station, 0.0, 0.0):.2f} J "
          f"(of {0.5 * m_ship * v_rel ** 2:.2f} J)")
    a, b = collide_1d(1.0, 2.0, 3.0, 0.0, e=1.0)
    print(f"elastic, 1 kg at 2 m/s on 3 kg at rest: {a:.2f} and {b:.2f} m/s; with e = 0.5: {collide_1d(1.0, 2.0, 3.0, 0.0, 0.5)}")

    # ---- the pellet rocket
    u = 300.0 * G0                       # Isp 300 s
    m0, mf = 10_000.0, 2_500.0           # mass ratio 4
    exact = tsiolkovsky_dv(u, m0, mf)
    print(f"\npellet rocket, m0 = {m0:.0f} kg, mf = {mf:.0f} kg (R = {m0 / mf:.0f}), u = {u:.1f} m/s (Isp 300 s); Tsiolkovsky {exact:.2f} m/s")
    ns = [1, 2, 4, 10, 100, 1000, 10000]
    errs = pellet_convergence(m0, mf, u, ns)
    for n, e in zip(ns, errs):
        print(f"  n = {n:5d}: v = {exact - e:8.2f} m/s, error {e:8.3f} m/s, n x error {n * e:8.2f}")
    print(f"  one pellet: u (m0 - mf)/m0 = {u * (m0 - mf) / m0:.2f} m/s; leading error term u (m0 - mf)(1/mf - 1/m0)/(2n) = "
          f"{u * (m0 - mf) * (1 / mf - 1 / m0) / 2:.2f}/n m/s")

    # ---- the continuous rocket: drag-free vertical ascent against the closed form
    mdot = 50.0                          # kg/s: T = 147 kN, T/W0 = 1.50
    T = mdot * u
    fl = simulate_rocket(m0, mf, mdot, u)
    vb, zb = burnout_closed_form(m0, mf, mdot, u)
    print(f"\nvertical ascent, mdot = {mdot} kg/s, thrust {T / 1e3:.1f} kN, T/W0 = {T / (m0 * G0):.2f}, burn {fl.burn_time:.0f} s:")
    print(f"  burnout v = {fl.v[-1]:.3f} m/s (closed form {vb:.3f}), z = {fl.z[-1] / 1e3:.4f} km (closed form {zb / 1e3:.4f})")
    ideal, grav, drag_loss = velocity_losses(fl, u)
    print(f"  ideal {ideal:.1f} m/s, gravity loss {grav:.1f} m/s ({100 * grav / ideal:.0f} %), drag loss {drag_loss:.2e} m/s")
    fs = simulate_rocket(m0, mf, mdot, u, g=0.0)
    P = fs.m * fs.v + fs.p_exhaust
    print(f"  in free space: burnout {fs.v[-1]:.3f} m/s; total momentum (rocket + exhaust) drifts by {np.ptp(P):.1e} kg m/s "
          f"(rocket momentum alone reaches {fs.m[-1] * fs.v[-1] / 1e3:.1f} t m/s)")
    k = 0.5 * 1.225 * 0.3 * math.pi * 0.5 ** 2     # rho0 1.225, CD 0.3, diameter 1 m
    fd = simulate_rocket(m0, mf, mdot, u, k_drag=k)
    i_q = int(np.argmax(k * np.exp(-fd.z / 8500.0) * fd.v ** 2))
    print(f"  with drag in an exponential atmosphere (k = {k:.3f} kg/m, H = 8.5 km): burnout {fd.v[-1]:.1f} m/s, drag loss "
          f"{velocity_losses(fd, u)[2]:.1f} m/s; maximum dynamic pressure {k * math.exp(-fd.z[i_q] / 8500.0) * fd.v[i_q] ** 2 / (0.3 * math.pi * 0.25) / 1e3:.1f} kPa "
          f"at t = {fd.t[i_q]:.0f} s, z = {fd.z[i_q] / 1e3:.1f} km")

    # ---- gravity loss against thrust-to-weight ratio
    print("\ngravity loss against initial thrust-to-weight ratio (same propellant, drag-free):")
    tw_list = [1.2, 1.5, 2.0, 3.0, 5.0]
    vbs = []
    for tw in tw_list:
        md = tw * m0 * G0 / u
        f2 = simulate_rocket(m0, mf, md, u, n_samples=11)
        vbs.append(f2.v[-1])
        print(f"  T/W0 = {tw:3.1f}: burn {f2.burn_time:6.1f} s, gravity loss {G0 * f2.burn_time:6.1f} m/s, burnout {f2.v[-1]:7.1f} m/s")

    # ---- staging
    payload, mprop, eps = 500.0, 20_000.0, 0.10
    ms = structure_mass(mprop, eps)
    dv1, _ = stack_dv(payload, [(mprop, ms, u)])
    f_half = stack_dv(payload, [(mprop / 2, structure_mass(mprop / 2, eps), u), (mprop / 2, structure_mass(mprop / 2, eps), u)])
    fbest, dvbest = best_two_stage_split(payload, mprop, eps, u)
    print(f"\nstaging: payload {payload:.0f} kg, {mprop / 1e3:.0f} t of propellant, structural coefficient {eps}, u = {u:.0f} m/s")
    print(f"  single stage: {dv1:.0f} m/s; two equal stages: {f_half[0]:.0f} m/s ({f_half[1][0]:.0f} + {f_half[1][1]:.0f}); "
          f"best split {100 * fbest:.1f} % / {100 * (1 - fbest):.1f} %: {dvbest:.0f} m/s; limit of the single stage as payload -> 0: "
          f"{u * math.log(1 / eps):.0f} m/s")
    print(f"  mass ratio for 9.4 km/s at u = {u:.0f} m/s: {mass_ratio_for_dv(9400.0, u):.1f}; at u = 4400 m/s (hydrogen): {mass_ratio_for_dv(9400.0, 4400.0):.1f}")

    # ---- figure
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    ax = axes[0]
    nn = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024])
    ax.loglog(nn, pellet_convergence(m0, mf, u, list(nn)), "o-", label="u ln(m0/mf) - v_n")
    ax.loglog(nn, u * (m0 - mf) * (1 / mf - 1 / m0) / 2 / nn, "k--", label="leading term / n")
    ax.set_xlabel("number of pellets n"); ax.set_ylabel("velocity error (m/s)"); ax.set_title("pellets -> Tsiolkovsky"); ax.legend()
    ax = axes[1]
    ax.plot(fl.t, fl.v, label="gravity, no drag")
    ax.plot(fs.t, fs.v, label="free space")
    ax.plot(fd.t, fd.v, label="gravity and drag")
    ax.set_xlabel("t (s)"); ax.set_ylabel("v (m/s)"); ax.set_title("vertical ascent"); ax.legend()
    ax = axes[2]
    fr = np.linspace(0.05, 0.95, 91)
    dv2 = [stack_dv(payload, [(f * mprop, structure_mass(f * mprop, eps), u), ((1 - f) * mprop, structure_mass((1 - f) * mprop, eps), u)])[0] for f in fr]
    ax.plot(100 * fr, dv2, label="two stages")
    ax.axhline(dv1, color="k", ls="--", label="one stage")
    ax.set_xlabel("propellant in stage 1 (%)"); ax.set_ylabel("dv (m/s)"); ax.set_title("staging"); ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "rocket_momentum.png", dpi=120)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'rocket_momentum.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
