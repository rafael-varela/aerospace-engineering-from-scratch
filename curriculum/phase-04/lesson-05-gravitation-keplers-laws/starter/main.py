"""P04.05 Gravitation & Kepler's Laws: starter code.

Complete every function marked NotImplementedError. When finished, running this file produces the
lesson's figure and prints the verification tables.

Newton's law of gravitation applied to satellites and to the general two-body problem:

* the inverse-square acceleration, g at altitude, circular speed, period, escape speed;
* specific orbital energy and angular momentum, the vis-viva relation, orbital elements (a, e)
  recovered from a state vector;
* a hand-written leapfrog (kick-drift-kick) integrator, kept here because the point of the
  lesson is conservation: its energy error stays bounded and its angular momentum is exact to
  rounding, which neither RK4 nor an adaptive library solver gives you;
* the full two-body problem (two masses, six first-order equations per body) reduced to the
  one-body problem about the centre of mass;
* Kepler's three laws checked on the simulated orbit; a Hohmann transfer from a low orbit to
  geostationary orbit as the engineering payoff.

Oracle: ``scipy.integrate.solve_ivp`` (DOP853 at tight tolerance) and the closed-form Kepler
relations.

Units: SI throughout (m, s, kg); mu = G M in m^3/s^2.

Run:  python main.py          (prints the verification tables, writes two_body.png)
Test: pytest curriculum/phase-04/lesson-05-gravitation-keplers-laws
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

G_NEWTON = 6.674e-11                 # m^3 kg^-1 s^-2
MU_EARTH = 3.986004418e14            # m^3/s^2
R_EARTH = 6.371e6                    # m (mean radius)
M_EARTH = 5.972e24                   # kg
M_MOON = 7.346e22                    # kg
D_EARTH_MOON = 3.844e8               # m (mean distance)
DAY_SIDEREAL = 86164.0905            # s


# ---------------------------------------------------------------- the law and its simplest consequences

def gravity_acceleration(r: np.ndarray, mu: float = MU_EARTH) -> np.ndarray:
    """Acceleration -mu r / |r|^3 of a test mass at position(s) ``r`` relative to a point mass.

    ``r`` may be one vector (shape (d,)) or a stack of vectors (shape (n, d)).
    """
    raise NotImplementedError


def gravitational_force(m1: float, m2: float, r1: np.ndarray, r2: np.ndarray, G: float = G_NEWTON) -> np.ndarray:
    """Force on mass 2 (at r2) due to mass 1 (at r1): attractive, inverse square, along the line of centres."""
    raise NotImplementedError


def g_at_altitude(h: float | np.ndarray, mu: float = MU_EARTH, R: float = R_EARTH) -> float | np.ndarray:
    """Gravitational field strength mu / (R + h)^2 [m/s^2] at altitude h above a sphere of radius R."""
    raise NotImplementedError


def circular_speed(r: float, mu: float = MU_EARTH) -> float:
    """Speed sqrt(mu / r) of a circular orbit of radius r."""
    raise NotImplementedError


def orbital_period(a: float, mu: float = MU_EARTH) -> float:
    """Kepler's third law: T = 2 pi sqrt(a^3 / mu) for semi-major axis a."""
    raise NotImplementedError


def escape_speed(r: float, mu: float = MU_EARTH) -> float:
    """Speed sqrt(2 mu / r) at which the specific energy is zero (parabolic escape)."""
    raise NotImplementedError


def vis_viva(r: float | np.ndarray, a: float, mu: float = MU_EARTH) -> float | np.ndarray:
    """Vis-viva: v^2 = mu (2/r - 1/a). Returns the speed v at radius r on an orbit of semi-major axis a."""
    raise NotImplementedError


def specific_energy(r: np.ndarray, v: np.ndarray, mu: float = MU_EARTH) -> float | np.ndarray:
    """Specific mechanical energy v^2/2 - mu/|r| [J/kg]; accepts single vectors or stacks of shape (n, d)."""
    raise NotImplementedError


def specific_angular_momentum(r: np.ndarray, v: np.ndarray) -> float | np.ndarray:
    """Out-of-plane component x vy - y vx of r x v for planar states [m^2/s]; accepts stacks of shape (n, 2)."""
    raise NotImplementedError


def orbital_elements(r: np.ndarray, v: np.ndarray, mu: float = MU_EARTH) -> tuple[float, float, float, float]:
    """Semi-major axis a, eccentricity e, periapsis radius r_p and apoapsis radius r_a of the orbit
    through the planar state (r, v). Uses eps = -mu / 2a and h^2 = mu a (1 - e^2).

    For an unbound state (eps >= 0) a is negative or infinite and r_a is +inf.
    """
    raise NotImplementedError


def perigee_state(a: float, e: float, mu: float = MU_EARTH) -> tuple[np.ndarray, np.ndarray]:
    """Planar position and velocity at periapsis of an orbit (a, e): r = a(1 - e) on the x axis, v along +y."""
    raise NotImplementedError


# ---------------------------------------------------------------- integrators

def two_body_rhs(t: float, y: np.ndarray, mu: float = MU_EARTH) -> np.ndarray:
    """Right-hand side for solve_ivp with state y = [x, y, vx, vy] of a test mass about a point mass."""
    raise NotImplementedError


def leapfrog(accel: Callable[[np.ndarray], np.ndarray], r0: np.ndarray, v0: np.ndarray, h: float,
             n_steps: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Kick-drift-kick leapfrog (velocity Verlet) with fixed step h.

    Returns (t, R, V) with R and V of shape (n_steps + 1, d). One acceleration evaluation per step
    (the kick at the end of one step is the kick at the start of the next). Second order, symplectic,
    time-reversible: the energy error oscillates but does not grow, and for a central force the
    angular momentum is conserved to rounding.
    """
    raise NotImplementedError


def propagate_leapfrog(a: float, e: float, n_orbits: float, steps_per_orbit: int,
                       mu: float = MU_EARTH) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Propagate the orbit (a, e) from periapsis for n_orbits with the leapfrog at a fixed step T/steps_per_orbit."""
    raise NotImplementedError


def propagate_solve_ivp(a: float, e: float, n_orbits: float, mu: float = MU_EARTH, rtol: float = 1e-12,
                        n_out: int = 2001) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Oracle: the same orbit with scipy's DOP853 at tight tolerance, sampled at n_out times."""
    raise NotImplementedError


def position_error(R: np.ndarray, R_ref: np.ndarray) -> float:
    """Largest distance between two sampled trajectories of the same shape [m]."""
    raise NotImplementedError


def convergence_order(hs: np.ndarray, errs: np.ndarray) -> float:
    """Least-squares slope of log(err) against log(h)."""
    raise NotImplementedError


def energy_drift(R: np.ndarray, V: np.ndarray, mu: float = MU_EARTH) -> np.ndarray:
    """Relative energy error (eps(t) - eps(0)) / |eps(0)| along a trajectory."""
    raise NotImplementedError


# ---------------------------------------------------------------- Kepler's laws on the simulated orbit

def areal_velocity(R: np.ndarray, dt: float) -> np.ndarray:
    """Area swept per unit time between successive samples, |r_k x r_{k+1}| / (2 dt) (Kepler's second law
    says this is constant and equal to |h| / 2). Returns an array of length n - 1."""
    raise NotImplementedError


def period_from_trajectory(t: np.ndarray, R: np.ndarray) -> float:
    """Period from the times of successive periapsis passages (minima of |r|), found by parabolic
    interpolation through the three samples around each minimum. Needs at least two full orbits."""
    raise NotImplementedError


def keplers_third_law_ratio(a: float | np.ndarray, T: float | np.ndarray) -> float | np.ndarray:
    """T^2 / a^3, which Kepler's third law says is the same for every orbit about the same primary (4 pi^2 / mu)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the full two-body problem

def two_body_full_rhs(t: float, y: np.ndarray, m1: float, m2: float, G: float = G_NEWTON) -> np.ndarray:
    """Both bodies move: y = [x1, y1, x2, y2, vx1, vy1, vx2, vy2]. Each feels the other's gravity."""
    raise NotImplementedError


def simulate_two_body(m1: float, m2: float, r_rel0: np.ndarray, v_rel0: np.ndarray, t_end: float,
                      n_out: int = 1001, G: float = G_NEWTON, rtol: float = 1e-11) -> dict:
    """Integrate two free masses started at rest in their centre-of-mass frame with the given
    relative position and velocity. Returns a dict with t, r1, r2, v1, v2 (shape (n, 2)), the
    centre of mass, the total momentum and the relative motion r2 - r1."""
    raise NotImplementedError


# ---------------------------------------------------------------- the engineering payoff: a transfer

def hohmann_transfer(r1: float, r2: float, mu: float = MU_EARTH) -> tuple[float, float, float]:
    """Hohmann transfer between circular orbits of radii r1 and r2: (dv1, dv2, transfer time).

    Both burns are positive magnitudes; the transfer ellipse has a = (r1 + r2)/2 and the time of
    flight is half its period. Works for r2 < r1 as well (both burns then slow the vehicle down).
    """
    raise NotImplementedError


# ---------------------------------------------------------------- driver

def main(out_dir: str | Path = ".") -> None:
    """Print the verification tables and write two_body.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    h_iss = 400e3
    r_iss = R_EARTH + h_iss
    print(f"g at the surface {g_at_altitude(0.0):.3f} m/s^2, at {h_iss/1e3:.0f} km {g_at_altitude(h_iss):.3f} m/s^2 "
          f"({100 * g_at_altitude(h_iss) / g_at_altitude(0.0):.1f} % of surface)")
    print(f"circular orbit at {h_iss/1e3:.0f} km: v = {circular_speed(r_iss):.1f} m/s, T = {orbital_period(r_iss)/60:.2f} min, "
          f"escape speed {escape_speed(r_iss):.1f} m/s = sqrt(2) x circular")
    r_geo = (MU_EARTH * (DAY_SIDEREAL / (2 * np.pi)) ** 2) ** (1 / 3)
    print(f"geostationary radius {r_geo/1e3:.1f} km (altitude {(r_geo - R_EARTH)/1e3:.1f} km), v = {circular_speed(r_geo):.1f} m/s")
    dv1, dv2, tof = hohmann_transfer(r_iss, r_geo)
    print(f"Hohmann {h_iss/1e3:.0f} km -> GEO: dv1 = {dv1:.1f} m/s, dv2 = {dv2:.1f} m/s, total {dv1 + dv2:.1f} m/s, "
          f"time of flight {tof/3600:.2f} h")
    a_t = 0.5 * (r_iss + r_geo)
    print(f"transfer ellipse: a = {a_t/1e3:.1f} km, e = {(r_geo - r_iss)/(r_geo + r_iss):.4f}, "
          f"energy -mu/2a = {-MU_EARTH/(2*a_t)/1e6:.3f} MJ/kg against {-MU_EARTH/(2*r_iss)/1e6:.3f} MJ/kg in the low orbit")

    # the test orbit
    a, e = 7500e3, 0.1
    T = orbital_period(a)
    r0, v0 = perigee_state(a, e)
    print(f"\ntest orbit a = {a/1e3:.0f} km, e = {e}: perigee altitude {(a*(1-e) - R_EARTH)/1e3:.0f} km, "
          f"apogee altitude {(a*(1+e) - R_EARTH)/1e3:.0f} km, T = {T:.1f} s, "
          f"v_p = {np.linalg.norm(v0):.1f} m/s, v_a = {float(vis_viva(a*(1+e), a)):.1f} m/s")
    el = orbital_elements(r0, v0)
    print(f"elements recovered from the perigee state: a = {el[0]/1e3:.1f} km, e = {el[1]:.6f}")

    print("\nleapfrog vs DOP853 (rtol 1e-12) after one orbit; energy and angular-momentum error over ten orbits")
    print(f"{'steps/orbit':>12} {'pos err 1 orbit (m)':>20} {'max |dE/E| 10 orbits':>22} {'max |dh/h| 10 orbits':>22}")
    hs, errs = [], []
    t_ref, R_ref, _ = propagate_solve_ivp(a, e, 1.0, n_out=2001)
    for spo in (250, 500, 1000, 2000):
        t, R, V = propagate_leapfrog(a, e, 1.0, spo)
        Ri = np.column_stack([np.interp(t_ref, t, R[:, 0]), np.interp(t_ref, t, R[:, 1])])
        err = position_error(Ri, R_ref)
        t10, R10, V10 = propagate_leapfrog(a, e, 10.0, spo)
        dE = np.max(np.abs(energy_drift(R10, V10)))
        hh = specific_angular_momentum(R10, V10)
        dh = np.max(np.abs(hh - hh[0]) / abs(hh[0]))
        print(f"{spo:>12d} {err:>20.3e} {dE:>22.3e} {dh:>22.3e}")
        hs.append(T / spo)
        errs.append(err)
    print(f"observed order of the leapfrog position error: {convergence_order(hs, errs):.2f}")

    spo = 1000
    t, R, V = propagate_leapfrog(a, e, 10.0, spo)
    dE = energy_drift(R, V)
    print(f"leapfrog at {spo} steps/orbit: |dE/E| after orbit 1: {np.max(np.abs(dE[:spo+1])):.3e}, "
          f"after orbit 10: {np.max(np.abs(dE)):.3e} (bounded, no secular growth)")
    Tsim = period_from_trajectory(t, R)
    print(f"period from the simulated periapsis passages {Tsim:.2f} s vs Kepler {T:.2f} s "
          f"(difference {Tsim - T:+.3f} s)")
    dA = areal_velocity(R, T / spo)
    print(f"areal velocity: mean {dA.mean():.6e} m^2/s, relative spread {(dA.max() - dA.min())/dA.mean():.2e}; "
          f"|h|/2 = {0.5*abs(specific_angular_momentum(r0, v0)):.6e} m^2/s")

    print("\nKepler's third law from simulated periods (T^2/a^3 should be 4 pi^2/mu = %.4e s^2/m^3)" % (4 * np.pi ** 2 / MU_EARTH))
    for ak in (7000e3, 10000e3, 20000e3, r_geo):
        tk, Rk, _ = propagate_leapfrog(ak, 0.2, 3.0, 1000)
        Tk = period_from_trajectory(tk, Rk)
        print(f"  a = {ak/1e3:>8.1f} km: T_sim = {Tk:>10.2f} s, T^2/a^3 = {float(keplers_third_law_ratio(ak, Tk)):.4e}")

    # the full two-body problem: Earth and Moon
    print("\nfull two-body problem: Earth and Moon from a circular relative orbit")
    mu_em = G_NEWTON * (M_EARTH + M_MOON)
    v_rel = circular_speed(D_EARTH_MOON, mu_em)
    T_em = orbital_period(D_EARTH_MOON, mu_em)
    sim = simulate_two_body(M_EARTH, M_MOON, np.array([D_EARTH_MOON, 0.0]), np.array([0.0, v_rel]), T_em)
    print(f"  period from G(M+m): {T_em/86400:.3f} days; Earth's orbit radius about the barycentre "
          f"{np.max(np.linalg.norm(sim['r1'], axis=1))/1e3:.0f} km (inside the Earth, radius {R_EARTH/1e3:.0f} km)")
    p_scale = M_MOON * v_rel
    print(f"  centre of mass moved {np.max(np.linalg.norm(sim['com'] - sim['com'][0], axis=1)):.2e} m; "
          f"total momentum changed by {np.max(np.linalg.norm(sim['momentum'] - sim['momentum'][0], axis=1))/p_scale:.2e} "
          f"of the Moon's momentum")
    eps_rel = specific_energy(sim["r_rel"], sim["v_rel"], mu_em)
    print(f"  relative-motion energy -mu/2a: {-mu_em/(2*D_EARTH_MOON):.4f} J/kg, simulated {eps_rel[0]:.4f} .. {eps_rel[-1]:.4f} J/kg")

    # figure
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.4))
    ax = axs[0]
    ax.plot(R[:spo + 1, 0] / 1e3, R[:spo + 1, 1] / 1e3, lw=1.2, label="leapfrog, 1000 steps/orbit")
    ax.plot(R_ref[:, 0] / 1e3, R_ref[:, 1] / 1e3, "k--", lw=0.8, label="DOP853 oracle")
    ax.add_patch(plt.Circle((0, 0), R_EARTH / 1e3, color="0.8"))
    ax.set_aspect("equal")
    ax.set_xlabel("x (km)")
    ax.set_ylabel("y (km)")
    ax.set_title(f"a = {a/1e3:.0f} km, e = {e}")
    ax.legend(fontsize=8)
    ax = axs[1]
    ax.plot(t / T, dE, lw=0.7)
    ax.set_xlabel("orbits")
    ax.set_ylabel("relative energy error")
    ax.set_title("leapfrog energy error: bounded")
    ax = axs[2]
    ax.loglog(np.array(hs), errs, "o-", label="leapfrog position error")
    ax.loglog(np.array(hs), errs[0] * (np.array(hs) / hs[0]) ** 2, "k:", label="slope 2")
    ax.set_xlabel("step h (s)")
    ax.set_ylabel("position error after one orbit (m)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "two_body.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
