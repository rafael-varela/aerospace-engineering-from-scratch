"""P03.07 Symplectic Integrators: starter code.

Fixed-step integrators for second-order systems q'' = a(q) (unit mass, so the momentum p equals
the velocity), written from scratch and compared on a two-body orbit:

* explicit Euler (``euler_step``): not symplectic, energy grows;
* symplectic Euler / Euler-Cromer (``symplectic_euler_step``): first order, symplectic;
* Stormer-Verlet / velocity Verlet / leapfrog (``stormer_verlet_step``): second order, symplectic,
  time-reversible, conserves angular momentum exactly for central forces;
* Yoshida's fourth-order composition of three Verlet steps (``yoshida4_step``): fourth order, symplectic;
* classical RK4 on the first-order system (``rk4_step``): fourth order, not symplectic.

Oracle: the closed-form two-body solution via Kepler's equation (``kepler_exact``) and
``scipy.integrate.solve_ivp`` (DOP853 at tight tolerance) for a short arc.

Conventions: states are (q, p) with q, p arrays of length 2 (planar motion, SI units); the
gravitational parameter mu is in m^3/s^2; energy is per unit mass (J/kg).

Run:  python main.py          (prints the verification tables, writes symplectic_orbit.png)
Test: pytest curriculum/phase-03/lesson-07-symplectic-integrators
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np

Accel = Callable[[np.ndarray], np.ndarray]
MU_EARTH = 3.986004418e14           # m^3/s^2
R_EARTH = 6.371e6                   # m

# Yoshida (1990) fourth-order composition coefficients: three Verlet steps of sizes w1 h, w0 h, w1 h.
_CBRT2 = 2.0 ** (1.0 / 3.0)
YOSHIDA_W1 = 1.0 / (2.0 - _CBRT2)
YOSHIDA_W0 = -_CBRT2 / (2.0 - _CBRT2)


# ---------------------------------------------------------------- forces and invariants

def kepler_acceleration(q: np.ndarray, mu: float = MU_EARTH) -> np.ndarray:
    """Two-body acceleration -mu q / |q|^3."""
    r2 = q[0] * q[0] + q[1] * q[1]
    return (-mu / (r2 * np.sqrt(r2))) * q


def energy(q: np.ndarray, p: np.ndarray, mu: float = MU_EARTH) -> np.ndarray:
    """Specific orbital energy |p|^2/2 - mu/|q| [J/kg]; works on arrays of shape (n, 2) too."""
    q, p = np.atleast_2d(q), np.atleast_2d(p)
    E = 0.5 * np.sum(p ** 2, axis=1) - mu / np.linalg.norm(q, axis=1)
    return E if E.size > 1 else float(E[0])


def angular_momentum(q: np.ndarray, p: np.ndarray) -> np.ndarray:
    """Specific angular momentum q x p (the z component, planar motion) [m^2/s]."""
    q, p = np.atleast_2d(q), np.atleast_2d(p)
    L = q[:, 0] * p[:, 1] - q[:, 1] * p[:, 0]
    return L if L.size > 1 else float(L[0])


def orbital_period(q: np.ndarray, p: np.ndarray, mu: float = MU_EARTH) -> float:
    """Period 2 pi sqrt(a^3 / mu) with a = -mu / (2 E); raises ValueError for an unbound orbit."""
    E = energy(q, p, mu)
    if E >= 0:
        raise ValueError("orbit is not bound (E >= 0)")
    a = -mu / (2.0 * E)
    return float(2.0 * np.pi * np.sqrt(a ** 3 / mu))


def perigee_state(a: float, e: float, mu: float = MU_EARTH) -> tuple[np.ndarray, np.ndarray]:
    """Position and velocity at perigee of an orbit with semi-major axis a and eccentricity e,
    perigee on the +x axis, motion counter-clockwise."""
    if not (0.0 <= e < 1.0):
        raise ValueError("need 0 <= e < 1")
    rp = a * (1.0 - e)
    vp = np.sqrt(mu * (1.0 + e) / rp)
    return np.array([rp, 0.0]), np.array([0.0, vp])


# ---------------------------------------------------------------- one-step maps (q, p) -> (q, p)

def euler_step(accel: Accel, q: np.ndarray, p: np.ndarray, h: float) -> tuple[np.ndarray, np.ndarray]:
    """Explicit Euler on the first-order system: q += h p, p += h a(q_old). Not symplectic."""
    raise NotImplementedError


def symplectic_euler_step(accel: Accel, q: np.ndarray, p: np.ndarray, h: float) -> tuple[np.ndarray, np.ndarray]:
    """Symplectic Euler (Euler-Cromer): p += h a(q_old), then q += h p_new. First order, symplectic."""
    raise NotImplementedError


def stormer_verlet_step(accel: Accel, q: np.ndarray, p: np.ndarray, h: float) -> tuple[np.ndarray, np.ndarray]:
    """Velocity Verlet (Stormer-Verlet, leapfrog): kick-drift-kick.

        p_half = p + (h/2) a(q);  q_new = q + h p_half;  p_new = p_half + (h/2) a(q_new).

    Second order, symplectic, symmetric (time-reversible). Two force evaluations as written; the
    second equals the first of the next step, so a careful implementation uses one per step.
    """
    raise NotImplementedError


def yoshida4_step(accel: Accel, q: np.ndarray, p: np.ndarray, h: float) -> tuple[np.ndarray, np.ndarray]:
    """Fourth-order symplectic step: Verlet with sub-steps w1 h, w0 h, w1 h (Yoshida's composition),
    w1 = 1/(2 - 2^(1/3)) = 1.3512, w0 = 1 - 2 w1 = -1.7024 (a backward step in the middle)."""
    raise NotImplementedError


def rk4_step(accel: Accel, q: np.ndarray, p: np.ndarray, h: float) -> tuple[np.ndarray, np.ndarray]:
    """Classical RK4 on y = (q, p), y' = (p, a(q)). Fourth order, not symplectic."""
    raise NotImplementedError


STEPPERS = {"euler": euler_step, "symplectic_euler": symplectic_euler_step, "verlet": stormer_verlet_step,
            "yoshida4": yoshida4_step, "rk4": rk4_step}
FORCE_EVALS = {"euler": 1, "symplectic_euler": 1, "verlet": 2, "yoshida4": 6, "rk4": 4}


# ---------------------------------------------------------------- drivers

def propagate(stepper: Callable, accel: Accel, q0: np.ndarray, p0: np.ndarray, h: float,
              n_steps: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Take ``n_steps`` steps of size h. Returns (t, q, p) with q, p of shape (n_steps + 1, 2)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the oracle: Kepler's equation

def solve_kepler(M: float, e: float, tol: float = 1e-14) -> float:
    """Eccentric anomaly E from M = E - e sin E by Newton's method (P03.02)."""
    raise NotImplementedError


def kepler_exact(t: np.ndarray | float, q0: np.ndarray, p0: np.ndarray, mu: float = MU_EARTH
                 ) -> tuple[np.ndarray, np.ndarray]:
    """Closed-form two-body position and velocity at times t from the state (q0, p0) at t = 0.

    Elements from the state (a from the energy, e and the eccentricity vector from h and v x h),
    the eccentric anomaly at t = 0 from the state, then E(t) from Kepler's equation and the
    position in the perifocal frame rotated by the argument of perigee. Planar, elliptic only.
    Returns arrays of shape (len(t), 2).
    """
    raise NotImplementedError


# ---------------------------------------------------------------- symplecticity and reversibility

def step_jacobian(stepper: Callable, accel: Accel, q: np.ndarray, p: np.ndarray, h: float,
                  eps: float = 1e-6) -> np.ndarray:
    """Central-difference Jacobian of the one-step map (q, p) -> (q', p'), shape (4, 4)."""
    raise NotImplementedError


def symplectic_defect(M: np.ndarray) -> float:
    """||M^T J M - J|| (Frobenius) with J = [[0, I], [-I, 0]]: zero iff the map's Jacobian is symplectic."""
    raise NotImplementedError


def reversibility_defect(stepper: Callable, accel: Accel, q: np.ndarray, p: np.ndarray, h: float) -> float:
    """|(q, p) - back| after one step of +h and one of -h: zero for a symmetric method."""
    raise NotImplementedError


def modified_energy_oscillator(q: np.ndarray, p: np.ndarray, omega: float, h: float) -> np.ndarray:
    """The quantity Stormer-Verlet conserves *exactly* on q'' = -omega^2 q:
    H_h = p^2/2 + (omega^2/2)(1 - (h omega)^2/4) q^2 (the modified or shadow Hamiltonian)."""
    raise NotImplementedError


# ---------------------------------------------------------------- studies

def energy_drift(stepper: Callable, q0: np.ndarray, p0: np.ndarray, n_orbits: int, steps_per_orbit: int,
                 mu: float = MU_EARTH) -> tuple[np.ndarray, np.ndarray]:
    """Relative energy error (E - E0)/|E0| at every step over ``n_orbits`` orbits. Returns (t/P, dE)."""
    raise NotImplementedError


def drift_slope(t: np.ndarray, dE: np.ndarray) -> float:
    """Least-squares slope of the relative energy error against time (per unit of t): the secular drift rate."""
    raise NotImplementedError


def convergence_order(hs: np.ndarray, errs: np.ndarray) -> float:
    """Least-squares slope of log(err) against log(h): the observed order of accuracy."""
    hs, errs = np.asarray(hs, dtype=float), np.asarray(errs, dtype=float)
    if np.any(errs <= 0):
        raise ValueError("errors must be positive to take logarithms")
    return float(np.polyfit(np.log(hs), np.log(errs), 1)[0])


def position_errors(stepper: Callable, q0: np.ndarray, p0: np.ndarray, t_end: float, n_steps_list: list[int],
                    mu: float = MU_EARTH) -> np.ndarray:
    """|q(t_end) - exact| for each step count: the data for an observed-order study."""
    raise NotImplementedError


def main(out_dir: str | Path = ".") -> None:
    """Run the lesson's experiments, print the tables and save ``symplectic_orbit.png``."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.integrate import solve_ivp

    out_dir = Path(out_dir)
    a, e = 7000e3, 0.1
    q0, p0 = perigee_state(a, e)
    P = orbital_period(q0, p0)
    print(f"orbit: a = {a / 1e3:.0f} km, e = {e}, perigee {np.linalg.norm(q0) / 1e3:.0f} km, v_p = {p0[1]:.1f} m/s, "
          f"P = {P:.1f} s = {P / 60:.2f} min, E = {energy(q0, p0):.4e} J/kg, h = {angular_momentum(q0, p0):.4e} m^2/s")

    # ---- the oracle against solve_ivp
    t_chk = np.array([0.25, 0.5, 1.0, 2.5]) * P
    q_ex, p_ex = kepler_exact(t_chk, q0, p0)
    f = lambda t, y: np.concatenate([y[2:], kepler_acceleration(y[:2])])
    sol = solve_ivp(f, (0.0, t_chk[-1]), np.concatenate([q0, p0]), method="DOP853", rtol=1e-13, atol=1e-6, t_eval=t_chk)
    print(f"Kepler closed form vs DOP853 (rtol 1e-13): max position difference {np.max(np.abs(sol.y[:2].T - q_ex)):.2e} m, "
          f"velocity {np.max(np.abs(sol.y[2:].T - p_ex)):.2e} m/s; after one period |q(P) - q0| = {np.linalg.norm(q_ex[2] - q0):.2e} m")

    # ---- observed orders over 1.5 orbits (perigee to apogee); and the full-period surprise
    ns = [300, 600, 1200, 2400]
    t_end = 1.5 * P
    print("\nobserved order of the position error after 1.5 orbits:")
    for name in ("symplectic_euler", "verlet", "rk4", "yoshida4"):
        errs = position_errors(STEPPERS[name], q0, p0, t_end, ns)
        print(f"  {name:17s} errors [m] {np.array2string(errs, precision=2)}  order {convergence_order(t_end / np.array(ns), errs):.2f}")
    errs = position_errors(symplectic_euler_step, q0, p0, P, [200, 400, 800, 1600])
    print(f"  symplectic Euler after exactly one orbit: errors {np.array2string(errs, precision=1)} m, "
          f"order {convergence_order(P / np.array([200, 400, 800, 1600]), errs):.2f} (the O(h) term cancels over a full period)")

    # ---- symplecticity, reversibility (in units mu = 1, r_p = 1 so that all entries of M are O(1))
    print("\none-step map of the Kepler problem in units mu = 1 (perigee of e = 0.1, h = 0.05):")
    acc1 = lambda q: kepler_acceleration(q, 1.0)
    q1, p1 = perigee_state(1.0 / (1.0 - 0.1), 0.1, mu=1.0)
    for name, st in STEPPERS.items():
        M = step_jacobian(st, acc1, q1, p1, 0.05)
        print(f"  {name:17s} det M - 1 = {np.linalg.det(M) - 1:+.2e}, symplectic defect {symplectic_defect(M):.2e}, "
              f"reversibility defect {reversibility_defect(st, acc1, q1, p1, 0.05):.2e}")
    omega = 2.0 * np.pi
    acc_osc = lambda q: -omega ** 2 * q
    for name in ("verlet", "rk4"):
        M = step_jacobian(STEPPERS[name], acc_osc, np.array([1.0, 0.0]), np.array([0.0, 0.0]), 1.0 / omega)
        print(f"  oscillator h omega = 1: {name:7s} det M = {np.linalg.det(M):.12f}")
    t_o, q_o, p_o = propagate(stormer_verlet_step, acc_osc, np.array([1.0, 0.0]), np.array([0.0, 0.0]), 0.5 / omega, 400)
    H = 0.5 * p_o[:, 0] ** 2 + 0.5 * omega ** 2 * q_o[:, 0] ** 2
    Hh = modified_energy_oscillator(q_o[:, 0], p_o[:, 0], omega, 0.5 / omega)
    print(f"  oscillator, Verlet, h omega = 0.5: true energy varies by {np.ptp(H) / H[0]:.3e} relative, "
          f"modified energy by {np.ptp(Hh) / Hh[0]:.3e}")

    # ---- long runs: energy and angular momentum
    n_orb = 60
    print(f"\n{n_orb} orbits (Euler: 5 orbits); force evaluations per orbit count one per Verlet step (reuse):")
    curves = {}
    for name, spo in (("euler", 200), ("symplectic_euler", 200), ("verlet", 200), ("verlet", 800), ("rk4", 200), ("yoshida4", 200)):
        n = n_orb if name != "euler" else 5
        tt, dE = energy_drift(STEPPERS[name], q0, p0, n, spo)
        if spo == 200:
            curves[name] = (tt, dE)
        _, qs, ps = propagate(STEPPERS[name], kepler_acceleration, q0, p0, P / spo, n * spo)
        L = angular_momentum(qs, ps)
        q_ex, _ = kepler_exact(tt[-1] * P, q0, p0)
        evals = spo * (1 if name == "verlet" else 3 if name == "yoshida4" else FORCE_EVALS[name])
        print(f"  {name:17s} {spo:4d} steps/orbit: max |dE/E| = {np.max(np.abs(dE)):.2e}, final dE/E = {dE[-1]:+.2e}, "
              f"slope {drift_slope(tt, dE):+.2e}/orbit, max |dL/L| = {np.max(np.abs(L - L[0]) / abs(L[0])):.1e}, "
              f"position error {np.linalg.norm(qs[-1] - q_ex[0]) / 1e3:.1f} km, {evals} force evals/orbit")

    # ---- figure
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    ax = axes[0, 0]
    for name in ("euler", "symplectic_euler", "verlet"):
        _, qs, _ = propagate(STEPPERS[name], kepler_acceleration, q0, p0, P / spo, 5 * spo)
        ax.plot(qs[:, 0] / 1e3, qs[:, 1] / 1e3, lw=0.8, label=name)
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(R_EARTH / 1e3 * np.cos(th), R_EARTH / 1e3 * np.sin(th), "k-", lw=0.5)
    ax.set_aspect("equal"); ax.set_xlabel("x [km]"); ax.set_ylabel("y [km]"); ax.set_title("5 orbits, 200 steps per orbit"); ax.legend(fontsize=8)
    ax = axes[0, 1]
    for name in ("symplectic_euler", "verlet", "rk4", "yoshida4"):
        tt, dE = curves[name]
        ax.plot(tt, dE, lw=0.6, label=name)
    ax.set_xlabel("orbits"); ax.set_ylabel("(E - E0)/|E0|"); ax.set_title(f"energy error over {n_orb} orbits, 200 steps per orbit"); ax.legend(fontsize=8)
    ax.set_yscale("symlog", linthresh=1e-9)
    ax = axes[1, 0]
    tt, dE = curves["verlet"]
    ax.plot(tt[: 3 * 200], dE[: 3 * 200], lw=0.8)
    ax.set_xlabel("orbits"); ax.set_ylabel("(E - E0)/|E0|"); ax.set_title("Verlet: the energy error oscillates with the orbit (first 3 orbits)")
    ax = axes[1, 1]
    for name in ("symplectic_euler", "verlet", "rk4", "yoshida4"):
        errs = position_errors(STEPPERS[name], q0, p0, t_end, ns)
        ax.loglog(t_end / np.array(ns), errs, "o-", label=f"{name} (order {convergence_order(t_end / np.array(ns), errs):.1f})")
    ax.set_xlabel("step h [s]"); ax.set_ylabel("|q(1.5 P) - exact| [m]"); ax.set_title("position error after 1.5 orbits"); ax.legend(fontsize=8)
    fig.suptitle("P03.07  Symplectic integration of a two-body orbit")
    fig.tight_layout()
    fig.savefig(out_dir / "symplectic_orbit.png", dpi=150)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'symplectic_orbit.png'}")


if __name__ == "__main__":
    main()
