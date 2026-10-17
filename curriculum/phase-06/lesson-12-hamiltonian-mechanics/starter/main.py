"""P06.12 Hamiltonian Mechanics: starter code.

Mechanics in phase space, where a state is a point (q, p) and a motion is a curve:

* the Legendre transform from a Lagrangian L(q, q') to a Hamiltonian H(q, p), Hamilton's equations
  q' = dH/dp, p' = -dH/dq and the Poisson bracket {f, g}, all in SymPy;
* the pendulum's phase portrait: libration, rotation and the separatrix; the exact period from the
  complete elliptic integral, checked against the derivative of the enclosed phase area, T = dA/dE;
* Liouville's theorem: the exact flow of a Hamiltonian system preserves phase-space area, which is
  measured by carrying the boundary of a blob of initial states and computing its area; with linear
  damping the area shrinks exactly as exp(-c t);
* the one-step maps of P03.07 (explicit Euler, RK4, Stormer-Verlet) as area maps: the determinant of
  each step's Jacobian, and what each does to the blob and to the energy;
* the Kepler problem in polar canonical coordinates: the effective potential, the conserved momentum
  p_theta, the Laplace-Runge-Lenz vector, the radial action J_r whose energy derivative is the period,
  the symplectic Jacobian of the flow, and the along-track growth of a semi-major-axis error.

Conventions: the pendulum uses q = theta (rad) and p = theta' (rad/s), i.e. H per unit m l^2, so
H = p^2/2 - w0^2 cos q; Kepler uses unit mass, so p is velocity (m/s) and H is specific energy (J/kg).

Run:  python main.py [out_dir]     (prints the tables, writes phase_space.png)
Test: pytest curriculum/phase-06/lesson-12-hamiltonian-mechanics
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable

import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.special import ellipk

MU_EARTH = 3.986004418e14       # m^3/s^2, Earth's gravitational parameter (as in P03.07)

Accel = Callable[[np.ndarray], np.ndarray]


# ---------------------------------------------------------------- symbolic: Legendre, Hamilton, Poisson

def legendre_transform(L: sp.Expr, qs: list[sp.Symbol], qdots: list[sp.Symbol],
                       ps: list[sp.Symbol]) -> tuple[sp.Expr, list[sp.Expr]]:
    """Hamiltonian H(q, p) = sum p_i q'_i - L from a Lagrangian L(q, q') written with plain symbols.
    The momenta are p_i = dL/dq'_i; solve those equations for the q'_i and substitute.
    Returns (H, [p_i as functions of q, q'])."""
    raise NotImplementedError


def hamilton_equations(H: sp.Expr, qs: list[sp.Symbol], ps: list[sp.Symbol]) -> tuple[list[sp.Expr], list[sp.Expr]]:
    """Hamilton's equations: returns ([dq_i/dt = dH/dp_i], [dp_i/dt = -dH/dq_i])."""
    raise NotImplementedError


def poisson_bracket(f: sp.Expr, g: sp.Expr, qs: list[sp.Symbol], ps: list[sp.Symbol]) -> sp.Expr:
    """{f, g} = sum_i (df/dq_i dg/dp_i - df/dp_i dg/dq_i). A quantity f with no explicit time
    dependence is conserved exactly when {f, H} = 0."""
    raise NotImplementedError


# ---------------------------------------------------------------- the pendulum in phase space

def pendulum_hamiltonian(q: np.ndarray, p: np.ndarray, w0: float) -> np.ndarray:
    """H = p^2/2 - w0^2 cos q (per unit m l^2; minimum -w0^2 at the bottom, w0^2 on the separatrix)."""
    raise NotImplementedError


def pendulum_flow(t: float, y: np.ndarray, w0: float, c: float = 0.0) -> np.ndarray:
    """Hamilton's equations of the pendulum, with optional linear damping c (1/s) on the momentum:
    q' = p, p' = -w0^2 sin q - c p. ``y`` may stack many states: y = [q_1..q_N, p_1..p_N]."""
    raise NotImplementedError


def separatrix_momentum(q: np.ndarray, w0: float) -> np.ndarray:
    """Upper branch of the separatrix H = w0^2: p = 2 w0 cos(q/2), for -pi <= q <= pi."""
    raise NotImplementedError


def pendulum_period(E: float, w0: float) -> float:
    """Exact period of the pendulum orbit of energy E (H value), from the complete elliptic integral
    K(m) (scipy.special.ellipk takes the parameter m = k^2). With k^2 = (E + w0^2)/(2 w0^2):
    libration (k < 1): T = 4 K(k^2)/w0; rotation (k > 1): one turn takes T = 2 K(1/k^2)/(k w0).
    Infinite on the separatrix (k = 1)."""
    raise NotImplementedError


def phase_area(E: float, w0: float) -> float:
    """Phase-space area of the pendulum orbit of energy E: for libration the area enclosed by the
    closed curve, 2 * integral of p dq between the turning points; for rotation the area under the
    upper curve over one turn, integral of p dq from -pi to pi. Integrated with scipy.integrate.quad."""
    raise NotImplementedError


def polygon_area(x: np.ndarray, y: np.ndarray) -> float:
    """Signed area of a closed polygon by the shoelace formula (positive when counter-clockwise) (given)."""
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def blob_boundary(q0: float, p0: float, radius_q: float, radius_p: float, n: int) -> tuple[np.ndarray, np.ndarray]:
    """n points counter-clockwise on an ellipse of initial states centred on (q0, p0) (given)."""
    s = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    return q0 + radius_q * np.cos(s), p0 + radius_p * np.sin(s)


def flow_blob(q: np.ndarray, p: np.ndarray, t_end: float, w0: float, c: float = 0.0,
              rtol: float = 1e-10, atol: float = 1e-12) -> tuple[np.ndarray, np.ndarray]:
    """Carry every boundary point of a blob along the exact flow for t_end with one ``solve_ivp`` call
    on the stacked state (DOP853). Returns the final (q, p) of each point."""
    raise NotImplementedError


# ---------------------------------------------------------------- one-step maps (from P03.07)

def euler_step(accel: Accel, q: np.ndarray, p: np.ndarray, h: float) -> tuple[np.ndarray, np.ndarray]:
    """Explicit Euler for q' = p, p' = a(q): both updates use the old values."""
    raise NotImplementedError


def rk4_step(accel: Accel, q: np.ndarray, p: np.ndarray, h: float) -> tuple[np.ndarray, np.ndarray]:
    """Classical fourth-order Runge-Kutta for q' = p, p' = a(q)."""
    raise NotImplementedError


def verlet_step(accel: Accel, q: np.ndarray, p: np.ndarray, h: float) -> tuple[np.ndarray, np.ndarray]:
    """Stormer-Verlet, kick-drift-kick: p += h/2 a(q); q += h p; p += h/2 a(q)."""
    raise NotImplementedError


def iterate(stepper: Callable, accel: Accel, q: np.ndarray, p: np.ndarray, h: float, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Apply a one-step map n times (given)."""
    for _ in range(n):
        q, p = stepper(accel, q, p, h)
    return q, p


def step_jacobian_det(stepper: Callable, accel: Accel, q: float, p: float, h: float, eps: float = 1e-6) -> float:
    """Determinant of the 2 x 2 Jacobian of one step of a 1-DOF map at (q, p), by central differences.
    It is the factor by which the step multiplies the area of a small patch: 1 for a symplectic map."""
    raise NotImplementedError


# ---------------------------------------------------------------- Kepler in phase space

def kepler_hamiltonian(q: np.ndarray, p: np.ndarray, mu: float = MU_EARTH) -> np.ndarray:
    """Specific energy H = |p|^2/2 - mu/|q| for Cartesian position q and velocity p (last axis of size 2)."""
    raise NotImplementedError


def kepler_polar_hamiltonian(r: np.ndarray, pr: np.ndarray, ptheta: float, mu: float = MU_EARTH) -> np.ndarray:
    """H = p_r^2/2 + p_theta^2/(2 r^2) - mu/r: the Kepler Hamiltonian in polar canonical coordinates.
    theta is cyclic (absent), so p_theta (the specific angular momentum) is conserved."""
    raise NotImplementedError


def effective_potential(r: np.ndarray, ptheta: float, mu: float = MU_EARTH) -> np.ndarray:
    """V_eff = p_theta^2/(2 r^2) - mu/r: the radial motion is a 1-DOF system in this potential."""
    raise NotImplementedError


def kepler_flow(t: float, y: np.ndarray, mu: float = MU_EARTH) -> np.ndarray:
    """Hamilton's equations in Cartesian coordinates, y = [x, y, px, py]: q' = p, p' = -mu q/|q|^3."""
    raise NotImplementedError


def lrl_vector(q: np.ndarray, p: np.ndarray, mu: float = MU_EARTH) -> np.ndarray:
    """Laplace-Runge-Lenz vector (per unit mass) in the plane: A = p x h - mu q/|q| with h = q x p
    along z, i.e. A = (p_y h, -p_x h) - mu q/r. It points to periapsis with magnitude mu e."""
    raise NotImplementedError


def radial_action(E: float, ptheta: float, mu: float = MU_EARTH) -> float:
    """Radial action J_r = closed integral of p_r dr = 2 * integral from r_p to r_a of
    sqrt(2 (E - V_eff(r))) dr, by quad, for a bound orbit (E < 0). Closed form for comparison:
    J_r = 2 pi (mu / sqrt(-2E) - p_theta)."""
    # turning points: roots of E r^2 + mu r - p_theta^2/2 = 0
    raise NotImplementedError


def flow_jacobian(y0: np.ndarray, t: float, mu: float = MU_EARTH, rel_eps: float = 1e-7) -> np.ndarray:
    """4 x 4 Jacobian d y(t) / d y(0) of the Kepler flow by central differences of DOP853 runs
    (rtol 1e-12), with a perturbation of rel_eps times each component's scale (|q| for positions,
    |p| for velocities). For a Hamiltonian flow M^T J M = J and det M = 1."""
    raise NotImplementedError


J4 = np.block([[np.zeros((2, 2)), np.eye(2)], [-np.eye(2), np.zeros((2, 2))]])


def symplectic_residual(M: np.ndarray) -> float:
    """max |M^T J M - J| for a 4 x 4 Jacobian in (q, p) ordering (given)."""
    return float(np.max(np.abs(M.T @ J4 @ M - J4)))


def periapsis_state(a: float, e: float, mu: float = MU_EARTH) -> np.ndarray:
    """Cartesian state [x, y, vx, vy] at periapsis of an orbit with semi-major axis a and eccentricity e (given)."""
    rp = a * (1.0 - e)
    vp = np.sqrt(mu * (1.0 + e) / rp)
    return np.array([rp, 0.0, 0.0, vp])


def along_track_drift(a: float, da: float, n_orbits: float, e: float = 0.0) -> float:
    """Along-track separation (m) at periapsis after n_orbits between two orbits that start together at
    periapsis with semi-major axes a and a + da, to first order in da/a. The period grows by
    dT = (3/2)(da/a) T, so after n orbits the second satellite arrives n dT late and is a distance
    v_p n dT behind; with v_p T = 2 pi a sqrt((1 + e)/(1 - e)) this is 3 pi da n sqrt((1 + e)/(1 - e))
    (3 pi da per orbit on a circular orbit)."""
    raise NotImplementedError


def energy_error_history(stepper: Callable, y0: np.ndarray, steps_per_orbit: int, n_orbits: int,
                         mu: float = MU_EARTH) -> tuple[np.ndarray, np.ndarray]:
    """Relative energy error (H - H0)/|H0| after every step of a fixed-step map, with the time in
    orbits (given)."""
    q, p = y0[:2].copy(), y0[2:].copy()
    E0 = float(kepler_hamiltonian(q, p, mu))
    a = -mu / (2 * E0)
    T = 2 * np.pi * np.sqrt(a ** 3 / mu)
    h = T / steps_per_orbit
    acc = lambda x: -mu * x / np.sqrt(x[0] * x[0] + x[1] * x[1]) ** 3
    n = steps_per_orbit * n_orbits
    errs = np.empty(n + 1)
    errs[0] = 0.0
    for k in range(n):
        q, p = stepper(acc, q, p, h)
        errs[k + 1] = (0.5 * (p[0] * p[0] + p[1] * p[1]) - mu / np.sqrt(q[0] * q[0] + q[1] * q[1]) - E0) / abs(E0)
    return np.arange(n + 1) / steps_per_orbit, errs


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save phase_space.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    # symbolic: pendulum and Kepler
    th, thd, pth = sp.symbols("theta thetadot p_theta")
    m, l, g = sp.symbols("m l g", positive=True)
    H, _ = legendre_transform(m * l ** 2 * thd ** 2 / 2 + m * g * l * sp.cos(th), [th], [thd], [pth])
    print("pendulum: H =", H, "  Hamilton:", hamilton_equations(H, [th], [pth]))
    r, rd, thetd, pr, pt, mu = sp.symbols("r rdot thetadot p_r p_theta mu", positive=True)
    Hk, _ = legendre_transform((rd ** 2 + r ** 2 * thetd ** 2) / 2 + mu / r, [r, th], [rd, thetd], [pr, pt])
    print("Kepler (polar): H =", Hk, "  {p_theta, H} =", poisson_bracket(pt, Hk, [r, th], [pr, pt]))

    # pendulum: period two ways
    w0 = 2 * np.pi / 2.0                         # a 1 m-class pendulum: 2 s small-angle period
    print(f"\npendulum w0 = {w0:.5f} rad/s (small-angle period 2.000 s)")
    print("   amplitude   T elliptic (s)   T = dA/dE (s)")
    for deg in (10.0, 60.0, 120.0, 170.0):
        E = -w0 ** 2 * np.cos(np.radians(deg))
        dE = 1e-6 * w0 ** 2
        TA = (phase_area(E + dE, w0) - phase_area(E - dE, w0)) / (2 * dE)
        print(f"   {deg:7.1f}     {pendulum_period(E, w0):.6f}         {TA:.6f}")

    # Liouville: a blob near the separatrix
    q0, p0, rq, rp = 0.0, 1.2 * w0, 0.3, 0.3 * w0
    qb, pb = blob_boundary(q0, p0, rq, rp, 4000)
    A0 = polygon_area(qb, pb)
    t_b = 6.0
    qe, pe = flow_blob(qb, pb, t_b, w0)
    qd, pdm = flow_blob(qb, pb, t_b, w0, c=0.3)
    print(f"\nLiouville: blob area {A0:.7f}; after {t_b:.0f} s exact flow {polygon_area(qe, pe):.7f} "
          f"(change {polygon_area(qe, pe) / A0 - 1:+.1e}, polygon error); with damping c = 0.3/s "
          f"{polygon_area(qd, pdm):.7f} (A0 exp(-c t) = {A0 * np.exp(-0.3 * t_b):.7f})")
    acc = lambda q: -w0 ** 2 * np.sin(q)
    h = 0.05
    for name, st, tt in (("Euler", euler_step, 1.0), ("RK4", rk4_step, t_b), ("Verlet", verlet_step, t_b)):
        qs, ps = iterate(st, acc, qb, pb, h, int(round(tt / h)))
        print(f"   {name:6s} h = {h}: area after {tt:.0f} s {polygon_area(qs, ps):.7f} (change {polygon_area(qs, ps) / A0 - 1:+.2e}); "
              f"one-step det at the blob centre {step_jacobian_det(st, acc, q0, p0, h):.9f}")

    # Kepler: invariants, action, Jacobian, along-track
    a, e = 7.0e6, 0.1
    y0 = periapsis_state(a, e)
    T = 2 * np.pi * np.sqrt(a ** 3 / MU_EARTH)
    E0 = float(kepler_hamiltonian(y0[:2], y0[2:]))
    h0 = y0[0] * y0[3] - y0[1] * y0[2]
    sol = solve_ivp(kepler_flow, (0, 5 * T), y0, method="DOP853", rtol=1e-12, atol=1e-6, t_eval=np.linspace(0, 5 * T, 2001))
    Y = sol.y.T
    Es = kepler_hamiltonian(Y[:, :2], Y[:, 2:])
    hs = Y[:, 0] * Y[:, 3] - Y[:, 1] * Y[:, 2]
    A = lrl_vector(Y[:, :2], Y[:, 2:])
    print(f"\nKepler a = {a / 1e3:.0f} km, e = {e}: period {T:.2f} s; over 5 orbits (DOP853, rtol 1e-12):")
    print(f"   energy {np.ptp(Es) / abs(E0):.1e}, angular momentum {np.ptp(hs) / h0:.1e}, "
          f"LRL direction {np.ptp(np.arctan2(A[:, 1], A[:, 0])):.1e} rad; |A|/mu = {np.linalg.norm(A[0]) / MU_EARTH:.6f}")
    Jr = radial_action(E0, h0)
    Jc = 2 * np.pi * (MU_EARTH / np.sqrt(-2 * E0) - h0)
    dE = 1e-6 * abs(E0)
    TJ = (radial_action(E0 + dE, h0) - radial_action(E0 - dE, h0)) / (2 * dE)
    print(f"   radial action J_r = {Jr:.6e} m^2/s (closed form {Jc:.6e}); dJ_r/dE = {TJ:.3f} s vs period {T:.3f} s")
    # flow Jacobian in canonical units (length a, time 1/n, mu = 1) so that its entries are of order 1
    yc = periapsis_state(1.0, e, 1.0)
    M = flow_jacobian(yc, 2 * np.pi, 1.0)
    print(f"   flow Jacobian over one orbit (canonical units): det {np.linalg.det(M):.10f}, "
          f"max|M^T J M - J| {symplectic_residual(M):.1e}")
    da = 10.0
    y1 = periapsis_state(a + da, e)
    s1 = solve_ivp(kepler_flow, (0, 20 * T), y1, method="DOP853", rtol=1e-12, atol=1e-6, dense_output=True)
    s0 = solve_ivp(kepler_flow, (0, 20 * T), y0, method="DOP853", rtol=1e-12, atol=1e-6, dense_output=True)
    for n in (1, 5, 20):
        d = s1.sol(n * T)[:2] - s0.sol(n * T)[:2]
        v = s0.sol(n * T)[2:]
        along = abs(d @ v) / np.linalg.norm(v)
        print(f"   da = {da:.0f} m: after {n:2d} orbits along-track {along:8.1f} m, total {np.linalg.norm(d):8.1f} m "
              f"(first order {along_track_drift(a, da, n, e):8.1f} m)")
    hist = {}
    for name, st in (("RK4", rk4_step), ("Verlet", verlet_step)):
        tt, err = energy_error_history(st, y0, 200, 100)
        hist[name] = (tt, err)
        first, last = np.max(np.abs(err[: 200 * 10 + 1])), np.max(np.abs(err[-200 * 10:]))
        print(f"   {name:6s} 200 steps/orbit: largest |energy error| in orbits 1-10 {first:.2e}, in orbits 91-100 {last:.2e}, "
              f"at the end {err[-1]:+.2e}")

    # figure
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    qq = np.linspace(-np.pi, np.pi, 401)
    for Ev in np.linspace(-0.9, 2.0, 12) * w0 ** 2:
        pp = np.sqrt(np.clip(2 * (Ev + w0 ** 2 * np.cos(qq)), 0, None))
        ax[0].plot(qq, pp, "0.6", lw=0.7)
        ax[0].plot(qq, -pp, "0.6", lw=0.7)
    ax[0].plot(qq, separatrix_momentum(qq, w0), "k", lw=1.5, label="separatrix")
    ax[0].plot(qq, -separatrix_momentum(qq, w0), "k", lw=1.5)
    ax[0].fill(qb, pb, alpha=0.5, label="blob t = 0")
    ax[0].fill(qe, pe, alpha=0.5, label="exact flow t = 6 s")
    ax[0].set_xlabel("q = theta (rad)")
    ax[0].set_ylabel("p = theta' (rad/s)")
    ax[0].set_title("pendulum: the flow keeps the blob's area, not its shape")
    ax[0].legend(fontsize=8)
    rr = np.linspace(0.75 * a * (1 - e), 1.3 * a * (1 + e), 400)
    ax[1].plot(rr / 1e3, effective_potential(rr, h0) / 1e6, label="V_eff")
    ax[1].axhline(E0 / 1e6, color="k", ls="--", label="orbit energy")
    ax[1].set_xlabel("r (km)")
    ax[1].set_ylabel("specific energy (MJ/kg)")
    ax[1].set_title("Kepler radial motion in V_eff")
    ax[1].legend()
    for name, (tt, err) in hist.items():
        ax[2].plot(tt, err, lw=0.6, label=name)
    ax[2].set_xlabel("orbits")
    ax[2].set_ylabel("relative energy error")
    ax[2].set_title("200 steps per orbit")
    ax[2].legend()
    fig.tight_layout()
    fig.savefig(out_dir / "phase_space.png", dpi=110)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'phase_space.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
