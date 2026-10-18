"""P06.18 Multibody Simulation: reference solution.

A landing-gear drop test as a small planar multibody system. Two bodies:

* body 1, the airframe (or the drop carriage of a test rig), mass M, guided so it moves only vertically
  and does not rotate;
* body 2, the strut's lower part with the axle and wheel, mass m ("unsprung mass"), which slides along
  the strut axis. The axis is raked by an angle beta from the vertical.

They are joined by an oleo-pneumatic shock strut: a polytropic gas spring (the pneumatic force),
quadratic orifice damping (the hydraulic force), and a stiff top-out stop that holds the strut fully
extended in the air. The tyre is a compression-only spring with light damping, on a frictionless
ground (no spin-up drag). A constant lift force on the airframe stands for wing lift.

The same system is written two ways:

1. reduced coordinates (y1, s) and Lagrange's equations: two equations, no constraint forces;
2. full planar coordinates q = (x1, y1, th1, x2, y2, th2) with four constraint equations and
   constrained Newton-Euler: M q'' = Q + Phi_q^T lambda, Phi_q q'' = gamma, solved as one linear system
   per step, with optional Baumgarte stabilisation. The multipliers give the constraint forces,
   including the side load on the strut bearings.

Checks: the two formulations agree; the constraint drift stays small; and the energy put in (initial
kinetic energy plus the work of gravity and lift) equals the energy stored (kinetic, gas, tyre, stop)
plus the energy dissipated (orifice and tyre damping) at every instant.

Units: SI (m, s, kg, N, Pa, J). Parameters are illustrative, of the size of a light business jet's main
gear; they are not data for any aircraft.

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes landing_gear.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-06/lesson-18-multibody-simulation
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

G0 = 9.80665
P_ATM = 101_325.0


@dataclass(frozen=True)
class Gear:
    """Parameters of the two-body landing-gear model (illustrative values)."""
    M: float = 1500.0          # kg, airframe mass carried by this gear (drop-carriage mass)
    m: float = 40.0            # kg, unsprung mass: piston, axle, brake, wheel and tyre
    beta: float = np.radians(8.0)   # rad, strut rake from the vertical
    ell0: float = 1.0          # m, attachment point to axle along the strut at full extension
    R: float = 0.30            # m, unloaded tyre radius (axle height when the tyre just touches)
    A_a: float = 7.854e-3      # m^2, pneumatic (piston) area, 100 mm diameter
    A_h: float = 7.854e-3      # m^2, hydraulic area that pushes oil through the orifice
    A_o: float = 4.0e-4        # m^2, orifice area (22.6 mm diameter)
    Cd: float = 0.7            # orifice discharge coefficient
    rho: float = 850.0         # kg/m^3, hydraulic oil density
    p0: float = 1.5e6          # Pa absolute, gas pressure at full extension
    V0: float = 3.1416e-3      # m^3, gas volume at full extension (0.40 m of piston area)
    n_poly: float = 1.1        # polytropic exponent of the gas compression
    stroke_max: float = 0.30   # m, available stroke
    k_stop: float = 1.0e7      # N/m, top-out stop (resists extension beyond s = 0)
    k_t: float = 8.0e5         # N/m, tyre stiffness
    c_t: float = 500.0         # N s/m, tyre damping
    lift_ratio: float = 1.0    # lift / total weight (drop tests commonly simulate lift = weight)
    I1: float = 2000.0         # kg m^2, airframe pitch inertia (constrained; enters only lambda)
    I2: float = 1.5            # kg m^2, wheel-assembly inertia (constrained)

    @property
    def W(self) -> float:
        """Total weight (M + m) g [N]."""
        return (self.M + self.m) * G0


# ---------------------------------------------------------------- force elements

def gas_force(s: np.ndarray | float, g: Gear) -> np.ndarray | float:
    """Pneumatic force [N] pushing the strut apart at stroke s (s > 0 compresses):
    A_a [p0 (V0 / (V0 - A_a s))^n - p_atm]."""
    return g.A_a * (g.p0 * (g.V0 / (g.V0 - g.A_a * np.asarray(s, float))) ** g.n_poly - P_ATM)


def gas_energy(s: np.ndarray | float, g: Gear) -> np.ndarray | float:
    """Work done compressing the strut from s = 0 to s against ``gas_force`` [J]:
    p0 V0/(n - 1) [(V0/(V0 - A_a s))^(n-1) - 1] - p_atm A_a s."""
    s = np.asarray(s, float)
    ratio = g.V0 / (g.V0 - g.A_a * s)
    return g.p0 * g.V0 / (g.n_poly - 1.0) * (ratio ** (g.n_poly - 1.0) - 1.0) - P_ATM * g.A_a * s


def orifice_force(sdot: np.ndarray | float, g: Gear) -> np.ndarray | float:
    """Hydraulic force [N] opposing the stroke rate: rho A_h^3 sdot |sdot| / (2 (Cd A_o)^2)."""
    sdot = np.asarray(sdot, float)
    return g.rho * g.A_h ** 3 * sdot * np.abs(sdot) / (2.0 * (g.Cd * g.A_o) ** 2)


def stop_force(s: np.ndarray | float, g: Gear) -> np.ndarray | float:
    """Top-out stop force [N], positive when it pushes the strut towards compression (s < 0 only)."""
    return g.k_stop * np.maximum(-np.asarray(s, float), 0.0)


def strut_force(s: float, sdot: float, g: Gear) -> float:
    """Net axial force [N] pushing the two bodies apart: gas + orifice - stop."""
    return float(gas_force(s, g) + orifice_force(sdot, g) - stop_force(s, g))


def tyre_force(y2: float, y2dot: float, g: Gear) -> tuple[float, float]:
    """Vertical tyre force [N] and its elastic part, for axle height y2 and rate y2dot. Deflection
    d = R - y2; force k_t d + c_t d' while in contact, never negative (the tyre cannot pull)."""
    d, ddot = g.R - y2, -y2dot
    if d <= 0.0:
        return 0.0, 0.0
    return max(g.k_t * d + g.c_t * ddot, 0.0), g.k_t * d


# ---------------------------------------------------------------- formulation 1: Lagrange, reduced

def mass_matrix_reduced(g: Gear) -> np.ndarray:
    """Mass matrix of T = M y1'^2/2 + m |r2'|^2/2 in coordinates (y1, s), with the axle at
    x2 = -(ell0 - s) sin(beta), y2 = y1 - (ell0 - s) cos(beta): [[M+m, m cos b], [m cos b, m]]."""
    cb = np.cos(g.beta)
    return np.array([[g.M + g.m, g.m * cb], [g.m * cb, g.m]])


def axle_position(y1: float, s: float, g: Gear) -> tuple[float, float]:
    """(x2, y2) of the axle for airframe height y1 and stroke s."""
    L = g.ell0 - s
    return -L * np.sin(g.beta), y1 - L * np.cos(g.beta)


def generalized_forces(y1: float, s: float, y1dot: float, sdot: float, g: Gear) -> np.ndarray:
    """(Q_y1, Q_s): gravity, lift and tyre on y1; the axle's share of gravity and tyre force along the
    axis, minus the strut force, on s (virtual work with dy2/dy1 = 1, dy2/ds = cos beta)."""
    cb = np.cos(g.beta)
    _, y2 = axle_position(y1, s, g)
    Ft, _ = tyre_force(y2, y1dot + sdot * cb, g)
    Fs = strut_force(s, sdot, g)
    Q1 = -g.W + g.lift_ratio * g.W + Ft
    Q2 = (Ft - g.m * G0) * cb - Fs
    return np.array([Q1, Q2])


def rhs_reduced(t: float, z: np.ndarray, g: Gear) -> np.ndarray:
    """State z = (y1, s, y1', s', D) with D the energy dissipated so far (orifice + tyre damping)."""
    y1, s, v1, vs = z[:4]
    acc = np.linalg.solve(mass_matrix_reduced(g), generalized_forces(y1, s, v1, vs, g))
    _, y2 = axle_position(y1, s, g)
    return np.array([v1, vs, acc[0], acc[1], dissipated_power(y2, v1 + vs * np.cos(g.beta), vs, g)])


def dissipated_power(y2: float, y2dot: float, sdot: float, g: Gear) -> float:
    """Power [W] dissipated by the orifice (F_h s') and by the tyre ((F_t - k_t d) d')."""
    Ft, Fe = tyre_force(y2, y2dot, g)
    return float(orifice_force(sdot, g) * sdot + (Ft - Fe) * (-y2dot))


def initial_stroke(g: Gear) -> float:
    """Stroke s0 < 0 at which the strut hangs before touchdown. Before contact the whole gear accelerates
    at -g (1 - lift_ratio) with s'' = 0, so the s-equation gives gas - stop = -lift_ratio m g cos(beta):
    the strut holds the unsprung weight when lift = weight, and carries nothing in free fall."""
    target = -g.lift_ratio * g.m * G0 * np.cos(g.beta)
    return float(brentq(lambda s: gas_force(s, g) - stop_force(s, g) - target, -0.05, 0.0, xtol=1e-15))


def initial_state(g: Gear, sink: float) -> np.ndarray:
    """Reduced state at touchdown: tyre just touching (y2 = R), strut on its stop, sinking at ``sink``."""
    s0 = initial_stroke(g)
    y1 = g.R + (g.ell0 - s0) * np.cos(g.beta)
    return np.array([y1, s0, -sink, 0.0, 0.0])


def simulate_reduced(g: Gear, sink: float, t_end: float = 0.5, n_out: int = 2001,
                     rtol: float = 1e-10, atol: float = 1e-12, max_step: float = 2e-3) -> dict[str, np.ndarray]:
    """Integrate the Lagrange model from touchdown. Returns arrays t, y1, s, y1dot, sdot, D and the
    derived axle height y2, tyre force Ft, strut force Fs (gas + orifice - stop)."""
    sol = solve_ivp(rhs_reduced, (0.0, t_end), initial_state(g, sink), args=(g,), method="DOP853",
                    t_eval=np.linspace(0.0, t_end, n_out), rtol=rtol, atol=atol, max_step=max_step)
    y1, s, v1, vs, D = sol.y
    cb = np.cos(g.beta)
    y2 = y1 - (g.ell0 - s) * cb
    Ft = np.array([tyre_force(a, b, g)[0] for a, b in zip(y2, v1 + vs * cb)])
    Fs = np.array([strut_force(a, b, g) for a, b in zip(s, vs)])
    return {"t": sol.t, "y1": y1, "s": s, "y1dot": v1, "sdot": vs, "D": D, "y2": y2, "Ft": Ft, "Fs": Fs}


# ---------------------------------------------------------------- energy accounting

def energy_terms(res: dict[str, np.ndarray], g: Gear, sink: float) -> dict[str, np.ndarray]:
    """Energy budget along a reduced-model run [J]. 'input' = initial kinetic energy - current kinetic
    energy + work of gravity and lift so far; 'stored' = gas + tyre + stop elastic energy gained since
    touchdown; 'dissipated' = D; 'residual' = input - stored - dissipated (zero for a correct model)."""
    cb, sb = np.cos(g.beta), np.sin(g.beta)
    v1, vs = res["y1dot"], res["sdot"]
    KE = 0.5 * g.M * v1 ** 2 + 0.5 * g.m * ((vs * sb) ** 2 + (v1 + vs * cb) ** 2)
    KE0 = 0.5 * (g.M + g.m) * sink ** 2
    y1_0, s0 = res["y1"][0], res["s"][0]
    y2_0 = res["y2"][0]
    W_grav = -g.M * G0 * (res["y1"] - y1_0) - g.m * G0 * (res["y2"] - y2_0)
    W_lift = g.lift_ratio * g.W * (res["y1"] - y1_0)
    d = np.maximum(g.R - res["y2"], 0.0)
    U = (gas_energy(res["s"], g) - gas_energy(s0, g) + 0.5 * g.k_t * d ** 2
         + 0.5 * g.k_stop * (np.maximum(-res["s"], 0) ** 2 - max(-s0, 0) ** 2))
    inp = KE0 - KE + W_grav + W_lift
    return {"input": inp, "stored": U, "dissipated": res["D"], "residual": inp - U - res["D"], "KE": KE}


def drop_metrics(res: dict[str, np.ndarray], g: Gear) -> dict[str, float]:
    """Headline numbers of a drop: peak tyre force and the gear load factor Ft_max / W, peak stroke,
    peak airframe deceleration in g, strut efficiency (work done by the strut force over the compression
    stroke divided by peak strut force times peak stroke), and whether the strut bottomed."""
    i_max = int(np.argmax(res["s"]))
    s, Fs = res["s"][: i_max + 1], res["Fs"][: i_max + 1]
    work = float(np.trapezoid(Fs, s))
    acc = np.gradient(res["y1dot"], res["t"])
    return {"Ft_max": float(res["Ft"].max()), "n_gear": float(res["Ft"].max() / g.W),
            "s_max": float(res["s"].max()), "n_airframe": float(acc.max() / G0),
            "Fs_max": float(Fs.max()), "efficiency": work / (float(Fs.max()) * (float(s.max()) - float(s[0]))),
            "bottomed": bool(res["s"].max() > g.stroke_max)}


# ---------------------------------------------------------------- formulation 2: constrained Newton-Euler

def _rot(th: float) -> np.ndarray:
    c, s = np.cos(th), np.sin(th)
    return np.array([[c, -s], [s, c]])


def _axis(g: Gear, th1: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Unit vector e along the strut (axle towards attachment), its normal n, and their th1-derivatives."""
    e0 = np.array([np.sin(g.beta), np.cos(g.beta)])
    n0 = np.array([np.cos(g.beta), -np.sin(g.beta)])
    R, Rp = _rot(th1), _rot(th1 + np.pi / 2)
    return R @ e0, R @ n0, Rp @ e0, Rp @ n0


def constraints(q: np.ndarray, g: Gear) -> np.ndarray:
    """Phi(q) = 0 for q = (x1, y1, th1, x2, y2, th2): x1 = 0 and th1 = 0 (the guide), th2 - th1 = 0 and
    (r2 - r1) . n(th1) = 0 (the translational joint along the strut)."""
    _, n, _, _ = _axis(g, q[2])
    d = q[3:5] - q[0:2]
    return np.array([q[0], q[2], q[5] - q[2], d @ n])


def constraint_jacobian(q: np.ndarray, g: Gear) -> np.ndarray:
    """Phi_q, 4 x 6."""
    _, n, _, n_p = _axis(g, q[2])
    d = q[3:5] - q[0:2]
    J = np.zeros((4, 6))
    J[0, 0] = 1.0
    J[1, 2] = 1.0
    J[2, 2], J[2, 5] = -1.0, 1.0
    J[3, 0:2], J[3, 2], J[3, 3:5] = -n, d @ n_p, n
    return J


def constraint_gamma(q: np.ndarray, qd: np.ndarray, g: Gear) -> np.ndarray:
    """Right-hand side gamma = -(Phi_q qd)_q qd of the acceleration constraint Phi_q q'' = gamma."""
    _, n, _, n_p = _axis(g, q[2])
    d, dd, w = q[3:5] - q[0:2], qd[3:5] - qd[0:2], qd[2]
    return np.array([0.0, 0.0, 0.0, -(2.0 * w * (dd @ n_p) - w * w * (d @ n))])


def stroke_full(q: np.ndarray, qd: np.ndarray, g: Gear) -> tuple[float, float]:
    """Stroke s = ell0 + (r2 - r1) . e and its rate, from full coordinates."""
    e, _, e_p, _ = _axis(g, q[2])
    d, dd = q[3:5] - q[0:2], qd[3:5] - qd[0:2]
    return float(g.ell0 + d @ e), float(dd @ e + qd[2] * (d @ e_p))


def applied_forces_full(q: np.ndarray, qd: np.ndarray, g: Gear) -> np.ndarray:
    """Generalized applied forces Q (6): gravity and lift on body 1, gravity and tyre on body 2, and the
    strut force F along the axis, +F e on body 1 and -F e on body 2, with its th1 term -F (d . e')."""
    e, _, e_p, _ = _axis(g, q[2])
    s, sd = stroke_full(q, qd, g)
    F = strut_force(s, sd, g)
    Ft, _ = tyre_force(q[4], qd[4], g)
    Q = np.zeros(6)
    Q[1] = -g.M * G0 + g.lift_ratio * g.W
    Q[4] = -g.m * G0 + Ft
    Q[0:2] += F * e
    Q[3:5] -= F * e
    Q[2] = -F * ((q[3:5] - q[0:2]) @ e_p)
    return Q


def accelerations_full(q: np.ndarray, qd: np.ndarray, g: Gear, alpha: float = 0.0,
                       beta_b: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """Solve [[Mq, -Phi_q^T], [Phi_q, 0]] [q''; lambda] = [Q; gamma - 2 alpha Phi' - beta_b^2 Phi].
    Returns (q'', lambda); the constraint forces on the bodies are Phi_q^T lambda."""
    Mq = np.diag([g.M, g.M, g.I1, g.m, g.m, g.I2])
    J = constraint_jacobian(q, g)
    rhs_c = constraint_gamma(q, qd, g) - 2.0 * alpha * (J @ qd) - beta_b ** 2 * constraints(q, g)
    A = np.zeros((10, 10))
    A[:6, :6], A[:6, 6:], A[6:, :6] = Mq, -J.T, J
    sol = np.linalg.solve(A, np.concatenate([applied_forces_full(q, qd, g), rhs_c]))
    return sol[:6], sol[6:]


def full_state_from_reduced(z: np.ndarray, g: Gear) -> np.ndarray:
    """Full state (q, qd, D) consistent with the constraints, from the reduced state (y1, s, y1', s', D)."""
    y1, s, v1, vs, D = z
    x2, y2 = axle_position(y1, s, g)
    q = np.array([0.0, y1, 0.0, x2, y2, 0.0])
    qd = np.array([0.0, v1, 0.0, vs * np.sin(g.beta), v1 + vs * np.cos(g.beta), 0.0])
    return np.concatenate([q, qd, [D]])


def simulate_full(g: Gear, sink: float, t_end: float = 0.5, n_out: int = 2001, alpha: float = 0.0,
                  beta_b: float = 0.0, rtol: float = 1e-10, atol: float = 1e-12, max_step: float = 2e-3,
                  dz0: np.ndarray | None = None) -> dict[str, np.ndarray]:
    """Integrate constrained Newton-Euler from touchdown (Baumgarte gains alpha, beta_b; zero = none).
    ``dz0`` (13 values) is added to the consistent initial state, to start with a constraint violation.
    Returns t, q (6 x N), qd, the multipliers lam (4 x N), the constraint violation phi (4 x N), the
    stroke s, the tyre force Ft and the bearing side load (the joint's normal force, lam[3])."""
    def rhs(t: float, z: np.ndarray) -> np.ndarray:
        q, qd = z[:6], z[6:12]
        qdd, _ = accelerations_full(q, qd, g, alpha, beta_b)
        s, sd = stroke_full(q, qd, g)
        return np.concatenate([qd, qdd, [dissipated_power(q[4], qd[4], sd, g)]])

    z0 = full_state_from_reduced(initial_state(g, sink), g)
    if dz0 is not None:
        z0 = z0 + np.asarray(dz0, float)
    sol = solve_ivp(rhs, (0.0, t_end), z0, method="DOP853", t_eval=np.linspace(0.0, t_end, n_out),
                    rtol=rtol, atol=atol, max_step=max_step)
    N = sol.t.size
    lam, phi, s = np.zeros((4, N)), np.zeros((4, N)), np.zeros(N)
    Ft = np.zeros(N)
    for k in range(N):
        q, qd = sol.y[:6, k], sol.y[6:12, k]
        lam[:, k] = accelerations_full(q, qd, g, alpha, beta_b)[1]
        phi[:, k] = constraints(q, g)
        s[k] = stroke_full(q, qd, g)[0]
        Ft[k] = tyre_force(q[4], qd[4], g)[0]
    return {"t": sol.t, "q": sol.y[:6], "qd": sol.y[6:12], "D": sol.y[12], "lam": lam, "phi": phi,
            "s": s, "Ft": Ft, "side_load": lam[3]}


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save landing_gear.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    g = Gear()
    sink = 3.05                                    # m/s (10 ft/s, a typical limit sink rate)
    print(f"gear: M = {g.M:.0f} kg, m = {g.m:.0f} kg, rake {np.degrees(g.beta):.0f} deg, W = {g.W:.0f} N, "
          f"lift = {g.lift_ratio:.1f} W, sink {sink} m/s, KE = {0.5 * (g.M + g.m) * sink ** 2:.0f} J")
    print(f"strut: gas force at full extension {gas_force(0.0, g):.0f} N, at {g.stroke_max} m {gas_force(g.stroke_max, g):.0f} N; "
          f"hanging stroke s0 = {initial_stroke(g) * 1e3:.4f} mm")
    s_static = brentq(lambda s: gas_force(s, g) - (g.M * G0) * np.cos(g.beta), 0.0, g.stroke_max)
    print(f"static stroke on the ground (no lift): {s_static * 1e3:.1f} mm")

    res = simulate_reduced(g, sink)
    met = drop_metrics(res, g)
    en = energy_terms(res, g, sink)
    print(f"\nLagrange (reduced) drop: peak tyre force {met['Ft_max']:.0f} N, gear load factor {met['n_gear']:.3f}, "
          f"peak stroke {met['s_max'] * 1e3:.1f} mm (bottomed: {met['bottomed']}), peak airframe deceleration "
          f"{met['n_airframe']:.3f} g, strut efficiency {met['efficiency']:.3f}")
    i = int(np.argmax(res["s"]))
    print(f"  at peak stroke (t = {res['t'][i]:.4f} s): input {en['input'][i]:.1f} J = stored {en['stored'][i]:.1f} J "
          f"+ dissipated {en['dissipated'][i]:.1f} J; residual {en['residual'][i]:.2e} J")
    print(f"  max |residual| over the run {np.max(np.abs(en['residual'])):.2e} J "
          f"({np.max(np.abs(en['residual'])) / en['input'].max():.1e} of the largest input)")
    print(f"  at t = {res['t'][-1]:.2f} s: stored {en['stored'][-1]:.1f} J, dissipated {en['dissipated'][-1]:.1f} J, "
          f"stroke {res['s'][-1] * 1e3:.1f} mm")

    full = simulate_full(g, sink, n_out=1001)
    dy1 = np.max(np.abs(full["q"][1] - res["y1"][::2]))
    ds = np.max(np.abs(full["s"] - res["s"][::2]))
    print(f"\nNewton-Euler (full, 6 coordinates, 4 constraints) against Lagrange: max |y1 diff| {dy1:.2e} m, "
          f"max |s diff| {ds:.2e} m, peak tyre force {full['Ft'].max():.0f} N")
    print(f"  constraint drift max |Phi| from a consistent start: {np.max(np.abs(full['phi'])):.2e} "
          f"(explicit Runge-Kutta preserves the constraints, linear along this motion, to round-off)")
    dz = np.zeros(13)
    dz[3], dz[9] = 1.0e-3, 0.05            # start the axle 1 mm and 0.05 m/s off the strut axis
    for a_b in (0.0, 50.0):
        pert = simulate_full(g, sink, t_end=0.2, n_out=201, alpha=a_b, beta_b=a_b, dz0=dz)
        print(f"  perturbed start, Baumgarte alpha = beta = {a_b:.0f}/s: |Phi_4| at 0, 0.1, 0.2 s = "
              + ", ".join(f"{abs(pert['phi'][3][j]) * 1e3:.3f}" for j in (0, 100, 200)) + " mm")
    k = int(np.argmax(full["Ft"]))
    print(f"  strut bearing side load at peak tyre force: {full['side_load'][k]:.0f} N "
          f"(tyre force x sin(beta) = {full['Ft'][k] * np.sin(g.beta):.0f} N); guide reaction "
          f"{full['lam'][0][k]:.0f} N")

    print("\nsink-rate sweep (Lagrange, first 0.3 s): sink [m/s], n_gear, peak stroke [mm], efficiency, max |residual| [J]")
    for v in (1.5, 2.5, 3.05, 4.0):
        r = simulate_reduced(g, v, t_end=0.3, n_out=601, rtol=1e-9, atol=1e-11)
        mm = drop_metrics(r, g)
        print(f"  {v:4.2f}   {mm['n_gear']:.3f}   {mm['s_max'] * 1e3:6.1f}{' (bottomed)' if mm['bottomed'] else ''}   "
              f"{mm['efficiency']:.3f}   {np.max(np.abs(energy_terms(r, g, v)['residual'])):.1e}")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    ax[0].plot(res["t"], res["Ft"] / 1e3, label="tyre force (Lagrange)")
    ax[0].plot(full["t"][::20], full["Ft"][::20] / 1e3, "o", ms=3, label="tyre force (Newton-Euler)")
    ax[0].plot(res["t"], res["Fs"] / 1e3, label="strut force")
    ax[0].plot(full["t"], full["side_load"] / 1e3, label="bearing side load (lambda_4)")
    ax[0].set_xlabel("time from touchdown (s)")
    ax[0].set_ylabel("force (kN)")
    ax[0].set_title(f"drop at {sink} m/s, lift = weight")
    ax[0].legend(fontsize=7)
    ax[1].plot(res["s"][: i + 1] * 1e3, res["Fs"][: i + 1] / 1e3, label="strut force")
    ax[1].plot(res["s"][: i + 1] * 1e3, gas_force(res["s"][: i + 1], g) / 1e3, "--", label="gas force alone")
    ax[1].set_xlabel("stroke (mm)")
    ax[1].set_ylabel("force (kN)")
    ax[1].set_title(f"compression stroke, efficiency {met['efficiency']:.2f}")
    ax[1].legend(fontsize=7)
    ax[2].plot(res["t"], en["input"] / 1e3, label="input (KE lost + work of gravity and lift)")
    ax[2].plot(res["t"], en["stored"] / 1e3, label="stored (gas, tyre, stop)")
    ax[2].plot(res["t"], en["dissipated"] / 1e3, label="dissipated (orifice, tyre)")
    ax[2].plot(res["t"], (en["stored"] + en["dissipated"]) / 1e3, "k:", label="stored + dissipated")
    ax[2].set_xlabel("time from touchdown (s)")
    ax[2].set_ylabel("energy (kJ)")
    ax[2].set_title("energy accounting")
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out_dir / "landing_gear.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'landing_gear.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
