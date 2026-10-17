"""P05.04 Frames, Machines & Friction: reference solution.

An elevator control run is a machine: a control stick pivoted on the cockpit floor drives a pushrod, the
pushrod turns a bellcrank, the bellcrank drives a second pushrod, and that pushrod turns the elevator by
its control horn. The pilot's question is how hard the stick pulls; the designer's is how the gearing
varies over the travel and how much friction adds. This module answers both:

* kinematics of the linkage: every pushrod has a fixed length, so each joint lies on the intersection of
  two circles; the stick angle fixes the bellcrank angle and the elevator angle;
* statics member by member (the method of frames and machines): pushrods are two-force members; the
  stick, the bellcrank and the elevator are multi-force members, each in moment equilibrium about its
  pivot; the result is the stick force for a given aerodynamic hinge moment;
* an independent check by virtual work: with no friction, F_grip * L_grip * dphi = -H * d(delta), with
  the gearing d(delta)/d(phi) measured by differencing the kinematics alone;
* friction in pins by the friction circle (moment r_f |R| against the rotation, r_f = r mu / sqrt(1 + mu^2)),
  solved by fixed-point iteration because the reactions depend on the friction; the force to drive the
  stick and the force to hold it against the air differ by a friction band;
* friction in cables: the capstan law T2 / T1 = exp(mu theta) for a cable sliding over a fixed fairlead
  (checked against a numerical integration), and the exact ratio for a pulley turning on a pin.

Coordinates: side view, x aft, y up. Angles anticlockwise positive. Stick angle phi > 0: grip forward.
Elevator angle delta > 0: trailing edge up. Units: SI (N, m, N m, rad).

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes control_linkage.png)
Test: pytest curriculum/phase-05/lesson-04-frames-machines-friction
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


def cross2(a: np.ndarray, b: np.ndarray) -> float:
    """z component of a x b for planar vectors."""
    return float(a[0] * b[1] - a[1] * b[0])


def unit(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def circle_intersections(c0: np.ndarray, r0: float, c1: np.ndarray, r1: float) -> np.ndarray:
    """The two intersection points of the circles (c0, r0) and (c1, r1), shape (2, 2). Raises ValueError
    when the circles do not meet (the linkage cannot be assembled in that position)."""
    c0, c1 = np.asarray(c0, float), np.asarray(c1, float)
    d = np.linalg.norm(c1 - c0)
    if d == 0 or d > r0 + r1 or d < abs(r0 - r1):
        raise ValueError("circles do not intersect: linkage cannot reach this position")
    a = (r0 ** 2 - r1 ** 2 + d ** 2) / (2 * d)
    h = np.sqrt(max(r0 ** 2 - a ** 2, 0.0))
    e = (c1 - c0) / d
    p = c0 + a * e
    n = np.array([-e[1], e[0]])
    return np.array([p + h * n, p - h * n])


# ---------------------------------------------------------------- the linkage

@dataclass
class Linkage:
    """Stick - pushrod 1 - bellcrank - pushrod 2 - horn. Lengths in m, angles in rad.

    Stick pivot S0; grip at Lg along the stick (up at phi = 0), pushrod 1 pinned a_s below the pivot.
    Bellcrank pivot C0; arm 1 (length b1, angle alpha1 at psi = 0) takes pushrod 1, arm 2 (b2, alpha2)
    takes pushrod 2. Elevator hinge E0; horn of length h pointing forward (-x) at delta = 0. Pushrod
    lengths are set so that phi = psi = delta = 0 is the neutral position.
    k_hinge: aerodynamic hinge-moment stiffness (N m/rad): H = -k_hinge * delta. e_cp: chordwise distance of
    the elevator's aerodynamic resultant aft of the hinge (m)."""
    S0: np.ndarray = field(default_factory=lambda: np.array([0.0, 0.0]))
    Lg: float = 0.55
    a_s: float = 0.08
    C0: np.ndarray = field(default_factory=lambda: np.array([1.20, 0.0]))
    b1: float = 0.08
    alpha1: float = -np.pi / 2
    b2: float = 0.10
    alpha2: float = 0.0
    E0: np.ndarray = field(default_factory=lambda: np.array([1.38, 0.60]))
    h: float = 0.08
    k_hinge: float = 82.7
    e_cp: float = 0.075

    def __post_init__(self) -> None:
        self.L1 = float(np.linalg.norm(self.bellcrank_points(0.0)[0] - self.stick_points(0.0)[1]))
        self.L2 = float(np.linalg.norm(self.horn_point(0.0) - self.bellcrank_points(0.0)[1]))

    def stick_points(self, phi: float) -> tuple[np.ndarray, np.ndarray]:
        """(grip, lower pushrod pin) at stick angle phi."""
        s, c = np.sin(phi), np.cos(phi)
        return self.S0 + self.Lg * np.array([-s, c]), self.S0 + self.a_s * np.array([s, -c])

    def bellcrank_points(self, psi: float) -> tuple[np.ndarray, np.ndarray]:
        """(arm-1 pin, arm-2 pin) at bellcrank angle psi."""
        return (self.C0 + self.b1 * np.array([np.cos(psi + self.alpha1), np.sin(psi + self.alpha1)]),
                self.C0 + self.b2 * np.array([np.cos(psi + self.alpha2), np.sin(psi + self.alpha2)]))

    def horn_point(self, delta: float) -> np.ndarray:
        """Horn pin at elevator angle delta."""
        return self.E0 + self.h * np.array([-np.cos(delta), -np.sin(delta)])


def _angle_on(center: np.ndarray, p: np.ndarray, offset: float) -> float:
    return float(np.arctan2(p[1] - center[1], p[0] - center[0]) - offset)


def solve_kinematics(lk: Linkage, phi: float) -> tuple[float, float]:
    """Bellcrank angle psi and elevator angle delta for stick angle phi, by two circle intersections.
    Of the two assembly solutions, the one nearer the neutral position is taken."""
    _, Ps = lk.stick_points(phi)
    cand = circle_intersections(lk.C0, lk.b1, Ps, lk.L1)
    B1 = cand[np.argmin(np.linalg.norm(cand - lk.bellcrank_points(0.0)[0], axis=1))]
    psi = _angle_on(lk.C0, B1, lk.alpha1)
    psi = (psi + np.pi) % (2 * np.pi) - np.pi
    _, B2 = lk.bellcrank_points(psi)
    cand = circle_intersections(lk.E0, lk.h, B2, lk.L2)
    H = cand[np.argmin(np.linalg.norm(cand - lk.horn_point(0.0), axis=1))]
    delta = _angle_on(lk.E0, H, np.pi)
    delta = (delta + np.pi) % (2 * np.pi) - np.pi
    return float(psi), float(delta)


def gearing(lk: Linkage, phi: float, dphi: float = 1e-6) -> tuple[float, float]:
    """(d psi / d phi, d delta / d phi) by central differences of ``solve_kinematics``: the velocity
    ratios of the bellcrank and the elevator, from the geometry alone."""
    p1, d1 = solve_kinematics(lk, phi + dphi)
    p0, d0 = solve_kinematics(lk, phi - dphi)
    return (p1 - p0) / (2 * dphi), (d1 - d0) / (2 * dphi)


def hinge_moment(lk: Linkage, delta: float) -> float:
    """Aerodynamic moment on the elevator about its hinge (N m, anticlockwise positive): -k_hinge delta,
    a restoring moment that grows with deflection (illustrative linear model)."""
    return -lk.k_hinge * delta


def friction_radius(r_pin: float, mu: float) -> float:
    """Radius of the friction circle of a pin of radius r_pin: r sin(atan mu) = r mu / sqrt(1 + mu^2).
    The pin reaction is tangent to this circle when the joint slips, so the friction moment is r_f |R|."""
    return r_pin * mu / np.sqrt(1.0 + mu ** 2)


def linkage_forces(lk: Linkage, phi: float, mu: float = 0.0, r_pin: float = 0.003175, direction: int = 1,
                   tol: float = 1e-12, max_iter: int = 200) -> dict:
    """Member-by-member statics of the linkage at stick angle phi.

    Free bodies, each in moment equilibrium about its own pivot:
      elevator:  H(delta) + P2 c_e + M_f,e = 0
      bellcrank: P1 c_b1 + P2 c_b2 + M_f,b = 0
      stick:     F c_g + P1 c_s + M_f,s = 0
    P1, P2 are pushrod tensions; c_* are the moments of unit forces along the pushrods (and of a unit grip
    force perpendicular to the stick, in the direction of increasing phi) about each pivot. Pin friction
    moments M_f = -sign(omega) r_f |R| oppose each body's rotation; ``direction`` is the sign of
    d(phi)/dt (+1: stick moving toward +phi). The reactions R depend on the forces, so the system is
    solved by fixed-point iteration from the frictionless answer.

    Returns P1, P2, F (grip force, N, positive toward +phi), the pin reactions R_stick, R_bell, R_hinge
    (vectors on each body), the friction moments, psi, delta and the number of iterations."""
    psi, delta = solve_kinematics(lk, phi)
    dpsi, ddelta = gearing(lk, phi)
    grip, Ps = lk.stick_points(phi)
    B1, B2 = lk.bellcrank_points(psi)
    H = lk.horn_point(delta)
    u1 = unit(B1 - Ps)                      # pushrod 1 tension pulls the stick pin toward B1
    u2 = unit(H - B2)                       # pushrod 2 tension pulls B2 toward H
    tg = np.array([-np.cos(phi), -np.sin(phi)])  # grip moves this way for increasing phi
    c_e = cross2(H - lk.E0, -u2)
    c_b1, c_b2 = cross2(B1 - lk.C0, -u1), cross2(B2 - lk.C0, u2)
    c_s, c_g = cross2(Ps - lk.S0, u1), cross2(grip - lk.S0, tg)
    Hm = hinge_moment(lk, delta)
    chord = np.array([np.cos(delta), np.sin(delta)])
    normal = np.array([-np.sin(delta), np.cos(delta)])
    F_aero = (Hm / lk.e_cp) * normal       # resultant whose moment about E0 is Hm
    rf = friction_radius(r_pin, mu)
    sgn = direction * np.sign([1.0, dpsi, ddelta])   # rotation sense of stick, bellcrank, elevator
    Mf = np.zeros(3)
    it = 0
    for it in range(1, max_iter + 1):
        P2 = -(Hm + Mf[2]) / c_e
        P1 = -(P2 * c_b2 + Mf[1]) / c_b1
        F = -(P1 * c_s + Mf[0]) / c_g
        R_hinge = -(F_aero - P2 * u2)
        R_bell = -(-P1 * u1 + P2 * u2)
        R_stick = -(F * tg + P1 * u1)
        new = -sgn * rf * np.array([np.linalg.norm(R_stick), np.linalg.norm(R_bell), np.linalg.norm(R_hinge)])
        if np.max(np.abs(new - Mf)) <= tol * max(1.0, abs(Hm)):
            Mf = new
            break
        Mf = new
    else:
        raise RuntimeError("friction iteration did not converge")
    P2 = -(Hm + Mf[2]) / c_e
    P1 = -(P2 * c_b2 + Mf[1]) / c_b1
    F = -(P1 * c_s + Mf[0]) / c_g
    return {"P1": float(P1), "P2": float(P2), "F": float(F), "R_stick": -(F * tg + P1 * u1),
            "R_bell": -(-P1 * u1 + P2 * u2), "R_hinge": -(F_aero - P2 * u2), "M_friction": Mf,
            "psi": psi, "delta": delta, "hinge_moment": Hm, "iterations": it,
            "_parts": {"grip": grip, "Ps": Ps, "B1": B1, "B2": B2, "H": H, "u1": u1, "u2": u2, "tg": tg,
                       "F_aero": F_aero}}


def stick_force_virtual_work(lk: Linkage, phi: float) -> float:
    """Frictionless grip force by virtual work: F Lg dphi + H d(delta) = 0, so
    F = -H(delta) (d delta / d phi) / Lg, with the gearing from the kinematics alone."""
    _, delta = solve_kinematics(lk, phi)
    _, ddelta = gearing(lk, phi)
    return -hinge_moment(lk, delta) * ddelta / lk.Lg


def friction_power_balance(lk: Linkage, phi: float, sol: dict, mu: float, r_pin: float = 0.003175,
                           direction: int = 1) -> float:
    """Energy check of a friction solution, per unit d(phi)/dt (N m):
    F Lg + H d(delta)/d(phi) - sum r_f |R_i| |d theta_i / d phi|, which must vanish (input power =
    power into the air + power lost in the pins)."""
    dpsi, ddelta = gearing(lk, phi)
    rf = friction_radius(r_pin, mu)
    loss = rf * (np.linalg.norm(sol["R_stick"]) * 1.0 + np.linalg.norm(sol["R_bell"]) * abs(dpsi)
                 + np.linalg.norm(sol["R_hinge"]) * abs(ddelta))
    return float(direction * sol["F"] * lk.Lg + direction * sol["hinge_moment"] * ddelta - loss)


def mechanical_advantage(lk: Linkage, phi: float) -> float:
    """Ideal mechanical advantage (H / h_grip-moment): elevator hinge moment resisted per unit of moment
    applied at the stick, = 1 / |d delta / d phi| for a frictionless linkage (dimensionless)."""
    return 1.0 / abs(gearing(lk, phi)[1])


def efficiency(lk: Linkage, phi: float, mu: float, r_pin: float = 0.003175) -> float:
    """Mechanical efficiency when the pilot drives the elevator against the air: F_ideal / F_with_friction."""
    f0 = linkage_forces(lk, phi, 0.0, r_pin)["F"]
    dirn = int(np.sign(phi)) if phi != 0 else 1
    f1 = linkage_forces(lk, phi, mu, r_pin, direction=dirn)["F"]
    return float(f0 / f1)


# ---------------------------------------------------------------- cables

def capstan_ratio(mu: float, wrap: float) -> float:
    """Tight-to-slack tension ratio of a cable slipping over a fixed curved surface: exp(mu theta)."""
    return float(np.exp(mu * wrap))


def capstan_numeric(mu: float, wrap: float) -> float:
    """The same ratio by integrating dT/d(theta) = mu T (from the equilibrium of a cable element: the normal
    force T d(theta) and the friction mu T d(theta)) with solve_ivp. An independent check."""
    sol = solve_ivp(lambda th, T: mu * T, (0.0, wrap), [1.0], rtol=1e-11, atol=1e-13)
    return float(sol.y[0, -1])


def pulley_ratio(mu_b: float, r_b: float, R: float, wrap: float) -> float:
    """Tight-to-slack tension ratio k of a cable turning a pulley of radius R on a pin of radius r_b.

    Moment equilibrium of the pulley about its centre: (k - 1) T R = r_f |T_tight + T_slack|, with
    |T_t + T_s|^2 = T_t^2 + T_s^2 - 2 T_t T_s cos(theta) for a wrap angle theta. Squared, this is the
    quadratic (1 - rho^2) k^2 - 2 (1 - rho^2 cos theta) k + (1 - rho^2) = 0 with rho = r_f / R; the root
    above 1 is returned. At theta = pi it reduces to (1 + rho) / (1 - rho)."""
    rho = friction_radius(r_b, mu_b) / R
    if rho >= 1.0:
        raise ValueError("friction circle larger than the pulley: it cannot turn")
    a = 1.0 - rho ** 2
    b = 1.0 - rho ** 2 * np.cos(wrap)
    return float((b + np.sqrt(max(b * b - a * a, 0.0))) / a)


def cable_run_ratio(elements: list[tuple]) -> float:
    """Overall tension ratio of a cable run through fairleads and pulleys in series:
    ('fairlead', mu, wrap) or ('pulley', mu_b, r_b, R, wrap). The run's efficiency is the inverse."""
    k = 1.0
    for e in elements:
        if e[0] == "fairlead":
            k *= capstan_ratio(e[1], e[2])
        elif e[0] == "pulley":
            k *= pulley_ratio(*e[1:])
        else:
            raise ValueError(f"unknown element {e[0]!r}")
    return k


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save control_linkage.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    lk = Linkage()
    q, Se, ce, chd = 0.5 * 1.225 * 50.0 ** 2, 1.2, 0.30, -0.15
    print(f"hinge stiffness: q = {q:.1f} Pa, S_e = {Se} m^2, c_e = {ce} m, Ch_delta = {chd}/rad -> "
          f"k = {-chd * q * Se * ce:.1f} N m/rad (illustrative)")
    print(f"pushrod lengths L1 = {lk.L1:.3f} m, L2 = {lk.L2:.3f} m")
    print(f"small-angle gearing -(a_s / b1)(b2 / h) = {-(lk.a_s / lk.b1) * (lk.b2 / lk.h):.4f}; "
          f"at neutral: d psi/d phi = {gearing(lk, 0.0)[0]:.4f}, d delta/d phi = {gearing(lk, 0.0)[1]:.4f}")

    phis = np.radians([-20.0, -10.0, 0.0, 10.0, 20.0])
    print("\nphi (deg), delta (deg), gearing, hinge moment (N m), F chain (N), F virtual work (N), P1, P2 (N)")
    for ph in phis:
        s = linkage_forces(lk, ph)
        print(f"  {np.degrees(ph):6.1f} {np.degrees(s['delta']):8.3f} {gearing(lk, ph)[1]:8.4f} {s['hinge_moment']:8.3f} "
              f"{s['F']:9.4f} {stick_force_virtual_work(lk, ph):9.4f} {s['P1']:9.2f} {s['P2']:9.2f}")

    mu, r_pin = 0.2, 0.003175
    ph = np.radians(20.0)
    s0 = linkage_forces(lk, ph)
    sd = linkage_forces(lk, ph, mu, r_pin, direction=+1)
    sb = linkage_forces(lk, ph, mu, r_pin, direction=-1)
    print(f"\npins mu = {mu}, r = {r_pin * 1e3:.3f} mm (friction radius {friction_radius(r_pin, mu) * 1e3:.4f} mm) at phi = 20 deg:")
    print(f"  frictionless F = {s0['F']:.3f} N; driving F = {sd['F']:.3f} N ({sd['iterations']} iterations); "
          f"holding/back-driven F = {sb['F']:.3f} N; band {sd['F'] - sb['F']:.3f} N; efficiency {s0['F'] / sd['F']:.4f}")
    print(f"  reactions |R| stick {np.linalg.norm(sd['R_stick']):.1f}, bellcrank {np.linalg.norm(sd['R_bell']):.1f}, "
          f"hinge {np.linalg.norm(sd['R_hinge']):.1f} N; friction moments {sd['M_friction'].round(4)} N m")
    print(f"  power balance residual {friction_power_balance(lk, ph, sd, mu, r_pin):.1e} N m")
    for m_ in (0.1, 0.2, 0.4, 0.8):
        print(f"  mu = {m_}: efficiency {efficiency(lk, ph, m_, r_pin):.4f}")
    big = 0.0127
    print(f"  with 1/2-in (12.7 mm) pins at mu = 0.2: efficiency {efficiency(lk, ph, 0.2, big / 2):.4f}")

    # ---- cables
    print(f"\ncapstan: mu = 0.2 over 180 deg: exp = {capstan_ratio(0.2, np.pi):.6f}, ODE = {capstan_numeric(0.2, np.pi):.6f}")
    print(f"fairlead, mu = 0.2, 3 deg deflection: {capstan_ratio(0.2, np.radians(3.0)):.5f}")
    for name, mub, rb in (("plain bushing", 0.15, 0.005), ("ball bearing (equivalent)", 0.005, 0.005)):
        print(f"pulley R = 38.1 mm, {name}: 90 deg {pulley_ratio(mub, rb, 0.0381, np.pi / 2):.5f}, "
              f"180 deg {pulley_ratio(mub, rb, 0.0381, np.pi):.5f}")
    run_plain = [("pulley", 0.15, 0.005, 0.0381, np.radians(a)) for a in (90, 45, 30, 90)] + \
                [("fairlead", 0.2, np.radians(3.0))] * 3
    run_ball = [("pulley", 0.005, 0.005, 0.0381, np.radians(a)) for a in (90, 45, 30, 90)] + \
               [("fairlead", 0.2, np.radians(3.0))] * 3
    for name, run in (("plain-bushing pulleys", run_plain), ("ball-bearing pulleys", run_ball)):
        k = cable_run_ratio(run)
        print(f"cable run (4 pulleys 90/45/30/90 deg + 3 fairleads at 3 deg), {name}: ratio {k:.4f}, efficiency {1 / k:.4f}")
    print(f"cable wrapped 180 deg on a seized 38.1 mm pulley (mu = 0.3): {capstan_ratio(0.3, np.pi):.3f}")

    # ---- sweep and figure
    sweep = np.radians(np.linspace(-20.0, 20.0, 41))
    dl, gr, F0, Fd, Fb = [], [], [], [], []
    for p in sweep:
        s = linkage_forces(lk, p)
        dl.append(np.degrees(s["delta"])); gr.append(gearing(lk, p)[1]); F0.append(s["F"])
        Fd.append(linkage_forces(lk, p, mu, r_pin, direction=+1)["F"])
        Fb.append(linkage_forces(lk, p, mu, r_pin, direction=-1)["F"])
    gr = np.array(gr)
    print(f"\nsweep phi -20..20 deg: delta {min(dl):.2f}..{max(dl):.2f} deg, gearing {gr.min():.4f}..{gr.max():.4f}")

    fig, axs = plt.subplots(1, 3, figsize=(14, 4.4))
    ax = axs[0]
    s = linkage_forces(lk, np.radians(15.0))
    pp = s["_parts"]
    for a, b in ((lk.S0, pp["grip"]), (lk.S0, pp["Ps"]), (pp["Ps"], pp["B1"]), (lk.C0, pp["B1"]), (lk.C0, pp["B2"]),
                 (pp["B2"], pp["H"]), (lk.E0, pp["H"])):
        ax.plot([a[0], b[0]], [a[1], b[1]], "k-", lw=2)
    ch = np.array([np.cos(s["delta"]), np.sin(s["delta"])])
    ax.plot([lk.E0[0], lk.E0[0] + 0.3 * ch[0]], [lk.E0[1], lk.E0[1] + 0.3 * ch[1]], "b-", lw=4)
    for p in (lk.S0, lk.C0, lk.E0):
        ax.plot(*p, "ro")
    ax.set_aspect("equal"); ax.set_title("linkage at phi = 15 deg"); ax.set_xlabel("x aft (m)")
    ax = axs[1]
    ax.plot(np.degrees(sweep), dl, label="delta (deg)")
    ax.plot(np.degrees(sweep), gr * 10, label="10 x gearing")
    ax.set_xlabel("stick angle phi (deg)"); ax.legend(); ax.set_title("kinematics")
    ax = axs[2]
    ax.plot(np.degrees(sweep), F0, label="frictionless")
    ax.plot(np.degrees(sweep), Fd, "--", label="driving, mu = 0.2")
    ax.plot(np.degrees(sweep), Fb, ":", label="back-driven, mu = 0.2")
    ax.set_xlabel("stick angle phi (deg)"); ax.set_ylabel("grip force (N)"); ax.legend(); ax.set_title("stick force")
    fig.tight_layout()
    fig.savefig(out_dir / "control_linkage.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
