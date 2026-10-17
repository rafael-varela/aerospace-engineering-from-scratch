"""P04.04 Rotation & Angular Momentum: starter code.

Rigid-body rotation from the moment of inertia to the full three-dimensional equations of motion:

* moments of inertia of simple shapes and the parallel-axis theorem (tensor form);
* Euler's equations I dw/dt + w x (I w) = tau in the body frame, with the attitude carried by a unit
  quaternion, integrated with `scipy.integrate.solve_ivp`;
* torque-free motion: conservation of the inertial angular momentum and of the kinetic energy, the
  body-frame nutation of a symmetric body, and the instability of rotation about the intermediate axis;
* a spacecraft with a reaction wheel (a gyrostat): the wheel stores angular momentum, the body slews,
  the total is conserved; the bang-bang slew in closed form;
* the heavy symmetric top on a fixed pivot: conservation of energy, of the vertical component of L and
  of the spin component, steady precession at the rate that solves
  Omega (I3 w3 - I1 Omega cos(theta)) = m g l, and the fast-top limit Omega = m g l / (I3 w3).

Conventions: body axes 1, 2, 3 are principal axes; the body-to-inertial rotation matrix R(q) maps body
components to inertial components; q = [w, x, y, z] with scalar first; the inertial z axis points up.

Run:  python main.py [out_dir]     (prints the tables, writes spinning_top.png)
Test: pytest curriculum/phase-04/lesson-04-rotation-angular-momentum
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

G0 = 9.80665


# ---------------------------------------------------------------- moments of inertia

def inertia_solid_cylinder(m: float, r: float, length: float) -> np.ndarray:
    """Inertia tensor (3 x 3, kg m^2) of a uniform solid cylinder about its centre of mass, axis along 3:
    I3 = m r^2 / 2, I1 = I2 = m (3 r^2 + L^2) / 12. A thin disk is L -> 0, a thin rod r -> 0."""
    raise NotImplementedError


def inertia_box(m: float, a: float, b: float, c: float) -> np.ndarray:
    """Inertia tensor of a uniform box with side lengths a, b, c along axes 1, 2, 3 about its centre:
    diag(m (b^2 + c^2), m (a^2 + c^2), m (a^2 + b^2)) / 12."""
    raise NotImplementedError


def parallel_axis(I_cm: np.ndarray, m: float, d: np.ndarray) -> np.ndarray:
    """Inertia tensor about a point displaced by d from the centre of mass: I = I_cm + m (|d|^2 E - d d^T)."""
    raise NotImplementedError


# ---------------------------------------------------------------- vectors and quaternions (given)

def cross3(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Cross product of two 3-vectors, written out (several times faster than np.cross on single vectors)."""
    return np.array([a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]])


def quat_multiply(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Hamilton product p * q of quaternions [w, x, y, z]."""
    w1, x1, y1, z1 = p
    w2, x2, y2, z2 = q
    return np.array([w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
                     w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                     w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
                     w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2])


def quat_to_matrix(q: np.ndarray) -> np.ndarray:
    """Rotation matrix R (body to inertial) of a quaternion [w, x, y, z]; q is normalised first."""
    w, x, y, z = np.asarray(q, dtype=float) / np.linalg.norm(q)
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
                     [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
                     [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]])


def quat_from_axis_angle(axis: np.ndarray, angle: float) -> np.ndarray:
    """Unit quaternion of a rotation by `angle` (rad) about `axis`."""
    a = np.asarray(axis, dtype=float)
    a = a / np.linalg.norm(a)
    return np.concatenate([[math.cos(angle / 2.0)], math.sin(angle / 2.0) * a])


# ---------------------------------------------------------------- energy and angular momentum

def rotational_energy(I: np.ndarray, w: np.ndarray) -> np.ndarray | float:
    """T = w . (I w) / 2 for one body rate (3,) or a stack (n, 3)."""
    raise NotImplementedError


def angular_momentum_inertial(q: np.ndarray, w: np.ndarray, I: np.ndarray, h_wheel: np.ndarray | None = None) -> np.ndarray:
    """Inertial components of the angular momentum, R(q) (I w + h_wheel), for one state (q (4,), w (3,))
    or stacks (n, 4) and (n, 3). h_wheel is the wheels' angular momentum in body axes (default zero)."""
    raise NotImplementedError


# ---------------------------------------------------------------- Euler's equations and the attitude

def euler_rates(w: np.ndarray, I: np.ndarray, torque: np.ndarray, h_wheel: np.ndarray | None = None,
                h_wheel_dot: np.ndarray | None = None) -> np.ndarray:
    """Body angular acceleration from Euler's equations for a rigid body carrying wheels:

        I dw/dt = tau - w x (I w + h) - dh/dt,

    with h the wheels' angular momentum in body axes and dh/dt = the torque the body applies to the wheels.
    Without wheels this is I dw/dt + w x (I w) = tau."""
    raise NotImplementedError


def quaternion_rate(q: np.ndarray, w: np.ndarray) -> np.ndarray:
    """dq/dt = q * [0, w] / 2 for body rates w (the attitude kinematics)."""
    raise NotImplementedError


TorqueFn = Callable[[float, np.ndarray, np.ndarray], np.ndarray]


@dataclass
class Rotation:
    t: np.ndarray
    q: np.ndarray            # (n, 4)
    w: np.ndarray            # (n, 3) body rates
    h: np.ndarray            # (n, 3) wheel momentum in body axes


def simulate_rigid_body(I: np.ndarray, w0: np.ndarray, q0: np.ndarray, t_end: float,
                        torque: TorqueFn | None = None, wheel_torque: TorqueFn | None = None,
                        h0: np.ndarray | None = None, rtol: float = 1e-11, n_samples: int = 2001,
                        max_step: float = np.inf) -> Rotation:
    """Integrate the state y = [q (4), w (3), h (3)] from (q0, w0, h0) to t_end with `solve_ivp`
    (DOP853, rtol given, atol 1e-12), sampled at n_samples equally spaced times:

        dq/dt = quaternion_rate(q, w),
        dw/dt = euler_rates(w, I, tau(t, q, w), h, tau_w(t, q, w)),
        dh/dt = tau_w(t, q, w)      (torque applied to the wheels by the body, body axes).

    `torque` is the external torque in body axes; both default to zero. The quaternion is renormalised
    on output (the integrator keeps |q| = 1 to about rtol).
    """
    raise NotImplementedError


# ---------------------------------------------------------------- torque-free motion: closed forms

def symmetric_nutation_rate(I1: float, I3: float, w3: float) -> float:
    """Body-frame nutation rate of a torque-free symmetric body (I1 = I2): w1 + i w2 rotates at
    lambda = (I3 - I1) w3 / I1 (positive: same sense as the spin, for an oblate body)."""
    raise NotImplementedError


def intermediate_axis_growth_rate(I: np.ndarray, w2: float) -> float:
    """Exponential growth rate of small perturbations of a spin w2 about the intermediate principal axis
    (I1 < I2 < I3, axes in that order): sigma = |w2| sqrt((I2 - I1)(I3 - I2) / (I1 I3))."""
    raise NotImplementedError


def principal_axis_wobble_frequency(I: np.ndarray, axis: int, w: float) -> float:
    """Frequency of small oscillations of the rate about a spin w about principal axis `axis` (0 or 2,
    the extreme axes): sqrt((Ia - Ij)(Ia - Ik) / (Ij Ik)) |w| with j, k the other two axes. Raises
    ValueError for the intermediate axis, where the motion is unstable."""
    raise NotImplementedError


# ---------------------------------------------------------------- reaction wheel slew (closed form)

def bang_bang_slew(I_axis: float, angle: float, duration: float) -> tuple[float, float, float]:
    """Rest-to-rest single-axis slew by `angle` in `duration` with a wheel torque of constant magnitude,
    reversed at mid-time. Returns (torque magnitude on the body 4 I angle / T^2, peak body rate 2 angle / T,
    peak wheel momentum I_axis x peak body rate)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the heavy top

@dataclass(frozen=True)
class Top:
    """A symmetric top on a fixed frictionless pivot: a disk of mass m and radius r whose centre is a
    distance l from the pivot along the symmetry axis (a light rod carries it)."""
    m: float = 0.10       # kg
    r: float = 0.030      # m
    l: float = 0.040      # m, pivot to centre of mass
    thickness: float = 0.004

    @property
    def I_pivot(self) -> np.ndarray:
        return parallel_axis(inertia_solid_cylinder(self.m, self.r, self.thickness), self.m, np.array([0.0, 0.0, self.l]))


def gravity_torque_body(q: np.ndarray, m: float, l: float, g: float = G0) -> np.ndarray:
    """Torque of gravity about the pivot in body axes: (l e3) x (m g_body), with g_body = R(q)^T (0, 0, -g)
    the gravitational acceleration expressed in body axes and the centre of mass at l e3."""
    raise NotImplementedError


def top_initial_state(theta: float, w3: float, Omega: float) -> tuple[np.ndarray, np.ndarray]:
    """Quaternion and body rates for a top tilted by theta from the vertical (rotation about the inertial x
    axis), spinning with body rate w3 about its axis and precessing at Omega about the vertical:
    w_body = Omega z_body + psi_dot e3 = (0, Omega sin(theta), w3)."""
    raise NotImplementedError


def steady_precession_rate(I1: float, I3: float, m: float, l: float, theta: float, w3: float, g: float = G0) -> float:
    """Slow steady-precession rate of a heavy symmetric top: the smaller root of
    I1 cos(theta) Omega^2 - I3 w3 Omega + m g l = 0 (the two cones of steady precession); at theta = 90 deg
    the equation is linear and Omega = m g l / (I3 w3). Raises ValueError if the top spins too slowly for
    steady precession at this tilt (negative discriminant)."""
    raise NotImplementedError


def fast_top_precession(I3: float, m: float, l: float, w3: float, g: float = G0) -> float:
    """Gyroscopic approximation Omega = m g l / (I3 w3): torque = rate of change of a spin angular momentum
    I3 w3 that dominates every other term."""
    raise NotImplementedError


def top_invariants(rot: Rotation, I: np.ndarray, m: float, l: float, g: float = G0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The three conserved quantities of the heavy symmetric top along a simulated motion:
    total energy w.(I w)/2 + m g l (R e3)_z, the vertical component of the inertial angular momentum,
    and the spin component (I w)_3."""
    raise NotImplementedError


def nutation_turning_tilt(I1: float, L3: float, m: float, l: float, theta0: float, g: float = G0) -> float:
    """Lower turning point of the axis of a heavy symmetric top released at tilt theta0 with no precession
    and no nodding (theta_dot = phi_dot = 0) and spin angular momentum L3 = I3 w3. Energy and L_z
    conservation at a turning point give

        cos(theta0) - cos(theta) = 2 I1 m g l sin(theta)^2 / L3^2,

    solved for theta in (theta0, pi) with `brentq`. For a fast top theta - theta0 ~ 2 I1 m g l sin(theta0) / L3^2."""
    raise NotImplementedError


def axis_tilt_and_azimuth(rot: Rotation) -> tuple[np.ndarray, np.ndarray]:
    """Tilt from the vertical (rad) and unwrapped azimuth (rad) of the symmetry axis R e3 along the motion."""
    raise NotImplementedError


def measured_precession_rate(rot: Rotation) -> float:
    """Mean precession rate (rad/s): the slope of a least-squares line through the azimuth of the axis."""
    raise NotImplementedError


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save spinning_top.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)

    # ---- the top
    top = Top()
    I = top.I_pivot
    I1, I3 = I[0, 0], I[2, 2]
    print(f"top: disk {top.m * 1e3:.0f} g, r = {top.r * 1e3:.0f} mm, l = {top.l * 1e3:.0f} mm; I1 (pivot) = {I1:.4e} kg m^2, I3 = {I3:.4e} kg m^2; "
          f"m g l = {top.m * G0 * top.l:.5f} N m")
    w3 = 2.0 * math.pi * 50.0                       # 50 rev/s
    theta = math.radians(30.0)
    Om_fast = fast_top_precession(I3, top.m, top.l, w3)
    Om = steady_precession_rate(I1, I3, top.m, top.l, theta, w3)
    print(f"spin 50 rev/s = {w3:.1f} rad/s, tilt 30 deg: L3 = I3 w3 = {I3 * w3:.4e} kg m^2/s; fast-top precession {Om_fast:.4f} rad/s "
          f"(period {2 * math.pi / Om_fast:.2f} s); steady precession {Om:.4f} rad/s; spin-to-precession ratio {w3 / Om:.0f}")
    q0, wb0 = top_initial_state(theta, w3, Om)
    tau_g = lambda t, q, w: gravity_torque_body(q, top.m, top.l)
    rot = simulate_rigid_body(I, wb0, q0, 1.5, torque=tau_g, n_samples=1501, rtol=1e-10)
    E, Lz, L3 = top_invariants(rot, I, top.m, top.l)
    tilt, az = axis_tilt_and_azimuth(rot)
    print(f"steady start, 1.5 s: precession measured {measured_precession_rate(rot):.5f} rad/s; tilt swing {math.degrees(np.ptp(tilt)):.2e} deg; "
          f"relative drift E {np.ptp(E) / abs(E[0]):.1e}, L_z {np.ptp(Lz) / abs(Lz[0]):.1e}, L_3 {np.ptp(L3) / abs(L3[0]):.1e}")
    q1, wb1 = top_initial_state(theta, w3, 0.0)
    rot1 = simulate_rigid_body(I, wb1, q1, 1.5, torque=tau_g, n_samples=6001, rtol=1e-10)
    tilt1, _ = axis_tilt_and_azimuth(rot1)
    E1, Lz1, L31 = top_invariants(rot1, I, top.m, top.l)
    nut_pred = 2 * I1 * top.m * G0 * top.l * math.sin(theta) / (I3 * w3) ** 2
    th_turn = nutation_turning_tilt(I1, I3 * w3, top.m, top.l, theta)
    print(f"released without precession: mean precession {measured_precession_rate(rot1):.4f} rad/s, tilt swing {math.degrees(np.ptp(tilt1)):.4f} deg "
          f"(turning point from the invariants {math.degrees(th_turn - theta):.4f} deg; fast-top estimate 2 I1 m g l sin(theta)/(I3 w3)^2 = {math.degrees(nut_pred):.3f} deg), "
          f"nutation frequency I3 w3/I1 = {I3 * w3 / I1:.0f} rad/s; "
          f"invariants drift {np.ptp(E1) / abs(E1[0]):.1e}, {np.ptp(Lz1) / abs(Lz1[0]):.1e}, {np.ptp(L31) / abs(L31[0]):.1e}")
    print("precession against spin rate (tilt 30 deg), steady start:")
    spins = [20.0, 35.0, 50.0, 100.0]
    meas = []
    for f in spins:
        ww = 2 * math.pi * f
        Oms = steady_precession_rate(I1, I3, top.m, top.l, theta, ww)
        qq, wwb = top_initial_state(theta, ww, Oms)
        r2 = simulate_rigid_body(I, wwb, qq, 1.0, torque=tau_g, n_samples=501, rtol=1e-9)
        meas.append(measured_precession_rate(r2))
        print(f"  {f:5.0f} rev/s: fast-top {fast_top_precession(I3, top.m, top.l, ww):.4f}, steady root {Oms:.4f}, simulated {meas[-1]:.4f} rad/s")
    w_min = 2 * math.sqrt(I1 * top.m * G0 * top.l * math.cos(theta)) / I3
    print(f"slowest spin for steady precession at 30 deg: (I3 w3)^2 = 4 I1 m g l cos(theta): {w_min:.1f} rad/s = {w_min / (2 * math.pi):.1f} rev/s")

    # ---- torque-free: symmetric nutation and the intermediate axis
    Isat = inertia_box(500.0, 1.0, 1.5, 2.0)          # a 500 kg box satellite
    d = np.diag(Isat)
    print(f"\nbox satellite 500 kg, 1.0 x 1.5 x 2.0 m: I = {d[0]:.1f}, {d[1]:.1f}, {d[2]:.1f} kg m^2")
    w_spin = 0.5
    for ax in (0, 1, 2):
        w0 = np.full(3, 1e-4); w0[ax] = w_spin
        r3 = simulate_rigid_body(Isat, w0, np.array([1.0, 0, 0, 0]), 120.0, n_samples=2401)
        L = angular_momentum_inertial(r3.q, r3.w, Isat)
        T = rotational_energy(Isat, r3.w)
        others = [i for i in range(3) if i != ax]
        dev = np.max(np.abs(r3.w[:, others]))
        print(f"  spin {w_spin} rad/s about axis {ax + 1}: max off-axis rate {dev:.2e} rad/s; |L| drift {np.ptp(np.linalg.norm(L, axis=1)) / np.linalg.norm(L[0]):.1e}, "
              f"L direction drift {np.max(np.linalg.norm(L - L[0], axis=1)) / np.linalg.norm(L[0]):.1e}, energy drift {np.ptp(T) / T[0]:.1e}")
    sig = intermediate_axis_growth_rate(Isat, w_spin)
    r4 = simulate_rigid_body(Isat, np.array([1e-6, w_spin, 1e-6]), np.array([1.0, 0, 0, 0]), 150.0, n_samples=3001)
    pert = np.hypot(r4.w[:, 0], r4.w[:, 2])
    first_big = np.argmax(pert > 1e-2)                 # fit only the linear (small-perturbation) stage
    win = (pert > 1e-5) & (np.arange(pert.size) < first_big)
    sig_meas = np.polyfit(r4.t[win], np.log(pert[win]), 1)[0]
    print(f"  intermediate axis: predicted growth rate {sig:.4f} 1/s, measured {sig_meas:.4f} 1/s; the spin first reverses at t = {r4.t[np.argmax(r4.w[:, 1] < 0)]:.1f} s (ln(0.5/1e-6)/sigma = {math.log(0.5 / 1e-6) / sig:.1f} s)")
    Ir = np.diag([800.0, 800.0, 1200.0])
    r5 = simulate_rigid_body(Ir, np.array([0.01, 0.0, 1.0]), np.array([1.0, 0, 0, 0]), 30.0, n_samples=6001)
    lam = symmetric_nutation_rate(800.0, 1200.0, 1.0)
    ph = np.unwrap(np.arctan2(r5.w[:, 1], r5.w[:, 0]))
    print(f"  symmetric spinner I1 = 800, I3 = 1200 kg m^2 at 1 rad/s: body nutation rate predicted {lam:.4f} rad/s, measured {np.polyfit(r5.t, ph, 1)[0]:.4f} rad/s")

    # ---- reaction wheel slew
    Iz = d[2]
    ang, Tsl = math.radians(30.0), 60.0
    tau, wmax, hmax = bang_bang_slew(Iz, ang, Tsl)
    wheel = lambda t, q, w: np.array([0.0, 0.0, -tau if t < Tsl / 2 else tau])     # the body feels the opposite torque
    r6 = simulate_rigid_body(Isat, np.zeros(3), np.array([1.0, 0, 0, 0]), Tsl, wheel_torque=wheel, n_samples=1201, max_step=0.5)
    yaw = 2 * np.arctan2(r6.q[:, 3], r6.q[:, 0])
    Ltot = angular_momentum_inertial(r6.q, r6.w, Isat, None) + np.array([quat_to_matrix(q) @ h for q, h in zip(r6.q, r6.h)])
    print(f"\nreaction-wheel slew of 30 deg about axis 3 in {Tsl:.0f} s: wheel torque {1e3 * tau:.2f} mN m, peak body rate {math.degrees(wmax):.3f} deg/s, "
          f"peak wheel momentum {hmax:.4f} N m s")
    print(f"  simulated: final angle {math.degrees(yaw[-1]):.4f} deg, final body rate {np.linalg.norm(r6.w[-1]):.1e} rad/s, "
          f"peak |h| {np.max(np.abs(r6.h[:, 2])):.4f} N m s, total angular momentum max {np.max(np.linalg.norm(Ltot, axis=1)):.1e} N m s")
    Iw = 0.02
    print(f"  a wheel with I_w = {Iw} kg m^2 must reach {hmax / Iw:.1f} rad/s = {hmax / Iw * 60 / (2 * math.pi):.0f} rpm relative to the body")

    # ---- figure
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    ax = axes[0]
    e3 = np.array([quat_to_matrix(q)[:, 2] for q in rot1.q])
    ax.plot(e3[:, 0], e3[:, 1], lw=0.5)
    ax.set_aspect("equal"); ax.set_xlabel("x of the axis tip"); ax.set_ylabel("y"); ax.set_title("top released without precession (1.5 s)")
    ax = axes[1]
    ff = np.linspace(18, 120, 50)
    ax.plot(ff, [fast_top_precession(I3, top.m, top.l, 2 * math.pi * x) for x in ff], "k--", label="fast top")
    ax.plot(ff, [steady_precession_rate(I1, I3, top.m, top.l, theta, 2 * math.pi * x) for x in ff], "C0", label="steady root")
    ax.plot(spins, meas, "o", color="C1", label="simulated")
    ax.set_xlabel("spin (rev/s)"); ax.set_ylabel("precession (rad/s)"); ax.set_title("precession against spin"); ax.legend()
    ax = axes[2]
    ax.plot(r4.t, r4.w[:, 0], label="w1"); ax.plot(r4.t, r4.w[:, 1], label="w2 (spin)"); ax.plot(r4.t, r4.w[:, 2], label="w3")
    ax.set_xlabel("t (s)"); ax.set_ylabel("body rate (rad/s)"); ax.set_title("spin about the intermediate axis"); ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "spinning_top.png", dpi=120)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'spinning_top.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
