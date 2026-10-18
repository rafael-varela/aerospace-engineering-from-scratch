"""P06.10 Gyroscopes & Spin Stability: reference solution.

A spin-stabilised satellite, from the ideal rigid spinner to the flat spin that ended Explorer 1's
planned attitude:

* torque-free nutation of an axisymmetric body: the body-frame nutation rate (I3 - It) w3 / It, the
  inertial coning rate H / It of the symmetry axis about the fixed angular momentum, and the constant
  half-cone angle with tan(theta) = It wt / (I3 w3);
* gyroscopic precession under a torque: an inertially fixed torque tilts H towards itself at the rate
  tau / H, and a torque kept perpendicular to H and to a fixed axis turns H about that axis at
  Omega = tau / (H sin(theta)) (the gyroscope);
* spin averaging: a torque fixed in a fast-spinning body (a thrust misalignment) only makes H wander on
  a small circle of angular diameter 2 tau / (I3 w3^2), so the pointing error falls as 1/spin^2;
* the major/minor-axis rule with energy dissipation: at fixed |H|, kinetic energy is lowest for spin about
  the axis of largest moment; the energy-sink argument gives d(theta)/dt from the dissipated power;
* a Kane damper (a sphere in a viscous fluid inside the body) as a concrete energy sink: the linearised
  nutation growth or decay rate in closed form (the eigenvalues of a 2 x 2 complex matrix), checked
  against the nonlinear simulation, which conserves the total angular momentum, loses energy
  monotonically, and ends in a flat spin with T = H^2 / (2 (It + J)).

Conventions as in P06.08/P06.09: body principal axes with axis 3 the symmetry (spin) axis, It = I1 = I2;
the quaternion q (scalar first) gives the passive C_b/i; H_i = C^T (I w); dq/dt = q * (0, w) / 2.

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes spin_stability.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-06/lesson-10-gyroscopes-spin-stability
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

TorqueFn = Callable[[float, np.ndarray, np.ndarray], np.ndarray]


# ---------------------------------------------------------------- helpers (given; from P06.08 and P06.09)

def quat_multiply(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Hamilton product p * q, scalar first."""
    p0, p1, p2, p3 = p
    q0, q1, q2, q3 = q
    return np.array([p0 * q0 - p1 * q1 - p2 * q2 - p3 * q3,
                     p0 * q1 + p1 * q0 + p2 * q3 - p3 * q2,
                     p0 * q2 - p1 * q3 + p2 * q0 + p3 * q1,
                     p0 * q3 + p1 * q2 - p2 * q1 + p3 * q0])


def dcm_from_quat(q: np.ndarray) -> np.ndarray:
    """Passive C_b/i of a quaternion (normalised first)."""
    q = np.asarray(q, dtype=float) / np.linalg.norm(q)
    q0, v = q[0], q[1:]
    vx = np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])
    return (q0 * q0 - v @ v) * np.eye(3) + 2.0 * np.outer(v, v) - 2.0 * q0 * vx


def cross3(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Cross product written out."""
    return np.array([a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]])


@dataclass
class Spin:
    """A simulation result: t (n,), q (n, 4), w (n, 3) body rates; wd (n, 3) damper rates if any."""
    t: np.ndarray
    q: np.ndarray
    w: np.ndarray
    wd: np.ndarray | None = None


def simulate_rigid(I: np.ndarray, w0: np.ndarray, t_end: float, torque: TorqueFn | None = None,
                   q0: np.ndarray | None = None, n_samples: int = 2001, rtol: float = 1e-11,
                   atol: float = 1e-12) -> Spin:
    """Euler's equations with quaternion attitude (the P06.09 simulator); torque(t, q, w) in body axes."""
    I = np.asarray(I, dtype=float)
    Iinv = np.linalg.inv(I)
    q0 = np.array([1.0, 0.0, 0.0, 0.0]) if q0 is None else np.asarray(q0, dtype=float)

    def rhs(t, y):
        q, w = y[:4], y[4:]
        tau = np.zeros(3) if torque is None else torque(t, q, w)
        return np.concatenate([0.5 * quat_multiply(q, np.array([0.0, *w])), Iinv @ (tau - cross3(w, I @ w))])

    sol = solve_ivp(rhs, (0.0, t_end), np.concatenate([q0, w0]), method="DOP853",
                    t_eval=np.linspace(0.0, t_end, n_samples), rtol=rtol, atol=atol)
    q = sol.y[:4].T
    return Spin(t=sol.t, q=q / np.linalg.norm(q, axis=1)[:, None], w=sol.y[4:].T)


def inertial_momentum(I: np.ndarray, res: Spin, J: float = 0.0) -> np.ndarray:
    """H_i(t) = C^T (I w + J wd) for every sample (J wd is the damper's share, if any)."""
    Hb = res.w @ np.asarray(I).T
    if res.wd is not None:
        Hb = Hb + J * res.wd
    return np.array([dcm_from_quat(q).T @ h for q, h in zip(res.q, Hb)])


def spin_axis_inertial(res: Spin) -> np.ndarray:
    """The body symmetry axis 3 in inertial components (row 3 of C_b/i) for every sample."""
    return np.array([dcm_from_quat(q)[2] for q in res.q])


def fit_growth_rate(t: np.ndarray, amp: np.ndarray, lo: float, hi: float) -> float:
    """Least-squares slope of ln(amp) against t over the samples with lo < amp < hi (an exponential's rate)."""
    m = (amp > lo) & (amp < hi)
    return float(np.polyfit(t[m], np.log(amp[m]), 1)[0])


# ---------------------------------------------------------------- the axisymmetric spinner

def nutation(It: float, I3: float, w3: float, wt: float) -> dict:
    """Torque-free nutation of an axisymmetric body spinning at w3 with transverse rate wt (magnitude).

    H = sqrt((I3 w3)^2 + (It wt)^2); body-frame nutation rate lam = (I3 - It) w3 / It (the transverse
    rate turns in the body at lam); inertial coning rate of the symmetry axis about H: H / It;
    half-cone angle theta between the axis and H: tan(theta) = It wt / (I3 w3).
    Returns {"H", "body_rate", "inertial_rate", "theta"}."""
    H = math.hypot(I3 * w3, It * wt)
    return {"H": H, "body_rate": (I3 - It) * w3 / It, "inertial_rate": H / It,
            "theta": math.atan2(It * wt, I3 * w3)}


def axis_coning_rate(res: Spin, H_dir: np.ndarray) -> float:
    """Measured rate (rad/s) at which the symmetry axis turns about the unit vector H_dir: the slope of the
    unwrapped azimuth of the axis in the plane perpendicular to H_dir."""
    h = np.asarray(H_dir, dtype=float) / np.linalg.norm(H_dir)
    a = np.cross(h, [1.0, 0.0, 0.0])
    if np.linalg.norm(a) < 1e-6:
        a = np.cross(h, [0.0, 1.0, 0.0])
    a /= np.linalg.norm(a)
    b = np.cross(h, a)
    e3 = spin_axis_inertial(res)
    az = np.unwrap(np.arctan2(e3 @ b, e3 @ a))
    return float(np.polyfit(res.t, az, 1)[0])


# ---------------------------------------------------------------- precession under torque

def precession_rate(tau: float, H: float, theta: float) -> float:
    """Gyroscope: a torque of size tau kept perpendicular to H and to a fixed axis z turns H about z at
    Omega = tau / (H sin(theta)), theta the angle between H and z (OpenStax 11.4's tau / L for theta = 90 deg)."""
    return tau / (H * math.sin(theta))


def gyroscope_torque(I: np.ndarray, tau0: float, z: np.ndarray | None = None) -> TorqueFn:
    """Body-axis torque function for the gyroscope: size tau0, along z x H_i in inertial axes (horizontal
    and perpendicular to H, like gravity's torque on a top), converted to body axes."""
    zz = np.array([0.0, 0.0, 1.0]) if z is None else np.asarray(z, dtype=float) / np.linalg.norm(z)
    I = np.asarray(I, dtype=float)

    def torque(t, q, w):
        C = dcm_from_quat(q)
        Hi = C.T @ (I @ w)
        d = np.cross(zz, Hi)
        return C @ (tau0 * d / np.linalg.norm(d))
    return torque


def tilt_under_fixed_torque(tau: float, H: float, t: np.ndarray) -> np.ndarray:
    """An inertially fixed torque perpendicular to H0: H(t) = H0 + tau t, so the tilt of H from H0 is
    atan(tau t / H). The axis moves towards the torque, at right angles to the force that makes it."""
    return np.arctan(tau * np.asarray(t, dtype=float) / H)


def spin_averaged_wander(tau: float, I3: float, w3: float) -> float:
    """Largest angle (rad) between H and its initial direction under a constant torque tau fixed in a body
    spinning fast at w3 about its symmetry axis: the torque turns with the body, H's tip runs round a circle
    of radius tau / w3 (N m s) through its start, so the largest tilt is about 2 tau / (I3 w3^2)."""
    return 2.0 * tau / (I3 * w3 * w3)


# ---------------------------------------------------------------- energy, dissipation and the major-axis rule

def energy_at_fixed_momentum(H: float, It: float, I3: float, theta: np.ndarray | float) -> np.ndarray | float:
    """Kinetic energy of an axisymmetric body with angular momentum H whose symmetry axis is at angle theta
    from H: T = (H^2 / 2)(sin^2(theta)/It + cos^2(theta)/I3)."""
    th = np.asarray(theta, dtype=float)
    return 0.5 * H * H * (np.sin(th) ** 2 / It + np.cos(th) ** 2 / I3)


def energy_sink_theta_rate(P: float, H: float, It: float, I3: float, theta: float) -> float:
    """Energy-sink estimate: if internal dissipation removes power P > 0 at fixed H, dT/dt = -P with T from
    `energy_at_fixed_momentum` gives d(theta)/dt = -P / (H^2 sin(theta) cos(theta) (1/It - 1/I3)).
    Positive (nutation grows) for a prolate spinner I3 < It; negative (decays) for an oblate one."""
    return -P / (H * H * math.sin(theta) * math.cos(theta) * (1.0 / It - 1.0 / I3))


def spin_stability(I_principal: np.ndarray, axis: int, dissipative: bool) -> str:
    """'stable' or 'unstable' for spin about principal axis `axis` (0, 1, 2): a rigid body is stable about
    the axes of largest and smallest moment; with any internal energy dissipation only the largest is."""
    d = np.diag(I_principal) if np.ndim(I_principal) == 2 else np.asarray(I_principal, dtype=float)
    if d[axis] == d.max() and np.sum(d == d.max()) == 1:
        return "stable"
    if not dissipative and d[axis] == d.min() and np.sum(d == d.min()) == 1:
        return "stable"
    return "unstable"


def flat_spin(H: float, I_max: float) -> dict:
    """The end state of dissipation at fixed H: spin about the axis of largest moment I_max at H/I_max,
    with kinetic energy H^2 / (2 I_max)."""
    return {"rate": H / I_max, "T": H * H / (2.0 * I_max)}


# ---------------------------------------------------------------- the Kane damper

def kane_linear_matrix(It: float, I3: float, J: float, c: float, Omega: float) -> np.ndarray:
    """Linearised transverse dynamics of a symmetric body (It, I3) spinning at Omega with a Kane damper
    (a sphere of moment of inertia J coupled by a viscous torque c (wd - w)), in complex form
    x = (w1 + i w2, wd1 + i wd2):  dx/dt = A x with
        It dw/dt  = i (I3 - It) Omega w + c (wd - w),
        J dwd/dt  = -c (wd - w) - i J Omega (wd - w)."""
    return np.array([[1j * (I3 - It) * Omega / It - c / It, c / It],
                     [c / J + 1j * Omega, -c / J - 1j * Omega]], dtype=complex)


def kane_nutation_mode(It: float, I3: float, J: float, c: float, Omega: float) -> complex:
    """The eigenvalue of `kane_linear_matrix` that continues the undamped body nutation mode
    i (I3 - It) Omega / It; its real part is the nutation growth (> 0) or decay (< 0) rate."""
    ev = np.linalg.eigvals(kane_linear_matrix(It, I3, J, c, Omega))
    target = 1j * (I3 - It) * Omega / It
    return complex(ev[np.argmin(np.abs(ev - target))])


def simulate_with_damper(I: np.ndarray, J: float, c: float, w0: np.ndarray, wd0: np.ndarray, t_end: float,
                         n_samples: int = 2001, rtol: float = 1e-10, atol: float = 1e-12) -> Spin:
    """Body with a Kane damper, states [q, w, wd] (wd the sphere's absolute rate in body axes):
        I dw/dt  = -w x (I w) + c (wd - w),
        J dwd/dt = -c (wd - w) - J w x wd.
    The total H = I w + J wd is conserved (in inertial axes); T = w.(I w)/2 + J |wd|^2/2 falls at
    c |wd - w|^2."""
    I = np.asarray(I, dtype=float)
    Iinv = np.linalg.inv(I)
    q0 = np.array([1.0, 0.0, 0.0, 0.0])

    def rhs(t, y):
        q, w, wd = y[:4], y[4:7], y[7:]
        rel = wd - w
        return np.concatenate([0.5 * quat_multiply(q, np.array([0.0, *w])),
                               Iinv @ (c * rel - cross3(w, I @ w)),
                               (-c * rel) / J - cross3(w, wd)])

    sol = solve_ivp(rhs, (0.0, t_end), np.concatenate([q0, w0, wd0]), method="DOP853",
                    t_eval=np.linspace(0.0, t_end, n_samples), rtol=rtol, atol=atol)
    q = sol.y[:4].T
    return Spin(t=sol.t, q=q / np.linalg.norm(q, axis=1)[:, None], w=sol.y[4:7].T, wd=sol.y[7:].T)


def damper_energy(I: np.ndarray, J: float, res: Spin) -> np.ndarray:
    """Total kinetic energy of body and damper at every sample."""
    return 0.5 * np.sum(res.w * (res.w @ np.asarray(I).T), axis=1) + 0.5 * J * np.sum(res.wd ** 2, axis=1)


# ---------------------------------------------------------------- the lesson's output

def explorer_like(m: float = 13.97, length: float = 2.03, diameter: float = 0.152) -> np.ndarray:
    """Principal inertia of a uniform solid cylinder with Explorer 1's approximate overall size and mass
    (illustrative: the real satellite was not uniform): I3 = m r^2/2, It = m (3 r^2 + L^2)/12."""
    r = 0.5 * diameter
    It = m * (3 * r * r + length * length) / 12.0
    return np.diag([It, It, 0.5 * m * r * r])


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save spin_stability.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    # an illustrative oblate spin-stabilised satellite
    I = np.diag([800.0, 800.0, 1200.0])
    w3, wt = 1.0, 0.05
    nu = nutation(800.0, 1200.0, w3, wt)
    w0 = np.array([wt, 0.0, w3])
    res = simulate_rigid(I, w0, 60.0, n_samples=3001)
    Hi = inertial_momentum(I, res)
    e3 = spin_axis_inertial(res)
    Hh = Hi[0] / np.linalg.norm(Hi[0])
    cone = np.degrees(np.arccos(np.clip(e3 @ Hh, -1, 1)))
    print(f"oblate satellite It = 800, I3 = 1200 kg m^2, spin 1 rad/s, transverse 0.05 rad/s: H = {nu['H']:.4f} N m s")
    print(f"  body nutation rate (I3 - It) w3 / It = {nu['body_rate']:.4f} rad/s; inertial coning rate H / It = {nu['inertial_rate']:.6f} rad/s; "
          f"half-cone angle {math.degrees(nu['theta']):.4f} deg")
    print(f"  simulated: coning rate {axis_coning_rate(res, Hh):.6f} rad/s, half-cone angle {cone.min():.4f} to {cone.max():.4f} deg, "
          f"H drift {np.abs(Hi - Hi[0]).max() / np.linalg.norm(Hi[0]):.1e}")
    # gyroscopic precession
    th = math.radians(30.0)
    H0 = 1200.0
    # start with H at 30 deg from inertial z: tilt the body about inertial x by 30 deg (pure spin)
    q_tilt = np.array([math.cos(th / 2), math.sin(th / 2), 0.0, 0.0])
    tau0 = 2.0
    resg = simulate_rigid(I, np.array([0.0, 0.0, 1.0]), 60.0, torque=gyroscope_torque(I, tau0), q0=q_tilt, n_samples=1201)
    Hg = inertial_momentum(I, resg)
    az = np.unwrap(np.arctan2(Hg[:, 1], Hg[:, 0]))
    tilt = np.degrees(np.arccos(Hg[:, 2] / np.linalg.norm(Hg, axis=1)))
    print(f"gyroscope: torque {tau0} N m kept perpendicular to H and to z, H = {H0} N m s at 30 deg: predicted precession "
          f"tau/(H sin theta) = {precession_rate(tau0, H0, th):.6f} rad/s; simulated {np.polyfit(resg.t, az, 1)[0]:.6f} rad/s; "
          f"tilt {tilt.min():.6f} to {tilt.max():.6f} deg; |H| change {np.ptp(np.linalg.norm(Hg, axis=1)):.1e} N m s")
    # spin averaging of a body-fixed torque (thrust misalignment)
    print("spin averaging: 1 N m torque fixed along body axis 1, H wander (deg) after 60 s:")
    for spin in [1.0, 2.0, 4.0]:
        rs = simulate_rigid(I, np.array([0.0, 0.0, spin]), 60.0, torque=lambda t, q, w: np.array([1.0, 0.0, 0.0]), n_samples=6001)
        Hs = inertial_momentum(I, rs)
        wander = np.degrees(np.arccos(np.clip(Hs @ Hs[0] / (np.linalg.norm(Hs, axis=1) * np.linalg.norm(Hs[0])), -1, 1)))
        print(f"  spin {spin:.1f} rad/s: simulated max {wander.max():.5f}, closed form 2 tau/(I3 w3^2) {math.degrees(spin_averaged_wander(1.0, 1200.0, spin)):.5f}; "
              f"(a non-spinning body would turn {math.degrees(0.5 * (1.0 / 800.0) * 60.0 ** 2):.0f} deg)")

    # Explorer 1, illustrative
    Ie = explorer_like()
    It, I3 = Ie[0, 0], Ie[2, 2]
    spin = 750 * 2 * math.pi / 60
    He = I3 * spin
    fs = flat_spin(He, It)
    print(f"Explorer-1-like cylinder (13.97 kg, 2.03 m x 0.152 m, uniform, illustrative): It = {It:.4f}, I3 = {I3:.5f} kg m^2; "
          f"spin 750 rpm = {spin:.3f} rad/s; H = {He:.4f} N m s; T0 = {0.5 * I3 * spin ** 2:.2f} J")
    print(f"  rigid model: minor-axis spin is {spin_stability(Ie, 2, False)}; with dissipation it is {spin_stability(Ie, 2, True)}; "
          f"flat spin {fs['rate']:.4f} rad/s ({fs['rate'] * 60 / (2 * math.pi):.2f} rpm), T = {fs['T']:.4f} J "
          f"({fs['T'] / (0.5 * I3 * spin ** 2) * 100:.2f} % of T0)")
    for J, c in [(4e-3, 1e-4), (4e-3, 3e-3), (2e-2, 1e-2)]:
        lam = kane_nutation_mode(It, I3, J, c, spin)
        print(f"  Kane damper J = {J}, c = {c}: nutation mode {lam.real:.3e} {lam.imag:+.3f} i 1/s; e-folding time {1 / lam.real:.1f} s "
              f"= {spin / (2 * math.pi) / lam.real:.0f} revolutions")

    # a prolate spinner with a strong damper, simulated to the flat spin (fast to run)
    It, I3, spin, J, c = 1500.0, 600.0, 3.0, 60.0, 100.0
    Ip = np.diag([It, It, I3])
    lam = kane_nutation_mode(It, I3, J, c, spin)
    rese = simulate_with_damper(Ip, J, c, np.array([1e-3, 0.0, spin]), np.array([0.0, 0.0, spin]), 600.0,
                                n_samples=6001, rtol=1e-10, atol=1e-12)
    Ht = inertial_momentum(Ip, rese, J)
    He = float(np.linalg.norm(Ht[0]))
    Et = damper_energy(Ip, J, rese)
    amp = np.hypot(rese.w[:, 0], rese.w[:, 1])
    sig = fit_growth_rate(rese.t, amp, 1e-2, 0.2)
    e3e = spin_axis_inertial(rese)
    nut = np.degrees(np.arccos(np.clip(e3e @ (Ht[0] / He), -1, 1)))
    Tf = He ** 2 / (2 * (It + J))
    print(f"prolate spinner It = {It}, I3 = {I3} kg m^2 at {spin} rad/s with a Kane damper J = {J} kg m^2, c = {c} N m s: "
          f"linear nutation mode {lam.real:.5f} {lam.imag:+.4f} i 1/s")
    print(f"  simulated 600 s: fitted growth rate {sig:.5f} 1/s; |H| and H direction drift {np.abs(Ht - Ht[0]).max() / He:.1e}; "
          f"energy never rises: {bool(np.all(np.diff(Et) <= 1e-9 * Et[0]))}; axis passes 89 deg from H at {rese.t[np.argmax(nut > 89)]:.1f} s")
    print(f"  final: axis {nut[-1]:.4f} deg from H; body rate {np.linalg.norm(rese.w[-1]):.6f} rad/s against H/(It + J) = {He / (It + J):.6f}; "
          f"T = {Et[-1]:.4f} J against H^2/(2 (It + J)) = {Tf:.4f} J (T0 = {Et[0]:.1f} J)")
    print("  energy-sink check (10 s averages; the co-rotating damper counted in the inertias It + J, I3 + J):")
    Pt = c * np.sum((rese.wd - rese.w) ** 2, axis=1)
    thr = np.radians(nut)
    for deg in (5.0, 20.0, 45.0, 70.0):
        k = int(np.argmax(nut > deg))
        win = (rese.t > rese.t[k] - 5.0) & (rese.t < rese.t[k] + 5.0)
        Pm = float(Pt[win].mean())
        pred = energy_sink_theta_rate(Pm, He, It + J, I3 + J, float(thr[k]))
        simr = float(np.polyfit(rese.t[win], thr[win], 1)[0])
        print(f"    nutation {deg:4.1f} deg at t = {rese.t[k]:6.1f} s: dissipated {Pm:7.3f} W; d(theta)/dt predicted {pred:.5f}, "
              f"simulated {simr:.5f} rad/s ({(pred / simr - 1) * 100:+.1f} %)")
    lo = kane_nutation_mode(800.0, 1200.0, 60.0, 100.0, 1.0)
    print(f"oblate satellite (It 800, I3 1200) with the same damper at 1 rad/s: nutation mode {lo.real:.5f} {lo.imag:+.4f} i 1/s (decays)")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    ax[0].plot(res.t, cone)
    ax[0].set_xlabel("t (s)"); ax[0].set_ylabel("deg"); ax[0].set_title("oblate spinner: axis-to-H angle")
    ax[1].semilogy(rese.t, amp, label="transverse rate |wt|")
    ax[1].semilogy(rese.t, 1e-3 * np.exp(lam.real * rese.t), "--", label="linear theory")
    ax[1].set_ylim(1e-4, 10); ax[1].set_xlabel("t (s)"); ax[1].legend(); ax[1].set_title("prolate spinner with a Kane damper")
    ax[2].plot(rese.t, nut, label="axis-to-H angle (deg)")
    ax2 = ax[2].twinx()
    ax2.plot(rese.t, Et / Et[0], "C1", label="T / T0")
    ax[2].set_xlabel("t (s)"); ax[2].set_title("to the flat spin")
    ax[2].legend(loc="center right"); ax2.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(out_dir / "spin_stability.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out_dir / 'spin_stability.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
