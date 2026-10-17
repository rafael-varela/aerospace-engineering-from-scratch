"""P06.09 3D Rigid Body Dynamics: Euler's Equations: starter code.

A tumbling-body simulator and the closed forms that check it:

* the inertia tensor's principal axes and moments (an eigenproblem), Euler's equations
  I dw/dt + w x (I w) = tau in any body axes (principal or not), with the attitude carried by a unit
  quaternion and integrated with `scipy.integrate.solve_ivp`;
* the conserved quantities of torque-free motion: the kinetic energy T, the magnitude |H| of the
  angular momentum, and the inertial angular momentum vector itself; the invariable plane w . H = 2T;
* the polhode, the curve the body rate traces on the energy ellipsoid, and its classification by the
  sign of H^2 - 2 T I2 (around the major axis, around the minor axis, or on the separatrix through the
  intermediate axis);
* the exact torque-free solution of an asymmetric body in Jacobi elliptic functions, with its period
  4 K(m)/lambda and the interval 2 K(m)/lambda between the flips of a body spun near its intermediate
  axis (the Dzhanibekov effect), checked against the simulation.

Conventions: body axes; I is the inertia tensor about the centre of mass in body axes (kg m^2); w the body
rate in body axes (rad/s); the attitude quaternion q = [q0, q1, q2, q3] (scalar first) is that of the body
frame relative to the inertial frame in the passive convention of P06.08, so C_b/i = dcm(q) maps inertial
components to body components, the inertial angular momentum is H_i = C^T (I w), and dq/dt = q * (0, w)/2.

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes tumbling_body.png)
Test: pytest curriculum/phase-06/lesson-09-euler-equations-rigid-body
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipj, ellipk, ellipkinc, ellipkm1

TorqueFn = Callable[[float, np.ndarray, np.ndarray], np.ndarray]   # (t, q, w) -> body torque


# ---------------------------------------------------------------- quaternion helpers (given, from P06.08)

def quat_multiply(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Hamilton product p * q, scalar first."""
    p0, p1, p2, p3 = p
    q0, q1, q2, q3 = q
    return np.array([p0 * q0 - p1 * q1 - p2 * q2 - p3 * q3,
                     p0 * q1 + p1 * q0 + p2 * q3 - p3 * q2,
                     p0 * q2 - p1 * q3 + p2 * q0 + p3 * q1,
                     p0 * q3 + p1 * q2 - p2 * q1 + p3 * q0])


def dcm_from_quat(q: np.ndarray) -> np.ndarray:
    """Passive C_b/i of a quaternion (normalised first): (q0^2 - v.v) E + 2 v v^T - 2 q0 [v x]."""
    q = np.asarray(q, dtype=float) / np.linalg.norm(q)
    q0, v = q[0], q[1:]
    vx = np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])
    return (q0 * q0 - v @ v) * np.eye(3) + 2.0 * np.outer(v, v) - 2.0 * q0 * vx


def cross3(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Cross product written out (faster than np.cross for single vectors inside an ODE right-hand side)."""
    return np.array([a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]])


@dataclass
class Tumble:
    """Result of `simulate`: t (n,), w (n, 3) body rates, q (n, 4) attitude quaternions (renormalised)."""
    t: np.ndarray
    w: np.ndarray
    q: np.ndarray


# ---------------------------------------------------------------- inertia

def principal_axes(I: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Principal moments (ascending) and the matrix R whose columns are the principal axes in the given
    body axes, a proper rotation (det R = +1): I = R diag(moments) R^T. The moments are the eigenvalues of
    the symmetric tensor; the axes its orthonormal eigenvectors."""
    raise NotImplementedError


def kinetic_energy(I: np.ndarray, w: np.ndarray) -> np.ndarray | float:
    """T = w . (I w) / 2 for one rate (3,) or a stack (n, 3)."""
    raise NotImplementedError


def angular_momentum_body(I: np.ndarray, w: np.ndarray) -> np.ndarray:
    """H_b = I w for one rate or a stack."""
    raise NotImplementedError


def angular_momentum_inertial(I: np.ndarray, w: np.ndarray, q: np.ndarray) -> np.ndarray:
    """H_i = C_b/i^T (I w) for one state or stacks (n, 3) and (n, 4)."""
    raise NotImplementedError


# ---------------------------------------------------------------- Euler's equations and the simulator

def euler_rhs(I: np.ndarray, w: np.ndarray, torque: np.ndarray | None = None) -> np.ndarray:
    """dw/dt from Euler's equations I dw/dt = tau - w x (I w), for any (symmetric, positive definite)
    inertia tensor in body axes. In principal axes: I1 dw1/dt = (I2 - I3) w2 w3 + tau1 and cyclic."""
    raise NotImplementedError


def simulate(I: np.ndarray, w0: np.ndarray, t_end: float, q0: np.ndarray | None = None,
             torque: TorqueFn | None = None, n_samples: int = 2001, rtol: float = 1e-11,
             atol: float = 1e-12, max_step: float = np.inf) -> Tumble:
    """Integrate Euler's equations and the quaternion kinematics together with solve_ivp (DOP853).
    States [q (4), w (3)]; torque(t, q, w) is in body axes (default zero). The output quaternion is
    renormalised."""
    raise NotImplementedError


def invariant_drift(I: np.ndarray, res: Tumble) -> dict:
    """Relative drifts over a run: max |T - T0|/T0, max ||H| - |H0||/|H0|, and the largest change of any
    inertial component of H relative to |H0|."""
    raise NotImplementedError


def invariable_plane_residual(I: np.ndarray, res: Tumble) -> float:
    """Poinsot: the inertial rate vector stays in the invariable plane w_i . H_i = 2T. Returns the largest
    relative deviation of w_i . H_i / |H| from 2T/|H| over the run."""
    raise NotImplementedError


# ---------------------------------------------------------------- polhodes and the elliptic solution

def polhode_type(I_principal: np.ndarray, w: np.ndarray, rel_tol: float = 1e-12) -> str:
    """Classify torque-free motion of a body with distinct principal moments (in principal axes):
    'major' if H^2 > 2 T I_mid (the rate circles the axis of largest moment), 'minor' if H^2 < 2 T I_mid,
    'separatrix' if they are equal. H^2 - 2 T I_mid = sum of I_k (I_k - I_mid) w_k^2 is summed term by term
    (no cancellation of large numbers), and 'equal' means within rel_tol of the sum of the terms' sizes."""
    raise NotImplementedError


def _sorted_frame(d: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Signed permutation Q (a proper rotation) with J = Q diag(d) Q^T sorted ascending: u = Q w."""
    order = np.argsort(d)
    Q = np.eye(3)[order]
    if np.linalg.det(Q) < 0:
        Q = -Q
    return Q, d[order]


def elliptic_parameters(I_principal: np.ndarray, w0: np.ndarray) -> dict:
    """Parameters of the exact torque-free solution for distinct principal moments.

    In sorted axes J1 < J2 < J3 with T and H^2 = |H|^2 from w0, and for H^2 > 2 T J2 (motion around the
    major axis 3):
        u1 = A1 cn(tau), u2 = A2 sn(tau), u3 = A3 dn(tau),   tau = lambda t + tau0,
        A1^2 = (2 T J3 - H^2) / (J1 (J3 - J1)),  A2^2 = (2 T J3 - H^2) / (J2 (J3 - J2)),
        A3^2 = (H^2 - 2 T J1) / (J3 (J3 - J1)),  lambda^2 = (J3 - J2)(H^2 - 2 T J1) / (J1 J2 J3),
        m = k^2 = (J2 - J1)(2 T J3 - H^2) / ((J3 - J2)(H^2 - 2 T J1)).
    For H^2 < 2 T J2 the roles of axes 1 and 3 swap (dn on axis 1, cn on axis 3). Signs and the direction
    of tau are fixed by w0 and by Euler's equations. Returns a dict with "period" = 4 K(m)/lambda (the
    period of the body-rate motion), "flip_interval" = 2 K(m)/lambda (the time between sign reversals of the
    intermediate-axis rate), "m", "lam", "type" and what `elliptic_solution` needs. Raises ValueError on
    the separatrix or for a symmetric body."""
    raise NotImplementedError


def elliptic_solution(I_principal: np.ndarray, w0: np.ndarray, t: np.ndarray) -> np.ndarray:
    """The exact torque-free body rates w(t) (n, 3) in the given principal axes, from
    `elliptic_parameters`: Jacobi sn, cn, dn of tau = tau0 + sigma lambda t, mapped back from sorted axes."""
    raise NotImplementedError


def small_wobble_frequency(I_principal: np.ndarray, axis: int, spin: float) -> float:
    """Linearised wobble frequency of spin about an extreme principal axis (P04.04):
    sqrt((Ia - Ij)(Ia - Ik) / (Ij Ik)) |spin|. Raises ValueError for the intermediate axis."""
    raise NotImplementedError


def intermediate_growth_rate(I_principal: np.ndarray, spin: float) -> float:
    """sigma = |spin| sqrt((I3 - I2)(I2 - I1) / (I1 I3)) for spin about the intermediate axis (P04.04)."""
    raise NotImplementedError


def flip_times(t: np.ndarray, w_mid: np.ndarray) -> np.ndarray:
    """Times at which the intermediate-axis rate changes sign (linear interpolation between samples)."""
    s = np.sign(w_mid)
    idx = np.nonzero(s[:-1] * s[1:] < 0)[0]
    return t[idx] - w_mid[idx] * (t[idx + 1] - t[idx]) / (w_mid[idx + 1] - w_mid[idx])


# ---------------------------------------------------------------- the lesson's output

T_HANDLE = np.diag([2.2e-6, 6.0e-6, 7.4e-6])     # illustrative T-handle (wing nut), kg m^2, principal axes
BOX_SAT = np.diag([135.4, 208.3, 260.4])         # the P04.04 box satellite, reordered I1 < I2 < I3


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save tumbling_body.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    # a tensor with products of inertia: the box satellite with an offset instrument boom
    R_true = dcm_from_quat(np.array([0.9, 0.1, -0.3, 0.2])).T
    I_full = R_true @ BOX_SAT @ R_true.T
    mom, R = principal_axes(I_full)
    print("box satellite in misaligned axes:")
    print(np.array2string(I_full, precision=3))
    print(f"  principal moments {np.round(mom, 3)} kg m^2; R^T I R diagonal to {np.abs(R.T @ I_full @ R - np.diag(mom)).max():.1e}")

    # Dzhanibekov: the T-handle spun about its intermediate axis with a small disturbance
    I = T_HANDLE
    spin = 2 * math.pi * 4.0                     # 4 rev/s about axis 2
    w0 = np.array([0.01 * spin, spin, 0.0])
    p = elliptic_parameters(I, w0)
    sig = intermediate_growth_rate(I, spin)
    print(f"T-handle I = {np.diag(I)} kg m^2, spin {spin:.3f} rad/s about axis 2, 1 % disturbance on axis 1:")
    print(f"  polhode type '{p['type']}', m = {p['m']:.8f}, lambda = {p['lam']:.4f} rad/s, K(m) = {ellipk(p['m']):.4f}; "
          f"period {p['period']:.4f} s, flip interval {p['flip_interval']:.4f} s; growth rate sigma = {sig:.4f} 1/s")
    res = simulate(I, w0, 3 * p["period"], n_samples=6001)
    dr = invariant_drift(I, res)
    fl = flip_times(res.t, res.w[:, 1])
    we = elliptic_solution(I, w0, res.t)
    print(f"  simulation {res.t[-1]:.2f} s: drifts T {dr['T']:.1e}, |H| {dr['H']:.1e}, inertial H {dr['H_inertial']:.1e}; "
          f"invariable plane {invariable_plane_residual(I, res):.1e}")
    print(f"  flips at {np.round(fl, 4)} s; mean interval {np.mean(np.diff(fl)):.5f} s against 2K/lambda = {p['flip_interval']:.5f} s")
    print(f"  max |w_sim - w_elliptic| / |w0| = {np.abs(res.w - we).max() / np.linalg.norm(w0):.1e}")
    print("  flip interval against the disturbance (closed form 2K/lambda; log estimate (2/lambda) ln(4/k'), k'^2 = 1 - m):")
    rows = []
    for eps in [1e-1, 1e-2, 1e-3, 1e-4, 1e-6]:
        pe = elliptic_parameters(I, np.array([eps * spin, spin, 0.0]))
        kp = math.sqrt(pe["m1"])
        est = (2 / pe["lam"]) * math.log(4 / kp)
        rows.append((eps, pe["flip_interval"], est))
        print(f"    eps {eps:7.0e}: 1 - m = {pe['m1']:.3e}, 2K/lambda = {pe['flip_interval']:.4f} s, log estimate {est:.4f} s")
    # the minor-axis side of the separatrix
    w0b = np.array([0.0, spin, 0.01 * spin])
    pb = elliptic_parameters(I, w0b)
    resb = simulate(I, w0b, 2 * pb["period"], n_samples=4001)
    eb = np.abs(resb.w - elliptic_solution(I, w0b, resb.t)).max() / np.linalg.norm(w0b)
    print(f"  disturbance on axis 3 instead: polhode type '{pb['type']}', period {pb['period']:.4f} s; "
          f"max |w_sim - w_elliptic| / |w0| = {eb:.1e}")
    # stable spin: period against the small-wobble frequency
    w0s = np.array([0.0, 0.01 * spin, spin])
    ps = elliptic_parameters(I, w0s)
    print(f"  spin about the major axis 3 with 1 % on axis 2: period {ps['period']:.5f} s; "
          f"2 pi / small-wobble frequency {2 * math.pi / small_wobble_frequency(I, 2, spin):.5f} s")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    for k, lab in enumerate(["w1", "w2", "w3"]):
        ax[0].plot(res.t, res.w[:, k], label=lab)
    ax[0].plot(res.t[::120], we[::120, 1], "k.", ms=3, label="w2 elliptic")
    ax[0].set_xlabel("t (s)"); ax[0].set_ylabel("rad/s"); ax[0].legend(); ax[0].set_title("T-handle: body rates")
    # polhodes on the energy ellipsoid seen along axis 1 (w2-w3 plane), for several H at fixed T
    T0 = kinetic_energy(I, w0)
    d = np.diag(I)
    for frac in np.linspace(0.05, 0.95, 9):
        # a family of initial rates on the energy ellipsoid from the axis-2-ish region to the poles
        a = frac * math.pi / 2
        for sgn in (1, -1):
            v = np.array([math.cos(a), sgn * 1e-3, math.sin(a)])
            v = v / math.sqrt(np.sum(d * v * v) / (2 * T0))
            tt = np.linspace(0, elliptic_parameters(I, v)["period"], 400)
            ww = elliptic_solution(I, v, tt)
            ax[1].plot(ww[:, 0], ww[:, 2], lw=0.8, color="C0" if polhode_type(I, v) == "major" else "C1")
    ax[1].set_xlabel("w1 (rad/s)"); ax[1].set_ylabel("w3 (rad/s)"); ax[1].set_title("polhodes (blue: major, orange: minor)")
    Hi = angular_momentum_inertial(I, res.w, res.q)
    ax[2].plot(res.t, (Hi - Hi[0]) / np.linalg.norm(Hi[0]))
    ax[2].set_xlabel("t (s)"); ax[2].set_title("inertial H drift (relative)")
    fig.tight_layout()
    fig.savefig(out_dir / "tumbling_body.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out_dir / 'tumbling_body.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
