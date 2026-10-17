"""P08.02 Fluid Kinematics: reference solution.

Can a smoke photograph of a wind-tunnel model be read as a picture of the streamlines when the
tunnel's gust generator swings the flow direction to and fro? This module describes a 2-D flow
field u(x, y, t) the way kinematics does, without asking what forces drive it:

* the three families of lines: streamlines (tangent to u at one frozen instant), pathlines (the
  track of one fluid particle over time) and streaklines (the locus at one instant of every
  particle that has passed a fixed point: what a smoke wand or a dye port shows);
* the velocity-gradient tensor at a point, split into the strain-rate tensor S (symmetric part:
  stretching and shearing) and the rotation tensor W (antisymmetric part), with the divergence
  (volume change) and the vorticity omega = dv/dx - du/dy (twice the local rotation rate);
* the material derivative Df/Dt = df/dt + u . grad f, the rate of change seen by a particle, and
  in particular the particle acceleration Du/Dt.

The independent checks: closed-form pathlines and streaklines of a uniform oscillating stream
against the RK4 integration; a steady flow in which all three lines coincide; a Taylor-Green
vortex that is unsteady but whose lines still coincide (only the speed changes, not the
direction); the material derivative against the second time derivative of a computed pathline;
and the area of a material loop, which an incompressible flow must keep.

Every velocity field here is a function ``vel(x, y, t) -> (u, v)`` that accepts arrays.
Units: SI (m, s); the dimensionless fields (Taylor-Green) use consistent units.

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes flow_lines.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-08/lesson-02-fluid-kinematics
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable

import numpy as np

Field = Callable[[np.ndarray, np.ndarray, float], tuple[np.ndarray, np.ndarray]]


# ---------------------------------------------------------------- example flow fields (given)

def oscillating_flow(U: float, V: float, omega: float) -> Field:
    """A uniform stream whose direction swings: u = U, v = V cos(omega t). It is the same at every
    point (no velocity gradient) but unsteady, which is enough to separate the three lines."""
    def vel(x, y, t):
        x = np.asarray(x, dtype=float)
        return np.full_like(x, U), np.full_like(x, V * np.cos(omega * t)) + 0.0 * np.asarray(y, dtype=float)
    return vel


def stagnation_flow(a: float) -> Field:
    """Plane stagnation-point flow u = a x, v = -a y (steady; streamlines are the hyperbolas xy = const)."""
    def vel(x, y, t):
        return a * np.asarray(x, dtype=float), -a * np.asarray(y, dtype=float)
    return vel


def shear_flow(gamma: float) -> Field:
    """Simple shear u = gamma y, v = 0 (the flow between sliding plates, `P08.06`)."""
    def vel(x, y, t):
        y = np.asarray(y, dtype=float)
        return gamma * y, 0.0 * y + 0.0 * np.asarray(x, dtype=float)
    return vel


def solid_body_rotation(Omega: float) -> Field:
    """Rigid rotation at angular velocity Omega: u = -Omega y, v = Omega x."""
    def vel(x, y, t):
        return -Omega * np.asarray(y, dtype=float), Omega * np.asarray(x, dtype=float)
    return vel


def taylor_green(nu: float, k: float = 1.0) -> Field:
    """Decaying Taylor-Green vortex array u = sin(kx) cos(ky) F(t), v = -cos(kx) sin(ky) F(t),
    F = exp(-2 nu k^2 t): unsteady, but the direction of u at each point never changes."""
    def vel(x, y, t):
        F = np.exp(-2.0 * nu * k ** 2 * t)
        x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        return np.sin(k * x) * np.cos(k * y) * F, -np.cos(k * x) * np.sin(k * y) * F
    return vel


def pulsating_stagnation_flow(a: float, eps: float, omega: float) -> Field:
    """A stagnation flow whose strength pulses, with an oscillating shear added:
    u = a(t) x + 0.3 a cos(omega t) y, v = -a(t) y, a(t) = a (1 + eps sin omega t).
    Non-uniform, unsteady and divergence-free: a field in which both the local and the convective
    parts of the material derivative are nonzero."""
    def vel(x, y, t):
        at = a * (1.0 + eps * np.sin(omega * t))
        x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        return at * x + 0.3 * a * y * np.cos(omega * t), -at * y
    return vel


# ---------------------------------------------------------------- integrating the lines

def rk4_step(vel: Field, x: np.ndarray, y: np.ndarray, t: float, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """One classical RK4 step of dx/dt = u(x, y, t), dy/dt = v(x, y, t) for arrays of particles."""
    u1, v1 = vel(x, y, t)
    u2, v2 = vel(x + dt / 2 * u1, y + dt / 2 * v1, t + dt / 2)
    u3, v3 = vel(x + dt / 2 * u2, y + dt / 2 * v2, t + dt / 2)
    u4, v4 = vel(x + dt * u3, y + dt * v3, t + dt)
    return x + dt / 6 * (u1 + 2 * u2 + 2 * u3 + u4), y + dt / 6 * (v1 + 2 * v2 + 2 * v3 + v4)


def pathline(vel: Field, x0: float, y0: float, t0: float, t1: float, n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Track of the particle that is at (x0, y0) at time t0, from t0 to t1 in n RK4 steps:
    dx/dt = u(x, t). Returns (t, x, y), each of length n + 1."""
    t = np.linspace(t0, t1, n + 1)
    dt = (t1 - t0) / n
    x, y = np.empty(n + 1), np.empty(n + 1)
    x[0], y[0] = x0, y0
    for i in range(n):
        x[i + 1], y[i + 1] = rk4_step(vel, x[i], y[i], t[i], dt)
    return t, x, y


def streakline(vel: Field, xs: float, ys: float, t_start: float, t_now: float, n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Streakline at time t_now from a source at (xs, ys) that has released a particle every
    dt = (t_now - t_start)/n since t_start (a smoke wand). All released particles are advanced
    together with RK4; a new one is added at the source after each step.
    Returns (tau, x, y): release times (oldest first) and the particles' positions at t_now."""
    dt = (t_now - t_start) / n
    x, y = np.array([float(xs)]), np.array([float(ys)])
    for i in range(n):
        t = t_start + i * dt
        x, y = rk4_step(vel, x, y, t, dt)
        x, y = np.append(x, xs), np.append(y, ys)
    tau = t_start + dt * np.arange(n + 1)
    return tau, x, y


def streamline(vel: Field, x0: float, y0: float, t: float, length: float, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Streamline through (x0, y0) at the frozen instant t, traced for an arc length ``length``
    in n RK4 steps of dx/ds = u/|u| (arc-length form, so the speed does not matter, only the
    direction). Returns (x, y) of length n + 1."""
    def unit(xx, yy, _s):
        u, v = vel(xx, yy, t)
        speed = np.hypot(u, v)
        if np.any(speed == 0):
            raise ValueError("streamline reached a stagnation point")
        return u / speed, v / speed
    ds = length / n
    x, y = np.empty(n + 1), np.empty(n + 1)
    x[0], y[0] = x0, y0
    for i in range(n):
        x[i + 1], y[i + 1] = rk4_step(unit, x[i], y[i], 0.0, ds)
    return x, y


# ---------------------------------------------------------------- closed forms for the oscillating stream

def oscillating_pathline_exact(U: float, V: float, omega: float, x0: float, y0: float, t0: float,
                               t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Exact pathline of u = U, v = V cos(omega t) from (x0, y0) at t0:
    x = x0 + U (t - t0), y = y0 + (V/omega)(sin omega t - sin omega t0)."""
    t = np.asarray(t, dtype=float)
    return x0 + U * (t - t0), y0 + V / omega * (np.sin(omega * t) - np.sin(omega * t0))


def oscillating_streakline_exact(U: float, V: float, omega: float, xs: float, ys: float, t_now: float,
                                 tau: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Exact streakline at t_now from a source at (xs, ys): the particle released at time tau is at
    x = xs + U (t_now - tau), y = ys + (V/omega)(sin omega t_now - sin omega tau). Eliminating tau,
    y - ys = (V/omega)[sin omega t_now - sin omega (t_now - (x - xs)/U)]: a sine wave of wavelength
    2 pi U/omega frozen into the smoke."""
    return oscillating_pathline_exact(U, V, omega, xs, ys, tau, t_now)


# ---------------------------------------------------------------- the velocity gradient

def velocity_gradient(vel: Field, x: float, y: float, t: float, h: float = 1e-5) -> np.ndarray:
    """Velocity-gradient tensor G_ij = d u_i / d x_j at (x, y, t) by central differences of step h:
    [[du/dx, du/dy], [dv/dx, dv/dy]]."""
    up, vp = vel(np.array(x + h), np.array(y), t)
    um, vm = vel(np.array(x - h), np.array(y), t)
    uq, vq = vel(np.array(x), np.array(y + h), t)
    ur, vr = vel(np.array(x), np.array(y - h), t)
    return np.array([[(up - um) / (2 * h), (uq - ur) / (2 * h)],
                     [(vp - vm) / (2 * h), (vq - vr) / (2 * h)]], dtype=float).reshape(2, 2)


def decompose_gradient(G: np.ndarray) -> dict[str, np.ndarray | float]:
    """Split G into the strain-rate tensor S = (G + G^T)/2 and the rotation tensor W = (G - G^T)/2,
    and return them with the divergence tr G = du/dx + dv/dy (rate of area change per unit area),
    the vorticity omega = dv/dx - du/dy = 2 W_21, and the principal strain rates (eigenvalues of S)."""
    G = np.asarray(G, dtype=float)
    S = 0.5 * (G + G.T)
    W = 0.5 * (G - G.T)
    return {"S": S, "W": W, "divergence": float(np.trace(G)), "vorticity": float(G[1, 0] - G[0, 1]),
            "principal_strain_rates": np.linalg.eigvalsh(S)}


def material_derivative(f: Callable[[np.ndarray, np.ndarray, float], np.ndarray], vel: Field,
                        x: float, y: float, t: float, h: float = 1e-5) -> float:
    """Df/Dt = df/dt + u df/dx + v df/dy of a scalar field f(x, y, t) at (x, y, t), by central
    differences of step h in x, y and t: the rate of change of f following the particle there."""
    u, v = vel(np.array(x), np.array(y), t)
    dfdt = (f(x, y, t + h) - f(x, y, t - h)) / (2 * h)
    dfdx = (f(x + h, y, t) - f(x - h, y, t)) / (2 * h)
    dfdy = (f(x, y + h, t) - f(x, y - h, t)) / (2 * h)
    return float(dfdt + float(u) * dfdx + float(v) * dfdy)


def particle_acceleration(vel: Field, x: float, y: float, t: float, h: float = 1e-5) -> np.ndarray:
    """Acceleration of the fluid particle at (x, y, t): Du/Dt = du/dt + (u . grad) u, both
    components by ``material_derivative``. Nonzero in a steady flow whenever the speed or the
    direction changes along the streamline (convective acceleration)."""
    fu = lambda xx, yy, tt: float(vel(np.array(xx), np.array(yy), tt)[0])
    fv = lambda xx, yy, tt: float(vel(np.array(xx), np.array(yy), tt)[1])
    return np.array([material_derivative(fu, vel, x, y, t, h), material_derivative(fv, vel, x, y, t, h)])


# ---------------------------------------------------------------- helpers (given)

def polygon_area(x: np.ndarray, y: np.ndarray) -> float:
    """Signed area of a closed polygon by the shoelace formula (counter-clockwise positive)."""
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def distance_to_curve(px: np.ndarray, py: np.ndarray, cx: np.ndarray, cy: np.ndarray) -> np.ndarray:
    """For each point (px, py), the distance to the nearest vertex of the densely sampled curve (cx, cy)."""
    d = np.hypot(np.asarray(px)[:, None] - np.asarray(cx)[None, :], np.asarray(py)[:, None] - np.asarray(cy)[None, :])
    return d.min(axis=1)


def observed_order(errors: np.ndarray, ratio: float = 2.0) -> np.ndarray:
    """Observed convergence orders log(e_k / e_{k+1}) / log(ratio) for steps refined by ``ratio``."""
    e = np.asarray(errors, dtype=float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


def smoke_error(U: float, V: float, omega: float, L: float, t_now: float, n: int = 400) -> float:
    """Largest lateral distance, over 0 <= x - xs <= L, between the smoke line (streakline) from a
    source and the instantaneous streamline through the same source, for the oscillating stream,
    divided by the largest lateral excursion of that streamline, V L / U. Uses the closed forms."""
    xi = np.linspace(0.0, L, n + 1)
    y_streak = V / omega * (np.sin(omega * t_now) - np.sin(omega * (t_now - xi / U)))
    y_stream = V * np.cos(omega * t_now) / U * xi
    return float(np.max(np.abs(y_streak - y_stream)) / (V * L / U))


# ---------------------------------------------------------------- the lesson's run

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save flow_lines.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    U, f_gust, ang = 10.0, 2.0, 5.0
    V, omega = U * np.tan(np.radians(ang)), 2 * np.pi * f_gust
    print(f"gust tunnel: U = {U} m/s, flow angle +/- {ang} deg (V = {V:.4f} m/s), {f_gust} Hz (omega = {omega:.3f} rad/s)")
    print(f"  streakline wavelength 2 pi U/omega = {2 * np.pi * U / omega:.3f} m, pathline amplitude V/omega = {V / omega * 1e3:.2f} mm")
    vel = oscillating_flow(U, V, omega)
    t_now = 1.0
    tau, xs, ys = streakline(vel, 0.0, 0.0, 0.0, t_now, 200)
    xe, ye = oscillating_streakline_exact(U, V, omega, 0.0, 0.0, t_now, tau)
    print(f"  streakline from RK4 (200 releases) vs closed form: max error {np.max(np.hypot(xs - xe, ys - ye)):.2e} m")
    tp, xp, yp = pathline(vel, 0.0, 0.0, 0.0, t_now, 200)
    xpe, ype = oscillating_pathline_exact(U, V, omega, 0.0, 0.0, 0.0, tp)
    print(f"  pathline from RK4 (200 steps) vs closed form: max error {np.max(np.hypot(xp - xpe, yp - ype)):.2e} m")
    errs = []
    for n in (10, 20, 40, 80):
        tp, xp, yp = pathline(vel, 0.0, 0.0, 0.0, 0.3, n)
        errs.append(abs(yp[-1] - oscillating_pathline_exact(U, V, omega, 0, 0, 0, 0.3)[1]))
    print(f"  pathline to t = 0.3 s, RK4 errors {np.array2string(np.array(errs), precision=2)}, orders {np.array2string(observed_order(errs), precision=2)}")
    print("  smoke-line error / streamline excursion over a model of length L (worst phase):")
    print("     L (m)   k = omega L / U   error")
    for L in (0.05, 0.15, 0.3, 1.5, 5.0):
        k = omega * L / U
        e = max(smoke_error(U, V, omega, L, tt) for tt in np.linspace(0, 1 / f_gust, 61))
        print(f"   {L:6.2f}      {k:7.3f}       {e:6.3f}")

    # steady stagnation flow: all three coincide
    stag = stagnation_flow(1.0)
    _, xs2, ys2 = streakline(stag, 0.2, 1.0, 0.0, 2.0, 200)
    _, xp2, yp2 = pathline(stag, 0.2, 1.0, 0.0, 2.0, 200)
    xl2, yl2 = streamline(stag, 0.2, 1.0, 0.0, 2.0, 400)
    print(f"\nsteady stagnation flow a = 1/s: x y along streakline {np.ptp(xs2 * ys2):.1e} spread, pathline {np.ptp(xp2 * yp2):.1e}, "
          f"streamline {np.ptp(xl2 * yl2):.1e} (all on x y = 0.2)")

    print("\nvelocity gradients (central differences):")
    for name, fld in (("stagnation a = 1/s", stagnation_flow(1.0)), ("shear gamma = 1/s", shear_flow(1.0)),
                      ("rotation Omega = 1/s", solid_body_rotation(1.0))):
        d = decompose_gradient(velocity_gradient(fld, 0.3, 0.4, 0.0))
        print(f"  {name:22s} div {d['divergence']:+.3f}  vorticity {d['vorticity']:+.3f}  principal strain rates "
              f"{np.array2string(d['principal_strain_rates'], precision=3)}")
    a_p = particle_acceleration(stagnation_flow(1.0), 2.0, 1.0, 0.0)
    print(f"  particle acceleration in steady stagnation flow at (2, 1): {np.array2string(a_p, precision=4)} m/s^2 (a^2 x, a^2 y)")

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    xl, yl = streamline(vel, 0.0, 0.0, t_now, 10.0, 200)
    ax[0].plot(xs, ys * 1e3, label="streakline (smoke from the origin)")
    tp, xp, yp = pathline(vel, 0.0, 0.0, 0.0, t_now, 200)
    ax[0].plot(xp, yp * 1e3, "--", label="pathline (particle released at t = 0)")
    ax[0].plot(xl, yl * 1e3, ":", label="streamline at t = 1 s")
    ax[0].set_xlabel("x (m)")
    ax[0].set_ylabel("y (mm)")
    ax[0].set_title("oscillating stream: three different lines")
    ax[0].legend(fontsize=8)
    for c in (0.05, 0.2, 0.5):
        ax[1].plot(*streamline(stag, c / 1.0 + 0.0, 1.0, 0.0, 3.0, 200), "k-", lw=0.8)
    ax[1].plot(xs2, ys2, "o", ms=2, label="streakline")
    ax[1].plot(xp2, yp2, "-", lw=2, alpha=0.5, label="pathline")
    ax[1].set_xlabel("x (m)")
    ax[1].set_ylabel("y (m)")
    ax[1].set_title("steady stagnation flow: lines coincide")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "flow_lines.png", dpi=110)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'flow_lines.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
