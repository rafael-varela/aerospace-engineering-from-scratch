"""P08.09 Vorticity & Circulation: starter code.

When the light aeroplane starts to fly, Kelvin's theorem says the circulation round a loop of air
particles that once enclosed the still wing cannot change: the wing's bound circulation must be
balanced by a starting vortex left behind on the runway. How strong is that vortex, how quickly would
viscosity alone spread and weaken it, and what downwash does the wing's vortex system induce at the
tailplane? This module builds the tools and checks each by an independent route:

* the Lamb-Oseen vortex, an exact solution of the Navier-Stokes equations: swirl
  u_theta = Gamma/(2 pi r) (1 - exp(-r^2/rc2)), vorticity Gamma/(pi rc2) exp(-r^2/rc2), with the
  squared core radius rc2 = rc0^2 + 4 nu t growing by diffusion (with nu = 0 it is a steady Euler
  solution with a finite core);
* circulation as a line integral round a closed polygon (midpoint rule) and as the area integral of
  vorticity over a disk (Gauss-Legendre in r, trapezoid in theta): Stokes's theorem;
* a material loop advected through the vortex two ways, by RK4 and exactly (every particle stays on
  its circle and turns through an angle known in closed form with the exponential integral E1), so
  Kelvin's theorem can be tested: constant circulation when nu = 0, and the Lamb-Oseen decay when
  nu > 0, where the circulation of the material loop equals the area integral of the current
  vorticity over the material region;
* the Biot-Savart law for a straight vortex segment, a horseshoe vortex for the wing, and plane
  point vortices for the bound and starting vortices.

Sign conventions: in this lesson Gamma is positive COUNTERCLOCKWISE (the mathematical sense used by
Stokes's and Kelvin's theorems), and vorticity omega_z = dv/dx - du/dy. The aerodynamic lift of a
wing in flow from left to right needs a clockwise circulation, i.e. Gamma < 0 here (P08.08 counted
clockwise as positive). In 3-D a filament's Gamma follows the right-hand rule about its direction.

Units: SI (m, s, m^2/s).

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes vorticity.png)
Test: pytest curriculum/phase-08/lesson-09-vorticity-circulation
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.special import exp1

G0 = 9.80665
RHO_3000 = 0.90913                                   # kg/m^3, ISA at 3000 m (P08.01)
MU_3000 = 1.458e-6 * 268.65 ** 1.5 / (268.65 + 110.4)  # Pa s, Sutherland at the ISA 3000 m temperature
NU_3000 = MU_3000 / RHO_3000                         # m^2/s
R_MAX_FACTOR = 1.1209                                # radius of peak swirl = 1.1209 sqrt(rc2)


# ---------------------------------------------------------------- the Lamb-Oseen vortex

def core_radius2(rc0: float, nu: float, t: float) -> float:
    """Squared core radius of a diffusing Lamb-Oseen vortex, rc^2 = rc0^2 + 4 nu t."""
    return rc0 ** 2 + 4.0 * nu * t


def lamb_oseen_utheta(r: np.ndarray, Gamma: float, rc2: float) -> np.ndarray:
    """Swirl speed u_theta = Gamma / (2 pi r) (1 - exp(-r^2 / rc2)); finite (zero) on the axis and
    tending to the point vortex Gamma/(2 pi r) for r >> rc."""
    raise NotImplementedError


def lamb_oseen_vorticity(r: np.ndarray, Gamma: float, rc2: float) -> np.ndarray:
    """Axial vorticity omega = Gamma / (pi rc2) exp(-r^2 / rc2), a Gaussian whose integral is Gamma."""
    raise NotImplementedError


def lamb_oseen_circulation(R: np.ndarray, Gamma: float, rc2: float) -> np.ndarray:
    """Circulation round a circle of radius R centred on the vortex, Gamma (1 - exp(-R^2 / rc2))."""
    raise NotImplementedError


def vortex_velocity(xy: np.ndarray, Gamma: float, rc2: float, center: tuple[float, float] = (0.0, 0.0)) -> np.ndarray:
    """Cartesian velocity (N, 2) of a Lamb-Oseen vortex centred at ``center`` at the points xy (N, 2)."""
    raise NotImplementedError


# ---------------------------------------------------------------- circulation two ways

def circulation_polygon(loop: np.ndarray, velocity: Callable[[np.ndarray], np.ndarray]) -> float:
    """Circulation round the closed polygon ``loop`` (N, 2), counterclockwise positive, by the
    midpoint rule: sum of u(midpoint) . (x_{k+1} - x_k). Second order in the segment length."""
    raise NotImplementedError


def disk_area_integral(f: Callable[[np.ndarray], np.ndarray], center: tuple[float, float], radius: float,
                       nr: int = 48, nt: int = 96) -> float:
    """Integral of f(xy) over a disk: Gauss-Legendre with nr points in r (with the Jacobian r) and the
    periodic trapezoid rule with nt points in theta; spectrally accurate for smooth f."""
    raise NotImplementedError


def circle_loop(center: tuple[float, float], radius: float, n: int) -> np.ndarray:
    """n points counterclockwise on a circle, as an (n, 2) array."""
    th = 2 * np.pi * np.arange(n) / n
    return np.column_stack([center[0] + radius * np.cos(th), center[1] + radius * np.sin(th)])


# ---------------------------------------------------------------- advecting a material loop

def turning_angle(r: np.ndarray, t: float, Gamma: float, rc0: float, nu: float) -> np.ndarray:
    """Angle a particle at radius r turns through in time t (it stays on its circle, since the flow is
    purely azimuthal): the integral of u_theta(r, s)/r ds from 0 to t. With nu = 0 it is
    Omega(r) t; with nu > 0, substituting w = r^2/rc2(s),

        theta = Gamma/(2 pi r^2) [t - r^2/(4 nu) (g(w_0) - g(w_t))],   g(w) = E1(w) - exp(-w)/w,

    where w_0 = r^2/rc0^2 and w_t = r^2/(rc0^2 + 4 nu t)."""
    raise NotImplementedError


def material_loop_exact(loop0: np.ndarray, t: float, Gamma: float, rc0: float, nu: float) -> np.ndarray:
    """Positions at time t of the particles that were at loop0 (N, 2) at t = 0 in a Lamb-Oseen vortex
    centred at the origin: each turns about the origin by ``turning_angle``."""
    raise NotImplementedError


def advect_rk4(xy0: np.ndarray, t0: float, t1: float, nt: int, Gamma: float, rc0: float, nu: float) -> np.ndarray:
    """Advect the points xy0 (N, 2) through the (possibly diffusing) Lamb-Oseen vortex at the origin
    from t0 to t1 with nt classical RK4 steps (fourth order in the step)."""
    raise NotImplementedError


def material_circulation_by_area(center0: tuple[float, float], radius0: float, t: float, Gamma: float,
                                 rc0: float, nu: float, nr: int = 64, nt: int = 128) -> float:
    """Circulation at time t of the material loop that was the circle (center0, radius0) at t = 0, by
    Stokes's theorem: the area integral of the CURRENT vorticity over the CURRENT material region. The
    flow turns each circle r = const rigidly, so the material region holds, at every r, the same arc
    lengths as the initial disk and the vorticity depends on r only: the integral over the current
    region equals the integral of omega(r, t) over the initial disk."""
    raise NotImplementedError


# ---------------------------------------------------------------- Biot-Savart

def biot_savart_segment(P: np.ndarray, A: np.ndarray, B: np.ndarray, Gamma: float) -> np.ndarray:
    """Velocity induced at the points P (N, 3) by a straight vortex filament from A to B of strength
    Gamma (right-hand rule about A -> B):

        u = Gamma/(4 pi) (r1 x r2)/|r1 x r2|^2 (r0 . (r1/|r1| - r2/|r2|)),  r0 = B - A, r1 = P - A, r2 = P - B.

    Points on the filament's line get zero velocity (the singular self-induction is dropped)."""
    raise NotImplementedError


def horseshoe_velocity(P: np.ndarray, b_prime: float, Gamma: float, length: float = 1e5) -> np.ndarray:
    """Velocity at P (N, 3) induced by a horseshoe vortex for a wing with the free stream along +x and
    lift along +z: trailing leg from (length, -b'/2, 0) to (0, -b'/2, 0), bound vortex to
    (0, +b'/2, 0) (Gamma along +y, as L = rho V x Gamma requires), trailing leg back to (length, +b'/2, 0).
    ``length`` stands for infinity."""
    raise NotImplementedError


def point_vortices_velocity(xy: np.ndarray, vortices: list[tuple[float, float, float]]) -> np.ndarray:
    """Plane velocity (N, 2) at xy of point vortices given as (x, y, Gamma), counterclockwise positive."""
    raise NotImplementedError


def observed_order(errors: np.ndarray, ratio: float = 2.0) -> np.ndarray:
    """Observed convergence orders log(e_k / e_{k+1}) / log(ratio) for refinement by ``ratio``."""
    e = np.asarray(errors, dtype=float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save vorticity.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    # the aeroplane
    mass, span, V, rho = 750.0, 10.0, 60.0, RHO_3000
    L = mass * G0
    b_prime = np.pi * span / 4                       # vortex span of an elliptic loading
    Gam0 = L / (rho * V * b_prime)
    print(f"light aeroplane: lift {L:.1f} N at {V} m/s, rho {rho} kg/m^3, span {span} m")
    print(f"  mean 2-D circulation L/(rho V b) = {L / (rho * V * span):.3f} m^2/s; "
          f"horseshoe (b' = pi b/4 = {b_prime:.4f} m): Gamma0 = {Gam0:.3f} m^2/s")
    print(f"  starting vortex: +{Gam0:.3f} m^2/s (counterclockwise seen from the left wing), bound vortex -{Gam0:.3f}")
    big = circle_loop((25.0, 0.0), 40.0, 4000)
    wing_only = circle_loop((0.0, 0.0), 2.0, 400)
    vort = [(0.0, 0.0, -Gam0), (50.0, 0.0, +Gam0)]
    vel = lambda p: point_vortices_velocity(p, vort)
    print(f"  Kelvin, plane model: circulation round both vortices {circulation_polygon(big, vel):.2e} m^2/s; "
          f"round the wing alone {circulation_polygon(wing_only, vel):.4f} m^2/s")
    w_wing = point_vortices_velocity(np.array([[0.0, 0.0]]), [(50.0, 0.0, Gam0)])[0, 1]
    print(f"  starting vortex 50 m behind induces {w_wing:.4f} m/s at the wing (2-D, Gamma/(2 pi d))")

    # viscous spreading of the starting vortex, laminar
    print(f"\nair at 3000 m: mu = {MU_3000:.4e} Pa s, nu = {NU_3000:.4e} m^2/s")
    for t in (1.0, 10.0, 60.0, 600.0):
        rc2 = core_radius2(0.0, NU_3000, t)
        rmax = R_MAX_FACTOR * np.sqrt(rc2)
        umax = float(lamb_oseen_utheta(np.array([rmax]), Gam0, rc2)[0])
        print(f"  t = {t:6.0f} s: core radius of peak swirl {rmax * 100:.2f} cm, peak swirl {umax:.2f} m/s, "
              f"circulation within 1 m {float(lamb_oseen_circulation(1.0, Gam0, rc2)) / Gam0 * 100:.4f} % of Gamma0")
    t_half = 1.0 / (4 * NU_3000 * np.log(2.0))
    print(f"  laminar diffusion takes {t_half:.0f} s ({t_half / 3600:.1f} h) to halve the circulation within 1 m")

    # downwash at the tail
    l_t = 4.5
    for x in (0.0, l_t, 50.0, 1e4):
        w = horseshoe_velocity(np.array([[x, 0.0, 0.0]]), b_prime, Gam0)[0, 2]
        print(f"  horseshoe downwash on the centreline {x:8.1f} m behind the bound vortex: w = {w:.4f} m/s, "
              f"angle {np.degrees(np.arctan(-w / V)):.3f} deg")
    print(f"  formulas: at the bound vortex -Gamma0/(pi b') = {-Gam0 / (np.pi * b_prime):.4f}; "
          f"far behind -2 Gamma0/(pi b') = {-2 * Gam0 / (np.pi * b_prime):.4f} m/s")
    CL = L / (0.5 * rho * V ** 2 * span * 1.5)
    AR = span / 1.5
    print(f"  C_L = {CL:.4f}, AR = {AR:.3f}: elliptic-wing estimate of the far downwash angle 2 C_L/(pi AR) = "
          f"{np.degrees(2 * CL / (np.pi * AR)):.3f} deg")

    # Kelvin's theorem: a material loop off the vortex centre (nondimensional: Gamma = 1, rc0 = 1)
    print("\nKelvin's theorem, material loop: circle radius 0.6 centred at (0.9, 0); Gamma = 1, rc0 = 1")
    c0, R0 = (0.9, 0.0), 0.6
    loop0 = circle_loop(c0, R0, 4000)
    G_init = circulation_polygon(loop0, lambda p: vortex_velocity(p, 1.0, 1.0))
    G_area = disk_area_integral(lambda p: lamb_oseen_vorticity(np.hypot(p[:, 0], p[:, 1]), 1.0, 1.0), c0, R0)
    print(f"  Stokes at t = 0: line integral {G_init:.10f}, area integral {G_area:.10f}")
    T_turn = 2 * np.pi / float(lamb_oseen_utheta(np.array([0.9]), 1.0, 1.0)[0] / 0.9)
    print(f"  turnover time at r = 0.9: {T_turn:.3f}")
    for nu in (0.0, 0.01):
        rows = []
        for t in (0.0, 1.0, 2.0, 4.0) * np.array(T_turn):
            lp = material_loop_exact(loop0, t, 1.0, 1.0, nu)
            rc2 = core_radius2(1.0, nu, t)
            G_line = circulation_polygon(lp, lambda p: vortex_velocity(p, 1.0, rc2))
            G_a = material_circulation_by_area(c0, R0, t, 1.0, 1.0, nu)
            length = np.sum(np.linalg.norm(np.roll(lp, -1, 0) - lp, axis=1))
            rows.append((t, G_line, G_a, length))
        print(f"  nu = {nu}:  t, line integral, area integral, loop length")
        for t, gl, ga, ln in rows:
            print(f"    {t:7.3f}  {gl:.8f}  {ga:.8f}  {ln:.3f}")
    errs = []
    t1 = T_turn
    exact = material_loop_exact(loop0[:200:10], t1, 1.0, 1.0, 0.0)
    for nt in (25, 50, 100, 200):
        errs.append(np.max(np.linalg.norm(advect_rk4(loop0[:200:10], 0.0, t1, nt, 1.0, 1.0, 0.0) - exact, axis=1)))
    print(f"  RK4 against the exact particle paths over one turnover: errors {np.array2string(np.array(errs), formatter={'float': lambda v: f'{v:.2e}'})}, "
          f"orders {np.array2string(observed_order(errs), precision=2)}")
    lp_rk = advect_rk4(loop0, 0.0, 2 * T_turn, 400, 1.0, 1.0, 0.0)
    G_rk = circulation_polygon(lp_rk, lambda p: vortex_velocity(p, 1.0, 1.0))
    print(f"  RK4-advected loop (400 steps, 2 turnovers): circulation {G_rk:.8f}, change {G_rk - G_init:.1e}")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    for t, col in zip((0.0, 1.0, 2.0, 4.0), ("k", "C0", "C1", "C2")):
        lp = material_loop_exact(loop0, t * T_turn, 1.0, 1.0, 0.0)
        ax[0].plot(*np.vstack([lp, lp[:1]]).T, col, lw=0.8, label=f"t = {t:g} turnovers")
    ax[0].plot(0, 0, "k+")
    ax[0].set_aspect("equal")
    ax[0].set_title("material loop in a vortex (nu = 0)")
    ax[0].legend(fontsize=7)
    ts = np.linspace(0, 4 * T_turn, 40)
    for nu, mk in ((0.0, "o"), (0.01, "s")):
        Gl = [circulation_polygon(material_loop_exact(loop0, t, 1.0, 1.0, nu),
                                  lambda p, t=t, nu=nu: vortex_velocity(p, 1.0, core_radius2(1.0, nu, t))) for t in ts[::4]]
        Ga = [material_circulation_by_area(c0, R0, t, 1.0, 1.0, nu) for t in ts]
        ax[1].plot(ts, Ga, "-", label=f"area integral, nu = {nu}")
        ax[1].plot(ts[::4], Gl, mk, ms=4, label=f"line integral, nu = {nu}")
    ax[1].set_xlabel("t (Gamma = 1, rc0 = 1 units)")
    ax[1].set_ylabel("circulation of the material loop")
    ax[1].set_title("Kelvin (nu = 0) and viscous decay")
    ax[1].legend(fontsize=7)
    xs = np.linspace(-2, 20, 300)
    w = horseshoe_velocity(np.column_stack([xs, np.zeros_like(xs), np.zeros_like(xs)]), b_prime, Gam0)[:, 2]
    ax[2].plot(xs, w)
    ax[2].axvline(l_t, color="k", lw=0.8, ls="--", label="tailplane")
    ax[2].axhline(-2 * Gam0 / (np.pi * b_prime), color="k", lw=0.8, ls=":", label="-2 Gamma/(pi b')")
    ax[2].set_xlabel("x behind the bound vortex (m)")
    ax[2].set_ylabel("w on the centreline (m/s)")
    ax[2].set_title("horseshoe vortex downwash")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "vorticity.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'vorticity.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
