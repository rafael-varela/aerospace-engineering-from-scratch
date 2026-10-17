"""P06.05 Rigid Body Kinematics in 2D: starter code.

A servo arm drives a control surface through a pushrod and a horn on the hinge. That is a four-bar
linkage: ground (servo shaft to hinge line), crank (servo arm), coupler (pushrod), rocker (horn). How far
does the surface deflect over the servo's travel, how linear is the gearing, how close does the linkage
come to locking up, and what servo torque does a hinge moment demand? This module analyses any four-bar:

* position: the loop closure r2 e^(i th2) + r3 e^(i th3) - r4 e^(i th4) - r1 = 0, two nonlinear equations
  for th3 and th4, solved by Newton's method, and independently by intersecting two circles;
* velocity and acceleration: the relative-velocity equations of rigid bodies, v_B = v_A + omega3 x r_BA
  and a_B = a_A + alpha3 x r_BA - omega3^2 r_BA, which are linear in the unknown rates;
* instant centres (Kennedy's theorem) as a third route to the velocity ratio;
* the Grashof criterion and the transmission angle.

The independent check: finite differences of the position solution, taken along a crank motion,
converge at second order to the analytic angular velocities and accelerations.

You write every function whose body is `raise NotImplementedError`. `_cross2`, `_intersect`,
`sweep_positions`, `central_rates`, `observed_order` and `main` are given.

Convention: O2 at the origin, O4 at (r1, 0); angles counterclockwise from +x; link vectors O2->A (r2),
A->B (r3), O4->B (r4). Units: m, rad, s (the lesson prints mm and degrees).

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes four_bar.png)
Test: pytest curriculum/phase-06/lesson-05-rigid-body-kinematics-2d
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# the lesson's control-surface linkage (illustrative dimensions for a small UAV)
SERVO_LINKS = (0.120, 0.020, 0.120, 0.025)   # m: ground r1, servo arm r2, pushrod r3, horn r4


def _cross2(w: float, r: np.ndarray) -> np.ndarray:
    """omega e_z x r for a planar vector r."""
    return w * np.array([-r[1], r[0]])


# ---------------------------------------------------------------- classification

def grashof(r1: float, r2: float, r3: float, r4: float) -> str:
    """Grashof classification. With s the shortest, l the longest and p, q the others:
    s + l < p + q (Grashof) gives a link that turns fully: 'crank-rocker' if the shortest is the input
    crank r2, 'rocker-crank' if it is the output r4, 'double-crank' if it is the ground r1, and
    'double-rocker' if it is the coupler r3. s + l = p + q is a 'change-point' linkage (it can fold flat);
    s + l > p + q is a non-Grashof 'triple-rocker' (no link turns fully)."""
    raise NotImplementedError


# ---------------------------------------------------------------- position

def loop_residual(th2: float, th3: float, th4: float, links) -> np.ndarray:
    """Loop-closure residual f = r2 e^(i th2) + r3 e^(i th3) - r4 e^(i th4) - r1 as (real, imaginary)."""
    raise NotImplementedError


def loop_jacobian(th3: float, th4: float, links) -> np.ndarray:
    """Jacobian of ``loop_residual`` with respect to (th3, th4):
    [[-r3 sin th3, r4 sin th4], [r3 cos th3, -r4 cos th4]]. Its determinant r3 r4 sin(th3 - th4) vanishes
    when the coupler and the rocker are in line (a toggle position)."""
    raise NotImplementedError


def solve_position_newton(th2: float, links, guess: tuple[float, float], tol: float = 1e-13,
                          maxit: int = 50) -> tuple[float, float, int]:
    """Newton's method on the loop closure: [th3, th4] <- [th3, th4] - J^-1 f until |f| < tol * r1.
    Returns (th3, th4, iterations). Raises RuntimeError if it does not converge (no assembly, or a
    guess on the wrong side of a toggle)."""
    raise NotImplementedError


def position_circle_intersection(th2: float, links, branch: int = 1) -> tuple[float, float]:
    """Closed-form position by geometry: B lies on the circle of radius r3 about A and on the circle of
    radius r4 about O4. With d = |O4 - A|, the chord's foot is a = (r3^2 - r4^2 + d^2) / (2d) along A->O4
    and its half-length h = sqrt(r3^2 - a^2); ``branch`` +1 puts B to the left of A->O4, -1 to the right.
    Raises ValueError if the circles do not meet (the linkage cannot be assembled at this th2)."""
    raise NotImplementedError


def joint_positions(th2: float, th3: float, th4: float, links) -> dict:
    """Positions of O2, A, B, O4 (planar arrays)."""
    raise NotImplementedError


def sweep_positions(th2: np.ndarray, links, guess: tuple[float, float]) -> tuple[np.ndarray, np.ndarray]:
    """Solve the position along a sequence of crank angles by continuation: each Newton solve starts
    from the previous answer, which keeps the linkage on one branch. Returns (th3, th4) arrays."""
    th3, th4 = np.empty(len(th2)), np.empty(len(th2))
    g = guess
    for i, t in enumerate(th2):
        a, b, _ = solve_position_newton(float(t), links, g)
        th3[i], th4[i] = a, b
        g = (a, b)
    return th3, th4


# ---------------------------------------------------------------- velocity and acceleration

def velocity_analysis(th2: float, th3: float, th4: float, w2: float, links) -> tuple[float, float]:
    """Angular velocities of coupler and rocker. The relative-velocity equation
    v_B = v_A + omega3 x r_BA = omega4 x r_B/O4, with v_A = omega2 x r_A, differentiated loop closure:
    J [w3, w4]^T = [r2 w2 sin th2, -r2 w2 cos th2]^T. Returns (w3, w4)."""
    raise NotImplementedError


def acceleration_analysis(th2: float, th3: float, th4: float, w2: float, al2: float, w3: float, w4: float,
                          links) -> tuple[float, float]:
    """Angular accelerations of coupler and rocker from the relative-acceleration equation
    a_B = a_A + alpha3 x r_BA - w3^2 r_BA = alpha4 x r_B/O4 - w4^2 r_B/O4, i.e. the loop closure
    differentiated twice: sum of r (i alpha - w^2) e^(i th) = 0 around the loop. Linear in
    (alpha3, alpha4) with the same Jacobian as the velocity. Returns (alpha3, alpha4)."""
    raise NotImplementedError


def coupler_point(th2: float, th3: float, w2: float, al2: float, w3: float, al3: float, links,
                  p: float, delta: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Position, velocity and acceleration of a point P on the coupler, a distance p from A at angle
    delta from A->B, by the rigid-body relations: r_P = r_A + r_PA, v_P = v_A + w3 x r_PA,
    a_P = a_A + al3 x r_PA - w3^2 r_PA, with v_A = w2 x r_A and a_A = al2 x r_A - w2^2 r_A."""
    raise NotImplementedError


# ---------------------------------------------------------------- instant centres

def _intersect(p1, p2, q1, q2) -> np.ndarray:
    """Intersection of the line through p1, p2 with the line through q1, q2 (inf if parallel)."""
    d1, d2 = p2 - p1, q2 - q1
    den = d1[0] * d2[1] - d1[1] * d2[0]
    if abs(den) < 1e-15 * np.linalg.norm(d1) * np.linalg.norm(d2):
        return np.array([np.inf, np.inf])
    s = ((q1[0] - p1[0]) * d2[1] - (q1[1] - p1[1]) * d2[0]) / den
    return p1 + s * d1


def instant_centres(th2: float, th3: float, th4: float, links) -> dict:
    """The six instant centres of a four-bar (1 = ground, 2 = crank, 3 = coupler, 4 = rocker). The four
    pins are I12 = O2, I23 = A, I34 = B, I14 = O4. Kennedy's theorem (the three instant centres of any
    three bodies lie on one line) places I13 on lines I12-I23 and I14-I34, and I24 on lines I12-I14 and
    I23-I34."""
    raise NotImplementedError


def velocity_ratio_ic(th2: float, th3: float, th4: float, links) -> tuple[float, float]:
    """Velocity ratios from instant centres alone. I24 is a point of both 2 and 4 with one velocity:
    w2 (I24 - I12) = w4 (I24 - I14), so w4/w2 = x24 / (x24 - r1) on the ground line. A = I23 moves as
    part of 2 about I12 and as part of 3 about I13: w3/w2 = (A - I12).(A - I13) / |A - I13|^2.
    Returns (w3/w2, w4/w2)."""
    raise NotImplementedError


def transmission_angle(th3: float, th4: float) -> float:
    """Acute angle between the coupler and the rocker, in [0, pi/2]. Near 0 the coupler pushes almost
    along the rocker and the linkage is close to toggle (it locks or binds); designers keep it above
    about 40 degrees as a rule of thumb."""
    raise NotImplementedError


def central_rates(t: np.ndarray, th: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """First and second derivatives of an angle history by central differences at the interior points.
    Returns (t[1:-1], dth/dt, d2th/dt2)."""
    h = t[1] - t[0]
    return t[1:-1], (th[2:] - th[:-2]) / (2 * h), (th[2:] - 2 * th[1:-1] + th[:-2]) / h ** 2


def observed_order(errors, ratio: float = 2.0) -> np.ndarray:
    """Observed convergence orders log(e_k / e_{k+1}) / log(ratio) on steps refined by ``ratio``."""
    e = np.asarray(errors, dtype=float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save four_bar.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    L = SERVO_LINKS
    mm = lambda v: v * 1e3
    print(f"servo linkage: ground {mm(L[0]):.0f} mm, servo arm {mm(L[1]):.0f} mm, pushrod {mm(L[2]):.0f} mm, "
          f"horn {mm(L[3]):.0f} mm: {grashof(*L)}")

    # neutral: servo arm vertical
    th2n = np.pi / 2
    th3n, th4n = position_circle_intersection(th2n, L, branch=1)
    a, b, it = solve_position_newton(th2n, L, (0.0, np.pi / 2))
    print(f"neutral (servo arm at 90 deg): pushrod {np.degrees(th3n):.4f} deg, horn {np.degrees(th4n):.4f} deg "
          f"(Newton from (0, 90 deg): {np.degrees(b):.4f} deg in {it} iterations)")
    # Newton's quadratic convergence
    x = np.array([0.2, 1.2])
    errs = []
    for _ in range(6):
        f = loop_residual(th2n, x[0], x[1], L)
        errs.append(np.linalg.norm(x - [th3n, th4n]))
        x = x - np.linalg.solve(loop_jacobian(x[0], x[1], L), f)
    print("Newton error per iteration from (11.5 deg, 68.8 deg): " + ", ".join(f"{e:.1e}" for e in errs))

    # servo travel +/- 45 deg
    travel = np.radians(np.linspace(-45.0, 45.0, 181))
    th2 = th2n + travel
    th3, th4 = sweep_positions(th2, L, (th3n, th4n))
    defl = th4 - th4n
    w42 = np.array([velocity_analysis(t2, t3, t4, 1.0, L)[1] for t2, t3, t4 in zip(th2, th3, th4)])
    mu = np.array([transmission_angle(t3, t4) for t3, t4 in zip(th3, th4)])
    print(f"\nservo travel -45 .. +45 deg: horn moves {np.degrees(defl[0]):.2f} .. {np.degrees(defl[-1]):.2f} deg "
          f"(asymmetry {np.degrees(defl[-1] + defl[0]):.2f} deg)")
    print(f"  gearing d(th4)/d(th2): neutral {w42[90]:.4f} (r2/r4 = {L[1] / L[3]:.4f}), "
          f"range {w42.min():.4f} .. {w42.max():.4f}")
    print(f"  transmission angle: neutral {np.degrees(mu[90]):.2f} deg, minimum {np.degrees(mu.min()):.2f} deg "
          f"at servo {np.degrees(travel[np.argmin(mu)]):+.0f} deg")
    lin = np.polyfit(travel, defl, 1)
    print(f"  best straight line {lin[0]:.4f} deg/deg; worst departure {np.degrees(np.max(np.abs(np.polyval(lin, travel) - defl))):.2f} deg")
    Mh = 0.30   # N m, illustrative hinge moment
    print(f"  hinge moment {Mh} N m needs servo torque {Mh * w42[90]:.3f} N m at neutral, "
          f"{Mh * w42[0]:.3f} and {Mh * w42[-1]:.3f} N m at -45 and +45 deg (virtual work: T2 w2 = M4 w4)")

    # instant centres agree with the velocity analysis
    k = 150
    w3, w4 = velocity_analysis(th2[k], th3[k], th4[k], 1.0, L)
    r32, r42 = velocity_ratio_ic(th2[k], th3[k], th4[k], L)
    ic = instant_centres(th2[k], th3[k], th4[k], L)
    print(f"\nat servo {np.degrees(travel[k]):+.0f} deg: w3/w2 = {w3:.6f}, w4/w2 = {w4:.6f} (loop closure); "
          f"{r32:.6f}, {r42:.6f} (instant centres)")
    print(f"  I13 = ({mm(ic['I13'][0]):.1f}, {mm(ic['I13'][1]):.1f}) mm, I24 = ({mm(ic['I24'][0]):.1f}, {mm(ic['I24'][1]):.1f}) mm")

    # independent check: finite differences of the position solution
    w2, al2 = 2.0, 3.0           # rad/s, rad/s^2 (a servo slewing and speeding up)
    t0 = 0.05
    th2f = lambda t: th2n - 0.5 + w2 * t + 0.5 * al2 * t ** 2
    t3a, t4a = sweep_positions(np.array([th2f(t0)]), L, (th3n, th4n))
    W3, W4 = velocity_analysis(th2f(t0), t3a[0], t4a[0], w2 + al2 * t0, L)
    A3, A4 = acceleration_analysis(th2f(t0), t3a[0], t4a[0], w2 + al2 * t0, al2, W3, W4, L)
    print(f"\nfinite differences of Newton positions, th2 = {np.degrees(th2f(t0)):.2f} deg, "
          f"w2 = {w2 + al2 * t0:.3f} rad/s, al2 = {al2} rad/s^2:")
    print(f"  analytic: w4 = {W4:.8f} rad/s, al4 = {A4:.8f} rad/s^2")
    print("      h (s)     |w4 FD - w4|   |al4 FD - al4|")
    ew, ea = [], []
    hs = [0.02, 0.01, 0.005, 0.0025]
    for h in hs:
        tt = np.array([t0 - h, t0, t0 + h])
        _, t4s = sweep_positions(th2f(tt), L, (t3a[0], t4a[0]))
        _, d1, d2 = central_rates(tt, t4s)
        ew.append(abs(d1[0] - W4))
        ea.append(abs(d2[0] - A4))
        print(f"   {h:8.4f}     {ew[-1]:.3e}      {ea[-1]:.3e}")
    print(f"  observed orders: velocity {np.array2string(observed_order(ew), precision=3)}, "
          f"acceleration {np.array2string(observed_order(ea), precision=3)}")

    # a horn-length sweep
    print("\nhorn length sweep (servo arm 20 mm, +/- 45 deg travel):")
    print("   horn (mm)  type            deflection range (deg)   asymmetry (deg)   min transmission (deg)")
    for r4 in (0.015, 0.020, 0.025, 0.030, 0.040):
        Ls = (L[0], L[1], L[2], r4)
        t3n, t4n = position_circle_intersection(th2n, Ls, branch=1)
        a3, a4 = sweep_positions(th2, Ls, (t3n, t4n))
        d = a4 - t4n
        m_ = min(transmission_angle(p, q) for p, q in zip(a3, a4))
        print(f"   {mm(r4):8.0f}  {grashof(*Ls):14s}  {np.degrees(d[0]):7.2f} .. {np.degrees(d[-1]):6.2f}"
              f"          {np.degrees(d[-1] + d[0]):6.2f}            {np.degrees(m_):6.2f}")

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    for kk, c in ((0, "C0"), (90, "C1"), (180, "C2")):
        P = joint_positions(th2[kk], th3[kk], th4[kk], L)
        pts = np.array([P["O2"], P["A"], P["B"], P["O4"]]) * 1e3
        ax[0].plot(pts[:, 0], pts[:, 1], "o-", color=c, label=f"servo {np.degrees(travel[kk]):+.0f} deg")
    ax[0].set_aspect("equal")
    ax[0].set_xlabel("x (mm)")
    ax[0].set_ylabel("y (mm)")
    ax[0].set_title("servo arm, pushrod, horn")
    ax[0].legend(fontsize=8)
    ax[1].plot(np.degrees(travel), np.degrees(defl), label="surface deflection")
    ax[1].plot(np.degrees(travel), np.degrees(np.polyval(lin, travel)), "k--", lw=0.8, label="best straight line")
    ax1b = ax[1].twinx()
    ax1b.plot(np.degrees(travel), np.degrees(mu), "C3", label="transmission angle")
    ax1b.set_ylabel("transmission angle (deg)")
    ax[1].set_xlabel("servo travel (deg)")
    ax[1].set_ylabel("horn rotation (deg)")
    ax[1].set_title("gearing and transmission")
    ax[1].legend(fontsize=8, loc="upper left")
    ax[2].loglog(hs, ew, "o-", label="w4")
    ax[2].loglog(hs, ea, "s-", label="alpha4")
    ax[2].loglog(hs, ea[0] * (np.array(hs) / hs[0]) ** 2, "k:", lw=0.8, label="slope 2")
    ax[2].set_xlabel("time step h (s)")
    ax[2].set_ylabel("|finite difference - analytic|")
    ax[2].set_title("independent check")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "four_bar.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'four_bar.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
