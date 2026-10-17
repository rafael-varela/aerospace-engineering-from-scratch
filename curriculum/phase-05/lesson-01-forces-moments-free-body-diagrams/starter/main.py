"""P05.01 Forces, Moments & Free-Body Diagrams: starter code.

A main landing-gear leg is the first structure in this phase. On touchdown the ground pushes up on the
tyre and, while the wheel spins up, drags it aft. The leg is held by a trunnion pin at its top and by a
drag brace, a strut that runs from the leg forward to the airframe. Which loads do the pin and the brace
carry? To answer, this module builds the tools a free-body diagram needs:

* forces as vectors, the moment of a force about a point (M = r x F) and about an axis (u . (r x F));
* the resultant of a force system, moved from one point to another (M_B = M_A + (A - B) x R), and its
  central axis, the line along which the system reduces to a force and a parallel couple (a wrench);
* the resultant of a distributed load and the point where it acts (uniform, triangular, elliptic lift);
* a 2-D free-body-diagram solver: known forces and couples plus unknown forces of known direction,
  assembled into the three equations sum Fx = 0, sum Fy = 0, sum Mz = 0 and solved;
* the landing gear: the trunnion reactions and the drag-brace force, by that solver and, independently,
  by one moment equation about the trunnion pin.

Units: SI (N, m, N m). Gear coordinates: x forward, y up, z to starboard (right-handed).

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes landing_gear_fbd.png)
Test: pytest curriculum/phase-05/lesson-01-forces-moments-free-body-diagrams
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

G0 = 9.80665                       # m/s^2, standard gravity


# ---------------------------------------------------------------- vectors and moments

def unit(v: np.ndarray) -> np.ndarray:
    """Unit vector along v. Raises ValueError for a zero vector."""
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if n == 0.0:
        raise ValueError("zero vector has no direction")
    return v / n


def force_vector(magnitude: float, direction: np.ndarray) -> np.ndarray:
    """Force of the given magnitude along ``direction`` (any non-zero vector, 2-D or 3-D)."""
    return magnitude * unit(direction)


def to3(v: np.ndarray) -> np.ndarray:
    """Pad 2-D vectors (shape (2,) or (n, 2)) with a zero z component; 3-D input is returned as float."""
    v = np.asarray(v, dtype=float)
    if v.shape[-1] == 3:
        return v
    if v.shape[-1] != 2:
        raise ValueError("vectors must have 2 or 3 components")
    pad = np.zeros(v.shape[:-1] + (1,))
    return np.concatenate([v, pad], axis=-1)


def moment_about_point(r: np.ndarray, F: np.ndarray, P: np.ndarray | None = None) -> np.ndarray:
    """Moment about P of a force F acting at r: M_P = (r - P) x F.

    r and F have shape (3,) or (n, 3) (2-D vectors are padded with z = 0); P defaults to the origin.
    The result has the shape of the broadcast inputs. Units: m x N = N m."""
    raise NotImplementedError


def moment_2d(r: np.ndarray, F: np.ndarray, P: np.ndarray | None = None) -> float | np.ndarray:
    """z component of the moment about P of a planar force F at r: (rx - Px) Fy - (ry - Py) Fx.
    Positive anticlockwise seen from +z."""
    raise NotImplementedError


def moment_about_axis(r: np.ndarray, F: np.ndarray, P: np.ndarray, u: np.ndarray) -> float:
    """Moment of F (acting at r) about the axis through P along u: u_hat . ((r - P) x F).

    It is the part of the moment that tries to turn a body about that axis, and it does not depend on
    which point P of the axis is used."""
    raise NotImplementedError


def resultant(points: np.ndarray, forces: np.ndarray, couples: np.ndarray | None = None,
              P: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Equivalent force-couple system at P of forces F_i at r_i plus free couples C_j.

    Returns (R, M_P) with R = sum F_i and M_P = sum (r_i - P) x F_i + sum C_j, both shape (3,).
    ``points`` and ``forces`` have shape (n, 2) or (n, 3); ``couples`` shape (m, 3) or (3,)."""
    raise NotImplementedError


def transport_moment(M_A: np.ndarray, R: np.ndarray, A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Moment about B of a system whose resultant R acts with moment M_A about A:
    M_B = M_A + (A - B) x R. The resultant itself does not change when the point moves."""
    raise NotImplementedError


def central_axis(R: np.ndarray, M_O: np.ndarray) -> tuple[np.ndarray, float]:
    """Central axis of a force system with resultant R and moment M_O about the origin.

    Returns (p, pitch): p = R x M_O / |R|^2 is the point of the axis closest to the origin, and the pitch
    is R . M_O / |R|^2 (metres). About any point of the axis the moment is pitch * R, parallel to R (a
    wrench). A planar system has R . M_O = 0, pitch 0: it reduces to a single force along that line."""
    raise NotImplementedError


# ---------------------------------------------------------------- distributed loads

def distributed_load_resultant(x: np.ndarray, w: np.ndarray) -> tuple[float, float]:
    """Resultant W = integral w dx and its point of action x_bar = integral x w dx / W of a line load
    w(x) (N/m) sampled at x, by the trapezoidal rule (error O(h^2) for a smooth w)."""
    raise NotImplementedError


def elliptic_lift(y: np.ndarray, semispan: float, lift_semispan: float) -> np.ndarray:
    """Elliptic spanwise lift per unit span on one wing half, w(y) = w0 sqrt(1 - (y/s)^2), scaled so that
    its integral over 0 <= y <= s is ``lift_semispan``: w0 = 4 L / (pi s). Its centroid is at 4 s / (3 pi)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the 2-D free-body diagram

def equilibrium_residual_2d(points: np.ndarray, forces: np.ndarray, couples: float = 0.0,
                            ref: np.ndarray | None = None) -> np.ndarray:
    """(sum Fx, sum Fy, sum Mz about ref) for planar forces F_i at r_i plus a free couple (N, N, N m).
    A body is in equilibrium when all three vanish."""
    raise NotImplementedError


def solve_fbd_2d(known_points: np.ndarray, known_forces: np.ndarray, unknown_points: np.ndarray,
                 unknown_dirs: np.ndarray, known_couple: float = 0.0, ref: np.ndarray | None = None,
                 cond_max: float = 1e10) -> np.ndarray:
    """Magnitudes of three unknown planar forces of known direction that put a body in equilibrium.

    Unknown k acts at ``unknown_points[k]`` along the unit vector of ``unknown_dirs[k]`` with signed
    magnitude x_k (a negative value means it acts the other way). Each unknown contributes a column
    [ux, uy, (r - ref) x u] to a 3 x 3 matrix A; the known loads give b = (sum Fx, sum Fy, sum Mz); the
    answer solves A x = -b. Raises ValueError when A is singular (cond > cond_max), which happens when the
    three unknown lines of action meet in one point or are all parallel: the supports cannot resist
    every load and the body is improperly constrained (lesson P05.02)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the landing gear

@dataclass
class GearGeometry:
    """Side view (x forward, y up) of a main-gear leg with a trunnion pin and a drag brace.

    contact: tyre contact point C, where the ground force acts.
    trunnion: pin A at the top of the leg (reacts a force in any direction, no moment).
    brace_leg: point B where the drag brace is pinned to the leg.
    brace_frame: point E where the drag brace is pinned to the airframe.
    The brace is pinned at both ends and carries no load between them, so it is a two-force member:
    its force on the leg acts along the line B-E (tension pulls B toward E)."""
    contact: np.ndarray = field(default_factory=lambda: np.array([0.0, 0.0]))
    trunnion: np.ndarray = field(default_factory=lambda: np.array([0.0, 1.20]))
    brace_leg: np.ndarray = field(default_factory=lambda: np.array([0.0, 0.70]))
    brace_frame: np.ndarray = field(default_factory=lambda: np.array([0.60, 1.10]))

    def brace_dir(self) -> np.ndarray:
        """Unit vector from B to E: the direction of a tension brace force on the leg."""
        return unit(np.asarray(self.brace_frame, float) - np.asarray(self.brace_leg, float))


def ground_force(V: float, D: float) -> np.ndarray:
    """Ground force on the tyre for a vertical reaction V (up, N) and a drag load D (aft, N): (-D, V)."""
    return np.array([-D, V])


def gear_reactions(geom: GearGeometry, V: float, D: float) -> dict[str, float]:
    """Trunnion reactions (Ax, Ay) and drag-brace force T (tension positive) on the gear leg, by
    ``solve_fbd_2d`` with the ground force at the contact point. Weight of the leg neglected."""
    raise NotImplementedError


def brace_force_by_moments(geom: GearGeometry, V: float, D: float) -> float:
    """Drag-brace force from one equation, sum of moments about the trunnion A = 0. The two trunnion
    reactions pass through A and drop out, so T = -M_A(ground force) / M_A(unit brace force).
    This is the independent check on ``gear_reactions``."""
    raise NotImplementedError


def brace_lever_arm(geom: GearGeometry) -> float:
    """Perpendicular distance from the trunnion A to the brace line B-E (m)."""
    raise NotImplementedError


def gear_fbd_forces(geom: GearGeometry, V: float, D: float) -> tuple[np.ndarray, np.ndarray]:
    """All forces on the leg and their points, ground first, then Ax, Ay and the brace force: arrays of
    shape (4, 2), ready for ``equilibrium_residual_2d``."""
    raise NotImplementedError


def landing_load_case(mass: float, main_share: float = 0.9, n_gear: float = 3.0,
                      drag_ratio: float = 0.25) -> tuple[float, float]:
    """Illustrative level-landing loads on ONE main gear: V = n_gear * (main_share / 2) * m g0 and
    D = drag_ratio * V. Not a certification load case; real cases come from the applicable regulations."""
    raise NotImplementedError


def brace_sweep(geom: GearGeometry, V: float, D: float, heights: np.ndarray) -> np.ndarray:
    """Brace force for a range of heights of the brace attachment B on the leg (E fixed): one value per
    height, NaN where the brace line passes through the trunnion."""
    raise NotImplementedError


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save landing_gear_fbd.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    np.set_printoptions(precision=4, suppress=True)

    # ---- moments: a worked 3-D example
    r, F, P = np.array([0.0, -1.2, 0.30]), np.array([-1000.0, 3000.0, 500.0]), np.zeros(3)
    print("moment of F = (-1000, 3000, 500) N at r = (0, -1.2, 0.3) m about the origin:",
          moment_about_point(r, F, P), "N m")
    print("  about the x axis (a fore-aft trunnion pin):", f"{moment_about_axis(r, F, P, [1, 0, 0]):.1f} N m;",
          "same axis through (5, 0, 0):", f"{moment_about_axis(r, F, [5, 0, 0], [1, 0, 0]):.1f} N m")

    # ---- the landing gear
    geom = GearGeometry()
    mass = 1100.0
    V, D = landing_load_case(mass)
    print(f"\nlevel-landing case (illustrative): m = {mass:.0f} kg, 90 % on the mains, n = 3.0 on the gear, "
          f"D/V = 0.25 -> V = {V:.0f} N, D = {D:.0f} N per main gear")
    rx = gear_reactions(geom, V, D)
    print(f"brace direction {geom.brace_dir()}, lever arm about A {brace_lever_arm(geom):.4f} m")
    print(f"system solve:   Ax = {rx['Ax']:.1f} N, Ay = {rx['Ay']:.1f} N, brace T = {rx['brace']:.1f} N "
          f"({'tension' if rx['brace'] > 0 else 'compression'})")
    T1 = brace_force_by_moments(geom, V, D)
    print(f"moments about A: T = {T1:.1f} N  (difference {abs(T1 - rx['brace']):.1e} N)")
    pts, Fs = gear_fbd_forces(geom, V, D)
    for ref in (geom.trunnion, geom.contact, np.array([3.0, -2.0])):
        print(f"  residual about {ref}: {equilibrium_residual_2d(pts, Fs, 0.0, ref)}")
    rv = gear_reactions(geom, V, 0.0)
    print(f"drag removed (D = 0): Ax = {rv['Ax']:.1f}, Ay = {rv['Ay']:.1f}, brace = {rv['brace']:.1f} N")
    print(f"trunnion resultant {np.hypot(rx['Ax'], rx['Ay']):.0f} N at "
          f"{np.degrees(np.arctan2(rx['Ay'], rx['Ax'])):.1f} deg from +x")

    # ---- distributed load: elliptic lift on a semispan
    s, L_half = 5.0, 0.5 * mass * G0
    for n in (11, 41, 161, 641):
        y = np.linspace(0.0, s, n)
        W, yb = distributed_load_resultant(y, elliptic_lift(y, s, L_half))
        print(f"elliptic lift, {n:4d} points: W = {W:.2f} N (exact {L_half:.2f}), y_bar = {yb:.5f} m "
              f"(exact 4s/3pi = {4 * s / (3 * np.pi):.5f})")

    # ---- the experiment: brace attachment height and drag ratio
    heights = np.linspace(0.05, 1.15, 23)
    Tb = brace_sweep(geom, V, D, heights)
    print("\nbrace attachment height h (m) vs brace force (kN):")
    for h, t in zip(heights[::2], Tb[::2]):
        print(f"  h = {h:.2f}: T = {t / 1e3:8.2f}")
    ratios = np.linspace(0.0, 0.8, 9)
    print("drag ratio D/V vs brace force and trunnion resultant (kN):")
    for q in ratios:
        rr = gear_reactions(geom, V, q * V)
        print(f"  D/V = {q:.1f}: T = {rr['brace'] / 1e3:6.2f}, |A| = {np.hypot(rr['Ax'], rr['Ay']) / 1e3:6.2f}")

    # sensitivity: brace attachment points moved by +-5 mm
    rng = np.random.default_rng(501)
    T_pert = []
    for _ in range(2000):
        g = GearGeometry(geom.contact, geom.trunnion, geom.brace_leg + rng.uniform(-0.005, 0.005, 2),
                         geom.brace_frame + rng.uniform(-0.005, 0.005, 2))
        T_pert.append(gear_reactions(g, V, D)["brace"])
    T_pert = np.array(T_pert)
    print(f"brace force with +-5 mm on B and E (2000 seeded samples): {T_pert.min():.0f} .. {T_pert.max():.0f} N, "
          f"spread {100 * (T_pert.max() - T_pert.min()) / rx['brace']:.2f} % of nominal")

    # ---- figure
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.6))
    ax = axs[0]
    for (a, b) in [(geom.contact, geom.trunnion), (geom.brace_leg, geom.brace_frame)]:
        ax.plot([a[0], b[0]], [a[1], b[1]], "k-", lw=3)
    sc = 1.0 / 40e3
    for p, f, lab in zip(pts, Fs, ["ground", "Ax", "Ay", "brace"]):
        ax.annotate("", xy=p + f * sc, xytext=p, arrowprops=dict(arrowstyle="->", color="tab:red", lw=1.8))
        ax.text(*(p + f * sc), f" {lab}\n {np.linalg.norm(f) / 1e3:.1f} kN", fontsize=8, color="tab:red")
    ax.plot(*geom.trunnion, "ko"); ax.plot(*geom.brace_frame, "k^")
    ax.set_aspect("equal"); ax.set_xlim(-0.6, 1.0); ax.set_ylim(-0.2, 1.6)
    ax.set_xlabel("x forward (m)"); ax.set_ylabel("y up (m)"); ax.set_title("free-body diagram of the leg")
    ax = axs[1]
    ax.plot(heights, Tb / 1e3, "o-")
    ax.axvline(geom.brace_leg[1], color="gray", lw=0.8)
    ax.set_ylim(0, 60); ax.set_xlabel("brace attachment height h (m)"); ax.set_ylabel("brace force (kN)")
    ax.set_title("brace force vs attachment height")
    ax = axs[2]
    y = np.linspace(0.0, s, 161)
    ax.plot(y, elliptic_lift(y, s, L_half), label="elliptic w(y)")
    ax.axvline(4 * s / (3 * np.pi), color="tab:red", ls="--", label="resultant line")
    ax.set_xlabel("y (m)"); ax.set_ylabel("lift per span (N/m)"); ax.set_title("distributed load and resultant")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "landing_gear_fbd.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
