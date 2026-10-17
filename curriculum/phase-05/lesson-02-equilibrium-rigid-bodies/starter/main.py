"""P05.02 Equilibrium of Rigid Bodies: starter code.

A rigid body is in equilibrium when the resultant force and the resultant moment of everything acting
on it vanish: six scalar equations in 3-D, three in a plane. The supports supply the unknown reactions.
This module turns that statement into a reaction solver:

* support models as lists of unknown scalars (a force along a known direction at a known point, or a
  couple about a known axis): roller, pin, fixed support and two-force link in 2-D; ball-and-socket,
  hinge, link and fixed support in 3-D;
* the equilibrium matrix A, one column per unknown ([u; (r - O) x u] for a force, [0; u] for a couple),
  and the load vector b, so that equilibrium is A x + b = 0;
* a classification by counting and by rank: determinate (rank = equations = unknowns), statically
  indeterminate (more unknowns than the rank: degree n - rank), partially constrained (fewer unknowns than
  equations) and improperly constrained (enough unknowns but rank-deficient: a mechanism);
* three aerospace bodies: a strut-braced wing in front view (2-D), an aileron on two hinges and a pushrod
  (3-D, checked by one moment equation about the hinge line), and an engine on a six-link mount (3-D,
  checked by virtual work with link-length changes measured from the geometry).

Units: SI (N, m, N m). 3-D bodies: x forward, y to starboard (spanwise), z up. The 2-D wing is drawn in
front view with x outboard along the span and y up.

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes equilibrium.png)
Test: pytest curriculum/phase-05/lesson-02-equilibrium-rigid-bodies
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

G0 = 9.80665                       # m/s^2


class SupportError(ValueError):
    """Raised when a support set cannot give a unique equilibrium solution for the applied load."""


# ---------------------------------------------------------------- unknowns and supports

@dataclass
class Unknown:
    """One unknown reaction scalar: a force of signed magnitude x along ``direction`` through ``point``
    (kind 'force'), or a couple of signed magnitude x about ``direction`` (kind 'couple')."""
    name: str
    point: np.ndarray
    direction: np.ndarray
    kind: str = "force"

    def __post_init__(self) -> None:
        p = np.asarray(self.point, dtype=float)
        d = np.asarray(self.direction, dtype=float)
        self.point = np.append(p, 0.0) if p.size == 2 else p
        d = np.append(d, 0.0) if d.size == 2 else d
        self.direction = d / np.linalg.norm(d)


def roller(name: str, point, normal) -> list[Unknown]:
    """2-D roller or rocker: one force along the surface normal."""
    return [Unknown(name, point, normal)]


def pin(name: str, point) -> list[Unknown]:
    """2-D pin: force components along x and y, no moment."""
    return [Unknown(name + "_x", point, [1.0, 0.0]), Unknown(name + "_y", point, [0.0, 1.0])]


def fixed_2d(name: str, point) -> list[Unknown]:
    """2-D fixed (built-in) support: two forces and a couple about z."""
    return pin(name, point) + [Unknown(name + "_M", point, [0.0, 0.0, 1.0], "couple")]


def link(name: str, body_point, frame_point) -> list[Unknown]:
    """Two-force member (strut, rod, cable) pinned at both ends and unloaded between them: one force along
    the line from the body point to the frame point. Positive = tension (pulls the body toward the frame)."""
    b, f = np.asarray(body_point, float), np.asarray(frame_point, float)
    return [Unknown(name, b, f - b)]


def ball_socket(name: str, point) -> list[Unknown]:
    """3-D ball-and-socket joint: three force components, no moment."""
    return [Unknown(f"{name}_{c}", point, e) for c, e in zip("xyz", np.eye(3))]


def hinge(name: str, point, axis, axial: bool = True, moments: bool = False) -> list[Unknown]:
    """3-D hinge (journal bearing) with its pin along ``axis``: two forces normal to the axis, plus an
    axial force if ``axial`` (a thrust bearing) and two couples normal to the axis if ``moments`` (a long
    bearing that resists tilting; for a body on two well-aligned hinges these are taken as zero)."""
    a = np.asarray(axis, float) / np.linalg.norm(axis)
    k = int(np.argmax(np.abs(a)))
    if abs(abs(a[k]) - 1.0) < 1e-12:                          # axis along x, y or z: name by coordinate
        (i, j) = [q for q in range(3) if q != k]
        t1, t2, l1, l2 = np.eye(3)[i], np.eye(3)[j], "xyz"[i], "xyz"[j]
    else:
        t = np.cross(a, [1.0, 0.0, 0.0]) if abs(a[0]) < 0.9 else np.cross(a, [0.0, 1.0, 0.0])
        t1 = t / np.linalg.norm(t)
        t2, l1, l2 = np.cross(a, t1), "n1", "n2"
    out = [Unknown(f"{name}_{l1}", point, t1), Unknown(f"{name}_{l2}", point, t2)]
    if axial:
        out.append(Unknown(name + "_a", point, a))
    if moments:
        out += [Unknown(f"{name}_M{l1}", point, t1, "couple"), Unknown(f"{name}_M{l2}", point, t2, "couple")]
    return out


def fixed_3d(name: str, point) -> list[Unknown]:
    """3-D fixed support: three forces and three couples."""
    return ball_socket(name, point) + [Unknown(name + "_M" + c, point, e, "couple") for c, e in zip("xyz", np.eye(3))]


# ---------------------------------------------------------------- the equilibrium equations

ROWS_2D = [0, 1, 5]                # Fx, Fy, Mz out of (Fx, Fy, Fz, Mx, My, Mz)


def equilibrium_matrix(unknowns: list[Unknown], dim: int = 3, ref=None) -> np.ndarray:
    """Equilibrium matrix: one column per unknown, [u; (r - ref) x u] for a force and [0; u] for a couple.
    Shape (6, n) for dim = 3, (3, n) with rows (Fx, Fy, Mz) for dim = 2."""
    raise NotImplementedError


def load_vector(points, forces, couples=None, dim: int = 3, ref=None) -> np.ndarray:
    """Resultant of the applied loads about ref, (sum F; sum (r - ref) x F + sum C): shape (6,) or, for
    dim = 2, (3,) = (Fx, Fy, Mz). Points and forces may be 2-D or 3-D."""
    raise NotImplementedError


def classify_supports(A: np.ndarray, rtol: float = 1e-9) -> dict:
    """Count and rank. Returns n_eq, n_unknowns, rank, degree of indeterminacy (n - rank: how many
    reactions equilibrium alone cannot find), degree of freedom (n_eq - rank: independent rigid motions
    the supports do not stop), the condition number of A, and a status:
    'determinate', 'indeterminate', 'partially constrained' or 'improperly constrained'.
    The rank uses singular values below rtol * s_max as zero."""
    raise NotImplementedError


def solve_reactions(unknowns: list[Unknown], points, forces, couples=None, dim: int = 3, ref=None) -> dict[str, float]:
    """Reactions that hold the body in equilibrium under the applied loads, as {name: value}.

    Solves A x = -b when the supports are determinate. Raises SupportError otherwise: an indeterminate
    body needs compatibility (deformation) to share the load (P05.16), and a partially or improperly
    constrained body is a mechanism for a general load."""
    raise NotImplementedError


def residual(unknowns: list[Unknown], values: dict[str, float], points, forces, couples=None, dim: int = 3,
             ref=None) -> np.ndarray:
    """A x + b about ref: the force and moment left over. Zero (to round-off) in equilibrium."""
    raise NotImplementedError


# ---------------------------------------------------------------- 1. strut-braced wing (2-D front view)

@dataclass
class StrutWing:
    """Front view of one wing half (x outboard along the span, y up), root pinned to the fuselage at the
    origin, braced by a strut from the fuselage point ``strut_frame`` to the wing at span ``y_strut``.
    lift_half: lift on this semispan (N, elliptic, acting at 4 s / 3 pi); wing_weight: weight of the
    half-wing (N, uniform, acting at s / 2); n: load factor applied to the wing's own weight (inertia)."""
    semispan: float = 5.0
    y_strut: float = 2.6
    strut_frame: tuple = (0.0, -1.0)
    lift_half: float = 0.0
    wing_weight: float = 0.0
    n: float = 1.0

    def unknowns(self) -> list[Unknown]:
        return pin("root", (0.0, 0.0)) + link("strut", (self.y_strut, 0.0), self.strut_frame)

    def loads(self) -> tuple[np.ndarray, np.ndarray]:
        ybar = 4.0 * self.semispan / (3.0 * np.pi)
        pts = np.array([[ybar, 0.0], [0.5 * self.semispan, 0.0]])
        F = np.array([[0.0, self.lift_half], [0.0, -self.n * self.wing_weight]])
        return pts, F


def strut_wing_reactions(w: StrutWing) -> dict[str, float]:
    """Root pin reactions and strut force (tension positive) by the general solver."""
    raise NotImplementedError


def strut_force_closed_form(w: StrutWing) -> float:
    """Strut force from moments about the root pin alone (the independent hand route):
    T * y_s * h / l_s = L ybar - n W s / 2, with h the depth of the strut's fuselage end below the root
    and l_s the strut length: the strut's lever arm about the root is y_s h / l_s."""
    raise NotImplementedError


# ---------------------------------------------------------------- 2. aileron on two hinges and a pushrod

@dataclass
class Aileron:
    """Aileron hinged along the y axis (x forward, z up). Inboard hinge at y = y1: ball-and-socket.
    Outboard hinge at y = y2: hinge free to slide along the axis (no axial force). Pushrod pinned to a horn
    at (0, y_horn, -horn) and running forward along +x to the airframe: a link. Load: the aerodynamic
    resultant (force N up at x_cp aft of the hinge line, spanwise station y_cp) and a chordwise force Fx."""
    y1: float = 0.10
    y2: float = 1.40
    y_horn: float = 0.75
    horn: float = 0.06
    N: float = 900.0
    x_cp: float = -0.08
    y_cp: float = 0.70
    Fx: float = -40.0

    def unknowns(self) -> list[Unknown]:
        hp = np.array([0.0, self.y_horn, -self.horn])
        return (ball_socket("hinge1", (0.0, self.y1, 0.0)) + hinge("hinge2", (0.0, self.y2, 0.0), (0, 1, 0), axial=False)
                + link("pushrod", hp, hp + np.array([0.4, 0.0, 0.0])))

    def loads(self) -> tuple[np.ndarray, np.ndarray]:
        return np.array([[self.x_cp, self.y_cp, 0.0]]), np.array([[self.Fx, 0.0, self.N]])


def aileron_reactions(a: Aileron) -> dict[str, float]:
    raise NotImplementedError


def hinge_moment(a: Aileron) -> float:
    """Moment of the aerodynamic load about the hinge line (the y axis), N m."""
    raise NotImplementedError


def pushrod_force_by_axis_moment(a: Aileron) -> float:
    """Pushrod force from one equation: moments about the hinge line. Both hinge reactions cut the line
    and drop out, so P * (moment of a unit pushrod force about y) = -(hinge moment)."""
    raise NotImplementedError


# ---------------------------------------------------------------- 3. engine on a six-link mount

def hexapod_mount(r_engine: float = 0.25, r_outer: float = 0.34, r_inner: float = 0.16, length_x: float = 0.40,
                  spread_deg: float = 40.0, extra_link: bool = False,
                  drop_link: int | None = None) -> tuple[list[Unknown], np.ndarray, np.ndarray]:
    """Six two-force links from three engine points to the firewall (the plane x = -length_x).

    Engine points sit at radius r_engine in the plane x = 0 at 90, 210 and 330 deg around the thrust axis
    (x). Each is joined to two firewall points: one at radius r_outer and angle theta + spread, one at
    radius r_inner and angle theta - spread. With spread 0 every link lies in a plane that contains the
    thrust axis, so no link has a moment about it: the mount cannot resist propeller torque.
    ``extra_link`` adds a seventh (redundant) link, ``drop_link`` removes one. Returns (unknowns,
    engine points, frame points), the points as arrays of shape (n_links, 3)."""
    th = np.radians([90.0, 210.0, 330.0])
    sp = np.radians(spread_deg)
    eng, frm = [], []
    for t in th:
        e = np.array([0.0, r_engine * np.cos(t), r_engine * np.sin(t)])
        for rad, s in ((r_outer, sp), (r_inner, -sp)):
            eng.append(e)
            frm.append(np.array([-length_x, rad * np.cos(t + s), rad * np.sin(t + s)]))
    if extra_link:
        eng.append(eng[0].copy())
        frm.append(np.array([-length_x, 0.0, 0.0]))
    if drop_link is not None:
        del eng[drop_link], frm[drop_link]
    eng, frm = np.array(eng), np.array(frm)
    unknowns = [Unknown(f"link{i + 1}", e, f - e) for i, (e, f) in enumerate(zip(eng, frm))]
    return unknowns, eng, frm


def engine_loads(mass: float, nz: float, thrust: float, torque: float, cg=(0.25, 0.0, 0.0)) -> tuple:
    """Loads on the engine: weight times load factor at the CG (down, -z), thrust along +x at the CG and
    the propeller-torque reaction as a couple about x. Returns (points, forces, couples) for 3-D."""
    pts = np.array([cg, cg], dtype=float)
    F = np.array([[0.0, 0.0, -nz * mass * G0], [thrust, 0.0, 0.0]])
    C = np.array([[torque, 0.0, 0.0]])
    return pts, F, C


def link_lengths(eng: np.ndarray, frm: np.ndarray, t: np.ndarray, theta: np.ndarray, ref=np.zeros(3)) -> np.ndarray:
    """Lengths of the links after the engine moves rigidly by translation t and rotation vector theta
    about ref (exact rotation by Rodrigues' formula), with the firewall fixed."""
    raise NotImplementedError


def virtual_work_residual(eng: np.ndarray, frm: np.ndarray, link_forces: np.ndarray, points, forces, couples,
                          t: np.ndarray, theta: np.ndarray, h: float = 1e-6) -> float:
    """Virtual work of every force on the engine for a small rigid motion (h t, h theta), divided by h.

    A tension T_i in link i does work -T_i * (change in its length): lengthening a link against its
    tension costs work. The link-length changes are measured from the moved geometry by a central
    difference, not from the equilibrium matrix, so this is an independent check: it vanishes for every
    rigid motion exactly when the body is in equilibrium (principle of virtual work)."""
    raise NotImplementedError


def spread_sweep(spreads_deg: np.ndarray, loads: tuple) -> tuple[np.ndarray, np.ndarray]:
    """Largest link force magnitude and condition number of A for a range of link spread angles."""
    raise NotImplementedError


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save equilibrium.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)

    # ---- 1. strut-braced wing
    m_ac, m_wing_half, n = 1100.0, 60.0, 3.8
    L_half = 0.5 * n * m_ac * G0
    w = StrutWing(lift_half=L_half, wing_weight=m_wing_half * G0, n=n)
    r = strut_wing_reactions(w)
    print(f"strut-braced wing, m = {m_ac:.0f} kg, n = {n}: lift per half {L_half:.0f} N at "
          f"{4 * w.semispan / (3 * np.pi):.4f} m, half-wing weight x n = {n * m_wing_half * G0:.0f} N at {w.semispan / 2:.2f} m")
    print(f"  solver: root_x = {r['root_x']:.0f} N, root_y = {r['root_y']:.0f} N, strut = {r['strut']:.0f} N "
          f"({'tension' if r['strut'] > 0 else 'compression'})")
    print(f"  moments about the root: strut = {strut_force_closed_form(w):.0f} N")
    wneg = StrutWing(lift_half=-0.5 * 1.52 * m_ac * G0, wing_weight=m_wing_half * G0, n=-1.52)
    print(f"  at n = -1.52 the strut carries {strut_wing_reactions(wneg)['strut']:.0f} N (compression: a buckling case, P05.15)")
    ys = np.linspace(1.0, 4.6, 10)
    print("  strut station y_s (m) vs strut force and root resultant (kN):")
    sweep_rows = []
    for y in ys:
        ww = StrutWing(y_strut=y, lift_half=L_half, wing_weight=m_wing_half * G0, n=n)
        rr = strut_wing_reactions(ww)
        sweep_rows.append((y, rr["strut"], np.hypot(rr["root_x"], rr["root_y"])))
        print(f"    {y:.1f}: {rr['strut'] / 1e3:7.2f}  {np.hypot(rr['root_x'], rr['root_y']) / 1e3:7.2f}")

    # ---- 2. aileron
    a = Aileron()
    ra = aileron_reactions(a)
    print(f"\naileron: N = {a.N:.0f} N at x = {a.x_cp} m, y = {a.y_cp} m; hinge moment {hinge_moment(a):.2f} N m")
    for k, v in ra.items():
        print(f"  {k:10s} {v:9.2f} N")
    print(f"  pushrod by moments about the hinge line: {pushrod_force_by_axis_moment(a):.2f} N")
    print(f"  classification: {classify_supports(equilibrium_matrix(a.unknowns()))['status']}")

    # ---- 3. engine mount
    U, eng, frm = hexapod_mount()
    loads = engine_loads(mass=120.0, nz=4.0, thrust=3000.0, torque=-450.0)
    info = classify_supports(equilibrium_matrix(U))
    re = solve_reactions(U, *loads)
    print(f"\nengine mount, six links: {info['status']}, rank {info['rank']}, cond {info['cond']:.2f}")
    for k, v in re.items():
        print(f"  {k}: {v:9.1f} N")
    x = np.array(list(re.values()))
    rng = np.random.default_rng(502)
    vw = [abs(virtual_work_residual(eng, frm, x, loads[0], loads[1], loads[2], rng.normal(size=3), rng.normal(size=3)))
          for _ in range(5)]
    print(f"  virtual-work residual over 5 random rigid motions: max {max(vw):.2e} N m (load scale {120 * 4 * G0:.0f} N)")
    for label, kw in [("seven links", dict(extra_link=True)), ("five links", dict(drop_link=5)),
                      ("spread 0 deg (all links cut the thrust axis)", dict(spread_deg=0.0))]:
        Ui, _, _ = hexapod_mount(**kw)
        ci = classify_supports(equilibrium_matrix(Ui))
        print(f"  {label}: {ci['status']}, rank {ci['rank']}, indeterminate degree {ci['indeterminate_degree']}, "
              f"freedom {ci['freedom_degree']}")
    spreads = np.array([40.0, 20.0, 10.0, 5.0, 2.0, 1.0])
    fmax, cond = spread_sweep(spreads, loads)
    print("  spread (deg), largest link force (N), cond(A):")
    for s, f, c in zip(spreads, fmax, cond):
        print(f"    {s:5.1f}  {f:10.0f}  {c:10.1f}")

    # ---- figure
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.4))
    ax = axs[0]
    sr = np.array(sweep_rows)
    ax.plot(sr[:, 0], sr[:, 1] / 1e3, "o-", label="strut force")
    ax.plot(sr[:, 0], sr[:, 2] / 1e3, "s-", label="root pin resultant")
    ax.set_xlabel("strut station y_s (m)"); ax.set_ylabel("kN"); ax.set_title("strut-braced wing, n = 3.8"); ax.legend()
    ax = axs[1]
    names = list(ra)
    ax.bar(range(len(names)), [ra[k] for k in names])
    ax.set_xticks(range(len(names)), names, rotation=45, fontsize=8); ax.set_ylabel("N"); ax.set_title("aileron reactions")
    ax = axs[2]
    sw = np.linspace(1.0, 40.0, 40)
    fm, cn = spread_sweep(sw, loads)
    ax.semilogy(sw, fm, label="largest link force (N)")
    ax.semilogy(sw, cn, label="cond(A)")
    ax.set_xlabel("link spread (deg)"); ax.set_title("engine mount near an improper constraint"); ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "equilibrium.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
