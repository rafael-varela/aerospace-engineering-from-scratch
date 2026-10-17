"""P05.03 Trusses: Joints & Sections: starter code.

A truss is a set of straight members joined at their ends by pins, loaded only at the joints. Each
member is then a two-force member: it carries one axial force, tension or compression. This module
solves trusses three independent ways and checks them against each other:

* the method of joints written as one linear system: d equations per joint (d = 2 or 3), unknowns the
  member forces and the support reactions; counting and rank say whether the truss is determinate;
* the method of sections: the external reactions from equilibrium of the whole truss as one rigid
  body, then a cut through at most three members (2-D) and three equations for the free body on one side;
* the direct stiffness method (P03.13): member stiffness EA/L, nodal displacements, member forces from
  elongations. For a determinate truss it must give the same forces whatever EA is; for an indeterminate
  truss it is the only one of the three that works.

Aerospace cases: a Pratt girder of 4130 steel tube (the geometry of truss spars and welded-tube
fuselage sides), and a 3-D engine-mount truss from a firewall to an engine ring.

Units: SI (N, m, Pa). 2-D: x along the span, y up. 3-D mount: x forward along the thrust axis, z up.

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes truss.png)
Test: pytest curriculum/phase-05/lesson-03-trusses-joints-sections
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

G0 = 9.80665                      # m/s^2
# Typical properties of normalized 4130 steel tube (see the lesson for the source and caveats).
E_4130 = 200e9                    # Pa
FTY_4130 = 517e6                  # Pa, tensile yield, typical
FTU_4130 = 655e6                  # Pa, tensile ultimate, typical


@dataclass
class Truss:
    """Pin-jointed truss.

    nodes: (n, d) coordinates. members: (m, 2) node index pairs. fixed: (n, d) booleans, True where a
    support stops that displacement component (each True is one reaction unknown). loads: (n, d) applied
    joint forces (N)."""
    nodes: np.ndarray
    members: np.ndarray
    fixed: np.ndarray
    loads: np.ndarray = field(default=None)

    def __post_init__(self) -> None:
        self.nodes = np.asarray(self.nodes, dtype=float)
        self.members = np.asarray(self.members, dtype=int)
        self.fixed = np.asarray(self.fixed, dtype=bool)
        self.loads = np.zeros_like(self.nodes) if self.loads is None else np.asarray(self.loads, dtype=float)

    @property
    def dim(self) -> int:
        return self.nodes.shape[1]


class TrussError(ValueError):
    """Raised when a truss is not statically determinate or a cut cannot be solved."""


# ---------------------------------------------------------------- geometry and the joint equations

def member_geometry(t: Truss) -> tuple[np.ndarray, np.ndarray]:
    """Lengths (m,) and unit vectors (m, d) pointing from the first node of each member to the second."""
    raise NotImplementedError


def joint_matrix(t: Truss) -> np.ndarray:
    """Equilibrium matrix of the method of joints, shape (n d, m + r).

    Row block i holds the d force equations of joint i. A tension T_k in member k = (i, j) pulls joint i
    toward j (+T_k u_k) and joint j toward i (-T_k u_k). Each fixed component adds a unit reaction column.
    Unknown order: the m member forces, then the r reactions in the order of np.argwhere(t.fixed)."""
    raise NotImplementedError


def determinacy(t: Truss, rtol: float = 1e-10) -> dict:
    """Counting (m + r against d n) and rank of the joint matrix. 'determinate' needs rank = d n = m + r;
    m + r > rank means redundant members or reactions (indeterminate); rank < d n means a mechanism."""
    raise NotImplementedError


def solve_joints(t: Truss) -> tuple[np.ndarray, np.ndarray]:
    """Method of joints as one linear system B x = -f. Returns (member forces (m,), tension positive;
    reactions (n, d), zero where a component is free). Raises TrussError unless determinate."""
    raise NotImplementedError


def joint_residual(t: Truss, forces: np.ndarray, reactions: np.ndarray) -> np.ndarray:
    """Net force on every joint, shape (n, d): loads + reactions + member pulls. Zero in equilibrium."""
    raise NotImplementedError


def zero_force_members(forces: np.ndarray, scale: float, rtol: float = 1e-10) -> np.ndarray:
    """Indices of members whose force is zero to round-off relative to ``scale`` (a typical force)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the method of sections (2-D)

def external_reactions_2d(t: Truss) -> np.ndarray:
    """Support reactions of a planar truss from the three equilibrium equations of the WHOLE truss as a
    rigid body (sum Fx = sum Fy = sum M_origin = 0), without any member force. Needs exactly three
    reaction components. Returns (n, 2) reactions."""
    raise NotImplementedError


def solve_section(t: Truss, cut: list[int], side: list[int], reactions: np.ndarray) -> np.ndarray:
    """Forces in the members ``cut`` by a section, from equilibrium of the free body made of the joints in
    ``side``. Each cut member pulls its ``side`` joint toward the joint across the cut (tension
    positive). The free body carries its joints' loads and reactions; members wholly inside it are
    internal and cancel. Three equations (sum Fx, sum Fy, sum M about the origin); up to three unknowns.
    Raises TrussError if a cut member does not cross the cut or the cut members are concurrent or
    parallel."""
    raise NotImplementedError


# ---------------------------------------------------------------- the direct stiffness method

def stiffness_solve(t: Truss, EA: np.ndarray | float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Direct stiffness method. Member k has stiffness EA_k / L_k along u_k; the global matrix
    K = sum (EA/L) [[uu^T, -uu^T], [-uu^T, uu^T]] is assembled, the fixed components are removed,
    K_ff q_f = f_f is solved, and each member force is EA/L times its elongation u . (q_j - q_i).
    Returns (displacements (n, d), member forces (m,), reactions (n, d)). Works for determinate and
    indeterminate trusses; fails (singular K) for a mechanism."""
    raise NotImplementedError


# ---------------------------------------------------------------- truss generators

def pratt_truss(n_panels: int = 6, span: float = 4.8, depth: float = 0.6, P: float = 2000.0,
                cross_diagonals: bool = False) -> Truss:
    """Planar Pratt girder with vertical end posts. Bottom joints 0..n at y = 0, top joints n+1..2n+1 at
    y = depth. Members: bottom chord, top chord, n + 1 verticals and one diagonal per panel sloping down
    toward mid-span (tension under downward loads). Pin at joint 0, roller (vertical) at joint n. A load
    -P (down) at every interior bottom joint. ``cross_diagonals`` adds the second diagonal in every panel,
    which makes the truss indeterminate of degree n_panels."""
    n = n_panels
    p = span / n
    bottom = [(i * p, 0.0) for i in range(n + 1)]
    top = [(i * p, depth) for i in range(n + 1)]
    nodes = np.array(bottom + top)
    T = lambda i: n + 1 + i
    members = [(i, i + 1) for i in range(n)] + [(T(i), T(i + 1)) for i in range(n)] + [(i, T(i)) for i in range(n + 1)]
    for i in range(n):
        if i < n / 2:
            members.append((T(i), i + 1))      # left half: from top outer to bottom inner
        else:
            members.append((i, T(i + 1)))      # right half: mirror image
        if cross_diagonals:
            members.append((i, T(i + 1)) if i < n / 2 else (T(i), i + 1))
    fixed = np.zeros_like(nodes, dtype=bool)
    fixed[0] = [True, True]
    fixed[n] = [False, True]
    loads = np.zeros_like(nodes)
    loads[1:n, 1] = -P
    return Truss(nodes, np.array(members), fixed, loads)


def pratt_members(n_panels: int) -> dict[str, list[int]]:
    """Member index groups of ``pratt_truss`` (without cross diagonals): bottom, top, verticals, diagonals."""
    n = n_panels
    return {"bottom": list(range(n)), "top": list(range(n, 2 * n)), "verticals": list(range(2 * n, 3 * n + 1)),
            "diagonals": list(range(3 * n + 1, 4 * n + 1))}


def pratt_section(n_panels: int, panel: int) -> tuple[list[int], list[int]]:
    """Cut through panel ``panel`` (0-based) of ``pratt_truss``: the bottom chord, top chord and diagonal
    of that panel, and the joints left of the cut. Returns (cut members, side joints)."""
    g = pratt_members(n_panels)
    cut = [g["bottom"][panel], g["top"][panel], g["diagonals"][panel]]
    side = list(range(panel + 1)) + [n_panels + 1 + i for i in range(panel + 1)]
    return cut, side


def beam_moment_shear(t: Truss, x: float) -> tuple[float, float]:
    """Bending moment M(x) and shear V(x) of the truss treated as a simply supported beam: the reactions
    and loads to the left of x. M is positive sagging (compression on top)."""
    raise NotImplementedError


def pratt_beam_analogy(t: Truss, n_panels: int, panel: int) -> tuple[float, float, float]:
    """(bottom chord, top chord, diagonal) forces of one Pratt panel from the beam analogy alone.

    Moments about the joint where the other two cut members meet: the bottom chord carries M / h at the
    top joint of the diagonal, the top chord -M / h at its bottom joint; the diagonal's vertical
    component carries the panel shear V. Exact for this truss (it is the method of sections written with
    the beam's M and V)."""
    raise NotImplementedError


def engine_mount_truss(mass: float = 120.0, nz: float = 4.0, thrust: float = 3000.0, torque: float = -450.0,
                       a: float = 0.25, b: float = 0.18, length: float = 0.55, cg_ahead: float = 0.30) -> Truss:
    """3-D engine-mount truss: four firewall joints (x = 0, corners (+-a, +-a) in y-z, fully fixed) and
    four engine-ring joints (x = length, corners (+-b, +-b)). Members: the four ring sides and two struts
    from each ring joint to the firewall corner behind it and to the next corner anticlockwise seen from
    ahead (12 members: determinate). Loads: the engine (weight times nz down, thrust forward, propeller
    torque about x) acts at its CG, ``cg_ahead`` in front of the ring, and reaches the four ring joints
    as the smallest nodal forces statically equivalent to it (``rigid_body_nodal_loads``)."""
    ang = np.radians([45.0, 135.0, 225.0, 315.0])
    fw = np.column_stack([np.zeros(4), a * np.sqrt(2) * np.cos(ang), a * np.sqrt(2) * np.sin(ang)])
    ring = np.column_stack([np.full(4, length), b * np.sqrt(2) * np.cos(ang), b * np.sqrt(2) * np.sin(ang)])
    nodes = np.vstack([fw, ring])
    members = [(4 + i, 4 + (i + 1) % 4) for i in range(4)]
    members += [(4 + i, i) for i in range(4)] + [(4 + i, (i + 1) % 4) for i in range(4)]
    fixed = np.zeros((8, 3), dtype=bool)
    fixed[:4] = True
    cg = np.array([length + cg_ahead, 0.0, 0.0])
    F = np.array([thrust, 0.0, -nz * mass * G0])
    M = np.array([torque, 0.0, 0.0])
    loads = np.zeros((8, 3))
    loads[4:] = rigid_body_nodal_loads(ring, cg, F, M)
    return Truss(nodes, np.array(members), fixed, loads)


def rigid_body_nodal_loads(points: np.ndarray, at: np.ndarray, F: np.ndarray, M: np.ndarray) -> np.ndarray:
    """Smallest (least-squares norm) set of nodal forces f_i at ``points`` whose resultant is F and whose
    moment about ``at`` is M: a statically equivalent load transfer through a rigid body. Shape (k, 3)."""
    k = len(points)
    A = np.zeros((6, 3 * k))
    for i, p in enumerate(points):
        r = p - at
        A[:3, 3 * i:3 * i + 3] = np.eye(3)
        A[3:, 3 * i:3 * i + 3] = np.array([[0, -r[2], r[1]], [r[2], 0, -r[0]], [-r[1], r[0], 0]])
    f = np.linalg.lstsq(A, np.concatenate([F, M]), rcond=None)[0]
    return f.reshape(k, 3)


# ---------------------------------------------------------------- sizing helpers

def tube_area(od: float, wall: float) -> float:
    """Cross-section area of a round tube, pi/4 (D^2 - (D - 2t)^2) (m^2)."""
    raise NotImplementedError


def tube_second_moment(od: float, wall: float) -> float:
    """Second moment of area of a round tube, pi/64 (D^4 - (D - 2t)^4) (m^4); derived in P05.05."""
    raise NotImplementedError


def euler_load_pinned(E: float, I: float, L: float) -> float:
    """Euler buckling load of a pin-ended strut, pi^2 E I / L^2 (N). A preview of P05.15: a truss
    member in compression must be checked against this as well as against yield."""
    raise NotImplementedError


def volume_index(t: Truss, forces: np.ndarray) -> float:
    """Sum of |F_k| L_k (N m). Dividing by an allowable stress gives the material volume of a truss whose
    every member is sized to that stress: a figure of merit for the geometry."""
    raise NotImplementedError


def depth_sweep(span: float, n_panels: int, P: float, depths: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Largest chord force and the volume index of the Pratt girder for each depth."""
    raise NotImplementedError


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save truss.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)

    # ---- Pratt girder
    n, span, h, P = 6, 4.8, 0.6, 2000.0
    t = pratt_truss(n, span, h, P)
    info = determinacy(t)
    print(f"Pratt girder: {n} panels, span {span} m, depth {h} m, {P:.0f} N at each of {n - 1} bottom joints")
    print(f"  m = {info['m']}, r = {info['r']}, 2j = {info['equations']}, rank {info['rank']}: {info['status']}")
    f, R = solve_joints(t)
    g = pratt_members(n)
    print(f"  reactions: joint 0 ({R[0, 0]:.1f}, {R[0, 1]:.1f}) N, joint {n} ({R[n, 0]:.1f}, {R[n, 1]:.1f}) N")
    for name in ("bottom", "top", "verticals", "diagonals"):
        print(f"  {name:9s}: " + " ".join(f"{f[k]:8.0f}" for k in g[name]))
    print(f"  joint residual max {np.abs(joint_residual(t, f, R)).max():.1e} N")
    zf = zero_force_members(f, P)
    print(f"  zero-force members: {[int(k) for k in zf]}")
    Rext = external_reactions_2d(t)
    print("  panel: (bottom, top, diagonal) by joints | by sections | by the beam analogy (N)")
    for panel in range(n):
        cut, side = pratt_section(n, panel)
        fs = solve_section(t, cut, side, Rext)
        ba = pratt_beam_analogy(t, n, panel)
        print(f"    {panel}: {np.array2string(f[cut], precision=1, suppress_small=True)} | "
              f"{np.array2string(fs, precision=1, suppress_small=True)} | "
              f"{np.array2string(np.array(ba), precision=1, suppress_small=True)}")
    # stiffness method: same forces whatever EA is
    od, wall = 0.75 * 0.0254, 0.035 * 0.0254
    A = tube_area(od, wall)
    _, f_k, R_k = stiffness_solve(t, E_4130 * A)
    rng = np.random.default_rng(503)
    _, f_k2, _ = stiffness_solve(t, E_4130 * A * rng.uniform(0.2, 5.0, len(t.members)))
    print(f"  stiffness method, uniform EA: max |difference from joints| {np.abs(f_k - f).max():.1e} N; "
          f"random EA (0.2x to 5x): {np.abs(f_k2 - f).max():.1e} N")
    Q, _, _ = stiffness_solve(t, E_4130 * A)
    print(f"  4130 tube {od * 1e3:.2f} x {wall * 1e3:.3f} mm: A = {A * 1e6:.2f} mm^2; largest tension "
          f"{f.max():.0f} N -> {f.max() / A / 1e6:.1f} MPa, largest compression {f.min():.0f} N -> "
          f"{f.min() / A / 1e6:.1f} MPa (typical Fty {FTY_4130 / 1e6:.0f} MPa); mid-span deflection {-Q[n // 2, 1] * 1e3:.3f} mm")
    Pcr = euler_load_pinned(E_4130, tube_second_moment(od, wall), span / n)
    print(f"  Euler load of a {span / n:.1f} m pin-ended bay of that tube: {Pcr:.0f} N (preview of P05.15): the top "
          f"chord's {-f.min():.0f} N is {-f.min() / Pcr:.2f} times that; the tube buckles at {Pcr / A / 1e6:.0f} MPa, "
          f"far below yield")
    tc = pratt_truss(n, span, h, P, cross_diagonals=True)
    ic = determinacy(tc)
    print(f"  with crossed diagonals: {ic['status']}, degree {ic['indeterminate_degree']}")
    _, fc1, _ = stiffness_solve(tc, E_4130 * A)
    EAc = np.full(len(tc.members), E_4130 * A); EAc[-1] *= 10
    _, fc2, _ = stiffness_solve(tc, EAc)
    print(f"    last diagonal pair, equal EA: {fc1[-2]:.0f}, {fc1[-1]:.0f} N; one diagonal 10x stiffer: "
          f"{fc2[-2]:.0f}, {fc2[-1]:.0f} N (load sharing depends on stiffness)")

    # ---- depth sweep
    ratios = np.array([0.05, 0.08, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5])
    fmax, vol = depth_sweep(span, n, P, ratios * span)
    print("\ndepth/span, largest chord force (N), volume index sum|F|L (N m):")
    for r_, fm, v in zip(ratios, fmax, vol):
        print(f"  {r_:5.3f}  {fm:9.0f}  {v:9.0f}")
    fine = np.linspace(0.05, 0.6, 221)
    _, vf = depth_sweep(span, n, P, fine * span)
    print(f"  minimum volume index at depth/span = {fine[np.argmin(vf)]:.3f}")

    # ---- 3-D engine mount
    em = engine_mount_truss()
    ie = determinacy(em)
    fe, Re = solve_joints(em)
    print(f"\nengine mount: m = {ie['m']}, r = {ie['r']}, 3j = {ie['equations']}, rank {ie['rank']}: {ie['status']}")
    print("  member forces (N): " + " ".join(f"{v:.0f}" for v in fe))
    print(f"  joint residual max {np.abs(joint_residual(em, fe, Re)).max():.1e} N")
    tot = em.loads.sum(axis=0) + Re.sum(axis=0)
    print(f"  sum of loads + reactions: {tot.round(9)} N")
    _, fe_k, _ = stiffness_solve(em, E_4130 * A)
    print(f"  stiffness method max |difference| {np.abs(fe_k - fe).max():.1e} N")
    print(f"  largest compression {fe.min():.0f} N = {fe.min() / A / 1e6:.0f} MPa in the 19.05 mm tube")

    # ---- figure
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.5))
    ax = axs[0]
    fm_abs = np.abs(f).max()
    for k, (i, j) in enumerate(t.members):
        c = "tab:blue" if f[k] > 1e-9 * P else ("tab:red" if f[k] < -1e-9 * P else "lightgray")
        ax.plot(*t.nodes[[i, j]].T, color=c, lw=1 + 4 * abs(f[k]) / fm_abs)
        ax.text(*t.nodes[[i, j]].mean(axis=0), f"{f[k] / 1e3:.1f}", fontsize=7, ha="center")
    ax.set_aspect("equal"); ax.set_title("Pratt girder: kN (blue tension, red compression)"); ax.set_xlabel("x (m)")
    ax = axs[1]
    ax.plot(ratios, fmax / 1e3, "o-", label="largest chord force (kN)")
    ax.plot(fine, vf / 1e3, "-", label="volume index (kN m)")
    ax.set_xlabel("depth / span"); ax.set_yscale("log"); ax.legend(); ax.set_title("depth sweep")
    fig.tight_layout()
    fig.savefig(out_dir / "truss.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
