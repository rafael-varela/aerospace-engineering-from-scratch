"""P10.20 Ground Effect & Interference: starter code.

How much does the ground help a wing in the landing flare, and how much does one lifting surface
disturb another? Both are interference problems, and the vortex lattice method (VLM) of P10.11
handles both with one idea: the velocity a vortex system induces at a point is the Biot-Savart
sum over its filaments, so anything that changes the vortex system near a surface changes that
surface's loads.

* Ground effect: a flat ground plane is a streamline. The method of images enforces it exactly by
  adding a mirror-image wing at the same distance below the plane with circulation of opposite
  sign: on the plane, the vertical velocities of each vortex and its image cancel. Solving the VLM
  with the images included gives the lift and the span load in ground effect; the Trefftz-plane
  wake (with its image) gives the induced drag.
* Wing-tail interference: the wing's horseshoe system induces a downwash epsilon at the tail;
  d(epsilon)/d(alpha) from the lattice is compared with the far-field estimate 2 CL_alpha/(pi AR).

Independent checks: h -> infinity recovers the free-air results; the span load is symmetric; the
image system leaves zero normal velocity on the ground plane (a residual); a single horseshoe's
Trefftz downwash reduction matches its closed form; the VLM induced-drag ratio is compared with
the classical empirical formula attributed to Wieselsberger.

Axes: x downstream, y to starboard, z up; ground plane z = 0; the wing's lattice at height h.
Angles in radians inside the code. The wing is a flat lifting surface (camber is ignored; angles
are measured from the zero-lift line), with nc x ns panels, bound vortices at panel quarter
chords and control points at three-quarter chords. Trailing legs run to +infinity along +x.

Units: SI.

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes ground_effect.png)
Test: pytest curriculum/phase-10/lesson-20-ground-effect-interference
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# Reference trainer (Phase 10 shared conventions): rectangular wing, b = 10 m, c = 1.2 m.
SPAN = 10.0
CHORD = 1.2
H_FLARE = 1.0          # m, assumed height of the wing above the runway in the flare (a low-wing layout)
TAIL_ARM = 4.5         # m behind the wing quarter chord
TAIL_AREA = 2.0        # m^2
TAIL_AR = 4.0          # assumed tail aspect ratio (span 2.83 m)
EPS_CORE = 1e-12          # a point closer to a filament line than 1e-6 x the segment length (1e-6 m for a trailing leg) gets nothing from it


# ---------------------------------------------------------------- given: Biot-Savart for horseshoes (from P10.11)

def segment_velocity(P: np.ndarray, A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Velocity induced at points P (..., 3) by straight vortex segments A -> B (..., 3) of unit
    circulation (right-hand rule along A -> B). Broadcasts over leading axes. Points on the line
    of a segment get zero (no core model needed: control points never lie on filaments)."""
    r1 = P - A
    r2 = P - B
    r0 = B - A
    c = np.cross(r1, r2)
    c2 = np.sum(c * c, axis=-1)
    n1 = np.linalg.norm(r1, axis=-1)
    n2 = np.linalg.norm(r2, axis=-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        k = np.sum(r0 * (r1 / n1[..., None] - r2 / n2[..., None]), axis=-1) / (4 * np.pi * c2)
    k = np.where(c2 > EPS_CORE * np.sum(r0 * r0, axis=-1) ** 2, k, 0.0)
    return c * k[..., None]


def semi_infinite_velocity(P: np.ndarray, A: np.ndarray, d: np.ndarray = np.array([1.0, 0.0, 0.0])) -> np.ndarray:
    """Velocity induced at P by a unit vortex starting at A and running to infinity along unit d."""
    r1 = P - A
    c = np.cross(d, r1)
    c2 = np.sum(c * c, axis=-1)
    n1 = np.linalg.norm(r1, axis=-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        k = (1 + np.sum(r1 * d, axis=-1) / n1) / (4 * np.pi * c2)
    k = np.where(c2 > EPS_CORE, k, 0.0)
    return c * k[..., None]


def horseshoe_velocity(P: np.ndarray, A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Velocity at points P (m, 3) from unit horseshoes with bound legs A -> B (n, 3) and trailing
    legs to +infinity along +x: returns (m, n, 3). A on the port side, B to starboard, so a
    positive circulation lifts (P10.11)."""
    P = P[:, None, :]
    A = A[None, :, :]
    B = B[None, :, :]
    return segment_velocity(P, A, B) + semi_infinite_velocity(P, B) - semi_infinite_velocity(P, A)


@dataclass
class Lattice:
    """A planar lattice: bound-vortex end points A, B, control points, panel areas and the spanwise
    edges of each strip."""
    A: np.ndarray
    B: np.ndarray
    cp: np.ndarray
    area: np.ndarray
    y_edges: np.ndarray
    nc: int
    ns: int
    S: float
    b: float


def rectangular_lattice(b: float, c: float, ns: int, nc: int = 4, z: float = 0.0,
                        spacing: str = "cosine") -> Lattice:
    """Given: an unswept rectangular wing of span b and chord c at height z, leading edge at x = 0,
    split into nc chordwise and ns spanwise panels (cosine spacing in span refines the tips).
    Panels are ordered strip by strip, chordwise fastest."""
    if spacing == "cosine":
        y = -0.5 * b * np.cos(np.linspace(0, np.pi, ns + 1))
    else:
        y = np.linspace(-b / 2, b / 2, ns + 1)
    xs = np.linspace(0, c, nc + 1)
    A, B, cp, area = [], [], [], []
    for j in range(ns):
        for i in range(nc):
            dx = xs[i + 1] - xs[i]
            xq = xs[i] + 0.25 * dx
            A.append([xq, y[j], z]); B.append([xq, y[j + 1], z])
            cp.append([xs[i] + 0.75 * dx, 0.5 * (y[j] + y[j + 1]), z])
            area.append(dx * (y[j + 1] - y[j]))
    return Lattice(np.array(A), np.array(B), np.array(cp), np.array(area), y, nc, ns, b * c, b)


# ---------------------------------------------------------------- the learner's functions

def mirror(points: np.ndarray, z_ground: float = 0.0) -> np.ndarray:
    """Reflect points in the horizontal plane z = z_ground."""
    raise NotImplementedError


def influence_matrix(lat: Lattice, ground: bool = False) -> np.ndarray:
    """Normal (z) velocity at each control point due to unit circulation on each horseshoe,
    including, when ground is True, the mirror-image horseshoe (reflected in z = 0) with circulation
    of opposite sign."""
    raise NotImplementedError


def solve_lattice(lat: Lattice, alpha: float, ground: bool = False, V: float = 1.0) -> np.ndarray:
    """Circulations satisfying flow tangency on the flat lattice at angle of attack alpha (small
    angle: the free stream's normal component is V sin(alpha))."""
    raise NotImplementedError


def strip_circulation(lat: Lattice, gamma: np.ndarray) -> np.ndarray:
    """Total circulation of each spanwise strip (sum over its chordwise panels): the span load Gamma(y)."""
    raise NotImplementedError


def lift_coefficient(lat: Lattice, gamma: np.ndarray, V: float = 1.0) -> float:
    """CL = 2 sum(Gamma_strip dy) / (V S) (Kutta-Joukowski on the bound vortices, small angles)."""
    raise NotImplementedError


def trefftz_downwash(y_edges: np.ndarray, G: np.ndarray, y_eval: np.ndarray, h: float | None = None) -> np.ndarray:
    """Downwash (positive down) in the Trefftz plane at span stations y_eval on the wake line, from
    the 2-D trailing vortices shed at y_edges and, when h is given, their images 2h below with
    opposite sign. The vortex at edge k points along +x with strength s_k = G_(k-1) - G_k (strip
    circulations G, zero beyond the tips); a +x vortex at y_k induces w = s_k (y - y_k)/(2 pi r^2)."""
    raise NotImplementedError


def induced_drag_coefficient(lat: Lattice, gamma: np.ndarray, h: float | None = None, V: float = 1.0) -> float:
    """Trefftz-plane induced drag, Di = (rho/2) int Gamma w dy with w the far-wake downwash (including
    the image wake if h is given) at the strip centres, so CDi = sum(Gamma w dy) / (V^2 S)."""
    raise NotImplementedError


def wing_in_ground_effect(b: float, c: float, h: float | None, alpha: float, ns: int = 40, nc: int = 4) -> dict:
    """Solve the rectangular wing at height h (None = free air) and return CL, CDi, the span
    efficiency e = CL^2/(pi AR CDi), the span load and the lattice."""
    raise NotImplementedError


def ground_plane_residual(lat: Lattice, gamma: np.ndarray, pts_xy: np.ndarray) -> float:
    """Largest |w| on the ground plane z = 0 at the given (x, y) points induced by the wing plus its
    image (it must vanish), relative to the largest |w| the wing alone induces there."""
    raise NotImplementedError


def single_horseshoe_trefftz_ratio(h_over_bv: np.ndarray) -> np.ndarray:
    """Closed form: ratio of the Trefftz-plane downwash at the centre of a single horseshoe of
    vortex span b_v in ground effect to that in free air, 1 - 1/(1 + (4h/b_v)^2) =
    (4h/b_v)^2 / (1 + (4h/b_v)^2) (the image pair, 2h below, cancels part of the downwash)."""
    raise NotImplementedError


def wieselsberger_ratio(h_over_b: np.ndarray) -> np.ndarray:
    """Classical empirical estimate of the induced-drag ratio in ground effect,
    CDi(h)/CDi(inf) = 1 - sigma with sigma = (1 - 1.32 h/b)/(1.05 + 7.4 h/b), usually quoted for
    roughly 0.033 < h/b < 0.25 (attributed to Wieselsberger, 1922)."""
    raise NotImplementedError


def tail_downwash(lat: Lattice, gamma: np.ndarray, x_tail: float, z_tail: float, b_tail: float,
                  ground: bool = False, n: int = 20, V: float = 1.0) -> float:
    """Average downwash angle epsilon [rad] (positive down) across a tail span b_tail at (x_tail, z_tail)
    induced by the wing's horseshoe system (and its image if ground is True)."""
    raise NotImplementedError


def downwash_gradient_far_field(CL_alpha: float, AR: float) -> float:
    """Far-field estimate d(epsilon)/d(alpha) = 2 CL_alpha / (pi AR) (elliptic wing, tail far behind:
    twice the downwash at the wing)."""
    raise NotImplementedError


# ---------------------------------------------------------------- given helpers

def observed_order(errors: np.ndarray, ratio: float = 2.0) -> np.ndarray:
    """Observed order p = log(e_i / e_{i+1}) / log(ratio)."""
    e = np.asarray(errors, dtype=float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


def main(out_dir: str | Path = ".", quick: bool = False) -> None:
    """Print the lesson's numbers and write ground_effect.png. quick=True (used by the tests) skips the
    80-strip lattices and thins the figure's downwash sweep."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    alpha = np.radians(5.0)

    free = wing_in_ground_effect(SPAN, CHORD, None, alpha)
    CLa_free = free["CL"] / alpha
    print(f"trainer wing, b = {SPAN} m, c = {CHORD} m, AR = {free['AR']:.3f}, flat VLM 4 x 40 panels, alpha = 5 deg from zero lift")
    print(f"free air: CL = {free['CL']:.4f} (CL_alpha = {CLa_free:.4f} /rad), CDi = {free['CDi']:.5f}, e = {free['e']:.4f}")

    # convergence of the free-air solution
    print("lattice convergence (free air, nc = 4): ns, CL, CDi")
    res = []
    for ns in (10, 20, 40) + (() if quick else (80,)):
        r = wing_in_ground_effect(SPAN, CHORD, None, alpha, ns=ns)
        res.append((r["CL"], r["CDi"]))
        print(f"  {ns:3d}  {r['CL']:.6f}  {r['CDi']:.7f}")
    res = np.array(res)
    print(f"  observed order CL: {np.round(observed_order(np.abs(np.diff(res[:, 0]))), 2)}"
          f"  CDi: {np.round(observed_order(np.abs(np.diff(res[:, 1]))), 2)}")

    print("\nground effect at alpha = 5 deg: h/b, CL/CL_free, CDi/CDi_free at fixed alpha,"
          " K = (CDi/CL^2)/(CDi/CL^2)_free, Wieselsberger, single-horseshoe (b_v = pi b/4)")
    rows = []
    for hb in [0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.5, 1.0, 50.0]:
        g = wing_in_ground_effect(SPAN, CHORD, hb * SPAN, alpha)
        K = (g["CDi"] / g["CL"] ** 2) / (free["CDi"] / free["CL"] ** 2)
        rows.append((hb, g["CL"] / free["CL"], g["CDi"] / free["CDi"], K))
        wr = f"{float(wieselsberger_ratio(hb)):.4f}" if 0.033 <= hb <= 0.25 else "  --  "
        print(f"  {hb:6.3f}  {g['CL']/free['CL']:.4f}  {g['CDi']/free['CDi']:.4f}  {K:.4f}  "
              f"{wr}  {float(single_horseshoe_trefftz_ratio(hb / (np.pi / 4))):.4f}")
    Ks = []
    for ns in (20, 40) + (() if quick else (80,)):
        f_ = wing_in_ground_effect(SPAN, CHORD, None, alpha, ns=ns)
        g_ = wing_in_ground_effect(SPAN, CHORD, H_FLARE, alpha, ns=ns)
        Ks.append((g_["CDi"] / g_["CL"] ** 2) / (f_["CDi"] / f_["CL"] ** 2))
    print(f"  K ratio at h/b = 0.1 for ns = 20, 40{'' if quick else ', 80'}: " + ", ".join(f"{x:.4f}" for x in Ks))

    # the reference trainer in the flare
    g = wing_in_ground_effect(SPAN, CHORD, H_FLARE, alpha)
    K = (g["CDi"] / g["CL"] ** 2) / (free["CDi"] / free["CL"] ** 2)
    print(f"\ntrainer in the flare, wing at h = {H_FLARE} m (h/b = {H_FLARE/SPAN}): CL_alpha up {100*(g['CL']/free['CL']-1):.1f} %,"
          f" induced-drag factor K ratio {K:.3f} ({100*(1-K):.1f} % less induced drag at the same CL)")
    m, V_app, rho = 750.0, 30.0, 1.225
    CL_app = m * 9.81 / (0.5 * rho * V_app ** 2 * SPAN * CHORD)
    Di_free = 0.5 * rho * V_app ** 2 * SPAN * CHORD * CL_app ** 2 * free["CDi"] / free["CL"] ** 2
    print(f"  at {V_app} m/s, CL = {CL_app:.3f}: induced drag {Di_free:.0f} N in free air, {K*Di_free:.0f} N at h = {H_FLARE} m;"
          f" angle of attack for that CL {np.degrees(CL_app/CLa_free):.2f} deg -> {np.degrees(CL_app/(g['CL']/alpha)):.2f} deg")
    print(f"  span efficiency e: free {free['e']:.4f}, effective in ground effect {g['CL']**2/(np.pi*free['AR']*g['CDi']):.4f}")
    if not quick:   # the 4 x 40 absolutes are not converged: Richardson on CL from ns = 20, 40, 80, CDi at ns = 80
        p_cl = float(observed_order(np.abs(np.diff(res[:, 0])))[-1])
        CL_R = res[2, 0] + (res[2, 0] - res[1, 0]) / (2 ** p_cl - 1)
        e_R = CL_R ** 2 / (np.pi * free["AR"] * res[2, 1])
        Di_R = 0.5 * rho * V_app ** 2 * SPAN * CHORD * CL_app ** 2 / (np.pi * free["AR"] * e_R)
        print(f"  Richardson-extrapolated free air (p = {p_cl:.2f}): CL_alpha = {CL_R/alpha:.3f} /rad, e = {e_R:.3f};"
              f" induced drag {Di_R:.0f} N in free air, {Ks[-1]*Di_R:.0f} N at h = {H_FLARE} m (K ratio at ns = 80)")

    # WIG craft: AR 3 wing, chord 4 m, span 12 m, at 1 m
    wig_free = wing_in_ground_effect(12.0, 4.0, None, alpha, ns=30, nc=6)
    for hh in [0.5, 1.0, 2.0]:
        wg = wing_in_ground_effect(12.0, 4.0, hh, alpha, ns=30, nc=6)
        print(f"WIG craft wing 12 x 4 m (AR 3) at h = {hh} m (h/c = {hh/4:.3f}, h/b = {hh/12:.3f}): CL x {wg['CL']/wig_free['CL']:.3f},"
              f" K x {(wg['CDi']/wg['CL']**2)/(wig_free['CDi']/wig_free['CL']**2):.3f},"
              f" L/Di at alpha = 5 deg: {wig_free['CL']/wig_free['CDi']:.1f} -> {wg['CL']/wg['CDi']:.1f}")

    # ground-plane residual
    rng = np.random.default_rng(20)
    pts = np.column_stack([rng.uniform(-5, 15, 200), rng.uniform(-8, 8, 200)])
    print(f"\nground-plane residual max|w|/max|w_wing| at 200 random points: {ground_plane_residual(g['lat'], g['gamma'], pts):.2e}")

    # wing-tail downwash
    b_t = np.sqrt(TAIL_AR * TAIL_AREA)
    x_t = 0.25 * CHORD + TAIL_ARM
    print(f"\nwing-tail interference: tail span {b_t:.2f} m at x = {x_t:.2f} m (4.5 m behind c/4)")
    lat = free["lat"]
    for z_t in [0.0, 0.5, 1.0]:
        eps = tail_downwash(lat, free["gamma"], x_t, z_t, b_t)
        print(f"  tail {z_t:.1f} m above the wing plane: d eps/d alpha = {eps/alpha:.4f}")
    print(f"  far-field 2 CL_alpha/(pi AR) = {downwash_gradient_far_field(CLa_free, free['AR']):.4f};"
          f" at the wing CL_alpha/(pi AR) = {CLa_free/(np.pi*free['AR']):.4f}")
    for xt in [2.0, 5.0, 10.0, 20.0, 100.0]:
        print(f"  tail at x = {xt:5.1f} m (z = 0): d eps/d alpha = {tail_downwash(lat, free['gamma'], xt, 0.0, b_t)/alpha:.4f}")
    eps_g = tail_downwash(g["lat"], g["gamma"], x_t, H_FLARE, b_t, ground=True)
    eps_f = tail_downwash(lat, free["gamma"], x_t, 0.0, b_t)
    print(f"  in the flare (h = {H_FLARE} m, tail at wing height): d eps/d alpha {eps_f/alpha:.4f} -> {eps_g/alpha:.4f};"
          f" at alpha = 5 deg the tail sees {np.degrees(eps_f - eps_g):.2f} deg less downwash")

    # figure
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    hbs = np.array([r[0] for r in rows[:-1]])
    ax[0].plot(hbs, [r[3] for r in rows[:-1]], "o-", label="VLM with images, K/K_free")
    hh = np.linspace(0.033, 0.25, 50)
    ax[0].plot(hh, wieselsberger_ratio(hh), "--", label="Wieselsberger (empirical)")
    hh2 = np.linspace(0.03, 1.0, 100)
    ax[0].plot(hh2, single_horseshoe_trefftz_ratio(hh2 / (np.pi / 4)), ":", label="single horseshoe, b_v = pi b/4")
    ax[0].set_xlabel("h / b"); ax[0].set_ylabel("induced drag / free-air value at the same CL")
    ax[0].set_title("ground effect on induced drag"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
    yc = 0.5 * (free["lat"].y_edges[1:] + free["lat"].y_edges[:-1])
    ax[1].plot(yc, free["G"], label="free air")
    for hb in (0.2, 0.1, 0.05):
        gg = wing_in_ground_effect(SPAN, CHORD, hb * SPAN, alpha)
        ax[1].plot(yc, gg["G"], label=f"h/b = {hb}")
    ax[1].set_xlabel("y [m]"); ax[1].set_ylabel("Gamma / V [m]"); ax[1].set_title("span load at alpha = 5 deg")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    xs = np.linspace(1.5, 40, 15 if quick else 60)
    ax[2].plot(xs, [tail_downwash(lat, free["gamma"], x, 0.0, b_t) / alpha for x in xs], label="VLM, tail in wing plane")
    ax[2].axhline(downwash_gradient_far_field(CLa_free, free["AR"]), color="k", ls="--", label="2 CL_alpha/(pi AR)")
    ax[2].axvline(x_t, color="0.5", ls=":", label="trainer tail")
    ax[2].set_xlabel("x behind the leading edge [m]"); ax[2].set_ylabel("d eps / d alpha")
    ax[2].set_title("downwash gradient behind the wing"); ax[2].legend(fontsize=8); ax[2].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "ground_effect.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
