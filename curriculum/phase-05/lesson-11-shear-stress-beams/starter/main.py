"""P05.11 Shear Stress in Beams: starter code.

A wing spar carries the lift to the fuselage as a bending moment (lesson P05.10) and as a shear
force. The shear force is carried mostly by the thin web, and it travels along the walls of a
thin-walled section as a *shear flow* q (force per unit length of wall). This module computes:

* the shear stress tau = V Q / (I t) in solid sections (rectangle, circle) and in an I-section,
  the "VQ/It" formula of elementary beam theory;
* the shear flow in any open thin-walled section made of straight wall segments, branched or not,
  for shear forces along both axes and sections with a product of inertia:

      q(s) = q_start - (Sx Ixx - Sy Ixy)/D * int t x ds - (Sy Iyy - Sx Ixy)/D * int t y ds,
      D = Ixx Iyy - Ixy^2  (x, y from the centroid);

* the resultant force and moment of a shear-flow distribution, and from them the shear centre: the
  point through which the shear force must act for the section to bend without twisting;
* an independent check: the shear flow obtained directly from longitudinal equilibrium of a slice of
  wall, as the integral of the change of the bending stress between two stations a distance dz apart
  (P05.10's unsymmetric bending formula, integrated numerically), with no VQ/It formula;
* fastener pitch from the shear flow at a cap-web joint.

Sign convention: x horizontal, y up, z along the span. A segment is directed from its start point to
its end point and q > 0 means the flow runs in that direction. (Sx, Sy) is the shear force acting on
the section face whose outward normal is +z; with P05.10's convention (M positive compressing the
top), dM/dz = -Sy and dM_h/dz = -Sx. The resultant of q is (Sx, Sy).

Units: SI (N, m, Pa, N/m).

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes shear_flow.png)
Test: pytest curriculum/phase-05/lesson-11-shear-stress-beams
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

G0 = 9.80665                       # m/s^2
F_SU_7075 = 331e6                  # Pa, typical ultimate shear strength of 7075-T6 (ASM Aerospace Specification Metals data sheet)
F_SU_2024 = 283e6                  # Pa, typical ultimate shear strength of 2024-T3 (ASM Aerospace Specification Metals data sheet)

Segment = tuple[tuple[float, float], tuple[float, float], float]   # ((x1, y1), (x2, y2), t)


# ---------------------------------------------------------------- solid sections: VQ/It

def vq_it(V: float, Q: float, I: float, t: float) -> float:
    """Average shear stress across a cut of width t: tau = V Q / (I t), Q the first moment about the
    neutral axis of the area beyond the cut."""
    raise NotImplementedError


def rect_shear_stress(V: float, b: float, h: float, y: np.ndarray) -> np.ndarray:
    """Shear stress at height y (from the centroid) in a solid b x h rectangle:
    tau = 3V/(2bh) (1 - 4 y^2/h^2), zero at the top and bottom, 1.5 V/A at the neutral axis."""
    raise NotImplementedError


def circle_max_shear(V: float, r: float) -> float:
    """Shear stress on the neutral axis of a solid circle of radius r: 4V/(3A)."""
    raise NotImplementedError


def i_section_shear_stress(V: float, h: float, b_f: float, t_f: float, t_w: float,
                           y: np.ndarray) -> np.ndarray:
    """Vertical shear stress VQ/(I b(y)) at height y in a doubly symmetric I-section (overall depth h,
    flanges b_f x t_f, web t_w), with b(y) = b_f in a flange and t_w in the web. In the flanges this
    'vertical' stress is small and not physically meaningful (the flange shear is horizontal, see
    the thin-walled shear flow); the formula is used there only to show the jump at the junction."""
    raise NotImplementedError


# ---------------------------------------------------------------- thin-walled open sections

def thin_wall_properties(segments: list[Segment]) -> dict:
    """Area, centroid and centroidal second moments of a thin-walled section given as straight
    centreline segments ((x1, y1), (x2, y2), t), using the thin-wall approximation (each wall's own
    bending about its centreline, of order t^3, is neglected). Along a straight segment
    int t y^2 ds = t L (y1^2 + y1 y2 + y2^2) / 3, and similarly for x^2 and x y."""
    raise NotImplementedError


def _incoming(segments: list[Segment], k: int, q_end: list[float], tol: float = 1e-12) -> float:
    """Sum of the end shear flows of the already processed segments that end where segment k starts."""
    p = np.array(segments[k][0])
    total = 0.0
    for j in range(k):
        if np.linalg.norm(np.array(segments[j][1]) - p) < tol:
            total += q_end[j]
    return total


def open_section_shear_flow(segments: list[Segment], Sx: float, Sy: float,
                            n: int = 51) -> list[tuple[np.ndarray, np.ndarray]]:
    """Shear flow in an open thin-walled section (shear force through the shear centre, so no twist).

    ``segments`` are listed so that every segment starts either at a free edge (q = 0) or at a
    junction where all the segments flowing into it have already been listed; the flow entering a
    junction leaves through the next segment that starts there. For the result to satisfy
    equilibrium the last segment must end at a free edge, where q comes back to 0 (a check).

    Returns, for each segment, (s, q): n points of arc length from its start and the shear flow
    (N/m) there, positive in the segment's direction, from the closed-form integral of the linear
    x(s), y(s) along the segment."""
    raise NotImplementedError


def shear_flow_resultant(segments: list[Segment], flows: list[tuple[np.ndarray, np.ndarray]],
                         about: tuple[float, float] = (0.0, 0.0)) -> tuple[float, float, float]:
    """Force (Fx, Fy) and moment Mz (anticlockwise, about the point ``about``) of a shear-flow
    distribution, by Simpson's rule along each segment (exact for the quadratic q of a straight wall)."""
    raise NotImplementedError


def shear_centre(segments: list[Segment]) -> tuple[float, float]:
    """Shear centre (xs, ys) of an open thin-walled section: the point where a shear force must act
    to produce the same moment as its shear flow. With a unit Sy, xs follows from Mz = xs Sy about the
    origin; with a unit Sx, ys from Mz = -ys Sx."""
    raise NotImplementedError


def channel_segments(h: float, b: float, t_f: float, t_w: float) -> list[Segment]:
    """Centreline model of a channel (C) section: web of height h on x = 0, flanges of length b
    pointing to +x. Directed from the top flange tip, down the web, out to the bottom flange tip."""
    return [((b, h / 2), (0.0, h / 2), t_f),
            ((0.0, h / 2), (0.0, -h / 2), t_w),
            ((0.0, -h / 2), (b, -h / 2), t_f)]


def channel_shear_centre_closed_form(h: float, b: float, t_f: float, t_w: float) -> float:
    """Distance e of the channel's shear centre from the web centreline, on the side away from the
    flanges: e = 3 b^2 t_f / (6 b t_f + h t_w) (centreline dimensions, thin walls)."""
    raise NotImplementedError


def i_beam_segments(h: float, b: float, t_f: float, t_w: float) -> list[Segment]:
    """Centreline model of an I-section (flanges of width b at y = +-h/2, web on x = 0). Directed so
    every junction's inflow is known before its outflow: both top half-flanges into the top junction,
    down the web, the bottom-left half-flange into the bottom junction, out along the bottom-right."""
    return [((-b / 2, h / 2), (0.0, h / 2), t_f),
            ((b / 2, h / 2), (0.0, h / 2), t_f),
            ((0.0, h / 2), (0.0, -h / 2), t_w),
            ((-b / 2, -h / 2), (0.0, -h / 2), t_f),
            ((0.0, -h / 2), (b / 2, -h / 2), t_f)]


def z_segments(h: float, b: float, t: float) -> list[Segment]:
    """Centreline model of a Z-section: top flange from (b, h/2) to the web, web down x = 0, bottom
    flange to (-b, -h/2). Antisymmetric, so its centroid and shear centre are at the origin."""
    return [((b, h / 2), (0.0, h / 2), t),
            ((0.0, h / 2), (0.0, -h / 2), t),
            ((0.0, -h / 2), (-b, -h / 2), t)]


# ---------------------------------------------------------------- independent check

def shear_flow_from_bending_gradient(segments: list[Segment], Sx: float, Sy: float,
                                     dz: float = 0.01, n: int = 2001) -> list[tuple[np.ndarray, np.ndarray]]:
    """Shear flow from longitudinal equilibrium of a slice of wall, without the VQ/It formula.

    Two stations dz apart carry moments that differ by dM = -Sy dz and dM_h = -Sx dz. The bending
    stress (unsymmetric formula of P05.10) changes by d_sigma at each wall point; a strip of wall
    between a free edge and arc length s is in equilibrium along z only if
    q(s) dz = q_start dz - integral(d_sigma t ds). The integral is done numerically (cumulative
    trapezoid on n points per segment) and flows are carried through junctions as in
    ``open_section_shear_flow``."""
    raise NotImplementedError


# ---------------------------------------------------------------- fasteners

def fastener_pitch(q: float, P_allow: float, factor: float = 1.0) -> float:
    """Largest spacing of single fasteners along a joint carrying shear flow q (N/m), each fastener
    allowed P_allow (N): factor * q * pitch = P_allow."""
    raise NotImplementedError


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save shear_flow.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    # the P05.10 spar root, limit load
    V = 3.8 * 750.0 * G0 / 2.0
    h, b, tf, tw = 0.180, 0.060, 0.010, 0.003
    hw = h - 2 * tf
    I = b * h ** 3 / 12 - (b - tw) * hw ** 3 / 12
    print(f"spar root shear V = {V:.1f} N (limit); I = {I * 1e12:.4g} mm^4")
    tau_na = float(i_section_shear_stress(V, h, b, tf, tw, 0.0))
    tau_j = float(i_section_shear_stress(V, h, b, tf, tw, hw / 2))
    print(f"VQ/It in the web: neutral axis {tau_na / 1e6:.2f} MPa, web top {tau_j / 1e6:.2f} MPa, "
          f"average V/(h_w t_w) {V / (hw * tw) / 1e6:.2f} MPa, V/(h t_w) {V / (h * tw) / 1e6:.2f} MPa")
    print(f"flange just above the junction (VQ/I b_f): {float(i_section_shear_stress(V, h, b, tf, tw, hw / 2 + 1e-9)) / 1e6:.2f} MPa")
    Q_cap = b * tf * (h - tf) / 2
    q_joint = V * Q_cap / I
    print(f"shear flow at the cap-web joint VQ_cap/I = {q_joint:.0f} N/m ({q_joint / 1e3:.1f} N/mm)")
    P_allow = 4000.0
    print(f"if the caps were fastened to the web with two rows of fasteners rated {P_allow / 1e3:.1f} kN each (assumed): "
          f"pitch at ultimate (1.5) {2 * fastener_pitch(q_joint, P_allow, 1.5) * 1e3:.1f} mm")
    y = np.linspace(-hw / 2, hw / 2, 2001)
    frac = np.trapezoid(i_section_shear_stress(V, h, b, tf, tw, y) * tw, y) / V
    print(f"fraction of V carried by the web (integral of tau t_w): {frac:.4f}")
    print(f"rectangle 1.5 V/A, circle 4V/(3A): {1.5:.2f}, {circle_max_shear(1.0, 1.0) * np.pi:.4f} times V/A")

    # thin-walled centreline model of the same spar
    segs = i_beam_segments(h - tf, b, tf, tw)
    tp = thin_wall_properties(segs)
    flows = open_section_shear_flow(segs, 0.0, V)
    Fx, Fy, Mz = shear_flow_resultant(segs, flows)
    q_web = flows[2][1]
    print(f"\nthin-wall I: I_xx = {tp['Ixx'] * 1e12:.4g} mm^4 (thick formula {I * 1e12:.4g}); "
          f"web q max {np.max(np.abs(q_web)):.0f} N/m -> tau {np.max(np.abs(q_web)) / tw / 1e6:.2f} MPa; "
          f"flange-tip-to-junction q {abs(flows[0][1][-1]):.0f} N/m each")
    print(f"resultant ({Fx:.3f}, {Fy:.3f}) N, moment about centroid {Mz:.2e} N m; "
          f"q at the last free edge {flows[-1][1][-1]:.2e} N/m")

    # channel shear centre
    hc, bc, tc = 0.100, 0.040, 0.002
    cs = channel_segments(hc, bc, tc, tc)
    xs, ys = shear_centre(cs)
    e = channel_shear_centre_closed_form(hc, bc, tc, tc)
    print(f"\nchannel {hc * 1e3:.0f} x {bc * 1e3:.0f} x {tc * 1e3:.0f} mm: shear centre at x = {xs * 1e3:.3f} mm, "
          f"y = {ys * 1e3:.2e} mm; closed form e = {e * 1e3:.3f} mm behind the web")
    cp = thin_wall_properties(cs)
    print(f"centroid at x = {cp['xc'] * 1e3:.3f} mm: a load through the centroid is "
          f"{(cp['xc'] - xs) * 1e3:.2f} mm from the shear centre and twists the channel by "
          f"{(cp['xc'] - xs):.4f} N m per N of shear")
    zs = z_segments(0.040, 0.019, 0.002)
    print(f"Z section shear centre: ({shear_centre(zs)[0] * 1e3:.2e}, {shear_centre(zs)[1] * 1e3:.2e}) mm")

    # independent check
    g = shear_flow_from_bending_gradient(segs, 0.0, V)
    d = max(np.max(np.abs(np.interp(fs, gs, gq) - fq)) for (fs, fq), (gs, gq) in zip(flows, g))
    print(f"bending-gradient route vs formula, I-beam: max |dq| = {d:.2e} N/m "
          f"({d / np.max(np.abs(q_web)):.1e} relative)")

    # figure
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    yy = np.linspace(-h / 2, h / 2, 801)
    ax[0].plot(i_section_shear_stress(V, h, b, tf, tw, yy) / 1e6, yy * 1e3)
    ax[0].set_xlabel("tau = VQ/(I b) (MPa)")
    ax[0].set_ylabel("y (mm)")
    ax[0].set_title("spar root: vertical shear stress")
    for (s, q), ((x1, y1), (x2, y2), _) in zip(flows, segs):
        L = np.hypot(x2 - x1, y2 - y1)
        X = x1 + (x2 - x1) * s / L
        Y = y1 + (y2 - y1) * s / L
        ax[1].plot(X * 1e3, Y * 1e3, "k", lw=0.8)
        sc = 0.02 / np.max(np.abs(q_web))
        nx, ny = -(y2 - y1) / L, (x2 - x1) / L
        ax[1].plot((X + sc * q * nx) * 1e3, (Y + sc * q * ny) * 1e3, "C0")
    ax[1].set_aspect("equal")
    ax[1].set_title("thin-wall shear flow (offset normal to wall)")
    cflows = open_section_shear_flow(cs, 0.0, 1.0)
    for (s, q), ((x1, y1), (x2, y2), _) in zip(cflows, cs):
        L = np.hypot(x2 - x1, y2 - y1)
        X = x1 + (x2 - x1) * s / L
        Y = y1 + (y2 - y1) * s / L
        ax[2].plot(X * 1e3, Y * 1e3, "k", lw=1.5)
    ax[2].plot([xs * 1e3], [ys * 1e3], "C3o", label="shear centre")
    ax[2].plot([cp["xc"] * 1e3], [cp["yc"] * 1e3], "C2s", label="centroid")
    ax[2].set_aspect("equal")
    ax[2].legend()
    ax[2].set_title("channel: shear centre outside the section")
    fig.tight_layout()
    fig.savefig(out_dir / "shear_flow.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'shear_flow.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
