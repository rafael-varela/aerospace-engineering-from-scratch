"""P10.11 Vortex Lattice Method: starter code.

What are the lift-curve slope, span efficiency and aerodynamic centre of the reference trainer's
wing, and of the swept and tapered wings that lifting-line theory cannot handle? This module builds
a vortex lattice method (VLM) from scratch:

* the planform (span, root chord, taper, quarter-chord sweep, dihedral) is cut into nx x ny panels;
* each panel carries a horseshoe vortex: its bound leg on the panel's quarter-chord line, two
  trailing legs to x = +infinity; the flow tangency condition is enforced at each panel's
  three-quarter-chord control point (the 1/4-3/4 rule, which makes a single panel reproduce the
  2-D lift slope 2 pi exactly);
* the influence-coefficient matrix (normal velocity at every control point from every unit
  horseshoe, by vectorised Biot-Savart) is solved with a dense direct solver (P03.08);
* lift and pitching moment from Kutta-Joukowski on the bound legs, induced drag in the Trefftz
  plane, the spanwise load, and the aerodynamic centre (neutral point of the wing alone)
  x_ac = x_ref - (dCm/dCL) c_ref.

Checks: lattice convergence with the observed order; VLM against a compact lifting-line solver for
unswept high-aspect-ratio wings; symmetry; the 2-D limit CL_alpha -> 2 pi as AR grows; swept wings
against the Helmbold-DATCOM formula CL_alpha = 2 pi A / (2 + sqrt(A^2 (1 + tan^2 Lambda) + 4)).

Axes: x downstream (along the root chord), y to starboard, z up; the wing's root leading edge at the
origin. alpha positive nose-up, Cm positive nose-up about x_ref. Thin flat lifting surface: camber
enters only as a zero-lift-angle shift alpha_eff = alpha - alpha_L0. V_inf = 1 inside the solver.

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes vlm.png)
Test: pytest curriculum/phase-10/lesson-11-vortex-lattice-method
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

RHO = 1.225
G = 9.81
TRAINER = dict(mass=750.0, S=12.0, b=10.0, c=1.2, V_cruise=55.0)
ALPHA_L0_2412 = np.radians(-2.0772)   # thin-airfoil zero-lift angle of the NACA 2412 mean line (P10.04)


# =============================================================== given helpers

def segment_velocity(P: np.ndarray, A: np.ndarray, B: np.ndarray, core: float = 1e-10) -> np.ndarray:
    """Given (from P10.09): velocity at points P (n, 3) from unit-strength straight vortex segments
    A -> B (m, 3), Biot-Savart: v = (1/4 pi) (r1 x r2)/|r1 x r2|^2 r0 . (r1/|r1| - r2/|r2|).
    Returns (n, m, 3)."""
    r0 = (B - A)[None, :, :]
    r1 = P[:, None, :] - A[None, :, :]
    r2 = P[:, None, :] - B[None, :, :]
    c = np.cross(r1, r2)
    c2 = np.sum(c * c, axis=-1)
    n1, n2 = np.linalg.norm(r1, axis=-1), np.linalg.norm(r2, axis=-1)
    small = c2 <= (core * np.sum(r0 * r0, axis=-1)) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        k = np.sum(r0 * (r1 / n1[..., None] - r2 / n2[..., None]), axis=-1) / c2
    return (np.where(small, 0.0, k) / (4 * np.pi))[..., None] * c


def semi_infinite_velocity(P: np.ndarray, A: np.ndarray, d: np.ndarray, core: float = 1e-10) -> np.ndarray:
    """Given (from P10.09): velocity at P (n, 3) from unit semi-infinite vortex lines starting at A
    (m, 3) with vector along the unit vector d: v = (1/4 pi) (d x r1)/|d x r1|^2 (1 + d.r1/|r1|)."""
    r1 = P[:, None, :] - A[None, :, :]
    c = np.cross(np.broadcast_to(d, r1.shape), r1)
    c2 = np.sum(c * c, axis=-1)
    n1 = np.linalg.norm(r1, axis=-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        k = (1 + (r1 @ d) / n1) / c2
    return (np.where(c2 <= core ** 2, 0.0, k) / (4 * np.pi))[..., None] * c


def lifting_line_CL_e(b: float, chord, alpha: float, a0: float = 2 * np.pi, N: int = 60) -> tuple[float, float]:
    """Given (compact P10.10 solver): Prandtl lifting line by Glauert collocation for a symmetric,
    flat, untwisted straight wing with chord(y); returns (CL, e) at alpha (rad)."""
    n = 2 * np.arange(N) + 1
    th = np.arange(1, N + 1) * np.pi / (2 * N)
    y = -0.5 * b * np.cos(th)
    mu = a0 * chord(y) / (4 * b)
    A = np.linalg.solve(np.sin(np.outer(th, n)) * (mu[:, None] * n + np.sin(th)[:, None]), mu * alpha * np.sin(th))
    t = np.linspace(0, np.pi, 4001)
    S = np.trapezoid(chord(-0.5 * b * np.cos(t)) * 0.5 * b * np.sin(t), t)
    AR = b ** 2 / S
    return float(np.pi * AR * A[0]), float(1 / (1 + np.sum(n[1:] * (A[1:] / A[0]) ** 2)))


def helmbold_datcom(AR: float, sweep: float, a0: float = 2 * np.pi) -> float:
    """Given: the Helmbold/DATCOM-type estimate of the lift slope (per rad) of a swept wing,
    CL_alpha = 2 pi A / (2 + sqrt(A^2 (1 + tan^2 Lambda) + 4)) for a0 = 2 pi (Lambda the mid-chord
    sweep; for an untapered wing all chord lines have the same sweep). An empirical fit, not a law."""
    k = a0 / (2 * np.pi)
    return float(a0 * AR / (2 + np.sqrt(AR ** 2 * (1 + np.tan(sweep) ** 2) / k ** 2 + 4)))


def mac(c_root: float, taper: float) -> float:
    """Given: mean aerodynamic chord of a straight-tapered wing, (2/3) c_r (1 + l + l^2)/(1 + l)."""
    return 2 / 3 * c_root * (1 + taper + taper ** 2) / (1 + taper)


# =============================================================== the lattice (learner)

def make_lattice(b: float, c_root: float, taper: float = 1.0, sweep_c4: float = 0.0, dihedral: float = 0.0,
                 nx: int = 4, ny: int = 20, spacing: str = "hough") -> dict:
    """Cut a straight-tapered wing (full span, root LE at the origin) into nx chordwise x ny spanwise
    panels. Spanwise strip edges: "uniform" (edges at the tips), "cosine" (fine at the tips), or
    "hough" (uniform strips of width d = b/(ny + 1/2) with the whole lattice inset by d/4 from each
    tip, Hough's rule, which places the tip vortices where a continuous load would put them and
    restores second-order convergence); chordwise edges uniform. Sweep is of the quarter-chord line (rad); dihedral (rad) raises z = |y| tan(dihedral).
    For every panel: the bound leg A -> B on its own quarter-chord line (left to right), the control
    point at its three-quarter chord on the strip's mid-span line, the unit normal, the area.
    Returns dict(A, B, cp, normal, area, strip (index of each panel's strip), y_edges, z_edges,
    nx, ny, S, b, c_root, taper, sweep_c4)."""
    raise NotImplementedError


def horseshoe_unit_velocity(P: np.ndarray, A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Velocity (n, m, 3) at points P from unit horseshoes: bound leg A -> B plus trailing legs along
    +x from B (vector +x) and into A (vector -x, i.e. a +x line of strength -1 from A)."""
    raise NotImplementedError


def solve_vlm(lat: dict, alpha: float, alpha_L0: float = 0.0) -> np.ndarray:
    """Circulations Gamma (V_inf = 1) of every horseshoe from flow tangency at the control points:
    sum_j AIC_ij Gamma_j = -V_inf . n_i, AIC_ij = (unit horseshoe j's velocity at cp_i) . n_i,
    with the free stream (cos a, 0, sin a) at a = alpha - alpha_L0. Dense LU solve."""
    raise NotImplementedError


def vlm_forces(lat: dict, gamma: np.ndarray, alpha: float, x_ref: float, c_ref: float, alpha_L0: float = 0.0) -> dict:
    """Coefficients from the circulations. Kutta-Joukowski on each bound leg with the free stream
    (F = rho V x Gamma l, per unit rho and V_inf = 1): CL (perpendicular to V), Cm about
    (x_ref, 0, 0) on c_ref (positive nose-up). Induced drag in the Trefftz plane: the trailing
    vortices at the strip edges (y_k, z_k), strength s_k = G_left - G_right of the strip
    circulations G (sum over each strip's chordwise panels), induce a 2-D velocity at each strip's
    centre; CDi = -(1/(q S)) (1/2) sum G_j (v . n_j) ds_j. Also the strip loads.
    Returns dict(CL, CDi, Cm, e, AR, y_strip, G_strip, cl_strip (local cl = 2 G / c))."""
    raise NotImplementedError


def analyse(lat: dict, alpha: float, x_ref: float = None, c_ref: float = None, alpha_L0: float = 0.0) -> dict:
    """solve_vlm + vlm_forces with defaults x_ref = root quarter chord, c_ref = MAC."""
    raise NotImplementedError


def lift_slope_and_ac(lat: dict, d_alpha: float = np.radians(2.0)) -> dict:
    """CL_alpha (per rad) by a central difference about alpha = 0 (the linear VLM makes it exact up
    to the small cos/sin of alpha), and the aerodynamic centre x_ac = x_ref - (dCm/dCL) c_ref with
    x_ref the root quarter chord, c_ref the MAC. Also x_ac as a fraction of the MAC behind the
    MAC's leading edge."""
    raise NotImplementedError


def observed_order(values: list[float], ratio: float = 2.0) -> float:
    """Observed order from three successive refinements by `ratio`:
    p = ln((f1 - f2)/(f2 - f3)) / ln(ratio) (no exact solution needed)."""
    raise NotImplementedError


# =============================================================== main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save vlm.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    m_, S, b, c, V = TRAINER["mass"], TRAINER["S"], TRAINER["b"], TRAINER["c"], TRAINER["V_cruise"]
    AR = b ** 2 / S

    # 1. the trainer
    lat = make_lattice(b, c, nx=6, ny=40)
    r = lift_slope_and_ac(lat)
    CLc = m_ * G / (0.5 * RHO * V ** 2 * S)
    a_c = CLc / r["CL_alpha"] + ALPHA_L0_2412
    rc = analyse(lat, a_c, alpha_L0=ALPHA_L0_2412)
    LL = lifting_line_CL_e(b, lambda y: np.full_like(y, c), 0.1)
    print(f"reference trainer wing (rectangular, AR {AR:.3f}), VLM nx = 6, ny = 40 (Hough spacing):")
    print(f"  CL_alpha = {r['CL_alpha']:.4f} /rad ({r['CL_alpha'] * np.pi / 180:.5f} /deg); lifting line {LL[0] / 0.1:.4f} "
          f"({(r['CL_alpha'] / (LL[0] / 0.1) - 1) * 100:+.2f} %)")
    print(f"  aerodynamic centre x_ac = {r['x_ac']:.4f} m = {r['x_ac_mac']:.4f} MAC behind the LE")
    print(f"  cruise CL {CLc:.4f} at alpha {np.degrees(a_c):.3f} deg: CDi = {rc['CDi']:.6f}, e = {rc['e']:.4f} "
          f"(lifting line e = {LL[1]:.4f}); induced drag {rc['CDi'] * 0.5 * RHO * V ** 2 * S:.1f} N")

    # 2. lattice convergence
    print("\nlattice convergence, trainer wing, alpha = 5 deg: CL for ny = 10, 20, 40, 80 (nx = 4)")
    a5 = np.radians(5.0)
    for sp in ("uniform", "cosine", "hough"):
        vals = [analyse(make_lattice(b, c, nx=4, ny=ny, spacing=sp), a5) for ny in (10, 20, 40, 80)]
        cl = [v["CL"] for v in vals]
        print(f"  {sp:8s}: CL {np.round(cl, 6)}, observed orders {[round(observed_order(cl[k:k + 3]), 2) for k in range(2)]}; "
              f"e {np.round([v['e'] for v in vals], 4)}")
    rowx = [analyse(make_lattice(b, c, nx=nx, ny=40), a5)["CL"] for nx in (1, 2, 4, 8, 16)]
    print(f"  hough ny 40: CL(nx = 1, 2, 4, 8, 16) = {np.round(rowx, 6)}, observed orders in nx "
          f"{[round(observed_order(rowx[k:k + 3]), 2) for k in range(3)]}")

    # 3. VLM vs lifting line, unswept
    print("\nVLM (nx 2, ny 100) vs lifting line (N = 200), unswept untapered wings, CL_alpha per rad:")
    for ARx in (4.0, 6.0, 8.33, 12.0, 20.0, 50.0):
        cx = b / ARx
        lx = lift_slope_and_ac(make_lattice(b, cx, nx=2, ny=100))["CL_alpha"]
        ll = lifting_line_CL_e(b, lambda y: np.full_like(y, cx), 0.1, N=200)[0] / 0.1
        print(f"  AR {ARx:6.2f}: VLM {lx:.4f} ({lx / (2 * np.pi):.4f} x 2 pi), LL {ll:.4f}, difference {(lx / ll - 1) * 100:+.2f} %, "
              f"AR x difference {(lx / ll - 1) * 100 * ARx:+.1f} %")

    # 4. swept wings vs Helmbold-DATCOM
    print("\nswept untapered wings, AR 6: VLM (nx 6, ny 40) vs Helmbold-DATCOM, and the aerodynamic centre")
    swept = []
    for L in (0.0, 15.0, 30.0, 45.0, 60.0):
        cx = b / 6.0
        rr = lift_slope_and_ac(make_lattice(b, cx, sweep_c4=np.radians(L), nx=6, ny=40))
        hd = helmbold_datcom(6.0, np.radians(L))
        swept.append((L, rr["CL_alpha"], hd))
        print(f"  sweep {L:4.1f} deg: VLM {rr['CL_alpha']:.4f}, formula {hd:.4f} ({(rr['CL_alpha'] / hd - 1) * 100:+.2f} %), "
              f"x_ac {rr['x_ac_mac']:.4f} MAC")
    print("taper 0.4, sweep 30 deg, AR 8.33 (an airliner-like planform), and with 5 deg dihedral:")
    for dih in (0.0, 5.0):
        lt = make_lattice(b, 2 * S / (b * 1.4), taper=0.4, sweep_c4=np.radians(30.0), dihedral=np.radians(dih), nx=6, ny=40)
        rr = lift_slope_and_ac(lt)
        aa = analyse(lt, a5)
        print(f"  dihedral {dih:.0f} deg: CL_alpha {rr['CL_alpha']:.4f}, x_ac {rr['x_ac_mac']:.4f} MAC behind the MAC LE, e {aa['e']:.4f}")

    # figure
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    lt = make_lattice(b, 2 * S / (b * 1.4), taper=0.4, sweep_c4=np.radians(30.0), nx=6, ny=20)
    for k in range(len(lt["A"])):
        ax[0].plot([lt["A"][k, 1], lt["B"][k, 1]], [lt["A"][k, 0], lt["B"][k, 0]], "b-", lw=0.6)
    ax[0].plot(lt["cp"][:, 1], lt["cp"][:, 0], "r.", ms=2)
    ax[0].invert_yaxis(); ax[0].set_aspect("equal"); ax[0].set_xlabel("y (m)"); ax[0].set_ylabel("x (m)")
    ax[0].set_title("bound legs (blue) and control points (red)")
    for L in (0.0, 30.0, 45.0):
        aa = analyse(make_lattice(b, b / 6.0, sweep_c4=np.radians(L), nx=6, ny=40), a5)
        ax[1].plot(aa["y_strip"], aa["cl_strip"] / aa["CL"], label=f"sweep {L:.0f} deg")
    ax[1].set_xlabel("y (m)"); ax[1].set_ylabel("cl / CL"); ax[1].set_title("swept AR 6 wings: tip loading"); ax[1].legend()
    Ls = np.array([s[0] for s in swept])
    ax[2].plot(Ls, [s[1] for s in swept], "o-", label="VLM")
    ax[2].plot(Ls, [s[2] for s in swept], "s--", label="Helmbold-DATCOM")
    ax[2].set_xlabel("sweep (deg)"); ax[2].set_ylabel("CL_alpha (/rad)"); ax[2].set_title("AR 6"); ax[2].legend()
    fig.tight_layout()
    fig.savefig(out_dir / "vlm.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
