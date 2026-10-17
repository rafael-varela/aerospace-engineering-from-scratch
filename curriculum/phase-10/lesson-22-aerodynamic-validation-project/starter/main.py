"""P10.22 Aerodynamic Validation Project: starter code.

Does your vortex lattice method (VLM) predict the lift-curve slope of real wings, and how sure can
you be? This module is a small verification-and-validation (V&V) study in the sense of P03.18:

* code verification: the VLM is solved correctly (the span load of a symmetric wing is symmetric,
  lift is linear in angle of attack, a very long unswept wing tends to the 2-D value, the
  half-wing solution with a symmetry image equals the full-wing solution, and the lattice
  converges with a measured observed order);
* solution verification: the discretisation error of each prediction is estimated from three
  lattices by Richardson extrapolation and Roache's grid convergence index (GCI);
* validation: the predictions are compared with the measured lift-curve slopes of twelve
  planforms from NACA TN-1669, TN-1671 and TN-1422 (the shared file
  curriculum/phase-10/data/finite_wing_cla.csv), with the experimental uncertainty (figure-reading
  uncertainty and the ~3 % run-to-run scatter the data README documents), the numerical
  uncertainty, and an input uncertainty from the section lift slope, which is taken from Ladson's
  NACA 0012 data (NASA TM-4074; curriculum/phase-10/data/ladson_naca0012_M015.csv).

Section-slope correction. A flat VLM panel has the thin-airfoil slope 2 pi. A real section has
a0 = kappa 2 pi with kappa < 1. In lifting-line theory a wing of chord c and section slope a0 behaves
exactly like a wing of chord kappa c with slope 2 pi (the chord appears only in the product a0 c),
so CL_alpha(AR, a0) = kappa * CL_alpha_2pi(AR / kappa) at the same span and quarter-chord sweep.
The same scaling (the DATCOM-type slope formula has it built in) is applied to the VLM here: run the
VLM on the planform with every chord multiplied by kappa about the quarter-chord line, multiply
the result by kappa. It is exact for lifting line and an approximation for the lattice.

Axes: x downstream, y to starboard, z up. Angles in radians inside the code; CL_alpha per radian
unless a name says per degree.

Run:  python main.py [out_dir]     (prints the V&V report numbers, writes vlm_validation.png)
Test: pytest curriculum/phase-10/lesson-22-aerodynamic-validation-project
"""
from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

EPS_CORE = 1e-12          # a point closer to a filament line than 1e-6 x the segment length (1e-6 m for a trailing leg) gets nothing from it
SCATTER_REL = 0.03          # run-to-run scatter between TN-1669 and TN-1671 on the same wing (data README)
GRIDS = ((2, 6), (4, 12), (8, 24))   # (chordwise, spanwise per half-wing) panels; refinement ratio 2


# ---------------------------------------------------------------- given: data location

def find_phase10_data(start: str | Path | None = None) -> Path:
    """Given: walk up from `start` (default this file) to the repository root that contains
    curriculum/phase-10/data, so the same code runs from curriculum/... and from solutions/...."""
    p = Path(start if start is not None else __file__).resolve()
    for parent in [p] + list(p.parents):
        cand = parent / "curriculum" / "phase-10" / "data"
        if (cand / "finite_wing_cla.csv").exists():
            return cand
    raise FileNotFoundError("curriculum/phase-10/data not found above " + str(p))


def _read_commented_csv(path: Path) -> list[dict]:
    return list(csv.DictReader(line for line in path.open() if not line.startswith("#")))


# ---------------------------------------------------------------- given: the VLM (compact version of P10.11)

def segment_velocity(P, A, B):
    """Given: velocity at P from unit vortex segments A -> B (Biot-Savart); zero on the segment's line."""
    r1, r2, r0 = P - A, P - B, B - A
    c = np.cross(r1, r2)
    c2 = np.sum(c * c, axis=-1)
    n1, n2 = np.linalg.norm(r1, axis=-1), np.linalg.norm(r2, axis=-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        k = np.sum(r0 * (r1 / n1[..., None] - r2 / n2[..., None]), axis=-1) / (4 * np.pi * c2)
    return c * np.where(c2 > EPS_CORE * np.sum(r0 * r0, axis=-1) ** 2, k, 0.0)[..., None]


def semi_infinite_velocity(P, A):
    """Given: velocity at P from a unit vortex running from A to +infinity along +x."""
    d = np.array([1.0, 0.0, 0.0])
    r1 = P - A
    c = np.cross(d, r1)
    c2 = np.sum(c * c, axis=-1)
    n1 = np.linalg.norm(r1, axis=-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        k = (1 + r1[..., 0] / n1) / (4 * np.pi * c2)
    return c * np.where(c2 > EPS_CORE, k, 0.0)[..., None]


def horseshoe_w(P: np.ndarray, A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Given: z-velocity at points P (m, 3) from unit horseshoes A -> B (n, 3), trailing to +x: (m, n)."""
    P, A, B = P[:, None, :], A[None], B[None]
    V = segment_velocity(P, A, B) + semi_infinite_velocity(P, B) - semi_infinite_velocity(P, A)
    return V[..., 2]


@dataclass
class Lattice:
    """Panels of the starboard half-wing (or the whole wing): bound-vortex ends A, B, control points
    cp, panel areas; reference area S (whole wing) and span b."""
    A: np.ndarray
    B: np.ndarray
    cp: np.ndarray
    area: np.ndarray
    S: float
    b: float
    half: bool


def trapezoidal_lattice(AR: float, taper: float, sweep_c4_deg: float, nc: int, ns: int,
                        chord_scale: float = 1.0, half: bool = True, b: float = 1.0) -> Lattice:
    """Given: a straight-tapered wing of aspect ratio AR, taper ratio (tip/root), quarter-chord sweep,
    span b, with streamwise tips; nc chordwise x ns spanwise panels per half (uniform spacing).
    chord_scale multiplies every chord about the quarter-chord line (the section-slope correction).
    S is the area of the unscaled planform times chord_scale."""
    S0 = b * b / AR
    cr = 2 * S0 / (b * (1 + taper))
    tl = np.tan(np.radians(sweep_c4_deg))
    ys = np.linspace(0, b / 2, ns + 1)
    if not half:
        ys = np.concatenate([-ys[::-1], ys[1:]])
    A, B, cp, area = [], [], [], []
    for j in range(len(ys) - 1):
        y0, y1 = ys[j], ys[j + 1]
        def chord(y):
            return chord_scale * cr * (1 - (1 - taper) * 2 * abs(y) / b)
        def xle(y):
            return abs(y) * tl - 0.25 * chord(y)
        for i in range(nc):
            f0, f1 = i / nc, (i + 1) / nc
            fq, fc = f0 + 0.25 / nc, f0 + 0.75 / nc
            A.append([xle(y0) + fq * chord(y0), y0, 0.0])
            B.append([xle(y1) + fq * chord(y1), y1, 0.0])
            ym = 0.5 * (y0 + y1)
            cp.append([xle(ym) + fc * chord(ym), ym, 0.0])
            area.append(0.5 * (chord(y0) + chord(y1)) * (y1 - y0) / nc)
    return Lattice(np.array(A), np.array(B), np.array(cp), np.array(area), S0 * chord_scale, b, half)


def solve_vlm(lat: Lattice, alpha: float) -> np.ndarray:
    """Given: horseshoe circulations (per unit V) for flow tangency at alpha on a flat lattice. With
    lat.half, the port half is the mirror image (A' = mirror(B), B' = mirror(A), same circulation)."""
    W = horseshoe_w(lat.cp, lat.A, lat.B)
    if lat.half:
        m = np.array([1.0, -1.0, 1.0])
        W = W + horseshoe_w(lat.cp, lat.B * m, lat.A * m)
    return np.linalg.solve(W, -np.sin(alpha) * np.ones(len(lat.cp)))


def lift_coefficient(lat: Lattice, gamma: np.ndarray) -> float:
    """Given: CL = 2 sum(Gamma * dy) / S (Kutta-Joukowski on the bound vortices, V = 1, small angles),
    doubled for a half-wing lattice."""
    dy = lat.B[:, 1] - lat.A[:, 1]
    return float((2.0 if lat.half else 1.0) * 2 * np.sum(gamma * dy) / lat.S)


# ---------------------------------------------------------------- the learner's functions

def vlm_cla(AR: float, taper: float, sweep_c4_deg: float, nc: int, ns: int, kappa: float = 1.0,
            half: bool = True) -> float:
    """Lift-curve slope per radian from the VLM with the section-slope correction: kappa times the
    2 pi VLM slope of the planform whose chords are scaled by kappa (span and sweep unchanged)."""
    raise NotImplementedError


def observed_order_3(f1: float, f2: float, f3: float, r: float = 2.0) -> float:
    """Observed order from three solutions on grids refined by r (f1 finest, f3 coarsest):
    p = ln((f3 - f2)/(f2 - f1)) / ln r. Requires monotone convergence."""
    raise NotImplementedError


def richardson(f1: float, f2: float, p: float, r: float = 2.0) -> float:
    """Richardson extrapolation of the two finest solutions: f1 + (f1 - f2)/(r^p - 1)."""
    raise NotImplementedError


def gci_fine(f1: float, f2: float, p: float, r: float = 2.0, Fs: float = 1.25) -> float:
    """Roache's grid convergence index on the fine grid, as a fraction of f1:
    GCI = Fs |(f2 - f1)/f1| / (r^p - 1). With Fs = 1.25 it is meant as a ~95 % band."""
    raise NotImplementedError


def grid_study(AR: float, taper: float, sweep_c4_deg: float, kappa: float = 1.0, grids=GRIDS) -> dict:
    """CL_alpha on three lattices refined by 2 in both directions, with the observed order (bounded to
    0.5 <= p <= 2 for the error estimate, as the GCI practice recommends), the Richardson value and the
    GCI. u_num is GCI / 2 (a standard uncertainty)."""
    raise NotImplementedError


def section_slope_from_ladson(reynolds: float = 2.0e6, lin_range: tuple = (-4.5, 8.5),
                              table: str = "I") -> tuple[float, float]:
    """Section lift slope a0 [per rad] and its standard error from a least-squares line through
    Ladson's (NASA TM-4074) cl(alpha) for one Reynolds-number block, inside lin_range [deg]."""
    raise NotImplementedError


def load_finite_wing_data() -> list[dict]:
    """The 15 rows of finite_wing_cla.csv with numbers parsed: aspect_ratio, taper, sweep (deg),
    cla_per_deg, and the stated reading uncertainty u_read_per_deg (TN-1669 +/-0.0005, TN-1671
    +/-0.0003, TN-1422 half its last printed digit, 0.0005)."""
    raise NotImplementedError


def validation_uncertainty(cla_meas_deg: float, u_read_deg: float, u_num_deg: float, u_input_deg: float,
                           scatter_rel: float = SCATTER_REL) -> float:
    """Combined standard uncertainty u_val of the comparison error E = S - D (the comparison-error approach of P03.18):
    sqrt(u_num^2 + u_input^2 + u_D^2), with the experimental part u_D^2 = (u_read/sqrt 3)^2 +
    (scatter_rel * D / sqrt 3)^2: both stated as half-widths of a uniform band."""
    raise NotImplementedError


def validation_table(kappa: float, u_kappa: float, grids=GRIDS) -> list[dict]:
    """For every row of the finite-wing file: the predicted CL_alpha per degree (fine lattice), the
    comparison error E = S - D, the relative error, u_num (from the GCI), u_input (from u_kappa by a
    central difference in kappa), u_val, and whether |E| <= 2 u_val. Planforms that occur more than
    once are computed once."""
    raise NotImplementedError


def validation_metrics(table: list[dict]) -> dict:
    """Mean and largest |relative error|, mean signed relative error, and the fraction of rows whose
    error lies within the k = 2 combined uncertainty."""
    raise NotImplementedError


def symmetry_residual(AR: float, taper: float, sweep_c4_deg: float, nc: int = 3, ns: int = 8) -> float:
    """Code verification: solve the full wing (no symmetry image) and return
    max |Gamma(y) - Gamma(-y)| / max |Gamma|."""
    raise NotImplementedError


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".", quick: bool = False) -> None:
    """Print the V&V report numbers and write vlm_validation.png. quick=True (used by the tests) skips
    the 16 x 96 lattice and the uncorrected (kappa = 1) validation table."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    print("CODE VERIFICATION")
    print(f"  symmetry residual, AR 2.61 taper 0.25 sweep 45 (full wing): {symmetry_residual(2.61, 0.25, 45):.2e}")
    full = vlm_cla(5.16, 1.0, 45.0, 3, 8, half=False)
    halfw = vlm_cla(5.16, 1.0, 45.0, 3, 8, half=True)
    print(f"  half-wing with image vs full wing, AR 5.16 sweep 45: {halfw:.10f} vs {full:.10f}")
    print(f"  2-D limit: AR 1000 unswept, CL_alpha = {vlm_cla(1000, 1.0, 0.0, 2, 40):.4f} /rad (2 pi = {2*np.pi:.4f})")
    lat = trapezoidal_lattice(5.16, 1.0, 0.0, 4, 12)
    c1 = lift_coefficient(lat, solve_vlm(lat, np.radians(2)))
    c2 = lift_coefficient(lat, solve_vlm(lat, np.radians(6)))
    print(f"  linearity: CL(6 deg)/CL(2 deg) = {c2/c1:.10f}, sin ratio {np.sin(np.radians(6))/np.sin(np.radians(2)):.10f}")
    print("  lattice convergence, AR 5.16 unswept, kappa = 1 (r = 2 in both directions):")
    seq = [(1, 3), (2, 6), (4, 12), (8, 24)] + ([] if quick else [(16, 48)])
    vals = [vlm_cla(5.16, 1.0, 0.0, nc, ns) for nc, ns in seq]
    for (nc, ns), v in zip(seq, vals):
        print(f"    {nc:2d} x {2*ns:3d}: CL_alpha = {v:.6f} /rad")
    p = [observed_order_3(vals[i + 2], vals[i + 1], vals[i]) for i in range(len(vals) - 2)]
    print(f"    observed orders: {np.round(p, 3)}")

    print("\nINPUT: section lift slope from Ladson (NASA TM-4074), Table I, -4.5 to 8.5 deg")
    a0_2, se_2 = section_slope_from_ladson(2.0e6)
    a0_6, se_6 = section_slope_from_ladson(5.97e6)
    kappa = a0_2 / (2 * np.pi)
    u_kappa = abs(a0_6 - a0_2) / (2 * np.pi)
    print(f"  a0 = {a0_2:.4f} +/- {se_2:.4f} /rad at Re 2.0e6; {a0_6:.4f} +/- {se_6:.4f} /rad at Re 5.97e6")
    print(f"  kappa = a0/(2 pi) = {kappa:.4f}; u_kappa = |a0(6e6) - a0(2e6)|/(2 pi) = {u_kappa:.4f} (Re extrapolation to ~1e6)")

    print("\nSOLUTION VERIFICATION and VALIDATION (fine lattice 8 x 48)")
    tab = validation_table(kappa, u_kappa)
    print("  report      AR   taper sweep  measured  predicted   E/D     p_obs  GCI%   u_num   u_input  u_val   |E|<=2u")
    for t in tab:
        print(f"  {t['report'][5:]:8s} {t['AR']:5.2f}  {t['taper']:4.2f}  {t['sweep']:4.0f}   {t['cla_per_deg']:.4f}   {t['pred_per_deg']:.4f}  "
              f"{100*t['rel_err']:+6.1f}%  {t['p_obs']:5.2f}  {100*t['gci']:4.2f}  {t['u_num']:.5f} {t['u_input']:.5f} {t['u_val']:.5f}  {'yes' if t['within'] else 'NO'}")
    m = validation_metrics(tab)
    print(f"  all {m['n']} rows: mean |E/D| = {100*m['mean_abs_rel']:.1f} %, max {100*m['max_abs_rel']:.1f} % ({m['worst']}),"
          f" mean signed {100*m['mean_rel']:+.1f} %, within k = 2: {100*m['fraction_within']:.0f} %")
    sub = [t for t in tab if t["AR"] > 2 and t["sweep"] < 50]
    ms = validation_metrics(sub)
    print(f"  excluding AR 1.34 and 60 deg sweep ({ms['n']} rows): mean |E/D| = {100*ms['mean_abs_rel']:.1f} %,"
          f" max {100*ms['max_abs_rel']:.1f} %, within: {100*ms['fraction_within']:.0f} %")
    if not quick:
        tab1 = validation_table(1.0, 0.0)
        m1 = validation_metrics(tab1)
        print(f"  without the section correction (kappa = 1): mean |E/D| = {100*m1['mean_abs_rel']:.1f} %, mean signed {100*m1['mean_rel']:+.1f} %,"
              f" within: {100*m1['fraction_within']:.0f} %")
    tn1422 = [t for t in tab if "1422" in t["report"]][0]
    print(f"  TN-1422 AR 9 taper 0.4: predicted {tn1422['pred_per_deg']:.4f}/deg (NACA 0012 a0 as a proxy for the 64-/65-210 sections),"
          f" measured 0.085-0.086/deg; report's own calculated value 0.083-0.085")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    D = np.array([t["cla_per_deg"] for t in tab]); S = np.array([t["pred_per_deg"] for t in tab])
    U = np.array([2 * t["u_val"] for t in tab])
    col = ["C3" if t["AR"] < 2 else ("C1" if t["sweep"] > 50 else "C0") for t in tab]
    for d, s, u, c in zip(D, S, U, col):
        ax[0].errorbar(d, s, yerr=u, fmt="o", color=c, ms=4, capsize=2)
    lim = [0.02, 0.095]
    ax[0].plot(lim, lim, "k--", lw=0.8)
    ax[0].set_xlabel("measured CL_alpha [/deg] (NACA TN-1669, TN-1671, TN-1422)")
    ax[0].set_ylabel("VLM with section correction [/deg]")
    ax[0].set_title("validation (bars: k = 2 u_val); red AR 1.34, orange 60 deg", fontsize=9)
    ax[0].grid(alpha=0.3)
    sw = [0, 45, 60]
    for i, ARv in enumerate((1.34, 2.61, 5.16)):
        meas = [t["cla_per_deg"] for t in tab if t["AR"] == ARv and t["report"].endswith("1669")]
        pred = [t["pred_per_deg"] for t in tab if t["AR"] == ARv and t["report"].endswith("1669")]
        ax[1].plot(sw, meas, "o", color=f"C{i}", label=f"TN-1669 AR {ARv}")
        ax[1].plot(sw, pred, "-", color=f"C{i}", lw=1)
    ax[1].set_xlabel("quarter-chord sweep [deg]"); ax[1].set_ylabel("CL_alpha [/deg]")
    ax[1].set_title("untapered wings: data (points) and VLM (lines)", fontsize=9)
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "vlm_validation.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
