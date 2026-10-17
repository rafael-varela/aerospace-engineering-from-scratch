"""P10.21 Using XFOIL, AVL & OpenVSP: starter code.

Free professional tools (XFOIL for airfoils, AVL for vortex-lattice wings and aircraft, OpenVSP
with VSPAERO for geometry and panel/VLM aerodynamics) are how practising engineers check a design
quickly. They are also how you check your own codes, and they can be wrong in the same ways your
code can: a wrong setting, an unconverged case, a polar file with no record of how it was made.
This module is the bookkeeping and the comparison, and it validates the learner's own codes
against published wind-tunnel data, which is the reference both your code and the tools answer to:

* `read_tool_polar` reads a polar stored as CSV and REFUSES any file whose header does not record
  the tool, its version and its settings (with the settings each tool needs: Reynolds number,
  Mach number, Ncrit, panel count and iteration limit for XFOIL, for example);
* `parse_xfoil_polar` reads the polar file XFOIL writes (its PACC accumulation file) and pulls the
  version, Re, Mach, Ncrit and transition settings from its header; `xfoil_polar_to_labelled_csv`
  adds what XFOIL does not write there (panel count, iteration limit, who ran it and when);
* `polar_metrics` and `compare_polars`: lift-curve slope, zero-lift angle, RMS cl difference over
  a linear range, minimum drag, between any two polars;
* the learner's own codes, as compact given helpers: a NACA 4-digit generator (P10.02), a
  linear-strength vortex panel method (P10.05) and a one-way Thwaites / Michel / turbulent
  momentum-integral / Squire-Young drag estimate (P10.06);
* validation against Ladson, NASA TM-4074 (NACA 0012 at M = 0.15, the shared file
  curriculum/phase-10/data/ladson_naca0012_M015.csv).

XFOIL, AVL and OpenVSP are NOT used, imported or run here; they were not available when this
lesson was written, and no tool output ships with it. Test fixtures in the tests are SYNTHETIC,
written by the tests themselves, and say so in their headers.

Run:  python main.py [out_dir]     (prints the validation, writes tool_validation.png)
Test: pytest curriculum/phase-10/lesson-21-xfoil-avl-openvsp
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import numpy as np

# Settings a polar file must record, by tool (lower-case keys). A file missing any is refused.
REQUIRED_SETTINGS = {
    "xfoil": {"re", "mach", "ncrit", "panels", "iter"},
    "avl": {"mach", "nchord", "nspan", "sref", "cref", "bref"},
    "openvsp": {"mach", "re", "method", "wake_iter", "tessellation"},
}
REQUIRED_META = ("tool", "version", "settings")
LADSON_FILE = "ladson_naca0012_M015.csv"


# ---------------------------------------------------------------- given: data location and loader

def find_phase10_data(start: str | Path | None = None) -> Path:
    """Given: walk up from `start` (default this file) to the repository root that contains
    curriculum/phase-10/data, so the same code runs from curriculum/... and from solutions/...."""
    p = Path(start if start is not None else __file__).resolve()
    for parent in [p] + list(p.parents):
        cand = parent / "curriculum" / "phase-10" / "data"
        if (cand / LADSON_FILE).exists():
            return cand
    raise FileNotFoundError("curriculum/phase-10/data not found above " + str(p))


def load_ladson(table: str = "I", reynolds: float = 5.97e6, path: str | Path | None = None) -> dict:
    """Given: read one Reynolds-number block of Ladson (NASA TM-4074) Table I (free transition) or
    Table VII (wraparound grit). Returns alpha_deg, cl, cd (NaN where the report prints none), cm.
    Values are used exactly as published; nothing is corrected here."""
    path = Path(path) if path else find_phase10_data() / LADSON_FILE
    rows = [r for r in csv.DictReader(line for line in path.open() if not line.startswith("#"))]
    sel = [r for r in rows if r["table"] == table and abs(float(r["reynolds"]) - reynolds) < 0.01 * reynolds]
    if not sel:
        raise ValueError(f"no rows for table {table} at Re = {reynolds:g}")
    f = lambda k: np.array([float(r[k]) if r[k] != "" else np.nan for r in sel])
    return {"alpha_deg": f("alpha_deg"), "cl": f("cl"), "cd": f("cd"), "cm": f("cm"),
            "reynolds": float(sel[0]["reynolds"]), "table": table,
            "source": "Ladson, NASA TM-4074 (1988), Table " + table}


# ---------------------------------------------------------------- given: the learner's own codes (compact)

def naca4(code: str, n: int = 80, closed_te: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """Given (from P10.02): NACA 4-digit coordinates, cosine spacing, ordered from the trailing edge
    along the lower surface to the leading edge and back along the upper surface (2n + 1 points, chord 1)."""
    m, p, t = int(code[0]) / 100, int(code[1]) / 10, int(code[2:]) / 100
    x = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n + 1)))
    a4 = -0.1036 if closed_te else -0.1015
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 + a4 * x ** 4)
    yc, dyc = np.zeros_like(x), np.zeros_like(x)
    if m > 0:
        f = x < p
        yc[f] = m / p ** 2 * (2 * p * x[f] - x[f] ** 2)
        dyc[f] = 2 * m / p ** 2 * (p - x[f])
        yc[~f] = m / (1 - p) ** 2 * (1 - 2 * p + 2 * p * x[~f] - x[~f] ** 2)
        dyc[~f] = 2 * m / (1 - p) ** 2 * (p - x[~f])
    th = np.arctan(dyc)
    xu, yu = x - yt * np.sin(th), yc + yt * np.cos(th)
    xl, yl = x + yt * np.sin(th), yc - yt * np.cos(th)
    return np.concatenate([xl[::-1], xu[1:]]), np.concatenate([yl[::-1], yu[1:]])


def panel_system(X: np.ndarray, Y: np.ndarray) -> dict:
    """Given (compact version of P10.05): the influence matrices of the linear-strength vortex panel
    method (Kuethe-Chow form) for nodes running clockwise from the trailing edge. They depend only on
    the geometry, so they are factorised once and reused for every angle of attack."""
    from scipy.linalg import lu_factor
    N = len(X) - 1
    xc, yc = 0.5 * (X[:-1] + X[1:]), 0.5 * (Y[:-1] + Y[1:])
    S = np.hypot(np.diff(X), np.diff(Y))
    th = np.arctan2(np.diff(Y), np.diff(X))
    xi, eta = xc[:, None] - X[None, :-1], yc[:, None] - Y[None, :-1]
    thi, thj, Sj = th[:, None], th[None, :], S[None, :]
    A = -xi * np.cos(thj) - eta * np.sin(thj)
    Bq = xi ** 2 + eta ** 2
    Cc, D = np.sin(thi - thj), np.cos(thi - thj)
    E = xi * np.sin(thj) - eta * np.cos(thj)
    F = np.log(1 + Sj * (Sj + 2 * A) / Bq)
    G = np.arctan2(E * Sj, Bq + A * Sj)
    P = xi * np.sin(thi - 2 * thj) + eta * np.cos(thi - 2 * thj)
    Q = xi * np.cos(thi - 2 * thj) - eta * np.sin(thi - 2 * thj)
    Cn2 = D + 0.5 * Q * F / Sj - (A * Cc + D * E) * G / Sj
    Cn1 = 0.5 * D * F + Cc * G - Cn2
    Ct2 = Cc + 0.5 * P * F / Sj + (A * D - Cc * E) * G / Sj
    Ct1 = 0.5 * Cc * F - D * G - Ct2
    i = np.arange(N)
    Cn1[i, i], Cn2[i, i], Ct1[i, i], Ct2[i, i] = -1.0, 1.0, 0.5 * np.pi, 0.5 * np.pi
    An = np.zeros((N + 1, N + 1)); At = np.zeros((N, N + 1))
    An[:N, :N] += Cn1; An[:N, 1:] += Cn2; At[:, :N] += Ct1; At[:, 1:] += Ct2
    An[N, 0] = An[N, N] = 1.0                      # Kutta condition: gamma_1 + gamma_N+1 = 0
    return {"lu": lu_factor(An), "At": At, "xc": xc, "yc": yc, "S": S, "theta": th,
            "chord": X.max() - X.min()}


def vortex_panel(X: np.ndarray, Y: np.ndarray, alpha: float, system: dict | None = None) -> dict:
    """Given (compact version of P10.05): linear-strength vortex panel solution at alpha [rad]. Returns
    panel midpoints xc, yc, lengths S, angles, tangential velocity Vt/V (positive along the panel
    direction), Cp and cl (chord 1). Pass `system` from panel_system to reuse the factorisation."""
    from scipy.linalg import lu_solve
    s = system if system is not None else panel_system(X, Y)
    th, S = s["theta"], s["S"]
    g = lu_solve(s["lu"], np.concatenate([np.sin(th - alpha), [0.0]]))
    Vt = np.cos(th - alpha) + s["At"] @ g
    Gamma = 2 * np.pi * np.sum(0.5 * (g[:-1] + g[1:]) * S)
    return {"xc": s["xc"], "yc": s["yc"], "S": S, "theta": th, "Vt": Vt, "Cp": 1 - Vt ** 2,
            "cl": 2 * Gamma / s["chord"]}


def _thwaites_H(lam: np.ndarray) -> np.ndarray:
    lam = np.clip(lam, -0.09, 0.1)
    return np.where(lam >= 0, 2.61 - 3.75 * lam + 5.24 * lam ** 2, 2.088 + 0.0731 / (lam + 0.14))


def _surface_bl(s, Ue, x, Re, x_tr):
    """March one surface from the stagnation point: Thwaites laminar; transition by Michel's
    criterion, laminar separation (lambda < -0.09) or a forced x_tr; then a turbulent momentum
    integral with fixed H = 1.4 and Ludwieg-Tillmann skin friction (a simplification of Head's method)."""
    nu = 1.0 / Re
    I = np.concatenate([[0.0], np.cumsum(0.5 * (Ue[1:] ** 5 + Ue[:-1] ** 5) * np.diff(s))])
    with np.errstate(divide="ignore", invalid="ignore"):
        th2 = 0.45 * nu * I / Ue ** 6
    th2[0] = th2[1]
    theta = np.sqrt(th2)
    dUe = np.gradient(Ue, s)
    lam = th2 * dUe / nu
    H = _thwaites_H(lam)
    itr, reason = len(s) - 1, "none"
    for i in range(2, len(s)):
        Rex, Ret = Ue[i] * s[i] / nu, Ue[i] * theta[i] / nu
        if x_tr is not None and x[i] >= x_tr:
            itr, reason = i, "forced"; break
        if Ret > 1.174 * (1 + 22400 / Rex) * Rex ** 0.46:
            itr, reason = i, "Michel"; break
        if lam[i] < -0.09:
            itr, reason = i, "laminar separation"; break

    def rhs(th, ue, due):
        cf = 0.246 * 10 ** (-0.678 * 1.4) * max(ue * th / nu, 50.0) ** -0.268
        return 0.5 * cf - 3.4 * th / ue * due

    for i in range(itr, len(s) - 1):
        h = s[i + 1] - s[i]
        k1 = rhs(theta[i], Ue[i], dUe[i])
        k2 = rhs(max(theta[i] + h * k1, 1e-12), Ue[i + 1], dUe[i + 1])
        theta[i + 1] = max(theta[i] + 0.5 * h * (k1 + k2), 1e-12)
        H[i + 1] = 1.4
    return {"theta": theta, "H": H, "x_tr": float(x[itr]), "reason": reason}


def bl_drag(panel: dict, Re: float, x_tr: tuple = (None, None), x_te: float = 0.99) -> dict:
    """Given (compact version of P10.06): split the panel solution at the stagnation point, march each
    surface, and apply Squire-Young, cd = sum 2 theta (Ue/V)^((H + 5)/2), at x/c = x_te. The default
    0.99 stops short of the last panel, where the inviscid velocity plunges toward a trailing-edge
    stagnation point that the real, viscous flow does not have. x_tr = (upper, lower) forced transition."""
    Vt, S, xc = panel["Vt"], panel["S"], panel["xc"]
    i0 = int(np.where((Vt[:-1] < 0) & (Vt[1:] >= 0))[0][0])
    sc = np.concatenate([[0.0], np.cumsum(0.5 * (S[:-1] + S[1:]))])
    s_stag = sc[i0] + (-Vt[i0] / (Vt[i0 + 1] - Vt[i0])) * (sc[i0 + 1] - sc[i0])
    out, cd = {}, 0.0
    for name, idx, xt in (("upper", np.arange(i0 + 1, len(Vt)), x_tr[0]), ("lower", np.arange(i0, -1, -1), x_tr[1])):
        s = np.concatenate([[0.0], np.abs(sc[idx] - s_stag)])
        Ue = np.concatenate([[0.0], np.abs(Vt[idx])])
        x = np.concatenate([[xc[idx[0]]], xc[idx]])
        keep = np.concatenate([[True], x[1:] <= x_te])
        r = _surface_bl(s[keep], Ue[keep], x[keep], Re, xt)
        cd += 2 * r["theta"][-1] * Ue[keep][-1] ** ((r["H"][-1] + 5) / 2)
        out[name] = r
    out["cd"] = float(cd)
    return out


def write_tool_polar(path: str | Path, meta: dict, alpha_deg, cl, cd=None, cm=None) -> Path:
    """Given: write a labelled polar CSV: '# key: value' header lines (tool, version, settings as
    'k=v; k=v', plus anything else in meta), then alpha_deg,cl,cd,cm columns (blank where absent)."""
    path = Path(path)
    n = len(alpha_deg)
    cd = [None] * n if cd is None else cd
    cm = [None] * n if cm is None else cm
    with path.open("w") as fh:
        for k, v in meta.items():
            if k == "settings" and isinstance(v, dict):
                v = "; ".join(f"{a}={b}" for a, b in v.items())
            fh.write(f"# {k}: {v}\n")
        fh.write("alpha_deg,cl,cd,cm\n")
        for row in zip(alpha_deg, cl, cd, cm):
            fh.write(",".join("" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{v:.6g}" for v in row) + "\n")
    return path


# ---------------------------------------------------------------- the learner's functions

def parse_header(lines: list[str]) -> dict:
    """Parse '# key: value' comment lines into a dict with lower-case keys. The 'settings' value
    'k1=v1; k2=v2' becomes a dict with lower-case keys."""
    raise NotImplementedError


def read_tool_polar(path: str | Path) -> dict:
    """Read a labelled polar CSV. Refuses (ValueError) a file without a header recording the tool,
    its version and its settings, or whose settings miss any that tool needs (REQUIRED_SETTINGS).
    Returns {'meta': dict, 'alpha_deg', 'cl', 'cd', 'cm': arrays (NaN where blank)}."""
    raise NotImplementedError


def parse_xfoil_polar(text: str) -> dict:
    """Parse the polar (PACC) file XFOIL writes: version, airfoil name, Mach, Re, Ncrit and the forced
    transition locations from the header, then the alpha, CL, CD, CDp, CM, Top_Xtr, Bot_Xtr table.
    Returns {'meta': {...}, 'alpha_deg', 'cl', 'cd', 'cdp', 'cm', 'xtr_top', 'xtr_bot'}."""
    raise NotImplementedError


def xfoil_polar_to_labelled_csv(text: str, path: str | Path, panels: int, iter_limit: int,
                                produced_by: str) -> Path:
    """Convert an XFOIL polar file to a labelled CSV. The panel count and iteration limit are not
    in XFOIL's file, so the person who ran it must supply them, with their name and the date."""
    raise NotImplementedError


def polar_metrics(alpha_deg, cl, cd=None, lin_range: tuple = (-4.0, 8.0)) -> dict:
    """Least-squares line cl = a (alpha - alpha_L0) over lin_range [deg]: cl_alpha per deg and per rad,
    alpha_L0 [deg], the RMS residual of the line; cd_min and the alpha where it occurs (if cd given);
    cl_max over the whole polar."""
    raise NotImplementedError


def compare_polars(p1: dict, p2: dict, lin_range: tuple = (-4.0, 8.0)) -> dict:
    """Compare polar p2 against reference p1 (dicts with alpha_deg, cl, cd): RMS of cl2 - cl1 over the
    reference's alphas inside lin_range and inside p2's alpha span (p2 interpolated linearly), and the
    differences in cl_alpha (relative), alpha_L0 (deg) and cd_min (relative)."""
    raise NotImplementedError


def own_polar(code: str, alpha_deg, Re: float, n: int = 100, x_tr: tuple = (None, None)) -> dict:
    """The learner's own codes as a polar: panel cl and the one-way boundary-layer cd at each alpha.
    Labelled with what produced it, so it can be compared like any tool polar."""
    raise NotImplementedError


def validate_against_ladson(Re: float = 5.97e6, n: int = 100, lin_range: tuple = (-4.5, 8.5),
                            table: str = "I", x_tr: tuple = (None, None)) -> dict:
    """Run own_polar at the published alphas of one Ladson block and compare: the metrics of each, the
    RMS cl difference over lin_range, and the cd ratio at each alpha where the report gives cd."""
    raise NotImplementedError


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".", quick: bool = False) -> None:
    """Validate the learner's own codes against Ladson's NACA 0012 data; write tool_validation.png.
    quick=True (used by the tests) skips the 400-panel case, the grit-strip comparison and the
    Squire-Young sensitivity."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    print("XFOIL, AVL and OpenVSP were not available when this lesson was written: no tool output ships.")
    print("Reference: Ladson, NASA TM-4074, NACA 0012 at M = 0.15 (authentic, public domain).\n")

    for Re in (2.0e6, 5.97e6):
        d = load_ladson("I", Re)
        m = polar_metrics(d["alpha_deg"], d["cl"], d["cd"], (-4.5, 8.5))
        print(f"Ladson Table I, Re = {Re:.3g}: cl_alpha = {m['cl_alpha_deg']:.5f}/deg ({m['cl_alpha_rad']:.3f}/rad), "
              f"alpha_L0 = {m['alpha_L0_deg']:+.3f} deg, line RMS {m['rms_line']:.4f}, cd_min = {m['cd_min']:.5f}, cl_max = {m['cl_max']:.3f}")

    # panel convergence of cl_alpha (own code, inviscid)
    print("\nown panel code, NACA 0012, cl at 4 deg vs panel count (inviscid):")
    cls = []
    for n in (25, 50, 100) + (() if quick else (200,)):
        X, Y = naca4("0012", n)
        cls.append(vortex_panel(X, Y, np.radians(4.0))["cl"])
        print(f"  {2*n:4d} panels: cl = {cls[-1]:.6f}")
    e = np.abs(np.diff(cls))
    print(f"  successive differences {', '.join(f'{x:.2e}' for x in e)}; observed order {np.round(np.log(e[:-1]/e[1:])/np.log(2), 2)}")

    results = {}
    for Re in (2.0e6, 5.97e6):
        v = validate_against_ladson(Re)
        results[Re] = v
        c = v["compare"]
        print(f"\nvalidation at Re = {Re:.3g} (Table I, free transition), linear range -4.5..8.5 deg, {c['n']} points:")
        print(f"  cl_alpha: own {c['m2']['cl_alpha_deg']:.5f}/deg vs data {c['m1']['cl_alpha_deg']:.5f}/deg "
              f"({100*c['rel_d_cl_alpha']:+.1f} %); alpha_L0 difference {c['d_alpha_L0_deg']:+.3f} deg")
        print(f"  RMS cl difference {c['rms_dcl']:.4f} (max {c['max_abs_dcl']:.4f}); cd_min own {c['m2']['cd_min']:.5f} vs data "
              f"{c['m1']['cd_min']:.5f} ({100*c['rel_d_cd_min']:+.1f} %)")
        print("  alpha, cd data, cd own, ratio, own x_tr(top):")
        for a, r in zip(v["alpha_cd"], v["cd_ratio"]):
            j = int(np.argmin(np.abs(v["data"]["alpha_deg"] - a)))
            print(f"   {a:6.2f}  {v['data']['cd'][j]:.5f}  {v['own']['cd'][j]:.5f}  {r:.3f}  {v['own']['xtr_top'][j]:.3f}")
    if not quick:
        vf = validate_against_ladson(6.0e6, table="VII", x_tr=(0.05, 0.05))
        c = vf["compare"]
        print(f"\nTable VII (wraparound grit) at Re = 6.0e6 against own code with transition forced at 0.05c:"
              f" cd_min own {c['m2']['cd_min']:.5f} vs data {c['m1']['cd_min']:.5f} ({100*c['rel_d_cd_min']:+.1f} %),"
              f" cl_alpha {100*c['rel_d_cl_alpha']:+.1f} %")

        # Sensitivity of the own-code cd to the Squire-Young station
        X, Y = naca4("0012", 100)
        p = vortex_panel(X, Y, 0.0)
        print("own cd at alpha = 0, Re = 5.97e6 vs Squire-Young station x/c: " +
              ", ".join(f"{xe}: {bl_drag(p, 5.97e6, x_te=xe)['cd']:.5f}" for xe in (0.95, 0.98, 0.99, 0.995, 1.0)))

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    for Re, col in ((2.0e6, "C0"), (5.97e6, "C1")):
        v = results[Re]
        ax[0].plot(v["data"]["alpha_deg"], v["data"]["cl"], "o", color=col, ms=4, label=f"Ladson TM-4074, Re {Re:.3g}")
        ax[0].plot(v["own"]["alpha_deg"], v["own"]["cl"], "-", color=col, lw=1, label=f"own panel code (inviscid)")
        ax[1].plot(v["data"]["cl"], v["data"]["cd"], "o", color=col, ms=4, label=f"Ladson, Re {Re:.3g}")
        ok = np.abs(v["own"]["alpha_deg"]) <= 10.5
        ax[1].plot(v["own"]["cl"][ok], v["own"]["cd"][ok], "-", color=col, lw=1, label="own panel + boundary layer")
    ax[0].set_xlabel("alpha [deg]"); ax[0].set_ylabel("cl"); ax[0].set_title("NACA 0012 lift: own code vs published data")
    ax[0].legend(fontsize=7); ax[0].grid(alpha=0.3)
    ax[1].set_xlabel("cl"); ax[1].set_ylabel("cd"); ax[1].set_title("drag polar"); ax[1].set_ylim(0, 0.03)
    ax[1].legend(fontsize=7); ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "tool_validation.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
