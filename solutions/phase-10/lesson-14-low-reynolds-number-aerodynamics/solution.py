"""P10.14 Low Reynolds Number Aerodynamics: reference solution.

At the reference UAV's cruise (V = 12 m/s, chord 0.20 m, Re = 1.6e5, section lift coefficient 0.45),
where does the laminar boundary layer on the upper surface separate, does natural transition happen
first, how long is the laminar separation bubble, and below which Reynolds number would the bubble
no longer close on the airfoil? Why does a thin cambered plate cope with Re ~ 5e4 better than a
12 % thick section? This module answers with a chain of simple, stated models:

* the inviscid surface speed from a compact Hess-Smith panel method (given; lesson P10.05 builds a
  panel method from scratch);
* Thwaites' integral method for the laminar momentum thickness, theta^2 = (0.45 nu / Ue^6) int Ue^5 ds,
  and its separation criterion lambda = theta^2 Ue' / nu = -0.09. The separation point does not depend
  on the Reynolds number; the boundary-layer thickness there scales as 1/sqrt(Re);
* Michel's transition criterion Re_theta > 1.174 (1 + 22400/Re_s) Re_s^0.46 (empirical);
* a laminar-separation-bubble length from the empirical correlation Ue_s l / nu ~ 5e4 (reported values
  span roughly 3e4 to 7e4; treated as uncertain by a factor of about 1.5), giving the Reynolds
  number below which a bubble cannot reattach before the trailing edge, and the pressure recovery
  that is left for a separated shear layer to make;
* the friction drag of a flat plate with laminar, transitional and turbulent boundary layers as a
  function of Re, with a thickness form factor, against the authentic NACA 0012 drag data (Ladson,
  NASA TM-4074) at the high-Re end, and cl_max against Re from the same data.

Units SI; Reynolds numbers based on chord; speeds normalised by V_inf where marked.

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes low_re.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-10/lesson-14-low-reynolds-number-aerodynamics
"""
from __future__ import annotations

import csv
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

RHO = 1.225
NU = 1.461e-5
G = 9.81
UAV_MASS, UAV_SPAN, UAV_CHORD, UAV_V = 1.2, 1.5, 0.20, 12.0
BUBBLE_RE = 5.0e4            # Ue_s * l / nu: empirical laminar-bubble length correlation (stated, uncertain)
THWAITES_SEP = -0.09


# ---------------------------------------------------------------- given helpers: data

def find_data_file(name: str, start: str | Path | None = None) -> Path:
    """Find curriculum/phase-10/data/<name> by searching upward from this file (works from both the
    curriculum/ starter and the solutions/ reference)."""
    here = Path(start or __file__).resolve()
    for parent in [here, *here.parents]:
        cand = parent / "curriculum" / "phase-10" / "data" / name
        if cand.exists():
            return cand
    raise FileNotFoundError(f"curriculum/phase-10/data/{name} not found above {here}")


def read_ladson(path: str | Path | None = None) -> list[dict]:
    """Rows of the authentic Ladson NACA 0012 data (NASA TM-4074, M = 0.15) as dicts of floats;
    cd and l_over_d are NaN where the report prints no value (post-stall)."""
    path = Path(path) if path else find_data_file("ladson_naca0012_M015.csv")
    with open(path) as fh:
        rows = list(csv.DictReader(line for line in fh if not line.startswith("#")))
    out = []
    for r in rows:
        d = {"table": r["table"], "transition": r["transition"]}
        for k in ("mach", "reynolds", "alpha_deg", "cd", "cl", "cm", "l_over_d"):
            d[k] = float(r[k]) if r[k].strip() else float("nan")
        out.append(d)
    return out


def ladson_blocks(rows: list[dict]) -> dict:
    """Group the Ladson rows into {(table, Re): dict of arrays alpha, cl, cd, cm}."""
    blocks: dict = {}
    for r in rows:
        key = (r["table"], r["reynolds"])
        blocks.setdefault(key, {"alpha": [], "cl": [], "cd": [], "cm": []})
        for k, src in (("alpha", "alpha_deg"), ("cl", "cl"), ("cd", "cd"), ("cm", "cm")):
            blocks[key][k].append(r[src])
    return {k: {kk: np.array(vv) for kk, vv in v.items()} for k, v in blocks.items()}


# ---------------------------------------------------------------- given helpers: geometry and panels

def naca4(code: str, n: int = 160) -> tuple[np.ndarray, np.ndarray]:
    """NACA 4-digit section with a closed trailing edge and cosine spacing (from P10.02): n + 1 points
    per surface; returned as one loop from the trailing edge along the lower surface to the leading
    edge and back along the upper surface (clockwise), chord 1."""
    m, p, t = int(code[0]) / 100, int(code[1]) / 10, int(code[2:]) / 100
    x = (1 - np.cos(np.linspace(0, np.pi, n + 1))) / 2
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1036 * x ** 4)
    if m > 0:
        yc = np.where(x < p, m / p ** 2 * (2 * p * x - x ** 2), m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x ** 2))
        dyc = np.where(x < p, 2 * m / p ** 2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))
    else:
        yc, dyc = np.zeros_like(x), np.zeros_like(x)
    th = np.arctan(dyc)
    xu, yu = x - yt * np.sin(th), yc + yt * np.cos(th)
    xl, yl = x + yt * np.sin(th), yc - yt * np.cos(th)
    return np.concatenate([xl[::-1], xu[1:]]), np.concatenate([yl[::-1], yu[1:]])


def hess_smith(X: np.ndarray, Y: np.ndarray, alpha: float) -> dict:
    """Hess-Smith panel method (constant sources per panel plus one uniform vortex strength, Kutta
    condition on the two trailing-edge panels); a compact version of the P10.05 build.
    Returns control points, panel lengths, tangential speed Vt / V_inf (signed along the clockwise
    panel direction), Cp and cl (from the circulation, cl = 2 Gamma / (V c))."""
    N = len(X) - 1
    xm, ym = (X[:-1] + X[1:]) / 2, (Y[:-1] + Y[1:]) / 2
    dx, dy = np.diff(X), np.diff(Y)
    l = np.hypot(dx, dy)
    th = np.arctan2(dy, dx)
    rx1, ry1 = xm[:, None] - X[None, :-1], ym[:, None] - Y[None, :-1]
    rx2, ry2 = xm[:, None] - X[None, 1:], ym[:, None] - Y[None, 1:]
    beta = np.arctan2(rx1 * ry2 - ry1 * rx2, rx1 * rx2 + ry1 * ry2)
    np.fill_diagonal(beta, np.pi)
    with np.errstate(divide="ignore"):
        lnr = np.log(np.hypot(rx2, ry2) / np.hypot(rx1, ry1))
    np.fill_diagonal(lnr, 0.0)
    d = th[:, None] - th[None, :]
    An = (np.sin(d) * lnr + np.cos(d) * beta) / (2 * np.pi)
    At = (np.sin(d) * beta - np.cos(d) * lnr) / (2 * np.pi)
    Bn = np.sum(np.cos(d) * lnr - np.sin(d) * beta, axis=1) / (2 * np.pi)
    Bt = np.sum(np.sin(d) * lnr + np.cos(d) * beta, axis=1) / (2 * np.pi)
    M = np.zeros((N + 1, N + 1))
    rhs = np.zeros(N + 1)
    M[:N, :N], M[:N, N], rhs[:N] = An, Bn, -np.sin(alpha - th)
    M[N, :N], M[N, N] = At[0] + At[-1], Bt[0] + Bt[-1]
    rhs[N] = -np.cos(alpha - th[0]) - np.cos(alpha - th[-1])
    sol = np.linalg.solve(M, rhs)
    Vt = At @ sol[:N] + Bt * sol[N] + np.cos(alpha - th)
    return {"x": xm, "y": ym, "l": l, "Vt": Vt, "Cp": 1 - Vt ** 2, "cl": 2 * sol[N] * np.sum(l)}


def alpha_for_cl(code: str, cl_target: float, n: int = 160) -> float:
    """Angle of attack (rad) at which the panel method gives cl_target (inviscid)."""
    X, Y = naca4(code, n)
    return brentq(lambda a: hess_smith(X, Y, a)["cl"] - cl_target, np.radians(-10), np.radians(15), xtol=1e-12)


def upper_surface_speed(code: str, alpha: float, n: int = 160) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Arc length s (from the stagnation point, chord units), x/c and Ue/V_inf along the upper surface
    from the stagnation point to the trailing edge. The stagnation point is where the signed panel
    speed changes sign; it is inserted as s = 0, Ue = 0."""
    X, Y = naca4(code, n)
    hs = hess_smith(X, Y, alpha)
    Vt = hs["Vt"]
    k = int(np.where((Vt[:-1] < 0) & (Vt[1:] >= 0))[0][-1])        # last sign change: lower -> upper
    f = -Vt[k] / (Vt[k + 1] - Vt[k])
    xs = hs["x"][k] + f * (hs["x"][k + 1] - hs["x"][k])
    ys = hs["y"][k] + f * (hs["y"][k + 1] - hs["y"][k])
    xu = np.concatenate([[xs], hs["x"][k + 1:]])
    yu = np.concatenate([[ys], hs["y"][k + 1:]])
    ue = np.concatenate([[0.0], Vt[k + 1:]])
    s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(xu), np.diff(yu)))])
    return s, xu, ue


@lru_cache(maxsize=64)
def surface_for_cl(code: str, cl: float, n: int = 160) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    """(alpha, s, x, ue) on the upper surface at inviscid lift coefficient cl (cached: the panel solve
    does not depend on the Reynolds number). Treat the returned arrays as read-only."""
    a = alpha_for_cl(code, cl, n)
    s, xu, ue = upper_surface_speed(code, a, n)
    return a, s, xu, ue


# ---------------------------------------------------------------- laminar boundary layer

def thwaites(s: np.ndarray, ue: np.ndarray, nu: float) -> tuple[np.ndarray, np.ndarray]:
    """Thwaites' method: theta^2 = (0.45 nu / ue^6) int_0^s ue^5 ds (exact for piecewise-linear ue) and
    lambda = theta^2 (d ue / ds) / nu. Where ue = 0 (a stagnation point) the limit is used:
    theta^2 = 0.075 nu / (d ue / ds), lambda = 0.075. Returns (theta, lambda)."""
    s = np.asarray(s, dtype=float)
    ue = np.asarray(ue, dtype=float)
    # ue is taken as linear on each interval, for which int ue^5 ds = h (u1^6 - u0^6) / (6 (u1 - u0))
    # exactly; the trapezoid rule would be badly wrong near a stagnation point, where ue^5 ~ s^5.
    u0, u1, h = ue[:-1], ue[1:], np.diff(s)
    du = u1 - u0
    small = np.abs(du) < 1e-12 * np.maximum(np.abs(u0) + np.abs(u1), 1e-300)
    with np.errstate(divide="ignore", invalid="ignore"):
        seg = np.where(small, h * u0 ** 5, h * (u1 ** 6 - u0 ** 6) / (6 * np.where(small, 1.0, du)))
    I = np.concatenate([[0.0], np.cumsum(seg)])
    due = np.gradient(ue, s)
    with np.errstate(divide="ignore", invalid="ignore"):
        th2 = np.where(ue > 0, 0.45 * nu * I / np.where(ue > 0, ue, 1.0) ** 6, 0.0)
    stag = ue <= 0
    if np.any(stag):
        th2[stag] = np.where(due[stag] > 0, 0.075 * nu / np.where(due[stag] > 0, due[stag], 1.0), 0.0)
    lam = th2 * due / nu
    return np.sqrt(th2), lam


def laminar_separation(s: np.ndarray, lam: np.ndarray, start: int = 1) -> float:
    """Arc length where lambda first falls through the Thwaites separation value -0.09 (linear
    interpolation); NaN if the laminar layer survives to the end."""
    for i in range(max(start, 1), len(s)):
        if lam[i] <= THWAITES_SEP < lam[i - 1]:
            f = (lam[i - 1] - THWAITES_SEP) / (lam[i - 1] - lam[i])
            return float(s[i - 1] + f * (s[i] - s[i - 1]))
    return float("nan")


def michel_transition(s: np.ndarray, ue: np.ndarray, theta: np.ndarray, nu: float) -> float:
    """Arc length where Michel's criterion Re_theta >= 1.174 (1 + 22400 / Re_s) Re_s^0.46 is first met
    (Re_s = ue s / nu, Re_theta = ue theta / nu); NaN if never met."""
    Re_s = ue * s / nu
    Re_t = ue * theta / nu
    with np.errstate(divide="ignore", invalid="ignore"):
        crit = 1.174 * (1 + 22400 / Re_s) * Re_s ** 0.46
    ok = (Re_s > 0) & (Re_t >= crit)
    if not np.any(ok):
        return float("nan")
    i = int(np.argmax(ok))
    g = Re_t - crit
    if i == 0 or not np.isfinite(g[i - 1]):
        return float(s[i])
    f = -g[i - 1] / (g[i] - g[i - 1])
    return float(s[i - 1] + f * (s[i] - s[i - 1]))


def bubble_length(ue_sep: float, nu: float, re_bubble: float = BUBBLE_RE) -> float:
    """Length (m, or chords if nu is per chord) from laminar separation to turbulent reattachment:
    l = re_bubble nu / ue_sep (empirical correlation; stated)."""
    return re_bubble * nu / ue_sep


def bubble_analysis(code: str, cl: float, Re: float, re_bubble: float = BUBBLE_RE, n: int = 160) -> dict:
    """Upper-surface laminar analysis of section `code` at inviscid lift coefficient cl and chord
    Reynolds number Re (chord 1, V 1, nu = 1/Re). Returns: alpha (deg), laminar separation x_s/c and
    arc length s_s, ue_s, natural-transition arc length (Michel), whether a bubble forms (separation
    before transition), its length and reattachment arc length, whether it closes before the trailing
    edge, the remaining pressure recovery dCp_rec = ue_s^2 - ue_TE^2, and the minimum Re for closure."""
    a, s, xu, ue = surface_for_cl(code, float(cl), n)
    nu = 1.0 / Re
    theta, lam = thwaites(s, ue, nu)
    i_peak = int(np.argmax(ue))
    s_sep = laminar_separation(s, lam, start=i_peak)
    s_tr = michel_transition(s, ue, theta, nu)
    s_te = float(s[-1])
    if np.isnan(s_sep):
        return {"alpha_deg": np.degrees(a), "s_sep": np.nan, "x_sep": np.nan, "ue_sep": np.nan, "s_tr": s_tr,
                "bubble": False, "length": 0.0, "closes": True, "dCp_rec": 0.0, "Re_min": 0.0, "ue_peak": ue[i_peak],
                "theta_sep_over_c": np.nan, "s_te": s_te}
    ue_s = float(np.interp(s_sep, s, ue))
    x_s = float(np.interp(s_sep, s, xu))
    bubble = not (np.isfinite(s_tr) and s_tr < s_sep)
    l = bubble_length(ue_s, nu, re_bubble)
    closes = s_sep + l < s_te
    return {"alpha_deg": np.degrees(a), "s_sep": s_sep, "x_sep": x_s, "ue_sep": ue_s, "s_tr": s_tr, "bubble": bubble,
            "length": l, "closes": bool(closes), "dCp_rec": ue_s ** 2 - ue[-1] ** 2,
            "Re_min": re_bubble / (ue_s * (s_te - s_sep)), "ue_peak": float(ue[i_peak]),
            "theta_sep_over_c": float(np.interp(s_sep, s, theta)), "s_te": s_te}


def natural_transition_re(code: str, cl: float, n: int = 160) -> float:
    """Chord Reynolds number above which Michel transition occurs before laminar separation (no
    bubble), found by Brent's method on log Re in [1e3, 1e9]; NaN if the regime does not change in
    that range."""
    a, s, xu, ue = surface_for_cl(code, float(cl), n)
    i_peak = int(np.argmax(ue))

    def margin(logRe):
        nu = 10 ** -logRe
        th, lam = thwaites(s, ue, nu)
        s_sep = laminar_separation(s, lam, start=i_peak)
        s_tr = michel_transition(s, ue, th, nu)
        s_sep = s[-1] if np.isnan(s_sep) else s_sep
        s_tr = 10.0 if np.isnan(s_tr) else s_tr
        return s_tr - s_sep

    if margin(3.0) * margin(9.0) > 0:
        return float("nan")
    return float(10 ** brentq(margin, 3.0, 9.0, xtol=1e-6))


# ---------------------------------------------------------------- friction drag vs Reynolds number

def cf_laminar(Re: np.ndarray) -> np.ndarray:
    """Blasius laminar flat-plate mean skin friction 1.328 / sqrt(Re)."""
    return 1.328 / np.sqrt(np.asarray(Re, dtype=float))


def cf_turbulent(Re: np.ndarray) -> np.ndarray:
    """Prandtl-Schlichting turbulent flat-plate mean skin friction 0.455 / (log10 Re)^2.58."""
    return 0.455 / np.log10(np.asarray(Re, dtype=float)) ** 2.58


def cf_mixed(Re: np.ndarray, Re_tr: float = 5e5) -> np.ndarray:
    """Mean skin friction of a plate with laminar flow up to Re_x = Re_tr and turbulent flow after it
    (the turbulent layer counted as if it had started at the leading edge, minus its share over the
    laminar part, plus the laminar friction there):
    Cf = Cf_t(Re) - (Re_tr / Re)(Cf_t(Re_tr) - Cf_l(Re_tr)) for Re > Re_tr, else Cf_l(Re)."""
    Re = np.asarray(Re, dtype=float)
    mixed = cf_turbulent(Re) - (Re_tr / Re) * (cf_turbulent(Re_tr) - cf_laminar(Re_tr))
    return np.where(Re > Re_tr, mixed, cf_laminar(Re))


def section_cd_friction(Re: np.ndarray, t_c: float = 0.12, Re_tr: float | None = 5e5) -> np.ndarray:
    """Profile drag at zero lift from friction on both surfaces with a thickness form factor:
    cd = 2 Cf (1 + 2 t/c + 60 (t/c)^4) (an empirical form factor; see P10.13). Re_tr None: fully turbulent."""
    cf = cf_turbulent(Re) if Re_tr is None else cf_mixed(Re, Re_tr)
    return 2 * cf * (1 + 2 * t_c + 60 * t_c ** 4)


def reynolds_number(V: float, c: float, nu: float = NU) -> float:
    """Chord Reynolds number V c / nu."""
    return V * c / nu


def clmax_vs_re(blocks: dict, table: str) -> tuple[np.ndarray, np.ndarray]:
    """Largest measured cl in each Reynolds-number block of a Ladson table: (Re, cl_max)."""
    keys = sorted(k for k in blocks if k[0] == table)
    return np.array([k[1] for k in keys]), np.array([np.nanmax(blocks[k]["cl"]) for k in keys])


def cd_at_zero_lift(blocks: dict, table: str) -> tuple[np.ndarray, np.ndarray]:
    """Measured cd at the angle of attack closest to zero in each block: (Re, cd)."""
    keys = sorted(k for k in blocks if k[0] == table)
    Re, cd = [], []
    for k in keys:
        b = blocks[k]
        i = int(np.argmin(np.abs(b["alpha"])))
        Re.append(k[1])
        cd.append(b["cd"][i])
    return np.array(Re), np.array(cd)


# ---------------------------------------------------------------- main

FLYERS = [  # illustrative round numbers (speed m/s, chord m), not measurements
    ("fruit fly", 2.0, 0.0008), ("bumblebee", 4.0, 0.003), ("pigeon", 15.0, 0.10),
    ("reference UAV", UAV_V, UAV_CHORD), ("model sailplane", 10.0, 0.25), ("sailplane", 30.0, 0.80),
    ("reference trainer, approach", 30.0, 1.2), ("reference trainer, cruise", 55.0, 1.2),
]


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save low_re.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    print("Reynolds numbers (sea-level air, illustrative speeds and chords):")
    for name, V, c in FLYERS:
        print(f"  {name:30s} V = {V:5.1f} m/s, c = {c:7.4f} m: Re = {reynolds_number(V, c):.3g}")
    S = UAV_SPAN * UAV_CHORD
    cl_uav = UAV_MASS * G / (0.5 * RHO * UAV_V ** 2 * S)
    Re_uav = reynolds_number(UAV_V, UAV_CHORD)
    print(f"UAV cruise lift coefficient {cl_uav:.3f}, Re = {Re_uav:.4g}")

    # Thwaites checks on classical flows
    x = np.linspace(0, 1, 2001)
    th, lam = thwaites(x, np.ones_like(x), 1e-6)
    print(f"flat plate: theta sqrt(U/(nu x)) at x = 1: {th[-1] / np.sqrt(1e-6):.4f} (Blasius 0.664)")
    th, lam = thwaites(x, 1 - x, 1e-6)
    print(f"Howarth linearly retarded flow Ue = 1 - x: Thwaites separation at x = {laminar_separation(x, lam):.4f} "
          f"(exact similarity-solution value about 0.120)")

    # the UAV at cruise
    print(f"upper-surface laminar analysis at cl = {cl_uav:.3f}, Re = {Re_uav:.4g} (bubble Re_l = {BUBBLE_RE:.0e}):")
    results = {}
    for code in ("0012", "2412", "4402"):
        r = bubble_analysis(code, cl_uav, Re_uav)
        results[code] = r
        Re_nat = natural_transition_re(code, cl_uav)
        print(f"  NACA {code}: alpha = {r['alpha_deg']:.2f} deg, peak Ue = {r['ue_peak']:.3f} V, laminar separation "
              f"x/c = {r['x_sep']:.4f} (Ue_s = {r['ue_sep']:.3f} V, theta_s = {r['theta_sep_over_c'] * 1e3:.3f}e-3 c), "
              f"Michel transition s = {r['s_tr']:.3f}, bubble {'yes' if r['bubble'] else 'no'}, length {r['length']:.3f} c, "
              f"closes {r['closes']}, recovery left dCp = {r['dCp_rec']:.3f}, Re_min = {r['Re_min']:.3g}, "
              f"no-bubble above Re = {Re_nat:.3g}")
    print("Reynolds-number sweep, NACA 0012 vs NACA 4402 at the cruise cl:")
    for Re in (2e4, 5e4, 1e5, 1.64e5, 3e5, 1e6):
        a = bubble_analysis("0012", cl_uav, Re)
        b = bubble_analysis("4402", cl_uav, Re)
        print(f"  Re = {Re:8.3g}: 0012 bubble {a['length']:.3f} c ({'closes' if a['closes'] else 'OPEN'}"
              f"{', none: transition first' if not a['bubble'] else ''}); 4402 bubble {b['length']:.3f} c "
              f"({'closes' if b['closes'] else 'OPEN'}{', none: transition first' if not b['bubble'] else ''})")
    print("sensitivity of Re_min to the bubble correlation (NACA 0012 / 4402):")
    for rb in (3e4, 5e4, 7e4):
        a = bubble_analysis("0012", cl_uav, Re_uav, re_bubble=rb)
        b = bubble_analysis("4402", cl_uav, Re_uav, re_bubble=rb)
        print(f"  Re_l = {rb:.0e}: Re_min = {a['Re_min']:.3g} / {b['Re_min']:.3g}")

    # friction drag model vs Re, with the authentic data at the high-Re end
    rows = read_ladson()
    blocks = ladson_blocks(rows)
    ReI, cdI = cd_at_zero_lift(blocks, "I")
    ReV, cdV = cd_at_zero_lift(blocks, "VII")
    print("NACA 0012 cd near alpha = 0, Ladson (NASA TM-4074), vs friction model 2 Cf FF:")
    for Re, cd in zip(ReI, cdI):
        print(f"  free transition Re = {Re:.3g}: measured {cd:.5f}; model Re_tr = 5e5 {float(section_cd_friction(Re)):.5f}, "
              f"Re_tr = 3e6 {float(section_cd_friction(Re, Re_tr=3e6)):.5f}, turbulent {float(section_cd_friction(Re, Re_tr=None)):.5f}")
    for Re, cd in zip(ReV, cdV):
        print(f"  wraparound grit Re = {Re:.3g}: measured {cd:.5f}; turbulent model {float(section_cd_friction(Re, Re_tr=None)):.5f}")
    for Re in (1e4, 5e4, 1.64e5, 5e5):
        print(f"  model at Re = {Re:.3g}: laminar {float(section_cd_friction(Re, Re_tr=1e12)):.5f}, "
              f"Re_tr 5e5 {float(section_cd_friction(Re)):.5f}, turbulent {float(section_cd_friction(Re, Re_tr=None)):.5f}")
    ReC, clI = clmax_vs_re(blocks, "I")
    ReF, clF = clmax_vs_re(blocks, "VII")
    pI = np.polyfit(np.log10(ReC), clI, 1)
    print("cl_max (largest measured cl per block), Ladson NACA 0012:")
    for Re, c in zip(ReC, clI):
        print(f"  free Re = {Re:.3g}: {c:.4f}")
    for Re, c in zip(ReF, clF):
        print(f"  grit Re = {Re:.3g}: {c:.4f}")
    print(f"  free transition: cl_max rises {pI[0]:.3f} per decade of Re (fit over {ReC.min():.2g}-{ReC.max():.2g}); "
          f"straight-line extrapolation to Re = {Re_uav:.3g} would give {np.polyval(pI, np.log10(Re_uav)):.3f} (not valid: "
          f"the bubble physics changes)")

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
    for code, col in (("0012", "C0"), ("4402", "C1")):
        a, s, xu, ue = surface_for_cl(code, float(cl_uav))
        ax[0].plot(xu, ue, col, label=f"NACA {code}, alpha {np.degrees(a):.1f} deg")
        r = results[code]
        ax[0].plot(r["x_sep"], r["ue_sep"], "o", color=col)
    ax[0].set_xlabel("x/c"); ax[0].set_ylabel("Ue / V (upper surface, inviscid)")
    ax[0].set_title(f"cl = {cl_uav:.2f}: dots = Thwaites laminar separation"); ax[0].legend(); ax[0].grid(alpha=.3)
    Res = np.logspace(3.5, 7.5, 200)
    ax[1].loglog(Res, section_cd_friction(Res, Re_tr=1e12), label="laminar (Blasius) x FF")
    ax[1].loglog(Res, section_cd_friction(Res), label="transition at Re_x = 5e5")
    ax[1].loglog(Res, section_cd_friction(Res, Re_tr=None), label="turbulent x FF")
    ax[1].loglog(ReI, cdI, "ks", label="Ladson 0012 free transition (data)")
    ax[1].loglog(ReV, cdV, "k^", mfc="none", label="Ladson 0012 wraparound grit (data)")
    ax[1].axvline(Re_uav, color="grey", ls="--", lw=1)
    ax[1].set_xlabel("Re"); ax[1].set_ylabel("cd at zero lift"); ax[1].legend(fontsize=7); ax[1].grid(alpha=.3, which="both")
    ax[2].semilogx(ReC, clI, "ks-", label="free transition (data)")
    ax[2].semilogx(ReF, clF, "k^--", mfc="none", label="wraparound grit (data)")
    ax[2].set_xlabel("Re"); ax[2].set_ylabel("largest measured cl"); ax[2].set_title("NACA 0012, Ladson NASA TM-4074")
    ax[2].legend(); ax[2].grid(alpha=.3, which="both")
    fig.tight_layout()
    fig.savefig(out_dir / "low_re.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
