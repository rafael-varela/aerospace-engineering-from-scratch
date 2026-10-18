"""P10.13 Drag Build-Up: starter code.

What is the reference trainer's zero-lift drag coefficient C_D0, what is its best lift-to-drag ratio
and the speed for it, and how much power does the airframe's drag absorb in the 55 m/s cruise? And
how good is the skin-friction-times-form-factor estimate on which the answer rests? This module
builds the answer component by component and then tests its key ingredient against authentic data:

* flat-plate skin friction: laminar (Blasius, 1.328/sqrt(Re)), turbulent (Prandtl-Schlichting,
  0.455/(log10 Re)^2.58) and mixed, with transition at a given fraction of the length;
* empirical form factors (wing and tail sections FF = 1 + 2 t/c + 60 (t/c)^4; bodies
  FF = 1 + 60/f^3 + f/400), stated interference factors Q, and wetted areas computed from the
  geometry (the NACA section perimeter and the fuselage's surface of revolution);
* C_D0 = sum Cf FF Q S_wet / S + the fixed landing gear (an assumed drag area) + a stated allowance
  for leakage and protuberances;
* the induced drag from a lifting-line span efficiency with a fuselage correction; the parabolic
  polar C_D = C_D0 + C_L^2 / (pi e AR), best L/D = 0.5 sqrt(pi e AR / C_D0);
* data-based investigation: the 2-D section estimate cd = 2 Cf (1 + 2 t/c + 60 (t/c)^4) against the
  authentic NACA 0012 drag near zero lift (Ladson, NASA TM-4074; free transition and a wraparound
  grit strip) over Re = 2e6 to 1.9e7.

Units SI; ISA sea level.

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes drag_build_up.png)
Test: pytest curriculum/phase-10/lesson-13-drag-build-up
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

RHO = 1.225
NU = 1.461e-5
G = 9.81

# reference trainer (shared Phase 10 conventions)
MASS = 750.0
S_REF = 12.0
SPAN = 10.0
CHORD = 1.2
V_CRUISE = 55.0
L_FUS = 7.0
D_FUS = 1.1
S_HT = 2.0
# stated assumptions for the parts the shared conventions do not fix
S_VT = 1.2              # m^2 vertical tail (exposed)
C_TAIL = 0.8            # m mean tail chord
TC_TAIL = 0.10          # tail sections NACA 0010
GEAR_DQ = 0.15          # m^2 drag area of a fixed tricycle gear with wheel fairings and struts
COOLING_DQ = 0.04       # m^2 engine cooling and induction drag area
K_VISC = 0.38           # lift-dependent profile drag factor in the Oswald estimate (empirical)
MISC_FRACTION = 0.05    # leakage, protuberances, antennas: +5 % of the component sum


# ---------------------------------------------------------------- given helpers

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


def naca4(code: str, n: int = 200) -> tuple[np.ndarray, np.ndarray]:
    """NACA 4-digit section, closed trailing edge, cosine spacing (from P10.02), as one loop from the
    trailing edge round the lower surface, the leading edge and back along the upper surface; chord 1."""
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


def trainer_fuselage_radius(x: np.ndarray) -> np.ndarray:
    """Equivalent-radius profile of the trainer's fuselage (from P10.16; assumed shape): elliptic nose
    0-1.8 m to r = 0.55 m, constant to 3.2 m, tail cone r = 0.55 - 0.47 ((x - 3.2)/3.8)^1.3 to 7.0 m."""
    x = np.asarray(x, dtype=float)
    Rm = D_FUS / 2
    tail = Rm - (Rm - 0.08) * (np.clip(x - 3.2, 0.0, None) / 3.8) ** 1.3
    return np.where(x < 1.8, Rm * np.sqrt(np.clip(1 - (1 - x / 1.8) ** 2, 0.0, None)), np.where(x <= 3.2, Rm, tail))


def lifting_line_e(AR: float, taper: float = 1.0, a0: float = 2 * np.pi, n_terms: int = 41) -> float:
    """Span efficiency e = 1/(1 + delta) of an untwisted trapezoidal wing by Prandtl's lifting line
    (Fourier collocation of the monoplane equation; from P10.10)."""
    N = n_terms
    theta = np.linspace(np.pi / (2 * N), np.pi / 2, N)                     # 0 < theta <= pi/2 (symmetric half)
    n = 2 * np.arange(N) + 1                                               # odd terms only (symmetric load)
    y_over_s = np.cos(theta)                                               # |y| / (b/2)
    c_over_cmean = 2 / (1 + taper) * (1 - (1 - taper) * y_over_s)
    mu = c_over_cmean * a0 / (4 * AR)                                      # c a0 / (4 b) with c = c_mean (c/c_mean)
    A = np.sin(np.outer(theta, n)) * (mu[:, None] * n[None, :] + np.sin(theta)[:, None])
    rhs = mu * np.sin(theta)
    An = np.linalg.solve(A, rhs)
    delta = np.sum(n[1:] * (An[1:] / An[0]) ** 2)
    return float(1 / (1 + delta))


# ---------------------------------------------------------------- skin friction

def cf_laminar(Re: np.ndarray) -> np.ndarray:
    """Blasius laminar flat-plate mean skin friction (both sides counted separately): 1.328 / sqrt(Re)."""
    raise NotImplementedError


def cf_turbulent(Re: np.ndarray) -> np.ndarray:
    """Prandtl-Schlichting turbulent flat-plate mean skin friction: 0.455 / (log10 Re)^2.58."""
    raise NotImplementedError


def cf_mixed(Re: float | np.ndarray, x_tr: float) -> np.ndarray:
    """Mean skin friction with laminar flow up to the fraction x_tr of the length and turbulent flow
    after it: Cf = Cf_t(Re) - x_tr (Cf_t(Re_tr) - Cf_l(Re_tr)), Re_tr = x_tr Re. x_tr = 0 is fully
    turbulent and x_tr = 1 fully laminar."""
    raise NotImplementedError


# ---------------------------------------------------------------- form factors, wetted areas

def form_factor_section(t_c: float) -> float:
    """Empirical section form factor FF = 1 + 2 t/c + 60 (t/c)^4 (velocity increase over the thickness
    plus a pressure-drag term that grows fast for thick sections)."""
    raise NotImplementedError


def form_factor_body(fineness: float) -> float:
    """Empirical body form factor FF = 1 + 60 / f^3 + f / 400, f = length / diameter."""
    raise NotImplementedError


def section_cd(Re: float | np.ndarray, t_c: float, x_tr: float = 0.0) -> np.ndarray:
    """2-D section profile drag at zero lift: cd = 2 Cf_mixed(Re, x_tr) FF(t/c)."""
    raise NotImplementedError


def perimeter_ratio(code: str, n: int = 400) -> float:
    """Perimeter of the NACA section divided by twice the chord: S_wet = 2 (ratio) S_exposed."""
    raise NotImplementedError


def body_wetted_area(x: np.ndarray, r: np.ndarray) -> float:
    """Surface of revolution through the profile points (chain of frustums)."""
    raise NotImplementedError


def component(name: str, L_ref: float, S_wet: float, FF: float, Q: float, V: float = V_CRUISE,
              x_tr: float = 0.0, nu: float = NU) -> dict:
    """One component's drag: Re = V L / nu, Cf = cf_mixed(Re, x_tr), D/q = Cf FF Q S_wet,
    Delta C_D = (D/q) / S_REF."""
    raise NotImplementedError


def trainer_components(V: float = V_CRUISE, x_tr_wing: float = 0.10) -> list[dict]:
    """The trainer's components (stated assumptions in the module constants): wing (NACA 2412, the
    1.32 m^2 inside the fuselage removed, laminar to x_tr_wing), fuselage, horizontal and vertical
    tails (NACA 0010, turbulent), and their interference factors (wing 1.0, fuselage 1.0, tails 1.05)."""
    raise NotImplementedError


def cd0_build_up(components: list[dict], gear_dq: float = GEAR_DQ, cooling_dq: float = COOLING_DQ,
                 misc: float = MISC_FRACTION) -> dict:
    """C_D0 = (sum of component C_D + (gear D/q + cooling D/q) / S) (1 + misc). Returns the parts and
    the total."""
    raise NotImplementedError


def oswald_e(CD0: float, AR: float = SPAN ** 2 / S_REF, d_over_b: float = D_FUS / SPAN,
             k_visc: float = K_VISC) -> float:
    """Oswald efficiency for the whole aeroplane (an empirical estimate in three factors):
    e_inv = e_LL (1 - 2 (d/b)^2), the lifting-line span efficiency reduced for the fuselage's break in
    the span load; then 1/e = 1/e_inv + k_visc pi AR C_D0 adds the profile drag that grows with C_L^2.
    k_visc = 0 gives the inviscid value."""
    raise NotImplementedError


def polar(CL: np.ndarray, CD0: float, e: float, AR: float = SPAN ** 2 / S_REF) -> np.ndarray:
    """Parabolic drag polar C_D = C_D0 + C_L^2 / (pi e AR)."""
    raise NotImplementedError


def best_ld(CD0: float, e: float, AR: float = SPAN ** 2 / S_REF, mass: float = MASS) -> dict:
    """Best L/D of the parabolic polar: C_L* = sqrt(pi e AR C_D0), (L/D)max = 0.5 sqrt(pi e AR / C_D0),
    at V* = sqrt(2 W / (rho S C_L*))."""
    raise NotImplementedError


def cruise_drag(V: float, CD0: float, e: float, AR: float = SPAN ** 2 / S_REF, mass: float = MASS) -> dict:
    """Level-flight C_L, C_D, drag (N) and drag power (W) at speed V."""
    raise NotImplementedError


# ---------------------------------------------------------------- data-based investigation

def ladson_cd_near_zero_lift(rows: list[dict], table: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """For each Reynolds-number block of a Ladson table, the cd at the angle of attack closest to 0:
    returns (Re, alpha_deg, cd) sorted by Re."""
    raise NotImplementedError


def effective_transition(Re: float, cd_measured: float, t_c: float = 0.12) -> float:
    """Transition location x_tr (fraction of chord, same on both surfaces) at which section_cd matches a
    measured cd; NaN if no x_tr in [0, 1) can (the measurement lies outside the estimate's range)."""
    raise NotImplementedError


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save drag_build_up.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    AR = SPAN ** 2 / S_REF
    print(f"flat plate at Re = 5e5: laminar Cf {float(cf_laminar(5e5)):.5f}, turbulent {float(cf_turbulent(5e5)):.5f}; "
          f"at Re = 4.52e6 laminar {float(cf_laminar(4.52e6)):.5f}, turbulent {float(cf_turbulent(4.52e6)):.5f}")
    print(f"perimeter / (2 c): NACA 2412 {perimeter_ratio('2412'):.4f}, NACA 0010 {perimeter_ratio('0010'):.4f}")
    comps = trainer_components()
    print("component build-up at 55 m/s:")
    print("  component         Re        Cf       FF      Q     S_wet m^2   D/q m^2   Delta C_D")
    for c in comps:
        print(f"  {c['name']:16s} {c['Re']:9.3g}  {c['Cf']:.5f}  {c['FF']:.4f}  {c['Q']:.2f}  {c['S_wet']:8.3f}  "
              f"{c['D_over_q']:8.4f}  {c['CD']:.5f}")
    cd0 = cd0_build_up(comps)
    print(f"  components {cd0['components']:.5f} + gear {cd0['gear']:.5f} (D/q = {GEAR_DQ} m^2) + cooling "
          f"{cd0['cooling']:.5f} (D/q = {COOLING_DQ} m^2) + misc "
          f"{cd0['misc']:.5f} ({MISC_FRACTION:.0%}) = C_D0 {cd0['CD0']:.5f}; equivalent flat-plate area "
          f"{cd0['CD0'] * S_REF:.3f} m^2")
    for c in comps:
        print(f"    share {c['name']}: {c['CD'] * (1 + MISC_FRACTION) / cd0['CD0']:.1%}")
    print(f"    share gear: {cd0['gear'] * (1 + MISC_FRACTION) / cd0['CD0']:.1%}; cooling "
          f"{cd0['cooling'] * (1 + MISC_FRACTION) / cd0['CD0']:.1%}; misc {cd0['misc'] / cd0['CD0']:.1%}")
    e_ll = lifting_line_e(AR)
    e_inv = oswald_e(cd0["CD0"], k_visc=0.0)
    e = oswald_e(cd0["CD0"])
    print(f"span efficiency: lifting line (rectangular, AR {AR:.3f}) {e_ll:.4f}; with fuselage factor {e_inv:.4f}; "
          f"Oswald e with lift-dependent profile drag {e:.4f}")
    bi = best_ld(cd0["CD0"], e_inv)
    print(f"  (with the inviscid e: best L/D {bi['LD']:.2f} at {bi['V']:.2f} m/s)")
    b = best_ld(cd0["CD0"], e)
    print(f"best L/D = {b['LD']:.2f} at C_L = {b['CL']:.3f}, V = {b['V']:.2f} m/s")
    cr = cruise_drag(V_CRUISE, cd0["CD0"], e)
    print(f"cruise 55 m/s: C_L = {cr['CL']:.4f}, C_D = {cr['CD']:.5f} (induced {cr['CDi']:.5f}), L/D = {cr['CL'] / cr['CD']:.2f}, "
          f"drag {cr['D']:.1f} N, drag power {cr['P'] / 1e3:.2f} kW")
    ap = cruise_drag(30.0, cd0["CD0"], e)
    print(f"approach 30 m/s: C_L = {ap['CL']:.3f}, drag {ap['D']:.1f} N, L/D = {ap['CL'] / ap['CD']:.2f}")
    print("sensitivity: wing transition location and gear drag area")
    for xt in (0.0, 0.10, 0.30):
        c2 = cd0_build_up(trainer_components(x_tr_wing=xt))
        print(f"  wing x_tr = {xt:.2f}: C_D0 = {c2['CD0']:.5f}, (L/D)max = {best_ld(c2['CD0'], oswald_e(c2['CD0']))['LD']:.2f}")
    for g in (0.0, 0.15, 0.30):
        c2 = cd0_build_up(comps, gear_dq=g)
        print(f"  gear D/q = {g:.2f} m^2: C_D0 = {c2['CD0']:.5f}, (L/D)max = {best_ld(c2['CD0'], oswald_e(c2['CD0']))['LD']:.2f}")

    # data-based investigation
    rows = read_ladson()
    ReF, aF, cdF = ladson_cd_near_zero_lift(rows, "I")
    ReG, aG, cdG = ladson_cd_near_zero_lift(rows, "VII")
    print("NACA 0012 near zero lift: estimate 2 Cf FF (FF = %.4f) vs Ladson NASA TM-4074" % form_factor_section(0.12))
    print("  free transition (Table I):")
    for Re, a, cd in zip(ReF, aF, cdF):
        turb = float(section_cd(Re, 0.12, 0.0))
        xt = effective_transition(Re, cd)
        print(f"    Re = {Re:9.3g} (alpha {a:+.2f}): measured {cd:.5f}; fully turbulent estimate {turb:.5f} "
              f"({(turb - cd) / cd:+.1%}); x_tr that matches: {xt:.3f}")
    print("  wraparound grit strip at 0.05c (Table VII):")
    for Re, a, cd in zip(ReG, aG, cdG):
        turb = float(section_cd(Re, 0.12, 0.0))
        t05 = float(section_cd(Re, 0.12, 0.05))
        print(f"    Re = {Re:9.3g} (alpha {a:+.2f}): measured {cd:.5f}; turbulent from LE {turb:.5f} ({(turb - cd) / cd:+.1%}); "
              f"turbulent from 0.05c {t05:.5f} ({(t05 - cd) / cd:+.1%}); grit increment {cd - t05:.5f}")
    rms_f = np.sqrt(np.mean(((section_cd(ReF, 0.12, 0.0) - cdF) / cdF) ** 2))
    rms_g = np.sqrt(np.mean(((section_cd(ReG, 0.12, 0.0) - cdG) / cdG) ** 2))
    print(f"  rms relative difference of the fully turbulent estimate: free {rms_f:.1%}, grit {rms_g:.1%}")

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    CL = np.linspace(0, 1.4, 100)
    ax[0].plot(polar(CL, cd0["CD0"], e), CL, label=f"C_D0 = {cd0['CD0']:.4f}, e = {e:.3f}")
    ax[0].plot([0, polar(b["CL"], cd0["CD0"], e) * 1.5], [0, b["CL"] * 1.5], "k:", label=f"(L/D)max = {b['LD']:.1f}")
    ax[0].plot(cr["CD"], cr["CL"], "o", label="cruise 55 m/s")
    ax[0].set_xlabel("C_D"); ax[0].set_ylabel("C_L"); ax[0].set_title("reference trainer drag polar (estimate)")
    ax[0].legend(); ax[0].grid(alpha=.3)
    Res = np.logspace(6.2, 7.4, 100)
    ax[1].semilogx(Res, section_cd(Res, 0.12, 0.0), label="2 Cf FF, turbulent from LE")
    ax[1].semilogx(Res, section_cd(Res, 0.12, 0.3), "--", label="2 Cf FF, transition at 0.3c")
    ax[1].semilogx(ReF, cdF, "ks", label="Ladson free transition (data)")
    ax[1].semilogx(ReG, cdG, "k^", mfc="none", label="Ladson wraparound grit (data)")
    ax[1].set_xlabel("Re"); ax[1].set_ylabel("cd near zero lift"); ax[1].set_title("NACA 0012: estimate vs NASA TM-4074")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=.3, which="both")
    fig.tight_layout()
    fig.savefig(out_dir / "drag_build_up.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
