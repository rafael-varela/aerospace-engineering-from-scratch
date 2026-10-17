"""P10.15 Stall & Post-Stall: starter code.

Where along the span does the reference trainer's wing stall first, and what is its maximum lift
coefficient, if the wing is rectangular, tapered, or twisted? The module answers with a section model
fitted to authentic data and a nonlinear lifting line:

* the Kirchhoff/Helmholtz separation-point model cl = a0 alpha ((1 + sqrt f)/2)^2, where f is the
  fraction of the chord over which the flow is attached (f = 1: thin-airfoil lift a0 alpha; f = 0:
  fully separated, a quarter of it), with f(alpha) = 1 / (1 + exp((|alpha| - alpha1) / S)) a sigmoid
  whose three parameters (a0, alpha1, S) are fitted by least squares to the authentic NACA 0012 lift
  curves of Ladson (NASA TM-4074, free transition, M = 0.15);
* a flat-plate post-stall limit for large angles: the normal force of a plate broadside to the flow is
  modelled as c_n = 2 sin(alpha) (a stated empirical approximation, close to the measured c_n of about
  2 for a long flat plate at 90 deg), so cl = c_n cos(alpha) = sin(2 alpha) and cd = 2 sin^2(alpha);
  the full-range section curve blends the fitted Kirchhoff curve into it above about 25 deg;
* a nonlinear lifting line: Prandtl's Fourier series for the circulation, with the section curve
  replacing the straight line, solved by relaxed fixed-point iteration on the Fourier coefficients;
  it reduces to the linear lifting line below stall.

Angles in radians inside the functions, degrees in printed output. Units SI.

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes stall.png)
Test: pytest curriculum/phase-10/lesson-15-stall-post-stall
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.optimize import least_squares

RHO = 1.225
NU = 1.461e-5
SPAN = 10.0
S_REF = 12.0
CHORD = 1.2


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


def read_ladson_block(table: str = "I", reynolds: float = 2.0e6, path: str | Path | None = None) -> dict:
    """alpha (deg), cl, cd, cm arrays of one Reynolds-number block of the authentic Ladson NACA 0012
    data (NASA TM-4074, M = 0.15); cd is NaN where the report prints none (post-stall rows)."""
    path = Path(path) if path else find_data_file("ladson_naca0012_M015.csv")
    with open(path) as fh:
        rows = [r for r in csv.DictReader(line for line in fh if not line.startswith("#"))
                if r["table"] == table and abs(float(r["reynolds"]) - reynolds) < 1e3]
    if not rows:
        raise ValueError(f"no block table={table} Re={reynolds}")
    f = lambda v: float(v) if v.strip() else float("nan")
    return {k: np.array([f(r[src]) for r in rows]) for k, src in
            (("alpha", "alpha_deg"), ("cl", "cl"), ("cd", "cd"), ("cm", "cm"))}


def ladson_reynolds_numbers(table: str = "I", path: str | Path | None = None) -> list[float]:
    """Reynolds numbers of the blocks of a Ladson table, ascending."""
    path = Path(path) if path else find_data_file("ladson_naca0012_M015.csv")
    with open(path) as fh:
        return sorted({float(r["reynolds"]) for r in csv.DictReader(line for line in fh if not line.startswith("#"))
                       if r["table"] == table})


# ---------------------------------------------------------------- given helpers: planforms

def rectangular_wing(n: int = 1) -> tuple[Callable, Callable]:
    """The trainer's rectangular wing: chord(y) = 1.2 m, no twist. Returns (chord(y), twist(y)) for
    y in [-b/2, b/2]."""
    return (lambda y: np.full_like(np.asarray(y, dtype=float), CHORD),
            lambda y: np.zeros_like(np.asarray(y, dtype=float)))


def tapered_wing(taper: float = 0.45) -> tuple[Callable, Callable]:
    """Same span and area, linear taper (tip chord / root chord = taper), no twist."""
    c_root = 2 * S_REF / (SPAN * (1 + taper))
    return (lambda y: c_root * (1 - (1 - taper) * np.abs(np.asarray(y, dtype=float)) / (SPAN / 2)),
            lambda y: np.zeros_like(np.asarray(y, dtype=float)))


def twisted_wing(washout_deg: float = 4.0) -> tuple[Callable, Callable]:
    """Rectangular wing with linear washout: the tip is washout_deg nose-down relative to the root."""
    w = np.radians(washout_deg)
    return (lambda y: np.full_like(np.asarray(y, dtype=float), CHORD),
            lambda y: -w * np.abs(np.asarray(y, dtype=float)) / (SPAN / 2))


# ---------------------------------------------------------------- section model

def kirchhoff_cl(alpha: np.ndarray, f: np.ndarray, a0: float = 2 * np.pi) -> np.ndarray:
    """Kirchhoff/Helmholtz separated-flow lift: cl = a0 alpha ((1 + sqrt f) / 2)^2, 0 <= f <= 1."""
    raise NotImplementedError


def separation_point(alpha: np.ndarray, alpha1: float, S: float) -> np.ndarray:
    """Attached-flow fraction f(alpha) = 1 / (1 + exp((|alpha| - alpha1) / S)): near 1 at small angles,
    1/2 at alpha1, tending to 0 beyond it over a width set by S (radians)."""
    raise NotImplementedError


def fit_kirchhoff(alpha_deg: np.ndarray, cl: np.ndarray) -> dict:
    """Least-squares fit of (a0, alpha1, S) in cl = kirchhoff_cl(alpha, separation_point(alpha, alpha1, S), a0)
    to measured (alpha, cl) points. Returns the parameters (alpha1 and S in radians), the residuals,
    the rms residual and the parameter standard errors from the Jacobian."""
    raise NotImplementedError


def flat_plate(alpha: np.ndarray, cd0: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """Post-stall flat-plate limit with normal force c_n = 2 sin(alpha) (stated empirical model):
    cl = c_n cos(alpha) = sin(2 alpha), cd = cd0 + c_n sin(alpha) = cd0 + 2 sin^2(alpha)."""
    raise NotImplementedError


def full_range_cl(alpha: np.ndarray, fit: dict, alpha_blend: float = np.radians(25.0),
                  width: float = np.radians(5.0)) -> np.ndarray:
    """Section lift over the full range: the fitted Kirchhoff curve blended into the flat-plate limit,
    cl = (1 - w) cl_K + w cl_plate, with a smoothstep weight w(t) = 3t^2 - 2t^3,
    t = (|alpha| - alpha_blend + width) / (2 width) clipped to [0, 1]: exactly the Kirchhoff curve below
    alpha_blend - width (20 deg) and exactly the flat plate above alpha_blend + width (30 deg)."""
    raise NotImplementedError


# ---------------------------------------------------------------- lifting line

def lifting_line_linear(alpha: float, chord: Callable, twist: Callable, a0: float, b: float = SPAN,
                        alpha_L0: float = 0.0, n_terms: int = 16) -> dict:
    """Prandtl's lifting line with a linear section cl = a0 (alpha_eff - alpha_L0), symmetric wing:
    collocation of sum_n A_n sin(n theta) (n mu + sin theta) = mu (alpha + twist - alpha_L0) sin theta,
    mu = c a0 / (4 b), odd n only, theta in (0, pi/2]. Returns A, theta, y, CL, CDi, e and cl(y)."""
    raise NotImplementedError


def _ll_grid(chord: Callable, b: float, n_terms: int):
    """Collocation angles on the half span, odd Fourier indices, stations y = (b/2) cos(theta), chords."""
    th = np.linspace(np.pi / (2 * n_terms), np.pi / 2, n_terms)
    n = 2 * np.arange(n_terms) + 1
    y = b / 2 * np.cos(th)
    return th, n, y, chord(y)


def _ll_out(A, th, n, y, c, b):
    """Lifting-line results from the Fourier coefficients: CL = pi AR A1, CDi = pi AR sum n A_n^2,
    e, induced angle and local cl at the collocation stations."""
    AR = b ** 2 / S_REF
    CL = np.pi * AR * A[0]
    delta = float(np.sum(n[1:] * (A[1:] / A[0]) ** 2)) if A[0] != 0 else 0.0
    CDi = np.pi * AR * float(np.sum(n * A ** 2))
    alpha_i = np.sin(np.outer(th, n)) @ (n * A) / np.sin(th)
    gamma_ratio = 2 * b * np.sin(np.outer(th, n)) @ A           # Gamma / V
    cl_local = 2 * gamma_ratio / c
    return {"A": A, "theta": th, "y": y, "c": c, "CL": CL, "CDi": CDi, "e": 1 / (1 + delta),
            "alpha_i": alpha_i, "cl": cl_local}


def lifting_line_nonlinear(alpha: float, chord: Callable, twist: Callable, section_cl: Callable,
                           a0: float = 2 * np.pi, b: float = SPAN, n_terms: int = 16, relax: float = 0.5,
                           tol: float = 1e-10, max_iter: int = 500, delta0: np.ndarray | None = None) -> dict:
    """Nonlinear lifting line: the section lift is cl = section_cl(alpha_eff), with
    alpha_eff = alpha + twist(y) - alpha_i(y), instead of a straight line.

    Iteration with relaxation on a local zero-lift-angle shift ("decambering"): each sweep solves the
    LINEAR lifting line with slope a0 and a spanwise zero-lift angle delta(y); that gives alpha_eff(y).
    The section curve says cl = section_cl(alpha_eff); the shift that would make the straight line pass
    through that point is delta_new = alpha_eff - section_cl(alpha_eff) / a0. Then
    delta <- delta + relax (delta_new - delta). At convergence the linear and nonlinear lift agree at
    every station, so the circulation satisfies the lifting-line equation with the true section curve.
    Solving the linear system every sweep keeps the stiff high harmonics implicit; a plain fixed point
    on the Fourier coefficients needs a relaxation factor below about 2 / (1 + a0 n_max / (4 AR)).
    Returns the lifting-line results plus 'alpha_eff', 'cl_section', 'delta', 'iterations', 'converged'
    and the residual history max |cl_linear - cl_section|.
    """
    raise NotImplementedError


def fixed_point_relaxation_limit(chord: Callable, a_slope: float, b: float = SPAN, n_terms: int = 16) -> float:
    """Largest stable relaxation factor of the Fourier fixed point for a section of constant slope a_slope:
    its iteration matrix J (A_new = J A + const) has real negative eigenvalues, and
    A <- A + w (J A - A) converges iff w < 2 / (1 - lambda_min(J))."""
    raise NotImplementedError


def lifting_line_fixed_point(alpha: float, chord: Callable, twist: Callable, section_cl: Callable,
                             b: float = SPAN, n_terms: int = 16, relax: float = 0.05, tol: float = 1e-10,
                             max_iter: int = 20000) -> dict:
    """The textbook iteration for comparison: relax the Fourier coefficients directly,
    A <- A + relax (A_new - A), with A_new from Gamma/V = c section_cl(alpha_eff) / 2. Explicit in the
    induced angle, so it is stable only for small relaxation factors."""
    raise NotImplementedError


def stall_onset(res: dict, cl_max: float) -> dict:
    """Where the wing is closest to stalling: the span station with the largest ratio of local cl to the
    section cl_max, the ratio and the spanwise position as a fraction of the semi-span 2|y|/b."""
    raise NotImplementedError


def section_stall_angle(fit: dict) -> tuple[float, float]:
    """Angle (rad) and value of the maximum of the fitted Kirchhoff section curve (dense search)."""
    raise NotImplementedError


def stall_onset_alpha(chord: Callable, twist: Callable, section_cl: Callable, a0: float, alpha_stall: float,
                      lo_deg: float = 0.0, hi_deg: float = 30.0, n_terms: int = 16, tol_deg: float = 1e-3) -> dict:
    """Wing angle of attack at which the first span station reaches the section stall angle
    (alpha_eff = alpha_stall), found by bisection with the nonlinear lifting line; also the station
    (2|y|/b), the wing CL there and the local-cl distribution. This is the 'critical section' estimate of
    the wing's maximum lift: CL at the first local stall."""
    raise NotImplementedError


def wing_cl_max(chord: Callable, twist: Callable, section_cl: Callable, alphas_deg: np.ndarray, a0: float,
                n_terms: int = 16, relax: float = 0.5) -> dict:
    """Sweep alpha upwards with the nonlinear lifting line, each solve started from the previous
    solution (continuation); return the CL curve, convergence flags, the largest converged CL and its
    alpha."""
    raise NotImplementedError


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save stall.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    print("Kirchhoff fits to Ladson NACA 0012, free transition (NASA TM-4074):")
    fits = {}
    for Re in ladson_reynolds_numbers("I"):
        d = read_ladson_block("I", Re)
        ft = fit_kirchhoff(d["alpha"], d["cl"])
        fits[Re] = (d, ft)
        print(f"  Re = {Re:.3g}: a0 = {ft['a0']:.3f} +- {ft['se'][0]:.3f} /rad, alpha1 = {np.degrees(ft['alpha1']):.2f} "
              f"+- {np.degrees(ft['se'][1]):.2f} deg, S = {np.degrees(ft['S']):.3f} deg, rms residual {ft['rms']:.4f}, "
              f"max |residual| {np.max(np.abs(ft['residuals'])):.4f}; model cl_max {ft['cl_max_model']:.3f}, "
              f"measured max cl {np.max(d['cl']):.3f}")
    Re_wing = 30.0 * CHORD / NU
    d, ft = fits[2.0e6]
    print(f"section for the wing: Re = 2.0e6 block (trainer approach Re = {Re_wing:.3g})")
    a = np.radians(np.array([4.0, 8.0, 12.0, 14.0]))
    for aa in a:
        f = separation_point(aa, ft["alpha1"], ft["S"])
        print(f"  alpha = {np.degrees(aa):4.1f} deg: f = {f:.4f}, cl = {float(kirchhoff_cl(aa, f, ft['a0'])):.4f}")
    cl_curve = lambda al: full_range_cl(al, ft)
    print("post-stall: Kirchhoff (f -> 0) and flat plate at large alpha:")
    for ad in (20, 30, 45, 60, 90):
        ar = np.radians(ad)
        clp, cdp = flat_plate(ar)
        print(f"  alpha = {ad:2d} deg: Kirchhoff {float(kirchhoff_cl(ar, separation_point(ar, ft['alpha1'], ft['S']), ft['a0'])):.3f}, "
              f"flat plate cl {float(clp):.3f}, cd {float(cdp):.3f}, blended {float(cl_curve(ar)):.3f}")

    # lifting line checks
    a0 = ft["a0"]
    chord, twist = rectangular_wing()
    lin = lifting_line_linear(np.radians(5.0), chord, twist, a0)
    nl = lifting_line_nonlinear(np.radians(5.0), chord, twist, lambda al: a0 * al, a0=a0)
    print(f"linear section, alpha = 5 deg: linear LL CL = {lin['CL']:.6f}, nonlinear iteration CL = {nl['CL']:.6f} "
          f"({nl['iterations']} iterations, converged {nl['converged']}), e = {lin['e']:.4f}")
    ell = lambda y: 4 * S_REF / (np.pi * SPAN) * np.sqrt(np.clip(1 - (2 * np.asarray(y) / SPAN) ** 2, 0, None))
    le = lifting_line_linear(np.radians(5.0), ell, twist, 2 * np.pi)
    AR = SPAN ** 2 / S_REF
    print(f"elliptic wing, 2 pi section: CL = {le['CL']:.6f}, closed form 2 pi alpha AR/(AR+2) = "
          f"{2 * np.pi * np.radians(5) * AR / (AR + 2):.6f}, e = {le['e']:.6f}")
    print("iteration at alpha = 12 deg (rectangular, fitted section, 16 terms):")
    for rl in (0.3, 0.5, 0.8, 1.0):
        r = lifting_line_nonlinear(np.radians(12.0), chord, twist, cl_curve, a0=a0, relax=rl)
        print(f"  decambering, relaxation {rl}: converged {r['converged']} in {r['iterations']} iterations, CL = {r['CL']:.6f}")
    for rl in (0.02, 0.05, 0.1, 0.2):
        r = lifting_line_fixed_point(np.radians(12.0), chord, twist, cl_curve, relax=rl)
        print(f"  Fourier fixed point, relaxation {rl}: converged {r['converged']} in {r['iterations']} iterations, "
              f"CL = {r['CL']:.6f}")
    print(f"  stability limit of the Fourier fixed point for a straight-line section of slope a0: relax < "
          f"{fixed_point_relaxation_limit(chord, a0):.4f}")
    for a_ in (15.5, 16.0):
        r1 = lifting_line_nonlinear(np.radians(a_), chord, twist, cl_curve, a0=2 * np.pi, relax=0.5)
        r2 = lifting_line_fixed_point(np.radians(a_), chord, twist, cl_curve, relax=0.05)
        print(f"  alpha = {a_} deg: decambering (from zero) CL = {r1['CL']:.6f} ({r1['iterations']} it), Fourier fixed "
              f"point (from zero) CL = {r2['CL']:.6f} ({r2['iterations']} it), both converged {r1['converged'] and r2['converged']}")
    print("Fourier-term convergence of CL at alpha = 12 deg:")
    prev = None
    for nt in (4, 8, 16, 32, 64):
        r = lifting_line_nonlinear(np.radians(12.0), chord, twist, cl_curve, a0=a0, n_terms=nt)
        print(f"  {nt:2d} terms: CL = {r['CL']:.7f}" + (f", change {r['CL'] - prev:+.2e}" if prev is not None else ""))
        prev = r["CL"]

    # stall pattern for three planforms
    a_st, cl_max_sec = section_stall_angle(ft)
    wings = {"rectangular": rectangular_wing(), "tapered 0.45": tapered_wing(0.45), "rectangular, 4 deg washout": twisted_wing(4.0)}
    alphas = np.arange(0.0, 24.01, 0.5)
    results = {}
    print(f"stall pattern (section: cl_max = {cl_max_sec:.3f} at {np.degrees(a_st):.2f} deg):")
    for name, (ch, tw) in wings.items():
        so = stall_onset_alpha(ch, tw, cl_curve, a0, a_st)
        w = wing_cl_max(ch, tw, cl_curve, alphas, a0)
        lin_ = lifting_line_linear(np.radians(5.0), ch, tw, a0)
        lin0 = lifting_line_linear(0.0, ch, tw, a0)
        print(f"  {name:28s}: e = {lin_['e']:.4f}, CL_alpha = {(lin_['CL'] - lin0['CL']) / np.radians(5.0):.3f}/rad; "
              f"first local stall at alpha = {so['alpha_deg']:.2f} deg, 2|y|/b = {so['eta']:.3f}, CL = {so['CL']:.3f} "
              f"(CL/cl_max = {so['CL'] / cl_max_sec:.3f}); largest converged CL in the sweep {w['CL_max']:.3f} at "
              f"{w['alpha_CL_max']:.1f} deg; not converged at {int(np.sum(~w['converged']))} of {len(alphas)} angles")
        results[name] = (w, so)

    # hysteresis probe: sweep up, then down, each solve started from the previous one
    ch, tw = wings["rectangular"]
    sweep = np.arange(12.0, 22.01, 0.5)
    up, dn, dl = [], [], None
    for a_ in sweep:
        r = lifting_line_nonlinear(np.radians(a_), ch, tw, cl_curve, a0=a0, relax=0.3, delta0=dl)
        up.append(r["CL"] if r["converged"] else np.nan)
        dl = r["delta"]
    for a_ in sweep[::-1]:
        r = lifting_line_nonlinear(np.radians(a_), ch, tw, cl_curve, a0=a0, relax=0.3, delta0=dl)
        dn.append(r["CL"] if r["converged"] else np.nan)
        dl = r["delta"]
    up, dn = np.array(up), np.array(dn[::-1])
    for a_, u, v in zip(sweep, up, dn):
        if abs(a_ - round(a_)) < 1e-9 and a_ % 2 == 0:
            print(f"  alpha = {a_:4.1f} deg: CL up {u:.4f}, CL down {v:.4f}")
    print(f"up/down sweep 12-22 deg (rectangular): largest |CL_up - CL_down| = {np.nanmax(np.abs(up - dn)):.3e}; "
          f"non-converged points up {int(np.sum(np.isnan(up)))}, down {int(np.sum(np.isnan(dn)))}")

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
    for Re, (dd, ff) in fits.items():
        p = ax[0].plot(dd["alpha"], dd["cl"], "o", ms=4, label=f"data Re {Re:.2g}")
        aa = np.radians(np.linspace(-5, 20, 300))
        ax[0].plot(np.degrees(aa), kirchhoff_cl(aa, separation_point(aa, ff["alpha1"], ff["S"]), ff["a0"]), color=p[0].get_color(), lw=1)
    ax[0].set_xlabel("alpha (deg)"); ax[0].set_ylabel("cl"); ax[0].set_title("NACA 0012 (Ladson) and Kirchhoff fits")
    ax[0].legend(fontsize=6); ax[0].grid(alpha=.3)
    aa = np.radians(np.linspace(0, 90, 500))
    ax[1].plot(np.degrees(aa), cl_curve(aa), label="fitted Kirchhoff + flat plate")
    ax[1].plot(np.degrees(aa), flat_plate(aa)[0], "--", label="flat plate sin 2 alpha")
    ax[1].plot(d["alpha"], d["cl"], "ko", ms=3, label="data Re 2e6")
    ax[1].set_xlabel("alpha (deg)"); ax[1].set_ylabel("cl"); ax[1].legend(); ax[1].grid(alpha=.3)
    for name, (w, so) in results.items():
        r = so["result"]
        order = np.argsort(r["y"])
        ax[2].plot(2 * r["y"][order] / SPAN, r["cl"][order] / cl_max_sec, label=f"{name}, alpha {so['alpha_deg']:.1f} deg")
    ax[2].axhline(1.0, color="k", lw=0.8)
    ax[2].set_xlabel("2y/b"); ax[2].set_ylabel("local cl / section cl_max"); ax[2].set_title("stall onset")
    ax[2].legend(fontsize=7); ax[2].grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(out_dir / "stall.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
