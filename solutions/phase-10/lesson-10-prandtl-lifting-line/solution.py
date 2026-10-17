"""P10.10 Prandtl Lifting-Line Theory: reference solution.

What lift-curve slope, span efficiency and induced drag does the reference trainer's rectangular
wing have, and how do taper, aspect ratio and twist change them? This module solves Prandtl's
lifting-line equation by Glauert's Fourier-series collocation:

* the span is y = -(b/2) cos(theta), and the circulation a sine series
  Gamma(theta) = 2 b V sum A_n sin(n theta);
* the monoplane equation at collocation stations theta_i,
  sum_n A_n sin(n theta_i) (mu_i n + sin theta_i) = mu_i (alpha_i - alpha_L0,i) sin theta_i,
  mu = a0 c / (4 b), with alpha_i = alpha + twist(y), is a square linear system for the A_n;
* CL = pi AR A_1, CDi = pi AR sum n A_n^2 = CL^2 (1 + delta)/(pi AR), e = 1/(1 + delta),
  delta = sum_{n>=2} n (A_n/A_1)^2, and the lift-curve slope defines Glauert's tau:
  CL_alpha = a0 / (1 + (a0/(pi AR))(1 + tau)).

Checks: an elliptic chord gives A_n = 0 for n > 1, e = 1 and CL = a0 alpha/(1 + a0/(pi AR)) exactly;
collocation convergence with the observed order; a symmetric wing has no even coefficients and a
symmetric span load.

Conventions (P10 brief): y to starboard, alpha and twist positive nose-up (washout = negative tip
twist), all angles in radians inside the functions, CL and CDi on the reference area S = int c dy.

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes lifting_line.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-10/lesson-10-prandtl-lifting-line
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np

RHO = 1.225                       # kg/m^3, ISA sea level
G = 9.81                          # m/s^2
TRAINER = dict(mass=750.0, S=12.0, b=10.0, c=1.2, V_cruise=55.0)
A0_THIN = 2 * np.pi               # thin-airfoil lift slope per rad (P10.04)
ALPHA_L0_2412 = np.radians(-2.0772)   # thin-airfoil zero-lift angle of the NACA 2412 mean line (P10.04, P10.07)


# =============================================================== given helpers

@dataclass
class Wing:
    """Given: a straight (unswept) wing for lifting-line analysis. chord, twist and alpha_L0 are
    functions of the spanwise station y in [-b/2, b/2] (twist and alpha_L0 in radians); a0 is the
    section lift slope per radian (a number, or a function of y)."""
    b: float
    chord: Callable[[np.ndarray], np.ndarray]
    twist: Callable[[np.ndarray], np.ndarray] = field(default=lambda y: np.zeros_like(np.asarray(y, dtype=float)))
    alpha_L0: Callable[[np.ndarray], np.ndarray] = field(default=lambda y: np.zeros_like(np.asarray(y, dtype=float)))
    a0: float | Callable[[np.ndarray], np.ndarray] = A0_THIN

    def area(self, n: int = 20001) -> float:
        """Reference area S = int c dy (Simpson-free: trapezoid on a fine cosine grid)."""
        th = np.linspace(0.0, np.pi, n)
        y = -0.5 * self.b * np.cos(th)
        return float(np.trapezoid(self.chord(y) * 0.5 * self.b * np.sin(th), th))

    def aspect_ratio(self) -> float:
        """AR = b^2 / S."""
        return self.b ** 2 / self.area()

    def a0_at(self, y: np.ndarray) -> np.ndarray:
        """Section lift slope at y."""
        return self.a0(y) if callable(self.a0) else np.full_like(np.asarray(y, dtype=float), self.a0)


def rectangular_wing(b: float, c: float, **kw) -> Wing:
    """Given: constant chord c."""
    return Wing(b, lambda y: np.full_like(np.asarray(y, dtype=float), c), **kw)


def tapered_wing(b: float, S: float, taper: float, **kw) -> Wing:
    """Given: straight-tapered planform of area S and taper ratio lambda = c_tip/c_root:
    c(y) = c_r (1 - (1 - lambda) |2y/b|), c_r = 2 S/(b (1 + lambda))."""
    cr = 2 * S / (b * (1 + taper))
    return Wing(b, lambda y: cr * (1 - (1 - taper) * np.abs(2 * np.asarray(y) / b)), **kw)


def elliptic_wing(b: float, S: float, **kw) -> Wing:
    """Given: elliptic planform of area S: c(y) = c0 sqrt(1 - (2y/b)^2), c0 = 4 S/(pi b)."""
    c0 = 4 * S / (np.pi * b)
    return Wing(b, lambda y: c0 * np.sqrt(np.clip(1 - (2 * np.asarray(y) / b) ** 2, 0, None)), **kw)


def linear_twist(b: float, tip_twist: float) -> Callable[[np.ndarray], np.ndarray]:
    """Given: twist growing linearly from 0 at the root to tip_twist (rad) at the tips."""
    return lambda y: tip_twist * np.abs(2 * np.asarray(y) / b)


# =============================================================== the solver (learner)

def collocation_points(N: int, symmetric: bool) -> np.ndarray:
    """Collocation angles. Symmetric case (odd n only, half-span): theta_i = i pi/(2N), i = 1..N, the
    last one at the root. General case (all n): theta_i = i pi/(N + 1), i = 1..N."""
    i = np.arange(1, N + 1)
    return i * np.pi / (2 * N) if symmetric else i * np.pi / (N + 1)


def solve_monoplane(wing: Wing, alpha: float, N: int = 40, symmetric: bool = True) -> dict:
    """Glauert's collocation solution of the monoplane equation. Unknowns A_n for n = 1, 3, 5, ...
    (N of them) if symmetric, else n = 1..N. Returns dict(n, A, theta (collocation))."""
    n = 2 * np.arange(N) + 1 if symmetric else np.arange(1, N + 1)
    th = collocation_points(N, symmetric)
    y = -0.5 * wing.b * np.cos(th)
    mu = wing.a0_at(y) * wing.chord(y) / (4 * wing.b)
    S = np.sin(np.outer(th, n))
    M = S * (mu[:, None] * n[None, :] + np.sin(th)[:, None])
    rhs = mu * (alpha + wing.twist(y) - wing.alpha_L0(y)) * np.sin(th)
    return dict(n=n, A=np.linalg.solve(M, rhs), theta=th)


def wing_coefficients(wing: Wing, sol: dict) -> dict:
    """CL = pi AR A_1, CDi = pi AR sum n A_n^2, delta = sum_{n>1} n (A_n/A_1)^2, e = 1/(1 + delta)."""
    AR = wing.aspect_ratio()
    n, A = sol["n"], sol["A"]
    A1 = A[n == 1][0]
    CL = np.pi * AR * A1
    CDi = np.pi * AR * np.sum(n * A ** 2)
    delta = float(np.sum(n[n > 1] * (A[n > 1] / A1) ** 2)) if A1 != 0 else np.nan
    return dict(CL=float(CL), CDi=float(CDi), delta=delta, e=1 / (1 + delta), AR=AR)


def span_loading(wing: Wing, sol: dict, y: np.ndarray, V: float = 1.0) -> dict:
    """Circulation Gamma(y) = 2 b V sum A_n sin(n theta), local lift coefficient cl = 2 Gamma/(V c),
    and induced angle alpha_i = sum n A_n sin(n theta)/sin(theta) at the stations y (|y| < b/2)."""
    th = np.arccos(np.clip(-2 * np.asarray(y, dtype=float) / wing.b, -1, 1))
    S = np.sin(np.outer(th, sol["n"]))
    Gam = 2 * wing.b * V * S @ sol["A"]
    with np.errstate(divide="ignore", invalid="ignore"):
        ai = (S * sol["n"][None, :]) @ sol["A"] / np.sin(th)
        cl = 2 * Gam / (V * wing.chord(y))
    return dict(gamma=Gam, cl=cl, alpha_i=ai)


def rolling_moment_coefficient(wing: Wing, sol: dict) -> float:
    """Rolling-moment coefficient Cl = (pi AR/4) A_2 (positive right wing down) from an
    unsymmetric solution: only the n = 2 term carries a net moment of lift about the x axis."""
    n, A = sol["n"], sol["A"]
    return float(np.pi * wing.aspect_ratio() / 4 * A[n == 2][0]) if np.any(n == 2) else 0.0


def lift_slope_and_tau(wing: Wing, N: int = 40) -> dict:
    """CL_alpha (per rad) from two solutions (the problem is linear in alpha) and Glauert's tau
    from CL_alpha = a0/(1 + (a0/(pi AR))(1 + tau)), with a0 the root section's slope."""
    c1 = wing_coefficients(wing, solve_monoplane(wing, 0.1, N))["CL"]
    c0 = wing_coefficients(wing, solve_monoplane(wing, 0.0, N))["CL"]
    CLa = (c1 - c0) / 0.1
    a0 = float(wing.a0_at(np.array([0.0]))[0])
    AR = wing.aspect_ratio()
    tau = (a0 / CLa - 1) * np.pi * AR / a0 - 1
    return dict(CL_alpha=float(CLa), tau=float(tau))


def elliptic_closed_form(a0: float, AR: float, alpha: float) -> dict:
    """Flat elliptic wing: CL = a0 alpha/(1 + a0/(pi AR)), CDi = CL^2/(pi AR), e = 1."""
    CL = a0 * alpha / (1 + a0 / (np.pi * AR))
    return dict(CL=float(CL), CDi=float(CL ** 2 / (np.pi * AR)), e=1.0)


def alpha_for_cl(wing: Wing, CL_target: float, N: int = 40) -> float:
    """Angle of attack (rad) for a target CL: the solution is affine in alpha."""
    c0 = wing_coefficients(wing, solve_monoplane(wing, 0.0, N))["CL"]
    c1 = wing_coefficients(wing, solve_monoplane(wing, 0.1, N))["CL"]
    return float((CL_target - c0) / ((c1 - c0) / 0.1))


def convergence_study(wing: Wing, alpha: float, Ns: list[int], N_ref: int = 1280) -> dict:
    """Errors of CL and CDi against a fine reference solution (N_ref terms), and the observed
    orders log2(e_k / e_{k+1}) for successive doublings of N."""
    ref = wing_coefficients(wing, solve_monoplane(wing, alpha, N_ref))
    eCL, eCD = [], []
    for N in Ns:
        r = wing_coefficients(wing, solve_monoplane(wing, alpha, N))
        eCL.append(abs(r["CL"] - ref["CL"])); eCD.append(abs(r["CDi"] - ref["CDi"]))
    oCL = np.log(np.array(eCL[:-1]) / np.array(eCL[1:])) / np.log(np.array(Ns[1:]) / np.array(Ns[:-1]))
    oCD = np.log(np.array(eCD[:-1]) / np.array(eCD[1:])) / np.log(np.array(Ns[1:]) / np.array(Ns[:-1]))
    return dict(err_CL=np.array(eCL), err_CDi=np.array(eCD), order_CL=oCL, order_CDi=oCD, ref=ref)


# =============================================================== main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save lifting_line.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    m_, S, b, c, V = TRAINER["mass"], TRAINER["S"], TRAINER["b"], TRAINER["c"], TRAINER["V_cruise"]

    # 1. elliptic check
    ew = elliptic_wing(b, S)
    a = np.radians(5.0)
    r = wing_coefficients(ew, solve_monoplane(ew, a, 10))
    cf = elliptic_closed_form(2 * np.pi, ew.aspect_ratio(), a)
    print(f"elliptic wing AR {ew.aspect_ratio():.4f}, alpha 5 deg: CL {r['CL']:.8f} (closed form {cf['CL']:.8f}), "
          f"e {r['e']:.10f}, max |A_n>1| = {np.max(np.abs(solve_monoplane(ew, a, 10)['A'][1:])):.1e}")

    # 2. the trainer
    tw = rectangular_wing(b, c, alpha_L0=lambda y: np.full_like(np.asarray(y, dtype=float), ALPHA_L0_2412))
    W = m_ * G
    CLc = W / (0.5 * RHO * V ** 2 * S)
    ls = lift_slope_and_tau(tw)
    ac = alpha_for_cl(tw, CLc)
    rc = wing_coefficients(tw, solve_monoplane(tw, ac))
    q = 0.5 * RHO * V ** 2
    print(f"\nreference trainer (rectangular, b {b} m, c {c} m, AR {tw.aspect_ratio():.3f}, a0 = 2 pi, alpha_L0 = -2.077 deg):")
    print(f"  CL_alpha = {ls['CL_alpha']:.4f} /rad = {ls['CL_alpha'] * np.pi / 180:.5f} /deg; elliptic would give "
          f"{2 * np.pi / (1 + 2 / tw.aspect_ratio()):.4f}; tau = {ls['tau']:.4f}")
    print(f"  cruise CL {CLc:.4f} at alpha = {np.degrees(ac):.3f} deg; e = {rc['e']:.4f} (delta {rc['delta']:.4f}); "
          f"CDi = {rc['CDi']:.6f} (elliptic {CLc ** 2 / (np.pi * tw.aspect_ratio()):.6f}); induced drag {rc['CDi'] * q * S:.1f} N")
    yy = np.array([0.0, 2.5, 4.0, 4.5, 4.9])
    sl = span_loading(tw, solve_monoplane(tw, ac), yy, V)
    print("  y (m)    Gamma (m^2/s)   cl/CL    alpha_i (deg)")
    for y, g_, cl_, ai in zip(yy, sl["gamma"], sl["cl"], sl["alpha_i"]):
        print(f"  {y:4.1f}     {g_:8.3f}       {cl_ / CLc:.4f}    {np.degrees(ai):.3f}")
    for k in (0.9, 1.0):
        twk = rectangular_wing(b, c, a0=k * 2 * np.pi)
        lk = lift_slope_and_tau(twk)
        print(f"  sensitivity: a0 = {k:.1f} x 2 pi -> CL_alpha {lk['CL_alpha']:.4f} /rad, "
              f"e {wing_coefficients(twk, solve_monoplane(twk, 0.1))['e']:.4f}")

    # 3. convergence
    print("\ncollocation convergence (alpha = 5 deg, flat wings, reference N = 1280):")
    Ns = [5, 10, 20, 40, 80, 160]
    for name, wg in (("rectangular AR 8.33", rectangular_wing(b, c)), ("taper 0.4, AR 8.33", tapered_wing(b, S, 0.4)),
                     ("elliptic + linear washout -4 deg", elliptic_wing(b, S, twist=linear_twist(b, np.radians(-4.0))))):
        cs = convergence_study(wg, a, Ns)
        print(f"  {name}: |CL error| {np.array2string(cs['err_CL'], precision=2)}")
        print(f"      observed orders CL {np.round(cs['order_CL'], 2)}, CDi {np.round(cs['order_CDi'], 2)}")

    # 4. delta and tau against taper and AR
    print("\ndelta (induced-drag factor) and tau (lift-slope factor) of flat, untwisted tapered wings, a0 = 2 pi:")
    tapers = [0.0, 0.2, 0.3, 0.35, 0.4, 0.5, 0.7, 1.0]
    print("  AR    " + "  ".join(f"l={t:4.2f}" for t in tapers))
    table = {}
    for AR in (4.0, 6.0, 8.33, 12.0):
        Sx = b ** 2 / AR
        row_d, row_t = [], []
        for t in tapers:
            wg = tapered_wing(b, Sx, t)
            row_d.append(wing_coefficients(wg, solve_monoplane(wg, 0.1, 80))["delta"])
            row_t.append(lift_slope_and_tau(wg, 80)["tau"])
        table[AR] = (row_d, row_t)
        print(f"  {AR:5.2f} d " + "  ".join(f"{d:6.4f}" for d in row_d))
        print(f"        t " + "  ".join(f"{t_:6.4f}" for t_ in row_t))
    dd = table[8.33][0]
    tf = np.linspace(0.0, 1.0, 101)
    dl = [wing_coefficients(tapered_wing(b, S, t), solve_monoplane(tapered_wing(b, S, t), 0.1, 80))["delta"] for t in tf]
    print(f"  AR 8.33: smallest delta {min(dl):.5f} at taper {tf[int(np.argmin(dl))]:.2f} (e = {1 / (1 + min(dl)):.4f})")

    # 5. symmetry and an asymmetric case (aileron-like twist)
    gw = rectangular_wing(b, c, twist=lambda y: np.radians(2.0) * np.sign(y) * (np.abs(y) > 3.0))
    sg = solve_monoplane(gw, a, 40, symmetric=False)
    print(f"\nanti-symmetric 'aileron' twist (+-2 deg outboard of 3 m) on the flat wing: A_2 = {sg['A'][1]:.5f}, "
          f"rolling-moment coefficient Cl = (pi AR/4) A_2 = {np.pi * gw.aspect_ratio() / 4 * sg['A'][1]:.5f} "
          f"(negative: right wing up), CL unchanged {wing_coefficients(gw, sg)['CL']:.5f} vs "
          f"{wing_coefficients(gw, solve_monoplane(rectangular_wing(b, c), a, 40))['CL']:.5f}")

    # figure
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    y = np.linspace(-b / 2 * 0.999, b / 2 * 0.999, 400)
    for name, wg in (("rectangular", rectangular_wing(b, c)), ("taper 0.4", tapered_wing(b, S, 0.4)), ("elliptic", ew)):
        aa = alpha_for_cl(wg, CLc)
        ax[0].plot(y, span_loading(wg, solve_monoplane(wg, aa), y, V)["cl"] / CLc, label=name)
    ax[0].set_xlabel("y (m)"); ax[0].set_ylabel("cl / CL"); ax[0].set_title("local lift coefficient at cruise CL"); ax[0].legend()
    for name, wg in (("rectangular", rectangular_wing(b, c)), ("taper 0.4", tapered_wing(b, S, 0.4))):
        cs = convergence_study(wg, a, Ns)
        ax[1].loglog(Ns, cs["err_CL"], "o-", label=f"CL, {name}")
    ax[1].set_xlabel("number of Fourier terms N"); ax[1].set_ylabel("|CL - CL_ref|"); ax[1].set_title("convergence"); ax[1].legend()
    ax[2].plot(tf, dl)
    ax[2].set_xlabel("taper ratio"); ax[2].set_ylabel("delta"); ax[2].set_title("induced-drag factor, AR 8.33")
    fig.tight_layout()
    fig.savefig(out_dir / "lifting_line.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
