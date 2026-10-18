"""P10.12 Wing Planform Effects: reference solution.

How should the reference trainer's wing be tapered and twisted, and what does each choice cost in
induced drag, root bending moment and stall behaviour? This module answers with the lifting line
(a compact P10.10 solver) and the vortex lattice (a compact P10.11 solver) as analysis engines and
scipy.optimize as the driver (P03.14, P03.16):

* span efficiency e(lambda) of flat, untwisted straight-tapered wings, computed (its maximum, the
  dip in the induced-drag factor delta, sits near lambda = 0.35);
* the local lift coefficient cl(y)/CL and the station where a wing of constant section cl_max
  first stalls; the effect of sweep on that station (VLM);
* the twist that minimises induced drag at a design CL: unconstrained (it recovers the elliptic
  load, e -> 1), with a linear twist (checked against the Lagrange-multiplier closed form), and
  with constraints on the root bending moment and on the tip's stall margin at a high-lift CL;
* the bending-moment trade and Prandtl's 1933 bell-shaped span load: with lift and the integral
  of lift x y^2 (a proxy for structural weight) fixed and the span free, the optimum span is 22 %
  longer and the induced drag 11 % lower than the elliptic wing's.

Conventions: y to starboard, eta = 2y/b, twist positive nose-up (washout negative), flat (zero
camber) sections with a0 = 2 pi; the camber of a real section only shifts alpha. Root bending
moment coefficient C_B = M_root/(q S b), M_root the moment of one wing's lift about the root.

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes planform.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-10/lesson-12-wing-planform-effects
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
from scipy.optimize import minimize, minimize_scalar

RHO = 1.225
G = 9.81
TRAINER = dict(mass=750.0, S=12.0, b=10.0, c=1.2, V_cruise=55.0, V_approach=30.0)
A0 = 2 * np.pi


# =============================================================== given helpers

def tapered_chord(b: float, S: float, taper: float) -> Callable[[np.ndarray], np.ndarray]:
    """Given: straight-tapered chord of area S, c(y) = c_r (1 - (1 - lambda)|2y/b|)."""
    cr = 2 * S / (b * (1 + taper))
    return lambda y: cr * (1 - (1 - taper) * np.abs(2 * np.asarray(y) / b))


def ll_solve(b: float, chord: Callable, twist: Callable, alpha: float, N: int = 40, a0: float = A0) -> np.ndarray:
    """Given (compact P10.10 solver): symmetric lifting line by Glauert collocation. Returns the odd
    coefficients A_1, A_3, ..., A_{2N-1} of Gamma = 2 b V sum A_n sin(n theta), y = -(b/2) cos(theta)."""
    n = 2 * np.arange(N) + 1
    th = np.arange(1, N + 1) * np.pi / (2 * N)
    y = -0.5 * b * np.cos(th)
    mu = a0 * chord(y) / (4 * b)
    M = np.sin(np.outer(th, n)) * (mu[:, None] * n + np.sin(th)[:, None])
    return np.linalg.solve(M, mu * (alpha + twist(y)) * np.sin(th))


def ll_odd_n(N: int) -> np.ndarray:
    """Given: the odd harmonics 1, 3, ..., 2N - 1 that ll_solve returns coefficients for."""
    return 2 * np.arange(N) + 1


def ll_gamma(b: float, A: np.ndarray, y: np.ndarray, V: float = 1.0) -> np.ndarray:
    """Given: Gamma(y) = 2 b V sum A_n sin(n theta) for the odd-harmonic coefficients A."""
    th = np.arccos(np.clip(-2 * np.asarray(y, dtype=float) / b, -1, 1))
    return 2 * b * V * np.sin(np.outer(th, ll_odd_n(len(A)))) @ A


def vlm_strip_cl(b: float, c_root: float, taper: float, sweep_c4: float, alpha: float = np.radians(5.0),
                 nx: int = 4, ny: int = 40) -> tuple[np.ndarray, np.ndarray, float]:
    """Given (compact P10.11 solver): flat VLM with Hough spanwise spacing, 1/4-3/4 rule, trailing
    legs along x. Returns (y_strip, cl_strip, CL)."""
    d = b / (ny + 0.5)
    ye = np.linspace(-0.5 * b + 0.25 * d, 0.5 * b - 0.25 * d, ny + 1)
    ch = lambda y: c_root * (1 - (1 - taper) * np.abs(2 * y / b))
    xle = lambda y: 0.25 * c_root + np.abs(y) * np.tan(sweep_c4) - 0.25 * ch(y)
    s = np.linspace(0.0, 1.0, nx + 1)
    A, B, P = [], [], []
    for j in range(ny):
        y0, y1 = ye[j], ye[j + 1]
        ym = 0.5 * (y0 + y1)
        for i in range(nx):
            fb, fc = s[i] + 0.25 / nx, s[i] + 0.75 / nx
            A.append([xle(y0) + fb * ch(y0), y0, 0.0]); B.append([xle(y1) + fb * ch(y1), y1, 0.0])
            P.append([xle(ym) + fc * ch(ym), ym, 0.0])
    A, B, P = np.array(A), np.array(B), np.array(P)

    def leg(P, A, B):
        r0 = (B - A)[None]; r1 = P[:, None] - A[None]; r2 = P[:, None] - B[None]
        c = np.cross(r1, r2); c2 = np.sum(c * c, -1)
        k = np.sum(r0 * (r1 / np.linalg.norm(r1, axis=-1)[..., None] - r2 / np.linalg.norm(r2, axis=-1)[..., None]), -1)
        return (np.where(c2 < 1e-20, 0.0, k / np.where(c2 < 1e-20, 1.0, c2)) / (4 * np.pi))[..., None] * c

    def semi(P, A):
        r1 = P[:, None] - A[None]; d_ = np.array([1.0, 0, 0])
        c = np.cross(np.broadcast_to(d_, r1.shape), r1); c2 = np.sum(c * c, -1)
        k = (1 + r1[..., 0] / np.linalg.norm(r1, axis=-1))
        return (np.where(c2 < 1e-20, 0.0, k / np.where(c2 < 1e-20, 1.0, c2)) / (4 * np.pi))[..., None] * c

    w = (leg(P, A, B) + semi(P, B) - semi(P, A))[..., 2]
    gam = np.linalg.solve(w, -np.full(len(P), np.sin(alpha)))
    Gs = gam.reshape(ny, nx).sum(axis=1)
    yc = 0.5 * (ye[:-1] + ye[1:])
    S = 0.5 * b * c_root * (1 + taper)
    CL = float(2 * np.sum(Gs * np.diff(ye)) / S)
    return yc, 2 * Gs / ch(yc), CL


# =============================================================== loads from coefficients (learner)

def wing_from_coefficients(b: float, S: float, chord: Callable, A: np.ndarray, y: np.ndarray | None = None) -> dict:
    """CL = pi AR A_1, CDi = pi AR sum n A_n^2, e = CL^2/(pi AR CDi), the root bending moment
    coefficient C_B = M_root/(q S b) = (2/(V S b)) int_0^{b/2} Gamma y dy (one wing), and the local
    cl(y) = 2 Gamma/(V c) at the stations y (V = 1)."""
    AR = b ** 2 / S
    n = ll_odd_n(len(A))
    CL = np.pi * AR * A[0]
    CDi = np.pi * AR * np.sum(n * A ** 2)
    th = np.pi / 2 + (np.arange(2000) + 0.5) * (np.pi / 2) / 2000           # right wing, midpoint rule
    yy = -0.5 * b * np.cos(th)
    Gam = 2 * b * np.sin(np.outer(th, n)) @ A
    MB = np.sum(Gam * yy * 0.5 * b * np.sin(th)) * (np.pi / 2) / 2000
    out = dict(CL=float(CL), CDi=float(CDi), e=float(CL ** 2 / (np.pi * AR * CDi)) if CDi > 0 else np.nan,
               CB=float(2 * MB / (S * b)))
    if y is not None:
        out["cl"] = 2 * ll_gamma(b, A, y) / chord(np.asarray(y, dtype=float))
    return out


def efficiency_vs_taper(AR: float, tapers: Sequence[float], b: float = 10.0, N: int = 60) -> np.ndarray:
    """Span efficiency e of flat, untwisted straight-tapered wings of aspect ratio AR."""
    S = b ** 2 / AR
    zero = lambda y: np.zeros_like(np.asarray(y, dtype=float))
    return np.array([wing_from_coefficients(b, S, tapered_chord(b, S, t), ll_solve(b, tapered_chord(b, S, t), zero, 0.1, N))["e"]
                     for t in tapers])


def cl_ratio(b: float, S: float, chord: Callable, twist: Callable, CL: float, y: np.ndarray, N: int = 40) -> np.ndarray:
    """Local cl(y)/CL of a wing flying at CL (alpha found from the affine dependence of CL on alpha)."""
    A0_ = ll_solve(b, chord, twist, 0.0, N)
    A1_ = ll_solve(b, chord, twist, 1.0, N)
    al = (CL / (np.pi * b ** 2 / S) - A0_[0]) / (A1_[0] - A0_[0])
    A = A0_ + al * (A1_ - A0_)
    return wing_from_coefficients(b, S, chord, A, y)["cl"] / CL


def stall_onset(b: float, S: float, chord: Callable, twist: Callable, CL: float, n_y: int = 401, N: int = 40) -> float:
    """Spanwise station eta = 2y/b (0 root, 1 tip) where cl/CL is largest: with the same section
    cl_max everywhere, the wing stalls there first as CL rises."""
    y = np.linspace(0.0, 0.5 * b * 0.995, n_y)
    return float(2 * y[int(np.argmax(cl_ratio(b, S, chord, twist, CL, y, N)))] / b)


def sweep_stall_station(b: float, S: float, taper: float, sweeps_deg: Sequence[float]) -> np.ndarray:
    """Station eta of the largest cl/CL (right wing) for each quarter-chord sweep, by the given VLM:
    aft sweep loads the tips (eta moves outboard), forward sweep the root."""
    cr = 2 * S / (b * (1 + taper))
    out = []
    for sw in sweeps_deg:
        yc, cl, _ = vlm_strip_cl(b, cr, taper, np.radians(sw))
        right = yc > 0
        out.append(2 * yc[right][np.argmax(cl[right])] / b)
    return np.array(out)


# =============================================================== twist optimisation (learner)

def hat_twists(eta_k: np.ndarray) -> list[Callable]:
    """Piecewise-linear 'hat' functions in |eta| on the stations eta_k (eta_k[0] = 0 is the root,
    whose twist is fixed at zero, so hats are returned for eta_k[1:] only)."""
    eta_k = np.asarray(eta_k, dtype=float)

    def make(k):
        e = np.zeros_like(eta_k); e[k] = 1.0
        return lambda y, b: np.interp(np.abs(2 * np.asarray(y, dtype=float) / b), eta_k, e)
    return [make(k) for k in range(1, len(eta_k))]


def twist_basis(b: float, chord: Callable, eta_k: np.ndarray, N: int = 40) -> np.ndarray:
    """Columns of Glauert coefficients: column 0 for alpha = 1 rad with no twist, column k for a unit
    (1 rad) twist hat at station eta_k[k] with alpha = 0. The lifting line is linear, so any
    (alpha, twists) gives A = basis @ [alpha, t_1, ..., t_K]."""
    zero = lambda y: np.zeros_like(np.asarray(y, dtype=float))
    cols = [ll_solve(b, chord, zero, 1.0, N)]
    for h in hat_twists(eta_k):
        cols.append(ll_solve(b, chord, lambda y, h=h: h(y, b), 0.0, N))
    return np.column_stack(cols)


def optimize_twist(b: float, S: float, chord: Callable, CL_design: float, eta_k: np.ndarray, N: int = 40,
                   CB_max: float | None = None, tip_margin: float | None = None, CL_high: float | None = None,
                   eta_tip: float = 0.7, washout_only: bool = False) -> dict:
    """Minimise CDi at CL_design over alpha and the twists at eta_k[1:] (root twist 0) with SLSQP.
    Constraints: CL = CL_design (equality); optionally C_B <= CB_max at CL_design; optionally, at
    CL_high (the same twist, alpha re-trimmed), cl(eta) <= (1 - tip_margin) CL_high at stations
    eta >= eta_tip: the outboard wing runs below the wing's average lift coefficient, so a wing of
    constant section cl_max stalls inboard first; optionally (washout_only) the twist may only
    decrease outboard, t_0 = 0 >= t_1 >= ... >= t_K, as a builder would make it. Returns dict(alpha, twist (rad, at eta_k),
    A, CL, CDi, e, CB, success)."""
    Bm = twist_basis(b, chord, eta_k, N)
    AR = b ** 2 / S
    n = ll_odd_n(N)
    Q = (np.pi * AR) * (Bm.T * n) @ Bm                                       # CDi = x' Q x (quadratic)
    f = lambda x: float(x @ Q @ x)
    g = lambda x: 2 * Q @ x
    cons = [dict(type="eq", fun=lambda x: np.pi * AR * (Bm[0] @ x) - CL_design,
                 jac=lambda x: np.pi * AR * Bm[0])]
    if CB_max is not None:
        cons.append(dict(type="ineq", fun=lambda x: CB_max - wing_from_coefficients(b, S, chord, Bm @ x)["CB"]))
    if tip_margin is not None:
        ys = 0.5 * b * np.linspace(eta_tip, 0.98, 9)
        th = np.arccos(-2 * ys / b)
        Ey = 2 * b * np.sin(np.outer(th, n)) @ Bm * (2 / chord(ys))[:, None]     # cl at ys = Ey @ x

        def tipc(x):
            a_hi = x.copy()
            a_hi[0] = x[0] + (CL_high - CL_design) / (np.pi * AR * Bm[0, 0])    # re-trim alpha at CL_high
            return (1 - tip_margin) * CL_high - Ey @ a_hi
        cons.append(dict(type="ineq", fun=tipc))
    if washout_only:
        K = Bm.shape[1] - 1
        D = np.zeros((K, K + 1))
        for k in range(K):
            D[k, k + 1] = -1.0
            if k > 0:
                D[k, k] = 1.0                                                 # t_k - t_{k+1} >= 0, t_0 = 0
        cons.append(dict(type="ineq", fun=lambda x: D @ x, jac=lambda x: D))
    x0 = np.zeros(Bm.shape[1]); x0[0] = CL_design / (np.pi * AR * Bm[0, 0])
    res = minimize(f, x0, jac=g, constraints=cons, method="SLSQP", options=dict(ftol=1e-14, maxiter=500))
    A = Bm @ res.x
    w = wing_from_coefficients(b, S, chord, A)
    return dict(alpha=float(res.x[0]), twist=np.concatenate([[0.0], res.x[1:]]), A=A, CL=w["CL"], CDi=w["CDi"],
                e=w["e"], CB=w["CB"], success=bool(res.success))


def linear_twist_optimum(b: float, S: float, chord: Callable, CL: float, N: int = 40) -> dict:
    """Closed-form (Lagrange) optimum of a linear twist t(eta) = t_tip |eta|. With A = alpha a + t_tip r
    and A_1 = CL/(pi AR) fixed, alpha = (A_1 - t r_1)/a_1, so A_n = p_n + t q_n with
    p = A_1 a/a_1, q = r - r_1 a/a_1, and sum n A_n^2 is minimised at t* = -sum n p_n q_n / sum n q_n^2
    (the n = 1 term is fixed). Returns dict(tip_twist, alpha, e, CDi)."""
    zero = lambda y: np.zeros_like(np.asarray(y, dtype=float))
    a = ll_solve(b, chord, zero, 1.0, N)
    r = ll_solve(b, chord, lambda y: np.abs(2 * np.asarray(y) / b), 0.0, N)
    AR = b ** 2 / S
    n = ll_odd_n(N)
    A1 = CL / (np.pi * AR)
    p = A1 * a / a[0]
    q = r - r[0] * a / a[0]
    t = -np.sum(n * p * q) / np.sum(n * q * q)
    A = p + t * q
    w = wing_from_coefficients(b, S, chord, A)
    return dict(tip_twist=float(t), alpha=float((A1 - t * r[0]) / a[0]), e=w["e"], CDi=w["CDi"])


def prandtl_bell(n_terms: int = 6) -> dict:
    """Prandtl (1933): least induced drag for given lift L and given I = int L'(y) y^2 dy (twice the
    integral of the bending moment along the span, a proxy for the weight of a spar sized by
    bending), span free. With Gamma = sum G_n sin(n theta) (odd n), y = -(b/2) cos theta:
    L ~ b G_1,  I ~ b^3 (G_1 + G_3),  D_i ~ sum n G_n^2 (independent of b). In units of the elliptic
    wing of span b0 (G_1 = 1): G_1 = 1/s, G_3 = 1/s^3 - 1/s with s = b/b0; higher harmonics only add
    drag, so D_i/D_i,ell = f(s) = 1/s^2 + 3 (1/s^3 - 1/s)^2. The lift must stay non-negative along the
    span (a spar sized to |M| gains nothing from download), which needs G_3/G_1 >= -1/3, i.e.
    s <= sqrt(3/2). f falls monotonically up to that bound, where f'(s) = 0 as well: the stationary
    point of the Lagrangian is the bell Gamma ~ (1 - eta^2)^(3/2). Also checked by SLSQP over
    (s, G_1, G_3, ..., G_{2 n_terms - 1}) with the two equality constraints and Gamma >= 0.
    Returns dict(span_ratio, drag_ratio, dfds, G3_over_G1, root_moment_ratio, G, slsqp_span, slsqp_drag)."""
    f = lambda s: 1 / s ** 2 + 3 * (1 / s ** 3 - 1 / s) ** 2
    s_max = np.sqrt(1.5)                                                      # G3/G1 = 1/s^2 - 1 >= -1/3
    r = minimize_scalar(f, bounds=(1.0, s_max), method="bounded", options=dict(xatol=1e-12))
    s = float(r.x)
    h = 1e-6
    dfds = (f(s) - f(s - h)) / h
    g3 = 1 / s ** 2 - 1                                                       # G3 / G1
    nn = 2 * np.arange(n_terms) + 1
    tq = np.linspace(0.01, np.pi / 2, 60)
    cons = [dict(type="eq", fun=lambda z: z[0] * z[1] - 1.0),
            dict(type="eq", fun=lambda z: z[0] ** 3 * (z[1] + z[2]) - 1.0),
            dict(type="ineq", fun=lambda z: np.sin(np.outer(tq, nn)) @ z[1:])]
    z0 = np.concatenate([[1.1], [1 / 1.1], np.zeros(n_terms - 1)])
    zz = minimize(lambda z: float(np.sum(nn * z[1:] ** 2)), z0, constraints=cons, method="SLSQP",
                  options=dict(ftol=1e-15, maxiter=1000)).x
    # root bending moment of one wing, int_0^{b/2} Gamma y dy ~ b^2 sum G_n J_n, relative to the elliptic wing
    th = np.pi / 2 + (np.arange(4000) + 0.5) * (np.pi / 2) / 4000
    Jn = lambda nk: float(np.sum(np.sin(nk * th) * (-np.cos(th)) * np.sin(th)) * (np.pi / 2) / 4000)
    M_bell = s ** 2 * (1 / s) * (Jn(1) + g3 * Jn(3))
    return dict(span_ratio=s, drag_ratio=float(f(s)), dfds=float(dfds), G3_over_G1=float(g3),
                root_moment_ratio=float(M_bell / Jn(1)), G=zz[1:] / zz[1], slsqp_span=float(zz[0]),
                slsqp_drag=float(np.sum(nn * zz[1:] ** 2)))


def bending_trade(b: float, S: float, chord: Callable, CL: float, eta_k: np.ndarray, ratios: Sequence[float],
                  N: int = 40) -> list[dict]:
    """At fixed span and CL, the minimum-CDi twist when the root bending moment may not exceed
    `ratio` x that of the elliptic load (C_B,ell = CL/(3 pi)). Returns one optimize_twist result per ratio."""
    out = []
    for r in ratios:
        out.append(optimize_twist(b, S, chord, CL, eta_k, N, CB_max=r * CL / (3 * np.pi)))
    return out


# =============================================================== main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save planform.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    m_, S, b, V = TRAINER["mass"], TRAINER["S"], TRAINER["b"], TRAINER["V_cruise"]
    AR = b ** 2 / S
    W = m_ * G
    CLc = W / (0.5 * RHO * V ** 2 * S)
    CLa = W / (0.5 * RHO * TRAINER["V_approach"] ** 2 * S)
    zero = lambda y: np.zeros_like(np.asarray(y, dtype=float))

    # 1. e vs taper
    tapers = np.round(np.linspace(0.0, 1.0, 101), 2)
    e = efficiency_vs_taper(AR, tapers)
    k = int(np.argmax(e))
    print(f"flat untwisted wings, AR {AR:.3f}: e at taper 0, 0.2, 0.35, 0.5, 1.0 = "
          f"{e[0]:.4f}, {e[20]:.4f}, {e[35]:.4f}, {e[50]:.4f}, {e[100]:.4f}; maximum e = {e[k]:.5f} at taper {tapers[k]:.2f}")
    for ARx in (6.0, 12.0, 20.0):
        ex = efficiency_vs_taper(ARx, tapers)
        print(f"  AR {ARx:4.1f}: maximum e = {ex.max():.5f} at taper {tapers[int(np.argmax(ex))]:.2f}; rectangular e = {ex[-1]:.4f}")

    # 2. local cl and stall onset (unswept, lifting line)
    print("\nstall-onset station eta (largest cl/CL), flat untwisted, AR 8.33, and the peak cl/CL:")
    yy = np.linspace(0, 0.5 * b * 0.995, 401)
    for t in (1.0, 0.7, 0.5, 0.35, 0.2):
        ch = tapered_chord(b, S, t)
        r = cl_ratio(b, S, ch, zero, CLa, yy)
        print(f"  taper {t:4.2f}: stall starts at eta = {stall_onset(b, S, ch, zero, CLa):.3f}, max cl/CL = {r.max():.4f}, "
              f"tip (eta 0.9) cl/CL = {np.interp(0.45 * b, yy, r):.4f}")
    ch35 = tapered_chord(b, S, 0.35)
    for tw in (-2.0, -4.0):
        twf = lambda y, tw=tw: np.radians(tw) * np.abs(2 * np.asarray(y) / b)
        print(f"  taper 0.35 with linear washout {tw:.0f} deg: stall starts at eta = {stall_onset(b, S, ch35, twf, CLa):.3f}")

    # sweep (VLM)
    print("sweep, VLM (AR 8.33, taper 0.35, flat, alpha 5 deg): station of max cl/CL")
    sws = (-20.0, 0.0, 20.0, 35.0)
    for sw, et in zip(sws, sweep_stall_station(b, S, 0.35, sws)):
        yc, cl, CLv = vlm_strip_cl(b, 2 * S / (b * 1.35), 0.35, np.radians(sw))
        print(f"  quarter-chord sweep {sw:+5.1f} deg: eta = {et:.3f}, max cl/CL = {cl.max() / CLv:.4f}")

    # 3. twist optimisation for the trainer and a tapered variant
    eta_k = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    print(f"\ntwist optimisation at the cruise CL = {CLc:.4f} (5 twist stations, root twist 0):")
    rect = lambda y: np.full_like(np.asarray(y, dtype=float), TRAINER["c"])
    r0 = wing_from_coefficients(b, S, rect, ll_solve(b, rect, zero, 0.1))
    ro = optimize_twist(b, S, rect, CLc, eta_k)
    print(f"  rectangular trainer: untwisted e = {r0['e']:.4f}; optimum e = {ro['e']:.5f}, twist (deg) at eta {eta_k} = "
          f"{np.round(np.degrees(ro['twist']), 3)}, alpha {np.degrees(ro['alpha']):.3f} deg, C_B {ro['CB']:.5f} "
          f"(elliptic {CLc / (3 * np.pi):.5f})")
    lin = linear_twist_optimum(b, S, rect, CLc)
    print(f"  linear twist, Lagrange closed form: tip twist {np.degrees(lin['tip_twist']):.4f} deg, e = {lin['e']:.5f}")
    lso = optimize_twist(b, S, rect, CLc, np.array([0.0, 1.0]))
    print(f"  linear twist, SLSQP: tip twist {np.degrees(lso['twist'][1]):.4f} deg, e = {lso['e']:.5f}")
    for K in (2, 4, 8, 16):
        ek = np.linspace(0, 1, K + 1)
        print(f"    {K:2d} twist stations: optimum e = {optimize_twist(b, S, rect, CLc, ek)['e']:.6f}")

    lam = 0.35
    ch = tapered_chord(b, S, lam)
    print(f"\ntapered variant (taper {lam}, same S and b), design CL {CLc:.4f}, high-lift CL {CLa:.4f}:")
    u = wing_from_coefficients(b, S, ch, ll_solve(b, ch, zero, 0.1))
    print(f"  untwisted: e = {u['e']:.5f}, stall onset at eta {stall_onset(b, S, ch, zero, CLa):.3f}")
    oo = optimize_twist(b, S, ch, CLc, eta_k)
    print(f"  min CDi, unconstrained: e = {oo['e']:.5f}, twist {np.round(np.degrees(oo['twist']), 3)} deg, "
          f"stall onset at CL_high: eta {stall_onset(b, S, ch, lambda y: np.interp(np.abs(2 * np.asarray(y) / b), eta_k, oo['twist']), CLa):.3f}")
    for m in (0.0, 0.05, 0.10):
        oc = optimize_twist(b, S, ch, CLc, eta_k, tip_margin=m, CL_high=CLa, washout_only=True)
        twf = lambda y: np.interp(np.abs(2 * np.asarray(y) / b), eta_k, oc["twist"])
        print(f"  tip margin {m:.2f}: e = {oc['e']:.5f} (CDi +{(oc['CDi'] / oo['CDi'] - 1) * 100:.2f} %), twist "
              f"{np.round(np.degrees(oc['twist']), 3)} deg, stall onset eta {stall_onset(b, S, ch, twf, CLa):.3f}, ok {oc['success']}")

    # 4. bending trade and Prandtl's bell
    print("\nbending-moment trade, rectangular trainer at cruise CL (fixed span): C_B limit / elliptic C_B -> e")
    trade = bending_trade(b, S, rect, CLc, eta_k, [1.00, 0.97, 0.94, 0.91, 0.88])
    for r, t in zip([1.00, 0.97, 0.94, 0.91, 0.88], trade):
        print(f"  {r:.2f}: e = {t['e']:.5f}, C_B/C_B,ell = {t['CB'] / (CLc / (3 * np.pi)):.4f}, tip twist {np.degrees(t['twist'][-1]):.3f} deg")
    pb = prandtl_bell()
    print(f"Prandtl 1933 (lift and int L y^2 dy fixed, span free, lift >= 0): span x {pb['span_ratio']:.4f} "
          f"(sqrt(3/2) = {np.sqrt(1.5):.4f}), df/ds there {pb['dfds']:.1e}, induced drag x {pb['drag_ratio']:.4f} "
          f"({(pb['drag_ratio'] - 1) * 100:.2f} %), G3/G1 = {pb['G3_over_G1']:.4f} (bell: -1/3), root moment x {pb['root_moment_ratio']:.4f}")
    print(f"  SLSQP over (s, G1, G3, ..., G11): span x {pb['slsqp_span']:.4f}, drag x {pb['slsqp_drag']:.4f}, "
          f"G/G1 = {np.round(pb['G'], 5)}")

    # figure
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    ax[0].plot(tapers, e)
    ax[0].set_xlabel("taper ratio"); ax[0].set_ylabel("span efficiency e"); ax[0].set_title(f"flat untwisted, AR {AR:.2f}")
    for t in (1.0, 0.35, 0.2):
        ax[1].plot(2 * yy / b, cl_ratio(b, S, tapered_chord(b, S, t), zero, CLa, yy), label=f"taper {t}")
    twf = lambda y: np.interp(np.abs(2 * np.asarray(y) / b), eta_k, optimize_twist(b, S, ch, CLc, eta_k, tip_margin=0.05, CL_high=CLa, washout_only=True)["twist"])
    ax[1].plot(2 * yy / b, cl_ratio(b, S, ch, twf, CLa, yy), "--", label="taper 0.35, optimised twist")
    ax[1].set_xlabel("eta = 2y/b"); ax[1].set_ylabel("cl / CL"); ax[1].set_title("local lift at the approach CL"); ax[1].legend(fontsize=8)
    th = np.linspace(0, np.pi, 300)
    ax[2].plot(-np.cos(th), np.sin(th), label="elliptic, span b0")
    s_ = pb["span_ratio"]
    ax[2].plot(-s_ * np.cos(th), (np.sin(th) + pb["G3_over_G1"] * np.sin(3 * th)) / s_, label=f"Prandtl bell, span {s_:.2f} b0")
    ax[2].set_xlabel("y / (b0/2)"); ax[2].set_ylabel("Gamma (same lift)"); ax[2].set_title("same lift and int L y^2"); ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "planform.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
