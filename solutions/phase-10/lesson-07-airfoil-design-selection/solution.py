"""P10.07 Airfoil Design & Selection: reference solution.

Which section should the reference trainer's wing use, and could we have designed it? This module
answers in two halves.

Inverse design. Thin-airfoil theory (P10.04) maps a camber line to a chordwise load. Run it
backwards: write the target load difference Delta_cp(x) = cp_lower - cp_upper as a Glauert sine
series, read off the coefficients A_n, and the camber-line slope is a cosine series with the same
coefficients, dz/dx = alpha_i + sum A_n cos(n theta), x/c = (1 - cos theta)/2. Closing the camber
line (z = 0 at both ends) fixes the ideal angle alpha_i. A uniform load gives the NACA a = 1.0 mean
line; a target (cl, cm_c/4) pair gives a cubic camber line. Thickness is then added and the result
checked with a panel method, which shows what the thin theory cannot see.

Selection. A weighted decision matrix scores candidate sections on criteria computed by the
lesson's own models (thickness, cl_ideal against the mission's design cl, |cm_c/4|, the inviscid
suction peak at the mission's high-lift cl) for two missions (the reference trainer and a sailplane),
and a weight-sensitivity study says how robust the winner is.

Conventions (P10 brief): alpha positive nose-up, x along the chord from the leading edge (x/c in
0..1), z up, cm about the quarter chord, positive nose-up, V_inf = 1 in the panel solver,
Delta_cp = cp_lower - cp_upper (positive when the section lifts).

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes airfoil_design.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-10/lesson-07-airfoil-design-selection
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable, Sequence

import numpy as np

RHO = 1.225                 # kg/m^3, ISA sea level
G = 9.81                    # m/s^2
TRAINER = dict(mass=750.0, S=12.0, V_cruise=55.0, V_approach=30.0)

# criteria of the decision matrix: (name, larger_is_better)
CRITERIA = [("t_c", True), ("cl_mismatch", False), ("abs_cm", False), ("peak_suction", False)]
WEIGHTS = {"trainer": np.array([0.30, 0.15, 0.20, 0.35]),
           "sailplane": np.array([0.15, 0.40, 0.20, 0.25])}


# =============================================================== given helpers

def naca4(code: str, n: int = 80, chord: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Given (from P10.02): NACA 4-digit section with cosine spacing and a closed trailing edge.
    Returns 2n + 1 nodes ordered clockwise from the trailing edge: lower surface to the leading
    edge, then upper surface back to the trailing edge (first node = last node)."""
    m, p, t = int(code[0]) / 100, int(code[1]) / 10, int(code[2:]) / 100
    x = 0.5 * (1 - np.cos(np.linspace(0.0, np.pi, n + 1)))
    yt = naca_thickness(x, t)
    yc, dyc = naca4_camber(code)(x), naca4_slope(code)(x)
    th = np.arctan(dyc)
    xu, yu = x - yt * np.sin(th), yc + yt * np.cos(th)
    xl, yl = x + yt * np.sin(th), yc - yt * np.cos(th)
    return np.concatenate([xl[::-1], xu[1:]]) * chord, np.concatenate([yl[::-1], yu[1:]]) * chord


def naca_thickness(x: np.ndarray, t: float) -> np.ndarray:
    """Given (from P10.02): NACA 00xx half-thickness at x/c for thickness ratio t (closed TE)."""
    x = np.asarray(x, dtype=float)
    return 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1036 * x ** 4)


def naca4_camber(code: str) -> Callable[[np.ndarray], np.ndarray]:
    """Given (from P10.02): the NACA 4-digit mean line z(x)/c as a function of x/c."""
    m, p = int(code[0]) / 100, int(code[1]) / 10
    if m == 0:
        return lambda x: np.zeros_like(np.asarray(x, dtype=float))
    return lambda x: np.where(np.asarray(x) < p, m / p ** 2 * (2 * p * np.asarray(x) - np.asarray(x) ** 2),
                              m / (1 - p) ** 2 * (1 - 2 * p + 2 * p * np.asarray(x) - np.asarray(x) ** 2))


def naca4_slope(code: str) -> Callable[[np.ndarray], np.ndarray]:
    """Given (from P10.02): the slope dz/dx of the NACA 4-digit mean line as a function of x/c."""
    m, p = int(code[0]) / 100, int(code[1]) / 10
    if m == 0:
        return lambda x: np.zeros_like(np.asarray(x, dtype=float))
    return lambda x: np.where(np.asarray(x) < p, 2 * m / p ** 2 * (p - np.asarray(x)),
                              2 * m / (1 - p) ** 2 * (p - np.asarray(x)))


def thin_airfoil(dzdx: Callable[[np.ndarray], np.ndarray], alpha: float, n_terms: int = 30,
                 m: int = 4000) -> dict:
    """Given (from P10.04): Glauert coefficients of a camber line with slope dzdx(x/c), by the
    midpoint rule in theta (x/c = (1 - cos theta)/2); alpha in radians.
    A0 = alpha - (1/pi) int dz/dx dtheta,  A_n = (2/pi) int dz/dx cos(n theta) dtheta.
    Returns dict(A0, An (n = 1..n_terms), cl = 2 pi (A0 + A1/2), cm_c4 = pi/4 (A2 - A1), alpha_L0)."""
    th = (np.arange(m) + 0.5) * np.pi / m
    s = dzdx(0.5 * (1 - np.cos(th)))
    n = np.arange(1, n_terms + 1)
    A0 = alpha - s.mean()
    An = 2 * (s[None, :] * np.cos(n[:, None] * th[None, :])).mean(axis=1)
    aL0 = -(s * (np.cos(th) - 1)).mean()
    return dict(A0=A0, An=An, cl=2 * np.pi * (A0 + An[0] / 2), cm_c4=np.pi / 4 * (An[1] - An[0]),
                alpha_L0=aL0)


def _panel_influence(xp, yp, X, Y, own: bool = False):
    """Velocity at points (xp, yp) from unit-strength constant source panels and clockwise constant
    vortex panels with nodes (X, Y). Returns (us, vs, uv, vv), each (points, panels). own=True when
    the points are these panels' own midpoints: the self-induced angle is then set to pi exactly."""
    dx, dy = np.diff(X), np.diff(Y)
    L = np.hypot(dx, dy)
    tx, ty = dx / L, dy / L
    rx, ry = xp[:, None] - X[None, :-1], yp[:, None] - Y[None, :-1]
    xi, eta = rx * tx + ry * ty, -rx * ty + ry * tx
    beta = np.arctan2(eta, xi - L) - np.arctan2(eta, xi)
    if own:
        np.fill_diagonal(beta, np.pi)                     # own panel, outer side
    lnr = np.log(np.hypot(xi, eta) / np.hypot(xi - L, eta)) / (2 * np.pi)
    b = beta / (2 * np.pi)
    return lnr * tx - b * ty, lnr * ty + b * tx, b * tx + lnr * ty, b * ty - lnr * tx


def panel_solve(elements: Sequence[tuple[np.ndarray, np.ndarray]], alpha_deg: float,
                x_ref: float = 0.25, chord: float = 1.0) -> dict:
    """Given (stands in for the P10.05 solver): a Hess-Smith panel method with constant-strength
    sources on every panel, one constant vortex strength per element and one Kutta condition per
    element (equal tangential speed on its two trailing-edge panels). Each element is (X, Y) nodes,
    clockwise from the trailing edge, closed. V_inf = 1. Circulation clockwise positive.
    Returns dict(cl, cm (about (x_ref, 0)), cl_kj = 2 Gamma/(V c), xc, yc, cp, parts=[per element
    dict(xc, yc, cp, cl, cm)])."""
    a = np.radians(alpha_deg)
    Xs = [np.asarray(e[0], dtype=float) for e in elements]
    Ys = [np.asarray(e[1], dtype=float) for e in elements]
    idx, start = [], 0
    for x in Xs:
        idx.append(np.arange(start, start + len(x) - 1))
        start += len(x) - 1
    N, K = start, len(Xs)
    xa = np.concatenate([x[:-1] for x in Xs]); xb = np.concatenate([x[1:] for x in Xs])
    ya = np.concatenate([y[:-1] for y in Ys]); yb = np.concatenate([y[1:] for y in Ys])
    L = np.hypot(xb - xa, yb - ya)
    tx, ty = (xb - xa) / L, (yb - ya) / L
    nx, ny = -ty, tx
    xc, yc = 0.5 * (xa + xb), 0.5 * (ya + yb)
    US = np.zeros((N, N)); VS = np.zeros((N, N)); UV = np.zeros((N, K)); VV = np.zeros((N, K))
    for k in range(K):
        us, vs, uv, vv = _panel_influence(xc, yc, Xs[k], Ys[k])
        j = idx[k]
        own = _panel_influence(xc[j], yc[j], Xs[k], Ys[k], own=True)
        us[j], vs[j], uv[j], vv[j] = own
        US[:, idx[k]], VS[:, idx[k]] = us, vs
        UV[:, k], VV[:, k] = uv.sum(axis=1), vv.sum(axis=1)
    An = np.hstack([US * nx[:, None] + VS * ny[:, None], UV * nx[:, None] + VV * ny[:, None]])
    At = np.hstack([US * tx[:, None] + VS * ty[:, None], UV * tx[:, None] + VV * ty[:, None]])
    A = np.zeros((N + K, N + K)); rhs = np.zeros(N + K)
    A[:N], rhs[:N] = An, -(np.cos(a) * nx + np.sin(a) * ny)
    vt_inf = np.cos(a) * tx + np.sin(a) * ty
    for k in range(K):
        i1, i2 = idx[k][0], idx[k][-1]
        A[N + k], rhs[N + k] = At[i1] + At[i2], -(vt_inf[i1] + vt_inf[i2])
    sol = np.linalg.solve(A, rhs)
    vt = At @ sol + vt_inf
    cp = 1 - vt ** 2
    fx, fy = -cp * L * nx, -cp * L * ny
    mz = -((xc - x_ref) * fy - yc * fx)
    parts, clt, cmt = [], 0.0, 0.0
    for k in range(K):
        j = idx[k]
        clk = (-fx[j].sum() * np.sin(a) + fy[j].sum() * np.cos(a)) / chord
        cmk = mz[j].sum() / chord ** 2
        parts.append(dict(xc=xc[j], yc=yc[j], cp=cp[j], cl=clk, cm=cmk))
        clt, cmt = clt + clk, cmt + cmk
    gam = sum(sol[N + k] * L[idx[k]].sum() for k in range(K))
    return dict(cl=clt, cm=cmt, cl_kj=2 * gam / chord, xc=xc, yc=yc, cp=cp, parts=parts)


def alpha_for_cl(X: np.ndarray, Y: np.ndarray, cl_target: float) -> tuple[float, dict]:
    """Given: the angle of attack (deg) at which the panel solution of a single section gives
    cl_target, by two solves and the exact inviscid relation cl = c1 cos(alpha) + c2 sin(alpha)."""
    c1 = panel_solve([(X, Y)], 0.0)["cl"]
    c2 = panel_solve([(X, Y)], 90.0)["cl"]
    R, phi = np.hypot(c1, c2), np.arctan2(c1, c2)          # cl = R sin(alpha + phi)
    alpha = np.degrees(np.arcsin(cl_target / R) - phi)
    return alpha, panel_solve([(X, Y)], alpha)


def mission_cls(mass: float, S: float, V_cruise: float, V_slow: float, rho: float = RHO) -> tuple[float, float]:
    """Given: the lift coefficients W/(q S) at the cruise and the slow (approach or thermalling)
    speed. Section cl is taken equal to the wing CL here (the 3-D correction is lessons 10-12)."""
    W = mass * G
    return W / (0.5 * rho * V_cruise ** 2 * S), W / (0.5 * rho * V_slow ** 2 * S)


# =============================================================== inverse design (learner)

def glauert_load(A0: float, An: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Chordwise load Delta_cp(x) = cp_lower - cp_upper of a thin airfoil from its Glauert
    coefficients: gamma/(2 V) = A0 (1 + cos th)/sin th + sum A_n sin(n th), Delta_cp = 2 gamma/V
    = 4 [A0 cot(th/2) + sum A_n sin(n th)], with x/c = (1 - cos th)/2."""
    x = np.asarray(x, dtype=float)
    th = np.arccos(1 - 2 * x)
    n = np.arange(1, len(An) + 1)
    s = np.sin(np.multiply.outer(th, n)) @ np.asarray(An)
    with np.errstate(divide="ignore", invalid="ignore"):
        lead = np.where(A0 == 0.0, 0.0, A0 / np.tan(th / 2))
    return 4 * (lead + s)


def load_coefficients(dcp: Callable[[np.ndarray], np.ndarray], n_terms: int = 200, m: int = 8000) -> np.ndarray:
    """Inverse step 1: the sine-series coefficients A_1..A_N of a target load Delta_cp(x) that is
    finite at the leading edge (so A0 = 0: the design point is the ideal angle, with no leading-edge
    singularity). A_n = (2/pi) int_0^pi (Delta_cp/4) sin(n th) dth, midpoint rule with m points."""
    th = (np.arange(m) + 0.5) * np.pi / m
    f = dcp(0.5 * (1 - np.cos(th))) / 4.0
    n = np.arange(1, n_terms + 1)
    return 2 * (np.sin(np.multiply.outer(n, th)) @ f) / m


def ideal_angle(An: np.ndarray) -> float:
    """Inverse step 2: the ideal angle alpha_i (rad) that closes the camber line, z(0) = z(1) = 0.
    int_0^pi dz/dx sin th dth = 0 with dz/dx = alpha_i + sum A_n cos(n th) gives
    alpha_i = -(1/2) sum A_n I_n, I_n = int_0^pi cos(n th) sin th dth = (1 + cos n pi)/(1 - n^2)
    (I_1 = 0)."""
    n = np.arange(1, len(An) + 1).astype(float)
    I = np.where(n == 1, 0.0, (1 + np.cos(n * np.pi)) / np.where(n == 1, 1.0, 1 - n ** 2))
    return -0.5 * float(np.dot(An, I))


def camber_from_coefficients(An: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Inverse step 3: the camber line z(x)/c and its slope from the load coefficients.
    dz/dx = alpha_i + sum A_n cos(n th); z(th) = (1/2) int_0^th dz/dx sin th' dth', term by term:
    int_0^th cos(n t) sin t dt = (1/2)[(1 - cos((n+1)th))/(n+1) - (1 - cos((n-1)th))/(n-1)]
    (= sin^2(th)/2 for n = 1). Returns (z, dzdx, alpha_i)."""
    x = np.asarray(x, dtype=float)
    th = np.arccos(np.clip(1 - 2 * x, -1, 1))
    An = np.asarray(An, dtype=float)
    ai = ideal_angle(An)
    n = np.arange(1, len(An) + 1)
    dzdx = ai + np.cos(np.multiply.outer(th, n)) @ An
    T, Nn = np.multiply.outer(th, np.ones_like(n)), n[None, :].astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        Jn = 0.5 * ((1 - np.cos((Nn + 1) * T)) / (Nn + 1) - (1 - np.cos((Nn - 1) * T)) / (Nn - 1))
    Jn[:, 0] = 0.5 * np.sin(th) ** 2
    z = 0.5 * (ai * (1 - np.cos(th)) + Jn @ An)
    return z, dzdx, ai


def naca_a1_camber(x: np.ndarray, cl_i: float) -> tuple[np.ndarray, np.ndarray]:
    """The NACA a = 1.0 mean line (uniform load Delta_cp = cl_i along the whole chord), closed form:
    z = -(cl_i/(4 pi)) [x ln x + (1 - x) ln(1 - x)],  dz/dx = (cl_i/(4 pi)) ln((1 - x)/x)."""
    x = np.asarray(x, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        xl = np.where(x > 0, x * np.log(np.where(x > 0, x, 1.0)), 0.0)
        ol = np.where(x < 1, (1 - x) * np.log(np.where(x < 1, 1 - x, 1.0)), 0.0)
        slope = cl_i / (4 * np.pi) * np.log((1 - x) / x)
    return -cl_i / (4 * np.pi) * (xl + ol), slope


def camber_from_cl_cm(cl_d: float, cm_d: float) -> tuple[np.ndarray, float]:
    """The simplest camber line with a design cl and cm_c/4 at its ideal angle: a two-term load
    (A0 = 0), cl = pi A1, cm_c/4 = (pi/4)(A2 - A1). Returns (An = [A1, A2], alpha_i = A2/3)."""
    A1 = cl_d / np.pi
    A2 = A1 + 4 * cm_d / np.pi
    An = np.array([A1, A2])
    return An, ideal_angle(An)


def load_cl_cm(dcp: Callable[[np.ndarray], np.ndarray], m: int = 8000) -> tuple[float, float]:
    """cl = int_0^1 Delta_cp dx and cm_c/4 = -int_0^1 Delta_cp (x - 1/4) dx of a load (midpoint rule
    in theta, which handles a load singular like 1/sqrt(x) at the leading edge)."""
    th = (np.arange(m) + 0.5) * np.pi / m
    x = 0.5 * (1 - np.cos(th))
    w = 0.5 * np.sin(th) * np.pi / m                      # dx = (1/2) sin th dth
    d = dcp(x)
    return float(np.sum(d * w)), float(-np.sum(d * (x - 0.25) * w))


def inverse_design(dcp: Callable[[np.ndarray], np.ndarray], n_terms: int = 200, n_x: int = 201) -> dict:
    """Thin-airfoil inverse design: target load -> camber line. Returns dict(x, z, dzdx, An,
    alpha_i (rad), cl, cm_c4, slope = callable dz/dx(x) for the forward check)."""
    An = load_coefficients(dcp, n_terms)
    x = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n_x)))
    z, dzdx, ai = camber_from_coefficients(An, x)
    cl, cm = load_cl_cm(dcp)
    return dict(x=x, z=z, dzdx=dzdx, An=An, alpha_i=ai, cl=cl, cm_c4=cm,
                slope=lambda xx: camber_from_coefficients(An, xx)[1])


def thicken(camber: Callable[[np.ndarray], np.ndarray], slope: Callable[[np.ndarray], np.ndarray],
            t: float, n: int = 80) -> tuple[np.ndarray, np.ndarray]:
    """Wrap the NACA 00xx thickness of ratio t around a camber line (perpendicular to it, as the
    NACA sections do), cosine spacing, closed trailing edge, nodes clockwise from the trailing edge
    (lower surface first) as the panel solver expects."""
    x = 0.5 * (1 - np.cos(np.linspace(0.0, np.pi, n + 1)))
    xs = np.clip(x, 1e-9, 1 - 1e-9)                       # the slope may be log-singular at the ends
    yt, yc, th = naca_thickness(x, t), camber(x), np.arctan(slope(xs))
    xu, yu = x - yt * np.sin(th), yc + yt * np.cos(th)
    xl, yl = x + yt * np.sin(th), yc - yt * np.cos(th)
    return np.concatenate([xl[::-1], xu[1:]]), np.concatenate([yl[::-1], yu[1:]])


def panel_load(X: np.ndarray, Y: np.ndarray, alpha_deg: float, x_eval: np.ndarray) -> np.ndarray:
    """Delta_cp(x) = cp_lower(x) - cp_upper(x) from a panel solution, by linear interpolation of the
    two surfaces' cp at the stations x_eval (the panel counterpart of glauert_load)."""
    r = panel_solve([(X, Y)], alpha_deg)
    n_low = int(np.argmin(r["xc"]))                       # control points: lower TE -> LE, then upper
    xl, cl_ = r["xc"][: n_low + 1][::-1], r["cp"][: n_low + 1][::-1]
    xu, cu = r["xc"][n_low + 1:], r["cp"][n_low + 1:]
    return np.interp(x_eval, xl, cl_) - np.interp(x_eval, xu, cu)


# =============================================================== selection (learner)

def section_metrics(X: np.ndarray, Y: np.ndarray, slope: Callable[[np.ndarray], np.ndarray],
                    cl_design: float, cl_high: float) -> dict:
    """The decision-matrix criteria of one section, all from the lesson's own models:
    t_c (largest vertical thickness of the panel nodes), alpha_L0_deg and cm0 (= cm_c/4, thin
    airfoil), cl_ideal (thin-airfoil cl at A0 = 0), cl_mismatch = |cl_design - cl_ideal|, abs_cm,
    and peak_suction = -min cp from the panel method at cl_high (an inviscid stall-severity proxy:
    a higher leading-edge suction peak means a steeper adverse pressure gradient behind it)."""
    ta = thin_airfoil(slope, 0.0)
    a_ideal = -(ta["A0"])                                   # alpha at which A0 = 0 (A0 = alpha - mean)
    ti = thin_airfoil(slope, a_ideal)
    xs = np.linspace(0.0, 1.0, 401)
    n_le = int(np.argmin(X))
    yl = np.interp(xs, X[: n_le + 1][::-1], Y[: n_le + 1][::-1])
    yu = np.interp(xs, X[n_le:], Y[n_le:])
    _, r = alpha_for_cl(X, Y, cl_high)
    return dict(t_c=float(np.max(yu - yl)), alpha_L0_deg=float(np.degrees(ta["alpha_L0"])),
                cm0=float(ta["cm_c4"]), cl_ideal=float(ti["cl"]),
                cl_mismatch=float(abs(cl_design - ti["cl"])), abs_cm=float(abs(ta["cm_c4"])),
                peak_suction=float(-r["cp"].min()))


def decision_matrix(raw: np.ndarray, larger_is_better: Sequence[bool], weights: np.ndarray
                    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Weighted decision matrix. raw: (candidates, criteria). Each criterion is min-max normalised
    over the candidates to [0, 1] with 1 the best (a criterion with no spread scores 1 for all).
    Weights are normalised to sum to 1. Returns (normalised scores, totals, ranking best-first)."""
    raw = np.asarray(raw, dtype=float)
    lo, hi = raw.min(axis=0), raw.max(axis=0)
    span = np.where(hi > lo, hi - lo, 1.0)
    s = (raw - lo) / span
    s = np.where(np.asarray(larger_is_better)[None, :], s, 1 - s)
    s = np.where((hi > lo)[None, :], s, 1.0)
    w = np.asarray(weights, dtype=float) / np.sum(weights)
    tot = s @ w
    return s, tot, np.argsort(-tot, kind="stable")


def weight_monte_carlo(norm: np.ndarray, weights: np.ndarray, n_samples: int = 20000,
                       concentration: float = 20.0, seed: int = 7) -> np.ndarray:
    """Weight uncertainty: draw weight vectors from a Dirichlet distribution with mean `weights`
    and the given concentration (larger = closer to the nominal weights), seeded, and return the
    fraction of draws in which each candidate wins."""
    rng = np.random.default_rng(seed)
    w = np.asarray(weights, dtype=float) / np.sum(weights)
    W = rng.dirichlet(concentration * w, size=n_samples)
    winners = np.argmax(W @ np.asarray(norm).T, axis=1)
    return np.bincount(winners, minlength=norm.shape[0]) / n_samples


def weight_flip_points(norm: np.ndarray, weights: np.ndarray, k: int) -> list[tuple[float, int]]:
    """Vary the weight of criterion k from 0 to 1, rescaling the others in proportion so the sum
    stays 1. Every total is then linear in w_k, T_a(w) = w s_ak + (1 - w) R_a with
    R_a = sum_{j != k} w_j s_aj / (1 - w_k). Returns the exact w_k values at which the winner changes,
    with the new winner, as [(w, winner), ...] in increasing w."""
    norm = np.asarray(norm, dtype=float)
    w0 = np.asarray(weights, dtype=float) / np.sum(weights)
    rest = np.delete(np.arange(norm.shape[1]), k)
    R = norm[:, rest] @ w0[rest] / (1 - w0[k])
    Sk = norm[:, k]
    T = lambda w: w * Sk + (1 - w) * R
    cands = sorted({float(w) for a in range(len(R)) for b in range(len(R))
                    if a != b and abs((Sk[a] - R[a]) - (Sk[b] - R[b])) > 1e-15
                    for w in [(R[b] - R[a]) / ((Sk[a] - R[a]) - (Sk[b] - R[b]))] if 0 < w < 1})
    flips, cur = [], int(np.argmax(T(0.0)))
    for w in cands:
        nxt = int(np.argmax(T(min(w + 1e-9, 1.0))))
        if nxt != cur:
            flips.append((w, nxt))
            cur = nxt
    return flips


def build_candidates(n: int = 80) -> dict:
    """The candidate sections: five NACA 4-digit sections and two inverse-designed ones (the a = 1.0
    mean line with cl_i = 0.4 and 0.8, thickened with NACA 0012 thickness). Returns
    {name: (X, Y, slope)}."""
    out = {}
    for code in ("0012", "2412", "4412", "4415", "6409"):
        X, Y = naca4(code, n)
        out["NACA " + code] = (X, Y, naca4_slope(code))
    for cli in (0.4, 0.8):
        cam = lambda x, c=cli: naca_a1_camber(x, c)[0]
        slo = lambda x, c=cli: naca_a1_camber(np.clip(x, 1e-9, 1 - 1e-9), c)[1]
        X, Y = thicken(cam, slo, 0.12, n)
        out[f"a=1.0 cl_i={cli} t=12%"] = (X, Y, slo)
    return out


def select_sections(mission: str, cl_design: float, cl_high: float, n: int = 80) -> dict:
    """Run the selection for one mission: metrics of every candidate, the decision matrix with the
    mission's weights, the Monte Carlo win fractions and the flip points for each criterion."""
    cands = build_candidates(n)
    names = list(cands)
    mets = [section_metrics(*cands[nm], cl_design, cl_high) for nm in names]
    raw = np.array([[m[c] for c, _ in CRITERIA] for m in mets])
    norm, tot, rank = decision_matrix(raw, [b for _, b in CRITERIA], WEIGHTS[mission])
    win = weight_monte_carlo(norm, WEIGHTS[mission])
    flips = {CRITERIA[k][0]: weight_flip_points(norm, WEIGHTS[mission], k) for k in range(len(CRITERIA))}
    return dict(names=names, metrics=mets, raw=raw, norm=norm, totals=tot, rank=rank, win=win, flips=flips)


# =============================================================== main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save airfoil_design.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    np.set_printoptions(precision=4, suppress=True)

    # 1. uniform load -> a = 1.0 mean line
    cl_i = 0.4
    des = inverse_design(lambda x: np.full_like(x, cl_i), n_terms=200)
    z_ex, _ = naca_a1_camber(des["x"], cl_i)
    print(f"uniform load Delta_cp = {cl_i}: alpha_i = {np.degrees(des['alpha_i']):.2e} deg, "
          f"max camber {des['z'].max():.5f} c at x = {des['x'][np.argmax(des['z'])]:.3f} c "
          f"(closed form {cl_i / (4 * np.pi) * np.log(2):.5f} c at 0.5 c)")
    print("  N terms   max |z - z_a1| / c")
    errs, Ns = [], [25, 50, 100, 200]
    for N in Ns:
        An = load_coefficients(lambda x: np.full_like(x, cl_i), N)
        z, _, _ = camber_from_coefficients(An, des["x"])
        errs.append(np.max(np.abs(z - z_ex)))
        print(f"  {N:5d}     {errs[-1]:.3e}")
    order = np.log(np.array(errs[:-1]) / np.array(errs[1:])) / np.log(2)
    print(f"  observed orders {order}")
    fwd = thin_airfoil(des["slope"], des["alpha_i"], n_terms=60, m=20000)
    print(f"forward thin-airfoil check of the recovered camber line at alpha_i: cl = {fwd['cl']:.4f} "
          f"(target {des['cl']:.4f}), cm_c/4 = {fwd['cm_c4']:.4f} (target {des['cm_c4']:.4f}), "
          f"A0 = {fwd['A0']:.1e}")
    xs = np.linspace(0.05, 0.95, 19)
    print(f"  forward load at x = 0.05..0.95: min {glauert_load(fwd['A0'], fwd['An'], xs).min():.4f}, "
          f"max {glauert_load(fwd['A0'], fwd['An'], xs).max():.4f} (target {cl_i})")

    # 2. a roof-top load: uniform to x = 0.6, then linear to zero at the trailing edge
    roof = lambda x, a=0.6, c=0.5: c * np.where(x < a, 1.0, (1 - x) / (1 - a)) / (1 - 0.5 * (1 - a))
    dr = inverse_design(roof, n_terms=200)
    fr = thin_airfoil(dr["slope"], dr["alpha_i"], n_terms=200, m=20000)
    xs2 = np.linspace(0.02, 0.98, 49)
    print(f"\nroof-top load (a = 0.6, cl = {dr['cl']:.3f}): alpha_i = {np.degrees(dr['alpha_i']):.3f} deg, "
          f"cm_c/4 = {dr['cm_c4']:.4f}, max camber {dr['z'].max():.4f} c at x = {dr['x'][np.argmax(dr['z'])]:.3f}")
    print(f"  forward check: cl = {fr['cl']:.4f}, cm_c/4 = {fr['cm_c4']:.4f}, "
          f"max |load error| on 0.02-0.98 = {np.max(np.abs(glauert_load(fr['A0'], fr['An'], xs2) - roof(xs2))):.4f}")

    # 3. target (cl, cm) -> cubic camber line, compared with NACA 2412
    t24 = thin_airfoil(naca4_slope("2412"), 0.0)
    An2, ai2 = camber_from_cl_cm(0.30, -0.053)
    xg = 0.5 * (1 - np.cos(np.linspace(0, np.pi, 201)))
    z2, s2, _ = camber_from_coefficients(An2, xg)
    chk = thin_airfoil(lambda x: camber_from_coefficients(An2, x)[1], ai2)
    print(f"\nNACA 2412 by thin airfoil: alpha_L0 = {np.degrees(t24['alpha_L0']):.3f} deg, cm_c/4 = {t24['cm_c4']:.4f}")
    print(f"target cl = 0.30, cm_c/4 = -0.053: A1 = {An2[0]:.4f}, A2 = {An2[1]:.4f}, alpha_i = {np.degrees(ai2):.3f} deg, "
          f"max camber {z2.max():.4f} c at x = {xg[np.argmax(z2)]:.3f}; forward cl = {chk['cl']:.4f}, cm = {chk['cm_c4']:.4f}")

    # 4. thickness: the panel method against the thin-airfoil design load
    print("\na = 1.0, cl_i = 0.4 with NACA 00tt thickness, panel method at alpha_i (thin theory cl = 0.4):")
    print("   t/c     cl_panel   mean |Delta_cp - 0.4| on 0.1-0.9")
    cam = lambda x: naca_a1_camber(x, cl_i)[0]
    slo = lambda x: naca_a1_camber(np.clip(x, 1e-9, 1 - 1e-9), cl_i)[1]
    xe = np.linspace(0.1, 0.9, 81)
    thick_rows = []
    for t in (0.02, 0.06, 0.12, 0.18):
        X, Y = thicken(cam, slo, t, 120)
        r = panel_solve([(X, Y)], 0.0)
        dev = np.mean(np.abs(panel_load(X, Y, 0.0, xe) - cl_i))
        thick_rows.append((t, r["cl"], dev))
        print(f"   {t:.2f}    {r['cl']:.4f}     {dev:.4f}")

    # 5. selection
    q = {}
    for mission, (cd, ch) in {"trainer": mission_cls(TRAINER["mass"], TRAINER["S"], TRAINER["V_cruise"],
                                                       TRAINER["V_approach"]),
                              "sailplane": (0.70, 1.20)}.items():
        res = select_sections(mission, cd, ch)
        q[mission] = res
        print(f"\nmission {mission}: cl_design = {cd:.3f}, cl_high = {ch:.3f}, weights {WEIGHTS[mission]}")
        print("  section                  t/c    aL0(deg)  cm0      cl_ideal  -cp_min@cl_high  total")
        for i, nm in enumerate(res["names"]):
            m = res["metrics"][i]
            print(f"  {nm:22s} {m['t_c']:.3f}  {m['alpha_L0_deg']:7.3f}  {m['cm0']:7.4f}  {m['cl_ideal']:7.4f}  "
                  f"{m['peak_suction']:8.3f}       {res['totals'][i]:.3f}")
        print("  ranking: " + " > ".join(res["names"][i] for i in res["rank"][:4]))
        print("  Monte Carlo win fractions (Dirichlet, concentration 20, seed 7): " +
              ", ".join(f"{res['names'][i]} {res['win'][i]:.3f}" for i in np.argsort(-res["win"]) if res["win"][i] > 0))
        for c, fl in res["flips"].items():
            print(f"  weight of {c:13s}: winner changes at " +
                  (", ".join(f"w = {w:.3f} -> {res['names'][j]}" for w, j in fl) if fl else "never"))

    # figure
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    ax[0].plot(des["x"], des["z"], label="inverse design, uniform load")
    ax[0].plot(des["x"], z_ex, "k--", lw=0.8, label="NACA a = 1.0 closed form")
    ax[0].plot(dr["x"], dr["z"], label="inverse design, roof-top a = 0.6")
    ax[0].plot(xg, z2, label="target cl = 0.3, cm = -0.053")
    ax[0].plot(xg, naca4_camber("2412")(xg), ":", label="NACA 2412 mean line")
    ax[0].set_xlabel("x/c"); ax[0].set_ylabel("z/c"); ax[0].set_title("camber lines"); ax[0].legend(fontsize=7)
    for t in (0.02, 0.12):
        X, Y = thicken(cam, slo, t, 120)
        ax[1].plot(xe, panel_load(X, Y, 0.0, xe), label=f"panel, t/c = {t}")
    ax[1].axhline(cl_i, color="k", ls="--", lw=0.8, label="design load (thin theory)")
    ax[1].set_xlabel("x/c"); ax[1].set_ylabel("cp_lower - cp_upper"); ax[1].set_title("design load with thickness")
    ax[1].legend(fontsize=7)
    for j, mission in enumerate(("trainer", "sailplane")):
        r = q[mission]
        ax[2].bar(np.arange(len(r["names"])) + 0.4 * j, r["totals"], width=0.4, label=mission)
    ax[2].set_xticks(np.arange(len(q["trainer"]["names"])) + 0.2)
    ax[2].set_xticklabels([n.replace(" t=12%", "") for n in q["trainer"]["names"]], rotation=35, fontsize=7, ha="right")
    ax[2].set_ylabel("weighted score"); ax[2].set_title("decision matrix"); ax[2].legend()
    fig.tight_layout()
    fig.savefig(out_dir / "airfoil_design.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
