"""P10.08 High-Lift Devices: starter code.

How slowly can the reference trainer approach, and what do flaps and slats buy? This module builds
the answer from three models of different fidelity:

* thin-airfoil flap theory (P10.04): a plain flap of chord fraction E = cf/c deflected by delta is a
  kinked camber line; its Glauert coefficients give the flap effectiveness
  tau = d(alpha_L0)/d(delta) = -[1 - (theta_h - sin theta_h)/pi], the lift increment
  2 pi tau delta, the moment increment -(1/2) sin theta_h (1 - cos theta_h) delta and the hinge
  moment, where cos theta_h = 2E - 1 places the hinge;
* a multi-element panel method (the given Hess-Smith solver standing in for P10.05) for a plain
  flap with thickness, a leading-edge slat and a slotted flap, which shows A.M.O. Smith's slat effect
  (the slat lowers the main element's leading-edge suction peak) and circulation effect (a flap
  raises the main element's lift);
* a suction-peak criterion for cl_max (stall when the highest suction peak reaches the value the
  clean section has at its stall), calibrated on an assumed clean cl_max, and the landing-speed
  payoff V_s = sqrt(2 W / (rho S CL_max)).

Conventions (P10 brief): alpha positive nose-up; flap deflection delta positive trailing edge down;
cm about the quarter chord, positive nose-up; hinge moment H positive when it would deflect the
flap trailing edge down, ch = H / (q cf^2) per unit span; chord c = 1 (the stowed chord) for every
coefficient; V_inf = 1 in the panel solver; ISA sea level for the trainer.

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes high_lift.png)
Test: pytest curriculum/phase-10/lesson-08-high-lift-devices
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
from scipy.optimize import brentq

RHO = 1.225                         # kg/m^3, ISA sea level
G = 9.81                            # m/s^2
TRAINER = dict(mass=750.0, S=12.0, V_approach=30.0)
CLMAX_CLEAN = 1.6                   # ASSUMED 2-D cl_max of the clean NACA 2412 (order of the NACA data at
                                    # Re of a few million); the calibration of the peak criterion. P10.15 fits stall to data.
K_3D = 0.90                         # ASSUMED CL_max(wing) / cl_max(section) for an unswept, untwisted wing
FLAPPED_SPAN = 0.60                 # ASSUMED fraction of the wing area ahead of the flaps (partial-span flaps)


# =============================================================== given helpers

def naca4(code: str, n: int = 80, chord: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Given (from P10.02): NACA 4-digit section, cosine spacing, closed trailing edge; nodes
    clockwise from the trailing edge (lower surface to the LE, then upper surface back)."""
    m, p, t = int(code[0]) / 100, int(code[1]) / 10, int(code[2:]) / 100
    x = 0.5 * (1 - np.cos(np.linspace(0.0, np.pi, n + 1)))
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1036 * x ** 4)
    if m > 0:
        yc = np.where(x < p, m / p ** 2 * (2 * p * x - x ** 2), m / (1 - p) ** 2 * (1 - 2 * p + 2 * p * x - x ** 2))
        dyc = np.where(x < p, 2 * m / p ** 2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))
    else:
        yc, dyc = np.zeros_like(x), np.zeros_like(x)
    th = np.arctan(dyc)
    xu, yu = x - yt * np.sin(th), yc + yt * np.cos(th)
    xl, yl = x + yt * np.sin(th), yc - yt * np.cos(th)
    return np.concatenate([xl[::-1], xu[1:]]) * chord, np.concatenate([yl[::-1], yu[1:]]) * chord


def naca4_camber(code: str) -> Callable[[np.ndarray], np.ndarray]:
    """Given (from P10.02): the NACA 4-digit mean line z(x)/c."""
    m, p = int(code[0]) / 100, int(code[1]) / 10
    if m == 0:
        return lambda x: np.zeros_like(np.asarray(x, dtype=float))
    return lambda x: np.where(np.asarray(x) < p, m / p ** 2 * (2 * p * np.asarray(x) - np.asarray(x) ** 2),
                              m / (1 - p) ** 2 * (1 - 2 * p + 2 * p * np.asarray(x) - np.asarray(x) ** 2))


def thin_airfoil(dzdx: Callable[[np.ndarray], np.ndarray], alpha: float, n_terms: int = 30,
                 m: int = 4000) -> dict:
    """Given (from P10.04): Glauert coefficients of a camber line with slope dzdx(x/c), midpoint rule
    in theta (x/c = (1 - cos theta)/2); alpha in radians. Returns dict(A0, An (n = 1..n_terms),
    cl = 2 pi (A0 + A1/2), cm_c4 = pi/4 (A2 - A1), alpha_L0)."""
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


def panel_basis(elements: Sequence[tuple[np.ndarray, np.ndarray]]) -> dict:
    """Given (stands in for the P10.05 solver): Hess-Smith panel method for one or more elements:
    constant sources on every panel, one constant vortex strength and one Kutta condition (equal
    tangential speed on its two trailing-edge panels) per element. Each element is (X, Y) nodes,
    clockwise from its trailing edge, closed. V_inf = 1, circulation clockwise positive.
    The flow is linear in the free stream, so the system is solved once for a stream along x and
    once along y; panel_at() then gives the solution at any alpha by superposition."""
    Xs = [np.asarray(e[0], dtype=float) for e in elements]
    Ys = [np.asarray(e[1], dtype=float) for e in elements]
    idx, upper, start = [], [], 0
    for x in Xs:
        npan = len(x) - 1
        idx.append(np.arange(start, start + npan))
        upper.append(np.arange(npan) >= int(np.argmin(x)))        # panels from the LE node to the TE
        start += npan
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
    A = np.zeros((N + K, N + K)); rhs = np.zeros((N + K, 2))
    A[:N], rhs[:N, 0], rhs[:N, 1] = An, -nx, -ny
    for k in range(K):
        i1, i2 = idx[k][0], idx[k][-1]
        A[N + k] = At[i1] + At[i2]
        rhs[N + k] = [-(tx[i1] + tx[i2]), -(ty[i1] + ty[i2])]
    sol = np.linalg.solve(A, rhs)
    vt = At @ sol + np.column_stack([tx, ty])                       # columns: stream along x, along y
    gam = np.array([[sol[N + k, j] * L[idx[k]].sum() for j in range(2)] for k in range(K)])
    return dict(vt=vt, L=L, nx=nx, ny=ny, xc=xc, yc=yc, idx=idx, upper=upper, gam=gam,
                elements=[(x, y) for x, y in zip(Xs, Ys)])


def panel_at(basis: dict, alpha_deg: float, x_ref: float = 0.25, chord: float = 1.0) -> dict:
    """Given: the panel solution at alpha from panel_basis(). Returns dict(cl, cm (about (x_ref, 0)),
    cl_kj = 2 Gamma/(V c), xc, yc, cp, parts=[dict(xc, yc, cp, cl, cm, upper)])."""
    a = np.radians(alpha_deg)
    b = basis
    vt = b["vt"] @ np.array([np.cos(a), np.sin(a)])
    cp = 1 - vt ** 2
    fx, fy = -cp * b["L"] * b["nx"], -cp * b["L"] * b["ny"]
    mz = -((b["xc"] - x_ref) * fy - b["yc"] * fx)
    parts, clt, cmt = [], 0.0, 0.0
    for j, up in zip(b["idx"], b["upper"]):
        clk = (-fx[j].sum() * np.sin(a) + fy[j].sum() * np.cos(a)) / chord
        cmk = mz[j].sum() / chord ** 2
        parts.append(dict(xc=b["xc"][j], yc=b["yc"][j], cp=cp[j], cl=clk, cm=cmk, upper=up))
        clt, cmt = clt + clk, cmt + cmk
    gam = float(b["gam"].sum(axis=0) @ np.array([np.cos(a), np.sin(a)]))
    return dict(cl=clt, cm=cmt, cl_kj=2 * gam / chord, xc=b["xc"], yc=b["yc"], cp=cp, parts=parts)


def panel_solve(elements: Sequence[tuple[np.ndarray, np.ndarray]], alpha_deg: float,
                x_ref: float = 0.25, chord: float = 1.0) -> dict:
    """Given: panel_at(panel_basis(elements), alpha_deg) in one call."""
    return panel_at(panel_basis(elements), alpha_deg, x_ref, chord)


def place_element(X: np.ndarray, Y: np.ndarray, chord: float, deflection_deg: float,
                  x_pivot: float, y_pivot: float, pivot: str = "te") -> tuple[np.ndarray, np.ndarray]:
    """Given: scale a unit-chord section to `chord`, rotate it by deflection_deg (positive = nose
    down for a slat about its trailing edge, trailing edge down for a flap about its leading edge;
    both are a counter-clockwise turn about the TE pivot or a clockwise turn about the LE pivot) and
    put the pivot (its 'te' or 'le') at (x_pivot, y_pivot)."""
    d = np.radians(deflection_deg)
    x, y = np.asarray(X) * chord, np.asarray(Y) * chord
    if pivot == "te":
        x = x - chord
        return x * np.cos(d) - y * np.sin(d) + x_pivot, x * np.sin(d) + y * np.cos(d) + y_pivot
    return x * np.cos(d) + y * np.sin(d) + x_pivot, -x * np.sin(d) + y * np.cos(d) + y_pivot


def slat_configuration(alpha_free: bool = True, n_main: int = 100, n_slat: int = 50,
                       slat_code: str = "6418", slat_chord: float = 0.15, slat_deflection: float = 30.0,
                       slat_te: tuple[float, float] = (0.03, 0.045)) -> list:
    """Given: the lesson's two-element configuration: a NACA 2412 main element (chord 1, LE at the
    origin) and a slat of `slat_code` shape and chord 0.15, deflected 30 deg nose down about its
    trailing edge, which sits at (0.03, 0.045), above the main element's nose (gap about 0.015 c).
    Returns [slat, main] as (X, Y) pairs."""
    Xs, Ys = naca4(slat_code, n_slat)
    return [place_element(Xs, Ys, slat_chord, slat_deflection, *slat_te, pivot="te"), naca4("2412", n_main)]


def slotted_flap_configuration(flap_deflection: float = 30.0, n_main: int = 100, n_flap: int = 50,
                               flap_chord: float = 0.30, flap_le: tuple[float, float] = (0.95, -0.03)) -> list:
    """Given: a NACA 2412 main element (chord 1) with a NACA 4412-shaped slotted flap of chord 0.30
    whose leading edge sits at (0.95, -0.03), under and just ahead of the main trailing edge (5 %
    overlap, 3 % gap), deflected trailing edge down about its leading edge. Returns [main, flap]."""
    Xf, Yf = naca4("4412", n_flap)
    return [naca4("2412", n_main), place_element(Xf, Yf, flap_chord, flap_deflection, *flap_le, pivot="le")]


# =============================================================== thin-airfoil flap theory (learner)

def hinge_angle(cf_c: float) -> float:
    """Glauert angle theta_h of the hinge at x_h/c = 1 - cf/c: x = (1 - cos theta)/2 gives
    cos theta_h = 1 - 2 x_h = 2 cf/c - 1."""
    raise NotImplementedError


def flap_effectiveness(cf_c: float) -> float:
    """tau = -d(alpha_L0)/d(delta) = 1 - (theta_h - sin theta_h)/pi (thin airfoil, plain flap)."""
    raise NotImplementedError


def flap_increments(cf_c: float, delta: float) -> dict:
    """Closed-form thin-airfoil increments of a plain flap deflected delta (rad, TE down):
    dcl = 2 pi tau delta, dalpha_L0 = -tau delta, dcm_c4 = -(1/2) sin th_h (1 - cos th_h) delta,
    from A0 += delta (pi - th_h)/pi, A1 = 2 delta sin th_h / pi, A2 = delta sin 2 th_h / pi."""
    raise NotImplementedError


def flap_slope(cf_c: float, delta: float, base: Callable[[np.ndarray], np.ndarray] | None = None
               ) -> Callable[[np.ndarray], np.ndarray]:
    """Camber-line slope of a section with a plain flap: the base slope (zero for a flat plate)
    plus -delta aft of the hinge x_h = 1 - cf/c (small-angle thin-airfoil model)."""
    raise NotImplementedError


def hinge_moment(cf_c: float, alpha: float, delta: float, n_terms: int = 300, m: int = 8000) -> float:
    """Thin-airfoil hinge-moment coefficient ch = H/(q cf^2), H positive TE down, of a flat plate
    with a plain flap at alpha (rad) and delta (rad):
    ch = -(1/cf^2) int_{x_h}^{1} Delta_cp (x - x_h) dx, with Delta_cp = 4[A0 cot(th/2) + sum A_n sin n th]
    from the Glauert coefficients (given thin_airfoil), integrated in theta (dx = sin th dth / 2)."""
    raise NotImplementedError


def flapped_section(code: str, cf_c: float, delta_deg: float, n: int = 100) -> tuple[np.ndarray, np.ndarray]:
    """A NACA 4-digit section with a plain flap: every node aft of the hinge x_h = 1 - cf/c is turned
    TE-down by delta about the hinge point on the mean line (x_h, z_c(x_h)). Lower-surface nodes that
    the rotation pushes forward of their unrotated neighbour are dropped so the contour does not fold.
    Nodes stay clockwise from the trailing edge."""
    raise NotImplementedError


def panel_flap_derivative(t: float, cf_c: float, d_delta_deg: float = 2.0, n: int = 120) -> float:
    """dcl/d(delta) (per rad) of a NACA 00tt section with a plain flap at alpha = 0, by the panel
    method and a central difference over +-d_delta_deg. As t -> 0 it should approach 2 pi tau."""
    raise NotImplementedError


# =============================================================== multi-element and cl_max (learner)

def element_peaks(basis: dict, alpha_deg: float, le_window: float = 0.2) -> dict:
    """Panel solution (from a panel_basis) at alpha with each element's cl and its leading-edge
    suction peak: the minimum cp on the element's upper surface within the forward le_window
    fraction of its own x-extent. Returns dict(cl, cls=[...], peaks=[...], result=panel_at dict)."""
    raise NotImplementedError


def critical_peak(clean: tuple[np.ndarray, np.ndarray], cl_max: float = CLMAX_CLEAN) -> tuple[float, float]:
    """Calibrate the suction-peak criterion: the angle (deg) at which the clean section's panel cl
    equals the assumed cl_max, and its leading-edge peak cp there (the critical peak).
    Returns (alpha_stall_deg, cp_crit)."""
    raise NotImplementedError


def clmax_by_peak(elements: list, cp_crit: float, checked: Sequence[int] | None = None,
                  alpha_lo: float = -5.0, alpha_hi: float = 30.0) -> dict:
    """Suction-peak estimate of cl_max: the smallest alpha in [alpha_lo, alpha_hi] at which the
    upper-surface leading-edge peak of any element listed in `checked` (default: all) reaches
    cp_crit (scanned in 0.5 deg steps, then refined by brentq). A slotted flap is normally left out
    of `checked`: its boundary layer starts fresh in the slot flow and tolerates a far higher peak
    than the one calibrated on the clean leading edge. Returns dict(alpha, cl, element (index of
    the limiting element), peaks (all elements' peaks at alpha))."""
    raise NotImplementedError


def stall_speed(mass: float, S: float, CL_max: float, rho: float = RHO) -> float:
    """V_s = sqrt(2 W / (rho S CL_max)), W = m g (1 g level flight at the stall)."""
    raise NotImplementedError


def wing_clmax(cl_max_2d_clean: float, dcl_max_2d: float, k3d: float = K_3D,
               flapped: float = FLAPPED_SPAN) -> float:
    """Wing CL_max from section values: k3d (cl_max_clean + dcl_max * S_flapped/S). The simple
    area-ratio rule for partial-span devices; k3d and the area ratio are assumptions."""
    raise NotImplementedError


# =============================================================== main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save high_lift.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)

    # 1. thin-airfoil flap theory
    print("thin-airfoil plain flap: E = cf/c, theta_h, tau, dcl/ddelta (per rad), dcm/ddelta, ch_alpha, ch_delta (per rad)")
    for E in (0.1, 0.2, 0.25, 0.3, 0.4, 0.5):
        inc = flap_increments(E, 1.0)
        cha = (hinge_moment(E, 0.01, 0.0) - hinge_moment(E, -0.01, 0.0)) / 0.02
        chd = (hinge_moment(E, 0.0, 0.01) - hinge_moment(E, 0.0, -0.01)) / 0.02
        print(f"  {E:.2f}  {np.degrees(hinge_angle(E)):6.2f} deg  {inc['tau']:.4f}  {inc['dcl']:.4f}  "
              f"{inc['dcm_c4']:.4f}  {cha:.4f}  {chd:.4f}")
    E, dd = 0.25, np.radians(20.0)
    inc = flap_increments(E, dd)
    num = thin_airfoil(flap_slope(E, dd), 0.0, n_terms=10, m=40000)
    print(f"E = 0.25, delta = 20 deg: dcl = {inc['dcl']:.4f} (quadrature {num['cl']:.4f}), dalpha_L0 = "
          f"{np.degrees(inc['dalpha_L0']):.3f} deg, dcm_c4 = {inc['dcm_c4']:.4f} (quadrature {num['cm_c4']:.4f})")

    # 2. independent route: panel method with a deflected plain flap as thickness -> 0
    tau = flap_effectiveness(E)
    print(f"\npanel dcl/ddelta (alpha = 0, +-2 deg) against thin theory 2 pi tau = {2 * np.pi * tau:.4f}:")
    rows = []
    for t in (0.02, 0.04, 0.08, 0.12):
        d = panel_flap_derivative(t, E)
        rows.append((t, d))
        print(f"  t/c = {t:.2f}: {d:.4f}  (ratio {d / (2 * np.pi * tau):.4f})")
    tt = np.array([r[0] for r in rows]); dv = np.array([r[1] for r in rows])
    fit = np.polyfit(tt, dv, 1)
    print(f"  linear fit through the four: intercept at t = 0 {fit[1]:.4f} ({(fit[1] / (2 * np.pi * tau) - 1) * 100:+.2f} %)")

    # 3. clean calibration and the slat effect
    clean = naca4("2412", 100)
    a_st, cp_crit = critical_peak(clean)
    print(f"\nclean NACA 2412: assumed cl_max = {CLMAX_CLEAN}, reached (inviscid) at alpha = {a_st:.2f} deg; "
          f"critical LE peak cp = {cp_crit:.3f}")
    slat = slat_configuration()
    print("slat effect (Smith): main-element LE peak cp, clean vs with slat, same alpha")
    print("  alpha   clean cl  clean peak   slat+main cl  main cl  main peak  slat peak")
    for a in (0.0, 5.0, 10.0, 15.0):
        c0 = element_peaks(panel_basis([clean]), a)
        c1 = element_peaks(panel_basis(slat), a)
        print(f"  {a:5.1f}   {c0['cl']:.3f}    {c0['peaks'][0]:7.3f}     {c1['cl']:.3f}        {c1['cls'][1]:.3f}    "
              f"{c1['peaks'][1]:7.3f}   {c1['peaks'][0]:7.3f}")

    # 4. slotted flap: circulation effect and dumping
    print("\nslotted flap (chord 0.30, 30 deg): main-element cl and its TE-region cp, alone vs with flap")
    fl = slotted_flap_configuration(30.0)
    for a in (0.0, 5.0):
        c0 = element_peaks(panel_basis([clean]), a)
        c1 = element_peaks(panel_basis(fl), a)
        te0 = c0["result"]["parts"][0]["cp"][[0, -1]].mean()
        te1 = c1["result"]["parts"][0]["cp"][[0, -1]].mean()
        print(f"  alpha {a:4.1f}: main cl {c0['cl']:.3f} -> {c1['cls'][0]:.3f} (+flap cl {c1['cls'][1]:.3f}, total {c1['cl']:.3f}); "
              f"main TE cp {te0:.3f} -> {te1:.3f}; main LE peak {c0['peaks'][0]:.3f} -> {c1['peaks'][0]:.3f}")

    # 5. cl_max by the peak criterion
    configs = {"clean": ([clean], None),
               "plain flap 0.25c, 10 deg": ([flapped_section("2412", 0.25, 10.0, 100)], None),
               "plain flap 0.25c, 20 deg": ([flapped_section("2412", 0.25, 20.0, 100)], None),
               "slat 0.15c, 30 deg": (slat, None),
               "slotted flap 0.30c, 20 deg": (slotted_flap_configuration(20.0), [0]),
               "slotted flap 0.30c, 30 deg": (slotted_flap_configuration(30.0), [0]),
               "slat + slotted flap 30 deg": ([slat[0], clean, slotted_flap_configuration(30.0)[1]], [0, 1])}
    print(f"\ncl_max by the suction-peak criterion (cp_crit = {cp_crit:.3f}); thin-airfoil lift shift for the plain flaps")
    res = {}
    for name, (el, chk) in configs.items():
        r = clmax_by_peak(el, cp_crit, chk)
        res[name] = r
        extra = ""
        if name.startswith("plain"):
            dl = float(name.split(",")[1].split()[0])
            extra = f"  (thin-airfoil shift 2 pi tau delta = {flap_increments(0.25, np.radians(dl))['dcl']:.3f})"
        print(f"  {name:28s} alpha_max {r['alpha']:6.2f} deg  cl_max {r['cl']:.3f}  dcl_max {r['cl'] - res['clean']['cl']:+.3f}  "
              f"limited by element {r['element']}, peaks {np.round(r['peaks'], 2)}{extra}")

    # 6. the trainer
    m_, S = TRAINER["mass"], TRAINER["S"]
    print(f"\nreference trainer, {m_:.0f} kg, S = {S} m^2, approach {TRAINER['V_approach']} m/s; "
          f"CL_max = {K_3D} (cl_max + dcl_max x {FLAPPED_SPAN})")
    need = 2 * m_ * G / (RHO * S * (TRAINER["V_approach"] / 1.3) ** 2)
    print(f"  CL_max needed for V_app = 1.3 V_s: {need:.3f}")
    for name in ("clean", "plain flap 0.25c, 20 deg", "slotted flap 0.30c, 30 deg", "slat + slotted flap 30 deg"):
        dcl = res[name]["cl"] - res["clean"]["cl"]
        CL = wing_clmax(CLMAX_CLEAN, dcl)
        Vs = stall_speed(m_, S, CL)
        print(f"  {name:28s} CL_max {CL:.3f}  V_s {Vs:.2f} m/s  V_app/V_s {TRAINER['V_approach'] / Vs:.3f}")
    print("  sensitivity to the assumed clean cl_max (the criterion is recalibrated each time):")
    sf = slotted_flap_configuration(30.0)
    for clm in (1.4, 1.6, 1.8):
        _, cpc = critical_peak(clean, clm)
        r = clmax_by_peak(sf, cpc, [0])
        for lab, CL in (("clean", K_3D * clm), ("slotted 30", wing_clmax(clm, r["cl"] - clm))):
            Vs = stall_speed(m_, S, CL)
            print(f"    cl_max {clm}: cp_crit {cpc:.3f}, {lab:10s} CL_max {CL:.3f}, V_s {Vs:.2f} m/s, V_app/V_s {TRAINER['V_approach'] / Vs:.3f}")

    # figure
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    Es = np.linspace(0.02, 0.98, 60)
    ax[0].plot(Es, [flap_effectiveness(e) for e in Es], label="tau, thin airfoil")
    ax[0].plot(E, dv[0] / (2 * np.pi), "o", label="panel, t/c = 0.02")
    ax[0].set_xlabel("flap chord ratio cf/c"); ax[0].set_ylabel("tau = -d alpha_L0 / d delta")
    ax[0].set_title("plain-flap effectiveness"); ax[0].legend()
    for a in (10.0,):
        c0 = element_peaks(panel_basis([clean]), a)["result"]["parts"][0]
        c1 = element_peaks(panel_basis(slat), a)["result"]["parts"]
        ax[1].plot(c0["xc"], c0["cp"], label=f"clean, alpha = {a:.0f}")
        ax[1].plot(c1[1]["xc"], c1[1]["cp"], label="main with slat")
        ax[1].plot(c1[0]["xc"], c1[0]["cp"], label="slat")
    ax[1].invert_yaxis(); ax[1].set_xlabel("x/c"); ax[1].set_ylabel("cp"); ax[1].set_title("Smith's slat effect")
    ax[1].legend(fontsize=8)
    for X, Y in slat + [slotted_flap_configuration(30.0)[1]]:
        ax[2].plot(X, Y, "-")
    ax[2].set_aspect("equal"); ax[2].set_title("slat + main + slotted flap (panel model)")
    fig.tight_layout()
    fig.savefig(out_dir / "high_lift.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
