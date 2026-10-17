"""P10.05 Vortex Panel Method: reference solution.

A 2-D linear-strength vortex panel method from scratch:

* the airfoil surface (Selig order: trailing edge -> upper -> leading edge -> lower -> trailing
  edge) is cut into N straight panels with N + 1 nodes;
* on each panel the vortex-sheet strength gamma varies linearly between its two end nodes; the
  N + 1 node values are the unknowns;
* flow tangency (zero normal velocity) at the N panel midpoints, plus the Kutta condition
  gamma_0 + gamma_N = 0 (equal speeds leaving the trailing edge on both surfaces), close the
  (N + 1) x (N + 1) linear system (`P01.13`, `P03.08`);
* the surface velocity, Cp, and cl two ways (from the circulation, L' = rho V Gamma, and from
  integrating Cp), cd from Cp (zero in exact inviscid flow: a measure of discretisation error),
  and cm about the quarter chord.

It is validated against the exact Karman-Trefftz solution of `P10.03` (copied here as given
helpers), with the observed order of convergence of cl and Cp, and against thin-airfoil theory
(`P10.04`, a given helper here) as the thickness goes to zero.

Sign convention (Phase 10): gamma and Gamma are positive CLOCKWISE, so L' = rho V Gamma > 0 for a
lifting airfoil in a left-to-right stream. A clockwise point vortex Gamma at the origin induces
u = Gamma y / (2 pi r^2), v = -Gamma x / (2 pi r^2).

Units: chord-scaled geometry, V = 1 unless stated; alpha in radians unless a name says _deg.

Run:  python solution.py [out_dir]   (prints the lesson's numbers, writes vortex_panel.png)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

TRAINER = dict(mass=750.0, S=12.0, c=1.2, V=55.0, rho=1.225, g=9.81)


# ---------------------------------------------------------------- given: NACA 4-digit generator (from P10.02)

def naca4(code: str, n: int = 80, closed_te: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """NACA 4-digit coordinates in Selig order, n panels per surface on cosine spacing, thickness
    perpendicular to the camber line, closed trailing edge by default (compact copy of P10.02)."""
    m, p, t = int(code[0]) / 100, int(code[1]) / 10, int(code[2:]) / 100
    return naca4_mpt(m, p, t, n, closed_te)


def naca4_mpt(m: float, p: float, t: float, n: int = 80, closed_te: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """As ``naca4`` but from the parameters (m, p, t), so t can be any value (thin sections)."""
    x = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n + 1)))
    a4 = -0.1036 if closed_te else -0.1015
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 + a4 * x ** 4)
    if m > 0:
        yc = np.where(x <= p, m / p ** 2 * (2 * p * x - x ** 2), m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x ** 2))
        dyc = np.where(x <= p, 2 * m / p ** 2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))
    else:
        yc, dyc = np.zeros_like(x), np.zeros_like(x)
    th = np.arctan(dyc)
    xu, yu, xl, yl = x - yt * np.sin(th), yc + yt * np.cos(th), x + yt * np.sin(th), yc - yt * np.cos(th)
    return np.concatenate([xu[::-1], xl[1:]]), np.concatenate([yu[::-1], yl[1:]])


# ---------------------------------------------------------------- given: exact Karman-Trefftz flow (from P10.03)

def circle_geometry(b: float, mu: complex) -> tuple[float, float]:
    """Radius a and angle beta of the circle through zeta = b centred at mu (b - mu = a e^{-i beta}) (P10.03)."""
    d = b - mu
    return float(abs(d)), float(-np.angle(d))


def karman_trefftz(zeta: np.ndarray, b: float, tau: float) -> np.ndarray:
    """Karman-Trefftz map z = n b (1 + r^n)/(1 - r^n), r = (zeta - b)/(zeta + b), n = 2 - tau/pi (P10.03)."""
    n = 2.0 - tau / np.pi
    r = (np.asarray(zeta, complex) - b) / (np.asarray(zeta, complex) + b)
    return n * b * (1 + r ** n) / (1 - r ** n)


def karman_trefftz_dzdzeta(zeta: np.ndarray, b: float, tau: float) -> np.ndarray:
    """dz/dzeta = 4 n^2 b^2 r^{n-1} / ((1 - r^n)^2 (zeta + b)^2) (P10.03)."""
    n = 2.0 - tau / np.pi
    zeta = np.asarray(zeta, complex)
    r = (zeta - b) / (zeta + b)
    return 4 * n ** 2 * b ** 2 * r ** (n - 1) / ((1 - r ** n) ** 2 * (zeta + b) ** 2)


def kt_exact(b: float, mu: complex, tau: float, alpha: float, theta: np.ndarray) -> dict:
    """Exact surface points, speed, Cp at circle angles theta, plus Gamma, chord and exact cl and
    cm_c/4 (Kutta circulation Gamma = 4 pi a sin(alpha + beta), V = 1) (P10.03). cm about the
    quarter-chord point on the TE-LE chord line, by Blasius' moment theorem evaluated numerically."""
    a, beta = circle_geometry(b, mu)
    G = 4 * np.pi * a * np.sin(alpha + beta)
    zeta = mu + a * np.exp(1j * np.asarray(theta, float))
    z = karman_trefftz(zeta, b, tau)
    wc = np.exp(-1j * alpha) - a ** 2 * np.exp(1j * alpha) / (zeta - mu) ** 2 + 1j * G / (2 * np.pi * (zeta - mu))
    q = np.abs(wc / karman_trefftz_dzdzeta(zeta, b, tau))
    # geometry: TE and LE (farthest point), on a fine grid
    tf = -beta + 2 * np.pi * (np.arange(8192) + 0.5) / 8192
    zf = karman_trefftz(mu + a * np.exp(1j * tf), b, tau)
    z_te = karman_trefftz(np.array([b + 0j]), b, tau)[0]
    z_le = zf[np.argmax(np.abs(zf - z_te))]
    chord = float(abs(z_te - z_le))
    # exact cm about the quarter chord: Blasius M = Re(-rho/2 oint z (dw/dz)^2 dz) (counter-clockwise
    # positive), evaluated on a circle of radius 2a in the zeta plane, mapped
    tt = 2 * np.pi * (np.arange(4096) + 0.5) / 4096
    zt = mu + 2 * a * np.exp(1j * tt)
    dz_ = karman_trefftz_dzdzeta(zt, b, tau)
    wz = (np.exp(-1j * alpha) - a ** 2 * np.exp(1j * alpha) / (zt - mu) ** 2 + 1j * G / (2 * np.pi * (zt - mu))) / dz_
    zq = z_le + 0.25 * (z_te - z_le)
    zz = karman_trefftz(zt, b, tau) - zq
    M_ccw = np.real(-0.5 * np.sum(zz * wz ** 2 * dz_ * 1j * 2 * a * np.exp(1j * tt)) * (2 * np.pi / 4096))
    return dict(z=z, q=q, cp=1 - q ** 2, Gamma=G, chord=chord, cl=2 * G / chord, cm=float(-M_ccw / (0.5 * chord ** 2)),
                z_te=z_te, z_le=z_le, a=a, beta=beta)


def kt_panel_nodes(b: float, mu: complex, tau: float, n_panels: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """N + 1 panel nodes on the exact Karman-Trefftz surface at equal circle-angle steps from the
    trailing edge (Selig order; the uniform circle angle clusters nodes at both edges). Returns x,
    y and the circle angles of the panel midpoints (for the exact Cp there) (given)."""
    a, beta = circle_geometry(b, mu)
    th = -beta + 2 * np.pi * np.arange(n_panels + 1) / n_panels
    z = karman_trefftz(mu + a * np.exp(1j * th), b, tau)
    z[0] = z[-1] = karman_trefftz(np.array([b + 0j]), b, tau)[0]
    return z.real, z.imag, 0.5 * (th[:-1] + th[1:])


def thin_airfoil_naca4(m: float, p: float, alpha: float) -> tuple[float, float, float]:
    """Thin-airfoil cl, cm_c/4 and alpha_L0 of the NACA 4-digit camber line, by Gauss-Legendre
    quadrature of the Glauert integrals split at the kink x = p (P10.04, given)."""
    if m == 0:
        return 2 * np.pi * alpha, 0.0, 0.0
    tp = np.arccos(1 - 2 * p)
    xg, wg = np.polynomial.legendre.leggauss(64)
    I = np.zeros(3)
    for lo, hi in ((0.0, tp), (tp, np.pi)):
        th = 0.5 * (hi - lo) * (xg + 1) + lo
        x = 0.5 * (1 - np.cos(th))
        dz = np.where(x <= p, 2 * m / p ** 2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))
        w = 0.5 * (hi - lo) * wg
        I += [np.sum(w * dz), np.sum(w * dz * np.cos(th)), np.sum(w * dz * np.cos(2 * th))]
    A0, A1, A2 = alpha - I[0] / np.pi, 2 * I[1] / np.pi, 2 * I[2] / np.pi
    aL0 = (I[0] - I[1]) / np.pi            # -(1/pi) int dz/dx (cos th - 1) d th
    return float(np.pi * (2 * A0 + A1)), float(np.pi / 4 * (A2 - A1)), float(aL0)


# ---------------------------------------------------------------- the panel method

def panel_geometry(x: np.ndarray, y: np.ndarray) -> dict:
    """Panels between consecutive nodes: midpoints (xc, yc), lengths, unit tangents (along the node
    order) and outward unit normals n = (t_y, -t_x) (outward for counter-clockwise Selig order)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    dx, dy = np.diff(x), np.diff(y)
    L = np.hypot(dx, dy)
    tx, ty = dx / L, dy / L
    return dict(xc=0.5 * (x[:-1] + x[1:]), yc=0.5 * (y[:-1] + y[1:]), L=L, tx=tx, ty=ty, nx=ty, ny=-tx,
                x0=x[:-1], y0=y[:-1])


def linear_vortex_influence(px: np.ndarray, py: np.ndarray, geo: dict,
                            outside_self: bool = False) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Velocity at points (px, py) induced by each panel j carrying a clockwise vortex sheet whose
    strength falls linearly from 1 at its start node to 0 at its end node (a), or rises from 0 to 1
    (b). Returns (ua, va, ub, vb), arrays (n_points, n_panels), global components.

    In panel coordinates (xi along the panel from its start, eta to the left), with
    D = (xi - s)^2 + eta^2, the clockwise element gamma(s) ds induces
    du = gamma eta / (2 pi D), dv = -gamma (xi - s) / (2 pi D), and with gamma = g_a (1 - s/L) + g_b s/L:

        I1 = int eta/D ds = beta (angle subtended),  I2 = int (xi - s)/D ds = ln(r1/r2),
        I3 = int s eta/D ds = eta ln(r2/r1) + xi beta,  I4 = int s (xi - s)/D ds = -L + eta beta + xi ln(r1/r2),
        u = [g_a (I1 - I3/L) + g_b I3/L]/(2 pi),  v = -[g_a (I2 - I4/L) + g_b I4/L]/(2 pi).

    If ``outside_self``, the points are the panel midpoints and the self-induced term is taken on the
    outer side (eta -> 0 from below, beta = -pi)."""
    px = np.asarray(px, float)[:, None]
    py = np.asarray(py, float)[:, None]
    tx, ty, L = geo["tx"][None, :], geo["ty"][None, :], geo["L"][None, :]
    rx, ry = px - geo["x0"][None, :], py - geo["y0"][None, :]
    xi = rx * tx + ry * ty
    eta = -rx * ty + ry * tx
    if outside_self:
        k = np.arange(min(px.shape[0], L.shape[1]))
        eta[k, k] = -0.0
    r1sq = xi ** 2 + eta ** 2
    r2sq = (xi - L) ** 2 + eta ** 2
    beta = np.arctan2(eta, xi - L) - np.arctan2(eta, xi)
    if outside_self:
        beta[k, k] = -np.pi
    lnr = 0.5 * np.log(r1sq / r2sq)                         # ln(r1/r2)
    I1, I2 = beta, lnr
    I3 = -eta * lnr + xi * beta
    I4 = -L + eta * beta + xi * lnr
    c = 1.0 / (2 * np.pi)
    ua_l, ub_l = c * (I1 - I3 / L), c * I3 / L
    va_l, vb_l = -c * (I2 - I4 / L), -c * I4 / L
    # back to global: local x along t, local y along the left normal (-ty, tx)
    g = lambda ul, vl: (ul * tx - vl * ty, ul * ty + vl * tx)
    ua, va = g(ua_l, va_l)
    ub, vb = g(ub_l, vb_l)
    return ua, va, ub, vb


def node_influence(ua: np.ndarray, va: np.ndarray, ub: np.ndarray, vb: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Assemble panel 'start' and 'end' influences into node influences: node j gets the (a) part of
    panel j and the (b) part of panel j - 1. Returns (U, V), (n_points, N + 1)."""
    npnt, N = ua.shape
    U = np.zeros((npnt, N + 1))
    Vv = np.zeros((npnt, N + 1))
    U[:, :N] += ua
    U[:, 1:] += ub
    Vv[:, :N] += va
    Vv[:, 1:] += vb
    return U, Vv


def vortex_panel(x: np.ndarray, y: np.ndarray, alpha: float, V: float = 1.0) -> dict:
    """Solve the linear-vortex panel method on the closed polygon (x, y) in Selig order at angle of
    attack alpha (free stream V (cos alpha, sin alpha)).

    Rows 0..N-1: normal velocity at each midpoint is zero; row N: gamma_0 + gamma_N = 0 (Kutta).
    Returns geometry, node gamma, the surface velocity Vt at the midpoints (along the node order,
    Vt = -(gamma_j + gamma_{j+1})/2: the velocity jump across the sheet, the fluid inside the body
    being at rest), cp = 1 - (Vt/V)^2, Gamma (clockwise), the largest residual of the linear system,
    and the coefficients from ``panel_forces``."""
    geo = panel_geometry(x, y)
    N = geo["L"].size
    ua, va, ub, vb = linear_vortex_influence(geo["xc"], geo["yc"], geo, outside_self=True)
    U, Vv = node_influence(ua, va, ub, vb)
    A = np.zeros((N + 1, N + 1))
    rhs = np.zeros(N + 1)
    A[:N] = U * geo["nx"][:, None] + Vv * geo["ny"][:, None]
    rhs[:N] = -V * (np.cos(alpha) * geo["nx"] + np.sin(alpha) * geo["ny"])
    A[N, 0] = A[N, N] = 1.0
    gamma = np.linalg.solve(A, rhs)
    # surface velocity: the jump across the sheet. With the fluid inside the body at rest, the
    # tangential velocity just outside, along the node order, is -gamma (gamma clockwise-positive).
    Vt = -0.5 * (gamma[:-1] + gamma[1:])
    cp = 1 - (Vt / V) ** 2
    Gamma = float(np.sum(geo["L"] * 0.5 * (gamma[:-1] + gamma[1:])))
    out = dict(geo=geo, gamma=gamma, Vt=Vt, cp=cp, Gamma=Gamma, alpha=alpha, V=V,
               residual=float(np.max(np.abs(A @ gamma - rhs))))
    out.update(panel_forces(geo, cp, alpha, x, y, Gamma, V))
    return out


def panel_forces(geo: dict, cp: np.ndarray, alpha: float, x: np.ndarray, y: np.ndarray, Gamma: float,
                 V: float = 1.0) -> dict:
    """Coefficients from a panel solution. The chord is the distance from the trailing edge (mean of
    the first and last node) to the farthest node (the leading edge); the quarter chord lies on that
    line. Pressure force per panel -cp n L (per unit q); c_x, c_y in body axes rotated to
    cl_p = c_y cos(alpha) - c_x sin(alpha), cd_p = c_x cos(alpha) + c_y sin(alpha);
    cl_gamma = 2 Gamma / (V c); cm_c/4 nose-up = sum[(y_i - y_q) F_x - (x_i - x_q) F_y]/c^2."""
    x_te, y_te = 0.5 * (x[0] + x[-1]), 0.5 * (y[0] + y[-1])
    i_le = int(np.argmax(np.hypot(x - x_te, y - y_te)))
    c = float(np.hypot(x[i_le] - x_te, y[i_le] - y_te))
    xq = x[i_le] + 0.25 * (x_te - x[i_le])
    yq = y[i_le] + 0.25 * (y_te - y[i_le])
    Fx = -cp * geo["nx"] * geo["L"]
    Fy = -cp * geo["ny"] * geo["L"]
    cx, cy = Fx.sum() / c, Fy.sum() / c
    cm = float(np.sum((geo["yc"] - yq) * Fx - (geo["xc"] - xq) * Fy) / c ** 2)
    return dict(chord=c, cl_p=float(cy * np.cos(alpha) - cx * np.sin(alpha)),
                cd_p=float(cx * np.cos(alpha) + cy * np.sin(alpha)), cl_gamma=2 * Gamma / (V * c), cm=cm)


def kt_convergence(b: float, mu: complex, tau: float, alpha: float, ns: list[int],
                   x_window: tuple[float, float] = (0.05, 0.95)) -> dict:
    """Panel solutions on the exact Karman-Trefftz surface for each N in ns: errors of cl (from Gamma
    and from Cp) against the exact cl, of cm_c/4, and the largest Cp error at midpoints whose chordwise
    position is inside x_window (the trailing-edge corner and the leading-edge peak converge more
    slowly and are reported separately by the full-surface maximum)."""
    ex = kt_exact(b, mu, tau, alpha, np.array([0.0]))
    out = dict(n=np.array(ns), e_clg=[], e_clp=[], e_cm=[], e_cp=[], e_cp_all=[])
    for n in ns:
        x, y, thm = kt_panel_nodes(b, mu, tau, n)
        s = vortex_panel(x, y, alpha)
        exm = kt_exact(b, mu, tau, alpha, thm)
        xs = np.real((exm["z"] - ex["z_le"]) * np.conj(ex["z_te"] - ex["z_le"])) / ex["chord"] ** 2
        win = (xs > x_window[0]) & (xs < x_window[1])
        out["e_clg"].append(abs(s["cl_gamma"] - ex["cl"]))
        out["e_clp"].append(abs(s["cl_p"] - ex["cl"]))
        out["e_cm"].append(abs(s["cm"] - ex["cm"]))
        out["e_cp"].append(float(np.max(np.abs(s["cp"][win] - exm["cp"][win]))))
        out["e_cp_all"].append(float(np.max(np.abs(s["cp"] - exm["cp"]))))
    for k in ("e_clg", "e_clp", "e_cm", "e_cp", "e_cp_all"):
        out[k] = np.array(out[k])
    out["cl_exact"], out["cm_exact"] = ex["cl"], ex["cm"]
    return out


def observed_order(errors: np.ndarray, ratio: float = 2.0) -> np.ndarray:
    """Observed orders log(e_k/e_{k+1})/log(ratio) for grids refined by ``ratio`` (given)."""
    e = np.asarray(errors, float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


def richardson(f_coarse: float, f_fine: float, order: float, ratio: float = 2.0) -> float:
    """Richardson extrapolation f_fine + (f_fine - f_coarse)/(ratio^order - 1) (given)."""
    return f_fine + (f_fine - f_coarse) / (ratio ** order - 1)


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".", n_max: int = 640) -> None:
    """Print the lesson's numbers and save vortex_panel.png (n_max: the finest panel count; the
    tests pass a smaller value to run quickly)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    tau = np.radians(10.0)
    mu = complex(-0.08, 0.06)
    al = np.radians(4.0)
    ns = [n for n in (20, 40, 80, 160, 320, 640) if n <= n_max]
    nn = [n for n in (40, 80, 160, 320, 640) if n <= n_max]
    cv = kt_convergence(1.0, mu, tau, al, ns)
    print(f"Karman-Trefftz test section (b = 1, mu = {mu}, tau = 10 deg), alpha = 4 deg: exact cl = {cv['cl_exact']:.6f}, "
          f"cm_c/4 = {cv['cm_exact']:.6f}")
    print("     N    |cl_G - cl|   |cl_p - cl|   |cm - cm_ex|  max|dCp| 5-95%   max|dCp| all")
    for i, n in enumerate(ns):
        print(f"  {n:4d}    {cv['e_clg'][i]:.3e}    {cv['e_clp'][i]:.3e}    {cv['e_cm'][i]:.3e}    {cv['e_cp'][i]:.3e}"
              f"      {cv['e_cp_all'][i]:.3e}")
    for k, name in (("e_clg", "cl from Gamma"), ("e_clp", "cl from Cp"), ("e_cm", "cm_c/4"), ("e_cp", "Cp 5-95 %"),
                    ("e_cp_all", "Cp all")):
        print(f"  observed order, {name:13s}: {np.array2string(observed_order(cv[k]), precision=2)}")

    # NACA 2412, the trainer's section
    t = TRAINER
    cl_need = t["mass"] * t["g"] / (0.5 * t["rho"] * t["V"] ** 2 * t["S"])
    res = {}
    for n in nn:
        x, y = naca4("2412", n // 2)
        s = vortex_panel(x, y, np.radians(2.0))
        res[n] = (s["cl_gamma"], s["cm"])
    print(f"\nNACA 2412 (closed TE), alpha = 2 deg: cl_Gamma by N = " +
          ", ".join(f"{n}: {v[0]:.5f}" for n, v in res.items()))
    diffs = np.abs(np.diff([res[n][0] for n in nn]))
    p_obs = observed_order(diffs)
    cl_inf = richardson(res[nn[-2]][0], res[nn[-1]][0], 2.0)
    print(f"  observed order from successive differences: {np.array2string(p_obs, precision=2)}; "
          f"Richardson (p = 2) cl = {cl_inf:.5f}")
    for n in nn:
        print(f"   N = {n:4d}: cl error vs Richardson {abs(res[n][0] - cl_inf) / cl_inf * 100:.3f} %, cm_c/4 = {res[n][1]:.5f}")
    cl_ta, cm_ta, aL0 = thin_airfoil_naca4(0.02, 0.4, np.radians(2.0))
    print(f"  thin-airfoil theory: cl = {cl_ta:.5f}, cm_c/4 = {cm_ta:.5f}, alpha_L0 = {np.degrees(aL0):.3f} deg")
    # alpha for the trainer's cruise cl (panel N = 160, linear in alpha)
    x, y = naca4("2412", 80)
    s0, s4 = vortex_panel(x, y, 0.0), vortex_panel(x, y, np.radians(4.0))
    slope = (s4["cl_gamma"] - s0["cl_gamma"]) / np.radians(4.0)
    print(f"  N = 160: lift slope {slope:.4f} /rad = {slope / (2 * np.pi):.4f} x 2 pi; alpha_L0 = "
          f"{-np.degrees(s0['cl_gamma'] / slope):.3f} deg; cruise cl {cl_need:.4f} at alpha = "
          f"{np.degrees((cl_need - s0['cl_gamma']) / slope):.3f} deg (thin airfoil "
          f"{np.degrees(cl_need / (2 * np.pi) + aL0):.3f} deg)")

    # thickness -> 0: NACA 24tt against thin-airfoil theory
    print("\nthin-airfoil limit, camber 2 % at 40 %, alpha = 4 deg, N = 320:")
    cl_ta, cm_ta, _ = thin_airfoil_naca4(0.02, 0.4, np.radians(4.0))
    thick = [0.12, 0.06, 0.03, 0.015, 0.0075]
    dcl, dcm = [], []
    for tt in thick:
        x, y = naca4_mpt(0.02, 0.4, tt, 160)
        s = vortex_panel(x, y, np.radians(4.0))
        dcl.append(s["cl_gamma"] - cl_ta)
        dcm.append(s["cm"] - cm_ta)
        print(f"   t/c = {tt:.4f}: cl {s['cl_gamma']:.5f} (thin {cl_ta:.5f}, ratio {s['cl_gamma'] / cl_ta:.4f}), "
              f"cm_c/4 {s['cm']:.5f} (thin {cm_ta:.5f})")
    print(f"   cl difference / t: {np.array2string(np.array(dcl) / np.array(thick), precision=4)}")

    # figure
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    x, y, thm = kt_panel_nodes(1.0, mu, tau, 80)
    s = vortex_panel(x, y, al)
    ex = kt_exact(1.0, mu, tau, al, np.linspace(-0.0, 2 * np.pi, 600) - circle_geometry(1.0, mu)[1])
    ax[0].plot(ex["z"].real, ex["cp"], "k-", lw=0.8, label="exact Karman-Trefftz")
    ax[0].plot(s["geo"]["xc"], s["cp"], "o", ms=3, mfc="none", label="panel, N = 80")
    ax[0].invert_yaxis(); ax[0].set_xlabel("x"); ax[0].set_ylabel("Cp"); ax[0].legend()
    ax[0].set_title("panel vs exact, alpha = 4 deg")
    h = 1.0 / np.array(ns)
    for k, lab in (("e_clg", "cl (Gamma)"), ("e_clp", "cl (Cp)"), ("e_cp", "max Cp error, 5-95 %")):
        ax[1].loglog(h, cv[k], "o-", label=lab)
    ax[1].loglog(h, cv["e_clg"][0] * (h / h[0]) ** 2, "k--", lw=0.8, label="slope 2")
    ax[1].set_xlabel("1/N"); ax[1].set_ylabel("error"); ax[1].legend(); ax[1].set_title("observed convergence")
    ax[2].plot(thick, np.array(dcl), "o-")
    ax[2].set_xlabel("t/c"); ax[2].set_ylabel("cl_panel - cl_thin"); ax[2].set_title("thin-airfoil limit (2 % camber, 4 deg)")
    fig.tight_layout()
    fig.savefig(out_dir / "vortex_panel.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'vortex_panel.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
