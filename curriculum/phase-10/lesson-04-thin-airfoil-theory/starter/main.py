"""P10.04 Thin Airfoil Theory: starter code.

A thin airfoil is replaced by a vortex sheet on its chord line whose strength makes the camber line
a streamline. With the Glauert substitution x = c (1 - cos theta)/2 the sheet strength is a Fourier
series (`P02.09`), and three of its coefficients give everything a designer needs:

    A0 = alpha - (1/pi) int_0^pi dz/dx d theta,
    An = (2/pi) int_0^pi dz/dx cos(n theta) d theta,
    cl = pi (2 A0 + A1) = 2 pi (alpha - alpha_L0),
    alpha_L0 = -(1/pi) int_0^pi dz/dx (cos theta - 1) d theta,
    cm_c/4 = (pi/4)(A2 - A1),     x_cp/c = 1/4 (1 + pi (A1 - A2)/cl).

This module evaluates those integrals by Gauss-Legendre quadrature for ANY camber-line slope,
split at the slope's kinks or jumps (the NACA 4-digit kink at x = p, the 5-digit kink at x = r, a
flap hinge), so each piece is smooth and the quadrature converges fast. A plain flap is just a
camber line with a slope jump at the hinge; its effectiveness tau = 1 - (theta_h - sin theta_h)/pi
follows.

Independent routes: closed forms (parabolic camber, flat plate, flap); a discrete-vortex (lumped
vortex) model of the camber line with its own convergence study; and a compact copy of the
`P10.05` vortex panel method, run on NACA sections of decreasing thickness.

Sign conventions (Phase 10): alpha nose-up positive, cm nose-up positive about the quarter chord,
x from the leading edge in chords, a flap deflection delta positive trailing-edge DOWN.
Circulation positive clockwise, L' = rho V Gamma.

Run:  python main.py [out_dir]   (prints the lesson's numbers, writes thin_airfoil.png)
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable

import numpy as np

TRAINER = dict(mass=750.0, S=12.0, c=1.2, V_cruise=55.0, V_approach=30.0, rho=1.225, g=9.81)

SlopeFunc = Callable[[np.ndarray], np.ndarray]


# ---------------------------------------------------------------- given: camber lines (from P10.02)

def naca4_camber_slope(m: float, p: float) -> tuple[SlopeFunc, list[float]]:
    """dz/dx of the NACA 4-digit camber line and its kink location(s) (from P10.02, given)."""
    if m == 0 or p == 0:
        return (lambda x: np.zeros_like(np.asarray(x, float))), []
    f = lambda x: np.where(np.asarray(x) <= p, 2 * m / p ** 2 * (p - np.asarray(x)),
                           2 * m / (1 - p) ** 2 * (p - np.asarray(x)))
    return f, [p]


def naca4_camber(x: np.ndarray, m: float, p: float) -> np.ndarray:
    """z(x) of the NACA 4-digit camber line (from P10.02, given)."""
    x = np.asarray(x, float)
    if m == 0 or p == 0:
        return np.zeros_like(x)
    return np.where(x <= p, m / p ** 2 * (2 * p * x - x ** 2), m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x ** 2))


def naca5_camber_slope(r: float, k1: float) -> tuple[SlopeFunc, list[float]]:
    """dz/dx of the standard NACA 5-digit camber line (cubic to x = r, then straight) (from P10.02, given).
    The 230 line: r = 0.2025, k1 = 15.957."""
    f = lambda x: np.where(np.asarray(x) <= r, k1 / 6 * (3 * np.asarray(x) ** 2 - 6 * r * np.asarray(x) + r ** 2 * (3 - r)),
                           -k1 * r ** 3 / 6)
    return f, [r]


# ---------------------------------------------------------------- thin-airfoil theory

def theta_of_x(x: float) -> float:
    """Glauert variable: x = (1 - cos theta)/2  =>  theta = arccos(1 - 2x)."""
    raise NotImplementedError


def glauert_integrals(dzdx: SlopeFunc, breaks: list[float] = (), n_max: int = 3,
                      n_gauss: int = 48) -> np.ndarray:
    """I_n = int_0^pi dz/dx(x(theta)) cos(n theta) d theta for n = 0..n_max, by Gauss-Legendre
    quadrature with n_gauss points on each sub-interval between the breaks (chordwise positions where
    the slope has a kink or a jump, converted to theta)."""
    raise NotImplementedError


def fourier_coefficients(dzdx: SlopeFunc, alpha: float, breaks: list[float] = (), n_max: int = 3,
                         n_gauss: int = 48) -> np.ndarray:
    """A0..A_{n_max}: A0 = alpha - I_0/pi, An = 2 I_n / pi."""
    raise NotImplementedError


def thin_airfoil(dzdx: SlopeFunc, alpha: float, breaks: list[float] = (), n_gauss: int = 48) -> dict:
    """Section coefficients of a thin cambered airfoil at alpha (radians):
    cl = pi (2 A0 + A1); alpha_L0 = -(I_1 - I_0)/pi; cm_c/4 = pi/4 (A2 - A1);
    cm_le = -(cl/4 + pi/4 (A1 - A2)); x_cp/c = -cm_le/cl (inf when cl = 0)."""
    raise NotImplementedError


def flap_camber_slope(E: float, delta: float) -> tuple[SlopeFunc, list[float]]:
    """A plain flap of chord E c deflected delta (trailing edge down, small angle) on a flat plate,
    as a camber line: dz/dx = 0 ahead of the hinge x_h = 1 - E and -delta behind it."""
    raise NotImplementedError


def flap_effectiveness(E: float) -> float:
    """Thin-airfoil flap effectiveness tau = d alpha_L0 / d delta (magnitude):
    tau = 1 - (theta_h - sin theta_h)/pi, theta_h = arccos(1 - 2 x_h), x_h = 1 - E."""
    raise NotImplementedError


def flap_cm_derivative(E: float) -> float:
    """d cm_c/4 / d delta = -(1/2) sin(theta_h) (1 - cos(theta_h)) for the plain flap."""
    raise NotImplementedError


def combine_slopes(*parts: tuple[SlopeFunc, list[float]]) -> tuple[SlopeFunc, list[float]]:
    """Sum of camber-line slopes (camber + flap), with the union of their breaks (given)."""
    fs = [p[0] for p in parts]
    br = sorted(set(b for p in parts for b in p[1]))
    return (lambda x: sum(f(x) for f in fs)), br


# ---------------------------------------------------------------- independent route 1: lumped vortices

def lumped_vortex(dzdx: SlopeFunc, alpha: float, n: int, spacing: str = "uniform") -> dict:
    """Discrete-vortex (lumped-vortex) model of the linearised camber line: n panels between nodes
    x_k on the chord, a clockwise point vortex Gamma_j at each panel's quarter point and a collocation
    point at its three-quarter point, both on the chord line (the thin-airfoil linearisation). Flow
    tangency at collocation point i, with the camber slope there:

        sum_j -Gamma_j / (2 pi (x_ci - x_vj)) = V (dz/dx(x_ci) - alpha),   V = 1.

    cl = 2 sum Gamma; cm_c/4 = -2 sum Gamma_j (x_vj - 1/4). Exact cl for a flat plate and for a
    parabolic camber line at any n."""
    raise NotImplementedError


def observed_order(errors: np.ndarray, ratio: float = 2.0) -> np.ndarray:
    """Observed orders log(e_k/e_{k+1})/log(ratio) (given)."""
    e = np.asarray(errors, float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


# ---------------------------------------------------------------- independent route 2: compact panel method (P10.05)

def naca4_coords(m: float, p: float, t: float, n: int = 100) -> tuple[np.ndarray, np.ndarray]:
    """Closed-trailing-edge NACA 4-digit coordinates in Selig order, n panels per surface (from P10.02, given)."""
    x = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n + 1)))
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1036 * x ** 4)
    yc = naca4_camber(x, m, p)
    th = np.arctan(naca4_camber_slope(m, p)[0](x))
    xu, yu, xl, yl = x - yt * np.sin(th), yc + yt * np.cos(th), x + yt * np.sin(th), yc - yt * np.cos(th)
    return np.concatenate([xu[::-1], xl[1:]]), np.concatenate([yu[::-1], yl[1:]])


def panel_cl_cm(x: np.ndarray, y: np.ndarray, alpha: float) -> tuple[float, float]:
    """Compact linear-strength vortex panel method of P10.05 (given): flow tangency at panel
    midpoints, Kutta condition gamma_0 + gamma_N = 0; returns (cl from Gamma, cm_c/4 from the
    surface pressure with the surface speed -gamma). Chord along x from 0 to 1."""
    dx, dy = np.diff(x), np.diff(y)
    L = np.hypot(dx, dy)
    tx, ty = dx / L, dy / L
    nx, ny = ty, -tx
    xc, yc = 0.5 * (x[:-1] + x[1:]), 0.5 * (y[:-1] + y[1:])
    N = L.size
    rx, ry = xc[:, None] - x[None, :-1], yc[:, None] - y[None, :-1]
    xi, eta = rx * tx + ry * ty, -rx * ty + ry * tx
    k = np.arange(N)
    eta[k, k] = -0.0
    beta = np.arctan2(eta, xi - L) - np.arctan2(eta, xi)
    beta[k, k] = -np.pi
    lnr = 0.5 * np.log((xi ** 2 + eta ** 2) / ((xi - L) ** 2 + eta ** 2))
    I3, I4 = -eta * lnr + xi * beta, -L + eta * beta + xi * lnr
    c = 1 / (2 * np.pi)
    ua, ub, va, vb = c * (beta - I3 / L), c * I3 / L, -c * (lnr - I4 / L), -c * I4 / L
    rot = lambda ul, vl: (ul * tx - vl * ty, ul * ty + vl * tx)
    (uA, vA), (uB, vB) = rot(ua, va), rot(ub, vb)
    U, V = np.zeros((N, N + 1)), np.zeros((N, N + 1))
    U[:, :N] += uA; U[:, 1:] += uB; V[:, :N] += vA; V[:, 1:] += vB
    A = np.zeros((N + 1, N + 1))
    A[:N] = U * nx[:, None] + V * ny[:, None]
    A[N, 0] = A[N, N] = 1.0
    rhs = np.concatenate([-(np.cos(alpha) * nx + np.sin(alpha) * ny), [0.0]])
    g = np.linalg.solve(A, rhs)
    cp = 1 - (0.5 * (g[:-1] + g[1:])) ** 2
    Fx, Fy = -cp * nx * L, -cp * ny * L
    cm = float(np.sum(yc * Fx - (xc - 0.25) * Fy))
    return float(2 * np.sum(L * 0.5 * (g[:-1] + g[1:]))), cm


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save thin_airfoil.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    t = TRAINER
    W = t["mass"] * t["g"]
    cl_cr = W / (0.5 * t["rho"] * t["V_cruise"] ** 2 * t["S"])
    cl_ap = W / (0.5 * t["rho"] * t["V_approach"] ** 2 * t["S"])

    f2412, br = naca4_camber_slope(0.02, 0.4)
    r = thin_airfoil(f2412, 0.0, br)
    print(f"NACA 2412 camber line: A1 = {r['A'][1]:.5f}, A2 = {r['A'][2]:.5f}, alpha_L0 = {np.degrees(r['alpha_L0']):.3f} deg, "
          f"cm_c/4 = {r['cm_c4']:.5f}")
    a_cr = cl_cr / (2 * np.pi) + r["alpha_L0"]
    rc = thin_airfoil(f2412, a_cr, br)
    print(f"trainer cruise cl = {cl_cr:.4f} at alpha = {np.degrees(a_cr):.3f} deg (section); x_cp/c = {rc['x_cp']:.4f}; "
          f"approach cl = {cl_ap:.4f} needs alpha = {np.degrees(cl_ap / (2 * np.pi) + r['alpha_L0']):.2f} deg without flap")
    for a_deg in (0.0, 2.0, 4.0, 8.0):
        ra = thin_airfoil(f2412, np.radians(a_deg), br)
        print(f"   alpha = {a_deg:4.1f} deg: cl = {ra['cl']:.4f}, cm_c/4 = {ra['cm_c4']:.5f}, x_cp/c = {ra['x_cp']:.4f}")
    # quadrature convergence (split vs unsplit)
    print("quadrature: |alpha_L0(n) - alpha_L0(200)| with and without the split at x = p:")
    ref = thin_airfoil(f2412, 0.0, br, 200)["alpha_L0"]
    for ng in (4, 8, 16, 32):
        e_s = abs(thin_airfoil(f2412, 0.0, br, ng)["alpha_L0"] - ref)
        e_u = abs(thin_airfoil(f2412, 0.0, [], 2 * ng)["alpha_L0"] - ref)
        print(f"   {ng:3d} points per piece: split {e_s:.2e}; unsplit with {2 * ng} points {e_u:.2e}")

    f230, br230 = naca5_camber_slope(0.2025, 15.957)
    r230 = thin_airfoil(f230, 0.0, br230)
    print(f"NACA 230 camber line: alpha_L0 = {np.degrees(r230['alpha_L0']):.3f} deg, cm_c/4 = {r230['cm_c4']:.5f}, "
          f"ideal cl (A0 = 0) = {np.pi * r230['A'][1]:.4f}")
    f4412, br4 = naca4_camber_slope(0.04, 0.4)
    r44 = thin_airfoil(f4412, 0.0, br4)
    print(f"NACA 4412 camber line: alpha_L0 = {np.degrees(r44['alpha_L0']):.3f} deg, cm_c/4 = {r44['cm_c4']:.5f}")

    # flaps
    print("\nplain flap, thin-airfoil theory: E = flap chord / c")
    for E in (0.1, 0.2, 0.25, 0.3, 0.4, 0.5):
        print(f"   E = {E:.2f}: tau = {flap_effectiveness(E):.4f}, dcm_c/4/d delta = {flap_cm_derivative(E):.4f} /rad")
    E = 0.25
    tau = flap_effectiveness(E)
    for a_deg in (4.0, 6.0):
        delta = (cl_ap / (2 * np.pi) - np.radians(a_deg) + r["alpha_L0"]) / tau
        sl, brs = combine_slopes((f2412, br), flap_camber_slope(E, delta))
        rf = thin_airfoil(sl, np.radians(a_deg), brs)
        print(f"approach cl {cl_ap:.4f} at alpha = {a_deg} deg with a 25 % plain flap: delta = {np.degrees(delta):.2f} deg "
              f"(check cl = {rf['cl']:.4f}); cm_c/4 = {rf['cm_c4']:.4f}")

    # lumped-vortex convergence (2412 camber, alpha = 4 deg)
    al = np.radians(4.0)
    ex = thin_airfoil(f2412, al, br)
    ns = [10, 20, 40, 80, 160, 320]
    e_cl, e_cm, e_clc = [], [], []
    for n in ns:
        lv = lumped_vortex(f2412, al, n)
        lc = lumped_vortex(f2412, al, n, "cosine")
        e_cl.append(abs(lv["cl"] - ex["cl"])); e_cm.append(abs(lv["cm_c4"] - ex["cm_c4"])); e_clc.append(abs(lc["cl"] - ex["cl"]))
    print(f"\nlumped vortex vs Glauert (NACA 2412 camber, alpha = 4 deg, cl = {ex['cl']:.5f}, cm_c/4 = {ex['cm_c4']:.5f}):")
    for n, a, b, c_ in zip(ns, e_cl, e_cm, e_clc):
        print(f"   n = {n:4d}: |dcl| {a:.3e}, |dcm| {b:.3e}  (cosine spacing |dcl| {c_:.3e})")
    print(f"   observed order cl: {np.array2string(observed_order(e_cl), precision=2)}; "
          f"cm: {np.array2string(observed_order(e_cm), precision=2)}; cosine cl: {np.array2string(observed_order(e_clc), precision=2)}")
    fp = lumped_vortex(lambda x: 0 * x, al, 7)
    pb = lumped_vortex(lambda x: 0.4 * (1 - 2 * x), al, 7)
    print(f"   parabolic camber z = 0.4 x (1 - x), n = 7: cl = {pb['cl']:.10f} (thin airfoil "
          f"{thin_airfoil(lambda x: 0.4 * (1 - 2 * x), al)['cl']:.10f})")
    print(f"   flat plate, n = 7: cl = {fp['cl']:.10f} (2 pi alpha = {2 * np.pi * al:.10f}), cm_c/4 = {fp['cm_c4']:.1e}")

    # panel method at decreasing thickness
    print("\npanel method (N = 240) on NACA 24tt, alpha = 4 deg, against thin-airfoil theory:")
    thick = [0.12, 0.06, 0.03, 0.015]
    pcl, pcm = [], []
    for tt in thick:
        x, y = naca4_coords(0.02, 0.4, tt, 120)
        c1, m1 = panel_cl_cm(x, y, al)
        pcl.append(c1); pcm.append(m1)
        print(f"   t/c = {tt:.3f}: cl {c1:.5f} (thin {ex['cl']:.5f}, +{(c1 / ex['cl'] - 1) * 100:.2f} %), "
              f"cm_c/4 {m1:.5f} (thin {ex['cm_c4']:.5f})")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    aa = np.linspace(-6, 12, 50)
    ax[0].plot(aa, [thin_airfoil(f2412, np.radians(a), br)["cl"] for a in aa], label="NACA 2412, thin airfoil")
    for d in (10, 20):
        sl, brs = combine_slopes((f2412, br), flap_camber_slope(0.25, np.radians(d)))
        ax[0].plot(aa, [thin_airfoil(sl, np.radians(a), brs)["cl"] for a in aa], label=f"+ 25 % flap, {d} deg")
    ax[0].axhline(cl_cr, color="k", lw=0.6, ls="--"); ax[0].axhline(cl_ap, color="k", lw=0.6, ls=":")
    ax[0].set_xlabel("alpha (deg)"); ax[0].set_ylabel("cl"); ax[0].legend(fontsize=8); ax[0].set_title("lift lines")
    Es = np.linspace(0.02, 0.98, 60)
    ax[1].plot(Es, [flap_effectiveness(e) for e in Es], label="tau")
    ax[1].plot(Es, [-flap_cm_derivative(e) for e in Es], label="-dcm_c/4/d delta")
    ax[1].set_xlabel("flap chord ratio E"); ax[1].legend(); ax[1].set_title("plain flap, thin-airfoil theory")
    h = 1 / np.array(ns)
    ax[2].loglog(h, e_cl, "o-", label="lumped vortex cl (uniform)")
    ax[2].loglog(h, e_cm, "s-", label="lumped vortex cm (uniform)")
    ax[2].loglog(h, e_cl[0] * (h / h[0]) ** 2, "k--", lw=0.8, label="slope 2")
    ax[2].set_xlabel("1/n"); ax[2].set_ylabel("error vs Glauert"); ax[2].legend(fontsize=8)
    ax[2].set_title("discrete-vortex convergence")
    fig.tight_layout()
    fig.savefig(out_dir / "thin_airfoil.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'thin_airfoil.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
