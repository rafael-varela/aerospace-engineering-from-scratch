"""P07.14 Heat Transfer: Radiation: reference solution.

How large must a spacecraft's radiator be to reject the heat of its electronics, how much does the
Sun or a neighbouring solar array that the radiator can "see" change the answer, and how do we know
the view factors that decide it are right? This module answers with:

* black-body radiation: Planck's spectral emissive power, its integral sigma T^4, the fraction of
  the energy in a wavelength band (why a white paint can absorb little sunlight yet emit well);
* view factors: closed forms for aligned parallel rectangles, perpendicular rectangles sharing an
  edge and coaxial parallel disks (Lienhard Table 10.3), two-dimensional strips (Table 10.2) and
  Hottel's crossed-string rule; a numerical double-area integration of the view-factor kernel
  cos(t1) cos(t2) / (pi s^2) that checks them; the six view factors of a rectangular box, whose
  rows must sum to one;
* the radiosity method for an enclosure of N gray, diffuse surfaces, each with either a known
  temperature or a known net heat flow;
* a spacecraft radiator: equilibrium temperature and required area in eclipse and in sunlight,
  and the penalty of a solar array that blocks part of its view of space.

Surface data: total emittances from Lienhard & Lienhard, A Heat Transfer Textbook 6th ed., Table 10.1;
solar absorptance and infrared emittance from Table 10.5; the solar irradiation 1361 W/m^2 and
the 5772 K effective solar temperature from section 10.6. Units: SI (K, W, m).

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes radiation.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-07/lesson-14-radiation
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

SIGMA = 5.670374419e-8        # W/(m^2 K^4), Stefan-Boltzmann constant (CODATA 2018)
H_PLANCK = 6.62607015e-34     # J s (exact, SI 2019)
C_LIGHT = 2.99792458e8        # m/s (exact)
K_B = 1.380649e-23            # J/K (exact, SI 2019)
G_SOLAR = 1361.0              # W/m^2, solar irradiation above the atmosphere at 1 AU (Lienhard section 10.6)
T_SUN = 5772.0                # K, black body giving the same irradiation (Lienhard section 10.6)
T_SPACE = 3.0                 # K, deep-space sink; its sigma T^4 is 5e-6 W/m^2, negligible
# Lienhard Table 10.5 (near 300 K): solar absorptance and infrared emittance
ZNO_WHITE_PAINT = {"alpha_s": (0.12, 0.18), "eps_ir": 0.93}
ACRYLIC_WHITE_PAINT = {"alpha_s": 0.26, "eps_ir": 0.90}
CARBON_BLACK_ACRYLIC = {"alpha_s": 0.94, "eps_ir": 0.83}
ALUMINIUM_PURE = {"alpha_s": 0.09, "eps_ir": 0.10}


# ---------------------------------------------------------------- black body

def blackbody_emissive_power(T):
    """Total black-body emissive power e_b = sigma T^4 (W/m^2)."""
    return SIGMA * np.asarray(T, dtype=float) ** 4


def planck_spectral(lam, T: float):
    """Planck's spectral emissive power e_b,lambda = 2 pi h c^2 / (lambda^5 (exp(h c / (lambda k T)) - 1)),
    in W/(m^2 m) with lambda in metres."""
    lam = np.asarray(lam, dtype=float)
    x = H_PLANCK * C_LIGHT / (lam * K_B * T)
    with np.errstate(over="ignore"):              # exp overflows to inf far below the peak: e_b -> 0
        return 2 * np.pi * H_PLANCK * C_LIGHT ** 2 / lam ** 5 / np.expm1(x)


def blackbody_fraction(lam1: float, lam2: float, T: float) -> float:
    """Fraction of sigma T^4 emitted between lam1 and lam2 (m; lam2 may be inf). With x = h c / (lambda k T)
    Planck's law integrates to F(0 -> lambda) = (15 / pi^4) int_x^inf t^3 / (e^t - 1) dt, a function of
    lambda T only; the integral is done by quadrature (the integrand is negligible beyond t = 700)."""
    c2 = H_PLANCK * C_LIGHT / K_B
    g = lambda t: t ** 3 / np.expm1(t)
    tail = lambda x: quad(g, x, max(700.0, x), limit=400, epsabs=0.0, epsrel=1e-12)[0] if x < 700 else 0.0
    x1 = c2 / (lam1 * T) if lam1 > 0 else np.inf
    x2 = 0.0 if np.isinf(lam2) else c2 / (lam2 * T)
    F_upto = lambda x: 15 / np.pi ** 4 * (tail(x) if np.isfinite(x) else 0.0)
    # F(0 -> lam2) - F(0 -> lam1); F(0 -> inf) = (15/pi^4) * pi^4/15 = 1
    F2 = 1.0 if x2 == 0.0 else F_upto(x2)
    return float(F2 - F_upto(x1))


# ---------------------------------------------------------------- closed-form view factors

def vf_parallel_rectangles(a: float, b: float, c: float) -> float:
    """Aligned, directly opposed, equal rectangles a x b a distance c apart (Lienhard Table 10.3, case 1):
    X = a/c, Y = b/c, F = (2/(pi X Y)) { ln[ ((1+X^2)(1+Y^2)/(1+X^2+Y^2))^(1/2) ] + X sqrt(1+Y^2) atan(X/sqrt(1+Y^2))
    + Y sqrt(1+X^2) atan(Y/sqrt(1+X^2)) - X atan X - Y atan Y }."""
    X, Y = a / c, b / c
    return 2 / (np.pi * X * Y) * (0.5 * np.log((1 + X ** 2) * (1 + Y ** 2) / (1 + X ** 2 + Y ** 2))
                                  + X * np.sqrt(1 + Y ** 2) * np.arctan(X / np.sqrt(1 + Y ** 2))
                                  + Y * np.sqrt(1 + X ** 2) * np.arctan(Y / np.sqrt(1 + X ** 2))
                                  - X * np.arctan(X) - Y * np.arctan(Y))


def vf_perpendicular_rectangles(l: float, w: float, h: float) -> float:
    """Two perpendicular rectangles sharing a common edge of length l (Lienhard Table 10.3, case 2): F from
    surface 1 (width w, normal to the shared edge) to surface 2 (height h). With W = w/l, H = h/l,
    F = (1/(pi W)) { W atan(1/W) + H atan(1/H) - sqrt(H^2+W^2) atan(1/sqrt(H^2+W^2))
         + (1/4) ln[ (1+W^2)(1+H^2)/(1+W^2+H^2) * (W^2(1+W^2+H^2)/((1+W^2)(W^2+H^2)))^(W^2)
                     * (H^2(1+H^2+W^2)/((1+H^2)(H^2+W^2)))^(H^2) ] }."""
    W, H = w / l, h / l
    S = np.sqrt(H ** 2 + W ** 2)
    log_term = (np.log((1 + W ** 2) * (1 + H ** 2) / (1 + W ** 2 + H ** 2))
                + W ** 2 * np.log(W ** 2 * (1 + W ** 2 + H ** 2) / ((1 + W ** 2) * (W ** 2 + H ** 2)))
                + H ** 2 * np.log(H ** 2 * (1 + H ** 2 + W ** 2) / ((1 + H ** 2) * (H ** 2 + W ** 2))))
    return (W * np.arctan(1 / W) + H * np.arctan(1 / H) - S * np.arctan(1 / S) + 0.25 * log_term) / (np.pi * W)


def vf_coaxial_disks(r1: float, r2: float, h: float) -> float:
    """Coaxial parallel disks of radii r1 (from) and r2 (to) a distance h apart (Lienhard Table 10.3, case 3):
    R1 = r1/h, R2 = r2/h, X = 1 + (1 + R2^2)/R1^2, F = (X - sqrt(X^2 - 4 (R2/R1)^2)) / 2."""
    R1, R2 = r1 / h, r2 / h
    X = 1 + (1 + R2 ** 2) / R1 ** 2
    return 0.5 * (X - np.sqrt(X ** 2 - 4 * (R2 / R1) ** 2))


def vf_parallel_strips(w: float, h: float) -> float:
    """Two directly opposed infinitely long strips of equal width w a distance h apart (Table 10.2, case 1):
    F = sqrt(1 + (h/w)^2) - h/w."""
    return np.sqrt(1 + (h / w) ** 2) - h / w


def vf_perpendicular_strips(w: float, h: float) -> float:
    """Infinitely long strip of width w to a perpendicular strip of width h sharing one edge (Table 10.2,
    case 3): F = (1 + h/w - sqrt(1 + (h/w)^2)) / 2."""
    return 0.5 * (1 + h / w - np.sqrt(1 + (h / w) ** 2))


def vf_crossed_strings(a1, b1, a2, b2) -> float:
    """Hottel's crossed-string rule for two-dimensional surfaces 1 (a1 to b1) and 2 (a2 to b2) given as
    (x, y) end points, each seeing the other unobstructed: F12 = [(crossed strings) - (uncrossed strings)]
    / (2 L1), with crossed = |a1 b2| + |b1 a2| and uncrossed = |a1 a2| + |b1 b2| (label the ends so the
    uncrossed strings do not cross)."""
    d = lambda p, q: float(np.hypot(p[0] - q[0], p[1] - q[1]))
    L1 = d(a1, b1)
    return ((d(a1, b2) + d(b1, a2)) - (d(a1, a2) + d(b1, b2))) / (2 * L1)


# ---------------------------------------------------------------- numerical view factor

def vf_numerical(o1, u1, v1, o2, u2, v2, n: int = 30) -> float:
    """View factor F12 between two planar parallelograms by double-area integration of the kernel
    cos(t1) cos(t2) / (pi s^2): surface i is the set o_i + s u_i + t v_i, 0 <= s, t <= 1. Each surface is
    divided into n x n cells and the kernel is evaluated between cell centres (midpoint rule):
    F12 = (1/A1) sum_i sum_j cos1 cos2 dA1 dA2 / (pi s^2). The normals are u x v, flipped if needed so
    that each faces the other surface's centroid; pairs facing away contribute zero (no blocking)."""
    o1, u1, v1, o2, u2, v2 = (np.asarray(z, dtype=float) for z in (o1, u1, v1, o2, u2, v2))
    g = (np.arange(n) + 0.5) / n

    def cells(o, u, v):
        S, T = np.meshgrid(g, g, indexing="ij")
        P = o + S.reshape(-1, 1) * u + T.reshape(-1, 1) * v
        nrm = np.cross(u, v)
        dA = np.linalg.norm(nrm) / n ** 2
        return P, nrm / np.linalg.norm(nrm), dA

    P1, n1, dA1 = cells(o1, u1, v1)
    P2, n2, dA2 = cells(o2, u2, v2)
    c1, c2 = o1 + 0.5 * (u1 + v1), o2 + 0.5 * (u2 + v2)
    if np.dot(n1, c2 - c1) < 0:
        n1 = -n1
    if np.dot(n2, c1 - c2) < 0:
        n2 = -n2
    R = P2[None, :, :] - P1[:, None, :]
    s2 = np.einsum("ijk,ijk->ij", R, R)
    cos1 = np.clip(R @ n1, 0.0, None) / np.sqrt(s2)
    cos2 = np.clip(-(R @ n2), 0.0, None) / np.sqrt(s2)
    A1 = dA1 * n ** 2
    return float(np.sum(cos1 * cos2 / (np.pi * s2)) * dA1 * dA2 / A1)


def box_view_factors(a: float, b: float, c: float) -> tuple[np.ndarray, np.ndarray]:
    """Areas and the 6 x 6 view-factor matrix of the inside of a rectangular box a (x) by b (y) by c (z),
    surfaces ordered x0, x1, y0, y1, z0, z1 (the faces normal to x, y, z at 0 and at the far side).
    Opposite faces use vf_parallel_rectangles, adjacent faces vf_perpendicular_rectangles; a face does not
    see itself. Returns (A, F)."""
    dims = {"x": (b, c, a), "y": (a, c, b), "z": (a, b, c)}       # (side 1, side 2, separation) per axis
    faces = [("x", 0), ("x", 1), ("y", 0), ("y", 1), ("z", 0), ("z", 1)]
    A = np.array([dims[ax][0] * dims[ax][1] for ax, _ in faces])
    span = {"x": a, "y": b, "z": c}
    F = np.zeros((6, 6))
    for i, (ai, si) in enumerate(faces):
        for j, (aj, sj) in enumerate(faces):
            if i == j:
                continue
            if ai == aj:
                p, q, sep = dims[ai]
                F[i, j] = vf_parallel_rectangles(p, q, sep)
            else:
                third = ({"x", "y", "z"} - {ai, aj}).pop()
                # common edge runs along the third axis; surface i extends along aj's axis, j along ai's
                F[i, j] = vf_perpendicular_rectangles(span[third], span[aj], span[ai])
    return A, F


# ---------------------------------------------------------------- enclosure (radiosity method)

def solve_enclosure(A, F, eps, T=None, Q=None) -> dict[str, np.ndarray]:
    """Radiosity method for an enclosure of N gray, diffuse, opaque surfaces with areas A, view factors F
    (rows sum to one) and emittances eps. For each surface give either its temperature T[i] (Q[i] None or
    nan) or its net heat flow out Q[i] (T[i] None or nan), in W. The equations are
        known T:  J_i - (1 - eps_i) sum_j F_ij J_j = eps_i sigma T_i^4
        known Q:  J_i - sum_j F_ij J_j = Q_i / A_i,
    and afterwards Q_i = A_i (J_i - sum_j F_ij J_j) and, for a known-Q surface,
    sigma T_i^4 = J_i + Q_i (1 - eps_i) / (eps_i A_i). Returns J, Q, T."""
    A = np.asarray(A, dtype=float)
    F = np.asarray(F, dtype=float)
    eps = np.asarray(eps, dtype=float)
    N = len(A)
    T = np.full(N, np.nan) if T is None else np.array([np.nan if t is None else t for t in T], dtype=float)
    Qk = np.full(N, np.nan) if Q is None else np.array([np.nan if q is None else q for q in Q], dtype=float)
    M = np.zeros((N, N))
    rhs = np.zeros(N)
    for i in range(N):
        if np.isfinite(T[i]):
            M[i] = -(1 - eps[i]) * F[i]
            M[i, i] += 1.0
            rhs[i] = eps[i] * SIGMA * T[i] ** 4
        elif np.isfinite(Qk[i]):
            M[i] = -F[i]
            M[i, i] += 1.0
            rhs[i] = Qk[i] / A[i]
        else:
            raise ValueError(f"surface {i}: give either T or Q")
    J = np.linalg.solve(M, rhs)
    Qn = A * (J - F @ J)
    Tout = T.copy()
    for i in range(N):
        if not np.isfinite(T[i]):
            Eb = J[i] + Qn[i] * (1 - eps[i]) / (eps[i] * A[i])
            Tout[i] = (Eb / SIGMA) ** 0.25
    return {"J": J, "Q": Qn, "T": Tout}


def two_surface_exchange(A1: float, A2: float, F12: float, eps1: float, eps2: float, T1: float, T2: float) -> float:
    """Net exchange from 1 to 2 in a two-surface enclosure, three resistances in series:
    Q12 = sigma (T1^4 - T2^4) / ((1-eps1)/(eps1 A1) + 1/(A1 F12) + (1-eps2)/(eps2 A2))."""
    return SIGMA * (T1 ** 4 - T2 ** 4) / ((1 - eps1) / (eps1 * A1) + 1 / (A1 * F12) + (1 - eps2) / (eps2 * A2))


def h_rad_linear(eps: float, T: float, Tsur: float) -> float:
    """Exact linearised coefficient eps sigma (T^2 + Tsur^2)(T + Tsur), so q = h_rad (T - Tsur)."""
    return eps * SIGMA * (T ** 2 + Tsur ** 2) * (T + Tsur)


def h_rad_mean(eps: float, T: float, Tsur: float) -> float:
    """Approximate coefficient 4 eps sigma T_m^3 with T_m = (T + Tsur)/2; its relative error is
    ((T - Tsur)/(T + Tsur))^2, low by that fraction."""
    return 4 * eps * SIGMA * (0.5 * (T + Tsur)) ** 3


# ---------------------------------------------------------------- spacecraft radiator

def radiator_temperature(Q: float, A: float, eps: float, alpha_s: float = 0.0, G_s: float = G_SOLAR,
                         cos_sun: float = 0.0, q_ir: float = 0.0, T_sink: float = T_SPACE) -> float:
    """Equilibrium temperature of a flat radiator of area A (one face radiating, back insulated) carrying
    Q watts from the electronics: eps sigma A (T^4 - T_sink^4) = Q + alpha_s G_s cos_sun A + eps q_ir A,
    where q_ir is any infrared irradiation from a planet (absorbed with absorptance eps, Kirchhoff)."""
    rhs = Q + alpha_s * G_s * max(cos_sun, 0.0) * A + eps * q_ir * A
    return float((rhs / (eps * SIGMA * A) + T_sink ** 4) ** 0.25)


def radiator_area(Q: float, T: float, eps: float, alpha_s: float = 0.0, G_s: float = G_SOLAR,
                  cos_sun: float = 0.0, q_ir: float = 0.0, T_sink: float = T_SPACE) -> float:
    """Area that rejects Q at temperature T: A = Q / (eps sigma (T^4 - T_sink^4) - alpha_s G_s cos_sun - eps q_ir);
    ValueError if the absorbed flux exceeds the emitted flux (no area is large enough)."""
    net = eps * SIGMA * (T ** 4 - T_sink ** 4) - alpha_s * G_s * max(cos_sun, 0.0) - eps * q_ir
    if net <= 0:
        raise ValueError("radiator_area: absorbed flux exceeds emission at this temperature")
    return Q / net


def radiator_with_array(Q: float, a: float, b: float, w_array: float, eps_r: float, eps_a: float,
                        T_array: float, T_sink: float = T_SPACE) -> dict[str, float]:
    """A radiator panel a x b (a along the shared edge) with a solar array of width w_array standing
    perpendicular to it along that edge (no sunlight on the radiator). Three-surface enclosure: radiator
    (gray, net heat Q given), array face (gray, at T_array) and space (black, at T_sink, closing the
    enclosure through the open sides, with a large area so that reciprocity makes F_space,i ~ 0).
    Returns the radiator temperature, the radiator-to-array view factor and the heat exchanged."""
    A_r, A_a = a * b, a * w_array
    F_ra = vf_perpendicular_rectangles(a, b, w_array)
    F_ar = A_r * F_ra / A_a
    A_s = 1e8
    F = np.array([[0.0, F_ra, 1 - F_ra],
                  [F_ar, 0.0, 1 - F_ar],
                  [A_r * (1 - F_ra) / A_s, A_a * (1 - F_ar) / A_s, 0.0]])
    F[2, 2] = 1 - F[2, 0] - F[2, 1]
    sol = solve_enclosure([A_r, A_a, A_s], F, [eps_r, eps_a, 1.0], T=[None, T_array, T_sink], Q=[Q, None, None])
    return {"T_rad": float(sol["T"][0]), "F_ra": float(F_ra), "Q_array": float(sol["Q"][1]),
            "Q_space": float(sol["Q"][2])}


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save radiation.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    # black body
    for T in (300.0, T_SUN):
        tot = quad(lambda lam: float(planck_spectral(lam, T)), 1e-8, 1e-3, limit=500, points=[2.9e-3 / T])[0]
        print(f"black body {T:.0f} K: sigma T^4 = {SIGMA * T ** 4:.4e} W/m^2, integral of Planck (10 nm - 1 mm) = {tot:.4e}")
    lam = np.logspace(-7, -3.5, 4000)
    print(f"Wien peak at 300 K: {lam[np.argmax(planck_spectral(lam, 300.0))] * 1e6:.2f} um (lambda T = "
          f"{lam[np.argmax(planck_spectral(lam, 300.0))] * 300 * 1e6:.0f} um K)")
    print(f"fraction below 3 um: sun (5772 K) {blackbody_fraction(1e-9, 3e-6, T_SUN):.3f}, 300 K body "
          f"{blackbody_fraction(1e-9, 3e-6, 300.0):.2e}")

    # view factors: closed forms against numerical integration
    a, b, c = 1.0, 0.8, 0.5
    Fp = vf_parallel_rectangles(a, b, c)
    print(f"\nparallel rectangles {a} x {b} m, gap {c} m: closed form F = {Fp:.6f}")
    for n in (5, 10, 20, 40):
        Fn = vf_numerical([0, 0, 0], [a, 0, 0], [0, b, 0], [0, 0, c], [a, 0, 0], [0, b, 0], n)
        print(f"  numerical n = {n:2d}: {Fn:.6f}, error {Fn - Fp:+.2e}")
    l, w, h = 1.0, 0.8, 1.2
    Fq = vf_perpendicular_rectangles(l, w, h)
    print(f"perpendicular rectangles, shared edge {l} m, widths {w} and {h} m: closed form F12 = {Fq:.6f}, "
          f"F21 = {vf_perpendicular_rectangles(l, h, w):.6f}; A1 F12 = {l * w * Fq:.6f}, A2 F21 = "
          f"{l * h * vf_perpendicular_rectangles(l, h, w):.6f}")
    for n in (5, 10, 20, 40):
        Fn = vf_numerical([0, 0, 0], [0, l, 0], [w, 0, 0], [0, 0, 0], [0, l, 0], [0, 0, h], n)
        print(f"  numerical n = {n:2d}: {Fn:.6f}, error {Fn - Fq:+.2e}")
    A, F = box_view_factors(1.0, 0.8, 0.5)
    print(f"box 1.0 x 0.8 x 0.5 m: row sums {np.array2string(F.sum(axis=1), precision=12)}; "
          f"max reciprocity residual {np.max(np.abs(A[:, None] * F - (A[:, None] * F).T)):.1e} m^2")
    print(f"crossed strings: opposed strips w = 1, h = 0.5: {vf_crossed_strings((0, 0), (1, 0), (0, 0.5), (1, 0.5)):.6f} "
          f"vs {vf_parallel_strips(1.0, 0.5):.6f}; perpendicular strips: "
          f"{vf_crossed_strings((1, 0), (0, 0), (0, 0.7), (0, 0)):.6f} vs {vf_perpendicular_strips(1.0, 0.7):.6f}")
    print(f"coaxial disks r1 = 0.1, r2 = 0.2, h = 0.3 m: F12 = {vf_coaxial_disks(0.1, 0.2, 0.3):.5f}, "
          f"A1 F12 = {np.pi * 0.01 * vf_coaxial_disks(0.1, 0.2, 0.3):.6f}, A2 F21 = {np.pi * 0.04 * vf_coaxial_disks(0.2, 0.1, 0.3):.6f}")

    # enclosure: a hot plate, a reradiating wall and a cold plate in a box
    Fpar = vf_parallel_rectangles(1.0, 1.0, 1.0)
    F3 = np.array([[0, Fpar, 1 - Fpar], [Fpar, 0, 1 - Fpar], [(1 - Fpar) / 4, (1 - Fpar) / 4, 0]])
    F3[2, 2] = 1 - F3[2, 0] - F3[2, 1]
    sol = solve_enclosure([1, 1, 4], F3, [0.8, 0.6, 0.5], T=[600.0, 300.0, None], Q=[None, None, 0.0])
    print(f"\nenclosure: hot plate 600 K, cold plate 300 K (1 m^2, 1 m apart), reradiating walls: Q = "
          f"{np.array2string(sol['Q'], precision=2)} W, sum {sol['Q'].sum():.1e} W, wall T = {sol['T'][2]:.1f} K")
    q2 = two_surface_exchange(1.0, 1.0, 1.0, 0.8, 0.6, 600.0, 300.0)
    print(f"infinite parallel plates 0.8/0.6, 600/300 K: {q2:.1f} W/m^2 (= sigma dT^4 / (1/e1 + 1/e2 - 1))")

    # spacecraft radiator
    Q, Tr = 200.0, 293.15
    eps = ZNO_WHITE_PAINT["eps_ir"]
    a_mid = float(np.mean(ZNO_WHITE_PAINT["alpha_s"]))
    A_ecl = radiator_area(Q, Tr, eps)
    A_sun = radiator_area(Q, Tr, eps, alpha_s=a_mid, cos_sun=1.0)
    A_sun45 = radiator_area(Q, Tr, eps, alpha_s=a_mid, cos_sun=np.cos(np.radians(60.0)))
    print(f"\nradiator {Q:.0f} W at {Tr - 273.15:.0f} C, ZnO white paint eps {eps}, alpha_s {a_mid}:")
    print(f"  area, no sun {A_ecl:.3f} m^2; sun at normal incidence {A_sun:.3f} m^2; sun 60 deg from normal {A_sun45:.3f} m^2")
    for name, al, ep in (("ZnO white, best", 0.12, 0.93), ("ZnO white, worst", 0.18, 0.93),
                         ("acrylic white", 0.26, 0.90), ("carbon black", 0.94, 0.83)):
        try:
            area = f"{radiator_area(Q, Tr, ep, alpha_s=al, cos_sun=1.0):.3f} m^2"
        except ValueError:
            area = "impossible (absorbs more than it emits at 20 C)"
        print(f"  {name:16s} alpha_s {al:.2f}, eps {ep:.2f}: area in normal sun {area}")
    A_design = 0.60
    print(f"  a {A_design} m^2 panel: eclipse {radiator_temperature(Q, A_design, eps) - 273.15:.1f} C, normal sun "
          f"{radiator_temperature(Q, A_design, eps, a_mid, cos_sun=1.0) - 273.15:.1f} C; with the electronics off, "
          f"holding it at -20 C takes a {eps * SIGMA * A_design * (253.15 ** 4 - T_SPACE ** 4):.0f} W heater")
    for w_arr in (0.0, 0.5, 1.0, 2.0):
        if w_arr == 0.0:
            print(f"  no array: radiator {radiator_temperature(Q, A_design, eps) - 273.15:.1f} C")
            continue
        r = radiator_with_array(Q, 1.0, A_design, w_arr, eps, 0.83, 330.0)
        print(f"  330 K array face {w_arr:.1f} m wide along the 1.0 m edge: F_ra {r['F_ra']:.3f}, radiator "
              f"{r['T_rad'] - 273.15:.1f} C; array face emits net {r['Q_array']:.0f} W, space receives {-r['Q_space']:.0f} W")
    A_need = brentq(lambda x: radiator_with_array(Q, 1.0, x, 1.0, eps, 0.83, 330.0)["T_rad"] - Tr, 0.3, 3.0)
    print(f"  area for 20 C beside a 1.0 m wide 330 K array: {A_need:.3f} m^2 (1.0 m edge, {A_need:.3f} m deep)")

    # figure
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    lam = np.logspace(-7, -4, 600)
    for T, ls in ((T_SUN, "-"), (300.0, "--")):
        e = planck_spectral(lam, T)
        ax[0].loglog(lam * 1e6, e / e.max(), ls, label=f"{T:.0f} K (normalised)")
    ax[0].set_xlabel("wavelength (um)"); ax[0].set_ylabel("e_b,lambda / max"); ax[0].set_ylim(1e-4, 2)
    ax[0].set_title("Planck: sunlight and room-temperature emission"); ax[0].legend(fontsize=8)
    ns = np.array([4, 8, 16, 32])
    ep = [abs(vf_numerical([0, 0, 0], [a, 0, 0], [0, b, 0], [0, 0, c], [a, 0, 0], [0, b, 0], int(n)) - Fp) for n in ns]
    eq = [abs(vf_numerical([0, 0, 0], [0, l, 0], [w, 0, 0], [0, 0, 0], [0, l, 0], [0, 0, h], int(n)) - Fq) for n in ns]
    ax[1].loglog(ns, ep, "o-", label="parallel rectangles")
    ax[1].loglog(ns, eq, "s-", label="perpendicular, shared edge")
    ax[1].loglog(ns, ep[0] * (ns[0] / ns) ** 2, "k--", lw=0.8, label="slope -2")
    ax[1].loglog(ns, eq[0] * (ns[0] / ns) ** 1, "k:", lw=0.8, label="slope -1")
    ax[1].set_xlabel("cells per side n"); ax[1].set_ylabel("|numerical - closed form|")
    ax[1].set_title("view factor by double-area integration"); ax[1].legend(fontsize=8)
    As = np.linspace(0.3, 2.0, 100)
    ax[2].plot(As, [radiator_temperature(Q, x, eps) - 273.15 for x in As], label="no sun")
    ax[2].plot(As, [radiator_temperature(Q, x, eps, a_mid, cos_sun=1.0) - 273.15 for x in As], label="normal sun, alpha 0.15")
    ax[2].plot(As, [radiator_temperature(Q, x, eps, 0.26, cos_sun=1.0) - 273.15 for x in As], label="normal sun, alpha 0.26")
    ax[2].axhline(20, color="k", lw=0.8, ls=":")
    ax[2].set_xlabel("radiator area (m^2)"); ax[2].set_ylabel("radiator temperature (C)")
    ax[2].set_title(f"radiator rejecting {Q:.0f} W"); ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "radiation.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'radiation.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
