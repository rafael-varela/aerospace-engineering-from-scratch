"""P10.09 Finite Wings & Induced Drag: reference solution.

How much of the reference trainer's cruise drag is the price of making lift with a finite wing, and
at what speed does it fly most efficiently? This module builds the answer from vortex filaments:

* the Biot-Savart law for a straight vortex segment and a semi-infinite one, vectorised over many
  points and segments;
* the horseshoe vortex (a bound segment and two trailing legs), checked against the infinite-line
  result Gamma/(2 pi h), Stokes' theorem (the circulation around a leg is its strength) and
  Kelvin/Helmholtz (no net circulation around the whole wake, a divergence-free field);
* a lifting line built from a stack of horseshoes: an elliptic span load gives a uniform downwash
  w = Gamma_0/(2b) and CDi = CL^2/(pi AR); induced drag computed two independent ways, near field
  (Kutta-Joukowski with the downwash at the bound vortex, 3-D Biot-Savart) and far field (Trefftz
  plane, 2-D vortices);
* the parabolic drag polar CD = CD0 + CL^2/(pi e AR), best L/D, the minimum-drag and minimum-power
  speeds of the trainer.

Axes: x downstream (along the free stream), y to starboard (right), z up. A bound vortex with its
vector along +y carries positive circulation and gives lift (F = rho V x Gamma). Downwash w is the
-z velocity, reported positive downward. SI units. ISA sea level.

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes induced_drag.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-10/lesson-09-finite-wings-induced-drag
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable

import numpy as np

RHO = 1.225                         # kg/m^3, ISA sea level
G = 9.81                            # m/s^2
TRAINER = dict(mass=750.0, S=12.0, b=10.0, V_cruise=55.0, V_approach=30.0)
CD0_ASSUMED = 0.030                 # ASSUMED parasite drag coefficient of the trainer (fixed gear; P10.13 builds it up)
E_ASSUMED = 0.80                    # ASSUMED Oswald efficiency of the whole aircraft (wing alone: P10.10-P10.11)


# =============================================================== given helpers

def cosine_layout(b: float, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Given: n strips across the span -b/2..b/2 laid out on the "semicircle": edges at
    y = -(b/2) cos(phi_k), phi_k = k pi/n, and control points at the phi-midpoints (not the
    y-midpoints). Fine at the tips; with this layout a sampled elliptic load has an exactly
    uniform Trefftz-plane downwash ratio and e = 1 to round-off. Returns (edges, y_control)."""
    phi = np.linspace(0.0, np.pi, n + 1)
    return -0.5 * b * np.cos(phi), -0.5 * b * np.cos(0.5 * (phi[:-1] + phi[1:]))


def elliptic_gamma(y: np.ndarray, b: float, gamma0: float) -> np.ndarray:
    """Given: the elliptic span load Gamma(y) = Gamma_0 sqrt(1 - (2y/b)^2)."""
    return gamma0 * np.sqrt(np.clip(1 - (2 * np.asarray(y) / b) ** 2, 0.0, None))


# =============================================================== Biot-Savart (learner)

def segment_velocity(P: np.ndarray, A: np.ndarray, B: np.ndarray, gamma: float | np.ndarray = 1.0,
                     core: float = 1e-10) -> np.ndarray:
    """Velocity induced at points P (n, 3) by straight vortex segments A -> B (m, 3) of strength
    gamma (scalar or (m,)), Biot-Savart for a finite segment:
        v = gamma/(4 pi) (r1 x r2)/|r1 x r2|^2 * r0 . (r1/|r1| - r2/|r2|),
    r0 = B - A, r1 = P - A, r2 = P - B. Points within `core` (relative to |r0|) of the segment's line
    get zero (the filament itself). Returns (n, m, 3): the contribution of every segment."""
    P = np.atleast_2d(P).astype(float)
    A = np.atleast_2d(A).astype(float); B = np.atleast_2d(B).astype(float)
    r0 = (B - A)[None, :, :]
    r1 = P[:, None, :] - A[None, :, :]
    r2 = P[:, None, :] - B[None, :, :]
    c = np.cross(r1, r2)
    c2 = np.sum(c * c, axis=-1)
    n1 = np.linalg.norm(r1, axis=-1); n2 = np.linalg.norm(r2, axis=-1)
    small = c2 <= (core * np.sum(r0 * r0, axis=-1)) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        k = np.sum(r0 * (r1 / n1[..., None] - r2 / n2[..., None]), axis=-1) / c2
    k = np.where(small, 0.0, k)
    g = np.broadcast_to(np.asarray(gamma, dtype=float), (A.shape[0],))[None, :]
    return (g / (4 * np.pi) * k)[..., None] * c


def semi_infinite_velocity(P: np.ndarray, A: np.ndarray, d: np.ndarray, gamma: float | np.ndarray = 1.0,
                           core: float = 1e-10) -> np.ndarray:
    """Velocity at points P (n, 3) from semi-infinite vortex lines that start at A (m, 3) and run to
    infinity along the unit vector d (3,) (the vortex vector points along d), strength gamma:
        v = gamma/(4 pi) (d x r1)/|d x r1|^2 (1 + d . r1/|r1|),   r1 = P - A
    (the finite-segment formula with B -> A + L d, L -> infinity). Returns (n, m, 3)."""
    P = np.atleast_2d(P).astype(float)
    A = np.atleast_2d(A).astype(float)
    d = np.asarray(d, dtype=float) / np.linalg.norm(d)
    r1 = P[:, None, :] - A[None, :, :]
    c = np.cross(np.broadcast_to(d, r1.shape), r1)
    c2 = np.sum(c * c, axis=-1)
    n1 = np.linalg.norm(r1, axis=-1)
    small = c2 <= core ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        k = (1 + (r1 @ d) / n1) / c2
    k = np.where(small, 0.0, k)
    g = np.broadcast_to(np.asarray(gamma, dtype=float), (A.shape[0],))[None, :]
    return (g / (4 * np.pi) * k)[..., None] * c


def horseshoe_velocity(P: np.ndarray, y_left: np.ndarray, y_right: np.ndarray, gamma: np.ndarray,
                       x_bound: float | np.ndarray = 0.0, z: float = 0.0) -> np.ndarray:
    """Total velocity (n, 3) at points P (n, 3) from horseshoe vortices: bound segments from
    (x_bound, y_left, z) to (x_bound, y_right, z) (vector along +y: positive gamma lifts) and
    trailing legs to x = +infinity: the right leg starts at the bound vortex's right end with vector
    +x, the left leg (vector -x, coming in from downstream) is a +x line of strength -gamma."""
    y_left, y_right = np.atleast_1d(y_left).astype(float), np.atleast_1d(y_right).astype(float)
    g = np.broadcast_to(np.asarray(gamma, dtype=float), y_left.shape)
    xb = np.broadcast_to(np.asarray(x_bound, dtype=float), y_left.shape)
    zz = np.full_like(y_left, z)
    A = np.column_stack([xb, y_left, zz]); B = np.column_stack([xb, y_right, zz])
    ex = np.array([1.0, 0.0, 0.0])
    v = segment_velocity(P, A, B, g).sum(axis=1)
    v += semi_infinite_velocity(P, B, ex, g).sum(axis=1)
    v += semi_infinite_velocity(P, A, ex, -g).sum(axis=1)
    return v


def loop_circulation(velocity: Callable[[np.ndarray], np.ndarray], center: np.ndarray, normal: np.ndarray,
                     radius: float, n: int = 400) -> float:
    """Circulation (line integral of v . dl) around a circle of the given centre, unit normal and
    radius, traversed counter-clockwise about the normal (right-hand rule), by the midpoint rule
    with n points (spectrally accurate for a smooth periodic integrand)."""
    nrm = np.asarray(normal, dtype=float) / np.linalg.norm(normal)
    a = np.cross(nrm, [1.0, 0.0, 0.0])
    if np.linalg.norm(a) < 1e-8:
        a = np.cross(nrm, [0.0, 1.0, 0.0])
    a /= np.linalg.norm(a)
    b = np.cross(nrm, a)
    t = (np.arange(n) + 0.5) * 2 * np.pi / n
    pts = np.asarray(center)[None, :] + radius * (np.cos(t)[:, None] * a + np.sin(t)[:, None] * b)
    dl = radius * (-np.sin(t)[:, None] * a + np.cos(t)[:, None] * b) * (2 * np.pi / n)
    return float(np.sum(velocity(pts) * dl))


# =============================================================== lifting line from horseshoes (learner)

def lifting_line_downwash(edges: np.ndarray, gamma_strips: np.ndarray, y_control: np.ndarray) -> np.ndarray:
    """Downwash w (positive down) at the control points (x = 0, y_control, z = 0) on the bound
    vortex line of a stack of horseshoes, one per strip between consecutive edges, by 3-D
    Biot-Savart. On the line the bound segments contribute nothing (collinear), so this is the
    effect of the semi-infinite trailing legs: half the Trefftz-plane value."""
    yc = np.asarray(y_control, dtype=float)
    P = np.column_stack([np.zeros_like(yc), yc, np.zeros_like(yc)])
    return -horseshoe_velocity(P, edges[:-1], edges[1:], gamma_strips)[:, 2]


def trefftz_downwash(edges: np.ndarray, gamma_strips: np.ndarray, y_control: np.ndarray) -> np.ndarray:
    """Far-field (Trefftz-plane) downwash at the control points. Far downstream the trailing legs
    are infinite 2-D vortices; the one at edge y_k has its vector along +x and strength
    s_k = Gamma_left(k) - Gamma_right(k) (zero circulation outside the span). Each induces
    v_z = s_k / (2 pi (y - y_k)), so w_T(y) = -sum_k s_k / (2 pi (y - y_k)) (positive down)."""
    G = np.concatenate([[0.0], gamma_strips, [0.0]])
    s = G[:-1] - G[1:]
    yc = np.asarray(y_control, dtype=float)
    return -np.sum(s[None, :] / (2 * np.pi * (yc[:, None] - edges[None, :])), axis=1)


def lift_and_drag(edges: np.ndarray, gamma_strips: np.ndarray, y_control: np.ndarray, V: float,
                  rho: float = RHO) -> dict:
    """Lift and induced drag of a discretised span load, Kutta-Joukowski strip by strip:
    L = rho V sum Gamma_i dy_i; near field D_i = rho sum Gamma_i w_i dy_i with w from
    lifting_line_downwash (3-D Biot-Savart); far field D_i = (rho/2) sum Gamma_i w_T,i dy_i with
    w_T from trefftz_downwash (2-D vortices). Returns dict(L, Di_near, Di_trefftz, w, wT)."""
    dy = np.diff(edges)
    w = lifting_line_downwash(edges, gamma_strips, y_control)
    wT = trefftz_downwash(edges, gamma_strips, y_control)
    return dict(L=float(rho * V * np.sum(gamma_strips * dy)),
                Di_near=float(rho * np.sum(gamma_strips * w * dy)),
                Di_trefftz=float(0.5 * rho * np.sum(gamma_strips * wT * dy)), w=w, wT=wT)


def span_efficiency(L: float, Di: float, b: float, V: float, rho: float = RHO) -> float:
    """e = L^2 / (pi q b^2 Di) = CL^2/(pi AR CDi) with q = rho V^2/2 (independent of S)."""
    return float(L ** 2 / (np.pi * 0.5 * rho * V ** 2 * b ** 2 * Di))


def elliptic_gamma0(CL: float, V: float, S: float, b: float) -> float:
    """Root circulation of the elliptic load giving CL: L = rho V Gamma_0 pi b/4 = q S CL, so
    Gamma_0 = 2 CL V S/(pi b)."""
    return 2 * CL * V * S / (np.pi * b)


# =============================================================== drag polar (learner)

def drag_polar(CL: np.ndarray, CD0: float, e: float, AR: float) -> np.ndarray:
    """Parabolic polar CD = CD0 + CL^2/(pi e AR)."""
    return CD0 + np.asarray(CL) ** 2 / (np.pi * e * AR)


def best_glide(CD0: float, e: float, AR: float) -> dict:
    """Maximum L/D of the parabolic polar: induced drag equals parasite drag at
    CL* = sqrt(pi e AR CD0), (L/D)max = (1/2) sqrt(pi e AR / CD0). Minimum power (max CL^1.5/CD)
    is at CL_mp = sqrt(3 pi e AR CD0) = sqrt(3) CL*, with CDi = 3 CD0 there."""
    k = 1 / (np.pi * e * AR)
    CLs = np.sqrt(CD0 / k)
    return dict(CL_md=float(CLs), LD_max=float(0.5 / np.sqrt(k * CD0)), CL_mp=float(np.sqrt(3 * CD0 / k)),
                LD_mp=float(np.sqrt(3 * CD0 / k) / (4 * CD0)))


def speed_for_cl(mass: float, S: float, CL: float, rho: float = RHO) -> float:
    """Level-flight speed at which the aircraft flies at CL: V = sqrt(2 W/(rho S CL))."""
    return float(np.sqrt(2 * mass * G / (rho * S * CL)))


def drag_and_power(V: np.ndarray, mass: float, S: float, CD0: float, e: float, AR: float,
                   rho: float = RHO) -> dict:
    """Level flight at speed V: CL = W/(q S), D = q S CD, P = D V; also the parasite and induced
    parts. Returns dict(CL, D, D0, Di, P) as arrays."""
    V = np.asarray(V, dtype=float)
    q = 0.5 * rho * V ** 2
    CL = mass * G / (q * S)
    D0 = q * S * CD0
    Di = q * S * CL ** 2 / (np.pi * e * AR)
    return dict(CL=CL, D=D0 + Di, D0=D0, Di=Di, P=(D0 + Di) * V)


# =============================================================== main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save induced_drag.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    m_, S, b, V = TRAINER["mass"], TRAINER["S"], TRAINER["b"], TRAINER["V_cruise"]
    AR = b ** 2 / S
    W = m_ * G

    # 1. Biot-Savart checks
    h = 0.5
    print("finite segment +-L at distance h = 0.5 m against the infinite line Gamma/(2 pi h):")
    for Lh in (10, 100, 1000):
        v = segment_velocity(np.array([[0.0, 0.0, h]]), np.array([[0.0, -Lh * h, 0.0]]), np.array([[0.0, Lh * h, 0.0]]))
        ratio = np.linalg.norm(v) / (1 / (2 * np.pi * h))
        print(f"  L/h = {Lh:5d}: v / (Gamma/(2 pi h)) = {ratio:.8f}  (1 - ratio = {1 - ratio:.3e}; "
              f"exact L/sqrt(L^2 + h^2) gives {1 - Lh / np.sqrt(Lh ** 2 + 1):.3e})")
    hs = lambda P: horseshoe_velocity(P, np.array([-1.0]), np.array([1.0]), np.array([1.0]))
    print(f"circulation round the right trailing leg {loop_circulation(hs, np.array([3.0, 1.0, 0.0]), np.array([1.0, 0, 0]), 0.2):.6f}, "
          f"left leg {loop_circulation(hs, np.array([3.0, -1.0, 0.0]), np.array([1.0, 0, 0]), 0.2):.6f}, "
          f"bound leg {loop_circulation(hs, np.array([0.0, 0.0, 0.0]), np.array([0, 1.0, 0]), 0.2):.6f}, "
          f"whole wake {loop_circulation(hs, np.array([5.0, 0.0, 0.0]), np.array([1.0, 0, 0]), 3.0):.2e}  (Gamma = 1)")

    # 2. elliptic wing from horseshoes
    CL = W / (0.5 * RHO * V ** 2 * S)
    g0 = elliptic_gamma0(CL, V, S, b)
    print(f"\nreference trainer cruise: CL = {CL:.4f}, AR = {AR:.3f}, elliptic Gamma_0 = {g0:.3f} m^2/s, "
          f"w = Gamma_0/(2b) = {g0 / (2 * b):.4f} m/s, alpha_i = CL/(pi AR) = {np.degrees(CL / (np.pi * AR)):.3f} deg")
    print("  n strips   max|w/w_exact - 1|   CDi near     CDi Trefftz   CL^2/(pi AR)   e (Trefftz)")
    q = 0.5 * RHO * V ** 2
    errs, ns = [], [10, 20, 40, 80, 160]
    for n in ns:
        ed, yc = cosine_layout(b, n)
        gs = elliptic_gamma(yc, b, g0)
        r = lift_and_drag(ed, gs, yc, V)
        err = np.max(np.abs(r["w"] / (g0 / (2 * b)) - 1))
        errs.append(abs(r["Di_trefftz"] / (q * S) - CL ** 2 / (np.pi * AR)))
        print(f"  {n:5d}      {err:.3e}         {r['Di_near'] / (q * S):.7f}    {r['Di_trefftz'] / (q * S):.7f}     "
              f"{CL ** 2 / (np.pi * AR):.7f}      {span_efficiency(r['L'], r['Di_trefftz'], b, V):.12f}")
    orders = np.log(np.array(errs[:-1]) / np.array(errs[1:])) / np.log(2)
    print(f"  observed order of the CDi error: {np.round(orders, 3)}")
    ed40 = cosine_layout(b, 40)[0]
    ym = 0.5 * (ed40[:-1] + ed40[1:])
    rm = lift_and_drag(ed40, elliptic_gamma(ym, b, g0), ym, V)
    print(f"  same 40 strips with control points at the y-midpoints instead: e = {span_efficiency(rm['L'], rm['Di_trefftz'], b, V):.4f}, "
          f"max|w/w_exact - 1| = {np.max(np.abs(rm['w'] / (g0 / (2 * b)) - 1)):.3e}")

    # a single horseshoe with the same lift
    gh = W / (RHO * V * b)
    w0 = -horseshoe_velocity(np.array([[0.0, 0.0, 0.0]]), np.array([-b / 2]), np.array([b / 2]), np.array([gh]))[0, 2]
    print(f"single horseshoe, same lift: Gamma = {gh:.3f}, centre downwash {w0:.4f} m/s (Gamma/(pi b) = {gh / (np.pi * b):.4f}) "
          f"vs elliptic {g0 / (2 * b):.4f}; its Trefftz drag is infinite (point vortices at the tips)")

    # 3. drag polar and speeds
    bg = best_glide(CD0_ASSUMED, E_ASSUMED, AR)
    Vmd = speed_for_cl(m_, S, bg["CL_md"]); Vmp = speed_for_cl(m_, S, bg["CL_mp"])
    print(f"\ndrag polar with ASSUMED CD0 = {CD0_ASSUMED}, e = {E_ASSUMED}: (L/D)max = {bg['LD_max']:.2f} at CL = {bg['CL_md']:.4f}, "
          f"V_md = {Vmd:.2f} m/s; minimum power at CL = {bg['CL_mp']:.4f}, V_mp = {Vmp:.2f} m/s (= V_md/3^0.25 = {Vmd / 3 ** 0.25:.2f}), L/D = {bg['LD_mp']:.2f}")
    for Vx in (V, TRAINER["V_approach"], Vmd, Vmp):
        d = drag_and_power(Vx, m_, S, CD0_ASSUMED, E_ASSUMED, AR)
        print(f"  V = {Vx:6.2f} m/s: CL {float(d['CL']):.4f}, D = {float(d['D']):7.1f} N (parasite {float(d['D0']):6.1f}, "
              f"induced {float(d['Di']):6.1f}, {float(d['Di'] / d['D']) * 100:4.1f} %), P = {float(d['P']) / 1e3:6.2f} kW, L/D {W / float(d['D']):.2f}")
    print("  sensitivity of (L/D)max and V_md to the assumptions:")
    for cd0, e in ((0.025, 0.80), (0.035, 0.80), (0.030, 0.70), (0.030, 0.90)):
        bb = best_glide(cd0, e, AR)
        print(f"    CD0 {cd0:.3f}, e {e:.2f}: (L/D)max {bb['LD_max']:.2f}, V_md {speed_for_cl(m_, S, bb['CL_md']):.2f} m/s")

    # figure
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    ed, yc = cosine_layout(b, 40)
    gs = elliptic_gamma(yc, b, g0)
    r = lift_and_drag(ed, gs, yc, V)
    ax[0].step(ed, np.concatenate([gs, [gs[-1]]]), where="post", label="40 horseshoes")
    yy = np.linspace(-b / 2, b / 2, 400)
    ax[0].plot(yy, elliptic_gamma(yy, b, g0), "k--", lw=0.8, label="elliptic")
    ax[0].set_xlabel("y (m)"); ax[0].set_ylabel("Gamma (m^2/s)"); ax[0].set_title("span load"); ax[0].legend()
    ax[1].plot(yc, r["w"], "o-", ms=3, label="at the lifting line (Biot-Savart)")
    ax[1].plot(yc, r["wT"], "s-", ms=3, label="Trefftz plane")
    ax[1].axhline(g0 / (2 * b), color="k", lw=0.8, ls="--", label="Gamma_0/(2b)")
    ax[1].set_xlabel("y (m)"); ax[1].set_ylabel("downwash w (m/s)"); ax[1].set_title("elliptic wing downwash"); ax[1].legend(fontsize=8)
    Vs = np.linspace(22, 70, 200)
    d = drag_and_power(Vs, m_, S, CD0_ASSUMED, E_ASSUMED, AR)
    ax[2].plot(Vs, d["D"], label="total"); ax[2].plot(Vs, d["D0"], label="parasite"); ax[2].plot(Vs, d["Di"], label="induced")
    ax[2].axvline(Vmd, color="k", lw=0.8, ls=":"); ax[2].axvline(Vmp, color="k", lw=0.8, ls="--")
    ax[2].set_xlabel("V (m/s)"); ax[2].set_ylabel("drag (N)"); ax[2].set_ylim(0, 2500)
    ax[2].set_title("trainer drag (assumed CD0, e)"); ax[2].legend()
    fig.tight_layout()
    fig.savefig(out_dir / "induced_drag.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
