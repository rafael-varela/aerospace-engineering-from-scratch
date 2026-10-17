"""P05.16 Energy Methods: reference solution.

An elastic structure stores the work done on it as strain energy. Castigliano's second theorem says
the deflection at a load, in the load's direction, is the derivative of that energy with respect to the
load; the unit-load (virtual-force) method is the same statement written as an integral of the real
internal forces times the internal forces of a unit load. For a statically indeterminate structure the
redundant force is the one that makes the energy stationary (dU/dX = 0), which is compatibility.

This module builds:

* strain energy in axial, bending and torsion members, and Castigliano by differentiating it;
* the unit-load method as a flexibility matrix F_jk = sum over members of
  integral (n_j n_k / EA + m_j m_k / EI + t_j t_k / GJ) dx, whose symmetry is Maxwell's reciprocal
  theorem and whose quadratic form is the energy;
* a cantilever spar under an elliptic lift distribution, checked against direct double integration of
  M / EI (the independent check);
* an out-of-plane L-bracket (bending plus torsion);
* trusses: deflection by unit loads, redundant bar forces by the force method (minimum complementary
  energy), both checked against the direct stiffness method;
* a strut-braced wing spar: the strut force as the redundant, with spar bending, spar axial and strut
  axial energy, checked by compatibility and by the propped-cantilever limit 3 q L / 8.

Units: SI throughout (N, m, Pa). Deflection w is positive upward; a bending moment is positive when it
bends the member concave up (E I w'' = M).

Run:  python main.py [out_dir]     (prints the tables, writes energy_methods.png)
Test: pytest curriculum/phase-05/lesson-16-energy-methods
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import cumulative_trapezoid, simpson


# ---------------------------------------------------------------- strain energy

def strain_energy_axial(N: float, L: float, EA: float) -> float:
    """U = N^2 L / (2 EA) for a bar of length L carrying a constant axial force N."""
    return N ** 2 * L / (2.0 * EA)


def strain_energy_torsion(T: float, L: float, GJ: float) -> float:
    """U = T^2 L / (2 GJ) for a shaft of length L carrying a constant torque T."""
    return T ** 2 * L / (2.0 * GJ)


def strain_energy_bending(x: np.ndarray, M: np.ndarray, EI: float | np.ndarray) -> float:
    """U = integral M^2 / (2 EI) dx by Simpson's rule on the samples (x, M)."""
    return float(simpson(np.asarray(M) ** 2 / (2.0 * np.asarray(EI)), x=x))


def castigliano_deflection(energy: Callable[[np.ndarray], float], loads: np.ndarray, j: int,
                           rel_step: float = 1e-3) -> float:
    """Castigliano's second theorem, numerically: delta_j = dU/dP_j by a central difference. For a
    linear structure U is quadratic in the loads, so the central difference is exact up to round-off."""
    P = np.asarray(loads, dtype=float)
    h = rel_step * max(np.max(np.abs(P)), 1.0)
    e = np.zeros_like(P)
    e[j] = h
    return (energy(P + e) - energy(P - e)) / (2.0 * h)


def unit_load_deflection(x: np.ndarray, M: np.ndarray, m: np.ndarray, EI: float | np.ndarray) -> float:
    """Unit-load method for bending: delta = integral M m / EI dx (Simpson), with M the real moment and
    m the moment caused by a unit load placed where, and in the direction in which, delta is wanted."""
    return float(simpson(np.asarray(M) * np.asarray(m) / np.asarray(EI), x=x))


# ---------------------------------------------------------------- members and the flexibility matrix

@dataclass
class Member:
    """A straight member of length L with stiffnesses EA, EI, GJ (np.inf for a term to ignore).

    actions(s) returns an array of shape (n_loads, 3, len(s)): the axial force, bending moment and
    torque at the stations s (0 <= s <= L) caused by a UNIT value of each load in turn."""
    L: float
    EA: float
    EI: float
    GJ: float
    actions: Callable[[np.ndarray], np.ndarray]


def flexibility_matrix(members: list[Member], n: int = 201) -> np.ndarray:
    """F_jk = sum over members of integral (n_j n_k / EA + m_j m_k / EI + t_j t_k / GJ) ds, by
    Simpson's rule with n stations per member. delta = F P; U = (1/2) P^T F P."""
    F = None
    for mb in members:
        s = np.linspace(0.0, mb.L, n)
        a = mb.actions(s)                                     # (n_loads, 3, n)
        inv = np.array([1.0 / mb.EA, 1.0 / mb.EI, 1.0 / mb.GJ])  # 1/inf = 0 drops the term
        integrand = np.einsum("jcs,kcs,c->jks", a, a, inv)
        Fm = simpson(integrand, x=s, axis=-1)
        F = Fm if F is None else F + Fm
    return F


def strain_energy_members(members: list[Member], loads: np.ndarray, n: int = 201) -> float:
    """Strain energy by integrating the actual N^2/2EA + M^2/2EI + T^2/2GJ along every member (the actual
    actions are the unit-load actions scaled by the loads and summed)."""
    P = np.asarray(loads, dtype=float)
    U = 0.0
    for mb in members:
        s = np.linspace(0.0, mb.L, n)
        act = np.einsum("j,jcs->cs", P, mb.actions(s))
        U += float(simpson(act[0] ** 2 / (2 * mb.EA) + act[1] ** 2 / (2 * mb.EI) + act[2] ** 2 / (2 * mb.GJ), x=s))
    return U


def l_bracket_members(a: float, b: float, EI: float, GJ: float) -> list[Member]:
    """An L-shaped bracket built in at its root: a leg of length b from the root to the corner, then an
    arm of length a at right angles in the same horizontal plane, both round tubes. Load 0 is a vertical
    (out-of-plane) force at the arm's tip, load 1 a vertical force at the corner.

    Arm (s measured from the tip): M = P0 s, no torque. Leg (s from the corner): M = P0 s + P1 s and a
    torque T = P0 a, because the tip force acts at the arm's length a from the leg's axis."""
    def arm(s):
        z = np.zeros_like(s)
        return np.array([[z, s, z], [z, z, z]])

    def leg(s):
        z = np.zeros_like(s)
        return np.array([[z, s, a + z], [z, s, z]])

    return [Member(a, np.inf, EI, GJ, arm), Member(b, np.inf, EI, GJ, leg)]


def l_bracket_tip_deflection(P: float, a: float, b: float, EI: float, GJ: float) -> float:
    """Closed form for the tip deflection of the L-bracket: P a^3/(3EI) + P b^3/(3EI) + P a^2 b/(GJ)."""
    return P * a ** 3 / (3 * EI) + P * b ** 3 / (3 * EI) + P * a ** 2 * b / GJ


# ---------------------------------------------------------------- the cantilever spar

def elliptic_lift(x: np.ndarray, L: float, total: float) -> np.ndarray:
    """Elliptic spanwise lift per unit length on a semi-span L: q0 sqrt(1 - (x/L)^2), q0 = 4 total/(pi L)."""
    q0 = 4.0 * total / (np.pi * L)
    return q0 * np.sqrt(np.clip(1.0 - (np.asarray(x) / L) ** 2, 0.0, None))


def cantilever_moment(x: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Bending moment of a cantilever (root at x[0], free at x[-1]) under an upward load q(x) per unit
    length: M(x) = integral from x to L of q(xi)(xi - x) d xi = S1(x) - x S0(x), with S0, S1 the
    integrals of q and xi q from x to the tip (cumulative trapezoid). Positive M bends the spar up."""
    x = np.asarray(x, dtype=float)
    q = np.asarray(q, dtype=float)
    c0 = cumulative_trapezoid(q, x, initial=0.0)
    c1 = cumulative_trapezoid(q * x, x, initial=0.0)
    S0 = c0[-1] - c0
    S1 = c1[-1] - c1
    return S1 - x * S0


def deflection_by_integration(x: np.ndarray, M: np.ndarray, EI: float | np.ndarray) -> np.ndarray:
    """Direct integration of E I w'' = M from a built-in root at x[0]: w(0) = w'(0) = 0, twice the
    cumulative trapezoid rule. The independent route to the unit-load deflection."""
    slope = cumulative_trapezoid(np.asarray(M) / np.asarray(EI), x, initial=0.0)
    return cumulative_trapezoid(slope, x, initial=0.0)


# ---------------------------------------------------------------- trusses

def truss_equilibrium_matrix(nodes: np.ndarray, bars: list[tuple[int, int]], free: np.ndarray) -> np.ndarray:
    """Equilibrium matrix C (n_free_dofs x n_bars) of a 2-D pin-jointed truss: C N = F, with N the bar
    tensions (positive in tension) and F the external loads on the free degrees of freedom (2 per node:
    x then y). At a node, a bar in tension pulls toward its other end, and the pulls balance the load:
    sum N e_toward_other_end + F = 0, so column k holds minus that unit vector at each end. Its transpose
    maps nodal displacements to bar elongations: (C^T u)_k = e . (u_j - u_i)."""
    nodes = np.asarray(nodes, dtype=float)
    C = np.zeros((2 * len(nodes), len(bars)))
    for k, (i, j) in enumerate(bars):
        e = nodes[j] - nodes[i]
        e = e / np.linalg.norm(e)        # unit vector from node i to node j
        C[2 * i:2 * i + 2, k] = -e
        C[2 * j:2 * j + 2, k] = e
    return C[np.asarray(free)]


def truss_lengths(nodes: np.ndarray, bars: list[tuple[int, int]]) -> np.ndarray:
    """Bar lengths."""
    nodes = np.asarray(nodes, dtype=float)
    return np.array([np.linalg.norm(nodes[j] - nodes[i]) for i, j in bars])


def truss_force_method(nodes: np.ndarray, bars: list[tuple[int, int]], EA: np.ndarray, loads: np.ndarray,
                       free: np.ndarray) -> np.ndarray:
    """Bar forces of a statically determinate or indeterminate truss by the force method.

    Any equilibrium solution is N = N0 + S X, with N0 a particular solution of C N = F (least squares)
    and the columns of S spanning the self-stress states (null space of C). The redundants X make the
    complementary energy sum N^2 L / (2 EA) stationary, which is Castigliano's dU/dX = 0:
    X = -(S^T D S)^-1 S^T D N0 with D = diag(L / EA)."""
    C = truss_equilibrium_matrix(nodes, bars, free)
    F = np.asarray(loads, dtype=float)
    N0 = np.linalg.lstsq(C, F, rcond=None)[0]
    u, sv, vt = np.linalg.svd(C)
    rank = int(np.sum(sv > 1e-12 * sv[0]))
    if rank < C.shape[0]:
        raise ValueError("truss is a mechanism: the loads cannot be equilibrated")
    S = vt[rank:].T
    if S.shape[1] == 0:
        return N0
    D = np.diag(truss_lengths(nodes, bars) / np.asarray(EA, dtype=float))
    X = -np.linalg.solve(S.T @ D @ S, S.T @ D @ N0)
    return N0 + S @ X


def truss_unit_load_displacement(nodes: np.ndarray, bars: list[tuple[int, int]], EA: np.ndarray, N: np.ndarray,
                                 free: np.ndarray, dof: int) -> float:
    """Displacement of free dof number `dof` (an index into `free`) by the unit-load method:
    delta = sum N n L / EA, with n ANY bar-force set in equilibrium with a unit load at that dof (the
    principle of virtual forces needs only static admissibility, so a least-squares particular solution
    serves even for an indeterminate truss)."""
    C = truss_equilibrium_matrix(nodes, bars, free)
    e = np.zeros(C.shape[0])
    e[dof] = 1.0
    n = np.linalg.lstsq(C, e, rcond=None)[0]
    return float(np.sum(np.asarray(N) * n * truss_lengths(nodes, bars) / np.asarray(EA, dtype=float)))


def truss_stiffness_method(nodes: np.ndarray, bars: list[tuple[int, int]], EA: np.ndarray, loads: np.ndarray,
                           free: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Direct stiffness method, the independent check: K = C diag(EA / L) C^T on the free dofs,
    u = K^-1 F, bar elongations C^T u, tensions N = (EA / L) C^T u. Returns (u_free, N)."""
    C = truss_equilibrium_matrix(nodes, bars, free)
    k = np.asarray(EA, dtype=float) / truss_lengths(nodes, bars)
    K = C @ np.diag(k) @ C.T
    u = np.linalg.solve(K, np.asarray(loads, dtype=float))
    return u, k * (C.T @ u)


# ---------------------------------------------------------------- the strut-braced spar

def braced_spar(L: float, a: float, h: float, EI: float, EA_spar: float, EA_strut: float,
                q: Callable[[np.ndarray], np.ndarray], n: int = 2001) -> dict:
    """A wing spar built in at the fuselage (x = 0), semi-span L, carrying upward lift q(x), braced by a
    strut from a fuselage point a height h below the root to the spar at x = a.

    Redundant: the strut tension S. With sin(th) = h / L_s and cos(th) = a / L_s (L_s the strut length),
    the strut pulls the spar down with S sin(th) and inboard with S cos(th), so
        M(x) = M_q(x) - S sin(th) (a - x) for x < a,   N_spar = -S cos(th) on 0 < x < a,   N_strut = S.
    dU/dS = 0 with U = int M^2/2EI + N_spar^2 a / 2EA_spar + S^2 L_s / 2EA_strut gives
        S = sin(th) int_0^a M_q (a - x)/EI dx / (sin^2(th) int_0^a (a - x)^2/EI dx + cos^2(th) a/EA_spar + L_s/EA_strut).
    Returns a dict with S, the vertical strut reaction R = S sin(th), x, q, the final M, the deflection w by
    direct integration, the root shear and moment, the tip deflection, and the compatibility residual
    (w(a) sin(th) - strut elongation - spar shortening x cos(th)) / w(a)."""
    x = np.linspace(0.0, L, n)
    qx = q(x)
    Mq = cantilever_moment(x, qx)
    Ls = np.hypot(a, h)
    sn, cs = h / Ls, a / Ls
    m = np.where(x < a, a - x, 0.0)                         # moment of a unit downward pull at x = a, sign flipped
    num = sn * simpson(Mq * m / EI, x=x)
    den = sn ** 2 * simpson(m ** 2 / EI, x=x) + cs ** 2 * a / EA_spar + Ls / EA_strut
    S = num / den
    M = Mq - S * sn * m
    w = deflection_by_integration(x, M, EI)
    wa = float(np.interp(a, x, w))
    strut_stretch = S * Ls / EA_strut
    spar_shortening = S * cs * a / EA_spar
    resid = (wa * sn - strut_stretch - spar_shortening * cs) / max(abs(wa), 1e-30)
    lift = float(simpson(qx, x=x))
    return {"S": S, "R": S * sn, "x": x, "q": qx, "M": M, "w": w, "lift": lift,
            "root_shear": lift - S * sn, "root_moment": float(M[0]), "tip": float(w[-1]),
            "w_strut": wa, "compatibility": resid, "L_strut": Ls}


# ---------------------------------------------------------------- driver

def main(out_dir: str | Path = ".") -> None:
    """Print the four worked cases and the experiment tables; save energy_methods.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)

    # ---- 1. cantilever spar, elliptic lift (2024-T3, I = 4e-6 m^4)
    E_al = 73.1e9
    L, I, total = 5.0, 4.0e-6, 18.6e3
    EI = E_al * I
    x = np.linspace(0.0, L, 4001)
    q = elliptic_lift(x, L, total)
    M = cantilever_moment(x, q)
    d_ul = unit_load_deflection(x, M, L - x, EI)
    w = deflection_by_integration(x, M, EI)
    print(f"cantilever spar L = {L} m, EI = {EI/1e3:.1f} kN m^2, elliptic lift {total/1e3:.1f} kN: root moment "
          f"{M[0]/1e3:.2f} kN m (4 W L / 3 pi = {4*total*L/(3*np.pi)/1e3:.2f}); tip deflection unit load "
          f"{d_ul*1e3:.2f} mm, direct integration {w[-1]*1e3:.2f} mm")
    print(f"  same total as a tip load: {total*L**3/(3*EI)*1e3:.1f} mm; uniform: {total*L**3/(8*EI)*1e3:.1f} mm")

    # ---- 2. L-bracket: Castigliano with bending and torsion (4130 tube 25.4 x 1.24 mm)
    E_st, G_st = 205e9, 80e9
    do, t = 0.0254, 0.00124
    I_t = np.pi / 64 * (do ** 4 - (do - 2 * t) ** 4)
    EIb, GJb = E_st * I_t, G_st * 2 * I_t
    a, b, P = 0.30, 0.50, 500.0
    mem = l_bracket_members(a, b, EIb, GJb)
    F = flexibility_matrix(mem)
    U = lambda Pv: strain_energy_members(mem, Pv)
    d_cast = castigliano_deflection(U, np.array([P, 0.0]), 0)
    print(f"\nL-bracket a = {a} m, b = {b} m, P = {P:.0f} N: closed form {l_bracket_tip_deflection(P, a, b, EIb, GJb)*1e3:.3f} mm, "
          f"flexibility {F[0,0]*P*1e3:.3f} mm, dU/dP {d_cast*1e3:.3f} mm")
    parts = [P * a ** 3 / (3 * EIb), P * b ** 3 / (3 * EIb), P * a ** 2 * b / GJb]
    print(f"  arm bending {parts[0]*1e3:.3f} mm, leg bending {parts[1]*1e3:.3f} mm, leg torsion {parts[2]*1e3:.3f} mm; "
          f"F12 = {F[0,1]*1e6:.4f}, F21 = {F[1,0]*1e6:.4f} mm/kN")

    # ---- 3. engine-mount truss: determinate and with one redundant bar
    nodes = np.array([[0.0, 0.0], [0.0, 0.30], [0.0, 0.60], [0.50, 0.25]])
    free = np.array([6, 7])
    load = np.array([2.0e3, -5.9e3])
    A_tube = np.pi / 4 * (0.01905 ** 2 - (0.01905 - 2 * 0.000889) ** 2)
    for bars, label in (([(0, 3), (2, 3)], "2 bars (determinate)"), ([(0, 3), (1, 3), (2, 3)], "3 bars (1 redundant)")):
        EA = np.full(len(bars), E_st * A_tube)
        Nf = truss_force_method(nodes, bars, EA, load, free)
        u, Ns = truss_stiffness_method(nodes, bars, EA, load, free)
        dv = truss_unit_load_displacement(nodes, bars, EA, Nf, free, 1)
        print(f"\nengine mount, {label}: bar forces force method " + ", ".join(f"{v/1e3:+.3f}" for v in Nf) +
              " kN; stiffness method " + ", ".join(f"{v/1e3:+.3f}" for v in Ns) + " kN")
        print(f"  vertical deflection: unit load {dv*1e3:.4f} mm, stiffness {u[1]*1e3:.4f} mm")

    # ---- 4. strut-braced spar
    a_s, h_s = 3.0, 1.2
    EA_spar, EA_strut = E_al * 1.2e-3, E_al * 3.0e-4
    qf = lambda xx: elliptic_lift(xx, L, total)
    r = braced_spar(L, a_s, h_s, EI, EA_spar, EA_strut, qf)
    print(f"\nbraced spar, strut at a = {a_s} m, h = {h_s} m: strut tension {r['S']/1e3:.2f} kN, vertical part {r['R']/1e3:.2f} kN; "
          f"root moment {r['root_moment']/1e3:.2f} kN m (unbraced {M[0]/1e3:.2f}); tip {r['tip']*1e3:.1f} mm (unbraced {w[-1]*1e3:.1f}); "
          f"compatibility residual {r['compatibility']:.1e}")
    rig = braced_spar(L, L, 1.0, EI, 1e30, 1e30, lambda xx: np.full_like(xx, 1.0e3))
    print(f"  rigid strut at the tip, uniform 1 kN/m: R = {rig['R']/1e3:.4f} kN against 3 q L / 8 = {3*1e3*L/8/1e3:.4f} kN")

    # ---- experiment: strut stiffness and strut position
    print("\nEA_strut / EA_ref   S [kN]    R [kN]   root M [kN m]   tip [mm]")
    for f in (0.0, 0.01, 0.1, 1.0, 10.0, 1e6):
        rr = braced_spar(L, a_s, h_s, EI, EA_spar, max(f, 1e-12) * EA_strut, qf)
        print(f"{f:18g}  {rr['S']/1e3:7.2f}  {rr['R']/1e3:7.2f}  {rr['root_moment']/1e3:14.2f}  {rr['tip']*1e3:8.1f}")
    print("\nstrut position a [m]   S [kN]   R [kN]   root M [kN m]   max |M| [kN m]   tip [mm]")
    for aa in (1.5, 2.0, 2.5, 3.0, 3.5, 4.0):
        rr = braced_spar(L, aa, h_s, EI, EA_spar, EA_strut, qf)
        print(f"{aa:20.1f}  {rr['S']/1e3:7.2f}  {rr['R']/1e3:7.2f}  {rr['root_moment']/1e3:14.2f}  "
              f"{np.max(np.abs(rr['M']))/1e3:14.2f}  {rr['tip']*1e3:8.1f}")

    # ---- figure
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    ax[0].plot(x, M / 1e3, label="unbraced")
    ax[0].plot(r["x"], r["M"] / 1e3, label="strut-braced")
    ax[0].set_xlabel("x [m]"); ax[0].set_ylabel("bending moment [kN m]"); ax[0].legend(); ax[0].grid(alpha=0.3)
    ax[0].set_title("spar under elliptic lift")
    ax[1].plot(x, w * 1e3, label="unbraced, direct integration")
    ax[1].plot(r["x"], r["w"] * 1e3, label="strut-braced")
    ax[1].plot([L], [d_ul * 1e3], "ko", label="unit-load tip deflection")
    ax[1].set_xlabel("x [m]"); ax[1].set_ylabel("deflection [mm]"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    ax[1].set_title("deflection")
    fig.tight_layout()
    fig.savefig(out_dir / "energy_methods.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'energy_methods.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
