"""P05.15 Column Buckling: starter code.

A slender member in compression does not fail by crushing; at a critical load the straight shape stops
being the only equilibrium and the member bows sideways. This module builds the column curve an
aerospace stress engineer uses for struts, longerons and landing-gear braces:

* the Euler load P_cr = pi^2 E I / (K L)^2 for the four classic end conditions, through the effective
  length factor K;
* the Johnson parabola for short and intermediate columns, which meets the Euler hyperbola tangentially
  at the transition slenderness sqrt(2 pi^2 E / F_cy) and caps the stress at F_cy;
* the tangent-modulus (Engesser) critical stress from a Ramberg-Osgood compression curve, the physics
  the Johnson parabola approximates;
* an independent check: the buckling load as the smallest eigenvalue of the finite-difference form of
  E I w'''' + P w'' = 0 with the boundary conditions of each end, converging at second order to the
  closed forms.

The worked case is a 4130 steel engine-mount strut; the margin is MS = P_cr / (1.5 P_limit) - 1.

Units: SI throughout (N, m, Pa).

Run:  python main.py [out_dir]     (prints the tables, writes column_buckling.png)
Test: pytest curriculum/phase-05/lesson-15-column-buckling
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eig
from scipy.optimize import brentq

# Theoretical effective-length factors (ideal ends). Design practice uses larger values for "fixed"
# ends because no real joint is perfectly rigid (for example 0.65 instead of 0.5, 0.8 instead of 0.699).
# fixed-pinned: K = pi / beta with beta = 4.4934... the smallest positive root of tan(beta) = beta.
BETA_FIXED_PINNED = 4.493409457909064
K_FACTORS = {"pinned-pinned": 1.0, "fixed-free": 2.0, "fixed-pinned": np.pi / BETA_FIXED_PINNED, "fixed-fixed": 0.5}

# Typical room-temperature compression properties (not design allowables), Pa.
E_4130, FCY_4130 = 205e9, 435e6
E_7075, FCY_7075 = 71.7e9, 503e6
E_2024, FCY_2024 = 73.1e9, 345e6


# ---------------------------------------------------------------- sections

def tube_properties(d_o: float, t: float) -> dict[str, float]:
    """Area A, second moment I and radius of gyration r = sqrt(I/A) of a round tube."""
    d_i = d_o - 2.0 * t
    A = np.pi / 4.0 * (d_o ** 2 - d_i ** 2)
    I = np.pi / 64.0 * (d_o ** 4 - d_i ** 4)
    return {"A": A, "I": I, "r": np.sqrt(I / A)}


# ---------------------------------------------------------------- Euler and Johnson

def effective_length_factor(ends: str) -> float:
    """K for the end condition name: 'pinned-pinned' 1, 'fixed-free' 2, 'fixed-pinned' 0.6992,
    'fixed-fixed' 0.5 (ideal ends)."""
    raise NotImplementedError


def euler_load(E: float, I: float, L: float, K: float = 1.0) -> float:
    """Euler critical load pi^2 E I / (K L)^2 [N]."""
    raise NotImplementedError


def euler_stress(E: float, slenderness: np.ndarray | float) -> np.ndarray | float:
    """Euler critical stress pi^2 E / (K L / r)^2 [Pa] as a function of the slenderness ratio KL/r."""
    raise NotImplementedError


def transition_slenderness(E: float, F_cy: float) -> float:
    """Slenderness at which the Johnson parabola meets the Euler curve: sqrt(2 pi^2 E / F_cy).
    There the Euler stress is F_cy / 2."""
    raise NotImplementedError


def johnson_stress(E: float, F_cy: float, slenderness: np.ndarray | float) -> np.ndarray | float:
    """Johnson parabola F_cy - (F_cy / (2 pi))^2 (KL/r)^2 / E [Pa]. Valid below the transition
    slenderness; equals F_cy at zero slenderness and is tangent to Euler at the transition."""
    raise NotImplementedError


def critical_stress(E: float, F_cy: float, slenderness: np.ndarray | float) -> np.ndarray | float:
    """Column curve: Johnson below the transition slenderness, Euler above it."""
    raise NotImplementedError


def strut_margin(P_cr: float, P_limit: float, fs: float = 1.5) -> float:
    """Margin of safety of a strut in compression: P_cr / (fs * P_limit) - 1."""
    raise NotImplementedError


# ---------------------------------------------------------------- tangent modulus

def ramberg_osgood_strain(sigma: np.ndarray | float, E: float, F_cy: float, n: float) -> np.ndarray | float:
    """Ramberg-Osgood compression curve in the yield-anchored form: epsilon = sigma/E + 0.002 (sigma/F_cy)^n.
    At sigma = F_cy the plastic strain is 0.2 %, the offset definition of yield."""
    raise NotImplementedError


def tangent_modulus(sigma: np.ndarray | float, E: float, F_cy: float, n: float) -> np.ndarray | float:
    """E_t = d sigma / d epsilon of the Ramberg-Osgood curve: 1 / (1/E + 0.002 n sigma^(n-1) / F_cy^n)."""
    raise NotImplementedError


def tangent_modulus_stress(slenderness: float, E: float, F_cy: float, n: float) -> float:
    """Engesser tangent-modulus critical stress: the root of sigma = pi^2 E_t(sigma) / (KL/r)^2."""
    raise NotImplementedError


# ---------------------------------------------------------------- finite-difference eigenproblem

def fd_buckling(L: float, EI: float, ends: str, n: int = 200) -> tuple[float, np.ndarray, np.ndarray]:
    """Smallest buckling load of a uniform column by finite differences.

    E I w'''' + P w'' = 0 on 0 <= x <= L, n intervals of h = L/n, nodes x_0..x_n, two ghost nodes beyond
    each end (n + 5 unknowns). The second-order central differences D4 and D2 are applied at every node
    x_0..x_n (n + 1 equations); the four boundary conditions close the system:
      pinned:  w = 0, w'' = 0          fixed: w = 0, w' = 0
      free:    w'' = 0, E I w''' + P w' = 0  (no moment, no transverse shear; the axial load stays axial)
    The result is a generalised eigenproblem A w = (P / E I) B w; boundary rows have B = 0 (infinite
    eigenvalues, discarded) except the free-end shear row, whose P w' term goes into B.
    'ends' names the x = 0 end first: 'pinned-pinned', 'fixed-free', 'fixed-pinned', 'fixed-fixed'.

    Returns (P_cr, x, w) with w the buckled shape at the nodes, normalised to max |w| = 1."""
    raise NotImplementedError


def euler_mode_residual(x: np.ndarray, w: np.ndarray, P: float, EI: float) -> float:
    """Equilibrium check of a pinned-pinned buckled shape: the bending moment E I w'' must balance the
    moment of the axial load about the deflected section, P w (M = -P w, so E I w'' + P w = 0).
    Returns max |E I w'' + P w| / max |P w| over interior nodes (second-order central differences)."""
    raise NotImplementedError


# ---------------------------------------------------------------- driver

def main(out_dir: str | Path = ".") -> None:
    """Print the worked strut, the FD convergence and the column-curve tables; save column_buckling.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    MPa = 1e6

    # ---- the engine-mount strut
    d_o, t, L = 0.0222, 0.000889, 0.50
    p = tube_properties(d_o, t)
    s = L / p["r"]
    sc = transition_slenderness(E_4130, FCY_4130)
    print(f"4130 strut {d_o*1e3:.1f} x {t*1e3:.3f} mm, L = {L*1e3:.0f} mm: A = {p['A']*1e6:.2f} mm^2, I = {p['I']*1e12:.0f} mm^4, "
          f"r = {p['r']*1e3:.2f} mm, L/r = {s:.1f}; transition slenderness {sc:.1f}")
    Pe = euler_load(E_4130, p["I"], L)
    Pj = johnson_stress(E_4130, FCY_4130, s) * p["A"]
    print(f"pinned ends: Euler {Pe/1e3:.2f} kN ({Pe/p['A']/MPa:.0f} MPa), Johnson {Pj/1e3:.2f} kN "
          f"({Pj/p['A']/MPa:.0f} MPa); Euler overestimates by {100*(Pe/Pj-1):.0f} %")
    P_limit = 10.0e3
    print(f"limit compression {P_limit/1e3:.1f} kN, ultimate {1.5*P_limit/1e3:.1f} kN: MS (Johnson) = {strut_margin(Pj, P_limit):+.3f}, "
          f"MS (Euler, wrong here) = {strut_margin(Pe, P_limit):+.3f}")
    Ptm = tangent_modulus_stress(s, E_4130, FCY_4130, 15.0) * p["A"]
    print(f"tangent modulus (Ramberg-Osgood n = 15): {Ptm/1e3:.2f} kN ({Ptm/p['A']/MPa:.0f} MPa)")

    # ---- FD eigenvalue check for every end condition
    EI = E_4130 * p["I"]
    print("\nends            K       closed form [kN]   FD n=50     FD n=100    FD n=200    observed order")
    for ends, K in K_FACTORS.items():
        exact = euler_load(E_4130, p["I"], L, K)
        errs, vals = [], []
        for n in (50, 100, 200):
            P, _, _ = fd_buckling(L, EI, ends, n)
            vals.append(P)
            errs.append(abs(P - exact) / exact)
        order = np.log(errs[1] / errs[2]) / np.log(2.0)
        print(f"{ends:14s} {K:6.4f}   {exact/1e3:12.4f}   " + "  ".join(f"{v/1e3:10.4f}" for v in vals) + f"   {order:.2f}")

    # ---- experiment: the column curve against length
    print("\nL [mm]   L/r     Euler [MPa]  Johnson [MPa]  tangent-modulus [MPa]  column curve [MPa]")
    for Lm in (0.15, 0.25, 0.35, 0.50, 0.70, 0.90, 1.20):
        sl = Lm / p["r"]
        eu = euler_stress(E_4130, sl)
        jo = johnson_stress(E_4130, FCY_4130, sl)
        tm = tangent_modulus_stress(sl, E_4130, FCY_4130, 15.0)
        cc = critical_stress(E_4130, FCY_4130, sl)
        jtxt = f"{jo/MPa:13.0f}" if sl < sc else f"{'(n/a)':>13s}"
        print(f"{Lm*1e3:6.0f}  {sl:6.1f}  {eu/MPa:11.0f}  {jtxt}  {tm/MPa:21.0f}  {cc/MPa:18.0f}")

    # ---- materials: transition slenderness
    print("\nmaterial   E [GPa]  F_cy [MPa]  transition KL/r")
    for nm, E, F in (("2024-T3", E_2024, FCY_2024), ("7075-T6", E_7075, FCY_7075), ("4130", E_4130, FCY_4130)):
        print(f"{nm:9s} {E/1e9:7.1f}  {F/MPa:10.0f}  {transition_slenderness(E, F):15.1f}")

    # ---- figure
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    sl = np.linspace(5, 200, 400)
    ax[0].plot(sl, euler_stress(E_4130, sl) / MPa, "--", label="Euler")
    ax[0].plot(sl, critical_stress(E_4130, FCY_4130, sl) / MPa, label="Johnson + Euler")
    ax[0].plot(sl, [tangent_modulus_stress(v, E_4130, FCY_4130, 15.0) / MPa for v in sl], ":", label="tangent modulus, n = 15")
    ax[0].axvline(sc, color="k", lw=0.6)
    ax[0].plot(s, Pj / p["A"] / MPa, "ko", label="engine-mount strut")
    ax[0].set_ylim(0, 1.3 * FCY_4130 / MPa); ax[0].set_xlabel("KL / r"); ax[0].set_ylabel("critical stress [MPa]")
    ax[0].set_title("4130 column curve"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
    for ends in K_FACTORS:
        _, x, w = fd_buckling(L, EI, ends, 100)
        ax[1].plot(x / L, w, label=ends)
    ax[1].set_xlabel("x / L"); ax[1].set_ylabel("buckled shape (normalised)")
    ax[1].set_title("first modes from the finite-difference eigenproblem"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "column_buckling.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'column_buckling.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
