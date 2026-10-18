"""P08.05 Dimensional Analysis & Similarity: starter code.

Can a 1:5 model of the light aeroplane's wing be tested in a wind tunnel at the Reynolds number
and the Mach number of the full-size wing in cruise, and if not, what would it take? This module
is a Buckingham Pi tool and a wind-tunnel scaling calculator:

* the dimension matrix D (rows: the base dimensions M, L, T, Theta; columns: the variables) of a
  physical relation; every dimensionless product of powers of the variables is an exponent vector
  e with D e = 0, so the dimensionless groups are the nullspace of D, of dimension n - rank D
  (the Pi theorem; `P01.14`);
* the nullspace found three ways: exactly in rational arithmetic with SymPy; in floating point
  from the singular value decomposition (`P01.18`), the right singular vectors whose singular
  values vanish; and by the classical method of repeating variables, which solves a small linear
  system for each group and gives the familiar groups (Reynolds, Mach, force coefficient);
* the conditions a model test needs: the speed for Reynolds-number similarity, the resulting
  Mach number, and the pressure or temperature of a tunnel that matches both.

The independent checks: each route's groups are verified dimensionless (D e = 0); the spans of the
SymPy and SVD bases are the same space; the Pi count is n - rank; Reynolds and Mach numbers do not
change when the inputs are converted to other units; and a pressurised-tunnel solution is
substituted back to check that both numbers match.

Units: SI (kg, m, s, K, Pa). Dimensions as exponents of (M, L, T, Theta).

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes similarity.png)
Test: pytest curriculum/phase-08/lesson-05-dimensional-analysis-similarity
"""
from __future__ import annotations

import sys
from fractions import Fraction
from math import gcd, lcm
from pathlib import Path

import numpy as np
import sympy

G0 = 9.80665
R_AIR = 287.05287
GAMMA = 1.4
T0, P0 = 288.15, 101325.0

BASE = ("M", "L", "T", "Theta")

# Dimensions of common variables as exponents of (M, L, T, Theta) (given).
DIMENSIONS: dict[str, tuple[int, int, int, int]] = {
    "F": (1, 1, -2, 0),        # force
    "D": (1, 1, -2, 0),        # drag force
    "p": (1, -1, -2, 0),       # pressure
    "dp": (1, -1, -2, 0),      # pressure difference
    "rho": (1, -3, 0, 0),      # density
    "V": (0, 1, -1, 0),        # speed
    "a": (0, 1, -1, 0),        # speed of sound
    "c": (0, 1, 0, 0),         # chord, a length
    "L": (0, 1, 0, 0),         # length
    "d": (0, 1, 0, 0),         # diameter
    "eps": (0, 1, 0, 0),       # roughness height
    "mu": (1, -1, -1, 0),      # dynamic viscosity
    "nu": (0, 2, -1, 0),       # kinematic viscosity
    "g": (0, 1, -2, 0),        # gravitational acceleration
    "omega": (0, 0, -1, 0),    # angular speed or frequency
    "f": (0, 0, -1, 0),        # frequency
    "sigma": (1, 0, -2, 0),    # surface tension
    "P": (1, 2, -3, 0),        # power
    "T": (0, 0, 0, 1),         # temperature
    "cp": (0, 2, -2, -1),      # specific heat
    "k": (1, 1, -3, -1),       # thermal conductivity
    "R": (0, 2, -2, -1),       # gas constant
}


# ---------------------------------------------------------------- the dimension matrix and its nullspace

def dimension_matrix(names: list[str], dims: dict[str, tuple[int, ...]] | None = None) -> np.ndarray:
    """Dimension matrix D (4 x n, float) with column j the (M, L, T, Theta) exponents of names[j]."""
    raise NotImplementedError


def pi_count(D: np.ndarray) -> int:
    """Number of independent dimensionless groups, n - rank(D) (the Pi theorem)."""
    raise NotImplementedError


def pi_groups_sympy(D: np.ndarray) -> list[list[Fraction]]:
    """A basis of the nullspace of D in exact rational arithmetic (SymPy ``Matrix.nullspace``), each
    vector scaled to the smallest integers with a positive first nonzero entry. Each vector holds the
    exponents of one dimensionless group."""
    raise NotImplementedError


def pi_groups_svd(D: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    """An orthonormal basis (columns) of the nullspace of D from the SVD D = U S V^T: the columns of
    V whose singular values are below tol times the largest. Floating point, not unique, not
    integer, but the same space as the SymPy basis."""
    raise NotImplementedError


def same_span(B1: np.ndarray, B2: np.ndarray, tol: float = 1e-9) -> bool:
    """True when the column spaces of B1 and B2 coincide: each has the other's rank and stacking
    them does not raise it."""
    raise NotImplementedError


def is_dimensionless(D: np.ndarray, e: np.ndarray, tol: float = 1e-12) -> bool:
    """True when the product of the variables raised to the exponents e has no dimension (D e = 0)."""
    raise NotImplementedError


def pi_groups_repeating(D: np.ndarray, names: list[str], repeating: list[str]) -> dict[str, np.ndarray]:
    """The method of repeating variables: for each variable q not in ``repeating``, find exponents
    x of the repeating variables such that q * prod(r_i^x_i) is dimensionless, by solving
    D_r x = -D_q (least squares on the base dimensions that occur). The repeating variables must be
    dimensionally independent (D_r of full column rank) and must span the dimensions of the others.
    Returns {q: exponent vector over all names}, with exponent 1 on q."""
    raise NotImplementedError


def format_group(names: list[str], e: np.ndarray) -> str:
    """A readable product such as 'rho^1 V^1 c^1 mu^-1' (zero exponents omitted; given)."""
    parts = []
    for n, x in zip(names, e):
        x = float(x)
        if abs(x) > 1e-12:
            parts.append(f"{n}^{Fraction(x).limit_denominator(12)}")
    return " ".join(parts)


# ---------------------------------------------------------------- air properties (self-contained helpers, given)

def isa_troposphere(h: float) -> tuple[float, float, float]:
    """ISA (T, p, rho) for 0 <= h <= 11 km, the closed form of `P08.01`."""
    T = T0 - 0.0065 * h
    p = P0 * (T / T0) ** (G0 / (0.0065 * R_AIR))
    return T, p, p / (R_AIR * T)


def sutherland(T: float) -> float:
    """Viscosity of air (Pa s), Sutherland's law of `P08.01`."""
    return 1.458e-6 * T ** 1.5 / (T + 110.4)


def speed_of_sound(T: float) -> float:
    """a = sqrt(gamma R T) for air (m/s) (`P08.17`)."""
    return float(np.sqrt(GAMMA * R_AIR * T))


# ---------------------------------------------------------------- similarity numbers and model tests

def reynolds(rho: float, V: float, L: float, mu: float) -> float:
    """Reynolds number rho V L / mu: inertia over viscous forces."""
    raise NotImplementedError


def mach(V: float, a: float) -> float:
    """Mach number V / a: flow speed over the speed of sound."""
    raise NotImplementedError


def full_scale_numbers(h: float, V: float, c: float) -> dict[str, float]:
    """Reynolds and Mach numbers of a wing of chord c at true airspeed V and ISA altitude h."""
    raise NotImplementedError


def speed_for_reynolds(Re: float, c_model: float, rho: float, mu: float) -> float:
    """Tunnel speed (m/s) that gives a model of chord c_model the Reynolds number Re."""
    raise NotImplementedError


def pressurised_tunnel(Re: float, M: float, c_model: float, T_tunnel: float) -> dict[str, float]:
    """Tunnel conditions that match both Re and M on a model of chord c_model at stagnation-free
    test-section temperature T_tunnel: Mach fixes V = M a(T); then Reynolds fixes the density,
    rho = Re mu(T) / (V c), and the ideal-gas law the pressure p = rho R T."""
    raise NotImplementedError


def cryogenic_tunnel(Re: float, M: float, c_model: float, p_tunnel: float) -> dict[str, float]:
    """Test-section temperature that matches both Re and M on a model of chord c_model at pressure
    p_tunnel: Re(T) = (p/(R T)) M a(T) c / mu(T) falls steadily with T, so bisect it on
    [20 K, 400 K] (scipy.optimize.brentq). Raises ValueError when no temperature in range works."""
    raise NotImplementedError


def stokes_drag(mu: float, V: float, d: float) -> float:
    """Drag on a small sphere in creeping flow (Re << 1), F = 3 pi mu V d (given, `P08.16`)."""
    return 3.0 * np.pi * mu * V * d


def drag_coefficient(F: float, rho: float, V: float, A: float) -> float:
    """C_D = F / (rho V^2 A / 2)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the lesson's run

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save similarity.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    names = ["F", "rho", "V", "c", "mu", "a"]
    D = dimension_matrix(names)
    print(f"wing force F = f(rho, V, c, mu, a): n = {len(names)}, rank = {np.linalg.matrix_rank(D)}, "
          f"{pi_count(D)} groups")
    sy = pi_groups_sympy(D)
    print("  SymPy nullspace (integers):")
    for e in sy:
        print(f"    {format_group(names, np.array(e, dtype=float)):40s} D e = {(D @ np.array(e, dtype=float)).astype(int)}")
    B = pi_groups_svd(D)
    print(f"  SVD nullspace: {B.shape[1]} orthonormal vectors, |D B| max {np.max(np.abs(D @ B)):.1e}, "
          f"same span as SymPy: {same_span(B, np.array(sy, dtype=float).T)}")
    rg = pi_groups_repeating(D, names, ["rho", "V", "c"])
    print("  repeating variables rho, V, c:")
    for q, e in rg.items():
        print(f"    Pi_{q}: {format_group(names, e)}")
    print(f"  same span as SymPy: {same_span(np.array(list(rg.values())).T, np.array(sy, dtype=float).T)}")

    names_p = ["dp", "rho", "V", "d", "mu", "eps", "L"]
    Dp = dimension_matrix(names_p)
    print(f"\npipe pressure drop dp = f(rho, V, d, mu, eps, L): {pi_count(Dp)} groups:")
    for q, e in pi_groups_repeating(Dp, names_p, ["rho", "V", "d"]).items():
        print(f"    Pi_{q}: {format_group(names_p, e)}")

    h, V, c = 3000.0, 60.0, 1.5
    fs = full_scale_numbers(h, V, c)
    print(f"\nlight aeroplane wing: chord {c} m, TAS {V} m/s at {h:.0f} m: Re = {fs['Re']:.4e}, M = {fs['M']:.4f}")
    rho0, mu0, a0 = 1.225, sutherland(T0), speed_of_sound(T0)
    print("model tests in a sea-level atmospheric tunnel (rho 1.225, mu %.4e, a %.2f):" % (mu0, a0))
    print("   scale   chord (m)   V for Re match (m/s)   tunnel Mach   Re at full-scale Mach")
    for k in (1.0, 2.0, 5.0, 10.0):
        cm = c / k
        Vr = speed_for_reynolds(fs["Re"], cm, rho0, mu0)
        Re_M = reynolds(rho0, fs["M"] * a0, cm, mu0)
        print(f"   1:{k:<4.0f}   {cm:6.3f}        {Vr:9.1f}             {mach(Vr, a0):7.3f}        {Re_M:.3e}")
    pt = pressurised_tunnel(fs["Re"], fs["M"], c / 5, T0)
    print(f"1:5 model, both matched, pressurised tunnel at {T0} K: V = {pt['V']:.2f} m/s, rho = {pt['rho']:.3f} kg/m^3, "
          f"p = {pt['p'] / 1e5:.3f} bar ({pt['p'] / P0:.2f} atm)")
    chk = reynolds(pt["rho"], pt["V"], c / 5, pt["mu"]), mach(pt["V"], speed_of_sound(T0))
    print(f"   substituted back: Re = {chk[0]:.6e}, M = {chk[1]:.6f}")
    ct = cryogenic_tunnel(fs["Re"], fs["M"], c / 5, P0)
    print(f"1:5 model, both matched, cryogenic tunnel at 1 atm: T = {ct['T']:.1f} K, V = {ct['V']:.2f} m/s, "
          f"rho = {ct['rho']:.3f} kg/m^3 (Sutherland's law extrapolated far below ordinary temperatures: indicative only)")
    nu_w = 1.002e-3 / 998.2
    print(f"1:5 model in a water tunnel (20 C): V = {fs['Re'] * nu_w / (c / 5):.2f} m/s for the Re match (Mach irrelevant)")

    ks = np.linspace(1, 10, 91)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].plot(ks, [mach(speed_for_reynolds(fs["Re"], c / k, rho0, mu0), a0) for k in ks], label="tunnel Mach at Re match")
    ax[0].axhline(fs["M"], color="k", ls="--", label="full-scale Mach")
    ax[0].axhline(1.0, color="r", ls=":", label="M = 1")
    ax[0].set_xlabel("model scale 1:k")
    ax[0].set_ylabel("Mach number")
    ax[0].set_title("Re and Mach cannot both match at 1 atm, 288 K")
    ax[0].legend(fontsize=8)
    ax[1].plot(ks, [pressurised_tunnel(fs["Re"], fs["M"], c / k, T0)["p"] / 1e5 for k in ks])
    ax[1].set_xlabel("model scale 1:k")
    ax[1].set_ylabel("tunnel pressure for both (bar)")
    ax[1].set_title("pressurised tunnel, 288.15 K")
    fig.tight_layout()
    fig.savefig(out_dir / "similarity.png", dpi=110)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'similarity.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
