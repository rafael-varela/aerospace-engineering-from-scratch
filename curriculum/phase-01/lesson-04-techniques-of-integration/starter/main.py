"""P01.04 Techniques of Integration: starter code.

Roll moment of inertia of a wing modelled as a tapered beam, for three planforms:
    linear taper     c(y) = c_r [1 - (1 - taper) y / s]          (polynomial: exact by Simpson)
    elliptic         c(y) = c_r sqrt(1 - (y/s)^2)                 (trigonometric substitution y = s sin(theta))
    exponential      c(y) = c_r exp(-y / ell)                     (integration by parts, twice; improper as s -> oo)
with mass per unit span lambda(y) = sigma c(y) (sigma = areal density, kg/m^2) and, for the wing pair,
    I_xx = 2 int_0^s y^2 sigma c(y) dy.
Also: the take-off time under thrust T and drag k v^2 (partial fractions), and the midpoint, trapezoid
and Simpson rules with their error bounds.

Complete every function marked NotImplementedError.

Run:  python main.py          (prints the inertia tables, writes beam_inertia_convergence.png)
Test: pytest curriculum/phase-01/lesson-04-techniques-of-integration
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
import sympy as sp

# The lesson's light-aircraft wing: span 11 m, root chord 1.6 m, taper 0.6, areal density 12 kg/m^2.
WING = {"semispan": 5.5, "root_chord": 1.6, "taper": 0.6, "sigma": 12.0}

# Shared symbols for the symbolic integrals.
y_s, s_s, ell_s, sigma_s, cr_s = sp.symbols("y s ell sigma c_r", positive=True)


def model_symbols() -> dict[str, sp.Symbol]:
    """Return the symbols: ``y`` (m), ``s`` semispan (m), ``ell`` (m), ``sigma`` (kg/m^2), ``cr`` root chord (m)."""
    return {"y": y_s, "s": s_s, "ell": ell_s, "sigma": sigma_s, "cr": cr_s}


# ---------------------------------------------------------------- planforms

def chord_linear(y: np.ndarray | float, c_r: float, taper: float, s: float) -> np.ndarray | float:
    """Return c(y) = c_r [1 - (1 - taper) y/s] for 0 <= y <= s (taper = tip chord / root chord)."""
    raise NotImplementedError


def chord_elliptic(y: np.ndarray | float, c_r: float, s: float) -> np.ndarray | float:
    """Return c(y) = c_r sqrt(1 - (y/s)^2) for 0 <= y <= s."""
    raise NotImplementedError


def chord_exponential(y: np.ndarray | float, c_r: float, ell: float) -> np.ndarray | float:
    """Return c(y) = c_r exp(-y / ell): the chord falls by a factor e every ``ell`` metres."""
    raise NotImplementedError


# ---------------------------------------------------------------- closed forms

def wing_mass_linear(sigma: float, c_r: float, taper: float, s: float) -> float:
    """Return the mass of the wing pair, 2 sigma int_0^s c dy = sigma c_r s (1 + taper), in kg."""
    raise NotImplementedError


def spanwise_centroid_linear(c_r: float, taper: float, s: float) -> float:
    """Return the spanwise centroid of one half-wing, y_bar = s (1 + 2 taper) / (3 (1 + taper)), in m."""
    raise NotImplementedError


def roll_inertia_linear(sigma: float, c_r: float, taper: float, s: float) -> float:
    """Return I_xx = 2 sigma c_r s^3 [1/3 - (1 - taper)/4] for the linearly tapered wing pair, in kg m^2."""
    raise NotImplementedError


def roll_inertia_elliptic(sigma: float, c_r: float, s: float) -> float:
    """Return I_xx = pi sigma c_r s^3 / 8 for the elliptic wing pair (y = s sin(theta): int sin^2 cos^2 = pi/16)."""
    raise NotImplementedError


def roll_inertia_exponential(sigma: float, c_r: float, ell: float, s: float) -> float:
    """Return I_xx = 2 sigma c_r [2 ell^3 - exp(-s/ell) (s^2 ell + 2 s ell^2 + 2 ell^3)] by parts (twice)."""
    raise NotImplementedError


def roll_inertia_exponential_infinite(sigma: float, c_r: float, ell: float) -> float:
    """Return the improper-integral limit s -> oo of ``roll_inertia_exponential``: 4 sigma c_r ell^3."""
    raise NotImplementedError


def symbolic_roll_inertia_exponential() -> sp.Expr:
    """Return 2 sigma c_r int_0^s y^2 exp(-y/ell) dy as a SymPy expression (SymPy does the parts for you)."""
    raise NotImplementedError


def symbolic_roll_inertia_elliptic() -> sp.Expr:
    """Return 2 sigma c_r int_0^s y^2 sqrt(1 - y^2/s^2) dy as a SymPy expression, after the substitution y = s sin(theta).

    With y = s sin(theta), dy = s cos(theta) dtheta and sqrt(1 - y^2/s^2) = cos(theta), so the integral becomes
    s^3 int_0^{pi/2} sin^2(theta) cos^2(theta) dtheta, a trigonometric integral SymPy evaluates at once
    (left in y, SymPy's general algorithm takes minutes).
    """
    raise NotImplementedError


# ---------------------------------------------------------------- numerical rules and their bounds

def midpoint(f: Callable[[np.ndarray], np.ndarray], a: float, b: float, n: int) -> float:
    """Return the midpoint rule with ``n`` panels. Raise ``ValueError`` if n < 1."""
    raise NotImplementedError


def trapezoid(f: Callable[[np.ndarray], np.ndarray], a: float, b: float, n: int) -> float:
    """Return the trapezoid rule with ``n`` panels: h [f_0/2 + f_1 + ... + f_{n-1} + f_n/2]. Raise ``ValueError`` if n < 1."""
    raise NotImplementedError


def simpson(f: Callable[[np.ndarray], np.ndarray], a: float, b: float, n: int) -> float:
    """Return Simpson's rule with ``n`` panels (n even): (h/3)[f_0 + 4 f_1 + 2 f_2 + ... + 4 f_{n-1} + f_n].

    Exact for every cubic. Raise ``ValueError`` if n < 2 or n is odd.
    """
    raise NotImplementedError


def trapezoid_error_bound(f2_max: float, a: float, b: float, n: int) -> float:
    """Return the bound (b - a)^3 max|f''| / (12 n^2) on the trapezoid error."""
    raise NotImplementedError


def simpson_error_bound(f4_max: float, a: float, b: float, n: int) -> float:
    """Return the bound (b - a)^5 max|f''''| / (180 n^4) on the Simpson error."""
    raise NotImplementedError


RULES: dict[str, Callable[..., float]] = {"midpoint": midpoint, "trapezoid": trapezoid, "simpson": simpson}


def roll_inertia_numeric(chord: Callable[[np.ndarray], np.ndarray], sigma: float, s: float, n: int, rule: str = "simpson") -> float:
    """Return 2 sigma int_0^s y^2 chord(y) dy by the named rule ('midpoint', 'trapezoid' or 'simpson')."""
    raise NotImplementedError


def roll_inertia_elliptic_substituted(sigma: float, c_r: float, s: float, n: int) -> float:
    """Return the elliptic I_xx with y = s sin(theta): 2 sigma c_r s^3 int_0^{pi/2} sin^2 cos^2 dtheta by the trapezoid rule.

    The substituted integrand is smooth and periodic, so the trapezoid rule converges far faster than any power of n.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- improper integrals

def improper_integral(f: Callable[[np.ndarray], np.ndarray], a: float, cutoffs: np.ndarray, n: int = 2000) -> np.ndarray:
    """Return int_a^R f dx for each R in ``cutoffs`` (Simpson, ``n`` panels): the sequence whose limit defines int_a^oo f."""
    raise NotImplementedError


# ---------------------------------------------------------------- partial fractions: take-off time

def takeoff_time(v: np.ndarray | float, thrust: float, k: float, mass: float) -> np.ndarray | float:
    """Return the time to accelerate from rest to speed ``v`` under m dv/dt = T - k v^2 (no rolling friction).

    Separating variables and splitting 1/(T - k v^2) by partial fractions gives
    t = (m / (2 sqrt(T k))) ln[(sqrt(T) + sqrt(k) v) / (sqrt(T) - sqrt(k) v)],  valid for v < sqrt(T/k).
    Raise ``ValueError`` if any v >= sqrt(T/k) or if T, k, m are not positive.
    """
    raise NotImplementedError


def takeoff_time_numeric(v: float, thrust: float, k: float, mass: float, n: int = 200) -> float:
    """Return the same time as ``takeoff_time`` by Simpson's rule on int_0^v m dv' / (T - k v'^2)."""
    raise NotImplementedError


# ---------------------------------------------------------------- driver

def main(out_dir: str | Path = ".") -> None:
    """Print mass, centroid and roll inertia for the three planforms with their checks, and save the convergence plot."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    s, c_r, taper, sigma = WING["semispan"], WING["root_chord"], WING["taper"], WING["sigma"]
    ell = s / np.log(1.0 / taper)        # exponential planform with the same tip chord
    lin = lambda y: chord_linear(y, c_r, taper, s)
    ell_fn = lambda y: chord_elliptic(y, c_r, s)
    exp_fn = lambda y: chord_exponential(y, c_r, ell)
    print(f"wing: b = {2 * s} m, c_r = {c_r} m, taper {taper}, sigma = {sigma} kg/m^2, exponential length ell = {ell:.3f} m")
    print(f"linear:  mass {wing_mass_linear(sigma, c_r, taper, s):.1f} kg, y_bar = {spanwise_centroid_linear(c_r, taper, s):.3f} m,"
          f" I_xx = {roll_inertia_linear(sigma, c_r, taper, s):.2f} kg m^2 (Simpson n = 2: {roll_inertia_numeric(lin, sigma, s, 2):.2f})")
    print(f"elliptic: I_xx = {roll_inertia_elliptic(sigma, c_r, s):.2f} kg m^2 (Simpson n = 1000: {roll_inertia_numeric(ell_fn, sigma, s, 1000):.2f},"
          f" substituted trapezoid n = 4: {roll_inertia_elliptic_substituted(sigma, c_r, s, 4):.4f})")
    print(f"exponential: I_xx = {roll_inertia_exponential(sigma, c_r, ell, s):.2f} kg m^2 (Simpson n = 100: {roll_inertia_numeric(exp_fn, sigma, s, 100):.2f});"
          f" infinite span limit {roll_inertia_exponential_infinite(sigma, c_r, ell):.1f} kg m^2")
    print("SymPy, exponential:", symbolic_roll_inertia_exponential())
    print("SymPy, elliptic:   ", symbolic_roll_inertia_elliptic())

    ns = np.array([4, 8, 16, 32, 64, 128, 256])
    exact_e = roll_inertia_elliptic(sigma, c_r, s)
    errs = {rule: np.array([abs(roll_inertia_numeric(ell_fn, sigma, s, int(n), rule) - exact_e) for n in ns]) for rule in RULES}
    errs["trapezoid, theta"] = np.array([abs(roll_inertia_elliptic_substituted(sigma, c_r, s, int(n)) - exact_e) for n in ns])
    print(f"\nelliptic planform, |I_xx error| in kg m^2:\n{'n':>5} " + " ".join(f"{k:>18}" for k in errs))
    for i, n in enumerate(ns):
        print(f"{n:5d} " + " ".join(f"{errs[k][i]:18.2e}" for k in errs))
    exact_x = roll_inertia_exponential(sigma, c_r, ell, s)
    f2 = 2.0 * sigma * c_r * 2.0  # |d^2/dy^2 (y^2 e^{-y/ell})| <= 2 near y = 0 (the maximum on [0, s])
    for n in (8, 32, 128):
        print(f"exponential, trapezoid n = {n}: error {abs(roll_inertia_numeric(exp_fn, sigma, s, n, 'trapezoid') - exact_x):.2e},"
              f" bound {trapezoid_error_bound(f2, 0.0, s, n):.2e}")

    cutoffs = np.array([1.0, 2.0, 3.0, 4.0, 6.0, 8.0]) * ell
    tail = improper_integral(lambda y: 2.0 * sigma * y**2 * exp_fn(y), 0.0, cutoffs)
    print("\nimproper integral: cutoff R/ell and I(R)/I(oo):", np.round(tail / roll_inertia_exponential_infinite(sigma, c_r, ell), 4))

    m, thrust, k = 1000.0, 2500.0, 0.5 * 1.225 * 0.04 * 15.0
    print(f"\ntake-off: m = {m} kg, T = {thrust} N, k = {k:.4f} kg/m, v_t = {np.sqrt(thrust / k):.1f} m/s")
    for v in (20.0, 30.0, 40.0):
        print(f"  to {v} m/s: {takeoff_time(v, thrust, k, m):.3f} s (Simpson {takeoff_time_numeric(v, thrust, k, m):.3f} s;"
              f" no drag {m * v / thrust:.3f} s)")

    fig, ax = plt.subplots(figsize=(6, 4))
    for key, marker in zip(errs, ("o", "s", "^", "x")):
        ax.loglog(ns, np.maximum(errs[key], 1e-16), marker=marker, label=key)
    ax.loglog(ns, errs["trapezoid"][0] * (ns / ns[0]) ** -1.5, "k:", lw=0.8, label=r"$n^{-1.5}$")
    ax.set_xlabel("panels n")
    ax.set_ylabel(r"$|I_{xx}$ error$|$ [kg m$^2$]")
    ax.set_title("Elliptic planform: raw rules vs trigonometric substitution")
    ax.grid(True, which="both", lw=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "beam_inertia_convergence.png", dpi=150)
    print(f"saved {out_dir / 'beam_inertia_convergence.png'}")


if __name__ == "__main__":
    main()
