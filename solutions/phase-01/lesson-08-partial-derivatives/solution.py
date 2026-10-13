"""P01.08 Multivariable Functions & Partial Derivatives: reference solution.

Lift of a finite wing as a function of two variables, angle of attack alpha and airspeed V:

    L(alpha, V) = 1/2 rho V**2 S C_L(alpha),
    C_L(alpha)  = a (alpha - alpha0) - k3 (alpha - alpha0)**3,     [alpha in rad]

with the cubic term modelling the roll-off of lift toward stall. The module computes the
partial derivatives dL/dalpha and dL/dV analytically, by central finite differences and by
SymPy, the tangent-plane (total differential) approximation, mixed partials (Clairaut's
theorem), the level curve L = W (the trim relation V(alpha)) with its implicit derivative, and
the chain-rule time derivative of lift along a manoeuvre alpha(t), V(t).

Run:  python solution.py          (prints the sensitivity table, writes lift_sensitivity.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-01/lesson-08-partial-derivatives
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import sympy as sp


@dataclass(frozen=True)
class Wing:
    """Wing and air parameters in SI units.

    rho : air density, kg/m**3
    S : wing area, m**2
    a : lift-curve slope, 1/rad
    alpha0 : zero-lift angle of attack, rad
    k3 : cubic stall coefficient, 1/rad**3 (0 gives the linear model)
    """

    rho: float = 1.225
    S: float = 16.0
    a: float = 5.0
    alpha0: float = np.radians(-2.0)
    k3: float = 12.0


DEFAULT_WING = Wing()


def lift_coefficient(alpha: np.ndarray | float, wing: Wing = DEFAULT_WING) -> np.ndarray | float:
    """Return C_L(alpha) = a (alpha - alpha0) - k3 (alpha - alpha0)**3, vectorized over ``alpha`` (rad)."""
    d = np.asarray(alpha, dtype=float) - wing.alpha0
    out = wing.a * d - wing.k3 * d**3
    return float(out) if out.ndim == 0 else out


def lift(alpha: np.ndarray | float, V: np.ndarray | float, wing: Wing = DEFAULT_WING) -> np.ndarray | float:
    """Return the lift L = 1/2 rho V**2 S C_L(alpha) in newtons; ``alpha`` and ``V`` broadcast together.

    Raise ``ValueError`` if any V is negative.
    """
    V = np.asarray(V, dtype=float)
    if np.any(V < 0.0):
        raise ValueError("airspeed must be non-negative")
    out = 0.5 * wing.rho * V**2 * wing.S * lift_coefficient(alpha, wing)
    return float(out) if np.ndim(out) == 0 else out


def lift_partials(alpha: np.ndarray | float, V: np.ndarray | float, wing: Wing = DEFAULT_WING) -> tuple[np.ndarray | float, np.ndarray | float]:
    """Return (dL/dalpha, dL/dV) analytically: (1/2 rho V**2 S C_L'(alpha), rho V S C_L(alpha)).

    Units: N/rad and N/(m/s).
    """
    alpha = np.asarray(alpha, dtype=float)
    V = np.asarray(V, dtype=float)
    d = alpha - wing.alpha0
    dcl = wing.a - 3.0 * wing.k3 * d**2
    dl_da = 0.5 * wing.rho * V**2 * wing.S * dcl
    dl_dv = wing.rho * V * wing.S * lift_coefficient(alpha, wing)
    if np.ndim(dl_da) == 0:
        return float(dl_da), float(dl_dv)
    return dl_da, dl_dv


def partial_derivative(f: Callable[[np.ndarray], float], x: np.ndarray, i: int, h: float) -> float:
    """Return the central-difference estimate of df/dx_i at the point ``x`` (array), step ``h``.

    Only coordinate ``i`` is perturbed; every other coordinate is held fixed, which is the
    definition of a partial derivative. Raise ``ValueError`` if h <= 0 or i is out of range.
    """
    x = np.asarray(x, dtype=float)
    if h <= 0.0:
        raise ValueError("h must be positive")
    if not 0 <= i < x.size:
        raise ValueError("coordinate index out of range")
    e = np.zeros_like(x)
    e[i] = h
    return float((f(x + e) - f(x - e)) / (2.0 * h))


def mixed_partial(f: Callable[[np.ndarray], float], x: np.ndarray, i: int, j: int, h: float) -> float:
    """Return the central-difference estimate of d2f/(dx_i dx_j) at ``x``: the partial in j of the partial in i."""
    return partial_derivative(lambda y: partial_derivative(f, y, i, h), x, j, h)


def lift_partials_sympy() -> tuple[sp.Expr, sp.Expr, dict[str, sp.Symbol]]:
    """Return (dL/dalpha, dL/dV) as SymPy expressions plus the symbol table.

    Symbols: alpha, V, rho, S, a, alpha0, k3 (all real; V, rho, S positive).
    """
    alpha, alpha0, k3 = sp.symbols("alpha alpha0 k3", real=True)
    V, rho, S, a = sp.symbols("V rho S a", positive=True)
    L = sp.Rational(1, 2) * rho * V**2 * S * (a * (alpha - alpha0) - k3 * (alpha - alpha0) ** 3)
    syms = {"alpha": alpha, "V": V, "rho": rho, "S": S, "a": a, "alpha0": alpha0, "k3": k3}
    return sp.diff(L, alpha), sp.diff(L, V), syms


def lift_linearized(alpha: np.ndarray | float, V: np.ndarray | float, alpha_ref: float, V_ref: float, wing: Wing = DEFAULT_WING) -> np.ndarray | float:
    """Return the tangent-plane approximation L_ref + L_alpha (alpha - alpha_ref) + L_V (V - V_ref)."""
    la, lv = lift_partials(alpha_ref, V_ref, wing)
    out = lift(alpha_ref, V_ref, wing) + la * (np.asarray(alpha, float) - alpha_ref) + lv * (np.asarray(V, float) - V_ref)
    return float(out) if np.ndim(out) == 0 else out


def relative_sensitivities(alpha: float, V: float, wing: Wing = DEFAULT_WING) -> tuple[float, float]:
    """Return the logarithmic sensitivities (dlnL/dln(alpha - alpha0), dlnL/dlnV).

    They say by what percentage L changes for a 1 % change in the effective angle of attack
    or in V: exactly 1 and 2 for the linear model, less than 1 in alpha as stall approaches.
    """
    la, lv = lift_partials(alpha, V, wing)
    L = lift(alpha, V, wing)
    return float(la * (alpha - wing.alpha0) / L), float(lv * V / L)


def level_curve_speed(alpha: np.ndarray | float, L_target: float, wing: Wing = DEFAULT_WING) -> np.ndarray | float:
    """Return V(alpha) on the level curve L(alpha, V) = L_target: V = sqrt(2 L_target / (rho S C_L(alpha))).

    This is the trim relation: the speed at which the wing carries exactly L_target. Raise
    ``ValueError`` where C_L(alpha) <= 0 (no positive speed can carry the load).
    """
    cl = np.asarray(lift_coefficient(alpha, wing), dtype=float)
    if np.any(cl <= 0.0):
        raise ValueError("C_L must be positive on the level curve")
    out = np.sqrt(2.0 * L_target / (wing.rho * wing.S * cl))
    return float(out) if out.ndim == 0 else out


def level_curve_slope(alpha: float, V: float, wing: Wing = DEFAULT_WING) -> float:
    """Return dV/dalpha along the level curve through (alpha, V) by implicit differentiation: -L_alpha / L_V.

    Raise ``ValueError`` if dL/dV vanishes (the level curve is vertical there).
    """
    la, lv = lift_partials(alpha, V, wing)
    if lv == 0.0:
        raise ValueError("dL/dV = 0: implicit derivative undefined")
    return float(-la / lv)


def lift_rate(alpha: float, V: float, dalpha_dt: float, dV_dt: float, wing: Wing = DEFAULT_WING) -> float:
    """Return dL/dt along a manoeuvre by the chain rule: L_alpha dalpha/dt + L_V dV/dt (N/s)."""
    la, lv = lift_partials(alpha, V, wing)
    return float(la * dalpha_dt + lv * dV_dt)


def stall_angle(wing: Wing = DEFAULT_WING) -> float:
    """Return the angle of attack (rad) of maximum C_L: where dC_L/dalpha = 0, alpha0 + sqrt(a / (3 k3)).

    Raise ``ValueError`` for the linear model (k3 = 0), which has no maximum.
    """
    if wing.k3 <= 0.0:
        raise ValueError("the linear model has no stall angle")
    return float(wing.alpha0 + np.sqrt(wing.a / (3.0 * wing.k3)))


def tangent_plane_error_study(alpha_ref: float, V_ref: float, scales: list[float], wing: Wing = DEFAULT_WING) -> np.ndarray:
    """Return max |L - L_linearized| over a square of half-width (scale * 1 deg, scale * 1 m/s) around the reference.

    Halving the scale should divide the error by four (the remainder of the tangent plane is quadratic).
    """
    errs = []
    for s in scales:
        da = np.linspace(-s * np.radians(1.0), s * np.radians(1.0), 11)
        dv = np.linspace(-s, s, 11)
        A, Vg = np.meshgrid(alpha_ref + da, V_ref + dv, indexing="ij")
        errs.append(np.max(np.abs(lift(A, Vg, wing) - lift_linearized(A, Vg, alpha_ref, V_ref, wing))))
    return np.array(errs)


def main(out_dir: str | Path = ".") -> None:
    """Print the lift sensitivity table and save lift_sensitivity.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    w = DEFAULT_WING
    alpha_ref, V_ref = np.radians(4.0), 55.0
    L0 = lift(alpha_ref, V_ref, w)
    la, lv = lift_partials(alpha_ref, V_ref, w)
    print(f"Wing: S = {w.S} m^2, a = {w.a} /rad, alpha0 = {np.degrees(w.alpha0):.1f} deg, k3 = {w.k3} /rad^3, rho = {w.rho} kg/m^3")
    print(f"Reference point: alpha = {np.degrees(alpha_ref):.1f} deg, V = {V_ref} m/s  ->  C_L = {lift_coefficient(alpha_ref, w):.4f}, L = {L0:.1f} N (= {L0/9.80665:.0f} kg)")
    print(f"  dL/dalpha = {la:.1f} N/rad = {la*np.pi/180:.1f} N/deg;  dL/dV = {lv:.1f} N/(m/s)")
    f = lambda x: lift(x[0], x[1], w)
    x0 = np.array([alpha_ref, V_ref])
    for h in (1e-2, 1e-3, 1e-4):
        fa = partial_derivative(f, x0, 0, h)
        fv = partial_derivative(f, x0, 1, h)
        print(f"  central differences h = {h:.0e}: dL/dalpha err {abs(fa - la):.2e} N/rad, dL/dV err {abs(fv - lv):.2e} N/(m/s)")
    dla, dlv, syms = lift_partials_sympy()
    subs = {syms["alpha"]: alpha_ref, syms["V"]: V_ref, syms["rho"]: w.rho, syms["S"]: w.S, syms["a"]: w.a, syms["alpha0"]: w.alpha0, syms["k3"]: w.k3}
    print(f"  sympy: dL/dalpha = {float(dla.subs(subs)):.1f}, dL/dV = {float(dlv.subs(subs)):.1f}")
    print(f"  mixed partials: L_alpha_V = {mixed_partial(f, x0, 0, 1, 1e-3):.3f}, L_V_alpha = {mixed_partial(f, x0, 1, 0, 1e-3):.3f} N/(rad m/s)")
    sa, sv = relative_sensitivities(alpha_ref, V_ref, w)
    print(f"  1 % more (alpha - alpha0) -> {sa:.3f} % more lift; 1 % more V -> {sv:.3f} % more lift")
    print(f"  stall angle (C_L max) = {np.degrees(stall_angle(w)):.2f} deg, C_L max = {lift_coefficient(stall_angle(w), w):.3f}")
    print(f"  trim: V(alpha) on L = {L0:.0f} N at alpha = 2, 4, 8 deg: {np.round(level_curve_speed(np.radians([2.0, 4.0, 8.0]), L0, w), 2)} m/s")
    print(f"  implicit dV/dalpha at reference = {level_curve_slope(alpha_ref, V_ref, w):.2f} (m/s)/rad = {level_curve_slope(alpha_ref, V_ref, w)*np.pi/180:.3f} (m/s)/deg")
    print(f"  gust: alpha + 5 deg at fixed V: exact dL = {lift(alpha_ref + np.radians(5), V_ref, w) - L0:+.0f} N, linear {la*np.radians(5):+.0f} N")
    print(f"  pull-up: dalpha/dt = 2 deg/s, dV/dt = -1 m/s^2: dL/dt = {lift_rate(alpha_ref, V_ref, np.radians(2.0), -1.0, w):.1f} N/s")
    errs = tangent_plane_error_study(alpha_ref, V_ref, [4.0, 2.0, 1.0, 0.5], w)
    print("  tangent-plane error over shrinking boxes:", np.round(errs, 3), "ratios", np.round(errs[:-1] / errs[1:], 2))

    A = np.radians(np.linspace(-2.0, 14.0, 161))
    Vg = np.linspace(30.0, 90.0, 121)
    AA, VV = np.meshgrid(A, Vg, indexing="ij")
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))
    cs = axs[0].contour(np.degrees(AA), VV, lift(AA, VV, w) / 1e3, levels=np.arange(0, 60, 5), cmap="viridis")
    axs[0].clabel(cs, fmt="%.0f kN", fontsize=7)
    axs[0].plot(np.degrees(A[lift_coefficient(A, w) > 0]), level_curve_speed(A[lift_coefficient(A, w) > 0], L0, w), "r-", label=f"level curve L = {L0/1e3:.1f} kN (trim)")
    axs[0].plot(np.degrees(alpha_ref), V_ref, "ko", label="reference")
    axs[0].axvline(np.degrees(stall_angle(w)), color="gray", ls="--", label="C_L max")
    axs[0].set_xlabel("angle of attack [deg]")
    axs[0].set_ylabel("airspeed [m/s]")
    axs[0].set_ylim(30, 90)
    axs[0].set_title("Lift L(alpha, V): level curves")
    axs[0].legend(fontsize=7, loc="upper right")
    da = np.radians(np.linspace(-6.0, 6.0, 121))
    axs[1].plot(np.degrees(da), lift(alpha_ref + da, V_ref, w) / 1e3, label="exact, V fixed")
    axs[1].plot(np.degrees(da), lift_linearized(alpha_ref + da, V_ref, alpha_ref, V_ref, w) / 1e3, "--", label="tangent plane")
    axs[1].set_xlabel("alpha - alpha_ref [deg]")
    axs[1].set_ylabel("lift [kN]")
    axs[1].set_title("Linearization at the reference point")
    axs[1].grid(True)
    axs[1].legend()
    fig.tight_layout()
    fig.savefig(out_dir / "lift_sensitivity.png", dpi=150)
    print(f"saved {out_dir / 'lift_sensitivity.png'}")


if __name__ == "__main__":
    main()
