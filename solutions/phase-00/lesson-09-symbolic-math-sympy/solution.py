"""P00.09 Symbolic Math with SymPy: reference solution.

Equation of motion of a point mass falling under gravity with quadratic drag, derived and
checked symbolically, then turned into fast NumPy functions with ``lambdify``.

Model (downward positive, object never moves upward):
    m dv/dt = m g - k v**2,     k = rho * C_D * A / 2        [k in kg/m]
Terminal velocity v_t = sqrt(m g / k); exact solution v(t) = v_t tanh(g t / v_t) for v(0) = 0.

Run:  python solution.py          (prints the derivation, writes fall_with_drag.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-00/lesson-09-symbolic-math-sympy
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
import sympy as sp

# One shared set of symbols so that expressions from different functions can be combined.
# ``positive=True`` lets SymPy simplify sqrt(g*m/k)**2 -> g*m/k without case analysis.
t, v, m, g, k = sp.symbols("t v m g k", positive=True)


def model_symbols() -> dict[str, sp.Symbol]:
    """Return the model's symbols by name: ``t`` (s), ``v`` (m/s), ``m`` (kg), ``g`` (m/s**2), ``k`` (kg/m)."""
    return {"t": t, "v": v, "m": m, "g": g, "k": k}


def drag_constant(rho: float, cd: float, area: float) -> float:
    """Return k = rho * C_D * A / 2 in kg/m. All inputs must be positive and finite."""
    for name, val in (("rho", rho), ("cd", cd), ("area", area)):
        if not np.isfinite(val) or val <= 0.0:
            raise ValueError(f"{name} must be positive and finite, got {val!r}")
    return 0.5 * rho * cd * area


def acceleration_expr() -> sp.Expr:
    """Return dv/dt = g - k v**2 / m as a SymPy expression in ``v``, ``m``, ``g``, ``k``.

    Derived from Newton's second law with weight m g downward and drag k v**2 upward,
    valid while v >= 0 (the object is falling).
    """
    return g - k * v**2 / m


def terminal_velocity_expr() -> sp.Expr:
    """Solve dv/dt = 0 for the positive root: v_t = sqrt(g m / k)."""
    roots = sp.solve(sp.Eq(acceleration_expr(), 0), v)
    (root,) = [r for r in roots if r.is_positive]
    return root


def velocity_expr() -> sp.Expr:
    """Return the exact solution v(t) = v_t tanh(g t / v_t) for v(0) = 0."""
    vt = terminal_velocity_expr()
    return vt * sp.tanh(g * t / vt)


def position_expr() -> sp.Expr:
    """Return x(t) = (v_t**2 / g) ln cosh(g t / v_t), the distance fallen for x(0) = 0.

    Obtained by integrating ``velocity_expr`` in ``t`` and fixing the constant so x(0) = 0.
    """
    vt = terminal_velocity_expr()
    return vt**2 / g * sp.log(sp.cosh(g * t / vt))


def ode_residual(v_candidate: sp.Expr) -> sp.Expr:
    """Return simplify(d v_candidate/dt - a(v_candidate)); zero if the candidate solves the ODE."""
    return sp.simplify(sp.diff(v_candidate, t) - acceleration_expr().subs(v, v_candidate))


def no_drag_limit() -> sp.Expr:
    """Return lim_{k -> 0} v(t). Must reduce to free fall, g t."""
    return sp.limit(velocity_expr(), k, 0)


def small_time_series(order: int = 4) -> sp.Expr:
    """Return the Taylor series of v(t) about t = 0 up to (not including) t**order, without the O() term."""
    if order < 2:
        raise ValueError("order must be at least 2 to show the leading term g t")
    return sp.series(velocity_expr(), t, 0, order).removeO()


def acceleration_derivative_expr() -> sp.Expr:
    """Return da/dv = -2 k v / m, the sensitivity of acceleration to speed."""
    return sp.diff(acceleration_expr(), v)


def lambdify_acceleration() -> Callable[..., np.ndarray]:
    """Return a NumPy function ``a(v, m, g, k)`` built from ``acceleration_expr``."""
    return sp.lambdify((v, m, g, k), acceleration_expr(), modules="numpy")


def lambdify_velocity() -> Callable[..., np.ndarray]:
    """Return a NumPy function ``v(t, m, g, k)`` built from ``velocity_expr``."""
    return sp.lambdify((t, m, g, k), velocity_expr(), modules="numpy")


def lambdify_position() -> Callable[..., np.ndarray]:
    """Return a NumPy function ``x(t, m, g, k)`` built from ``position_expr``."""
    return sp.lambdify((t, m, g, k), position_expr(), modules="numpy")


def lambdify_acceleration_derivative() -> Callable[..., np.ndarray]:
    """Return a NumPy function ``da_dv(v, m, g, k)`` built from ``acceleration_derivative_expr``."""
    return sp.lambdify((v, m, g, k), acceleration_derivative_expr(), modules="numpy")


def central_difference(f: Callable[[np.ndarray], np.ndarray], x: np.ndarray, h: float) -> np.ndarray:
    """Return the central-difference estimate (f(x+h) - f(x-h)) / (2h) of f'(x). ``h`` must be positive."""
    if h <= 0.0:
        raise ValueError("h must be positive")
    x = np.asarray(x, dtype=float)
    return (f(x + h) - f(x - h)) / (2.0 * h)


def derivative_check(v0: float, m_val: float, g_val: float, k_val: float, h: float) -> float:
    """Return |da/dv (symbolic) - da/dv (central difference)| at speed ``v0``, in 1/s.

    For the quadratic acceleration the central difference is exact up to rounding, so the
    returned number measures floating-point error, not truncation error.
    """
    a = lambdify_acceleration()
    dadv = lambdify_acceleration_derivative()
    fd = central_difference(lambda vv: a(vv, m_val, g_val, k_val), np.asarray(v0, dtype=float), h)
    return float(abs(dadv(v0, m_val, g_val, k_val) - fd))


def velocity_derivative_check(t0: float, m_val: float, g_val: float, k_val: float, h: float) -> float:
    """Return |a(v(t0)) - dv/dt by central difference| in m/s**2: the ODE checked numerically.

    Here the error does have a truncation part proportional to h**2, so halving ``h``
    should divide the error by about four until rounding takes over.
    """
    vel = lambdify_velocity()
    a = lambdify_acceleration()
    fd = central_difference(lambda tt: vel(tt, m_val, g_val, k_val), np.asarray(t0, dtype=float), h)
    return float(abs(a(vel(t0, m_val, g_val, k_val), m_val, g_val, k_val) - fd))


def main(out_dir: str | Path = ".") -> None:
    """Print the derivation and compare the exact solution with a SciPy integration."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.integrate import solve_ivp

    out_dir = Path(out_dir)
    print("equation of motion   dv/dt =", acceleration_expr())
    print("terminal velocity      v_t =", terminal_velocity_expr())
    print("exact solution        v(t) =", velocity_expr())
    print("distance fallen       x(t) =", position_expr())
    print("ODE residual (should be 0) :", ode_residual(velocity_expr()))
    print("dx/dt - v (should be 0)    :", sp.simplify(sp.diff(position_expr(), t) - velocity_expr()))
    print("limit k -> 0               :", no_drag_limit())
    print("series about t = 0         :", small_time_series(6))

    # A 1.2 kg instrument package, C_D = 1.0, frontal area 0.01 m**2, ISA sea-level density.
    m_val, g_val = 1.2, 9.80665
    k_val = drag_constant(1.225, 1.0, 0.01)
    vt = float(terminal_velocity_expr().subs({m: m_val, g: g_val, k: k_val}))
    print(f"\nexample: m = {m_val} kg, k = {k_val:.5f} kg/m  ->  v_t = {vt:.2f} m/s")
    for h in (1e-1, 1e-2, 1e-3, 1e-4):
        print(f"  h = {h:7.0e}: |da/dv symbolic - FD| = {derivative_check(20.0, m_val, g_val, k_val, h):.2e} 1/s,"
              f"  |a - dv/dt(FD)| = {velocity_derivative_check(3.0, m_val, g_val, k_val, h):.2e} m/s^2")

    vel = lambdify_velocity()
    a = lambdify_acceleration()
    ts = np.linspace(0.0, 15.0, 301)
    v_exact = vel(ts, m_val, g_val, k_val)
    sol = solve_ivp(lambda _t, y: [a(y[0], m_val, g_val, k_val)], (0.0, 15.0), [0.0],
                    t_eval=ts, rtol=1e-10, atol=1e-12)
    print(f"max |exact - solve_ivp| = {np.max(np.abs(v_exact - sol.y[0])):.2e} m/s")

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(ts, v_exact, label="exact: $v_t \\tanh(g t / v_t)$")
    ax.plot(ts, sol.y[0], "--", label="scipy solve_ivp")
    ax.plot(ts, g_val * ts, ":", label="no drag: $g t$")
    ax.axhline(vt, color="gray", lw=0.8, label=f"$v_t$ = {vt:.1f} m/s")
    ax.set_ylim(0, 1.3 * vt)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("speed [m/s]")
    ax.set_title("Fall with quadratic drag")
    ax.grid(True)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "fall_with_drag.png", dpi=150)
    print(f"saved {out_dir / 'fall_with_drag.png'}")


if __name__ == "__main__":
    main()
