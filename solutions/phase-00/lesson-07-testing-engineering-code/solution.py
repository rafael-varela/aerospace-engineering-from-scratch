"""P00.07 Testing Engineering Code: reference solution.

A small lift-and-drag coefficient calculator built to be tested. Every function states its
units, checks its inputs, and has a limiting case that a test can pin down. The module also
holds two deliberate demonstrations of floating-point behaviour (``is_close`` and
``rocket_mass_stepwise``) so that the lesson's tests can show why exact equality fails.

Run:  python solution.py          (prints the verification table, writes drag_polar.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-00/lesson-07-testing-engineering-code
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

ArrayLike = float | np.ndarray

# Default parabolic drag polar used throughout the lesson: a clean light aircraft.
CD0_DEFAULT = 0.020   # zero-lift drag coefficient, dimensionless
K_DEFAULT = 0.050     # induced-drag factor k = 1 / (pi * e * AR), dimensionless


def _check_positive(name: str, value: ArrayLike) -> None:
    """Raise ``ValueError`` if ``value`` is not strictly positive and finite (element-wise)."""
    arr = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite, got {value!r}")
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} must be strictly positive, got {value!r}")


def _check_nonnegative(name: str, value: ArrayLike) -> None:
    """Raise ``ValueError`` if ``value`` is negative or not finite (element-wise)."""
    arr = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite, got {value!r}")
    if np.any(arr < 0.0):
        raise ValueError(f"{name} must be non-negative, got {value!r}")


def dynamic_pressure(rho: ArrayLike, v: ArrayLike) -> ArrayLike:
    """Return the dynamic pressure q = rho * v**2 / 2 in Pa.

    Parameters
    ----------
    rho : air density in kg/m**3, strictly positive.
    v : true airspeed in m/s, non-negative (speed, not velocity component).
    """
    _check_positive("rho", rho)
    _check_nonnegative("v", v)
    return 0.5 * np.asarray(rho, dtype=float) * np.asarray(v, dtype=float) ** 2


def lift_coefficient(lift: ArrayLike, rho: ArrayLike, v: ArrayLike, area: ArrayLike) -> ArrayLike:
    """Return C_L = L / (q S), dimensionless.

    Parameters
    ----------
    lift : lift force in N (any sign).
    rho : air density in kg/m**3, strictly positive.
    v : true airspeed in m/s, strictly positive (q = 0 would divide by zero).
    area : reference wing area in m**2, strictly positive.
    """
    _check_positive("v", v)
    _check_positive("area", area)
    q = dynamic_pressure(rho, v)
    return np.asarray(lift, dtype=float) / (q * np.asarray(area, dtype=float))


def drag_coefficient(drag: ArrayLike, rho: ArrayLike, v: ArrayLike, area: ArrayLike) -> ArrayLike:
    """Return C_D = D / (q S), dimensionless. Same input rules as ``lift_coefficient``."""
    _check_nonnegative("drag", drag)
    return lift_coefficient(drag, rho, v, area)


def drag_polar(cl: ArrayLike, cd0: float = CD0_DEFAULT, k: float = K_DEFAULT) -> ArrayLike:
    """Return C_D from the parabolic drag polar C_D = C_D0 + k * C_L**2.

    Parameters
    ----------
    cl : lift coefficient (scalar or array), dimensionless.
    cd0 : zero-lift drag coefficient, non-negative.
    k : induced-drag factor 1 / (pi * e * AR), non-negative.
    """
    _check_nonnegative("cd0", cd0)
    _check_nonnegative("k", k)
    cl_arr = np.asarray(cl, dtype=float)
    if not np.all(np.isfinite(cl_arr)):
        raise ValueError(f"cl must be finite, got {cl!r}")
    return cd0 + k * cl_arr ** 2


def lift_to_drag(cl: ArrayLike, cd0: float = CD0_DEFAULT, k: float = K_DEFAULT) -> ArrayLike:
    """Return L/D = C_L / C_D along the parabolic drag polar."""
    cd = drag_polar(cl, cd0, k)
    if np.any(np.asarray(cd) <= 0.0):
        raise ValueError("drag coefficient must be positive to form L/D")
    return np.asarray(cl, dtype=float) / cd


def max_lift_to_drag(cd0: float = CD0_DEFAULT, k: float = K_DEFAULT) -> tuple[float, float]:
    """Return ``(cl_star, ld_max)``: the C_L of best L/D and that best L/D.

    From d(C_L/C_D)/dC_L = 0 on the parabolic polar: C_L* = sqrt(C_D0 / k) and
    (L/D)_max = 1 / (2 sqrt(C_D0 k)). Both ``cd0`` and ``k`` must be strictly positive.
    """
    _check_positive("cd0", cd0)
    _check_positive("k", k)
    cl_star = math.sqrt(cd0 / k)
    ld_max = 1.0 / (2.0 * math.sqrt(cd0 * k))
    return cl_star, ld_max


def polar_table(cls: np.ndarray, cd0: float = CD0_DEFAULT, k: float = K_DEFAULT) -> np.ndarray:
    """Return an (n, 3) array with columns C_L, C_D, L/D for each lift coefficient in ``cls``.

    This is the table the regression test compares against a stored reference.
    """
    cls = np.asarray(cls, dtype=float)
    cd = drag_polar(cls, cd0, k)
    return np.column_stack([cls, cd, cls / cd])


def is_close(a: float, b: float, rel_tol: float = 1e-9, abs_tol: float = 0.0) -> bool:
    """Return True if ``a`` and ``b`` agree within a relative or an absolute tolerance.

    Implements |a - b| <= max(rel_tol * max(|a|, |b|), abs_tol), the same rule as
    ``math.isclose``. ``rel_tol`` and ``abs_tol`` must be non-negative.
    """
    if rel_tol < 0.0 or abs_tol < 0.0:
        raise ValueError("tolerances must be non-negative")
    if math.isnan(a) or math.isnan(b):
        return False
    if a == b:  # handles equal infinities and exact ties
        return True
    if math.isinf(a) or math.isinf(b):
        return False
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def rocket_mass_closed_form(m0: float, mdot: float, t: ArrayLike) -> ArrayLike:
    """Return m(t) = m0 - mdot * t in kg for constant propellant mass flow ``mdot`` (kg/s)."""
    _check_positive("m0", m0)
    _check_nonnegative("mdot", mdot)
    _check_nonnegative("t", t)
    m = m0 - mdot * np.asarray(t, dtype=float)
    if np.any(np.asarray(m) < 0.0):
        raise ValueError("propellant exhausted: m(t) would be negative")
    return m


def rocket_mass_stepwise(m0: float, mdot: float, dt: float, n_steps: int) -> np.ndarray:
    """Return the mass after each of ``n_steps`` repeated subtractions of ``mdot * dt``.

    Returns an array of length ``n_steps + 1`` starting at ``m0``. This mimics how a
    simulation loop depletes mass, and it accumulates floating-point rounding error that
    the closed form does not have. Raises ``ValueError`` if the mass would go negative.
    """
    _check_positive("m0", m0)
    _check_nonnegative("mdot", mdot)
    _check_positive("dt", dt)
    if n_steps < 0:
        raise ValueError("n_steps must be non-negative")
    masses = np.empty(n_steps + 1, dtype=float)
    m = float(m0)
    masses[0] = m
    for i in range(1, n_steps + 1):
        m = m - mdot * dt
        if m < 0.0:
            raise ValueError(f"propellant exhausted at step {i}")
        masses[i] = m
    return masses


def main(out_dir: str | Path = ".") -> None:
    """Print the verification table and save the drag-polar figure as ``drag_polar.png``."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    print("Floating point: 0.1 + 0.2 == 0.3 ->", 0.1 + 0.2 == 0.3)
    print("                is_close(0.1 + 0.2, 0.3) ->", is_close(0.1 + 0.2, 0.3))
    m_step = float(rocket_mass_stepwise(50.0, 0.1, 1.0, 100)[-1])
    m_exact = float(rocket_mass_closed_form(50.0, 0.1, 100.0))
    print(f"Rocket mass after 100 steps: stepwise {m_step!r}, closed form {m_exact!r}, "
          f"exact equality -> {m_step == m_exact}, difference {m_step - m_exact:.3e} kg")

    cls = np.linspace(0.0, 1.4, 15)
    table = polar_table(cls)
    print("\n   C_L      C_D      L/D")
    for cl, cd, ld in table:
        print(f"{cl:6.3f}  {cd:7.4f}  {ld:7.3f}")
    cl_star, ld_max = max_lift_to_drag()
    print(f"\nbest L/D = {ld_max:.3f} at C_L* = {cl_star:.3f}")
    print(f"table maximum of L/D = {table[:, 2].max():.3f} (should not exceed the analytic value)")

    fig, ax = plt.subplots(figsize=(6, 4))
    fine = np.linspace(0.0, 1.4, 200)
    ax.plot(drag_polar(fine), fine, label=r"$C_D = C_{D0} + k C_L^2$")
    ax.plot(drag_polar(cl_star), cl_star, "o", label="best L/D")
    ax.set_xlabel(r"$C_D$ (dimensionless)")
    ax.set_ylabel(r"$C_L$ (dimensionless)")
    ax.set_title(f"Parabolic drag polar, $C_{{D0}}$={CD0_DEFAULT}, k={K_DEFAULT}")
    ax.grid(True)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "drag_polar.png", dpi=150)
    print(f"saved {out_dir / 'drag_polar.png'}")


if __name__ == "__main__":
    main()
