"""P00.07 Testing Engineering Code: starter code.

A small lift-and-drag coefficient calculator built to be tested. Complete every function
marked NotImplementedError. Every function states its units, checks its inputs, and has a
limiting case that a test can pin down. Two functions (``is_close`` and
``rocket_mass_stepwise``) exist so the tests can show why exact float equality fails.

Run:  python main.py            (prints the verification table, writes drag_polar.png)
Test: pytest curriculum/phase-00/lesson-07-testing-engineering-code
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

ArrayLike = float | np.ndarray

# Default parabolic drag polar used throughout the lesson: a clean light aircraft.
CD0_DEFAULT = 0.020   # zero-lift drag coefficient, dimensionless
K_DEFAULT = 0.050     # induced-drag factor k = 1 / (pi * e * AR), dimensionless


def dynamic_pressure(rho: ArrayLike, v: ArrayLike) -> ArrayLike:
    """Return the dynamic pressure q = rho * v**2 / 2 in Pa.

    Parameters
    ----------
    rho : air density in kg/m**3, strictly positive (raise ``ValueError`` otherwise).
    v : true airspeed in m/s, non-negative (raise ``ValueError`` otherwise).
    """
    raise NotImplementedError


def lift_coefficient(lift: ArrayLike, rho: ArrayLike, v: ArrayLike, area: ArrayLike) -> ArrayLike:
    """Return C_L = L / (q S), dimensionless.

    Parameters
    ----------
    lift : lift force in N (any sign).
    rho : air density in kg/m**3, strictly positive.
    v : true airspeed in m/s, strictly positive (q = 0 would divide by zero).
    area : reference wing area in m**2, strictly positive.

    Raise ``ValueError`` on any non-positive or non-finite ``rho``, ``v`` or ``area``.
    """
    raise NotImplementedError


def drag_coefficient(drag: ArrayLike, rho: ArrayLike, v: ArrayLike, area: ArrayLike) -> ArrayLike:
    """Return C_D = D / (q S), dimensionless. Same input rules as ``lift_coefficient``;
    additionally ``drag`` must be non-negative."""
    raise NotImplementedError


def drag_polar(cl: ArrayLike, cd0: float = CD0_DEFAULT, k: float = K_DEFAULT) -> ArrayLike:
    """Return C_D from the parabolic drag polar C_D = C_D0 + k * C_L**2.

    Parameters
    ----------
    cl : lift coefficient (scalar or array), dimensionless, finite.
    cd0 : zero-lift drag coefficient, non-negative.
    k : induced-drag factor 1 / (pi * e * AR), non-negative.

    Raise ``ValueError`` for negative ``cd0`` or ``k`` or a non-finite ``cl``.
    """
    raise NotImplementedError


def lift_to_drag(cl: ArrayLike, cd0: float = CD0_DEFAULT, k: float = K_DEFAULT) -> ArrayLike:
    """Return L/D = C_L / C_D along the parabolic drag polar."""
    raise NotImplementedError


def max_lift_to_drag(cd0: float = CD0_DEFAULT, k: float = K_DEFAULT) -> tuple[float, float]:
    """Return ``(cl_star, ld_max)``: the C_L of best L/D and that best L/D.

    From d(C_L/C_D)/dC_L = 0 on the parabolic polar: C_L* = sqrt(C_D0 / k) and
    (L/D)_max = 1 / (2 sqrt(C_D0 k)). Both ``cd0`` and ``k`` must be strictly positive.
    """
    raise NotImplementedError


def polar_table(cls: np.ndarray, cd0: float = CD0_DEFAULT, k: float = K_DEFAULT) -> np.ndarray:
    """Return an (n, 3) array with columns C_L, C_D, L/D for each lift coefficient in ``cls``.

    This is the table the regression test compares against a stored reference.
    """
    raise NotImplementedError


def is_close(a: float, b: float, rel_tol: float = 1e-9, abs_tol: float = 0.0) -> bool:
    """Return True if ``a`` and ``b`` agree within a relative or an absolute tolerance.

    Implement |a - b| <= max(rel_tol * max(|a|, |b|), abs_tol), the same rule as
    ``math.isclose``. NaN is never close to anything; equal infinities are close.
    ``rel_tol`` and ``abs_tol`` must be non-negative (raise ``ValueError`` otherwise).
    """
    raise NotImplementedError


def rocket_mass_closed_form(m0: float, mdot: float, t: ArrayLike) -> ArrayLike:
    """Return m(t) = m0 - mdot * t in kg for constant propellant mass flow ``mdot`` (kg/s).

    Raise ``ValueError`` if ``m0 <= 0``, ``mdot < 0``, ``t < 0`` or the mass would go negative.
    """
    raise NotImplementedError


def rocket_mass_stepwise(m0: float, mdot: float, dt: float, n_steps: int) -> np.ndarray:
    """Return the mass after each of ``n_steps`` repeated subtractions of ``mdot * dt``.

    Returns an array of length ``n_steps + 1`` starting at ``m0``. This mimics how a
    simulation loop depletes mass, and it accumulates floating-point rounding error that
    the closed form does not have. Raises ``ValueError`` if the mass would go negative.
    """
    raise NotImplementedError


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
