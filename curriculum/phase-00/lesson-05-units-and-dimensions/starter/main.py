"""P00.05 Units & Dimensional Consistency: starter code.

Complete every function marked NotImplementedError. When finished, running this file
prints a unit-check report and a re-enactment of the Mars Climate Orbiter (MCO) unit
mix-up: an impulse computed in lbf·s and read as N·s.

Run:  python main.py
Test: pytest tests/
"""
from __future__ import annotations

import numpy as np
import pint
from pint import DimensionalityError, Quantity, UnitRegistry

#: One registry for the whole module. Quantities from different registries cannot be
#: mixed, so every function here (and every test) uses this one.
ureg: UnitRegistry = UnitRegistry()
Q_ = ureg.Quantity

#: Exact definition of the pound-force in newtons (SI): 0.45359237 kg × 9.80665 m/s².
LBF_IN_N: float = 0.45359237 * 9.80665


# --------------------------------------------------------------------------- basics

def parse_quantity(text: str) -> Quantity:
    """Parse a string such as ``"35 000 ft"`` or ``"250 knot"`` into a pint Quantity."""
    raise NotImplementedError


def to_si(q: Quantity) -> Quantity:
    """Convert ``q`` to SI base units (m, kg, s, ...). See ``Quantity.to_base_units``."""
    raise NotImplementedError


def magnitude_in(q: Quantity, unit: str) -> float:
    """Return the plain number of ``q`` expressed in ``unit``.

    This is the only sanctioned way to hand a value to NumPy or Matplotlib: the unit is
    named at the boundary. Raises ``DimensionalityError`` if ``unit`` has the wrong
    dimension (so you cannot ask for a pressure "in metres").
    """
    raise NotImplementedError


def round_trip(value: float, unit_a: str, unit_b: str) -> float:
    """Convert ``value unit_a`` to ``unit_b`` and back; return the magnitude in ``unit_a``.

    A correct conversion table round-trips to within floating-point error.
    """
    raise NotImplementedError


def add_checked(a: Quantity, b: Quantity) -> Quantity:
    """Add two quantities. pint raises ``DimensionalityError`` if their dimensions differ.

    The function exists to make that behaviour explicit and testable: 1 m + 1 s is an
    error, 1 m + 1 ft is 1.3048 m.
    """
    raise NotImplementedError


def has_dimension(q: Quantity, dimension: str) -> bool:
    """True if ``q`` has the given dimensionality, written as e.g. ``"[length] / [time]"``.

    See ``Quantity.check``.
    """
    raise NotImplementedError


# --------------------------------------------------------------------------- aerospace formulas

def dynamic_pressure(rho: Quantity, v: Quantity) -> Quantity:
    """Dynamic pressure ``q = 1/2 rho v^2`` in pascal.

    Raises ``DimensionalityError`` if ``rho`` is not a density or ``v`` not a speed:
    check with ``rho.check("[mass] / [length] ** 3")`` and raise
    ``DimensionalityError(rho.units, "kg / m ** 3")`` yourself; pint would otherwise
    happily multiply a speed by a density squared and return a nonsense unit.
    """
    raise NotImplementedError


def total_impulse(thrust: Quantity, burn_time: Quantity) -> Quantity:
    """Total impulse ``I = F * t`` in newton-second, for any force and time units.

    Raise ``DimensionalityError`` if ``thrust`` is not a force or ``burn_time`` not a time.
    """
    raise NotImplementedError


def delta_v(impulse: Quantity, mass: Quantity) -> Quantity:
    """Velocity change ``dv = I / m`` in m/s for a short burn (mass treated as constant).

    Raise ``DimensionalityError`` if ``impulse`` is not ``[force] * [time]`` or ``mass``
    is not a mass.
    """
    raise NotImplementedError


# --------------------------------------------------------------------------- Mars Climate Orbiter

def lbf_s_to_N_s_factor() -> float:
    """The number of newton-seconds in one pound-force-second (about 4.448).

    This is the factor by which MCO's small-forces data were off: the ground software
    wrote impulses in lbf·s, the navigation software read them as N·s, so every
    thruster impulse was under-estimated by this factor.
    """
    raise NotImplementedError


def mco_misread_impulse(impulse_lbf_s: float) -> tuple[Quantity, Quantity]:
    """Re-enact the MCO interface error for one thruster firing.

    ``impulse_lbf_s`` is the number the ground software wrote to the file (its true
    unit was lbf·s). Returns ``(true_impulse, misread_impulse)``: the value the file
    meant, and the value the trajectory software used after attaching N·s to the bare
    number. Both are returned in N·s so they can be compared directly.
    """
    raise NotImplementedError


def mco_underestimate_ratio(impulse_lbf_s: float) -> float:
    """Ratio true/misread for a firing; equals ``lbf_s_to_N_s_factor()`` for any size."""
    raise NotImplementedError


def accumulated_delta_v_error(
    impulses_lbf_s: np.ndarray, spacecraft_mass: Quantity
) -> Quantity:
    """Total velocity change *missed* by the navigation model over many firings.

    Each firing's true impulse is ``x lbf·s``; the model booked ``x N·s``. The missed
    impulse is the difference, summed over all firings, divided by the mass. Return
    a Quantity in m/s (reuse :func:`delta_v`).
    """
    raise NotImplementedError


# --------------------------------------------------------------------------- report

def main() -> None:
    """Print a short unit-check report and the MCO re-enactment table."""
    print("Round trips (should be exact to ~1e-12):")
    for v, a, b in ((35_000.0, "ft", "m"), (250.0, "knot", "m/s"), (14.7, "psi", "Pa")):
        rt = round_trip(v, a, b)
        print(f"  {v:>10g} {a:<5s} -> {b:<4s} -> {a:<5s} = {rt:.12g}   (error {rt - v:+.2e})")

    print("\nA cruise condition, unit-checked:")
    rho = Q_(0.3796, "kg / m**3")        # ISA density near 11 km
    v = parse_quantity("250 knot")
    q = dynamic_pressure(rho, v)
    print(f"  rho = {rho:~P}, V = {v:~P} = {v.to('m/s'):~.1fP}")
    print(f"  q   = {q:~.0fP} = {q.to('psi'):~.3fP}")

    print("\nA mismatched addition is rejected:")
    try:
        add_checked(Q_(1.0, "m"), Q_(1.0, "s"))
    except DimensionalityError as err:
        print(f"  DimensionalityError: {err}")

    print("\nMars Climate Orbiter small-forces re-enactment (illustrative numbers):")
    print(f"  1 lbf*s = {lbf_s_to_N_s_factor():.6f} N*s")
    true, misread = mco_misread_impulse(1.0)
    print(f"  file value 1.000 (meant {true:~.3fP}, read as {misread:~.3fP}); "
          f"ratio {mco_underestimate_ratio(1.0):.3f}")
    rng = np.random.default_rng(0)
    firings = rng.uniform(0.2, 0.6, size=400)  # illustrative lbf*s per firing, not flight data
    mass = Q_(600.0, "kg")                     # illustrative spacecraft mass
    dv_missed = accumulated_delta_v_error(firings, mass)
    print(f"  {firings.size} illustrative firings, {mass:~P}: missed dv = {dv_missed:~.3fP}")


if __name__ == "__main__":
    main()
