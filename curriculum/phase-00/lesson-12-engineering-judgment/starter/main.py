"""P00.12 Engineering Judgment & Sanity Checks: starter code.

An order-of-magnitude checker for an airliner in cruise. Given wing area, mass, cruise speed and
altitude it computes the lift coefficient, wing loading, Mach number and the thrust-to-weight ratio
needed for level flight, and flags any value outside the range an airliner can plausibly have,
with the reason. The ISA relations are written here in full (not imported) so the file stands alone.

Complete every function marked NotImplementedError. When finished, running this file prints the
checker's verdict on three airliner cases.

Run:  python main.py
Test: pytest curriculum/phase-00/lesson-12-engineering-judgment
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

# --- International Standard Atmosphere, troposphere and lower stratosphere (0 to 20 km) -----------
G0 = 9.80665            # m/s², standard gravity
R_AIR = 287.05287       # J/(kg·K), specific gas constant of air
GAMMA = 1.4             # ratio of specific heats
T0 = 288.15             # K, sea-level temperature
P0 = 101325.0           # Pa, sea-level pressure
LAPSE = -0.0065         # K/m, tropospheric temperature gradient
H_TROPOPAUSE = 11000.0  # m
H_MAX = 20000.0         # m, top of the isothermal layer modelled here
T_STRAT = T0 + LAPSE * H_TROPOPAUSE  # 216.65 K


def isa(altitude_m: float) -> tuple[float, float, float]:
    """Return ``(temperature_K, pressure_Pa, density_kg_m3)`` of the ISA at a geopotential altitude.

    Troposphere (0–11 km): ``T = T0 + L h``, ``p = p0 (T/T0)^(-g0/(L R))``.
    Lower stratosphere (11–20 km): ``T = 216.65 K``, ``p = p11 exp(-g0 (h - 11000)/(R T))``.
    Density from the ideal-gas law ``rho = p/(R T)``. Raises ``ValueError`` outside -500 m..20 km.
    """
    raise NotImplementedError


def speed_of_sound(temperature_k: float) -> float:
    """``a = sqrt(gamma R T)`` in m/s."""
    raise NotImplementedError


def dynamic_pressure(speed_m_s: float, altitude_m: float) -> float:
    """``q = 1/2 rho V²`` in Pa, with ISA density at ``altitude_m``."""
    raise NotImplementedError


def wing_loading(mass_kg: float, wing_area_m2: float) -> float:
    """Weight per unit wing area ``W/S`` in N/m²."""
    raise NotImplementedError


def lift_coefficient(mass_kg: float, wing_area_m2: float, speed_m_s: float, altitude_m: float) -> float:
    """Lift coefficient for steady level flight, ``C_L = W / (q S)``."""
    raise NotImplementedError


def mach_number(speed_m_s: float, altitude_m: float) -> float:
    """Flight Mach number ``V/a`` with ISA temperature at ``altitude_m``."""
    raise NotImplementedError


def lift_to_drag(cl: float, cd0: float = 0.02, oswald: float = 0.8, aspect_ratio: float = 9.0) -> float:
    """``L/D`` from the parabolic drag polar ``C_D = C_D0 + C_L² / (pi e AR)``."""
    raise NotImplementedError


def thrust_to_weight_required(
    mass_kg: float,
    wing_area_m2: float,
    speed_m_s: float,
    altitude_m: float,
    cd0: float = 0.02,
    oswald: float = 0.8,
    aspect_ratio: float = 9.0,
) -> float:
    """Thrust-to-weight ratio for steady level flight: ``T/W = D/L = 1/(L/D)``."""
    raise NotImplementedError


def order_of_magnitude(x: float) -> int:
    """Nearest power of ten: ``order_of_magnitude(950) == 3``, ``order_of_magnitude(0.002) == -3``."""
    raise NotImplementedError


# --- plausibility ranges for a subsonic jet airliner in cruise ---------------------------------------
# (quantity: (low, high, unit, reason if too low, reason if too high))
PLAUSIBLE: dict[str, tuple[float, float, str, str, str]] = {
    "mass": (
        20_000.0, 600_000.0, "kg",
        "lighter than any jet airliner; a regional turboprop or business jet, or mass given in tonnes",
        "heavier than any aircraft ever certified; check for kg given as lb or as grams",
    ),
    "altitude": (
        6_000.0, 14_000.0, "m",
        "below normal jet cruise altitude; fuel burn would be far above design",
        "above the ceiling of subsonic airliners; the cabin and engines are not designed for it",
    ),
    "mach": (
        0.60, 0.92, "",
        "slower than jet cruise; is the speed in m/s and not km/h or knots?",
        "at or beyond the transonic drag rise; a conventional airliner cannot cruise here",
    ),
    "wing_loading": (
        3_000.0, 8_000.0, "N/m²",
        "wing far too large for the mass; glider or light-aircraft territory",
        "wing far too small for the mass; approach speeds would be unflyable",
    ),
    "lift_coefficient": (
        0.30, 0.80, "",
        "very low C_L: flying far too fast for the wing or the wing is far too big; drag would be wasteful",
        "C_L above the cruise range and approaching a clean wing's maximum; too slow, too heavy, or wing too small",
    ),
    "thrust_to_weight": (
        0.04, 0.12, "",
        "L/D above 25 is beyond any airliner; check inputs",
        "more than 12 % of weight as thrust in cruise means L/D below 8; the aircraft would be a fuel sieve",
    ),
}


@dataclass(frozen=True)
class Flag:
    """One checked quantity: its value, the plausible range, and the reason if it falls outside."""

    quantity: str
    value: float
    unit: str
    low: float
    high: float
    ok: bool
    reason: str

    def __str__(self) -> str:
        status = "ok  " if self.ok else "FLAG"
        return f"{status} {self.quantity:<17} {self.value:10.4g} {self.unit:<5} (plausible {self.low:g}–{self.high:g})" + (
            "" if self.ok else f"  <- {self.reason}"
        )


@dataclass
class SanityReport:
    """The result of ``check_airliner``: every quantity checked, in order."""

    flags: list[Flag] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return all(f.ok for f in self.flags)

    @property
    def problems(self) -> list[Flag]:
        return [f for f in self.flags if not f.ok]

    def summary(self) -> str:
        head = "PASS: every quantity is inside the airliner-class range." if self.ok else f"FLAGGED: {len(self.problems)} quantity(ies) outside plausible range."
        return "\n".join([head, *map(str, self.flags)])


def _flag(quantity: str, value: float) -> Flag:
    raise NotImplementedError


def check_airliner(
    mass_kg: float,
    wing_area_m2: float,
    speed_m_s: float,
    altitude_m: float,
    cd0: float = 0.02,
    oswald: float = 0.8,
    aspect_ratio: float = 9.0,
) -> SanityReport:
    """Compute cruise quantities for an airliner and flag each one that is outside ``PLAUSIBLE``.

    Checks, in order: mass, altitude, Mach number, wing loading, lift coefficient, thrust-to-weight
    required. Raises ``ValueError`` for non-positive inputs or an altitude outside the ISA range
    modelled here, because those are not implausible values but impossible ones.
    """
    raise NotImplementedError


def main() -> None:
    """Run the checker on an A320-class airliner and on the same aircraft with a units mistake."""
    cases = {
        "A320-class, correct inputs": dict(mass_kg=70_000, wing_area_m2=122.6, speed_m_s=230.0, altitude_m=11_000),
        "same aircraft, mass typed in lb": dict(mass_kg=154_000, wing_area_m2=122.6, speed_m_s=230.0, altitude_m=11_000),
        "same aircraft, speed typed in km/h": dict(mass_kg=70_000, wing_area_m2=122.6, speed_m_s=830.0, altitude_m=11_000),
    }
    for name, kw in cases.items():
        print(f"\n== {name} ==")
        print(check_airliner(**kw).summary())


if __name__ == "__main__":
    main()
