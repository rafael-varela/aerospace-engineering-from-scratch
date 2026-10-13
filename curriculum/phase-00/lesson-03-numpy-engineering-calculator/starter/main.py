"""P00.03 NumPy as an Engineering Calculator: starter code.

Build a vectorized International Standard Atmosphere (ISA) from 0 to 86 km geometric
altitude: temperature, pressure, density and speed of sound. Complete every function marked
NotImplementedError. When finished, running this file prints the ISA table every 2 km and the
sea-level and tropopause verification lines.

Run:  python main.py
Test: pytest curriculum/phase-00/lesson-03-numpy-engineering-calculator   (from the repo root)
"""
from __future__ import annotations

import numpy as np

# --- Defining constants of the ISA (SI units) -------------------------------------------
T0: float = 288.15          # sea-level temperature, K
P0: float = 101325.0        # sea-level pressure, Pa
G0: float = 9.80665         # standard gravitational acceleration, m/s^2
R_AIR: float = 287.05287    # specific gas constant of dry air, J/(kg K)
GAMMA: float = 1.4          # ratio of specific heats for air
R_EARTH: float = 6356766.0  # nominal Earth radius used in the geopotential conversion, m
Z_TOP: float = 86000.0      # top of the model, geometric altitude, m

# Layer bases in geopotential altitude (m) and the lapse rate (K/m) of the layer that
# starts at each base. Seven layers; the model top (about 84 852 m geopotential) is H_TOP.
H_BASE: np.ndarray = np.array([0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0])
LAPSE: np.ndarray = np.array([-6.5e-3, 0.0, 1.0e-3, 2.8e-3, 0.0, -2.8e-3, -2.0e-3])
H_TOP: float = R_EARTH * Z_TOP / (R_EARTH + Z_TOP)  # geopotential altitude of Z_TOP, m


def geometric_to_geopotential(z: np.ndarray | float) -> np.ndarray | float:
    """Convert geometric altitude ``z`` (m above mean sea level) to geopotential altitude (m).

    The ISA defines its layers in geopotential altitude, which absorbs the weakening of
    gravity with height so that g can be treated as the constant ``G0`` in the hydrostatic
    equation. For a spherical Earth of radius ``R_EARTH``: h = R_EARTH * z / (R_EARTH + z).
    Return a float for scalar input and an array of the same shape otherwise.
    """
    raise NotImplementedError


def layer_bases() -> tuple[np.ndarray, np.ndarray]:
    """Return (T_base, p_base), the temperature (K) and pressure (Pa) at each layer base.

    Derive them by integrating the hydrostatic equation layer by layer from sea level, so the
    numbers follow from the defining constants: nothing is copied from a table.
    """
    raise NotImplementedError


def temperature(h: np.ndarray | float) -> np.ndarray | float:
    """ISA temperature (K) at geopotential altitude ``h`` (m), 0 <= h <= H_TOP.

    Piecewise linear: T = T_b + L_b (h - h_b) inside the layer that starts at base b.
    Accepts a scalar or any array shape; returns the same shape (a float for a scalar).
    Raise ValueError if any altitude is below 0 or above H_TOP.
    """
    raise NotImplementedError


def pressure(h: np.ndarray | float) -> np.ndarray | float:
    """ISA pressure (Pa) at geopotential altitude ``h`` (m).

    Hydrostatic equilibrium plus the ideal gas law gives, inside one layer,
    p = p_b (T/T_b)^(-g0/(R L))        if L != 0 (gradient layer)
    p = p_b exp(-g0 (h - h_b)/(R T_b))  if L == 0 (isothermal layer)
    Evaluate both branches on the whole array and combine them with ``np.where``.
    Same shape and range rules as ``temperature``.
    """
    raise NotImplementedError


def density(h: np.ndarray | float) -> np.ndarray | float:
    """ISA density (kg/m^3) at geopotential altitude ``h`` (m), from the ideal gas law rho = p/(R T)."""
    raise NotImplementedError


def speed_of_sound(t: np.ndarray | float) -> np.ndarray | float:
    """Speed of sound (m/s) in air at temperature ``t`` (K): a = sqrt(gamma R T).

    Raise ValueError if any temperature is not positive.
    """
    raise NotImplementedError


def isa(h: np.ndarray | float) -> dict[str, np.ndarray | float]:
    """All four ISA properties at geopotential altitude ``h`` (m).

    Returns a dict with keys ``"T"`` (K), ``"p"`` (Pa), ``"rho"`` (kg/m^3), ``"a"`` (m/s).
    """
    raise NotImplementedError


def isa_table(z_geometric: np.ndarray) -> np.ndarray:
    """Tabulate the ISA at geometric altitudes ``z_geometric`` (m, 1-D array).

    Returns an array of shape (n, 6) with columns
    [z geometric (m), h geopotential (m), T (K), p (Pa), rho (kg/m^3), a (m/s)].
    Raise ValueError if the input is not one-dimensional.
    """
    raise NotImplementedError


def main() -> None:
    """Print the ISA table every 2 km from 0 to 86 km and the sea-level / tropopause checks."""
    z = np.arange(0.0, 86000.0 + 1.0, 2000.0)
    table = isa_table(z)
    print(f"{'z km':>7} {'h km':>8} {'T K':>8} {'p Pa':>12} {'rho kg/m3':>12} {'a m/s':>8}")
    for row in table:
        print(f"{row[0]/1000:7.1f} {row[1]/1000:8.3f} {row[2]:8.2f} {row[3]:12.3f} {row[4]:12.6f} {row[5]:8.2f}")
    sl = isa(0.0)
    tp = isa(11000.0)
    print()
    print("Verification against the defining values")
    print(f"  sea level : T={sl['T']:.2f} K (288.15)  p={sl['p']:.0f} Pa (101325)  "
          f"rho={sl['rho']:.4f} kg/m^3 (1.2250)  a={sl['a']:.2f} m/s (340.29)")
    print(f"  11 km     : T={tp['T']:.2f} K (216.65)  p={tp['p']:.1f} Pa (22632)  "
          f"rho={tp['rho']:.4f} kg/m^3 (0.3639)  a={tp['a']:.2f} m/s (295.07)")


if __name__ == "__main__":
    main()
