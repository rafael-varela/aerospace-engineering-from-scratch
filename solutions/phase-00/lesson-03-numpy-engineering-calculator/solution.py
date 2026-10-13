"""P00.03 NumPy as an Engineering Calculator: reference solution.

A vectorized International Standard Atmosphere (ISA) from 0 to 86 km geometric altitude:
temperature, pressure, density and speed of sound, with the layer bases derived from the
defining constants rather than typed in from a table.

Run:  python solution.py
Test: AEFS_IMPL=solution pytest curriculum/phase-00/lesson-03-numpy-engineering-calculator
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
# starts at each base. Seven layers, eight bases; the last base is the model top.
H_BASE: np.ndarray = np.array([0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0])
LAPSE: np.ndarray = np.array([-6.5e-3, 0.0, 1.0e-3, 2.8e-3, 0.0, -2.8e-3, -2.0e-3])


def geometric_to_geopotential(z: np.ndarray | float) -> np.ndarray | float:
    """Convert geometric altitude ``z`` (m above mean sea level) to geopotential altitude (m).

    The ISA defines its layers in geopotential altitude, which absorbs the weakening of
    gravity with height so that g can be treated as the constant ``G0`` in the hydrostatic
    equation. For a spherical Earth of radius ``R_EARTH``: h = R_EARTH * z / (R_EARTH + z).
    """
    z = np.asarray(z, dtype=float)
    h = R_EARTH * z / (R_EARTH + z)
    return float(h) if h.ndim == 0 else h


H_TOP: float = float(geometric_to_geopotential(Z_TOP))  # about 84 852 m


def layer_bases() -> tuple[np.ndarray, np.ndarray]:
    """Return (T_base, p_base), the temperature (K) and pressure (Pa) at each layer base.

    Derived by integrating the hydrostatic equation layer by layer from sea level, so the
    numbers follow from the defining constants: nothing is copied from a table.
    """
    n = len(H_BASE)
    t_base = np.empty(n)
    p_base = np.empty(n)
    t_base[0], p_base[0] = T0, P0
    for i in range(n - 1):
        dh = H_BASE[i + 1] - H_BASE[i]
        t_base[i + 1] = t_base[i] + LAPSE[i] * dh
        if LAPSE[i] == 0.0:
            p_base[i + 1] = p_base[i] * np.exp(-G0 * dh / (R_AIR * t_base[i]))
        else:
            p_base[i + 1] = p_base[i] * (t_base[i + 1] / t_base[i]) ** (-G0 / (R_AIR * LAPSE[i]))
    return t_base, p_base


def _prepare(h: np.ndarray | float) -> tuple[np.ndarray, np.ndarray, bool]:
    """Validate ``h`` and return (h_array, layer_index, was_scalar)."""
    arr = np.asarray(h, dtype=float)
    scalar = arr.ndim == 0
    arr = np.atleast_1d(arr)
    if np.any(np.isnan(arr)):
        raise ValueError("altitude contains NaN")
    if np.any(arr < 0.0) or np.any(arr > H_TOP + 1e-6):
        raise ValueError(
            f"geopotential altitude must lie in [0, {H_TOP:.1f}] m (0 to 86 km geometric)"
        )
    idx = np.searchsorted(H_BASE, arr, side="right") - 1
    idx = np.clip(idx, 0, len(H_BASE) - 1)
    return arr, idx, scalar


def _finish(values: np.ndarray, scalar: bool) -> np.ndarray | float:
    return float(values[0]) if scalar else values


def temperature(h: np.ndarray | float) -> np.ndarray | float:
    """ISA temperature (K) at geopotential altitude ``h`` (m), 0 <= h <= H_TOP.

    Piecewise linear: T = T_b + L_b (h - h_b) inside the layer that starts at base b.
    Accepts a scalar or any array shape; returns the same shape (a float for a scalar).
    """
    arr, idx, scalar = _prepare(h)
    t_base, _ = layer_bases()
    t = t_base[idx] + LAPSE[idx] * (arr - H_BASE[idx])
    return _finish(t, scalar)


def pressure(h: np.ndarray | float) -> np.ndarray | float:
    """ISA pressure (Pa) at geopotential altitude ``h`` (m).

    Hydrostatic equilibrium plus the ideal gas law gives, inside one layer,
    p = p_b (T/T_b)^(-g0/(R L))        if L != 0 (gradient layer)
    p = p_b exp(-g0 (h - h_b)/(R T_b))  if L == 0 (isothermal layer)
    Both branches are evaluated on the whole array and combined with ``np.where``.
    """
    arr, idx, scalar = _prepare(h)
    t_base, p_base = layer_bases()
    lapse = LAPSE[idx]
    t = t_base[idx] + lapse * (arr - H_BASE[idx])
    iso = lapse == 0.0
    lapse_safe = np.where(iso, 1.0, lapse)  # avoid 0 in a denominator we will discard
    p_gradient = p_base[idx] * (t / t_base[idx]) ** (-G0 / (R_AIR * lapse_safe))
    p_isothermal = p_base[idx] * np.exp(-G0 * (arr - H_BASE[idx]) / (R_AIR * t_base[idx]))
    p = np.where(iso, p_isothermal, p_gradient)
    return _finish(p, scalar)


def density(h: np.ndarray | float) -> np.ndarray | float:
    """ISA density (kg/m^3) at geopotential altitude ``h`` (m), from the ideal gas law rho = p/(R T)."""
    arr, _, scalar = _prepare(h)
    rho = np.asarray(pressure(arr)) / (R_AIR * np.asarray(temperature(arr)))
    return _finish(np.atleast_1d(rho), scalar)


def speed_of_sound(t: np.ndarray | float) -> np.ndarray | float:
    """Speed of sound (m/s) in air at temperature ``t`` (K): a = sqrt(gamma R T)."""
    arr = np.asarray(t, dtype=float)
    if np.any(arr <= 0.0):
        raise ValueError("temperature must be positive (kelvin)")
    a = np.sqrt(GAMMA * R_AIR * arr)
    return float(a) if a.ndim == 0 else a


def isa(h: np.ndarray | float) -> dict[str, np.ndarray | float]:
    """All four ISA properties at geopotential altitude ``h`` (m).

    Returns a dict with keys ``"T"`` (K), ``"p"`` (Pa), ``"rho"`` (kg/m^3), ``"a"`` (m/s).
    """
    t = temperature(h)
    p = pressure(h)
    rho = np.asarray(p) / (R_AIR * np.asarray(t))
    a = speed_of_sound(t)
    if np.ndim(h) == 0:
        rho = float(rho)
    return {"T": t, "p": p, "rho": rho, "a": a}


def isa_table(z_geometric: np.ndarray) -> np.ndarray:
    """Tabulate the ISA at geometric altitudes ``z_geometric`` (m, 1-D array).

    Returns an array of shape (n, 6) with columns
    [z geometric (m), h geopotential (m), T (K), p (Pa), rho (kg/m^3), a (m/s)].
    """
    z = np.atleast_1d(np.asarray(z_geometric, dtype=float))
    if z.ndim != 1:
        raise ValueError("z_geometric must be one-dimensional")
    h = np.asarray(geometric_to_geopotential(z))
    out = isa(h)
    return np.column_stack([z, h, out["T"], out["p"], out["rho"], out["a"]])


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
