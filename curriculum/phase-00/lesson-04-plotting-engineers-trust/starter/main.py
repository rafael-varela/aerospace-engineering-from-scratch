"""P00.04 Plotting That Engineers Trust: starter code.

Complete every function marked NotImplementedError. When finished, running this file
builds the three lesson figures, checks them against the house style, and saves them
as SVG files in ``figures/``.

Run:  python main.py
Test: pytest tests/
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless backend: must come before importing pyplot

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.axes import Axes  # noqa: E402
from matplotlib.container import ErrorbarContainer  # noqa: E402,F401
from matplotlib.figure import Figure  # noqa: E402

# --------------------------------------------------------------------------- house style

#: One color per *role* a curve can play. Using roles instead of "blue, orange, ..."
#: means every figure in a report reads the same way: measured data is always the
#: same color, the model is always the same color. Keep these keys; you may change
#: the hex values.
HOUSE_COLORS: dict[str, str] = {
    "measured": "#1f5fa8",   # blue: things you measured
    "model": "#d1495b",      # red: things you computed from a model
    "reference": "#4d4d4d",  # grey: published or tabulated values
    "fit": "#2a9d8f",        # teal: a curve fitted to the measured data
    "limit": "#8a6d3b",      # brown: a limit or envelope line
}

#: File suffixes the house style allows for figures that go into a report.
VECTOR_SUFFIXES: tuple[str, ...] = (".svg", ".pdf")


def axis_label(quantity: str, unit: str) -> str:
    """Return a house-style axis label: ``"Quantity [unit]"``.

    Dimensionless quantities use ``"-"`` as the unit, so a reader is never left
    guessing whether the unit was forgotten. An empty ``unit`` means dimensionless.
    Raise ``ValueError`` if ``quantity`` is empty.

    >>> axis_label("Pressure", "Pa")
    'Pressure [Pa]'
    >>> axis_label("Lift coefficient", "-")
    'Lift coefficient [-]'
    """
    raise NotImplementedError


def apply_house_style() -> None:
    """Set Matplotlib rcParams for the house style (fonts, grid, line widths, vector output).

    Use ``plt.rcParams.update({...})``. At minimum set a figure size, a font size,
    ``axes.grid``, and ``savefig.format`` to ``"svg"``.
    """
    raise NotImplementedError


def new_figure(
    xquantity: str, xunit: str, yquantity: str, yunit: str, *, title: str | None = None
) -> tuple[Figure, Axes]:
    """Create a Figure with one Axes whose x and y labels carry units.

    Returns ``(fig, ax)``. Use the explicit object-oriented API (``plt.subplots``)
    so the caller always holds the objects it later saves or inspects.
    """
    raise NotImplementedError


def plot_series(
    ax: Axes,
    x: np.ndarray,
    y: np.ndarray,
    *,
    role: str,
    label: str,
    yerr: np.ndarray | float | None = None,
    xerr: np.ndarray | float | None = None,
):
    """Draw one data series in its house color.

    ``role`` selects the color from :data:`HOUSE_COLORS`; raise ``KeyError`` for an
    unknown role and ``ValueError`` if ``x`` and ``y`` differ in shape. If ``yerr``
    or ``xerr`` is given the series is drawn with error bars (``ax.errorbar``,
    markers only); otherwise as a line (``ax.plot``). Return the Matplotlib artist:
    the ``Line2D`` from ``ax.plot`` or the ``ErrorbarContainer`` from ``ax.errorbar``.
    """
    raise NotImplementedError


def set_log_axes(ax: Axes, *, x: bool = False, y: bool = False) -> None:
    """Switch the requested axes to a base-10 logarithmic scale."""
    raise NotImplementedError


def finish_axes(ax: Axes) -> None:
    """Add a legend when more than one labelled series is present."""
    raise NotImplementedError


def save_vector(fig: Figure, path: str | Path) -> Path:
    """Save ``fig`` as a vector file (.svg or .pdf) and return the path.

    Raise ``ValueError`` for raster suffixes such as .png: a raster figure cannot
    be re-scaled in a report without becoming blurry. Create parent folders as needed.
    """
    raise NotImplementedError


def check_figure(fig: Figure) -> list[str]:
    """Return a list of house-style violations found in ``fig`` (empty if none).

    Check, for every Axes: an x label and a y label exist (message contains
    "x label" / "y label"); each label ends in a ``[unit]`` tag (message contains
    "unit"); a legend exists when more than one labelled series is drawn (message
    contains "legend"); and no line is drawn in a color outside :data:`HOUSE_COLORS`
    (message contains "color"). Compare colors with ``matplotlib.colors.to_hex``.
    """
    raise NotImplementedError


# --------------------------------------------------------------------------- example data

def isa_pressure(h_m: np.ndarray) -> np.ndarray:
    """ISA pressure in Pa for geopotential altitude ``h_m`` in m, 0–20 km.

    Troposphere (0–11 km): linear temperature lapse 6.5 K/km. Lower stratosphere
    (11–20 km): isothermal at 216.65 K. Beyond 20 km the model is not valid here.
    Constants: g0 = 9.80665 m/s^2, R = 287.05287 J/(kg K), T0 = 288.15 K, p0 = 101325 Pa.
    """
    raise NotImplementedError


def drag_polar(cl: np.ndarray, *, cd0: float, k: float) -> np.ndarray:
    """Parabolic drag polar ``CD = CD0 + k * CL**2`` (dimensionless)."""
    raise NotImplementedError


def synthetic_thrust_curve(seed: int = 0, n: int = 25) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """A *synthetic* model-rocket thrust curve with simulated noise.

    Returns ``(t_s, thrust_N, sigma_N)`` with ``t`` from 0 to 2 s, the model
    ``12 t exp(-t/0.35)`` plus Gaussian noise of standard deviation ``sigma``, and
    ``sigma = 0.08 * model + 0.15`` (always positive). Use
    ``np.random.default_rng(seed)`` so the same seed gives the same samples.
    This is a simulated dataset generated in code, not a measurement.
    """
    raise NotImplementedError


# --------------------------------------------------------------------------- figures

def figure_atmosphere() -> Figure:
    """Pressure vs altitude in km, log y, ISA model line plus tabulated reference
    points with error bars, legend. Return the Figure."""
    raise NotImplementedError


def figure_drag_polar() -> Figure:
    """Drag polar for two aircraft configurations (clean, gear down) on linear axes,
    dimensionless labels, legend. Return the Figure."""
    raise NotImplementedError


def figure_thrust_curve() -> Figure:
    """Synthetic thrust samples with error bars plus the generating model line,
    legend. Return the Figure."""
    raise NotImplementedError


def main(out_dir: str | Path = "figures") -> list[Path]:
    """Build the three lesson figures, lint them, and save them as SVG. Returns the paths."""
    out = Path(out_dir)
    written: list[Path] = []
    for name, builder in (
        ("isa_pressure", figure_atmosphere),
        ("drag_polar", figure_drag_polar),
        ("thrust_curve", figure_thrust_curve),
    ):
        fig = builder()
        problems = check_figure(fig)
        status = "ok" if not problems else "; ".join(problems)
        print(f"{name:<14s} style check: {status}")
        written.append(save_vector(fig, out / f"{name}.svg"))
        plt.close(fig)
    for p in written:
        print("wrote", p)
    return written


if __name__ == "__main__":
    main()
