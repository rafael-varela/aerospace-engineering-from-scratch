"""P00.04 Plotting That Engineers Trust: reference solution.

A small plotting module with a house style. Every function works on Matplotlib's
object model (Figure, Axes, Line2D) so that tests can inspect what was drawn
instead of looking at pixels.

Run:  python solution.py            (writes figures/*.svg)
Test: AEFS_IMPL=solution pytest curriculum/phase-00/lesson-04-plotting-engineers-trust
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless backend: must come before importing pyplot

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.axes import Axes  # noqa: E402
from matplotlib.container import ErrorbarContainer  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

# --------------------------------------------------------------------------- house style

#: One color per *role* a curve can play. Using roles instead of "blue, orange, ..."
#: means every figure in a report reads the same way: measured data is always the
#: same color, the model is always the same color.
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
    guessing whether the unit was forgotten.

    >>> axis_label("Pressure", "Pa")
    'Pressure [Pa]'
    >>> axis_label("Lift coefficient", "-")
    'Lift coefficient [-]'
    """
    quantity = quantity.strip()
    unit = unit.strip() or "-"
    if not quantity:
        raise ValueError("quantity name must not be empty")
    return f"{quantity} [{unit}]"


def apply_house_style() -> None:
    """Set Matplotlib rcParams for the house style (fonts, grid, line widths, vector output)."""
    plt.rcParams.update({
        "figure.figsize": (6.0, 4.0),
        "figure.dpi": 100,
        "font.size": 10,
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 1.6,
        "lines.markersize": 4,
        "legend.frameon": False,
        "savefig.format": "svg",
        "savefig.bbox": "tight",
    })


def new_figure(
    xquantity: str, xunit: str, yquantity: str, yunit: str, *, title: str | None = None
) -> tuple[Figure, Axes]:
    """Create a Figure with one Axes whose x and y labels carry units.

    Returns ``(fig, ax)``. The explicit object-oriented API is used so the caller
    always holds the objects it later saves or inspects.
    """
    apply_house_style()
    fig, ax = plt.subplots()
    ax.set_xlabel(axis_label(xquantity, xunit))
    ax.set_ylabel(axis_label(yquantity, yunit))
    if title:
        ax.set_title(title)
    return fig, ax


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

    ``role`` selects the color from :data:`HOUSE_COLORS`. If ``yerr`` or ``xerr`` is
    given the series is drawn with error bars (``ax.errorbar``) and markers only;
    otherwise as a line (``ax.plot``). Returns the Matplotlib artist container.
    """
    if role not in HOUSE_COLORS:
        raise KeyError(f"unknown role {role!r}; choose one of {sorted(HOUSE_COLORS)}")
    color = HOUSE_COLORS[role]
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.shape != y.shape:
        raise ValueError(f"x and y must have the same shape, got {x.shape} and {y.shape}")
    if yerr is not None or xerr is not None:
        return ax.errorbar(
            x, y, yerr=yerr, xerr=xerr, fmt="o", color=color, ecolor=color,
            capsize=2.5, elinewidth=1.0, label=label,
        )
    (line,) = ax.plot(x, y, color=color, label=label)
    return line


def set_log_axes(ax: Axes, *, x: bool = False, y: bool = False) -> None:
    """Switch the requested axes to a base-10 logarithmic scale."""
    if x:
        ax.set_xscale("log")
    if y:
        ax.set_yscale("log")


def finish_axes(ax: Axes) -> None:
    """Add a legend when more than one labelled series is present."""
    handles, labels = ax.get_legend_handles_labels()
    if len(labels) > 1:
        ax.legend()


def save_vector(fig: Figure, path: str | Path) -> Path:
    """Save ``fig`` as a vector file (.svg or .pdf) and return the path.

    Raises ``ValueError`` for raster suffixes such as .png: a raster figure cannot
    be re-scaled in a report without becoming blurry.
    """
    path = Path(path)
    if path.suffix.lower() not in VECTOR_SUFFIXES:
        raise ValueError(f"house style requires a vector format {VECTOR_SUFFIXES}, got {path.suffix!r}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, format=path.suffix.lstrip(".").lower())
    return path


def check_figure(fig: Figure) -> list[str]:
    """Return a list of house-style violations found in ``fig`` (empty if none).

    Checks, for every Axes: an x label and a y label exist; each label ends in a
    ``[unit]`` tag; a legend exists when more than one labelled series is drawn;
    and no data series is drawn in a color outside :data:`HOUSE_COLORS`.
    """
    problems: list[str] = []
    allowed = {matplotlib.colors.to_hex(c) for c in HOUSE_COLORS.values()}
    for i, ax in enumerate(fig.get_axes()):
        for name, text in (("x", ax.get_xlabel()), ("y", ax.get_ylabel())):
            if not text:
                problems.append(f"axes {i}: missing {name} label")
            elif not (text.rstrip().endswith("]") and "[" in text):
                problems.append(f"axes {i}: {name} label {text!r} has no [unit] tag")
        handles, labels = ax.get_legend_handles_labels()
        if len(labels) > 1 and ax.get_legend() is None:
            problems.append(f"axes {i}: {len(labels)} labelled series but no legend")
        for line in ax.get_lines():
            if matplotlib.colors.to_hex(line.get_color()) not in allowed:
                problems.append(f"axes {i}: line color {line.get_color()!r} is not a house color")
    return problems


# --------------------------------------------------------------------------- example data

def isa_pressure(h_m: np.ndarray) -> np.ndarray:
    """ISA pressure in Pa for geopotential altitude ``h_m`` in m, 0–20 km.

    Troposphere (0–11 km): linear temperature lapse 6.5 K/km. Lower stratosphere
    (11–20 km): isothermal at 216.65 K. Beyond 20 km the model is not valid here.
    """
    h = np.asarray(h_m, dtype=float)
    g0, R, T0, p0, L = 9.80665, 287.05287, 288.15, 101325.0, -0.0065
    h11 = 11_000.0
    T11 = T0 + L * h11
    p11 = p0 * (T11 / T0) ** (-g0 / (L * R))
    p_tropo = p0 * ((T0 + L * np.minimum(h, h11)) / T0) ** (-g0 / (L * R))
    p_strato = p11 * np.exp(-g0 * (h - h11) / (R * T11))
    return np.where(h <= h11, p_tropo, p_strato)


def drag_polar(cl: np.ndarray, *, cd0: float, k: float) -> np.ndarray:
    """Parabolic drag polar ``CD = CD0 + k * CL**2`` (dimensionless)."""
    return cd0 + k * np.asarray(cl, dtype=float) ** 2


def synthetic_thrust_curve(seed: int = 0, n: int = 25) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """A *synthetic* model-rocket thrust curve with simulated noise.

    Returns ``(t_s, thrust_N, sigma_N)``. This is a simulated dataset generated in
    code, not a measurement; its purpose is to practise plotting with error bars.
    """
    rng = np.random.default_rng(seed)
    t = np.linspace(0.0, 2.0, n)
    model = 12.0 * t * np.exp(-t / 0.35)
    sigma = 0.08 * model + 0.15
    return t, model + rng.normal(0.0, sigma), sigma


# --------------------------------------------------------------------------- figures

def figure_atmosphere() -> Figure:
    """Pressure vs altitude, log y, model vs tabulated reference points."""
    fig, ax = new_figure("Altitude", "km", "Pressure", "Pa", title="ISA pressure (0–20 km)")
    h = np.linspace(0.0, 20_000.0, 201)
    plot_series(ax, h / 1000.0, isa_pressure(h), role="model", label="ISA model")
    h_ref = np.array([0.0, 5_000.0, 10_000.0, 15_000.0, 20_000.0])
    p_ref = np.array([101_325.0, 54_020.0, 26_436.0, 12_045.0, 5_475.0])
    plot_series(ax, h_ref / 1000.0, p_ref, role="reference", label="tabulated", yerr=0.005 * p_ref)
    set_log_axes(ax, y=True)
    finish_axes(ax)
    return fig


def figure_drag_polar() -> Figure:
    """Drag polar for two aircraft configurations on linear axes."""
    fig, ax = new_figure("Lift coefficient", "-", "Drag coefficient", "-", title="Parabolic drag polar")
    cl = np.linspace(-0.2, 1.4, 81)
    plot_series(ax, cl, drag_polar(cl, cd0=0.020, k=0.045), role="model", label="clean")
    plot_series(ax, cl, drag_polar(cl, cd0=0.060, k=0.045), role="limit", label="gear down")
    finish_axes(ax)
    return fig


def figure_thrust_curve() -> Figure:
    """Synthetic thrust curve with error bars and the generating model."""
    fig, ax = new_figure("Time", "s", "Thrust", "N", title="Synthetic motor thrust curve")
    t, f, sigma = synthetic_thrust_curve()
    plot_series(ax, t, f, role="measured", label="synthetic samples", yerr=sigma)
    tt = np.linspace(0.0, 2.0, 200)
    plot_series(ax, tt, 12.0 * tt * np.exp(-tt / 0.35), role="model", label="generating model")
    finish_axes(ax)
    return fig


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
