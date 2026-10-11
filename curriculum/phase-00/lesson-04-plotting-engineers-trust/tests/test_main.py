"""Public tests for P00.04 Plotting That Engineers Trust.

They inspect figures through Matplotlib's object model (labels, scales, lines,
error-bar containers, colors), never pixels. Run `pytest <lesson dir>` against your
starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>` against the reference
solution. Hidden validation cases live in assessments/phase-00/parameterized-cases/.
"""
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402
from matplotlib.container import ErrorbarContainer  # noqa: E402

from conftest import load_impl  # noqa: E402

impl = load_impl(__file__)


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


# ------------------------------------------------------------------ labels and units

def test_axis_label_has_unit_in_brackets():
    assert impl.axis_label("Pressure", "Pa") == "Pressure [Pa]"
    assert impl.axis_label("Lift coefficient", "-") == "Lift coefficient [-]"
    # a blank unit is never silently dropped: dimensionless is written explicitly
    assert impl.axis_label("Mach number", "") == "Mach number [-]"


def test_new_figure_labels_carry_units():
    fig, ax = impl.new_figure("Altitude", "km", "Pressure", "Pa")
    assert "[km]" in ax.get_xlabel()
    assert "[Pa]" in ax.get_ylabel()
    assert ax.get_xlabel().startswith("Altitude")


def test_backend_is_agg():
    assert matplotlib.get_backend().lower() == "agg"


# ------------------------------------------------------------------ series and colors

def test_plot_series_uses_house_color_and_counts_lines():
    fig, ax = impl.new_figure("Time", "s", "Thrust", "N")
    t = np.linspace(0.0, 1.0, 11)
    impl.plot_series(ax, t, t**2, role="measured", label="a")
    impl.plot_series(ax, t, 2 * t, role="model", label="b")
    lines = ax.get_lines()
    assert len(lines) == 2
    colors = {matplotlib.colors.to_hex(ln.get_color()) for ln in lines}
    expected = {matplotlib.colors.to_hex(impl.HOUSE_COLORS["measured"]),
                matplotlib.colors.to_hex(impl.HOUSE_COLORS["model"])}
    assert colors == expected


def test_plot_series_data_is_what_was_passed():
    fig, ax = impl.new_figure("x", "m", "y", "m")
    x = np.array([0.0, 1.0, 2.0])
    y = np.array([1.0, 3.0, 5.0])
    line = impl.plot_series(ax, x, y, role="model", label="line")
    xd, yd = line.get_data()
    assert np.allclose(xd, x) and np.allclose(yd, y)


def test_plot_series_rejects_unknown_role_and_shape_mismatch():
    fig, ax = impl.new_figure("x", "m", "y", "m")
    with pytest.raises(KeyError):
        impl.plot_series(ax, [0, 1], [0, 1], role="purple-ish", label="bad")
    with pytest.raises(ValueError):
        impl.plot_series(ax, [0, 1, 2], [0, 1], role="model", label="bad")


def test_error_bars_create_an_errorbar_container():
    fig, ax = impl.new_figure("Time", "s", "Thrust", "N")
    t, f, sigma = impl.synthetic_thrust_curve(seed=1, n=10)
    impl.plot_series(ax, t, f, role="measured", label="data", yerr=sigma)
    containers = [c for c in ax.containers if isinstance(c, ErrorbarContainer)]
    assert len(containers) == 1
    assert containers[0].has_yerr and not containers[0].has_xerr


# ------------------------------------------------------------------ log scales

def test_set_log_axes_sets_only_requested_axis():
    fig, ax = impl.new_figure("Altitude", "km", "Pressure", "Pa")
    impl.set_log_axes(ax, y=True)
    assert ax.get_yscale() == "log"
    assert ax.get_xscale() == "linear"


def test_atmosphere_figure_is_log_pressure_vs_linear_altitude():
    fig = impl.figure_atmosphere()
    (ax,) = fig.get_axes()
    assert ax.get_yscale() == "log"
    assert ax.get_xscale() == "linear"
    assert "[Pa]" in ax.get_ylabel() and "[km]" in ax.get_xlabel()
    # one model line plus a reference series with error bars
    assert len(ax.get_lines()) >= 1
    assert any(isinstance(c, ErrorbarContainer) for c in ax.containers)
    assert ax.get_legend() is not None


# ------------------------------------------------------------------ data generators

def test_isa_pressure_limits_and_monotonicity():
    h = np.linspace(0.0, 20_000.0, 401)
    p = impl.isa_pressure(h)
    assert np.isclose(p[0], 101_325.0)
    assert np.all(np.diff(p) < 0), "pressure must fall monotonically with altitude"
    # tropopause value from the standard-atmosphere derivation: 22 632 Pa (within 0.1 %)
    assert np.isclose(impl.isa_pressure(np.array([11_000.0]))[0], 22_632.0, rtol=1e-3)
    # continuity across the tropopause
    assert np.isclose(impl.isa_pressure(np.array([10_999.0]))[0],
                      impl.isa_pressure(np.array([11_001.0]))[0], rtol=1e-3)


def test_drag_polar_is_quadratic_and_symmetric():
    cl = np.array([-0.5, 0.0, 0.5])
    cd = impl.drag_polar(cl, cd0=0.02, k=0.05)
    assert np.isclose(cd[1], 0.02)
    assert np.isclose(cd[0], cd[2])
    assert np.isclose(cd[2] - cd[1], 0.05 * 0.25)


def test_synthetic_thrust_curve_is_reproducible_and_shaped():
    t1, f1, s1 = impl.synthetic_thrust_curve(seed=3, n=20)
    t2, f2, s2 = impl.synthetic_thrust_curve(seed=3, n=20)
    assert t1.shape == f1.shape == s1.shape == (20,)
    assert np.array_equal(f1, f2), "same seed must give the same samples"
    assert np.all(s1 > 0)


# ------------------------------------------------------------------ linter and saving

def test_check_figure_passes_house_figures():
    for builder in (impl.figure_atmosphere, impl.figure_drag_polar, impl.figure_thrust_curve):
        fig = builder()
        assert impl.check_figure(fig) == [], builder.__name__


def test_check_figure_flags_missing_units_and_legend():
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1], label="a", color="magenta")
    ax.plot([0, 1], [1, 0], label="b", color="magenta")
    ax.set_xlabel("time")  # no unit tag, no y label, no legend, off-palette color
    problems = impl.check_figure(fig)
    joined = " ".join(problems).lower()
    assert "unit" in joined
    assert "y label" in joined
    assert "legend" in joined
    assert "color" in joined


def test_save_vector_writes_svg_and_rejects_png(tmp_path):
    fig = impl.figure_drag_polar()
    path = impl.save_vector(fig, tmp_path / "polar.svg")
    assert path.exists()
    head = path.read_text(encoding="utf-8")[:300].lower()
    assert "<svg" in head or "<?xml" in head
    with pytest.raises(ValueError):
        impl.save_vector(fig, tmp_path / "polar.png")


def test_main_writes_three_vector_figures(tmp_path):
    written = impl.main(tmp_path / "figs")
    assert len(written) == 3
    for p in written:
        assert p.exists() and p.suffix == ".svg" and p.stat().st_size > 1000
