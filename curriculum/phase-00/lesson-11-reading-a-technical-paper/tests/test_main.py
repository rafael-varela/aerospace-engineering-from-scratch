"""Public tests for P00.11 Reading a Technical Paper.

They assert geometry and reading rules, not one magic number: the camber line peaks where the
designation says, the thickness polynomial has its known shape, the lift-curve fit recovers a
known slope, the drag-polar reader finds the analytic L/D maximum, and the notes validator
reports exactly what is missing.

Run `pytest curriculum/phase-00/lesson-11-reading-a-technical-paper` against your starter/main.py,
or `AEFS_IMPL=solution pytest ...` against the reference solution.
"""
import numpy as np
import pytest

from conftest import load_impl

impl = load_impl(__file__)


# --- designation parser --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "designation, m, p, t",
    [("2412", 0.02, 0.4, 0.12), ("0012", 0.0, 0.0, 0.12), ("4415", 0.04, 0.4, 0.15), ("6409", 0.06, 0.4, 0.09)],
)
def test_parse_naca4_values(designation, m, p, t):
    s = impl.parse_naca4(designation)
    assert s.max_camber == pytest.approx(m)
    assert s.camber_position == pytest.approx(p)
    assert s.thickness == pytest.approx(t)
    assert s.designation == designation
    assert s.is_symmetric == (m == 0.0)


@pytest.mark.parametrize("spelling", ["NACA 2412", "naca2412", "  2412 ", "Naca  2412"])
def test_parse_accepts_common_spellings(spelling):
    assert impl.parse_naca4(spelling) == impl.parse_naca4("2412")


@pytest.mark.parametrize("bad", ["", "241", "24120", "23012", "24a2", "2400", "2012", "NACA 64-212"])
def test_parse_rejects_malformed_designations(bad):
    with pytest.raises(ValueError):
        impl.parse_naca4(bad)


# --- geometry ------------------------------------------------------------------------------------

def test_symmetric_camber_line_is_zero():
    x = np.linspace(0, 1, 51)
    assert np.all(impl.camber_line(x, 0.0, 0.0) == 0.0)
    assert np.all(impl.camber_slope(x, 0.0, 0.0) == 0.0)


@pytest.mark.parametrize("m, p", [(0.02, 0.4), (0.04, 0.4), (0.06, 0.3), (0.03, 0.6)])
def test_camber_line_peaks_at_p_with_value_m(m, p):
    x = np.linspace(0, 1, 1001)
    yc = impl.camber_line(x, m, p)
    assert yc[0] == pytest.approx(0.0, abs=1e-12)
    assert yc[-1] == pytest.approx(0.0, abs=1e-12)
    assert np.max(yc) == pytest.approx(m, rel=1e-6)
    assert x[np.argmax(yc)] == pytest.approx(p, abs=1e-3)
    # continuous value and zero slope where the two arcs meet
    eps = 1e-9
    assert impl.camber_line(np.array([p - eps]), m, p)[0] == pytest.approx(impl.camber_line(np.array([p + eps]), m, p)[0], abs=1e-8)
    assert impl.camber_slope(np.array([p]), m, p)[0] == pytest.approx(0.0, abs=1e-8)


def test_camber_slope_is_derivative_of_camber_line():
    x = np.linspace(0.01, 0.99, 400)
    m, p = 0.04, 0.4
    h = 1e-6
    numeric = (impl.camber_line(x + h, m, p) - impl.camber_line(x - h, m, p)) / (2 * h)
    assert np.allclose(numeric, impl.camber_slope(x, m, p), atol=1e-6)


def test_thickness_distribution_shape():
    x = np.linspace(0, 1, 2001)
    for t in (0.09, 0.12, 0.15):
        yt = impl.thickness_distribution(x, t)
        assert yt[0] == pytest.approx(0.0, abs=1e-12)
        assert np.max(yt) == pytest.approx(t / 2, rel=0.01)           # max thickness = t
        assert 0.25 < x[np.argmax(yt)] < 0.35                          # at about 30 % chord
        assert 0 < yt[-1] < 0.02 * t                                   # small, finite trailing edge
        assert np.all(yt >= 0)
    # linear in t
    assert np.allclose(impl.thickness_distribution(x, 0.24), 2 * impl.thickness_distribution(x, 0.12))


def test_airfoil_surfaces_symmetric_section_mirrors():
    xu, yu, xl, yl = impl.airfoil_surfaces("0012", 101)
    assert np.allclose(xu, xl)
    assert np.allclose(yu, -yl)
    assert xu[0] == pytest.approx(0.0) and xu[-1] == pytest.approx(1.0)
    assert np.max(yu) == pytest.approx(0.06, rel=0.02)


def test_airfoil_surfaces_cambered_section_upper_above_lower():
    xu, yu, xl, yl = impl.airfoil_surfaces("4415", 201)
    assert np.all(yu[1:-1] > yl[1:-1])
    # cambered: mean of upper and lower near the max-camber station sits near m
    i = np.argmin(np.abs(xu - 0.4))
    assert 0.5 * (yu[i] + yl[i]) == pytest.approx(0.04, abs=0.005)
    assert len(xu) == len(yu) == len(xl) == len(yl) == 201


# --- reading curves --------------------------------------------------------------------------------

def test_lift_curve_fit_recovers_slope_and_zero_lift_angle():
    alpha = np.arange(-8.0, 17.0, 1.0)
    a0 = 0.11
    cl = a0 * (alpha + 2.0)
    cl = np.where(alpha > 10, a0 * 12 - 0.02 * (alpha - 10), cl)      # stall flattening
    slope, alpha0 = impl.lift_curve_fit(alpha, cl)
    assert slope == pytest.approx(a0, rel=1e-6)
    assert alpha0 == pytest.approx(-2.0, abs=1e-6)
    # thin-airfoil theory gives 2*pi per radian = 0.1097 per degree
    assert abs(slope - 2 * np.pi * np.pi / 180) < 0.01


def test_lift_curve_fit_needs_two_points_in_range():
    with pytest.raises(ValueError):
        impl.lift_curve_fit(np.array([12.0, 14.0]), np.array([1.2, 1.3]))


def test_drag_polar_summary_matches_parabolic_polar():
    cd0, k = 0.008, 0.04
    cl = np.linspace(-0.4, 1.4, 1801)
    cd = cd0 + k * cl**2
    s = impl.drag_polar_summary(cl, cd)
    assert s["cd_min"] == pytest.approx(cd0, rel=1e-9)
    assert s["cl_at_cd_min"] == pytest.approx(0.0, abs=1e-3)
    assert s["cl_at_ld_max"] == pytest.approx(np.sqrt(cd0 / k), abs=2e-3)   # analytic optimum
    assert s["ld_max"] == pytest.approx(1 / (2 * np.sqrt(cd0 * k)), rel=1e-4)


def test_drag_polar_summary_ignores_negative_lift_for_ld():
    cl = np.array([-0.5, -0.2, 0.1, 0.3])
    cd = np.array([0.001, 0.009, 0.010, 0.012])
    s = impl.drag_polar_summary(cl, cd)
    assert s["cl_at_ld_max"] == pytest.approx(0.3)


# --- reading notes --------------------------------------------------------------------------------

def test_example_notes_are_complete():
    assert impl.validate_reading_notes(impl.EXAMPLE_NOTES) == []


def test_missing_field_is_reported():
    text = impl.EXAMPLE_NOTES.replace("## Limitations", "## Caveats")
    problems = impl.validate_reading_notes(text)
    assert problems == ["missing field: Limitations"]


@pytest.mark.parametrize("placeholder", ["", "...", "TODO", "tbd"])
def test_placeholder_body_is_reported_as_empty(placeholder):
    lines = impl.EXAMPLE_NOTES.split("## Result\n")
    text = lines[0] + "## Result\n" + placeholder + "\n\n## Limitations" + lines[1].split("## Limitations", 1)[1]
    problems = impl.validate_reading_notes(text)
    assert problems == ["empty field: Result"]


def test_notes_headings_are_case_insensitive_and_any_level():
    text = impl.EXAMPLE_NOTES.replace("## Key figure", "### KEY FIGURE")
    assert impl.validate_reading_notes(text) == []


# --- end to end -----------------------------------------------------------------------------------

def test_main_writes_plot(tmp_path, capsys):
    impl.main(tmp_path)
    assert (tmp_path / "naca_airfoils.png").stat().st_size > 1000
    out = capsys.readouterr().out
    assert "NACA 2412" in out and "symmetric" in out and "complete" in out
