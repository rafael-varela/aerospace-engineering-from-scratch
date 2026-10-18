"""Public tests for P10.02 Airfoil Geometry. They assert geometry, not one magic number: cosine
spacing clusters at both ends and is symmetric; the thickness form peaks at 30 % chord with the
designated thickness, closes when asked, and has the leading-edge radius its sqrt term implies; the
camber lines are continuous with zero slope at their maximum; the 5-digit constants put the camber
maximum at p and give the design ideal lift coefficient (and reproduce the quoted 230 values); the
Selig ordering is closed and counter-clockwise; the shoelace area matches the exact integral and
converges at order 2 with cosine spacing and only 1.5 with uniform spacing; and symmetric sections
are symmetric.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)


def test_cosine_spacing_endpoints_symmetry_and_clustering():
    x = impl.cosine_spacing(40)
    assert x.size == 41 and x[0] == 0.0 and x[-1] == pytest.approx(1.0)
    npt.assert_allclose(x + x[::-1], 1.0, atol=1e-14)                 # symmetric about mid-chord
    d = np.diff(x)
    assert d[0] < d[20] / 10 and d[-1] < d[20] / 10                    # clustered at both edges
    assert np.all(d > 0)


def test_thickness_form_maximum_is_t_at_30_percent():
    for t in (0.06, 0.12, 0.21):
        xm, tm = impl.max_thickness(t)
        assert xm == pytest.approx(0.30, abs=0.002)
        assert tm == pytest.approx(t, rel=5e-4)
    xm, tm = impl.max_thickness(0.12, closed_te=True)
    assert xm == pytest.approx(0.30, abs=0.002) and tm == pytest.approx(0.12, rel=5e-4)


def test_trailing_edge_open_and_closed():
    assert 2 * impl.naca4_thickness(1.0, 0.12) == pytest.approx(0.021 * 0.12, rel=1e-6)
    assert abs(impl.naca4_thickness(1.0, 0.12, closed_te=True)) < 1e-12              # coefficients sum to zero
    assert impl.naca4_thickness(0.0, 0.12) == 0.0


def test_leading_edge_radius_from_the_sqrt_term():
    t = 0.15
    x = np.array([1e-8, 1e-9])
    # near x = 0, y_t^2 / (2 x) -> r
    npt.assert_allclose(impl.naca4_thickness(x, t) ** 2 / (2 * x), impl.leading_edge_radius(t), rtol=1e-3)
    assert impl.leading_edge_radius(0.12) == pytest.approx(1.1019 * 0.12 ** 2, rel=1e-4)


def test_thickness_area_coefficient_is_the_exact_integral():
    from scipy import integrate
    for closed in (False, True):
        exact, _ = integrate.quad(lambda x: 2 * impl.naca4_thickness(x, 1.0, closed), 0, 1, limit=200)
        assert impl.thickness_area_coefficient(closed) == pytest.approx(exact, rel=1e-8)
    assert impl.thickness_area_coefficient() == pytest.approx(0.685, abs=5e-4)


def test_naca4_camber_line_is_smooth_at_p_and_peaks_there():
    m, p = 0.04, 0.4
    yc, dyc = impl.naca4_camber(np.array([p - 1e-9, p, p + 1e-9]), m, p)
    npt.assert_allclose(yc, m, atol=1e-9)
    npt.assert_allclose(dyc, 0.0, atol=1e-7)
    x = np.linspace(0, 1, 2001)
    yc, dyc = impl.naca4_camber(x, m, p)
    assert yc[0] == 0.0 and abs(yc[-1]) < 1e-15
    npt.assert_allclose(np.gradient(yc, x)[5:-5], dyc[5:-5], atol=2e-4)  # slope is the derivative
    y0, d0 = impl.naca4_camber(x, 0.0, 0.0)
    assert np.all(y0 == 0) and np.all(d0 == 0)


@pytest.mark.parametrize("p", [0.05, 0.10, 0.15, 0.20, 0.25])
def test_naca5_constants_satisfy_their_defining_conditions(p):
    r, k1 = impl.naca5_constants(p, 0.3)
    assert p < r < 1
    yc, dyc = impl.naca5_camber(np.array([p]), r, k1)
    assert abs(dyc[0]) < 1e-10                                         # camber maximum at x = p
    assert impl.naca5_ideal_cl(r, k1) == pytest.approx(0.3, rel=1e-8)  # design ideal cl
    yr, _ = impl.naca5_camber(np.array([r - 1e-12, r + 1e-12]), r, k1)
    assert yr[0] == pytest.approx(yr[1], rel=1e-9)                     # continuous at r
    # k1 scales with the design lift coefficient (the line is linear in k1)
    assert impl.naca5_constants(p, 0.6)[1] == pytest.approx(2 * k1)


def test_quoted_230_constants_are_reproduced():
    r, k1 = impl.naca5_constants(0.15)
    assert r == pytest.approx(0.2025, abs=5e-4)
    assert k1 == pytest.approx(15.957, rel=5e-3)
    assert impl.naca5_ideal_cl(0.2025, 15.957) == pytest.approx(0.3, abs=5e-4)


def test_parse_naca():
    d = impl.parse_naca("2412")
    assert (d["kind"], d["m"], d["p"], d["t"]) == (4, 0.02, 0.4, 0.12)
    d = impl.parse_naca("23012")
    assert d["kind"] == 5 and d["cl_design"] == pytest.approx(0.3) and d["p"] == pytest.approx(0.15)
    with pytest.raises(ValueError):
        impl.parse_naca("23112")


def test_selig_order_closed_counter_clockwise():
    X, Y = impl.naca_coordinates("2412", n=60)
    assert X.size == 121
    assert X[0] == pytest.approx(1.0, abs=1e-3) and X[-1] == pytest.approx(1.0, abs=1e-3)
    assert Y[0] > Y[-1]                                                 # starts on the upper surface
    i_le = int(np.argmin(X))
    assert i_le == 60 and np.all(Y[1:i_le] > Y[-2:i_le:-1])             # upper above lower at each station
    assert impl.polygon_area(X, Y) > 0                                   # counter-clockwise
    Xc, Yc = impl.naca_coordinates("2412", n=60, closed_te=True)
    assert Yc[0] == pytest.approx(Yc[-1], abs=1e-12)


def test_symmetric_section_is_symmetric_and_camber_is_recovered():
    X, Y = impl.naca_coordinates("0015", n=50)
    npt.assert_allclose(X, X[::-1], atol=1e-14)
    npt.assert_allclose(Y, -Y[::-1], atol=1e-14)
    X, Y = impl.naca_coordinates("4412", n=300)
    xs, th, cam = impl.surface_thickness_camber(X, Y)
    assert cam.max() == pytest.approx(0.04, rel=2e-3)
    assert xs[np.argmax(cam)] == pytest.approx(0.4, abs=0.01)
    assert th.max() == pytest.approx(0.12, rel=5e-3)                     # vertical vs perpendicular thickness


def test_area_matches_exact_and_scales_with_chord_squared():
    X, Y = impl.naca_coordinates("0012", n=400)
    assert impl.polygon_area(X, Y) == pytest.approx(impl.thickness_area_coefficient() * 0.12, rel=1e-4)
    X2, Y2 = impl.naca_coordinates("0012", n=400, chord=2.0)
    assert impl.polygon_area(X2, Y2) == pytest.approx(4 * impl.polygon_area(X, Y))
    # camber barely changes the area (thickness is laid perpendicular to a gently curved line)
    X4, Y4 = impl.naca_coordinates("2412", n=400)
    assert impl.polygon_area(X4, Y4) == pytest.approx(impl.polygon_area(X, Y), rel=3e-3)


def test_area_convergence_orders():
    ns = [20, 40, 80, 160]
    o_cos = impl.observed_order(impl.area_convergence("0012", ns, "cosine"))
    o_uni = impl.observed_order(impl.area_convergence("0012", ns, "uniform"))
    assert np.all(np.abs(o_cos - 2.0) < 0.1)
    assert np.all(np.abs(o_uni - 1.5) < 0.1)                             # sqrt(x) at the leading edge


def test_spar_box_area_is_part_of_the_total():
    total = impl.polygon_area(*impl.naca_coordinates("2412", n=400))
    box = impl.area_between("2412", 0.15, 0.60)
    assert 0.5 * total < box < total
    assert impl.area_between("2412", 0.15, 0.60, chord=1.2) == pytest.approx(1.44 * box)
    assert impl.area_between("2412", 0.0, 1.0) == pytest.approx(total, rel=2e-3)


def test_main_writes_outputs(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "airfoil_geometry.png").stat().st_size > 1000
    lines = (tmp_path / "naca2412.dat").read_text().splitlines()
    assert lines[0].startswith("NACA 2412") and len(lines) == 162
