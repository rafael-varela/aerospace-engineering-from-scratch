"""Public tests for P05.17 CAD & Technical Drawings I. The code checks the drawing, not the CAD file,
and nothing here imports a CAD package. The tests assert properties: the Green's-theorem polygon
formulas reproduce rectangles and triangles, ignore vertex order and the starting vertex; circles
approximated by n-gons converge at second order to the exact formulas; the bracket's mass properties
agree with a raster integration; the inertia tensor obeys the perpendicular-axis and thin-plate limits;
worst case bounds RSS and a seeded Monte Carlo reproduces the RSS statistics; the lesson's dimension
scheme is fully constrained without redundancy, an extra dimension is reported redundant and a missing
one leaves a degree of freedom, and the geometry is rebuilt from the dimensions alone.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)


# ---------------------------------------------------------------- polygons

def test_rectangle_properties_match_closed_forms():
    b, h = 30.0, 12.0
    P = np.array([[0, 0], [b, 0], [b, h], [0, h]], dtype=float)
    s = impl.polygon_properties(P)
    assert s["A"] == pytest.approx(b * h)
    assert s["Qy"] / s["A"] == pytest.approx(b / 2) and s["Qx"] / s["A"] == pytest.approx(h / 2)
    assert s["Ixx"] == pytest.approx(b * h ** 3 / 3)                 # about the base
    assert s["Iyy"] == pytest.approx(h * b ** 3 / 3)
    assert s["Ixy"] == pytest.approx(b ** 2 * h ** 2 / 4)


def test_triangle_and_order_independence():
    P = np.array([[0.0, 0.0], [6.0, 0.0], [0.0, 9.0]])
    s = impl.polygon_properties(P)
    assert s["A"] == pytest.approx(27.0)
    assert s["Qy"] / s["A"] == pytest.approx(2.0) and s["Qx"] / s["A"] == pytest.approx(3.0)
    assert s["Ixx"] == pytest.approx(6.0 * 9.0 ** 3 / 12)              # b h^3 / 12 about the base
    for Q in (P[::-1], np.roll(P, 1, axis=0)):                        # clockwise, other start vertex
        t = impl.polygon_properties(Q)
        for k in s:
            assert t[k] == pytest.approx(s[k])


def test_polygon_circle_converges_to_exact_at_second_order():
    exact = impl.circle_properties(3.0, -2.0, 5.0)
    errs = []
    for n in (16, 32, 64, 128):
        p = impl.polygon_properties(impl.regular_polygon(3.0, -2.0, 5.0, n))
        errs.append(abs(p["A"] - exact["A"]) / exact["A"])
        assert p["Qy"] / p["A"] == pytest.approx(3.0) and p["Qx"] / p["A"] == pytest.approx(-2.0)
    orders = np.log(np.array(errs[:-1]) / np.array(errs[1:])) / np.log(2)
    npt.assert_allclose(orders, 2.0, atol=0.05)


def test_bracket_area_properties_agree_with_raster_integration():
    b = impl.Bracket()
    s = impl.section_properties(b)
    r = impl.raster_properties(b, 0.1)
    # boundary cells give an O(h) error: about 2e-4 in area and 5e-4 in second moments at h = 0.1 mm
    assert r["A"] == pytest.approx(s["A"], rel=1e-3)
    assert r["xc"] == pytest.approx(s["xc"], abs=0.01) and r["yc"] == pytest.approx(s["yc"], abs=0.01)
    for k in ("Ixx_c", "Iyy_c", "Ixy_c"):
        assert r[k] == pytest.approx(s[k], rel=2e-3)
    s64 = impl.section_properties(b, 256)
    assert s64["A"] == pytest.approx(s["A"], rel=1e-4)


def test_net_area_is_outline_minus_holes_and_centroid_moves_away_from_holes():
    b = impl.Bracket()
    s = impl.section_properties(b)
    holes = sum(np.pi * r ** 2 for _, _, r in b.holes)
    assert s["A"] == pytest.approx(impl.polygon_properties(b.outline)["A"] - holes)
    solid = impl.section_properties(impl.Bracket(holes=[]))
    big = impl.section_properties(impl.Bracket(holes=[(10.0, 10.0, 5.0)]))
    assert big["xc"] > solid["xc"] and big["yc"] > solid["yc"]        # a hole near the corner pushes it away


def test_principal_axes_by_eigh():
    b = impl.Bracket()
    s = impl.section_properties(b)
    I1, I2, ang = impl.principal_axes(s["Ixx_c"], s["Iyy_c"], s["Ixy_c"])
    assert I1 + I2 == pytest.approx(s["Ixx_c"] + s["Iyy_c"])
    c, si = np.cos(ang), np.sin(ang)
    # I about an axis at angle ang: Ixx c^2 + Iyy s^2 - 2 Ixy s c
    assert s["Ixx_c"] * c ** 2 + s["Iyy_c"] * si ** 2 - 2 * s["Ixy_c"] * si * c == pytest.approx(I1, rel=1e-10)


# ---------------------------------------------------------------- mass properties

def test_mass_properties_scale_and_limits():
    b = impl.Bracket()
    m = impl.mass_properties(b)
    assert m["V"] == pytest.approx(impl.section_properties(b)["A"] * b.t)
    assert m["m"] == pytest.approx(m["V"] * b.rho * 1e-9)
    J = m["J"]
    npt.assert_allclose(J, J.T)
    assert np.all(np.linalg.eigvalsh(J) > 0)
    # thick-plate correction: Jzz = Jxx + Jyy - m t^2 / 6 (perpendicular-axis theorem plus the thickness)
    assert J[2, 2] == pytest.approx(J[0, 0] + J[1, 1] - m["m"] * b.t ** 2 / 6, rel=1e-12)
    m2 = impl.mass_properties(impl.Bracket(t=8.0))
    assert m2["m"] == pytest.approx(2 * m["m"])
    assert m2["com"][2] == pytest.approx(4.0)


def test_rectangular_plate_inertia_matches_the_box_formula():
    a, c, t, rho = 50.0, 20.0, 4.0, 2780.0
    b = impl.Bracket(outline=np.array([[0, 0], [a, 0], [a, c], [0, c]], dtype=float), holes=[], t=t, rho=rho)
    m = impl.mass_properties(b)
    mass = rho * 1e-9 * a * c * t
    npt.assert_allclose(np.diag(m["J"]), [mass * (c ** 2 + t ** 2) / 12, mass * (a ** 2 + t ** 2) / 12,
                                          mass * (a ** 2 + c ** 2) / 12], rtol=1e-12)
    assert abs(m["J"][0, 1]) < 1e-12


# ---------------------------------------------------------------- tolerances

def test_worst_case_bounds_rss_and_single_term_equal():
    rng = np.random.default_rng(1701)
    for _ in range(10):
        a = rng.normal(size=4)
        t = rng.uniform(0.01, 0.2, 4)
        assert impl.worst_case_stack(a, t) >= impl.rss_stack(a, t)
    assert impl.worst_case_stack([1.0], [0.1]) == pytest.approx(impl.rss_stack([1.0], [0.1]))
    assert impl.rss_stack([1, -1], [0.1, 0.1]) == pytest.approx(0.1 * np.sqrt(2))


def test_monte_carlo_reproduces_rss_statistics():
    nom = np.array([62.0, 60.0, 0.5, 0.5])
    a = np.array([1.0, -1.0, -1.0, -1.0])
    t = np.array([0.10, 0.05, 0.02, 0.02])
    y = impl.monte_carlo_stack(nom, a, t, n=200_000, seed=3)
    rss = impl.rss_stack(a, t)
    # with 2e5 samples the standard error of the mean is 2e-4 sigma and of the std about 0.16 %
    assert np.mean(y) == pytest.approx(nom @ a, abs=5 * rss / 3 / np.sqrt(2e5))
    assert 3 * np.std(y) == pytest.approx(rss, rel=0.01)
    frac = np.mean(np.abs(y - nom @ a) <= rss)
    assert frac == pytest.approx(0.9973, abs=0.001)                   # 3-sigma coverage, binomial sd 1.2e-4
    yu = impl.monte_carlo_stack(nom, a, t, n=200_000, seed=3, dist="uniform")
    assert np.max(np.abs(yu - nom @ a)) <= impl.worst_case_stack(a, t)
    assert np.std(yu) == pytest.approx(rss / np.sqrt(3), rel=0.01)    # uniform: sigma_i = t_i / sqrt(3)


def test_floating_fastener():
    assert impl.floating_fastener_position_tolerance(6.6, 6.0) == pytest.approx(0.3)
    assert impl.floating_fastener_position_tolerance(6.0, 6.0) == 0.0


# ---------------------------------------------------------------- the dimension scheme

def test_lesson_scheme_is_fully_constrained_without_redundancy():
    b = impl.Bracket()
    z, cons, n = impl.bracket_sketch(b)
    st = impl.constraint_status(z, cons, n)
    assert st["dof"] == 0 and st["redundant"] == 0
    npt.assert_allclose(impl.constraint_residuals(z, cons, n), 0.0, atol=1e-12)


def test_extra_dimension_is_redundant_and_missing_one_leaves_a_dof():
    b = impl.Bracket()
    z, cons, n = impl.bracket_sketch(b)
    st = impl.constraint_status(z, cons + [("dy", 2, 3, 25.0)], n)
    assert st["redundant"] == 1 and st["dof"] == 0
    st = impl.constraint_status(z, cons[:-1], n)                     # drop the last hole's radius
    assert st["dof"] == 1 and st["redundant"] == 0
    st = impl.constraint_status(z, [c for c in cons if c[0] != "fixx"], n)
    assert st["dof"] == 1                                            # the whole part can slide in x


def test_sketch_is_rebuilt_from_the_dimensions_alone():
    b = impl.Bracket()
    z, cons, n = impl.bracket_sketch(b)
    rng = np.random.default_rng(1702)
    zs = impl.solve_sketch(z + rng.normal(0.0, 1.0, z.size), cons, n)
    npt.assert_allclose(zs, z, atol=1e-8)
    outline = zs[:10].reshape(5, 2)
    assert impl.polygon_properties(outline)["A"] == pytest.approx(impl.polygon_properties(b.outline)["A"])


def test_distance_and_angle_constraints():
    z = np.array([0.0, 0.0, 3.0, 4.0, 3.0, 0.0])
    r = impl.constraint_residuals(z, [("dist", 0, 1, 5.0), ("angle", 1, 2, 0, np.pi / 2)], 3)
    npt.assert_allclose(r, 0.0, atol=1e-12)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "bracket_check.png").stat().st_size > 10_000
