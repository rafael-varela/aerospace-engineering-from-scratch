"""Public tests for P05.05 Centroids & Moments of Inertia. They assert geometry, not magic numbers: the
polygon integrals give bh^3/12 for a rectangle and the triangle formulas, ignore vertex order and
starting vertex, and move correctly with a translation; inscribed n-gons converge to the circle at second
order; composite parallel-axis sections equal polygon (Green's theorem) sections for the I spar, the Z
stringer and a plate with a hole; rotation preserves the two invariants; the principal axes from Mohr's
circle equal numpy.linalg.eigh and make the product of inertia vanish; symmetric sections have Ixy = 0;
the thin-walled error falls like (t/b)^2; the lumped boom idealisation beats booms alone; the centroid
of masses is a weighted mean.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)


def _close(s1, s2, rtol=1e-10, scale=None):
    a, b = np.array(s1.as_tuple()), np.array(s2.as_tuple())
    scale = np.abs(a).max() if scale is None else scale
    npt.assert_allclose(a, b, rtol=rtol, atol=rtol * scale)


# ---------------------------------------------------------------- polygon integrals

def test_rectangle_and_triangle_closed_forms():
    b, h = 0.04, 0.10
    s = impl.polygon_section(impl.rect_polygon(b, h, 0.3, -0.2))
    assert s.A == pytest.approx(b * h)
    assert (s.cx, s.cy) == pytest.approx((0.3, -0.2))
    assert s.Ix == pytest.approx(b * h ** 3 / 12) and s.Iy == pytest.approx(h * b ** 3 / 12)
    assert abs(s.Ixy) < 1e-12 * s.A * (0.3 ** 2 + 0.2 ** 2)      # round-off of the parallel-axis shift
    # right triangle with legs b (x) and h (y): centroid (b/3, h/3), Ix = b h^3/36, Ixy = -b^2 h^2/72
    tri = impl.polygon_section(np.array([[0, 0], [b, 0], [0, h]]))
    assert (tri.A, tri.cx, tri.cy) == pytest.approx((b * h / 2, b / 3, h / 3))
    assert tri.Ix == pytest.approx(b * h ** 3 / 36)
    assert tri.Ixy == pytest.approx(-b ** 2 * h ** 2 / 72)


def test_polygon_is_independent_of_vertex_order_start_and_translation():
    P = impl.z_section_polygon(0.025, 0.015, 0.0016)
    s = impl.polygon_section(P)
    _close(impl.polygon_section(P[::-1]), s)
    _close(impl.polygon_section(np.roll(P, 3, axis=0)), s)
    shifted = impl.polygon_section(P + np.array([1.5, -0.7]))
    npt.assert_allclose([shifted.cx - 1.5, shifted.cy + 0.7], [s.cx, s.cy], atol=1e-12)
    npt.assert_allclose([shifted.Ix, shifted.Iy, shifted.Ixy], [s.Ix, s.Iy, s.Ixy], rtol=1e-7)


def test_parallel_axis_theorem_on_raw_integrals():
    P = impl.angle_polygon(0.03, 0.02, 0.002)
    A, Qx, Qy, Ixx, Iyy, Ixy = impl.polygon_integrals(P)
    s = impl.polygon_section(P)
    assert Ixx == pytest.approx(s.Ix + A * s.cy ** 2)
    assert Iyy == pytest.approx(s.Iy + A * s.cx ** 2)
    assert Ixy == pytest.approx(s.Ixy + A * s.cx * s.cy)
    assert s.Ix < Ixx and s.Iy < Iyy                         # the centroidal axis gives the least I


def test_inscribed_polygons_converge_to_the_circle_at_second_order():
    r = 0.01
    errs = []
    for n in (16, 32, 64, 128):
        s = impl.polygon_section(impl.circle_polygon(r, n))
        errs.append(abs(s.Ix - np.pi * r ** 4 / 4) / (np.pi * r ** 4 / 4))
        assert s.Ix == pytest.approx(s.Iy, rel=1e-12) and abs(s.Ixy) < 1e-12 * s.Ix
    orders = np.diff(np.log(errs)) / np.log(0.5)
    assert np.all((orders > 1.9) & (orders < 2.1))
    c = impl.circle(r, 0.0, 0.0)
    assert c.Ix == pytest.approx(np.pi * r ** 4 / 4) and c.A == pytest.approx(np.pi * r ** 2)


# ---------------------------------------------------------------- composite against polygon

def test_i_spar_composite_equals_polygon_and_textbook_formula():
    h, b, tf, tw = 0.150, 0.050, 0.006, 0.003
    sc = impl.composite_section(impl.i_section_parts(h, b, tf, tw))
    sp = impl.polygon_section(impl.i_section_polygon(h, b, tf, tw))
    _close(sc, sp, rtol=1e-10)
    hw = h - 2 * tf
    assert sc.Ix == pytest.approx((b * h ** 3 - (b - tw) * hw ** 3) / 12)       # box minus the two cut-outs
    assert sc.Ixy == pytest.approx(0.0, abs=1e-20)                             # symmetric: no product


def test_z_stringer_composite_equals_polygon():
    for h, b, t in ((0.025, 0.015, 0.0016), (0.040, 0.020, 0.001), (0.030, 0.010, 0.003)):
        _close(impl.composite_section(impl.z_section_parts(h, b, t)), impl.polygon_section(impl.z_section_polygon(h, b, t)))


def test_holes_by_negative_parts_and_negative_polygons_agree():
    # 60 x 40 mm plate with an off-centre 10 mm radius hole: composite (exact circle) vs 2048-gon hole
    plate = impl.rect(0.06, 0.04, 0.0, 0.0)
    hole = impl.circle(0.01, 0.012, 0.005, sign=-1.0)
    sc = impl.composite_section([plate, hole])
    sp = impl.polygon_section(impl.rect_polygon(0.06, 0.04), [impl.circle_polygon(0.01, 2048, 0.012, 0.005)])
    _close(sc, sp, rtol=2e-5)
    assert sc.cx < 0 and sc.cy < 0                           # the centroid moves away from the hole
    tube = impl.composite_section([impl.circle(0.01, 0, 0), impl.circle(0.009, 0, 0, sign=-1.0)])
    assert tube.Ix == pytest.approx(np.pi / 4 * (0.01 ** 4 - 0.009 ** 4))


# ---------------------------------------------------------------- rotation and principal axes

def test_rotation_preserves_invariants_and_quarter_turn_swaps_axes():
    s = impl.polygon_section(impl.angle_polygon(0.03, 0.02, 0.002))
    for th in np.linspace(-1.5, 1.5, 7):
        Iu, Iv, Iuv = impl.rotate_axes(s.Ix, s.Iy, s.Ixy, th)
        assert Iu + Iv == pytest.approx(s.Ix + s.Iy)
        assert Iu * Iv - Iuv ** 2 == pytest.approx(s.Ix * s.Iy - s.Ixy ** 2)
    Iu, Iv, Iuv = impl.rotate_axes(s.Ix, s.Iy, s.Ixy, np.pi / 2)
    assert (Iu, Iv, Iuv) == pytest.approx((s.Iy, s.Ix, -s.Ixy))


def test_rotate_axes_matches_direct_integration_of_the_rotated_polygon():
    P = impl.z_section_polygon(0.025, 0.015, 0.0016)
    s = impl.polygon_section(P)
    th = 0.37
    c, sn = np.cos(th), np.sin(th)
    # coordinates in the rotated frame: u = x c + y s, v = -x s + y c
    Q = (P - [s.cx, s.cy]) @ np.array([[c, -sn], [sn, c]])
    r = impl.polygon_section(Q)
    npt.assert_allclose(impl.rotate_axes(s.Ix, s.Iy, s.Ixy, th), (r.Ix, r.Iy, r.Ixy), rtol=1e-9, atol=1e-9 * s.Ix)


def test_principal_axes_match_eigh_and_zero_the_product():
    rng = np.random.default_rng(5051)
    shapes = [impl.z_section_polygon(0.025, 0.015, 0.0016), impl.angle_polygon(0.03, 0.02, 0.002),
              impl.z_section_polygon(0.025, 0.005, 0.0016)]
    for P in shapes:
        s = impl.polygon_section(P)
        I1, I2, th = impl.principal_axes(s.Ix, s.Iy, s.Ixy)
        w, V = impl.principal_axes_eig(s.Ix, s.Iy, s.Ixy)
        npt.assert_allclose([I2, I1], w, rtol=1e-12)
        Iu, Iv, Iuv = impl.rotate_axes(s.Ix, s.Iy, s.Ixy, th)
        assert Iu == pytest.approx(I1) and Iv == pytest.approx(I2)
        assert abs(Iuv) < 1e-12 * I1
        n_max = V[:, 1]                                       # normal of the I_max axis
        e_max = np.array([n_max[1], -n_max[0]])
        assert abs(abs(np.dot(e_max, [np.cos(th), np.sin(th)])) - 1.0) < 1e-10
        # no axis has more than I1 or less than I2
        for t in rng.uniform(-np.pi, np.pi, 20):
            Iu = impl.rotate_axes(s.Ix, s.Iy, s.Ixy, t)[0]
            assert I2 * (1 - 1e-12) <= Iu <= I1 * (1 + 1e-12)


def test_symmetric_sections_have_principal_axes_along_x_and_y():
    s = impl.polygon_section(impl.i_section_polygon(0.15, 0.05, 0.006, 0.003))
    I1, I2, th = impl.principal_axes(s.Ix, s.Iy, s.Ixy)
    assert th == pytest.approx(0.0, abs=1e-12) and I1 == pytest.approx(s.Ix) and I2 == pytest.approx(s.Iy)


# ---------------------------------------------------------------- approximations and helpers

def test_thin_walled_error_falls_like_t_squared():
    h, b = 0.025, 0.015
    errs = []
    for t in (0.0002, 0.0004, 0.0008):
        exact = impl.polygon_section(impl.z_section_polygon(h, b, t))
        thin = impl.z_section_thin(h, b, t)
        assert thin.A == pytest.approx(exact.A, rel=1e-12)
        assert np.sign(thin.Ixy) == np.sign(exact.Ixy)
        errs.append(abs(thin.Ix / exact.Ix - 1))
    ratios = np.array(errs[1:]) / np.array(errs[:-1])
    npt.assert_allclose(ratios, 4.0, rtol=0.05)


def test_wing_box_idealisations():
    wb = dict(width=0.40, height=0.12, t_skin=0.0016, t_web=0.002, cap_area=300e-6, n_stringers=4, stringer_area=60e-6)
    full = impl.composite_section(impl.wing_box_parts(**wb))
    booms = impl.composite_section(impl.wing_box_parts(**wb, idealisation="booms"))
    lump = impl.composite_section(impl.wing_box_parts(**wb, idealisation="lumped"))
    assert booms.Ix < full.Ix
    assert abs(lump.Ix / full.Ix - 1) < 0.05 < abs(booms.Ix / full.Ix - 1)
    assert (full.cx, full.cy) == pytest.approx((0.0, 0.0), abs=1e-15)
    assert abs(full.Ixy) < 1e-12 * full.Ix
    with pytest.raises(ValueError):
        impl.wing_box_parts(**wb, idealisation="magic")


def test_section_moduli_and_radius_of_gyration():
    b, h = 0.02, 0.06
    s = impl.polygon_section(impl.rect_polygon(b, h))
    St, Sb = impl.section_moduli(s, impl.rect_polygon(b, h))
    assert St == pytest.approx(b * h ** 2 / 6) and Sb == pytest.approx(St)
    assert impl.radius_of_gyration(s.Ix, s.A) == pytest.approx(h / np.sqrt(12))
    T = np.array([[0, 0], [b, 0], [0, h]])
    st = impl.polygon_section(T)
    St, Sb = impl.section_moduli(st, T)
    assert Sb == pytest.approx(2 * St)                       # the base is twice as close to the centroid


def test_centroid_of_masses_and_percent_mac():
    m = np.array([700.0, 80.0, 20.0])
    x = np.array([2.1, 2.0, 3.0])
    xcg = impl.composite_centroid(m, x)
    assert xcg == pytest.approx(np.dot(m, x) / m.sum())
    assert impl.composite_centroid(2 * m, x) == pytest.approx(xcg)                   # scale-free
    P = np.array([[0.0, 1.0], [2.0, 3.0]])
    npt.assert_allclose(impl.composite_centroid([1.0, 1.0], P), [1.0, 2.0])
    assert impl.cg_percent_mac(1.8 + 0.375, 1.8, 1.5) == pytest.approx(25.0)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "section_properties.png").exists()
