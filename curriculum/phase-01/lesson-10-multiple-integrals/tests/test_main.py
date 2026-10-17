"""Public tests for P01.10 Multiple Integrals. They assert closed forms, symmetry, scaling and convergence, not magic numbers.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-01/parameterized-cases/.
"""
import functools

import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)


def _box(**kw):
    return impl.WingBox(**kw)


@functools.lru_cache(maxsize=None)
def _analytic(**kw):
    """Cache the SymPy result (about a second) so several tests can share it."""
    return impl.mass_properties_analytic(_box(**kw))


# ---------------------------------------------------------------- geometry and the change of variables

def test_chord_is_linear_and_validated():
    b = _box()
    assert impl.chord(0.0, b) == b.c_root and impl.chord(b.L, b) == b.c_tip
    assert impl.chord(b.L / 2, b) == pytest.approx(0.5 * (b.c_root + b.c_tip))
    y = np.linspace(0, b.L, 7)
    npt.assert_allclose(np.diff(impl.chord(y, b)), (b.c_tip - b.c_root) / 6, rtol=1e-12)
    with pytest.raises(ValueError):
        impl.chord(-0.1, b)
    with pytest.raises(ValueError):
        impl.chord(0.0, _box(f0=0.7, f1=0.6))


def test_map_to_box_corners_and_jacobian():
    b = _box()
    x, y, z, jac = impl.map_to_box(0.0, 0.0, 0.5, b)
    assert (x, y, z) == pytest.approx((b.f0 * b.c_root, 0.0, 0.0))                   # root front wall, mid-depth
    x, y, z, jac = impl.map_to_box(1.0, 1.0, 1.0, b)
    assert x == pytest.approx(b.L * np.tan(b.sweep) + b.f1 * b.c_tip)
    assert y == pytest.approx(b.L) and z == pytest.approx(0.5 * b.tau * b.c_tip)
    # Jacobian = product of the three edge lengths of the mapped cell: dx/dxi * dy/deta * dz/dzeta.
    assert jac == pytest.approx((b.f1 - b.f0) * b.c_tip * b.L * b.tau * b.c_tip)
    # Check the Jacobian against finite differences of the map at an interior point.
    h = 1e-6
    p = np.array([0.3, 0.6, 0.7])
    J = np.empty((3, 3))
    for j in range(3):
        e = np.zeros(3)
        e[j] = h
        plus = np.array(impl.map_to_box(*(p + e), b)[:3])
        minus = np.array(impl.map_to_box(*(p - e), b)[:3])
        J[:, j] = (plus - minus) / (2 * h)
    assert abs(np.linalg.det(J)) == pytest.approx(impl.map_to_box(*p, b)[3], rel=1e-6)


def test_midpoint_3d_integrates_constants_and_the_volume():
    b = _box()
    vol = impl.mass_closed_form(b) / b.rho
    assert impl.midpoint_3d(lambda x, y, z: np.ones_like(x), 1, b) == pytest.approx(
        (b.f1 - b.f0) * b.tau * b.L * ((b.c_root + b.c_tip) / 2) ** 2, rel=1e-12)   # n = 1: chord at mid-span
    assert impl.midpoint_3d(lambda x, y, z: np.ones_like(x), 64, b) == pytest.approx(vol, rel=1e-4)
    with pytest.raises(ValueError):
        impl.midpoint_3d(lambda x, y, z: np.ones_like(x), 0, b)


# ---------------------------------------------------------------- mass and centre of mass

def test_mass_three_ways():
    b = _box()
    m_a, cg_a, _ = _analytic()
    assert m_a == pytest.approx(impl.mass_closed_form(b), rel=1e-12)
    m_q, _ = impl.tplquad_mass_and_ixx(b)
    assert m_q == pytest.approx(m_a, rel=1e-8)
    assert cg_a[1] == pytest.approx(impl.spanwise_cg_closed_form(b), rel=1e-12)
    # Mass is linear in density and cubic in a uniform scale of all lengths.
    assert impl.mass_closed_form(_box(rho=2 * b.rho)) == pytest.approx(2 * m_a, rel=1e-12)
    k = 1.7
    assert impl.mass_closed_form(_box(L=k * b.L, c_root=k * b.c_root, c_tip=k * b.c_tip)) == pytest.approx(k**3 * m_a, rel=1e-12)


def test_midpoint_mass_converges_at_second_order():
    errs = impl.convergence_study([2, 4, 8, 16, 32], _box())
    assert np.all(errs[1:] < errs[:-1])
    rates = np.log2(errs[:-1] / errs[1:])
    npt.assert_allclose(rates, 2.0, atol=0.05)
    assert errs[-1] < 1e-4


def test_centre_of_mass_lies_inside_the_box_and_obeys_symmetry():
    b = _box()
    m, cg, _ = _analytic()
    assert 0.0 < cg[1] < b.L and cg[1] < b.L / 2                                       # tapered: inboard of mid-span
    xle = cg[1] * np.tan(b.sweep)
    c = impl.chord(cg[1], b)
    assert xle + b.f0 * c < cg[0] < xle + b.f1 * c                                      # within the box at that station
    assert cg[2] == pytest.approx(0.0, abs=1e-12)                                       # z-symmetric
    # Untapered: y_cg = L/2 exactly; the closed form reduces to it.
    m2, cg2, _ = _analytic(c_tip=b.c_root)
    assert cg2[1] == pytest.approx(b.L / 2, rel=1e-12)
    assert impl.spanwise_cg_closed_form(_box(c_tip=b.c_root)) == pytest.approx(b.L / 2, rel=1e-12)
    # The midpoint CG converges to the analytic one.
    m16, cg16, _ = impl.mass_properties_midpoint(16, b)
    npt.assert_allclose(cg16, cg, atol=5e-3)


# ---------------------------------------------------------------- the inertia tensor

def test_inertia_tensor_is_symmetric_positive_definite_and_obeys_triangle_inequalities():
    m, cg, I = _analytic()
    npt.assert_allclose(I, I.T, rtol=1e-12)
    w = np.linalg.eigvalsh(I)
    assert np.all(w > 0.0)
    Ixx, Iyy, Izz = np.diag(I)
    assert Ixx + Iyy >= Izz and Iyy + Izz >= Ixx and Izz + Ixx >= Iyy                   # true for any mass distribution
    assert I[0, 2] == pytest.approx(0.0, abs=1e-9) and I[1, 2] == pytest.approx(0.0, abs=1e-9)   # z-symmetry kills I_xz, I_yz
    assert I[0, 1] < 0.0                                                                 # -int x y rho dV with x, y > 0
    # Thin box: I_zz ~ I_xx + I_yy - 2 I_z-part ... simply: Izz > Ixx > Iyy for a long, slender span.
    assert Izz > Ixx > Iyy


def test_inertia_matches_midpoint_and_tplquad():
    b = _box()
    m, cg, I = _analytic()
    m32, cg32, I32 = impl.mass_properties_midpoint(32, b)
    npt.assert_allclose(I32, I, rtol=2e-4, atol=1e-6)
    _, ixx_q = impl.tplquad_mass_and_ixx(b)
    assert ixx_q == pytest.approx(I[0, 0], rel=1e-8)


def test_parallel_axis_theorem_round_trip_and_minimum_at_cg():
    m, cg, I_o = _analytic()
    I_cg = impl.inertia_about_cg(I_o, m, cg)
    npt.assert_allclose(impl.parallel_axis(I_cg, m, cg), I_o, rtol=1e-12)
    # The moments about any axis through the CG are the smallest among parallel axes.
    assert np.all(np.diag(I_cg) < np.diag(I_o))
    assert np.all(np.linalg.eigvalsh(I_cg) > 0.0)
    # Shifting along z by d adds m d**2 to I_xx and I_yy and nothing to I_zz.
    d = np.array([0.0, 0.0, 0.3])
    I_s = impl.parallel_axis(I_cg, m, d)
    assert I_s[0, 0] == pytest.approx(I_cg[0, 0] + m * 0.09) and I_s[2, 2] == pytest.approx(I_cg[2, 2])
    # Untapered, unswept box: symmetric about its mid-span and mid-chord, so I_xy about the CG vanishes.
    m2, cg2, I2 = _analytic(c_tip=1.6, sweep=0.0)
    I2cg = impl.inertia_about_cg(I2, m2, cg2)
    assert I2cg[0, 1] == pytest.approx(0.0, abs=1e-9)


def test_inertia_scales_as_density_and_length_to_the_fifth():
    m, cg, I = _analytic()
    m2, cg2, I2 = _analytic(rho=360.0)
    npt.assert_allclose(I2, 2 * I, rtol=1e-12)
    npt.assert_allclose(cg2, cg, rtol=1e-12)
    k = 2.0
    m3, cg3, I3 = _analytic(L=k * 5.0, c_root=k * 1.6, c_tip=k * 0.8)
    npt.assert_allclose(I3, k**5 * I, rtol=1e-12)
    npt.assert_allclose(cg3, k * cg, rtol=1e-12)


# ---------------------------------------------------------------- polar double integral

def test_ring_inertia_polar_matches_closed_form_and_limits():
    m, I = impl.ring_inertia_polar(0.95, 1.0, 0.05, 2700.0, 64)
    mc, Ic = impl.ring_inertia_closed_form(0.95, 1.0, 0.05, 2700.0)
    assert m == pytest.approx(mc, rel=1e-12)                                            # midpoint is exact for r dr
    assert I == pytest.approx(Ic, rel=1e-6)
    # Solid disc: I = 1/2 m R**2.  Thin shell: I -> m R**2.
    ms, Is = impl.ring_inertia_closed_form(0.0, 1.0, 0.05, 2700.0)
    assert Is == pytest.approx(0.5 * ms * 1.0**2, rel=1e-12)
    mt, It = impl.ring_inertia_polar(0.999, 1.0, 0.05, 2700.0, 16)
    assert It / (mt * 1.0**2) == pytest.approx(1.0, abs=2e-3)
    errs = [abs(impl.ring_inertia_polar(0.5, 1.0, 1.0, 1.0, n)[1] - impl.ring_inertia_closed_form(0.5, 1.0, 1.0, 1.0)[1]) for n in (4, 8, 16)]
    assert np.log2(errs[0] / errs[1]) > 1.9 and np.log2(errs[1] / errs[2]) > 1.9
    with pytest.raises(ValueError):
        impl.ring_inertia_polar(1.0, 0.5, 1.0, 1.0, 8)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "wing_box_mass.png").stat().st_size > 1000
