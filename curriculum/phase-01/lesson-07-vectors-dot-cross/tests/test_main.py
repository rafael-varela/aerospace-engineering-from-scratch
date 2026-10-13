"""Public tests for P01.07 Vectors, Dot & Cross Products. They assert identities and physics, not magic numbers.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-01/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

MU = 3.986004418e14
RNG = np.random.default_rng(7)
A = RNG.normal(size=(20, 3))
B = RNG.normal(size=(20, 3))
C = RNG.normal(size=(20, 3))
R_SAT = np.array([3500.0e3, -5000.0e3, 3000.0e3])     # m
V_SAT = np.array([6.020e3, 4.670e3, 0.730e3])         # m/s


# ---------------------------------------------------------------- against NumPy as the oracle

def test_dot_cross_norm_match_numpy_on_single_and_stacked_vectors():
    for a, b in zip(A, B):
        assert impl.dot(a, b) == pytest.approx(np.dot(a, b), rel=1e-14)
        npt.assert_allclose(impl.cross(a, b), np.cross(a, b), rtol=1e-14, atol=1e-15)
        assert impl.norm(a) == pytest.approx(np.linalg.norm(a), rel=1e-14)
    npt.assert_allclose(impl.dot(A, B), np.einsum("ij,ij->i", A, B), rtol=1e-14)
    npt.assert_allclose(impl.cross(A, B), np.cross(A, B), rtol=1e-14, atol=1e-15)
    npt.assert_allclose(impl.norm(A), np.linalg.norm(A, axis=1), rtol=1e-14)
    assert impl.cross(A, B).shape == (20, 3) and np.shape(impl.dot(A, B)) == (20,)
    assert isinstance(impl.dot(A[0], B[0]), float) and impl.cross(A[0], B[0]).shape == (3,)


def test_shape_validation():
    with pytest.raises(ValueError):
        impl.dot(np.ones(2), np.ones(2))
    with pytest.raises(ValueError):
        impl.cross(np.ones(4), np.ones(4))
    with pytest.raises(ValueError):
        impl.unit(np.zeros(3))
    with pytest.raises(ValueError):
        impl.angle_between(np.zeros(3), np.ones(3))


# ---------------------------------------------------------------- algebraic identities

def test_dot_is_symmetric_bilinear_and_gives_the_norm():
    npt.assert_allclose(impl.dot(A, B), impl.dot(B, A), rtol=1e-14)
    npt.assert_allclose(impl.dot(2.0 * A + C, B), 2.0 * impl.dot(A, B) + impl.dot(C, B), rtol=1e-12)
    npt.assert_allclose(np.sqrt(impl.dot(A, A)), impl.norm(A), rtol=1e-14)
    e1, e2 = np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])
    assert impl.dot(e1, e2) == 0.0 and impl.dot(e1, e1) == 1.0


def test_cross_is_antisymmetric_perpendicular_and_right_handed():
    npt.assert_allclose(impl.cross(A, B), -impl.cross(B, A), rtol=1e-14)
    npt.assert_allclose(impl.cross(A, A), 0.0, atol=1e-15)
    c = impl.cross(A, B)
    npt.assert_allclose(impl.dot(c, A), 0.0, atol=1e-12)                       # a x b is perpendicular to a ...
    npt.assert_allclose(impl.dot(c, B), 0.0, atol=1e-12)                       # ... and to b
    e1, e2, e3 = np.eye(3)
    npt.assert_allclose(impl.cross(e1, e2), e3)                                 # i x j = k
    npt.assert_allclose(impl.cross(e2, e3), e1)                                 # j x k = i
    npt.assert_allclose(impl.cross(e3, e1), e2)                                 # k x i = j


def test_lagrange_identity_and_angle():
    """|a x b|**2 + (a . b)**2 = |a|**2 |b|**2, i.e. sin**2 + cos**2 = 1 for the angle between them."""
    lhs = impl.norm(impl.cross(A, B)) ** 2 + impl.dot(A, B) ** 2
    npt.assert_allclose(lhs, impl.norm(A) ** 2 * impl.norm(B) ** 2, rtol=1e-12)
    th = impl.angle_between(A, B)
    npt.assert_allclose(impl.dot(A, B), impl.norm(A) * impl.norm(B) * np.cos(th), rtol=1e-12)
    npt.assert_allclose(impl.norm(impl.cross(A, B)), impl.norm(A) * impl.norm(B) * np.sin(th), rtol=1e-12)
    assert np.all(th >= 0.0) and np.all(th <= np.pi)


def test_angle_between_is_robust_for_nearly_parallel_vectors():
    a = np.array([1.0, 0.0, 0.0])
    tiny = 1e-9
    b = np.array([1.0, tiny, 0.0])
    assert impl.angle_between(a, b) == pytest.approx(tiny, rel=1e-6)            # arccos would return ~1e-8 or 0
    assert impl.angle_between(a, -a) == pytest.approx(np.pi)
    assert impl.angle_between(a, 3.0 * a) == 0.0
    assert impl.angle_between(a, np.array([0.0, 2.0, 0.0])) == pytest.approx(np.pi / 2)


def test_scalar_triple_product_is_the_determinant_and_cyclic():
    det = np.array([np.linalg.det(np.vstack([a, b, c])) for a, b, c in zip(A, B, C)])
    npt.assert_allclose(impl.scalar_triple(A, B, C), det, rtol=1e-12)
    npt.assert_allclose(impl.scalar_triple(A, B, C), impl.scalar_triple(B, C, A), rtol=1e-12)
    npt.assert_allclose(impl.scalar_triple(A, B, C), -impl.scalar_triple(B, A, C), rtol=1e-12)
    npt.assert_allclose(impl.scalar_triple(A, A, C), 0.0, atol=1e-12)          # coplanar: zero volume


def test_vector_triple_product_bac_cab():
    lhs = impl.cross(A, impl.cross(B, C))
    rhs = B * impl.dot(A, C)[:, None] - C * impl.dot(A, B)[:, None]
    npt.assert_allclose(lhs, rhs, rtol=1e-12, atol=1e-12)
    # Jacobi identity: a x (b x c) + b x (c x a) + c x (a x b) = 0.
    jac = impl.cross(A, impl.cross(B, C)) + impl.cross(B, impl.cross(C, A)) + impl.cross(C, impl.cross(A, B))
    npt.assert_allclose(jac, 0.0, atol=1e-12)


def test_projection_and_rejection_decompose_the_vector():
    p = impl.projection(A, B)
    q = impl.rejection(A, B)
    npt.assert_allclose(p + q, A, rtol=1e-12, atol=1e-14)
    npt.assert_allclose(impl.dot(q, B), 0.0, atol=1e-12)                       # rejection is perpendicular to b
    npt.assert_allclose(impl.norm(impl.cross(p, B)), 0.0, atol=1e-12)          # projection is parallel to b
    npt.assert_allclose(impl.norm(p) ** 2 + impl.norm(q) ** 2, impl.norm(A) ** 2, rtol=1e-12)   # Pythagoras
    npt.assert_allclose(impl.projection(A, B), impl.projection(A, 5.0 * B), rtol=1e-12)         # scale-free


# ---------------------------------------------------------------- mechanics: moments

def test_moment_is_unchanged_by_sliding_the_force_along_its_line():
    r = np.array([2.0, 1.0, 0.5])
    f = np.array([10.0, -3.0, 4.0])
    m = impl.moment(r, f)
    npt.assert_allclose(impl.moment(r + 0.7 * f, f), m, rtol=1e-12)            # principle of transmissibility
    npt.assert_allclose(impl.dot(m, f), 0.0, atol=1e-12)
    npt.assert_allclose(impl.dot(m, r), 0.0, atol=1e-12)
    # Lever arm: |M| = |F| x perpendicular distance from origin to the line of action.
    assert impl.norm(m) == pytest.approx(impl.norm(f) * impl.point_line_distance(np.zeros(3), r, f), rel=1e-12)


def test_resultant_obeys_varignon_and_shifts_correctly():
    pts, fs = RNG.normal(size=(6, 3)), RNG.normal(size=(6, 3))
    F, M = impl.resultant(pts, fs)
    npt.assert_allclose(F, fs.sum(axis=0), rtol=1e-12)
    npt.assert_allclose(M, impl.cross(pts, fs).sum(axis=0), rtol=1e-12)        # Varignon: moment of sum = sum of moments
    about = np.array([0.3, -1.2, 0.8])
    F2, M2 = impl.resultant(pts, fs, about)
    npt.assert_allclose(F2, F, rtol=1e-12)
    npt.assert_allclose(M2, M - impl.cross(about, F), rtol=1e-12, atol=1e-12)  # M_B = M_A - r_BA x F
    with pytest.raises(ValueError):
        impl.resultant(pts, fs[:4])


def test_wing_root_loads_match_the_elliptic_closed_form():
    L, s, sweep = 15.5e3, 8.0, np.radians(25.0)
    y = np.linspace(0.0, s, 2001)
    # The distribution integrates to L/2 over the half-span (trapezoid with a fine grid).
    half = np.trapezoid(impl.elliptic_lift_distribution(y, L, s), y)
    assert half == pytest.approx(L / 2, rel=1e-4)
    F, M = impl.wing_root_loads(L, s, sweep, 400)
    assert F[2] == pytest.approx(L / 2, rel=2e-4) and abs(F[0]) < 1e-9 and abs(F[1]) < 1e-9
    bending = L / 2 * 4 * s / (3 * np.pi)                                       # int y l(y) dy = 2 L s / (3 pi)
    assert M[0] == pytest.approx(bending, rel=5e-4)
    assert M[1] == pytest.approx(-np.tan(sweep) * bending, rel=5e-4)           # torsion from the swept arm
    assert abs(M[2]) < 1e-9
    F0, M0 = impl.wing_root_loads(L, s, 0.0, 400)
    assert abs(M0[1]) < 1e-9                                                   # no sweep: no root torsion
    assert M0[0] == pytest.approx(M[0], rel=1e-12)                             # bending does not depend on sweep
    with pytest.raises(ValueError):
        impl.elliptic_lift_distribution(np.array([1.5 * s]), L, s)


# ---------------------------------------------------------------- mechanics: the satellite state vector

def test_angular_momentum_defines_the_orbital_plane():
    h = impl.specific_angular_momentum(R_SAT, V_SAT)
    npt.assert_allclose(h, np.cross(R_SAT, V_SAT), rtol=1e-14)
    assert impl.dot(h, R_SAT) == pytest.approx(0.0, abs=1e-3 * impl.norm(h))
    assert impl.dot(h, V_SAT) == pytest.approx(0.0, abs=1e-3 * impl.norm(h))
    n, inc = impl.orbital_plane(R_SAT, V_SAT)
    assert impl.norm(n) == pytest.approx(1.0, rel=1e-14)
    assert inc == pytest.approx(np.arccos(n[2]), rel=1e-12) and 0.0 <= inc <= np.pi
    # Equatorial prograde orbit: normal is +z, inclination 0; polar orbit: inclination 90 deg.
    r_eq, v_eq = np.array([7.0e6, 0.0, 0.0]), np.array([0.0, 7.5e3, 0.0])
    n_eq, i_eq = impl.orbital_plane(r_eq, v_eq)
    npt.assert_allclose(n_eq, [0.0, 0.0, 1.0]) and i_eq == pytest.approx(0.0)
    _, i_pol = impl.orbital_plane(r_eq, np.array([0.0, 0.0, 7.5e3]))
    assert i_pol == pytest.approx(np.pi / 2)
    with pytest.raises(ValueError):
        impl.orbital_plane(r_eq, 2.0 * r_eq)


def test_angular_momentum_is_conserved_along_a_conic_from_p01_06():
    """r = p/(1 + e cos th) with v from (12) of P01.06: r x v must be the same vector at every th."""
    p, e = 6797.3e3, 0.05
    h0 = np.sqrt(MU * p)
    hs = []
    for th in np.linspace(0.0, 2 * np.pi, 9):
        r = p / (1 + e * np.cos(th))
        rv = np.array([r * np.cos(th), r * np.sin(th), 0.0])
        vr, vt = h0 / p * e * np.sin(th), h0 / p * (1 + e * np.cos(th))
        vv = vr * np.array([np.cos(th), np.sin(th), 0.0]) + vt * np.array([-np.sin(th), np.cos(th), 0.0])
        hs.append(impl.specific_angular_momentum(rv, vv))
        assert impl.semi_latus_rectum(rv, vv, MU) == pytest.approx(p, rel=1e-12)
    hs = np.array(hs)
    npt.assert_allclose(hs, np.broadcast_to(hs[0], hs.shape), rtol=1e-12, atol=1e-6 * h0)
    npt.assert_allclose(hs[:, 2], h0, rtol=1e-12)


def test_flight_path_angle_from_state_vector():
    r = np.array([7.0e6, 0.0, 0.0])
    v_h = np.array([0.0, 7.5e3, 0.0])
    assert impl.flight_path_angle(r, v_h) == 0.0
    v_up = v_h + np.array([1.0e3, 0.0, 0.0])
    assert impl.flight_path_angle(r, v_up) == pytest.approx(np.arctan2(1.0e3, 7.5e3), rel=1e-12)
    assert impl.flight_path_angle(r, v_h - np.array([1.0e3, 0.0, 0.0])) == pytest.approx(-np.arctan2(1.0e3, 7.5e3))
    # Rotating the whole picture does not change gamma (it depends on r and v only through their angle).
    rng = np.random.default_rng(3)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    assert impl.flight_path_angle(q @ r, q @ v_up) == pytest.approx(impl.flight_path_angle(r, v_up), rel=1e-12)


def test_tangential_and_normal_components():
    v = np.array([3.0, 4.0, 0.0])
    a = np.array([0.0, 0.0, 2.0])                                             # perpendicular: purely normal
    at, an = impl.tangential_normal_components(v, a)
    assert at == pytest.approx(0.0, abs=1e-14) and an == pytest.approx(2.0)
    at, an = impl.tangential_normal_components(v, 2.0 * v)                   # parallel: purely tangential
    assert at == pytest.approx(10.0) and an == pytest.approx(0.0, abs=1e-14)
    a = np.array([1.0, -2.0, 0.7])
    at, an = impl.tangential_normal_components(v, a)
    assert at**2 + an**2 == pytest.approx(impl.dot(a, a), rel=1e-12)          # Pythagoras
    # Circular orbit: gravity is purely normal and v**2 / a_n is the radius.
    r = np.array([7.0e6, 0.0, 0.0])
    vc = np.sqrt(MU / 7.0e6) * np.array([0.0, 1.0, 0.0])
    g = -MU / 7.0e6**3 * r
    at, an = impl.tangential_normal_components(vc, g)
    assert at == pytest.approx(0.0, abs=1e-9) and 7.0e6 * an == pytest.approx(impl.dot(vc, vc), rel=1e-12)
    with pytest.raises(ValueError):
        impl.tangential_normal_components(np.zeros(3), a)


# ---------------------------------------------------------------- lines and planes

def test_point_line_and_plane_distances():
    assert impl.point_line_distance([0.0, 3.0, 4.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]) == pytest.approx(5.0)
    assert impl.point_line_distance([2.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]) == pytest.approx(0.0, abs=1e-14)
    assert impl.point_line_distance([0.0, 1.0, 0.0], [5.0, 0.0, 0.0], [3.0, 0.0, 0.0]) == pytest.approx(1.0)   # direction scale-free
    p1, p2, p3 = np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0])
    n, d = impl.plane_from_points(p1, p2, p3)
    assert impl.norm(n) == pytest.approx(1.0)
    for p in (p1, p2, p3):
        assert impl.point_plane_distance(p, n, d) == pytest.approx(0.0, abs=1e-14)
    assert abs(impl.point_plane_distance(np.zeros(3), n, d)) == pytest.approx(1.0 / np.sqrt(3.0))
    assert impl.point_plane_distance(np.zeros(3), n, d) == -impl.point_plane_distance(np.zeros(3), -n, -d)
    with pytest.raises(ValueError):
        impl.plane_from_points(p1, 2.0 * p1, 3.0 * p1)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "vector_mechanics.png").stat().st_size > 1000
