"""Public tests for P01.11 Vector Fields, Line & Surface Integrals. They assert properties, not magic numbers.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-01/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

GAMMA = 7.5           # m**2/s, an arbitrary counter-clockwise circulation
Q2D = 3.0             # m**2/s, 2D source strength
MU = 3.986004418e14   # m**3/s**2, Earth's gravitational parameter


def vortex(x, y):
    return impl.vortex_velocity(x, y, GAMMA)


def source(x, y):
    return impl.source_velocity(x, y, Q2D)


# ---------------------------------------------------------------- fields

def test_vortex_velocity_is_tangential_with_the_right_magnitude_and_sense():
    th = np.linspace(0.0, 2 * np.pi, 37)
    r = 2.0
    x, y = r * np.cos(th), r * np.sin(th)
    u, v = vortex(x, y)
    npt.assert_allclose(u * x + v * y, 0.0, atol=1e-12)                      # no radial component
    npt.assert_allclose(np.hypot(u, v), GAMMA / (2 * np.pi * r), rtol=1e-12)  # 1/r decay
    assert impl.vortex_velocity(1.0, 0.0, GAMMA)[1] > 0.0                     # counter-clockwise at (1, 0)
    assert impl.vortex_velocity(1.0, 0.0, -GAMMA)[1] < 0.0


def test_source_velocity_is_radial_and_decays_as_one_over_r():
    x, y = np.array([1.0, 0.0, -3.0]), np.array([0.0, 2.0, 4.0])
    u, v = source(x, y)
    npt.assert_allclose(u * y - v * x, 0.0, atol=1e-12)                      # no tangential component
    npt.assert_allclose(np.hypot(u, v) * 2 * np.pi * np.hypot(x, y), Q2D, rtol=1e-12)


def test_uniform_velocity_broadcasts_and_respects_the_angle():
    x = np.zeros((3, 4))
    u, v = impl.uniform_velocity(x, x, 10.0, np.pi / 6)
    assert u.shape == (3, 4) and v.shape == (3, 4)
    npt.assert_allclose(u, 10.0 * np.cos(np.pi / 6))
    npt.assert_allclose(v, 5.0)


def test_cylinder_flow_has_no_normal_velocity_on_the_surface_and_recovers_the_stream_far_away():
    v_inf, a, gamma = 50.0, 0.5, -20.0
    th = np.linspace(0.0, 2 * np.pi, 73)
    u, v = impl.cylinder_flow_velocity(a * np.cos(th), a * np.sin(th), v_inf, a, gamma)
    npt.assert_allclose(u * np.cos(th) + v * np.sin(th), 0.0, atol=1e-10)
    u_far, v_far = impl.cylinder_flow_velocity(np.array([1e4, -1e4]), np.array([1e4, 1e4]), v_inf, a, gamma)
    npt.assert_allclose(u_far, v_inf, rtol=1e-4)
    npt.assert_allclose(v_far, 0.0, atol=1e-3)
    # Without circulation the flow is symmetric top to bottom: same speed at +-y.
    u1, v1 = impl.cylinder_flow_velocity(0.0, 1.0, v_inf, a, 0.0)
    u2, v2 = impl.cylinder_flow_velocity(0.0, -1.0, v_inf, a, 0.0)
    assert np.hypot(u1, v1) == pytest.approx(np.hypot(u2, v2))


def test_rankine_vortex_is_continuous_at_the_core_and_solid_body_inside():
    a = 1.5
    for r in (0.3, 0.9):
        u, v = impl.rankine_vortex_velocity(r, 0.0, GAMMA, a)
        assert v == pytest.approx(GAMMA * r / (2 * np.pi * a**2)) and u == pytest.approx(0.0, abs=1e-15)
    inside = impl.rankine_vortex_velocity(a * (1 - 1e-9), 0.0, GAMMA, a)[1]
    outside = impl.rankine_vortex_velocity(a * (1 + 1e-9), 0.0, GAMMA, a)[1]
    assert inside == pytest.approx(outside, rel=1e-6)
    u0, v0 = impl.rankine_vortex_velocity(0.0, 0.0, GAMMA, a)
    assert np.isfinite(u0) and np.isfinite(v0) and u0 == 0.0 and v0 == 0.0


def test_gravity_field_is_minus_the_gradient_of_the_potential():
    x, y = 6.771e6, 2.0e6
    h = 10.0
    gx = -(impl.gravity_potential(x + h, y, MU) - impl.gravity_potential(x - h, y, MU)) / (2 * h)
    gy = -(impl.gravity_potential(x, y + h, MU) - impl.gravity_potential(x, y - h, MU)) / (2 * h)
    fx, fy = impl.gravity_field(x, y, MU)
    npt.assert_allclose([fx, fy], [gx, gy], rtol=1e-6)
    assert np.hypot(fx, fy) == pytest.approx(MU / (x**2 + y**2))


# ---------------------------------------------------------------- paths

def test_path_helpers_are_counter_clockwise_with_the_right_counts():
    c = impl.circle_path(2.0, 50, (1.0, -1.0))
    assert c.shape == (50, 2)
    npt.assert_allclose(np.hypot(c[:, 0] - 1.0, c[:, 1] + 1.0), 2.0, rtol=1e-12)
    r = impl.rectangle_path(-1.0, 3.0, 0.0, 2.0, 7)
    assert r.shape == (28, 2) and len(np.unique(r, axis=0)) == 28
    for pts in (c, r):   # shoelace area positive <=> counter-clockwise
        x, y = pts[:, 0], pts[:, 1]
        assert 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y) > 0.0
    with pytest.raises(ValueError):
        impl.circle_path(1.0, 2)
    with pytest.raises(ValueError):
        impl.rectangle_path(1.0, 0.0, 0.0, 1.0, 4)


# ---------------------------------------------------------------- circulation and flux

def test_circulation_equals_gamma_for_any_enclosing_path_and_zero_otherwise():
    for path in (impl.circle_path(1.0, 2000), impl.circle_path(4.0, 2000, (0.5, -0.5)),
                 impl.rectangle_path(-1.0, 2.0, -0.5, 1.5, 500)):
        assert impl.circulation(vortex, path) == pytest.approx(GAMMA, rel=1e-5)
    for path in (impl.circle_path(0.5, 2000, (3.0, 0.0)), impl.rectangle_path(1.0, 2.0, 1.0, 3.0, 500)):
        assert impl.circulation(vortex, path) == pytest.approx(0.0, abs=1e-5 * GAMMA)


def test_circulation_is_linear_in_gamma_and_flips_with_orientation():
    path = impl.circle_path(1.0, 500)
    c1 = impl.circulation(lambda x, y: impl.vortex_velocity(x, y, 1.0), path)
    c3 = impl.circulation(lambda x, y: impl.vortex_velocity(x, y, 3.0), path)
    assert c3 == pytest.approx(3 * c1, rel=1e-12)
    assert impl.circulation(vortex, path[::-1]) == pytest.approx(-impl.circulation(vortex, path), rel=1e-12)


def test_circulation_converges_at_second_order_on_circles_and_squares():
    ns = [32, 64, 128, 256]
    errs = impl.circulation_convergence(vortex, GAMMA, ns)
    assert errs.shape == (4, 2) and np.all(errs > 0.0)
    rates = np.log2(errs[:-1] / errs[1:])
    assert np.all(rates > 1.8), f"expected second-order convergence, got {rates}"
    # The circle error has a closed form: n tan(pi/n)/pi - 1 times Gamma.
    n = np.array(ns, dtype=float)
    npt.assert_allclose(errs[:, 0], GAMMA * (n * np.tan(np.pi / n) / np.pi - 1.0), rtol=1e-8)


def test_source_has_flux_q_and_no_circulation_vortex_has_circulation_and_no_flux():
    enclosing = impl.rectangle_path(-1.0, 1.5, -2.0, 0.5, 400)
    outside = impl.circle_path(0.4, 400, (2.0, 2.0))
    assert impl.flux_2d(source, enclosing) == pytest.approx(Q2D, rel=1e-5)
    assert impl.flux_2d(source, outside) == pytest.approx(0.0, abs=1e-9)
    assert impl.circulation(source, enclosing) == pytest.approx(0.0, abs=1e-5)
    assert impl.flux_2d(vortex, enclosing) == pytest.approx(0.0, abs=1e-5)
    assert impl.flux_2d(source, impl.circle_path(0.1, 400)) == pytest.approx(Q2D, rel=1e-4)   # any radius


def test_uniform_stream_has_zero_circulation_and_zero_net_flux():
    stream = lambda x, y: impl.uniform_velocity(x, y, 30.0, 0.3)  # noqa: E731
    path = impl.rectangle_path(0.0, 2.0, 0.0, 1.0, 10)
    assert impl.circulation(stream, path) == pytest.approx(0.0, abs=1e-12)
    assert impl.flux_2d(stream, path) == pytest.approx(0.0, abs=1e-12)


def test_cylinder_circulation_is_gamma_on_every_path_outside_the_body():
    v_inf, a, gamma = 50.0, 0.5, -20.0
    cyl = lambda x, y: impl.cylinder_flow_velocity(x, y, v_inf, a, gamma)  # noqa: E731
    for path in (impl.circle_path(0.6, 4000), impl.circle_path(3.0, 1000), impl.rectangle_path(-2, 2, -1, 1, 500)):
        assert impl.circulation(cyl, path) == pytest.approx(gamma, rel=1e-4)
    assert impl.circulation(cyl, impl.circle_path(0.2, 1000, (3.0, 0.0))) == pytest.approx(0.0, abs=1e-6)
    assert impl.flux_2d(cyl, impl.circle_path(0.6, 4000)) == pytest.approx(0.0, abs=1e-6)   # no mass crosses


def test_rankine_circulation_grows_with_enclosed_vorticity_then_saturates():
    a = 1.5
    rank = lambda x, y: impl.rankine_vortex_velocity(x, y, GAMMA, a)  # noqa: E731
    for r in (0.5, 1.0):
        assert impl.circulation(rank, impl.circle_path(r, 2000)) == pytest.approx(GAMMA * (r / a) ** 2, rel=1e-5)
    for r in (1.5, 4.0):
        assert impl.circulation(rank, impl.circle_path(r, 2000)) == pytest.approx(GAMMA, rel=1e-5)


def test_kutta_joukowski_sign_and_magnitude():
    rho, v_inf, chord, cl = 1.225, 50.0, 1.0, 0.8
    gamma_cw = -0.5 * cl * chord * v_inf             # clockwise circulation for an airfoil lifting upward
    lift = impl.kutta_joukowski_lift(rho, v_inf, gamma_cw)
    assert lift == pytest.approx(0.5 * rho * v_inf**2 * chord * cl)   # q c c_l
    assert impl.kutta_joukowski_lift(rho, v_inf, -gamma_cw) == pytest.approx(-lift)


# ---------------------------------------------------------------- conservative fields and work

def test_gravity_work_is_path_independent_and_equals_the_potential_drop():
    grav = lambda x, y: impl.gravity_field(x, y, MU)  # noqa: E731
    p1, p2 = np.array([6.771e6, 0.0]), np.array([0.0, 7.171e6])
    straight = np.linspace(p1, p2, 4001)
    th = np.linspace(0.0, np.pi / 2, 4001)
    rr = np.linspace(p1[0], p2[1], 4001)
    arc = np.column_stack([rr * np.cos(th), rr * np.sin(th)])
    w_straight, w_arc = impl.work_along_path(grav, straight), impl.work_along_path(grav, arc)
    expected = impl.gravity_potential(*p1, MU) - impl.gravity_potential(*p2, MU)
    assert w_straight == pytest.approx(expected, rel=1e-6)
    assert w_arc == pytest.approx(expected, rel=1e-6)
    assert w_straight < 0.0                                   # climbing against gravity costs energy
    loop = impl.circle_path(7.0e6, 3000, (0.5e6, 0.0))
    assert impl.circulation(grav, loop) == pytest.approx(0.0, abs=1e-6 * abs(expected))


def test_open_line_integral_of_a_gradient_field_is_the_potential_difference():
    phi = lambda x, y: x**2 * y + np.sin(y)  # noqa: E731
    grad = lambda x, y: (2 * x * y, x**2 + np.cos(y))  # noqa: E731
    t = np.linspace(0.0, 1.0, 2001)
    wiggle = np.column_stack([t * 2.0, np.sin(3 * t) + t])      # an arbitrary curve from (0,0) to (2, sin3 + 1)
    assert impl.line_integral(grad, wiggle, closed=False) == pytest.approx(phi(2.0, np.sin(3) + 1.0) - phi(0.0, 0.0), rel=1e-6)
    with pytest.raises(ValueError):
        impl.line_integral(grad, np.zeros((3, 3)))


# ---------------------------------------------------------------- surface integrals

def test_surface_normal_of_the_sphere_points_outward_and_has_area_r2_sin():
    R = 2.0
    sph = impl.sphere_surface(R)
    s, t = np.array([0.3, 1.2, 2.5]), np.array([0.1, 2.0, 4.0])
    n = impl.surface_normal(sph, s, t)
    x, y, z = sph(s, t)
    assert np.all(n[:, 0] * x + n[:, 1] * y + n[:, 2] * z > 0.0)
    npt.assert_allclose(np.linalg.norm(n, axis=1), R**2 * np.sin(s), rtol=1e-6)


def test_source_flux_through_a_sphere_is_q_for_any_radius_and_centre_offset():
    q = 2.0
    src = lambda x, y, z: impl.point_source_3d(x, y, z, q)  # noqa: E731
    for R, c in ((0.5, (0, 0, 0)), (3.0, (0, 0, 0)), (2.0, (0.5, -0.3, 0.2))):
        flux = impl.surface_flux(src, impl.sphere_surface(R, c), (0.0, np.pi), (0.0, 2 * np.pi), 80, 80)
        assert flux == pytest.approx(q, rel=2e-4)
    far = impl.surface_flux(src, impl.sphere_surface(0.5, (5.0, 0.0, 0.0)), (0.0, np.pi), (0.0, 2 * np.pi), 80, 80)
    assert far == pytest.approx(0.0, abs=1e-6)


def test_uniform_flow_has_zero_net_flux_through_a_closed_sphere_and_rho_v_a_through_a_disc():
    rho, v, R = 0.4135, 250.0, 1.0
    stream = lambda x, y, z: (np.zeros_like(x), np.zeros_like(x), np.full_like(x, rho * v))  # noqa: E731
    closed = impl.surface_flux(stream, impl.sphere_surface(1.7), (0.0, np.pi), (0.0, 2 * np.pi), 40, 40)
    assert closed == pytest.approx(0.0, abs=1e-6 * rho * v)
    mdot = impl.surface_flux(stream, impl.disc_surface(R), (0.0, R), (0.0, 2 * np.pi), 40, 40)
    assert mdot == pytest.approx(impl.mass_flow_rate(rho, v, R), rel=1e-6)
    assert impl.mass_flow_rate(rho, v, R) == pytest.approx(rho * v * np.pi * R**2)
    with pytest.raises(ValueError):
        impl.mass_flow_rate(-1.0, v, R)


def test_surface_flux_converges_at_second_order():
    pos = lambda x, y, z: (x, y, z)  # noqa: E731
    R = 1.3
    exact = 4 * np.pi * R**3                                   # position field through a sphere
    errs = np.array([abs(impl.surface_flux(pos, impl.sphere_surface(R), (0.0, np.pi), (0.0, 2 * np.pi), n, n) - exact)
                     for n in (8, 16, 32)])
    rates = np.log2(errs[:-1] / errs[1:])
    assert np.all(rates > 1.8), f"expected second-order convergence, got {rates}"
    with pytest.raises(ValueError):
        impl.surface_flux(pos, impl.sphere_surface(R), (0.0, np.pi), (0.0, 2 * np.pi), 0, 4)
