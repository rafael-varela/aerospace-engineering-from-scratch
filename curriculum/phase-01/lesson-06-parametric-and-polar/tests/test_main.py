"""Public tests for P01.06 Parametric Curves & Polar Coordinates. They assert geometry and physics, not magic numbers.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-01/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

MU = 3.986004418e14                      # m**3/s**2
R_P, R_A = 6678.137e3, 42164.0e3         # geostationary transfer ellipse, metres
E = (R_A - R_P) / (R_A + R_P)              # eccentricity from the apsides
P = R_P * (1.0 + E)                       # semi-latus rectum, metres
A = P / (1.0 - E**2)                      # semi-major axis, metres


# ---------------------------------------------------------------- the polar conic

def test_conic_from_apsides_and_back():
    p, e = impl.conic_from_apsides(R_P, R_A)
    assert p == pytest.approx(P, rel=1e-12) and e == pytest.approx(E, rel=1e-12)
    assert impl.semi_major_axis(P, E) == pytest.approx(A, rel=1e-12)
    rp, ra = impl.periapsis_apoapsis(P, E)
    assert rp == pytest.approx(R_P, rel=1e-12) and ra == pytest.approx(R_A, rel=1e-12)
    assert 0.0 < E < 1.0 and P == pytest.approx(A * (1.0 - E**2), rel=1e-12)
    with pytest.raises(ValueError):
        impl.conic_from_apsides(R_A, R_P)          # periapsis above apoapsis


def test_radius_extremes_and_symmetry():
    assert impl.conic_radius(0.0, P, E) == pytest.approx(R_P, rel=1e-12)
    assert impl.conic_radius(np.pi, P, E) == pytest.approx(R_A, rel=1e-12)
    assert impl.conic_radius(np.pi / 2, P, E) == pytest.approx(P, rel=1e-12)   # r = p at theta = 90 deg
    th = np.linspace(-np.pi, np.pi, 41)
    r = impl.conic_radius(th, P, E)
    npt.assert_allclose(r, r[::-1], rtol=1e-12)                                # symmetric about the apse line
    assert np.all(r >= R_P - 1e-6) and np.all(r <= R_A + 1e-6)


def test_circle_is_the_e_zero_limit():
    th = np.linspace(0.0, 2 * np.pi, 13)
    npt.assert_allclose(impl.conic_radius(th, P, 0.0), P, rtol=1e-14)
    npt.assert_allclose(impl.conic_radius_derivative(th, P, 0.0), 0.0, atol=1e-14)
    npt.assert_allclose(impl.flight_path_angle(th, 0.0), 0.0, atol=1e-14)


def test_radius_derivative_matches_central_difference():
    th = np.array([0.3, 1.0, 2.0, 2.9])
    h = 1e-5
    fd = (impl.conic_radius(th + h, P, E) - impl.conic_radius(th - h, P, E)) / (2 * h)
    npt.assert_allclose(impl.conic_radius_derivative(th, P, E), fd, rtol=1e-7)
    assert impl.conic_radius_derivative(0.0, P, E) == 0.0 and abs(impl.conic_radius_derivative(np.pi, P, E)) < 1e-6


def test_open_conics_have_an_asymptote():
    assert impl.true_anomaly_limit(0.5) == pytest.approx(np.pi)
    assert impl.true_anomaly_limit(1.0) == pytest.approx(np.pi)
    assert impl.true_anomaly_limit(2.0) == pytest.approx(2 * np.pi / 3)     # cos(theta_inf) = -1/2
    with pytest.raises(ValueError):
        impl.conic_radius(np.pi, P, 1.5)                                    # beyond the asymptote
    th = 0.999 * impl.true_anomaly_limit(1.5)
    assert impl.conic_radius(th, P, 1.5) > 100.0 * P                        # radius blows up near it
    assert impl.semi_major_axis(P, 1.5) < 0.0 and np.isinf(impl.semi_major_axis(P, 1.0))
    assert np.isinf(impl.periapsis_apoapsis(P, 1.0)[1])


def test_input_validation():
    for bad in [(-1.0, 0.5), (0.0, 0.5), (P, -0.1), (float("nan"), 0.5)]:
        with pytest.raises(ValueError):
            impl.conic_radius(0.0, *bad)
    with pytest.raises(ValueError):
        impl.orbital_period(P, 1.2)
    with pytest.raises(ValueError):
        impl.simpson(np.sin, 0.0, 1.0, 5)                                   # odd panel count


# ---------------------------------------------------------------- parametric vs polar descriptions

def test_parametric_ellipse_lies_on_the_polar_conic():
    Es = np.linspace(0.0, 2 * np.pi, 25)
    x, y = impl.ellipse_parametric(A, E, Es)
    r = np.hypot(x, y)
    th = np.arctan2(y, x)
    npt.assert_allclose(r, impl.conic_radius(th, P, E), rtol=1e-12)
    # Centre is offset from the focus by a e along -x; semi-axes a and b.
    assert x.max() == pytest.approx(A * (1 - E)) and x.min() == pytest.approx(-A * (1 + E))
    assert y.max() == pytest.approx(A * np.sqrt(1 - E**2), rel=1e-6)
    with pytest.raises(ValueError):
        impl.ellipse_parametric(A, 1.0, 0.0)


def test_orbit_points_round_trip():
    x, y = impl.orbit_points(P, E, 101)
    assert x.shape == (101,) and y.shape == (101,)
    th = np.arctan2(y, x)
    npt.assert_allclose(np.hypot(x, y), impl.conic_radius(th, P, E), rtol=1e-12)
    assert x[50] == pytest.approx(R_P) and abs(y[50]) < 1e-6                 # middle point is periapsis
    x, y = impl.orbit_points(P, 1.5, 51)
    assert np.all(np.isfinite(x)) and np.all(np.isfinite(y))


# ---------------------------------------------------------------- integrals: area, arc length, time

def test_simpson_is_exact_for_cubics_and_fourth_order():
    assert impl.simpson(lambda x: x**3 - 2 * x + 1, 0.0, 2.0, 2) == pytest.approx(4.0 - 4.0 + 2.0, abs=1e-12)
    errs = np.array([abs(impl.simpson(np.exp, 0.0, 1.0, n) - (np.e - 1.0)) for n in (4, 8, 16, 32)])
    rates = np.log2(errs[:-1] / errs[1:])
    assert np.all(rates > 3.8), rates


def test_full_sector_area_is_pi_a_b():
    b = A * np.sqrt(1 - E**2)
    assert impl.sector_area(P, E, -np.pi, np.pi, 400) == pytest.approx(np.pi * A * b, rel=1e-8)
    assert impl.sector_area(P, 0.0, 0.0, 2 * np.pi, 8) == pytest.approx(np.pi * P**2, rel=1e-12)
    # Additivity and sign: A(0, pi) + A(pi, 2 pi) = A(0, 2 pi); reversing the limits flips the sign.
    a1 = impl.sector_area(P, E, 0.0, np.pi, 200)
    a2 = impl.sector_area(P, E, np.pi, 2 * np.pi, 200)
    assert a1 + a2 == pytest.approx(impl.sector_area(P, E, 0.0, 2 * np.pi, 400), rel=1e-10)
    assert impl.sector_area(P, E, np.pi, 0.0, 200) == pytest.approx(-a1)


def test_arc_length_of_circle_and_ellipse_perimeter():
    from scipy.special import ellipe

    assert impl.polar_arc_length(P, 0.0, 0.0, 2 * np.pi, 8) == pytest.approx(2 * np.pi * P, rel=1e-12)
    assert impl.polar_arc_length(P, E, -np.pi, np.pi, 400) == pytest.approx(4 * A * ellipe(E**2), rel=1e-8)
    # Arc length is at least the chord between the end points.
    th0, th1 = 0.2, 1.4
    x, y = impl.polar_to_cartesian(impl.conic_radius(np.array([th0, th1]), P, E), np.array([th0, th1]))
    chord = np.hypot(x[1] - x[0], y[1] - y[0])
    assert impl.polar_arc_length(P, E, th0, th1, 64) >= chord


def test_arc_length_converges_at_fourth_order():
    errs = impl.arc_length_convergence(P, E, [8, 16, 32, 64])
    assert np.all(errs[1:] < errs[:-1])
    rates = np.log2(errs[:-1] / errs[1:])
    assert np.all(rates > 3.6), rates
    assert errs[-1] < 1e-6 * impl.polar_arc_length(P, E, 0.0, np.pi / 2, 64)


def test_keplers_second_and_third_laws():
    h = impl.specific_angular_momentum(P, MU)
    assert h == pytest.approx(np.sqrt(MU * P))
    T = impl.orbital_period(P, E, MU)
    # Time for the full revolution from the swept area equals Kepler's third-law period.
    assert impl.time_of_flight(P, E, 0.0, 2 * np.pi, MU, 400) == pytest.approx(T, rel=1e-8)
    # Equal areas in equal times: the two halves of the orbit take the same time by symmetry ...
    assert impl.time_of_flight(P, E, 0.0, np.pi, MU, 200) == pytest.approx(T / 2, rel=1e-8)
    # ... but the apoapsis half (90 deg to 270 deg) takes much longer than the periapsis half.
    t_peri = impl.time_of_flight(P, E, -np.pi / 2, np.pi / 2, MU, 200)
    t_apo = impl.time_of_flight(P, E, np.pi / 2, 3 * np.pi / 2, MU, 200)
    assert t_apo > 5.0 * t_peri and t_peri + t_apo == pytest.approx(T, rel=1e-8)
    # Circular orbit: T = 2 pi r / v exactly.
    v_circ = np.sqrt(MU / P)
    assert impl.orbital_period(P, 0.0, MU) == pytest.approx(2 * np.pi * P / v_circ, rel=1e-12)


def test_period_scales_as_a_to_the_three_halves():
    T1 = impl.orbital_period(P, E, MU)
    T2 = impl.orbital_period(4.0 * P, E, MU)       # a -> 4 a at fixed e
    assert T2 / T1 == pytest.approx(8.0, rel=1e-12)


# ---------------------------------------------------------------- calculus of the parametric curve

def test_parametric_speed_equals_vis_viva():
    th = np.linspace(-np.pi, np.pi, 37)
    r = impl.conic_radius(th, P, E)
    npt.assert_allclose(impl.orbital_speed(th, P, E, MU), impl.vis_viva_speed(r, A, MU), rtol=1e-12)
    # Parabola: v = sqrt(2 mu / r) (escape speed) everywhere.
    th = np.linspace(-2.0, 2.0, 9)
    r = impl.conic_radius(th, P, 1.0)
    npt.assert_allclose(impl.orbital_speed(th, P, 1.0, MU), np.sqrt(2 * MU / r), rtol=1e-12)
    # Periapsis is the fastest point, apoapsis the slowest, and r v is the same at both (gamma = 0 there).
    vp, va = impl.orbital_speed(0.0, P, E, MU), impl.orbital_speed(np.pi, P, E, MU)
    assert vp > va and R_P * vp == pytest.approx(R_A * va, rel=1e-12) == pytest.approx(impl.specific_angular_momentum(P, MU), rel=1e-12)


def test_speed_from_finite_difference_of_position_in_time():
    """Differentiate x(t), y(t) numerically: theta advances by h/r**2 dt, so dx/dt follows from the chain rule."""
    h = impl.specific_angular_momentum(P, MU)
    th0 = 1.1
    r0 = impl.conic_radius(th0, P, E)
    dt = 0.5
    dth = h / r0**2 * dt                                    # Kepler II, first order in dt
    ths = np.array([th0 - dth, th0 + dth])
    x, y = impl.polar_to_cartesian(impl.conic_radius(ths, P, E), ths)
    v_fd = np.hypot(x[1] - x[0], y[1] - y[0]) / (2 * dt)
    assert v_fd == pytest.approx(impl.orbital_speed(th0, P, E, MU), rel=1e-5)


def test_flight_path_angle_properties():
    assert impl.flight_path_angle(0.0, E) == 0.0 and abs(impl.flight_path_angle(np.pi, E)) < 1e-12
    th = np.linspace(0.05, np.pi - 0.05, 30)
    g = impl.flight_path_angle(th, E)
    assert np.all(g > 0.0)                                   # climbing from periapsis to apoapsis
    npt.assert_allclose(impl.flight_path_angle(-th, E), -g, atol=1e-14)
    # gamma is the angle between the velocity and the local horizontal: tan(gamma) = r' / (r theta').
    tan_g = impl.conic_radius_derivative(th, P, E) / impl.conic_radius(th, P, E)
    npt.assert_allclose(np.tan(g), tan_g, rtol=1e-12)
    # At theta = 90 deg on this ellipse tan(gamma) = e.
    assert np.tan(impl.flight_path_angle(np.pi / 2, E)) == pytest.approx(E)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "conic_orbits.png").stat().st_size > 1000
