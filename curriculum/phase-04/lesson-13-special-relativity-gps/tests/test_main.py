"""Public tests for P04.13 Special Relativity for GPS. They assert the non-relativistic limits, the
invariance of the spacetime interval, the group property of boosts, the Doppler symmetries, the
energy-momentum relation, the weak-field gravitational rate, and the GPS clock budget (velocity term
about -7 us/day, gravitational term about +46 us/day, net about +38.6 us/day) against an independent
oracle, plus the eccentricity, Sagnac and accumulated-error terms.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-04/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)
C = impl.C
DAY = 86400.0


# ---------------------------------------------------------------- Lorentz factor and its limits

def test_gamma_values_and_monotonicity():
    assert impl.gamma(0.0) == 1.0
    assert impl.gamma(0.6 * C) == pytest.approx(1.25)
    assert impl.gamma(0.8 * C) == pytest.approx(5.0 / 3.0)
    v = np.linspace(0.0, 0.99 * C, 50)
    assert np.all(np.diff(impl.gamma(v)) > 0)
    for bad in (C, 1.5 * C, -C):
        with pytest.raises(ValueError):
            impl.gamma(bad)


def test_non_relativistic_limit_of_gamma_and_time_dilation():
    # gamma - 1 = v^2/(2c^2) + 3 v^4/(8 c^4) + ...; at v = 4 km/s the next term is 1e-21, so the agreement
    # must be to the v^4 term, which the cancellation-free form reaches while gamma(v) - 1 would not.
    v = 3873.8
    b2 = (v / C) ** 2
    assert impl.gamma_minus_one(v) == pytest.approx(b2 / 2 + 3 * b2 ** 2 / 8, rel=1e-12)
    assert impl.time_dilation_rate(v) == pytest.approx(-b2 / 2, rel=1e-9)
    assert impl.time_dilation_rate(v) < 0
    # an airliner at 250 m/s loses about 30 ns per day
    assert impl.time_dilation_rate(250.0) * DAY == pytest.approx(-30.0e-9, rel=0.01)


def test_length_contraction_limits():
    assert impl.length_contraction(2.0, 0.0) == 2.0
    assert impl.length_contraction(2.0, 0.6 * C) == pytest.approx(1.6)
    assert impl.length_contraction(2.0, 3873.8) < 2.0


# ---------------------------------------------------------------- Lorentz transformation

def test_interval_is_invariant_under_boosts():
    rng = np.random.default_rng(413)
    t = rng.uniform(-1.0, 1.0, 20)
    x = rng.uniform(-1.0, 1.0, 20) * C
    for v in rng.uniform(-0.9, 0.9, 5) * C:
        tp, xp = impl.lorentz_transform(t, x, v)
        npt.assert_allclose(impl.spacetime_interval(tp, xp), impl.spacetime_interval(t, x), rtol=1e-9, atol=1e-3)
        t2, x2 = impl.lorentz_inverse(tp, xp, v)
        npt.assert_allclose(t2, t, rtol=1e-12, atol=1e-12)
        npt.assert_allclose(x2, x, rtol=1e-12, atol=1e-6)


def test_boost_at_zero_speed_is_identity_and_light_stays_light():
    t, x = np.array([0.0, 1.0, 2.0]), np.array([0.0, C, 2 * C])      # a light ray x = c t
    tp, xp = impl.lorentz_transform(t, x, 0.0)
    npt.assert_array_equal(tp, t); npt.assert_array_equal(xp, x)
    tp, xp = impl.lorentz_transform(t, x, 0.5 * C)
    npt.assert_allclose(xp, C * tp, rtol=1e-12, atol=1e-6)


def test_two_collinear_boosts_compose_by_velocity_addition():
    t, x = 1.3, 0.4 * C
    v1, v2 = 0.3 * C, 0.5 * C
    t1, x1 = impl.lorentz_transform(t, x, v1)
    t12, x12 = impl.lorentz_transform(t1, x1, v2)
    tc, xc = impl.lorentz_transform(t, x, impl.velocity_addition(v1, v2))
    assert t12 == pytest.approx(tc, rel=1e-12)
    assert x12 == pytest.approx(xc, rel=1e-12)


def test_velocity_addition_limits():
    assert impl.velocity_addition(0.5 * C, 0.5 * C) == pytest.approx(0.8 * C)
    assert impl.velocity_addition(C, 0.9 * C) == pytest.approx(C)                # c plus anything is c
    assert impl.velocity_addition(100.0, 200.0) == pytest.approx(300.0, rel=1e-12)  # Galilean limit
    assert abs(impl.velocity_addition(0.99 * C, 0.99 * C)) < C


# ---------------------------------------------------------------- Doppler

def test_doppler_symmetries_and_non_relativistic_limit():
    f0, v = 1575.42e6, 0.3 * C
    fa, fr = impl.relativistic_doppler(f0, v, 0.0), impl.relativistic_doppler(f0, v, np.pi)
    assert fa > f0 > fr
    assert fa * fr == pytest.approx(f0 ** 2, rel=1e-12)                      # approaching x receding = f0^2
    assert impl.relativistic_doppler(f0, v, np.pi / 2) == pytest.approx(f0 / impl.gamma(v))   # transverse
    assert impl.relativistic_doppler(f0, 800.0, 0.0) == pytest.approx(f0 * (1 + 800.0 / C), rel=1e-11)
    assert impl.relativistic_doppler(f0, 0.0, 0.7) == f0


def test_pass_doppler_is_antisymmetric_to_first_order_with_a_constant_second_order_offset():
    f0, vr = 1575.42e6, np.array([-800.0, 0.0, 800.0])
    f = impl.pass_doppler(f0, vr, 3873.8)
    assert f[2] - f0 == pytest.approx(4204.0, abs=1.0)
    assert (f[2] - f[1]) == pytest.approx(-(f[0] - f[1]), rel=1e-9)
    assert f[1] - f0 == pytest.approx(-f0 * 3873.8 ** 2 / (2 * C ** 2), rel=1e-4)   # -0.13 Hz (round-off: f0 is 1.6e9)


# ---------------------------------------------------------------- dynamics

def test_energy_momentum_relation_and_classical_limit():
    m = 1.0
    for v in (0.0, 1e3, 0.5 * C, 0.99 * C):
        E, p = impl.relativistic_energy(m, v), impl.momentum(m, v)
        assert E ** 2 - (p * C) ** 2 == pytest.approx((m * C ** 2) ** 2, rel=1e-9)
        assert impl.kinetic_energy(m, v) == pytest.approx(E - m * C ** 2, rel=1e-6, abs=1e-3)
    assert impl.kinetic_energy(1.0, 7700.0) == pytest.approx(0.5 * 7700.0 ** 2, rel=1e-9)
    assert impl.kinetic_energy(1.0, 0.9 * C) > 0.5 * (0.9 * C) ** 2


# ---------------------------------------------------------------- gravity

def test_weak_field_gravitational_rate():
    assert impl.gravitational_rate(impl.R_EARTH) == pytest.approx(-6.9535e-10, rel=1e-3)
    assert impl.gravitational_rate(np.inf) == 0.0
    r = np.array([1.0, 2.0, 4.0]) * impl.R_EARTH
    g = impl.gravitational_rate(r)
    assert np.all(np.diff(g) > 0)                                             # higher is faster
    npt.assert_allclose(g * r, g[0] * r[0])                                   # 1/r law


def test_geoid_potential_is_dominated_by_gm_over_r():
    phi = impl.geoid_potential_ratio()
    assert phi < 0
    assert phi == pytest.approx(-impl.GM_EARTH / (impl.R_EARTH * C ** 2), rel=3e-3)
    assert phi < -impl.GM_EARTH / (impl.R_EARTH * C ** 2)                      # rotation and J2 deepen it


# ---------------------------------------------------------------- the GPS budget

def test_gps_clock_budget_terms_and_net():
    b = impl.gps_clock_budget()
    u = b.us_per_day()
    assert -7.5 < u["velocity"] < -6.9
    assert 45.0 < u["gravity"] < 46.3
    assert 0.0 < u["ground_rotation"] + u["ground_j2"] < 0.3
    assert 38.0 < u["net"] < 39.0                                             # the 38 us/day figure
    assert b.net == pytest.approx(4.4647e-10, rel=2e-4)                       # Ashby's published value


def test_gps_budget_matches_independent_oracle():
    for a in (impl.R_EARTH + 400e3, impl.GPS_SEMI_MAJOR_AXIS, 42164e3):
        assert impl.gps_clock_budget(a).net == pytest.approx(impl.gps_frequency_offset_oracle(a), rel=1e-12)


def test_budget_scales_with_mu_and_vanishes_without_gravity():
    # dimensional sanity: every term is proportional to GM (through v^2 or the potential) except the rotation term
    b1, b2 = impl.gps_clock_budget(mu=impl.GM_EARTH), impl.gps_clock_budget(mu=2 * impl.GM_EARTH)
    assert b2.velocity == pytest.approx(2 * b1.velocity) and b2.gravity == pytest.approx(2 * b1.gravity)
    assert b2.ground_rotation == b1.ground_rotation
    b0 = impl.gps_clock_budget(omega=0.0, J2=0.0)
    assert b0.ground_rotation == 0.0 and b0.ground_j2 == 0.0


def test_crossover_altitude_is_half_an_earth_radius():
    a = impl.crossover_radius()
    assert a == pytest.approx(1.5 * impl.R_EARTH, rel=1e-6)
    b = impl.gps_clock_budget(a, omega=0.0, J2=0.0)
    assert b.velocity + b.gravity == pytest.approx(0.0, abs=1e-18)
    iss, geo = impl.gps_clock_budget(impl.R_EARTH + 400e3), impl.gps_clock_budget(42164e3)
    assert iss.net < 0 < impl.gps_clock_budget().net < geo.net             # ISS clocks lose, GEO gains more than GPS


def test_factory_frequency():
    f = impl.factory_clock_frequency()
    assert f < impl.GPS_CLOCK_HZ
    assert impl.GPS_CLOCK_HZ - f == pytest.approx(4.567e-3, rel=2e-3)        # 10 229 999.995 43 Hz
    assert impl.factory_clock_frequency(10.23e6, 0.0) == 10.23e6


def test_eccentricity_correction_properties():
    E = np.linspace(0, 2 * np.pi, 9)
    npt.assert_array_equal(impl.eccentricity_correction(0.0, E), 0.0)
    d = impl.eccentricity_correction(0.02, E)
    assert np.max(np.abs(d)) == pytest.approx(45.8e-9, rel=0.01)
    assert impl.eccentricity_correction(0.02, 0.0) == pytest.approx(0.0, abs=1e-20)
    assert impl.eccentricity_correction(0.02, np.pi / 2) < 0                 # near perigee the clock lagged
    assert impl.eccentricity_correction(0.04, 1.0) == pytest.approx(2 * impl.eccentricity_correction(0.02, 1.0))
    assert d.mean() == pytest.approx(0.0, abs=1e-20)                         # averages to zero over an orbit


def test_sagnac_term_geometry():
    a, R = impl.GPS_SEMI_MAJOR_AXIS, impl.R_EARTH
    rec = np.array([R, 0.0, 0.0])
    assert impl.sagnac_correction(np.array([a, 0.0, 0.0]), rec) == 0.0        # satellite overhead
    west = np.array([a * np.cos(1.0), -a * np.sin(1.0), 0.0])
    east = np.array([a * np.cos(1.0), a * np.sin(1.0), 0.0])
    assert impl.sagnac_correction(west, rec) > 0                              # receiver moving away
    assert impl.sagnac_correction(east, rec) == pytest.approx(-impl.sagnac_correction(west, rec))
    assert impl.sagnac_correction(west, rec, omega=2 * impl.OMEGA_EARTH) == pytest.approx(2 * impl.sagnac_correction(west, rec))
    assert impl.max_sagnac() == pytest.approx(133.4e-9, rel=0.01)
    theta = np.arccos(R / a)
    horizon = np.array([a * np.cos(theta), -a * np.sin(theta), 0.0])
    assert impl.sagnac_correction(horizon, rec) == pytest.approx(impl.max_sagnac(), rel=1e-9)


def test_position_error_grows_linearly():
    t = np.array([1.0, 60.0, 3600.0, 86400.0])
    e = impl.position_error_without_correction(t)
    npt.assert_allclose(e / t, e[0])
    assert e[-1] == pytest.approx(11.56e3, rel=0.01)                          # 11.6 km per day
    assert impl.position_error_without_correction(10.0, 0.0) == 0.0


def test_orbit_speed_is_keplerian():
    a = impl.GPS_SEMI_MAJOR_AXIS
    v = impl.circular_orbit_speed(a)
    assert v ** 2 * a == pytest.approx(impl.GM_EARTH, rel=1e-12)
    assert 2 * np.pi * a / v == pytest.approx(43082.0, rel=1e-3)             # half a sidereal day


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "gps_relativity.png").exists()
