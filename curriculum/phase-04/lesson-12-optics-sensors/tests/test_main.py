"""Public tests for P04.12 Optics & Sensors. They assert physics, not magic numbers: Snell's law, its
reversibility and total internal reflection beyond the critical angle; the thin-lens equation and its
limits (object at infinity, at 2f, inside f); pixel scale, field of view and their small-angle limit;
the Airy radius scaling with lambda N; star counts and photon counts scaling as they must; the
pixel-integrated Gaussian conserving flux; an unbiased centroid for a centred spot; photon-limited
centroid error falling as 1/sqrt(N) in SYNTHETIC images (seeded); pixel locking for a narrow spot;
Bernoulli's pitot formula as the low-Mach limit of the compressible one; the Wheatstone bridge's
balance, its linear limit and the full-bridge factor 4; the thermistor model's inversion; the
accelerometer's 1/omega^2; and ADC resolution in physical units.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-04/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

ARCSEC = np.pi / 648000


# ---------------------------------------------------------------- refraction

def test_snells_law_and_reversibility():
    th1 = np.radians(np.array([0.0, 10.0, 30.0, 60.0, 89.0]))
    th2 = impl.refraction_angle(1.0, 1.5, th1)
    npt.assert_allclose(1.0 * np.sin(th1), 1.5 * np.sin(th2), atol=1e-14)
    assert th2[0] == 0.0                                                          # normal incidence undeviated
    npt.assert_allclose(impl.refraction_angle(1.5, 1.0, th2), th1, atol=1e-12)    # the path is reversible
    assert np.all(th2 <= th1)                                                     # bends toward the normal


def test_total_internal_reflection_beyond_the_critical_angle():
    tc = impl.critical_angle(1.5, 1.0)
    assert np.sin(tc) == pytest.approx(1 / 1.5)
    assert np.degrees(tc) == pytest.approx(41.81, abs=0.01)
    assert np.isfinite(impl.refraction_angle(1.5, 1.0, tc - 1e-6))
    assert np.isnan(impl.refraction_angle(1.5, 1.0, tc + 1e-6))
    assert impl.refraction_angle(1.5, 1.0, tc - 1e-9) == pytest.approx(np.pi / 2, abs=1e-3)   # grazes the surface
    with pytest.raises(ValueError):
        impl.critical_angle(1.0, 1.5)
    assert impl.fresnel_reflectance_normal(1.0, 1.5) == pytest.approx(0.04)
    assert impl.fresnel_reflectance_normal(1.5, 1.0) == pytest.approx(impl.fresnel_reflectance_normal(1.0, 1.5))
    assert impl.fresnel_reflectance_normal(1.5, 1.5) == 0.0


# ---------------------------------------------------------------- thin lens and camera

def test_thin_lens_equation_and_limits():
    f = 0.025
    d_o = np.array([0.03, 0.05, 1.0, 100.0])
    d_i = impl.thin_lens_image_distance(f, d_o)
    npt.assert_allclose(1 / d_o + 1 / d_i, 1 / f, rtol=1e-12)
    assert impl.thin_lens_image_distance(f, np.inf) == pytest.approx(f)              # stars focus at f
    assert impl.thin_lens_image_distance(f, 2 * f) == pytest.approx(2 * f)          # 2f -> 2f
    assert impl.magnification(2 * f, 2 * f) == pytest.approx(-1.0)
    assert impl.thin_lens_image_distance(f, 0.5 * f) < 0                             # virtual image inside f
    assert d_i[-1] == pytest.approx(f, rel=1e-3)


def test_pixel_scale_fov_and_f_number():
    assert impl.f_number(0.05, 0.025) == pytest.approx(2.0)
    s = impl.pixel_scale(5.5e-6, 0.025)
    assert s == pytest.approx(2.2e-4)
    fov = impl.field_of_view(1024, 5.5e-6, 0.025)
    assert fov == pytest.approx(2 * np.arctan(1024 * 5.5e-6 / 0.05))
    assert fov < 1024 * s                                                            # arctan is below its argument
    assert impl.field_of_view(10, 5.5e-6, 0.025) == pytest.approx(10 * s, rel=1e-6)  # small-angle limit
    # the solid angle of a small square field tends to fov^2; a hemisphere-sized one to 2 pi
    assert impl.fov_solid_angle(0.01, 0.01) == pytest.approx(1e-4, rel=1e-4)
    assert impl.fov_solid_angle(np.pi, np.pi) == pytest.approx(2 * np.pi)


def test_airy_radius_scales_with_wavelength_and_f_number():
    r = impl.airy_radius(600e-9, 1.2)
    assert r == pytest.approx(1.22 * 600e-9 * 1.2)
    assert impl.airy_radius(600e-9, 2.4) == pytest.approx(2 * r)
    assert impl.airy_radius(1200e-9, 1.2) == pytest.approx(2 * r)
    assert impl.airy_gaussian_sigma(600e-9, 1.2) < r
    st = impl.StarTracker()
    assert r < st.pitch                                                              # the focused star is sub-pixel


# ---------------------------------------------------------------- stars and photons

def test_star_and_photon_counts_scale_correctly():
    assert impl.stars_brighter_than(7.0) / impl.stars_brighter_than(6.0) == pytest.approx(np.exp(impl.STAR_COUNT_B))
    Om = impl.fov_solid_angle(0.2, 0.2)
    assert impl.mean_stars_in_fov(6.0, 2 * Om) == pytest.approx(2 * impl.mean_stars_in_fov(6.0, Om))
    assert impl.mean_stars_in_fov(6.0, 4 * np.pi) == pytest.approx(impl.stars_brighter_than(6.0))
    e5 = impl.star_electrons(5.0, 0.02, 0.1)
    assert impl.star_electrons(0.0, 0.02, 0.1) / e5 == pytest.approx(100.0)          # 5 magnitudes = factor 100
    assert impl.star_electrons(5.0, 0.04, 0.1) == pytest.approx(4 * e5)               # aperture area
    assert impl.star_electrons(5.0, 0.02, 0.2) == pytest.approx(2 * e5)               # exposure


# ---------------------------------------------------------------- centroiding

def test_pixel_integrated_spot_conserves_flux_and_centred_centroid_is_exact():
    img = impl.gaussian_spot(0.0, 0.0, 1.0, 1000.0, 11)
    assert img.shape == (11, 11)
    assert img.sum() == pytest.approx(1000.0, rel=1e-6)                             # 5.5 sigma half-width
    npt.assert_allclose(img, img.T, atol=1e-10)
    npt.assert_allclose(img, img[::-1, ::-1], atol=1e-10)
    assert impl.centroid(img) == pytest.approx((0.0, 0.0), abs=1e-12)
    # a wide spot is centroided without bias at any sub-pixel position
    for x0, y0 in ((0.3, -0.2), (-0.45, 0.1)):
        cx, cy = impl.centroid(impl.gaussian_spot(x0, y0, 1.0, 1.0, 11))
        assert cx == pytest.approx(x0, abs=1e-4) and cy == pytest.approx(y0, abs=1e-4)
    # background subtraction
    assert impl.centroid(impl.gaussian_spot(0.2, 0.0, 1.0, 100.0, 11) + 5.0, background=5.0)[0] == pytest.approx(0.2, abs=1e-4)


def test_photon_limited_centroid_error_falls_as_one_over_sqrt_n():
    fluxes = np.array([200.0, 800.0, 3200.0])
    rms = np.array([np.std(impl.simulate_centroid_errors(F, 1.0, 0.0, 600, size=9, seed=21)) for F in fluxes])
    slope = np.polyfit(np.log(fluxes), np.log(rms), 1)[0]
    assert slope == pytest.approx(-0.5, abs=0.07)                     # 600 trials: rms known to about 3 % each
    # level: sqrt(sigma^2 + 1/12) / sqrt(N) (the pixel adds its own 1/12 px^2 variance)
    pred = impl.centroid_sigma_photon_limit(np.sqrt(1.0 + 1.0 / 12.0), fluxes)
    npt.assert_allclose(rms, pred, rtol=0.12)
    errs = impl.simulate_centroid_errors(800.0, 1.0, 0.0, 600, size=9, seed=22)
    assert abs(errs.mean()) < 4 * errs.std() / np.sqrt(errs.size)     # unbiased within 4 standard errors
    # reproducible with the seed
    npt.assert_array_equal(impl.simulate_centroid_errors(500.0, 1.0, 2.0, 20, seed=5),
                           impl.simulate_centroid_errors(500.0, 1.0, 2.0, 20, seed=5))


def test_read_noise_penalises_large_windows():
    e5 = impl.simulate_centroid_errors(543.0, 1.0, 10.0, 600, size=5, seed=31)
    e9 = impl.simulate_centroid_errors(543.0, 1.0, 10.0, 600, size=9, seed=31)
    assert np.sqrt(np.mean(e9 ** 2)) > 1.5 * np.sqrt(np.mean(e5 ** 2))


def test_pixel_locking_is_large_for_a_sharp_spot_and_vanishes_for_a_wide_one():
    x0, b_sharp = impl.pixel_locking_bias(0.2)
    _, b_wide = impl.pixel_locking_bias(1.0)
    assert np.max(np.abs(b_sharp)) > 0.1
    assert np.max(np.abs(b_wide)) < 1e-3
    npt.assert_allclose(b_sharp, -b_sharp[::-1], atol=1e-12)          # odd in x0
    assert b_sharp[np.argmin(np.abs(x0))] == pytest.approx(0.0, abs=1e-12)
    # attitude: n independent stars average down as 1/sqrt(n)
    a1 = impl.attitude_error_cross_boresight(0.1, 2.2e-4, 1)
    assert impl.attitude_error_cross_boresight(0.1, 2.2e-4, 16) == pytest.approx(a1 / 4)


# ---------------------------------------------------------------- sensors

def test_pitot_incompressible_is_the_low_mach_limit_of_the_compressible_formula():
    assert impl.pitot_airspeed(0.5 * 1.225 * 50.0 ** 2) == pytest.approx(50.0)
    a0 = np.sqrt(1.4 * 101325.0 / 1.225)
    for M in (0.05, 0.3, 0.6, 0.9):
        qc = impl.pitot_delta_p_compressible(M)
        assert impl.calibrated_airspeed(qc) == pytest.approx(M * a0, rel=1e-10)    # sea level: CAS = TAS
    V_inc = impl.pitot_airspeed(impl.pitot_delta_p_compressible(0.05))
    assert V_inc == pytest.approx(0.05 * a0, rel=1e-3)                                # same at low Mach
    V_inc6 = impl.pitot_airspeed(impl.pitot_delta_p_compressible(0.6))
    assert V_inc6 > 1.03 * 0.6 * a0                                                   # overreads at Mach 0.6


def test_wheatstone_bridge_balance_linearity_and_full_bridge():
    assert impl.bridge_output(5.0, 350, 350, 350, 350) == pytest.approx(0.0)
    assert impl.bridge_output(5.0, 100, 200, 400, 200) == pytest.approx(0.0)        # R1 R3 = R2 R4
    eps = np.array([-2000e-6, -100e-6, 0.0, 100e-6, 2000e-6])
    q = impl.quarter_bridge(3.3, eps, 2.0)
    R = 350.0
    exact = [impl.bridge_output(3.3, R, R, R * (1 + 2.0 * e), R) for e in eps]
    npt.assert_allclose(q, exact, atol=1e-15)
    npt.assert_allclose(q[1:4], 3.3 * 2.0 * eps[1:4] / 4, rtol=2e-4)                # linear for small strain
    assert abs(q[-1] - 3.3 * 2.0 * eps[-1] / 4) > 1e-6                              # nonlinearity visible at 2000 ue
    npt.assert_allclose(impl.full_bridge(3.3, eps, 2.0), 4 * 3.3 * 2.0 * eps / 4)


def test_thermistor_model_inverts_and_divider_is_monotone():
    T = np.array([253.15, 273.15, 298.15, 333.15])
    R = impl.thermistor_resistance(T)
    assert R[2] == pytest.approx(10e3)
    assert np.all(np.diff(R) < 0)                                                    # NTC
    npt.assert_allclose(impl.thermistor_temperature(R), T, rtol=1e-12)
    v = impl.divider_voltage(3.3, 10e3, R)
    assert v[2] == pytest.approx(1.65)
    assert np.all(np.diff(v) < 0)
    assert impl.divider_voltage(3.3, 1.0, 1e12) == pytest.approx(3.3)


def test_accelerometer_and_adc_resolution():
    x = impl.accelerometer_deflection(9.80665, 5e3)
    assert x == pytest.approx(9.80665 / (2 * np.pi * 5e3) ** 2)
    assert impl.accelerometer_deflection(9.80665, 10e3) == pytest.approx(x / 4)
    assert impl.adc_resolution(32.0, 12) == pytest.approx(32.0 / 4096)
    assert impl.bits_needed(32.0, 1e-3) == 15
    assert impl.adc_resolution(32.0, impl.bits_needed(32.0, 1e-3)) <= 1e-3
    assert impl.adc_resolution(32.0, impl.bits_needed(32.0, 1e-3) - 1) > 1e-3


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "star_tracker.png").exists()
