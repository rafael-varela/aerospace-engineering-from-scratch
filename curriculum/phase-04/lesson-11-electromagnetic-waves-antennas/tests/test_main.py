"""Public tests for P04.11 Electromagnetic Waves & Antennas. They assert physics, not magic numbers:
c = 1/sqrt(mu0 eps0); B = E/c and I = E^2/(2 eta0) consistently; the inverse-square law, and the
power through any sphere around a Hertzian dipole equal to its radiated power; radiation pressure
doubles on a mirror; the short and half-wave dipole patterns, their symmetry and nulls, directivities
1.5 and 1.641 checked against ``scipy.integrate.quad``; the half-wave pattern equal to the sum of
Hertzian elements; R_rad = 73.1 ohm; A_e = G lambda^2 / (4 pi); polarisation factors; the Friis
equation's 20 dB per decade of distance and frequency scaling; slant range at zenith and horizon;
and the link budget's internal consistency.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-04/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import quad

from conftest import load_impl

impl = load_impl(__file__)

C = 299_792_458.0
ETA0 = 4e-7 * np.pi * C


# ---------------------------------------------------------------- plane waves

def test_wave_speed_from_maxwell_is_c_and_scales_in_a_medium():
    assert impl.wave_speed() == pytest.approx(C, rel=1e-9)
    assert impl.wave_speed(eps=4 * 8.8541878128e-12) == pytest.approx(C / 2, rel=1e-9)   # n = 2


def test_plane_wave_fields_and_intensity_are_consistent():
    E0 = 1000.0
    B0 = impl.b_from_e(E0)
    assert B0 == pytest.approx(E0 / C)
    I = impl.plane_wave_intensity(E0)
    assert I == pytest.approx(E0 * B0 / (2 * 4e-7 * np.pi), rel=1e-9)           # <S> = E0 B0 / (2 mu0)
    assert I == pytest.approx(0.5 * C * 8.8541878128e-12 * E0 ** 2, rel=1e-6)
    assert impl.e_from_intensity(I) == pytest.approx(E0)
    assert impl.plane_wave_intensity(2 * E0) == pytest.approx(4 * I)           # intensity ~ E^2


def test_inverse_square_law():
    r = np.array([1.0, 2.0, 10.0, 1e6])
    I = impl.isotropic_intensity(100.0, r)
    npt.assert_allclose(I * 4 * np.pi * r ** 2, 100.0)
    assert I[0] / I[1] == pytest.approx(4.0)
    assert I[0] / I[2] == pytest.approx(100.0)


def test_power_through_any_sphere_equals_the_radiated_power():
    lam = C / 437e6
    P = impl.hertzian_radiated_power(1.0, 0.05 * lam, lam)
    for r in (5.0, 1e3, 1e6):
        assert impl.power_through_sphere(1.0, 0.05 * lam, lam, r) == pytest.approx(P, rel=1e-6)
    # far field falls as 1/r and as sin(theta)
    E1 = impl.hertzian_far_field(1.0, 0.1, lam, 1e3, np.pi / 2)
    assert impl.hertzian_far_field(1.0, 0.1, lam, 2e3, np.pi / 2) == pytest.approx(E1 / 2)
    assert impl.hertzian_far_field(1.0, 0.1, lam, 1e3, np.pi / 6) == pytest.approx(E1 / 2)
    assert impl.hertzian_radiated_power(2.0, 0.05 * lam, lam) == pytest.approx(4 * P)


def test_radiation_pressure_and_sail():
    I = 1361.0
    assert impl.radiation_pressure(I, 0.0) == pytest.approx(I / C)
    assert impl.radiation_pressure(I, 1.0) == pytest.approx(2 * I / C)
    a1 = impl.sail_acceleration(32.0, 5.0)
    assert a1 == pytest.approx(2 * I * 32.0 / (5.0 * C))
    assert impl.sail_acceleration(32.0, 5.0, r_au=2.0) == pytest.approx(a1 / 4)              # inverse square
    assert impl.sail_acceleration(32.0, 5.0, cone_deg=60.0) == pytest.approx(a1 / 4)         # cos^2
    assert impl.sail_acceleration(64.0, 10.0) == pytest.approx(a1)                           # area / mass


# ---------------------------------------------------------------- dipole patterns

def test_patterns_are_symmetric_with_nulls_on_the_axis_and_maximum_broadside():
    t = np.linspace(0.05, np.pi / 2, 50)
    for f in (impl.hertzian_pattern, impl.halfwave_pattern):
        npt.assert_allclose(f(t), f(np.pi - t), atol=1e-14)                    # symmetric about broadside
        assert f(0.0) == pytest.approx(0.0, abs=1e-12) and f(np.pi) == pytest.approx(0.0, abs=1e-12)
        assert f(np.pi / 2) == pytest.approx(1.0)
        assert np.all(np.diff(f(t)) >= -1e-15)                                 # rises toward broadside
    # the half-wave pattern is narrower than the short dipole's
    assert impl.halfwave_pattern(np.pi / 4) < impl.hertzian_pattern(np.pi / 4)
    # a short dipole of any length tends to sin^2 once normalised
    small = impl.dipole_pattern(t, 0.01) / impl.dipole_pattern(np.pi / 2, 0.01)
    npt.assert_allclose(small, np.sin(t) ** 2, rtol=1e-3)


def test_directivities_match_closed_forms_and_scipy_quad():
    assert impl.directivity(impl.hertzian_pattern) == pytest.approx(1.5, rel=1e-6)
    # independent quad oracle for the half-wave dipole
    U = lambda t: (np.cos(np.pi / 2 * np.cos(t)) / np.sin(t)) ** 2 if np.sin(t) > 1e-12 else 0.0
    P = 2 * np.pi * quad(lambda t: U(t) * np.sin(t), 0, np.pi)[0]
    D_ref = 4 * np.pi / P
    assert D_ref == pytest.approx(1.641, abs=1e-3)
    assert impl.directivity(impl.halfwave_pattern) == pytest.approx(D_ref, rel=1e-6)
    assert impl.directivity_quad(impl.halfwave_pattern) == pytest.approx(D_ref, rel=1e-9)
    assert impl.to_db(D_ref) == pytest.approx(2.15, abs=0.01)


def test_beamwidths():
    assert impl.half_power_beamwidth(impl.hertzian_pattern) == pytest.approx(90.0, abs=0.01)
    assert impl.half_power_beamwidth(impl.halfwave_pattern) == pytest.approx(78.1, abs=0.1)


def test_halfwave_pattern_is_the_sum_of_hertzian_elements():
    t = np.linspace(0.1, np.pi - 0.1, 13)
    for Lw in (0.5, 1.0, 1.25):
        npt.assert_allclose(impl.dipole_pattern_numeric(t, Lw, n_seg=400), impl.dipole_pattern(t, Lw), atol=1e-4 * impl.dipole_pattern(t, Lw).max())
    # midpoint sum converges at second order in the segment length
    e1 = np.max(np.abs(impl.dipole_pattern_numeric(t, 0.5, 50) - impl.dipole_pattern(t, 0.5)))
    e2 = np.max(np.abs(impl.dipole_pattern_numeric(t, 0.5, 100) - impl.dipole_pattern(t, 0.5)))
    assert np.log2(e1 / e2) == pytest.approx(2.0, abs=0.2)


def test_radiation_resistance():
    assert impl.radiation_resistance_dipole(0.5) == pytest.approx(73.08, abs=0.02)
    assert impl.radiation_resistance_short(0.05) == pytest.approx(80 * np.pi ** 2 * 0.0025, rel=1e-3)
    assert impl.radiation_resistance_short(0.1) == pytest.approx(4 * impl.radiation_resistance_short(0.05))
    # consistency: U_max = eta0 I^2 / (8 pi^2) and P = I^2 R / 2 give D = eta0 / (pi R) for the half-wave dipole
    D = impl.directivity(impl.halfwave_pattern)
    assert ETA0 / (np.pi * impl.radiation_resistance_dipole(0.5)) == pytest.approx(D, rel=1e-4)


# ---------------------------------------------------------------- gain, aperture, polarisation

def test_gain_aperture_and_dish():
    lam = 0.1
    assert impl.effective_aperture(1.0, lam) == pytest.approx(lam ** 2 / (4 * np.pi))
    assert impl.gain(1.64, 0.8) == pytest.approx(1.312)
    D = 2.0
    G = impl.dish_gain(D, lam, efficiency=0.6)
    assert impl.effective_aperture(G, lam) == pytest.approx(0.6 * np.pi * D ** 2 / 4)     # aperture = eta x area
    assert impl.dish_gain(D, lam / 2) == pytest.approx(4 * impl.dish_gain(D, lam))
    assert impl.from_db(impl.to_db(123.4)) == pytest.approx(123.4)
    assert impl.to_db(2.0) == pytest.approx(3.0103, abs=1e-4)


def test_polarization_loss_factors():
    assert impl.polarization_loss("linear", "linear", 0.0) == pytest.approx(1.0)
    assert impl.polarization_loss("linear", "linear", 90.0) == pytest.approx(0.0, abs=1e-15)
    assert impl.polarization_loss("linear", "linear", 45.0) == pytest.approx(0.5)
    for ang in (0.0, 37.0, 90.0):
        assert impl.polarization_loss("linear", "rhcp", ang) == pytest.approx(0.5)
    assert impl.polarization_loss("rhcp", "rhcp") == 1.0
    assert impl.polarization_loss("rhcp", "lhcp") == 0.0
    with pytest.raises(ValueError):
        impl.polarization_loss("rhcp", "elliptical")


# ---------------------------------------------------------------- Friis and the link

def test_friis_scaling_and_reciprocity():
    lam = C / 437e6
    P1 = impl.friis_received_power(1.0, 1.64, 25.0, 1e6, lam)
    assert impl.friis_received_power(1.0, 25.0, 1.64, 1e6, lam) == pytest.approx(P1)        # reciprocity
    assert impl.friis_received_power(1.0, 1.64, 25.0, 1e7, lam) == pytest.approx(P1 / 100)  # 20 dB per decade
    # Friis = isotropic intensity x effective aperture of the receiver
    S = impl.isotropic_intensity(1.0 * 1.64, 1e6)
    assert P1 == pytest.approx(S * impl.effective_aperture(25.0, lam))
    assert impl.to_db(impl.free_space_path_loss(1e6, lam)) == pytest.approx(
        20 * np.log10(4 * np.pi * 1e6 / lam))


def test_slant_range_limits():
    h = 500e3
    assert impl.slant_range(h, 90.0) == pytest.approx(h)
    horizon = np.sqrt((impl.R_EARTH + h) ** 2 - impl.R_EARTH ** 2)
    assert impl.slant_range(h, 0.0) == pytest.approx(horizon)
    d = impl.slant_range(h, np.array([5.0, 10.0, 30.0, 60.0, 90.0]))
    assert np.all(np.diff(d) < 0)
    # law of cosines check at 30 degrees: (R+h)^2 = R^2 + d^2 + 2 R d sin(e)
    R, e = impl.R_EARTH, np.radians(30.0)
    dd = impl.slant_range(h, 30.0)
    assert (R + h) ** 2 == pytest.approx(R ** 2 + dd ** 2 + 2 * R * dd * np.sin(e))


def test_link_budget_is_internally_consistent():
    link = impl.Link()
    b = link.budget()
    Pr = impl.friis_received_power(link.P_t, impl.from_db(link.G_t_dBi), impl.from_db(link.G_r_dBi),
                                   link.distance, link.wavelength) / impl.from_db(link.other_losses_dB)
    assert b["Pr_dBW"] == pytest.approx(10 * np.log10(Pr), abs=1e-9)
    assert b["Pr_dBm"] == pytest.approx(b["Pr_dBW"] + 30)
    EbN0 = Pr / (1.380649e-23 * link.T_sys * link.bit_rate)
    assert b["EbN0_dB"] == pytest.approx(10 * np.log10(EbN0), abs=1e-9)
    # doubling the bit rate costs 3 dB; going from 10 deg to zenith gains 20 log10(d10/h)
    b2 = impl.Link(bit_rate=2 * link.bit_rate).budget()
    assert b["margin_dB"] - b2["margin_dB"] == pytest.approx(3.0103, abs=1e-3)
    bz = impl.Link(elevation_deg=90.0).budget()
    assert bz["margin_dB"] - b["margin_dB"] == pytest.approx(20 * np.log10(link.distance / link.altitude), abs=1e-9)
    assert 0.0 < b["margin_dB"] < 30.0                                          # a closing but not absurd link


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "dipole_pattern.png").exists()
