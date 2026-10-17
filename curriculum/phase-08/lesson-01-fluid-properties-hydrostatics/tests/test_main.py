"""Public tests for P08.01 Fluid Properties & Hydrostatics. They assert physics, not magic numbers:
the closed-form ISA reproduces the published anchor values (labelled published data) and is
continuous at the layer boundaries; it satisfies dp/dh = -rho g0 when substituted back; an
independent RK4 integration of the hydrostatic equation agrees with it and converges at fourth
order; limiting cases (an isothermal layer, small heights, zero acceleration, free fall) hold;
Sutherland's law has the published sea-level value and the right trend; and Archimedes' principle is
recovered by integrating surface pressure, for a liquid and for the compressible atmosphere.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

# ISA anchor values, published data (geopotential altitude in m: pressure in Pa)
ANCHORS_P = {11000.0: 22632.06, 20000.0: 5474.889, 32000.0: 868.019, 47000.0: 110.906}


# ---------------------------------------------------------------- the standard atmosphere

def test_temperature_profile_matches_the_layer_definition():
    h = np.array([0.0, 11000.0, 15000.0, 20000.0, 32000.0, 47000.0])
    npt.assert_allclose(impl.isa_temperature(h), [288.15, 216.65, 216.65, 216.65, 228.65, 270.65], atol=1e-9)
    # lapse rates by finite difference inside each layer
    for h_mid, lapse in ((5000.0, -6.5e-3), (15000.0, 0.0), (25000.0, 1.0e-3), (40000.0, 2.8e-3)):
        d = (impl.isa_temperature(h_mid + 1.0) - impl.isa_temperature(h_mid - 1.0)) / 2.0
        assert d == pytest.approx(lapse, abs=1e-12)


def test_pressure_reproduces_published_isa_anchors():
    assert float(impl.isa_pressure(0.0)) == pytest.approx(101325.0, rel=1e-12)
    for h, p in ANCHORS_P.items():
        assert float(impl.isa_pressure(h)) == pytest.approx(p, rel=1e-5)        # published data


def test_density_reproduces_published_values():
    assert float(impl.isa_density(0.0)) == pytest.approx(1.225, rel=1e-4)
    assert float(impl.isa_density(11000.0)) == pytest.approx(0.36392, rel=1e-4)


def test_pressure_and_temperature_are_continuous_at_layer_boundaries():
    for hb in (11000.0, 20000.0, 32000.0):
        for f in (impl.isa_pressure, impl.isa_temperature):
            lo, hi = float(f(hb - 1e-6)), float(f(hb + 1e-6))
            assert lo == pytest.approx(hi, rel=1e-9)


def test_closed_form_satisfies_the_hydrostatic_equation():
    h = np.array([1000.0, 6000.0, 14000.0, 19000.0, 25000.0, 31000.0, 36000.0, 45000.0])
    dh = 0.5
    dpdh = (impl.isa_pressure(h + dh) - impl.isa_pressure(h - dh)) / (2 * dh)
    npt.assert_allclose(dpdh, -impl.isa_density(h) * impl.G0, rtol=1e-7)


def test_ideal_gas_law_is_obeyed():
    h = np.linspace(0.0, 47000.0, 48)
    npt.assert_allclose(impl.isa_pressure(h), impl.isa_density(h) * impl.R_AIR * impl.isa_temperature(h), rtol=1e-12)


def test_isothermal_layer_is_exponential_with_the_scale_height():
    H = impl.scale_height(216.65)
    p11 = float(impl.isa_pressure(11000.0))
    for dz in (1000.0, 5000.0, 9000.0):
        assert float(impl.isa_pressure(11000.0 + dz)) == pytest.approx(p11 * np.exp(-dz / H), rel=1e-12)


def test_vectorised_and_rejects_out_of_range():
    h = np.linspace(0.0, 47000.0, 11).reshape(11)
    assert np.shape(impl.isa_pressure(h)) == (11,)
    assert np.all(np.diff(impl.isa_pressure(h)) < 0)                    # pressure falls monotonically
    with pytest.raises(ValueError):
        impl.isa_pressure(50000.0)


# ---------------------------------------------------------------- the independent route: RK4

def test_rk4_integration_matches_closed_forms_and_anchors():
    h, p = impl.integrate_hydrostatic(47000.0, 100.0)
    assert h[0] == 0.0 and p[0] == 101325.0 and h[-1] == pytest.approx(47000.0)
    npt.assert_allclose(p, impl.isa_pressure(h), rtol=1e-8)
    for hh, pp in ANCHORS_P.items():
        assert p[int(round(hh / 100.0))] == pytest.approx(pp, rel=1e-5)          # published data


def test_rk4_converges_at_fourth_order():
    exact = float(impl.isa_pressure(47000.0))
    errs = [abs(impl.integrate_hydrostatic(47000.0, dh)[1][-1] - exact) for dh in (1000.0, 500.0, 250.0)]
    orders = impl.observed_order(errs)
    assert np.all(orders > 3.7) and np.all(orders < 4.3)


def test_rk4_with_constant_temperature_gives_an_exponential():
    h, p = impl.integrate_hydrostatic(8000.0, 100.0, T_of_h=lambda z: 250.0)
    npt.assert_allclose(p, 101325.0 * np.exp(-h / impl.scale_height(250.0)), rtol=1e-9)


def test_incompressible_estimate_is_first_order_accurate():
    # p0 - rho0 g h is the first two Taylor terms of the ISA: the error grows like h^2
    e1 = abs(float(impl.incompressible_pressure(250.0)) - float(impl.isa_pressure(250.0)))
    e2 = abs(float(impl.incompressible_pressure(500.0)) - float(impl.isa_pressure(500.0)))
    assert e2 / e1 == pytest.approx(4.0, rel=0.02)
    assert float(impl.incompressible_pressure(0.0)) == 101325.0


def test_geopotential_altitude_is_below_geometric_and_equal_near_the_ground():
    z = np.array([10.0, 11000.0, 47000.0])
    h = impl.geopotential_altitude(z)
    assert np.all(h < z)
    assert h[0] == pytest.approx(10.0, rel=1e-5)
    assert h[2] == pytest.approx(47000.0 * 6356766.0 / (6356766.0 + 47000.0))


# ---------------------------------------------------------------- viscosity

def test_sutherland_viscosity_sea_level_trend_and_scaling():
    assert float(impl.sutherland_viscosity(288.15)) == pytest.approx(1.7894e-5, rel=1e-4)   # published
    T = np.array([150.0, 216.65, 288.15, 400.0])
    mu = impl.sutherland_viscosity(T)
    assert np.all(np.diff(mu) > 0)                                       # a gas: mu rises with T
    # far above S the law tends to C1 sqrt(T): mu(4T)/mu(T) -> 2
    assert float(impl.sutherland_viscosity(4e6) / impl.sutherland_viscosity(1e6)) == pytest.approx(2.0, rel=1e-4)
    with pytest.raises(ValueError):
        impl.sutherland_viscosity(-1.0)


def test_kinematic_viscosity_grows_with_altitude():
    nu = impl.kinematic_viscosity(np.array([0.0, 3000.0, 11000.0, 20000.0]))
    assert np.all(np.diff(nu) > 0)
    assert float(nu[0]) == pytest.approx(1.7894e-5 / 1.225, rel=2e-4)


# ---------------------------------------------------------------- manometers, tanks, buoyancy

def test_manometer_limits_and_inclination():
    assert float(impl.manometer_delta_p(0.1, rho_m=1000.0, rho_f=1000.0)) == 0.0     # same fluid: no reading
    dp = float(impl.manometer_delta_p(0.02, rho_m=998.2, rho_f=0.0))
    assert dp == pytest.approx(998.2 * 9.80665 * 0.02)
    # an inclined tube spreads the same pressure over 1/sin(angle) more length
    assert float(impl.manometer_delta_p(0.02 / np.sin(np.radians(15.0)), angle_deg=15.0)) == pytest.approx(
        float(impl.manometer_delta_p(0.02)), rel=1e-12)


def test_tank_pressure_under_acceleration():
    p_pad = impl.tank_bottom_pressure(3e5, 1141.0, 2.0)
    assert p_pad == pytest.approx(3e5 + 1141.0 * 9.80665 * 2.0)
    assert impl.tank_bottom_pressure(3e5, 1141.0, 2.0, accel=-9.80665) == pytest.approx(3e5)    # free fall
    # the hydrostatic head is linear in the effective gravity
    d1 = impl.tank_bottom_pressure(0.0, 1141.0, 2.0, accel=9.80665)
    assert d1 == pytest.approx(2 * (p_pad - 3e5))


def test_archimedes_by_surface_pressure_in_a_liquid():
    rho, R = 998.2, 0.3
    weight = rho * 9.80665 * 4.0 / 3.0 * np.pi * R ** 3
    for depth in (1.0, 5.0):                                             # buoyancy does not depend on depth
        F = impl.pressure_force_on_sphere(lambda z: 101325.0 + rho * 9.80665 * (-z), -depth, R, n=400)
        assert F == pytest.approx(weight, rel=1e-4)
    # a uniform pressure exerts no net force on a closed surface
    assert abs(impl.pressure_force_on_sphere(lambda z: np.full_like(z, 5e5), 0.0, R)) < 1e-6 * 5e5 * R ** 2


def test_archimedes_in_the_compressible_atmosphere():
    R, zc = 2.0, 3000.0
    F = impl.pressure_force_on_sphere(lambda z: impl.isa_pressure(z), zc, R, n=400)
    W = impl.displaced_weight(lambda z: impl.isa_density(z), zc, R)
    assert F == pytest.approx(W, rel=1e-4)
    # midpoint rule: error falls as n^-2
    errs = [abs(impl.pressure_force_on_sphere(lambda z: impl.isa_pressure(z), zc, R, n=n) - W) for n in (25, 50, 100)]
    assert np.all(impl.observed_order(errs) > 1.8)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "isa_atmosphere.png").stat().st_size > 10_000
