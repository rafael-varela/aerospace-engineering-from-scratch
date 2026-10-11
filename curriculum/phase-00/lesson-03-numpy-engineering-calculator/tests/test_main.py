"""Public tests for P00.03. They assert physics, not magic numbers.

Run `pytest curriculum/phase-00/lesson-03-numpy-engineering-calculator` against your
starter/main.py, or with `AEFS_IMPL=solution` against the reference solution.
The only fixed numbers used are the ISA's defining sea-level values and the 11 km
tropopause values; everything else is a property: continuity, monotonicity, the ideal gas law,
the hydrostatic equation, scaling laws and shape preservation.
"""
import numpy as np
import pytest

from conftest import load_impl

impl = load_impl(__file__)

R_AIR = 287.05287
G0 = 9.80665


# --- Defining values -----------------------------------------------------------------------
def test_sea_level_values():
    out = impl.isa(0.0)
    assert out["T"] == pytest.approx(288.15, abs=1e-9)
    assert out["p"] == pytest.approx(101325.0, abs=1e-6)
    assert out["rho"] == pytest.approx(1.225, rel=1e-6)
    assert out["a"] == pytest.approx(340.294, abs=0.01)


def test_tropopause_values():
    out = impl.isa(11000.0)
    assert out["T"] == pytest.approx(216.65, abs=1e-9)
    assert out["p"] == pytest.approx(22632.0, rel=2e-4)
    assert out["rho"] == pytest.approx(0.3639, rel=5e-4)
    assert out["a"] == pytest.approx(295.07, abs=0.02)


def test_layer_bases_start_at_sea_level_and_descend_in_pressure():
    t_base, p_base = impl.layer_bases()
    assert t_base[0] == pytest.approx(288.15) and p_base[0] == pytest.approx(101325.0)
    assert np.all(np.diff(p_base) < 0), "pressure at each higher layer base must be lower"
    assert t_base[1] == pytest.approx(216.65, abs=1e-9)
    assert t_base[2] == pytest.approx(216.65, abs=1e-9), "11-20 km is isothermal"


# --- Continuity and monotonicity --------------------------------------------------------
def test_properties_are_continuous_across_layer_boundaries():
    eps = 1e-6
    for hb in impl.H_BASE[1:]:
        below, above = impl.isa(hb - eps), impl.isa(hb + eps)
        assert above["T"] == pytest.approx(below["T"], abs=1e-6)
        assert above["p"] == pytest.approx(below["p"], rel=1e-9)
        assert above["rho"] == pytest.approx(below["rho"], rel=1e-9)


def test_pressure_and_density_decrease_monotonically():
    h = np.linspace(0.0, impl.H_TOP, 5000)
    assert np.all(np.diff(impl.pressure(h)) < 0)
    assert np.all(np.diff(impl.density(h)) < 0)


def test_temperature_is_constant_in_isothermal_layers():
    h1 = np.linspace(11000.0, 20000.0, 50)
    h2 = np.linspace(47000.0, 51000.0, 50)
    assert np.ptp(impl.temperature(h1)) < 1e-9
    assert np.ptp(impl.temperature(h2)) < 1e-9


def test_temperature_gradient_matches_lapse_rate_in_each_layer():
    bases = np.append(impl.H_BASE, impl.H_TOP)
    for i, lapse in enumerate(impl.LAPSE):
        a, b = bases[i] + 10.0, bases[i + 1] - 10.0
        slope = (impl.temperature(b) - impl.temperature(a)) / (b - a)
        assert slope == pytest.approx(lapse, abs=1e-12)


# --- Physical laws -----------------------------------------------------------------------
def test_ideal_gas_law_holds_everywhere():
    h = np.linspace(0.0, impl.H_TOP, 1001)
    out = impl.isa(h)
    assert np.allclose(out["rho"], out["p"] / (R_AIR * out["T"]), rtol=1e-12)


def test_hydrostatic_equation_residual_is_small():
    """dp/dh = -rho g0 must hold to the accuracy of a central finite difference."""
    h = np.linspace(500.0, impl.H_TOP - 500.0, 2000)
    dh = 1.0
    dpdh = (impl.pressure(h + dh) - impl.pressure(h - dh)) / (2 * dh)
    rho_g = impl.density(h) * G0
    residual = np.abs(dpdh + rho_g) / rho_g
    assert np.max(residual) < 1e-5, f"max relative hydrostatic residual {np.max(residual):.2e}"


def test_speed_of_sound_scales_with_sqrt_temperature():
    t = np.array([200.0, 250.0, 300.0])
    assert np.allclose(impl.speed_of_sound(4 * t), 2 * impl.speed_of_sound(t), rtol=1e-12)
    with pytest.raises(ValueError):
        impl.speed_of_sound(-5.0)


def test_isothermal_layer_pressure_is_exponential():
    """In 11-20 km, ln p must be linear in h with slope -g0/(R T)."""
    h = np.linspace(11000.0, 20000.0, 200)
    logp = np.log(impl.pressure(h))
    slope = np.polyfit(h, logp, 1)[0]
    assert slope == pytest.approx(-G0 / (R_AIR * 216.65), rel=1e-9)


# --- Geometric to geopotential -----------------------------------------------------------
def test_geopotential_conversion_properties():
    assert impl.geometric_to_geopotential(0.0) == pytest.approx(0.0, abs=1e-12)
    z = np.linspace(0.0, 86000.0, 100)
    h = impl.geometric_to_geopotential(z)
    assert np.all(h[1:] <= z[1:]), "geopotential altitude never exceeds geometric altitude"
    assert np.all(np.diff(h) > 0)
    assert impl.geometric_to_geopotential(86000.0) == pytest.approx(84852.0, rel=1e-4)
    assert impl.H_TOP == pytest.approx(impl.geometric_to_geopotential(86000.0), rel=1e-9)


# --- Vectorization and interface ----------------------------------------------------------
def test_scalar_in_float_out_and_shape_preserved():
    assert isinstance(impl.temperature(5000.0), float)
    assert isinstance(impl.pressure(5000.0), float)
    h2d = np.linspace(0.0, 80000.0, 12).reshape(3, 4)
    for f in (impl.temperature, impl.pressure, impl.density):
        assert f(h2d).shape == (3, 4)


def test_vectorized_agrees_with_scalar_loop():
    h = np.linspace(0.0, impl.H_TOP, 97)
    vec = impl.pressure(h)
    loop = np.array([impl.pressure(float(x)) for x in h])
    assert np.allclose(vec, loop, rtol=1e-12, atol=0.0)


def test_out_of_range_raises():
    with pytest.raises(ValueError):
        impl.temperature(-1.0)
    with pytest.raises(ValueError):
        impl.pressure(np.array([0.0, 90000.0]))
    with pytest.raises(ValueError):
        impl.density(impl.H_TOP + 10.0)


def test_isa_table_shape_columns_and_consistency():
    z = np.array([0.0, 11000.0, 50000.0, 86000.0])
    table = impl.isa_table(z)
    assert table.shape == (4, 6)
    assert np.allclose(table[:, 0], z)
    assert np.allclose(table[:, 1], impl.geometric_to_geopotential(z))
    assert np.allclose(table[:, 2], impl.temperature(table[:, 1]))
    assert np.allclose(table[:, 4], table[:, 3] / (R_AIR * table[:, 2]), rtol=1e-12)
    assert np.allclose(table[:, 5], impl.speed_of_sound(table[:, 2]))
    with pytest.raises(ValueError):
        impl.isa_table(z.reshape(2, 2))
