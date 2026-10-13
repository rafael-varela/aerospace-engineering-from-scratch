"""Public tests for P00.05 Units & Dimensional Consistency.

They assert unit round-trips, dimension checks, and the lbf·s vs N·s factor behind the
Mars Climate Orbiter loss. Run `pytest <lesson dir>` against your starter/main.py, or
`AEFS_IMPL=solution pytest <lesson dir>` against the reference solution. Hidden
validation cases live in assessments/phase-00/parameterized-cases/.
"""
import numpy as np
import pint
import pytest

from conftest import load_impl

impl = load_impl(__file__)
Q_ = impl.ureg.Quantity

LBF_IN_N = 0.45359237 * 9.80665  # exact SI definition of the pound-force


# ------------------------------------------------------------------ round trips

@pytest.mark.parametrize("value,a,b", [
    (35_000.0, "ft", "m"),
    (250.0, "knot", "m/s"),
    (14.696, "psi", "Pa"),
    (1.0, "lbf", "N"),
    (3.0, "slug / ft**3", "kg / m**3"),
])
def test_round_trip_returns_to_start(value, a, b):
    assert impl.round_trip(value, a, b) == pytest.approx(value, rel=1e-12)


def test_to_si_gives_base_units():
    q = impl.to_si(Q_(1.0, "knot"))
    assert q.check("[length] / [time]")
    assert str(q.units) == "meter / second"
    assert q.magnitude == pytest.approx(1852.0 / 3600.0, rel=1e-12)


def test_parse_quantity_and_unknown_unit():
    q = impl.parse_quantity("35000 ft")
    assert q.check("[length]")
    assert impl.magnitude_in(q, "m") == pytest.approx(10_668.0, rel=1e-12)
    with pytest.raises(pint.UndefinedUnitError):
        impl.parse_quantity("3 furlongs_per_fortnite")


def test_magnitude_in_refuses_wrong_dimension():
    with pytest.raises(pint.DimensionalityError):
        impl.magnitude_in(Q_(101_325.0, "Pa"), "m")


# ------------------------------------------------------------------ dimension checks

def test_mismatched_addition_raises_dimensionality_error():
    with pytest.raises(pint.DimensionalityError):
        impl.add_checked(Q_(1.0, "m"), Q_(1.0, "s"))


def test_compatible_addition_converts():
    total = impl.add_checked(Q_(1.0, "m"), Q_(1.0, "ft"))
    assert impl.magnitude_in(total, "m") == pytest.approx(1.3048, rel=1e-12)


def test_add_checked_refuses_bare_numbers():
    with pytest.raises(TypeError):
        impl.add_checked(1.0, Q_(1.0, "m"))


def test_has_dimension():
    assert impl.has_dimension(Q_(3.0, "N * s"), "[force] * [time]")
    assert impl.has_dimension(Q_(3.0, "N * s"), "[mass] * [length] / [time]")
    assert not impl.has_dimension(Q_(3.0, "N"), "[energy]")


# ------------------------------------------------------------------ aerospace formulas

def test_dynamic_pressure_value_and_units():
    q = impl.dynamic_pressure(Q_(1.225, "kg / m**3"), Q_(100.0, "m/s"))
    assert q.check("[pressure]")
    assert impl.magnitude_in(q, "Pa") == pytest.approx(0.5 * 1.225 * 100.0**2, rel=1e-12)
    # same physics in imperial inputs gives the same pressure
    q2 = impl.dynamic_pressure(Q_(1.225, "kg / m**3").to("slug / ft**3"), Q_(100.0, "m/s").to("knot"))
    assert impl.magnitude_in(q2, "Pa") == pytest.approx(impl.magnitude_in(q, "Pa"), rel=1e-10)


def test_dynamic_pressure_rejects_swapped_arguments():
    with pytest.raises(pint.DimensionalityError):
        impl.dynamic_pressure(Q_(100.0, "m/s"), Q_(1.225, "kg / m**3"))


def test_total_impulse_and_delta_v_dimensions():
    impulse = impl.total_impulse(Q_(10.0, "N"), Q_(2.0, "s"))
    assert impl.magnitude_in(impulse, "N * s") == pytest.approx(20.0)
    dv = impl.delta_v(impulse, Q_(4.0, "kg"))
    assert dv.check("[length] / [time]")
    assert impl.magnitude_in(dv, "m/s") == pytest.approx(5.0)
    with pytest.raises(pint.DimensionalityError):
        impl.delta_v(Q_(10.0, "N"), Q_(4.0, "kg"))  # a force is not an impulse


def test_total_impulse_scales_linearly():
    i1 = impl.total_impulse(Q_(10.0, "N"), Q_(2.0, "s"))
    i2 = impl.total_impulse(Q_(20.0, "N"), Q_(2.0, "s"))
    assert impl.magnitude_in(i2, "N*s") == pytest.approx(2 * impl.magnitude_in(i1, "N*s"))


# ------------------------------------------------------------------ Mars Climate Orbiter

def test_lbf_s_to_N_s_factor_matches_definition():
    f = impl.lbf_s_to_N_s_factor()
    assert f == pytest.approx(LBF_IN_N, rel=1e-12)
    assert 4.44 < f < 4.45


def test_mco_misread_impulse_is_underestimated_by_the_factor():
    true, misread = impl.mco_misread_impulse(0.5)
    assert true.check("[force] * [time]") and misread.check("[force] * [time]")
    assert impl.magnitude_in(true, "N*s") == pytest.approx(0.5 * LBF_IN_N, rel=1e-12)
    assert impl.magnitude_in(misread, "N*s") == pytest.approx(0.5)
    assert impl.mco_underestimate_ratio(0.5) == pytest.approx(LBF_IN_N, rel=1e-12)
    # the ratio is independent of the size of the firing
    assert impl.mco_underestimate_ratio(7.3) == pytest.approx(impl.mco_underestimate_ratio(0.01))


def test_accumulated_delta_v_error_is_additive_and_has_speed_units():
    mass = Q_(600.0, "kg")
    one = impl.accumulated_delta_v_error(np.array([1.0]), mass)
    many = impl.accumulated_delta_v_error(np.full(10, 1.0), mass)
    assert one.check("[length] / [time]")
    assert impl.magnitude_in(many, "m/s") == pytest.approx(10 * impl.magnitude_in(one, "m/s"), rel=1e-12)
    expected = (LBF_IN_N - 1.0) / 600.0
    assert impl.magnitude_in(one, "m/s") == pytest.approx(expected, rel=1e-12)
    # doubling the mass halves the missed velocity change
    half = impl.accumulated_delta_v_error(np.array([1.0]), 2 * mass)
    assert impl.magnitude_in(half, "m/s") == pytest.approx(0.5 * impl.magnitude_in(one, "m/s"))


def test_main_runs(capsys):
    impl.main()
    out = capsys.readouterr().out
    assert "DimensionalityError" in out
    assert "lbf" in out
