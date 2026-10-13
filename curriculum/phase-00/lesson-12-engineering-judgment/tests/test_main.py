"""Public tests for P00.12 Engineering Judgment & Sanity Checks.

They assert physics and limiting cases: the ISA reproduces its defining constants and is continuous
at the tropopause, the lift coefficient scales as 1/V², the thrust-to-weight minimum sits at the
analytic optimum of the drag polar, airliner-class inputs pass, and absurd inputs are flagged with
a reason.

Run `pytest curriculum/phase-00/lesson-12-engineering-judgment` against your starter/main.py, or
`AEFS_IMPL=solution pytest ...` against the reference solution.
"""
import math

import numpy as np
import pytest

from conftest import load_impl

impl = load_impl(__file__)

A320 = dict(mass_kg=70_000.0, wing_area_m2=122.6, speed_m_s=230.0, altitude_m=11_000.0)
B777 = dict(mass_kg=250_000.0, wing_area_m2=427.8, speed_m_s=250.0, altitude_m=10_700.0)
B747 = dict(mass_kg=350_000.0, wing_area_m2=541.0, speed_m_s=250.0, altitude_m=10_700.0)


# --- ISA ---------------------------------------------------------------------------------------------

def test_isa_sea_level_constants():
    t, p, rho = impl.isa(0.0)
    assert t == pytest.approx(288.15)
    assert p == pytest.approx(101325.0)
    assert rho == pytest.approx(1.225, rel=1e-3)


def test_isa_tropopause_values():
    t, p, rho = impl.isa(11_000.0)
    assert t == pytest.approx(216.65, abs=1e-6)
    assert p == pytest.approx(22632.0, rel=2e-3)
    assert rho == pytest.approx(0.3639, rel=2e-3)


def test_isa_isothermal_above_tropopause_and_continuous():
    t15, p15, _ = impl.isa(15_000.0)
    assert t15 == pytest.approx(216.65, abs=1e-6)
    below = impl.isa(11_000.0 - 1e-6)
    above = impl.isa(11_000.0 + 1e-6)
    assert below[1] == pytest.approx(above[1], rel=1e-8)
    assert below[2] == pytest.approx(above[2], rel=1e-8)
    assert p15 == pytest.approx(12045.0, rel=3e-3)


def test_isa_pressure_and_density_decrease_monotonically():
    hs = np.linspace(0, 20_000, 81)
    ps = [impl.isa(h)[1] for h in hs]
    rhos = [impl.isa(h)[2] for h in hs]
    assert all(a > b for a, b in zip(ps, ps[1:]))
    assert all(a > b for a, b in zip(rhos, rhos[1:]))


def test_isa_obeys_ideal_gas_law():
    for h in (0.0, 5_000.0, 11_000.0, 18_000.0):
        t, p, rho = impl.isa(h)
        assert p == pytest.approx(rho * 287.05287 * t, rel=1e-9)


@pytest.mark.parametrize("h", [-1000.0, 20_001.0, 30_000.0])
def test_isa_rejects_altitude_outside_model(h):
    with pytest.raises(ValueError):
        impl.isa(h)


def test_speed_of_sound_sea_level():
    assert impl.speed_of_sound(288.15) == pytest.approx(340.3, abs=0.1)
    assert impl.speed_of_sound(4 * 288.15) == pytest.approx(2 * 340.3, abs=0.2)   # a ~ sqrt(T)


# --- cruise quantities ---------------------------------------------------------------------------------

def test_wing_loading_is_weight_per_area():
    assert impl.wing_loading(70_000.0, 122.6) == pytest.approx(70_000 * 9.80665 / 122.6)
    with pytest.raises(ValueError):
        impl.wing_loading(70_000.0, 0.0)


def test_lift_coefficient_scales_as_inverse_speed_squared():
    cl1 = impl.lift_coefficient(70_000.0, 122.6, 200.0, 10_000.0)
    cl2 = impl.lift_coefficient(70_000.0, 122.6, 400.0, 10_000.0)
    assert cl1 / cl2 == pytest.approx(4.0, rel=1e-9)
    # and linearly with mass
    assert impl.lift_coefficient(140_000.0, 122.6, 200.0, 10_000.0) == pytest.approx(2 * cl1)


def test_lift_coefficient_is_airliner_like_for_a320_class():
    cl = impl.lift_coefficient(**A320)
    assert 0.4 < cl < 0.7


def test_mach_number_a320_cruise():
    assert 0.75 < impl.mach_number(230.0, 11_000.0) < 0.82
    assert impl.mach_number(340.3, 0.0) == pytest.approx(1.0, abs=2e-3)


def test_thrust_to_weight_minimum_at_analytic_optimum():
    cd0, e, ar = 0.02, 0.8, 9.0
    speeds = np.linspace(120.0, 400.0, 2801)
    tw = np.array([impl.thrust_to_weight_required(70_000.0, 122.6, v, 11_000.0, cd0, e, ar) for v in speeds])
    i = int(np.argmin(tw))
    cl_opt = impl.lift_coefficient(70_000.0, 122.6, speeds[i], 11_000.0)
    assert cl_opt == pytest.approx(math.sqrt(math.pi * e * ar * cd0), rel=5e-3)
    assert tw[i] == pytest.approx(2 * math.sqrt(cd0 / (math.pi * e * ar)), rel=1e-4)
    # T/W is 1/(L/D) by construction
    assert tw[i] * impl.lift_to_drag(cl_opt, cd0, e, ar) == pytest.approx(1.0)


def test_lift_to_drag_limiting_cases():
    assert impl.lift_to_drag(0.0) == 0.0
    big = impl.lift_to_drag(1e6)
    assert big < 1e-3                                  # induced drag dominates
    assert impl.lift_to_drag(0.5, cd0=0.02, oswald=0.8, aspect_ratio=9.0) > impl.lift_to_drag(0.5, cd0=0.03, oswald=0.8, aspect_ratio=9.0)


@pytest.mark.parametrize("x, expected", [(950, 3), (0.002, -3), (1.0, 0), (20, 1), (3.0, 0), (4.0, 1), (-700, 3)])
def test_order_of_magnitude(x, expected):
    assert impl.order_of_magnitude(x) == expected


def test_order_of_magnitude_zero_is_an_error():
    with pytest.raises(ValueError):
        impl.order_of_magnitude(0.0)


# --- the sanity checker ------------------------------------------------------------------------------------

@pytest.mark.parametrize("case", [A320, B777, B747], ids=["A320", "B777", "B747"])
def test_airliner_class_inputs_pass(case):
    report = impl.check_airliner(**case)
    assert report.ok, report.summary()
    assert report.problems == []
    assert {f.quantity for f in report.flags} == {"mass", "altitude", "mach", "wing_loading", "lift_coefficient", "thrust_to_weight"}


def _flagged(report):
    return {f.quantity: f for f in report.problems}


def test_mass_in_pounds_is_flagged_with_reason():
    report = impl.check_airliner(mass_kg=154_000.0, wing_area_m2=122.6, speed_m_s=230.0, altitude_m=11_000.0)
    bad = _flagged(report)
    assert "lift_coefficient" in bad and "wing_loading" in bad
    assert all(f.reason for f in bad.values())
    assert not report.ok


def test_speed_in_km_h_is_flagged_as_mach():
    report = impl.check_airliner(mass_kg=70_000.0, wing_area_m2=122.6, speed_m_s=830.0, altitude_m=11_000.0)
    bad = _flagged(report)
    assert "mach" in bad and bad["mach"].value > 1.0
    assert "lift_coefficient" in bad and bad["lift_coefficient"].value < 0.3


def test_tiny_wing_is_flagged():
    report = impl.check_airliner(mass_kg=70_000.0, wing_area_m2=12.26, speed_m_s=230.0, altitude_m=11_000.0)
    bad = _flagged(report)
    assert {"wing_loading", "lift_coefficient", "thrust_to_weight"} <= set(bad)
    assert bad["lift_coefficient"].value > 0.8
    assert bad["thrust_to_weight"].value > 0.12 and "L/D" in bad["thrust_to_weight"].reason


def test_low_slow_flight_is_flagged():
    report = impl.check_airliner(mass_kg=70_000.0, wing_area_m2=122.6, speed_m_s=80.0, altitude_m=1_000.0)
    bad = _flagged(report)
    assert {"altitude", "mach", "lift_coefficient"} <= set(bad)
    # C_L ~ 1.6 puts T/W near 0.08: inside the band. A sanity check catches gross errors, not every error.
    assert "thrust_to_weight" not in bad


def test_impossible_inputs_raise_rather_than_flag():
    with pytest.raises(ValueError):
        impl.check_airliner(mass_kg=70_000.0, wing_area_m2=122.6, speed_m_s=230.0, altitude_m=30_000.0)
    with pytest.raises(ValueError):
        impl.check_airliner(mass_kg=-70_000.0, wing_area_m2=122.6, speed_m_s=230.0, altitude_m=11_000.0)


def test_every_flag_has_range_and_reason_when_failed():
    report = impl.check_airliner(**A320)
    for f in report.flags:
        assert f.low < f.high
        assert f.low <= f.value <= f.high
        assert f.reason == ""
    report = impl.check_airliner(mass_kg=1_000_000.0, wing_area_m2=122.6, speed_m_s=230.0, altitude_m=11_000.0)
    for f in report.problems:
        assert f.reason
        assert not (f.low <= f.value <= f.high)
    assert "FLAG" in report.summary() and "PASS" not in report.summary()


def test_main_runs(capsys):
    impl.main()
    out = capsys.readouterr().out
    assert "PASS" in out and "FLAG" in out
