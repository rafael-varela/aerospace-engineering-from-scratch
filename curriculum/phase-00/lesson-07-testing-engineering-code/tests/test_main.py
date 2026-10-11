"""Public tests for P00.07 Testing Engineering Code. They assert physics, not magic numbers.

This file is also the lesson's worked example of a test suite: tolerance-based assertions
with ``pytest.approx`` and ``numpy.testing``, limiting cases, raising on invalid inputs, a
parametrized test, and a regression test against a stored reference table.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-00/parameterized-cases/.
"""
import math

import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

RHO_SL = 1.225  # kg/m**3, ISA sea-level density


# ---------------------------------------------------------------- 1. why exact equality fails

def test_exact_float_equality_fails_for_decimal_fractions():
    """0.1, 0.2 and 0.3 are not representable in binary; their sum is off by one ulp."""
    assert 0.1 + 0.2 != 0.3                       # the trap
    assert 0.1 + 0.2 == pytest.approx(0.3)        # the fix
    assert impl.is_close(0.1 + 0.2, 0.3)          # your own version of the fix


def test_is_close_matches_math_isclose():
    cases = [(1.0, 1.0 + 1e-12), (1.0, 1.0 + 1e-6), (1e-300, 2e-300), (0.0, 1e-12), (-5.0, -5.0)]
    for a, b in cases:
        for rel, ab in [(1e-9, 0.0), (1e-5, 0.0), (0.0, 1e-9)]:
            assert impl.is_close(a, b, rel, ab) == math.isclose(a, b, rel_tol=rel, abs_tol=ab), (a, b, rel, ab)


def test_is_close_edge_cases():
    assert impl.is_close(math.inf, math.inf)
    assert not impl.is_close(math.inf, -math.inf)
    assert not impl.is_close(math.nan, math.nan)
    assert not impl.is_close(0.0, 1e-12)                  # relative tolerance alone cannot pass near zero
    assert impl.is_close(0.0, 1e-12, abs_tol=1e-9)        # an absolute tolerance is needed there
    with pytest.raises(ValueError):
        impl.is_close(1.0, 1.0, rel_tol=-1e-9)


def test_rocket_mass_stepwise_drifts_but_stays_within_tolerance():
    """Repeated subtraction accumulates rounding error; the closed form does not."""
    m0, mdot, dt, n = 50.0, 0.1, 1.0, 100
    stepwise = impl.rocket_mass_stepwise(m0, mdot, dt, n)
    exact = impl.rocket_mass_closed_form(m0, mdot, dt * np.arange(n + 1))
    assert stepwise.shape == (n + 1,)
    assert stepwise[0] == m0
    # Exact equality at the last step is NOT guaranteed; a tolerance is.
    assert stepwise[-1] == pytest.approx(exact[-1], rel=1e-12)
    npt.assert_allclose(stepwise, exact, rtol=1e-12, atol=0.0)
    # The drift is real but tiny: bounded by n ulps of m0.
    drift = np.max(np.abs(stepwise - exact))
    assert drift <= n * np.finfo(float).eps * m0
    # Mass is monotonically non-increasing (a physical property).
    assert np.all(np.diff(stepwise) <= 0.0)


def test_rocket_mass_raises_when_propellant_exhausted():
    with pytest.raises(ValueError):
        impl.rocket_mass_stepwise(1.0, 0.5, 1.0, 3)
    with pytest.raises(ValueError):
        impl.rocket_mass_closed_form(1.0, 0.5, 3.0)


# ---------------------------------------------------------------- 2. the calculator: analytic and limiting cases

def test_dynamic_pressure_analytic_value():
    # ISA sea level, 50 m/s: q = 0.5 * 1.225 * 2500 = 1531.25 Pa exactly in decimal, not in binary.
    assert impl.dynamic_pressure(RHO_SL, 50.0) == pytest.approx(1531.25, rel=1e-12)


def test_dynamic_pressure_limits_and_scaling():
    assert impl.dynamic_pressure(RHO_SL, 0.0) == 0.0                     # v -> 0 gives q = 0 exactly
    q1 = impl.dynamic_pressure(RHO_SL, 30.0)
    q2 = impl.dynamic_pressure(RHO_SL, 60.0)
    assert q2 / q1 == pytest.approx(4.0)                                 # q scales with v**2
    assert impl.dynamic_pressure(2 * RHO_SL, 30.0) == pytest.approx(2 * q1)  # and linearly with rho


def test_coefficients_round_trip():
    """C_L from a lift force, then the lift force back from C_L: must recover the input."""
    rho, v, s = 0.9, 80.0, 16.2
    lift = 12_000.0
    cl = impl.lift_coefficient(lift, rho, v, s)
    q = impl.dynamic_pressure(rho, v)
    assert cl * q * s == pytest.approx(lift, rel=1e-12)
    cd = impl.drag_coefficient(600.0, rho, v, s)
    assert cl / cd == pytest.approx(lift / 600.0, rel=1e-12)             # L/D is the force ratio


def test_coefficients_vectorized():
    v = np.array([40.0, 50.0, 60.0])
    cl = impl.lift_coefficient(5000.0, RHO_SL, v, 10.0)
    assert cl.shape == v.shape
    assert np.all(np.diff(cl) < 0.0)                                      # same lift at higher speed needs less C_L
    npt.assert_allclose(cl * impl.dynamic_pressure(RHO_SL, v) * 10.0, 5000.0, rtol=1e-12)


def test_drag_polar_limiting_cases():
    assert impl.drag_polar(0.0, cd0=0.02, k=0.05) == pytest.approx(0.02)   # no lift: only C_D0
    assert impl.drag_polar(0.7, cd0=0.02, k=0.0) == pytest.approx(0.02)    # k = 0: C_D independent of C_L
    assert impl.drag_polar(0.7, cd0=0.0, k=0.05) == pytest.approx(0.05 * 0.49)
    cls = np.linspace(-1.0, 1.0, 21)
    cd = impl.drag_polar(cls, cd0=0.02, k=0.05)
    npt.assert_allclose(cd, cd[::-1], rtol=1e-14)                          # even in C_L
    assert np.all(cd >= 0.02)


# ---------------------------------------------------------------- 3. parametrized test

@pytest.mark.parametrize(
    "cd0, k",
    [(0.020, 0.050), (0.015, 0.040), (0.030, 0.080), (0.008, 0.020)],
    ids=["light-aircraft", "glider-ish", "draggy", "sailplane"],
)
def test_max_lift_to_drag_is_the_maximum(cd0, k):
    """The analytic optimum must beat every sampled point on the polar and sit at its peak."""
    cl_star, ld_max = impl.max_lift_to_drag(cd0, k)
    cls = np.linspace(0.05, 2.0, 400)
    ld = impl.lift_to_drag(cls, cd0, k)
    assert ld_max >= np.max(ld)                                   # analytic max bounds the sampled curve
    assert ld_max == pytest.approx(np.max(ld), rel=1e-4)          # and the sampled curve nearly reaches it
    assert impl.lift_to_drag(cl_star, cd0, k) == pytest.approx(ld_max, rel=1e-12)
    # At the optimum, C_D0 equals the induced drag k C_L**2 (a classic result).
    assert k * cl_star ** 2 == pytest.approx(cd0, rel=1e-12)


@pytest.mark.parametrize(
    "bad_kwargs",
    [dict(rho=0.0, v=50.0, area=10.0), dict(rho=-1.2, v=50.0, area=10.0),
     dict(rho=1.2, v=0.0, area=10.0), dict(rho=1.2, v=50.0, area=0.0),
     dict(rho=float("nan"), v=50.0, area=10.0), dict(rho=1.2, v=float("inf"), area=10.0)],
)
def test_lift_coefficient_rejects_invalid_inputs(bad_kwargs):
    with pytest.raises(ValueError):
        impl.lift_coefficient(1000.0, **bad_kwargs)


def test_other_invalid_inputs_raise():
    with pytest.raises(ValueError):
        impl.dynamic_pressure(1.225, -10.0)             # a negative speed is a sign error upstream
    with pytest.raises(ValueError):
        impl.drag_coefficient(-5.0, 1.225, 50.0, 10.0)  # negative drag is unphysical
    with pytest.raises(ValueError):
        impl.drag_polar(0.5, cd0=-0.01)
    with pytest.raises(ValueError):
        impl.drag_polar(float("nan"))
    with pytest.raises(ValueError):
        impl.max_lift_to_drag(cd0=0.02, k=0.0)          # k = 0 would make (L/D)max infinite


# ---------------------------------------------------------------- 4. regression test against a stored table

# Generated once from a reviewed implementation with cd0 = 0.020, k = 0.050 and stored to
# 8 significant figures. If this test fails after a code change, either the change broke the
# calculator or the table is deliberately being updated: never "fix" the table to make it pass.
REFERENCE_POLAR = np.array([
    # C_L,   C_D,        L/D
    (0.0, 0.02, 0.0),
    (0.2, 0.022, 9.0909091),
    (0.4, 0.028, 14.285714),
    (0.6, 0.038, 15.789474),
    (0.8, 0.052, 15.384615),
    (1.0, 0.07, 14.285714),
    (1.2, 0.092, 13.043478),
    (1.4, 0.118, 11.864407),
])


def test_polar_table_matches_stored_reference():
    table = impl.polar_table(REFERENCE_POLAR[:, 0], cd0=0.020, k=0.050)
    assert table.shape == REFERENCE_POLAR.shape
    # The stored values carry 8 significant figures, so the tolerance must be looser than
    # that rounding (1e-7 relative) but far tighter than any physical effect.
    npt.assert_allclose(table, REFERENCE_POLAR, rtol=2e-7, atol=1e-12)


def test_polar_table_columns_are_consistent():
    cls = np.linspace(0.1, 1.5, 15)
    table = impl.polar_table(cls, cd0=0.025, k=0.045)
    npt.assert_allclose(table[:, 0], cls)
    npt.assert_allclose(table[:, 1], impl.drag_polar(cls, 0.025, 0.045))
    npt.assert_allclose(table[:, 2], table[:, 0] / table[:, 1], rtol=1e-14)
