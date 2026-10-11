"""Public tests for P00.09 Symbolic Math with SymPy. They assert identities and physics, not magic numbers.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-00/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
import sympy as sp

from conftest import load_impl

impl = load_impl(__file__)

S = impl.model_symbols()
t, v, m, g, k = S["t"], S["v"], S["m"], S["g"], S["k"]

# A 1.2 kg instrument package, C_D = 1.0, A = 0.01 m**2 at ISA sea level.
M, G, RHO, CD, AREA = 1.2, 9.80665, 1.225, 1.0, 0.01
K = 0.5 * RHO * CD * AREA
VT = np.sqrt(M * G / K)


# ---------------------------------------------------------------- symbolic identities

def test_drag_constant_value_and_validation():
    assert impl.drag_constant(RHO, CD, AREA) == pytest.approx(K)
    assert impl.drag_constant(2 * RHO, CD, AREA) == pytest.approx(2 * K)   # linear in rho
    for bad in [(0.0, CD, AREA), (RHO, -1.0, AREA), (RHO, CD, float("nan"))]:
        with pytest.raises(ValueError):
            impl.drag_constant(*bad)


def test_acceleration_expression_structure():
    a = impl.acceleration_expr()
    assert sp.simplify(a - (g - k * v**2 / m)) == 0
    assert a.subs(v, 0) == g                                   # at rest only gravity acts
    assert sp.simplify(a.subs(k, 0) - g) == 0                  # no drag: free fall
    assert sp.diff(a, v).subs(v, 0) == 0                       # drag is quadratic, no linear term


def test_terminal_velocity_zeroes_acceleration_and_scales_correctly():
    vt = impl.terminal_velocity_expr()
    assert sp.simplify(impl.acceleration_expr().subs(v, vt)) == 0
    assert sp.simplify(vt - sp.sqrt(g * m / k)) == 0
    # Doubling the mass multiplies v_t by sqrt(2); doubling k divides by sqrt(2).
    assert sp.simplify(vt.subs(m, 2 * m) / vt - sp.sqrt(2)) == 0
    assert sp.simplify(vt.subs(k, 2 * k) / vt - 1 / sp.sqrt(2)) == 0


def test_velocity_solves_the_ode_with_zero_initial_speed():
    vel = impl.velocity_expr()
    assert vel.subs(t, 0) == 0
    assert impl.ode_residual(vel) == 0
    # Late-time limit is the terminal velocity.
    assert sp.simplify(sp.limit(vel, t, sp.oo) - impl.terminal_velocity_expr()) == 0


def test_ode_residual_is_nonzero_for_a_wrong_candidate():
    """A test that only ever sees zero proves nothing: free fall must fail the drag ODE."""
    assert impl.ode_residual(g * t) != 0


def test_position_is_the_integral_of_velocity():
    x = impl.position_expr()
    assert x.subs(t, 0) == 0
    assert sp.simplify(sp.diff(x, t) - impl.velocity_expr()) == 0


def test_no_drag_limit_is_free_fall():
    assert sp.simplify(impl.no_drag_limit() - g * t) == 0


def test_small_time_series_leading_terms():
    ser = impl.small_time_series(4)
    poly = sp.Poly(sp.expand(ser), t)
    assert sp.simplify(poly.coeff_monomial(t) - g) == 0
    assert poly.coeff_monomial(t**2) == 0
    assert sp.simplify(poly.coeff_monomial(t**3) + g**2 * k / (3 * m)) == 0  # drag slows the fall
    with pytest.raises(ValueError):
        impl.small_time_series(1)


def test_acceleration_derivative_expression():
    assert sp.simplify(impl.acceleration_derivative_expr() + 2 * k * v / m) == 0


# ---------------------------------------------------------------- lambdified functions vs known values

def test_lambdified_acceleration_matches_hand_values():
    a = impl.lambdify_acceleration()
    assert a(0.0, M, G, K) == pytest.approx(G)
    assert a(VT, M, G, K) == pytest.approx(0.0, abs=1e-12)
    vs = np.array([0.0, 10.0, 20.0, 30.0])
    npt.assert_allclose(a(vs, M, G, K), G - K * vs**2 / M, rtol=1e-14)
    assert np.all(np.diff(a(vs, M, G, K)) < 0.0)                # acceleration falls as speed rises


def test_lambdified_velocity_matches_hand_values_and_limits():
    vel = impl.lambdify_velocity()
    ts = np.linspace(0.0, 30.0, 61)
    out = vel(ts, M, G, K)
    assert out.shape == ts.shape
    npt.assert_allclose(out, VT * np.tanh(G * ts / VT), rtol=1e-13)
    assert out[0] == 0.0
    assert np.all(np.diff(out) > 0.0) and np.all(out < VT)      # monotone, bounded by v_t
    assert out[-1] == pytest.approx(VT, rel=1e-5)                # 30 s is ~7 time constants
    # Early on, drag has not yet mattered: v ~ g t within 1 %.
    assert vel(0.2, M, G, K) == pytest.approx(G * 0.2, rel=1e-2)
    # Vanishing drag reproduces free fall to high accuracy.
    assert vel(3.0, M, G, 1e-9) == pytest.approx(G * 3.0, rel=1e-6)


def test_lambdified_position_is_consistent_with_velocity():
    x = impl.lambdify_position()
    vel = impl.lambdify_velocity()
    assert x(0.0, M, G, K) == 0.0
    # Trapezoid integral of v(t) against the closed form x(t): error ~ dt**2.
    ts = np.linspace(0.0, 10.0, 4001)
    x_num = np.concatenate([[0.0], np.cumsum(0.5 * (vel(ts[1:], M, G, K) + vel(ts[:-1], M, G, K)) * np.diff(ts))])
    npt.assert_allclose(x_num, x(ts, M, G, K), atol=1e-4, rtol=0.0)
    assert np.all(np.diff(x(ts, M, G, K)) > 0.0) and np.all(x(ts[1:], M, G, K) < VT * ts[1:])  # falls, slower than v_t
    # No-drag check: x -> g t**2 / 2 for tiny k.
    assert x(2.0, M, G, 1e-9) == pytest.approx(0.5 * G * 4.0, rel=1e-6)


# ---------------------------------------------------------------- symbolic derivative vs finite difference

def test_central_difference_on_a_known_function():
    fd = impl.central_difference(np.sin, np.array([0.0, 1.0]), 1e-5)
    npt.assert_allclose(fd, np.cos([0.0, 1.0]), atol=1e-9)
    with pytest.raises(ValueError):
        impl.central_difference(np.sin, np.array([0.0]), 0.0)


def test_symbolic_derivative_agrees_with_finite_difference():
    # a(v) is quadratic, so the central difference is exact: only rounding remains.
    for h in (1e-1, 1e-2, 1e-3):
        assert impl.derivative_check(20.0, M, G, K, h) < 1e-10


def test_velocity_derivative_check_converges_at_second_order():
    errs = np.array([impl.velocity_derivative_check(3.0, M, G, K, h) for h in (4e-2, 2e-2, 1e-2)])
    assert np.all(errs < 1e-4)
    rates = np.log2(errs[:-1] / errs[1:])
    assert np.all(rates > 1.8), f"expected second-order convergence, got {rates}"
    # And the tightest check is essentially zero against rounding.
    assert impl.velocity_derivative_check(3.0, M, G, K, 1e-4) < 1e-8
