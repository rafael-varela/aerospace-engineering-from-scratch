"""Public tests for P01.08 Multivariable Functions & Partial Derivatives. They assert calculus identities, not magic numbers.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-01/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
import sympy as sp

from conftest import load_impl

impl = load_impl(__file__)

ALPHA, V = np.radians(4.0), 55.0


def _wing(**kw):
    return impl.Wing(**kw)


# ---------------------------------------------------------------- the function of two variables

def test_lift_is_quadratic_in_speed_and_broadcasts():
    w = _wing()
    L1 = impl.lift(ALPHA, V, w)
    assert impl.lift(ALPHA, 2 * V, w) == pytest.approx(4 * L1, rel=1e-12)
    assert impl.lift(ALPHA, 0.0, w) == 0.0
    assert impl.lift(w.alpha0, V, w) == pytest.approx(0.0, abs=1e-9)             # zero-lift angle
    A = np.radians(np.linspace(0, 10, 5))[:, None]
    Vs = np.linspace(40, 80, 3)[None, :]
    out = impl.lift(A, Vs, w)
    assert out.shape == (5, 3)
    npt.assert_allclose(out[:, 0], impl.lift(A[:, 0], 40.0, w), rtol=1e-12)
    assert isinstance(L1, float)
    with pytest.raises(ValueError):
        impl.lift(ALPHA, -1.0, w)


def test_linear_model_is_the_k3_zero_limit():
    w = _wing(k3=0.0)
    npt.assert_allclose(impl.lift_coefficient(np.radians([0.0, 4.0, 8.0]), w),
                        w.a * (np.radians([0.0, 4.0, 8.0]) - w.alpha0), rtol=1e-12)
    la, lv = impl.lift_partials(ALPHA, V, w)
    assert la == pytest.approx(0.5 * w.rho * V**2 * w.S * w.a, rel=1e-12)
    assert lv == pytest.approx(2 * impl.lift(ALPHA, V, w) / V, rel=1e-12)       # d(V**2)/dV = 2 V
    sa, sv = impl.relative_sensitivities(ALPHA, V, w)
    assert sa == pytest.approx(1.0, rel=1e-12) and sv == pytest.approx(2.0, rel=1e-12)
    with pytest.raises(ValueError):
        impl.stall_angle(w)


# ---------------------------------------------------------------- partial derivatives three ways

def test_analytic_partials_match_central_differences_at_second_order():
    w = _wing()
    f = lambda x: impl.lift(x[0], x[1], w)
    x0 = np.array([ALPHA, V])
    la, lv = impl.lift_partials(ALPHA, V, w)
    errs = np.array([abs(impl.partial_derivative(f, x0, 0, h) - la) for h in (4e-2, 2e-2, 1e-2)])
    rates = np.log2(errs[:-1] / errs[1:])
    assert np.all(rates > 1.9), rates                                            # L is cubic in alpha: h**2 truncation
    assert impl.partial_derivative(f, x0, 0, 1e-4) == pytest.approx(la, rel=1e-7)
    # L is quadratic in V, so the central difference in V is exact up to rounding.
    assert impl.partial_derivative(f, x0, 1, 1e-2) == pytest.approx(lv, rel=1e-10)
    with pytest.raises(ValueError):
        impl.partial_derivative(f, x0, 0, 0.0)
    with pytest.raises(ValueError):
        impl.partial_derivative(f, x0, 2, 1e-3)


def test_partial_derivative_holds_the_other_variable_fixed():
    """df/dx of f(x, y) = x**2 y must give 2 x y, not the total derivative along some path."""
    f = lambda x: x[0] ** 2 * x[1]
    assert impl.partial_derivative(f, np.array([3.0, 2.0]), 0, 1e-3) == pytest.approx(12.0, rel=1e-8)
    assert impl.partial_derivative(f, np.array([3.0, 2.0]), 1, 1e-3) == pytest.approx(9.0, rel=1e-10)


def test_sympy_partials_agree_with_analytic_ones():
    w = _wing()
    dla, dlv, s = impl.lift_partials_sympy()
    subs = {s["alpha"]: ALPHA, s["V"]: V, s["rho"]: w.rho, s["S"]: w.S, s["a"]: w.a, s["alpha0"]: w.alpha0, s["k3"]: w.k3}
    la, lv = impl.lift_partials(ALPHA, V, w)
    assert float(dla.subs(subs)) == pytest.approx(la, rel=1e-12)
    assert float(dlv.subs(subs)) == pytest.approx(lv, rel=1e-12)
    # Structure: dL/dV is 2 L / V for any C_L, and the mixed partials are equal (Clairaut).
    L = sp.Rational(1, 2) * s["rho"] * s["V"] ** 2 * s["S"] * (s["a"] * (s["alpha"] - s["alpha0"]) - s["k3"] * (s["alpha"] - s["alpha0"]) ** 3)
    assert sp.simplify(dlv - 2 * L / s["V"]) == 0
    assert sp.simplify(sp.diff(dla, s["V"]) - sp.diff(dlv, s["alpha"])) == 0


def test_mixed_partials_are_equal_clairaut():
    w = _wing()
    f = lambda x: impl.lift(x[0], x[1], w)
    x0 = np.array([ALPHA, V])
    fav = impl.mixed_partial(f, x0, 0, 1, 1e-3)
    fva = impl.mixed_partial(f, x0, 1, 0, 1e-3)
    assert fav == pytest.approx(fva, rel=1e-6)
    # Analytic: d/dV (1/2 rho V**2 S C_L') = rho V S C_L'(alpha).
    dcl = w.a - 3 * w.k3 * (ALPHA - w.alpha0) ** 2
    assert fav == pytest.approx(w.rho * V * w.S * dcl, rel=1e-5)
    g = lambda x: np.sin(x[0]) * np.exp(x[1])
    assert impl.mixed_partial(g, np.array([0.7, 0.2]), 0, 1, 1e-3) == pytest.approx(np.cos(0.7) * np.exp(0.2), rel=1e-5)


# ---------------------------------------------------------------- linearization

def test_tangent_plane_error_is_second_order():
    w = _wing()
    errs = impl.tangent_plane_error_study(ALPHA, V, [4.0, 2.0, 1.0, 0.5], w)
    ratios = errs[:-1] / errs[1:]
    assert np.all(ratios > 3.8) and np.all(ratios < 4.2), ratios
    assert impl.lift_linearized(ALPHA, V, ALPHA, V, w) == pytest.approx(impl.lift(ALPHA, V, w), rel=1e-14)
    # The tangent plane is exact along V for the linear-in-alpha part? No: L is quadratic in V, so the
    # error at fixed alpha is exactly L_VV dV**2 / 2 = 1/2 rho S C_L dV**2.
    dv = 5.0
    err = impl.lift(ALPHA, V + dv, w) - impl.lift_linearized(ALPHA, V + dv, ALPHA, V, w)
    assert err == pytest.approx(0.5 * w.rho * w.S * impl.lift_coefficient(ALPHA, w) * dv**2, rel=1e-10)


def test_relative_sensitivity_drops_toward_stall():
    w = _wing()
    sa_low, _ = impl.relative_sensitivities(np.radians(2.0), V, w)
    sa_high, _ = impl.relative_sensitivities(np.radians(14.0), V, w)
    assert 0.9 < sa_low <= 1.0 and sa_high < sa_low
    sa_stall, sv = impl.relative_sensitivities(impl.stall_angle(w), V, w)
    assert sa_stall == pytest.approx(0.0, abs=1e-12)                              # dC_L/dalpha = 0 at C_L max
    assert sv == pytest.approx(2.0, rel=1e-12)


def test_stall_angle_is_the_maximum_of_cl():
    w = _wing()
    a_s = impl.stall_angle(w)
    assert w.alpha0 < a_s < np.radians(25.0)
    la, _ = impl.lift_partials(a_s, V, w)
    assert la == pytest.approx(0.0, abs=1e-9)
    al = np.radians(np.linspace(-2.0, 25.0, 2701))
    cl = impl.lift_coefficient(al, w)
    assert abs(al[np.argmax(cl)] - a_s) < np.radians(0.02)


# ---------------------------------------------------------------- level curves and the chain rule

def test_level_curve_carries_constant_lift_and_implicit_slope_matches():
    w = _wing()
    L0 = impl.lift(ALPHA, V, w)
    al = np.radians(np.linspace(0.0, 12.0, 25))
    Vs = impl.level_curve_speed(al, L0, w)
    npt.assert_allclose(impl.lift(al, Vs, w), L0, rtol=1e-12)
    assert np.all(np.diff(Vs) < 0.0)                                              # more alpha, less speed needed
    assert impl.level_curve_speed(ALPHA, L0, w) == pytest.approx(V, rel=1e-12)
    h = 1e-5
    fd = (impl.level_curve_speed(ALPHA + h, L0, w) - impl.level_curve_speed(ALPHA - h, L0, w)) / (2 * h)
    assert impl.level_curve_slope(ALPHA, V, w) == pytest.approx(fd, rel=1e-7)
    with pytest.raises(ValueError):
        impl.level_curve_speed(w.alpha0, L0, w)                                   # C_L = 0: no speed carries L0


def test_chain_rule_rate_matches_finite_difference_along_a_manoeuvre():
    w = _wing()
    alpha_t = lambda t: ALPHA + np.radians(2.0) * t + np.radians(0.3) * t**2
    V_t = lambda t: V - 1.0 * t + 0.05 * t**2
    t0, h = 1.5, 1e-4
    dalpha = (alpha_t(t0 + h) - alpha_t(t0 - h)) / (2 * h)
    dV = (V_t(t0 + h) - V_t(t0 - h)) / (2 * h)
    dL_fd = (impl.lift(alpha_t(t0 + h), V_t(t0 + h), w) - impl.lift(alpha_t(t0 - h), V_t(t0 - h), w)) / (2 * h)
    assert impl.lift_rate(alpha_t(t0), V_t(t0), dalpha, dV, w) == pytest.approx(dL_fd, rel=1e-7)
    # Along the level curve the lift does not change: L_alpha + L_V dV/dalpha = 0.
    slope = impl.level_curve_slope(ALPHA, V, w)
    assert impl.lift_rate(ALPHA, V, 1.0, slope, w) == pytest.approx(0.0, abs=1e-6 * impl.lift(ALPHA, V, w))


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "lift_sensitivity.png").stat().st_size > 1000
