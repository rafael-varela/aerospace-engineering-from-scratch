"""Public tests for P01.05 Sequences, Series & Taylor Expansions. They assert properties, not magic numbers.

Run `pytest curriculum/phase-01/lesson-05-series-and-taylor` against your starter/main.py, or with
`AEFS_IMPL=solution` against the reference solution. Hidden validation cases live in
assessments/phase-01/parameterized-cases/.
"""
from math import factorial

import numpy as np
import numpy.testing as npt
import pytest
import sympy as sp

from conftest import load_impl

impl = load_impl(__file__)

x = impl.model_symbols()["x"]


# ---------------------------------------------------------------- sequences and series

def test_geometric_partial_sums_approach_the_closed_form_and_diverge_outside_the_disc():
    for r in (0.5, -0.5, 0.9, 0.1):
        s = impl.partial_sums(impl.geometric_terms(r, 60))
        assert s.shape == (60,) and s[0] == 1.0
        npt.assert_allclose(s[:5], [(1 - r ** (k + 1)) / (1 - r) for k in range(5)], rtol=1e-14)   # finite geometric sum
        assert abs(s[-1] - impl.geometric_sum(r)) <= abs(r) ** 60 / (1 - abs(r)) * (1 + 1e-9) + 1e-15   # the tail bound, plus rounding
    assert np.all(np.diff(impl.partial_sums(impl.geometric_terms(0.5, 20))) > 0)                  # monotone for r > 0
    s = impl.partial_sums(impl.geometric_terms(1.5, 30))
    assert s[-1] > 1e4 and np.all(np.diff(s) > 0)                                                  # divergent: unbounded
    for r in (1.0, -1.0, 2.0):
        with pytest.raises(ValueError):
            impl.geometric_sum(r)
    with pytest.raises(ValueError):
        impl.geometric_terms(0.5, 0)


def test_harmonic_series_grows_like_log_n_without_bound():
    h = impl.harmonic_partial_sums(100000)
    assert h[0] == 1.0 and h[1] == 1.5 and np.all(np.diff(h) > 0)
    n = np.array([10, 100, 1000, 10000, 100000])
    gamma = 0.5772156649
    npt.assert_allclose(h[n - 1], np.log(n) + gamma, atol=0.06)                                    # H_n ~ ln n + gamma
    assert h[-1] - h[9999] == pytest.approx(np.log(10.0), abs=1e-3)                                # each decade adds ln 10


def test_ratio_test_estimate_and_radius_of_convergence():
    assert impl.ratio_test_estimate(impl.geometric_terms(0.7, 20)) == pytest.approx(0.7, rel=1e-12)
    fact = np.array([1.0 / factorial(k) for k in range(25)])                                       # exp: ratio -> 0
    assert impl.ratio_test_estimate(fact) < 0.05
    assert impl.ratio_test_estimate(1.0 / np.arange(1, 40.0)) == pytest.approx(1.0, abs=0.03)      # harmonic: test silent
    assert impl.radius_of_convergence(impl.taylor_coefficients(1 / (1 - 2 * x), 0, 15)) == pytest.approx(0.5, rel=1e-12)
    assert impl.radius_of_convergence(impl.taylor_coefficients(1 / (1 + x**2), 0, 30)) == pytest.approx(1.0, rel=1e-12)
    assert impl.radius_of_convergence(impl.taylor_coefficients(sp.exp(x), 0, 30)) > 25.0          # grows without bound
    r_pg = impl.radius_of_convergence(impl.binomial_coefficients(-0.5, 200))
    assert 1.0 < r_pg < 1.01                                                                       # Prandtl-Glauert: R -> 1 from above
    with pytest.raises(ValueError):
        impl.ratio_test_estimate(np.array([1.0, 0.0, 0.0]))
    with pytest.raises(ValueError):
        impl.radius_of_convergence(np.array([0.0, 3.0]))


def test_alternating_series_bound_is_honoured_by_sin_and_cos_partial_sums():
    for xv in (0.3, 1.0, 2.5):
        for kind, exact in (("sin", np.sin(xv)), ("cos", np.cos(xv))):
            t = impl.series_terms(xv, kind, 12)
            for k in range(3, 10):
                err = abs(np.sum(t[: k + 1]) - exact)
                assert err <= impl.alternating_series_bound(t, k) * (1 + 1e-12) + 4e-16                   # 4e-16: rounding of the sum
    with pytest.raises(ValueError):
        impl.alternating_series_bound(np.array([1.0, -0.5]), 1)


# ---------------------------------------------------------------- Taylor polynomials

def test_taylor_polynomial_matches_sympy_and_has_the_right_contact():
    p = impl.taylor_polynomial(sp.sin(x), 0, 7)
    assert sp.simplify(p - (x - x**3 / 6 + x**5 / 120 - x**7 / 5040)) == 0
    assert sp.simplify(impl.taylor_polynomial(sp.exp(x), 0, 4) - sp.series(sp.exp(x), x, 0, 5).removeO()) == 0
    q = impl.taylor_polynomial(sp.cos(x), np.pi / 3, 2)                                            # centre recognised as pi/3
    assert sp.simplify(q.subs(x, sp.pi / 3) - sp.Rational(1, 2)) == 0
    assert sp.simplify(sp.diff(q, x).subs(x, sp.pi / 3) + sp.sqrt(3) / 2) == 0
    assert sp.simplify(sp.diff(q, x, 2).subs(x, sp.pi / 3) + sp.Rational(1, 2)) == 0             # same value, slope, curvature
    f = sp.log(x)
    for n in (1, 2, 3):
        t = impl.taylor_polynomial(f, 1, n)
        for k in range(n + 1):
            assert sp.simplify(sp.diff(t, x, k).subs(x, 1) - sp.diff(f, x, k).subs(x, 1)) == 0
    with pytest.raises(ValueError):
        impl.taylor_polynomial(sp.sin(x), 0, -1)


def test_taylor_coefficients_match_the_known_series():
    npt.assert_allclose(impl.taylor_coefficients(sp.exp(x), 0, 8), [1.0 / factorial(k) for k in range(9)], rtol=1e-12)
    c = impl.taylor_coefficients(sp.sin(x), 0, 7)
    npt.assert_allclose(c[::2], 0.0, atol=1e-15)
    npt.assert_allclose(c[1::2], [1.0, -1.0 / 6, 1.0 / 120, -1.0 / 5040], rtol=1e-12)
    npt.assert_allclose(impl.taylor_coefficients(1 / sp.sqrt(1 - x), 0, 4), impl.binomial_coefficients(-0.5, 4) * (-1.0) ** np.arange(5), rtol=1e-12)
    npt.assert_allclose(impl.taylor_coefficients(sp.cos(x), np.pi / 3, 1), [0.5, -np.sqrt(3) / 2], rtol=1e-12)


def test_maclaurin_partial_sums_converge_to_numpy_everywhere_but_cancel_badly_for_large_x():
    xs = np.array([-3.0, -0.5, 0.0, 0.7, 2.0, 5.0])
    for kind, f in (("sin", np.sin), ("cos", np.cos), ("exp", np.exp)):
        t = impl.series_terms(xs, kind, 30)
        assert t.shape == (30, xs.size)
        npt.assert_allclose(impl.maclaurin_partial_sum(xs, kind, 30), f(xs), rtol=1e-10, atol=1e-12)
        npt.assert_allclose(np.sum(t, axis=0), impl.maclaurin_partial_sum(xs, kind, 30), rtol=1e-14)
    s1 = impl.maclaurin_partial_sum(0.2, "sin", 1)
    assert s1 == pytest.approx(0.2, rel=1e-15) and abs(s1 - np.sin(0.2)) < 0.2**3 / 6 * (1 + 1e-9)
    big = impl.maclaurin_partial_sum(40.0, "sin", 60)                                              # converges in exact arithmetic...
    assert abs(big - np.sin(40.0)) > 1e-3                                                          # ...but cancellation ruins it in float64
    assert abs(impl.maclaurin_partial_sum(-40.0, "exp", 120) - np.exp(-40.0)) > 1e-10 * np.exp(-40.0)
    with pytest.raises(ValueError):
        impl.series_terms(1.0, "tan", 3)
    with pytest.raises(ValueError):
        impl.series_terms(1.0, "sin", 0)


def test_lagrange_remainder_bound_holds_and_terms_needed_is_minimal():
    for xv in (0.3, 1.2, 3.0):
        for n in (1, 3, 5, 7):
            err = abs(impl.maclaurin_partial_sum(xv, "exp", n + 1) - np.exp(xv))
            assert err <= impl.lagrange_remainder_bound(np.exp(xv), xv, 0.0, n) * (1 + 1e-12)       # M = e^x on [0, x]
            err_sin = abs(np.sum(impl.series_terms(xv, "sin", 20)[: (n + 1) // 2 + 1]) - np.sin(xv))
            assert err_sin <= impl.lagrange_remainder_bound(1.0, xv, 0.0, 2 * ((n + 1) // 2) + 2) * (1 + 1e-12)   # |sin^(k)| <= 1
    assert impl.lagrange_remainder_bound(1.0, 0.5, 0.0, 3) == pytest.approx(0.5**4 / 24, rel=1e-14)
    for xv, tol in ((0.5, 1e-10), (2.0, 1e-10), (10.0, 1e-10), (1.0, 1e-3)):
        k = impl.terms_needed(xv, "sin", tol)
        t = impl.series_terms(xv, "sin", k + 2)
        assert abs(np.sum(t[:k]) - np.sin(xv)) < tol                                              # k terms suffice...
        assert k == 1 or abs(t[k - 1]) >= tol or abs(t[k - 2]) < abs(t[k - 1])                      # ...and k-1 would not be guaranteed
    assert impl.terms_needed(10.0, "sin", 1e-10) > impl.terms_needed(2.0, "sin", 1e-10) > impl.terms_needed(0.5, "sin", 1e-10)
    assert abs(impl.maclaurin_partial_sum(1.0, "exp", impl.terms_needed(1.0, "exp", 1e-8)) - np.e) < 1e-8
    with pytest.raises(ValueError):
        impl.lagrange_remainder_bound(1.0, 0.5, 0.0, -1)
    with pytest.raises(ValueError):
        impl.terms_needed(1.0, "sin", 0.0)


# ---------------------------------------------------------------- small-angle approximations

def test_small_angle_errors_follow_their_leading_taylor_term():
    th = np.radians(np.array([0.5, 1.0, 2.0, 4.0]))
    for a in impl.APPROXIMATIONS:
        err = impl.small_angle_error(th, a)
        lead = impl.small_angle_error_leading(th, a)
        assert err.shape == th.shape and np.all(err > 0) and np.all(np.diff(err) > 0)
        npt.assert_allclose(err, lead, rtol=2e-2)                                                  # next Taylor term rules at small angles
        order = np.polyfit(np.log(th), np.log(err), 1)[0]
        assert order == pytest.approx(4.0 if a == "cos2" else 2.0, abs=0.05)
    assert impl.small_angle_error(0.0, "sin") == 0.0
    assert impl.small_angle_error(np.radians(10.0), "cos2") < impl.small_angle_error(np.radians(10.0), "sin") < impl.small_angle_error(np.radians(10.0), "tan") < impl.small_angle_error(np.radians(10.0), "cos")
    with pytest.raises(ValueError):
        impl.small_angle_error(0.1, "sec")
    with pytest.raises(ValueError):
        impl.small_angle_error_leading(0.1, "sec")


def test_max_angle_for_error_inverts_the_error_and_orders_the_approximations():
    for a in impl.APPROXIMATIONS:
        for tol in (1e-2, 1e-3):
            th = impl.max_angle_for_error(a, tol)
            assert 0.0 < th < 0.5 * np.pi
            assert impl.small_angle_error(th, a) == pytest.approx(tol, rel=1e-6)
            assert impl.small_angle_error(th * 0.9, a) < tol < impl.small_angle_error(th * 1.1, a)
        assert impl.max_angle_for_error(a, 1e-3) < impl.max_angle_for_error(a, 1e-2)
    deg = {a: np.degrees(impl.max_angle_for_error(a, 1e-2)) for a in impl.APPROXIMATIONS}
    assert deg["cos"] < deg["tan"] < deg["sin"] < deg["cos2"]
    assert deg["sin"] == pytest.approx(np.degrees(np.sqrt(0.06)), rel=0.03)                        # theta^2/6 = 1 %
    assert impl.max_angle_for_error("sin", 0.9) == pytest.approx(0.5 * np.pi, abs=1e-6)           # sin ~ x is 57 % off at pi/2: never reaches 90 %
    assert impl.max_angle_for_error("cos2", 0.5) < 0.5 * np.pi                                     # cos -> 0 at pi/2: relative error explodes
    with pytest.raises(ValueError):
        impl.max_angle_for_error("sin", 1.5)


# ---------------------------------------------------------------- Prandtl-Glauert and the pendulum

def test_binomial_series_reproduces_prandtl_glauert_inside_the_unit_disc():
    npt.assert_allclose(impl.binomial_coefficients(3, 5), [1, 3, 3, 1, 0, 0], atol=1e-14)         # integer exponent terminates
    npt.assert_allclose(impl.binomial_coefficients(-0.5, 3), [1, -0.5, 0.375, -0.3125], rtol=1e-14)
    m = np.array([0.0, 0.3, 0.6, 0.8])
    npt.assert_allclose(impl.prandtl_glauert_series(m, 2), 1 + m**2 / 2, rtol=1e-14)
    npt.assert_allclose(impl.prandtl_glauert_series(m, 3), 1 + m**2 / 2 + 3 * m**4 / 8, rtol=1e-14)
    npt.assert_allclose(impl.prandtl_glauert_series(m, 80), impl.prandtl_glauert(m), rtol=1e-6)
    errs = [abs(impl.prandtl_glauert_series(0.6, k) - impl.prandtl_glauert(0.6)) for k in (2, 4, 8, 16)]
    assert all(e1 > e2 for e1, e2 in zip(errs, errs[1:]))                                          # monotone convergence
    assert abs(impl.prandtl_glauert_series(0.95, 10) - impl.prandtl_glauert(0.95)) > 0.3            # slow near M = 1
    assert np.all(np.diff(impl.prandtl_glauert(np.linspace(0, 0.99, 50))) > 0)
    with pytest.raises(ValueError):
        impl.prandtl_glauert(1.0)
    with pytest.raises(ValueError):
        impl.binomial_coefficients(0.5, -1)


def test_pendulum_period_series_converges_to_the_elliptic_integral_and_to_one_at_small_amplitude():
    for deg in (5.0, 20.0, 45.0, 90.0):
        th = np.radians(deg)
        exact = impl.pendulum_period_ratio_exact(th)
        assert exact > 1.0
        assert impl.pendulum_period_ratio(th, 1) == 1.0
        errs = np.array([abs(impl.pendulum_period_ratio(th, k) - exact) for k in (1, 2, 4, 8, 40)])
        assert np.all(np.diff(errs) <= 0) and errs[1] < errs[0] and errs[-1] < 1e-9 * (1 if deg < 60 else 1e3)
    th = np.radians(5.0)
    assert impl.pendulum_period_ratio(th, 2) - 1.0 == pytest.approx(np.sin(th / 2) ** 2 / 4, rel=1e-14)   # 1 + k^2/4
    assert impl.pendulum_period_ratio_exact(np.radians(90.0)) == pytest.approx(1.18034, rel=1e-5)
    assert abs(impl.pendulum_period_ratio_exact(np.radians(10.0)) - 1.0 - np.radians(10.0) ** 2 / 16) < 1e-5  # T/T0 ~ 1 + theta0^2/16
    with pytest.raises(ValueError):
        impl.pendulum_period_ratio(0.1, 0)


def test_euler_formula_holds_term_by_term():
    for th in (0.3, 1.0, 2.0, -1.5):
        assert impl.euler_formula_residual(th, 12) < 1e-12
    t = impl.series_terms(0.7, "exp", 8)
    npt.assert_allclose(t[::2], np.abs(impl.series_terms(0.7, "cos", 4)), rtol=1e-14)             # even powers are the cos terms
    npt.assert_allclose(t[1::2], np.abs(impl.series_terms(0.7, "sin", 4)), rtol=1e-14)
