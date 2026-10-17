"""Public tests for P01.04 Techniques of Integration. They assert exactness, bounds, limits and convergence, not magic numbers.

Run `pytest curriculum/phase-01/lesson-04-techniques-of-integration` against your starter/main.py, or with
`AEFS_IMPL=solution` against the reference solution. Hidden validation cases live in
assessments/phase-01/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
import sympy as sp

from conftest import load_impl

impl = load_impl(__file__)

S_, CR, TAPER, SIGMA = 5.5, 1.6, 0.6, 12.0
ELL = S_ / np.log(1.0 / TAPER)
SYM = impl.model_symbols()
y, s, ell, sigma, cr = SYM["y"], SYM["s"], SYM["ell"], SYM["sigma"], SYM["cr"]


# ---------------------------------------------------------------- planforms and closed forms

def test_chord_functions_at_root_and_tip():
    assert impl.chord_linear(0.0, CR, TAPER, S_) == pytest.approx(CR)
    assert impl.chord_linear(S_, CR, TAPER, S_) == pytest.approx(TAPER * CR)
    assert impl.chord_elliptic(0.0, CR, S_) == pytest.approx(CR)
    assert impl.chord_elliptic(S_, CR, S_) == pytest.approx(0.0, abs=1e-12)
    assert impl.chord_exponential(0.0, CR, ELL) == pytest.approx(CR)
    assert impl.chord_exponential(S_, CR, ELL) == pytest.approx(TAPER * CR)          # ell chosen for the same tip chord
    yy = np.linspace(0.0, S_, 7)
    for fn in (lambda q: impl.chord_linear(q, CR, TAPER, S_), lambda q: impl.chord_elliptic(q, CR, S_), lambda q: impl.chord_exponential(q, CR, ELL)):
        assert np.all(np.diff(fn(yy)) < 0.0)                                            # chord shrinks outboard


def test_linear_closed_forms_agree_with_numerical_integrals_and_scale():
    lin = lambda q: impl.chord_linear(q, CR, TAPER, S_)
    mass = impl.wing_mass_linear(SIGMA, CR, TAPER, S_)
    assert mass == pytest.approx(2.0 * SIGMA * impl.simpson(lin, 0.0, S_, 2))           # Simpson exact for a line
    ybar = impl.spanwise_centroid_linear(CR, TAPER, S_)
    assert ybar == pytest.approx(impl.simpson(lambda q: q * lin(q), 0.0, S_, 2) / impl.simpson(lin, 0.0, S_, 2))
    assert S_ / 3.0 < ybar < S_ / 2.0                                                  # between a triangle's and a rectangle's centroid
    assert impl.spanwise_centroid_linear(CR, 1.0, S_) == pytest.approx(S_ / 2.0)        # rectangular wing
    assert impl.spanwise_centroid_linear(CR, 0.0, S_) == pytest.approx(S_ / 3.0)        # triangular wing
    i_lin = impl.roll_inertia_linear(SIGMA, CR, TAPER, S_)
    assert i_lin == pytest.approx(impl.roll_inertia_numeric(lin, SIGMA, S_, 2, "simpson"), rel=1e-13)   # y^2 c(y) is a cubic
    assert impl.roll_inertia_linear(SIGMA, CR, TAPER, 2.0 * S_) == pytest.approx(8.0 * i_lin)            # ~ s^3
    assert impl.roll_inertia_linear(2.0 * SIGMA, CR, TAPER, S_) == pytest.approx(2.0 * i_lin)
    assert impl.roll_inertia_linear(SIGMA, CR, 1.0, S_) == pytest.approx(2.0 * SIGMA * CR * S_**3 / 3.0)  # rectangular: 2 (1/3) m' s^3


def test_elliptic_closed_form_matches_numeric_and_symbolic():
    ell_fn = lambda q: impl.chord_elliptic(q, CR, S_)
    exact = impl.roll_inertia_elliptic(SIGMA, CR, S_)
    assert exact == pytest.approx(impl.roll_inertia_numeric(ell_fn, SIGMA, S_, 4000, "simpson"), rel=1e-5)   # order 1.5 at the tip: 2e-6 at n = 4000
    assert exact == pytest.approx(impl.roll_inertia_elliptic_substituted(SIGMA, CR, S_, 4), rel=1e-14)    # periodic integrand: exact already
    sym = impl.symbolic_roll_inertia_elliptic()
    assert sp.simplify(sym - sp.pi * sigma * cr * s**3 / 8) == 0
    assert exact < impl.roll_inertia_linear(SIGMA, CR, 1.0, S_)                           # less than the rectangular wing


def test_exponential_closed_form_by_parts_matches_symbolic_numeric_and_limit():
    exact = impl.roll_inertia_exponential(SIGMA, CR, ELL, S_)
    sym = impl.symbolic_roll_inertia_exponential()
    assert float(sym.subs({sigma: SIGMA, cr: CR, ell: ELL, s: S_})) == pytest.approx(exact, rel=1e-12)
    hand = 2 * sigma * cr * (2 * ell**3 - sp.exp(-s / ell) * (s**2 * ell + 2 * s * ell**2 + 2 * ell**3))
    assert sp.simplify(sym - hand) == 0
    exp_fn = lambda q: impl.chord_exponential(q, CR, ELL)
    assert exact == pytest.approx(impl.roll_inertia_numeric(exp_fn, SIGMA, S_, 200, "simpson"), rel=1e-8)
    inf = impl.roll_inertia_exponential_infinite(SIGMA, CR, ELL)
    assert inf > exact                                                                    # the tail beyond s is positive
    assert impl.roll_inertia_exponential(SIGMA, CR, ELL, 40.0 * ELL) == pytest.approx(inf, rel=1e-12)     # s >> ell
    assert impl.roll_inertia_exponential(SIGMA, CR, ELL, 1e-6) == pytest.approx(0.0, abs=1e-12)            # s -> 0


# ---------------------------------------------------------------- rules, exactness and error bounds

def test_rules_exactness_on_polynomials_and_input_checks():
    cubic = lambda x: 2.0 * x**3 - x**2 + 3.0 * x - 1.0
    exact = 2.0 * 16.0 / 4.0 - 8.0 / 3.0 + 3.0 * 4.0 / 2.0 - 2.0                           # on [0, 2]
    assert impl.simpson(cubic, 0.0, 2.0, 2) == pytest.approx(exact, rel=1e-14)             # exact for cubics with 2 panels
    assert impl.trapezoid(lambda x: 3.0 * x + 1.0, 0.0, 2.0, 1) == pytest.approx(8.0)       # exact for lines
    assert impl.midpoint(lambda x: 3.0 * x + 1.0, 0.0, 2.0, 1) == pytest.approx(8.0)
    assert abs(impl.trapezoid(cubic, 0.0, 2.0, 2) - exact) > 1e-3                           # but not for cubics
    for bad in (0, -2):
        with pytest.raises(ValueError):
            impl.trapezoid(cubic, 0.0, 2.0, bad)
    for bad in (0, 1, 3):
        with pytest.raises(ValueError):
            impl.simpson(cubic, 0.0, 2.0, bad)
    with pytest.raises(ValueError):
        impl.roll_inertia_numeric(cubic, 1.0, 1.0, 4, "gauss")


def test_trapezoid_is_second_order_and_simpson_fourth_on_a_smooth_integrand():
    f, exact = np.exp, np.e - 1.0
    ns = np.array([4, 8, 16, 32])
    et = np.array([abs(impl.trapezoid(f, 0.0, 1.0, int(n)) - exact) for n in ns])
    es = np.array([abs(impl.simpson(f, 0.0, 1.0, int(n)) - exact) for n in ns])
    npt.assert_allclose(np.log2(et[:-1] / et[1:]), 2.0, atol=0.05)
    npt.assert_allclose(np.log2(es[:-1] / es[1:]), 4.0, atol=0.1)
    # Midpoint error is about half the trapezoid error and of opposite sign for a convex function.
    em = impl.midpoint(f, 0.0, 1.0, 16) - exact
    assert em < 0.0 < impl.trapezoid(f, 0.0, 1.0, 16) - exact
    assert abs(em) == pytest.approx(0.5 * et[2], rel=0.02)


def test_error_bounds_are_honoured_and_not_absurdly_loose():
    f, exact = np.sin, 1.0 - np.cos(2.0)
    for n in (4, 16, 64):
        et = abs(impl.trapezoid(f, 0.0, 2.0, n) - exact)
        bt = impl.trapezoid_error_bound(1.0, 0.0, 2.0, n)                                  # |f''| <= 1
        assert et <= bt < 20.0 * et
        es = abs(impl.simpson(f, 0.0, 2.0, n) - exact)
        bs = impl.simpson_error_bound(1.0, 0.0, 2.0, n)
        assert es <= bs < 20.0 * es
    assert impl.simpson_error_bound(1.0, 0.0, 2.0, 8) == pytest.approx(impl.simpson_error_bound(1.0, 0.0, 2.0, 4) / 16.0)


def test_tip_singularity_slows_simpson_but_the_substitution_removes_it():
    ell_fn = lambda q: impl.chord_elliptic(q, CR, S_)
    exact = impl.roll_inertia_elliptic(SIGMA, CR, S_)
    ns = np.array([16, 32, 64, 128, 256])
    es = np.array([abs(impl.roll_inertia_numeric(ell_fn, SIGMA, S_, int(n), "simpson") - exact) for n in ns])
    rates = np.log2(es[:-1] / es[1:])
    assert np.all(rates > 1.3) and np.all(rates < 1.7)                                     # order 1.5, not 4: sqrt at the tip
    sub = np.array([abs(impl.roll_inertia_elliptic_substituted(SIGMA, CR, S_, int(n)) - exact) for n in (2, 4, 8)])
    assert np.all(sub[1:] < 1e-10 * exact)                                                 # trapezoid on a periodic integrand


# ---------------------------------------------------------------- improper integrals

def test_improper_integral_converges_to_the_closed_form_limit():
    exp_fn = lambda q: 2.0 * SIGMA * q**2 * impl.chord_exponential(q, CR, ELL)
    cutoffs = ELL * np.array([1.0, 2.0, 4.0, 8.0, 16.0])
    seq = impl.improper_integral(exp_fn, 0.0, cutoffs)
    inf = impl.roll_inertia_exponential_infinite(SIGMA, CR, ELL)
    assert np.all(np.diff(seq) > 0.0) and np.all(seq < inf)                                # monotone from below
    assert seq[-1] == pytest.approx(inf, rel=1e-4)                                          # tail beyond 16 ell is 1.6e-5 of the total
    assert inf - seq[-1] == pytest.approx(inf - impl.roll_inertia_exponential(SIGMA, CR, ELL, cutoffs[-1]), rel=1e-3)   # and equals the closed-form tail
    npt.assert_allclose(seq, [impl.roll_inertia_exponential(SIGMA, CR, ELL, r) for r in cutoffs], rtol=1e-8)
    # A divergent tail: int_1^R dx/x grows without bound (as log R).
    div = impl.improper_integral(lambda x: 1.0 / x, 1.0, np.array([10.0, 100.0, 1000.0]))
    npt.assert_allclose(div[:2], np.log([10.0, 100.0]), rtol=1e-6)
    npt.assert_allclose(np.diff(div), np.log(10.0), rtol=1e-3)     # each decade adds ln 10: no limit (Simpson with h = 0.5 near x = 1 is only 2e-4 accurate)
    with pytest.raises(ValueError):
        impl.improper_integral(exp_fn, 0.0, np.array([-1.0]))


# ---------------------------------------------------------------- partial fractions: take-off time

def test_takeoff_time_closed_form_matches_simpson_and_its_limits():
    m, thrust, k = 1000.0, 2500.0, 0.5 * 1.225 * 0.04 * 15.0
    vt = np.sqrt(thrust / k)
    for v in (10.0, 30.0, 60.0):
        assert impl.takeoff_time(v, thrust, k, m) == pytest.approx(impl.takeoff_time_numeric(v, thrust, k, m, 400), rel=1e-8)
    assert impl.takeoff_time(0.0, thrust, k, m) == 0.0
    assert impl.takeoff_time(30.0, thrust, k, m) > m * 30.0 / thrust                       # drag lengthens the roll
    assert impl.takeoff_time(30.0, thrust, 1e-9, m) == pytest.approx(m * 30.0 / thrust, rel=1e-6)   # no drag: t = m v / T
    assert impl.takeoff_time(0.999 * vt, thrust, k, m) > 5.0 * impl.takeoff_time(0.5 * vt, thrust, k, m)    # log divergence at v_t (ratio 6.9)
    assert np.all(np.diff(impl.takeoff_time(np.array([5.0, 10.0, 20.0, 40.0]), thrust, k, m)) > 0.0)
    for bad in (vt, 1.1 * vt, -1.0):
        with pytest.raises(ValueError):
            impl.takeoff_time(bad, thrust, k, m)
    with pytest.raises(ValueError):
        impl.takeoff_time(10.0, 0.0, k, m)
