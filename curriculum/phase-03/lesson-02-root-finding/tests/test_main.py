"""Public tests for P03.02 Root Finding. They assert residuals, convergence orders and orbital invariants, with scipy.optimize as the oracle.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import math

import numpy as np
import numpy.testing as npt
import pytest
from scipy import optimize

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(302)
GOLDEN = (1.0 + math.sqrt(5.0)) / 2.0


def kepler_ref(M, e):
    return optimize.brentq(lambda E: E - e * math.sin(E) - M, M - e - 1e-9, M + e + 1e-9, xtol=1e-15, rtol=4 * 2.0 ** -52)


# ---------------------------------------------------------------- the four iterations on textbook functions

def test_bisection_halves_the_bracket_and_meets_its_bound():
    tol = 1e-10
    res = impl.bisection(math.cos, 0.0, 2.0, tol=tol)
    assert res.root == pytest.approx(math.pi / 2, abs=tol)
    assert res.iterates.size == impl.bisection_steps_needed(0.0, 2.0, tol)
    # Every midpoint error is within the shrinking bound (b - a) / 2^k.
    k = np.arange(1, res.iterates.size + 1)
    assert np.all(np.abs(res.iterates - math.pi / 2) <= 2.0 / 2.0 ** k + 1e-15)
    assert res.n_evals == res.iterates.size + 2
    with pytest.raises(ValueError):
        impl.bisection(math.cos, 0.0, 1.0)                     # no sign change
    # An endpoint that is an exact root is returned at once.
    assert impl.bisection(lambda x: x - 1.0, 1.0, 3.0).root == 1.0
    with pytest.raises(ValueError):
        impl.bisection_steps_needed(1.0, 0.0, 1e-3)


def test_newton_converges_quadratically_to_sqrt2_and_matches_scipy():
    f, df = (lambda x: x * x - 2.0), (lambda x: 2.0 * x)
    res = impl.newton(f, df, 1.0)
    assert res.root == pytest.approx(math.sqrt(2.0), abs=1e-15)
    assert res.iterates.size - 1 <= 6                        # quadratic convergence from 1.0
    assert res.root == pytest.approx(optimize.newton(f, 1.0, fprime=df, tol=1e-14), abs=1e-14)
    errs = impl.errors_from_iterates(res.iterates, math.sqrt(2.0))
    q = impl.observed_convergence_order(errs)
    assert q.size >= 2 and abs(q[-1] - 2.0) < 0.1
    assert res.n_evals == 2 * (res.iterates.size - 1)


def test_newton_failure_modes_raise():
    with pytest.raises(ZeroDivisionError):
        impl.newton(lambda x: x * x - 1.0, lambda x: 2.0 * x, 0.0)       # f'(0) = 0
    with pytest.raises((RuntimeError, ZeroDivisionError)):
        impl.newton(math.atan, lambda x: 1.0 / (1.0 + x * x), 1.5)     # diverges: |x0| > 1.39
    with pytest.raises(RuntimeError):
        impl.newton(lambda x: x ** 3 - 2 * x + 2, lambda x: 3 * x * x - 2, 0.0, max_iter=20)   # the 0 -> 1 -> 0 cycle


def test_newton_on_a_double_root_is_only_linear_with_rate_one_half():
    res = impl.newton(lambda x: (x - 1.0) ** 2, lambda x: 2.0 * (x - 1.0), 2.0, tol=1e-10)
    errs = impl.errors_from_iterates(res.iterates, 1.0)
    npt.assert_allclose(errs[1:8] / errs[:7], 0.5, rtol=1e-9)
    assert res.iterates.size - 1 > 25                        # ~34 steps instead of ~6


def test_secant_matches_brentq_and_raises_on_a_flat_slope():
    f = lambda x: x ** 3 - 2.0                              # noqa: E731
    res = impl.secant(f, 1.0, 2.0)
    assert res.root == pytest.approx(2.0 ** (1.0 / 3.0), abs=1e-14)
    assert res.root == pytest.approx(optimize.brentq(f, 1.0, 2.0, xtol=1e-15), abs=1e-14)
    assert res.n_evals == res.iterates.size - 1             # two to start, one per step; the last iterate is not evaluated
    with pytest.raises(ZeroDivisionError):
        impl.secant(lambda x: x * x, -1.0, 1.0)              # f(x0) = f(x1), both nonzero


def test_fixed_point_converges_at_the_rate_g_prime_and_diverges_when_g_prime_exceeds_one():
    res = impl.fixed_point(math.cos, 1.0, tol=1e-13)
    root = optimize.brentq(lambda x: math.cos(x) - x, 0.0, 1.0, xtol=1e-15)
    assert res.root == pytest.approx(root, abs=1e-12)
    rate = impl.observed_linear_rate(impl.errors_from_iterates(res.iterates, root), floor=1e-10)
    assert rate == pytest.approx(abs(math.sin(root)), abs=0.01)   # |g'(x*)| = sin(0.739) = 0.674
    with pytest.raises(RuntimeError):
        impl.fixed_point(lambda x: 2.0 * x, 0.1)              # g' = 2: diverges to inf


# ---------------------------------------------------------------- measuring convergence

def test_order_estimators_recover_synthetic_sequences():
    e_quad = np.array([1e-1, 1e-2, 1e-4, 1e-8])              # e_{k+1} = e_k^2
    npt.assert_allclose(impl.observed_convergence_order(e_quad), 2.0, atol=1e-12)
    assert impl.fitted_convergence_order(np.array([1e-1, 1e-2, 1e-4, 1e-8, 1e-16]), floor=1e-20) == pytest.approx(2.0, abs=1e-12)
    e_lin = 0.5 ** np.arange(1, 12)
    npt.assert_allclose(impl.observed_convergence_order(e_lin), 1.0, atol=1e-12)
    assert impl.observed_linear_rate(e_lin) == pytest.approx(0.5)
    # Entries below the floor are ignored; zeros are dropped by errors_from_iterates.
    assert impl.observed_convergence_order(np.array([1e-14, 1e-15, 1e-16])).size == 0
    assert impl.errors_from_iterates(np.array([1.0, 2.0, 2.0]), 2.0).size == 1
    with pytest.raises(ValueError):
        impl.observed_linear_rate(np.array([1e-20]))


def test_observed_orders_on_kepler_are_2_for_newton_1p6_for_secant_and_1_for_fixed_point():
    e, M = 0.9, math.radians(10.0)
    ref = kepler_ref(M, e)
    errs_n = impl.errors_from_iterates(impl.solve_kepler(M, e, "newton").iterates, ref)
    errs_s = impl.errors_from_iterates(impl.solve_kepler(M, e, "secant").iterates, ref)
    errs_f = impl.errors_from_iterates(impl.solve_kepler(M, e, "fixed_point").iterates, ref)
    # Fitted slope of log e_{k+1} vs log e_k: theory 2, (1+sqrt5)/2 = 1.618, 1. Tolerances reflect
    # the short sequences (5 to 12 usable errors) and the pre-asymptotic first steps.
    assert impl.fitted_convergence_order(errs_n) == pytest.approx(2.0, abs=0.2)
    assert impl.fitted_convergence_order(errs_s) == pytest.approx(GOLDEN, abs=0.2)
    assert impl.fitted_convergence_order(errs_f) == pytest.approx(1.0, abs=0.05)
    # The last two per-step Newton estimates are within 0.1 of 2.
    q = impl.observed_convergence_order(errs_n)
    assert np.all(np.abs(q[-2:] - 2.0) < 0.1)
    # Fixed-point rate is |g'(E*)| = e |cos E*|.
    assert impl.observed_linear_rate(errs_f) == pytest.approx(e * abs(math.cos(ref)), abs=0.01)
    # Cost ordering for the same tolerance: Newton fewest steps, then secant; fixed point at rate
    # e |cos E*| = 0.89 here is even slower than bisection's rate 0.5.
    steps = {m: impl.solve_kepler(M, e, m).iterates.size for m in ("newton", "secant", "fixed_point", "bisection")}
    assert steps["newton"] < steps["secant"] < steps["bisection"] < steps["fixed_point"]


# ---------------------------------------------------------------- Kepler's equation

@pytest.mark.parametrize("e", [0.0, 0.1, 0.5, 0.9, 0.99])
def test_all_methods_solve_kepler_to_the_oracle(e):
    for M in np.linspace(0.05, 2 * math.pi - 0.05, 9):
        ref = kepler_ref(M, e)
        for method in ("newton", "secant", "bisection", "fixed_point"):
            res = impl.solve_kepler(M, e, method, tol=1e-14)
            assert abs(res.root - ref) < 1e-12, (method, M)
            assert abs(impl.kepler_residual(res.root, M, e)) < 1e-13
    if e == 0.0:
        assert impl.solve_kepler(1.234, 0.0).root == pytest.approx(1.234, abs=1e-15)   # circular: E = M


def test_kepler_handles_revolutions_and_rejects_bad_input():
    M, e = 1.0, 0.6
    base = impl.solve_kepler(M, e).root
    assert impl.solve_kepler(M + 6 * math.pi, e).root == pytest.approx(base + 6 * math.pi, abs=1e-12)
    assert impl.solve_kepler(M - 2 * math.pi, e).root == pytest.approx(base - 2 * math.pi, abs=1e-12)
    assert 1.0 - e <= impl.kepler_derivative(2.0, e) <= 1.0 + e
    with pytest.raises(ValueError):
        impl.solve_kepler(M, 1.0)
    with pytest.raises(ValueError):
        impl.solve_kepler(M, -0.1)
    with pytest.raises(ValueError):
        impl.solve_kepler(M, e, method="halley")


def test_true_anomaly_limits_and_round_trip():
    for E in np.linspace(0.1, 2 * math.pi - 0.1, 7):
        assert math.remainder(impl.true_anomaly(E, 0.0) - E, 2 * math.pi) == pytest.approx(0.0, abs=1e-13)   # circular: nu = E
    assert impl.true_anomaly(math.pi, 0.7) == pytest.approx(math.pi, abs=1e-12)
    assert impl.true_anomaly(0.0, 0.7) == 0.0
    e = 0.73
    for E in np.linspace(0.2, math.pi - 0.2, 5):
        nu = impl.true_anomaly(E, e)
        assert nu > E                                        # the satellite is ahead in angle near perigee
        assert math.cos(E) == pytest.approx((e + math.cos(nu)) / (1 + e * math.cos(nu)), abs=1e-13)   # inverse relation


def test_orbit_state_conserves_energy_and_angular_momentum_and_returns_after_one_period():
    a, e = 24_396.1e3, 0.7283                                # GTO
    mu = impl.MU_EARTH
    n = impl.mean_motion(a, mu)
    T = 2 * math.pi / n
    h_expected = math.sqrt(mu * a * (1 - e * e))
    for t in RNG.uniform(0.0, 2 * T, 12):
        r, v = impl.orbit_state(a, e, t, mu)
        rn = np.linalg.norm(r)
        assert np.dot(v, v) == pytest.approx(mu * (2.0 / rn - 1.0 / a), rel=1e-12)     # vis-viva
        assert abs(r[0] * v[1] - r[1] * v[0]) == pytest.approx(h_expected, rel=1e-12)   # |r x v|
        assert impl.specific_energy(r, v, mu) == pytest.approx(-mu / (2 * a), rel=1e-12)
        assert a * (1 - e) - 1e-6 <= rn <= a * (1 + e) + 1e-6
    r0, v0 = impl.orbit_state(a, e, 0.0, mu)
    npt.assert_allclose(r0, [a * (1 - e), 0.0], atol=1e-6)
    assert np.linalg.norm(v0) == pytest.approx(math.sqrt(mu / a * (1 + e) / (1 - e)), rel=1e-12)
    rT, vT = impl.orbit_state(a, e, T, mu)
    npt.assert_allclose(rT, r0, atol=1e-3)                   # back to perigee after one period (mm)
    npt.assert_allclose(vT, v0, atol=1e-9)
    assert impl.mean_motion(6778e3) == pytest.approx(1.1314e-3, rel=1e-4)   # the P02.05 value
    with pytest.raises(ValueError):
        impl.mean_motion(-1.0)


# ---------------------------------------------------------------- Rayleigh pitot

def test_rayleigh_pitot_inverts_round_trip_and_is_monotone():
    assert impl.rayleigh_pitot_ratio(1.0) == pytest.approx(1.8929, abs=1e-4)    # the normal-shock M = 1 limit
    ms = np.linspace(1.0, 8.0, 50)
    ratios = np.array([impl.rayleigh_pitot_ratio(m) for m in ms])
    assert np.all(np.diff(ratios) > 0.0)
    for M in (1.2, 2.0, 3.5, 6.0):
        res = impl.mach_from_pitot_ratio(impl.rayleigh_pitot_ratio(M))
        assert res.root == pytest.approx(M, abs=1e-9)
    # SciPy reaches the same answer.
    target = impl.rayleigh_pitot_ratio(2.5)
    assert optimize.brentq(lambda m: impl.rayleigh_pitot_ratio(m) - target, 1.0, 10.0, xtol=1e-13) == pytest.approx(2.5, abs=1e-9)
    with pytest.raises(ValueError):
        impl.rayleigh_pitot_ratio(0.8)
    with pytest.raises(ValueError):
        impl.mach_from_pitot_ratio(1.5)


def test_main_runs_and_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "root_finding.png").exists()
