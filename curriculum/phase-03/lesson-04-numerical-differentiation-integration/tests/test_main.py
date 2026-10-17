"""Public tests for P03.04 Numerical Differentiation & Integration. They assert exactness classes, error bounds,
convergence orders and the closed-form re-entry results, with scipy.integrate and numpy.polynomial.legendre as oracles.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import math

import numpy as np
import numpy.testing as npt
import pytest
from numpy.polynomial import legendre
from scipy import integrate

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(304)


# ---------------------------------------------------------------- finite differences

def test_stencils_have_their_theoretical_orders_and_leading_constants():
    f, x = math.sin, 1.0
    hs = np.array([1e-1, 1e-2, 1e-3])
    ef = [abs(impl.forward_diff(f, x, h) - math.cos(x)) for h in hs]
    ec = [abs(impl.central_diff(f, x, h) - math.cos(x)) for h in hs]
    e2 = [abs(impl.second_diff(f, x, h) + math.sin(x)) for h in hs]
    assert impl.observed_order(hs, ef) == pytest.approx(1.0, abs=0.02)
    assert impl.observed_order(hs, ec) == pytest.approx(2.0, abs=0.02)
    assert impl.observed_order(hs, e2) == pytest.approx(2.0, abs=0.02)
    # Leading error constants: h f''/2, h^2 f'''/6, h^2 f''''/12 (f'' = -sin, f''' = -cos, f'''' = sin at 1).
    h = 1e-3
    assert ef[2] == pytest.approx(h * math.sin(x) / 2, rel=1e-2)
    assert ec[2] == pytest.approx(h * h * math.cos(x) / 6, rel=1e-2)
    assert e2[2] == pytest.approx(h * h * math.sin(x) / 12, rel=1e-2)


def test_stencils_are_exact_on_their_polynomial_classes():
    lin = lambda t: 3.0 * t - 2.0                            # noqa: E731
    quad = lambda t: 1.5 * t * t - t + 4.0                    # noqa: E731
    for h in (0.5, 0.1):
        assert impl.forward_diff(lin, 0.7, h) == pytest.approx(3.0, abs=1e-12)
        assert impl.central_diff(quad, 0.7, h) == pytest.approx(3.0 * 0.7 - 1.0, abs=1e-12)   # central is exact for quadratics
        assert impl.second_diff(quad, 0.7, h) == pytest.approx(3.0, abs=1e-10)


def test_richardson_raises_the_order_and_recovers_the_exact_derivative_of_a_cubic():
    f, x = math.sin, 1.0
    hs = np.array([1e-1, 2e-2])
    er = [abs(impl.richardson(lambda s: impl.central_diff(f, x, s), h, 2) - math.cos(x)) for h in hs]
    assert impl.observed_order(hs, er) == pytest.approx(4.0, abs=0.1)
    ef = [abs(impl.richardson(lambda s: impl.forward_diff(f, x, s), h, 1) - math.cos(x)) for h in hs]
    assert impl.observed_order(hs, ef) == pytest.approx(2.0, abs=0.1)
    cubic = lambda t: t ** 3 - 2 * t                          # noqa: E731
    assert impl.richardson(lambda s: impl.central_diff(cubic, 0.5, s), 0.3, 2) == pytest.approx(3 * 0.25 - 2, abs=1e-12)
    assert impl.richardson(lambda s: s * s * 7.0 + 1.0, 0.1, 2) == pytest.approx(1.0, abs=1e-14)   # pure h^2 error is removed exactly
    with pytest.raises(ValueError):
        impl.richardson(lambda s: s, 0.1, 0)


def test_central_difference_optimal_step_and_rounding_floor():
    hs = 10.0 ** -np.arange(1, 13, dtype=float)
    ec = np.array([abs(impl.central_diff(math.sin, 1.0, h) - math.cos(1.0)) for h in hs])
    h_star = impl.optimal_step_central(1.0, 1.0)
    assert h_star == pytest.approx((3 * 2.0 ** -52) ** (1 / 3))
    best = hs[np.argmin(ec)]
    assert 0.1 * h_star <= best <= 10 * h_star
    assert ec.min() < 1e-10                                   # two thirds of the digits
    assert ec[-1] > ec.min() * 100                            # rounding regime at h = 1e-12
    with pytest.raises(ValueError):
        impl.optimal_step_central(1.0, 0.0)


# ---------------------------------------------------------------- Newton-Cotes and Romberg

def test_trapezoid_and_simpson_exactness_classes_and_orders():
    lin = lambda x: 2.0 * x + 1.0                             # noqa: E731
    cubic = lambda x: x ** 3 - x + 0.5                        # noqa: E731
    assert impl.trapezoid(lin, 0.0, 3.0, 1) == pytest.approx(12.0, abs=1e-12)
    assert impl.simpson(cubic, -1.0, 2.0, 2) == pytest.approx(15.0 / 4 - 3.0 / 2 + 1.5, abs=1e-12)   # exact for cubics with 2 panels
    assert impl.simpson(lambda x: x ** 4, 0.0, 1.0, 2) != pytest.approx(0.2, abs=1e-6)            # not for quartics
    f = lambda x: np.exp(-x) * np.cos(3.0 * x)                # noqa: E731
    exact = (1.0 - math.exp(-2.0) * (math.cos(6.0) - 3.0 * math.sin(6.0))) / 10.0
    ns = np.array([8, 16, 32, 64, 128])
    et = [abs(impl.trapezoid(f, 0.0, 2.0, int(n)) - exact) for n in ns]
    es = [abs(impl.simpson(f, 0.0, 2.0, int(n)) - exact) for n in ns]
    assert impl.observed_order(1.0 / ns, et) == pytest.approx(2.0, abs=0.1)
    assert impl.observed_order(1.0 / ns, es) == pytest.approx(4.0, abs=0.2)
    # Euler-Maclaurin: the trapezoid error is -h^2 (f'(b) - f'(a)) / 12 + O(h^4), with f' = -e^-x (cos 3x + 3 sin 3x).
    h = 2.0 / 128
    fp = lambda x: -math.exp(-x) * (math.cos(3 * x) + 3 * math.sin(3 * x))   # noqa: E731
    assert et[-1] == pytest.approx(abs(h * h * (fp(2.0) - fp(0.0)) / 12), rel=0.02)
    # Against SciPy's sample-based rules on the same grid.
    x = np.linspace(0.0, 2.0, 65)
    assert impl.trapezoid(f, 0.0, 2.0, 64) == pytest.approx(integrate.trapezoid(f(x), x), rel=1e-13)
    assert impl.simpson(f, 0.0, 2.0, 64) == pytest.approx(integrate.simpson(f(x), x=x), rel=1e-13)
    with pytest.raises(ValueError):
        impl.simpson(f, 0.0, 1.0, 3)
    with pytest.raises(ValueError):
        impl.trapezoid(f, 0.0, 1.0, 0)


def test_trapezoid_is_spectrally_accurate_on_a_periodic_integrand():
    # int_0^{2pi} exp(cos x) dx = 2 pi I_0(1); the trapezoid error falls exponentially with n, not as h^2.
    from scipy.special import i0
    exact = 2 * math.pi * i0(1.0)
    f = lambda x: np.exp(np.cos(x))                           # noqa: E731
    e8 = abs(impl.trapezoid(f, 0.0, 2 * math.pi, 8) - exact)
    e16 = abs(impl.trapezoid(f, 0.0, 2 * math.pi, 16) - exact)
    assert e8 < 1e-5 and e16 < 1e-14 and e16 < e8 / 1e6      # far beyond the factor 4 of a second-order rule


def test_romberg_table_reproduces_simpson_and_converges_fast():
    f = lambda x: 1.0 / (1.0 + x * x)                         # noqa: E731
    R = impl.romberg(f, 0.0, 1.0, 6)
    exact = math.pi / 4
    npt.assert_allclose(R[:, 0], [impl.trapezoid(f, 0.0, 1.0, 2 ** i) for i in range(6)], rtol=1e-14)
    for i in range(1, 6):
        assert R[i, 1] == pytest.approx(impl.simpson(f, 0.0, 1.0, 2 ** i), rel=1e-13)
    assert abs(R[5, 5] - exact) < 1e-10 and abs(R[5, 0] - exact) > 1e-5
    assert np.all(np.abs(np.diag(R)[1:] - exact) <= np.abs(np.diag(R)[:-1] - exact) + 1e-16)
    with pytest.raises(ValueError):
        impl.romberg(f, 0.0, 1.0, 0)


# ---------------------------------------------------------------- Gauss-Legendre

@pytest.mark.parametrize("n", [1, 2, 3, 5, 8, 16])
def test_gauss_legendre_nodes_match_numpy_and_integrate_degree_2n_minus_1_exactly(n):
    x, w = impl.gauss_legendre_nodes(n)
    x_ref, w_ref = legendre.leggauss(n)
    npt.assert_allclose(x, x_ref, atol=1e-14)
    npt.assert_allclose(w, w_ref, atol=1e-14)
    assert w.sum() == pytest.approx(2.0, abs=1e-14) and np.all(w > 0)
    npt.assert_allclose(x, -x[::-1], atol=1e-14)             # symmetric nodes
    # Exact for every monomial up to degree 2n - 1 on [-1, 1], not for degree 2n.
    for k in range(2 * n):
        exact = 0.0 if k % 2 else 2.0 / (k + 1)
        assert impl.gauss_legendre(lambda t, k=k: t ** k, -1.0, 1.0, n) == pytest.approx(exact, abs=1e-13)
    k = 2 * n
    assert abs(impl.gauss_legendre(lambda t: t ** k, -1.0, 1.0, n) - 2.0 / (k + 1)) > 1e-11


def test_gauss_legendre_converges_spectrally_and_maps_intervals():
    f = lambda x: np.exp(-x) * np.cos(3.0 * x)                # noqa: E731
    exact = (1.0 - math.exp(-2.0) * (math.cos(6.0) - 3.0 * math.sin(6.0))) / 10.0
    e4 = abs(impl.gauss_legendre(f, 0.0, 2.0, 4) - exact)
    e8 = abs(impl.gauss_legendre(f, 0.0, 2.0, 8) - exact)
    e16 = abs(impl.gauss_legendre(f, 0.0, 2.0, 16) - exact)
    assert e8 < e4 / 1e5 and e16 < 1e-14                      # doubling n gains > 5 digits: spectral, not algebraic
    assert impl.gauss_legendre(lambda x: x ** 2, 1.0, 4.0, 2) == pytest.approx(21.0, abs=1e-12)   # mapped interval, exact
    with pytest.raises(ValueError):
        impl.gauss_legendre_nodes(0)


# ---------------------------------------------------------------- adaptive Simpson

def test_adaptive_simpson_meets_its_tolerance_and_spends_evaluations_where_needed():
    f = lambda x: math.exp(-x) * math.cos(3.0 * x)            # noqa: E731
    exact = (1.0 - math.exp(-2.0) * (math.cos(6.0) - 3.0 * math.sin(6.0))) / 10.0
    for tol in (1e-6, 1e-10):
        q, n = impl.adaptive_simpson(f, 0.0, 2.0, tol)
        assert abs(q - exact) < tol
    q6, n6 = impl.adaptive_simpson(f, 0.0, 2.0, 1e-6)
    q10, n10 = impl.adaptive_simpson(f, 0.0, 2.0, 1e-10)
    assert n10 > n6                                           # tighter tolerance costs more
    # A nearly singular peak: adaptivity concentrates evaluations and still matches quad.
    g = lambda x: 1.0 / (1e-3 + (x - 0.3) ** 2)               # noqa: E731
    ref = integrate.quad(g, 0.0, 1.0, epsabs=1e-12, epsrel=1e-12)[0]
    q, n = impl.adaptive_simpson(g, 0.0, 1.0, 1e-8)
    assert abs(q - ref) < 1e-6 * ref and n < 5000
    with pytest.raises(RuntimeError):
        impl.adaptive_simpson(lambda x: abs(x - 0.5) ** -0.9 if x != 0.5 else 0.0, 0.0, 1.0, 1e-14, max_depth=8)


def test_observed_order_recovers_a_power_law():
    hs = np.array([0.4, 0.2, 0.1, 0.05])
    assert impl.observed_order(hs, 2.0 * hs ** 4) == pytest.approx(4.0, abs=1e-12)
    with pytest.raises(ValueError):
        impl.observed_order(hs, -hs)


# ---------------------------------------------------------------- the re-entry model

def test_entry_profile_limits_and_closed_form_peaks():
    p = impl.EntryParams()
    B = impl.entry_B(p)
    assert B == pytest.approx(p.rho0 * p.H / (2 * p.beta * math.sin(p.gamma)))
    assert impl.entry_velocity(300e3, p) == pytest.approx(p.V_e, rel=1e-12)     # no air, no slowing
    assert impl.entry_velocity(0.0, p) == pytest.approx(p.V_e * math.exp(-B))   # essentially stopped at sea level
    h_a, a_max = impl.peak_deceleration(p)
    assert h_a == pytest.approx(p.H * math.log(2 * B))
    assert a_max == pytest.approx(p.V_e ** 2 * math.sin(p.gamma) / (2 * math.e * p.H))
    assert impl.deceleration_analytic(h_a, p) == pytest.approx(a_max, rel=1e-12)
    # The analytic peak really is the maximum of the analytic profile.
    hh = np.linspace(10e3, 100e3, 9001)
    assert hh[np.argmax(impl.deceleration_analytic(hh, p))] == pytest.approx(h_a, abs=20.0)
    assert hh[np.argmax(impl.heating_rate(hh, p))] == pytest.approx(impl.peak_heating_altitude(p), abs=20.0)
    assert 5.0 < a_max / impl.G0 < 12.0                       # a ballistic capsule-class deceleration


def test_numerical_deceleration_matches_the_analytic_one_and_converges_at_second_order():
    p = impl.EntryParams()
    h = 50e3
    a_ref = impl.deceleration_analytic(h, p)
    errs = [abs(impl.deceleration_from_velocity(h, p, dh) - a_ref) for dh in (400.0, 200.0, 100.0)]
    assert impl.observed_order([400.0, 200.0, 100.0], errs) == pytest.approx(2.0, abs=0.1)
    assert abs(impl.deceleration_from_velocity(h, p, 10.0) - a_ref) < 1e-6 * a_ref


def test_heat_load_numeric_matches_the_erf_closed_form_with_the_right_orders():
    p = impl.EntryParams()
    Q = impl.heat_load_analytic(120e3, p)
    # Independent oracle: scipy.integrate.quad on the same integrand.
    q_ref = integrate.quad(lambda h: impl.heat_load_integrand(h, p), 0.0, 120e3, epsrel=1e-12)[0]
    assert Q == pytest.approx(q_ref, rel=1e-9)
    ns = np.array([32, 64, 128, 256])
    et = [abs(impl.heat_load_numeric(120e3, p, "trapezoid", int(n)) - Q) / Q for n in ns]
    es = [abs(impl.heat_load_numeric(120e3, p, "simpson", int(n)) - Q) / Q for n in ns]
    assert impl.observed_order(1.0 / ns, et) == pytest.approx(2.0, abs=0.15)
    # The integrand is a bump of width ~H = 7.2 km on a 120 km range, so Simpson is still pre-asymptotic at
    # 32-128 panels (its h^4 f'''' term is not yet the leading one) and its apparent order exceeds 4 before
    # the error reaches the rounding floor at 256 panels: assert at least fourth order, not exactly fourth.
    assert impl.observed_order(1.0 / ns[:3], es[:3]) > 3.5
    assert es[-1] < 1e-9 and et[-1] < 1e-5
    assert abs(impl.heat_load_numeric(120e3, p, "gauss", 64) - Q) / Q < 1e-10
    # Heat load is monotone in the starting altitude and dominated by the band below 80 km.
    assert impl.heat_load_analytic(80e3, p) > 0.9 * Q and impl.heat_load_analytic(120e3, p) > impl.heat_load_analytic(80e3, p)
    with pytest.raises(ValueError):
        impl.heat_load_numeric(120e3, p, "monte-carlo", 8)


def test_main_runs_and_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "reentry_heating.png").exists()
