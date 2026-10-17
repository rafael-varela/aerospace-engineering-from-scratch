"""Public tests for P03.03 Interpolation & Splines. They assert reproduction of data, continuity, error bounds and
convergence orders, with scipy.interpolate as the oracle and the analytic ISA as the truth.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy import interpolate as si

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(303)


# ---------------------------------------------------------------- the ISA reference (given)

def test_isa_reproduces_the_defining_values():
    T, p, rho = impl.isa_atmosphere(0.0)
    assert (T, p) == (288.15, 101325.0) and rho == pytest.approx(1.2250, abs=5e-5)
    T, p, rho = impl.isa_atmosphere(11000.0)
    assert T == pytest.approx(216.65) and p == pytest.approx(22632.0, rel=2e-5) and rho == pytest.approx(0.3639, abs=5e-5)
    T, p, _ = impl.isa_atmosphere(20000.0)
    assert T == pytest.approx(216.65) and p == pytest.approx(5474.9, rel=2e-4)
    T, p, _ = impl.isa_atmosphere(47000.0)
    assert T == pytest.approx(270.65) and p == pytest.approx(110.91, rel=2e-4)
    table = impl.isa_table(2000.0)
    assert table.shape == (44, 4) and table[0, 0] == 0.0 and table[-1, 0] == 86000.0
    with pytest.raises(ValueError):
        impl.isa_atmosphere(-1.0)


# ---------------------------------------------------------------- polynomial interpolation

def test_lagrange_reproduces_nodes_matches_scipy_and_is_exact_for_polynomials():
    xs = np.array([0.0, 0.5, 1.2, 2.0, 3.1])
    ys = RNG.standard_normal(5)
    npt.assert_allclose(impl.lagrange_interpolate(xs, ys, xs), ys, atol=1e-13)
    x = np.linspace(0.0, 3.1, 50)
    ref = si.BarycentricInterpolator(xs, ys)(x)
    npt.assert_allclose(impl.lagrange_interpolate(xs, ys, x), ref, atol=1e-12)
    # Degree-4 polynomial through 5 nodes is recovered exactly.
    poly = lambda t: 2 - t + 0.5 * t ** 2 - 0.1 * t ** 4                         # noqa: E731
    npt.assert_allclose(impl.lagrange_interpolate(xs, poly(xs), x), poly(x), atol=1e-12)
    assert impl.lagrange_interpolate(xs, ys, 1.0) == pytest.approx(float(si.BarycentricInterpolator(xs, ys)(1.0)))
    with pytest.raises(ValueError):
        impl.lagrange_interpolate(xs, ys, 3.5)                                     # no extrapolation
    with pytest.raises(ValueError):
        impl.lagrange_interpolate(np.array([0.0, 1.0, 1.0]), np.ones(3), 0.5)      # repeated node


def test_newton_form_equals_lagrange_form_and_adds_nodes_incrementally():
    xs = np.array([-1.0, -0.3, 0.2, 0.9, 1.5, 2.2])
    ys = np.sin(xs)
    c = impl.newton_divided_differences(xs, ys)
    x = np.linspace(-1.0, 2.2, 40)
    npt.assert_allclose(impl.newton_evaluate(xs, c, x), impl.lagrange_interpolate(xs, ys, x), atol=1e-12)
    # The first k coefficients do not change when more nodes are appended.
    c5 = impl.newton_divided_differences(xs[:5], ys[:5])
    npt.assert_allclose(c[:5], c5, atol=1e-12)
    assert c[0] == ys[0] and c[1] == pytest.approx((ys[1] - ys[0]) / (xs[1] - xs[0]))


def test_runge_phenomenon_grows_on_equispaced_nodes_and_vanishes_on_chebyshev_nodes():
    xf = np.linspace(-1.0, 1.0, 801)
    errs_eq, errs_ch = [], []
    for n in (5, 9, 13, 17, 21):
        xe = np.linspace(-1.0, 1.0, n)
        errs_eq.append(np.max(np.abs(impl.lagrange_interpolate(xe, impl.runge(xe), xf) - impl.runge(xf))))
        xc = impl.chebyshev_nodes(n, -1.0, 1.0)
        inside = xf[(xf >= xc[0]) & (xf <= xc[-1])]
        errs_ch.append(np.max(np.abs(impl.lagrange_interpolate(xc, impl.runge(xc), inside) - impl.runge(inside))))
    assert np.all(np.diff(errs_eq) > 0) and errs_eq[-1] > 10.0           # diverges: error 60 at n = 21
    assert np.all(np.diff(errs_ch) < 0) and errs_ch[-1] < 0.02           # converges
    # Chebyshev nodes: sorted, inside (a, b), symmetric, the right count and formula.
    xc = impl.chebyshev_nodes(7, 2.0, 4.0)
    assert xc.size == 7 and np.all(np.diff(xc) > 0) and 2.0 < xc[0] and xc[-1] < 4.0
    npt.assert_allclose(xc + xc[::-1], 6.0, atol=1e-12)
    npt.assert_allclose(xc, np.sort(3.0 + np.cos((2 * np.arange(7) + 1) * np.pi / 14)), atol=1e-12)
    with pytest.raises(ValueError):
        impl.chebyshev_nodes(0, 0.0, 1.0)


# ---------------------------------------------------------------- piecewise interpolation

def test_piecewise_linear_matches_numpy_interp_and_log_linear_is_exact_for_exponentials():
    xs = np.sort(RNG.uniform(0.0, 10.0, 12))
    ys = RNG.standard_normal(12)
    x = np.linspace(xs[0], xs[-1], 200)
    npt.assert_allclose(impl.piecewise_linear(xs, ys, x), np.interp(x, xs, ys), atol=1e-13)
    npt.assert_allclose(impl.piecewise_linear(xs, ys, xs), ys, atol=1e-13)
    # Log-linear is exact for A exp(B x): the ISA pressure in an isothermal layer.
    h = np.linspace(11000.0, 20000.0, 10)
    p = impl.isa_atmosphere(h)[1]
    hf = np.linspace(11000.0, 20000.0, 300)
    npt.assert_allclose(impl.log_linear_interpolate(h, p, hf), impl.isa_atmosphere(hf)[1], rtol=1e-12)
    assert np.max(np.abs(impl.piecewise_linear(h, p, hf) - impl.isa_atmosphere(hf)[1])) > 10.0   # linear is not (Pa)
    with pytest.raises(ValueError):
        impl.log_linear_interpolate(xs, ys, x)                                     # negative values
    with pytest.raises(ValueError):
        impl.piecewise_linear(xs, ys, xs[-1] + 1.0)


@pytest.mark.parametrize("bc", ["natural", "clamped", "not-a-knot"])
def test_cubic_spline_matches_scipy_and_is_c2(bc):
    xs = np.sort(RNG.uniform(0.0, 5.0, 9))
    ys = np.cos(xs) + 0.2 * RNG.standard_normal(9)
    slopes = (0.7, -1.3)
    sp = impl.cubic_spline(xs, ys, bc, slopes)
    scipy_bc = {"natural": "natural", "clamped": ((1, slopes[0]), (1, slopes[1])), "not-a-knot": "not-a-knot"}[bc]
    ref = si.CubicSpline(xs, ys, bc_type=scipy_bc)
    x = np.linspace(xs[0], xs[-1], 400)
    npt.assert_allclose(impl.spline_evaluate(sp, x), ref(x), atol=1e-10)
    npt.assert_allclose(impl.spline_evaluate(sp, x, 1), ref(x, 1), atol=1e-9)
    npt.assert_allclose(impl.spline_evaluate(sp, x, 2), ref(x, 2), atol=1e-8)
    npt.assert_allclose(impl.spline_evaluate(sp, xs), ys, atol=1e-12)              # interpolates
    # C^2: first and second derivatives agree from the left and the right at every interior knot.
    d = 1e-9
    for xk in xs[1:-1]:
        for deriv in (1, 2):
            left = impl.spline_evaluate(sp, xk - d, deriv)
            right = impl.spline_evaluate(sp, xk + d, deriv)
            assert left == pytest.approx(right, rel=1e-6, abs=1e-6)
    if bc == "natural":
        assert sp.M[0] == 0.0 and sp.M[-1] == 0.0
        assert impl.spline_evaluate(sp, xs[0], 2) == pytest.approx(0.0, abs=1e-12)
    if bc == "clamped":
        assert impl.spline_evaluate(sp, xs[0], 1) == pytest.approx(slopes[0], abs=1e-10)
        assert impl.spline_evaluate(sp, xs[-1], 1) == pytest.approx(slopes[1], abs=1e-10)
    if bc == "not-a-knot":
        # Third derivative continuous at x_1: the first two pieces are one cubic.
        h0, h1 = xs[1] - xs[0], xs[2] - xs[1]
        assert (sp.M[1] - sp.M[0]) / h0 == pytest.approx((sp.M[2] - sp.M[1]) / h1, abs=1e-8)
    with pytest.raises(ValueError):
        impl.spline_evaluate(sp, xs[-1] + 0.1)
    with pytest.raises(ValueError):
        impl.cubic_spline(xs, ys, "periodic")


def test_spline_reproduces_a_cubic_exactly_and_clamped_beats_natural():
    xs = np.linspace(0.0, 2.0, 6)
    cubic = lambda t: 1.0 - 2.0 * t + 0.5 * t ** 2 + 0.3 * t ** 3                 # noqa: E731
    dcubic = lambda t: -2.0 + t + 0.9 * t ** 2                                     # noqa: E731
    x = np.linspace(0.0, 2.0, 101)
    sp_c = impl.cubic_spline(xs, cubic(xs), "clamped", (dcubic(0.0), dcubic(2.0)))
    sp_k = impl.cubic_spline(xs, cubic(xs), "not-a-knot")
    sp_n = impl.cubic_spline(xs, cubic(xs), "natural")
    npt.assert_allclose(impl.spline_evaluate(sp_c, x), cubic(x), atol=1e-12)        # exact
    npt.assert_allclose(impl.spline_evaluate(sp_k, x), cubic(x), atol=1e-12)        # exact
    assert np.max(np.abs(impl.spline_evaluate(sp_n, x) - cubic(x))) > 1e-3          # natural is not (f'' != 0 at ends)


def test_error_bounds_hold_for_linear_and_clamped_spline_on_a_smooth_function():
    f = np.sin
    for n in (10, 20, 40):
        h = np.pi / n
        e_lin = impl.max_interpolation_error(f, 0.0, np.pi, n, "linear")
        e_cl = impl.max_interpolation_error(f, 0.0, np.pi, n, "clamped")
        assert e_lin <= impl.linear_error_bound(h, 1.0) and e_lin > 0.5 * impl.linear_error_bound(h, 1.0)
        assert e_cl <= impl.spline_error_bound(h, 1.0)


def test_observed_convergence_orders_linear_2_natural_2_at_the_ends_not_a_knot_and_clamped_4():
    f = lambda t: np.exp(-t / 2.0) * np.cos(2.0 * t)                               # noqa: E731
    ns = np.array([20, 40, 80, 160])
    hs = 3.0 / ns
    orders = {}
    for m in ("linear", "natural", "not-a-knot", "clamped"):
        errs = [impl.max_interpolation_error(f, 0.0, 3.0, int(n), m) for n in ns]
        orders[m] = impl.observed_order(hs, errs)
    assert orders["linear"] == pytest.approx(2.0, abs=0.15)
    assert orders["natural"] == pytest.approx(2.0, abs=0.3)      # the free ends force f'' = 0 where it is not
    assert orders["not-a-knot"] == pytest.approx(4.0, abs=0.4)   # the end pieces converge a little slower in the max norm
    assert orders["clamped"] == pytest.approx(4.0, abs=0.3)
    # A function with f'' = 0 at both ends lets the natural spline reach order 4 too.
    errs = [impl.max_interpolation_error(np.sin, 0.0, np.pi, int(n), "natural") for n in ns]
    assert impl.observed_order(np.pi / ns, errs) == pytest.approx(4.0, abs=0.3)
    assert impl.observed_order([1.0, 0.5, 0.25], [1.0, 0.25, 0.0625]) == pytest.approx(2.0, abs=1e-12)
    with pytest.raises(ValueError):
        impl.observed_order([1.0], [1.0])


# ---------------------------------------------------------------- the atmosphere table

def test_isa_table_interpolation_errors_fall_with_spacing_and_log_linear_wins_for_pressure():
    hf = np.linspace(500.0, 85500.0, 1701)
    _, p_true, rho_true = impl.isa_atmosphere(hf)
    for col, truth in (("p", p_true), ("rho", rho_true)):
        e2 = {m: np.max(np.abs(impl.interpolate_isa(impl.isa_table(2000.0), hf, col, m) - truth) / truth)
              for m in ("linear", "loglinear", "natural")}
        e1 = {m: np.max(np.abs(impl.interpolate_isa(impl.isa_table(1000.0), hf, col, m) - truth) / truth)
              for m in ("linear", "loglinear", "natural")}
        for m in e2:
            assert e1[m] < e2[m] / 2.5                                  # at least ~order 1.3 improvement
        if col == "p":
            assert e2["loglinear"] < e2["linear"] / 4                   # pressure is near-exponential everywhere
    # Density in the troposphere interior (no layer boundary inside a table interval): log-linear beats linear.
    ht = np.linspace(500.0, 9500.0, 901)
    rho_t = impl.isa_atmosphere(ht)[2]
    e_lin = np.max(np.abs(impl.interpolate_isa(impl.isa_table(2000.0), ht, "rho", "linear") - rho_t) / rho_t)
    e_log = np.max(np.abs(impl.interpolate_isa(impl.isa_table(2000.0), ht, "rho", "loglinear") - rho_t) / rho_t)
    assert e_log < e_lin / 2
    # Pressure from the 1 km table: log-linear is exact inside the isothermal layer 11-20 km.
    hi = np.linspace(11000.0, 20000.0, 91)
    npt.assert_allclose(impl.interpolate_isa(impl.isa_table(1000.0), hi, "p", "loglinear"), impl.isa_atmosphere(hi)[1], rtol=1e-12)
    # Temperature is piecewise linear with knots on the table: linear interpolation is exact.
    npt.assert_allclose(impl.interpolate_isa(impl.isa_table(1000.0), hf, "T", "linear"), impl.isa_atmosphere(hf)[0], atol=1e-10)
    with pytest.raises(ValueError):
        impl.interpolate_isa(impl.isa_table(1000.0), 5000.0, "rho", "quintic")


def test_spline_overshoots_at_the_tropopause_kink_while_linear_does_not():
    table = impl.isa_table(1000.0)
    hk = np.linspace(9000.0, 13000.0, 401)
    T_true = impl.isa_atmosphere(hk)[0]
    err_sp = np.max(np.abs(impl.interpolate_isa(table, hk, "T", "natural") - T_true))
    err_li = np.max(np.abs(impl.interpolate_isa(table, hk, "T", "linear") - T_true))
    assert 0.2 < err_sp < 2.0 and err_li < 1e-9                        # the spline rings around the kink by ~0.5 K
    # The ringing decays away from the kink: error at 5 km is far smaller than at 10.5 km.
    e5 = abs(impl.interpolate_isa(table, 5500.0, "T", "natural") - impl.isa_atmosphere(5500.0)[0])
    e105 = abs(impl.interpolate_isa(table, 10500.0, "T", "natural") - impl.isa_atmosphere(10500.0)[0])
    assert e5 < 0.05 * e105


def test_main_runs_and_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "atmosphere_interpolation.png").exists()
