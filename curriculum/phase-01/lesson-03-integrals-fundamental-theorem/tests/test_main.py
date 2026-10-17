"""Public tests for P01.03 Integrals & the Fundamental Theorem. They assert properties, not magic numbers.

Run `pytest curriculum/phase-01/lesson-03-integrals-fundamental-theorem` against your starter/main.py, or with
`AEFS_IMPL=solution` against the reference solution. Hidden validation cases live in
assessments/phase-01/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
import sympy as sp

from conftest import load_impl

impl = load_impl(__file__)

T, C = 0.12, 1.0
S = impl.model_symbols()
t_s, c_s, x_s = S["t"], S["c"], S["x"]


# ---------------------------------------------------------------- geometry

def test_thickness_distribution_shape():
    assert impl.naca4_thickness(0.0, T, C) == 0.0                                   # sharp leading edge point
    x = np.linspace(0.0, C, 2001)
    yt = impl.naca4_thickness(x, T, C)
    assert yt.shape == x.shape and np.all(yt >= 0.0)
    i = int(np.argmax(yt))
    assert 2.0 * yt[i] == pytest.approx(T, rel=2e-3)                                # max thickness is t c ...
    assert x[i] == pytest.approx(0.30 * C, abs=0.01)                                # ... at 30 % chord
    assert impl.naca4_thickness(C, T, C, closed_te=True) == pytest.approx(0.0, abs=1e-12)
    assert impl.naca4_thickness(C, T, C) > 0.0                                      # open trailing edge
    npt.assert_allclose(impl.naca4_thickness(x * 2.0, T, 2.0 * C), 2.0 * yt, rtol=1e-12)   # scales with chord
    npt.assert_allclose(impl.naca4_thickness(x, 2.0 * T, C), 2.0 * yt, rtol=1e-12)         # linear in t
    for bad in ((-0.1, T, C), (1.1 * C, T, C), (0.5, 0.0, C)):
        with pytest.raises(ValueError):
            impl.naca4_thickness(*bad)


def test_camber_line_properties():
    x = np.linspace(0.0, C, 1001)
    assert np.all(impl.naca4_camber(x, 0.0, 0.4, C) == 0.0)
    yc = impl.naca4_camber(x, 0.02, 0.4, C)
    assert yc[0] == pytest.approx(0.0) and yc[-1] == pytest.approx(0.0, abs=1e-12)
    assert np.max(yc) == pytest.approx(0.02 * C, rel=1e-6)                           # max camber m c ...
    assert x[int(np.argmax(yc))] == pytest.approx(0.4 * C, abs=1e-3)                 # ... at p c
    left = impl.naca4_camber(0.4 * C - 1e-9, 0.02, 0.4, C)
    right = impl.naca4_camber(0.4 * C + 1e-9, 0.02, 0.4, C)
    assert left == pytest.approx(right, abs=1e-8)                                    # continuous at the join
    with pytest.raises(ValueError):
        impl.naca4_camber(0.3, 0.02, 1.0, C)


# ---------------------------------------------------------------- Riemann sums

def test_riemann_sum_bounds_and_exactness():
    f = np.exp
    exact = np.e - 1.0
    left, right, mid = (impl.riemann_sum(f, 0.0, 1.0, 50, r) for r in ("left", "right", "midpoint"))
    assert left < exact < right                                                      # increasing f: left under, right over
    assert abs(mid - exact) < abs(left - exact) and abs(mid - exact) < abs(right - exact)
    assert impl.riemann_sum(lambda x: 3.0 * x + 1.0, 0.0, 2.0, 7, "midpoint") == pytest.approx(8.0)   # exact for linear
    assert impl.riemann_sum(f, 1.0, 0.0, 50, "midpoint") == pytest.approx(-mid)      # reversed limits flip the sign
    with pytest.raises(ValueError):
        impl.riemann_sum(f, 0.0, 1.0, 0)
    with pytest.raises(ValueError):
        impl.riemann_sum(f, 0.0, 1.0, 10, "simpson")


def test_riemann_convergence_orders_on_a_smooth_integrand():
    ns = np.array([20, 40, 80, 160, 320])
    exact = np.sin(1.0)
    assert impl.observed_order(ns, impl.convergence_study(np.cos, exact, 0.0, 1.0, ns, "left")) == pytest.approx(1.0, abs=0.05)
    assert impl.observed_order(ns, impl.convergence_study(np.cos, exact, 0.0, 1.0, ns, "midpoint")) == pytest.approx(2.0, abs=0.05)


def test_square_root_singularity_degrades_midpoint_to_order_one_and_a_half():
    ns = np.array([40, 80, 160, 320, 640, 1280])
    f = lambda x: 2.0 * impl.naca4_thickness(x, T, C)
    errs = impl.convergence_study(f, impl.section_area_exact(T, C), 0.0, C, ns, "midpoint")
    assert np.all(np.diff(errs) < 0.0)
    assert impl.observed_order(ns, errs) == pytest.approx(1.5, abs=0.1)
    # The substitution x = c u^2 restores second order.
    g = lambda u: 2.0 * impl.naca4_thickness(C * u**2, T, C) * 2.0 * C * u
    errs_u = impl.convergence_study(g, impl.section_area_exact(T, C), 0.0, 1.0, ns, "midpoint")
    assert impl.observed_order(ns, errs_u) == pytest.approx(2.0, abs=0.1)
    assert np.all(errs_u < errs)


def test_observed_order_recovers_synthetic_data_and_rejects_bad_input():
    ns = np.array([10, 20, 40, 80])
    assert impl.observed_order(ns, 2.0 * ns ** -1.5) == pytest.approx(1.5, abs=1e-9)
    with pytest.raises(ValueError):
        impl.observed_order(ns, np.array([1.0, 0.0, 1.0, 1.0]))


# ---------------------------------------------------------------- the Fundamental Theorem

def test_ftc_part1_derivative_of_accumulated_area_is_the_integrand():
    f = lambda x: 2.0 * impl.naca4_thickness(x, T, C)
    big_f = lambda x: impl.cumulative_integral(f, 0.0, x, n=2000)
    assert big_f(np.array([0.0]))[0] == 0.0
    xs = np.array([0.2, 0.5, 0.8])
    for h in (1e-2, 5e-3):
        d = (big_f(xs + h) - big_f(xs - h)) / (2.0 * h)
        npt.assert_allclose(d, f(xs), rtol=2e-3)
    assert np.all(np.diff(big_f(np.linspace(0.0, C, 11))) > 0.0)                    # f > 0: F increasing


def test_ftc_part2_antiderivative_matches_riemann_sum_and_symbolic_area():
    poly = x_s**3 - 2 * x_s + 1
    exact = impl.definite_from_antiderivative(poly, x_s, 0.0, 2.0)
    assert exact == pytest.approx(4.0 - 4.0 + 2.0)
    assert impl.riemann_sum(lambda x: x**3 - 2 * x + 1, 0.0, 2.0, 4000, "midpoint") == pytest.approx(exact, rel=1e-6)
    a_sym = impl.symbolic_section_area()
    assert sp.simplify(a_sym / (t_s * c_s**2)).is_number                             # A = (number) t c^2
    assert float(a_sym.subs({t_s: T, c_s: C})) == pytest.approx(impl.section_area_exact(T, C), rel=1e-12)
    assert float(impl.symbolic_section_area(True).subs({t_s: T, c_s: C})) == pytest.approx(impl.section_area_exact(T, C, True), rel=1e-12)


def test_area_numeric_agrees_with_closed_form_and_scales():
    a_exact = impl.section_area_exact(T, C)
    assert 0.08 < a_exact / (T * C**2) * T < 0.09                                   # NACA 0012: about 0.082 c^2
    assert impl.section_area(T, C, n=20000) == pytest.approx(a_exact, rel=1e-6)
    assert impl.section_area_substituted(T, C, n=4000) == pytest.approx(a_exact, rel=1e-7)
    assert impl.section_area_exact(T, 2.0 * C) == pytest.approx(4.0 * a_exact)        # ~ c^2
    assert impl.section_area_exact(2.0 * T, C) == pytest.approx(2.0 * a_exact)        # ~ t
    assert impl.section_area_exact(T, C, closed_te=True) < a_exact                    # closing the TE removes area


def test_centroid_of_symmetric_and_cambered_sections():
    xb, yb = impl.section_centroid(T, C, 0.0, 0.4, n=4000)
    assert yb == pytest.approx(0.0, abs=1e-12)                                       # symmetric about the chord
    assert xb == pytest.approx(impl.section_centroid_x_exact(T, C), rel=1e-4)
    assert 0.40 * C < xb < 0.45 * C                                                   # ahead of mid-chord: the nose is fat
    assert impl.section_centroid_x_exact(T, 3.0 * C) == pytest.approx(3.0 * impl.section_centroid_x_exact(T, C))
    xb2, yb2 = impl.section_centroid(T, C, 0.02, 0.4, n=4000)
    assert yb2 > 0.0 and yb2 < 0.02 * C                                              # lifted by the camber, below max camber
    assert xb2 == pytest.approx(xb, rel=1e-12)                                       # vertical thickness: x_bar unchanged
    xb3, yb3 = impl.section_centroid(T, C, 0.04, 0.4, n=4000)
    assert yb3 == pytest.approx(2.0 * yb2, rel=1e-12)                                # y_bar linear in camber


# ---------------------------------------------------------------- net change

def test_net_change_and_total_variation():
    v = lambda s: 2.0 * np.sin(np.pi * s / 60.0) - 0.8
    net = impl.net_change(v, 0.0, 120.0, n=4000)
    assert net == pytest.approx(-0.8 * 120.0, rel=1e-9)                              # the sine integrates to zero over a period
    tv = impl.total_variation(v, 0.0, 120.0, n=4000)
    assert tv > abs(net)
    assert impl.total_variation(lambda s: np.ones_like(s), 0.0, 5.0) == pytest.approx(5.0)
    assert impl.net_change(np.cos, 0.0, 2.0 * np.pi, n=4000) == pytest.approx(0.0, abs=1e-9)
