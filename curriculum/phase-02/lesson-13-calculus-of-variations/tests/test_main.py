"""Public tests for P02.13 Calculus of Variations. They assert properties: Euler-Lagrange equations checked
symbolically, the cycloid meets its end conditions and its closed-form time, every competing path is slower,
the first variation vanishes on the extremal (second-order time increase), Snell's invariant is constant,
and the direct (Ritz) minimisation approaches the cycloid from above.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-02/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
import sympy as sp

from conftest import load_impl

impl = load_impl(__file__)

C = impl.chute()
G = C["g"]
RNG = np.random.default_rng(213)


# ---------------------------------------------------------------- symbolic

def test_euler_lagrange_reproduces_classical_results():
    x = sp.symbols("x")
    y = sp.Function("y")(x)
    yp = sp.Derivative(y, x)
    # Shortest path: the Euler-Lagrange expression vanishes iff y'' = 0.
    el = impl.euler_lagrange(sp.sqrt(1 + yp ** 2), y, x)
    assert sp.simplify(el * (1 + yp ** 2) ** sp.Rational(3, 2) + sp.Derivative(y, (x, 2))) == 0
    # Harmonic oscillator action: F = y'^2 / 2 - w^2 y^2 / 2 gives y'' + w^2 y = 0.
    w = sp.symbols("omega", positive=True)
    el = impl.euler_lagrange(yp ** 2 / 2 - w ** 2 * y ** 2 / 2, y, x)
    assert sp.simplify(el + sp.Derivative(y, (x, 2)) + w ** 2 * y) == 0
    # A quadratic integrand with explicit x: F = y'^2 + 2 x y gives y'' = x.
    el = impl.euler_lagrange(yp ** 2 + 2 * x * y, y, x)
    assert sp.simplify(el - (2 * x - 2 * sp.Derivative(y, (x, 2)))) == 0
    # An extremal substituted back satisfies the equation: y = x^3 / 6 for y'' = x.
    assert sp.simplify(el.subs(y, x ** 3 / 6).doit()) == 0


def test_beltrami_identity_gives_the_cycloid_equation_and_rejects_explicit_x():
    x = sp.symbols("x")
    y = sp.Function("y")(x)
    yp = sp.Derivative(y, x)
    g = sp.symbols("g", positive=True)
    belt = impl.beltrami_first_integral(sp.sqrt((1 + yp ** 2) / (2 * g * y)), y, x)
    # 1 / (2 g belt^2) must be y (1 + y'^2).
    assert sp.simplify(1 / (2 * g * belt ** 2) - y * (1 + yp ** 2)) == 0
    # The cycloid satisfies y (1 + y'^2) = 2k: check with the parametric form.
    th, k = sp.symbols("theta k", positive=True)
    X, Y = k * (th - sp.sin(th)), k * (1 - sp.cos(th))
    yprime = sp.diff(Y, th) / sp.diff(X, th)
    assert sp.simplify(Y * (1 + yprime ** 2) - 2 * k) == 0
    # Shortest path: F - y' F_y' = 1 / sqrt(1 + y'^2) constant means y' constant.
    belt = impl.beltrami_first_integral(sp.sqrt(1 + yp ** 2), y, x)
    assert sp.simplify(belt - 1 / sp.sqrt(1 + yp ** 2)) == 0
    with pytest.raises(ValueError):
        impl.beltrami_first_integral(yp ** 2 + x * y, y, x)


def test_hamiltons_principle_gives_the_pendulum_and_spring_mass_equations():
    t = sp.symbols("t")
    g, m, l, k = sp.symbols("g m l k", positive=True)
    th = sp.Function("theta")(t)
    eom = impl.lagrangian_equation(sp.Rational(1, 2) * m * l ** 2 * sp.Derivative(th, t) ** 2, -m * g * l * sp.cos(th), th, t)
    assert sp.simplify(eom / (m * l ** 2) - (sp.Derivative(th, (t, 2)) + g / l * sp.sin(th))) == 0
    q = sp.Function("q")(t)
    eom = impl.lagrangian_equation(sp.Rational(1, 2) * m * sp.Derivative(q, t) ** 2, sp.Rational(1, 2) * k * q ** 2, q, t)
    assert sp.simplify(eom - (m * sp.Derivative(q, (t, 2)) + k * q)) == 0
    # Free particle: no force.
    assert sp.simplify(impl.lagrangian_equation(sp.Rational(1, 2) * m * sp.Derivative(q, t) ** 2, 0, q, t) - m * sp.Derivative(q, (t, 2))) == 0


# ---------------------------------------------------------------- path times

def test_path_time_is_exact_for_the_ramp_and_second_order_for_the_cycloid():
    x_f, y_f = C["x_f"], C["y_f"]
    x = np.linspace(0.0, x_f, 5)
    assert impl.travel_time(x, y_f * x / x_f, G) == pytest.approx(impl.straight_line_time(x_f, y_f, G), rel=1e-12)   # exact, any n
    k, th = impl.cycloid_parameter(x_f, y_f)
    T = impl.brachistochrone_time(k, th, G)
    errs = []
    for n in (21, 41, 81, 161):
        xb, yb, _, _ = impl.brachistochrone(x_f, y_f, n)
        errs.append(impl.travel_time(xb, yb, G) - T)
    errs = np.array(errs)
    assert np.all(errs > 0.0)                                                     # chords are slower than the curve
    npt.assert_allclose(errs[:-1] / errs[1:], 4.0, rtol=0.1)                      # O(h^2)
    # Constant speed: time = length / v.
    assert impl.path_time([0, 3, 3], [0, 4, 4], lambda x, y: 2.0 + 0 * x) == pytest.approx(2.5)
    with pytest.raises(ValueError):
        impl.path_time([0, 1], [0, 0], lambda x, y: 0 * x)                        # zero speed at both ends
    with pytest.raises(ValueError):
        impl.path_time([0, 1], [0, 1], lambda x, y: -1.0 + 0 * x)
    with pytest.raises(ValueError):
        impl.straight_line_time(1.0, 0.0, G)


def test_straight_line_time_scales_as_expected():
    T = impl.straight_line_time(4.0, 2.0, G)
    assert impl.straight_line_time(8.0, 4.0, G) == pytest.approx(np.sqrt(2) * T)   # lengths x2 -> time x sqrt 2
    assert impl.straight_line_time(4.0, 2.0, 4 * G) == pytest.approx(T / 2)         # g x4 -> time / 2
    assert impl.straight_line_time(0.0, 2.0, G) == pytest.approx(np.sqrt(2 * 2.0 / G))   # free fall


# ---------------------------------------------------------------- the cycloid

def test_gravity_speed_and_cycloid_scaling():
    v = impl.gravity_speed(G)
    npt.assert_allclose(v(np.zeros(3), np.array([0.0, 1.0, 4.0])), np.sqrt(2 * G) * np.array([0.0, 1.0, 2.0]))
    assert v(0.0, -1.0) == 0.0                                                    # above the start: no speed, not NaN
    # Geometric similarity: scaling both coordinates scales k and keeps theta_f; time scales as sqrt(size).
    k1, th1 = impl.cycloid_parameter(4.0, 2.0)
    k2, th2 = impl.cycloid_parameter(8.0, 4.0)
    assert k2 == pytest.approx(2 * k1) and th2 == pytest.approx(th1)
    assert impl.brachistochrone_time(k2, th2, G) == pytest.approx(np.sqrt(2) * impl.brachistochrone_time(k1, th1, G))
    # Nearly vertical drop: theta_f -> 0 and the time -> free-fall time sqrt(2 y_f / g).
    k, th = impl.cycloid_parameter(1e-3, 2.0)
    assert th < 0.1
    assert impl.brachistochrone_time(k, th, G) == pytest.approx(np.sqrt(2 * 2.0 / G), rel=1e-3)


def test_cycloid_meets_the_end_points_and_the_closed_form_time():
    for x_f, y_f in ((4.0, 2.0), (1.0, 1.0), (np.pi, 2.0), (4.0, 0.5), (1.0, 4.0)):
        xb, yb, k, th = impl.brachistochrone(x_f, y_f, 301)
        assert xb[0] == 0.0 and yb[0] == 0.0
        assert xb[-1] == pytest.approx(x_f, abs=1e-10) and yb[-1] == pytest.approx(y_f, abs=1e-10)
        assert 0.0 < th < 2 * np.pi and k > 0.0
        assert impl.travel_time(xb, yb, G) == pytest.approx(impl.brachistochrone_time(k, th, G), rel=1e-4)
        assert impl.brachistochrone_time(k, th, G) < impl.straight_line_time(x_f, y_f, G)
        # y (1 + y'^2) = 2k along the curve (Beltrami), away from the vertical start.
        yp = np.gradient(yb, xb)
        npt.assert_allclose((yb * (1 + yp ** 2))[40:-1], 2 * k, rtol=1e-2)      # finite-difference y' is crude near the vertical start
    # The end point at the bottom of the cycloid: x_f / y_f = pi / 2 gives theta_f = pi exactly.
    k, th = impl.cycloid_parameter(np.pi, 2.0)
    assert th == pytest.approx(np.pi, abs=1e-9) and k == pytest.approx(1.0, abs=1e-9)
    with pytest.raises(ValueError):
        impl.cycloid_parameter(-1.0, 1.0)


def test_cycloid_beats_every_competing_path():
    x_f, y_f = C["x_f"], C["y_f"]
    xb, yb, k, th = impl.brachistochrone(x_f, y_f, 801)
    T = impl.brachistochrone_time(k, th, G)
    x = np.linspace(0.0, x_f, 801)
    for p in (0.5, 1.0, 1.5, 2.0, 3.0):
        assert impl.travel_time(x, y_f * (x / x_f) ** p, G) > T                     # power-law paths
    for _ in range(20):                                                           # random admissible perturbations of the ramp
        a = RNG.normal(scale=0.3, size=3)
        y = y_f * x / x_f + sum(ak * np.sin((j + 1) * np.pi * x / x_f) for j, ak in enumerate(a))
        if np.all(y[1:] > 0.0):
            assert impl.travel_time(x, y, G) > T


def test_first_variation_vanishes_on_the_cycloid_but_not_on_the_ramp():
    x_f, y_f = C["x_f"], C["y_f"]
    xb, yb, _, _ = impl.brachistochrone(x_f, y_f, 801)
    eps = np.array([0.1, 0.05, 0.025, 0.0125])
    dT = np.array([impl.perturbed_time_increase(xb, yb, e, G) for e in eps])
    assert np.all(dT > 0.0)
    slopes = np.log(dT[:-1] / dT[1:]) / np.log(2.0)
    npt.assert_allclose(slopes, 2.0, atol=0.1)                                   # second order: delta T ~ eps^2
    x = np.linspace(0.0, x_f, 801)
    dTl = np.array([abs(impl.perturbed_time_increase(x, y_f * x / x_f, e, G)) for e in eps])
    slopes_l = np.log(dTl[:-1] / dTl[1:]) / np.log(2.0)
    npt.assert_allclose(slopes_l, 1.0, atol=0.15)                                # first order on a non-extremal


def test_snell_invariant_is_constant_on_the_cycloid_and_not_on_the_ramp():
    x_f, y_f = C["x_f"], C["y_f"]
    xb, yb, k, _ = impl.brachistochrone(x_f, y_f, 801)
    sn = impl.snell_invariant(xb, yb, impl.gravity_speed(G))[5:]
    npt.assert_allclose(sn, 1.0 / np.sqrt(4 * G * k), rtol=2e-3)
    x = np.linspace(0.0, x_f, 801)
    snl = impl.snell_invariant(x, y_f * x / x_f, impl.gravity_speed(G))[5:]
    assert np.ptp(snl) > 0.5 * np.max(snl)
    # Straight segment in a uniform medium: sin(theta)/v is the direction cosine over the speed.
    npt.assert_allclose(impl.snell_invariant([0, 1, 2], [0, 1, 2], lambda x, y: 3.0 + 0 * x), np.sqrt(0.5) / 3.0)


# ---------------------------------------------------------------- Ritz

def test_ritz_minimum_approaches_the_cycloid_time_from_above():
    x_f, y_f = C["x_f"], C["y_f"]
    k, th = impl.cycloid_parameter(x_f, y_f)
    T = impl.brachistochrone_time(k, th, G)
    T_line = impl.straight_line_time(x_f, y_f, G)
    prev = T_line
    for n_modes in (1, 2, 4):
        x, y, T_r, a = impl.ritz_minimum_time(x_f, y_f, n_modes, G)
        assert T < T_r < prev                                                     # above the true minimum, improving with modes
        assert y[0] == pytest.approx(0.0) and y[-1] == pytest.approx(y_f)
        assert a.size == n_modes and a[0] > 0.0                                   # the path sags below the ramp
        prev = T_r
    assert (T_r - T) / T < 0.02
    # The Ritz path with zero coefficients is the ramp.
    xr, yr = impl.ritz_path(np.zeros(3), x_f, y_f)
    npt.assert_allclose(yr, y_f * xr / x_f)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "brachistochrone.png").stat().st_size > 1000
