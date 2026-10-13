"""Public tests for P01.09 Gradients, Directional Derivatives & Optimization. They assert calculus and physics, not magic numbers.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-01/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
import sympy as sp

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(11)


def _ac(**kw):
    return impl.Aircraft(**kw)


# ---------------------------------------------------------------- the drag function

def test_drag_has_parasite_and_induced_parts_that_scale_correctly():
    ac = _ac()
    V, AR = 50.0, 8.0
    D = impl.drag(V, AR, ac)
    q = 0.5 * ac.rho * V**2
    parasite = q * ac.S * ac.cd0
    induced = impl.induced_factor(AR, ac) * impl.weight(AR, ac) ** 2 / (q * ac.S)
    assert D == pytest.approx(parasite + induced, rel=1e-12)
    # No induced drag when the weight is zero; no parasite drag when C_D0 = 0.
    assert impl.drag(V, AR, _ac(W0=0.0, c_w=0.0)) == pytest.approx(parasite, rel=1e-12)
    assert impl.drag(V, AR, _ac(cd0=0.0)) == pytest.approx(induced, rel=1e-12)
    # Parasite part grows as V**2, induced part falls as V**-2.
    assert impl.drag(2 * V, AR, _ac(W0=0.0, c_w=0.0)) == pytest.approx(4 * parasite, rel=1e-12)
    assert impl.drag(2 * V, AR, _ac(cd0=0.0)) == pytest.approx(induced / 4, rel=1e-12)
    assert impl.induced_factor(2 * AR, ac) == pytest.approx(impl.induced_factor(AR, ac) / 2)
    grid = impl.drag(np.linspace(30, 60, 4)[:, None], np.linspace(5, 15, 3)[None, :], ac)
    assert grid.shape == (4, 3)
    for bad in [(0.0, AR), (-10.0, AR), (V, 0.0)]:
        with pytest.raises(ValueError):
            impl.drag(*bad, ac)


# ---------------------------------------------------------------- gradient, three ways

def test_analytic_gradient_matches_central_differences_at_second_order():
    ac = _ac()
    f = lambda x: impl.drag(x[0], x[1], ac)
    for x0 in (np.array([35.0, 6.0]), np.array([55.0, 12.0]), np.array([46.0, 10.0])):
        g = impl.drag_gradient(x0[0], x0[1], ac)
        errs = np.array([np.linalg.norm(impl.gradient_fd(f, x0, h) - g) for h in (0.4, 0.2, 0.1)])
        rates = np.log2(errs[:-1] / errs[1:])
        assert np.all(rates > 1.9), rates
        npt.assert_allclose(impl.gradient_fd(f, x0, np.array([1e-3, 1e-4])), g, rtol=1e-7)
    with pytest.raises(ValueError):
        impl.gradient_fd(f, np.array([35.0, 6.0]), 0.0)


def test_gradient_fd_on_a_known_function():
    f = lambda x: x[0] ** 2 * x[1] + np.sin(x[2])
    x = np.array([1.5, -2.0, 0.3])
    npt.assert_allclose(impl.gradient_fd(f, x, 1e-4), [2 * 1.5 * -2.0, 1.5**2, np.cos(0.3)], rtol=1e-7)
    H = impl.hessian_fd(f, x, 1e-3)
    H_exact = np.array([[2 * -2.0, 2 * 1.5, 0.0], [2 * 1.5, 0.0, 0.0], [0.0, 0.0, -np.sin(0.3)]])
    npt.assert_allclose(H, H_exact, atol=1e-5)
    npt.assert_allclose(H, H.T)                                                  # symmetric (Clairaut)


def test_gradient_vanishes_at_the_minimum_drag_speed_for_fixed_ar():
    ac = _ac()
    for AR in (6.0, 10.0, 14.0):
        W = impl.weight(AR, ac)
        V_md = impl.minimum_drag_speed(W, AR, ac)
        assert impl.drag_gradient(V_md, AR, ac)[0] == pytest.approx(0.0, abs=1e-9)
        # At V_md the parasite and induced drags are equal and D = 2 W sqrt(K C_D0).
        q = 0.5 * ac.rho * V_md**2
        assert q * ac.S * ac.cd0 == pytest.approx(impl.induced_factor(AR, ac) * W**2 / (q * ac.S), rel=1e-12)
        assert impl.drag(V_md, AR, ac) == pytest.approx(2 * W * np.sqrt(impl.induced_factor(AR, ac) * ac.cd0), rel=1e-12)
        # Below V_md the slope is negative, above it positive.
        assert impl.drag_gradient(0.8 * V_md, AR, ac)[0] < 0.0 < impl.drag_gradient(1.2 * V_md, AR, ac)[0]
    assert impl.minimum_drag_speed(4 * 15000.0, 10.0, ac) == pytest.approx(2 * impl.minimum_drag_speed(15000.0, 10.0, ac), rel=1e-12)
    with pytest.raises(ValueError):
        impl.minimum_drag_speed(-1.0, 10.0, ac)


# ---------------------------------------------------------------- directional derivatives

def test_directional_derivative_properties():
    g = np.array([-48.0, -115.0])
    assert impl.directional_derivative(g, [1.0, 0.0]) == pytest.approx(g[0])
    assert impl.directional_derivative(g, [0.0, 3.0]) == pytest.approx(g[1])          # direction is normalised
    # Steepest ascent along the gradient equals |grad|; steepest descent is its negative; level direction gives 0.
    assert impl.directional_derivative(g, g) == pytest.approx(np.linalg.norm(g))
    assert impl.directional_derivative(g, -g) == pytest.approx(-np.linalg.norm(g))
    assert impl.directional_derivative(g, [g[1], -g[0]]) == pytest.approx(0.0, abs=1e-12)
    # No direction beats the gradient: random unit vectors give smaller slopes.
    for d in RNG.normal(size=(50, 2)):
        assert impl.directional_derivative(g, d) <= np.linalg.norm(g) + 1e-12
    with pytest.raises(ValueError):
        impl.directional_derivative(g, [0.0, 0.0])


def test_directional_derivative_matches_a_finite_difference_along_the_line():
    ac = _ac()
    x0 = np.array([40.0, 8.0])
    g = impl.drag_gradient(*x0, ac)
    for d in ([1.0, 0.0], [0.0, 1.0], [1.0, 0.1], [-2.0, 0.3]):
        u = np.asarray(d) / np.linalg.norm(d)
        h = 1e-4
        fd = (impl.drag(*(x0 + h * u), ac) - impl.drag(*(x0 - h * u), ac)) / (2 * h)
        assert impl.directional_derivative(g, d) == pytest.approx(fd, rel=1e-6)


# ---------------------------------------------------------------- the optimum: analytic, searched, classified

def test_analytic_optimum_is_a_critical_point_and_a_minimum():
    ac = _ac()
    V, AR, D = impl.optimum_analytic(ac)
    npt.assert_allclose(impl.drag_gradient(V, AR, ac), 0.0, atol=1e-8)
    assert D == pytest.approx(impl.drag(V, AR, ac), rel=1e-12)
    H = impl.hessian_fd(lambda x: impl.drag(x[0], x[1], ac), np.array([V, AR]), np.array([1e-2, 1e-3]))
    assert impl.classify_critical_point(H) == "minimum"
    assert np.all(np.linalg.eigvalsh(H) > 0.0)
    # Random nearby points all have more drag.
    for dv, da in RNG.normal(size=(50, 2)):
        assert impl.drag(V + 3 * dv, AR + 0.5 * da, ac) >= D
    # AR* = sqrt(W0 / (3 c_w)) scales as expected; V* is V_md at that AR.
    assert AR == pytest.approx(np.sqrt(ac.W0 / (3 * ac.c_w)), rel=1e-12)
    assert V == pytest.approx(impl.minimum_drag_speed(impl.weight(AR, ac), AR, ac), rel=1e-12)
    _, AR4, _ = impl.optimum_analytic(_ac(c_w=4 * ac.c_w))
    assert AR4 == pytest.approx(AR / 2, rel=1e-12)


def test_classify_critical_point_on_known_hessians():
    assert impl.classify_critical_point(np.diag([2.0, 3.0])) == "minimum"
    assert impl.classify_critical_point(np.diag([-2.0, -3.0])) == "maximum"
    assert impl.classify_critical_point(np.array([[1.0, 0.0], [0.0, -1.0]])) == "saddle"
    assert impl.classify_critical_point(np.array([[1.0, 2.0], [2.0, 1.0]])) == "saddle"     # eigenvalues 3, -1
    assert impl.classify_critical_point(np.diag([1.0, 0.0])) == "degenerate"


def test_scaled_gradient_descent_converges_to_the_analytic_optimum():
    ac = _ac()
    V, AR, D = impl.optimum_analytic(ac)
    f, g = impl.scaled_drag(ac)
    z_opt, path = impl.gradient_descent(f, g, np.array([0.7, 0.6]), step=0.2, tol=1e-9, max_iter=20000)
    assert 50 * z_opt[0] == pytest.approx(V, rel=1e-5)
    assert 10 * z_opt[1] == pytest.approx(AR, rel=1e-5)
    assert path.shape[1] == 2 and path.shape[0] >= 2
    npt.assert_allclose(path[0], [0.7, 0.6])
    fs = np.array([f(z) for z in path])
    assert np.all(np.diff(fs) < 0.0)                                             # monotone decrease
    with pytest.raises(ValueError):
        impl.gradient_descent(f, g, np.array([0.7, 0.6]), step=0.0)


def test_gradient_descent_on_a_quadratic_bowl():
    A = np.array([[3.0, 1.0], [1.0, 2.0]])
    b = np.array([1.0, -2.0])
    f = lambda x: 0.5 * x @ A @ x - b @ x
    g = lambda x: A @ x - b
    x_star = np.linalg.solve(A, b)
    x, path = impl.gradient_descent(f, g, np.array([5.0, 5.0]), step=0.3, tol=1e-10)
    npt.assert_allclose(x, x_star, atol=1e-8)
    xn, pn = impl.newton_minimize(f, g, lambda x: A, np.array([5.0, 5.0]))
    npt.assert_allclose(xn, x_star, atol=1e-12)
    assert len(pn) <= 3                                                           # Newton is exact on a quadratic in one step


def test_newton_converges_quadratically_on_the_scaled_drag():
    ac = _ac()
    V, AR, _ = impl.optimum_analytic(ac)
    f, g = impl.scaled_drag(ac)
    h = lambda z: impl.hessian_fd(f, z, 1e-4)
    z, path = impl.newton_minimize(f, g, h, np.array([0.7, 0.6]), tol=1e-12)
    assert 50 * z[0] == pytest.approx(V, rel=1e-8) and 10 * z[1] == pytest.approx(AR, rel=1e-8)
    assert len(path) <= 12
    errs = np.linalg.norm(path - np.array([V / 50, AR / 10]), axis=1)
    errs = errs[errs > 1e-13]
    # Error roughly squares each step once in the basin: log-error differences grow.
    ratios = errs[1:] / errs[:-1]
    assert ratios[-1] < 1e-2


# ---------------------------------------------------------------- Lagrange multipliers

def test_lagrange_solution_reproduces_minimum_drag_speed():
    ac = _ac()
    V_s, AR_s, _ = impl.optimum_analytic(ac)
    W = impl.weight(AR_s, ac)
    V, cl, lam = impl.lagrange_min_drag(W, ac)
    K = impl.induced_factor(AR_s, ac)
    assert cl == pytest.approx(np.sqrt(ac.cd0 / K), rel=1e-12)
    assert V == pytest.approx(impl.minimum_drag_speed(W, AR_s, ac), rel=1e-12)
    # The constraint holds and the stationarity conditions grad D = lambda grad L hold in (V, C_L).
    q = 0.5 * ac.rho * V**2
    assert q * ac.S * cl == pytest.approx(W, rel=1e-12)
    dD_dV, dD_dcl = ac.rho * V * ac.S * (ac.cd0 + K * cl**2), 2 * q * ac.S * K * cl
    dL_dV, dL_dcl = ac.rho * V * ac.S * cl, q * ac.S
    assert dD_dV == pytest.approx(lam * dL_dV, rel=1e-12) and dD_dcl == pytest.approx(lam * dL_dcl, rel=1e-12)
    # lambda = dD/dW: the marginal drag of carrying one more newton equals 1/(L/D) at the optimum.
    dW = 10.0

    def drag_at_optimum(Wt):
        Vt, clt, _ = impl.lagrange_min_drag(Wt, ac)
        return 0.5 * ac.rho * Vt**2 * ac.S * (ac.cd0 + K * clt**2)

    assert (drag_at_optimum(W + dW) - drag_at_optimum(W - dW)) / (2 * dW) == pytest.approx(lam, rel=1e-6)
    assert lam == pytest.approx(drag_at_optimum(W) / W, rel=1e-12)
    assert drag_at_optimum(W) == pytest.approx(impl.drag(V, AR_s, ac), rel=1e-12)


def test_lagrange_sympy_matches_closed_forms():
    sol = impl.lagrange_sympy()
    K, cd0, W, rho, S = sp.symbols("K C_D0 W rho S", positive=True)
    assert sp.simplify(sol["C_L"] - sp.sqrt(cd0 / K)) == 0
    assert sp.simplify(sol["lambda"] - 2 * sp.sqrt(K * cd0)) == 0
    assert sp.simplify(sol["V"] ** 2 - 2 * W / (rho * S) * sp.sqrt(K / cd0)) == 0


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "drag_landscape.png").stat().st_size > 1000
