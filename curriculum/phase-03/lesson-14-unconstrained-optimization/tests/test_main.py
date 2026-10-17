"""Tests for P03.14 Unconstrained Optimization (BFGS from scratch)."""
import numpy as np
import numpy.testing as npt
import pytest
from conftest import load_impl

impl = load_impl(__file__)

X_STAR = np.array([1.0, 1.0])
X0 = np.array([-1.2, 1.0])


def _spd_quadratic():
    Q = np.array([[10.0, 2.0, 0.0], [2.0, 4.0, 1.0], [0.0, 1.0, 1.0]])
    b = np.array([1.0, -2.0, 0.5])
    return Q, b


# ---------------------------------------------------------------- test functions and derivatives

def test_rosenbrock_values_and_gradient():
    assert impl.rosenbrock(X_STAR) == 0.0
    npt.assert_allclose(impl.rosenbrock_grad(X_STAR), 0.0)
    assert impl.rosenbrock(X0) == pytest.approx(24.2)
    npt.assert_allclose(impl.rosenbrock_grad(X0), [-215.6, -88.0])
    H = impl.rosenbrock_hess(X_STAR)
    npt.assert_allclose(H, H.T)
    assert np.all(np.linalg.eigvalsh(H) > 0)


def test_gradient_fd_observed_order_two():
    """Central differences: halving h divides the error by four until rounding takes over."""
    x = np.array([0.3, -0.7])
    exact = impl.rosenbrock_grad(x)
    hs = np.array([1e-2, 5e-3, 2.5e-3, 1.25e-3])
    errs = np.array([np.max(np.abs(impl.gradient_fd(impl.rosenbrock, x, h) - exact)) for h in hs])
    orders = np.log(errs[:-1] / errs[1:]) / np.log(2.0)
    npt.assert_allclose(orders, 2.0, atol=0.05)
    assert np.max(np.abs(impl.gradient_fd(impl.rosenbrock, x, 1e-6) - exact)) < 1e-6
    # too small a step: rounding error grows like eps / h
    assert np.max(np.abs(impl.gradient_fd(impl.rosenbrock, x, 1e-10) - exact)) > 1e-6


def test_hessian_fd_matches_analytic():
    x = np.array([0.3, -0.7])
    H = impl.hessian_fd(impl.rosenbrock, x, 1e-4)
    npt.assert_allclose(H, impl.rosenbrock_hess(x), rtol=1e-5, atol=1e-4)
    npt.assert_allclose(H, H.T)


# ---------------------------------------------------------------- line searches

def test_backtracking_satisfies_armijo_and_rejects_ascent():
    x = X0
    g = impl.rosenbrock_grad(x)
    a = impl.line_search_backtracking(impl.rosenbrock, x, -g, g)
    assert 0.0 < a <= 1.0
    assert impl.rosenbrock(x - a * g) <= impl.rosenbrock(x) + 1e-4 * a * (g @ -g)
    with pytest.raises(ValueError):
        impl.line_search_backtracking(impl.rosenbrock, x, g, g)


def test_wolfe_conditions_hold():
    c1, c2 = 1e-4, 0.9
    x = X0
    g = impl.rosenbrock_grad(x)
    p = -g / np.linalg.norm(g)
    a = impl.line_search_wolfe(impl.rosenbrock, impl.rosenbrock_grad, x, p, c1=c1, c2=c2)
    phi0, dphi0 = impl.rosenbrock(x), g @ p
    assert impl.rosenbrock(x + a * p) <= phi0 + c1 * a * dphi0
    assert abs(impl.rosenbrock_grad(x + a * p) @ p) <= c2 * abs(dphi0)
    # on a quadratic with alpha0 = exact step the search accepts alpha0 at once
    Q, b = _spd_quadratic()
    f, grad = impl.quadratic(Q, b)
    x = np.zeros(3)
    g = grad(x)
    a_exact = (g @ g) / (g @ Q @ g)
    assert impl.line_search_wolfe(f, grad, x, -g, alpha0=a_exact) == pytest.approx(a_exact)
    with pytest.raises(ValueError):
        impl.line_search_wolfe(impl.rosenbrock, impl.rosenbrock_grad, X0, impl.rosenbrock_grad(X0))


# ---------------------------------------------------------------- BFGS update algebra

def test_bfgs_update_secant_and_positive_definite():
    rng = np.random.default_rng(3)
    H = np.eye(4)
    for _ in range(5):
        s = rng.standard_normal(4)
        y = s + 0.3 * rng.standard_normal(4)
        if y @ s <= 0:
            y = -y
        H = impl.bfgs_update(H, s, y)
        npt.assert_allclose(H @ y, s, atol=1e-12)          # secant equation
        npt.assert_allclose(H, H.T, atol=1e-12)
        assert np.all(np.linalg.eigvalsh(H) > 0)
    with pytest.raises(ValueError):
        impl.bfgs_update(np.eye(2), np.array([1.0, 0.0]), np.array([-1.0, 0.0]))


def test_bfgs_quadratic_terminates_in_n_steps_with_exact_inverse():
    Q, b = _spd_quadratic()
    f, grad = impl.quadratic(Q, b)
    res = impl.bfgs(f, grad, np.zeros(3), tol=1e-12, exact_quadratic=Q)
    assert res["converged"]
    assert res["iterations"] <= 3
    npt.assert_allclose(res["x"], np.linalg.solve(Q, b), atol=1e-10)
    npt.assert_allclose(res["H"], np.linalg.inv(Q), atol=1e-10)


# ---------------------------------------------------------------- convergence behaviour

def test_steepest_descent_linear_rate_within_bound():
    Q, b = _spd_quadratic()
    f, grad = impl.quadratic(Q, b)
    x_star = np.linalg.solve(Q, b)
    res = impl.steepest_descent(f, grad, np.zeros(3), tol=1e-10, exact_quadratic=Q)
    assert res["converged"]
    fs = np.array([f(x) for x in res["path"]]) - f(x_star)
    fs = fs[fs > 1e-12]
    ratios = fs[1:] / fs[:-1]
    bound = impl.linear_rate_bound(Q)
    assert bound == pytest.approx(((np.linalg.cond(Q) - 1) / (np.linalg.cond(Q) + 1)) ** 2)
    assert np.all(ratios <= bound + 1e-9)
    assert np.max(ratios[5:]) > 0.3                   # genuinely linear, not superlinear
    assert res["iterations"] > 50


def test_newton_converges_quadratically_on_rosenbrock():
    res = impl.newton_method(impl.rosenbrock, impl.rosenbrock_grad, impl.rosenbrock_hess, X0, tol=1e-12)
    assert res["converged"]
    npt.assert_allclose(res["x"], X_STAR, atol=1e-10)
    assert res["iterations"] < 40
    e = impl.error_sequence(res["path"], X_STAR)
    # once full Newton steps are taken (error below 2e-3) each step roughly squares the error
    tail = e[e > 1e-13]
    k = np.where(tail < 2e-3)[0][0]
    assert tail[k + 1] < 20.0 * tail[k] ** 2
    assert tail[k + 2] < 20.0 * tail[k + 1] ** 2
    assert tail[k + 2] < 1e-3 * tail[k + 1]


def test_bfgs_rosenbrock_superlinear_and_matches_scipy():
    from scipy.optimize import minimize

    res = impl.bfgs(impl.rosenbrock, impl.rosenbrock_grad, X0, tol=1e-10)
    assert res["converged"]
    npt.assert_allclose(res["x"], X_STAR, atol=1e-8)
    assert res["iterations"] < 80
    e = impl.error_sequence(res["path"], X_STAR)
    e = e[e > 1e-13]
    ratios = e[1:] / e[:-1]
    assert np.mean(ratios[-4:]) < 0.2                 # superlinear: ratios tend to zero
    assert np.mean(ratios[:5]) > 0.5                  # ... but only once close
    sp = minimize(impl.rosenbrock, X0, jac=impl.rosenbrock_grad, method="BFGS", options={"gtol": 1e-10})
    npt.assert_allclose(res["x"], sp.x, atol=1e-6)
    npt.assert_allclose(res["H"], np.linalg.inv(impl.rosenbrock_hess(X_STAR)), atol=1e-3)


def test_convergence_order_estimator():
    lin = 0.5 ** np.arange(12)                        # linear: order 1
    npt.assert_allclose(impl.convergence_order(lin), 1.0, atol=1e-12)
    quad = np.array([0.5 ** (2 ** k) for k in range(5)])     # quadratic: order 2
    npt.assert_allclose(impl.convergence_order(quad), 2.0, atol=1e-12)
    assert impl.convergence_order(np.array([1.0, 1e-20])).size == 0


def test_bfgs_with_fd_gradient_and_custom_h0():
    f = impl.rosenbrock
    g = lambda x: impl.gradient_fd(f, x, 1e-6)        # noqa: E731
    res = impl.bfgs(f, g, X0, tol=1e-6)
    assert res["converged"]
    npt.assert_allclose(res["x"], X_STAR, atol=1e-5)
    res2 = impl.bfgs(impl.rosenbrock, impl.rosenbrock_grad, X0, tol=1e-10, H0=0.01 * np.eye(2))
    assert res2["converged"]
    npt.assert_allclose(res2["x"], X_STAR, atol=1e-8)


# ---------------------------------------------------------------- cruise problem

def test_isa_reference_points():
    T, p, rho, a = impl.isa(0.0)
    assert (T, p) == pytest.approx((288.15, 101325.0))
    assert rho == pytest.approx(1.225, rel=1e-3)
    assert a == pytest.approx(340.3, rel=1e-3)
    T, p, rho, a = impl.isa(11000.0)
    assert T == pytest.approx(216.65)
    assert p == pytest.approx(22632.0, rel=1e-3)
    T2, p2, _, _ = impl.isa(20000.0)
    assert T2 == pytest.approx(216.65)
    assert p2 == pytest.approx(5474.9, rel=2e-3)


def test_fuel_per_km_physics():
    ac = impl.cruise_aircraft()
    x = np.array([0.78, 1.0])
    bd = impl.cruise_breakdown(x, ac)
    T, p, rho, a = impl.isa(10000.0)
    V = 0.78 * a
    q = 0.5 * rho * V ** 2
    CL = ac["W"] / (q * ac["S"])
    assert bd["CL"] == pytest.approx(CL)
    assert bd["V_mps"] == pytest.approx(V)
    assert bd["drag_N"] == pytest.approx(ac["W"] / bd["L_over_D"])
    assert impl.fuel_per_km(x, ac) == pytest.approx(ac["cT0"] * np.sqrt(T / 288.15) * bd["drag_N"] / V * 1000.0)
    assert impl.fuel_per_km(x, ac) == pytest.approx(bd["fuel_per_km"])
    # below drag divergence there is no wave drag: fuel/km is independent of k_wave
    ac2 = dict(ac, k_wave=0.0)
    slow = np.array([0.5, 0.8])
    assert impl.fuel_per_km(slow, ac) == pytest.approx(impl.fuel_per_km(slow, ac2))
    with pytest.raises(ValueError):
        impl.fuel_per_km(np.array([-0.5, 1.0]), ac)
    with pytest.raises(ValueError):
        impl.fuel_per_km(np.array([0.5, -1.0]), ac)


def test_optimal_cruise_matches_scipy_and_is_a_minimum():
    from scipy.optimize import minimize

    ac = impl.cruise_aircraft()
    res = impl.optimal_cruise(ac)
    assert res["converged"]
    f = lambda x: impl.fuel_per_km(x, ac)             # noqa: E731
    sp = minimize(f, np.array([0.70, 0.80]), method="BFGS", options={"gtol": 1e-9})
    assert abs(res["x"][0] - sp.x[0]) < 1e-4          # Mach
    assert abs(res["x"][1] - sp.x[1]) * 1e4 < 10.0    # altitude within 10 m
    assert res["f"] == pytest.approx(sp.fun, rel=1e-6)
    H = impl.hessian_fd(f, res["x"], 1e-3)
    assert np.all(np.linalg.eigvalsh(H) > 0)
    # the optimum beats a grid of nearby points
    for dM in (-0.02, 0.0, 0.02):
        for dh in (-0.05, 0.0, 0.05):
            assert f(res["x"] + [dM, dh]) >= res["f"] - 1e-9
    assert 0.70 < res["x"][0] < 0.85
    assert 9000.0 < res["x"][1] * 1e4 < 13500.0


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "bfgs_cruise.png").stat().st_size > 10_000
