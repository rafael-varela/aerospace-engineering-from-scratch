"""Public tests for P03.15 Constrained Optimization. They assert KKT properties, agreement with the
analytic vertex and with scipy.optimize as the oracle, shadow prices and observed convergence orders.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-03/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.optimize import NonlinearConstraint, minimize

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(315)
PARAMS = impl.StrutParams()
X0 = np.array([1.0, 2.0])          # r = 50 mm, t = 2 mm: feasible and heavy


def _fd_jacobian(fun, x, h=1e-6):
    x = np.asarray(x, dtype=float)
    f0 = np.atleast_1d(fun(x))
    J = np.zeros((f0.size, x.size))
    for j in range(x.size):
        e = np.zeros(x.size); e[j] = h
        J[:, j] = (np.atleast_1d(fun(x + e)) - np.atleast_1d(fun(x - e))) / (2 * h)
    return J


# ---------------------------------------------------------------- the strut model

def test_strut_mass_and_gradient_are_consistent_and_physical():
    x = np.array([1.2, 0.8])
    r, t = impl.unscale(x, PARAMS)
    assert r == pytest.approx(0.06) and t == pytest.approx(0.8e-3)
    assert impl.strut_mass(x, PARAMS) == pytest.approx(2700 * 2.5 * 2 * np.pi * 0.06 * 0.8e-3)
    npt.assert_allclose(impl.strut_mass_grad(x, PARAMS), _fd_jacobian(lambda z: impl.strut_mass(z, PARAMS), x)[0], rtol=1e-6)
    # Mass is linear in each variable: doubling the thickness doubles the mass.
    assert impl.strut_mass(x * np.array([1, 2]), PARAMS) == pytest.approx(2 * impl.strut_mass(x, PARAMS))


def test_strut_constraints_match_the_physics_and_their_jacobian():
    x = np.array([1.1, 0.7])
    r, t = impl.unscale(x, PARAMS)
    g = impl.strut_constraints(x, PARAMS)
    assert g.shape == (5,)
    sigma = PARAMS.P / (2 * np.pi * r * t)
    assert g[0] == pytest.approx(1 - sigma / PARAMS.sigma_allow)
    P_euler = np.pi ** 2 * PARAMS.E * (np.pi * r ** 3 * t) / PARAMS.L ** 2      # pi^2 E I / L^2
    assert g[1] == pytest.approx(1 - PARAMS.P / P_euler)
    assert g[2] == pytest.approx(1 - sigma / (PARAMS.k_local * PARAMS.E * t / r))
    assert g[3] == pytest.approx(t / PARAMS.t_min - 1) and g[4] == pytest.approx(r / PARAMS.r_min - 1)
    npt.assert_allclose(impl.strut_jacobian(x, PARAMS), _fd_jacobian(lambda z: impl.strut_constraints(z, PARAMS), x), rtol=1e-5, atol=1e-8)
    # The Lagrangian Hessian by finite differences of the Lagrangian gradient.
    lam = np.array([0.1, 0.5, 0.3, 0.0, 0.2])
    grad_L = lambda z: impl.strut_mass_grad(z, PARAMS) - impl.strut_jacobian(z, PARAMS).T @ lam
    H = impl.strut_hessian_lagrangian(x, lam, PARAMS)
    npt.assert_allclose(H, H.T, atol=1e-9)
    npt.assert_allclose(H, _fd_jacobian(grad_L, x, 1e-5), rtol=1e-4, atol=1e-6)


def test_analytic_optimum_is_a_feasible_vertex_with_zero_active_constraints():
    x_star, active, mass = impl.analytic_optimum(PARAMS)
    g = impl.strut_constraints(x_star, PARAMS)
    assert np.all(g >= -1e-10)
    npt.assert_allclose(g[list(active)], 0.0, atol=1e-10)
    assert mass == pytest.approx(impl.strut_mass(x_star, PARAMS))
    assert set(active) == {1, 2}                               # Euler and local buckling govern at 40 kN
    # Hand check of the vertex: t^2 = P / (2 pi k E), r^3 t = P L^2 / (pi^3 E).
    r, t = impl.unscale(x_star, PARAMS)
    assert t == pytest.approx(np.sqrt(PARAMS.P / (2 * np.pi * PARAMS.k_local * PARAMS.E)), rel=1e-10)
    assert r ** 3 * t == pytest.approx(PARAMS.P * PARAMS.L ** 2 / (np.pi ** 3 * PARAMS.E), rel=1e-10)
    # It really is the minimum: random feasible points are all heavier.
    for _ in range(200):
        x = np.array([RNG.uniform(0.3, 2.5), RNG.uniform(0.5, 3.0)])
        if np.all(impl.strut_constraints(x, PARAMS) >= 0):
            assert impl.strut_mass(x, PARAMS) >= mass - 1e-12


def test_active_set_changes_with_the_load():
    chart = impl.design_chart(np.array([5e3, 40e3, 160e3]), PARAMS)
    masses = [m for _, m, _ in chart]
    assert masses == sorted(masses)                            # heavier load, heavier strut
    assert set(chart[0][2]) == {1, 3}                          # light load: Euler + minimum gauge
    assert set(chart[2][2]) == {0, 1}                          # heavy load: stress + Euler
    # When stress and Euler govern, the mass is proportional to the load.
    chart2 = impl.design_chart(np.array([160e3, 320e3]), PARAMS)
    assert chart2[1][1] == pytest.approx(2 * chart2[0][1], rel=1e-10)


# ---------------------------------------------------------------- KKT machinery

def test_kkt_residual_is_zero_at_the_solution_and_positive_elsewhere():
    prob = impl.textbook_problem()
    assert impl.kkt_residual(prob, np.array([0.5, 0.5]), np.array([1.0, 0.0])) < 1e-14
    assert impl.kkt_residual(prob, np.array([0.5, 0.5]), np.array([0.0, 0.0])) == pytest.approx(1.0)   # stationarity fails
    assert impl.kkt_residual(prob, np.array([0.2, 0.2]), np.array([1.0, 0.0])) > 0.5                   # infeasible
    assert impl.kkt_residual(prob, np.array([0.5, 0.5]), np.array([1.0, -0.3])) >= 0.3                     # dual infeasibility
    # Least-squares multipliers recover (1, 0) at the solution.
    lam = impl.multipliers_least_squares(prob, np.array([0.5, 0.5]), [0])
    npt.assert_allclose(lam, [1.0, 0.0], atol=1e-12)


def test_qp_solver_matches_hand_solutions_and_scipy():
    H, c = np.eye(2), np.zeros(2)
    # Project the origin onto x + y >= 1: p = (1/2, 1/2); stationarity p = lam (1, 1) gives lam = 1/2.
    p, lam = impl.solve_qp(H, c, np.array([[1.0, 1.0]]), np.array([1.0]))
    npt.assert_allclose(p, [0.5, 0.5], atol=1e-12)
    npt.assert_allclose(lam, [0.5], atol=1e-12)
    # An inactive constraint gets a zero multiplier and the unconstrained minimiser is returned.
    p, lam = impl.solve_qp(2 * np.eye(2), np.array([-2.0, -4.0]), np.array([[1.0, 0.0]]), np.array([-5.0]))
    npt.assert_allclose(p, [1.0, 2.0], atol=1e-12)
    assert lam[0] == 0.0
    # Random convex QPs against scipy's trust-constr.
    for _ in range(5):
        M = RNG.standard_normal((3, 3)); H = M @ M.T + 0.5 * np.eye(3)
        c = RNG.standard_normal(3); A = RNG.standard_normal((4, 3)); b = RNG.standard_normal(4) - 2.0
        p, lam = impl.solve_qp(H, c, A, b)
        assert np.all(A @ p - b >= -1e-8) and np.all(lam >= 0)
        npt.assert_allclose(lam * (A @ p - b), 0.0, atol=1e-8)          # complementarity
        npt.assert_allclose(H @ p + c - A.T @ lam, 0.0, atol=1e-8)      # stationarity
        ref = minimize(lambda z: 0.5 * z @ H @ z + c @ z, np.zeros(3), jac=lambda z: H @ z + c, method="trust-constr",
                       constraints=[NonlinearConstraint(lambda z: A @ z - b, 0, np.inf, jac=lambda z: A)],
                       options={"gtol": 1e-10, "xtol": 1e-12})
        npt.assert_allclose(p, ref.x, atol=1e-5)
    with pytest.raises(ValueError):
        impl.solve_qp(-np.eye(2), np.zeros(2), np.eye(2), np.zeros(2))        # not positive definite
    with pytest.raises(RuntimeError):
        impl.solve_qp(np.eye(1), np.zeros(1), np.array([[1.0], [-1.0]]), np.array([1.0, 1.0]))   # x >= 1 and x <= -1


def test_damped_bfgs_keeps_positive_definiteness_and_satisfies_the_secant_equation():
    B = np.eye(3)
    for _ in range(20):
        s = RNG.standard_normal(3); y = RNG.standard_normal(3)
        B_new = impl.damped_bfgs_update(B, s, y)
        assert np.min(np.linalg.eigvalsh(B_new)) > 0                       # stays SPD even for s.y < 0
        if s @ y >= 0.2 * s @ B @ s:
            npt.assert_allclose(B_new @ s, y, atol=1e-10)                  # undamped: secant equation holds
        B = B_new
    assert impl.merit_l1(impl.textbook_problem(), np.array([0.0, 0.0]), 3.0) == pytest.approx(3.0)   # f = 0, violation 1


# ---------------------------------------------------------------- SQP

@pytest.mark.parametrize("hessian", ["bfgs", "exact"])
def test_sqp_solves_the_textbook_problem_from_several_starts(hessian):
    prob = impl.textbook_problem()
    for x0 in ([3.0, -1.0], [0.0, 0.0], [-0.5, 2.0]):
        res = impl.sqp(prob, np.array(x0), hessian=hessian, tol=1e-10)
        assert res.converged
        npt.assert_allclose(res.x, [0.5, 0.5], atol=1e-8)
        npt.assert_allclose(res.lam, [1.0, 0.0], atol=1e-8)
        assert impl.kkt_residual(prob, res.x, res.lam) < 1e-10
        assert res.history[0].tolist() == x0 and len(res.history) == res.n_iter + 1


@pytest.mark.parametrize("hessian", ["bfgs", "exact"])
def test_sqp_sizes_the_strut_to_the_analytic_vertex(hessian):
    prob = impl.strut_problem(PARAMS)
    x_star, active, m_star = impl.analytic_optimum(PARAMS)
    res = impl.sqp(prob, X0, hessian=hessian, tol=1e-10)
    assert res.converged and res.n_iter <= 30
    npt.assert_allclose(res.x, x_star, rtol=1e-8)
    assert res.fun == pytest.approx(m_star, rel=1e-10)
    lam_ref = impl.multipliers_least_squares(prob, x_star, list(active))
    npt.assert_allclose(res.lam, lam_ref, atol=1e-7)
    assert np.all(res.lam >= 0) and res.lam[0] == 0.0 and res.lam[3] == 0.0     # stress and gauge inactive
    # From an infeasible start the l1 merit function still brings SQP home.
    res2 = impl.sqp(prob, np.array([0.5, 0.4]), hessian=hessian, tol=1e-10)
    assert res2.converged
    npt.assert_allclose(res2.x, x_star, rtol=1e-8)


def test_sqp_converges_superlinearly_once_close():
    """Observed order test: the last three nonzero errors of the SQP iteration give an order
    q = ln(e3/e2)/ln(e2/e1) near 2 (quadratic for the exact Hessian, superlinear for BFGS), far above
    the q = 1 of a penalty or augmented-Lagrangian iteration."""
    prob = impl.strut_problem(PARAMS)
    x_star, _, _ = impl.analytic_optimum(PARAMS)
    for hessian in ("bfgs", "exact"):
        res = impl.sqp(prob, X0, hessian=hessian, tol=1e-12)
        err = np.array([np.linalg.norm(h - x_star) for h in res.history])
        err = err[err > 1e-12]                                     # stop before rounding floors the error
        q = impl.observed_order(err[-3:])
        assert q > 1.5, f"{hessian}: observed order {q:.2f}, errors {err[-3:]}"
        assert err[-1] / err[-2] < 0.05                            # the last ratio is tiny: superlinear


def test_sqp_agrees_with_scipy_slsqp_and_trust_constr():
    prob = impl.strut_problem(PARAMS)
    res = impl.sqp(prob, X0, hessian="bfgs", tol=1e-10)
    cons = [{"type": "ineq", "fun": prob.cons, "jac": prob.jac}]
    ref = minimize(prob.f, X0, jac=prob.grad, constraints=cons, method="SLSQP", options={"ftol": 1e-14, "maxiter": 200})
    npt.assert_allclose(res.x, ref.x, rtol=1e-7)
    ref2 = minimize(prob.f, X0, jac=prob.grad, method="trust-constr",
                    constraints=[NonlinearConstraint(prob.cons, 0, np.inf, jac=prob.jac)],
                    options={"gtol": 1e-10, "xtol": 1e-12, "maxiter": 1000})
    npt.assert_allclose(res.x, ref2.x, rtol=1e-4)
    # And on Rosenbrock inside a disk (a nonconvex objective, one curved constraint).
    rb = impl.rosenbrock_disk_problem(1.5)
    res_rb = impl.sqp(rb, np.array([0.0, 0.0]), hessian="bfgs", tol=1e-9, max_iter=300)
    ref_rb = minimize(rb.f, [0.0, 0.0], jac=rb.grad, constraints=[{"type": "ineq", "fun": rb.cons, "jac": rb.jac}],
                      method="SLSQP", options={"ftol": 1e-14, "maxiter": 500})
    assert res_rb.converged
    npt.assert_allclose(res_rb.x, ref_rb.x, atol=1e-6)
    assert rb.cons(res_rb.x)[0] == pytest.approx(0.0, abs=1e-8)          # the disk constraint is active
    assert res_rb.lam[0] > 0


# ---------------------------------------------------------------- sensitivities, penalty, augmented Lagrangian

def test_multipliers_are_shadow_prices():
    pred, fd = impl.shadow_price_check(PARAMS)
    assert pred == pytest.approx(fd, rel=1e-4)
    assert pred > 0                                                        # more load, more mass
    # Doubling the load does not double the mass here (Euler + local: t ~ P^(1/2), r ~ P^(1/6), mass ~ P^(2/3)).
    chart = impl.design_chart(np.array([40e3, 80e3]), PARAMS)
    assert 1.5 < chart[1][1] / chart[0][1] < 2.0


def test_quadratic_penalty_error_is_first_order_in_one_over_mu():
    """Observed order test: |x_mu - x*| ~ C / mu, so the log-log slope against 1/mu is 1."""
    prob = impl.strut_problem(PARAMS)
    x_star, _, _ = impl.analytic_optimum(PARAMS)
    mus = np.array([1e2, 1e3, 1e4])
    errs = np.array([np.linalg.norm(impl.penalty_minimize(prob, mu, X0) - x_star) for mu in mus])
    assert np.all(np.diff(errs) < 0)
    order = impl.observed_order_from_resolutions(errs, 1.0 / mus)
    assert order == pytest.approx(1.0, abs=0.1)
    # The penalty minimiser is always slightly infeasible (it sits outside the feasible set).
    xm = impl.penalty_minimize(prob, 1e3, X0)
    assert np.min(prob.cons(xm)) < 0
    assert impl.observed_order(np.array([1e-1, 1e-2, 1e-4])) == pytest.approx(2.0)
    with pytest.raises(ValueError):
        impl.observed_order(np.array([1.0, 0.1]))


def test_augmented_lagrangian_reaches_the_sqp_solution_with_a_fixed_penalty():
    prob = impl.strut_problem(PARAMS)
    x_star, active, _ = impl.analytic_optimum(PARAMS)
    x, lam, k = impl.augmented_lagrangian(prob, X0, mu=10.0, tol=1e-8)
    npt.assert_allclose(x, x_star, rtol=1e-6)
    npt.assert_allclose(lam, impl.multipliers_least_squares(prob, x_star, list(active)), atol=1e-6)
    assert k <= 30
    assert impl.kkt_residual(prob, x, lam) < 1e-8


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "strut_sqp.png").stat().st_size > 1000
