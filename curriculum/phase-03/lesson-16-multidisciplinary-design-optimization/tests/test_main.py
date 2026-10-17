"""Public tests for P03.16 Multidisciplinary Design Optimization. They assert the physics of the two
disciplines, the fixed point and its convergence rates, agreement of the three derivative routes with
each other and with finite differences (observed order 2), and agreement of the MDF and IDF optima,
with scipy.optimize (root, SLSQP) as the oracle. No OpenMDAO is imported anywhere.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-03/parameterized-cases/.
"""
import dataclasses

import numpy as np
import numpy.testing as npt
import pytest
from scipy.optimize import root

from conftest import load_impl

impl = load_impl(__file__)

P = impl.WingParams()
X0 = np.array([12.0, 16.8])


def _fd(fun, z, h=1e-6):
    z = np.asarray(z, dtype=float)
    f0 = np.atleast_1d(fun(z))
    J = np.zeros((f0.size, z.size))
    for j in range(z.size):
        e = np.zeros(z.size); e[j] = h * max(1.0, abs(z[j]))
        J[:, j] = (np.atleast_1d(fun(z + e)) - np.atleast_1d(fun(z - e))) / (2 * e[j])
    return J


# ---------------------------------------------------------------- the disciplines

def test_aero_mission_reproduces_the_formulas_and_their_limits():
    W_wing = 2500.0
    a = impl.aero_mission(X0, W_wing, P)
    b, S = X0
    L = P.W0 + W_wing
    assert a["L"] == pytest.approx(L)
    D_i = L ** 2 / (P.q * np.pi * P.e * b ** 2)
    assert a["D"] == pytest.approx(P.q * (P.f_rest + P.c_dw * S) + D_i)
    assert a["LD"] == pytest.approx(L / a["D"])
    assert a["W_fuel"] == pytest.approx(L * (np.exp(P.kappa * a["D"] / L) - 1.0))
    assert a["W_TO"] == pytest.approx(L + a["W_fuel"])
    assert a["M_root"] == pytest.approx(P.n_ult * a["W_TO"] * b / (3 * np.pi))
    # Limits: zero range burns no fuel; doubling the span quarters the induced drag.
    p0 = dataclasses.replace(P, R=0.0)
    assert impl.aero_mission(X0, W_wing, p0)["W_fuel"] == pytest.approx(0.0, abs=1e-9)
    a2 = impl.aero_mission(np.array([2 * b, S]), W_wing, P)
    assert (a["D"] - a2["D"]) == pytest.approx(0.75 * D_i)
    # More range or more empty weight always costs more fuel.
    assert impl.aero_mission(X0, W_wing, dataclasses.replace(P, R=2 * P.R))["W_fuel"] > a["W_fuel"]
    assert impl.aero_mission(X0, W_wing, dataclasses.replace(P, W0=1.2 * P.W0))["W_fuel"] > a["W_fuel"]


def test_aero_partials_match_finite_differences():
    for x, W_wing in ((X0, 2500.0), (np.array([16.0, 20.0]), 4000.0)):
        d = impl.aero_mission_partials(x, W_wing, P)
        for key in ("M_root", "W_fuel", "W_TO"):
            fd = _fd(lambda z: impl.aero_mission(z[:2], z[2], P)[key], np.array([x[0], x[1], W_wing]))[0]
            npt.assert_allclose(d[key], fd, rtol=1e-6, atol=1e-9)


def test_structures_reproduce_the_formulas_and_the_load_free_deflection():
    M = 1.5e5
    s = impl.structures(X0, M, P)
    b, S = X0
    h = P.tau * S / b
    assert s["A_cap"] == pytest.approx(M / (h * P.sigma_a))
    assert s["W_caps"] == pytest.approx(P.rho_m * impl.G0 * (2 / 3) * s["A_cap"] * b)
    assert s["W_wing"] == pytest.approx(s["W_caps"] + P.w_skin * S)
    assert s["defl_ratio"] == pytest.approx(P.sigma_a * b / (3 * P.E * h))
    assert impl.structures(X0, 2 * M, P)["defl_ratio"] == s["defl_ratio"]        # fully stressed: load-free
    assert impl.structures(X0, 2 * M, P)["W_caps"] == pytest.approx(2 * s["W_caps"])   # caps linear in M
    fd = _fd(lambda z: impl.structures(z[:2], z[2], P)["W_wing"], np.array([b, S, M]))[0]
    npt.assert_allclose(impl.structures_partials(X0, M, P), fd, rtol=1e-6)
    # A zero-moment wing still weighs its skins.
    assert impl.structures(X0, 0.0, P)["W_wing"] == pytest.approx(P.w_skin * S)


# ---------------------------------------------------------------- the coupled system

def test_coupled_residual_and_jacobians_are_consistent():
    y = np.array([1.5e5, 2.5e3])
    r = impl.coupled_residual(X0, y, P)
    assert r.shape == (2,)
    assert r[0] == pytest.approx(y[0] - impl.aero_mission(X0, y[1], P)["M_root"])
    assert r[1] == pytest.approx(y[1] - impl.structures(X0, y[0], P)["W_wing"])
    dRdy, dRdx = impl.coupled_jacobians(X0, y, P)
    npt.assert_allclose(dRdy, _fd(lambda yy: impl.coupled_residual(X0, yy, P), y), rtol=1e-6, atol=1e-8)
    npt.assert_allclose(dRdx, _fd(lambda xx: impl.coupled_residual(xx, y, P), X0), rtol=1e-6, atol=1e-8)
    npt.assert_allclose(np.diag(dRdy), 1.0)                       # R = y - G with G_i independent of y_i


def test_gauss_seidel_converges_linearly_at_the_coupling_strength():
    """Observed rate test: a contraction y <- G(y) converges linearly with ratio |k| where k is the
    product of the two off-diagonal coupling derivatives; the measured ratio of successive residuals
    must match it to 3 % once the iteration is in its asymptotic regime."""
    res = impl.mda_gauss_seidel(X0, P, tol=1e-11)
    assert res.converged and res.n_iter < 30
    npt.assert_allclose(impl.coupled_residual(X0, res.y, P), 0.0, atol=1e-4)     # 1e-9 of M_root ~ 1e5 N m
    r = np.array(res.residuals)
    assert np.all(np.diff(r) < 0)
    k = impl.coupling_strength(X0, res.y, P)
    assert 0 < k < 0.2
    ratios = r[3:7] / r[2:6]
    npt.assert_allclose(ratios, k, rtol=0.03)
    # The coupling strengthens with span (caps ~ b^3 / S) and weakens with area.
    y18 = impl.mda_newton(np.array([18.0, 16.8]), P).y
    assert impl.coupling_strength(np.array([18.0, 16.8]), y18, P) > 1.5 * k


def test_newton_mda_converges_quadratically_and_matches_gauss_seidel_and_scipy_root():
    nw = impl.mda_newton(X0, P, tol=1e-13)
    gs = impl.mda_gauss_seidel(X0, P, tol=1e-12)
    assert nw.converged and nw.n_iter <= 6
    npt.assert_allclose(nw.y, gs.y, rtol=1e-9)
    # Observed order test: from a poor start, three residuals above the rounding floor (1e-9 of
    # M_root ~ 1e5) give q = ln(r3/r2) / ln(r2/r1) near 2.
    r = np.array(impl.mda_newton(X0, P, y0=np.array([0.0, 0.0]), tol=1e-13).residuals)
    r = r[r > 1e-9]
    assert r.size >= 3
    q = np.log(r[-1] / r[-2]) / np.log(r[-2] / r[-3])
    assert q > 1.7, f"observed order {q:.2f} from residuals {r}"
    ref = root(lambda y: impl.coupled_residual(X0, y, P), np.array([1e5, 2e3]), tol=1e-14)
    assert ref.success
    npt.assert_allclose(nw.y, ref.x, rtol=1e-9)
    # From a poor start Newton still converges (the problem is mildly nonlinear in y).
    nw2 = impl.mda_newton(X0, P, y0=np.array([0.0, 0.0]))
    assert nw2.converged
    npt.assert_allclose(nw2.y, nw.y, rtol=1e-9)


# ---------------------------------------------------------------- derivatives

def test_direct_adjoint_and_finite_difference_total_derivatives_agree():
    F_d, dF_d = impl.total_derivatives_direct(X0, P)
    F_a, dF_a, psi = impl.total_derivatives_adjoint(X0, P)
    dF_fd = impl.total_derivatives_fd(X0, P, h_rel=1e-6)
    npt.assert_allclose(F_d, F_a)
    npt.assert_allclose(dF_d, dF_a, rtol=1e-10)                   # same linear algebra, transposed
    npt.assert_allclose(dF_fd, dF_a, rtol=1e-5)
    assert dF_a.shape == (2, 2) and psi.shape == (2, 2)
    # The adjoint variables solve (dR/dy)^T psi_i = (dF_i/dy)^T.
    y = impl.mda_newton(X0, P).y
    dRdy, _ = impl.coupled_jacobians(X0, y, P)
    _, dFdy = impl.functions_partials(X0, y, P)
    npt.assert_allclose(dRdy.T @ psi.T, dFdy.T, rtol=1e-10, atol=1e-14)
    # Signs: more span saves fuel at the baseline, more area costs fuel, more area helps the stall margin.
    assert dF_a[0, 0] < 0 and dF_a[0, 1] > 0 and dF_a[1, 1] > 0
    # The partial derivative with frozen coupling is not the total: the weight feedback matters.
    dFdx, _ = impl.functions_partials(X0, y, P)
    assert abs(dFdx[0, 0] - dF_a[0, 0]) > 0.05 * abs(dF_a[0, 0])


def test_functions_partials_match_finite_differences():
    y = np.array([1.5e5, 2.5e3])
    dFdx, dFdy = impl.functions_partials(X0, y, P)
    npt.assert_allclose(dFdx, _fd(lambda xx: impl.functions_of_interest(xx, y, P), X0), rtol=1e-6, atol=1e-10)
    npt.assert_allclose(dFdy, _fd(lambda yy: impl.functions_of_interest(X0, yy, P), y), rtol=1e-6, atol=1e-12)
    assert dFdy[0, 0] == 0.0 and dFdy[1, 0] == 0.0                # neither function sees M_root directly


def test_finite_difference_error_is_second_order_in_the_step():
    """Observed order test: central differences of the converged MDA against the adjoint derivative
    have error ~ h^2; the log-log slope over h = 1e-1, 1e-2, 1e-3 must be 2 within 0.1."""
    hs = np.array([1e-1, 1e-2, 1e-3])
    errs = impl.fd_convergence_study(X0, P, hs)
    assert np.all(np.diff(errs) < 0)
    slope = np.polyfit(np.log(hs), np.log(errs), 1)[0]
    assert slope == pytest.approx(2.0, abs=0.1)
    npt.assert_allclose(errs[1:] / errs[:-1], 0.01, rtol=0.05)


# ---------------------------------------------------------------- optimisation

def test_mdf_optimum_is_feasible_better_than_baseline_and_sits_on_the_deflection_limit():
    x, f, info = impl.mdf_optimize(P, X0)
    assert info["result"].success
    y = impl.mda_newton(x, P).y
    g = np.concatenate([[impl.functions_of_interest(x, y, P)[1]], impl.design_constraints_x_only(x, P)])
    assert np.all(g >= -1e-8)
    assert abs(g[1]) < 1e-6                                       # deflection (aspect-ratio) limit active
    base = impl.aero_mission(X0, impl.mda_newton(X0, P).y[1], P)["W_fuel"]
    assert f < 0.95 * base
    assert f == pytest.approx(impl.functions_of_interest(x, y, P)[0], rel=1e-10)
    assert x[0] > X0[0]                                           # the optimiser bought span ...
    assert x[0] ** 2 / x[1] == pytest.approx(P.defl_max * 3 * P.E * P.tau / P.sigma_a, rel=1e-6)   # ... up to the AR limit
    # Same optimum from a different start.
    x2, f2, _ = impl.mdf_optimize(P, np.array([10.0, 14.0]))
    npt.assert_allclose(x2, x, rtol=1e-5)
    assert f2 == pytest.approx(f, rel=1e-8)


def test_idf_and_mdf_agree_and_adjoint_gradients_save_analyses():
    x_m, f_m, info_m = impl.mdf_optimize(P, X0, use_adjoint=True)
    x_fd, f_fd, info_fd = impl.mdf_optimize(P, X0, use_adjoint=False)
    x_i, y_i, f_i, info_i = impl.idf_optimize(P, X0)
    npt.assert_allclose(x_fd, x_m, rtol=1e-5)
    npt.assert_allclose(x_i, x_m, rtol=1e-5)
    assert f_i == pytest.approx(f_m, rel=1e-7) and f_fd == pytest.approx(f_m, rel=1e-7)
    # IDF's consistency constraints hold at its optimum, so y_i is the MDA solution there.
    npt.assert_allclose(impl.coupled_residual(x_i, y_i, P), 0.0, atol=1e-6)
    npt.assert_allclose(y_i, impl.mda_newton(x_i, P).y, rtol=1e-8)
    assert info_fd["n_mda"] > info_m["n_mda"]                     # FD gradients cost extra MDA solves


def test_optimum_responds_to_parameters_in_the_physical_direction():
    x, f, _ = impl.mdf_optimize(P, X0)
    # A stiffer spar material relaxes the deflection limit: more span, less fuel.
    x_e, f_e, _ = impl.mdf_optimize(dataclasses.replace(P, E=1.5 * P.E), X0)
    assert x_e[0] > x[0] and f_e < f
    # A longer range costs more fuel.
    _, f_r, _ = impl.mdf_optimize(dataclasses.replace(P, R=1.2 * P.R), X0)
    assert f_r > f


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "wing_mdo.png").stat().st_size > 1000
