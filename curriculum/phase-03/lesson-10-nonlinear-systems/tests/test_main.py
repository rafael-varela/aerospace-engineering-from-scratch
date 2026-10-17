"""Public tests for P03.10 Nonlinear Systems. They assert Jacobian identities, convergence to the
``scipy.optimize`` solutions, observed convergence orders (Newton 2, Picard 1, forward differences 1),
the failure and rescue of undamped Newton, and the continuation sweep of trim conditions.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-03/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.optimize import root, newton_krylov

from conftest import load_impl

impl = load_impl(__file__)

AC = impl.Aircraft()
H, V = 5000.0, 200.0
RHO = impl.isa_density(H)


def trim_f(x):
    return impl.trim_residual(x, AC, V, RHO)


def trim_j(x):
    return impl.trim_jacobian(x, AC, V, RHO)


# ---------------------------------------------------------------- the model

def test_atmosphere_and_dynamic_pressure():
    assert impl.isa_density(0.0) == pytest.approx(1.225, rel=1e-3)
    assert impl.isa_density(5000.0) == pytest.approx(0.7364, rel=1e-3)
    assert impl.dynamic_pressure(1.225, 100.0) == pytest.approx(6125.0)


def test_trim_residual_is_zero_at_a_hand_built_state():
    # choose alpha and delta_e, then set T from the drag equation: F2 vanishes exactly, F1 and F3 are what they are
    a, de = 0.02, 0.03
    qS = impl.dynamic_pressure(RHO, V) * AC.S
    CL = AC.CL0 + AC.CLa * a + AC.CLde * de
    CD = AC.CD0 + AC.k * CL ** 2
    T = qS * CD / np.cos(a)
    F = impl.trim_residual(np.array([a, de, T / AC.W]), AC, V, RHO)
    assert F[1] == pytest.approx(0.0, abs=1e-14)
    assert F[0] == pytest.approx((qS * CL + T * np.sin(a) - AC.W) / AC.W)
    assert F[2] == pytest.approx((qS * AC.c * (AC.Cm0 + AC.Cma * a + AC.Cmde * de) + T * AC.zT) / (AC.W * AC.c))


def test_analytic_jacobian_matches_finite_differences_and_complex_step():
    x = np.array([0.03, -0.02, 0.08])
    J = trim_j(x)
    assert J.shape == (3, 3)
    npt.assert_allclose(impl.fd_jacobian(trim_f, x), J, rtol=1e-6, atol=1e-7)
    npt.assert_allclose(impl.complex_step_jacobian(trim_f, x), J, rtol=1e-13, atol=1e-14)


def test_forward_difference_jacobian_is_first_order_in_h():
    f = lambda x: np.array([x[0] ** 3 + x[1], np.exp(x[0] * x[1])])
    x = np.array([0.7, -0.4])
    J_exact = np.array([[3 * x[0] ** 2, 1.0], [x[1] * np.exp(x[0] * x[1]), x[0] * np.exp(x[0] * x[1])]])
    hs = np.array([1e-2, 1e-3, 1e-4, 1e-5])
    errs = np.array([np.max(np.abs(impl.fd_jacobian(f, x, h) - J_exact)) for h in hs])
    slope = np.polyfit(np.log(hs), np.log(errs), 1)[0]
    assert 0.9 < slope < 1.1
    assert np.max(np.abs(impl.fd_jacobian(f, x) - J_exact)) < 1e-6          # default step ~ sqrt(eps)
    assert np.max(np.abs(impl.complex_step_jacobian(f, x) - J_exact)) < 1e-14


# ---------------------------------------------------------------- Newton on scalar problems

def test_newton_is_quadratically_convergent_on_sqrt2():
    f = lambda x: np.array([x[0] ** 2 - 2.0])
    j = lambda x: np.array([[2.0 * x[0]]])
    res = impl.newton(f, np.array([3.0]), j, tol=1e-14)
    assert res.converged
    assert res.x[0] == pytest.approx(np.sqrt(2.0), rel=1e-14)
    errs = impl.iterate_errors(res.xs, np.array([np.sqrt(2.0)]))
    assert 1.8 < impl.convergence_order(errs) < 2.2
    # Newton for sqrt is the Babylonian iteration: x_{k+1} = (x_k + 2/x_k)/2
    assert res.xs[1][0] == pytest.approx((3.0 + 2.0 / 3.0) / 2)


def test_newton_fd_jacobian_counts_evaluations():
    f = lambda x: np.array([x[0] ** 2 - 2.0])
    res = impl.newton(f, np.array([3.0]), None, tol=1e-12)
    assert res.converged
    assert res.nfev == 1 + res.iterations * 2                                # one residual + one FD column per step


def test_undamped_newton_diverges_on_arctan_and_damping_rescues_it():
    for x0 in (1.5, 3.0):
        ru = impl.newton(impl.arctan_residual, np.array([x0]), impl.arctan_jacobian, tol=1e-12, max_iter=20)
        assert not ru.converged
        assert abs(ru.x[0]) > 10.0
        rd = impl.damped_newton(impl.arctan_residual, np.array([x0]), impl.arctan_jacobian, tol=1e-12)
        assert rd.converged
        assert abs(rd.x[0]) < 1e-12
        assert rd.step_lengths[0] < 1.0                                      # the first step was shortened
        assert all(l == 1.0 for l in rd.step_lengths[1:])                    # then full Newton steps
    ru = impl.newton(impl.arctan_residual, np.array([1.0]), impl.arctan_jacobian, tol=1e-12)
    assert ru.converged                                                      # inside the basin, no damping needed


def test_convergence_order_estimator():
    quadratic = 0.5 ** (2.0 ** np.arange(6))                                 # e_{k+1} = e_k^2
    assert impl.convergence_order(quadratic) == pytest.approx(2.0, abs=1e-8)
    linear = 0.3 * 0.5 ** np.arange(12)
    assert impl.convergence_order(linear) == pytest.approx(1.0, abs=1e-8)
    with pytest.raises(ValueError):
        impl.convergence_order(np.array([0.1, 1e-20]))


# ---------------------------------------------------------------- trim

def test_newton_trim_matches_scipy_and_converges_quadratically():
    res = impl.solve_trim(AC, V, H, "newton", crude=True)
    assert res.converged and res.iterations <= 6
    sol = root(trim_f, np.zeros(3), jac=trim_j, method="hybr", tol=1e-13)
    assert sol.success
    npt.assert_allclose(res.x, sol.x, rtol=1e-10, atol=1e-12)
    assert np.linalg.norm(trim_f(res.x)) < 1e-12
    errs = impl.iterate_errors(res.xs, res.x)
    assert 1.6 < impl.convergence_order(errs) < 2.7
    # physics: lift ~ weight, thrust ~ drag, elevator a few degrees
    qS = impl.dynamic_pressure(RHO, V) * AC.S
    CL = AC.CL0 + AC.CLa * res.x[0] + AC.CLde * res.x[1]
    assert CL == pytest.approx(AC.W / qS, rel=5e-3)
    assert res.x[2] * AC.W == pytest.approx(qS * (AC.CD0 + AC.k * CL ** 2), rel=1e-2)
    assert 0.0 < np.degrees(res.x[1]) < 5.0


def test_fd_newton_broyden_and_damped_reach_the_same_trim():
    ref = impl.solve_trim(AC, V, H, "newton", crude=True)
    fd = impl.solve_trim(AC, V, H, "newton-fd", crude=True)
    br = impl.solve_trim(AC, V, H, "broyden", crude=True)
    da = impl.solve_trim(AC, V, H, "damped", crude=True)
    for r in (fd, br, da):
        assert r.converged
        npt.assert_allclose(r.x, ref.x, rtol=1e-9, atol=1e-11)
    assert fd.nfev == 1 + 4 * fd.iterations
    assert br.nfev < fd.nfev                                                  # one evaluation per step after the first Jacobian
    assert br.iterations > ref.iterations                                    # but more steps than Newton
    e = impl.iterate_errors(br.xs, ref.x)
    ratios = e[1:] / e[:-1]
    assert ratios[-1] < ratios[1]                                            # superlinear: the contraction improves
    assert all(l == 1.0 for l in da.step_lengths)                             # the crude start is inside Newton's basin
    with pytest.raises(ValueError):
        impl.solve_trim(AC, V, H, "simplex")


def test_initial_guesses():
    x0 = impl.trim_initial_guess(AC, V, RHO)
    qS = impl.dynamic_pressure(RHO, V) * AC.S
    assert x0[1] == 0.0
    assert AC.CL0 + AC.CLa * x0[0] == pytest.approx(AC.W / qS)
    npt.assert_array_equal(impl.trim_initial_guess(AC, V, RHO, crude=True), np.zeros(3))
    res = impl.solve_trim(AC, V, H, "newton")                                 # the physical guess converges faster
    assert res.converged and res.iterations <= impl.solve_trim(AC, V, H, "newton", crude=True).iterations


def test_trim_sweep_continuation_matches_cold_starts_and_is_physical():
    Vs = np.arange(120.0, 251.0, 10.0)
    Xw, itw = impl.trim_sweep(AC, Vs, H, warm_start=True)
    Xc, itc = impl.trim_sweep(AC, Vs, H, warm_start=False)
    assert Xw.shape == (Vs.size, 3)
    npt.assert_allclose(Xw, Xc, rtol=1e-8, atol=1e-10)
    assert itw.sum() <= itc.sum()
    assert np.all(np.diff(Xw[:, 0]) < 0)                                     # alpha falls with speed
    assert np.all(np.diff(Xw[:, 1]) > 0)                                     # elevator moves trailing-edge-down with speed
    thrust = Xw[:, 2] * AC.W
    V_md = np.sqrt(2 * AC.W / (RHO * AC.S)) * (AC.k / AC.CD0) ** 0.25        # minimum-drag speed, CD0 = k CL^2
    assert abs(Vs[np.argmin(thrust)] - V_md) <= 10.0
    assert impl.stall_speed(AC, H, 1.6) == pytest.approx(np.sqrt(2 * AC.W / (RHO * AC.S * 1.6)))


# ---------------------------------------------------------------- the hot plate

def test_conduction_matrix_reduces_to_the_kronecker_laplacian_and_is_spd():
    n, h, k = 5, 0.04, 0.34
    A = impl.conduction_matrix(np.full((n + 2, n + 2), k), n, h)
    T = np.diag(2 * np.ones(n)) + np.diag(-np.ones(n - 1), 1) + np.diag(-np.ones(n - 1), -1)
    ref = (k / h ** 2) * (np.kron(np.eye(n), T) + np.kron(T, np.eye(n)))
    npt.assert_allclose(A, ref, atol=1e-12)
    rng = np.random.default_rng(310)
    kf = rng.uniform(0.2, 1.0, (n + 2, n + 2))
    A = impl.conduction_matrix(kf, n, h)
    npt.assert_allclose(A, A.T, atol=1e-15)
    assert np.all(np.linalg.eigvalsh(A) > 0)
    b = impl.conduction_boundary_rhs(kf, n, h, 300.0)
    npt.assert_allclose(np.linalg.solve(A, b), 300.0, rtol=1e-12)              # no source: wall temperature everywhere


def test_plate_residual_vanishes_without_source_and_linear_case_needs_one_step():
    pl = impl.HotPlate(n=6, q=0.0)
    npt.assert_allclose(impl.plate_residual(np.full(36, pl.T_wall), pl), 0.0, atol=1e-9)
    lin = impl.HotPlate(n=6, beta=0.0)
    rp = impl.picard_plate(lin)
    assert rp.converged and rp.iterations == 1
    rn = impl.newton_plate(lin)
    assert rn.converged and rn.iterations == 1


def test_picard_is_first_order_and_newton_second_order_on_the_plate():
    pl = impl.HotPlate(n=8)
    rp = impl.picard_plate(pl, tol=1e-9)
    rn = impl.newton_plate(pl, tol=1e-9)
    assert rp.converged and rn.converged
    assert rn.iterations <= 6 and rp.iterations > rn.iterations
    npt.assert_allclose(rp.x, rn.x, rtol=1e-9)
    ep = impl.iterate_errors(rp.xs, rn.x) / 100.0
    en = impl.iterate_errors(rn.xs, rn.x) / 100.0
    assert 0.85 < impl.convergence_order(ep) < 1.15
    assert impl.convergence_order(en) > 1.6
    # the conductance falls with temperature, so the plate runs hotter than the linear-k0 model
    lin = impl.newton_plate(impl.HotPlate(n=8, beta=0.0))
    assert rn.x.max() > lin.x.max()
    assert rn.x.max() - pl.T_wall == pytest.approx(47.0, abs=3.0)


def test_plate_matches_scipy_newton_krylov():
    pl = impl.HotPlate(n=8)
    rn = impl.newton_plate(pl, tol=1e-10)
    Tk = newton_krylov(lambda T: impl.plate_residual(T, pl), np.full(64, pl.T_wall), f_tol=1e-6)
    npt.assert_allclose(rn.x, Tk, rtol=1e-8)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "nonlinear_systems.png").exists()
