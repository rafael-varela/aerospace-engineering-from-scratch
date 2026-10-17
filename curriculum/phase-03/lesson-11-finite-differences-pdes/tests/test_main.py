"""Public tests for P03.11 Finite Differences for PDEs. They assert properties, with SciPy as the oracle.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-03/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import solve_ivp
from scipy.linalg import expm, solve_banded

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(311)
L, ALPHA = 0.01, 2.8e-6


def _grid(n):
    return np.linspace(0.0, L, n + 1), L / n


# ---------------------------------------------------------------- exact solutions

def test_exact_mode_and_step_series_satisfy_the_pde_and_the_boundary_conditions():
    x, dx = _grid(400)
    t = 2.0
    # Finite-difference residual of the PDE for both exact solutions, with a central time difference.
    for f in (lambda tt: impl.heat_exact_mode(x, tt, L, ALPHA, 2),
              lambda tt: impl.heat_exact_step(x, tt, L, ALPHA, 300.0, 700.0, 300.0)):
        h = 1e-3
        u_t = (f(t + h) - f(t - h)) / (2 * h)
        u_xx = impl.second_difference(f(t), dx)
        scale = np.max(np.abs(u_t[1:-1])) + 1e-30
        assert np.max(np.abs(u_t[1:-1] - ALPHA * u_xx)) / scale < 1e-3
    T = impl.heat_exact_step(x, t, L, ALPHA, 300.0, 700.0, 300.0)
    assert T[0] == pytest.approx(700.0) and T[-1] == pytest.approx(300.0)
    assert np.all(np.diff(T) <= 1e-9)                                     # monotone from hot to cold
    # Steady state is the linear profile; initial state is T0 away from the ends.
    npt.assert_allclose(impl.heat_exact_step(x, 1e4, L, ALPHA, 300.0, 700.0, 300.0), 700.0 - 400.0 * x / L, atol=1e-8)
    T0 = impl.heat_exact_step(x[20:-20], 0.0, L, ALPHA, 300.0, 700.0, 300.0, n_terms=20000)
    npt.assert_allclose(T0, 300.0, atol=0.5)
    # Mode decay rate.
    u1 = impl.heat_exact_mode(x, 1.0, L, ALPHA, 1)
    u2 = impl.heat_exact_mode(x, 2.0, L, ALPHA, 1)
    k = np.pi / L
    assert np.max(u2) / np.max(u1) == pytest.approx(np.exp(-ALPHA * k * k), rel=1e-10)


def test_cold_face_flux_matches_a_finite_difference_of_the_series_and_its_steady_limit():
    k_cond = 7.0
    q = impl.cold_face_flux_exact(5.0, L, ALPHA, k_cond, 300.0, 700.0, 300.0)
    h = 1e-7
    xs = np.array([L - h, L])
    T = impl.heat_exact_step(xs, 5.0, L, ALPHA, 300.0, 700.0, 300.0, n_terms=400)
    assert q == pytest.approx(-k_cond * (T[1] - T[0]) / h, rel=1e-4)
    assert q > 0.0
    assert impl.cold_face_flux_exact(1e4, L, ALPHA, k_cond, 300.0, 700.0, 300.0) == pytest.approx(k_cond * 400.0 / L, rel=1e-8)
    assert impl.cold_face_flux_exact(1e-3, L, ALPHA, k_cond, 300.0, 700.0, 300.0) == pytest.approx(0.0, abs=1e-6)


def test_thermal_diffusivity_and_the_case_values():
    assert impl.thermal_diffusivity(7.0, 4430.0, 560.0) == pytest.approx(7.0 / (4430.0 * 560.0))
    c = impl.titanium_bulkhead()
    a = impl.thermal_diffusivity(c["k"], c["rho"], c["cp"])
    assert 1e-6 < a < 1e-5                                                 # titanium: a poor conductor among metals
    with pytest.raises(ValueError):
        impl.thermal_diffusivity(-1.0, 1.0, 1.0)


# ---------------------------------------------------------------- stencils and the tridiagonal solver

def test_second_difference_is_second_order_accurate():
    f = lambda x: np.sin(3 * x) + x ** 3                                    # noqa: E731
    f2 = lambda x: -9 * np.sin(3 * x) + 6 * x                               # noqa: E731
    errs = []
    for n in (20, 40, 80, 160):
        x = np.linspace(0.0, 1.0, n + 1)
        errs.append(np.max(np.abs(impl.second_difference(f(x), x[1] - x[0]) - f2(x[1:-1]))))
    orders = impl.observed_order(np.array(errs))
    npt.assert_allclose(orders, 2.0, atol=0.05)
    # Exact for quadratics (the truncation error is proportional to the fourth derivative).
    x = np.linspace(0.0, 2.0, 11)
    npt.assert_allclose(impl.second_difference(3 * x ** 2 - x + 1, 0.2), 6.0, atol=1e-10)
    with pytest.raises(ValueError):
        impl.second_difference(np.ones(2), 0.1)


def test_thomas_algorithm_matches_scipy_solve_banded():
    for n in (1, 2, 5, 50):
        diag = 4.0 + RNG.uniform(0, 1, n)
        lower = RNG.uniform(-1, 1, max(n - 1, 0))
        upper = RNG.uniform(-1, 1, max(n - 1, 0))
        rhs = RNG.standard_normal(n)
        x = impl.solve_tridiagonal(lower, diag, upper, rhs)
        A = np.diag(diag) + np.diag(lower, -1) + np.diag(upper, 1)
        npt.assert_allclose(A @ x, rhs, atol=1e-12)
        if n > 1:
            ab = np.zeros((3, n))
            ab[0, 1:] = upper
            ab[1] = diag
            ab[2, :-1] = lower
            npt.assert_allclose(x, solve_banded((1, 1), ab, rhs), rtol=1e-12, atol=1e-12)
    with pytest.raises(ValueError):
        impl.solve_tridiagonal(np.ones(3), np.ones(3), np.ones(2), np.ones(3))
    with pytest.raises(ValueError):
        impl.solve_tridiagonal(np.ones(1), np.array([0.0, 1.0]), np.ones(1), np.ones(2))


# ---------------------------------------------------------------- one step of each scheme

def test_steps_hold_the_dirichlet_values_and_agree_to_first_order_in_dt():
    x, dx = _grid(50)
    u = 300.0 + 400.0 * np.exp(-((x - 0.3 * L) / (0.1 * L)) ** 2)
    u[0], u[-1] = 500.0, 300.0
    r = 0.01
    for step in (impl.step_ftcs, impl.step_btcs, impl.step_crank_nicolson):
        v = step(u, r)
        assert v[0] == 500.0 and v[-1] == 300.0 and v.shape == u.shape
        # All three agree with u + r D2 u to O(r^2).
        explicit = impl.step_ftcs(u, r)
        assert np.max(np.abs(v - explicit)) < 5 * r ** 2 * np.max(np.abs(np.diff(u, 2)))
    # r = 0 is the identity.
    for step in (impl.step_ftcs, impl.step_btcs, impl.step_crank_nicolson):
        npt.assert_allclose(step(u, 0.0), u, atol=1e-12)
        with pytest.raises(ValueError):
            step(u, -0.1)
    # Implicit steps satisfy their defining linear equations exactly.
    v = impl.step_btcs(u, 3.0)
    npt.assert_allclose(v[1:-1] - 3.0 * (v[:-2] - 2 * v[1:-1] + v[2:]), u[1:-1], atol=1e-9)
    w = impl.step_crank_nicolson(u, 3.0)
    lhs = w[1:-1] - 1.5 * (w[:-2] - 2 * w[1:-1] + w[2:])
    rhs = u[1:-1] + 1.5 * (u[:-2] - 2 * u[1:-1] + u[2:])
    npt.assert_allclose(lhs, rhs, atol=1e-9)


def test_crank_nicolson_is_the_average_of_ftcs_and_btcs_operators():
    # On the interior unknowns: (I - rA/2) G_cn = (I + rA/2), so G_cn differs from (G_ftcs + G_btcs)/2 at O(r^2) only
    # for small r, but the exact identity is that CN is the trapezoidal rule. Check it via the update matrices.
    n, r = 9, 0.3
    A = -2 * np.eye(n) + np.eye(n, k=1) + np.eye(n, k=-1)
    G = impl.update_matrix("cn", n, r)
    npt.assert_allclose((np.eye(n) - 0.5 * r * A) @ G, np.eye(n) + 0.5 * r * A, atol=1e-12)
    npt.assert_allclose(impl.update_matrix("ftcs", n, r), np.eye(n) + r * A, atol=1e-12)
    npt.assert_allclose(impl.update_matrix("btcs", n, r) @ (np.eye(n) - r * A), np.eye(n), atol=1e-12)
    with pytest.raises(ValueError):
        impl.update_matrix("leapfrog", n, r)


# ---------------------------------------------------------------- stability

def test_von_neumann_factors_and_matrix_eigenvalues_agree_and_give_the_half_rule():
    th = np.linspace(0.0, np.pi, 501)
    for r in (0.1, 0.5, 0.75, 4.0):
        g_f = impl.amplification_factor("ftcs", r, th)
        g_b = impl.amplification_factor("btcs", r, th)
        g_c = impl.amplification_factor("cn", r, th)
        assert g_f[0] == pytest.approx(1.0) and g_b[0] == pytest.approx(1.0) and g_c[0] == pytest.approx(1.0)
        assert np.all(np.abs(g_b) <= 1.0 + 1e-12) and np.all(np.abs(g_c) <= 1.0 + 1e-12)
        assert (np.max(np.abs(g_f)) <= 1.0 + 1e-12) == (r <= 0.5)
        # Matrix stability: the eigenvalues of G on n interior nodes are g(j pi / (n+1)).
        n = 30
        thj = np.arange(1, n + 1) * np.pi / (n + 1)
        for scheme in ("ftcs", "btcs", "cn"):
            eig = np.sort(np.linalg.eigvals(impl.update_matrix(scheme, n, r)).real)
            npt.assert_allclose(eig, np.sort(impl.amplification_factor(scheme, r, thj)), atol=1e-9)
            assert impl.spectral_radius(scheme, n, r) == pytest.approx(np.max(np.abs(impl.amplification_factor(scheme, r, thj))), rel=1e-9)
    assert impl.max_stable_dt("ftcs", 1e-4, ALPHA) == pytest.approx(1e-8 / (2 * ALPHA))
    assert impl.max_stable_dt("btcs", 1e-4, ALPHA) == np.inf and impl.max_stable_dt("cn", 1e-4, ALPHA) == np.inf
    # Crank-Nicolson's grid-scale mode tends to -1 as r grows: bounded but not damped.
    assert impl.amplification_factor("cn", 100.0, np.pi) == pytest.approx(-1.0, abs=0.02)
    assert abs(impl.amplification_factor("btcs", 100.0, np.pi)) < 0.01
    with pytest.raises(ValueError):
        impl.amplification_factor("dufort", 0.1, th)


def test_ftcs_blows_up_past_the_limit_and_march_refuses_it():
    x, dx = _grid(40)
    u0 = impl.heat_exact_mode(x, 0.0, L, ALPHA, 1) + 0.01 * np.sin(np.arange(x.size) * np.pi)  # seed the grid mode
    r_bad = 0.6
    u = u0.copy()
    for _ in range(200):
        u = impl.step_ftcs(u, r_bad)
    assert np.max(np.abs(u)) > 1e6                                         # (4 r - 1)^200 growth of the grid mode
    with pytest.raises(ValueError):
        impl.march(u0, dx, r_bad * dx ** 2 / ALPHA, ALPHA, 10, "ftcs")
    with pytest.raises(ValueError):
        impl.march(u0, dx, 1.0, ALPHA, 10, "upwind")
    # The implicit schemes take the same step without complaint and stay bounded.
    for scheme in ("btcs", "cn"):
        v = impl.march(u0, dx, 50 * dx ** 2 / ALPHA, ALPHA, 200, scheme)
        assert np.max(np.abs(v)) <= np.max(np.abs(u0)) + 1e-9


def test_implicit_schemes_obey_the_maximum_principle_differently_on_a_step():
    # A step initial condition at large r: BTCS stays monotone, Crank-Nicolson rings (grid-scale oscillation).
    n = 100
    x, dx = _grid(n)
    u0 = np.full(n + 1, 300.0)
    u0[0] = 700.0
    r = 100.0
    dt = r * dx ** 2 / ALPHA
    ub = impl.march(u0, dx, dt, ALPHA, 2, "btcs")
    uc = impl.march(u0, dx, dt, ALPHA, 2, "cn")
    assert np.all(np.diff(ub) <= 1e-9)                                     # monotone
    assert np.min(ub) >= 300.0 - 1e-9 and np.max(ub) <= 700.0 + 1e-9
    assert np.max(np.abs(np.diff(uc, 2))) > 20.0                           # ringing of tens of kelvin
    assert np.any(np.diff(uc) > 1e-6)                                      # not monotone


# ---------------------------------------------------------------- accuracy and convergence

def test_march_matches_the_exact_mode_the_semi_discrete_solution_and_scipy():
    n = 100
    x, dx = _grid(n)
    t_final = L ** 2 / (ALPHA * np.pi ** 2)
    u0 = impl.heat_exact_mode(x, 0.0, L, ALPHA, 1)
    exact = impl.heat_exact_mode(x, t_final, L, ALPHA, 1)
    semi = impl.semi_discrete_exact(x, t_final, L, ALPHA, 1)
    # The semi-discrete solution equals expm of the discrete operator applied to u0 (scipy oracle).
    A = ALPHA / dx ** 2 * (-2 * np.eye(n - 1) + np.eye(n - 1, k=1) + np.eye(n - 1, k=-1))
    npt.assert_allclose(expm(A * t_final) @ u0[1:-1], semi[1:-1], rtol=1e-10, atol=1e-12)
    # And it is what scipy's stiff integrator returns for the method of lines.
    mol = impl.method_of_lines_oracle(u0, dx, ALPHA, t_final)
    npt.assert_allclose(mol, semi, atol=1e-7)
    # Each marcher with a small step is within the time error of the semi-discrete solution and within
    # the spatial error of the exact one.
    dt = 0.25 * dx ** 2 / ALPHA
    steps = int(round(t_final / dt))
    for scheme in ("ftcs", "btcs", "cn"):
        u = impl.march(u0, dx, dt, ALPHA, steps, scheme)
        assert np.max(np.abs(u - semi)) < 2e-3
        assert np.max(np.abs(u - exact)) < 3e-3
    # Spatial error of the grid itself: (k dx)^2 / 12 relative to the decay exponent.
    k = np.pi / L
    predicted = np.exp(-1.0) * (1.0 - np.exp(-ALPHA * ((4 / dx ** 2) * np.sin(0.5 * k * dx) ** 2 - k * k) * t_final))
    assert np.max(np.abs(semi - exact)) == pytest.approx(abs(predicted), rel=1e-6)


def test_observed_spatial_convergence_order_is_two_at_fixed_r_and_four_at_r_one_sixth():
    t_final = L ** 2 / (ALPHA * np.pi ** 2)
    ns = [10, 20, 40, 80]
    e_f = impl.spatial_convergence("ftcs", ns, L, ALPHA, t_final, 0.4)
    e_b = impl.spatial_convergence("btcs", ns, L, ALPHA, t_final, 0.4)
    e_c = impl.spatial_convergence("cn", [20, 40, 80, 160], L, ALPHA, t_final, 0.5)
    for e in (e_f, e_b, e_c):
        assert np.all(np.diff(e) < 0)
        npt.assert_allclose(impl.observed_order(e), 2.0, atol=0.1)         # second order: error / 4 per halving
    e_magic = impl.spatial_convergence("ftcs", ns, L, ALPHA, t_final, 1.0 / 6.0)
    npt.assert_allclose(impl.observed_order(e_magic), 4.0, atol=0.1)       # the dx^2 and dt terms cancel at r = 1/6
    assert e_magic[0] < 1e-2 * e_f[0]


def test_observed_temporal_convergence_order_is_one_for_btcs_and_two_for_crank_nicolson():
    t_final = L ** 2 / (ALPHA * np.pi ** 2)
    dts = [t_final / 5, t_final / 10, t_final / 20, t_final / 40]
    e_b = impl.temporal_convergence("btcs", dts, 200, L, ALPHA, t_final)
    e_c = impl.temporal_convergence("cn", dts, 200, L, ALPHA, t_final)
    npt.assert_allclose(impl.observed_order(e_b), 1.0, atol=0.1)
    npt.assert_allclose(impl.observed_order(e_c), 2.0, atol=0.1)
    assert np.all(e_c < e_b)
    # BTCS decays each mode too slowly (1/(1 + z) > e^{-z}), FTCS too fast (1 - z < e^{-z}): check the signs.
    x, dx = _grid(200)
    semi = impl.semi_discrete_exact(x, t_final, L, ALPHA, 1)
    u_b = impl.march(impl.heat_exact_mode(x, 0.0, L, ALPHA, 1), dx, dts[0], ALPHA, 5, "btcs")
    assert np.max(u_b) > np.max(semi)
    dt_f = 0.4 * dx ** 2 / ALPHA
    u_f = impl.march(impl.heat_exact_mode(x, 0.0, L, ALPHA, 1), dx, dt_f, ALPHA, int(round(t_final / dt_f)), "ftcs")
    assert np.max(u_f) < np.max(impl.semi_discrete_exact(x, int(round(t_final / dt_f)) * dt_f, L, ALPHA, 1))
    with pytest.raises(ValueError):
        impl.temporal_convergence("cn", [t_final / 3.3], 50, L, ALPHA, t_final)
    orders = impl.observed_order(np.array([1.0, 0.25, 0.0625]))
    npt.assert_allclose(orders, [2.0, 2.0])


# ---------------------------------------------------------------- the engineering answer

def test_time_to_limit_agrees_across_schemes_and_with_the_series():
    c = impl.titanium_bulkhead()
    t_exact = impl.time_to_limit_exact(c)
    assert 1.0 < t_exact < 60.0
    alpha = impl.thermal_diffusivity(c["k"], c["rho"], c["cp"])
    n = 50
    dx = c["L"] / n
    t_f, steps_f = impl.time_to_limit(c, "ftcs", n, 0.45 * dx ** 2 / alpha)
    t_c, steps_c = impl.time_to_limit(c, "cn", n, 0.05)
    t_b, steps_b = impl.time_to_limit(c, "btcs", n, 0.05)
    assert t_f == pytest.approx(t_exact, rel=5e-3)
    assert t_c == pytest.approx(t_exact, rel=5e-3)
    assert t_b == pytest.approx(t_exact, rel=3e-2)                         # first order in dt: visibly worse
    assert steps_c < steps_f / 5                                           # the point of an implicit scheme
    assert t_b > t_c                                                       # BTCS decays the transient too slowly, so the limit arrives later
    # Mid-plane steady temperature is (T_hot + T_cold)/2 = 500 K > T_limit, so the limit is reached;
    # raise it above 500 K and it never is.
    c_safe = dict(c, T_limit=520.0)
    with pytest.raises(RuntimeError):
        impl.time_to_limit(c_safe, "cn", n, 0.5, t_max=50.0)
    with pytest.raises(ValueError):
        impl.time_to_limit(c, "cn", n, 0.5, x_probe=0.37 * dx)
    # A probe nearer the hot face reaches the limit sooner.
    t_near, _ = impl.time_to_limit(c, "cn", n, 0.05, x_probe=0.3 * c["L"])
    assert t_near < t_c


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "heat_fd_schemes.png").stat().st_size > 1000
