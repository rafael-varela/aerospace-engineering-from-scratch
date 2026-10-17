"""Public tests for P02.11 PDEs: Heat, Wave & Laplace Equations. They assert properties: PDE residuals,
boundary conditions, agreement between the series and an independent method, conservation and convergence.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-02/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

TILE = impl.reentry_tile()
CABLE = impl.sensor_cable()


def triangle(L: float, a: float, h: float):
    return lambda x: np.where(x < a, h * x / a, h * (L - x) / (L - a))


# ---------------------------------------------------------------- sine coefficients

def test_sine_coefficients_reproduce_known_series():
    L = 2.0
    b = impl.sine_coefficients(lambda x: np.ones_like(x), L, 8)
    n = np.arange(1, 9)
    npt.assert_allclose(b, np.where(n % 2 == 1, 4.0 / (n * np.pi), 0.0), atol=1e-6)   # the square wave (trapezoid quadrature)
    b = impl.sine_coefficients(lambda x: np.sin(3 * np.pi * x / L), L, 6)
    npt.assert_allclose(b, [0, 0, 1, 0, 0, 0], atol=1e-8)                              # orthogonality
    # Plucked string: closed form equals quadrature.
    Lc, a, h = CABLE["L"], CABLE["a"], CABLE["h"]
    npt.assert_allclose(impl.plucked_string_coefficients(Lc, a, h, 30),
                        impl.sine_coefficients(triangle(Lc, a, h), Lc, 30), atol=1e-9)
    with pytest.raises(ValueError):
        impl.sine_coefficients(np.ones_like, -1.0, 3)
    with pytest.raises(ValueError):
        impl.plucked_string_coefficients(1.0, 1.5, 0.01, 3)


# ---------------------------------------------------------------- the heat equation

def test_heat_series_satisfies_pde_boundary_and_initial_conditions():
    L, alpha = 0.5, 1e-5
    b = impl.sine_coefficients(lambda x: x * (L - x), L, 40)
    x = np.linspace(0.0, L, 101)
    npt.assert_allclose(impl.heat_series(x, 0.0, L, alpha, b), x * (L - x), atol=1e-5)   # initial data
    t = 2000.0
    assert abs(impl.heat_series(np.array([0.0, L]), t, L, alpha, b)).max() < 1e-12       # u = 0 at the ends
    # PDE residual by central differences in x and t on the series.
    h, k = 1e-3, 1.0
    xi = np.linspace(0.05, 0.45, 9)
    u = lambda xx, tt: impl.heat_series(xx, tt, L, alpha, b)
    u_t = (u(xi, t + k) - u(xi, t - k)) / (2 * k)
    u_xx = (u(xi + h, t) - 2 * u(xi, t) + u(xi - h, t)) / h ** 2
    npt.assert_allclose(u_t, alpha * u_xx, atol=1e-9)
    # Long-time decay is governed by the first mode: u ~ b1 sin(pi x / L) exp(-alpha (pi/L)^2 t).
    t_long = 3.0 * L ** 2 / (alpha * np.pi ** 2)
    ratio = impl.heat_series(np.array([L / 2]), t_long, L, alpha, b)[0] / (b[0] * np.exp(-alpha * (np.pi / L) ** 2 * t_long))
    assert ratio == pytest.approx(1.0, rel=1e-6)
    with pytest.raises(ValueError):
        impl.heat_series(x, -1.0, L, alpha, b)


def test_tile_series_limits_and_the_one_term_approximation():
    alpha, L, T0, Ts = TILE["alpha"], TILE["L"], TILE["T0"], TILE["Ts"]
    x = np.linspace(0.0, L, 51)
    # t = 0: the interior is at T0 (Gibbs only at the hot face itself).
    npt.assert_allclose(impl.tile_temperature(x[2:], 0.0, L, alpha, T0, Ts, n_terms=400), T0, atol=15.0)
    # Hot face always at Ts; insulated back face has zero slope.
    for t in (10.0, 600.0, 1200.0):
        T = impl.tile_temperature(x, t, L, alpha, T0, Ts)
        assert T[0] == pytest.approx(Ts, abs=1e-9)
        h = 1e-6
        slope = (impl.tile_temperature(np.array([L]), t, L, alpha, T0, Ts)[0]
                 - impl.tile_temperature(np.array([L - h]), t, L, alpha, T0, Ts)[0]) / h
        assert abs(slope) < 1e-3 * (Ts - T0) / L
        assert np.all(np.diff(T) <= 1e-9)                                              # monotone from hot to back
        assert np.all((T >= T0 - 1e-9) & (T <= Ts + 1e-9))                            # maximum principle
    # Very long time: everything at Ts.
    npt.assert_allclose(impl.tile_temperature(x, 50 * L ** 2 / alpha, L, alpha, T0, Ts), Ts, atol=1e-6)
    # Once the Fourier number exceeds ~0.2 the first term carries the back face to 0.5 %.
    t = TILE["t_pulse"]
    one = Ts + (T0 - Ts) * 4 / np.pi * np.exp(-alpha * (np.pi / (2 * L)) ** 2 * t)
    assert impl.tile_temperature(np.array([L]), t, L, alpha, T0, Ts)[0] == pytest.approx(one, rel=5e-3)
    # Scaling: the solution depends on x/L and alpha t / L^2 only.
    T_a = impl.tile_temperature(np.array([0.4 * L]), 500.0, L, alpha, T0, Ts)[0]
    T_b = impl.tile_temperature(np.array([0.8 * L]), 2000.0, 2 * L, alpha, T0, Ts)[0]
    assert T_a == pytest.approx(T_b, rel=1e-10)


def test_tile_thickness_for_limit_brackets_and_inverts():
    alpha, T0, Ts, t = TILE["alpha"], TILE["T0"], TILE["Ts"], TILE["t_pulse"]
    L_req = impl.tile_thickness_for_limit(t, alpha, T0, Ts, TILE["T_limit"])
    assert impl.tile_temperature(np.array([L_req]), t, L_req, alpha, T0, Ts)[0] == pytest.approx(TILE["T_limit"], abs=1e-6)
    assert L_req > TILE["L"]                                                           # the 50 mm tile is too thin here
    # A hotter limit needs less tile; four times the time needs twice the thickness (Fourier scaling).
    assert impl.tile_thickness_for_limit(t, alpha, T0, Ts, 600.0) < L_req
    assert impl.tile_thickness_for_limit(4 * t, alpha, T0, Ts, TILE["T_limit"]) == pytest.approx(2 * L_req, rel=1e-6)
    with pytest.raises(ValueError):
        impl.tile_thickness_for_limit(t, alpha, T0, Ts, 2000.0)


def test_ftcs_matches_the_series_converges_at_second_order_and_refuses_unstable_steps():
    alpha, L, T0, Ts, t = TILE["alpha"], TILE["L"], TILE["T0"], TILE["Ts"], 600.0
    errs = []
    for nx in (26, 51, 101):
        x = np.linspace(0.0, L, nx)
        dx = x[1] - x[0]
        n_steps = int(np.ceil(t / (0.25 * dx ** 2 / alpha)))
        dt = t / n_steps
        u = impl.heat_ftcs(np.full(nx, T0), dx, dt, alpha, n_steps, left=("dirichlet", Ts), right=("neumann", 0.0))
        errs.append(np.max(np.abs(u - impl.tile_temperature(x, t, L, alpha, T0, Ts))))
    assert errs[-1] < 0.5                                                              # kelvin, on a 1000 K range
    assert errs[0] / errs[1] > 3.0 and errs[1] / errs[2] > 3.0                         # ~ second order in dx (dt ~ dx^2)
    # Dirichlet-Dirichlet slab against heat_series.
    Ld = 0.5
    b = impl.sine_coefficients(lambda xx: np.sin(np.pi * xx / Ld), Ld, 5)
    x = np.linspace(0.0, Ld, 41)
    dx = x[1] - x[0]
    dt = 0.4 * dx ** 2 / 1e-4
    u = impl.heat_ftcs(np.sin(np.pi * x / Ld), dx, dt, 1e-4, 200)
    npt.assert_allclose(u, impl.heat_series(x, 200 * dt, Ld, 1e-4, b), atol=2e-3)
    with pytest.raises(ValueError):
        impl.heat_ftcs(np.full(11, T0), dx, 0.6 * dx ** 2 / alpha, alpha, 1)
    with pytest.raises(ValueError):
        impl.heat_ftcs(np.full(11, T0), dx, 0.1 * dx ** 2 / alpha, alpha, 1, left=("robin", 0.0))
    # Zero steps returns the initial data; the sum (heat content) of an insulated-insulated slab is conserved.
    npt.assert_array_equal(impl.heat_ftcs(np.full(11, T0), dx, dt, alpha, 0), np.full(11, T0))
    u0 = np.linspace(300.0, 400.0, 41)
    u = impl.heat_ftcs(u0, dx, dt, 1e-4, 300, left=("neumann", 0.0), right=("neumann", 0.0))
    assert np.trapezoid(u, dx=dx) == pytest.approx(np.trapezoid(u0, dx=dx), rel=1e-6)


# ---------------------------------------------------------------- the wave equation

def test_mode_frequencies_scale_with_n_tension_and_length():
    L, T, mu = CABLE["L"], CABLE["T"], CABLE["mu"]
    f1 = impl.string_mode_frequency(1, L, T, mu)
    assert f1 == pytest.approx(np.sqrt(T / mu) / (2 * L))
    assert impl.string_mode_frequency(3, L, T, mu) == pytest.approx(3 * f1)
    assert impl.string_mode_frequency(1, L, 4 * T, mu) == pytest.approx(2 * f1)
    assert impl.string_mode_frequency(1, 2 * L, T, mu) == pytest.approx(f1 / 2)
    with pytest.raises(ValueError):
        impl.string_mode_frequency(0, L, T, mu)


def test_wave_series_equals_dalembert_and_satisfies_the_pde():
    L, T, mu, a, h = CABLE["L"], CABLE["T"], CABLE["mu"], CABLE["a"], CABLE["h"]
    c = np.sqrt(T / mu)
    f = triangle(L, a, h)
    b = impl.plucked_string_coefficients(L, a, h, 200)
    x = np.linspace(0.0, L, 301)
    for t in (0.0, 0.7e-3, 2.3e-3, 11.0e-3):
        npt.assert_allclose(impl.wave_series(x, t, L, c, b), impl.dalembert(x, t, c, L, f), atol=2e-5)   # 0.4 % of h
    # Period 2L/c: the string returns to its shape; at half a period it is inverted.
    P = 2 * L / c
    npt.assert_allclose(impl.dalembert(x, P, c, L, f), f(x), atol=1e-12)
    npt.assert_allclose(impl.dalembert(x, P / 2, c, L, f), -f(L - x), atol=1e-12)
    # Fixed ends.
    assert abs(impl.wave_series(np.array([0.0, L]), 1e-3, L, c, b)).max() < 1e-12
    # PDE residual on a smooth single-mode solution.
    b1 = np.array([1e-3])
    xi, t, dx, dt = np.linspace(0.1, 1.1, 7), 2e-3, 1e-4, 1e-7
    u = lambda xx, tt: impl.wave_series(xx, tt, L, c, b1)
    u_tt = (u(xi, t + dt) - 2 * u(xi, t) + u(xi, t - dt)) / dt ** 2
    u_xx = (u(xi + dx, t) - 2 * u(xi, t) + u(xi - dx, t)) / dx ** 2
    npt.assert_allclose(u_tt, c ** 2 * u_xx, rtol=1e-4)
    # Initial velocity data: u_t(x, 0) = sum d_n sin, displacement zero.
    d = np.array([0.0, 0.2, 0.0])
    npt.assert_allclose(impl.wave_series(x, 0.0, L, c, np.zeros(3), d), 0.0, atol=1e-15)
    w2 = 2 * np.pi * c / L
    npt.assert_allclose(impl.wave_series(x, 1e-4, L, c, np.zeros(3), d), 0.2 / w2 * np.sin(w2 * 1e-4) * np.sin(2 * np.pi * x / L), atol=1e-14)


def test_odd_periodic_extension_has_the_right_symmetries():
    L = 1.2
    F = impl.odd_periodic_extension(triangle(L, 0.3, 0.005), L)
    x = np.linspace(-3 * L, 3 * L, 601)
    npt.assert_allclose(F(-x), -F(x), atol=1e-15)
    npt.assert_allclose(F(x + 2 * L), F(x), atol=1e-15)
    assert F(np.array([0.0, L, 2 * L]))[0] == 0.0 and abs(F(np.array([L]))[0]) < 1e-15


def test_leapfrog_tracks_dalembert_conserves_energy_and_enforces_cfl():
    L, T, mu, a, h = CABLE["L"], CABLE["T"], CABLE["mu"], CABLE["a"], CABLE["h"]
    c = np.sqrt(T / mu)
    f = triangle(L, a, h)
    nx = 241
    x = np.linspace(0.0, L, nx)
    dx = x[1] - x[0]
    dt = 0.8 * dx / c
    E0 = impl.string_energy(f(x), np.zeros(nx), dx, T, mu)
    for n_steps in (50, 300):
        u, ut = impl.wave_leapfrog(f(x), np.zeros(nx), dx, dt, c, n_steps)
        npt.assert_allclose(u, impl.dalembert(x, n_steps * dt, c, L, f), atol=1e-4)     # 2 % of h
        assert impl.string_energy(u, ut, dx, T, mu) == pytest.approx(E0, rel=5e-3)
        assert u[0] == 0.0 and u[-1] == 0.0
    # Courant number exactly 1 reproduces the exact solution on the grid (a classic property of leapfrog).
    dt1 = dx / c
    u, _ = impl.wave_leapfrog(f(x), np.zeros(nx), dx, dt1, c, 100)
    npt.assert_allclose(u, impl.dalembert(x, 100 * dt1, c, L, f), atol=1e-12)
    with pytest.raises(ValueError):
        impl.wave_leapfrog(f(x), np.zeros(nx), dx, 1.01 * dx / c, c, 10)
    # Energy of the triangle: exact (1/2) T h^2 (1/a + 1/(L-a)) within the gradient's kink error.
    assert E0 == pytest.approx(0.5 * T * h ** 2 * (1 / a + 1 / (L - a)), rel=0.02)


# ---------------------------------------------------------------- Laplace's equation

def test_laplace_series_is_harmonic_and_meets_its_boundary_conditions():
    a, b = 0.3, 0.2
    f_top = lambda x: 60.0 * 4 * x * (a - x) / a ** 2
    bt = impl.sine_coefficients(f_top, a, 40)
    x = np.linspace(0.0, a, 61)
    y = np.linspace(0.0, b, 41)
    X, Y = np.meshgrid(x, y, indexing="ij")
    u = impl.laplace_rectangle(X, Y, a, b, bt)
    npt.assert_allclose(u[:, -1], f_top(x), atol=0.05)                                 # top edge (40 terms)
    npt.assert_allclose(u[:, 0], 0.0, atol=1e-12)
    npt.assert_allclose(u[0, :], 0.0, atol=1e-12)
    npt.assert_allclose(u[-1, :], 0.0, atol=1e-12)
    # Harmonic: the five-point Laplacian on a fine grid is small compared with u / h^2.
    xf = np.linspace(0.05, 0.25, 201)
    yf = np.linspace(0.02, 0.18, 161)
    Xf, Yf = np.meshgrid(xf, yf, indexing="ij")
    res = impl.laplacian_residual(impl.laplace_rectangle(Xf, Yf, a, b, bt), xf[1] - xf[0], yf[1] - yf[0])
    assert np.max(np.abs(res)) < 1e-3 * 60.0 / (xf[1] - xf[0]) ** 2
    # Maximum principle and mean-value property at the centre.
    assert 0.0 <= u.min() and u.max() <= 60.0 + 1e-9
    r = 0.05
    th = np.linspace(0, 2 * np.pi, 721)[:-1]
    ring = impl.laplace_rectangle(a / 2 + r * np.cos(th), b / 2 + r * np.sin(th), a, b, bt)
    assert np.mean(ring) == pytest.approx(impl.laplace_rectangle(a / 2, b / 2, a, b, bt), rel=1e-6)
    # Each term alone is harmonic: sin(kx) sinh(ky) exactly.
    one = impl.laplace_rectangle(Xf, Yf, a, b, np.array([1.0]))
    npt.assert_allclose(one, np.sin(np.pi * Xf / a) * np.sinh(np.pi * Yf / a) / np.sinh(np.pi * b / a), rtol=1e-12)


def test_jacobi_relaxation_converges_to_the_series():
    a, b = 0.3, 0.2
    f_top = lambda x: 60.0 * 4 * x * (a - x) / a ** 2
    bt = impl.sine_coefficients(f_top, a, 40)
    errs, iters = [], []
    for nx, ny in ((16, 11), (31, 21)):
        u, it = impl.laplace_jacobi(nx, ny, a, b, f_top, tol=1e-10)
        X, Y = np.meshgrid(np.linspace(0, a, nx), np.linspace(0, b, ny), indexing="ij")
        errs.append(np.max(np.abs(u - impl.laplace_rectangle(X, Y, a, b, bt))))
        iters.append(it)
        assert np.max(np.abs(impl.laplacian_residual(u, a / (nx - 1), b / (ny - 1)))) < 1e-6 * 60.0 / (a / (nx - 1)) ** 2
    assert errs[1] < errs[0] / 3.0                                                     # second-order grid error
    assert iters[1] > 2.5 * iters[0]                                                   # Jacobi cost grows ~ n^2
    with pytest.raises(RuntimeError):
        impl.laplace_jacobi(31, 21, a, b, f_top, tol=1e-12, max_iter=5)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "pde_solutions.png").stat().st_size > 1000
