"""Public tests for P08.18 Incompressible CFD from Scratch (lid-driven cavity).

They assert numerics and physics, not one magic number: the Poisson solver converges at second order
to a manufactured solution and is exact for harmonic quadratics; Thom's wall vorticity converges at
first order; the convection operator is exact on quadratics; the discrete divergence of the velocity
field vanishes to round-off and the net flux through a centreline falls as h^2 (mass conservation);
Newton and pseudo-time marching reach the same discrete steady state; Newton converges in a handful of
iterations; the Stokes limit is symmetric and the Re = 100 vortex moves downstream; the solution obeys
Reynolds-number similarity; the centre values converge with observed order near 2 on 17, 33 and 65
grids; and the centreline profiles approach the published data of Ghia, Ghia & Shin (1982), Re = 100.

Grids are small (N <= 65) so the whole file runs in a few seconds.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

RE = 100.0


@pytest.fixture(scope="module")
def sols():
    """Newton solutions at Re = 100 on the three test grids (computed once for the module)."""
    return {N: impl.solve_cavity_newton(N, RE) for N in (17, 33, 65)}


# ---------------------------------------------------------------- the Poisson operator

def test_laplacian_matrix_is_symmetric_negative_definite_with_five_point_sparsity():
    N = 9
    A = impl.laplacian_matrix(N)
    n = N - 2
    assert A.shape == (n * n, n * n)
    assert A.nnz == 5 * n * n - 4 * n
    Ad = A.toarray()
    npt.assert_allclose(Ad, Ad.T, atol=0)
    ev = np.linalg.eigvalsh(Ad)
    h = 1.0 / (N - 1)
    # the extreme eigenvalues of the five-point Laplacian are -(8/h^2) sin^2(pi h/2) and -(8/h^2) cos^2(pi h/2)
    npt.assert_allclose(ev.max(), -8 / h ** 2 * np.sin(np.pi * h / 2) ** 2, rtol=1e-10)
    npt.assert_allclose(ev.min(), -8 / h ** 2 * np.cos(np.pi * h / 2) ** 2, rtol=1e-10)


def test_poisson_converges_at_second_order_to_a_manufactured_solution():
    errs = []
    for N in (9, 17, 33, 65):
        x, _ = impl.grid(N)
        X, Y = np.meshgrid(x, x, indexing="ij")
        exact = np.sin(np.pi * X) * np.sin(2 * np.pi * Y) + np.exp(X) * np.cos(Y)   # non-zero edges too
        f = -5 * np.pi ** 2 * np.sin(np.pi * X) * np.sin(2 * np.pi * Y)              # lap(e^x cos y) = 0
        phi = impl.solve_poisson(f, boundary=exact)
        errs.append(np.max(np.abs(phi - exact)))
    orders = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert np.all((orders > 1.9) & (orders < 2.1)), orders


def test_poisson_is_exact_for_a_harmonic_quadratic():
    x, _ = impl.grid(17)
    X, Y = np.meshgrid(x, x, indexing="ij")
    exact = X ** 2 - Y ** 2 + 3 * X * Y + 0.5          # harmonic; the five-point stencil is exact on it
    phi = impl.solve_poisson(np.zeros_like(X), boundary=exact)
    npt.assert_allclose(phi, exact, atol=1e-12)


# ---------------------------------------------------------------- velocities, walls, convection

def test_velocities_divergence_free_to_round_off_for_any_stream_function():
    rng = np.random.default_rng(18)
    psi = rng.standard_normal((21, 21))                 # even a random psi: the operators commute,
    psi[0, :] = psi[-1, :] = psi[:, 0] = psi[:, -1] = 0.0   # as long as psi = 0 on the walls (no through-flow)
    u, v = impl.velocities(psi)
    assert np.max(np.abs(impl.divergence(u, v))) < 1e-10 * np.max(np.abs(psi)) * 20 ** 2


def test_velocities_exact_on_quadratic_stream_function_with_wall_values():
    x, _ = impl.grid(11)
    X, Y = np.meshgrid(x, x, indexing="ij")
    psi = 2 * X ** 2 - X * Y + 3 * Y ** 2
    u, v = impl.velocities(psi, U=1.7)
    npt.assert_allclose(u[1:-1, 1:-1], (-X + 6 * Y)[1:-1, 1:-1], atol=1e-12)
    npt.assert_allclose(v[1:-1, 1:-1], -(4 * X - Y)[1:-1, 1:-1], atol=1e-12)
    assert np.all(u[1:-1, -1] == 1.7) and np.all(u[:, 0] == 0) and np.all(u[0, :-1] == 0)
    assert np.all(v[0, :] == 0) and np.all(v[-1, :] == 0) and np.all(v[:, -1] == 0)


def test_thom_formula_lid_term_and_first_order_wall_vorticity():
    # a resting fluid under a moving lid: omega_lid = -2U/h exactly
    N = 17
    w = impl.wall_vorticity(np.zeros((N, N)), U=2.0)
    npt.assert_allclose(w[1:-1, -1], -2 * 2.0 * (N - 1), rtol=1e-14)
    assert np.all(w[1:-1, 0] == 0) and np.all(w[1:-1, 1:-1] == 0)
    # psi = x^2(1-x)^2 y^2(1-y)^2 has psi = dpsi/dn = 0 on every wall; exact wall vorticity -lap(psi)
    errs = []
    for N in (17, 33, 65, 129):
        x, _ = impl.grid(N)
        X, Y = np.meshgrid(x, x, indexing="ij")
        gx, gy = X ** 2 * (1 - X) ** 2, Y ** 2 * (1 - Y) ** 2
        psi = gx * gy
        exact_bottom = -(2 * x ** 2 * (1 - x) ** 2)          # -psi_yy at y = 0 (psi_xx = 0 on the wall)
        wb = impl.wall_vorticity(psi, U=0.0)[1:-1, 0]
        errs.append(np.max(np.abs(wb - exact_bottom[1:-1])))
    orders = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert np.all((orders > 0.9) & (orders < 1.1)), orders


def test_convection_operator_exact_on_quadratics():
    x, _ = impl.grid(13)
    X, Y = np.meshgrid(x, x, indexing="ij")
    psi = X ** 2 + 2 * X * Y - Y ** 2                    # u = 2x - 2y, v = -(2x + 2y)
    omega = 3 * X ** 2 - X * Y + Y ** 2                  # omega_x = 6x - y, omega_y = -x + 2y
    exact = (2 * X - 2 * Y) * (6 * X - Y) - (2 * X + 2 * Y) * (-X + 2 * Y)
    npt.assert_allclose(impl.convection(psi, omega), exact[1:-1, 1:-1], atol=1e-11)


# ---------------------------------------------------------------- the solvers

def test_newton_converges_fast_and_satisfies_the_discrete_equations(sols):
    s = sols[33]
    assert s["iterations"] <= 6
    r1, r2 = impl.cavity_residual(s["psi"], s["omega"], RE)
    assert r1 < 1e-9 and r2 < 1e-9
    # the boundary conditions: psi = 0 on every wall
    assert np.all(s["psi"][0, :] == 0) and np.all(s["psi"][-1, :] == 0)
    assert np.all(s["psi"][:, 0] == 0) and np.all(s["psi"][:, -1] == 0)


def test_marching_and_newton_reach_the_same_steady_state(sols):
    m = impl.march_cavity(17, RE)
    n = sols[17]
    assert np.max(np.abs(m["psi"] - n["psi"])) < 1e-7
    assert np.max(np.abs(m["omega"] - n["omega"])) < 1e-5 * np.max(np.abs(n["omega"]))


def test_mass_conservation_divergence_and_centreline_flux(sols):
    imbalance = []
    for N in (17, 33, 65):
        s = sols[N]
        u, v = impl.velocities(s["psi"])
        assert np.max(np.abs(impl.divergence(u, v))) < 1e-12        # discrete continuity, round-off
        y, uc, x, vc = impl.centreline_profiles(s)
        imbalance.append(abs(np.trapezoid(uc, y)))                   # net flow through x = 0.5: zero exactly
    orders = np.log2(np.array(imbalance[:-1]) / np.array(imbalance[1:]))
    assert np.all(orders > 1.8), orders                              # the trapezoid flux error falls as h^2
    assert imbalance[-1] < 2e-3


def test_stokes_limit_is_symmetric_and_re100_vortex_moves_downstream(sols):
    st = impl.solve_cavity_newton(33, 1e-6)
    assert np.max(np.abs(st["psi"] - st["psi"][::-1, :])) < 1e-9     # mirror symmetry about x = 0.5
    _, xs, _ = impl.vortex_centre(st)
    assert abs(xs - 0.5) < 1e-6
    pmin, xc, yc = impl.vortex_centre(sols[65])
    assert 0.58 < xc < 0.66 and 0.70 < yc < 0.78                     # downstream of the centreline


def test_reynolds_similarity_scaling_with_lid_speed():
    # dimensional sanity: with Re fixed, U*L/nu is what matters; doubling U and halving the code's
    # 1/nu parameter is the same flow with psi and omega scaled by 2
    a = impl.solve_cavity_newton(17, 100.0, U=1.0)
    b = impl.solve_cavity_newton(17, 50.0, U=2.0)
    npt.assert_allclose(b["psi"], 2 * a["psi"], atol=1e-10)
    npt.assert_allclose(b["omega"], 2 * a["omega"], rtol=1e-8, atol=1e-8)


def test_grid_convergence_with_observed_order_near_two(sols):
    psis = [impl.centre_values(sols[N])[0] for N in (17, 33, 65)]
    us = [impl.centre_values(sols[N])[1] for N in (17, 33, 65)]
    p_psi, ext_psi = impl.observed_order(*psis)
    p_u, ext_u = impl.observed_order(*us)
    assert 1.7 < p_psi < 2.3 and 1.7 < p_u < 2.3
    # the extrapolate is closer to the fine-grid value than the fine grid is to the middle one
    assert abs(ext_u - us[2]) < abs(us[2] - us[1])
    # the vortex-centre stream function converges monotonically towards about -0.1035
    pm = [impl.vortex_centre(sols[N])[0] for N in (17, 33, 65)]
    assert pm[0] > pm[1] > pm[2] > -0.1040


def test_centreline_profiles_approach_ghia_published_data(sols):
    # PUBLISHED DATA: Ghia, Ghia & Shin (1982), Tables I and II, Re = 100, a 129 x 129 solution
    g17, g33, g65 = (impl.compare_with_ghia(sols[N]) for N in (17, 33, 65))
    assert g17["max_du"] > g33["max_du"] > g65["max_du"]
    assert g33["max_du"] < 0.012 and g33["max_dv"] < 0.012          # about 1 % of the lid speed
    assert g65["max_du"] < 0.005 and g65["max_dv"] < 0.008
    pmin, _, _ = impl.vortex_centre(sols[65])
    assert pmin == pytest.approx(impl.GHIA_RE100_PSI_MIN, rel=0.01)  # Table V: -0.103423


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path, grids=(17, 33, 65))
    assert (tmp_path / "cavity.png").stat().st_size > 10_000
