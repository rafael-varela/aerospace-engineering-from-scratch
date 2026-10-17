"""Public tests for P08.10 Boundary Layers: Laminar. They assert physics, not magic numbers: the
RK4 integrator is fourth-order; shooting satisfies all three boundary conditions and agrees with an
independent collocation solve (scipy.integrate.solve_bvp); the converged profile satisfies the
Blasius equation when substituted back; the wall curvature matches the published value
f''(0) = 0.332057 (labelled published); the momentum integral theta = 2 f''(0) holds; the
dimensional layer scales as sqrt(x) and is invariant to a change of units; and a marching solution
of the boundary-layer PDEs, which never assumes similarity, collapses onto the Blasius profile.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

FPP0_PUBLISHED = 0.332057          # published Blasius value (Howarth 1938; many later sources)


@pytest.fixture(scope="module")
def shot():
    return impl.shoot_blasius(10.0, 1000)


@pytest.fixture(scope="module")
def bvp():
    return impl.blasius_bvp(10.0)


# ---------------------------------------------------------------- the ODE and RK4

def test_rhs_is_the_blasius_system_and_vectorised():
    Y = np.array([0.3, 0.5, 0.2])
    npt.assert_allclose(impl.blasius_rhs(0.0, Y), [0.5, 0.2, -0.5 * 0.3 * 0.2])
    Ym = np.random.default_rng(1).random((3, 7))
    out = impl.blasius_rhs(0.0, Ym)
    assert out.shape == (3, 7)
    npt.assert_allclose(out[2], -0.5 * Ym[0] * Ym[2])


def test_rk4_starts_at_the_wall_and_has_the_right_shape():
    eta, Y = impl.rk4_blasius(0.33, 5.0, 50)
    assert eta.shape == (51,) and Y.shape == (51, 3)
    npt.assert_allclose(Y[0], [0.0, 0.0, 0.33])
    assert eta[-1] == pytest.approx(5.0)


def test_rk4_is_fourth_order():
    # integrate a fixed initial value problem; compare with a fine run (no root finding involved)
    ref = impl.rk4_blasius(0.332, 6.0, 3200)[1][-1]
    errs = [np.max(np.abs(impl.rk4_blasius(0.332, 6.0, n)[1][-1] - ref)) for n in (20, 40, 80)]
    orders = impl.observed_order(errs)
    assert np.all(orders > 3.5) and np.all(orders < 4.5), orders


def test_zero_wall_curvature_gives_no_flow():
    # f''(0) = 0 with f = f' = 0 is the trivial solution: the wall never sets the fluid moving
    _, Y = impl.rk4_blasius(0.0, 8.0, 100)
    assert np.all(Y == 0.0)


# ---------------------------------------------------------------- shooting against solve_bvp

def test_shooting_meets_all_boundary_conditions(shot):
    fpp0, eta, Y = shot
    assert Y[0, 0] == 0.0 and Y[0, 1] == 0.0
    assert Y[-1, 1] == pytest.approx(1.0, abs=1e-10)
    assert Y[0, 2] == pytest.approx(fpp0)


def test_shooting_matches_published_and_collocation(shot, bvp):
    assert shot[0] == pytest.approx(FPP0_PUBLISHED, abs=2e-6)          # published value, 6 figures
    assert shot[0] == pytest.approx(bvp[0], abs=1e-7)                  # two independent methods


def test_profiles_from_the_two_methods_agree(shot, bvp):
    _, eta, Y = shot
    _, eb, Yb = bvp
    for k in range(3):
        npt.assert_allclose(np.interp(eb, eta, Y[:, k]), Yb[:, k], atol=2e-5)   # onto the bvp mesh


def test_solution_satisfies_the_blasius_equation_when_substituted(shot):
    # f''' from a finite difference of f'' must equal -f f''/2 everywhere
    _, eta, Y = shot
    h = eta[1] - eta[0]
    fppp = (Y[2:, 2] - Y[:-2, 2]) / (2 * h)
    resid = fppp + 0.5 * Y[1:-1, 0] * Y[1:-1, 2]
    assert np.max(np.abs(resid)) < 1e-5


def test_profile_is_monotone_and_shear_decays(shot):
    _, eta, Y = shot
    assert np.all(np.diff(Y[:, 1]) > -1e-14)           # u/U rises monotonically to 1
    assert np.all(Y[:, 2] > 0) and Y[-1, 2] < 1e-6 * Y[0, 2]
    assert np.all(np.diff(Y[:, 2]) <= 1e-14)           # f'' decreases: no inflection (favourable to stability)


def test_domain_truncation_converges_from_above():
    vals = [impl.shoot_blasius(em, int(100 * em))[0] for em in (5.0, 6.0, 8.0)]
    assert vals[0] > vals[1] > vals[2] > FPP0_PUBLISHED - 1e-6       # too short a domain forces a steeper wall layer
    assert vals[2] - FPP0_PUBLISHED < 1e-5


# ---------------------------------------------------------------- integral thicknesses

def test_integral_thicknesses_and_momentum_integral(shot, bvp):
    c = impl.blasius_integrals(shot[1], shot[2])
    assert c["delta99"] == pytest.approx(4.91, abs=0.01)
    assert c["delta_star"] == pytest.approx(1.7208, abs=2e-4)
    assert c["H"] == pytest.approx(2.591, abs=1e-3)
    # momentum integral for a flat plate: d(theta)/dx = c_f / 2 with theta ~ sqrt(x) gives theta = 2 f''(0)
    assert c["theta"] == pytest.approx(2 * c["fpp0"], rel=1e-4)
    # outward displacement velocity at the edge: (eta f' - f)/2 -> delta*/2
    assert c["v_edge"] == pytest.approx(c["delta_star"] / 2, rel=1e-4)
    cb = impl.blasius_integrals(bvp[1], bvp[2])
    assert cb["theta"] == pytest.approx(c["theta"], rel=1e-3) and cb["H"] == pytest.approx(c["H"], rel=1e-3)


def test_dimensional_layer_scales_as_root_x_and_is_unit_invariant(shot):
    c = impl.blasius_integrals(shot[1], shot[2])
    U, nu, rho = 60.0, 1.863e-5, 0.9091
    x = np.array([0.1, 0.4, 1.6])
    bl = impl.plate_boundary_layer(x, U, nu, rho, c)
    npt.assert_allclose(bl["delta99"][1:] / bl["delta99"][:-1], 2.0, rtol=1e-12)   # 4x the distance, 2x the thickness
    npt.assert_allclose(bl["cf"][1:] / bl["cf"][:-1], 0.5, rtol=1e-12)
    npt.assert_allclose(bl["cf"] * np.sqrt(bl["Re_x"]), 2 * c["fpp0"], rtol=1e-12)
    npt.assert_allclose(bl["tau_w"], 0.5 * rho * U ** 2 * bl["cf"], rtol=1e-12)
    # the same plate in millimetres and milliseconds: lengths scale by 1000, Re and c_f do not change
    bl2 = impl.plate_boundary_layer(x * 1e3, U * 1e3 / 1e3, nu * 1e6 / 1e3, rho, c)
    npt.assert_allclose(bl2["Re_x"], bl["Re_x"], rtol=1e-12)
    npt.assert_allclose(bl2["delta99"], bl["delta99"] * 1e3, rtol=1e-12)
    npt.assert_allclose(bl2["cf"], bl["cf"], rtol=1e-12)


def test_integrated_drag_equals_integral_of_local_shear(shot):
    c = impl.blasius_integrals(shot[1], shot[2])
    U, nu, rho, L, b = 60.0, 1.863e-5, 0.9091, 1.5, 10.0
    D, CF = impl.friction_drag(U, L, b, rho, nu, c["fpp0"], sides=2)
    # integrate tau_w along the chord numerically (substitute x = L s^2 to remove the 1/sqrt(x) singularity)
    s = np.linspace(0, 1, 2001)[1:]
    x = L * s ** 2
    tau = impl.plate_boundary_layer(x, U, nu, rho, c)["tau_w"]
    integrand = tau * 2 * L * s                                     # tau dx = tau 2 L s ds; finite at s -> 0
    s0 = np.concatenate([[0.0], s])
    i0 = 0.5 * rho * U ** 2 * 2 * c["fpp0"] * 2 * np.sqrt(L * nu / U)  # limit of the integrand at s = 0
    D_num = 2 * b * np.trapezoid(np.concatenate([[i0], integrand]), s0)
    assert D == pytest.approx(D_num, rel=1e-6)
    assert CF == pytest.approx(2 * impl.plate_boundary_layer(L, U, nu, rho, c)["cf"], rel=1e-12)
    assert impl.friction_drag(U, L, b, rho, nu, c["fpp0"], sides=1)[0] == pytest.approx(D / 2)


# ---------------------------------------------------------------- independent check: marching the PDEs

def test_marched_pde_collapses_onto_blasius(shot):
    _, eta_b, Y = shot
    U, nu, L = 10.0, 1.4607e-5, 1.0
    d99 = 5.0 * L / np.sqrt(U * L / nu)
    x, y, u, v = impl.march_boundary_layer(U, nu, L, 3 * d99, ny=200, nx=800)
    assert u.shape == (801, 201) and np.all(u[:, 0] == 0) and np.all(u[:, -1] == U)
    for xf in (0.3, 1.0):
        j = int(np.argmin(abs(x - xf)))
        eta = y * np.sqrt(U / (nu * x[j]))
        npt.assert_allclose(u[j] / U, np.interp(eta, eta_b, Y[:, 1], right=1.0), atol=6e-3)
        cf = impl.wall_cf_from_profile(y, u[j], U, nu)
        assert cf * np.sqrt(U * x[j] / nu) == pytest.approx(2 * FPP0_PUBLISHED, rel=6e-3)


def test_marching_error_is_first_order_in_x():
    U, nu, L = 10.0, 1.4607e-5, 1.0
    d99 = 5.0 * L / np.sqrt(U * L / nu)
    errs = []
    for nx in (100, 200, 400):
        x, y, u, v = impl.march_boundary_layer(U, nu, L, 3 * d99, ny=300, nx=nx)
        errs.append(abs(impl.wall_cf_from_profile(y, u[-1], U, nu) * np.sqrt(U * L / nu) - 2 * FPP0_PUBLISHED))
    orders = impl.observed_order(errs)
    assert np.all(orders > 0.8) and np.all(orders < 1.3), orders


def test_marched_flow_conserves_mass():
    # continuity: the mass flux deficit int (U - u) dy grows exactly as U v at the top of the domain requires
    U, nu, L = 10.0, 1.4607e-5, 1.0
    d99 = 5.0 * L / np.sqrt(U * L / nu)
    x, y, u, v = impl.march_boundary_layer(U, nu, L, 3 * d99, ny=200, nx=400)
    deficit = np.trapezoid(U - u, y, axis=1)
    # d/dx int u dy = -v_top  (u = U at the top)  =>  d(deficit)/dx = v_top / U * U
    lhs = np.diff(deficit) / np.diff(x)
    npt.assert_allclose(lhs, v[1:, -1], rtol=1e-8, atol=1e-12)
    assert np.all(v[1:, -1] > 0)                           # the layer pushes fluid outward


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "blasius.png").stat().st_size > 10_000
