"""Public tests for P08.06 Navier-Stokes Equations. They assert fluid mechanics, not magic numbers:
each exact solution is substituted back into the full Navier-Stokes equations by SymPy and must
leave a zero residual (while a wrong profile must not); each satisfies its boundary conditions; its
flow rate equals the integral of its profile; the wall shear stresses close the momentum balance of
the gap; the finite-difference solvers reproduce the parabolic profiles to round-off (their
truncation error is proportional to u'''', which vanishes), converge at second order on a
manufactured problem and on the two unsteady problems, and relax to the steady state.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
import sympy as sp
from scipy.integrate import quad

from conftest import load_impl

impl = load_impl(__file__)

X, Y, T = impl.X, impl.Y, impl.T
h_, U_, G_, mu_, rho_ = 1.0e-3, 0.4, 2.0e3, 1.2e-2, 850.0


# ---------------------------------------------------------------- substitution into the equations

def test_couette_poiseuille_satisfies_navier_stokes_and_matches_the_function():
    h, U, G, mu, rho = sp.symbols("h U G mu rho", positive=True)
    u = U * Y / h + G * Y * (h - Y) / (2 * mu)
    res = impl.ns_residuals(u, sp.Integer(0), -G * X, rho, mu)
    assert all(sp.simplify(r) == 0 for r in res)
    # the boundary conditions: no slip on both walls
    assert sp.simplify(u.subs(Y, 0)) == 0 and sp.simplify(u.subs(Y, h) - U) == 0
    # the numerical function is the same profile
    y = np.linspace(0.0, h_, 7)
    f = sp.lambdify(Y, u.subs({h: h_, U: U_, G: G_, mu: mu_}), "numpy")
    npt.assert_allclose(impl.couette_poiseuille(y, h_, U_, G_, mu_), f(y), rtol=1e-12, atol=1e-15)


def test_residual_detects_a_wrong_profile_and_accepts_taylor_green():
    mu, rho, G = sp.symbols("mu rho G", positive=True)
    wrong = G * Y * (1 - Y ** 2) / mu                       # a cubic: mu u'' is not constant
    assert sp.simplify(impl.ns_residuals(wrong, sp.Integer(0), -G * X, rho, mu)[1]) != 0
    # Taylor-Green vortex: an unsteady exact solution that exercises the nonlinear terms
    nu = mu / rho
    F = sp.exp(-2 * nu * T)
    u = -sp.cos(X) * sp.sin(Y) * F
    v = sp.sin(X) * sp.cos(Y) * F
    p = -rho / 4 * (sp.cos(2 * X) + sp.cos(2 * Y)) * F ** 2
    assert all(sp.simplify(r) == 0 for r in impl.ns_residuals(u, v, p, rho, mu))
    # drop the pressure and the x-momentum residual must not vanish
    assert sp.simplify(impl.ns_residuals(u, v, sp.Integer(0), rho, mu)[1]) != 0


def test_hydrostatic_state_balances_gravity():
    rho, mu, g = sp.symbols("rho mu g", positive=True)
    res = impl.ns_residuals(sp.Integer(0), sp.Integer(0), -rho * g * Y, rho, mu, gy=-g)
    assert all(r == 0 for r in res)


def test_stokes_first_problem_satisfies_navier_stokes_and_conditions():
    U, nu, rho = sp.symbols("U nu rho", positive=True)
    u = U * sp.erfc(Y / (2 * sp.sqrt(nu * T)))
    res = impl.ns_residuals(u, sp.Integer(0), sp.Integer(0), rho, rho * nu)
    assert all(sp.simplify(r) == 0 for r in res)
    y = np.array([0.0, 1e-3, 5e-3, 1.0])
    u_num = impl.stokes_first_problem(y, 2.0, 1.5, 1e-6)
    assert u_num[0] == pytest.approx(1.5) and u_num[-1] < 1e-12             # wall speed; fluid far away at rest
    # self-similar in y / sqrt(nu t): quadrupling t and doubling y gives the same speed
    assert impl.stokes_first_problem(2e-3, 8.0, 1.5, 1e-6) == pytest.approx(float(u_num[1]))
    # the 99 % thickness is 3.64 sqrt(nu t)
    assert impl.stokes_first_problem(3.64 * np.sqrt(1e-6 * 2.0), 2.0, 1.0, 1e-6) == pytest.approx(0.01, abs=2e-4)


def test_pipe_poiseuille_satisfies_axial_momentum_and_its_flow_rate():
    G, mu, R = sp.symbols("G mu R", positive=True)
    w = G * (R ** 2 - impl.R_SYM ** 2) / (4 * mu)
    assert impl.pipe_axial_residual(w, -G * impl.Z_SYM, mu) == 0
    assert impl.pipe_axial_residual(w, -2 * G * impl.Z_SYM, mu) != 0
    Rn, Gn = 2e-3, 5e3
    assert impl.pipe_poiseuille(Rn, Rn, Gn, mu_) == pytest.approx(0.0, abs=1e-15)
    Q = quad(lambda r: 2 * np.pi * r * impl.pipe_poiseuille(r, Rn, Gn, mu_), 0, Rn)[0]
    assert impl.pipe_flow_rate(Rn, Gn, mu_) == pytest.approx(Q, rel=1e-10)
    # mean speed is half the centreline speed
    assert Q / (np.pi * Rn ** 2) == pytest.approx(0.5 * float(impl.pipe_poiseuille(0.0, Rn, Gn, mu_)))


# ---------------------------------------------------------------- flow rate, shear and the momentum balance

def test_plane_flow_rate_is_the_integral_of_the_profile():
    for U, G in ((U_, 0.0), (0.0, G_), (U_, G_), (-U_, G_)):
        q = quad(lambda y: impl.couette_poiseuille(y, h_, U, G, mu_), 0, h_)[0]
        assert impl.plane_flow_rate(h_, U, G, mu_) == pytest.approx(q, rel=1e-10, abs=1e-18)


def test_wall_shear_closes_the_momentum_balance():
    tau0, tauh = impl.wall_shear_stress(h_, U_, G_, mu_)
    # a slab of fluid of length dx: pressure force G h dx is balanced by the wall shear
    assert tauh - tau0 == pytest.approx(-G_ * h_, rel=1e-12)
    # and each equals mu du/dy from the profile (central difference)
    e = 1e-9
    du = lambda y: (impl.couette_poiseuille(y + e, h_, U_, G_, mu_) - impl.couette_poiseuille(y - e, h_, U_, G_, mu_)) / (2 * e)
    assert mu_ * du(0.0) == pytest.approx(tau0, rel=1e-5)
    assert mu_ * du(h_) == pytest.approx(tauh, rel=1e-5)
    # pure Couette flow: the same stress everywhere
    a, b = impl.wall_shear_stress(h_, U_, 0.0, mu_)
    assert a == pytest.approx(b) == pytest.approx(mu_ * U_ / h_)


def test_superposition_and_zero_limits():
    y = np.linspace(0, h_, 11)
    npt.assert_allclose(impl.couette_poiseuille(y, h_, U_, G_, mu_),
                        impl.couette_poiseuille(y, h_, U_, 0.0, mu_) + impl.couette_poiseuille(y, h_, 0.0, G_, mu_),
                        rtol=1e-13, atol=1e-16)
    assert np.all(impl.couette_poiseuille(y, h_, 0.0, 0.0, mu_) == 0.0)
    assert impl.plane_flow_rate(h_, 0.0, 0.0, mu_) == 0.0


def test_annulus_formula_matches_its_profile_and_the_thin_gap_limit():
    Ri, Ro, G = 5e-3, 8e-3, 1e4
    lnk = np.log(Ro / Ri)
    w = lambda r: G / (4 * mu_) * (Ro ** 2 - r ** 2 + (Ro ** 2 - Ri ** 2) * np.log(r / Ro) / lnk)
    assert abs(w(Ri)) < 1e-12 * w(0.5 * (Ri + Ro)) and abs(w(Ro)) < 1e-15
    Q = quad(lambda r: 2 * np.pi * r * w(r), Ri, Ro)[0]
    assert impl.annulus_flow_rate(Ri, Ro, G, mu_) == pytest.approx(Q, rel=1e-9)
    # a thin gap is a plane slot of width pi (Ri + Ro), with an error of order (gap/R)^2
    for gap in (1e-3, 1e-4, 1e-5):
        Ri = 12.5e-3
        slot = np.pi * (2 * Ri + gap) * gap ** 3 * G / (12 * mu_)
        assert impl.annulus_flow_rate(Ri, Ri + gap, G, mu_) / slot - 1 == pytest.approx(0.0, abs=(gap / Ri) ** 2)


# ---------------------------------------------------------------- finite differences

def test_fd_channel_is_exact_for_the_parabola():
    for n in (2, 5, 32):
        y, u = impl.fd_steady_channel(n, h_, U_, G_, mu_)
        assert y.size == n + 1 and y[0] == 0.0 and y[-1] == pytest.approx(h_)
        ex = impl.couette_poiseuille(y, h_, U_, G_, mu_)
        assert np.max(np.abs(u - ex)) < 1e-12 * np.max(np.abs(ex))


def test_fd_channel_converges_at_second_order_on_a_manufactured_problem():
    # u = sin(pi y / h) solves mu u'' = -f with f = mu (pi/h)^2 sin(pi y / h): u'''' != 0
    f = lambda y: mu_ * (np.pi / h_) ** 2 * np.sin(np.pi * y / h_)
    errs = []
    for n in (10, 20, 40, 80):
        y, u = impl.fd_steady_channel(n, h_, 0.0, f, mu_)
        errs.append(np.max(np.abs(u - np.sin(np.pi * y / h_))))
    orders = impl.observed_order(errs)
    assert np.all(np.abs(orders - 2.0) < 0.05)


def test_fd_pipe_is_exact_and_conserves_the_flow_rate():
    R, G = 3e-3, 4e3
    r, w = impl.fd_steady_pipe(16, R, G, mu_)
    ex = impl.pipe_poiseuille(r, R, G, mu_)
    assert np.max(np.abs(w - ex)) < 1e-11 * ex[0]
    assert w[-1] == 0.0
    r, w = impl.fd_steady_pipe(400, R, G, mu_)
    Q = np.trapezoid(2 * np.pi * r * w, r)
    assert Q == pytest.approx(impl.pipe_flow_rate(R, G, mu_), rel=1e-4)


def test_crank_nicolson_converges_at_second_order_on_stokes_first_problem():
    nu, t0, t1, Yd = 1.0, 0.01, 0.05, 2.0
    errs = []
    for n in (20, 40, 80):
        y = np.linspace(0, Yd, n + 1)
        u = impl.crank_nicolson(y, impl.stokes_first_problem(y, t0, 1.0, nu), t0, t1, n // 2, nu,
                                u_lo=1.0, u_hi=lambda t: float(impl.stokes_first_problem(Yd, t, 1.0, nu)))
        errs.append(np.max(np.abs(u - impl.stokes_first_problem(y, t1, 1.0, nu))))
    orders = impl.observed_order(errs)
    assert np.all(np.abs(orders - 2.0) < 0.1)
    assert errs[-1] < 1e-3


def test_startup_flow_series_and_finite_differences_agree_and_relax_to_poiseuille():
    h, G, rho, mu = 1.0, 1.0, 1.0, 1.0
    y = np.linspace(0, h, 41)
    steady = impl.couette_poiseuille(y, h, 0.0, G, mu)
    assert np.max(np.abs(impl.startup_channel_series(y, 0.0, h, G, rho, mu))) < 1e-5 * steady.max()
    npt.assert_allclose(impl.startup_channel_series(y, 5.0, h, G, rho, mu), steady, atol=1e-15)
    errs = []
    for n in (20, 40, 80):
        yy = np.linspace(0, h, n + 1)
        u = impl.crank_nicolson(yy, np.zeros_like(yy), 0.0, 0.05, n, mu / rho, source=G / rho)
        errs.append(np.max(np.abs(u - impl.startup_channel_series(yy, 0.05, h, G, rho, mu))))
    assert np.all(impl.observed_order(errs) > 1.8)
    u = impl.crank_nicolson(y, np.zeros_like(y), 0.0, 3.0, 300, mu / rho, source=G / rho)
    npt.assert_allclose(u, steady, atol=1e-6 * steady.max())


def test_startup_time_reaches_99_percent_of_the_centreline_speed():
    h, nu = 2e-5, 1.4e-5
    t99 = impl.startup_time(h, nu)
    uc = impl.startup_channel_series(h / 2, t99, h, 1.0, 1.0, nu)
    assert float(uc) / (h ** 2 / (8 * nu)) == pytest.approx(0.99, abs=1e-6)
    assert impl.startup_time(2 * h, nu) == pytest.approx(4 * t99)          # diffusion time ~ h^2/nu


# ---------------------------------------------------------------- the actuator

def test_actuator_leakage_scaling_and_force_balance():
    D, h, L, dp, mu, rho = 25e-3, 10e-6, 12e-3, 20.7e6, impl.MU_HYD, impl.RHO_HYD
    a = impl.actuator_leakage(D, h, L, dp, 0.0, mu, rho)
    assert a["Q"] == pytest.approx(np.pi * D * h ** 3 * dp / (12 * mu * L))
    assert impl.actuator_leakage(D, 2 * h, L, dp, 0.0, mu, rho)["Q"] == pytest.approx(8 * a["Q"])  # h^3
    assert a["Re_gap"] < 10                                                  # laminar by a wide margin
    # with the piston at rest the pressure force on the gap fluid is shared equally by the two walls
    assert 2 * a["force"] == pytest.approx(dp * np.pi * D * h, rel=1e-12)
    b = impl.actuator_leakage(D, h, L, dp, 0.1, mu, rho)
    assert b["force"] == pytest.approx(b["force_pressure"] + b["force_couette"])
    assert b["Q"] > a["Q"] and b["force"] < a["force"]                      # moving with the leak
    assert a["t99"] < 1e-4                                                   # steady in microseconds


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "navier_stokes.png").stat().st_size > 10_000
