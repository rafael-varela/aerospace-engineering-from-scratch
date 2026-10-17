"""Public tests for P06.12 Hamiltonian Mechanics. They assert mechanics, not magic numbers: the Legendre
transform of the pendulum and of Kepler in polar coordinates; Hamilton's equations reproduce Lagrange's;
canonical Poisson brackets and the conserved angular momentum and Laplace-Runge-Lenz vector; the
pendulum's elliptic-integral period equals the derivative of its phase area (two independent routes) and
a simulated period; Liouville's theorem for the exact flow, exponential area contraction with damping;
the area factor of each one-step map; Kepler's invariants, the radial action and its energy derivative,
the symplectic flow Jacobian; the along-track drift of a semi-major-axis error; bounded Verlet energy
error against secular RK4 drift.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
import sympy as sp
from scipy.integrate import solve_ivp

from conftest import load_impl

impl = load_impl(__file__)

W0 = np.pi                       # a 2 s small-angle pendulum
MU = impl.MU_EARTH


# ---------------------------------------------------------------- symbolic

def test_legendre_transform_of_the_pendulum_and_hamilton_equations():
    th, thd, p = sp.symbols("theta thetadot p")
    m, l, g = sp.symbols("m l g", positive=True)
    L = m * l ** 2 * thd ** 2 / 2 + m * g * l * sp.cos(th)
    H, mom = impl.legendre_transform(L, [th], [thd], [p])
    assert sp.simplify(H - (p ** 2 / (2 * m * l ** 2) - m * g * l * sp.cos(th))) == 0
    assert sp.simplify(mom[0] - m * l ** 2 * thd) == 0
    (qd,), (pd,) = impl.hamilton_equations(H, [th], [p])
    # eliminating p gives Lagrange's equation m l^2 th'' = -m g l sin th
    assert sp.simplify(qd - p / (m * l ** 2)) == 0 and sp.simplify(pd + m * g * l * sp.sin(th)) == 0


def test_legendre_transform_of_kepler_in_polar_coordinates():
    r, th, rd, thd, pr, pt, mu = sp.symbols("r theta rdot thetadot p_r p_theta mu", positive=True)
    H, _ = impl.legendre_transform((rd ** 2 + r ** 2 * thd ** 2) / 2 + mu / r, [r, th], [rd, thd], [pr, pt])
    assert sp.simplify(H - (pr ** 2 / 2 + pt ** 2 / (2 * r ** 2) - mu / r)) == 0
    assert impl.poisson_bracket(pt, H, [r, th], [pr, pt]) == 0               # theta cyclic: p_theta conserved
    assert sp.simplify(impl.poisson_bracket(pr, H, [r, th], [pr, pt]) - (pt ** 2 / r ** 3 - mu / r ** 2)) == 0


def test_poisson_brackets_canonical_and_kepler_invariants():
    x, y, px, py, mu = sp.symbols("x y p_x p_y mu", positive=True)
    qs, ps = [x, y], [px, py]
    assert impl.poisson_bracket(x, px, qs, ps) == 1 and impl.poisson_bracket(x, py, qs, ps) == 0
    assert impl.poisson_bracket(px, x, qs, ps) == -1
    H = (px ** 2 + py ** 2) / 2 - mu / sp.sqrt(x ** 2 + y ** 2)
    h = x * py - y * px
    assert impl.poisson_bracket(h, H, qs, ps) == 0
    Ax = py * h - mu * x / sp.sqrt(x ** 2 + y ** 2)
    assert sp.simplify(impl.poisson_bracket(Ax, H, qs, ps)) == 0                # LRL is conserved


# ---------------------------------------------------------------- pendulum

def test_separatrix_and_hamiltonian():
    q = np.linspace(-np.pi, np.pi, 41)
    npt.assert_allclose(impl.pendulum_hamiltonian(q, impl.separatrix_momentum(q, W0), W0), W0 ** 2, rtol=1e-12)
    assert impl.pendulum_hamiltonian(0.0, 0.0, W0) == pytest.approx(-W0 ** 2)


@pytest.mark.parametrize("E_over", [-0.999, -0.5, 0.3, 0.9, 1.5, 3.0])
def test_period_equals_derivative_of_phase_area(E_over):
    E = E_over * W0 ** 2
    dE = 1e-6 * W0 ** 2
    TA = (impl.phase_area(E + dE, W0) - impl.phase_area(E - dE, W0)) / (2 * dE)
    assert TA == pytest.approx(impl.pendulum_period(E, W0), rel=1e-6)


def test_period_limits_and_simulation():
    assert impl.pendulum_period(-W0 ** 2 * np.cos(1e-4), W0) == pytest.approx(2 * np.pi / W0, rel=1e-8)
    E = -W0 ** 2 * np.cos(np.radians(120.0))
    assert impl.pendulum_period(E, W0) > impl.pendulum_period(-0.5 * W0 ** 2, W0)       # slower when larger
    # simulated: time between successive upward zero crossings of q
    sol = solve_ivp(impl.pendulum_flow, (0, 20), [np.radians(120.0), 0.0], args=(W0,), rtol=1e-11, atol=1e-12,
                    events=lambda t, y, w: y[0], dense_output=True)
    ev = sol.t_events[0]
    assert np.mean(np.diff(ev[::2])) == pytest.approx(impl.pendulum_period(E, W0), rel=1e-8)
    # small-oscillation phase area is an ellipse: A = 2 pi (E - E_min)/w0
    Es = -W0 ** 2 + 1e-4
    assert impl.phase_area(Es, W0) == pytest.approx(2 * np.pi * 1e-4 / W0, rel=1e-4)
    with pytest.raises(ValueError):
        impl.pendulum_period(-2 * W0 ** 2, W0)


# ---------------------------------------------------------------- Liouville

def _blob(n=2000):
    return impl.blob_boundary(0.0, 1.2 * W0, 0.3, 0.3 * W0, n)


def test_liouville_exact_flow_keeps_area_damping_shrinks_it_exponentially():
    qb, pb = _blob()
    A0 = impl.polygon_area(qb, pb)
    assert A0 == pytest.approx(np.pi * 0.3 * 0.3 * W0, rel=1e-5)
    qe, pe = impl.flow_blob(qb, pb, 6.0, W0)
    assert impl.polygon_area(qe, pe) == pytest.approx(A0, rel=2e-5)     # polygon chord error ~5e-6 at n = 2000
    assert np.max(np.abs(qe - qb)) > 0.5                                # it did move and shear
    c = 0.3
    qd, pd = impl.flow_blob(qb, pb, 6.0, W0, c=c)
    # div f = -c everywhere, so the area decays exactly as exp(-c t)
    assert impl.polygon_area(qd, pd) == pytest.approx(A0 * np.exp(-c * 6.0), rel=2e-5)


def test_step_area_factors():
    acc = lambda q: -W0 ** 2 * np.sin(q)
    rng = np.random.default_rng(2)
    h = 0.05
    for q, p in rng.uniform(-2, 2, (5, 2)):
        assert impl.step_jacobian_det(impl.verlet_step, acc, q, p, h) == pytest.approx(1.0, abs=1e-9)
        # explicit Euler: J = [[1, h], [h a'(q), 1]], det = 1 + h^2 w0^2 cos q
        assert impl.step_jacobian_det(impl.euler_step, acc, q, p, h) == pytest.approx(1 + h * h * W0 ** 2 * np.cos(q), abs=1e-9)
    d = impl.step_jacobian_det(impl.rk4_step, acc, 0.0, 0.0, h)
    assert 0 < 1 - d < 1e-5                                             # RK4 contracts, slightly


def test_blob_under_one_step_maps():
    acc = lambda q: -W0 ** 2 * np.sin(q)
    qb, pb = _blob()
    A0 = impl.polygon_area(qb, pb)
    qe, pe = impl.flow_blob(qb, pb, 6.0, W0)
    Ae = impl.polygon_area(qe, pe)                                      # same polygon error as the maps
    qv, pv = impl.iterate(impl.verlet_step, acc, qb, pb, 0.05, 120)
    qr, pr = impl.iterate(impl.rk4_step, acc, qb, pb, 0.05, 120)
    qx, px = impl.iterate(impl.euler_step, acc, qb, pb, 0.05, 20)
    Av, Ar, Ax = impl.polygon_area(qv, pv), impl.polygon_area(qr, pr), impl.polygon_area(qx, px)
    assert abs(Av - Ae) < 1e-6 * A0                                     # symplectic: area of the exact flow
    assert abs(Ar - Ae) > 5 * abs(Av - Ae) and Ar < Ae                  # RK4 loses area
    assert Ax > 1.15 * A0                                               # Euler gains area: (1 + h^2 w0^2 cos q)^n


# ---------------------------------------------------------------- Kepler

def test_kepler_invariants_and_lrl_vector():
    a, e = 7.0e6, 0.1
    y0 = impl.periapsis_state(a, e)
    T = 2 * np.pi * np.sqrt(a ** 3 / MU)
    sol = solve_ivp(impl.kepler_flow, (0, 2 * T), y0, method="DOP853", rtol=1e-12, atol=1e-6,
                    t_eval=np.linspace(0, 2 * T, 401))
    Y = sol.y.T
    E = impl.kepler_hamiltonian(Y[:, :2], Y[:, 2:])
    assert E[0] == pytest.approx(-MU / (2 * a), rel=1e-12)
    assert np.ptp(E) < 1e-10 * abs(E[0])
    A = impl.lrl_vector(Y[:, :2], Y[:, 2:])
    npt.assert_allclose(np.linalg.norm(A, axis=1), MU * e, rtol=1e-8)     # |A| = mu e
    assert np.max(np.abs(A[:, 1])) < 1e-8 * MU                          # points at periapsis (+x)
    # polar and Cartesian Hamiltonians agree for the same states
    r = np.linalg.norm(Y[:, :2], axis=1)
    prr = np.sum(Y[:, :2] * Y[:, 2:], axis=1) / r
    pt = Y[0, 0] * Y[0, 3] - Y[0, 1] * Y[0, 2]
    npt.assert_allclose(impl.kepler_polar_hamiltonian(r, prr, pt), E, rtol=1e-9)


def test_effective_potential_minimum_is_the_circular_orbit():
    h = np.sqrt(MU * 7.0e6)
    r = np.linspace(5e6, 9e6, 40001)
    V = impl.effective_potential(r, h)
    assert r[np.argmin(V)] == pytest.approx(h ** 2 / MU, rel=1e-4)
    assert V.min() == pytest.approx(-MU ** 2 / (2 * h ** 2), rel=1e-9)


def test_radial_action_and_its_energy_derivative_is_the_period():
    a, e = 7.0e6, 0.3
    y0 = impl.periapsis_state(a, e)
    E = float(impl.kepler_hamiltonian(y0[:2], y0[2:]))
    h = y0[0] * y0[3]
    Jr = impl.radial_action(E, h)
    assert Jr == pytest.approx(2 * np.pi * (MU / np.sqrt(-2 * E) - h), rel=1e-6)
    dE = 1e-6 * abs(E)
    T = (impl.radial_action(E + dE, h) - impl.radial_action(E - dE, h)) / (2 * dE)
    assert T == pytest.approx(2 * np.pi * np.sqrt(a ** 3 / MU), rel=1e-5)


def test_kepler_flow_jacobian_is_symplectic():
    y0 = impl.periapsis_state(1.0, 0.3, 1.0)                          # canonical units, mu = 1
    M = impl.flow_jacobian(y0, 2 * np.pi * 1.3, 1.0)
    assert np.linalg.det(M) == pytest.approx(1.0, abs=1e-6)
    assert impl.symplectic_residual(M) < 1e-5
    assert np.max(np.abs(M)) > 5                                         # not trivially the identity: it shears


def test_along_track_drift_of_a_semi_major_axis_error():
    a, e, da = 7.0e6, 0.1, 10.0
    T = 2 * np.pi * np.sqrt(a ** 3 / MU)
    kw = dict(method="DOP853", rtol=1e-12, atol=1e-6)
    n = 10
    s0 = solve_ivp(impl.kepler_flow, (0, n * T), impl.periapsis_state(a, e), **kw).y[:, -1]
    s1 = solve_ivp(impl.kepler_flow, (0, n * T), impl.periapsis_state(a + da, e), **kw).y[:, -1]
    d, v = s1[:2] - s0[:2], s0[2:]
    along = abs(d @ v) / np.linalg.norm(v)
    assert along == pytest.approx(impl.along_track_drift(a, da, n, e), rel=1e-3)
    assert impl.along_track_drift(a, da, 2 * n, 0.0) == pytest.approx(2 * impl.along_track_drift(a, da, n, 0.0))


def test_verlet_energy_error_bounded_rk4_drifts():
    y0 = impl.periapsis_state(7.0e6, 0.1)
    _, ev = impl.energy_error_history(impl.verlet_step, y0, 200, 40)
    _, er = impl.energy_error_history(impl.rk4_step, y0, 200, 40)
    k = 200 * 5
    assert np.max(np.abs(ev[-k:])) == pytest.approx(np.max(np.abs(ev[:k])), rel=0.05)   # bounded swing
    assert np.max(np.abs(er[-k:])) > 4 * np.max(np.abs(er[:k]))                          # secular drift
    npt.assert_allclose(impl.kepler_hamiltonian(y0[:2], y0[2:]), -MU / (2 * 7.0e6), rtol=1e-12)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "phase_space.png").stat().st_size > 10_000
