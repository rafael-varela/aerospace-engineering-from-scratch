"""Public tests for P06.11 Lagrangian Mechanics. They assert mechanics, not magic numbers: Lagrange's
equations reproduce the spring-mass and the simple pendulum; the SymPy double pendulum has the hand
mass matrix and forcing; the Lagrange accelerations equal an independent Newton-Euler solution at random
states; the link tensions have the right static and rigid-rotation limits; the integrated motion
conserves energy and is time-reversible; the small-angle modes solve the closed-form frequency equation,
are M-orthogonal, tend to two independent pendulums when the lower mass vanishes, and are what a
small-amplitude simulation actually does; a large-amplitude release is chaotic and a small one is not.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
import sympy as sp

from conftest import load_impl

impl = load_impl(__file__)

P_LAB = impl.DoublePendulum(1.0, 1.0, 1.0, 1.0)
P_ODD = impl.DoublePendulum(m1=2.5, m2=0.7, l1=0.8, l2=1.9, g=9.81)      # unequal everything
P_SLUNG = impl.DoublePendulum(m1=80.0, m2=1500.0, l1=12.0, l2=4.0)


def _random_states(n, seed):
    rng = np.random.default_rng(seed)
    return [np.concatenate([rng.uniform(-np.pi, np.pi, 2), rng.uniform(-4.0, 4.0, 2)]) for _ in range(n)]


# ---------------------------------------------------------------- Lagrange's equations, symbolically

def test_euler_lagrange_spring_mass_and_simple_pendulum():
    t = sp.Symbol("t")
    m, k, g, l = sp.symbols("m k g l", positive=True)
    x = sp.Function("x")(t)
    (eq,) = impl.euler_lagrange(m * sp.diff(x, t) ** 2 / 2 - k * x ** 2 / 2, [x], t)
    assert sp.simplify(eq - (m * sp.diff(x, t, 2) + k * x)) == 0
    th = sp.Function("th")(t)
    L = m * l ** 2 * sp.diff(th, t) ** 2 / 2 + m * g * l * sp.cos(th)
    (eq,) = impl.euler_lagrange(L, [th], t)
    assert sp.simplify(eq - (m * l ** 2 * sp.diff(th, t, 2) + m * g * l * sp.sin(th))) == 0


def test_sympy_mass_matrix_and_forcing_match_the_hand_forms():
    M, f, syms = impl.derive_equations()
    Mf = sp.lambdify(syms, M, "numpy")
    ff = sp.lambdify(syms, f, "numpy")
    for p in (P_LAB, P_ODD):
        for y in _random_states(10, 3):
            args = (*y, p.m1, p.m2, p.l1, p.l2, p.g)
            npt.assert_allclose(np.array(Mf(*args), dtype=float), impl.mass_matrix(y[:2], p), rtol=1e-12)
            npt.assert_allclose(np.array(ff(*args), dtype=float).ravel(), impl.generalized_forces(y, p),
                                rtol=1e-12, atol=1e-12)


def test_lower_mass_zero_reduces_to_a_simple_pendulum():
    M, f, syms = impl.derive_equations()
    a1, a2, w1, w2, m1, m2, l1, l2, g = syms
    # with m2 = 0 the first equation decouples and is the simple pendulum th1'' = -(g/l1) sin th1
    assert sp.simplify(M[0, 1].subs(m2, 0)) == 0
    assert sp.simplify(f[0].subs(m2, 0) / M[0, 0].subs(m2, 0) + g / l1 * sp.sin(a1)) == 0
    # and the SymPy accelerations at a tiny m2 agree numerically
    (d1, _), _ = impl.derive_accelerations()
    val = float(d1.subs({a1: 0.7, a2: -0.4, w1: 0.3, w2: 1.1, m1: 2.0, m2: 1e-12, l1: 1.5, l2: 0.5, g: 9.81}))
    assert val == pytest.approx(-9.81 / 1.5 * np.sin(0.7), rel=1e-9)


def test_mass_matrix_symmetric_positive_definite():
    for y in _random_states(20, 4):
        M = impl.mass_matrix(y[:2], P_ODD)
        npt.assert_allclose(M, M.T, atol=0)
        assert np.all(np.linalg.eigvalsh(M) > 0)


# ---------------------------------------------------------------- the independent check: Newton-Euler

@pytest.mark.parametrize("p", [P_LAB, P_ODD, P_SLUNG])
def test_lagrange_equals_newton_euler_at_random_states(p):
    for y in _random_states(50, 11):
        aL = impl.lagrange_accelerations(y, p)
        aN, _ = impl.newton_euler(y, p)
        npt.assert_allclose(aL, aN, rtol=1e-10, atol=1e-10 * (1 + np.max(np.abs(aN))))
        npt.assert_allclose(impl.rhs(0.0, y, p)[2:], aN, rtol=1e-10, atol=1e-10 * (1 + np.max(np.abs(aN))))
        npt.assert_allclose(impl.rhs(0.0, y, p)[:2], y[2:])


def test_tensions_static_and_rigid_rotation_limits():
    p = P_ODD
    acc, (T1, T2) = impl.newton_euler(np.zeros(4), p)
    npt.assert_allclose(acc, 0.0, atol=1e-12)
    assert T2 == pytest.approx(p.m2 * p.g) and T1 == pytest.approx((p.m1 + p.m2) * p.g)
    w = 1.7                                          # both links vertical, turning together at w
    acc, (T1, T2) = impl.newton_euler(np.array([0.0, 0.0, w, w]), p)
    npt.assert_allclose(acc, 0.0, atol=1e-12)
    assert T2 == pytest.approx(p.m2 * (p.g + (p.l1 + p.l2) * w ** 2))
    assert T1 == pytest.approx(T2 + p.m1 * (p.g + p.l1 * w ** 2))


# ---------------------------------------------------------------- integration and energy

def test_energy_is_conserved_and_matches_cartesian_kinetic_energy():
    t = np.linspace(0.0, 10.0, 5001)
    y0 = np.array([2.0, -1.0, 0.5, 1.5])                    # far from small angles
    Y = impl.simulate(y0, t, P_ODD, rtol=1e-11, atol=1e-11)
    E = impl.energy(Y, P_ODD)
    assert np.max(np.abs(E - E[0])) < 1e-8 * abs(E[0]) + 1e-8      # DOP853 at rtol 1e-11
    # energy from the Cartesian positions (velocities by central differences): an independent formula
    x1, y1, x2, y2 = impl.positions(Y, P_ODD)
    dt = t[1] - t[0]
    v = lambda s: (s[:-4] - 8 * s[1:-3] + 8 * s[3:-1] - s[4:]) / (12 * dt)   # fourth-order differences
    p = P_ODD
    T = 0.5 * p.m1 * (v(x1) ** 2 + v(y1) ** 2) + 0.5 * p.m2 * (v(x2) ** 2 + v(y2) ** 2)
    V = p.g * (p.m1 * y1 + p.m2 * y2)[2:-2]
    npt.assert_allclose(T + V, E[2:-2], atol=1e-5 * np.max(np.abs(E)))     # O(dt^4) differences


def test_time_reversibility():
    t = np.linspace(0.0, 3.0, 301)
    y0 = np.array([1.2, 0.4, -0.3, 0.8])
    yf = impl.simulate(y0, t, P_LAB)[:, -1]
    yb = impl.simulate(yf * np.array([1, 1, -1, -1]), t, P_LAB)[:, -1]
    npt.assert_allclose(yb * np.array([1, 1, -1, -1]), y0, atol=1e-7)


# ---------------------------------------------------------------- small-angle normal modes

@pytest.mark.parametrize("p", [P_LAB, P_ODD, P_SLUNG])
def test_normal_modes_solve_the_frequency_equation_and_are_orthogonal(p):
    w, phi = impl.normal_modes(p)
    npt.assert_allclose(w, impl.normal_mode_frequencies_closed_form(p), rtol=1e-10)
    assert w[0] < w[1]
    M0, K0 = impl.linearized_matrices(p)
    for i in range(2):
        assert abs(np.linalg.det(K0 - w[i] ** 2 * M0)) < 1e-8 * np.linalg.norm(K0) ** 2
        npt.assert_allclose(K0 @ phi[:, i], w[i] ** 2 * M0 @ phi[:, i], rtol=1e-10)
    assert abs(phi[:, 0] @ M0 @ phi[:, 1]) < 1e-10 * np.linalg.norm(M0) * np.linalg.norm(phi) ** 2
    assert phi[1, 0] > 0 > phi[1, 1]                 # in phase, then out of phase


def test_linearisation_of_the_sympy_equations_is_minus_Minv_K():
    p = P_ODD
    h = 1e-6
    J = np.zeros((2, 2))
    for j in range(2):
        e = np.zeros(4)
        e[j] = h
        J[:, j] = (impl.lagrange_accelerations(e, p) - impl.lagrange_accelerations(-e, p)) / (2 * h)
    M0, K0 = impl.linearized_matrices(p)
    npt.assert_allclose(J, -np.linalg.solve(M0, K0), rtol=1e-7)


def test_vanishing_lower_mass_gives_two_independent_pendulums():
    p = impl.DoublePendulum(m1=1.0, m2=1e-9, l1=1.0, l2=0.25)
    npt.assert_allclose(impl.normal_mode_frequencies_closed_form(p), np.sqrt(p.g / np.array([1.0, 0.25])), rtol=1e-6)


@pytest.mark.parametrize("mode", [0, 1])
def test_small_amplitude_motion_is_the_linear_mode(mode):
    p = P_SLUNG
    w, phi = impl.normal_modes(p)
    T = 2 * np.pi / w[mode]
    t = np.linspace(0.0, 20 * T, 4001)
    y0 = np.concatenate([1e-3 * phi[:, mode], [0.0, 0.0]])
    Y = impl.simulate(y0, t, p)
    assert impl.dominant_frequency(t, Y[0]) == pytest.approx(w[mode], rel=2e-3)   # FFT bin over 20 periods
    # the shape is held: th2 = (phi2/phi1) th1 throughout
    assert np.max(np.abs(Y[1] - phi[1, mode] * Y[0])) < 1e-4 * np.max(np.abs(Y[1]))


def test_large_swing_is_slower_than_the_linear_mode():
    p = P_SLUNG
    w, phi = impl.normal_modes(p)
    t = np.linspace(0.0, 10 * 2 * np.pi / w[0], 4001)
    a = np.radians(20.0)
    Y = impl.simulate(np.concatenate([a * phi[:, 0], [0.0, 0.0]]), t, p)
    ws = impl.dominant_frequency(t, Y[0])
    # single-pendulum amplitude correction 1 - a^2/16 as the estimate (the mode is almost rigid)
    assert ws < w[0] and ws == pytest.approx(w[0] * (1 - a ** 2 / 16), rel=2e-3)


# ---------------------------------------------------------------- chaos

def test_large_release_is_chaotic_small_release_is_not():
    lam_big, tk, run = impl.lyapunov_estimate(np.radians([120.0, 120.0, 0.0, 0.0]), P_LAB, 20.0)
    lam_small, _, _ = impl.lyapunov_estimate(np.radians([5.0, 5.0, 0.0, 0.0]), P_LAB, 20.0)
    assert len(tk) == len(run) == 40
    assert lam_big > 0.5            # separation grows by more than e^10 in 20 s
    assert lam_small < 0.1          # regular motion: separation grows at most algebraically
    assert lam_big > 10 * lam_small


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "double_pendulum.png").stat().st_size > 10_000
