"""Tests for P02.05 State transition matrix (matrix exponential).

Run against the starter (default) or the reference solution (AEFS_IMPL=solution).
"""
from __future__ import annotations

import numpy as np
import numpy.testing as npt
import pytest
from conftest import load_impl
from scipy.integrate import solve_ivp
from scipy.linalg import expm

impl = load_impl(__file__)

RNG = np.random.default_rng(20250205)
SP = np.array([[-0.8, 1.0], [-4.0, -1.2]])          # the short-period matrix of P02.04
OSC = np.array([[0.0, 1.0], [-4.0, -0.4]])          # a damped oscillator, omega_n = 2


def random_matrix(n: int, scale: float = 1.0) -> np.ndarray:
    return scale * RNG.standard_normal((n, n))


# ---------------------------------------------------------------- the exponential itself

def test_taylor_series_is_exact_for_nilpotent_and_matches_the_scalar_exponential():
    N = np.array([[0.0, 1.0], [0.0, 0.0]])
    npt.assert_allclose(impl.expm_taylor(N * 2.5, 5), np.array([[1.0, 2.5], [0.0, 1.0]]), atol=1e-15)
    npt.assert_allclose(impl.expm_taylor(np.array([[0.3]]), 30), np.array([[np.exp(0.3)]]), rtol=1e-14)
    npt.assert_allclose(impl.expm_taylor(np.zeros((3, 3))), np.eye(3), atol=0.0)


def test_taylor_series_converges_for_small_norm_but_fails_for_large():
    A = random_matrix(3, 0.3)
    npt.assert_allclose(impl.expm_taylor(A, 30), expm(A), rtol=1e-12, atol=1e-13)
    big = SP * 30.0                                     # ||A|| ~ 150: 30 terms are hopeless
    assert np.abs(impl.expm_taylor(big, 30) - expm(big)).max() > 1.0


@pytest.mark.parametrize("n,scale", [(2, 1.0), (3, 5.0), (4, 40.0), (6, 2.0)])
def test_scaling_and_squaring_matches_scipy_on_random_matrices(n, scale):
    A = random_matrix(n, scale)
    E = impl.expm_scaling_squaring(A)
    ref = expm(A)
    npt.assert_allclose(E, ref, rtol=1e-9, atol=1e-9 * np.abs(ref).max())


def test_eigen_method_matches_and_refuses_defective_matrices():
    A = random_matrix(3, 2.0)
    npt.assert_allclose(impl.expm_eig(A), expm(A), rtol=1e-8, atol=1e-10)
    npt.assert_allclose(impl.expm_eig(SP), expm(SP), rtol=1e-10)          # complex eigenvalues, real result
    assert np.isrealobj(impl.expm_eig(SP))
    with pytest.raises(ValueError):
        impl.expm_eig(impl.cw_matrix(1e-3))
    with pytest.raises(ValueError):
        impl.expm_eig(np.array([[1.0, 1.0], [0.0, 1.0]]))


def test_commuting_matrices_split_the_exponential_and_non_commuting_ones_do_not():
    A = SP
    B = 2.0 * SP + 0.5 * np.eye(2)                     # commutes with A
    npt.assert_allclose(impl.expm_scaling_squaring(A + B), impl.expm_scaling_squaring(A) @ impl.expm_scaling_squaring(B), rtol=1e-9)
    C = np.array([[0.0, 1.0], [0.0, 0.0]])
    D = np.array([[0.0, 0.0], [1.0, 0.0]])
    assert np.abs(impl.expm_scaling_squaring(C + D) - impl.expm_scaling_squaring(C) @ impl.expm_scaling_squaring(D)).max() > 0.1


# ---------------------------------------------------------------- state transition matrix properties

@pytest.mark.parametrize("A", [SP, OSC, random_matrix(3, 1.5), random_matrix(4, 0.7)])
def test_state_transition_group_properties(A):
    I = np.eye(A.shape[0])
    npt.assert_allclose(impl.state_transition(A, 0.0), I, atol=1e-15)
    P1, P2, P3 = impl.state_transition(A, 0.7), impl.state_transition(A, 1.9), impl.state_transition(A, 2.6)
    npt.assert_allclose(P1 @ P2, P3, rtol=1e-9, atol=1e-12)                        # Phi(t+s) = Phi(t) Phi(s)
    npt.assert_allclose(P2 @ P1, P3, rtol=1e-9, atol=1e-12)                        # ... and they commute
    npt.assert_allclose(impl.state_transition(A, -0.7) @ P1, I, atol=1e-10)         # Phi(-t) = Phi(t)^-1
    npt.assert_allclose(np.linalg.det(P1), np.exp(np.trace(A) * 0.7), rtol=1e-10)  # Liouville


@pytest.mark.parametrize("A", [SP, OSC, random_matrix(3, 1.0)])
def test_state_transition_satisfies_the_matrix_differential_equation(A):
    t, h = 1.3, 1e-5
    dPhi = (impl.state_transition(A, t + h) - impl.state_transition(A, t - h)) / (2 * h)
    npt.assert_allclose(dPhi, A @ impl.state_transition(A, t), rtol=1e-7, atol=1e-8)


def test_propagate_matches_solve_ivp_and_the_closed_form_of_P02_04():
    x0 = np.array([np.radians(5.0), 0.0])
    t = np.linspace(0.0, 8.0, 81)
    x = impl.propagate(SP, x0, t)
    assert x.shape == (81, 2)
    sol = solve_ivp(lambda tt, y: SP @ y, (0.0, 8.0), x0, t_eval=t, rtol=1e-11, atol=1e-13)
    npt.assert_allclose(x, sol.y.T, atol=1e-9)
    # eigen-decomposition closed form, as in P02.04
    lam, V, Vinv = impl.modal_decomposition(SP)
    coef = Vinv @ x0
    xc = np.real((np.exp(np.outer(t, lam)) * coef) @ V.T)
    npt.assert_allclose(x, xc, atol=1e-10)
    npt.assert_allclose(V @ np.diag(lam) @ Vinv, SP, atol=1e-12)


def test_stability_follows_the_eigenvalues():
    assert impl.is_stable(SP)
    assert impl.is_stable(OSC)
    assert not impl.is_stable(np.array([[-0.8, 1.0], [2.0, -1.2]]))       # the saddle of P02.04
    assert not impl.is_stable(impl.cw_matrix(1e-3))                       # marginal: eigenvalues on the axis
    assert not impl.is_stable(np.array([[0.0, 1.0], [-4.0, 0.0]]))        # undamped: marginal, not asymptotically stable


# ---------------------------------------------------------------- forced response and discretisation

def test_discretize_reproduces_the_exact_step_and_response_integral():
    B = np.array([[0.0], [1.0]])
    dt = 0.05
    Ad, Bd = impl.discretize(OSC, B, dt)
    npt.assert_allclose(Ad, expm(OSC * dt), rtol=1e-10)
    # B_d = int_0^dt e^{A s} ds B by fine quadrature
    s = np.linspace(0.0, dt, 2001)
    integrand = np.array([expm(OSC * si) @ B for si in s])
    Bd_quad = np.trapezoid(integrand, s, axis=0)
    npt.assert_allclose(Bd, Bd_quad, rtol=1e-7)
    # small-dt limit: B_d -> dt B
    Ad2, Bd2 = impl.discretize(OSC, B, 1e-5)
    npt.assert_allclose(Bd2, 1e-5 * B, rtol=1e-4, atol=1e-10)
    npt.assert_allclose(Ad2, np.eye(2) + 1e-5 * OSC, atol=1e-9)


def test_zero_order_hold_response_is_exact_for_piecewise_constant_input():
    B = np.array([[0.0], [1.0]])
    dt = 0.1
    u = RNG.standard_normal(40)                        # a random piecewise-constant force
    x0 = np.array([0.1, 0.0])
    xs = impl.forced_response_zoh(OSC, B, x0, u, dt)
    assert xs.shape == (41, 2)
    npt.assert_allclose(xs[0], x0)
    # integrate the same piecewise-constant input tightly with solve_ivp, one hold at a time
    y = x0.copy()
    for k in range(40):
        sol = solve_ivp(lambda tt, yy: OSC @ yy + (B * u[k]).ravel(), (0.0, dt), y, rtol=1e-12, atol=1e-14)
        y = sol.y[:, -1]
        npt.assert_allclose(xs[k + 1], y, atol=1e-9)


def test_constant_input_reaches_the_static_deflection():
    B = np.array([[0.0], [1.0]])
    u = np.full(10000, 2.0)                            # constant force 2 on x'' + 0.4 x' + 4 x = F, 100 s (transient ~ e^-20)
    xs = impl.forced_response_zoh(OSC, B, np.zeros(2), u, 0.01)
    npt.assert_allclose(xs[-1], [2.0 / 4.0, 0.0], atol=1e-6)


# ---------------------------------------------------------------- Clohessy-Wiltshire

def test_cw_matrix_closed_form_stm_agrees_with_our_exponential_and_scipy():
    n = impl.iss_orbit()["n"]
    A = impl.cw_matrix(n)
    assert n == pytest.approx(np.sqrt(impl.MU_EARTH / 6778e3 ** 3))
    lam = np.linalg.eigvals(A)
    npt.assert_allclose(lam.real, 0.0, atol=1e-12)
    npt.assert_allclose(np.sort(lam.imag), [-n, 0.0, 0.0, n], atol=1e-12)
    for t in (10.0, 600.0, np.pi / n, 2 * np.pi / n, 5000.0):
        ref = impl.cw_stm_closed_form(n, t)
        npt.assert_allclose(impl.state_transition(A, t), ref, rtol=1e-8, atol=1e-8)
        npt.assert_allclose(expm(A * t), ref, rtol=1e-8, atol=1e-8)
        npt.assert_allclose(np.linalg.det(ref), 1.0, rtol=1e-9)          # tr A = 0


def test_cw_closed_form_satisfies_the_equations_and_shows_the_drift():
    n = 1.1e-3
    x0 = np.array([-50.0, -1000.0, 0.0, 0.0])
    P = 2 * np.pi / n
    sol = solve_ivp(impl.cw_rhs, (0.0, P), x0, args=(n,), rtol=1e-12, atol=1e-12)
    npt.assert_allclose(impl.cw_stm_closed_form(n, P) @ x0, sol.y[:, -1], atol=1e-6)
    # after one orbit the radial offset repeats and the along-track position has drifted by -3 pi (-2 ... ) : y(P) - y0 = 6*(-nP) x0 ... = -12 pi x0
    xP = impl.cw_stm_closed_form(n, P) @ x0
    assert xP[0] == pytest.approx(x0[0], abs=1e-9)
    assert xP[1] - x0[1] == pytest.approx(-12 * np.pi * x0[0], rel=1e-9)
    # a chaser on the same orbit, ahead of the target, stays put
    xs = impl.propagate(impl.cw_matrix(n), np.array([0.0, 500.0, 0.0, 0.0]), np.linspace(0.0, P, 7))
    npt.assert_allclose(xs, np.tile([0.0, 500.0, 0.0, 0.0], (7, 1)), atol=1e-9)


@pytest.mark.parametrize("frac", [0.25, 0.5, 0.75, 1.3])
def test_two_impulse_rendezvous_arrives_at_rest(frac):
    n = impl.iss_orbit()["n"]
    T = frac * 2 * np.pi / n
    r0, v0 = np.array([-50.0, -1000.0]), np.array([0.01, -0.05])
    dv1, dv2 = impl.rendezvous_delta_v(n, r0, v0, T)
    x0 = np.concatenate([r0, v0 + dv1])
    xT = impl.cw_stm_closed_form(n, T) @ x0
    npt.assert_allclose(xT[:2], 0.0, atol=1e-8)
    npt.assert_allclose(xT[2:] + dv2, 0.0, atol=1e-12)
    sol = solve_ivp(impl.cw_rhs, (0.0, T), x0, args=(n,), rtol=1e-12, atol=1e-12)
    npt.assert_allclose(sol.y[:2, -1], 0.0, atol=1e-5)


def test_rendezvous_is_singular_after_a_whole_orbit():
    n = impl.iss_orbit()["n"]
    with pytest.raises(ValueError):
        impl.rendezvous_delta_v(n, np.array([-50.0, -1000.0]), np.zeros(2), 2 * np.pi / n)


def test_thrust_arc_by_zoh_matches_solve_ivp():
    n = impl.iss_orbit()["n"]
    A = impl.cw_matrix(n)
    B = np.array([[0.0, 0.0], [0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    x0 = np.array([-50.0, -1000.0, 0.0, 0.0])
    u = np.tile([[0.0], [2e-4]], (1, 60))          # 10 min of 0.2 mm/s^2 along-track, 10 s holds
    xs = impl.forced_response_zoh(A, B, x0, u, 10.0)
    sol = solve_ivp(impl.cw_rhs, (0.0, 600.0), x0, args=(n, np.array([0.0, 2e-4])), rtol=1e-12, atol=1e-12)
    npt.assert_allclose(xs[-1], sol.y[:, -1], atol=1e-7)


# ---------------------------------------------------------------- figure

def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "state_transition.png").stat().st_size > 1000
