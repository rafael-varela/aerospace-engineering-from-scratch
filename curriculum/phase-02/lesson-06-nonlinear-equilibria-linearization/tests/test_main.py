"""Public tests for P02.06 Nonlinear ODEs, Equilibria & Linearization. They assert properties, with SciPy as the oracle.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-02/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import solve_ivp

from conftest import load_impl

impl = load_impl(__file__)

CASE = impl.uav_swing_test()
# Reference values computed here, not through impl, so the tests collect against an empty starter.
L_EFF = (CASE["I_cg"] + CASE["m"] * CASE["d"] ** 2) / (CASE["m"] * CASE["d"])
W0 = np.sqrt(9.80665 / L_EFF)
T0 = 2.0 * np.pi / W0


# ---------------------------------------------------------------- the compound pendulum

def test_effective_length_and_inertia_round_trip():
    # L_eff = (I_cg + m d^2) / (m d) reduces to d for a point mass (I_cg -> 0) and grows with I_cg.
    assert impl.effective_length(1e-9, 1.0, 0.5) == pytest.approx(0.5, rel=1e-6)
    assert impl.effective_length(0.5, 1.0, 0.5) > 0.5
    # Small-angle period -> inertia recovers I_cg exactly (the formula is its own inverse).
    I_back = impl.inertia_from_period(T0, CASE["m"], CASE["d"])
    assert I_back == pytest.approx(CASE["I_cg"], rel=1e-12)
    with pytest.raises(ValueError):
        impl.effective_length(0.2, -1.0, 0.3)
    with pytest.raises(ValueError):
        impl.inertia_from_period(0.0, 1.0, 1.0)


# ---------------------------------------------------------------- right-hand sides, equilibria, Jacobian

def test_rhs_vanishes_at_both_equilibria_and_matches_linear_model_for_small_angles():
    for eq in impl.pendulum_equilibria():
        npt.assert_allclose(impl.pendulum_rhs(0.0, eq, W0, 0.3), 0.0, atol=1e-12)
    # sin(theta) - theta is O(theta^3): the two models agree to that order.
    for th in (1e-3, 1e-2, 1e-1):
        y = np.array([th, 0.0])
        diff = impl.pendulum_rhs(0.0, y, W0, 0.0) - impl.linear_pendulum_rhs(0.0, y, W0, 0.0)
        assert abs(diff[1]) == pytest.approx(W0 ** 2 * th ** 3 / 6.0, rel=0.01)
    # The inverted equilibrium is at theta = pi, the hanging one at 0.
    eqs = impl.pendulum_equilibria()
    assert len(eqs) == 2
    assert sorted(float(e[0]) for e in eqs) == pytest.approx([0.0, np.pi])


def test_analytic_jacobian_matches_central_differences_and_sympy():
    rng = np.random.default_rng(206)
    for _ in range(5):
        th, zeta = rng.uniform(-np.pi, np.pi), rng.uniform(0.0, 1.5)
        J = impl.pendulum_jacobian(th, W0, zeta)
        Jn = impl.jacobian_numeric(lambda x: impl.pendulum_rhs(0.0, x, W0, zeta), np.array([th, 0.3]))
        npt.assert_allclose(J, Jn, atol=1e-6)                        # central differences: O(eps^2) + rounding
        npt.assert_allclose(J, impl.pendulum_jacobian_sympy(th, W0, zeta), rtol=1e-12, atol=1e-12)
    # The Jacobian of a linear field is the field's matrix, exactly up to rounding.
    A = np.array([[1.0, -2.0], [3.0, 0.5]])
    npt.assert_allclose(impl.jacobian_numeric(lambda x: A @ x, np.array([0.7, -0.2])), A, atol=1e-9)


def test_classification_follows_the_trace_determinant_plane():
    # Hanging equilibrium: center without damping, stable spiral with light damping, stable node when overdamped.
    assert impl.classify_equilibrium(impl.pendulum_jacobian(0.0, W0, 0.0)) == "center"
    assert impl.classify_equilibrium(impl.pendulum_jacobian(0.0, W0, 0.05)) == "stable spiral"
    assert impl.classify_equilibrium(impl.pendulum_jacobian(0.0, W0, 1.5)) == "stable node"
    # Inverted equilibrium: a saddle for every damping (det < 0).
    for zeta in (0.0, 0.05, 1.5):
        assert impl.classify_equilibrium(impl.pendulum_jacobian(np.pi, W0, zeta)) == "saddle"
    # Negative damping flips the verdicts; a singular Jacobian is degenerate.
    assert impl.classify_equilibrium(np.array([[0.0, 1.0], [-4.0, 0.3]])) == "unstable spiral"
    assert impl.classify_equilibrium(np.array([[1.0, 0.0], [0.0, 2.0]])) == "unstable node"
    assert impl.classify_equilibrium(np.array([[0.0, 1.0], [0.0, 0.0]])) == "degenerate"
    with pytest.raises(ValueError):
        impl.classify_equilibrium(np.eye(3))


def test_eigenvalues_agree_with_the_damped_oscillator_formula_and_with_the_saddle():
    zeta = 0.2
    lam = np.sort_complex(impl.pendulum_eigenvalues(0.0, W0, zeta))
    expected = np.sort_complex(np.array([-zeta * W0 - 1j * W0 * np.sqrt(1 - zeta ** 2),
                                         -zeta * W0 + 1j * W0 * np.sqrt(1 - zeta ** 2)]))
    npt.assert_allclose(lam, expected, rtol=1e-12)
    lam_inv = impl.pendulum_eigenvalues(np.pi, W0, 0.0)
    npt.assert_allclose(np.sort(lam_inv.real), [-W0, W0], rtol=1e-12)   # growth rate omega0 at the saddle
    npt.assert_allclose(lam_inv.imag, 0.0, atol=1e-12)


# ---------------------------------------------------------------- RK4 from scratch

def test_rk4_is_fourth_order_against_the_closed_form():
    theta0, zeta = 0.1, 0.05
    f = lambda t, y: impl.linear_pendulum_rhs(t, y, W0, zeta)
    hs = np.array([0.04, 0.02, 0.01, 0.005])
    errs = []
    for h in hs:
        t, Y = impl.integrate_rk4(f, [theta0, 0.0], 2.0 * T0, h)
        errs.append(np.max(np.abs(Y[:, 0] - impl.linear_response(t, theta0, W0, zeta))))
    rates = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert np.all(rates > 3.7), f"expected fourth order, got rates {rates}"
    assert t[0] == 0.0 and t[-1] == pytest.approx(2.0 * T0)
    with pytest.raises(ValueError):
        impl.integrate_rk4(f, [0.1, 0.0], 1.0, -0.01)


def test_rk4_agrees_with_solve_ivp_on_the_nonlinear_pendulum():
    theta0 = np.radians(60.0)
    f = lambda t, y: impl.pendulum_rhs(t, y, W0, 0.02)
    t, Y = impl.integrate_rk4(f, [theta0, 0.0], 4.0 * T0, 0.002)
    sol = solve_ivp(f, (0.0, t[-1]), [theta0, 0.0], t_eval=t, rtol=1e-10, atol=1e-12)
    npt.assert_allclose(Y[:, 0], sol.y[0], atol=1e-7)
    npt.assert_allclose(Y[:, 1], sol.y[1], atol=1e-6)


def test_undamped_energy_is_conserved_and_damped_energy_decays():
    theta0 = np.radians(120.0)
    t, Y = impl.integrate_rk4(lambda t, y: impl.pendulum_rhs(t, y, W0, 0.0), [theta0, 0.0], 5.0 * T0, 0.002)
    E = impl.energy(Y[:, 0], Y[:, 1], W0)
    assert np.max(np.abs(E - E[0])) / E[0] < 1e-9                       # RK4 at h = 2 ms: drift far below 1e-9
    assert E[0] == pytest.approx(W0 ** 2 * (1.0 - np.cos(theta0)))      # released from rest: all potential
    t, Yd = impl.integrate_rk4(lambda t, y: impl.pendulum_rhs(t, y, W0, 0.05), [theta0, 0.0], 5.0 * T0, 0.002)
    Ed = impl.energy(Yd[:, 0], Yd[:, 1], W0)
    assert Ed[-1] < 0.5 * Ed[0]
    # Energy never increases with positive damping (sampled at the turning points, where it is cleanest).
    turning = np.where(np.abs(Yd[:, 1]) < 1e-3)[0]
    assert np.all(np.diff(Ed[turning]) <= 1e-12)


def test_linear_response_solves_the_linear_ode():
    t = np.linspace(0.0, 3.0, 3001)
    th = impl.linear_response(t, 0.2, W0, 0.1)
    dt = t[1] - t[0]
    th_dd = (th[2:] - 2 * th[1:-1] + th[:-2]) / dt ** 2
    th_d = (th[2:] - th[:-2]) / (2 * dt)
    residual = th_dd + 2 * 0.1 * W0 * th_d + W0 ** 2 * th[1:-1]
    assert np.max(np.abs(residual)) < 1e-3 * W0 ** 2 * 0.2            # second-order finite differences at dt = 1 ms
    assert th[0] == pytest.approx(0.2)
    assert abs(th[1] - th[0]) / dt < W0 ** 2 * 0.2 * dt                # released from rest: slope is O(dt)
    assert impl.effective_length(CASE["I_cg"], CASE["m"], CASE["d"]) == pytest.approx(L_EFF)
    assert impl.natural_frequency(L_EFF) == pytest.approx(W0)
    with pytest.raises(ValueError):
        impl.linear_response(t, 0.2, W0, 1.2)


# ---------------------------------------------------------------- period, separatrix

def test_exact_period_limits_series_and_measured_period_agree():
    # Vanishing amplitude: the exact period is T0; the series matches to O(theta^6).
    assert impl.period_exact(1e-6, W0) == pytest.approx(T0, rel=1e-10)
    for th in (0.1, 0.3, 0.6):
        assert impl.period_exact(th, W0) == pytest.approx(impl.period_series(th, W0), rel=0.6 * th ** 6)
    # Monotone in amplitude and 18.03 % longer at 90 degrees (4 K(1/2) / (2 pi) = 1.18034).
    amps = np.linspace(0.0, 3.0, 50)
    T = np.array([impl.period_exact(a, W0) for a in amps])
    assert np.all(np.diff(T) > 0.0)
    assert impl.period_exact(np.pi / 2, W0) / T0 == pytest.approx(1.18034, rel=1e-5)
    # The period measured from an RK4 trajectory at 60 degrees matches the elliptic integral.
    theta0 = np.radians(60.0)
    t, Y = impl.integrate_rk4(lambda t, y: impl.pendulum_rhs(t, y, W0, 0.0), [theta0, 0.0], 6.0 * T0, 0.001)
    assert impl.period_from_trajectory(t, Y[:, 0]) == pytest.approx(impl.period_exact(theta0, W0), rel=1e-5)
    with pytest.raises(ValueError):
        impl.period_exact(np.pi, W0)
    with pytest.raises(ValueError):
        impl.period_from_trajectory(np.linspace(0, 1, 11), np.linspace(-1, 1, 11))


def test_period_scales_as_sqrt_length_and_is_independent_of_mass():
    # omega0 = sqrt(g/L): doubling L multiplies the period by sqrt 2, at any amplitude.
    L2 = 2.0 * L_EFF
    w2 = impl.natural_frequency(L2)
    assert impl.period_exact(1.0, w2) / impl.period_exact(1.0, W0) == pytest.approx(np.sqrt(2.0), rel=1e-12)
    assert impl.period_small_angle(W0) == pytest.approx(T0)


def test_separatrix_divides_swinging_from_rotating():
    assert impl.separatrix_energy(W0) == pytest.approx(2.0 * W0 ** 2)
    # Just below the separatrix the pendulum turns back; just above it goes over the top.
    v_sep = 2.0 * W0
    assert not impl.is_rotating(0.0, 0.999 * v_sep, W0)
    assert impl.is_rotating(0.0, 1.001 * v_sep, W0)
    for v, rotates in ((0.98 * v_sep, False), (1.02 * v_sep, True)):
        t, Y = impl.integrate_rk4(lambda t, y: impl.pendulum_rhs(t, y, W0, 0.0), [0.0, v], 3.0 * T0, 0.001)
        assert (np.max(Y[:, 0]) > np.pi) == rotates


def test_swing_test_error_grows_quadratically_with_amplitude_and_is_amplified_by_the_lever():
    lever = (CASE["I_cg"] + CASE["m"] * CASE["d"] ** 2) / CASE["I_cg"]
    r1 = impl.swing_test_error(np.radians(5.0), CASE)
    r2 = impl.swing_test_error(np.radians(10.0), CASE)
    assert r2["rel_error"] / r1["rel_error"] == pytest.approx(4.0, rel=0.02)       # quadratic in amplitude
    assert r2["rel_error"] / (r2["T_exact"] / r2["T0"] - 1.0) == pytest.approx(2.0 * lever, rel=0.02)
    assert r2["rel_error"] > 0.0                                                   # the linear model always overestimates I
    a_T = impl.max_amplitude_for_period_error(0.01)
    assert impl.period_exact(a_T, W0) / T0 - 1.0 == pytest.approx(0.01, abs=1e-9)
    assert np.degrees(a_T) == pytest.approx(22.8, abs=0.1)                          # P01.05's number, reproduced
    a_I = impl.max_amplitude_for_inertia_error(0.01, CASE)
    assert a_I < a_T
    assert impl.swing_test_error(a_I, CASE)["rel_error"] == pytest.approx(0.01, abs=1e-9)
    with pytest.raises(ValueError):
        impl.max_amplitude_for_period_error(-0.01)


# ---------------------------------------------------------------- the limit cycle

def test_van_der_pol_limit_cycle_is_reached_from_inside_and_outside():
    J = impl.jacobian_numeric(lambda x: impl.van_der_pol_rhs(0.0, x, 1.0), np.zeros(2))
    assert impl.classify_equilibrium(J) == "unstable spiral"          # the only equilibrium repels
    a_in = impl.limit_cycle_amplitude(1.0, [0.1, 0.0])
    a_out = impl.limit_cycle_amplitude(1.0, [4.0, 0.0])
    assert a_in == pytest.approx(a_out, rel=1e-3)                       # same orbit from either side
    assert a_in == pytest.approx(2.0, rel=0.01)                         # amplitude 2.0 + O(mu^2)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "pendulum_linearization.png").stat().st_size > 1000
