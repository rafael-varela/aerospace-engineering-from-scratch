"""Public tests for P02.08 Laplace Transform. They assert properties, with scipy.signal, solve_ivp and SymPy as the oracles.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-02/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import solve_ivp
from scipy.signal import lti

from conftest import load_impl

impl = load_impl(__file__)

ACT = impl.elevator_actuator()
WN, ZETA = ACT["wn"], ACT["zeta"]
T = np.linspace(0.0, 0.6, 3001)


# ---------------------------------------------------------------- the transform itself

def test_quadrature_reproduces_the_transform_table():
    s = 1.7
    cases = [("step", lambda t: 1.0, {}), ("ramp", lambda t: t, {}), ("exp", lambda t: np.exp(-2.0 * t), {"a": 2.0}),
             ("cos", lambda t: np.cos(3.0 * t), {"w": 3.0}), ("sin", lambda t: np.sin(3.0 * t), {"w": 3.0}),
             ("damped_sin", lambda t: np.exp(-0.5 * t) * np.sin(4.0 * t), {"a": 0.5, "w": 4.0}),
             ("damped_cos", lambda t: np.exp(-0.5 * t) * np.cos(4.0 * t), {"a": 0.5, "w": 4.0}),
             ("delayed_step", lambda t: 1.0 if t >= 0.8 else 0.0, {"c": 0.8})]
    for name, f, p in cases:
        assert impl.laplace_numeric(f, s) == pytest.approx(impl.transform_table(name, s, **p), rel=1e-7, abs=1e-9)
    with pytest.raises(ValueError):
        impl.laplace_numeric(lambda t: 1.0, 0.0)
    with pytest.raises(KeyError):
        impl.transform_table("no-such-entry", s)


def test_linearity_derivative_and_shift_rules_hold_numerically():
    s = 2.3
    # Linearity: L{3 e^{-t} + 2 cos 4t} = 3/(s+1) + 2 s/(s^2+16).
    assert impl.laplace_numeric(lambda t: 3 * np.exp(-t) + 2 * np.cos(4 * t), s) == pytest.approx(
        3 * impl.transform_table("exp", s, a=1.0) + 2 * impl.transform_table("cos", s, w=4.0), rel=1e-8)
    # Derivative rule: L{f'} = s F(s) - f(0) for f = e^{-2t} cos 3t (f(0) = 1, f' = -2 e^{-2t} cos 3t - 3 e^{-2t} sin 3t).
    F = impl.transform_table("damped_cos", s, a=2.0, w=3.0)
    fprime = lambda t: -2 * np.exp(-2 * t) * np.cos(3 * t) - 3 * np.exp(-2 * t) * np.sin(3 * t)
    assert impl.laplace_numeric(fprime, s) == pytest.approx(s * F - 1.0, rel=1e-8)
    # Shift rule: L{f(t - c) u_c(t)} = e^{-cs} F(s) for f = sin 3t, c = 0.4.
    c = 0.4
    shifted = lambda t: np.sin(3 * (t - c)) if t >= c else 0.0
    assert impl.laplace_numeric(shifted, s) == pytest.approx(np.exp(-c * s) * impl.transform_table("sin", s, w=3.0), rel=1e-7)
    # s-shift: L{e^{-at} f(t)} = F(s + a).
    assert impl.transform_table("damped_sin", s, a=0.7, w=3.0) == pytest.approx(impl.transform_table("sin", s + 0.7, w=3.0))


# ---------------------------------------------------------------- transfer functions, poles, residues

def test_transfer_function_static_gain_poles_and_stability():
    num, den = impl.transfer_function(WN, ZETA)
    assert impl.evaluate(num, den, 0.0) == pytest.approx(1.0)                       # unit static gain
    p = impl.poles(den)
    npt.assert_allclose(np.sort_complex(p), np.sort_complex(np.array([-ZETA * WN - 1j * WN * np.sqrt(1 - ZETA ** 2),
                                                                      -ZETA * WN + 1j * WN * np.sqrt(1 - ZETA ** 2)])), rtol=1e-12)
    assert impl.is_stable(den)
    assert not impl.is_stable(np.array([1.0, -2 * ZETA * WN, WN ** 2]))             # negative damping: unstable
    assert not impl.is_stable(np.array([1.0, 0.0, WN ** 2]))                         # zero damping: marginal, not stable
    # Overdamped: two real negative poles whose product is wn^2 and sum -2 zeta wn.
    p2 = impl.poles(impl.transfer_function(WN, 2.0)[1])
    npt.assert_allclose(p2.imag, 0.0, atol=1e-12)
    assert np.prod(p2).real == pytest.approx(WN ** 2) and np.sum(p2).real == pytest.approx(-4.0 * WN)
    with pytest.raises(ValueError):
        impl.transfer_function(-1.0, 0.5)


def test_residues_match_the_cover_up_rule_and_reconstruct_the_function():
    # 1 / ((s + 1)(s + 3)) = (1/2)/(s + 1) - (1/2)/(s + 3): the cover-up rule by hand.
    p, r = impl.residues(np.array([1.0]), np.array([1.0, 4.0, 3.0]))
    order = np.argsort(p.real)
    npt.assert_allclose(p[order], [-3.0, -1.0], atol=1e-12)
    npt.assert_allclose(r[order], [-0.5, 0.5], atol=1e-12)
    # A general strictly proper function: the partial-fraction sum reproduces N/D at random points.
    rng = np.random.default_rng(208)
    num, den = rng.normal(size=2), np.array([1.0, 6.0, 25.0, 30.0])
    p, r = impl.residues(num, den)
    for s in rng.normal(size=5) + 1j * rng.normal(size=5):
        assert np.sum(r / (s - p)) == pytest.approx(impl.evaluate(num, den, s), rel=1e-10)
    with pytest.raises(ValueError):
        impl.residues(np.array([1.0, 0.0, 0.0]), np.array([1.0, 1.0, 1.0]))          # not strictly proper
    with pytest.raises(ValueError):
        impl.residues(np.array([1.0]), np.array([1.0, 2.0, 1.0]))                    # repeated pole at -1


def test_inverse_transform_reproduces_table_entries_and_is_causal():
    t = np.linspace(-0.2, 2.0, 1101)
    # 1/(s+2) -> e^{-2t}; 3/(s^2+9) -> sin 3t; (s+1)/((s+1)^2+4) -> e^{-t} cos 2t.
    f1 = impl.inverse_laplace(np.array([1.0]), np.array([1.0, 2.0]), t)
    f2 = impl.inverse_laplace(np.array([3.0]), np.array([1.0, 0.0, 9.0]), t)
    f3 = impl.inverse_laplace(np.array([1.0, 1.0]), np.array([1.0, 2.0, 5.0]), t)
    pos = t >= 0.0
    npt.assert_allclose(f1[pos], np.exp(-2 * t[pos]), atol=1e-12)
    npt.assert_allclose(f2[pos], np.sin(3 * t[pos]), atol=1e-12)
    npt.assert_allclose(f3[pos], np.exp(-t[pos]) * np.cos(2 * t[pos]), atol=1e-12)
    for f in (f1, f2, f3):
        npt.assert_allclose(f[~pos], 0.0)                                            # zero before t = 0
        assert f.dtype == float
    # Round trip: the quadrature transform of the inverse equals the rational function.
    s = 1.5
    assert impl.laplace_numeric(lambda tt: impl.inverse_laplace(np.array([1.0, 1.0]), np.array([1.0, 2.0, 5.0]), np.array([tt]))[0], s) \
        == pytest.approx(impl.evaluate(np.array([1.0, 1.0]), np.array([1.0, 2.0, 5.0]), s), rel=1e-7)


# ---------------------------------------------------------------- responses of the actuator

def test_step_response_matches_closed_form_scipy_and_solve_ivp():
    y = impl.step_response(WN, ZETA, T)
    npt.assert_allclose(y, impl.step_response_closed_form(WN, ZETA, T), atol=1e-12)
    num, den = impl.transfer_function(WN, ZETA)
    _, y_lti = lti(num, den).step(T=T)
    npt.assert_allclose(y, y_lti, atol=1e-9)
    sol = solve_ivp(impl.actuator_rhs, (0.0, T[-1]), [0.0, 0.0], t_eval=T, args=(WN, ZETA, lambda t: 1.0), rtol=1e-10, atol=1e-12)
    npt.assert_allclose(y, sol.y[0], atol=1e-8)
    assert y[0] == pytest.approx(0.0, abs=1e-12) and y[-1] == pytest.approx(1.0, abs=2e-4)   # from rest to the static gain (e^{-9} left at 0.6 s)
    with pytest.raises(ValueError):
        impl.step_response_closed_form(WN, 1.2, T)


def test_initial_conditions_enter_through_the_derivative_rule():
    num, den = impl.transfer_function(WN, ZETA)
    y0, yd0 = 0.4, -3.0
    Y_num, Y_den = impl.response_with_initial_conditions(num, den, np.array([1.0]), np.array([1.0, 0.0]), y0, yd0)
    y = impl.inverse_laplace(Y_num, Y_den, T)
    sol = solve_ivp(impl.actuator_rhs, (0.0, T[-1]), [y0, yd0], t_eval=T, args=(WN, ZETA, lambda t: 1.0), rtol=1e-10, atol=1e-12)
    npt.assert_allclose(y, sol.y[0], atol=1e-8)
    assert y[0] == pytest.approx(y0, abs=1e-12)
    assert (y[1] - y[0]) / (T[1] - T[0]) == pytest.approx(yd0, rel=0.05)              # initial slope on a 0.2 ms grid
    # Zero input with initial conditions: the free response, which decays to zero.
    Yn, Yd = impl.response_with_initial_conditions(num, den, np.array([0.0]), np.array([1.0]), y0, yd0)
    y_free = impl.inverse_laplace(Yn, Yd, T)
    assert abs(y_free[-1]) < 1e-3 * abs(y0)
    with pytest.raises(ValueError):
        impl.response_with_initial_conditions(num, np.array([2.0, 1.0, 1.0]), np.array([1.0]), np.array([1.0, 0.0]), 0.0, 0.0)


def test_impulse_response_is_the_derivative_of_the_step_and_integrates_to_the_static_gain():
    t = np.linspace(0.0, 1.0, 20001)
    y = impl.step_response(WN, ZETA, t)
    g = impl.impulse_response(WN, ZETA, t)
    dt = t[1] - t[0]
    npt.assert_allclose(g[2:-2], np.gradient(y, dt)[2:-2], atol=2e-3 * WN)
    assert np.trapezoid(g, t) == pytest.approx(1.0, abs=1e-5)                       # G(0) = 1
    _, g_lti = lti(*impl.transfer_function(WN, ZETA)).impulse(T=t)
    npt.assert_allclose(g, g_lti, atol=1e-8 * WN)
    # Convolution with a step reproduces the step response (first-order rule: error O(dt)).
    npt.assert_allclose(impl.convolve_causal(g, np.ones_like(t), dt), y, atol=0.01)


def test_delayed_step_is_the_shifted_step_and_superposes_into_a_pulse():
    t = np.linspace(0.0, 1.0, 2001)
    c = 0.3
    yd = impl.delayed_step_response(WN, ZETA, t, c)
    y = impl.step_response(WN, ZETA, t)
    npt.assert_allclose(yd[t < c], 0.0, atol=1e-12)                                  # nothing happens before the step
    npt.assert_allclose(yd[t >= c], y[: np.sum(t >= c)], atol=1e-9)                   # same shape, shifted by c
    # A 0.1 s pulse = step minus delayed step; integrate the ODE with that input and compare.
    pulse = y - impl.delayed_step_response(WN, ZETA, t, 0.1)
    sol = solve_ivp(impl.actuator_rhs, (0.0, 1.0), [0.0, 0.0], t_eval=t, args=(WN, ZETA, lambda tt: 1.0 if tt < 0.1 else 0.0),
                    rtol=1e-10, atol=1e-12, max_step=1e-3)
    npt.assert_allclose(pulse, sol.y[0], atol=1e-6)
    with pytest.raises(ValueError):
        impl.delayed_step_response(WN, ZETA, t, -0.1)


def test_overshoot_peak_time_and_settling_time_agree_with_the_response():
    t = np.linspace(0.0, 1.0, 100001)
    for zeta in (0.2, ZETA, 0.9):
        y = impl.step_response(WN, zeta, t)
        assert y.max() - 1.0 == pytest.approx(impl.overshoot(zeta), rel=1e-4)
        assert t[np.argmax(y)] == pytest.approx(impl.peak_time(WN, zeta), rel=1e-3)
        ts = impl.settling_time(t, y, 1.0)
        assert np.all(np.abs(y[t >= ts] - 1.0) <= 0.02 + 1e-12)
        assert ts < 4.0 / (zeta * WN) * 1.2 and ts > 0.0                              # the 4 / (zeta wn) rule is a loose upper estimate
    assert impl.overshoot(0.0) == pytest.approx(1.0)                                 # undamped: doubles
    assert impl.overshoot(0.7) < impl.overshoot(0.3)
    with pytest.raises(ValueError):
        impl.overshoot(1.0)
    with pytest.raises(ValueError):
        impl.settling_time(np.array([0.0, 1.0]), np.array([0.0, 0.5]), 1.0)


# ---------------------------------------------------------------- frequency response and the symbolic check

def test_frequency_response_is_the_complex_gain_of_the_previous_lesson():
    num, den = impl.transfer_function(WN, ZETA)
    omega = np.array([0.0, 5.0, WN, 100.0])
    H = impl.frequency_response(num, den, omega)
    r = omega / WN
    npt.assert_allclose(H, 1.0 / (1.0 - r ** 2 + 2j * ZETA * r), rtol=1e-12)
    assert H[0] == pytest.approx(1.0)
    assert abs(H[2]) == pytest.approx(1.0 / (2 * ZETA)) and np.angle(H[2]) == pytest.approx(-np.pi / 2)
    # Steady-state sinusoidal response from solve_ivp has the amplitude |H| (after the transient).
    w = 20.0
    t = np.linspace(0.0, 2.0, 20001)
    sol = solve_ivp(impl.actuator_rhs, (0.0, 2.0), [0.0, 0.0], t_eval=t, args=(WN, ZETA, lambda tt: np.cos(w * tt)), rtol=1e-10, atol=1e-12)
    tail = t > 1.0
    amp = 0.5 * (sol.y[0][tail].max() - sol.y[0][tail].min())
    assert amp == pytest.approx(abs(impl.frequency_response(num, den, np.array([w]))[0]), rel=1e-3)


def test_sympy_step_response_agrees_with_the_residue_arithmetic():
    import sympy as sp
    expr = impl.laplace_symbolic_step_response(WN, ZETA)
    tsym = sp.Symbol("t", positive=True)
    f = sp.lambdify(tsym, expr, "numpy")
    tt = np.linspace(0.01, 0.5, 50)
    npt.assert_allclose(f(tt), impl.step_response(WN, ZETA, tt), atol=1e-10)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "laplace_transform.png").stat().st_size > 1000
