"""Public tests for P02.01 First-Order ODEs. They assert physics and convergence, with scipy.integrate.solve_ivp as the oracle.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-02/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import solve_ivp

from conftest import load_impl

impl = load_impl(__file__)

P = impl.sounding_rocket()
TB = (P["m0"] - P["m_dry"]) / P["mdot"]


def fd_derivative(f, t, dt=1e-4):
    """Central difference of a scalar-in, array-out function of time."""
    return (f(t + dt) - f(t - dt)) / (2.0 * dt)


# ---------------------------------------------------------------- separable equations

def test_mass_depletion_is_linear_and_stops_at_the_dry_mass():
    assert impl.burn_time(P["m0"], P["m_dry"], P["mdot"]) == pytest.approx(TB)
    t = np.linspace(0.0, 2.0 * TB, 41)
    m = impl.mass_at(t, P["m0"], P["m_dry"], P["mdot"])
    assert m[0] == pytest.approx(P["m0"])
    # dm/dt = -mdot during the burn, 0 after it.
    burn = t < TB
    npt.assert_allclose(np.diff(m[burn]) / np.diff(t[burn]), -P["mdot"], rtol=1e-12)
    npt.assert_allclose(m[t >= TB], P["m_dry"])
    with pytest.raises(ValueError):
        impl.mass_at(-1.0, P["m0"], P["m_dry"], P["mdot"])
    with pytest.raises(ValueError):
        impl.burn_time(P["m0"], P["m0"] + 1.0, P["mdot"])


def test_tsiolkovsky_is_logarithmic_in_the_mass_ratio():
    dv = impl.tsiolkovsky_delta_v(P["ve"], P["m0"], P["m_dry"])
    assert dv == pytest.approx(P["ve"] * np.log(P["m0"] / P["m_dry"]))
    assert impl.tsiolkovsky_delta_v(P["ve"], P["m0"], P["m0"]) == 0.0
    # Two burns in sequence add: ln(a/b) + ln(b/c) = ln(a/c).
    mid = 160.0
    assert (impl.tsiolkovsky_delta_v(P["ve"], P["m0"], mid) + impl.tsiolkovsky_delta_v(P["ve"], mid, P["m_dry"])
            == pytest.approx(dv))
    with pytest.raises(ValueError):
        impl.tsiolkovsky_delta_v(P["ve"], 100.0, 200.0)


def test_vacuum_velocity_satisfies_the_rocket_equation():
    """dv/dt = mdot ve / m(t) - g, checked by central differences at several times, and v(0) = 0."""
    assert float(impl.vacuum_velocity(0.0, P)) == pytest.approx(0.0, abs=1e-12)
    for t in (1.0, 5.0, 12.0, 19.0):
        lhs = fd_derivative(lambda s: impl.vacuum_velocity(s, P), t)
        rhs = P["mdot"] * P["ve"] / (P["m0"] - P["mdot"] * t) - P["g"]
        assert lhs == pytest.approx(rhs, rel=1e-7)
    # Burnout velocity = Tsiolkovsky delta-v minus the gravity loss g t_b.
    assert float(impl.vacuum_velocity(TB, P)) == pytest.approx(
        impl.tsiolkovsky_delta_v(P["ve"], P["m0"], P["m_dry"]) - P["g"] * TB)
    with pytest.raises(ValueError):
        impl.vacuum_velocity(TB + 1.0, P)


def test_vacuum_altitude_is_the_integral_of_velocity():
    assert float(impl.vacuum_altitude(0.0, P)) == pytest.approx(0.0, abs=1e-9)
    for t in (2.0, 10.0, 18.0):
        assert fd_derivative(lambda s: impl.vacuum_altitude(s, P), t) == pytest.approx(float(impl.vacuum_velocity(t, P)), rel=1e-7)
    # Independent check: trapezoid rule on a fine grid of v(t) reproduces h(t_b).
    t = np.linspace(0.0, TB, 20001)
    h_trap = np.trapezoid(impl.vacuum_velocity(t, P), t)
    assert h_trap == pytest.approx(float(impl.vacuum_altitude(TB, P)), rel=1e-7)


def test_apogee_without_drag_is_energy_conservation():
    assert impl.apogee_from_burnout(1000.0, 0.0, 9.81) == pytest.approx(1000.0)
    assert impl.apogee_from_burnout(0.0, 100.0, 10.0) == pytest.approx(500.0)
    # Gravity-free check: doubling the velocity quadruples the height gained.
    g1 = impl.apogee_from_burnout(0.0, 50.0, 9.81)
    assert impl.apogee_from_burnout(0.0, 100.0, 9.81) == pytest.approx(4.0 * g1)


# ---------------------------------------------------------------- steppers

def exp_rhs(t, y):
    return -2.0 * y


def test_single_steps_match_the_taylor_series_of_exp():
    """One Euler step of y' = -2y matches e^{-2h} to O(h^2); one RK4 step matches it to O(h^5)."""
    y0 = np.array([1.0])
    for h in (0.1, 0.05):
        exact = np.exp(-2.0 * h)
        e_euler = abs(impl.euler_step(exp_rhs, 0.0, y0, h)[0] - exact)
        e_rk4 = abs(impl.rk4_step(exp_rhs, 0.0, y0, h)[0] - exact)
        assert e_euler == pytest.approx(0.5 * (2.0 * h) ** 2, rel=0.3)        # leading term of the Euler error
        assert e_rk4 < (2.0 * h) ** 5 / 120.0 * 1.5                            # first missing Taylor term


def test_euler_is_first_order_and_rk4_fourth_order_on_the_vacuum_burn():
    ns = [10, 20, 40, 80]
    e_euler = impl.convergence_study("euler", ns, P)
    e_rk4 = impl.convergence_study("rk4", ns, P)
    r_euler = np.log2(e_euler[:-1] / e_euler[1:])
    r_rk4 = np.log2(e_rk4[:-1] / e_rk4[1:])
    assert np.all(np.abs(r_euler - 1.0) < 0.1), r_euler
    assert np.all(np.abs(r_rk4 - 4.0) < 0.15), r_rk4
    assert e_rk4[0] < 1e-3 * e_euler[0]                                         # RK4 wins by orders of magnitude at equal cost


def test_fixed_step_integrator_validates_its_inputs_and_returns_the_grid():
    ts, ys = impl.integrate_fixed_step(exp_rhs, 0.0, np.array([1.0]), 1.0, 10, "rk4")
    assert ts.shape == (11,) and ys.shape == (11, 1)
    npt.assert_allclose(ts, np.linspace(0.0, 1.0, 11))
    assert ys[-1, 0] == pytest.approx(np.exp(-2.0), rel=1e-4)                  # h = 0.1: RK4 error ~ h^4
    with pytest.raises(ValueError):
        impl.integrate_fixed_step(exp_rhs, 0.0, np.array([1.0]), 1.0, 10, "rk3")
    with pytest.raises(ValueError):
        impl.integrate_fixed_step(exp_rhs, 0.0, np.array([1.0]), 1.0, 0)


def test_drag_free_simulation_reproduces_the_closed_form_during_the_burn():
    ts, ys = impl.simulate_rocket(P, TB, 400, "rk4", drag=False)
    npt.assert_allclose(ys[:, 1], impl.vacuum_velocity(ts, P), atol=1e-5)
    npt.assert_allclose(ys[:, 0], impl.vacuum_altitude(ts, P), atol=1e-4)
    npt.assert_allclose(ys[:, 2], impl.mass_at(ts, P["m0"], P["m_dry"], P["mdot"]), rtol=1e-12)


def test_full_model_agrees_with_solve_ivp_and_drag_costs_velocity():
    t_end = 60.0
    ts, ys = impl.simulate_rocket(P, t_end, 1200, "rk4", drag=True)
    i_b = int(np.argmin(np.abs(ts - TB)))
    assert ts[i_b] == pytest.approx(TB)                                         # the switch-off sits on the grid
    burn = solve_ivp(lambda t, y: impl.rocket_rhs(t, y, P, True, True), (0.0, TB), [0.0, 0.0, P["m0"]],
                     rtol=1e-10, atol=1e-8)
    coast = solve_ivp(lambda t, y: impl.rocket_rhs(t, y, P, False, True), (TB, t_end), burn.y[:, -1],
                      rtol=1e-10, atol=1e-8, t_eval=ts[ts >= TB])
    assert ys[i_b, 1] == pytest.approx(burn.y[1, -1], abs=1e-4)
    npt.assert_allclose(ys[ts >= TB, 0], coast.y[0], atol=1e-2)
    # Drag only ever removes momentum: burnout velocity and altitude are below the vacuum values.
    assert ys[i_b, 1] < float(impl.vacuum_velocity(TB, P))
    assert ys[i_b, 0] < float(impl.vacuum_altitude(TB, P))
    # Mass is piecewise linear and never below the dry mass.
    assert ys[:, 2].min() == pytest.approx(P["m_dry"])


def test_rhs_components_have_the_right_signs_and_limits():
    y = np.array([0.0, 0.0, P["m0"]])
    d = impl.rocket_rhs(0.0, y, P, True, True)
    assert d[0] == 0.0                                                          # dh/dt = v
    assert d[1] == pytest.approx(P["mdot"] * P["ve"] / P["m0"] - P["g"])         # no drag at rest
    assert d[2] == pytest.approx(-P["mdot"])
    d_coast = impl.rocket_rhs(0.0, np.array([0.0, 100.0, P["m_dry"]]), P, False, True)
    assert d_coast[2] == 0.0 and d_coast[1] < -P["g"]                           # drag adds to gravity when rising
    d_fall = impl.rocket_rhs(0.0, np.array([0.0, -100.0, P["m_dry"]]), P, False, True)
    assert -P["g"] < d_fall[1] < 0.0                                             # drag opposes a fall
    assert impl.rocket_rhs(0.0, np.array([0.0, 100.0, P["m_dry"]]), P, False, False)[1] == pytest.approx(-P["g"])
    # Exponential atmosphere: sea-level value and one scale height.
    assert float(impl.air_density(0.0)) == pytest.approx(1.225)
    assert float(impl.air_density(8500.0)) == pytest.approx(1.225 / np.e)


# ---------------------------------------------------------------- coast closed forms

def test_quadratic_drag_coast_matches_rk4_and_reduces_to_vacuum():
    v_b, k, g = 800.0, 2e-4, 9.81
    t_apex, h_gain = impl.quadratic_drag_coast(v_b, k, g)
    # Closed form satisfies dv/dt = -g - k v^2 and v(t_apex) = 0.
    for t in (0.1 * t_apex, 0.3 * t_apex, 0.9 * t_apex):
        lhs = fd_derivative(lambda s: impl.quadratic_drag_velocity(s, v_b, k, g), t)
        v = float(impl.quadratic_drag_velocity(t, v_b, k, g))
        assert lhs == pytest.approx(-g - k * v * v, rel=1e-6)
    assert float(impl.quadratic_drag_velocity(t_apex, v_b, k, g)) == pytest.approx(0.0, abs=1e-8)
    # RK4 on (h, v) reproduces the height gained.
    rhs = lambda t, y: np.array([y[1], -g - k * y[1] ** 2])
    _, ys = impl.integrate_fixed_step(rhs, 0.0, np.array([0.0, v_b]), t_apex, 2000, "rk4")
    assert ys[-1, 0] == pytest.approx(h_gain, rel=1e-8)
    assert abs(ys[-1, 1]) < 1e-6
    # k -> 0 recovers the vacuum values continuously.
    t0, h0 = impl.quadratic_drag_coast(v_b, 0.0, g)
    t1, h1 = impl.quadratic_drag_coast(v_b, 1e-9, g)
    assert t0 == pytest.approx(v_b / g) and h0 == pytest.approx(v_b ** 2 / (2 * g))
    assert t1 == pytest.approx(t0, rel=1e-4) and h1 == pytest.approx(h0, rel=1e-4)
    assert t_apex < t0 and h_gain < h0                                          # drag shortens both


def test_linear_drag_solution_satisfies_its_ode_and_tends_to_terminal_velocity():
    v0, c, m, g = 300.0, 2.0, 100.0, 9.81
    assert float(impl.linear_drag_velocity(0.0, v0, c, m, g)) == pytest.approx(v0)
    for t in (1.0, 20.0, 100.0):
        lhs = fd_derivative(lambda s: impl.linear_drag_velocity(s, v0, c, m, g), t)
        v = float(impl.linear_drag_velocity(t, v0, c, m, g))
        assert lhs + (c / m) * v == pytest.approx(-g, rel=1e-6)
    # After many time constants m/c the velocity is the terminal value -m g / c, whatever v0.
    t_long = 20.0 * m / c
    for v_start in (-500.0, 0.0, 300.0):
        assert float(impl.linear_drag_velocity(t_long, v_start, c, m, g)) == pytest.approx(-m * g / c, rel=1e-6)
    npt.assert_allclose(impl.linear_drag_velocity(np.array([0.0, 3.0]), v0, 0.0, m, g), [v0, v0 - 3 * g])
    with pytest.raises(ValueError):
        impl.linear_drag_velocity(1.0, v0, -1.0, m, g)


def test_terminal_velocity_is_a_stable_equilibrium_of_the_autonomous_equation():
    vt = impl.terminal_velocity(P["m_dry"], 1.225, P["cd"], P["area"], P["g"])
    assert vt == pytest.approx(np.sqrt(2 * P["m_dry"] * P["g"] / (1.225 * P["cd"] * P["area"])))
    k = 1.225 * P["cd"] * P["area"] / (2.0 * P["m_dry"])
    f = lambda v: -P["g"] - k * v * abs(v)
    eq = impl.equilibria(f, np.linspace(-600.0, 600.0, 121))
    assert len(eq) == 1
    v_star, stable = eq[0]
    assert v_star == pytest.approx(-vt, rel=1e-8)
    assert stable
    # Scaling: four times the mass doubles the terminal velocity.
    assert impl.terminal_velocity(4 * P["m_dry"], 1.225, P["cd"], P["area"], P["g"]) == pytest.approx(2 * vt)
    with pytest.raises(ValueError):
        impl.terminal_velocity(-1.0, 1.225, P["cd"], P["area"], P["g"])


def test_equilibria_classifies_the_logistic_equation():
    """y' = y (1 - y): y = 0 unstable, y = 1 stable (the phase line of the logistic equation)."""
    eq = impl.equilibria(lambda y: y * (1.0 - y), np.linspace(-0.5, 1.5, 41))
    roots = sorted(eq)
    assert len(roots) == 2
    assert roots[0][0] == pytest.approx(0.0, abs=1e-10) and roots[0][1] is False
    assert roots[1][0] == pytest.approx(1.0, abs=1e-10) and roots[1][1] is True


def test_direction_field_arrows_are_unit_length_and_follow_the_slope():
    tg, yg = np.linspace(0.0, 1.0, 5), np.linspace(-1.0, 1.0, 7)
    U, V = impl.direction_field(lambda t, y: t - y, tg, yg)
    assert U.shape == V.shape == (7, 5)
    npt.assert_allclose(U ** 2 + V ** 2, 1.0, atol=1e-12)
    T, Y = np.meshgrid(tg, yg)
    npt.assert_allclose(V / U, T - Y, atol=1e-12)
    assert np.all(U > 0.0)                                                      # time always runs forward


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "rocket_ascent.png").stat().st_size > 1000
