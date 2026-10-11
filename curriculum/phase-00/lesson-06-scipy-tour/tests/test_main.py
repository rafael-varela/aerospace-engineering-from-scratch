"""Public tests for P00.06 SciPy Tour.

Every check compares a SciPy result with an analytic value derived in the lesson:
the exact integral of the thrust model, its exact peak, the closed-form ISA pressure,
and a balance system whose answer was constructed. Run `pytest <lesson dir>` against
your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>` against the
reference solution. Hidden validation cases live in assessments/phase-00/parameterized-cases/.
"""
import numpy as np
import pytest

from conftest import load_impl

impl = load_impl(__file__)

A, TAU, T_END = 60.0, 0.3, 2.0


# ------------------------------------------------------------------ model and analytic reference

def test_thrust_model_limits():
    t = np.array([-1.0, 0.0, TAU])
    f = impl.thrust_model(t, A, TAU)
    assert f[0] == 0.0 and f[1] == 0.0
    assert f[2] == pytest.approx(A * TAU / np.e)


def test_analytic_impulse_limits():
    assert impl.analytic_impulse(A, TAU, 0.0) == pytest.approx(0.0, abs=1e-15)
    # as t_end -> infinity the impulse tends to a * tau**2
    assert impl.analytic_impulse(A, TAU, 50.0 * TAU) == pytest.approx(A * TAU**2, rel=1e-12)
    assert impl.analytic_impulse(2 * A, TAU, T_END) == pytest.approx(2 * impl.analytic_impulse(A, TAU, T_END))


def test_synthesize_thrust_data_is_reproducible_and_exact_without_noise():
    t, f = impl.synthesize_thrust_data(A, TAU, T_END, n=41, noise=0.0)
    assert t.shape == f.shape == (41,)
    assert np.allclose(f, impl.thrust_model(t, A, TAU))
    t1, f1 = impl.synthesize_thrust_data(A, TAU, T_END, n=41, noise=0.5, seed=7)
    t2, f2 = impl.synthesize_thrust_data(A, TAU, T_END, n=41, noise=0.5, seed=7)
    assert np.array_equal(f1, f2)
    assert not np.allclose(f1, impl.thrust_model(t1, A, TAU))


# ------------------------------------------------------------------ integrate

def test_quad_matches_analytic_impulse():
    assert impl.total_impulse_model(A, TAU, T_END) == pytest.approx(impl.analytic_impulse(A, TAU, T_END), rel=1e-8)


def test_trapezoid_matches_analytic_on_fine_grid():
    t, f = impl.synthesize_thrust_data(A, TAU, T_END, n=2001, noise=0.0)
    assert impl.total_impulse_samples(t, f) == pytest.approx(impl.analytic_impulse(A, TAU, T_END), rel=1e-5)


def test_trapezoid_is_second_order():
    errs = impl.trapezoid_convergence(A, TAU, T_END, [11, 21, 41, 81, 161])
    rates = np.log2(errs[:-1] / errs[1:])
    assert np.all(rates > 1.8), f"expected second-order convergence, got {rates}"


def test_total_impulse_samples_rejects_mismatched_input():
    with pytest.raises(ValueError):
        impl.total_impulse_samples(np.array([0.0, 1.0, 2.0]), np.array([1.0, 2.0]))


def test_motor_class_boundaries_are_monotone():
    assert impl.motor_class(2.0) == "A"
    assert impl.motor_class(2.5) == "A"
    assert impl.motor_class(2.6) == "B"
    assert impl.motor_class(7.0) == "C"
    assert impl.motor_class(15.0) == "D"
    assert impl.motor_class(1.0) == "micro"
    with pytest.raises(ValueError):
        impl.motor_class(0.0)
    letters = [impl.motor_class(x) for x in np.geomspace(1.3, 1000.0, 40)]
    assert letters == sorted(letters), "class letter must never decrease with impulse"


# ------------------------------------------------------------------ optimize

def test_fit_recovers_exact_parameters_without_noise():
    t, f = impl.synthesize_thrust_data(A, TAU, T_END, n=81, noise=0.0)
    a_fit, tau_fit = impl.fit_thrust_model(t, f)
    assert a_fit == pytest.approx(A, rel=1e-6)
    assert tau_fit == pytest.approx(TAU, rel=1e-6)


def test_fit_is_close_with_noise():
    t, f = impl.synthesize_thrust_data(A, TAU, T_END, n=161, noise=0.3, seed=1)
    a_fit, tau_fit = impl.fit_thrust_model(t, f)
    assert a_fit == pytest.approx(A, rel=0.1)
    assert tau_fit == pytest.approx(TAU, rel=0.1)
    # the fitted impulse is within a few percent of the true one
    assert impl.analytic_impulse(a_fit, tau_fit, T_END) == pytest.approx(impl.analytic_impulse(A, TAU, T_END), rel=0.05)


def test_peak_thrust_matches_analytic():
    t_pk, f_pk = impl.peak_thrust(A, TAU, T_END)
    assert t_pk == pytest.approx(TAU, abs=1e-6)
    assert f_pk == pytest.approx(A * TAU / np.e, rel=1e-8)


def test_burnout_time_is_a_root_after_the_peak():
    t_b = impl.burnout_time(A, TAU, 1.0, T_END)
    assert t_b > TAU
    assert impl.thrust_model(t_b, A, TAU) == pytest.approx(1.0, abs=1e-8)
    with pytest.raises(ValueError):
        impl.burnout_time(A, TAU, 100.0, T_END)  # above the peak
    with pytest.raises(ValueError):
        impl.burnout_time(A, TAU, 1e-9, 0.5)     # not yet reached


# ------------------------------------------------------------------ interpolate

def test_isa_exact_reference_values():
    p = impl.isa_pressure_exact(np.array([0.0, 11_000.0, 20_000.0]))
    assert p[0] == pytest.approx(101_325.0)
    assert p[1] == pytest.approx(22_632.0, rel=1e-3)
    assert p[2] == pytest.approx(5_474.9, rel=1e-3)


@pytest.mark.parametrize("kind", ["linear", "log-linear", "cubic"])
def test_interpolation_reproduces_table_nodes(kind):
    p = impl.interpolate_pressure(impl.ISA_TABLE_H, kind)
    assert np.allclose(p, impl.ISA_TABLE_P, rtol=1e-10)


def test_interpolation_is_monotone_and_log_linear_beats_linear():
    h = np.linspace(0.0, 20_000.0, 401)
    exact = impl.isa_pressure_exact(h)
    err = {}
    for kind in ("linear", "log-linear", "cubic"):
        p = impl.interpolate_pressure(h, kind)
        assert np.all(np.diff(p) <= 0.0), f"{kind}: pressure must not rise with altitude"
        err[kind] = np.max(np.abs(p - exact) / exact)
    assert err["log-linear"] < err["linear"]
    assert err["log-linear"] < 0.01
    assert err["cubic"] < err["log-linear"]


def test_interpolation_refuses_to_extrapolate():
    with pytest.raises(ValueError):
        impl.interpolate_pressure(np.array([25_000.0]))
    with pytest.raises(ValueError):
        impl.interpolate_pressure(np.array([-10.0]))


# ------------------------------------------------------------------ linalg

def test_balance_round_trip_and_residual():
    truth = np.array([120.0, 8.0, 3.0])
    readings = impl.balance_matrix() @ truth
    x = impl.resolve_balance(readings)
    assert np.allclose(x, truth, rtol=1e-12, atol=1e-12)
    assert np.allclose(impl.balance_matrix() @ x - readings, 0.0, atol=1e-12)


def test_balance_hand_check():
    # equal cell readings mean zero moment; lift is their sum
    x = impl.resolve_balance(np.array([50.0, 50.0, 2.0]))
    assert x[0] == pytest.approx(100.0)
    assert x[2] == pytest.approx(0.0, abs=1e-12)
    # drag cell reads D + 0.02 L -> D = 2 - 0.02*100 = 0
    assert x[1] == pytest.approx(0.0, abs=1e-12)


def test_balance_is_linear():
    r1 = np.array([70.0, 50.0, 10.4])
    r2 = np.array([10.0, -5.0, 1.0])
    assert np.allclose(impl.resolve_balance(r1 + r2), impl.resolve_balance(r1) + impl.resolve_balance(r2))


def test_balance_condition_number_grows_as_arm_shrinks():
    assert impl.balance_condition_number() >= 1.0
    assert impl.balance_condition_number(arm=0.015) > impl.balance_condition_number(arm=0.15)
    with pytest.raises(ValueError):
        impl.resolve_balance(np.array([1.0, 2.0]))


def test_main_runs_and_writes_figure(tmp_path, capsys):
    impl.main(tmp_path)
    out = capsys.readouterr().out
    assert "integrate" in out and "linalg" in out
    assert (tmp_path / "thrust_fit.svg").exists()
