"""Public tests for P03.18 Verification & Validation. They verify the tridiagonal solver against
scipy.linalg.solve_banded and by its backward error, the heat solver against exact discrete, separable
and manufactured solutions (observed order two for Crank-Nicolson, one for backward Euler, and the
planted bugs that the order test must catch), the Richardson/GCI routine on an exact sequence and on a
real functional, and the validation metric and calibration against the lesson's SYNTHETIC thermocouple
data (no measurement anywhere in this lesson).

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-03/parameterized-cases/.
"""
import math

import numpy as np
import numpy.testing as npt
import pytest
import sympy as sp
from scipy.linalg import solve_banded

from conftest import load_impl

impl = load_impl(__file__)

P = impl.HeatSetup()


def _banded(lower, diag, upper, b):
    n = diag.size
    ab = np.zeros((3, n))
    ab[0, 1:] = upper[:-1]
    ab[1] = diag
    ab[2, :-1] = lower[1:]
    return solve_banded((1, 1), ab, b)


def _mms():
    x, t = sp.symbols("x t")
    tau = 200.0
    u = 40 * sp.sin(sp.pi * x / P.L) * (1 - sp.exp(-t / tau)) + 10 * sp.sin(2 * sp.pi * x / P.L) * sp.exp(-t / tau)
    return impl.manufactured_source(u, x, t, P.alpha)


def _mms_error(nx, theta=0.5, bug=None, t_end=300.0):
    u_m, q_m = _mms()
    xg = np.linspace(0.0, P.L, nx + 1)
    U = impl.heat_theta(xg, u_m(xg, 0.0), t_end, nx, P.alpha, source=q_m, theta=theta, bug=bug)
    return impl.l2_error(U[-1], u_m(xg, t_end), xg[1] - xg[0])


# ---------------------------------------------------------------- linear solver

def test_thomas_matches_scipy_solve_banded():
    rng = np.random.default_rng(18)
    n = 150
    lower, upper = rng.standard_normal(n), rng.standard_normal(n)
    diag = np.abs(lower) + np.abs(upper) + 1.0 + rng.random(n)        # diagonally dominant
    b = rng.standard_normal(n)
    x = impl.thomas_solve(lower, diag, upper, b)
    npt.assert_allclose(x, _banded(lower, diag, upper, b), rtol=1e-12, atol=1e-13)


def test_thomas_is_backward_stable_even_when_ill_conditioned():
    n = 200
    lower = np.full(n, -1.0)
    upper = np.full(n, -1.0)
    diag = np.full(n, 2.0 - 1e-6)                                     # condition number ~ 1.6e4
    b = np.random.default_rng(1).standard_normal(n)
    x = impl.thomas_solve(lower, diag, upper, b)
    be = impl.backward_error(lower, diag, upper, x, b)
    assert be < 20 * np.finfo(float).eps                              # a few ulps whatever the conditioning
    npt.assert_allclose(x, _banded(lower, diag, upper, b), rtol=1e-9)


def test_thomas_raises_on_zero_pivot():
    with pytest.raises(ValueError):
        impl.thomas_solve(np.array([0.0, 1.0]), np.array([0.0, 1.0]), np.array([1.0, 0.0]), np.ones(2))


# ---------------------------------------------------------------- exact solutions of the discrete equations

def test_constant_solution_is_reproduced_to_rounding():
    x = np.linspace(0.0, P.L, 31)
    U = impl.heat_theta(x, np.full(x.size, 7.0), 500.0, 20, P.alpha, bc=(lambda t: 7.0, lambda t: 7.0))
    npt.assert_allclose(U, 7.0, rtol=0, atol=1e-12)


def test_polynomial_solution_is_an_exact_discrete_solution():
    # u = 2 t x (L - x): central differences are exact for quadratics, the trapezoidal time average
    # for a source linear in t, so Crank-Nicolson reproduces it to rounding on any grid.
    x = np.linspace(0.0, P.L, 21)
    u = lambda xx, t: 2.0 * t * xx * (P.L - xx)
    q = lambda xx, t: 2.0 * xx * (P.L - xx) + 4.0 * P.alpha * t
    U = impl.heat_theta(x, u(x, 0.0), 200.0, 8, P.alpha, source=q)
    npt.assert_allclose(U[-1], u(x, 200.0), rtol=0, atol=1e-11)


def test_separable_solution_accuracy_and_second_order():
    def run(nx):
        xg = np.linspace(0.0, P.L, nx + 1)
        U = impl.heat_theta(xg, impl.separable_solution(xg, 0.0, P.alpha, P.L, 50.0), 300.0, nx, P.alpha)
        return impl.l2_error(U[-1], impl.separable_solution(xg, 300.0, P.alpha, P.L, 50.0), xg[1] - xg[0])
    err = impl.convergence_study(run, [25, 50, 100])
    assert err[-1] < 1e-3                                             # 50 K amplitude: relative 2e-5
    npt.assert_allclose(impl.observed_order(err), 2.0, atol=0.05)


# ---------------------------------------------------------------- manufactured solutions

def test_manufactured_source_matches_hand_formula():
    x, t = sp.symbols("x t")
    tau, beta = 150.0, 1e-3
    u_fn, q_fn = impl.manufactured_source(30 * sp.sin(sp.pi * x / P.L) * sp.exp(-t / tau), x, t, P.alpha, beta)
    xg = np.linspace(0.0, P.L, 11)
    u = 30 * np.sin(np.pi * xg / P.L) * math.exp(-40.0 / tau)
    q_hand = (-1.0 / tau + P.alpha * (math.pi / P.L) ** 2 + beta) * u
    npt.assert_allclose(u_fn(xg, 40.0), u, rtol=1e-12)
    npt.assert_allclose(q_fn(xg, 40.0), q_hand, rtol=1e-10, atol=1e-14)
    assert q_fn(xg, 40.0).shape == xg.shape and abs(u_fn(np.array([0.0, P.L]), 40.0)).max() < 1e-12


def test_observed_order_crank_nicolson_is_two():
    # dt = t_end/nx so both discretisation errors shrink together: the observed order must be 2 +- 0.1
    err = impl.convergence_study(_mms_error, [25, 50, 100, 200])
    orders = impl.observed_order(err)
    assert np.all(np.diff(err) < 0)
    npt.assert_allclose(orders, 2.0, atol=0.1)


def test_backward_euler_is_first_order():
    err = impl.convergence_study(lambda n: _mms_error(n, theta=1.0), [25, 50, 100])
    npt.assert_allclose(impl.observed_order(err), 1.0, atol=0.1)
    assert err[-1] > 10 * _mms_error(100)                           # and far less accurate than CN


def test_planted_bug_source_time_drops_the_order_to_one():
    err = impl.convergence_study(lambda n: _mms_error(n, bug="source_time"), [25, 50, 100])
    orders = impl.observed_order(err)
    assert np.all((orders > 0.9) & (orders < 1.3))


def test_planted_bug_diffusivity_converges_to_the_wrong_answer():
    err = impl.convergence_study(lambda n: _mms_error(n, bug="diffusivity"), [25, 50, 100])
    assert err[-1] > 0.1                                              # error plateau (0.15 K here) ...
    assert abs(impl.observed_order(err)[-1]) < 0.3                    # ... with observed order ~ 0


def test_unknown_bug_flag_raises():
    x = np.linspace(0.0, P.L, 5)
    with pytest.raises(ValueError):
        impl.heat_theta(x, np.zeros(5), 1.0, 1, P.alpha, bug="typo")


def test_temporal_order_from_point_values_without_an_exact_solution():
    # fixed fine grid, time step halved twice; Richardson on the point value gives the time order even
    # though the spatial error floor would spoil a comparison with the exact solution.
    u_m, q_m = _mms()
    nx = 200
    xg = np.linspace(0.0, P.L, nx + 1)
    vals = {}
    for theta in (0.5, 1.0):
        v = [float(impl.heat_theta(xg, u_m(xg, 0.0), 300.0, nt, P.alpha, source=q_m, theta=theta)[-1][60])
             for nt in (10, 20, 40)]
        vals[theta] = impl.richardson(v[2], v[1], v[0])["p"]
    assert abs(vals[0.5] - 2.0) < 0.1 and abs(vals[1.0] - 1.0) < 0.1


# ---------------------------------------------------------------- order and Richardson tools

def test_observed_order_recovers_the_exponent_of_an_exact_error_law():
    h = 0.4 / 2.0 ** np.arange(5)
    npt.assert_allclose(impl.observed_order(3.0 * h ** 1.7), 1.7, rtol=1e-12)
    npt.assert_allclose(impl.observed_order(3.0 * h ** 2, r=2.0), 2.0, rtol=1e-12)


def test_richardson_on_an_exact_sequence():
    f = lambda h: 10.0 + 3.0 * h ** 2
    res = impl.richardson(f(0.1), f(0.2), f(0.4))
    npt.assert_allclose(res["p"], 2.0, rtol=1e-10)
    npt.assert_allclose(res["f_extrap"], 10.0, rtol=1e-10)
    npt.assert_allclose(res["gci_fine"], 1.25 * abs(f(0.2) - f(0.1)) / f(0.1) / 3.0, rtol=1e-10)
    npt.assert_allclose(res["asymptotic_ratio"], f(0.1) / f(0.2), rtol=1e-10)   # tautology with the observed p
    res2 = impl.richardson(f(0.1), f(0.2), f(0.4), p_formal=2.0)
    npt.assert_allclose(res2["f_extrap"], 10.0, rtol=1e-10)
    # a sequence converging at order 1 judged with a formal order 2: the ratio exposes the mismatch
    g = lambda h: 10.0 + 3.0 * h
    assert abs(impl.richardson(g(0.1), g(0.2), g(0.4), p_formal=2.0)["asymptotic_ratio"] - 1.0) > 0.4


def test_richardson_raises_on_oscillatory_convergence():
    with pytest.raises(ValueError):
        impl.richardson(1.0, 1.1, 1.05)


def test_richardson_on_the_heater_functional_bounds_the_true_error():
    def at_heater(nx):
        xg = np.linspace(0.0, P.L, nx + 1)
        U = impl.heat_theta(xg, np.zeros(nx + 1), P.t_end, nx, P.alpha, source=impl.heater_source(P))
        return float(U[-1][int(round(nx * P.x0 / P.L))])
    f3, f2, f1, f_ref = at_heater(25), at_heater(50), at_heater(100), at_heater(400)
    res = impl.richardson(f1, f2, f3)
    assert 1.9 < res["p"] < 2.2
    assert abs(res["f_extrap"] - f_ref) < 0.2 * abs(f1 - f_ref)      # extrapolation beats the fine grid
    assert abs(f1 - f_ref) / abs(f1) <= res["gci_fine"]               # the GCI (Fs = 1.25) covers the error
    assert 1.5 < impl.required_refinement(res["gci_fine"], res["p"], 1e-4) < 2.0


def test_l2_error_is_the_discrete_l2_norm():
    assert impl.l2_error(np.zeros(5), np.zeros(5), 0.1) == 0.0
    npt.assert_allclose(impl.l2_error(np.ones(11), np.zeros(11), 0.1), math.sqrt(1.1), rtol=1e-12)


# ---------------------------------------------------------------- validation (synthetic data)

def _model_factory(nx=100, nt=300):
    xv = np.linspace(0.0, P.L, nx + 1)
    tv = np.linspace(0.0, P.t_end, nt + 1)

    def model(beta):
        U = impl.heat_theta(xv, np.zeros(nx + 1), P.t_end, nt, P.alpha, beta=beta, source=impl.heater_source(P))
        return impl.sample_solution(xv, tv, U, P.x_sensors, P.t_sensors)
    return model


def test_synthetic_data_is_labelled_and_has_the_stated_noise():
    meas = impl.synthetic_measurements(P, nx=100)
    assert meas["source"] == "synthetic"
    rms = math.sqrt(np.mean((meas["data"] - meas["truth"]) ** 2))
    assert 0.2 < rms < 0.4                                            # 24 draws of sd 0.3: sd of the rms ~ 0.045


def test_validation_metric_flags_the_missing_loss_term_and_passes_the_true_model():
    meas = impl.synthetic_measurements(P, nx=100)
    model = _model_factory()
    bad = impl.validation_comparison(model(0.0), meas["data"], meas["u_data"], 0.02)
    good = impl.validation_comparison(model(P.beta_true), meas["data"], meas["u_data"], 0.02)
    assert bad["n_within_2"] <= 10 and bad["chi2_reduced"] > 10 and bad["max_abs_z"] > 8
    assert np.all(bad["E"][-1] > 2.0)                                 # the model overpredicts late temperatures
    assert good["n_within_2"] >= 20 and 0.5 < good["chi2_reduced"] < 2.0
    npt.assert_allclose(good["u_val"], math.sqrt(0.3 ** 2 + 0.02 ** 2))


def test_fit_recovers_the_loss_coefficient_within_its_uncertainty():
    meas = impl.synthetic_measurements(P, nx=100)
    fit = impl.fit_loss_coefficient(_model_factory(), meas["data"], meas["u_data"])
    assert abs(fit["beta"] - P.beta_true) < 0.1 * P.beta_true         # 24 points at 0.3 K give ~ 3 % on beta
    assert fit["chi2"] < 2.0 * meas["data"].size and fit["n_eval"] < 40
