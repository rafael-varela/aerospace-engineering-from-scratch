"""Public tests for P06.03 Work-Energy & Impulse-Momentum. They assert mechanics, not magic numbers: the
work-energy solution u = v^2(x) solves its ODE and has the right limits (brake, spring, damper); the
stopping distance and the peak force (at engagement or at the stop) come from work-energy alone;
the stopping time comes from integral dx / v; the constant-force design really is constant and has
efficiency 1; and the independent check: `solve_ivp` reproduces the stopping distance, time and peak,
its force impulse equals m v0, and its energy audit closes to the kinetic energy.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

M, V0 = 15000.0, 60.0
CASES = [  # (F0, k, c): mixed, damper + brake, spring + damper, spring only, brake only
    (5.0e4, 3000.0, 80.0),
    (1.0e5, 0.0, 150.0),
    (0.0, 6000.0, 60.0),
    (0.0, 8000.0, 0.0),
    (2.0e5, 0.0, 0.0),
]


@pytest.mark.parametrize("F0,k,c", CASES)
def test_work_energy_solution_solves_its_ode(F0, k, c):
    xs = impl.stopping_distance(V0, M, F0, k, c)
    x = np.linspace(0.0, 0.95 * xs, 50)
    h = 1e-4 * xs
    u = impl.v_squared(x, V0, M, F0, k, c)
    du = (impl.v_squared(x + h, V0, M, F0, k, c) - impl.v_squared(x - h, V0, M, F0, k, c)) / (2 * h)
    # d(m u / 2)/dx = -F  ->  u' = -(2/m)(F0 + k x + c u)
    npt.assert_allclose(du, -2.0 / M * (F0 + k * x + c * u), rtol=1e-6, atol=1e-8)
    assert float(impl.v_squared(0.0, V0, M, F0, k, c)) == pytest.approx(V0 ** 2)
    assert abs(float(impl.v_squared(xs, V0, M, F0, k, c))) < 1e-8 * V0 ** 2


def test_limits_brake_spring_and_damper():
    F0 = 3.0e5
    assert impl.stopping_distance(V0, M, F0, 0.0, 0.0) == pytest.approx(M * V0 ** 2 / (2 * F0))
    k = 9000.0
    assert impl.stopping_distance(V0, M, 0.0, k, 0.0) == pytest.approx(V0 * np.sqrt(M / k))
    # a pure v^2 damper never stops the aircraft: v = v0 exp(-c x / m)
    assert np.isinf(impl.stopping_distance(V0, M, 0.0, 0.0, 100.0))
    x = np.array([0.0, 50.0, 200.0])
    npt.assert_allclose(np.sqrt(impl.v_squared(x, V0, M, 0.0, 0.0, 100.0)), V0 * np.exp(-100.0 * x / M), rtol=1e-12)
    assert np.isinf(impl.stopping_time(V0, M, 0.0, 0.0, 100.0))


def test_small_damper_is_continuous_with_no_damper():
    # the closed form must not lose accuracy as c -> 0 (catastrophic cancellation in the naive form)
    k = 6000.0
    x0 = impl.stopping_distance(V0, M, 0.0, k, 0.0)
    for c in (1e-9, 1e-6, 1e-3):
        assert impl.stopping_distance(V0, M, 0.0, k, c) == pytest.approx(x0, rel=1e-5)


def test_more_force_stops_sooner():
    base = impl.stopping_distance(V0, M, 5e4, 3000.0, 80.0)
    assert impl.stopping_distance(V0, M, 1e5, 3000.0, 80.0) < base
    assert impl.stopping_distance(V0, M, 5e4, 6000.0, 80.0) < base
    assert impl.stopping_distance(V0, M, 5e4, 3000.0, 160.0) < base
    assert impl.stopping_distance(V0, 2 * M, 5e4, 3000.0, 80.0) > base      # heavier runs further


@pytest.mark.parametrize("F0,k,c", CASES)
def test_peak_force_is_at_an_end_of_the_stroke(F0, k, c):
    xs = impl.stopping_distance(V0, M, F0, k, c)
    x = np.linspace(0.0, xs, 2001)
    F = impl.force_along_stroke(x, V0, M, F0, k, c)
    d = np.diff(F)
    assert np.all(d >= -1e-6 * F.max()) or np.all(d <= 1e-6 * F.max())        # monotonic
    Fp, xp = impl.peak_force(V0, M, F0, k, c)
    assert Fp == pytest.approx(F.max(), rel=1e-9)
    assert xp in (0.0, pytest.approx(xs))
    # never below either average
    assert Fp >= impl.energy_average_force(V0, M, xs) * (1 - 1e-12)
    assert Fp >= impl.momentum_average_force(V0, M, impl.stopping_time(V0, M, F0, k, c)) * (1 - 1e-12)


def test_constant_force_design():
    c, k = impl.constant_force_design(V0, M, 90.0)
    assert impl.stopping_distance(V0, M, 0.0, k, c) == pytest.approx(90.0, rel=1e-10)
    x = np.linspace(0.0, 90.0, 101)
    F = impl.force_along_stroke(x, V0, M, 0.0, k, c)
    npt.assert_allclose(F, M * V0 ** 2 / (2 * 90.0), rtol=1e-9)
    Fp, _ = impl.peak_force(V0, M, 0.0, k, c)
    assert impl.efficiency(V0, M, 90.0, Fp) == pytest.approx(1.0, rel=1e-9)
    # with a constant force, time = 2 x / v0 and the two averages coincide
    t = impl.stopping_time(V0, M, 0.0, k, c)
    assert t == pytest.approx(2 * 90.0 / V0, rel=1e-8)
    assert impl.momentum_average_force(V0, M, t) == pytest.approx(impl.energy_average_force(V0, M, 90.0), rel=1e-8)
    with pytest.raises(ValueError):
        impl.constant_force_design(V0, M, 90.0, F0=1e9)


def test_spring_only_time_is_a_quarter_period_and_efficiency_one_half():
    k = 8000.0
    assert impl.stopping_time(V0, M, 0.0, k, 0.0) == pytest.approx(0.5 * np.pi * np.sqrt(M / k), rel=1e-9)
    xs = impl.stopping_distance(V0, M, 0.0, k, 0.0)
    Fp, xp = impl.peak_force(V0, M, 0.0, k, 0.0)
    assert xp == pytest.approx(xs) and impl.efficiency(V0, M, xs, Fp) == pytest.approx(0.5)


# ---------------------------------------------------------------- independent check: time integration

@pytest.mark.parametrize("F0,k,c", CASES)
def test_time_integration_agrees_with_both_routes(F0, k, c):
    sol = impl.simulate_arrestment(V0, M, F0, k, c)
    assert sol.status == 1                                                    # stopped by the event
    xs, ts = impl.stopping_distance(V0, M, F0, k, c), impl.stopping_time(V0, M, F0, k, c)
    assert sol.y_events[0][0][0] == pytest.approx(xs, rel=1e-7)
    assert sol.t_events[0][0] == pytest.approx(ts, rel=1e-7)
    au = impl.energy_audit(sol, M, F0, k, c)
    assert au["impulse"] == pytest.approx(M * V0, rel=1e-5)                   # impulse-momentum
    assert au["W_F0"] + au["E_k"] + au["W_damper"] == pytest.approx(au["KE0"], rel=1e-5)   # work-energy
    assert au["F_peak"] == pytest.approx(impl.peak_force(V0, M, F0, k, c)[0], rel=1e-6)


def test_speed_along_the_stroke_matches_simulation():
    F0, k, c = CASES[0]
    sol = impl.simulate_arrestment(V0, M, F0, k, c)
    t = np.linspace(0.0, 0.9 * sol.t[-1], 200)
    x, v = sol.sol(t)
    npt.assert_allclose(v ** 2, impl.v_squared(x, V0, M, F0, k, c), rtol=1e-7, atol=1e-7)


def test_arrestor_force_vectorised():
    F = impl.arrestor_force(np.array([0.0, 10.0]), np.array([50.0, 40.0]), 1e4, 2e3, 50.0)
    npt.assert_allclose(F, [1e4 + 50.0 * 2500, 1e4 + 2e4 + 50.0 * 1600])


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "arrestment.png").stat().st_size > 10_000
