"""Public tests for P09.01 Signals & Systems.

They assert properties, not magic numbers: the elementary signals and their integral relations; energy
and power against closed forms; the even/odd split is exact and orthogonal; shift and scale compose in
the right order. The independent check of the property testers is the ANALYTIC verdict table below,
derived by hand for each system (a delay, a squarer, an accumulator, a time-varying gain and more): the
numerical tests must reproduce it. LTI systems are also checked against convolution with their impulse
response and against the BIBO bound |y| <= sum|h| max|x|.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
from __future__ import annotations

import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

# Hand-derived verdicts (linear, time-invariant, causal, memoryless, BIBO stable), independent of the code.
ANALYTIC = {
    "delay 2":          (True,  True,  True,  False, True),    # y[n] = x[n-2]
    "squarer":          (False, True,  True,  True,  True),    # (x1 + x2)^2 != x1^2 + x2^2
    "accumulator":      (True,  True,  True,  False, False),   # a step gives a ramp
    "gain n x[n]":      (True,  False, True,  True,  False),   # the gain depends on n and is unbounded
    "bias correction":  (False, True,  True,  True,  True),    # affine: S(0) = -b
    "saturation 1":     (False, True,  True,  True,  True),    # linear only below the limit
    "moving average 4": (True,  True,  True,  False, True),
    "low-pass a=0.9":   (True,  True,  True,  False, True),    # pole at 0.9 inside the unit circle
    "recursion a=1.05": (True,  True,  True,  False, False),   # pole at 1.05 outside
    "centred smoother": (True,  True,  False, False, True),    # uses x[n+1]
    "decimator":        (True,  False, False, False, True),    # y[n] = x[2n] reads ahead for n > 0
}


# ---------------------------------------------------------------- signals

def test_elementary_signals_and_their_relations():
    t = np.linspace(-1.0, 3.0, 4001)
    dt = t[1] - t[0]
    u = impl.unit_step(t)
    assert set(np.unique(u)) <= {0.0, 1.0} and u[t < 0].sum() == 0 and u[t >= 0].min() == 1.0
    # the ramp is the running integral of the step (rectangle rule, error <= dt)
    npt.assert_allclose(np.cumsum(u) * dt - dt, impl.ramp(t), atol=2 * dt)
    # rect(t, w) = u(t) - u(t - w); its area is w
    r = impl.rect_pulse(t, 1.5)
    npt.assert_array_equal(r, impl.unit_step(t) - impl.unit_step(t - 1.5))
    assert impl.energy(r, dt) == pytest.approx(1.5, abs=2 * dt)
    d = impl.unit_impulse(8, 3)
    npt.assert_array_equal(d, [0, 0, 0, 1, 0, 0, 0, 0])


@pytest.mark.parametrize("a", [0.5, 2.0, 10.0])
def test_energy_of_a_decaying_exponential_converges_to_1_over_2a(a):
    for dt in (1e-2, 1e-4):                                # left rectangle rule: relative error ~ a dt
        t = np.arange(0.0, 40.0 / a, dt)
        assert impl.energy(np.exp(-a * t), dt) == pytest.approx(1.0 / (2.0 * a), rel=1.2 * a * dt)


def test_power_of_a_sinusoid_is_half_the_square_of_its_amplitude_and_energy_grows_with_time():
    dt = 1e-4
    for amp, f, ph in ((2.0, 5.0, 0.0), (0.7, 13.0, 1.1), (3.0, 0.5, -0.4)):
        t = np.arange(0.0, 10.0, dt)                       # a whole number of periods
        x = impl.sinusoid(t, amp, f, ph)
        assert impl.average_power(x) == pytest.approx(amp ** 2 / 2, rel=1e-6)
    t1, t2 = np.arange(0.0, 10.0, dt), np.arange(0.0, 20.0, dt)
    e1 = impl.energy(impl.sinusoid(t1, 1.0, 5.0), dt)
    e2 = impl.energy(impl.sinusoid(t2, 1.0, 5.0), dt)
    assert e2 == pytest.approx(2 * e1, rel=1e-6)           # a power signal has infinite energy


def test_even_odd_split_is_exact_symmetric_and_orthogonal():
    t = np.linspace(-3.0, 3.0, 6001)
    dt = t[1] - t[0]
    x = np.exp(-t) * impl.unit_step(t) + 0.3 * np.sin(4 * t) + 0.1 * t ** 2
    e, o = impl.even_part(x), impl.odd_part(x)
    npt.assert_allclose(e + o, x, atol=1e-15)
    npt.assert_allclose(e, e[::-1], atol=1e-15)
    npt.assert_allclose(o, -o[::-1], atol=1e-15)
    assert np.sum(e * o) * dt == pytest.approx(0.0, abs=1e-12)
    assert impl.energy(x, dt) == pytest.approx(impl.energy(e, dt) + impl.energy(o, dt), rel=1e-12)
    npt.assert_allclose(impl.even_part(np.cos(t)), np.cos(t), atol=1e-15)
    npt.assert_allclose(impl.odd_part(np.cos(t)), 0.0, atol=1e-15)


def test_shift_and_scale_compose_in_the_right_order():
    f = lambda t: impl.rect_pulse(t, 1.0) * (1.0 + t)      # an asymmetric pulse on [0, 1)
    t = np.linspace(-2.0, 5.0, 7001)
    target = f(2.0 * t - 3.0)                              # x(2t - 3) occupies [1.5, 2)
    npt.assert_array_equal(impl.time_scale(impl.time_shift(f, 3.0), 2.0)(t), target)   # shift by 3, then scale
    npt.assert_array_equal(impl.time_shift(impl.time_scale(f, 2.0), 1.5)(t), target)   # scale, then shift by 3/2
    assert not np.array_equal(impl.time_shift(impl.time_scale(f, 2.0), 3.0)(t), target)
    npt.assert_array_equal(impl.time_reverse(impl.time_reverse(f))(t), f(t))
    g = impl.time_shift(f, 0.5)(t)
    assert t[np.argmax(g > 0)] == pytest.approx(0.5, abs=1e-3)   # a positive t0 delays the signal


# ---------------------------------------------------------------- property tests: the analytic verdicts

@pytest.mark.parametrize("name", list(ANALYTIC))
def test_property_tests_reproduce_the_analytic_verdict(name):
    S = impl.gyro_path_zoo()[name]
    v = impl.classify(S)
    got = tuple(v[p] for p in ("linear", "time_invariant", "causal", "memoryless", "bibo_stable"))
    assert got == ANALYTIC[name], f"{name}: got {got}, analytic {ANALYTIC[name]}"


def test_saturation_is_linear_in_its_linear_region_and_not_beyond():
    sat = impl.saturation(1.0)
    rng = np.random.default_rng(3)
    x1, x2 = rng.uniform(-1, 1, 300), rng.uniform(-1, 1, 300)
    # |0.4 x1 + 0.5 x2| <= 0.9 < 1: every signal stays inside the limits, superposition holds exactly
    assert impl.superposition_error(sat, x1, x2, 0.4, 0.5) < 1e-15
    assert impl.superposition_error(sat, 5 * x1, 5 * x2, 0.4, 0.5) > 0.05
    assert impl.is_linear(sat, amplitude=0.01) is True     # the falsification test cannot see the limit
    assert impl.is_linear(sat, amplitude=10.0) is False


def test_linear_systems_map_zero_to_zero_and_affine_ones_do_not():
    z = np.zeros(64)
    for name in ("delay 2", "accumulator", "moving average 4", "low-pass a=0.9", "centred smoother"):
        npt.assert_array_equal(impl.gyro_path_zoo()[name](z), 0.0)
    assert np.all(impl.gyro_path_zoo()["bias correction"](z) == -0.05)


# ---------------------------------------------------------------- LTI systems and the impulse response

@pytest.mark.parametrize("a,b", [(0.9, 0.1), (0.5, 1.0), (-0.7, 2.0)])
def test_impulse_response_of_the_recursion_is_geometric(a, b):
    h = impl.impulse_response(impl.first_order_recursion(a, b), 60)
    npt.assert_allclose(h, b * a ** np.arange(60), rtol=1e-12)
    assert impl.abs_sum(h) == pytest.approx(abs(b) * (1 - abs(a) ** 60) / (1 - abs(a)), rel=1e-12)


def test_causal_lti_output_is_the_convolution_with_the_impulse_response():
    rng = np.random.default_rng(7)
    x = rng.standard_normal(300)
    for name in ("delay 2", "accumulator", "moving average 4", "low-pass a=0.9"):
        S = impl.gyro_path_zoo()[name]
        h = impl.impulse_response(S, 300)
        npt.assert_allclose(S(x), np.convolve(x, h)[:300], atol=1e-10)


def test_bibo_bound_holds_and_is_reached_by_the_worst_case_input():
    S = impl.first_order_recursion(-0.8, 1.0)              # alternating impulse response
    n = 400
    h = impl.impulse_response(S, n)
    bound = impl.abs_sum(h)
    assert bound == pytest.approx(1 / (1 - 0.8), rel=1e-9)  # sum |(-0.8)^k| = 5
    rng = np.random.default_rng(11)
    for _ in range(5):
        x = rng.uniform(-1, 1, n)
        assert np.max(np.abs(S(x))) <= bound * np.max(np.abs(x)) + 1e-12
    x_worst = np.sign(h[::-1])                             # x[m] = sign h[n - 1 - m] makes y[n-1] = sum |h|
    assert S(x_worst)[-1] == pytest.approx(bound, rel=1e-9)


def test_accumulator_impulse_response_is_not_absolutely_summable():
    acc = impl.accumulator()
    sums = [impl.abs_sum(impl.impulse_response(acc, n)) for n in (100, 1000, 10000)]
    assert sums == pytest.approx([100.0, 1000.0, 10000.0])
    # and an integrated gyro bias grows without bound: 0.01 deg/s for one hour at 100 Hz -> 36 deg
    angle = impl.accumulator(0.01)(np.full(360000, 0.01))
    assert angle[-1] == pytest.approx(36.0, rel=1e-9)


def test_stability_verdict_needs_a_record_longer_than_the_time_constant():
    lp = impl.first_order_recursion(0.999, 0.001)          # time constant 1000 samples
    assert impl.is_bibo_stable(lp, n=50) is False          # a false alarm on a short record
    assert impl.is_bibo_stable(lp, n=5000) is True
    assert impl.is_bibo_stable(impl.first_order_recursion(1.001, 1.0), n=5000) is False


# ---------------------------------------------------------------- figure

def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "signals_and_systems.png").stat().st_size > 1000
