"""Public tests for P00.15. They assert the formulas, their limiting cases, and that a synthetic
dataset consistent with free fall is judged consistent.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. The only file written goes to tmp_path.
"""
import math

import numpy as np
import pytest

from conftest import load_impl

impl = load_impl(__file__)

G = impl.G_STANDARD


# ---------- prediction ----------

def test_fall_time_formula():
    t, _ = impl.predict_fall_time(1.2, 0.0, g=9.81, u_g=0.0)
    assert t == pytest.approx(math.sqrt(2 * 1.2 / 9.81))


def test_fall_time_scales_as_sqrt_height():
    t1, _ = impl.predict_fall_time(0.5, 0.0, u_g=0.0)
    t4, _ = impl.predict_fall_time(2.0, 0.0, u_g=0.0)
    assert t4 / t1 == pytest.approx(2.0)


def test_zero_input_uncertainty_gives_zero_prediction_uncertainty():
    _, u = impl.predict_fall_time(1.0, 0.0, u_g=0.0)
    assert u == 0.0


def test_prediction_uncertainty_matches_finite_difference():
    h, uh, g, ug = 1.5, 0.01, 9.81, 0.02
    t, u = impl.predict_fall_time(h, uh, g, ug)
    f = lambda hh, gg: math.sqrt(2 * hh / gg)
    eps = 1e-6
    dh = (f(h + eps, g) - f(h - eps, g)) / (2 * eps)
    dg = (f(h, g + eps) - f(h, g - eps)) / (2 * eps)
    assert u == pytest.approx(math.hypot(dh * uh, dg * ug), rel=1e-6)


def test_prediction_uncertainty_is_linear_in_height_uncertainty():
    _, u1 = impl.predict_fall_time(1.0, 0.01, u_g=0.0)
    _, u2 = impl.predict_fall_time(1.0, 0.02, u_g=0.0)
    assert u2 == pytest.approx(2 * u1)


def test_predict_rejects_nonpositive_height():
    with pytest.raises(ValueError):
        impl.predict_fall_time(0.0, 0.01)


# ---------- drag as a limiting case ----------

def test_drag_free_limit_recovers_free_fall():
    t0 = math.sqrt(2 * 1.8 / G)
    assert impl.fall_time_with_drag(1.8, math.inf) == pytest.approx(t0)
    assert impl.fall_time_with_drag(1.8, 1e6) == pytest.approx(t0, rel=1e-9)


def test_drag_lengthens_fall_and_more_so_for_lower_terminal_speed():
    t0 = math.sqrt(2 * 1.8 / G)
    t_fast = impl.fall_time_with_drag(1.8, 50.0)
    t_slow = impl.fall_time_with_drag(1.8, 10.0)
    assert t0 < t_fast < t_slow


def test_drag_correction_small_for_tennis_ball_heights():
    # v_t ~ 22 m/s; the correction at 2 m should be a few percent at most
    t0 = math.sqrt(2 * 2.0 / G)
    assert 1.0 < impl.fall_time_with_drag(2.0, 22.0) / t0 < 1.05


def test_drag_solution_satisfies_the_distance_law():
    v_t, h = 15.0, 1.5
    t = impl.fall_time_with_drag(h, v_t)
    x = (v_t**2 / G) * math.log(math.cosh(G * t / v_t))
    assert x == pytest.approx(h, rel=1e-9)


# ---------- measurement reduction ----------

def test_frame_timing_value_and_uncertainty_scale_with_fps():
    t30, u30 = impl.frame_timing(15, 30.0)
    t60, u60 = impl.frame_timing(30, 60.0)
    assert t30 == pytest.approx(0.5) and t60 == pytest.approx(0.5)
    assert u30 == pytest.approx(2 * u60)
    assert u30 == pytest.approx(math.sqrt(2) * 0.5 / (math.sqrt(3) * 30.0))
    with pytest.raises(ValueError):
        impl.frame_timing(10, 0.0)


def test_reduce_timings_type_a():
    mean, u = impl.reduce_timings([0.50, 0.52, 0.48])
    assert mean == pytest.approx(0.50)
    assert u == pytest.approx(np.std([0.50, 0.52, 0.48], ddof=1) / math.sqrt(3))


def test_reduce_timings_adds_type_b_in_quadrature_and_it_does_not_average_down():
    times = np.full(50, 0.6)
    mean, u = impl.reduce_timings(times, u_single=0.01)
    assert mean == pytest.approx(0.6)
    assert u == pytest.approx(0.01)
    with pytest.raises(ValueError):
        impl.reduce_timings([0.6])


# ---------- comparison metric ----------

def test_normalized_error_properties():
    assert impl.normalized_error(1.0, 0.1, 1.0, 0.1) == 0.0
    assert impl.normalized_error(1.0, 0.03, 1.05, 0.04) == pytest.approx(1.0)
    assert impl.normalized_error(1.0, 0.03, 1.05, 0.04) == impl.normalized_error(1.05, 0.04, 1.0, 0.03)
    assert impl.normalized_error(1.0, 0.0, 1.0, 0.0) == 0.0
    assert math.isinf(impl.normalized_error(1.0, 0.0, 1.1, 0.0))


def test_normalized_error_grows_when_uncertainty_shrinks():
    loose = impl.normalized_error(1.0, 0.05, 1.1, 0.05)
    tight = impl.normalized_error(1.0, 0.01, 1.1, 0.01)
    assert tight > loose


# ---------- end to end on a consistent synthetic dataset ----------

HEIGHTS = [0.6, 0.9, 1.2, 1.5, 1.8, 2.0]


@pytest.mark.parametrize("seed", [0, 1, 2, 3])
def test_consistent_synthetic_video_dataset_has_small_normalized_error(seed):
    fps = 30.0
    timings = impl.make_synthetic_timings(HEIGHTS, n_rep=6, fps=fps, seed=seed)
    _, u_frame = impl.frame_timing(1, fps)
    comp = impl.compare(HEIGHTS, 0.005, timings, u_single=u_frame)
    assert comp.e_n.shape == (len(HEIGHTS),)
    assert np.all(comp.e_n < 2.0), comp.e_n
    assert np.all(comp.u_meas > 0) and np.all(comp.u_pred > 0)
    assert np.all(np.diff(comp.t_pred) > 0) and np.all(np.diff(comp.t_meas) > 0)


def test_biased_dataset_is_flagged():
    # A stopwatch started late by 0.15 s every time: wrong by far more than the budget.
    fps = 1e6
    timings = [np.asarray(t) - 0.15 for t in impl.make_synthetic_timings(HEIGHTS, 5, fps, seed=7, reaction_sd=0.02)]
    comp = impl.compare(HEIGHTS, 0.005, timings)
    assert np.all(comp.e_n > 2.0)


def test_compare_rejects_mismatched_lengths():
    with pytest.raises(ValueError):
        impl.compare([1.0, 2.0], 0.01, [[0.4, 0.5]])


def test_fit_g_recovers_g_from_consistent_data():
    timings = impl.make_synthetic_timings(HEIGHTS, n_rep=8, fps=60.0, seed=11)
    t_meas = [np.mean(t) for t in timings]
    g_fit, u_g = impl.fit_g(HEIGHTS, t_meas)
    assert abs(g_fit - G) < 0.05 * G
    assert abs(g_fit - G) < 4 * u_g + 0.02 * G     # within scatter plus frame quantisation
    with pytest.raises(ValueError):
        impl.fit_g([1.0, 2.0], [0.45, 0.64])


def test_plot_is_written(tmp_path):
    timings = impl.make_synthetic_timings(HEIGHTS[:3], 4, 30.0, seed=5)
    comp = impl.compare(HEIGHTS[:3], 0.005, timings, u_single=impl.frame_timing(1, 30.0)[1])
    out = tmp_path / "cmp.png"
    impl.plot_comparison(comp, str(out))
    assert out.exists() and out.stat().st_size > 1000
