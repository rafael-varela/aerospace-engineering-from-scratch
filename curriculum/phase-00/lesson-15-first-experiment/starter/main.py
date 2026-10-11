"""P00.15 First Experiment: Predict, Measure, Compare: starter code.

Predict the fall time of a dropped object with its uncertainty, reduce repeated timings to a
mean and a Type A uncertainty, compare prediction and measurement with a normalized error,
and plot both with error bars.

Complete every function marked NotImplementedError. When finished, running this file prints
the comparison table for the illustrative example and saves the plot.

Run:  python main.py [output.png]
Test: pytest curriculum/phase-00/lesson-15-first-experiment
"""
from __future__ import annotations

import math
import sys
from collections.abc import Sequence
from dataclasses import dataclass

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

G_STANDARD = 9.80665
"""Standard acceleration of free fall, m/s^2."""

U_G_LOCAL = 0.03
"""Standard uncertainty to use for g when you do not know your local value, m/s^2.

Local g varies between about 9.78 and 9.83 m/s^2 over the Earth's surface; treating the
true value as uniformly spread over ±0.03 m/s^2 about 9.807 gives u = 0.03/sqrt(3) = 0.017,
and 0.03 is a deliberately generous round number.
"""


def predict_fall_time(h: float, u_h: float, g: float = G_STANDARD, u_g: float = U_G_LOCAL) -> tuple[float, float]:
    """Predicted fall time t = sqrt(2 h / g) and its standard uncertainty.

    Propagation (law of propagation of uncertainty, independent inputs):
    dt/dh = t / (2 h) and dt/dg = -t / (2 g), so
    u_t = t/2 * sqrt((u_h/h)^2 + (u_g/g)^2).
    ``h`` must be positive (ValueError otherwise).
    """
    raise NotImplementedError


def fall_time_with_drag(h: float, v_terminal: float, g: float = G_STANDARD) -> float:
    """Fall time from rest through height ``h`` with quadratic air drag.

    With drag proportional to v^2 the object approaches a terminal speed ``v_terminal``
    and the distance fallen is x(t) = (v_t^2/g) ln cosh(g t / v_t). Inverting:
    t = (v_t / g) arccosh(exp(g h / v_t^2)). As v_terminal -> infinity this tends to
    sqrt(2 h / g). Pass ``math.inf`` for the drag-free case (return sqrt(2h/g) directly).
    Numerics: with a = g h / v_t^2, write arccosh(e^a) = log1p(x + sqrt(x (x + 2))) where
    x = expm1(a); this keeps full precision when a is tiny (weak drag), where acosh(exp(a))
    cancels badly. For a > 300 (exp would overflow) use arccosh(e^a) ~ a + ln 2.
    """
    raise NotImplementedError


def frame_timing(n_frames: float, fps: float) -> tuple[float, float]:
    """Fall time from a frame count and its Type B uncertainty.

    t = n_frames / fps. The release and the impact are each located to within half a frame
    with a rectangular distribution, so each end contributes (0.5/fps)/sqrt(3) and together
    u_t = sqrt(2) * 0.5 / (sqrt(3) * fps). The frame rate itself is taken as exact.
    ``fps`` must be positive (ValueError otherwise).
    """
    raise NotImplementedError


def reduce_timings(times: Sequence[float] | np.ndarray, u_single: float = 0.0) -> tuple[float, float]:
    """Mean of repeated timings and its standard uncertainty.

    The Type A part is s / sqrt(n). ``u_single`` is a Type B standard uncertainty that
    applies to every reading in the same way (for example the frame-quantisation term from
    ``frame_timing``); it does not average down and is added in quadrature.
    Needs at least two timings (ValueError otherwise).
    """
    raise NotImplementedError


def normalized_error(pred: float, u_pred: float, meas: float, u_meas: float) -> float:
    """|pred - meas| / sqrt(u_pred^2 + u_meas^2).

    Below about 2 the prediction and measurement agree within their combined uncertainty
    (k = 2). Larger values mean something is wrong with the model, the measurement, or the
    uncertainty budget of either. If both uncertainties are zero, return 0.0 when the values
    are equal and ``math.inf`` otherwise.
    """
    raise NotImplementedError


@dataclass(frozen=True)
class Comparison:
    """Prediction and measurement for each drop height."""

    heights: np.ndarray
    t_pred: np.ndarray
    u_pred: np.ndarray
    t_meas: np.ndarray
    u_meas: np.ndarray
    e_n: np.ndarray


def compare(
    heights: Sequence[float],
    u_h: float,
    timings: Sequence[Sequence[float]],
    g: float = G_STANDARD,
    u_g: float = U_G_LOCAL,
    u_single: float = 0.0,
) -> Comparison:
    """Predict, reduce, and compare for several heights.

    ``timings[i]`` holds the repeated fall times measured at ``heights[i]``
    (ValueError if the lengths differ).
    """
    raise NotImplementedError


def fit_g(heights: Sequence[float], t_meas: Sequence[float]) -> tuple[float, float]:
    """Estimate g from h = (g/2) t^2 by a least-squares line through the origin in (t^2, h).

    Returns ``(g, u_g)``; the uncertainty comes from the residual scatter and needs at
    least three heights (ValueError otherwise).
    Slope through the origin: b = sum(t2*h) / sum(t2*t2); residual variance s^2 = sum(r^2)/(n-1);
    u_b = sqrt(s^2 / sum(t2*t2)); g = 2 b, u_g = 2 u_b.
    """
    raise NotImplementedError


def make_synthetic_timings(
    heights: Sequence[float], n_rep: int, fps: float, seed: int, g: float = G_STANDARD,
    reaction_sd: float = 0.0,
) -> list[np.ndarray]:
    """Synthetic frame-counted timings that are CONSISTENT with free fall, for testing the pipeline.

    This is simulation, not data: it exists so the tests can check the comparison logic on a
    dataset whose true behaviour is known. Each timing is the free-fall time plus an optional
    normal reaction-time error (``reaction_sd``, in seconds, for stopwatch runs), quantised
    to whole frames at ``fps``: frames = floor(t * fps + phase) with phase uniform on [0, 1),
    returned as frames / fps. Use ``numpy.random.default_rng(seed)``.
    """
    raise NotImplementedError


def plot_comparison(comp: Comparison, path: str) -> None:
    """Save a figure of predicted and measured fall time against height, with error bars.

    Upper panel: the curve sqrt(2h/g), the predictions and the measurements with k = 2 error
    bars. Lower panel: E_n per height with a dashed line at 2. Save to ``path`` and close the figure.
    """
    raise NotImplementedError


# Illustrative example for main(): frame counts at 30 fps for a tennis ball, made up to show the
# output format. They are NOT a measurement; replace them with your own data/raw/ files.
EXAMPLE_HEIGHTS_M = [0.60, 0.90, 1.20, 1.50, 1.80]
EXAMPLE_U_H_M = 0.005
EXAMPLE_FPS = 30.0
EXAMPLE_FRAMES = [
    [10, 11, 10, 11, 10],
    [13, 13, 12, 13, 13],
    [15, 15, 14, 15, 15],
    [16, 17, 17, 16, 17],
    [18, 18, 19, 18, 18],
]


def main() -> None:
    """Run the comparison on the illustrative example and save the plot."""
    out = sys.argv[1] if len(sys.argv) > 1 else "fall_time_comparison.png"
    print("Illustrative example (replace with your own data from projects/p00-first-experiment/)")
    _, u_frame = frame_timing(1, EXAMPLE_FPS)
    timings = [np.asarray(f) / EXAMPLE_FPS for f in EXAMPLE_FRAMES]
    comp = compare(EXAMPLE_HEIGHTS_M, EXAMPLE_U_H_M, timings, u_single=u_frame)
    print(f"{'h (m)':>6} {'t_pred (s)':>11} {'u_pred':>7} {'t_meas (s)':>11} {'u_meas':>7} {'E_n':>5}")
    for h, tp, up, tm, um, en in zip(comp.heights, comp.t_pred, comp.u_pred, comp.t_meas, comp.u_meas, comp.e_n):
        print(f"{h:6.2f} {tp:11.3f} {up:7.3f} {tm:11.3f} {um:7.3f} {en:5.2f}")
    g_fit, u_gfit = fit_g(comp.heights, comp.t_meas)
    print(f"g from the fit of h against t^2: {g_fit:.2f} ± {u_gfit:.2f} m/s^2 (standard 9.81)")
    t_drag = fall_time_with_drag(1.80, v_terminal=22.0)
    print(f"drag check at 1.80 m for a tennis ball (v_t ~ 22 m/s): {t_drag:.4f} s vs "
          f"{predict_fall_time(1.80, 0.0, u_g=0.0)[0]:.4f} s free fall")
    plot_comparison(comp, out)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
