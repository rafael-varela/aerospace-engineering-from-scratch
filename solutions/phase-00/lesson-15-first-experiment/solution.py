"""P00.15 First Experiment: Predict, Measure, Compare: reference solution.

Predict the fall time of a dropped object with its uncertainty, reduce repeated timings to a
mean and a Type A uncertainty, compare prediction and measurement with a normalized error,
and plot both with error bars.

Run:  python solution.py [output.png]
Test: AEFS_IMPL=solution pytest curriculum/phase-00/lesson-15-first-experiment
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
    ``h`` must be positive.
    """
    if h <= 0:
        raise ValueError("height must be positive")
    t = math.sqrt(2.0 * h / g)
    u_t = 0.5 * t * math.sqrt((u_h / h) ** 2 + (u_g / g) ** 2)
    return t, u_t


def fall_time_with_drag(h: float, v_terminal: float, g: float = G_STANDARD) -> float:
    """Fall time from rest through height ``h`` with quadratic air drag.

    With drag proportional to v^2 the object approaches a terminal speed ``v_terminal``
    and the distance fallen is x(t) = (v_t^2/g) ln cosh(g t / v_t). Inverting:
    t = (v_t / g) arccosh(exp(g h / v_t^2)). As v_terminal -> infinity this tends to
    sqrt(2 h / g). Pass ``math.inf`` for the drag-free case.
    """
    if h <= 0:
        raise ValueError("height must be positive")
    if math.isinf(v_terminal):
        return math.sqrt(2.0 * h / g)
    a = g * h / v_terminal**2
    if a > 300:  # exp would overflow; arccosh(e^a) ~ a + ln 2 for large a
        return (v_terminal / g) * (a + math.log(2.0))
    # arccosh(e^a) = ln(e^a + sqrt(e^2a - 1)). Written with x = e^a - 1 = expm1(a) and log1p it
    # keeps full precision when a is tiny (weak drag), where acosh(exp(a)) would cancel badly.
    x = math.expm1(a)
    return (v_terminal / g) * math.log1p(x + math.sqrt(x * (x + 2.0)))


def frame_timing(n_frames: float, fps: float) -> tuple[float, float]:
    """Fall time from a frame count and its Type B uncertainty.

    t = n_frames / fps. The release and the impact are each located to within half a frame
    with a rectangular distribution, so each end contributes (0.5/fps)/sqrt(3) and together
    u_t = sqrt(2) * 0.5 / (sqrt(3) * fps). The frame rate itself is taken as exact.
    """
    if fps <= 0:
        raise ValueError("frame rate must be positive")
    t = n_frames / fps
    u_t = math.sqrt(2.0) * 0.5 / (math.sqrt(3.0) * fps)
    return t, u_t


def reduce_timings(times: Sequence[float] | np.ndarray, u_single: float = 0.0) -> tuple[float, float]:
    """Mean of repeated timings and its standard uncertainty.

    The Type A part is s / sqrt(n). ``u_single`` is a Type B standard uncertainty that
    applies to every reading in the same way (for example the frame-quantisation term from
    ``frame_timing``); it does not average down and is added in quadrature.
    """
    x = np.asarray(times, dtype=float).ravel()
    n = x.size
    if n < 2:
        raise ValueError("need at least two timings for a Type A evaluation")
    mean = float(x.mean())
    u_a = float(x.std(ddof=1) / math.sqrt(n))
    return mean, math.hypot(u_a, u_single)


def normalized_error(pred: float, u_pred: float, meas: float, u_meas: float) -> float:
    """|pred - meas| / sqrt(u_pred^2 + u_meas^2).

    Below about 2 the prediction and measurement agree within their combined uncertainty
    (k = 2). Larger values mean something is wrong with the model, the measurement, or the
    uncertainty budget of either.
    """
    denom = math.hypot(u_pred, u_meas)
    if denom == 0.0:
        return 0.0 if pred == meas else math.inf
    return abs(pred - meas) / denom


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

    ``timings[i]`` holds the repeated fall times measured at ``heights[i]``.
    """
    hs = np.asarray(heights, dtype=float)
    if len(timings) != hs.size:
        raise ValueError("one list of timings per height")
    tp, up, tm, um, en = [], [], [], [], []
    for h, ts in zip(hs, timings):
        t_p, u_p = predict_fall_time(h, u_h, g, u_g)
        t_m, u_m = reduce_timings(ts, u_single)
        tp.append(t_p); up.append(u_p); tm.append(t_m); um.append(u_m)
        en.append(normalized_error(t_p, u_p, t_m, u_m))
    return Comparison(hs, np.array(tp), np.array(up), np.array(tm), np.array(um), np.array(en))


def fit_g(heights: Sequence[float], t_meas: Sequence[float]) -> tuple[float, float]:
    """Estimate g from h = (g/2) t^2 by a least-squares line through the origin in (t^2, h).

    Returns ``(g, u_g)``; the uncertainty comes from the residual scatter and needs at
    least three heights.
    """
    hs = np.asarray(heights, dtype=float)
    t2 = np.asarray(t_meas, dtype=float) ** 2
    if hs.size < 3 or t2.size != hs.size:
        raise ValueError("need at least three (height, time) pairs")
    slope = float(np.sum(t2 * hs) / np.sum(t2 * t2))
    resid = hs - slope * t2
    s2 = float(np.sum(resid**2) / (hs.size - 1))
    u_slope = math.sqrt(s2 / float(np.sum(t2 * t2)))
    return 2.0 * slope, 2.0 * u_slope


def make_synthetic_timings(
    heights: Sequence[float], n_rep: int, fps: float, seed: int, g: float = G_STANDARD,
    reaction_sd: float = 0.0,
) -> list[np.ndarray]:
    """Synthetic frame-counted timings that are CONSISTENT with free fall, for testing the pipeline.

    This is simulation, not data: it exists so the tests can check the comparison logic on a
    dataset whose true behaviour is known. Each timing is the free-fall time plus an optional
    normal reaction-time error (``reaction_sd``, in seconds, for stopwatch runs), quantised
    to whole frames at ``fps``.
    """
    rng = np.random.default_rng(seed)
    out = []
    for h in heights:
        t_true = math.sqrt(2.0 * h / g)
        t = t_true + rng.normal(0.0, reaction_sd, n_rep) if reaction_sd > 0 else np.full(n_rep, t_true)
        phase = rng.uniform(0.0, 1.0, n_rep)          # where in a frame the release happened
        frames = np.floor(t * fps + phase)
        out.append(frames / fps)
    return out


def plot_comparison(comp: Comparison, path: str) -> None:
    """Save a figure of predicted and measured fall time against height, with error bars."""
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(6, 6), sharex=True, height_ratios=[3, 1])
    hs = np.linspace(0.8 * comp.heights.min(), 1.1 * comp.heights.max(), 100)
    ax.plot(hs, np.sqrt(2 * hs / G_STANDARD), "-", color="0.6", label="t = sqrt(2h/g)")
    ax.errorbar(comp.heights, comp.t_pred, yerr=2 * comp.u_pred, fmt="s", capsize=3, label="prediction (k = 2)")
    ax.errorbar(comp.heights, comp.t_meas, yerr=2 * comp.u_meas, fmt="o", capsize=3, label="measurement (k = 2)")
    ax.set_ylabel("fall time (s)")
    ax.legend()
    ax.set_title("Falling-object timing: prediction versus measurement")
    ax2.axhline(2.0, color="0.5", ls="--", lw=0.8)
    ax2.plot(comp.heights, comp.e_n, "o")
    ax2.set_ylim(bottom=0)
    ax2.set_xlabel("drop height (m)")
    ax2.set_ylabel("E_n")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


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
