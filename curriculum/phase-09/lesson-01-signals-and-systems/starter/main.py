"""P09.01 Signals & Systems: starter code.

A signal toolkit and a system-property tester. Signals: the unit step, ramp, rectangular pulse, discrete
impulse and sinusoid; energy and average power; even and odd parts; time shift, time scale and time
reversal as operations on callables. Systems: sampled discrete-time operators y = S(x), each a Python
callable that maps an array x[0..N-1] (with x[n] = 0 for n < 0, "initial rest") to an array y of the same
length. Five properties are tested numerically by trying to falsify them: linearity (superposition),
time invariance (shift in, shift out), causality (the past of y ignores the future of x), memory (y[n]
depends only on x[n]) and bounded-input bounded-output stability (bounded test inputs on a long and a
four-times-longer record). The verdicts are checked against analytic verdicts for known systems.

The example is the signal path of a MEMS rate gyro in a flight computer: a delay, a bias correction, a
saturating converter, an integrator to angle, moving averages, a recursive low-pass, a centred
(look-ahead) smoother and a decimator.

Conventions:
    * Discrete signals are 1-D float arrays indexed n = 0, 1, ..., N-1; a system is a callable
      ``S(x: np.ndarray) -> np.ndarray`` returning an array of the same length.
    * Continuous signals are sampled on a grid t with spacing dt; energy and power are Riemann sums.

Complete every function marked NotImplementedError.

Run:  python main.py          (prints the property table and the examples, writes signals_and_systems.png)
Test: pytest curriculum/phase-09/lesson-01-signals-and-systems
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np

System = Callable[[np.ndarray], np.ndarray]
PROPERTIES = ("linear", "time_invariant", "causal", "memoryless", "bibo_stable")


# ---------------------------------------------------------------- elementary signals

def unit_step(t: np.ndarray) -> np.ndarray:
    """Return u(t) = 1 for t >= 0 and 0 for t < 0 (the value at t = 0 is taken as 1)."""
    raise NotImplementedError


def ramp(t: np.ndarray) -> np.ndarray:
    """Return r(t) = t u(t), the running integral of the unit step."""
    raise NotImplementedError


def rect_pulse(t: np.ndarray, width: float) -> np.ndarray:
    """Return a unit pulse that is 1 for 0 <= t < width and 0 elsewhere, i.e. u(t) - u(t - width)."""
    raise NotImplementedError


def unit_impulse(n: int, k: int = 0) -> np.ndarray:
    """Return the length-n discrete impulse delta[m - k]: 1 at index k, 0 elsewhere."""
    raise NotImplementedError


def sinusoid(t: np.ndarray, amplitude: float, freq_hz: float, phase: float = 0.0) -> np.ndarray:
    """Return amplitude * cos(2 pi freq_hz t + phase) [phase in rad]."""
    return amplitude * np.cos(2.0 * np.pi * freq_hz * np.asarray(t, dtype=float) + phase)


def energy(x: np.ndarray, dt: float) -> float:
    """Return the signal energy sum |x|^2 dt (a Riemann sum for integral |x(t)|^2 dt) [units of x^2 times s]."""
    raise NotImplementedError


def average_power(x: np.ndarray) -> float:
    """Return the average power mean |x|^2 over the record (energy divided by duration) [units of x^2]."""
    raise NotImplementedError


def even_part(x: np.ndarray) -> np.ndarray:
    """Return (x[n] + x[-n]) / 2 for x sampled on a grid symmetric about zero (t = -T ... T, odd length)."""
    raise NotImplementedError


def odd_part(x: np.ndarray) -> np.ndarray:
    """Return (x[n] - x[-n]) / 2 for x sampled on a grid symmetric about zero, so x = even + odd."""
    raise NotImplementedError


def time_shift(f: Callable[[np.ndarray], np.ndarray], t0: float) -> Callable[[np.ndarray], np.ndarray]:
    """Return the callable t -> f(t - t0): the signal delayed by t0 (moved right for t0 > 0)."""
    raise NotImplementedError


def time_scale(f: Callable[[np.ndarray], np.ndarray], a: float) -> Callable[[np.ndarray], np.ndarray]:
    """Return the callable t -> f(a t): compressed in time for |a| > 1, stretched for |a| < 1."""
    raise NotImplementedError


def time_reverse(f: Callable[[np.ndarray], np.ndarray]) -> Callable[[np.ndarray], np.ndarray]:
    """Return the callable t -> f(-t)."""
    raise NotImplementedError


# ---------------------------------------------------------------- property tests on sampled operators

def superposition_error(system: System, x1: np.ndarray, x2: np.ndarray, a: float, b: float) -> float:
    """Return max |S(a x1 + b x2) - (a S(x1) + b S(x2))| divided by (1 + max |a S(x1) + b S(x2)|).

    Zero (to rounding) for a linear system; a scale-free measure of how badly superposition fails otherwise.
    """
    raise NotImplementedError


def is_linear(system: System, n: int = 256, trials: int = 20, amplitude: float = 10.0, seed: int = 0,
              tol: float = 1e-9) -> bool:
    """Return True if superposition holds for ``trials`` random input pairs and random weights.

    Inputs are Gaussian with standard deviation ``amplitude``; the weights a, b are uniform in
    [-amplitude, amplitude]. A single failure above ``tol`` proves non-linearity; passing every trial is
    evidence, not proof (a saturation is linear for inputs below its limit).
    """
    raise NotImplementedError


def shift_signal(x: np.ndarray, k: int) -> np.ndarray:
    """Return x delayed by k >= 0 samples with zeros shifted in: y[n] = x[n - k] (initial rest)."""
    x = np.asarray(x, dtype=float)
    if k == 0:
        return x.copy()
    return np.concatenate([np.zeros(k), x[:-k]])


def is_time_invariant(system: System, n: int = 256, shifts: tuple[int, ...] = (1, 3, 17), trials: int = 5,
                      seed: int = 0, tol: float = 1e-9) -> bool:
    """Return True if delaying the input by k samples delays the output by k samples, for every k tried.

    Test inputs are random on the first half of the record and zero after it, and k <= n/4, so a finite
    look-ahead or the end of the record cannot fake a failure. Compares S(shift(x, k))[k:] with S(x)[:n-k].
    """
    raise NotImplementedError


def is_causal(system: System, n: int = 256, trials: int = 5, seed: int = 0, tol: float = 1e-12) -> bool:
    """Return True if changing the input from index k onward never changes the output before index k.

    Tries several k across the record for each random input; one change in y[:k] proves non-causality.
    """
    raise NotImplementedError


def is_memoryless(system: System, n: int = 256, trials: int = 5, seed: int = 0, tol: float = 1e-12) -> bool:
    """Return True if perturbing the single input sample x[k] changes the output only at index k."""
    raise NotImplementedError


def bounded_test_inputs(n: int, seed: int = 0) -> list[np.ndarray]:
    """Return bounded test inputs with |x| <= 1: a step, an alternating +-1 sequence, random signs, a slow sine."""
    rng = np.random.default_rng(seed)
    m = np.arange(n)
    return [np.ones(n), (-1.0) ** m, rng.choice([-1.0, 1.0], n), np.sin(2 * np.pi * m / 50.0)]


def is_bibo_stable(system: System, n: int = 2000, growth: float = 1.5, seed: int = 0) -> bool:
    """Return False if some bounded input (|x| <= 1) gives an output whose peak keeps growing.

    Each input of ``bounded_test_inputs`` is run for 4n samples; the system is declared unstable if the peak
    |y| over the whole record exceeds ``growth`` times the peak over the first n samples, or is not finite.
    The verdict is only as good as the record: a stable system whose time constant is longer than n
    samples looks unstable, and an unstable mode that no test input excites is missed.
    """
    raise NotImplementedError


def impulse_response(system: System, n: int) -> np.ndarray:
    """Return h = S(delta), the response to the unit impulse at n = 0, of length n."""
    raise NotImplementedError


def abs_sum(h: np.ndarray) -> float:
    """Return sum |h[n]|. An LTI system is BIBO stable if and only if this is finite (and |y| <= sum|h| max|x|)."""
    raise NotImplementedError


def classify(system: System, n: int = 256) -> dict[str, bool]:
    """Return the five property verdicts {linear, time_invariant, causal, memoryless, bibo_stable}."""
    raise NotImplementedError


# ---------------------------------------------------------------- a zoo of gyro-path systems (given)

def delay(k: int) -> System:
    """Return the k-sample delay y[n] = x[n - k] (a transport or processing delay)."""
    return lambda x: shift_signal(x, k)


def squarer() -> System:
    """Return y[n] = x[n]^2 (a power detector, as in a vibration-level monitor)."""
    return lambda x: np.asarray(x, dtype=float) ** 2


def accumulator(dt: float = 1.0) -> System:
    """Return y[n] = dt * sum_{m <= n} x[m]: integration of a rate to an angle (rectangle rule)."""
    return lambda x: dt * np.cumsum(np.asarray(x, dtype=float))


def time_varying_gain() -> System:
    """Return y[n] = n x[n]: a gain that grows with time (a scale factor schedule gone wrong)."""
    return lambda x: np.arange(len(x)) * np.asarray(x, dtype=float)


def bias_correction(bias: float) -> System:
    """Return y[n] = x[n] - bias: subtracting a gyro bias. Affine, not linear: S(0) = -bias != 0."""
    return lambda x: np.asarray(x, dtype=float) - bias


def saturation(limit: float) -> System:
    """Return y[n] = clip(x[n], -limit, limit): a converter or actuator that saturates."""
    return lambda x: np.clip(np.asarray(x, dtype=float), -limit, limit)


def moving_average(M: int) -> System:
    """Return the causal M-point moving average y[n] = (x[n] + ... + x[n - M + 1]) / M."""
    def S(x: np.ndarray) -> np.ndarray:
        c = np.cumsum(np.concatenate([np.zeros(M), np.asarray(x, dtype=float)]))
        return (c[M:] - c[:-M]) / M
    return S


def first_order_recursion(a: float, b: float = 1.0) -> System:
    """Return y[n] = a y[n - 1] + b x[n] from rest: a recursive low-pass for 0 < a < 1, unstable for |a| > 1."""
    def S(x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        y = np.empty_like(x)
        prev = 0.0
        for i, xi in enumerate(x):
            prev = a * prev + b * xi
            y[i] = prev
        return y
    return S


def centred_smoother() -> System:
    """Return y[n] = (x[n - 1] + x[n] + x[n + 1]) / 3 (uses one future sample; x beyond the record is 0)."""
    def S(x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        p = np.concatenate([[0.0], x, [0.0]])
        return (p[:-2] + p[1:-1] + p[2:]) / 3.0
    return S


def decimator() -> System:
    """Return y[n] = x[2n] (zero once 2n runs past the record): keeping every second sample."""
    def S(x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        y = np.zeros_like(x)
        m = (len(x) + 1) // 2
        y[:m] = x[::2][:m]
        return y
    return S


def gyro_path_zoo() -> dict[str, System]:
    """Return the named systems of the lesson's gyro signal path."""
    return {
        "delay 2": delay(2),
        "squarer": squarer(),
        "accumulator": accumulator(),
        "gain n x[n]": time_varying_gain(),
        "bias correction": bias_correction(0.05),
        "saturation 1": saturation(1.0),
        "moving average 4": moving_average(4),
        "low-pass a=0.9": first_order_recursion(0.9, 0.1),
        "recursion a=1.05": first_order_recursion(1.05),
        "centred smoother": centred_smoother(),
        "decimator": decimator(),
    }


# ---------------------------------------------------------------- the example

def main(out_dir: str | Path = ".") -> None:
    """Print the property table of the gyro-path zoo and the signal examples; save signals_and_systems.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("system              linear  TI     causal memless stable  sum|h| (first 2000 samples)")
    for name, S in gyro_path_zoo().items():
        v = classify(S)
        h = abs_sum(impulse_response(S, 2000)) if v["linear"] and v["time_invariant"] else float("nan")
        flags = "  ".join(f"{str(v[p]):6s}" for p in PROPERTIES)
        print(f"{name:18s}  {flags}  {h:.4g}")

    # energy and power
    dt = 1e-4
    t = np.arange(0.0, 20.0, dt)
    a = 2.0
    print(f"\nenergy of e^(-{a:g} t) u(t): {energy(np.exp(-a * t), dt):.5f} (exact 1/(2a) = {1 / (2 * a):.5f})")
    tp = np.arange(0.0, 10.0, dt)
    print(f"average power of 2 cos(2 pi 5 t) over 10 s: {average_power(sinusoid(tp, 2.0, 5.0)):.5f} (exact A^2/2 = 2)")
    ts = np.linspace(-2.0, 2.0, 40001)
    dts = ts[1] - ts[0]
    doublet = rect_pulse(ts + 1.0, 1.0) - rect_pulse(ts, 1.0)        # +1 on [-1, 0), -1 on [0, 1): an elevator doublet
    print(f"elevator doublet (+1 deg for 1 s, -1 deg for 1 s): energy {energy(doublet, dts):.4f} deg^2 s,"
          f" energy of even part {energy(even_part(doublet), dts):.1e}, of odd part {energy(odd_part(doublet), dts):.4f}")

    # integrating a gyro bias
    fs = 100.0
    bias_deg_s = 0.01
    hour = np.full(int(3600 * fs), bias_deg_s)
    angle = accumulator(1.0 / fs)(hour)
    print(f"\ngyro bias {bias_deg_s} deg/s integrated at {fs:.0f} Hz for one hour: attitude error {angle[-1]:.2f} deg"
          f" (accumulator: bounded input, unbounded output)")

    # saturation looks linear for small inputs
    print("\nsuperposition error of a +-1 saturation against input amplitude (x1, x2 Gaussian, a = 0.6, b = 0.3):")
    rng = np.random.default_rng(1)
    x1, x2 = rng.standard_normal(512), rng.standard_normal(512)
    sat = saturation(1.0)
    for amp in (0.1, 0.3, 0.5, 1.0, 3.0):
        e = superposition_error(sat, amp * x1, amp * x2, 0.6, 0.3)
        print(f"  amplitude {amp:4.1f}: error {e:.3e}")

    # the stability verdict depends on the record length
    print("\nBIBO verdict for the stable low-pass y = a y + (1 - a) x against record length n:")
    for a_lp in (0.9, 0.99, 0.999):
        verdicts = [is_bibo_stable(first_order_recursion(a_lp, 1 - a_lp), n=n) for n in (50, 500, 5000)]
        print(f"  a = {a_lp}: time constant {1 / (1 - a_lp):6.0f} samples; n = 50, 500, 5000 -> {verdicts}")

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    tt = np.linspace(-1.0, 3.0, 2001)
    axes[0].plot(tt, unit_step(tt), label="u(t)")
    axes[0].plot(tt, ramp(tt), label="r(t)")
    axes[0].plot(tt, rect_pulse(tt, 1.5) * 0.5, label="0.5 rect, width 1.5")
    axes[0].plot(tt, np.exp(-tt) * unit_step(tt) * np.cos(2 * np.pi * tt), label="e^{-t} cos 2 pi t u(t)")
    axes[0].set_ylim(-1.1, 2.2); axes[0].legend(fontsize=7); axes[0].grid(lw=0.4); axes[0].set_xlabel("t [s]")
    axes[0].set_title("elementary signals", fontsize=9)
    x_eo = np.exp(-ts) * unit_step(ts)
    axes[1].plot(ts, x_eo, label="x = e^{-t} u(t)")
    axes[1].plot(ts, even_part(x_eo), "--", label="even part")
    axes[1].plot(ts, odd_part(x_eo), ":", label="odd part")
    axes[1].legend(fontsize=7); axes[1].grid(lw=0.4); axes[1].set_xlabel("t [s]")
    axes[1].set_title("even and odd parts", fontsize=9)
    n_ax = np.arange(200)
    step = np.ones(200)
    for name in ("accumulator", "low-pass a=0.9", "moving average 4", "saturation 1"):
        axes[2].plot(n_ax, gyro_path_zoo()[name](0.5 * step), label=name)
    axes[2].set_ylim(-0.1, 3.0); axes[2].legend(fontsize=7); axes[2].grid(lw=0.4); axes[2].set_xlabel("n")
    axes[2].set_title("response to x[n] = 0.5 u[n]", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "signals_and_systems.png", dpi=150)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'signals_and_systems.png'}")


if __name__ == "__main__":
    import sys

    main(sys.argv[1] if len(sys.argv) > 1 else ".")
