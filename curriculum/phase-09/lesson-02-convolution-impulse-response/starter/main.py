"""P09.02 Convolution & Impulse Response: starter code.

Discrete convolution from the defining sum, checked against FFT multiplication with zero padding; circular
convolution and the wrap-around error that zero padding prevents; continuous convolution by a Riemann sum
(rectangle and trapezoid rules) against closed forms; the impulse and step response of a first-order
system, used twice: the RC anti-aliasing front end of a MEMS rate gyro and a first-order elevator servo.
The output of any LTI system is the input convolved with the impulse response, so a doublet command, a
step or a synthetic gyro record can be pushed through either system by one convolution.

Conventions:
    * Discrete sequences start at n = 0; ``convolve_direct(x, h)`` returns all len(x) + len(h) - 1 samples.
    * Continuous signals are sampled at t_k = k dt, k = 0 .. N-1, and are zero for t < 0 (causal).
    * A first-order system has transfer function K/(tau s + 1): impulse response (K/tau) e^{-t/tau},
      step response K (1 - e^{-t/tau}).

Complete every function marked NotImplementedError.

Run:  python main.py          (prints the checks and the two examples, writes convolution.png)
Test: pytest curriculum/phase-09/lesson-02-convolution-impulse-response
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


# ---------------------------------------------------------------- discrete convolution

def convolve_direct(x: np.ndarray, h: np.ndarray) -> np.ndarray:
    """Return y[n] = sum_k x[k] h[n - k] for n = 0 .. len(x) + len(h) - 2 (full linear convolution).

    Built from the definition: each input sample x[k] launches a copy of h scaled by x[k] and delayed by k,
    and the copies add (superposition plus time invariance). Cost len(x) * len(h) multiply-adds.
    """
    raise NotImplementedError


def next_power_of_two(n: int) -> int:
    """Return the smallest power of two >= n."""
    return 1 << max(0, int(n - 1).bit_length())


def convolve_fft(x: np.ndarray, h: np.ndarray) -> np.ndarray:
    """Return the full linear convolution by FFT multiplication, zero-padded to a power of two >= N + M - 1.

    The DFT multiplies *circular* convolutions; padding both sequences with zeros to at least N + M - 1
    samples makes the circular and the linear convolution identical (no wrap-around). Real inputs give a
    real result (the tiny imaginary rounding part is dropped).
    """
    raise NotImplementedError


def circular_convolve(x: np.ndarray, h: np.ndarray, L: int) -> np.ndarray:
    """Return the L-point circular convolution y[n] = sum_k x[k] h[(n - k) mod L], n = 0 .. L-1.

    Sequences shorter than L are zero-padded; longer ones are an error. Equal to the linear convolution
    when L >= len(x) + len(h) - 1; otherwise the tail of the linear result wraps around onto its start.
    """
    raise NotImplementedError


def causal_lti_output(h: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Return the first len(x) samples of x * h: the output of the causal LTI system with impulse response h."""
    raise NotImplementedError


def operation_counts(N: int, M: int) -> tuple[int, int]:
    """Return (direct, fft) real multiply counts for convolving lengths N and M.

    Direct: N * M. FFT: three length-L complex FFTs (two forward, one inverse) at (L/2) log2 L complex
    multiplies each, plus L for the product, with 4 real multiplies per complex one; L = next power of two
    >= N + M - 1. An estimate: it ignores additions, memory traffic and real-input savings.
    """
    L = next_power_of_two(N + M - 1)
    fft_complex = 3 * (L // 2) * int(np.log2(L)) + L
    return N * M, 4 * fft_complex


# ---------------------------------------------------------------- continuous convolution by Riemann sums

def convolve_continuous(f: np.ndarray, g: np.ndarray, dt: float, rule: str = "rect") -> np.ndarray:
    """Return y(t_k) = int_0^{t_k} f(tau) g(t_k - tau) d tau at t_k = k dt, k = 0 .. len(f) - 1.

    ``f`` and ``g`` are samples of causal signals on the same grid. ``rule="rect"`` is dt times the
    discrete convolution (first-order accurate); ``rule="trap"`` removes half of each end-point term
    (the trapezoid rule, second-order accurate for smooth integrands).
    """
    raise NotImplementedError


def exp_conv_exp(t: np.ndarray, a: float, b: float) -> np.ndarray:
    """Return the closed form of (e^{-at} u) * (e^{-bt} u) = (e^{-at} - e^{-bt}) / (b - a), or t e^{-at} if a == b."""
    raise NotImplementedError


def rect_conv_rect(t: np.ndarray, width: float) -> np.ndarray:
    """Return the closed form of rect * rect for two unit pulses on [0, width): a triangle of height width on [0, 2 width]."""
    raise NotImplementedError


# ---------------------------------------------------------------- first-order systems

def rc_time_constant(R: float, C: float) -> float:
    """Return tau = R C [s] of an RC low-pass, with corner frequency 1/(2 pi tau)."""
    return R * C


def first_order_impulse_response(t: np.ndarray, K: float, tau: float) -> np.ndarray:
    """Return h(t) = (K / tau) e^{-t/tau} u(t), the impulse response of K/(tau s + 1)."""
    raise NotImplementedError


def first_order_step_response(t: np.ndarray, K: float, tau: float) -> np.ndarray:
    """Return s(t) = K (1 - e^{-t/tau}) u(t), the step response of K/(tau s + 1)."""
    raise NotImplementedError


def step_from_impulse(h: np.ndarray, dt: float) -> np.ndarray:
    """Return the running integral s(t_k) = int_0^{t_k} h dt by the trapezoid rule (s[0] = 0)."""
    raise NotImplementedError


def impulse_from_step(s: np.ndarray, dt: float) -> np.ndarray:
    """Return h = ds/dt by central differences (one-sided at the ends): the impulse response is the derivative of the step response."""
    raise NotImplementedError


def impulse_response_centroid(h: np.ndarray, dt: float) -> float:
    """Return int t h dt / int h dt: the mean delay of a low-pass system (equal to tau for K/(tau s + 1))."""
    raise NotImplementedError


def sampled_first_order(K: float, tau: float, dt: float, n: int) -> np.ndarray:
    """Return h[n] of the difference equation y[n] = a y[n-1] + (1 - a) K x[n], a = e^{-dt/tau}: the
    sampled-data version of the first-order system (unit DC gain times K; h[n] = (1 - a) K a^n)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the examples (given)

def doublet(t: np.ndarray, amplitude: float, half_period: float, start: float = 0.0) -> np.ndarray:
    """Return a flight-test doublet: +amplitude on [start, start + T), -amplitude on [start + T, start + 2T)."""
    t = np.asarray(t, dtype=float)
    up = (t >= start) & (t < start + half_period)
    down = (t >= start + half_period) & (t < start + 2 * half_period)
    return amplitude * (up.astype(float) - down.astype(float))


def doublet_response_closed_form(t: np.ndarray, K: float, tau: float, amplitude: float, half_period: float,
                                 start: float = 0.0) -> np.ndarray:
    """Return the first-order response to a doublet as a sum of shifted step responses: A[s(t - t0) - 2 s(t - t0 - T) + s(t - t0 - 2T)]."""
    s = lambda tt: first_order_step_response(tt, K, tau)
    T = half_period
    return amplitude * (s(t - start) - 2 * s(t - start - T) + s(t - start - 2 * T))


def synthetic_gyro_record(fs: float = 20000.0, duration: float = 0.2, seed: int = 9) -> tuple[np.ndarray, np.ndarray]:
    """Return (t, omega) [s, deg/s], a SYNTHETIC rate-gyro signal: a 10 deg/s pitch-rate step at 20 ms, a
    3 deg/s vibration line at 1200 Hz and 0.2 deg/s rms white noise (seeded). Not a measurement."""
    rng = np.random.default_rng(seed)
    t = np.arange(0.0, duration, 1.0 / fs)
    omega = 10.0 * (t >= 0.02) + 3.0 * np.sin(2 * np.pi * 1200.0 * t) + 0.2 * rng.standard_normal(len(t))
    return t, omega


def main(out_dir: str | Path = ".") -> None:
    """Print the convolution checks and the RC and servo examples; save convolution.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(2)

    x, h = rng.standard_normal(1000), rng.standard_normal(200)
    yd, yf = convolve_direct(x, h), convolve_fft(x, h)
    print(f"direct vs FFT (zero-padded to {next_power_of_two(1199)}): max difference {np.max(np.abs(yd - yf)):.2e}")
    yc = circular_convolve(x, h, 1000)
    print(f"circular convolution with L = 1000 (no padding): max error against linear {np.max(np.abs(yc - yd[:1000])):.2f};"
          f" wrap-around region first {len(h) - 1} samples; error after it {np.max(np.abs(yc[199:] - yd[199:1000])):.1e}")
    for N, M in ((100, 10), (1000, 100), (10000, 1000), (100000, 1000)):
        d, f = operation_counts(N, M)
        print(f"  multiplies N = {N:6d}, M = {M:4d}: direct {d:.2e}, FFT {f:.2e}")

    # continuous convolution: convergence
    print("\n(e^-t u) * (e^-3t u) on [0, 8] s by Riemann sum: max error against (e^-t - e^-3t)/2")
    for dt in (0.1, 0.05, 0.025, 0.0125):
        t = np.arange(0.0, 8.0, dt)
        exact = exp_conv_exp(t, 1.0, 3.0)
        er = np.max(np.abs(convolve_continuous(np.exp(-t), np.exp(-3 * t), dt, "rect") - exact))
        et = np.max(np.abs(convolve_continuous(np.exp(-t), np.exp(-3 * t), dt, "trap") - exact))
        print(f"  dt = {dt:7.4f}: rectangle {er:.3e}, trapezoid {et:.3e}")

    # RC front end of a gyro
    R, C = 1.6e3, 1.0e-6
    tau = rc_time_constant(R, C)
    fs = 20000.0
    dt = 1.0 / fs
    t = np.arange(0.0, 0.02, dt)
    hrc = first_order_impulse_response(t, 1.0, tau)
    s_conv = convolve_continuous(hrc, np.ones_like(t), dt, "trap")
    print(f"\nRC front end R = {R / 1e3:.1f} kohm, C = {C * 1e6:.1f} uF: tau = {tau * 1e3:.2f} ms, corner {1 / (2 * np.pi * tau):.1f} Hz")
    print(f"  step by convolution vs 1 - e^(-t/tau): max error {np.max(np.abs(s_conv - first_order_step_response(t, 1.0, tau))):.2e}"
          f" (dt = {dt * 1e6:.0f} us)")
    print(f"  impulse-response centroid {impulse_response_centroid(first_order_impulse_response(np.arange(0, 0.05, dt), 1.0, tau), dt) * 1e3:.3f} ms (tau)")
    tg, wg = synthetic_gyro_record(fs)
    wf = convolve_continuous(wg, first_order_impulse_response(tg, 1.0, tau), dt, "trap")
    late = tg > 0.1
    vib = lambda sig: np.sqrt(2) * np.std(sig[late] - np.mean(sig[late]))
    gain_1200 = 1 / np.sqrt(1 + (2 * np.pi * 1200 * tau) ** 2)
    print(f"  SYNTHETIC gyro record: 1200 Hz line + noise amplitude ~ {vib(wg):.2f} deg/s before, {vib(wf):.2f} deg/s after the RC"
          f" (|H(1200 Hz)| = {gain_1200:.3f}); mean rate after the step {np.mean(wf[late]):.3f} deg/s")

    # first-order elevator servo
    K, tau_s = 1.0, 0.05
    dts = 1e-3
    ts = np.arange(0.0, 2.0, dts)
    cmd = doublet(ts, 2.0, 0.5, start=0.1)
    hs = first_order_impulse_response(ts, K, tau_s)
    ys = convolve_continuous(cmd, hs, dts, "trap")
    yc_s = doublet_response_closed_form(ts, K, tau_s, 2.0, 0.5, start=0.1)
    ys_r = convolve_continuous(cmd, hs, dts, "rect")
    print(f"\nelevator servo K/(tau s + 1), tau = {tau_s * 1e3:.0f} ms: 2 deg doublet, 0.5 s per half")
    print(f"  convolution vs sum of shifted steps (dt = 1 ms): max error rectangle {np.max(np.abs(ys_r - yc_s)):.3e} deg,"
          f" trapezoid {np.max(np.abs(ys - yc_s)):.3e} deg")
    for d in (2e-3, 5e-4):
        tt = np.arange(0.0, 2.0, d)
        e_tr = np.max(np.abs(convolve_continuous(doublet(tt, 2.0, 0.5, 0.1), first_order_impulse_response(tt, K, tau_s), d, "trap")
                             - doublet_response_closed_form(tt, K, tau_s, 2.0, 0.5, 0.1)))
        print(f"    trapezoid at dt = {d * 1e3:.1f} ms: max error {e_tr:.3e} deg")
    print(f"  peak surface deflection {yc_s.max():.4f} deg at t = {ts[np.argmax(yc_s)]:.3f} s; minimum {yc_s.min():.4f} deg")
    print(f"  sampled difference equation at dt = 1 ms: h[0] = {sampled_first_order(K, tau_s, dts, 3)[0]:.5f},"
          f" continuous h(0) dt = {K / tau_s * dts:.5f}")

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    axes[0].plot(ts, cmd, "k--", lw=0.8, label="command")
    axes[0].plot(ts, ys, label="convolution")
    axes[0].plot(ts, yc_s, ":", label="shifted steps")
    axes[0].set_xlabel("t [s]"); axes[0].set_ylabel("elevator [deg]"); axes[0].legend(fontsize=7); axes[0].grid(lw=0.4)
    axes[0].set_title("servo doublet response", fontsize=9)
    m = tg < 0.04
    axes[1].plot(tg[m] * 1e3, wg[m], lw=0.5, label="SYNTHETIC gyro input")
    axes[1].plot(tg[m] * 1e3, wf[m], label="after RC (convolution)")
    axes[1].set_xlabel("t [ms]"); axes[1].set_ylabel("rate [deg/s]"); axes[1].legend(fontsize=7); axes[1].grid(lw=0.4)
    axes[1].set_title("RC front end, SYNTHETIC record", fontsize=9)
    dts_list = np.array([0.1, 0.05, 0.025, 0.0125, 0.00625])
    er, et = [], []
    for d in dts_list:
        tt = np.arange(0.0, 8.0, d)
        ex = exp_conv_exp(tt, 1.0, 3.0)
        er.append(np.max(np.abs(convolve_continuous(np.exp(-tt), np.exp(-3 * tt), d, "rect") - ex)))
        et.append(np.max(np.abs(convolve_continuous(np.exp(-tt), np.exp(-3 * tt), d, "trap") - ex)))
    axes[2].loglog(dts_list, er, "o-", label="rectangle")
    axes[2].loglog(dts_list, et, "s-", label="trapezoid")
    axes[2].set_xlabel("dt [s]"); axes[2].set_ylabel("max error"); axes[2].legend(fontsize=7); axes[2].grid(True, which="both", lw=0.4)
    axes[2].set_title("Riemann-sum convergence", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "convolution.png", dpi=150)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'convolution.png'}")


if __name__ == "__main__":
    import sys

    main(sys.argv[1] if len(sys.argv) > 1 else ".")
