"""P09.21 Home Lab: Servo Feedback Control: reference solution.

Analysis code for the home lab in projects/p09-servo-control/:

* identification: a first-order-plus-delay model K e^{-Ls}/(tau s + 1) (Option A, a hobby servo with its
  own loop treated as part of the plant, command angle -> potentiometer angle) by the log-linear method,
  and a first-order-plus-integrator model K e^{-Ls}/(s (tau s + 1)) (Option B, a DC gear motor, duty ->
  output angle) from the straight tail of the step response and a log-linear fit of what is left; both
  refined by Levenberg-Marquardt least squares with parameter standard errors; the motor's dead zone
  from steps of several sizes;
* prediction: the closed loop under P or PI control by transfer-function algebra, its step response by an
  RK4 map, and its rise time (10-90 %), overshoot and 2 % settling time; the second-order formulas;
* design: P gain for an overshoot target on the integrator plant, PI gain by pole-zero cancellation for a
  rise-time target on the first-order plant;
* simulation: the sampled loop as the firmware runs it, with a duty (or command) limit, clamping
  anti-windup, the one-sample firmware delay, a plant dead time, an input dead zone, sensor quantisation
  and noise;
* comparison: predicted against measured (or SYNTHETIC) metrics.

A SYNTHETIC dataset (a seeded simulation of a "bench" whose true parameters are written below, labelled
synthetic in every file name, header line and row) lets a learner without the kit run the whole analysis.
It is a simulation, not a measurement.

Run:  python solution.py [out_dir]     (prints the analysis, writes the synthetic CSVs and the figure)
Test: AEFS_IMPL=solution pytest curriculum/phase-09/lesson-21-home-lab-servo-feedback
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Callable

import numpy as np

SEED = 2026
T_PRE = 0.1                     # s of logged samples before every step

# ---------------------------------------------------------------- the SYNTHETIC bench (truth)
# These numbers define the simulated bench that produces the SYNTHETIC dataset. The analysis never reads
# them; it has to recover them from the step tests, as it would on real hardware.
# Option B: small DC gear motor on a motor-driver module, duty in [-1, 1] -> output-shaft angle in deg,
# read by a 10 kOhm potentiometer (300 deg electrical travel) on a 12-bit ADC.
MOTOR_TRUE = {"K": 540.0, "tau": 0.045, "dz_pos": 0.10, "dz_neg": 0.11}   # deg/s per unit duty, s, duty
MOTOR_TS = 0.01                 # s, 100 Hz loop
MOTOR_DEG_PER_COUNT = 300.0 / 4096
MOTOR_NOISE_DEG = 0.08          # rms wiper noise before quantisation
MOTOR_DUTY_LIMIT = 0.6          # software duty limit in the firmware
MOTOR_OPEN_DUTIES = (0.3, 0.45, 0.6)
MOTOR_OPEN_T_END = 0.5
MOTOR_KP = 0.029                # duty per degree: the firmware gain designed in the lesson, rounded
MOTOR_KP_REDESIGN = 0.022       # duty per degree: redesigned on the sampled model (lesson part 5), rounded
MOTOR_DZ_COMP = 0.09            # duty added in the direction of the error (dead-zone compensation)
MOTOR_CLOSED_RUNS = ((15.0, MOTOR_KP), (60.0, MOTOR_KP), (15.0, MOTOR_KP_REDESIGN))   # (step deg, Kp)
MOTOR_CLOSED_T_END = 1.0
# Option A: hobby micro servo with its feedback potentiometer wiper brought out; command angle (deg) ->
# wiper angle (deg). The servo's own loop is part of the plant: a first-order lag, a frame latency, a
# slew-rate limit and a small internal dead band.
SERVO_TRUE = {"K": 0.95, "tau": 0.06, "delay": 0.012, "slew": 450.0, "deadband": 0.5}
SERVO_TS = 0.02                 # s, one loop per 50 Hz servo frame
SERVO_DEG_PER_COUNT = 0.066
SERVO_NOISE_DEG = 0.06
SERVO_CMD_LIMIT = 60.0          # software command limit, deg either side of centre
SERVO_OPEN_STEPS = (20.0, 40.0)
SERVO_OPEN_T_END = 0.6
SERVO_KP = 0.33                 # command deg per deg of error (designed in the lesson, rounded)
SERVO_KI = 5.9                  # command deg per (deg s)
SERVO_RISE_TARGET = 0.40        # s
SERVO_REF_STEPS = (20.0, 45.0)
SERVO_CLOSED_T_END = 1.2
MOTOR_OVERSHOOT_TARGET = 10.0   # percent

SYNTHETIC_FILES = ("synthetic_motor_open_loop.csv", "synthetic_motor_closed_loop.csv",
                   "synthetic_servo_open_loop.csv", "synthetic_servo_closed_loop.csv")
MEASURED_FILES = ("motor_open_loop.csv", "motor_closed_loop.csv", "servo_open_loop.csv", "servo_closed_loop.csv")


# ---------------------------------------------------------------- given helpers

def tf_to_ss(num: np.ndarray, den: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Controllable canonical state-space (A, B, C, D) of a proper transfer function num(s)/den(s).

    Coefficients are highest power first, as in numpy.polyval.
    """
    num = np.atleast_1d(np.asarray(num, float))
    den = np.atleast_1d(np.asarray(den, float))
    den = np.trim_zeros(den, "f")
    num = np.trim_zeros(num, "f") if np.any(num) else np.array([0.0])
    if len(num) > len(den):
        raise ValueError("transfer function must be proper")
    num, den = num / den[0], den / den[0]
    n = len(den) - 1
    num = np.concatenate([np.zeros(n + 1 - len(num)), num])
    D = float(num[0])
    if n == 0:
        return np.zeros((0, 0)), np.zeros((0, 1)), np.zeros((1, 0)), D
    b = num[1:] - D * den[1:]
    A = np.zeros((n, n))
    A[0, :] = -den[1:]
    A[1:, :-1] = np.eye(n - 1)
    B = np.zeros((n, 1))
    B[0, 0] = 1.0
    C = b.reshape(1, n)
    return A, B, C, D


def rk4_map(A: np.ndarray, B: np.ndarray, h: float) -> tuple[np.ndarray, np.ndarray]:
    """The classical RK4 step for x' = A x + B u with u held constant, as a linear map x -> M x + N u."""
    n = A.shape[0]
    Ah = A * h
    M = np.eye(n) + Ah + Ah @ Ah / 2 + Ah @ Ah @ Ah / 6 + Ah @ Ah @ Ah @ Ah / 24
    N = (np.eye(n) + Ah / 2 + Ah @ Ah / 6 + Ah @ Ah @ Ah / 24) @ B * h
    return M, N


def quantize(y: np.ndarray | float, q: float) -> np.ndarray | float:
    """Round to the nearest multiple of q (q = 0: no quantisation)."""
    if q <= 0:
        return y
    return q * np.round(np.asarray(y) / q)


def dead_zone(u: float, dz_pos: float, dz_neg: float | None = None) -> float:
    """Input dead zone: 0 for -dz_neg <= u <= dz_pos, else u minus the threshold (a friction/driver model)."""
    dz_neg = dz_pos if dz_neg is None else dz_neg
    if u > dz_pos:
        return u - dz_pos
    if u < -dz_neg:
        return u + dz_neg
    return 0.0


def pi_update(integ: float, e: float, kp: float, ki: float, Ts: float, u_min: float, u_max: float,
              dz_comp: float = 0.0) -> tuple[float, float]:
    """One firmware control update: P or PI with clamping anti-windup and dead-zone compensation.

    u = kp e + ki I + dz_comp sign(e), clipped to [u_min, u_max]; then I += e Ts (forward Euler: this
    sample's error reaches the integral at the next sample) unless the output is saturated and the error
    would drive it further into saturation. Returns (u, new I). This is the same law as the firmware listing
    in projects/p09-servo-control/procedure.md.
    """
    i_new = integ + e * Ts
    u = kp * e + ki * integ + dz_comp * float(np.sign(e))
    if u > u_max:
        if e > 0:
            i_new = integ
        u = u_max
    elif u < u_min:
        if e < 0:
            i_new = integ
        u = u_min
    return float(u), float(i_new)


# ---------------------------------------------------------------- identification

def first_order_step_response(t: np.ndarray, K: float, tau: float, delay: float, u: float = 1.0) -> np.ndarray:
    """Response of K e^{-delay s}/(tau s + 1) to a step of size u at t = 0, from rest."""
    s = np.clip(np.asarray(t, float) - delay, 0.0, None)
    return K * u * (1.0 - np.exp(-s / tau))


def integrator_step_response(t: np.ndarray, K: float, tau: float, delay: float, u: float = 1.0) -> np.ndarray:
    """Response of K e^{-delay s}/(s (tau s + 1)) to a step of size u at t = 0, from rest.

    y = K u [s - tau (1 - e^{-s/tau})], s = max(t - delay, 0): a ramp of slope K u that lags by tau + delay.
    """
    s = np.clip(np.asarray(t, float) - delay, 0.0, None)
    return K * u * (s - tau * (1.0 - np.exp(-s / tau)))


def fit_first_order_loglinear(t: np.ndarray, y: np.ndarray, u: float, t_step: float = 0.0,
                              tail_fraction: float = 0.2, window: tuple[float, float] = (0.1, 0.9)) -> dict:
    """Log-linear fit of K e^{-Ls}/(tau s + 1) to a step test.

    y0 = mean before the step, y_ss = mean of the last ``tail_fraction`` of the post-step samples,
    K = (y_ss - y0)/u. The unfinished fraction z = 1 - (y - y0)/(y_ss - y0) obeys ln z = L/tau - (t - t_step)/tau,
    so a straight line through ln z for the samples with window[0] < z < window[1] gives tau = -1/slope and
    L = intercept * tau (clipped at 0). Returns dict(K, tau, delay, y0, y_ss).
    """
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    pre = t < t_step
    post = ~pre
    y0 = float(np.mean(y[pre])) if pre.any() else float(y[0])
    yp = y[post]
    n_tail = max(3, int(round(tail_fraction * yp.size)))
    yss = float(np.mean(yp[-n_tail:]))
    dy = yss - y0
    if dy == 0 or u == 0:
        raise ValueError("no response to the step")
    z = 1.0 - (y - y0) / dy
    sel = post & (z > window[0]) & (z < window[1])
    if sel.sum() < 3:
        raise ValueError("too few samples in the log-linear window; sample faster or use a longer step")
    slope, icpt = np.polyfit(t[sel] - t_step, np.log(z[sel]), 1)
    tau = -1.0 / slope
    return {"K": dy / u, "tau": float(tau), "delay": float(max(0.0, icpt * tau)), "y0": y0, "y_ss": yss}


def fit_integrator_loglinear(t: np.ndarray, y: np.ndarray, u: float, t_step: float = 0.0,
                             tail_fraction: float = 0.4, window: tuple[float, float] = (0.1, 0.6)) -> dict:
    """Fit K e^{-Ls}/(s (tau s + 1)) to a step test from the tail line and a log-linear fit of the remainder.

    The tail (last ``tail_fraction`` of the post-step samples) is the straight line y - y0 = m (s - T0) with
    m = K u and T0 = tau + L. What is left, r = (y - y0) - m (s - T0), equals K u tau e^{-(s - L)/tau} for
    s > L: a line through ln r for window[0] r_ref < r < window[1] r_ref, r_ref = m T0, before the tail,
    gives tau; then L = max(T0 - tau, 0). Assumes L is small compared with tau. Returns dict(K, tau, delay,
    y0, slope).
    """
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    pre = t < t_step
    post = ~pre
    y0 = float(np.mean(y[pre])) if pre.any() else float(y[0])
    s = t - t_step
    idx = np.flatnonzero(post)
    n_tail = max(3, int(round(tail_fraction * idx.size)))
    tail = idx[-n_tail:]
    m, c = np.polyfit(s[tail], y[tail] - y0, 1)
    if m == 0 or u == 0:
        raise ValueError("no response to the step")
    T0 = -c / m
    r = (y - y0) - m * (s - T0)
    r_ref = m * T0
    before_tail = np.zeros_like(post)
    before_tail[idx[:-n_tail]] = True
    sel = before_tail & (r / r_ref > window[0]) & (r / r_ref < window[1])
    if sel.sum() < 3:
        raise ValueError("too few samples in the log-linear window; sample faster or use a longer step")
    slope, _ = np.polyfit(s[sel], np.log(r[sel] / r_ref), 1)
    tau = -1.0 / slope
    return {"K": float(m / u), "tau": float(tau), "delay": float(max(0.0, T0 - tau)), "y0": y0,
            "slope": float(m)}


def refine_least_squares(model: Callable[[np.ndarray, np.ndarray], np.ndarray], t: np.ndarray, y: np.ndarray,
                         p0: np.ndarray, lower: np.ndarray | None = None, max_iter: int = 200,
                         tol: float = 1e-12) -> dict:
    """Levenberg-Marquardt least squares: minimise sum (y - model(t, p))^2 from p0.

    Central-difference Jacobian; damping lambda halved on success and multiplied by 10 on failure; optional
    lower bounds by clipping. Returns dict(p, se, rms, cov, n_iter): se = sqrt(diag(s^2 (J^T J)^-1)) with
    s^2 = SSR/(n - p), the usual linearised standard errors.
    """
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    p = np.asarray(p0, float).copy()
    lo = None if lower is None else np.asarray(lower, float)
    if lo is not None:
        p = np.maximum(p, lo)

    def resid(pp: np.ndarray) -> np.ndarray:
        return y - model(t, pp)

    def jac(pp: np.ndarray) -> np.ndarray:
        J = np.empty((y.size, pp.size))
        for j in range(pp.size):
            h = 1e-6 * max(abs(pp[j]), 1e-3)
            dp = np.zeros_like(pp)
            dp[j] = h
            J[:, j] = (model(t, pp + dp) - model(t, pp - dp)) / (2 * h)
        return J

    r = resid(p)
    cost = float(r @ r)
    lam = 1e-3
    it = 0
    for it in range(1, max_iter + 1):
        J = jac(p)
        g = J.T @ r
        H = J.T @ J
        step = np.linalg.solve(H + lam * np.diag(np.diag(H) + 1e-30), g)
        p_new = p + step
        if lo is not None:
            p_new = np.maximum(p_new, lo)
        r_new = resid(p_new)
        cost_new = float(r_new @ r_new)
        if cost_new < cost:
            converged = cost - cost_new <= tol * max(cost, 1e-300) and np.all(np.abs(step) <= 1e-10 + 1e-8 * np.abs(p))
            p, r, cost = p_new, r_new, cost_new
            lam = max(lam / 2, 1e-12)
            if converged:
                break
        else:
            lam *= 10
            if lam > 1e12:
                break
    J = jac(p)
    dof = max(y.size - p.size, 1)
    s2 = cost / dof
    cov = s2 * np.linalg.pinv(J.T @ J)
    return {"p": p, "se": np.sqrt(np.abs(np.diag(cov))), "rms": float(np.sqrt(cost / y.size)), "cov": cov,
            "n_iter": it}


def fit_dead_zone(duties: np.ndarray, slopes: np.ndarray) -> tuple[float, float]:
    """Gain and dead zone from steps of several sizes: slope = K (u - u_dz) for u > u_dz, by a straight line.

    ``slopes`` are the steady output rates (deg/s) reached at each step size u. Returns (K, u_dz).
    """
    duties = np.asarray(duties, float)
    slopes = np.asarray(slopes, float)
    if duties.size < 2:
        raise ValueError("need at least two step sizes")
    a, b = np.polyfit(duties, slopes, 1)
    return float(a), float(-b / a)


# ---------------------------------------------------------------- metrics and second-order formulas

def step_metrics(t: np.ndarray, y: np.ndarray, t_step: float = 0.0, y_initial: float | None = None,
                 y_final: float | None = None, ref: float | None = None, band: float = 0.02,
                 tail_fraction: float = 0.2) -> dict:
    """Rise time (10-90 %), percent overshoot, peak time and settling time (band, default 2 %) of a sampled step.

    y_initial defaults to the mean before t_step (or y[0]); y_final to the mean of the last ``tail_fraction``
    of the post-step samples. Crossing times are interpolated linearly between samples. Times are measured
    from t_step. The settling time is when |y - y_final| last leaves the band (nan if it never settles
    inside the record). ss_error = ref - y_final when ref is given. Returns a dict.
    """
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    pre = t < t_step
    post = ~pre
    if y_initial is None:
        y_initial = float(np.mean(y[pre])) if pre.any() else float(y[0])
    tp, yp = t[post], y[post]
    if y_final is None:
        n_tail = max(3, int(round(tail_fraction * yp.size)))
        y_final = float(np.mean(yp[-n_tail:]))
    dy = y_final - y_initial
    if dy == 0:
        raise ValueError("no change between initial and final value")
    z = (yp - y_initial) / dy

    def first_crossing(level: float) -> float:
        k = np.flatnonzero(z >= level)
        if k.size == 0:
            return np.nan
        k = int(k[0])
        if k == 0:
            return float(tp[0])
        return float(tp[k - 1] + (level - z[k - 1]) * (tp[k] - tp[k - 1]) / (z[k] - z[k - 1]))

    t10, t90 = first_crossing(0.1), first_crossing(0.9)
    kpk = int(np.argmax(z))
    overshoot = max(0.0, float(z[kpk] - 1.0)) * 100.0
    dev = np.abs(z - 1.0) - band
    out_idx = np.flatnonzero(dev > 0)
    if out_idx.size == 0:
        ts = float(tp[0] - t_step)
    elif out_idx[-1] == z.size - 1:
        ts = np.nan
    else:
        k = int(out_idx[-1])
        ts = float(tp[k] + dev[k] * (tp[k + 1] - tp[k]) / (dev[k] - dev[k + 1])) - t_step
    return {"rise_time": t90 - t10, "t10": t10 - t_step, "t90": t90 - t_step, "overshoot_pct": overshoot,
            "peak_time": float(tp[kpk] - t_step), "settling_time": ts, "y_initial": y_initial,
            "y_final": y_final, "ss_error": (np.nan if ref is None else float(ref - y_final))}


def second_order_step(t: np.ndarray, zeta: float, wn: float) -> np.ndarray:
    """Unit-step response of wn^2/(s^2 + 2 zeta wn s + wn^2) from rest, for any zeta >= 0 (closed form)."""
    t = np.asarray(t, float)
    if zeta < 1:
        wd = wn * np.sqrt(1 - zeta ** 2)
        return 1 - np.exp(-zeta * wn * t) * (np.cos(wd * t) + zeta / np.sqrt(1 - zeta ** 2) * np.sin(wd * t))
    if zeta == 1:
        return 1 - np.exp(-wn * t) * (1 + wn * t)
    r = np.sqrt(zeta ** 2 - 1)
    p1, p2 = -wn * (zeta - r), -wn * (zeta + r)
    return 1 + (p2 * np.exp(p1 * t) - p1 * np.exp(p2 * t)) / (p1 - p2)


def _bisect(f: Callable[[float], float], a: float, b: float, n: int = 100) -> float:
    fa = f(a)
    for _ in range(n):
        m = 0.5 * (a + b)
        fm = f(m)
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    return 0.5 * (a + b)


def second_order_metrics(zeta: float, wn: float, band: float = 0.02) -> dict:
    """Closed-form and exact metrics of the standard second-order step (0 < zeta < 1 for overshoot).

    overshoot_pct = 100 exp(-pi zeta/sqrt(1 - zeta^2)); peak_time = pi/(wn sqrt(1 - zeta^2));
    settling_approx = 4/(zeta wn) (the envelope estimate for a 2 % band); rise_time and settling_time are the
    exact 10-90 % rise and band settling of the closed form, found by bisection on its crossings.
    """
    f = lambda tt: float(second_order_step(np.array([tt]), zeta, wn)[0])  # noqa: E731
    if zeta < 1:
        wd = wn * np.sqrt(1 - zeta ** 2)
        mp = 100 * np.exp(-np.pi * zeta / np.sqrt(1 - zeta ** 2))
        tpk = np.pi / wd
    else:
        mp, tpk = 0.0, np.inf
    t_first = tpk if zeta < 1 else 50.0 / wn
    t10 = _bisect(lambda tt: f(tt) - 0.1, 0.0, t_first)
    t90 = _bisect(lambda tt: f(tt) - 0.9, 0.0, t_first)
    # exact settling: last time |y - 1| = band, from a dense scan then bisection
    tg = np.linspace(0, 20.0 / (min(zeta, 1.0) * wn), 200001)
    dev = np.abs(second_order_step(tg, zeta, wn) - 1) - band
    k = int(np.flatnonzero(dev > 0)[-1])
    ts = _bisect(lambda tt: abs(f(tt) - 1) - band, tg[k], tg[k + 1])
    return {"overshoot_pct": float(mp), "peak_time": float(tpk), "rise_time": t90 - t10,
            "settling_time": ts, "settling_approx": 4.0 / (zeta * wn)}


def zeta_from_overshoot(overshoot_pct: float) -> float:
    """Damping ratio of the standard second-order system with the given percent overshoot."""
    if not 0 < overshoot_pct < 100:
        raise ValueError("overshoot must be between 0 and 100 %")
    L = np.log(overshoot_pct / 100.0)
    return float(-L / np.sqrt(np.pi ** 2 + L ** 2))


def overshoot_from_zeta(zeta: float) -> float:
    """Percent overshoot 100 exp(-pi zeta/sqrt(1 - zeta^2)) of the standard second-order step (0 for zeta >= 1)."""
    if zeta >= 1:
        return 0.0
    return float(100 * np.exp(-np.pi * zeta / np.sqrt(1 - zeta ** 2)))


# ---------------------------------------------------------------- design and prediction

def design_p_position(K: float, tau: float, overshoot_pct: float) -> dict:
    """P gain on K/(s (tau s + 1)) for a percent-overshoot target.

    Closed loop Kp K/(tau s^2 + s + Kp K): wn = sqrt(Kp K/tau), zeta = 1/(2 sqrt(tau Kp K)), so
    Kp = 1/(4 tau zeta^2 K). Returns dict(kp, zeta, wn).
    """
    zeta = zeta_from_overshoot(overshoot_pct)
    kpk = 1.0 / (4 * tau * zeta ** 2)
    return {"kp": kpk / K, "zeta": zeta, "wn": float(np.sqrt(kpk / tau))}


def design_pi_first_order(K: float, tau: float, rise_time: float) -> dict:
    """PI gains on K/(tau s + 1) by pole-zero cancellation for a 10-90 % rise-time target.

    The PI zero at -ki/kp cancels the plant pole (ki = kp/tau), leaving the loop K kp/(tau s) and the
    closed loop 1/(tau_cl s + 1), tau_cl = tau/(K kp), whose rise time is tau_cl ln 9. Returns
    dict(kp, ki, tau_cl).
    """
    tau_cl = rise_time / np.log(9.0)
    kp = tau / (K * tau_cl)
    return {"kp": float(kp), "ki": float(kp / tau), "tau_cl": float(tau_cl)}


def _padd(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    n = max(len(a), len(b))
    return np.concatenate([np.zeros(n - len(a)), a]) + np.concatenate([np.zeros(n - len(b)), b])


def closed_loop_tf(num: np.ndarray, den: np.ndarray, kp: float, ki: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """Unity-feedback closed loop T = C G/(1 + C G) with C = kp (ki = 0) or kp + ki/s, as (num, den).

    Polynomials highest power first; leading zeros removed and the denominator made monic.
    """
    num = np.asarray(num, float)
    den = np.asarray(den, float)
    cn, cd = (np.array([kp, ki]), np.array([1.0, 0.0])) if ki != 0 else (np.array([kp]), np.array([1.0]))
    ln = np.convolve(cn, num)
    ld = np.convolve(cd, den)
    dcl = np.trim_zeros(_padd(ld, ln), "f")
    ncl = np.trim_zeros(ln, "f")
    return ncl / dcl[0], dcl / dcl[0]


def predict_step(num: np.ndarray, den: np.ndarray, amplitude: float, t_end: float,
                 dt: float = 1e-4) -> tuple[np.ndarray, np.ndarray]:
    """Step response of num/den to a step of ``amplitude`` at t = 0, from rest, by the RK4 map with step dt."""
    A, B, C, D = tf_to_ss(num, den)
    n = int(round(t_end / dt))
    t = np.arange(n + 1) * dt
    y = np.empty(n + 1)
    M, N = rk4_map(A, B, dt)
    x = np.zeros((A.shape[0], 1))
    for k in range(n + 1):
        y[k] = (C @ x)[0, 0] + D * amplitude
        x = M @ x + N * amplitude
    return t, y


def predict_closed_loop(num: np.ndarray, den: np.ndarray, kp: float, ki: float, r_step: float, t_end: float,
                        dt: float = 1e-4) -> dict:
    """Model-predicted closed-loop step: the continuous LTI loop, no delay, no saturation, no quantisation.

    Returns dict(t, y, num_cl, den_cl, poles, metrics), with metrics from step_metrics using the exact
    initial value 0 and the final value r_step T(0) from the final-value theorem.
    """
    ncl, dcl = closed_loop_tf(num, den, kp, ki)
    t, y = predict_step(ncl, dcl, r_step, t_end, dt)
    y_final = r_step * ncl[-1] / dcl[-1]
    m = step_metrics(t, y, 0.0, y_initial=0.0, y_final=y_final, ref=r_step)
    return {"t": t, "y": y, "num_cl": ncl, "den_cl": dcl, "poles": np.roots(dcl), "metrics": m}


# ---------------------------------------------------------------- the sampled loop

def simulate_closed_loop(num: np.ndarray, den: np.ndarray, kp: float, ki: float, r_step: float, Ts: float,
                         t_end: float, *, u_min: float = -np.inf, u_max: float = np.inf, delay_samples: int = 1,
                         plant_delay: float = 0.0, dz: float = 0.0, dz_comp: float = 0.0, quant: float = 0.0,
                         noise_sd: float = 0.0, substeps: int = 20, t_pre: float = 0.0, seed: int = 0) -> dict:
    """Simulate the firmware loop on the model plant num/den (output from rest at 0).

    Every Ts the firmware reads y_meas = quantize(y + noise, quant), forms e = r - y_meas (r = r_step for
    t >= 0, else 0) and calls ``pi_update`` (limit [u_min, u_max], clamping anti-windup, dead-zone
    compensation dz_comp). The value computed at tick k is applied (zero-order hold) over tick
    k + delay_samples (delay_samples = 1 is the firmware's write-last-then-compute structure). The plant sees
    the applied value ``plant_delay`` seconds later, through a symmetric input dead zone dz. The plant is
    integrated by the RK4 map with ``substeps`` steps per sample; plant_delay is rounded to whole substeps.
    Returns dict(t, ref, y, y_true, u) with u the applied value over [t_k, t_k + Ts).
    """
    rng = np.random.default_rng(seed)
    A, B, C, D = tf_to_ss(num, den)
    h = Ts / substeps
    M, N = rk4_map(A, B, h)
    n_pre = int(round(t_pre / Ts))
    n = n_pre + int(round(t_end / Ts)) + 1
    t = (np.arange(n) - n_pre) * Ts
    n_pd = int(round(plant_delay / h))
    x = np.zeros((A.shape[0], 1))
    u_comp = np.zeros(n)            # computed at tick k
    u_app = np.zeros(n)             # applied over [t_k, t_k + Ts)
    u_sub = np.zeros(n * substeps)  # applied value on the substep grid (before the plant delay)
    y_true = np.zeros(n)
    y_meas = np.zeros(n)
    ref = np.where(t >= 0, r_step, 0.0)
    integ = 0.0
    for k in range(n):
        y_true[k] = (C @ x)[0, 0]
        yn = y_true[k] + (rng.normal(0.0, noise_sd) if noise_sd > 0 else 0.0)
        y_meas[k] = quantize(yn, quant)
        u_comp[k], integ = pi_update(integ, ref[k] - y_meas[k], kp, ki, Ts, u_min, u_max, dz_comp)
        u_app[k] = u_comp[k - delay_samples] if k >= delay_samples else 0.0
        for j in range(substeps):
            g = k * substeps + j
            u_sub[g] = u_app[k]
            ud = u_sub[g - n_pd] if g >= n_pd else 0.0
            x = M @ x + N * dead_zone(ud, dz)
    return {"t": t, "ref": ref, "y": y_meas, "y_true": y_true, "u": u_app}


def simulate_open_loop(num: np.ndarray, den: np.ndarray, u_step: float, Ts: float, t_end: float, *,
                       plant_delay: float = 0.0, dz: float = 0.0, quant: float = 0.0, noise_sd: float = 0.0,
                       substeps: int = 20, t_pre: float = 0.0, seed: int = 0) -> dict:
    """Simulate an open-loop step test on the model plant: u = u_step from t = 0, sampled every Ts.

    Same plant, dead zone, plant delay, quantisation and noise conventions as simulate_closed_loop.
    Returns dict(t, y, y_true, u).
    """
    rng = np.random.default_rng(seed)
    A, B, C, D = tf_to_ss(num, den)
    h = Ts / substeps
    M, N = rk4_map(A, B, h)
    n_pre = int(round(t_pre / Ts))
    n = n_pre + int(round(t_end / Ts)) + 1
    t = (np.arange(n) - n_pre) * Ts
    n_pd = int(round(plant_delay / h))
    u_app = np.where(t >= 0, u_step, 0.0)
    u_sub = np.repeat(u_app, substeps)
    x = np.zeros((A.shape[0], 1))
    y_true = np.zeros(n)
    y = np.zeros(n)
    for k in range(n):
        y_true[k] = (C @ x)[0, 0]
        y[k] = quantize(y_true[k] + (rng.normal(0.0, noise_sd) if noise_sd > 0 else 0.0), quant)
        for j in range(substeps):
            g = k * substeps + j
            ud = u_sub[g - n_pd] if g >= n_pd else 0.0
            x = M @ x + N * dead_zone(ud, dz)
    return {"t": t, "y": y, "y_true": y_true, "u": u_app}


def compare_metrics(pred: dict, meas: dict,
                    keys: tuple[str, ...] = ("rise_time", "overshoot_pct", "settling_time", "ss_error")) -> dict:
    """Side-by-side table: key -> (predicted, measured, measured - predicted, relative difference or nan)."""
    out = {}
    for k in keys:
        p, m = float(pred[k]), float(meas[k])
        rel = (m - p) / p if p not in (0.0,) and np.isfinite(p) and abs(p) > 1e-12 else np.nan
        out[k] = (p, m, m - p, rel)
    return out


def design_p_sampled(num: np.ndarray, den: np.ndarray, r_step: float, Ts: float, overshoot_pct: float,
                     kp_lo: float, kp_hi: float, t_end: float = 1.0, n_iter: int = 40, **sim_kwargs) -> float:
    """P gain giving ``overshoot_pct`` on the sampled loop (simulate_closed_loop, noise-free), by bisection.

    Assumes the overshoot rises with kp on [kp_lo, kp_hi] and brackets the target; extra keyword arguments
    (u_min, u_max, delay_samples, plant_delay, dz, dz_comp, quant, ...) go to simulate_closed_loop.
    """
    def os_of(kp: float) -> float:
        sim = simulate_closed_loop(num, den, kp, 0.0, r_step, Ts, t_end, **sim_kwargs)
        return step_metrics(sim["t"], sim["y"])["overshoot_pct"] - overshoot_pct

    f_lo, f_hi = os_of(kp_lo), os_of(kp_hi)
    if f_lo > 0 or f_hi < 0:
        raise ValueError("the gain interval does not bracket the overshoot target")
    for _ in range(n_iter):
        mid = 0.5 * (kp_lo + kp_hi)
        if os_of(mid) > 0:
            kp_hi = mid
        else:
            kp_lo = mid
    return 0.5 * (kp_lo + kp_hi)


# ---------------------------------------------------------------- the SYNTHETIC bench (given)

def _motor_world(duty_step: float | None, r_step: float | None, kp: float, rng: np.random.Generator) -> dict:
    """The simulated motor bench: exact ZOH integration of the true plant with its asymmetric dead zone.

    Open loop when duty_step is given (duty applied from t = 0), closed loop (gain kp, MOTOR_DZ_COMP, duty
    limit, one-sample firmware delay) when r_step is given.
    """
    K, tau = MOTOR_TRUE["K"], MOTOR_TRUE["tau"]
    Ts = MOTOR_TS
    t_end = MOTOR_OPEN_T_END if duty_step is not None else MOTOR_CLOSED_T_END
    n_pre = int(round(T_PRE / Ts))
    n = n_pre + int(round(t_end / Ts)) + 1
    t = (np.arange(n) - n_pre) * Ts
    e = np.exp(-Ts / tau)
    theta = omega = 0.0
    u_pending = integ = 0.0
    rec = {"t": t, "theta": np.zeros(n), "u": np.zeros(n), "ref": np.zeros(n), "kp": np.full(n, kp),
           "ki": np.zeros(n)}
    for k in range(n):
        y = float(quantize(theta + rng.normal(0.0, MOTOR_NOISE_DEG), MOTOR_DEG_PER_COUNT))
        if duty_step is not None:
            u = duty_step if t[k] >= 0 else 0.0
        else:
            u = u_pending
            r = r_step if t[k] >= 0 else 0.0
            rec["ref"][k] = r
            u_pending, integ = pi_update(integ, r - y, kp, 0.0, Ts, -MOTOR_DUTY_LIMIT, MOTOR_DUTY_LIMIT,
                                         MOTOR_DZ_COMP)
        rec["theta"][k], rec["u"][k] = y, u
        ku = K * dead_zone(u, MOTOR_TRUE["dz_pos"], MOTOR_TRUE["dz_neg"])
        theta, omega = theta + ku * Ts + (omega - ku) * tau * (1 - e), omega * e + ku * (1 - e)
    return rec


def _servo_world(cmd_step: float | None, r_step: float | None, rng: np.random.Generator) -> dict:
    """The simulated servo bench: the servo's internal loop as a slew-limited first-order lag with a dead band
    and a frame latency, integrated by Euler steps of 0.5 ms (the definition of the synthetic truth).

    Open loop when cmd_step is given, closed loop (SERVO_KP, SERVO_KI, command limit, one-sample firmware
    delay) when r_step is given.
    """
    K, tau, L = SERVO_TRUE["K"], SERVO_TRUE["tau"], SERVO_TRUE["delay"]
    Ts, h = SERVO_TS, 0.0005
    nsub = int(round(Ts / h))
    jL = int(round(L / h))
    t_end = SERVO_OPEN_T_END if cmd_step is not None else SERVO_CLOSED_T_END
    n_pre = int(round(T_PRE / Ts))
    n = n_pre + int(round(t_end / Ts)) + 1
    t = (np.arange(n) - n_pre) * Ts
    theta = 0.0
    c_prev = u_pending = integ = 0.0
    rec = {"t": t, "theta": np.zeros(n), "u": np.zeros(n), "ref": np.zeros(n), "kp": np.full(n, SERVO_KP),
           "ki": np.full(n, SERVO_KI)}
    for k in range(n):
        y = float(quantize(theta + rng.normal(0.0, SERVO_NOISE_DEG), SERVO_DEG_PER_COUNT))
        if cmd_step is not None:
            c = cmd_step if t[k] >= 0 else 0.0
        else:
            c = u_pending
            r = r_step if t[k] >= 0 else 0.0
            rec["ref"][k] = r
            u_pending, integ = pi_update(integ, r - y, SERVO_KP, SERVO_KI, Ts, -SERVO_CMD_LIMIT, SERVO_CMD_LIMIT)
        rec["theta"][k], rec["u"][k] = y, c
        for j in range(nsub):
            target = K * (c if j >= jL else c_prev)
            err = target - theta
            if abs(err) > SERVO_TRUE["deadband"]:
                theta += float(np.clip(err / tau, -SERVO_TRUE["slew"], SERVO_TRUE["slew"])) * h
        c_prev = c
    return rec


def make_synthetic_dataset(seed: int = SEED) -> dict[str, list[dict]]:
    """Generate the SYNTHETIC dataset from the simulated bench: four lists of runs (dicts of arrays).

    Keys: motor_open (one run per duty in MOTOR_OPEN_DUTIES), motor_closed (MOTOR_CLOSED_RUNS), servo_open
    (SERVO_OPEN_STEPS), servo_closed (SERVO_REF_STEPS). A simulation, not a measurement.
    """
    rng = np.random.default_rng(seed)
    out: dict[str, list[dict]] = {"motor_open": [], "motor_closed": [], "servo_open": [], "servo_closed": []}
    for d in MOTOR_OPEN_DUTIES:
        out["motor_open"].append(_motor_world(d, None, 0.0, rng))
    for r, kp in MOTOR_CLOSED_RUNS:
        out["motor_closed"].append(_motor_world(None, r, kp, rng))
    for c in SERVO_OPEN_STEPS:
        out["servo_open"].append(_servo_world(c, None, rng))
    for r in SERVO_REF_STEPS:
        out["servo_closed"].append(_servo_world(None, r, rng))
    return out


def write_synthetic_csvs(out_dir: str | Path, seed: int = SEED) -> list[Path]:
    """Write the four SYNTHETIC CSV files (header line, column names, one row per sample, note = synthetic)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ds = make_synthetic_dataset(seed)
    head = (f"# SYNTHETIC DATA generated by make_synthetic_dataset(seed={seed}) in P09.21; "
            "a simulation, not a measurement")
    specs = [
        (SYNTHETIC_FILES[0], "motor_open", ["run", "t_s", "duty_cmd", "theta_deg", "note"], ("u", "theta")),
        (SYNTHETIC_FILES[1], "motor_closed", ["run", "t_s", "ref_deg", "theta_deg", "duty_cmd", "kp", "ki", "note"],
         ("ref", "theta", "u", "kp", "ki")),
        (SYNTHETIC_FILES[2], "servo_open", ["run", "t_s", "cmd_deg", "theta_deg", "note"], ("u", "theta")),
        (SYNTHETIC_FILES[3], "servo_closed", ["run", "t_s", "ref_deg", "theta_deg", "cmd_deg", "kp", "ki", "note"],
         ("ref", "theta", "u", "kp", "ki")),
    ]
    paths = []
    for fname, key, header, cols in specs:
        p = out_dir / fname
        with open(p, "w", newline="") as fh:
            fh.write(head + "\n")
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(header)
            for i, run in enumerate(ds[key], start=1):
                for k in range(run["t"].size):
                    w.writerow([i, f"{run['t'][k]:.3f}"] + [f"{run[c][k]:.4f}" for c in cols] + ["synthetic"])
        paths.append(p)
    return paths


def read_runs_csv(path: str | Path) -> list[dict[str, np.ndarray]]:
    """Read a lab CSV (comment lines starting with '#' skipped) into a list of runs, one dict of columns each."""
    with open(path) as fh:
        rows = [r for r in csv.reader(line for line in fh if not line.startswith("#")) if r]
    header, body = rows[0], rows[1:]
    cols: dict[str, np.ndarray] = {}
    for j, name in enumerate(header):
        col = [r[j] for r in body]
        try:
            cols[name] = np.array([float(v) for v in col])
        except ValueError:
            cols[name] = np.array(col)
    runs = []
    for rid in np.unique(cols["run"]):
        sel = cols["run"] == rid
        runs.append({k: v[sel] for k, v in cols.items()})
    return runs


# ---------------------------------------------------------------- the analysis (given orchestration)

def _closed_loop_runs(runs: list[dict], num: np.ndarray, den: np.ndarray, Ts: float, t_end: float,
                      u_lim: float, plant_delay: float, dz: float, dz_comp: float, quant: float,
                      out_col: str) -> list[dict]:
    """Predicted (continuous), measured (from the CSV) and sampled-model metrics for each closed-loop run."""
    out = []
    for run in runs:
        r, kp, ki = float(run["ref_deg"][-1]), float(run["kp"][0]), float(run["ki"][0])
        pred = predict_closed_loop(num, den, kp, ki, r, t_end)
        theta = 1.5 * Ts + plant_delay           # firmware sample + half a sample of hold + dead time
        pred_lag = predict_closed_loop(num, np.convolve(den, [theta, 1.0]), kp, ki, r, t_end)
        meas = step_metrics(run["t_s"], run["theta_deg"], ref=r)
        sim = simulate_closed_loop(num, den, kp, ki, r, Ts, t_end, u_min=-u_lim, u_max=u_lim, delay_samples=1,
                                   plant_delay=plant_delay, dz=dz, dz_comp=dz_comp, quant=quant, t_pre=T_PRE)
        simm = step_metrics(sim["t"], sim["y"], ref=r)
        out.append({"r": r, "kp": kp, "ki": ki, "data": run, "pred": pred, "pred_lag": pred_lag, "theta": theta,
                    "cmp_lag": compare_metrics(pred_lag["metrics"], meas), "meas": meas, "sim": sim,
                    "sim_metrics": simm, "cmp_pred": compare_metrics(pred["metrics"], meas),
                    "cmp_sim": compare_metrics(simm, meas), "peak_u": float(np.max(np.abs(run[out_col])))})
    return out


def analyse_dataset(data_dir: str | Path, files: tuple[str, str, str, str] = SYNTHETIC_FILES) -> dict:
    """Identify, design, predict, and compare for both options from the lab CSVs in data_dir.

    ``files`` names the motor open-loop, motor closed-loop, servo open-loop and servo closed-loop files
    (SYNTHETIC_FILES by default; MEASURED_FILES, or your own names, for your data). An option whose
    open-loop file is missing is skipped; a missing closed-loop file gives identification and design only.
    The firmware constants (Ts, limits, dead-zone compensation, quantisation, targets) are this module's;
    edit them to match your firmware. Uses only the learner-written functions and the CSVs: the true
    synthetic-bench numbers are never read. Returns a nested dict used by main() and by the tests.
    """
    data_dir = Path(data_dir)
    opt = lambda i: read_runs_csv(data_dir / files[i]) if (data_dir / files[i]).exists() else []  # noqa: E731
    mo, mc, so, sc = opt(0), opt(1), opt(2), opt(3)
    res: dict = {"motor": {}, "servo": {}}
    if mo:
        _analyse_motor(res["motor"], mo, mc)
    if so:
        _analyse_servo(res["servo"], so, sc)
    return res


def _analyse_motor(out: dict, mo: list[dict], mc: list[dict]) -> None:
    """Option B: fit every open-loop run, the dead zone from the slopes, then design, predict and compare."""
    fits = []
    for run in mo:
        u = float(run["duty_cmd"][-1])
        g = fit_integrator_loglinear(run["t_s"], run["theta_deg"], u)
        model = lambda tt, p, y0=g["y0"], u=u: y0 + integrator_step_response(tt, p[0], p[1], p[2], u)  # noqa: E731
        ls = refine_least_squares(model, run["t_s"], run["theta_deg"], np.array([g["K"], g["tau"], g["delay"]]),
                                  lower=np.array([0.0, 1e-4, 0.0]))
        fits.append({"u": u, "loglin": g, "ls": ls})
    duties = np.array([f["u"] for f in fits])
    slopes = np.array([f["ls"]["p"][0] * f["u"] for f in fits])
    K, udz = fit_dead_zone(duties, slopes)
    tau = float(np.mean([f["ls"]["p"][1] for f in fits]))
    L = float(np.mean([f["ls"]["p"][2] for f in fits]))
    num, den = np.array([K]), np.array([tau, 1.0, 0.0])
    des = design_p_position(K, tau, MOTOR_OVERSHOOT_TARGET)
    kp_sampled = design_p_sampled(num, den, MOTOR_CLOSED_RUNS[0][0], MOTOR_TS, MOTOR_OVERSHOOT_TARGET, 0.005,
                                  0.05, t_end=MOTOR_CLOSED_T_END, u_min=-MOTOR_DUTY_LIMIT, u_max=MOTOR_DUTY_LIMIT,
                                  delay_samples=1, plant_delay=L, dz=udz, dz_comp=MOTOR_DZ_COMP,
                                  quant=MOTOR_DEG_PER_COUNT)
    out.update(fits=fits, K=K, dead_zone=udz, tau=tau, delay=L, num=num, den=den, design=des,
               kp_sampled=kp_sampled,
               runs=_closed_loop_runs(mc, num, den, MOTOR_TS, MOTOR_CLOSED_T_END, MOTOR_DUTY_LIMIT, L, udz,
                                      MOTOR_DZ_COMP, MOTOR_DEG_PER_COUNT, "duty_cmd"))


def _analyse_servo(out: dict, so: list[dict], sc: list[dict]) -> None:
    """Option A: fit every open-loop step (the model is the first, small one; a large one shows the slew
    limit), then design, predict and compare."""
    sfits = []
    for run in so:
        c = float(run["cmd_deg"][-1])
        g = fit_first_order_loglinear(run["t_s"], run["theta_deg"], c)
        model = lambda tt, p, y0=g["y0"], c=c: y0 + first_order_step_response(tt, p[0], p[1], p[2], c)  # noqa: E731
        ls = refine_least_squares(model, run["t_s"], run["theta_deg"], np.array([g["K"], g["tau"], g["delay"]]),
                                  lower=np.array([0.0, 1e-4, 0.0]))
        sfits.append({"u": c, "loglin": g, "ls": ls})
    Ks, taus, Ls = (float(v) for v in sfits[0]["ls"]["p"])
    nums, dens = np.array([Ks]), np.array([taus, 1.0])
    sdes = design_pi_first_order(Ks, taus, SERVO_RISE_TARGET)
    out.update(fits=sfits, K=Ks, tau=taus, delay=Ls, num=nums, den=dens, design=sdes,
               runs=_closed_loop_runs(sc, nums, dens, SERVO_TS, SERVO_CLOSED_T_END, SERVO_CMD_LIMIT, Ls, 0.0, 0.0,
                                      SERVO_DEG_PER_COUNT, "cmd_deg"))


def _fmt_cmp(cmp: dict) -> list[str]:
    names = {"rise_time": "rise time 10-90 % (s)", "overshoot_pct": "overshoot (%)",
             "settling_time": "settling 2 % (s)", "ss_error": "steady-state error (deg)"}
    return [f"    {names[k]:26s} {p:9.4f} {m:9.4f} {d:+9.4f}" for k, (p, m, d, _) in cmp.items()]


def _print_runs(runs: list[dict], what: str, limit: float, unit: str) -> None:
    for run in runs:
        tag = f"step {run['r']:.0f} deg, Kp {run['kp']}" + (f", Ki {run['ki']}" if run["ki"] else "")
        print(f"  {tag}: predicted (continuous linear model) | SYNTHETIC | difference")
        print("\n".join(_fmt_cmp(run["cmp_pred"])))
        print(f"  {tag}: delay as a lag 1/(theta s + 1), theta = {run['theta'] * 1e3:.1f} ms | SYNTHETIC | difference")
        print("\n".join(_fmt_cmp(run["cmp_lag"])))
        print(f"  {tag}: sampled model (Ts, 1-sample delay, dead time, limit, dead zone, quantisation) | SYNTHETIC "
              "| difference")
        print("\n".join(_fmt_cmp(run["cmp_sim"])))
        print(f"    peak {what} in the SYNTHETIC run {run['peak_u']:.3f}{unit} (software limit {limit}{unit})")


def main(out_dir: str | Path = ".") -> dict:
    """Write the SYNTHETIC CSVs to out_dir, analyse them, print the report and write servo_feedback_lab.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    write_synthetic_csvs(out_dir, SEED)
    res = analyse_dataset(out_dir)
    m, s = res["motor"], res["servo"]
    print("SYNTHETIC DATASET (seeded simulation of the bench, not a measurement)")
    print("OPTION B: DC GEAR MOTOR, duty -> angle, model K e^{-Ls}/(s(tau s + 1)) with a dead zone")
    for f in m["fits"]:
        g, p, se = f["loglin"], f["ls"]["p"], f["ls"]["se"]
        print(f"  step duty {f['u']:.2f}: log-linear K_app {g['K']:.1f} deg/s, tau {g['tau'] * 1e3:.1f} ms, "
              f"L {g['delay'] * 1e3:.1f} ms | LM K_app {p[0]:.1f} +- {se[0]:.1f}, tau {p[1] * 1e3:.1f} +- "
              f"{se[1] * 1e3:.1f} ms, L {p[2] * 1e3:.1f} +- {se[2] * 1e3:.1f} ms, rms {f['ls']['rms']:.3f} deg")
    print(f"  dead-zone fit: K = {m['K']:.1f} deg/s per unit duty, dead zone {m['dead_zone']:.4f}; "
          f"mean tau {m['tau'] * 1e3:.2f} ms, mean L {m['delay'] * 1e3:.2f} ms")
    d = m["design"]
    print(f"  P design for {MOTOR_OVERSHOOT_TARGET:.0f} % overshoot: zeta {d['zeta']:.4f}, wn {d['wn']:.2f} rad/s, "
          f"Kp {d['kp']:.5f} duty/deg -> firmware Kp {MOTOR_KP}, dead-zone compensation {MOTOR_DZ_COMP}")
    kpk = MOTOR_KP * m["K"]
    print(f"  with the firmware Kp: zeta {1 / (2 * np.sqrt(m['tau'] * kpk)):.4f}, wn {np.sqrt(kpk / m['tau']):.2f} "
          f"rad/s, closed-loop poles {np.round(m['runs'][0]['pred']['poles'], 3)}")
    print(f"  redesign on the sampled model for {MOTOR_OVERSHOOT_TARGET:.0f} % overshoot: Kp {m['kp_sampled']:.5f} "
          f"-> firmware Kp {MOTOR_KP_REDESIGN} (run 3)")
    _print_runs(m["runs"], "duty", MOTOR_DUTY_LIMIT, "")
    print("  uncompensated dead-zone band (dz - compensation)/Kp, SYNTHETIC bench truth "
          f"(dz {MOTOR_TRUE['dz_pos']:.2f} forward, {MOTOR_TRUE['dz_neg']:.2f} back): " + "; ".join(
              f"Kp {kp}: {(MOTOR_TRUE['dz_pos'] - MOTOR_DZ_COMP) / kp:.2f} deg forward, "
              f"{(MOTOR_TRUE['dz_neg'] - MOTOR_DZ_COMP) / kp:.2f} deg back" for kp in (MOTOR_KP, MOTOR_KP_REDESIGN)))
    print("OPTION A: HOBBY SERVO, command -> wiper angle, model K e^{-Ls}/(tau s + 1)")
    for f in s["fits"]:
        g, p, se = f["loglin"], f["ls"]["p"], f["ls"]["se"]
        print(f"  step {f['u']:.0f} deg: log-linear K {g['K']:.4f}, tau {g['tau'] * 1e3:.1f} ms, L {g['delay'] * 1e3:.1f} ms"
              f" | LM K {p[0]:.4f} +- {se[0]:.4f}, tau {p[1] * 1e3:.1f} +- {se[1] * 1e3:.1f} ms, "
              f"L {p[2] * 1e3:.1f} +- {se[2] * 1e3:.1f} ms, rms {f['ls']['rms']:.3f} deg")
    d = s["design"]
    print(f"  PI design for rise time {SERVO_RISE_TARGET} s: tau_cl {d['tau_cl'] * 1e3:.1f} ms, Kp {d['kp']:.4f}, "
          f"Ki {d['ki']:.3f} 1/s -> firmware Kp {SERVO_KP}, Ki {SERVO_KI}")
    print(f"  closed-loop poles {np.round(s['runs'][0]['pred']['poles'], 3)}, "
          f"zero {-SERVO_KI / SERVO_KP:.3f}")
    _print_runs(s["runs"], "command", SERVO_CMD_LIMIT, " deg")

    # Section 5, "Is the remaining difference real?": ten seeds, each identified afresh (SYNTHETIC)
    import tempfile
    seeds = range(15, 25)
    os_bench, os_model, os_sym = [], [], []
    with tempfile.TemporaryDirectory() as tmp:
        for sd in seeds:
            write_synthetic_csvs(tmp, sd)
            run0 = analyse_dataset(tmp)["motor"]["runs"][0]
            os_bench.append(run0["meas"]["overshoot_pct"])
            os_model.append(run0["sim_metrics"]["overshoot_pct"])
    dz_neg = MOTOR_TRUE["dz_neg"]
    try:
        MOTOR_TRUE["dz_neg"] = MOTOR_TRUE["dz_pos"]          # the falsification test: a symmetric dead zone
        for sd in seeds:
            rec = make_synthetic_dataset(sd)["motor_closed"][0]
            os_sym.append(step_metrics(rec["t"], rec["theta"], ref=MOTOR_CLOSED_RUNS[0][0])["overshoot_pct"])
    finally:
        MOTOR_TRUE["dz_neg"] = dz_neg
    ms = lambda v: f"{np.mean(v):.2f} +- {np.std(v, ddof=1):.2f} %"  # noqa: E731
    print(f"SEEDS 15-24 (SYNTHETIC), Option B {MOTOR_CLOSED_RUNS[0][0]:.0f} deg step, Kp {MOTOR_KP}: bench overshoot "
          f"{ms(os_bench)}, sampled-model prediction {ms(os_model)}, mean difference "
          f"{np.mean(np.subtract(os_model, os_bench)):+.2f} points; same bench and noise with a symmetric dead zone "
          f"{MOTOR_TRUE['dz_pos']:.2f}: {ms(os_sym)}")

    fig, ax = plt.subplots(2, 2, figsize=(13, 8))
    for f, run in zip(m["fits"], read_runs_csv(out_dir / SYNTHETIC_FILES[0])):
        tt = run["t_s"]
        ax[0, 0].plot(tt, run["theta_deg"], ".", ms=3, label=f"SYNTHETIC duty {f['u']:.2f}")
        p = f["ls"]["p"]
        ax[0, 0].plot(tt, f["loglin"]["y0"] + integrator_step_response(tt, p[0], p[1], p[2], f["u"]), "k-", lw=0.8)
    ax[0, 0].set_title("Option B open-loop steps (SYNTHETIC) and fitted K/(s(tau s+1))")
    ax[0, 0].set_xlabel("t (s)"); ax[0, 0].set_ylabel("angle (deg)"); ax[0, 0].legend(fontsize=8)
    for f, run in zip(s["fits"], read_runs_csv(out_dir / SYNTHETIC_FILES[2])):
        tt = run["t_s"]
        ax[1, 0].plot(tt, run["theta_deg"], ".", ms=3, label=f"SYNTHETIC command {f['u']:.0f} deg")
        p = f["ls"]["p"]
        ax[1, 0].plot(tt, f["loglin"]["y0"] + first_order_step_response(tt, p[0], p[1], p[2], f["u"]), "k-", lw=0.8)
    ax[1, 0].set_title("Option A open-loop steps (SYNTHETIC) and fitted K e^{-Ls}/(tau s+1)")
    ax[1, 0].set_xlabel("t (s)"); ax[1, 0].set_ylabel("angle (deg)"); ax[1, 0].legend(fontsize=8)
    for axx, runs in ((ax[0, 1], m["runs"][:1] + m["runs"][2:]), (ax[1, 1], s["runs"][:1])):
        for run in runs:
            lab = f"Kp {run['kp']}"
            axx.plot(run["data"]["t_s"], run["data"]["theta_deg"], ".", ms=3, label=f"SYNTHETIC closed loop, {lab}")
            axx.plot(run["pred"]["t"], run["pred"]["y"], "-", label=f"predicted (continuous model), {lab}")
            axx.step(run["sim"]["t"], run["sim"]["y"], where="post", lw=0.8, label=f"sampled model, {lab}")
        r = runs[0]["r"]
        axx.axhline(r * 1.02, color="gray", lw=0.5); axx.axhline(r * 0.98, color="gray", lw=0.5)
        axx.set_xlabel("t (s)"); axx.set_ylabel("angle (deg)"); axx.legend(fontsize=7)
    ax[0, 1].set_title(f"Option B closed loop, {MOTOR_CLOSED_RUNS[0][0]:.0f} deg step (SYNTHETIC)")
    ax[1, 1].set_title(f"Option A closed loop, {SERVO_REF_STEPS[0]:.0f} deg step, PI (SYNTHETIC)")
    fig.tight_layout()
    fig.savefig(out_dir / "servo_feedback_lab.png", dpi=110)
    plt.close(fig)
    return res


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
