"""P09.04 Time Response & Poles: reference solution.

A step-response analyzer. From a sampled response it measures the steady-state value, the 10-90 % rise
time, the peak and peak time, the percentage overshoot (and undershoot) and the 2 % settling time, each
interpolated between samples. Independently it computes the same quantities from the pole locations: the
second-order zeta-omega_n formulas, the exact second-order closed form, and the step response of any TF
with simple poles by partial fractions (residues), which is checked against simulation. Then it asks when
a higher-order system is "second order enough": dominant poles, the effect of a third pole, and the
effect of a zero (left half plane: more overshoot; right half plane: undershoot).

The example is a control-surface position servo: a DC motor and gearbox K_m/(s (tau s + 1)) under
proportional feedback, the plant used again in the home lab (P09.21), and the pitch-attitude loop of
P09.03.

Conventions:
    * TF = (num, den) arrays, highest power first; time in s; steps are unit steps from rest.
    * Second order: omega_n^2 / (s^2 + 2 zeta omega_n s + omega_n^2), unit DC gain.

Run:  python solution.py          (prints the servo design and the dominance studies, writes time_response.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-09/lesson-04-time-response-poles
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np

TF = tuple[np.ndarray, np.ndarray]


# ---------------------------------------------------------------- given: TF realisation and simulation (P09.03)

def tf(num, den) -> TF:
    """Return (num, den) as float arrays with den monic (leading zeros of num removed)."""
    num = np.trim_zeros(np.atleast_1d(np.asarray(num, dtype=float)), "f")
    num = num if num.size else np.array([0.0])
    den = np.trim_zeros(np.atleast_1d(np.asarray(den, dtype=float)), "f")
    return num / den[0], den / den[0]


def tf_to_ss(G: TF) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Return (A, B, C, D) in controllable canonical form for a proper TF (as built in P09.03)."""
    num, den = tf(*G)
    n = len(den) - 1
    num = np.concatenate([np.zeros(n + 1 - len(num)), num])
    D = num[0]
    rem = num - D * den
    A = np.zeros((n, n))
    if n:
        A[:-1, 1:] = np.eye(n - 1)
        A[-1, :] = -den[1:][::-1]
    B = np.zeros(n)
    if n:
        B[-1] = 1.0
    return A, B, rem[1:][::-1].copy(), float(D)


def simulate_tf(G: TF, u: Callable[[float], float], t: np.ndarray, substeps: int = 4) -> np.ndarray:
    """Return y(t_k) for G driven by u(t) from rest, by fixed-step RK4 (as built in P09.03)."""
    A, B, C, D = tf_to_ss(G)
    x = np.zeros(A.shape[0])
    y = np.empty(len(t))
    y[0] = C @ x + D * u(t[0])
    f = lambda tt, xx: A @ xx + B * u(tt)
    for k in range(len(t) - 1):
        h = (t[k + 1] - t[k]) / substeps
        tt = t[k]
        for _ in range(substeps):
            k1 = f(tt, x)
            k2 = f(tt + h / 2, x + h / 2 * k1)
            k3 = f(tt + h / 2, x + h / 2 * k2)
            k4 = f(tt + h, x + h * k3)
            x = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            tt += h
        y[k + 1] = C @ x + D * u(t[k + 1])
    return y


def simulate_step(G: TF, t: np.ndarray) -> np.ndarray:
    """Return the unit-step response of G on the grid t by RK4."""
    return simulate_tf(G, lambda tt: 1.0, t)


def second_order_tf(zeta: float, wn: float) -> TF:
    """Return wn^2 / (s^2 + 2 zeta wn s + wn^2)."""
    return tf([wn * wn], [1.0, 2 * zeta * wn, wn * wn])


# ---------------------------------------------------------------- measuring a sampled step response

def _crossing(t: np.ndarray, y: np.ndarray, level: float) -> float:
    """Return the first time y crosses ``level`` upward, by linear interpolation between samples."""
    i = int(np.argmax(y >= level))
    if y[i] < level:
        return float("nan")
    if i == 0:
        return float(t[0])
    return float(t[i - 1] + (level - y[i - 1]) * (t[i] - t[i - 1]) / (y[i] - y[i - 1]))


def step_metrics(t: np.ndarray, y: np.ndarray, y_final: float | None = None, band: float = 0.02) -> dict[str, float]:
    """Return the step-response metrics of a sampled response y(t) to a step applied at t[0].

    Keys: ``steady_state`` (``y_final`` if given, else the last sample), ``rise_time`` (10 % to 90 % of the
    steady state, each crossing interpolated), ``peak`` and ``peak_time`` (largest sample), ``overshoot``
    (percent above the steady state, 0 if none), ``undershoot`` (percent of the steady state that y goes
    below zero, 0 if none) and ``settling_time`` (the time after which |y - y_final| stays within
    ``band`` |y_final|, interpolated on the last exit from the band; nan if never inside at the end).
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    yf = float(y[-1]) if y_final is None else float(y_final)
    sgn = 1.0 if yf >= 0 else -1.0
    ys, yfs = sgn * y, sgn * yf
    rise = _crossing(t, ys, 0.9 * yfs) - _crossing(t, ys, 0.1 * yfs)
    i_pk = int(np.argmax(ys))
    overshoot = max(0.0, 100.0 * (ys[i_pk] - yfs) / yfs)
    undershoot = max(0.0, -100.0 * ys.min() / yfs)
    tol = band * abs(yf)
    out = np.abs(y - yf) > tol
    if out[-1]:
        settle = float("nan")
    elif not out.any():
        settle = float(t[0])
    else:
        j = int(np.nonzero(out)[0][-1])               # last sample outside the band
        e0, e1 = abs(y[j] - yf), abs(y[j + 1] - yf)
        settle = float(t[j] + (e0 - tol) * (t[j + 1] - t[j]) / (e0 - e1))
    return {"steady_state": yf, "rise_time": float(rise), "peak": float(y[i_pk]), "peak_time": float(t[i_pk]),
            "overshoot": float(overshoot), "undershoot": float(undershoot), "settling_time": settle}


# ---------------------------------------------------------------- second-order formulas

def second_order_params(den: np.ndarray) -> tuple[float, float]:
    """Return (zeta, omega_n) of a quadratic den = a s^2 + b s + c with a, c > 0: omega_n = sqrt(c/a), zeta = b/(2 sqrt(a c))."""
    a, b, c = (float(v) for v in den)
    return b / (2.0 * np.sqrt(a * c)), float(np.sqrt(c / a))


def second_order_step(t: np.ndarray, zeta: float, wn: float) -> np.ndarray:
    """Return the exact unit-step response of wn^2/(s^2 + 2 zeta wn s + wn^2) for zeta >= 0 (all three cases)."""
    t = np.asarray(t, dtype=float)
    if zeta < 1.0:
        wd = wn * np.sqrt(1.0 - zeta * zeta)
        phi = np.arccos(zeta)
        return 1.0 - np.exp(-zeta * wn * t) * np.sin(wd * t + phi) / np.sqrt(1.0 - zeta * zeta)
    if zeta == 1.0:
        return 1.0 - np.exp(-wn * t) * (1.0 + wn * t)
    r = wn * np.sqrt(zeta * zeta - 1.0)
    p1, p2 = -zeta * wn + r, -zeta * wn - r           # both negative, p1 the slower
    return 1.0 + (p2 * np.exp(p1 * t) - p1 * np.exp(p2 * t)) / (p1 - p2)


def overshoot_formula(zeta: float) -> float:
    """Return the percentage overshoot 100 exp(-pi zeta / sqrt(1 - zeta^2)) of the standard second order (0 for zeta >= 1)."""
    if zeta >= 1.0:
        return 0.0
    return float(100.0 * np.exp(-np.pi * zeta / np.sqrt(1.0 - zeta * zeta)))


def peak_time_formula(zeta: float, wn: float) -> float:
    """Return the peak time pi / omega_d with omega_d = omega_n sqrt(1 - zeta^2) (inf for zeta >= 1)."""
    if zeta >= 1.0:
        return float("inf")
    return float(np.pi / (wn * np.sqrt(1.0 - zeta * zeta)))


def settling_time_envelope(zeta: float, wn: float, band: float = 0.02) -> float:
    """Return the time at which the decay envelope e^{-zeta wn t}/sqrt(1 - zeta^2) falls to ``band``: -ln(band sqrt(1 - zeta^2))/(zeta wn).

    An upper bound on the true settling time for 0 < zeta < 1; the familiar 4/(zeta wn) is this with the
    sqrt(1 - zeta^2) dropped and ln 50 = 3.9 rounded to 4.
    """
    return float(-np.log(band * np.sqrt(1.0 - zeta * zeta)) / (zeta * wn))


def rise_time_exact(zeta: float, wn: float) -> float:
    """Return the exact 10-90 % rise time of the standard second order by bisection on ``second_order_step``."""
    def cross(level: float) -> float:
        hi = 1.0 / wn
        while second_order_step(np.array([hi]), zeta, wn)[0] < level:
            hi *= 2.0
        lo = 0.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if second_order_step(np.array([mid]), zeta, wn)[0] < level:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)
    return cross(0.9) - cross(0.1)


def zeta_from_overshoot(overshoot_pct: float) -> float:
    """Return the damping ratio that gives the stated percentage overshoot: zeta = -ln M / sqrt(pi^2 + ln^2 M), M = overshoot/100."""
    L = np.log(overshoot_pct / 100.0)
    return float(-L / np.sqrt(np.pi ** 2 + L * L))


# ---------------------------------------------------------------- partial fractions

def step_residues(G: TF, tol: float = 1e-6) -> tuple[np.ndarray, np.ndarray]:
    """Return (p, r): the roots of D(s) = s den(s) and the residues r_k = N(p_k)/D'(p_k) of Y(s) = G(s)/s.

    The residue at p = 0 is the steady-state value G(0); the others say how much of each mode the step
    excites. Raise ``ValueError`` if two roots coincide within ``tol`` (repeated poles need the derivative
    formula, not covered here).
    """
    num, den = tf(*G)
    D = np.polymul(den, [1.0, 0.0])
    p = np.roots(D)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if abs(p[i] - p[j]) <= tol * max(1.0, abs(p[i])):
                raise ValueError("repeated poles: partial fractions with simple poles do not apply")
    r = np.polyval(num, p) / np.polyval(np.polyder(D), p)
    return p, r


def residue_step_response(G: TF, t: np.ndarray, tol: float = 1e-6) -> np.ndarray:
    """Return the unit-step response of a proper G with simple poles by partial fractions.

    Y(s) = G(s)/s = N(s)/D(s) with D(s) = s den(s); for simple roots p_k of D, y(t) = sum_k r_k e^{p_k t}
    with r_k = N(p_k)/D'(p_k) from ``step_residues``. (For a biproper G the direct term is included
    automatically: it sits in the residue at s = 0 and in the others.) Raise ``ValueError`` on repeated poles.
    """
    p, r = step_residues(G, tol)
    t = np.asarray(t, dtype=float)
    return np.real(np.exp(np.outer(t, p)) @ r)


# ---------------------------------------------------------------- dominance, zeros, the servo

def dominant_poles(G: TF) -> tuple[np.ndarray, float]:
    """Return (the dominant pole or conjugate pair, dominance ratio).

    The dominant poles are those with the largest real part (closest to the imaginary axis); the ratio is
    Re(next pole) / Re(dominant), so a ratio of 5 or more is the usual "second order is enough" rule.
    Assumes a stable G; zeros near the dominant poles can still spoil the approximation.
    """
    p = np.roots(tf(*G)[1])
    p = p[np.argsort(-p.real, kind="stable")]
    lead = p[0]
    k = 2 if abs(lead.imag) > 1e-12 * max(1.0, abs(lead)) else 1
    dom = p[:k]
    ratio = float(p[k].real / dom[0].real) if len(p) > k else float("inf")
    return dom, ratio


def with_third_pole(zeta: float, wn: float, p3: float) -> TF:
    """Return wn^2 p3 / ((s + p3)(s^2 + 2 zeta wn s + wn^2)): the standard second order times a unit-DC-gain real pole at -p3."""
    return tf([wn * wn * p3], np.polymul([1.0, p3], [1.0, 2 * zeta * wn, wn * wn]))


def with_zero(zeta: float, wn: float, z: float) -> TF:
    """Return wn^2 (1 - s/z) / (s^2 + 2 zeta wn s + wn^2): a zero at s = z (z < 0 left half plane, z > 0 right half plane), unit DC gain."""
    return tf([-wn * wn / z, wn * wn], [1.0, 2 * zeta * wn, wn * wn])


def zero_step_closed_form(t: np.ndarray, zeta: float, wn: float, z: float) -> np.ndarray:
    """Return the step response of ``with_zero``: y = y2 - (1/z) dy2/dt, with dy2/dt the impulse response of the standard second order (zeta < 1)."""
    t = np.asarray(t, dtype=float)
    wd = wn * np.sqrt(1.0 - zeta * zeta)
    dy2 = wn / np.sqrt(1.0 - zeta * zeta) * np.exp(-zeta * wn * t) * np.sin(wd * t)
    return second_order_step(t, zeta, wn) - dy2 / z


def servo_closed_loop(K: float, tau: float) -> TF:
    """Return the unity-feedback position loop around K/(s (tau s + 1)): K / (tau s^2 + s + K)."""
    return tf([K], [tau, 1.0, K])


def servo_gain_for_overshoot(overshoot_pct: float, tau: float) -> float:
    """Return the loop gain K [1/s] that gives the stated overshoot: zeta = 1/(2 sqrt(K tau)), so K = 1/(4 zeta^2 tau)."""
    z = zeta_from_overshoot(overshoot_pct)
    return 1.0 / (4.0 * z * z * tau)


def pitch_closed_loop() -> TF:
    """Return the pitch-attitude closed loop theta/theta_c of P09.03 (K_q = 0.3, K_theta = 3), coefficients as reduced there."""
    return tf([720.0, 73080.0, 108000.0], [1.0, 123.2, 2400.0, 16240.0, 115880.0, 108000.0])


# ---------------------------------------------------------------- the example

def main(out_dir: str | Path = ".") -> None:
    """Print the servo design, the formula checks and the dominance studies; save time_response.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tau = 0.05
    print(f"elevator position servo K/(s (tau s + 1)), tau = {tau * 1e3:.0f} ms, proportional feedback")
    t = np.linspace(0.0, 1.5, 3001)
    for os_target in (5.0, 10.0, 20.0):
        K = servo_gain_for_overshoot(os_target, tau)
        G = servo_closed_loop(K, tau)
        z, wn = second_order_params(G[1])
        m = step_metrics(t, simulate_step(G, t), 1.0)
        print(f"  target {os_target:4.1f} %: K = {K:6.2f} 1/s, zeta = {z:.4f}, wn = {wn:.3f} rad/s, zeta wn = {z * wn:.2f} 1/s")
        print(f"    measured: overshoot {m['overshoot']:.3f} %, peak time {m['peak_time']:.4f} s, rise {m['rise_time']:.4f} s,"
              f" settling {m['settling_time']:.4f} s")
        print(f"    formulas: overshoot {overshoot_formula(z):.3f} %, peak time {peak_time_formula(z, wn):.4f} s,"
              f" rise (exact) {rise_time_exact(z, wn):.4f} s, envelope settling {settling_time_envelope(z, wn):.4f} s, 4/(zeta wn) {4 / (z * wn):.4f} s")
    Ks = np.linspace(5.0, 60.0, 551)
    tf_ = np.linspace(0.0, 1.5, 30001)
    ts_all = np.array([step_metrics(tf_, second_order_step(tf_, *second_order_params(servo_closed_loop(K, tau)[1])), 1.0)["settling_time"] for K in Ks])
    i = int(np.nanargmin(ts_all))
    ok = Ks[ts_all <= 0.3]
    print(f"  K with 2 % settling <= 0.3 s: {ok.min():.1f} to {ok.max():.1f} 1/s (overshoot"
          f" {overshoot_formula(second_order_params(servo_closed_loop(ok.min(), tau)[1])[0]):.2f} to"
          f" {overshoot_formula(second_order_params(servo_closed_loop(ok.max(), tau)[1])[0]):.2f} %)")
    print(f"  best 2 % settling over K in [5, 60]: {ts_all[i]:.4f} s at K = {Ks[i]:.1f} (overshoot "
          f"{overshoot_formula(second_order_params(servo_closed_loop(Ks[i], tau)[1])[0]):.2f} %): Re(poles) = -1/(2 tau) = {-1 / (2 * tau):.0f} 1/s for every K above {1 / (4 * tau):.0f}")

    # partial fractions against simulation on the pitch loop of P09.03
    Gp = pitch_closed_loop()
    tp = np.linspace(0.0, 6.0, 1201)
    y_pf, y_sim = residue_step_response(Gp, tp), simulate_step(Gp, tp)
    mp = step_metrics(tp, y_pf, 1.0)
    dom, ratio = dominant_poles(Gp)
    pp, rr = step_residues(Gp)
    print(f"\npitch loop of P09.03: partial fractions vs RK4: max difference {np.max(np.abs(y_pf - y_sim)):.2e}")
    for pk, rk in sorted(zip(pp, rr), key=lambda a: -a[0].real):
        print(f"  pole {pk.real:+9.3f} {pk.imag:+8.3f}j   residue |r| = {abs(rk):.4f}")
    print(f"  slowest pole {np.round(dom, 3)} (ratio to the next {ratio:.3f}); its residue is comparable with the pair's:"
          f" no single pole or pair dominates (the zero at -1.5 sits between them)")
    print(f"  measured overshoot {mp['overshoot']:.2f} %, peak time {mp['peak_time']:.3f} s, settling {mp['settling_time']:.3f} s")
    pair = pp[np.argmax(np.abs(pp.imag))]
    z_d, wn_d = second_order_params(np.real(np.poly([pair, np.conj(pair)])))
    print(f"  complex pair alone: zeta = {z_d:.3f}, wn = {wn_d:.3f}: overshoot {overshoot_formula(z_d):.2f} %, peak time {peak_time_formula(z_d, wn_d):.3f} s,"
          f" envelope settling {settling_time_envelope(z_d, wn_d):.3f} s")

    # third pole
    zeta, wn = 0.5, 1.0
    print(f"\nthird pole at -p3 added to zeta = {zeta}, wn = {wn} (zeta wn = {zeta * wn}): overshoot and settling against p3/(zeta wn)")
    t3 = np.linspace(0.0, 30.0, 30001)
    base = step_metrics(t3, second_order_step(t3, zeta, wn), 1.0)
    print(f"  second order alone: overshoot {base['overshoot']:.2f} %, settling {base['settling_time']:.3f} s")
    for ratio3 in (0.5, 1.0, 2.0, 5.0, 10.0):
        m3 = step_metrics(t3, residue_step_response(with_third_pole(zeta, wn, ratio3 * zeta * wn), t3), 1.0)
        print(f"  ratio {ratio3:4.1f}: overshoot {m3['overshoot']:6.2f} %, settling {m3['settling_time']:6.3f} s")

    # zero
    print(f"\nzero at s = z added to zeta = {zeta}, wn = {wn}:")
    for zz in (-0.5, -1.0, -2.0, -5.0, 2.0, 0.5):
        y = zero_step_closed_form(t3, zeta, wn, zz)
        mz = step_metrics(t3, y, 1.0)
        print(f"  z = {zz:+5.1f}: overshoot {mz['overshoot']:6.2f} %, undershoot {mz['undershoot']:6.2f} %, peak time {mz['peak_time']:.3f} s")

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    for z_ in (0.0, 0.2, 0.5, 0.7, 1.0, 2.0):
        axes[0].plot(t3[t3 < 15], second_order_step(t3[t3 < 15], z_, 1.0), label=f"zeta = {z_}")
    axes[0].set_xlabel("omega_n t"); axes[0].set_ylabel("y"); axes[0].legend(fontsize=7); axes[0].grid(lw=0.4)
    axes[0].set_title("second-order step responses", fontsize=9)
    for ratio3 in (0.5, 1.0, 5.0):
        axes[1].plot(t3[t3 < 15], residue_step_response(with_third_pole(zeta, wn, ratio3 * zeta * wn), t3[t3 < 15]), label=f"p3 = {ratio3} zeta wn")
    axes[1].plot(t3[t3 < 15], second_order_step(t3[t3 < 15], zeta, wn), "k--", label="no third pole")
    axes[1].set_xlabel("t [s]"); axes[1].legend(fontsize=7); axes[1].grid(lw=0.4)
    axes[1].set_title("a third real pole", fontsize=9)
    for zz in (-0.5, -2.0, 2.0, 0.5):
        axes[2].plot(t3[t3 < 15], zero_step_closed_form(t3[t3 < 15], zeta, wn, zz), label=f"zero at {zz:+}")
    axes[2].plot(t3[t3 < 15], second_order_step(t3[t3 < 15], zeta, wn), "k--", label="no zero")
    axes[2].set_xlabel("t [s]"); axes[2].legend(fontsize=7); axes[2].grid(lw=0.4)
    axes[2].set_title("a zero (RHP zeros undershoot)", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "time_response.png", dpi=150)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'time_response.png'}")


if __name__ == "__main__":
    import sys

    main(sys.argv[1] if len(sys.argv) > 1 else ".")
