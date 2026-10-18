"""P09.05 Frequency Response & Bode Plots: reference solution.

A Bode plotter from scratch. The exact magnitude (dB) and phase (deg) come from G(j omega); the phase is
unwrapped from scratch and, independently, built as the sum of the angles of the individual pole and zero
factors, which is continuous by construction. The asymptotic (straight-line) Bode diagram is built from
the Bode form K_B (j omega)^(-type) prod(1 + j omega/omega_b)^(+-1) ..., with the classic corrections at
the corners (-3.01 dB at a real pole, -20 log10(2 zeta) at a complex pair). The exact curves are checked
against a simulated sine test: drive the model with a sinusoid, wait for the transient to die, fit the
steady-state amplitude and phase by least squares. ``scipy.signal.bode`` is used only as an oracle in the
tests.

The examples are the elevator position servo of P09.04 (closed loop, 10 % overshoot design), the RC
anti-aliasing front end of a MEMS gyro from P09.02, and the pitch-attitude loop of P09.03.

Conventions:
    * TF = (num, den), highest power first; omega in rad/s, f = omega / (2 pi) in Hz.
    * Magnitude in dB = 20 log10 |G|; phase in degrees, unwrapped so that it starts at the physical
      low-frequency value (for example -90 deg for a type-1 system).
    * The sine test is a SIMULATION of a frequency-response test, not a measurement.

Run:  python solution.py          (prints the servo, gyro filter and pitch-loop results, writes bode.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-09/lesson-05-frequency-response-bode
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np

TF = tuple[np.ndarray, np.ndarray]


# ---------------------------------------------------------------- given: TF basics and simulation (P09.03)

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


def simulate_tf(G: TF, u: Callable[[np.ndarray], np.ndarray], t: np.ndarray) -> np.ndarray:
    """Return y(t_k) for G driven by u(t) from rest, by fixed-step RK4 with the grid step (as built in P09.03)."""
    A, B, C, D = tf_to_ss(G)
    x = np.zeros(A.shape[0])
    y = np.empty(len(t))
    ut = u(np.asarray(t, dtype=float))
    y[0] = C @ x + D * ut[0]
    for k in range(len(t) - 1):
        h = t[k + 1] - t[k]
        um = u(np.array([t[k] + h / 2]))[0]
        k1 = A @ x + B * ut[k]
        k2 = A @ (x + h / 2 * k1) + B * um
        k3 = A @ (x + h / 2 * k2) + B * um
        k4 = A @ (x + h * k3) + B * ut[k + 1]
        x = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        y[k + 1] = C @ x + D * ut[k + 1]
    return y


def servo_closed_loop(K: float = 14.31, tau: float = 0.05) -> TF:
    """Return the P09.04 elevator servo loop K/(tau s^2 + s + K) (K = 14.31 1/s gives 10 % overshoot)."""
    return tf([K], [tau, 1.0, K])


def rc_filter(tau: float = 1.6e-3, stages: int = 1) -> TF:
    """Return 1/(tau s + 1)^stages: the gyro RC front end of P09.02 (buffered stages)."""
    den = np.array([1.0])
    for _ in range(stages):
        den = np.polymul(den, [tau, 1.0])
    return tf([1.0], den)


def pitch_open_loop() -> TF:
    """Return the attitude-loop return ratio L(s) = K_theta theta/u of P09.03 (K_q = 0.3, K_theta = 3)."""
    return tf([720.0, 73080.0, 108000.0], [1.0, 123.2, 2400.0, 15520.0, 42800.0, 0.0])


# ---------------------------------------------------------------- exact frequency response

def freq_response(G: TF, w: np.ndarray) -> np.ndarray:
    """Return G(j omega) at the frequencies w [rad/s] by evaluating num and den on the imaginary axis."""
    s = 1j * np.asarray(w, dtype=float)
    return np.polyval(G[0], s) / np.polyval(G[1], s)


def mag_db(H: np.ndarray) -> np.ndarray:
    """Return 20 log10 |H| [dB]."""
    return 20.0 * np.log10(np.abs(H))


def unwrap_deg(phase_deg: np.ndarray) -> np.ndarray:
    """Return the phase with every jump larger than 180 deg between neighbours removed by adding a multiple of 360 deg."""
    p = np.asarray(phase_deg, dtype=float).copy()
    d = np.diff(p)
    corr = -360.0 * np.round(d / 360.0)
    p[1:] += np.cumsum(corr)
    return p


def phase_from_factors(G: TF, w: np.ndarray, tol: float = 1e-9) -> np.ndarray:
    """Return the phase [deg] as the sum of the angles of the Bode-form factors of G.

    Write G(s) = C s^(n0) prod(1 - s/z_i) / prod(1 - s/p_i) over the non-zero zeros z and poles p, with n0
    = (zeros at 0) - (poles at 0) and the real constant C = k prod(-z_i)/prod(-p_i). The phase is
    (180 if C < 0 else 0) + 90 n0 + sum angle(1 - j w/z_i) - sum angle(1 - j w/p_i). Each factor starts at
    0 deg at omega = 0 and moves continuously (no factor can cross the negative real axis unless its root
    is on the imaginary axis), so the sum is the physical, unwrapped phase with no 360-degree ambiguity.
    """
    num, den = tf(*G)
    s = 1j * np.asarray(w, dtype=float)
    z = np.roots(num) if len(num) > 1 else np.array([], dtype=complex)
    p = np.roots(den) if len(den) > 1 else np.array([], dtype=complex)
    C = complex(num[0] / den[0])
    n0 = 0
    ph = np.zeros(s.shape)
    for roots, sgn in ((z, +1), (p, -1)):
        for r in roots:
            if abs(r) <= tol:
                n0 += sgn
                continue
            C *= (-r) ** sgn
            ph += sgn * np.degrees(np.angle(1.0 - s / r))
    return ph + 90.0 * n0 + (180.0 if C.real < 0 else 0.0)


def bode(G: TF, w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (magnitude [dB], phase [deg]) of G at w: the phase of G(j w) unwrapped, then shifted by the
    multiple of 360 deg that makes its first value agree with ``phase_from_factors``."""
    H = freq_response(G, w)
    ph = unwrap_deg(np.degrees(np.angle(H)))
    ref = phase_from_factors(G, np.asarray(w, dtype=float)[:1])[0]
    ph += 360.0 * np.round((ref - ph[0]) / 360.0)
    return mag_db(H), ph


# ---------------------------------------------------------------- asymptotic (straight-line) Bode

def bode_form(G: TF, tol: float = 1e-9) -> dict:
    """Return the Bode form of G: K_B, type, and the corner list.

    G(s) = K_B s^(-type) prod (1 -+ s/w_b)^(+-1) prod (1 + 2 zeta s/w_n + s^2/w_n^2)^(+-1). Keys:
    ``K`` (K_B, the gain of the low-frequency asymptote at omega = 1 rad/s when type != 0, the DC gain
    when type = 0), ``type`` (poles minus zeros at s = 0), ``real`` (list of (w_b, +1 zero / -1 pole,
    lhp: bool)), ``complex`` (list of (w_n, zeta, +1 / -1)). Complex pairs with zeta <= 0 are reported with
    zeta = |Re| / w_n sign included (negative for a right-half-plane pair).
    """
    num, den = tf(*G)
    z = np.roots(num) if len(num) > 1 else np.array([], dtype=complex)
    p = np.roots(den) if len(den) > 1 else np.array([], dtype=complex)
    out = {"type": 0, "real": [], "complex": []}
    K = num[0] / den[0]
    for roots, sgn in ((z, +1), (p, -1)):
        done = np.zeros(len(roots), dtype=bool)
        for i, r in enumerate(roots):
            if done[i]:
                continue
            done[i] = True
            mag = abs(r)
            if mag <= tol:
                out["type"] -= sgn                   # a pole at 0 raises the type, a zero lowers it
                continue
            if abs(r.imag) <= tol * max(1.0, mag):
                out["real"].append((mag, sgn, bool(r.real < 0)))
                K *= (-r.real) ** sgn                # (s - r) = -r (1 - s/r)
            else:
                j = next(k for k in range(len(roots)) if not done[k] and abs(roots[k] - np.conj(r)) <= 1e-6 * mag)
                done[j] = True
                out["complex"].append((mag, float(-r.real / mag), sgn))
                K *= mag ** (2 * sgn)                # (s - r)(s - r*) = w_n^2 (1 + 2 zeta s/w_n + s^2/w_n^2)
    out["K"] = float(np.real(K))
    return out


def asymptotic_mag_db(G: TF, w: np.ndarray) -> np.ndarray:
    """Return the straight-line magnitude [dB]: 20 log|K_B| - 20 type log w, then +-20 dB/decade above each
    real corner and +-40 dB/decade above each complex-pair corner (zeros +, poles -; RHP same as LHP)."""
    bf = bode_form(G)
    lw = np.log10(np.asarray(w, dtype=float))
    m = 20.0 * np.log10(abs(bf["K"])) - 20.0 * bf["type"] * lw
    for wb, sgn, _ in bf["real"]:
        m += sgn * 20.0 * np.maximum(0.0, lw - np.log10(wb))
    for wn, _, sgn in bf["complex"]:
        m += sgn * 40.0 * np.maximum(0.0, lw - np.log10(wn))
    return m


def asymptotic_phase_deg(G: TF, w: np.ndarray) -> np.ndarray:
    """Return the straight-line phase [deg].

    Start at -90 type (+180 if K_B < 0). A real LHP zero adds 0 below w_b/10, +90 above 10 w_b and a
    straight line in log w between (a LHP pole the same with -90; a RHP zero -90 and a RHP pole +90). A
    complex pair adds +-180 by a straight line in log w from w_n 10^(-zeta) to w_n 10^(zeta) (a step when
    zeta -> 0); the sign flips for a RHP pair.
    """
    bf = bode_form(G)
    lw = np.log10(np.asarray(w, dtype=float))
    ph = np.full(lw.shape, -90.0 * bf["type"] + (180.0 if bf["K"] < 0 else 0.0))
    for wb, sgn, lhp in bf["real"]:
        frac = np.clip((lw - (np.log10(wb) - 1.0)) / 2.0, 0.0, 1.0)
        ph += sgn * (1 if lhp else -1) * 90.0 * frac
    for wn, zeta, sgn in bf["complex"]:
        width = max(abs(zeta), 1e-9)
        frac = np.clip((lw - (np.log10(wn) - width)) / (2 * width), 0.0, 1.0)
        ph += sgn * (1 if zeta > 0 else -1) * 180.0 * frac
    return ph


# ---------------------------------------------------------------- features of the response

def bandwidth(G: TF, drop_db: float = 3.0103, w_max: float = 1e7) -> float:
    """Return the lowest omega at which |G(j w)| first falls ``drop_db`` below |G(0)| (log grid, then bisection)."""
    ref = mag_db(freq_response(G, np.array([1e-12])))[0] - drop_db
    w = np.logspace(-6, np.log10(w_max), 4001)
    m = mag_db(freq_response(G, w))
    i = int(np.argmax(m < ref))
    if m[i] >= ref:
        return float("inf")
    lo, hi = w[i - 1], w[i]
    for _ in range(100):
        mid = np.sqrt(lo * hi)
        if mag_db(freq_response(G, np.array([mid])))[0] >= ref:
            lo = mid
        else:
            hi = mid
    return float(np.sqrt(lo * hi))


def gain_crossover(G: TF, w_min: float = 1e-4, w_max: float = 1e6) -> tuple[float, float]:
    """Return (omega_c, phase at omega_c [deg]) where |G(j w)| first falls through 1 (0 dB), for a loop gain G."""
    w = np.logspace(np.log10(w_min), np.log10(w_max), 8001)
    m = mag_db(freq_response(G, w))
    i = int(np.argmax(m < 0.0))
    if i == 0:
        raise ValueError("no 0 dB crossing from above in the range")
    lo, hi = w[i - 1], w[i]
    for _ in range(100):
        mid = np.sqrt(lo * hi)
        if mag_db(freq_response(G, np.array([mid])))[0] >= 0.0:
            lo = mid
        else:
            hi = mid
    wc = float(np.sqrt(lo * hi))
    return wc, float(phase_from_factors(G, np.array([wc]))[0])


def resonance(zeta: float, wn: float) -> tuple[float, float]:
    """Return (omega_r, M_r) of the standard second order: omega_n sqrt(1 - 2 zeta^2) and 1/(2 zeta sqrt(1 - zeta^2)), or (0, 1) if zeta >= 1/sqrt 2."""
    if zeta >= 1.0 / np.sqrt(2.0):
        return 0.0, 1.0
    return float(wn * np.sqrt(1.0 - 2.0 * zeta * zeta)), float(1.0 / (2.0 * zeta * np.sqrt(1.0 - zeta * zeta)))


# ---------------------------------------------------------------- the simulated sine test

def fit_sinusoid(t: np.ndarray, y: np.ndarray, w: float) -> tuple[float, float]:
    """Return (amplitude, phase [deg]) of the least-squares fit y ~ a cos(w t) + b sin(w t) + c.

    With y = A cos(w t + phi): a = A cos phi, b = -A sin phi, so A = hypot(a, b), phi = atan2(-b, a).
    """
    M = np.column_stack([np.cos(w * t), np.sin(w * t), np.ones_like(t)])
    (a, b, _), *_ = np.linalg.lstsq(M, y, rcond=None)
    return float(np.hypot(a, b)), float(np.degrees(np.arctan2(-b, a)))


def sine_test(G: TF, w: float, amplitude: float = 1.0, settle: float | None = None, n_periods: int = 5,
              pts_per_period: int = 200) -> tuple[float, float]:
    """Return the SIMULATED steady-state (gain, phase [deg]) of G at omega = w.

    Drive G from rest with u = amplitude cos(w t), integrate by RK4 for ``settle`` seconds (default: 10 /
    the slowest decay rate |Re p|, so the transient has decayed by e^-10) plus ``n_periods`` periods, and fit
    the last ``n_periods`` with ``fit_sinusoid``. Gain = fitted amplitude / input amplitude; phase relative
    to the input, wrapped to (-180, 180].
    """
    p = np.roots(tf(*G)[1])
    if settle is None:
        slow = np.min(np.abs(p.real)) if len(p) else 1.0
        settle = 10.0 / slow
    T = 2 * np.pi / w
    h = T / pts_per_period
    fast = np.max(np.abs(p)) if len(p) else 0.0
    if fast > 0:
        h = min(h, 0.5 / fast)                       # keep RK4 accurate on the fastest mode
    n = int(np.ceil((settle + n_periods * T) / h))
    t = np.arange(n + 1) * h
    y = simulate_tf(G, lambda tt: amplitude * np.cos(w * tt), t)
    keep = t >= t[-1] - n_periods * T
    A, ph = fit_sinusoid(t[keep], y[keep], w)
    return A / amplitude, ph


# ---------------------------------------------------------------- the example

def main(out_dir: str | Path = ".") -> None:
    """Print the Bode results for the servo, the gyro filter and the pitch loop; save bode.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    G = servo_closed_loop()
    zeta, wn = 1 / (2 * np.sqrt(14.31 * 0.05)), np.sqrt(14.31 / 0.05)
    wr, Mr = resonance(zeta, wn)
    wb = bandwidth(G)
    print(f"elevator servo closed loop 14.31/(0.05 s^2 + s + 14.31): zeta = {zeta:.4f}, wn = {wn:.3f} rad/s")
    print(f"  resonance: omega_r = {wr:.3f} rad/s = {wr / (2 * np.pi):.2f} Hz, M_r = {Mr:.4f} ({20 * np.log10(Mr):.3f} dB);"
          f" -3 dB bandwidth {wb:.3f} rad/s = {wb / (2 * np.pi):.3f} Hz")
    print("  f [Hz]   |G| exact [dB]  phase exact [deg]   SIMULATED sine test: |G| [dB], phase [deg]")
    for f in (0.5, 1.0, 2.0, 4.0, 8.0):
        w = 2 * np.pi * f
        m, ph = bode(G, np.array([w]))
        g_s, ph_s = sine_test(G, w)
        print(f"  {f:5.1f}   {m[0]:+12.4f}   {ph[0]:+14.3f}        {20 * np.log10(g_s):+9.4f}, {ph_s:+9.3f}")
    print(f"  phase lag at 1 Hz is a time delay of {-bode(G, np.array([2 * np.pi]))[1][0] / 360.0 * 1e3:.1f} ms")

    print("\ngyro RC front end, tau = 1.6 ms (corner 99.5 Hz):")
    for n in (1, 2):
        H = rc_filter(1.6e-3, n)
        for f in (10.0, 99.47, 1200.0):
            m, ph = bode(H, np.array([2 * np.pi * f]))
            ma = asymptotic_mag_db(H, np.array([2 * np.pi * f]))[0]
            print(f"  {n} stage(s), {f:7.2f} Hz: exact {m[0]:+8.3f} dB, {ph[0]:+8.2f} deg; straight line {ma:+8.3f} dB")
    g_s, ph_s = sine_test(rc_filter(1.6e-3, 1), 2 * np.pi * 1200.0)
    print(f"  SIMULATED 1200 Hz sine test on one stage: {20 * np.log10(g_s):+.3f} dB, {ph_s:+.2f} deg")

    L = pitch_open_loop()
    bf = bode_form(L)
    wc, phc = gain_crossover(L)
    print(f"\npitch attitude loop L(s) of P09.03: Bode gain K_B = {bf['K']:.4f}, type {bf['type']}")
    print("  corners: " + ", ".join(f"{'zero' if s > 0 else 'pole'} {w_:.3f}" for w_, s, _ in sorted(bf["real"]))
          + "; " + ", ".join(f"{'zero' if s > 0 else 'pole'} pair wn = {w_:.3f} zeta = {z_:.3f}" for w_, z_, s in bf["complex"]))
    print(f"  gain crossover {wc:.3f} rad/s, phase there {phc:.2f} deg (180 + phase = {180 + phc:.2f} deg, P09.10 calls this the phase margin)")
    ww = np.logspace(-1, 3, 2001)
    me, pe = bode(L, ww)
    ma, pa = asymptotic_mag_db(L, ww), asymptotic_phase_deg(L, ww)
    print(f"  straight line vs exact over 0.1-1000 rad/s: max magnitude error {np.max(np.abs(ma - me)):.2f} dB,"
          f" max phase error {np.max(np.abs(pa - pe)):.1f} deg")

    num_l, den_l = tf(*L)
    Gcl = (num_l, np.polyadd(den_l, num_l))          # closed pitch loop L/(1 + L) (Problem 4)
    H25 = freq_response(Gcl, np.array([2.5]))
    m25, p25 = mag_db(H25)[0], float(np.angle(H25[0], deg=True))
    print(f"\nclosed pitch loop at 2.5 rad/s: exact {m25:+.4f} dB, {p25:+.3f} deg;"
          f" slowest pole {-np.min(np.abs(np.roots(Gcl[1]).real)):.3f} 1/s")
    for st in (0.0, 1.0, 3.0, 10.0):
        g_s, ph_s = sine_test(Gcl, 2.5, settle=st, n_periods=3)
        print(f"  SIMULATED sine test, settle {st:4.1f} s, 3 periods: gain error {20 * np.log10(g_s) - m25:+.2e} dB,"
              f" phase error {ph_s - p25:+.2e} deg")
    w_sweep = np.logspace(-1, 1, 20)
    print(f"  sweep 0.1-10 rad/s, 20 points, 3 s settle + 3 periods each: {np.sum(3.0 + 3 * 2 * np.pi / w_sweep):.0f} s")

    fig, axes = plt.subplots(2, 2, figsize=(11, 6.5), sharex="col")
    w1 = np.logspace(0, 3, 600)
    m1, p1 = bode(G, w1)
    axes[0, 0].semilogx(w1, m1, label="exact")
    axes[0, 0].semilogx(w1, asymptotic_mag_db(G, w1), "--", label="straight line")
    fs = np.array([0.5, 1, 2, 4, 8]) * 2 * np.pi
    tests = [sine_test(G, w_) for w_ in fs]
    axes[0, 0].semilogx(fs, [20 * np.log10(g) for g, _ in tests], "ko", ms=4, label="SIMULATED sine test")
    axes[1, 0].semilogx(w1, p1)
    axes[1, 0].semilogx(w1, asymptotic_phase_deg(G, w1), "--")
    axes[1, 0].semilogx(fs, [p for _, p in tests], "ko", ms=4)
    axes[0, 0].set_title("elevator servo, closed loop", fontsize=9)
    axes[0, 0].legend(fontsize=7)
    m2, p2 = bode(L, ww)
    axes[0, 1].semilogx(ww, m2, label="exact")
    axes[0, 1].semilogx(ww, ma, "--", label="straight line")
    axes[0, 1].axhline(0, color="k", lw=0.5)
    axes[1, 1].semilogx(ww, p2)
    axes[1, 1].semilogx(ww, pa, "--")
    axes[1, 1].axhline(-180, color="k", lw=0.5)
    axes[0, 1].set_title("pitch attitude loop L(j omega)", fontsize=9)
    axes[0, 1].legend(fontsize=7)
    for ax in axes.ravel():
        ax.grid(True, which="both", lw=0.3)
    axes[0, 0].set_ylabel("|G| [dB]"); axes[1, 0].set_ylabel("phase [deg]")
    axes[1, 0].set_xlabel("omega [rad/s]"); axes[1, 1].set_xlabel("omega [rad/s]")
    fig.tight_layout()
    fig.savefig(out_dir / "bode.png", dpi=150)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'bode.png'}")


if __name__ == "__main__":
    import sys

    main(sys.argv[1] if len(sys.argv) > 1 else ".")
