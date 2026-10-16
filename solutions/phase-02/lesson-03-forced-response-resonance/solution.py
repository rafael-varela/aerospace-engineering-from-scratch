"""P02.03 Forced Response & Resonance: reference solution.

Steady-state and transient response of m x'' + c x' + k x = F(t) to sinusoidal, step and near-resonant
forcing, applied to an avionics bracket driven by engine vibration: undetermined coefficients, the
exponential response formula and complex gain, magnification and phase lag, half-power bandwidth and
Q, pure resonance (secular growth) and beats, step response and overshoot, base-excitation
transmissibility. Closed forms are checked against an RK4 stepper written from scratch and against
``scipy.integrate.solve_ivp`` as the oracle.

Conventions:
    * Forcing F(t) = F0 cos(omega t) [N]; steady state x_p = X cos(omega t - phi) with X = F0 |H|,
      H(omega) = 1 / (k - m omega^2 + i c omega) the complex gain [m/N], phi = -arg H in [0, pi].
    * Magnification M = k |H| = X / (F0 / k) is the amplitude relative to the static deflection.
    * omega_n = sqrt(k/m), zeta = c / (2 sqrt(k m)), r = omega / omega_n.

Run:  python solution.py          (prints the bracket example, writes resonance_explorer.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-02/lesson-03-forced-response-resonance
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
from scipy.optimize import brentq


def avionics_bracket() -> dict[str, float]:
    """Return the model of a 2.5 kg avionics box on its bracket (SI units).

    Natural frequency 55 Hz (k = m (2 pi 55)^2 = 298.5 kN/m), damping ratio 0.02 (a bolted aluminium
    bracket), driven by an engine-induced force of amplitude 20 N whose frequency sweeps 50-60 Hz.
    """
    m, fn, zeta = 2.5, 55.0, 0.02
    k = m * (2.0 * np.pi * fn) ** 2
    return {"m": m, "k": k, "c": 2.0 * zeta * np.sqrt(k * m), "F0": 20.0}


# ---------------------------------------------------------------- parameters

def natural_frequency(m: float, k: float) -> float:
    """Return omega_n = sqrt(k/m) [rad/s]."""
    if m <= 0.0 or k <= 0.0:
        raise ValueError("mass and stiffness must be positive")
    return float(np.sqrt(k / m))


def damping_ratio(m: float, c: float, k: float) -> float:
    """Return zeta = c / (2 sqrt(k m))."""
    if c < 0.0:
        raise ValueError("damping must be non-negative")
    return float(c / (2.0 * np.sqrt(k * m)))


# ---------------------------------------------------------------- steady state: complex gain and undetermined coefficients

def complex_gain(omega: np.ndarray | float, m: float, c: float, k: float) -> np.ndarray:
    """Return H(omega) = 1 / (k - m omega^2 + i c omega) [m/N], the complex gain from force to displacement.

    From the exponential response formula: the particular solution of m x'' + c x' + k x = F0 e^{i omega t}
    is F0 H(omega) e^{i omega t}. Raise ``ValueError`` if c = 0 and omega = omega_n (the gain is infinite).
    """
    omega = np.asarray(omega, dtype=float)
    denom = k - m * omega ** 2 + 1j * c * omega
    if np.any(denom == 0.0):
        raise ValueError("undamped system driven exactly at resonance: no bounded steady state")
    return 1.0 / denom


def magnification(omega: np.ndarray | float, m: float, c: float, k: float) -> np.ndarray:
    """Return M = k |H| = 1 / sqrt((1 - r^2)^2 + (2 zeta r)^2): steady amplitude over the static deflection F0/k."""
    return k * np.abs(complex_gain(omega, m, c, k))


def phase_lag(omega: np.ndarray | float, m: float, c: float, k: float) -> np.ndarray:
    """Return phi(omega) in [0, pi] [rad]: the response x lags the force by phi, x_p = |H| F0 cos(omega t - phi)."""
    omega = np.asarray(omega, dtype=float)
    return np.arctan2(c * omega, k - m * omega ** 2)


def undetermined_coefficients(omega: float, F0: float, m: float, c: float, k: float) -> tuple[float, float]:
    """Return (A, B) such that x_p = A cos(omega t) + B sin(omega t) solves m x'' + c x' + k x = F0 cos(omega t).

    Substituting and matching the cos and sin coefficients gives the 2x2 system
        (k - m omega^2) A + c omega B = F0,   -c omega A + (k - m omega^2) B = 0,
    solved here with ``numpy.linalg.solve`` so that the algebra can be checked against ``complex_gain``.
    """
    a = k - m * omega ** 2
    b = c * omega
    A, B = np.linalg.solve(np.array([[a, b], [-b, a]]), np.array([F0, 0.0]))
    return float(A), float(B)


def steady_state_response(t: np.ndarray | float, omega: float, F0: float, m: float, c: float,
                          k: float) -> np.ndarray:
    """Return x_p(t) = Re[F0 H(omega) e^{i omega t}] = F0 |H| cos(omega t - phi): the steady-state response."""
    t = np.asarray(t, dtype=float)
    H = complex_gain(omega, m, c, k)
    return np.real(F0 * H * np.exp(1j * omega * t))


# ---------------------------------------------------------------- the resonance peak

def peak_frequency(m: float, c: float, k: float) -> float:
    """Return the forcing frequency [rad/s] of maximum magnification, omega_n sqrt(1 - 2 zeta^2); 0.0 if zeta >= 1/sqrt(2) (no peak)."""
    z = damping_ratio(m, c, k)
    if z >= 1.0 / np.sqrt(2.0):
        return 0.0
    return float(natural_frequency(m, k) * np.sqrt(1.0 - 2.0 * z * z))


def peak_magnification(m: float, c: float, k: float) -> float:
    """Return M at the peak, 1 / (2 zeta sqrt(1 - zeta^2)) for zeta < 1/sqrt(2); 1.0 otherwise (the maximum is at omega = 0)."""
    z = damping_ratio(m, c, k)
    if z >= 1.0 / np.sqrt(2.0):
        return 1.0
    return float(1.0 / (2.0 * z * np.sqrt(1.0 - z * z)))


def half_power_bandwidth(m: float, c: float, k: float) -> tuple[float, float]:
    """Return (omega_lo, omega_hi) [rad/s] where M^2 falls to half its peak value, found with ``brentq`` on each side of the peak.

    For light damping omega_hi - omega_lo ~ 2 zeta omega_n. Raise ``ValueError`` if there is no peak (zeta >= 1/sqrt(2)).
    """
    wp = peak_frequency(m, c, k)
    if wp == 0.0:
        raise ValueError("no resonance peak for zeta >= 1/sqrt(2)")
    target = 0.5 * peak_magnification(m, c, k) ** 2
    f = lambda w: float(magnification(w, m, c, k)) ** 2 - target
    wn = natural_frequency(m, k)
    lo = brentq(f, 1e-9 * wn, wp)
    hi = brentq(f, wp, 10.0 * wn)
    return float(lo), float(hi)


def quality_factor(m: float, c: float, k: float) -> float:
    """Return Q = omega_peak / (omega_hi - omega_lo), the resonance sharpness; ~ 1 / (2 zeta) for light damping."""
    lo, hi = half_power_bandwidth(m, c, k)
    return float(peak_frequency(m, c, k) / (hi - lo))


# ---------------------------------------------------------------- transients: full response, resonance, beats, step

def free_response(t: np.ndarray, m: float, c: float, k: float, x0: float, v0: float) -> np.ndarray:
    """Return the underdamped homogeneous solution with x(0) = x0, x'(0) = v0 (P02.02, eq. 5). Raise ``ValueError`` if zeta >= 1."""
    t = np.asarray(t, dtype=float)
    wn, z = natural_frequency(m, k), damping_ratio(m, c, k)
    if z >= 1.0:
        raise ValueError("free_response here covers the underdamped case only")
    wd, sig = wn * np.sqrt(1.0 - z * z), z * wn
    return np.exp(-sig * t) * (x0 * np.cos(wd * t) + (v0 + sig * x0) / wd * np.sin(wd * t))


def full_response(t: np.ndarray | float, omega: float, F0: float, m: float, c: float, k: float,
                  x0: float, v0: float) -> np.ndarray:
    """Return x(t) = x_h(t) + x_p(t) for F0 cos(omega t) forcing with x(0) = x0, x'(0) = v0 (underdamped, c > 0).

    The homogeneous part is chosen so that the *sum* meets the initial conditions: x_h(0) = x0 - x_p(0),
    x_h'(0) = v0 - x_p'(0). The transient dies as e^{-zeta omega_n t} and only x_p remains.
    """
    t = np.asarray(t, dtype=float)
    H = complex_gain(omega, m, c, k)
    xp0 = float(np.real(F0 * H))
    vp0 = float(np.real(F0 * H * 1j * omega))
    return free_response(t, m, c, k, x0 - xp0, v0 - vp0) + steady_state_response(t, omega, F0, m, c, k)


def undamped_resonance(t: np.ndarray | float, F0: float, m: float, k: float) -> np.ndarray:
    """Return x(t) = (F0 / (2 m omega_n)) t sin(omega_n t): the response from rest of an undamped system driven at omega_n.

    The exceptional case of undetermined coefficients: the trial A cos + B sin is itself a homogeneous solution,
    so the particular solution carries a factor t and the amplitude grows without bound.
    """
    t = np.asarray(t, dtype=float)
    wn = natural_frequency(m, k)
    return F0 / (2.0 * m * wn) * t * np.sin(wn * t)


def beats(t: np.ndarray | float, omega: float, F0: float, m: float, k: float) -> np.ndarray:
    """Return the undamped response from rest to F0 cos(omega t), omega != omega_n:
    x = F0 / (m (omega_n^2 - omega^2)) (cos(omega t) - cos(omega_n t)) = 2 F0 / (m (omega_n^2 - omega^2)) sin((omega_n - omega) t / 2) sin((omega_n + omega) t / 2),
    a fast oscillation inside a slow envelope of frequency |omega_n - omega| / 2 (beats). Raise ``ValueError`` at omega = omega_n.
    """
    t = np.asarray(t, dtype=float)
    wn = natural_frequency(m, k)
    if omega == wn:
        raise ValueError("omega = omega_n: use undamped_resonance")
    return F0 / (m * (wn ** 2 - omega ** 2)) * (np.cos(omega * t) - np.cos(wn * t))


def step_response(t: np.ndarray | float, F0: float, m: float, c: float, k: float) -> np.ndarray:
    """Return the response from rest to a constant force F0 applied at t = 0 (underdamped):
    x = (F0/k) [1 - e^{-zeta wn t} (cos(wd t) + (zeta wn / wd) sin(wd t))]."""
    t = np.asarray(t, dtype=float)
    wn, z = natural_frequency(m, k), damping_ratio(m, c, k)
    if z >= 1.0:
        raise ValueError("step_response here covers the underdamped case only")
    wd, sig = wn * np.sqrt(1.0 - z * z), z * wn
    return F0 / k * (1.0 - np.exp(-sig * t) * (np.cos(wd * t) + sig / wd * np.sin(wd * t)))


def step_overshoot(zeta: float) -> float:
    """Return the fractional overshoot of the step response, exp(-pi zeta / sqrt(1 - zeta^2)), for 0 <= zeta < 1."""
    if not 0.0 <= zeta < 1.0:
        raise ValueError("overshoot is defined for 0 <= zeta < 1")
    return float(np.exp(-np.pi * zeta / np.sqrt(1.0 - zeta * zeta)))


def transient_decay_time(m: float, c: float, k: float, fraction: float = 0.01) -> float:
    """Return the time [s] for the transient envelope e^{-zeta omega_n t} to fall to ``fraction``: -ln(fraction) / (zeta omega_n)."""
    sig = damping_ratio(m, c, k) * natural_frequency(m, k)
    if sig <= 0.0:
        return float("inf")
    return float(-np.log(fraction) / sig)


# ---------------------------------------------------------------- base excitation

def transmissibility(omega: np.ndarray | float, m: float, c: float, k: float) -> np.ndarray:
    """Return T = |x| / |y| for base motion y = Y cos(omega t): T = sqrt((1 + (2 zeta r)^2) / ((1 - r^2)^2 + (2 zeta r)^2)).

    All curves cross T = 1 at r = sqrt(2) whatever the damping; above it more damping transmits *more*.
    """
    omega = np.asarray(omega, dtype=float)
    return np.abs((k + 1j * c * omega) * complex_gain(omega, m, c, k))


# ---------------------------------------------------------------- numerical oracle

def rk4_forced(m: float, c: float, k: float, force: Callable[[float], float], x0: float, v0: float,
               t_end: float, n_steps: int) -> tuple[np.ndarray, np.ndarray]:
    """Integrate m x'' + c x' + k x = force(t) with classical RK4 on a uniform grid; return (times, states (n+1, 2))."""
    if n_steps < 1:
        raise ValueError("n_steps must be at least 1")

    def f(t: float, y: np.ndarray) -> np.ndarray:
        return np.array([y[1], (force(t) - c * y[1] - k * y[0]) / m])

    ts = np.linspace(0.0, t_end, n_steps + 1)
    h = t_end / n_steps
    ys = np.zeros((n_steps + 1, 2))
    ys[0] = [x0, v0]
    for i in range(n_steps):
        t, y = ts[i], ys[i]
        k1 = f(t, y)
        k2 = f(t + 0.5 * h, y + 0.5 * h * k1)
        k3 = f(t + 0.5 * h, y + 0.5 * h * k2)
        k4 = f(t + h, y + h * k3)
        ys[i + 1] = y + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    return ts, ys


# ---------------------------------------------------------------- the example

def main(out_dir: str | Path = ".") -> None:
    """Print the avionics-bracket resonance numbers and save resonance_explorer.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.integrate import solve_ivp

    out_dir = Path(out_dir)
    p = avionics_bracket()
    m, c, k, F0 = p["m"], p["c"], p["k"], p["F0"]
    wn, z = natural_frequency(m, k), damping_ratio(m, c, k)
    print(f"bracket: k = {k / 1e3:.1f} kN/m, c = {c:.2f} N s/m, f_n = {wn / 2 / np.pi:.2f} Hz, zeta = {z:.3f},"
          f" static deflection F0/k = {F0 / k * 1e6:.1f} um")
    wp = peak_frequency(m, c, k)
    lo, hi = half_power_bandwidth(m, c, k)
    print(f"peak at {wp / 2 / np.pi:.3f} Hz, M_max = {peak_magnification(m, c, k):.2f} (1/(2 zeta) = {1 / (2 * z):.2f}),"
          f" half-power band {lo / 2 / np.pi:.2f}-{hi / 2 / np.pi:.2f} Hz, Q = {quality_factor(m, c, k):.1f}")
    print(f"{'f [Hz]':>7} {'M':>7} {'X [um]':>8} {'phase lag [deg]':>16} {'T (base)':>9}")
    for f in (40.0, 50.0, 53.0, 55.0, 57.0, 60.0, 80.0, 110.0):
        w = 2 * np.pi * f
        print(f"{f:7.1f} {float(magnification(w, m, c, k)):7.2f} {float(magnification(w, m, c, k)) * F0 / k * 1e6:8.1f}"
              f" {np.degrees(float(phase_lag(w, m, c, k))):16.1f} {float(transmissibility(w, m, c, k)):9.2f}")
    print("design options at a 55 Hz drive: stiffen to f_n = 80 Hz ->",
          f"M = {float(magnification(2 * np.pi * 55, m, c, m * (2 * np.pi * 80) ** 2)):.2f};",
          f"add damping to zeta = 0.10 -> M = {float(magnification(wn, m, 2 * 0.1 * np.sqrt(k * m), k)):.2f}")

    # Transient to steady state at resonance, closed form vs RK4 vs solve_ivp.
    T_dec = transient_decay_time(m, c, k)
    t = np.linspace(0.0, 1.0, 20001)
    x_full = full_response(t, wn, F0, m, c, k, 0.0, 0.0)
    ts, ys = rk4_forced(m, c, k, lambda tt: F0 * np.cos(wn * tt), 0.0, 0.0, 1.0, 20000)
    sol = solve_ivp(lambda tt, y: [y[1], (F0 * np.cos(wn * tt) - c * y[1] - k * y[0]) / m], (0.0, 1.0), [0.0, 0.0],
                    rtol=1e-10, atol=1e-13, t_eval=t)
    print(f"\ntransient decays to 1 % in {T_dec:.3f} s ({T_dec * wn / 2 / np.pi:.0f} cycles);"
          f" at resonance the amplitude after 0.05, 0.1, 0.2, 0.5 s is "
          + ", ".join(f"{np.max(np.abs(x_full[(t >= tt - 0.02) & (t <= tt)])) * 1e6:.0f}" for tt in (0.05, 0.1, 0.2, 0.5))
          + f" um (steady {F0 * abs(complex_gain(wn, m, c, k)) * 1e6:.0f} um)")
    print(f"max |RK4 - closed form| = {np.max(np.abs(ys[:, 0] - x_full)):.2e} m; max |solve_ivp - closed form| = {np.max(np.abs(sol.y[0] - x_full)):.2e} m")
    A, B = undetermined_coefficients(wn, F0, m, c, k)
    H = complex_gain(wn, m, c, k)
    print(f"undetermined coefficients at omega_n: A = {A * 1e6:.3f} um, B = {B * 1e6:.3f} um; from H: Re = {F0 * H.real * 1e6:.3f}, -Im = {-F0 * H.imag * 1e6:.3f}")
    print(f"undamped drive at omega_n: amplitude after 0.1 s would be {F0 / (2 * m * wn) * 0.1 * 1e6:.0f} um and growing {F0 / (2 * m * wn) * 1e6:.0f} um per second")
    print(f"step response: overshoot {step_overshoot(z) * 100:.1f} %, final F0/k = {F0 / k * 1e6:.1f} um")

    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    fr = np.linspace(20.0, 120.0, 2000)
    w = 2 * np.pi * fr
    for zz in (0.02, 0.05, 0.1, 0.2, 0.5):
        cc = 2 * zz * np.sqrt(k * m)
        axes[0, 0].semilogy(fr, magnification(w, m, cc, k), label=f"zeta = {zz}")
        axes[0, 1].plot(fr, np.degrees(phase_lag(w, m, cc, k)))
    axes[0, 0].axvspan(50, 60, color="0.9"); axes[0, 0].set_ylabel("magnification M"); axes[0, 0].legend(fontsize=8)
    axes[0, 0].set_xlabel("forcing frequency [Hz]"); axes[0, 0].grid(True, which="both", lw=0.4)
    axes[0, 1].set_ylabel("phase lag [deg]"); axes[0, 1].set_xlabel("forcing frequency [Hz]"); axes[0, 1].grid(lw=0.4)
    axes[1, 0].plot(t[t <= 0.6], x_full[t <= 0.6] * 1e6, lw=0.6)
    axes[1, 0].set_xlabel("t [s]"); axes[1, 0].set_ylabel("x [um], drive at f_n"); axes[1, 0].grid(lw=0.4)
    tb = np.linspace(0.0, 1.0, 20001)
    axes[1, 1].plot(tb, beats(tb, 2 * np.pi * 52.0, F0, m, k) * 1e6, lw=0.5, label="undamped, 52 Hz drive: beats")
    axes[1, 1].plot(tb, undamped_resonance(tb, F0, m, k) * 1e6, lw=0.5, label="undamped, 55 Hz drive: t sin")
    axes[1, 1].set_xlabel("t [s]"); axes[1, 1].set_ylabel("x [um]"); axes[1, 1].legend(fontsize=8); axes[1, 1].grid(lw=0.4)
    fig.suptitle("Avionics bracket, f_n = 55 Hz: resonance explorer", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_dir / "resonance_explorer.png", dpi=150)
    print(f"saved {out_dir / 'resonance_explorer.png'}")


if __name__ == "__main__":
    main()
