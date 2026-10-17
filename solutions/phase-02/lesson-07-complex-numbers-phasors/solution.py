"""P02.07 Complex Numbers & Phasors: reference solution.

A phasor calculator for same-frequency vibration: complex arithmetic and Euler's formula from the
power series, roots of unity, sinusoids as the real part of a rotating complex number, the sum of
two propeller excitations as a function of their phase difference, the complex gain of a
mass-spring-damper, and the transmissibility of an instrument isolation mount, each checked
against a time-domain integration with scipy.integrate.solve_ivp and a least-squares sinusoid fit.

Conventions:
    * A sinusoid is written x(t) = A cos(omega t + phi) = Re(Z e^{i omega t}) with phasor Z = A e^{i phi}.
      Amplitudes A >= 0, phases in radians, returned in (-pi, pi].
    * Complex gain of m x'' + c x' + k x = F0 cos(omega t):  X / F0 = 1 / (k - m omega^2 + i c omega).
      In dimensionless form with r = omega / omega_n:  H(r) = 1 / (1 - r^2 + 2 i zeta r), x = (F0/k) Re(H e^{i omega t}).
    * Transmissibility of a base-excited mount: T(r) = |1 + 2 i zeta r| / |1 - r^2 + 2 i zeta r|.

Run:  python solution.py          (prints the cabin-vibration example, writes phasor_calculator.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-02/lesson-07-complex-numbers-phasors
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


def cabin_vibration_case() -> dict[str, float]:
    """Return the nominal case: a twin-propeller cabin floor and an instrument tray on an isolation mount.

    Both propellers turn at 2400 rpm (40 Hz once-per-revolution). Their floor accelerations at the tray
    position are a1 = 1.2 m/s^2 and a2 = 0.9 m/s^2 at the same frequency; the phase between them, set by the
    synchrophaser, is the variable. The tray (m = 2.5 kg) sits on a mount of natural frequency 30 Hz with
    damping ratio 0.08.
    """
    return {"f_drive": 40.0, "a1": 1.2, "a2": 0.9, "m_tray": 2.5, "f_mount": 30.0, "zeta": 0.08}


# ---------------------------------------------------------------- complex arithmetic and Euler's formula

def to_polar(z: complex) -> tuple[float, float]:
    """Return (modulus, argument) of z with the argument in (-pi, pi] (numpy.angle convention)."""
    z = complex(z)
    return abs(z), float(np.angle(z))


def from_polar(r: float, theta: float) -> complex:
    """Return r e^{i theta} = r (cos theta + i sin theta). Raise ``ValueError`` for a negative modulus."""
    if r < 0.0:
        raise ValueError("modulus must be non-negative; put the sign into the angle")
    return complex(r * np.cos(theta), r * np.sin(theta))


def exp_series(z: complex, n_terms: int = 40) -> complex:
    """Return e^z from its power series sum_{k<n_terms} z^k / k!, computed term by term in complex arithmetic.

    Written from scratch so that Euler's formula can be *observed*: exp_series(1j * theta) returns
    cos(theta) + i sin(theta) to rounding for |theta| of a few radians. Converges for every z; the number
    of terms needed grows with |z| (40 terms suffice to 1e-15 for |z| <= 5).
    """
    z = complex(z)
    term = complex(1.0, 0.0)
    total = complex(1.0, 0.0)
    for k in range(1, n_terms):
        term *= z / k
        total += term
    return total


def euler(theta: float) -> complex:
    """Return cos(theta) + i sin(theta), the unit phasor at angle theta."""
    return complex(np.cos(theta), np.sin(theta))


def de_moivre(theta: float, n: int) -> complex:
    """Return (cos theta + i sin theta)^n by repeated complex multiplication; de Moivre says it equals euler(n theta)."""
    if int(n) != n or n < 0:
        raise ValueError("n must be a non-negative integer")
    w = euler(theta)
    out = complex(1.0, 0.0)
    for _ in range(int(n)):
        out *= w
    return out


def nth_roots(z: complex, n: int) -> np.ndarray:
    """Return the n complex n-th roots of z, equally spaced by 2 pi / n on the circle of radius |z|^{1/n}.

    Raise ``ValueError`` for n < 1.
    """
    if n < 1:
        raise ValueError("n must be a positive integer")
    r, th = to_polar(z)
    k = np.arange(n)
    return r ** (1.0 / n) * np.exp(1j * (th + 2.0 * np.pi * k) / n)


def rotate_2d(v: np.ndarray, theta: float) -> np.ndarray:
    """Rotate the plane vector v = (x, y) anticlockwise by theta using complex multiplication by e^{i theta}.

    The same operation as the P01.15 rotation matrix [[cos, -sin], [sin, cos]] applied to v.
    """
    z = complex(v[0], v[1]) * euler(theta)
    return np.array([z.real, z.imag])


# ---------------------------------------------------------------- phasors and sinusoids

def phasor(amplitude: float, phase: float) -> complex:
    """Return the phasor Z = A e^{i phi} of the sinusoid A cos(omega t + phi). Raise ``ValueError`` for A < 0."""
    if amplitude < 0.0:
        raise ValueError("amplitude must be non-negative; a negative amplitude is a phase shift of pi")
    return from_polar(amplitude, phase)


def sinusoid(t: np.ndarray, amplitude: float, omega: float, phase: float) -> np.ndarray:
    """Return A cos(omega t + phi) sampled at t."""
    return amplitude * np.cos(omega * np.asarray(t, dtype=float) + phase)


def from_phasor(Z: complex, omega: float, t: np.ndarray) -> np.ndarray:
    """Return Re(Z e^{i omega t}): the time signal of a phasor at angular frequency omega."""
    return np.real(Z * np.exp(1j * omega * np.asarray(t, dtype=float)))


def phasor_derivative(Z: complex, omega: float) -> complex:
    """Return i omega Z, the phasor of d/dt Re(Z e^{i omega t}). Differentiation is a multiplication."""
    return 1j * omega * Z


def add_sinusoids(components: list[tuple[float, float]]) -> tuple[float, float]:
    """Return (A, phi) of the sum of same-frequency sinusoids given as a list of (amplitude, phase) pairs.

    The phasors add as complex numbers; the result is converted back with ``to_polar``.
    """
    Z = sum(phasor(a, p) for a, p in components)
    return to_polar(Z)


def combined_amplitude(a1: float, a2: float, delta_phase: float) -> float:
    """Return the amplitude of a1 cos(wt) + a2 cos(wt + delta): sqrt(a1^2 + a2^2 + 2 a1 a2 cos delta)."""
    return float(np.sqrt(a1 * a1 + a2 * a2 + 2.0 * a1 * a2 * np.cos(delta_phase)))


def amplitude_range(a1: float, a2: float) -> tuple[float, float]:
    """Return (|a1 - a2|, a1 + a2): the least and greatest amplitude two same-frequency sinusoids can sum to."""
    return abs(a1 - a2), a1 + a2


def fit_sinusoid(t: np.ndarray, x: np.ndarray, omega: float) -> tuple[float, float]:
    """Fit x(t) ~ a cos(omega t) + b sin(omega t) by least squares (normal equations) and return (A, phi).

    A cos(omega t + phi) = a cos(omega t) + b sin(omega t) with a = A cos(phi), b = -A sin(phi), so the
    phasor is Z = a - i b. This is how an amplitude and a phase are read off a simulated or measured record.
    """
    t = np.asarray(t, dtype=float)
    x = np.asarray(x, dtype=float)
    M = np.column_stack([np.cos(omega * t), np.sin(omega * t)])
    a, b = np.linalg.solve(M.T @ M, M.T @ x)
    return to_polar(complex(a, -b))


# ---------------------------------------------------------------- complex gain of a second-order system

def complex_gain(r: float | np.ndarray, zeta: float) -> complex | np.ndarray:
    """Return H(r) = 1 / (1 - r^2 + 2 i zeta r), the displacement per static deflection of m x'' + c x' + k x = F0 cos(omega t).

    r = omega / omega_n. x(t) = (F0 / k) Re(H e^{i omega t}); |H| is the dynamic magnification and -arg H
    the phase lag of the displacement behind the force.
    """
    r = np.asarray(r, dtype=float)
    H = 1.0 / (1.0 - r * r + 2j * zeta * r)
    return complex(H) if H.ndim == 0 else H


def gain_and_phase_lag(r: float, zeta: float) -> tuple[float, float]:
    """Return (|H|, phase lag in [0, pi]) of the complex gain at frequency ratio r."""
    H = complex_gain(r, zeta)
    lag = -np.angle(H)
    if lag < 0.0:
        lag += 2.0 * np.pi
    return abs(H), float(lag)


def steady_state_displacement(t: np.ndarray, F0: float, omega: float, m: float, c: float, k: float) -> np.ndarray:
    """Return the steady-state x(t) of m x'' + c x' + k x = F0 cos(omega t) as Re(X e^{i omega t}), X = F0 / (k - m omega^2 + i c omega)."""
    X = F0 / (k - m * omega * omega + 1j * c * omega)
    return from_phasor(X, omega, t)


def forced_oscillator_rhs(t: float, y: np.ndarray, F0: float, omega: float, m: float, c: float, k: float) -> np.ndarray:
    """Return dy/dt for y = (x, x') under the force F0 cos(omega t): the time-domain check of the phasor solution."""
    x, v = y
    return np.array([v, (F0 * np.cos(omega * t) - c * v - k * x) / m])


def transmissibility(r: float | np.ndarray, zeta: float) -> float | np.ndarray:
    """Return |X / Y| = |1 + 2 i zeta r| / |1 - r^2 + 2 i zeta r| for a mass on a spring-damper whose base moves as Y cos(omega t).

    Above r = sqrt(2) the mount isolates (T < 1); below it, it amplifies, however much damping it has.
    """
    r = np.asarray(r, dtype=float)
    T = np.abs(1.0 + 2j * zeta * r) / np.abs(1.0 - r * r + 2j * zeta * r)
    return float(T) if T.ndim == 0 else T


def mount_frequency_for_transmissibility(T_target: float, f_drive: float, zeta: float) -> float:
    """Return the mount natural frequency [Hz] at which the transmissibility at ``f_drive`` equals ``T_target`` < 1.

    Bisection on r in [sqrt(2), 50] where T(r) is monotone decreasing; f_mount = f_drive / r.
    Raise ``ValueError`` unless 0 < T_target < 1.
    """
    if not 0.0 < T_target < 1.0:
        raise ValueError("isolation needs 0 < T_target < 1")
    lo, hi = np.sqrt(2.0), 50.0
    if transmissibility(hi, zeta) > T_target:
        raise ValueError("target not reachable below r = 50")
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if transmissibility(mid, zeta) > T_target:
            lo = mid
        else:
            hi = mid
    return f_drive / (0.5 * (lo + hi))


def main(out_dir: str | Path = ".") -> None:
    """Print the cabin-vibration phasor analysis and save the phasor diagram and transmissibility figure."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.integrate import solve_ivp

    out_dir = Path(out_dir)
    case = cabin_vibration_case()
    omega = 2.0 * np.pi * case["f_drive"]
    a1, a2 = case["a1"], case["a2"]
    print(f"Euler's formula from the series: exp_series(i pi) = {exp_series(1j * np.pi):.3e} (should be -1)")
    print(f"cube roots of 8: {np.round(nth_roots(8.0, 3), 6)}")

    print(f"\ntwo propellers at {case['f_drive']} Hz, floor accelerations {a1} and {a2} m/s^2:")
    print("  phase difference [deg]   combined amplitude [m/s^2]   phase of the sum [deg]")
    for deg in (0, 45, 90, 135, 180):
        A, phi = add_sinusoids([(a1, 0.0), (a2, np.radians(deg))])
        print(f"  {deg:8d}                 {A:8.4f}                      {np.degrees(phi):7.2f}")
    lo, hi = amplitude_range(a1, a2)
    print(f"  range: {lo:.2f} (anti-phase) to {hi:.2f} (in phase) m/s^2: a synchrophaser is worth a factor {hi / lo:.1f}")

    m, fn, zeta = case["m_tray"], case["f_mount"], case["zeta"]
    wn = 2.0 * np.pi * fn
    k = m * wn * wn
    c = 2.0 * zeta * m * wn
    r = omega / wn
    H = complex_gain(r, zeta)
    mag, lag = gain_and_phase_lag(r, zeta)
    print(f"\ntray on its mount: m = {m} kg, f_n = {fn} Hz, k = {k:.0f} N/m, c = {c:.2f} N s/m, r = {r:.4f}")
    print(f"  complex gain H = {H:.4f}: magnification {mag:.4f}, phase lag {np.degrees(lag):.2f} deg")
    F0 = 5.0
    t = np.linspace(0.0, 1.0, 20001)
    x_ss = steady_state_displacement(t, F0, omega, m, c, k)
    sol = solve_ivp(forced_oscillator_rhs, (0.0, 1.0), [0.0, 0.0], t_eval=t, args=(F0, omega, m, c, k), rtol=1e-10, atol=1e-13)
    tail = t > 0.6
    A_fit, phi_fit = fit_sinusoid(t[tail], sol.y[0][tail], omega)
    A_ph, phi_ph = to_polar(F0 / (k - m * omega ** 2 + 1j * c * omega))
    print(f"  force {F0} N at {case['f_drive']} Hz: phasor amplitude {1e6 * A_ph:.3f} um, phase {np.degrees(phi_ph):.2f} deg;"
          f" solve_ivp after 0.6 s: {1e6 * A_fit:.3f} um, {np.degrees(phi_fit):.2f} deg;"
          f" max |x_ode - x_phasor| over the last 0.4 s = {1e6 * np.max(np.abs(sol.y[0][tail] - x_ss[tail])):.4f} um")
    T = transmissibility(r, zeta)
    y_base = hi / omega ** 2
    print(f"  transmissibility at r = {r:.3f}: {T:.4f}; floor {hi:.2f} m/s^2 = {1e6 * y_base:.2f} um -> tray {1e6 * T * y_base:.2f} um, {T * hi:.3f} m/s^2")
    for Tt in (0.5, 0.25, 0.1):
        print(f"  for T = {Tt}: mount natural frequency {mount_frequency_for_transmissibility(Tt, case['f_drive'], zeta):.2f} Hz")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2))
    for deg, col in ((0, "tab:green"), (90, "tab:blue"), (150, "tab:red")):
        Z1, Z2 = phasor(a1, 0.0), phasor(a2, np.radians(deg))
        ax1.annotate("", xy=(Z1.real, Z1.imag), xytext=(0, 0), arrowprops=dict(arrowstyle="->", color="0.4"))
        ax1.annotate("", xy=((Z1 + Z2).real, (Z1 + Z2).imag), xytext=(Z1.real, Z1.imag), arrowprops=dict(arrowstyle="->", color=col))
        ax1.annotate("", xy=((Z1 + Z2).real, (Z1 + Z2).imag), xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=col, lw=2))
        ax1.text((Z1 + Z2).real, (Z1 + Z2).imag, f" {deg} deg: {abs(Z1 + Z2):.2f}", color=col, fontsize=8)
    ax1.set_aspect("equal")
    ax1.set_xlim(-0.5, 2.4)
    ax1.set_ylim(-0.3, 1.3)
    ax1.grid(True, lw=0.4)
    ax1.set_xlabel("Re [m/s^2]")
    ax1.set_ylabel("Im [m/s^2]")
    ax1.set_title("phasor sum of the two propeller excitations")
    rr = np.linspace(0.05, 4.0, 800)
    for z in (0.02, 0.08, 0.3):
        ax2.semilogy(rr, transmissibility(rr, z), label=f"zeta = {z}")
    ax2.axvline(np.sqrt(2.0), color="0.5", lw=0.8, ls="--")
    ax2.axhline(1.0, color="0.5", lw=0.8)
    ax2.plot([r], [T], "ko")
    ax2.set_xlabel("r = f_drive / f_mount")
    ax2.set_ylabel("transmissibility |X/Y|")
    ax2.set_title("isolation mount: amplifies below r = sqrt 2")
    ax2.grid(True, which="both", lw=0.4)
    ax2.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "phasor_calculator.png", dpi=150)
    print(f"saved {out_dir / 'phasor_calculator.png'}")


if __name__ == "__main__":
    main()
