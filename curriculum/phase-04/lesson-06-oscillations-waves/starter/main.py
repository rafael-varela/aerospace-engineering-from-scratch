"""P04.06 Oscillations & Waves: starter code.

Complete every function marked NotImplementedError. When finished, running this file produces the
lesson's figure and prints the verification tables.

Simple harmonic motion, damping, resonance, the pendulum beyond the small-angle limit, and the
one-dimensional wave equation on a string with fixed ends:

* closed forms for the free, damped and driven oscillator, checked against solve_ivp;
* the pendulum period: small-angle 2 pi sqrt(L/g), the exact elliptic-integral result, the series
  T0 (1 + theta0^2/16 + 11 theta0^4/3072), and a simulated period from solve_ivp with an event;
* the resonance curve (steady-state amplitude and phase), quality factor and half-power bandwidth;
* a finite-difference wave propagation visualizer (second order in x and t, fixed ends) checked
  against d'Alembert's solution, the standing-wave modes and energy conservation, with the CFL limit;
* the reduction of pendulum timing data: period from repeated timings, g with its uncertainty, and
  a power-law fit of T against L. ``make_synthetic_pendulum_data`` is SYNTHETIC: it exists so the
  tests can exercise the reduction on data whose truth is known. It is not a measurement.

Units: SI (m, s, kg, rad). Angles in radians unless a name says deg.

Run:  python main.py          (prints the verification tables, writes oscillations_waves.png)
Test: pytest curriculum/phase-04/lesson-06-oscillations-waves
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipk

G_STANDARD = 9.80665                 # m/s^2


# ---------------------------------------------------------------- simple harmonic motion and damping

def shm(t: np.ndarray | float, x0: float, v0: float, omega: float) -> np.ndarray | float:
    """x(t) = x0 cos(omega t) + (v0 / omega) sin(omega t): the solution of x'' + omega^2 x = 0."""
    raise NotImplementedError


def damped_oscillator_rhs(t: float, y: np.ndarray, omega0: float, zeta: float) -> np.ndarray:
    """State y = [x, v] of x'' + 2 zeta omega0 x' + omega0^2 x = 0, for solve_ivp."""
    raise NotImplementedError


def damped_response(t: np.ndarray | float, x0: float, v0: float, omega0: float, zeta: float) -> np.ndarray | float:
    """Closed-form free response of x'' + 2 zeta omega0 x' + omega0^2 x = 0 for any zeta >= 0.

    zeta < 1: e^{-zeta omega0 t} (x0 cos(wd t) + (v0 + zeta omega0 x0)/wd sin(wd t)), wd = omega0 sqrt(1 - zeta^2)
    zeta = 1: e^{-omega0 t} (x0 + (v0 + omega0 x0) t)
    zeta > 1: e^{-zeta omega0 t} (x0 cosh(s t) + (v0 + zeta omega0 x0)/s sinh(s t)), s = omega0 sqrt(zeta^2 - 1)
    """
    raise NotImplementedError


def log_decrement_to_zeta(delta: float) -> float:
    """Damping ratio from the logarithmic decrement delta = ln(x_n / x_{n+1}) of successive peaks."""
    raise NotImplementedError


def shm_energy(x: np.ndarray, v: np.ndarray, m: float, k: float) -> np.ndarray:
    """Total mechanical energy m v^2 / 2 + k x^2 / 2 along a trajectory."""
    raise NotImplementedError


# ---------------------------------------------------------------- the pendulum

def pendulum_rhs(t: float, y: np.ndarray, L: float, g: float = G_STANDARD, b: float = 0.0) -> np.ndarray:
    """State y = [theta, theta_dot] of theta'' + b theta' + (g / L) sin(theta) = 0 (b: damping per unit time)."""
    raise NotImplementedError


def pendulum_period_small(L: float, g: float = G_STANDARD) -> float:
    """Small-angle period T0 = 2 pi sqrt(L / g)."""
    raise NotImplementedError


def pendulum_period_exact(L: float, theta0: float, g: float = G_STANDARD) -> float:
    """Exact period of the undamped pendulum released from rest at amplitude theta0 [rad]:
    T = 4 sqrt(L / g) K(k^2), k = sin(theta0 / 2), K the complete elliptic integral of the first kind."""
    raise NotImplementedError


def pendulum_period_series(L: float, theta0: float, g: float = G_STANDARD) -> float:
    """Series T0 (1 + theta0^2 / 16 + 11 theta0^4 / 3072) for the amplitude dependence of the period."""
    raise NotImplementedError


def pendulum_period_simulated(L: float, theta0: float, g: float = G_STANDARD, rtol: float = 1e-10) -> float:
    """Period from solve_ivp: the time between the first and third zero crossings of theta
    (released from rest at theta0, so the first crossing is at T/4 and the third at 3T/4)."""
    raise NotImplementedError


def amplitude_correction(theta0: float) -> float:
    """T / T0 for amplitude theta0 [rad], from the exact result (independent of L and g)."""
    raise NotImplementedError


# ---------------------------------------------------------------- forced response and resonance

def steady_state_amplitude(omega: np.ndarray | float, omega0: float, zeta: float, f0: float) -> np.ndarray | float:
    """Steady-state amplitude of x'' + 2 zeta omega0 x' + omega0^2 x = f0 cos(omega t):
    X = f0 / sqrt((omega0^2 - omega^2)^2 + (2 zeta omega0 omega)^2)."""
    raise NotImplementedError


def steady_state_phase(omega: np.ndarray | float, omega0: float, zeta: float) -> np.ndarray | float:
    """Phase lag of the response behind the drive, in [0, pi]: atan2(2 zeta omega0 omega, omega0^2 - omega^2)."""
    raise NotImplementedError


def quality_factor(zeta: float) -> float:
    """Q = 1 / (2 zeta)."""
    raise NotImplementedError


def half_power_bandwidth(omega0: float, zeta: float) -> float:
    """Width 2 zeta omega0 (= omega0 / Q) between the two frequencies where the amplitude is 1/sqrt(2) of
    its peak, exact as zeta -> 0 and within a few percent for zeta < 0.1."""
    raise NotImplementedError


def resonance_peak_frequency(omega0: float, zeta: float) -> float:
    """Frequency omega0 sqrt(1 - 2 zeta^2) at which the steady-state amplitude peaks (0 if zeta >= 1/sqrt(2))."""
    raise NotImplementedError


# ---------------------------------------------------------------- waves on a string

def string_wave_speed(tension: float, mu: float) -> float:
    """c = sqrt(T / mu) for a string under tension T [N] with mass per length mu [kg/m]."""
    raise NotImplementedError


def standing_wave_frequencies(L: float, c: float, n_max: int) -> np.ndarray:
    """f_n = n c / (2 L), n = 1..n_max, for a string fixed at both ends."""
    raise NotImplementedError


def standing_wave(x: np.ndarray, t: float, n: int, L: float, c: float, A: float = 1.0) -> np.ndarray:
    """Mode n of the fixed-fixed string released from rest: A sin(n pi x / L) cos(n pi c t / L)."""
    raise NotImplementedError


def odd_periodic_extension(f: Callable[[np.ndarray], np.ndarray], L: float) -> Callable[[np.ndarray], np.ndarray]:
    """Extend a profile f on [0, L] as an odd function about 0 and L (period 2L), as d'Alembert's solution
    for fixed ends requires."""
    raise NotImplementedError


def dalembert(x: np.ndarray, t: float, c: float, f: Callable[[np.ndarray], np.ndarray], L: float) -> np.ndarray:
    """d'Alembert's solution u(x, t) = [F(x - c t) + F(x + c t)] / 2 for a string fixed at 0 and L,
    released from rest with profile f (F is the odd periodic extension of f)."""
    raise NotImplementedError


def wave_equation_fd(c: float, L: float, nx: int, t_end: float, courant: float,
                     u0: Callable[[np.ndarray], np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Solve u_tt = c^2 u_xx on [0, L] with u(0) = u(L) = 0, released from rest with u(x, 0) = u0(x),
    by central differences (second order in x and t). ``courant`` = c dt / dx sets the step; the
    scheme is stable for courant <= 1 and reproduces d'Alembert exactly on the grid at courant = 1.

    Returns (x, t, U) with U of shape (n_t + 1, nx + 1).
    """
    raise NotImplementedError


def wave_energy(U: np.ndarray, dt: float, dx: float, c: float) -> np.ndarray:
    """Discrete energy per unit mu at the half-steps: sum over x of [ (u_t)^2 + c^2 (u_x)^2 ] dx / 2,
    with u_t from the centred difference between levels n and n+1 and u_x averaged over the two
    levels. Returns an array of length n_t (one per interval); it is conserved by the scheme."""
    raise NotImplementedError


def is_cfl_stable(courant: float) -> bool:
    """The explicit scheme is stable iff c dt / dx <= 1."""
    raise NotImplementedError


# ---------------------------------------------------------------- reducing pendulum data

def period_from_timings(times_for_n_swings: np.ndarray, n_swings: int, u_single: float = 0.0) -> tuple[float, float]:
    """Mean period and its standard uncertainty from repeated timings of n_swings full swings.

    Type A: sample standard deviation of the per-swing periods over sqrt(count). Type B: ``u_single``
    (for example the frame-quantisation term 0.41 / fps of a video timing) divided by n_swings, added in
    quadrature after the averaging so it does not shrink with the number of repeats."""
    raise NotImplementedError


def g_from_pendulum(L: float, T: float, u_L: float, u_T: float) -> tuple[float, float]:
    """g = 4 pi^2 L / T^2 and its uncertainty by the power rule: u_g / g = sqrt((u_L / L)^2 + (2 u_T / T)^2)."""
    raise NotImplementedError


def fit_power_law(L: np.ndarray, T: np.ndarray) -> tuple[float, float]:
    """Least-squares fit of T = A L^p in log-log form. Returns (A, p); the model predicts p = 1/2 and
    A = 2 pi / sqrt(g)."""
    raise NotImplementedError


def make_synthetic_pendulum_data(seed: int, lengths: np.ndarray | None = None, g: float = 9.81,
                                 theta0: float = 0.15, n_swings: int = 10, repeats: int = 5,
                                 fps: float = 30.0) -> tuple[np.ndarray, np.ndarray]:
    """SYNTHETIC timings, for the tests only: for each length, ``repeats`` timings of ``n_swings`` swings,
    equal to the exact period at amplitude theta0 times n_swings, quantised to whole frames at fps with a
    random phase, plus a small release scatter. Not a measurement; nothing from it belongs in a report.

    Returns (lengths, times) with times of shape (len(lengths), repeats)."""
    raise NotImplementedError


# ---------------------------------------------------------------- driver

def main(out_dir: str | Path = ".") -> None:
    """Print the verification tables and write oscillations_waves.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    g = G_STANDARD

    # damped oscillator against solve_ivp
    w0, z = 2.0 * np.pi * 1.5, 0.05
    t = np.linspace(0.0, 4.0, 401)
    sol = solve_ivp(damped_oscillator_rhs, (0.0, 4.0), [0.02, 0.0], args=(w0, z), t_eval=t, method="DOP853",
                    rtol=1e-11, atol=1e-14)
    x = damped_response(t, 0.02, 0.0, w0, z)
    print(f"damped oscillator f0 = 1.5 Hz, zeta = {z}: closed form vs DOP853 max |dx| = {np.max(np.abs(x - sol.y[0])):.2e} m; "
          f"Q = {quality_factor(z):.1f}, half-power bandwidth {half_power_bandwidth(w0, z)/(2*np.pi):.3f} Hz, "
          f"amplitude at resonance / static = {steady_state_amplitude(w0, w0, z, 1.0) / steady_state_amplitude(0.0, w0, z, 1.0):.1f}")
    peaks = x[np.r_[False, (x[1:-1] > x[:-2]) & (x[1:-1] > x[2:]), False]]
    delta = float(np.mean(np.log(peaks[:-1] / peaks[1:])))
    print(f"  log decrement from the peaks {delta:.4f} -> zeta = {log_decrement_to_zeta(delta):.4f}")

    # the pendulum
    L = 1.000
    T0 = pendulum_period_small(L, g)
    print(f"\npendulum L = {L:.3f} m: T0 = {T0:.4f} s")
    print(f"{'theta0 (deg)':>12} {'T exact (s)':>12} {'T series (s)':>13} {'T simulated (s)':>16} {'T/T0 - 1':>10}")
    for deg in (5, 10, 20, 30, 45, 60, 90):
        th = np.radians(deg)
        Te, Ts, Tsim = pendulum_period_exact(L, th, g), pendulum_period_series(L, th, g), pendulum_period_simulated(L, th, g)
        print(f"{deg:>12d} {Te:>12.5f} {Ts:>13.5f} {Tsim:>16.5f} {Te/T0 - 1:>10.4%}")
    u_frame = np.sqrt(2.0) * 0.5 / (np.sqrt(3.0) * 30.0)
    print(f"  video at 30 fps: frame term {u_frame*1e3:.1f} ms on one timing, {u_frame/10*1e3:.2f} ms per period over 10 swings "
          f"({u_frame/10/T0:.2%} of T0); the 30 deg correction is {amplitude_correction(np.radians(30)) - 1:.2%}")

    # synthetic data through the reduction (SYNTHETIC: labelled as such in the lesson)
    lengths, times = make_synthetic_pendulum_data(seed=4061)
    print("\nSYNTHETIC pendulum dataset (seed 4061, g = 9.81 m/s^2, theta0 = 0.15 rad, 10 swings x 5, 30 fps) through the reduction:")
    print(f"{'L (m)':>7} {'T (s)':>8} {'u_T (ms)':>9} {'g (m/s^2)':>10} {'u_g':>6}")
    Ts, uTs = [], []
    for Li, row in zip(lengths, times):
        T, uT = period_from_timings(row, 10, u_single=u_frame)
        gi, ugi = g_from_pendulum(Li, T, 0.002, uT)
        Ts.append(T)
        uTs.append(uT)
        print(f"{Li:>7.3f} {T:>8.4f} {uT*1e3:>9.2f} {gi:>10.3f} {ugi:>6.3f}")
    A_fit, p_fit = fit_power_law(lengths, np.array(Ts))
    print(f"  power-law fit T = A L^p: p = {p_fit:.4f} (model: 0.5), g from A: {4*np.pi**2/A_fit**2:.3f} m/s^2 "
          f"(uncorrected for amplitude; x {amplitude_correction(0.15)**2:.4f} after correction: {4*np.pi**2/A_fit**2*amplitude_correction(0.15)**2:.3f})")

    # the string
    tension, mu_s, Ls = 50.0, 2.0e-3, 0.65
    c = string_wave_speed(tension, mu_s)
    fn = standing_wave_frequencies(Ls, c, 4)
    print(f"\nstring: T = {tension} N, mu = {mu_s*1e3:.1f} g/m, L = {Ls} m: c = {c:.1f} m/s, modes {np.round(fn, 1)} Hz")
    nx = 200
    pulse = lambda xx: np.exp(-((xx - 0.5 * Ls) / 0.03) ** 2)
    x, tt, U = wave_equation_fd(c, Ls, nx, 2.0 * Ls / c, 1.0, pulse)            # one round trip
    err = np.max(np.abs(U - np.array([dalembert(x, ti, c, pulse, Ls) for ti in tt])))
    Em = wave_energy(U, tt[1] - tt[0], x[1] - x[0], c)
    print(f"  FD at Courant 1, {nx} cells, one round trip: max |u - d'Alembert| = {err:.2e}; "
          f"energy spread {(Em.max() - Em.min())/Em.mean():.2e}")
    x, tt, U = wave_equation_fd(c, Ls, nx, 2.0 * Ls / c, 0.8, pulse)
    err8 = np.max(np.abs(U - np.array([dalembert(x, ti, c, pulse, Ls) for ti in tt])))
    Em8 = wave_energy(U, tt[1] - tt[0], x[1] - x[0], c)
    print(f"  FD at Courant 0.8: max |u - d'Alembert| = {err8:.2e} (dispersion); energy spread {(Em8.max() - Em8.min())/Em8.mean():.2e}")
    x, tt, U3 = wave_equation_fd(c, Ls, nx, 3.0 / fn[2], 0.9, lambda xx: standing_wave(xx, 0.0, 3, Ls, c))
    err_mode = np.max(np.abs(U3 - np.array([standing_wave(x, ti, 3, Ls, c) for ti in tt])))
    print(f"  mode 3 at Courant 0.9 over three periods: max error {err_mode:.2e} (phase error of the scheme)")
    x, tt, Ubad = wave_equation_fd(c, Ls, nx, 0.5 * Ls / c, 1.05, pulse)
    print(f"  Courant 1.05: max |u| after half a transit = {np.max(np.abs(Ubad)):.2e} (unstable)")

    # figure
    fig, axs = plt.subplots(2, 2, figsize=(11, 7.5))
    ax = axs[0, 0]
    ws = np.linspace(0.0, 3.0, 601)
    for zz in (0.05, 0.1, 0.25, 0.5):
        ax.plot(ws, steady_state_amplitude(ws, 1.0, zz, 1.0), label=f"zeta = {zz}")
    ax.set_xlabel("omega / omega0")
    ax.set_ylabel("X / X_static")
    ax.set_title("resonance curves")
    ax.legend(fontsize=8)
    ax = axs[0, 1]
    degs = np.linspace(0.0, 120.0, 121)
    ax.plot(degs, [amplitude_correction(np.radians(d)) for d in degs], label="exact (elliptic K)")
    ax.plot(degs, [pendulum_period_series(1.0, np.radians(d), 1.0) / pendulum_period_small(1.0, 1.0) for d in degs],
            "--", label="series to theta0^4")
    ax.set_xlabel("amplitude (deg)")
    ax.set_ylabel("T / T0")
    ax.set_title("pendulum period against amplitude")
    ax.legend(fontsize=8)
    ax = axs[1, 0]
    x, tt, U = wave_equation_fd(c, Ls, nx, 2.0 * Ls / c, 1.0, pulse)
    for frac, lab in ((0.0, "t = 0"), (0.2, "0.2 round trips"), (0.5, "0.5 (reflected, inverted)"), (1.0, "1.0 (back)")):
        i = int(round(frac * (len(tt) - 1)))
        ax.plot(x, U[i] + 0.0, label=lab)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("u")
    ax.set_title("pulse on a fixed-fixed string (FD, Courant 1)")
    ax.legend(fontsize=8)
    ax = axs[1, 1]
    ax.errorbar(lengths, Ts, yerr=2 * np.array(uTs), fmt="o", label="SYNTHETIC timings (k = 2)")
    LL = np.linspace(0.3, 1.6, 100)
    ax.plot(LL, 2 * np.pi * np.sqrt(LL / 9.81), "k--", label="2 pi sqrt(L / g), g = 9.81")
    ax.set_xlabel("L (m)")
    ax.set_ylabel("T (s)")
    ax.set_title("period against length (synthetic data)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "oscillations_waves.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
