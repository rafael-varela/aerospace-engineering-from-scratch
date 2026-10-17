"""P02.02 Second-Order Linear ODEs: reference solution.

Free response of the damped spring-mass system m x'' + c x' + k x = 0, applied to an aircraft
landing-gear strut at touchdown: characteristic roots, the three damping regimes, closed-form
solutions, the logarithmic decrement, energy decay, the Wronskian, and a reduction to a first-order
system integrated by an RK4 stepper written from scratch, with ``scipy.integrate.solve_ivp`` and
``scipy.signal.find_peaks`` as oracles.

Conventions:
    * x is measured from the *static* equilibrium (spring force balances weight), positive in compression,
      so the free equation is homogeneous. ``touchdown_response`` adds the static deflection back.
    * omega_n = sqrt(k/m) [rad/s], zeta = c / (2 sqrt(k m)), omega_d = omega_n sqrt(1 - zeta^2).

Run:  python solution.py          (prints the landing-gear table, writes landing_gear_response.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-02/lesson-02-second-order-linear-odes
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

CRITICAL_TOL = 1e-9   # |zeta - 1| below this is treated as critical damping


def landing_gear() -> dict[str, float]:
    """Return the main-gear leg model of a light twin at touchdown (SI units).

    1200 kg carried by the leg, strut stiffness 150 kN/m (natural frequency 1.78 Hz), touchdown sink
    rate 3.0 m/s, g = 9.81 m/s^2. The damper coefficient is chosen per case with ``damper_for_ratio``.
    """
    return {"m": 1200.0, "k": 150e3, "sink_rate": 3.0, "g": 9.81}


# ---------------------------------------------------------------- parameters of the equation

def natural_frequency(m: float, k: float) -> float:
    """Return omega_n = sqrt(k / m) [rad/s]. Raise ``ValueError`` unless m, k > 0."""
    if m <= 0.0 or k <= 0.0:
        raise ValueError("mass and stiffness must be positive")
    return float(np.sqrt(k / m))


def damping_ratio(m: float, c: float, k: float) -> float:
    """Return zeta = c / (2 sqrt(k m)) (dimensionless). Raise ``ValueError`` for c < 0."""
    if c < 0.0:
        raise ValueError("damping coefficient must be non-negative")
    return float(c / (2.0 * np.sqrt(k * m)))


def damper_for_ratio(m: float, k: float, zeta: float) -> float:
    """Return the damper coefficient c = 2 zeta sqrt(k m) [N s/m] that gives damping ratio ``zeta``."""
    if zeta < 0.0:
        raise ValueError("damping ratio must be non-negative")
    return float(2.0 * zeta * np.sqrt(k * m))


def damped_frequency(m: float, c: float, k: float) -> float:
    """Return omega_d = omega_n sqrt(1 - zeta^2) [rad/s] for zeta < 1, and 0.0 for zeta >= 1 (no oscillation)."""
    z = damping_ratio(m, c, k)
    return float(natural_frequency(m, k) * np.sqrt(max(0.0, 1.0 - z * z)))


def characteristic_roots(m: float, c: float, k: float) -> np.ndarray:
    """Return the two roots of m r^2 + c r + k = 0 as a complex array, sorted by imaginary part then real part."""
    if m <= 0.0 or k <= 0.0 or c < 0.0:
        raise ValueError("need m, k > 0 and c >= 0")
    disc = complex(c * c - 4.0 * m * k)
    sq = np.sqrt(disc)
    roots = np.array([(-c - sq) / (2.0 * m), (-c + sq) / (2.0 * m)])
    return roots[np.lexsort((roots.real, roots.imag))]


def classify_damping(m: float, c: float, k: float) -> str:
    """Return 'undamped', 'underdamped', 'critical' or 'overdamped' according to zeta = 0, < 1, = 1 (to CRITICAL_TOL), > 1."""
    z = damping_ratio(m, c, k)
    if z == 0.0:
        return "undamped"
    if abs(z - 1.0) < CRITICAL_TOL:
        return "critical"
    return "underdamped" if z < 1.0 else "overdamped"


# ---------------------------------------------------------------- the closed-form free response

def free_response(t: np.ndarray | float, m: float, c: float, k: float, x0: float,
                  v0: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (x(t), x'(t)) solving m x'' + c x' + k x = 0 with x(0) = x0, x'(0) = v0, in every damping regime.

    Underdamped: x = e^{-zeta wn t} [x0 cos(wd t) + (v0 + zeta wn x0)/wd sin(wd t)].
    Critical:    x = e^{-wn t} [x0 + (v0 + wn x0) t].
    Overdamped:  x = A e^{r1 t} + B e^{r2 t} with real r1 < r2 < 0, A + B = x0, A r1 + B r2 = v0.
    The three formulas agree in the limit zeta -> 1 from either side.
    """
    t = np.asarray(t, dtype=float)
    wn = natural_frequency(m, k)
    z = damping_ratio(m, c, k)
    if abs(z - 1.0) < CRITICAL_TOL:
        a = v0 + wn * x0
        e = np.exp(-wn * t)
        x = e * (x0 + a * t)
        v = e * (a - wn * (x0 + a * t))
        return x, v
    if z < 1.0:
        wd = wn * np.sqrt(1.0 - z * z)
        sig = z * wn
        b = (v0 + sig * x0) / wd
        e = np.exp(-sig * t)
        cs, sn = np.cos(wd * t), np.sin(wd * t)
        x = e * (x0 * cs + b * sn)
        v = e * ((-sig * x0 + b * wd) * cs + (-sig * b - x0 * wd) * sn)
        return x, v
    r1, r2 = characteristic_roots(m, c, k).real
    bcoef = (v0 - r1 * x0) / (r2 - r1)
    acoef = x0 - bcoef
    e1, e2 = np.exp(r1 * t), np.exp(r2 * t)
    return acoef * e1 + bcoef * e2, acoef * r1 * e1 + bcoef * r2 * e2


def mechanical_energy(x: np.ndarray, v: np.ndarray, m: float, k: float) -> np.ndarray:
    """Return E = m v^2 / 2 + k x^2 / 2 [J], the energy stored in motion and in the spring."""
    return 0.5 * m * np.asarray(v, dtype=float) ** 2 + 0.5 * k * np.asarray(x, dtype=float) ** 2


def wronskian(m: float, c: float, k: float, t: np.ndarray | float) -> np.ndarray:
    """Return the Wronskian W(t) = x1 x2' - x1' x2 of the fundamental pair x1 (x0 = 1, v0 = 0) and x2 (x0 = 0, v0 = 1).

    Abel's identity for x'' + (c/m) x' + (k/m) x = 0 gives W(t) = W(0) exp(-c t / m) with W(0) = 1: the two
    solutions are independent for all t, and the decay rate of W is the damping c/m.
    """
    x1, v1 = free_response(t, m, c, k, 1.0, 0.0)
    x2, v2 = free_response(t, m, c, k, 0.0, 1.0)
    return x1 * v2 - v1 * x2


# ---------------------------------------------------------------- reduction to first order and RK4

def to_first_order(m: float, c: float, k: float) -> np.ndarray:
    """Return the 2x2 state matrix A of d/dt (x, v) = A (x, v): [[0, 1], [-k/m, -c/m]]."""
    return np.array([[0.0, 1.0], [-k / m, -c / m]])


def rk4_step(f: Callable[[float, np.ndarray], np.ndarray], t: float, y: np.ndarray, h: float) -> np.ndarray:
    """Return one classical Runge-Kutta step of dy/dt = f(t, y)."""
    k1 = f(t, y)
    k2 = f(t + 0.5 * h, y + 0.5 * h * k1)
    k3 = f(t + 0.5 * h, y + 0.5 * h * k2)
    k4 = f(t + h, y + h * k3)
    return y + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0


def rk4_free_response(m: float, c: float, k: float, x0: float, v0: float, t_end: float,
                      n_steps: int) -> tuple[np.ndarray, np.ndarray]:
    """Integrate the first-order system with RK4 from (x0, v0); return (times, states of shape (n+1, 2))."""
    if n_steps < 1:
        raise ValueError("n_steps must be at least 1")
    A = to_first_order(m, c, k)
    f = lambda t, y: A @ y
    ts = np.linspace(0.0, t_end, n_steps + 1)
    h = t_end / n_steps
    ys = np.zeros((n_steps + 1, 2))
    ys[0] = [x0, v0]
    for i in range(n_steps):
        ys[i + 1] = rk4_step(f, ts[i], ys[i], h)
    return ts, ys


# ---------------------------------------------------------------- reading a decaying record

def find_peaks_simple(t: np.ndarray, x: np.ndarray, min_separation: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """Return (times, values) of the strict local maxima of the sampled signal x(t), in time order.

    With ``min_separation`` > 0, candidate maxima are taken in order of decreasing height and a candidate
    closer than ``min_separation`` (in t) to an already accepted peak is dropped, so noise ripples riding on
    an oscillation are rejected; use a little more than half the expected period. Interior points only.
    """
    t = np.asarray(t, dtype=float)
    x = np.asarray(x, dtype=float)
    idx = np.where((x[1:-1] > x[:-2]) & (x[1:-1] > x[2:]))[0] + 1
    if min_separation > 0.0 and idx.size:
        kept: list[int] = []
        for i in idx[np.argsort(-x[idx], kind="stable")]:
            if all(abs(t[i] - t[j]) >= min_separation for j in kept):
                kept.append(int(i))
        idx = np.array(sorted(kept), dtype=int)
    return t[idx], x[idx]


def logarithmic_decrement(peaks: np.ndarray) -> float:
    """Return delta = ln(x_i / x_{i+n}) / n averaged over the successive positive peaks ``peaks``.

    Uses the first and last peak: delta = ln(peaks[0] / peaks[-1]) / (len(peaks) - 1). Raise ``ValueError``
    with fewer than two peaks or a non-positive peak.
    """
    p = np.asarray(peaks, dtype=float)
    if p.size < 2 or np.any(p <= 0.0):
        raise ValueError("need at least two positive peaks")
    return float(np.log(p[0] / p[-1]) / (p.size - 1))


def zeta_from_decrement(delta: float) -> float:
    """Return zeta = delta / sqrt(4 pi^2 + delta^2): the exact inversion of delta = 2 pi zeta / sqrt(1 - zeta^2)."""
    return float(delta / np.sqrt(4.0 * np.pi ** 2 + delta ** 2))


def settling_time(t: np.ndarray, x: np.ndarray, band: float) -> float:
    """Return the first time after which |x| stays within ``band`` for the rest of the record (inf if never)."""
    t = np.asarray(t, dtype=float)
    out = np.abs(np.asarray(x, dtype=float)) > band
    if not np.any(out):
        return float(t[0])
    last = int(np.max(np.where(out)[0]))
    return float("inf") if last == t.size - 1 else float(t[last + 1])


def synthetic_decay_record(fn: float, zeta: float, a0: float, duration: float, fs: float,
                           noise_std: float, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Return (t, a) for a SYNTHETIC accelerometer record of one decaying mode plus Gaussian noise.

    a(t) = a0 exp(-zeta wn t) cos(wd t) + N(0, noise_std^2), wn = 2 pi fn, sampled at ``fs`` for ``duration``
    seconds, with ``numpy.random.default_rng(seed)``. This is simulated data, not a measurement.
    """
    rng = np.random.default_rng(seed)
    t = np.arange(0.0, duration, 1.0 / fs)
    wn = 2.0 * np.pi * fn
    wd = wn * np.sqrt(1.0 - zeta * zeta)
    return t, a0 * np.exp(-zeta * wn * t) * np.cos(wd * t) + rng.normal(0.0, noise_std, t.size)


def estimate_damping_from_record(t: np.ndarray, x: np.ndarray, f_guess: float,
                                 height_threshold: float) -> dict[str, float]:
    """Estimate (f_d, zeta) of a decaying oscillation from its positive peaks above ``height_threshold``.

    Peaks are found with ``find_peaks_simple`` using a minimum separation of 0.7 / f_guess. Returns
    {'n_peaks', 'f_d' (from the median peak spacing), 'zeta_decrement' (first and last peak),
    'zeta_fit' (least-squares slope of ln(peak) against time, divided by -2 pi f_d / sqrt(1 - zeta^2) ~ -wn)}.
    Raise ``ValueError`` with fewer than three usable peaks.
    """
    tp, xp = find_peaks_simple(t, x, min_separation=0.7 / f_guess)
    keep = xp > height_threshold
    tp, xp = tp[keep], xp[keep]
    if xp.size < 3:
        raise ValueError("fewer than three peaks above the threshold")
    f_d = 1.0 / float(np.median(np.diff(tp)))
    z_dec = zeta_from_decrement(logarithmic_decrement(xp))
    A = np.vstack([tp, np.ones_like(tp)]).T
    slope = np.linalg.lstsq(A, np.log(xp), rcond=None)[0][0]
    sigma = -slope                                    # decay rate zeta * wn
    wd = 2.0 * np.pi * f_d
    z_fit = sigma / np.sqrt(sigma ** 2 + wd ** 2)     # zeta wn / wn, with wn^2 = sigma^2 + wd^2
    return {"n_peaks": float(xp.size), "f_d": f_d, "zeta_decrement": z_dec, "zeta_fit": float(z_fit)}


# ---------------------------------------------------------------- the landing gear

def touchdown_response(t: np.ndarray | float, m: float, k: float, zeta: float, sink_rate: float,
                       g: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (compression x(t) from the unloaded strut, x'(t)) after touchdown at the given sink rate.

    The strut is unloaded at t = 0 (x = 0) and compressing at ``sink_rate``. About the static deflection
    x_s = m g / k the motion is free, so x(t) = x_s + y(t) with y(0) = -x_s, y'(0) = sink_rate.
    """
    c = damper_for_ratio(m, k, zeta)
    xs = m * g / k
    y, v = free_response(t, m, c, k, -xs, sink_rate)
    return xs + y, v


def max_stroke(m: float, k: float, zeta: float, sink_rate: float, g: float) -> tuple[float, float]:
    """Return (time of maximum compression, maximum compression) of ``touchdown_response``.

    The maximum is the first zero of the velocity, bracketed between 0 and the first time at which the
    closed-form velocity is negative and located with ``scipy.optimize.brentq``. If the velocity never
    changes sign within 200 / omega_n (a heavily overdamped strut creeping up to the static deflection)
    the pair (inf, m g / k) is returned: the compression approaches the static value from below.
    """
    wn = natural_frequency(m, k)
    vel = lambda tt: float(touchdown_response(tt, m, k, zeta, sink_rate, g)[1])
    t_hi = 0.1 / wn
    while vel(t_hi) > 0.0:
        t_hi *= 1.5
        if t_hi > 200.0 / wn:
            return float("inf"), float(m * g / k)
    t_star = brentq(vel, 0.0, t_hi, xtol=1e-14)
    return float(t_star), float(touchdown_response(t_star, m, k, zeta, sink_rate, g)[0])


def strut_force(t: np.ndarray | float, m: float, k: float, zeta: float, sink_rate: float,
                g: float) -> np.ndarray:
    """Return the force the strut applies to the airframe, F = k x + c x' [N], along ``touchdown_response``.

    At t = 0 the spring force is zero and the damper force is c * sink_rate: a stiff damper hits the
    airframe hardest at the instant of touchdown, which is the design penalty of heavy damping.
    """
    c = damper_for_ratio(m, k, zeta)
    x, v = touchdown_response(t, m, k, zeta, sink_rate, g)
    return k * x + c * v


def loses_contact(x: np.ndarray) -> bool:
    """Return True if the compression ever becomes negative: the linear model is then asking the strut to pull, i.e. the wheel bounces."""
    return bool(np.any(np.asarray(x, dtype=float) < 0.0))


# ---------------------------------------------------------------- the example

def main(out_dir: str | Path = ".") -> None:
    """Print the landing-gear touchdown table for several damping ratios and save landing_gear_response.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.integrate import solve_ivp

    out_dir = Path(out_dir)
    p = landing_gear()
    m, k, v0, g = p["m"], p["k"], p["sink_rate"], p["g"]
    wn = natural_frequency(m, k)
    xs = m * g / k
    print(f"omega_n = {wn:.3f} rad/s ({wn / 2 / np.pi:.3f} Hz), period {2 * np.pi / wn:.3f} s, static deflection {xs * 1e3:.1f} mm")
    print(f"touchdown kinetic energy {0.5 * m * v0 ** 2 / 1e3:.2f} kJ; undamped stroke estimate from energy"
          f" sqrt(x_s^2 + m v0^2 / k) + x_s = {(np.sqrt(xs ** 2 + m * v0 ** 2 / k) + xs) * 1e3:.0f} mm")
    t = np.linspace(0.0, 3.0, 3001)
    print(f"{'zeta':>5} {'c [kN s/m]':>11} {'regime':>12} {'roots [1/s]':>28} {'max stroke [mm]':>16} {'at t [s]':>9}"
          f" {'bounce?':>8} {'settle 2% [s]':>14} {'peak force [kN]':>16}")
    rows = []
    for z in (0.0, 0.2, 0.5, 0.7, 1.0, 2.0):
        c = damper_for_ratio(m, k, z)
        r = characteristic_roots(m, c, k)
        t_s, x_s = max_stroke(m, k, z, v0, g)
        x, v = touchdown_response(t, m, k, z, v0, g)
        settle = settling_time(t, x - xs, 0.02 * xs)
        rows.append((z, x))
        root_str = f"{r[0]:.2f}, {r[1]:.2f}" if np.iscomplexobj(r) and abs(r.imag).max() > 0 else f"{r[0].real:.2f}, {r[1].real:.2f}"
        f_peak = np.max(strut_force(t, m, k, z, v0, g))
        print(f"{z:5.2f} {c / 1e3:11.2f} {classify_damping(m, c, k):>12} {root_str:>28} {x_s * 1e3:16.1f} {t_s:9.3f}"
              f" {str(loses_contact(x)):>8} {settle:14.3f} {f_peak / 1e3:16.1f}")

    # Decrement from the zeta = 0.2 record, RK4 and solve_ivp against the closed form.
    z = 0.2
    c = damper_for_ratio(m, k, z)
    x, v = free_response(t, m, c, k, -xs, v0)
    tp, xp = find_peaks_simple(t, x)
    delta = logarithmic_decrement(xp)
    print(f"\nzeta = 0.2 record: {xp.size} positive peaks, log decrement {delta:.4f}, zeta recovered {zeta_from_decrement(delta):.4f},"
          f" damped period from peaks {np.mean(np.diff(tp)):.4f} s vs 2 pi / omega_d = {2 * np.pi / damped_frequency(m, c, k):.4f} s")
    for n in (300, 600, 1200, 2400):
        ts, ys = rk4_free_response(m, c, k, -xs, v0, 3.0, n)
        xe, _ = free_response(ts, m, c, k, -xs, v0)
        print(f"RK4 with {n:5d} steps: max |x - x_exact| = {np.max(np.abs(ys[:, 0] - xe)):.3e} m")
    sol = solve_ivp(lambda tt, y: to_first_order(m, c, k) @ y, (0.0, 3.0), [-xs, v0], rtol=1e-10, atol=1e-12, t_eval=t)
    print(f"solve_ivp (rtol 1e-10): max |x - x_exact| = {np.max(np.abs(sol.y[0] - x)):.3e} m")
    E = mechanical_energy(x, v, m, k)
    print(f"energy: {E[0]:.1f} J at touchdown, {E[-1]:.3f} J after 3 s; dissipated fraction {1 - E[-1] / E[0]:.6f}")
    print(f"Wronskian at t = 1 s: {float(wronskian(m, c, k, 1.0)):.6f} vs exp(-c t / m) = {np.exp(-c / m):.6f}")

    print("\nSYNTHETIC panel record: 42 Hz, zeta = 0.015, 25 m/s^2 initial, 2000 Hz, 2 s (simulated, not measured)")
    for noise in (0.1, 0.5, 1.0):
        for thr in (5.0, 20.0):
            tt, acc = synthetic_decay_record(42.0, 0.015, 25.0, 2.0, 2000.0, noise, seed=202)
            try:
                est = estimate_damping_from_record(tt, acc, 42.0, thr * noise)
                print(f"  noise {noise:4.1f} m/s^2, threshold {thr:4.1f} sigma: {int(est['n_peaks']):3d} peaks, f_d = {est['f_d']:.2f} Hz,"
                      f" zeta (first/last) = {est['zeta_decrement']:.4f}, zeta (fit) = {est['zeta_fit']:.4f}")
            except ValueError as exc:
                print(f"  noise {noise:4.1f} m/s^2, threshold {thr:4.1f} sigma: {exc}")

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    for z, x in rows:
        axes[0].plot(t, x * 1e3, label=f"zeta = {z}")
    axes[0].axhline(xs * 1e3, color="k", lw=0.6, ls=":")
    axes[0].axhline(0.0, color="k", lw=0.6)
    axes[0].set_xlabel("t [s]"); axes[0].set_ylabel("strut compression [mm]"); axes[0].set_xlim(0, 2.0)
    axes[0].legend(fontsize=8); axes[0].grid(lw=0.4)
    axes[0].set_title("touchdown at 3 m/s, 1200 kg, 150 kN/m", fontsize=9)
    zs = np.linspace(0.05, 3.0, 300)
    strokes = np.array([max_stroke(m, k, zz, v0, g)[1] for zz in zs])
    settles = np.array([settling_time(t, touchdown_response(t, m, k, zz, v0, g)[0] - xs, 0.02 * xs) for zz in zs])
    axes[1].plot(zs, strokes * 1e3, label="max stroke [mm]")
    axes[1].plot(zs, settles * 100, label="2 % settling time [x 0.01 s]")
    forces = np.array([np.max(strut_force(t, m, k, zz, v0, g)) for zz in zs])
    axes[1].plot(zs, forces / 1e3, label="peak strut force [kN]")
    axes[1].set_xlabel("damping ratio zeta"); axes[1].legend(fontsize=8); axes[1].grid(lw=0.4)
    fig.tight_layout()
    fig.savefig(out_dir / "landing_gear_response.png", dpi=150)
    print(f"saved {out_dir / 'landing_gear_response.png'}")


if __name__ == "__main__":
    main()
