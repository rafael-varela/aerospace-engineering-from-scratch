"""P09.10 Nyquist & Stability Margins: starter code.

The Nyquist criterion counts the unstable closed-loop poles of 1 + L(s) from the open-loop frequency
response alone: map a contour enclosing the right half plane through L, count how many times the image
encircles -1 clockwise (N), add the number of open-loop poles in the right half plane (P), and
Z = N + P closed-loop poles are unstable. Because it needs only L(j omega), it works for loops that
contain a pure time delay, where there is no characteristic polynomial to factor. The distances of the
curve from -1 are the gain, phase and delay margins.

Three representative aerospace loops (typical numbers, not one vehicle's data):

* a control-surface servo, L = K / (s (tau_m s + 1)(tau_f s + 1)), tau_m = 0.05 s (motor and
  linkage), tau_f = 0.01 s (driver and sensor filter): one pole on the imaginary axis (s = 0);
* a reaction-wheel satellite, single axis J theta'' = torque, with a PD controller and a derivative
  filter, L = Kp (Td s + 1) e^{-sT} / (J s^2 (tau s + 1)): a double pole at s = 0 and an optional
  computation and actuation delay T;
* a statically unstable airframe (pitch, with a relaxed-stability margin) under PD control,
  L = K (s + 3) / ((s - 2)(s + 10)): P = 1, so stability needs one counterclockwise encirclement.

This module builds, from scratch:

* the Nyquist D-contour with small semicircular indentations to the right of imaginary-axis poles;
* the winding number of the image about -1 (sum of wrapped angle increments) and Z = N + P;
* the gain margin (at the phase crossover, L real and negative), the phase margin (at the gain
  crossover, |L| = 1) and the delay margin PM / omega_gc;
* a Pade approximant of the delay, used only in the tests as an independent route to the closed-loop
  poles of the delayed loop.

Units: s, rad/s; the servo gain K in 1/s; satellite J in kg m^2, Kp in N m/rad.

Run:  python main.py [out_dir]     (prints the lesson's tables, writes nyquist_margins.png)
Test: pytest curriculum/phase-09/lesson-10-nyquist-stability-margins
Complete every function whose body is `raise NotImplementedError`; every other function
and main() are given and complete.
"""
from __future__ import annotations

import sys
from math import factorial
from pathlib import Path

import numpy as np

TAU_M, TAU_F = 0.05, 0.01           # s, servo motor and filter
J_SAT, TD_SAT, TAU_SAT = 10.0, 10.0, 1.0   # kg m^2, s, s (satellite PD with derivative filter)
KP_SAT = 0.3162                      # N m / rad, puts the gain crossover near 0.316 rad/s


# ---------------------------------------------------------------- given helpers

def poly_add(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Sum of two polynomials (highest power first)."""
    p, q = np.atleast_1d(np.asarray(p, dtype=float)), np.atleast_1d(np.asarray(q, dtype=float))
    n = max(p.size, q.size)
    return np.concatenate([np.zeros(n - p.size), p]) + np.concatenate([np.zeros(n - q.size), q])


def tf_eval(num: np.ndarray, den: np.ndarray, s: np.ndarray | complex) -> np.ndarray:
    """num(s)/den(s) at complex s."""
    return np.polyval(num, s) / np.polyval(den, s)


def servo_loop(K: float, tau_m: float = TAU_M, tau_f: float = TAU_F) -> tuple[np.ndarray, np.ndarray]:
    """(num, den) of L = K / (s (tau_m s + 1)(tau_f s + 1))."""
    return np.array([K]), np.convolve(np.convolve([1.0, 0.0], [tau_m, 1.0]), [tau_f, 1.0])


def satellite_loop(kp: float = KP_SAT, td: float = TD_SAT, tau: float = TAU_SAT, J: float = J_SAT) -> tuple[np.ndarray, np.ndarray]:
    """(num, den) of L = kp (td s + 1) / (J s^2 (tau s + 1)) (without the delay)."""
    return kp * np.array([td, 1.0]), np.convolve([J, 0.0, 0.0], [tau, 1.0])


def unstable_airframe_loop(K: float) -> tuple[np.ndarray, np.ndarray]:
    """(num, den) of L = K (s + 3) / ((s - 2)(s + 10))."""
    return K * np.array([1.0, 3.0]), np.convolve([1.0, -2.0], [1.0, 10.0])


def pade_delay(T: float, order: int = 6) -> tuple[np.ndarray, np.ndarray]:
    """(num, den) of the [order/order] Pade approximant of e^{-sT} (highest power first).

    den = sum_k c_k (sT)^k, num = sum_k c_k (-sT)^k, c_k = (2n - k)! n! / ((2n)! k! (n - k)!).
    """
    n = order
    c = np.array([factorial(2 * n - k) * factorial(n) / (factorial(2 * n) * factorial(k) * factorial(n - k)) for k in range(n + 1)])
    den = (c * T ** np.arange(n + 1))[::-1]
    num = (c * (-T) ** np.arange(n + 1))[::-1]
    return num, den


# ---------------------------------------------------------------- the contour and its image

def imag_axis_poles(den: np.ndarray, tol: float = 1e-8) -> np.ndarray:
    """Sorted distinct imaginary parts omega_k of the poles of 1/den on the imaginary axis (|Re| < tol)."""
    raise NotImplementedError


def nyquist_contour(den: np.ndarray, R: float = 1e4, eps: float = 1e-4, n_axis: int = 40000, n_arc: int = 400) -> np.ndarray:
    """Points of the clockwise Nyquist D-contour: up the imaginary axis from -jR to +jR, detouring
    around each imaginary-axis pole j omega_k along s = j omega_k + eps e^{j phi}, phi from -90 to +90
    deg (to the right, so the pole stays outside), then back along s = R e^{j phi}, phi from +90 to -90 deg.

    The axis points are spaced uniformly in asinh(omega) (logarithmic for large |omega|), with extra
    points packed logarithmically near each indentation.
    """
    raise NotImplementedError


def loop_on_contour(num: np.ndarray, den: np.ndarray, s: np.ndarray, delay: float = 0.0) -> np.ndarray:
    """L(s) e^{-s delay} at the contour points."""
    raise NotImplementedError


def winding_number(z: np.ndarray, point: complex = -1.0 + 0j) -> float:
    """Net number of counterclockwise turns of the closed curve z (last point joined to the first)
    about `point`: the sum of the wrapped angle increments divided by 2 pi."""
    raise NotImplementedError


def rhp_poles(den: np.ndarray) -> int:
    """P: number of open-loop poles with positive real part."""
    raise NotImplementedError


def nyquist_count(num: np.ndarray, den: np.ndarray, delay: float = 0.0, **contour_kw) -> tuple[int, int, int]:
    """Return (N, P, Z): N clockwise encirclements of -1 by L(s) e^{-s delay} on the D-contour,
    P open-loop right-half-plane poles, Z = N + P closed-loop right-half-plane poles."""
    raise NotImplementedError


# ---------------------------------------------------------------- margins

def _freq_grid(w_min: float = 1e-4, w_max: float = 1e4, n: int = 20000) -> np.ndarray:
    return np.geomspace(w_min, w_max, n)


def phase_margin(num: np.ndarray, den: np.ndarray, delay: float = 0.0, w: np.ndarray | None = None) -> tuple[float, float]:
    """Return (PM in degrees, omega_gc) at the lowest-frequency gain crossover |L(j omega)| = 1 (refined by
    bisection on log |L|). PM = 180 deg + arg L(j omega_gc), wrapped to (-180, 180]. (nan, nan) if |L|
    never crosses 1 on the grid."""
    raise NotImplementedError


def phase_crossovers(num: np.ndarray, den: np.ndarray, delay: float = 0.0, w: np.ndarray | None = None,
                     min_mag: float = 1e-3) -> list[tuple[float, float]]:
    """Frequencies omega > 0 where L(j omega) crosses the negative real axis, with |L| there:
    [(omega_pc, |L(j omega_pc)|), ...], found by sign changes of Im L with Re L < 0 and bisection.
    Crossings where |L| < min_mag (a gain margin above 60 dB) are ignored: they cannot matter, and a
    delayed loop spirals through thousands of them."""
    raise NotImplementedError


def gain_margin(num: np.ndarray, den: np.ndarray, delay: float = 0.0, w: np.ndarray | None = None) -> tuple[float, float]:
    """Return (GM as a factor, omega_pc): 1/|L| at the phase crossover where |L| is largest (the first
    point on the negative real axis the curve would push through -1 as the gain rises).
    (inf, nan) if the curve never crosses the negative real axis for omega > 0."""
    raise NotImplementedError


def delay_margin(pm_deg: float, w_gc: float) -> float:
    """Largest extra pure delay [s] the loop tolerates: PM (in rad) / omega_gc."""
    raise NotImplementedError


def margins(num: np.ndarray, den: np.ndarray, delay: float = 0.0) -> dict[str, float]:
    """All three margins: GM (factor), GM_db, w_pc, PM (deg), w_gc, DM (s)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the lesson's numbers

def closed_loop_rhp_roots(num: np.ndarray, den: np.ndarray, delay: float = 0.0, pade_order: int = 8) -> int:
    """Closed-loop right-half-plane roots from the characteristic polynomial (the delay replaced by
    its Pade approximant). The independent route against which Z = N + P is checked."""
    if delay > 0:
        pn, pd = pade_delay(delay, pade_order)
        num, den = np.convolve(num, pn), np.convolve(den, pd)
    return int(np.sum(np.roots(poly_add(den, num)).real > 1e-9))


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's tables and save nyquist_margins.png in out_dir."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    Kc = (TAU_M + TAU_F) / (TAU_M * TAU_F)
    print(f"servo L = K / (s ({TAU_M} s + 1)({TAU_F} s + 1)); Routh: stable for K < (tau_m + tau_f)/(tau_m tau_f) = {Kc:.1f} 1/s, "
          f"phase crossover 1/sqrt(tau_m tau_f) = {1 / np.sqrt(TAU_M * TAU_F):.3f} rad/s")
    for K in (5.0, 20.0, 60.0, 119.0, 121.0, 200.0):
        num, den = servo_loop(K)
        N, P, Z = nyquist_count(num, den)
        m = margins(num, den)
        print(f"  K = {K:6.1f}: N = {N:+d}, P = {P}, Z = {Z}, roots in RHP {closed_loop_rhp_roots(num, den)}; "
              f"GM = {m['GM']:.4f} ({m['GM_db']:+.2f} dB) at {m['w_pc']:.3f} rad/s, PM = {m['PM']:+.2f} deg at {m['w_gc']:.3f} rad/s, "
              f"DM = {1e3 * m['DM']:.2f} ms")

    num, den = satellite_loop()
    m = margins(num, den)
    print(f"\nsatellite L = {KP_SAT} ({TD_SAT} s + 1) / ({J_SAT} s^2 ({TAU_SAT} s + 1)), double pole at 0: "
          f"GM = {m['GM']}, PM = {m['PM']:.2f} deg at {m['w_gc']:.4f} rad/s, DM = {m['DM']:.3f} s")
    for T in (0.0, 0.5, 0.9 * m["DM"], 1.1 * m["DM"], 2 * m["DM"]):
        N, P, Z = nyquist_count(num, den, delay=T)
        mT = margins(num, den, delay=T)
        print(f"  delay {T:6.3f} s: N = {N:+d}, P = {P}, Z = {Z}; Pade(8) closed loop RHP roots {closed_loop_rhp_roots(num, den, T)}; "
              f"PM {mT['PM']:+7.2f} deg")
    s_gc = 1j * m["w_gc"]
    print(f"  |1 + L(j w_gc) e^(-j w_gc DM)| = {abs(1 + tf_eval(num, den, s_gc) * np.exp(-s_gc * m['DM'])):.2e} (the curve passes through -1)")

    print("\nunstable airframe L = K (s + 3)/((s - 2)(s + 10)), P = 1; closed loop s^2 + (8 + K) s + (3K - 20): stable for K > 20/3 = 6.667")
    for K in (3.0, 6.0, 7.0, 20.0, 100.0):
        num, den = unstable_airframe_loop(K)
        N, P, Z = nyquist_count(num, den)
        pm, wgc = phase_margin(num, den)
        f0 = -1 / tf_eval(num, den, 0.0).real
        verdict = (f"gain may fall to {f0:.3f} times this value ({20 * np.log10(f0):+.2f} dB) before instability" if f0 < 1
                   else f"gain must rise by a factor {f0:.3f} for stability")
        print(f"  K = {K:6.1f}: L(0) = {tf_eval(num, den, 0.0).real:+.3f}, N = {N:+d}, P = {P}, Z = {Z}, roots in RHP "
              f"{closed_loop_rhp_roots(num, den)}; PM = {pm:+.2f} deg at {wgc:.3f} rad/s; {verdict}")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, (name, (num, den), delay, lim) in zip(axes, [
        ("servo, K = 20", servo_loop(20.0), 0.0, 3.0),
        (f"satellite, delay 0 and {m['DM']:.2f} s", satellite_loop(), None, 3.0),
        ("unstable airframe, K = 20", unstable_airframe_loop(20.0), 0.0, 3.0)]):
        s = nyquist_contour(den)
        for T, style in (((0.0, "-"), (m["DM"], "--")) if delay is None else ((delay, "-"),)):
            Lc = loop_on_contour(num, den, s, T)
            Lc = np.where(np.abs(Lc) > 50, np.nan, Lc)
            ax.plot(Lc.real, Lc.imag, style, lw=1, label=f"delay {T:.2f} s")
        ax.plot([-1], [0], "r+", ms=12, mew=2)
        uc = np.exp(1j * np.linspace(0, 2 * np.pi, 200))
        ax.plot(uc.real, uc.imag, "k:", lw=0.5)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
        ax.set_xlabel("Re L"); ax.set_ylabel("Im L"); ax.set_title(name, fontsize=9); ax.legend(fontsize=7); ax.grid(lw=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "nyquist_margins.png", dpi=130)
    plt.close(fig)
    print(f"saved {out_dir / 'nyquist_margins.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
