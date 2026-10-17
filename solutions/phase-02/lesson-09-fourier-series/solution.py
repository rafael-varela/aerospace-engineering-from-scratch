"""P02.09 Fourier Series: reference solution.

Square-wave synthesis: Fourier coefficients computed from scratch by the trapezoid rule on one period,
the closed-form coefficients of the square and triangle waves, partial sums and the Gibbs overshoot,
Parseval's identity, the convergence rate of the truncation error, the complex form, and the
steady-state response of a lightly damped structural mode to a periodic (on-off thruster) torque,
harmonic by harmonic with the complex gain of P02.07. Oracles: scipy.integrate.quad for the
coefficients, SymPy for one symbolic integral, scipy.integrate.solve_ivp for the response.

Conventions:
    * f(t) periodic with period T, written  f = a0/2 + sum_{n>=1} [a_n cos(n w0 t) + b_n sin(n w0 t)],  w0 = 2 pi / T,
      a_n = (2/T) int_0^T f cos(n w0 t) dt,  b_n = (2/T) int_0^T f sin(n w0 t) dt.
    * Coefficient arrays are indexed from 1: an[0] is a_1. a0 is returned separately.
    * Complex form c_n = (a_n - i b_n)/2 for n >= 1, c_0 = a0/2, c_{-n} = conj(c_n).
    * The structural mode is  q'' + 2 zeta wn q' + wn^2 q = f(t)  (unit modal mass), so the static
      deflection for unit forcing is 1 / wn^2.

Run:  python solution.py          (prints the thruster-pulse example, writes fourier_series.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-02/lesson-09-fourier-series
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np


def thruster_boom_case() -> dict[str, float]:
    """Return the nominal case: an on-off thruster pulse train driving a solar-array boom mode.

    The thruster alternates sign every half period (a square wave of period T = 1.0 s, normalised
    amplitude 1 in the modal equation); the boom's first bending mode is at 5.0 Hz with damping ratio 0.01,
    i.e. exactly on the fifth harmonic of the pulse train.
    """
    return {"T": 1.0, "f_mode": 5.0, "zeta": 0.01}


# ---------------------------------------------------------------- the waveforms

def square_wave(t: np.ndarray, T: float, amplitude: float = 1.0) -> np.ndarray:
    """Return the odd square wave: +A on (0, T/2), -A on (T/2, T), repeated; zero exactly at the jumps."""
    phase = np.mod(np.asarray(t, dtype=float), T) / T
    return amplitude * np.where(phase == 0.0, 0.0, np.where(phase < 0.5, 1.0, np.where(phase == 0.5, 0.0, -1.0)))


def triangle_wave(t: np.ndarray, T: float, amplitude: float = 1.0) -> np.ndarray:
    """Return the even triangle wave with value +A at t = 0 and -A at t = T/2 (continuous, slope +-4A/T)."""
    phase = np.mod(np.asarray(t, dtype=float) + 0.5 * T, T) / T       # 0 at t = -T/2, 0.5 at t = 0
    return amplitude * (1.0 - 4.0 * np.abs(phase - 0.5))


def square_wave_coefficients(n_max: int, amplitude: float = 1.0) -> tuple[float, np.ndarray, np.ndarray]:
    """Return (a0, an, bn) of the odd square wave up to harmonic n_max: a_n = 0, b_n = 4A/(pi n) for odd n, 0 for even n."""
    n = np.arange(1, n_max + 1)
    bn = np.where(n % 2 == 1, 4.0 * amplitude / (np.pi * n), 0.0)
    return 0.0, np.zeros(n_max), bn


def triangle_wave_coefficients(n_max: int, amplitude: float = 1.0) -> tuple[float, np.ndarray, np.ndarray]:
    """Return (a0, an, bn) of the even triangle wave: b_n = 0, a_n = 8A/(pi^2 n^2) for odd n, 0 for even n."""
    n = np.arange(1, n_max + 1)
    an = np.where(n % 2 == 1, 8.0 * amplitude / (np.pi ** 2 * n ** 2), 0.0)
    return 0.0, an, np.zeros(n_max)


# ---------------------------------------------------------------- coefficients from scratch

def fourier_coefficients(f: Callable[[np.ndarray], np.ndarray], T: float, n_max: int,
                         n_samples: int = 4096) -> tuple[float, np.ndarray, np.ndarray]:
    """Return (a0, an, bn) up to harmonic n_max by the trapezoid rule on n_samples equal subintervals of one period.

    For a periodic integrand the trapezoid rule on a full period is spectrally accurate when f is smooth,
    and first-order accurate (error ~ 1/n_samples) when f has jumps. Sample points are the midpoints of
    the subintervals, so a jump at t = 0 or T/2 is never sampled exactly. Raise ``ValueError`` if
    n_samples <= 2 n_max (fewer than two samples per cycle of the highest harmonic: aliasing).
    """
    if n_samples <= 2 * n_max:
        raise ValueError("need more than 2 n_max samples per period to resolve harmonic n_max")
    t = (np.arange(n_samples) + 0.5) * T / n_samples
    w0 = 2.0 * np.pi / T
    ft = np.asarray(f(t), dtype=float)
    a0 = 2.0 * np.mean(ft)
    n = np.arange(1, n_max + 1)
    cos_nt = np.cos(np.outer(n, w0 * t))
    sin_nt = np.sin(np.outer(n, w0 * t))
    an = 2.0 * (cos_nt @ ft) / n_samples
    bn = 2.0 * (sin_nt @ ft) / n_samples
    return float(a0), an, bn


def partial_sum(t: np.ndarray, a0: float, an: np.ndarray, bn: np.ndarray, T: float, n_terms: int | None = None) -> np.ndarray:
    """Return S_N(t) = a0/2 + sum_{n=1}^{N} [a_n cos(n w0 t) + b_n sin(n w0 t)], N = n_terms (default: all coefficients given)."""
    t = np.asarray(t, dtype=float)
    N = len(an) if n_terms is None else int(n_terms)
    if N > len(an) or N < 0:
        raise ValueError("n_terms exceeds the number of coefficients supplied")
    w0 = 2.0 * np.pi / T
    n = np.arange(1, N + 1)
    arg = np.outer(n, w0 * t)
    return 0.5 * a0 + an[:N] @ np.cos(arg) + bn[:N] @ np.sin(arg)


def complex_coefficients(a0: float, an: np.ndarray, bn: np.ndarray) -> np.ndarray:
    """Return c_n for n = -N..N as an array of length 2N + 1 (index N is c_0): c_0 = a0/2, c_n = (a_n - i b_n)/2, c_{-n} = conj(c_n)."""
    an = np.asarray(an, dtype=float)
    bn = np.asarray(bn, dtype=float)
    pos = 0.5 * (an - 1j * bn)
    return np.concatenate([np.conj(pos[::-1]), [0.5 * a0 + 0j], pos])


def from_complex_coefficients(c: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """Invert complex_coefficients: a0 = 2 c_0, a_n = 2 Re c_n, b_n = -2 Im c_n. Raise ``ValueError`` for even length."""
    c = np.asarray(c, dtype=complex)
    if c.size % 2 == 0:
        raise ValueError("expected an odd-length array c_{-N..N}")
    N = c.size // 2
    pos = c[N + 1:]
    return float(2.0 * c[N].real), 2.0 * pos.real, -2.0 * pos.imag


# ---------------------------------------------------------------- Gibbs, Parseval, convergence

def gibbs_overshoot(n_terms: int, T: float = 1.0) -> float:
    """Return the maximum of the N-term partial sum of the unit square wave on (0, T/2), N = n_terms.

    Tends to (2/pi) Si(pi) = 1.17898 as N grows, i.e. an overshoot of 0.179 = 8.95 % of the jump height 2: the
    overshoot does not shrink as N grows, it only narrows.
    Evaluated on a grid of 200 points per half-cycle of the highest harmonic near the jump at t = 0.
    """
    a0, an, bn = square_wave_coefficients(n_terms)
    t = np.linspace(0.0, T / 4.0, 200 * max(n_terms, 1))
    return float(np.max(partial_sum(t, a0, an, bn, T)))


def parseval_power(a0: float, an: np.ndarray, bn: np.ndarray) -> float:
    """Return the mean square (1/T) int_0^T f^2 dt computed from the coefficients: a0^2/4 + (1/2) sum (a_n^2 + b_n^2)."""
    return 0.25 * a0 * a0 + 0.5 * float(np.sum(np.asarray(an) ** 2 + np.asarray(bn) ** 2))


def mean_square(f: Callable[[np.ndarray], np.ndarray], T: float, n_samples: int = 4096) -> float:
    """Return (1/T) int_0^T f^2 dt by the midpoint rule with n_samples points."""
    t = (np.arange(n_samples) + 0.5) * T / n_samples
    return float(np.mean(np.asarray(f(t), dtype=float) ** 2))


def truncation_error_rms(f: Callable[[np.ndarray], np.ndarray], T: float, a0: float, an: np.ndarray, bn: np.ndarray,
                         n_terms: int, n_samples: int = 8192) -> float:
    """Return the rms of f - S_N over one period (midpoint rule). By Parseval it equals sqrt((1/2) sum_{n>N} (a_n^2 + b_n^2))."""
    t = (np.arange(n_samples) + 0.5) * T / n_samples
    return float(np.sqrt(np.mean((np.asarray(f(t), dtype=float) - partial_sum(t, a0, an, bn, T, n_terms)) ** 2)))


def convergence_rate(errors: np.ndarray, n_terms: np.ndarray) -> float:
    """Return the slope p of a least-squares fit log(error) = const - p log(N): error ~ N^{-p}."""
    x = np.log(np.asarray(n_terms, dtype=float))
    y = np.log(np.asarray(errors, dtype=float))
    A = np.column_stack([np.ones_like(x), -x])
    coef = np.linalg.lstsq(A, y, rcond=None)[0]
    return float(coef[1])


# ---------------------------------------------------------------- response to periodic forcing

def complex_gain(r: np.ndarray, zeta: float) -> np.ndarray:
    """Return H(r) = 1 / (1 - r^2 + 2 i zeta r) (P02.07): displacement per static deflection at frequency ratio r."""
    r = np.asarray(r, dtype=float)
    return 1.0 / (1.0 - r * r + 2j * zeta * r)


def response_coefficients(a0: float, an: np.ndarray, bn: np.ndarray, T: float, wn: float, zeta: float) -> tuple[float, np.ndarray, np.ndarray]:
    """Return the Fourier coefficients (A0, An, Bn) of the steady-state q(t) of q'' + 2 zeta wn q' + wn^2 q = f(t).

    Each harmonic is a phasor problem: input phasor Z_n = a_n - i b_n at omega = n w0, output phasor
    H(n w0 / wn) Z_n / wn^2, so A_n = Re(.), B_n = -Im(.); the mean responds statically, A0 = a0 / wn^2.
    """
    an = np.asarray(an, dtype=float)
    bn = np.asarray(bn, dtype=float)
    n = np.arange(1, an.size + 1)
    r = n * (2.0 * np.pi / T) / wn
    out = complex_gain(r, zeta) * (an - 1j * bn) / wn ** 2
    return a0 / wn ** 2, out.real, -out.imag


def harmonic_amplitudes(an: np.ndarray, bn: np.ndarray) -> np.ndarray:
    """Return sqrt(a_n^2 + b_n^2), the amplitude of each harmonic."""
    return np.sqrt(np.asarray(an, dtype=float) ** 2 + np.asarray(bn, dtype=float) ** 2)


def dominant_harmonic(an: np.ndarray, bn: np.ndarray) -> int:
    """Return the harmonic number (1-based) with the largest amplitude."""
    return int(np.argmax(harmonic_amplitudes(an, bn))) + 1


def forced_mode_rhs(t: float, y: np.ndarray, f: Callable[[float], float], wn: float, zeta: float) -> np.ndarray:
    """Return dy/dt for y = (q, q') of q'' + 2 zeta wn q' + wn^2 q = f(t): the solve_ivp check of the harmonic sum."""
    return np.array([y[1], f(t) - 2.0 * zeta * wn * y[1] - wn * wn * y[0]])


def sympy_square_wave_coefficient(n: int, T: float = 1.0):
    """Return b_n of the unit odd square wave by SymPy's exact integration, (2/T)[int_0^{T/2} sin - int_{T/2}^T sin].

    An independent symbolic check of 4/(pi n) for odd n and 0 for even n; returns a SymPy expression.
    """
    import sympy as sp

    t = sp.symbols("t", real=True)
    Ts = sp.nsimplify(T)
    w0 = 2 * sp.pi / Ts
    expr = (2 / Ts) * (sp.integrate(sp.sin(n * w0 * t), (t, 0, Ts / 2)) - sp.integrate(sp.sin(n * w0 * t), (t, Ts / 2, Ts)))
    return sp.simplify(expr)


def main(out_dir: str | Path = ".") -> None:
    """Print the square-wave synthesis and thruster-boom analysis; save fourier_series.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.integrate import quad, solve_ivp

    out_dir = Path(out_dir)
    case = thruster_boom_case()
    T, wn, zeta = case["T"], 2.0 * np.pi * case["f_mode"], case["zeta"]
    sq = lambda t: square_wave(t, T)
    a0, an, bn = fourier_coefficients(sq, T, 51)
    _, an_c, bn_c = square_wave_coefficients(51)
    print("square wave, T = 1 s: b_n by trapezoid (4096 samples) vs 4/(pi n):")
    for n in (1, 2, 3, 5, 7, 15):
        q, _ = quad(lambda t: sq(np.array([t]))[0] * np.sin(2 * np.pi * n * t / T), 0.0, T, points=[T / 2], limit=200)
        sym = sympy_square_wave_coefficient(n, T) if n <= 3 else "(see n <= 3)"
        print(f"  n = {n:2d}: trapezoid {bn[n - 1]:+.6f}, closed form {bn_c[n - 1]:+.6f}, quad {2 * q / T:+.6f}, sympy {sym}")
    print(f"  a_n and a0 are zero to {max(abs(a0), np.max(np.abs(an))):.1e} (odd function)")
    print("Gibbs overshoot of the N-term partial sum:", {N: round(gibbs_overshoot(N, T), 5) for N in (1, 3, 9, 29, 99)},
          "-> (2/pi) Si(pi) = 1.17898, i.e. 8.95 % of the jump of 2")
    ms_c = parseval_power(0.0, an_c, bn_c)
    print(f"Parseval: mean square from 15 harmonics {ms_c:.5f}, from 99 harmonics {parseval_power(*square_wave_coefficients(99)):.5f}, exact 1 (8/pi^2 sum 1/n^2 odd = 1)")
    Ns = np.array([5, 11, 21, 41, 81, 161])
    err_sq = np.array([truncation_error_rms(sq, T, *square_wave_coefficients(N), N) for N in Ns])
    tri = lambda t: triangle_wave(t, T)
    err_tri = np.array([truncation_error_rms(tri, T, *triangle_wave_coefficients(N), N) for N in Ns])
    print(f"rms truncation error ~ N^-p: square p = {convergence_rate(err_sq, Ns):.2f} (jump: coefficients ~ 1/n),"
          f" triangle p = {convergence_rate(err_tri, Ns):.2f} (kink: coefficients ~ 1/n^2)")

    A0, An, Bn = response_coefficients(0.0, an_c, bn_c, T, wn, zeta)
    amp = harmonic_amplitudes(An, Bn) * wn ** 2
    k = dominant_harmonic(An, Bn)
    print(f"\nboom mode at {case['f_mode']} Hz, zeta = {zeta}, driven by the 1 Hz square wave:")
    print("  harmonic  forcing amp  response amp / static deflection")
    for n in (1, 3, 5, 7, 9):
        print(f"  {n:5d}     {bn_c[n - 1]:8.4f}       {amp[n - 1]:8.4f}")
    print(f"  dominant harmonic {k}: the mode answers the {k}th harmonic {amp[k - 1] / amp[0]:.1f} times more than the fundamental")
    print(f"  total response amplitude (peak of the harmonic sum) {np.max(np.abs(partial_sum(np.linspace(0, T, 4001), A0, An, Bn, T))) * wn ** 2:.3f} static deflections,"
          f" against {abs(complex_gain(2 * np.pi / T / wn, zeta)):.3f} for a unit sinusoid at 1 Hz")
    t = np.linspace(0.0, 40.0, 40001)
    sol = solve_ivp(forced_mode_rhs, (0.0, 40.0), [0.0, 0.0], t_eval=t, args=(lambda tt: float(square_wave(np.array([tt]), T)[0]), wn, zeta),
                    rtol=1e-9, atol=1e-12, max_step=2e-3)
    tail = t >= 35.0
    q_series = partial_sum(t[tail], A0, An, Bn, T)
    print(f"  solve_ivp after 35 s (transient e^{{-zeta wn t}} = {np.exp(-zeta * wn * 35):.1e}) vs 51-harmonic series: max |dq| = {np.max(np.abs(sol.y[0][tail] - q_series)) * wn ** 2:.2e} static deflections")

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))
    tt = np.linspace(-0.05 * T, 1.05 * T, 2001)
    axes[0].plot(tt, square_wave(tt, T), "k", lw=0.8, label="square wave")
    for N in (1, 3, 9, 29):
        a0_, an_, bn_ = square_wave_coefficients(N)
        axes[0].plot(tt, partial_sum(tt, a0_, an_, bn_, T), lw=1, label=f"N = {N}")
    axes[0].set_xlabel("t / T")
    axes[0].set_title("partial sums: the Gibbs ears do not shrink")
    axes[0].legend(fontsize=7)
    axes[0].grid(True, lw=0.4)
    axes[1].loglog(Ns, err_sq, "o-", label="square (jump)")
    axes[1].loglog(Ns, err_tri, "s-", label="triangle (kink)")
    axes[1].loglog(Ns, err_sq[0] * (Ns / Ns[0]) ** -0.5, "k--", lw=0.7, label=r"$N^{-1/2}$")
    axes[1].loglog(Ns, err_tri[0] * (Ns / Ns[0]) ** -1.5, "k:", lw=0.7, label=r"$N^{-3/2}$")
    axes[1].set_xlabel("harmonics kept N")
    axes[1].set_ylabel("rms truncation error")
    axes[1].set_title("smoothness sets the convergence rate")
    axes[1].legend(fontsize=7)
    axes[1].grid(True, which="both", lw=0.4)
    n = np.arange(1, 16)
    axes[2].stem(n, bn_c[:15], linefmt="0.6", markerfmt="ko", basefmt=" ", label="forcing b_n")
    axes[2].stem(n + 0.25, amp[:15], linefmt="r-", markerfmt="ro", basefmt=" ", label="response / static")
    axes[2].set_xlabel("harmonic n (of 1 Hz)")
    axes[2].set_title("boom at 5 Hz picks the 5th harmonic")
    axes[2].legend(fontsize=7)
    axes[2].grid(True, lw=0.4)
    fig.tight_layout()
    fig.savefig(out_dir / "fourier_series.png", dpi=150)
    print(f"saved {out_dir / 'fourier_series.png'}")


if __name__ == "__main__":
    main()
