"""P10.19 Unsteady Aerodynamics: starter code.

How much lift does a wing section really make while it oscillates, or in the moments after a
sudden change of angle of attack? Quasi-steady aerodynamics says "2 pi times the angle, now".
The wake disagrees: every change of circulation sheds vorticity into the wake, and that shed
vorticity induces a downwash that holds the lift back. This module computes the classical
incompressible thin-airfoil answers and checks them against each other:

* Theodorsen's function C(k) = H1(k) / (H1(k) + i H0(k)) with Hankel functions of the second
  kind (scipy.special.hankel2), and the lift and pitching moment of a section pitching and
  plunging harmonically at reduced frequency k = omega b / V;
* the Wagner function phi(s) (lift after a step change in angle of attack) in R.T. Jones's
  two-exponential approximation, and the Kussner function psi(s) (lift on entering a
  sharp-edged gust) in Sears's approximation;
* Duhamel's integral, which turns an indicial (step) response into the response to any motion;
* independent checks: Theodorsen tends to the quasi-steady 2 pi alpha_3/4 as k -> 0; the
  Laplace/Fourier transform of Jones's Wagner approximation gives back an approximate C(k) that
  is compared with the exact one; and a time-domain Duhamel simulation of harmonic pitch, after
  its transient dies out, reproduces Theodorsen's amplitude and phase.

Conventions (Theodorsen's): b is the semi-chord (c = 2b); the pivot (elastic axis) sits a
semi-chords aft of mid-chord (a = -1/2 is the quarter chord, a = 0 mid-chord); alpha positive
nose-up; plunge h positive DOWN (Theodorsen's sign, so that a downward plunge velocity raises the
effective angle of attack); lift positive up; moment about the pivot positive nose-up.
Reduced time s = V t / b (semi-chords travelled). Harmonic motion is alpha(t) = Re(alpha_hat
e^{i omega t}). Coefficients: cl = L'/(q c), cm = M'/(q c^2), q = rho V^2 / 2.

Units: SI; angles in radians inside the code.

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes unsteady_aero.png)
Test: pytest curriculum/phase-10/lesson-19-unsteady-aerodynamics
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.signal import fftconvolve
from scipy.special import hankel2

# Reference trainer (Phase 10 shared conventions): chord 1.2 m, cruise 55 m/s at sea level.
CHORD = 1.2
B_SEMI = CHORD / 2
V_CRUISE = 55.0
RHO = 1.225
F_TORSION = 5.0          # Hz: an illustrative wing torsion frequency used in the lesson (assumption)

# R.T. Jones's approximation of the Wagner function and Sears's of the Kussner function:
# phi(s) = 1 - sum A_i exp(-beta_i s)
WAGNER_A = (0.165, 0.335)
WAGNER_B = (0.0455, 0.300)
KUSSNER_A = (0.5, 0.5)
KUSSNER_B = (0.130, 1.000)


# ---------------------------------------------------------------- frequency domain

def theodorsen(k: np.ndarray | float) -> np.ndarray:
    """Theodorsen's function C(k) = H1(k) / (H1(k) + i H0(k)), H_n = Hankel function of the second
    kind. Complex; C(0) = 1 (the quasi-steady limit, returned exactly for k = 0), C -> 1/2 as k -> inf."""
    raise NotImplementedError


def theodorsen_from_indicial(k: np.ndarray | float, A=WAGNER_A, beta=WAGNER_B) -> np.ndarray:
    """C(k) implied by an exponential indicial function phi(s) = 1 - sum A_i exp(-beta_i s).

    For harmonic input alpha = e^{iks}, Duhamel gives cl = 2 pi alpha * (ik) * Phi(ik), where Phi(p)
    is the Laplace transform of phi. Phi(p) = 1/p - sum A_i/(p + beta_i), so
    C(k) = 1 - sum A_i ik/(ik + beta_i)."""
    raise NotImplementedError


def theodorsen_lift(k: float, alpha_hat: complex, h_hat_over_b: complex, a: float,
                    C: complex | None = None) -> tuple[complex, complex]:
    """Complex lift-coefficient amplitudes (circulatory, non-circulatory) for harmonic pitch
    alpha_hat about the pivot a and plunge h_hat/b (positive down) at reduced frequency k.

    L' = pi rho b^2 (hdd + V alphad - b a alphadd) + 2 pi rho V b C(k) (hd + V alpha + b(1/2 - a) alphad)
    divided by q c = rho V^2 b, with d/dt -> i omega:
    cl_nc = pi (-k^2 h/b + i k alpha + a k^2 alpha), cl_c = 2 pi C (i k h/b + alpha + (1/2 - a) i k alpha).
    C defaults to theodorsen(k); pass another value (e.g. Jones's) to compare."""
    raise NotImplementedError


def theodorsen_moment(k: float, alpha_hat: complex, h_hat_over_b: complex, a: float,
                      C: complex | None = None) -> complex:
    """Complex pitching-moment coefficient amplitude about the pivot (nose-up positive), cm = M'/(q c^2).

    M' = pi rho b^2 [b a hdd - V b (1/2 - a) alphad - b^2 (1/8 + a^2) alphadd]
         + 2 pi rho V b^2 (a + 1/2) C(k) [hd + V alpha + b (1/2 - a) alphad],
    divided by q c^2 = 2 rho V^2 b^2."""
    raise NotImplementedError


def quasi_steady_lift(k: float, alpha_hat: complex, h_hat_over_b: complex, a: float) -> complex:
    """Quasi-steady lift 2 pi alpha_3/4: the angle of attack seen at the three-quarter-chord point,
    alpha + hd/V + b (1/2 - a) alphad / V, with no wake lag and no apparent mass."""
    raise NotImplementedError


def aero_work_per_cycle(x_hat: complex, f_hat: complex) -> float:
    """Work done by a harmonic generalized force f on a harmonic displacement x over one cycle,
    oint f dx = pi Im(conj(x_hat) f_hat). Positive: the air feeds energy into the motion."""
    raise NotImplementedError


# ---------------------------------------------------------------- time domain

def wagner_jones(s: np.ndarray | float) -> np.ndarray:
    """Wagner function (lift build-up after a step in angle of attack, as a fraction of steady lift),
    R.T. Jones's approximation phi(s) = 1 - 0.165 e^{-0.0455 s} - 0.335 e^{-0.3 s}, s >= 0."""
    raise NotImplementedError


def kussner_sears(s: np.ndarray | float) -> np.ndarray:
    """Kussner function (lift build-up as the section penetrates a sharp-edged gust), Sears's
    approximation psi(s) = 1 - 0.5 e^{-0.13 s} - 0.5 e^{-s}, s >= 0 from the gust front at the
    leading edge. psi(0) = 0: no lift until the gust reaches the section."""
    raise NotImplementedError


def duhamel(s: np.ndarray, x: np.ndarray, indicial) -> np.ndarray:
    """Duhamel's integral y(s) = x(0) phi(s) + int_0^s phi(s - sigma) dx/dsigma dsigma on a uniform
    grid s_0 = 0 .. s_N, by the midpoint rule on each interval (second order in the step):
    y_n = x_0 phi(s_n) + sum_{j=1..n} (x_j - x_{j-1}) phi(s_n - s_{j-1/2})."""
    raise NotImplementedError


def time_domain_lift(s: np.ndarray, alpha: np.ndarray, h_over_b: np.ndarray, a: float,
                     indicial=wagner_jones) -> tuple[np.ndarray, np.ndarray]:
    """Lift coefficient (circulatory, non-circulatory) for an arbitrary motion alpha(s), h(s)/b
    (h positive down) on a uniform reduced-time grid. Circulatory: 2 pi times the Duhamel integral
    of alpha_3/4 = alpha + h'(s)/b + (1/2 - a) alpha'(s) with the Wagner function. Apparent mass:
    pi (h''/b + alpha' - a alpha''). Derivatives by second-order finite differences."""
    raise NotImplementedError


def fit_harmonic(s: np.ndarray, y: np.ndarray, k: float) -> complex:
    """Least-squares complex amplitude Y of y(s) ~ c0 + Re(Y e^{iks}) (a constant is fitted too)."""
    M = np.column_stack([np.ones_like(s), np.cos(k * s), np.sin(k * s)])
    c, *_ = np.linalg.lstsq(M, y, rcond=None)
    return complex(c[1], -c[2])


def simulate_harmonic_pitch(k: float, alpha_amp: float, a: float, n_cycles: int = 12,
                            pts_per_cycle: int = 400, indicial=wagner_jones) -> dict:
    """Start a pitch oscillation alpha = alpha_amp sin(k s) from rest, integrate Duhamel's integral,
    and fit the last two cycles (after the wake transient has decayed). Returns the complex
    amplitudes of the circulatory, non-circulatory and total lift relative to alpha_hat, where
    alpha = Re(alpha_hat e^{iks}) with alpha_hat = -i alpha_amp."""
    raise NotImplementedError


def rise_distance(level: float, indicial=wagner_jones) -> float:
    """Reduced time s (semi-chords) at which an indicial function first reaches `level`."""
    raise NotImplementedError


def gust_response(s: np.ndarray, w_over_V: np.ndarray, indicial=kussner_sears) -> np.ndarray:
    """Lift coefficient on a section flying through a vertical gust w(s)/V (s = position of the gust
    front relative to the leading edge, in semi-chords): cl = 2 pi Duhamel(w/V, Kussner)."""
    raise NotImplementedError


def lift_loop(k: float, alpha_mean: float, alpha_amp: float, a: float, n: int = 361) -> tuple[np.ndarray, np.ndarray]:
    """One cycle of (alpha, cl) for alpha = alpha_mean + alpha_amp sin(theta) from Theodorsen:
    the hysteresis ellipse. Returns alpha [rad] and cl."""
    raise NotImplementedError


def loop_area(x: np.ndarray, y: np.ndarray) -> float:
    """oint y dx around a closed sampled loop (trapezoid rule)."""
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def pitch_flutter_boundary(a: float, k_lo: float = 1e-3, k_hi: float = 1.0) -> float | None:
    """Reduced frequency below which a single-degree-of-freedom pitch oscillation about the pivot a
    is fed energy by the air (aerodynamic work per cycle > 0), or None if it is damped throughout."""
    raise NotImplementedError


# ---------------------------------------------------------------- given helpers

def observed_order(errors: np.ndarray, ratio: float = 2.0) -> np.ndarray:
    """Observed order p = log(e_i / e_{i+1}) / log(ratio) for successive refinements."""
    e = np.asarray(errors, dtype=float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


def reduced_frequency(f_hz: float, b: float, V: float) -> float:
    """k = omega b / V."""
    return 2 * np.pi * f_hz * b / V


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and write unsteady_aero.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    print("Theodorsen's function C(k) = F + iG")
    print("     k        F         G        |C|    phase [deg]   Jones-approx |dC|")
    for k in [0.0, 0.01, 0.05, 0.1, 0.2, 0.3427, 0.5, 1.0, 2.0, 10.0]:
        C = theodorsen(k)
        Cj = theodorsen_from_indicial(k)
        print(f"  {k:6.4f}  {C.real:8.4f}  {C.imag:8.4f}  {abs(C):7.4f}  {np.degrees(np.angle(C)):9.3f}      {abs(Cj - C):.4f}")
    ks = np.linspace(1e-4, 5.0, 20001)
    err = np.abs(theodorsen_from_indicial(ks) - theodorsen(ks))
    print(f"max |C_Jones - C| over 0 < k <= 5: {err.max():.4f} at k = {ks[err.argmax()]:.3f}"
          f" ({100 * (err / np.abs(theodorsen(ks))).max():.2f} % of |C| at most)")
    for k in [1e-2, 1e-3, 1e-4]:
        C = theodorsen(k)
        print(f"  k = {k:.0e}: |1 - C| = {abs(1 - C):.4f}")

    # The engineering question: the trainer's wing section pitching about the quarter chord at 5 Hz
    k = reduced_frequency(F_TORSION, B_SEMI, V_CRUISE)
    a = -0.5
    c, nc = theodorsen_lift(k, 1.0, 0.0, a)
    qs = quasi_steady_lift(k, 1.0, 0.0, a)
    tot = c + nc
    print(f"\ntrainer section, c = {CHORD} m at {V_CRUISE} m/s, pitch about c/4 at {F_TORSION} Hz: k = {k:.4f}")
    print(f"  |C| = {abs(theodorsen(k)):.4f}, phase {np.degrees(np.angle(theodorsen(k))):.2f} deg")
    print(f"  cl/alpha_hat: steady 2 pi = {2*np.pi:.4f}; quasi-steady |2 pi alpha_3/4| = {abs(qs):.4f} "
          f"(phase {np.degrees(np.angle(qs)):+.2f} deg)")
    print(f"  circulatory: |{abs(c):.4f}| phase {np.degrees(np.angle(c)):+.2f} deg; apparent mass |{abs(nc):.4f}| "
          f"phase {np.degrees(np.angle(nc)):+.2f} deg")
    print(f"  total Theodorsen: |{abs(tot):.4f}| = {abs(tot)/(2*np.pi):.3f} x steady, phase {np.degrees(np.angle(tot)):+.2f} deg"
          f" (time {np.degrees(np.angle(tot))/360/F_TORSION*1e3:+.2f} ms)")
    cl_ac = theodorsen_lift(k, 1.0, 0.0, a)[0] / (1 + 1j * k)
    print(f"  wake factor on the circulatory lift: |C| = {abs(cl_ac)/(2*np.pi):.4f}, lag {-np.degrees(np.angle(cl_ac)):.2f} deg")

    # Step response
    for lev in (0.75, 0.9, 0.95):
        sr = rise_distance(lev)
        print(f"Wagner (Jones): {lev*100:.0f} % of steady lift at s = {sr:.2f} semi-chords"
              f" = {sr*B_SEMI:.2f} m = {sr*B_SEMI/V_CRUISE*1e3:.1f} ms at {V_CRUISE} m/s")
    for lev in (0.5, 0.9):
        sr = rise_distance(lev, kussner_sears)
        print(f"Kussner (Sears): {lev*100:.0f} % at s = {sr:.2f} semi-chords = {sr*B_SEMI:.2f} m")

    # Independent check: time-domain Duhamel vs Theodorsen
    print("\nharmonic pitch about c/4, Duhamel with Jones's Wagner function vs frequency domain")
    print("     k    |cl|/alpha TD   with C_Jones   with exact C    phase TD   phase exact   |TD-Jones|/|.|  |TD-exact|/|.|")
    rows = []
    for kk in [0.05, 0.1, 0.2, 0.3427, 0.5, 1.0]:
        sim = simulate_harmonic_pitch(kk, np.radians(2.0), a)
        cj = sum(theodorsen_lift(kk, 1.0, 0.0, a, C=theodorsen_from_indicial(kk)))
        ce = sum(theodorsen_lift(kk, 1.0, 0.0, a))
        td = sim["total"]
        rows.append((kk, td, cj, ce))
        print(f"  {kk:6.4f}  {abs(td):10.4f}  {abs(cj):12.4f}  {abs(ce):12.4f}  {np.degrees(np.angle(td)):9.3f}  "
              f"{np.degrees(np.angle(ce)):10.3f}   {abs(td-cj)/abs(cj):.2e}      {abs(td-ce)/abs(ce):.4f}")

    # Duhamel convergence against the closed-form convolution for alpha = sin(k s)
    kk = 0.3
    errs = []
    for n in [100, 200, 400, 800]:
        s = np.linspace(0, 40, n + 1)
        y = duhamel(s, np.sin(kk * s), wagner_jones)
        errs.append(np.max(np.abs(y - _exact_jones_sine(s, kk))))
    print(f"\nDuhamel midpoint rule vs exact convolution (alpha = sin 0.3 s, 0 <= s <= 40): errors "
          + ", ".join(f"{e:.2e}" for e in errs) + f"; observed orders {np.round(observed_order(errs), 3)}")

    # Lift loops and energy
    print("\nlift loops, pitch about c/4, alpha = 4 +/- 2 deg")
    for kk in [0.05, 0.1, 0.3427, 1.0]:
        al, cl = lift_loop(kk, np.radians(4), np.radians(2), a)
        print(f"  k = {kk:6.4f}: cl from {cl.min():.4f} to {cl.max():.4f}, loop area oint cl dalpha = {loop_area(al, cl):+.5f}")
    print("work per cycle on pure plunge (h_hat/b = 1), oint cl d(h/b) with h down:")
    for kk in [0.05, 0.3, 1.0]:
        c, nc = theodorsen_lift(kk, 0.0, 1.0, a)
        print(f"  k = {kk}: {aero_work_per_cycle(1.0, -(c + nc)):+.4f} (negative = damped)")

    print("single-degree-of-freedom pitch: aerodynamic work per cycle (alpha_hat = 1 rad) about pivot a")
    print("     a   x_pivot/c     k=0.02     k=0.05     k=0.1      k=0.3    k below which energy is fed in")
    for aa in [-0.5, -0.6, -0.7, -0.8, -1.0, 0.0]:
        w = [aero_work_per_cycle(1.0, theodorsen_moment(kk, 1.0, 0.0, aa)) for kk in (0.02, 0.05, 0.1, 0.3)]
        kb = pitch_flutter_boundary(aa)
        print(f"  {aa:5.2f}   {(1 + aa)/2:6.3f}  " + "  ".join(f"{x:+8.4f}" for x in w)
              + f"    {('%.4f' % kb) if kb else 'none (damped)'}")

    # Gusts on the trainer
    w_g = 5.0
    s = np.linspace(0, 60, 6001)
    cl_sharp = gust_response(s, np.full_like(s, w_g / V_CRUISE))
    print(f"\nsharp-edged {w_g} m/s gust at {V_CRUISE} m/s: quasi-steady dcl = {2*np.pi*w_g/V_CRUISE:.4f};"
          f" Kussner: {cl_sharp[np.searchsorted(s, 2.0)]:.4f} after 1 chord, {cl_sharp[np.searchsorted(s, 10.0)]:.4f} after 5 chords")
    print("1-cos gust of gradient length H (gust length 2H), peak dcl / quasi-steady peak:")
    for H in [2.0, 5.0, 12.5, 25.0, 50.0]:          # semi-chords: 1, 2.5, 6.25, 12.5, 25 chords
        ss = np.linspace(0, 2 * H + 30, int((2 * H + 30) * 100) + 1)
        wv = np.where(ss <= 2 * H, 0.5 * (1 - np.cos(np.pi * ss / H)), 0.0) * w_g / V_CRUISE
        r = gust_response(ss, wv).max() / (2 * np.pi * w_g / V_CRUISE)
        print(f"  H = {H/2:5.2f} chords ({H*B_SEMI:5.2f} m): {r:.4f}")

    # Figure
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
    kp = np.concatenate([[0], np.logspace(-3, 1.3, 300)])
    C = theodorsen(kp); Cj = theodorsen_from_indicial(kp)
    ax[0, 0].plot(C.real, C.imag, label="exact C(k) (Hankel functions)")
    ax[0, 0].plot(Cj.real, Cj.imag, "--", label="from Jones's Wagner approximation")
    for kk in [0.05, 0.1, 0.2, 0.5, 1.0]:
        z = theodorsen(kk); ax[0, 0].plot(z.real, z.imag, "k.")
        ax[0, 0].annotate(f"k={kk}", (z.real, z.imag), fontsize=7, xytext=(3, -9), textcoords="offset points")
    ax[0, 0].set_xlabel("F = Re C"); ax[0, 0].set_ylabel("G = Im C"); ax[0, 0].set_title("Theodorsen's function")
    ax[0, 0].legend(fontsize=8); ax[0, 0].grid(alpha=0.3)
    sp = np.linspace(0, 40, 400)
    ax[0, 1].plot(sp, wagner_jones(sp), label="Wagner phi(s) (Jones)")
    ax[0, 1].plot(sp, kussner_sears(sp), label="Kussner psi(s) (Sears)")
    ax[0, 1].axhline(1, color="k", lw=0.6)
    ax[0, 1].set_xlabel("s = V t / b (semi-chords)"); ax[0, 1].set_ylabel("lift / steady lift")
    ax[0, 1].set_title("indicial lift"); ax[0, 1].legend(fontsize=8); ax[0, 1].grid(alpha=0.3)
    for kk in [0.05, 0.1, 0.3427, 1.0]:
        al, cl = lift_loop(kk, np.radians(4), np.radians(2), a)
        ax[1, 0].plot(np.degrees(al), cl, label=f"k = {kk}")
    al = np.radians(np.linspace(2, 6, 3)); ax[1, 0].plot(np.degrees(al), 2 * np.pi * al, "k--", label="steady 2 pi alpha")
    ax[1, 0].set_xlabel("alpha [deg]"); ax[1, 0].set_ylabel("cl"); ax[1, 0].set_title("lift loops, pitch about c/4")
    ax[1, 0].legend(fontsize=8); ax[1, 0].grid(alpha=0.3)
    sim = simulate_harmonic_pitch(k, np.radians(2.0), a, n_cycles=6)
    ax[1, 1].plot(sim["s"], sim["cl"], label="Duhamel (time domain)")
    ce = sum(theodorsen_lift(k, 1.0, 0.0, a))
    ax[1, 1].plot(sim["s"], np.real(ce * (-1j) * np.radians(2.0) * np.exp(1j * k * sim["s"])), "--",
                  label="Theodorsen (periodic)")
    ax[1, 1].set_xlabel("s (semi-chords)"); ax[1, 1].set_ylabel("cl"); ax[1, 1].set_title(f"pitch from rest, k = {k:.3f}")
    ax[1, 1].legend(fontsize=8); ax[1, 1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "unsteady_aero.png", dpi=120)
    plt.close(fig)


def _exact_jones_sine(s: np.ndarray, k: float) -> np.ndarray:
    """Closed-form Duhamel response of Jones's Wagner function to x(s) = sin(k s), x(0) = 0."""
    y = np.sin(k * s)
    for Ai, bi in zip(WAGNER_A, WAGNER_B):
        y = y - Ai * k * (bi * np.cos(k * s) + k * np.sin(k * s) - bi * np.exp(-bi * s)) / (bi ** 2 + k ** 2)
    return y


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
