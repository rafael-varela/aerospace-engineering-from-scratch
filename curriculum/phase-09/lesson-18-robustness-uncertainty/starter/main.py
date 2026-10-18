"""P09.18 Robustness & Uncertainty: starter code.

How much can the real plant differ from the model before the loop goes unstable? This module answers
with the sensitivity functions of a feedback loop and a first taste of H-infinity analysis:

    * the four closed-loop maps S = 1/(1 + L), T = L/(1 + L), KS and GS, as polynomial transfer functions;
    * their peaks Ms = ||S||inf and Mt = ||T||inf, and the disk-margin guarantees that Ms implies:
      gain margin in [Ms/(Ms + 1), Ms/(Ms - 1)] and phase margin >= 2 arcsin(1/(2 Ms));
    * the H-infinity norm of a stable state-space system two independent ways: a frequency sweep with local
      refinement, and bisection on the imaginary-axis eigenvalues of the Hamiltonian matrix;
    * multiplicative uncertainty G = G0 (1 + W Delta), |Delta| <= 1, and the small-gain robust-stability
      test ||W T||inf < 1, including the destabilising perturbation that exists when the test fails;
    * the Bode sensitivity integral  int_0^inf ln|S(jw)| dw = pi sum Re(p_i)  over open-loop RHP poles.

The example is a control-surface servo (DC-motor position loop) with a lead compensator, whose real
actuator has an unmodelled gear-compliance mode and a small transport delay.

Conventions:
    * Transfer functions are (num, den) coefficient arrays in descending powers of s.
    * Loop: u = C(s) (r - y), y = G(s) u; L = C G; negative feedback.
    * Frequencies w in rad/s; margins as ratios (not dB) unless a name says _db; phases in degrees.

Complete every function marked NotImplementedError.

Run:  python main.py          (prints the servo robustness study, writes robustness.png)
Test: pytest curriculum/phase-09/lesson-18-robustness-uncertainty
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


# ---------------------------------------------------------------- given helpers

def freqresp(num: np.ndarray, den: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Return G(jw) = num(jw)/den(jw)."""
    s = 1j * np.asarray(w, dtype=float)
    return np.polyval(num, s) / np.polyval(den, s)


def pade(tau: float, order: int = 3) -> tuple[np.ndarray, np.ndarray]:
    """Return (num, den) of the [order/order] Pade approximant of the delay e^{-s tau}."""
    from math import factorial
    c = [factorial(2 * order - k) * factorial(order) / (factorial(2 * order) * factorial(k) * factorial(order - k))
         for k in range(order + 1)]
    den = np.array([c[k] * tau ** k for k in range(order, -1, -1)])
    num = np.array([c[k] * (-tau) ** k for k in range(order, -1, -1)])
    return num / den[0], den / den[0]


def servo_plant() -> tuple[np.ndarray, np.ndarray]:
    """Return the nominal servo model G0(s) = Km / (s (tm s + 1)), Km = 30 rad/(V s), tm = 0.04 s
    (a representative small control-surface actuator: output shaft angle per volt)."""
    return np.array([30.0]), np.array([0.04, 1.0, 0.0])


def true_servo_plant(pade_order: int = 3) -> tuple[np.ndarray, np.ndarray]:
    """Return the 'true' servo: G0 times an unmodelled gear-compliance mode wa^2/(s^2 + 2 za wa s + wa^2),
    wa = 150 rad/s, za = 0.3, times a 5 ms transport delay (Pade approximant). Representative numbers."""
    n0, d0 = servo_plant()
    wa, za = 150.0, 0.3
    npd, dpd = pade(0.005, pade_order)
    num = np.polymul(n0 * wa ** 2, npd)
    den = np.polymul(np.polymul(d0, [1.0, 2 * za * wa, wa ** 2]), dpd)
    return num, den


def true_servo_relative_error(w: np.ndarray) -> np.ndarray:
    """Return Delta_true(jw) = G_true(jw)/G0(jw) - 1 with the exact delay e^{-jw 0.005}."""
    s = 1j * np.asarray(w, dtype=float)
    wa, za = 150.0, 0.3
    return wa ** 2 / (s ** 2 + 2 * za * wa * s + wa ** 2) * np.exp(-0.005 * s) - 1.0


def uncertainty_weight(r0: float = 0.05, r_inf: float = 4.0, wb: float = 40.0) -> tuple[np.ndarray, np.ndarray]:
    """Return W(s) = (s/wb + r0) / (s/(wb r_inf) + 1): relative uncertainty r0 at low frequency, rising
    through 100 % near wb and levelling at r_inf at high frequency."""
    return np.array([1.0 / wb, r0]), np.array([1.0 / (wb * r_inf), 1.0])


def lead(Kc: float, z: float, p: float) -> tuple[np.ndarray, np.ndarray]:
    """Return C(s) = Kc (s/z + 1)/(s/p + 1)."""
    return np.array([Kc / z, Kc]), np.array([1.0 / p, 1.0])


# ---------------------------------------------------------------- closed-loop maps

def closed_loop_tfs(nC: np.ndarray, dC: np.ndarray, nG: np.ndarray, dG: np.ndarray) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Return {"L", "S", "T", "KS", "GS"} as (num, den) polynomial pairs.

    With L = nC nG / (dC dG) and characteristic polynomial chi = dC dG + nC nG:
    S = dC dG / chi, T = nC nG / chi, KS = nC dG / chi, GS = nG dC / chi.
    """
    raise NotImplementedError


def is_closed_loop_stable(nC: np.ndarray, dC: np.ndarray, nG: np.ndarray, dG: np.ndarray) -> bool:
    """Return True if every root of dC dG + nC nG has negative real part (internal stability, assuming no
    unstable pole-zero cancellation between C and G)."""
    raise NotImplementedError


def tf_to_ss(num: np.ndarray, den: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (A, B, C, D) of a proper SISO transfer function in controllable canonical form.

    Normalise den to be monic of degree n, pad num to length n + 1; D = num[0]; with the strictly proper
    remainder r(s) = num(s) - D den(s) = r_1 s^{n-1} + ... + r_n:
    A = [[-a_1, ..., -a_n], [1, 0, ...], ..., [0, ..., 1, 0]], B = e_1, C = [r_1, ..., r_n].
    Raise ``ValueError`` if improper.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- H-infinity norm two ways

def hinf_norm_sweep(A: np.ndarray, B: np.ndarray, C: np.ndarray, D: np.ndarray,
                    w: np.ndarray | None = None) -> tuple[float, float]:
    """Return (peak gain, frequency): max over w of sigma_max(C (jwI - A)^{-1} B + D) on a log grid,
    refined by golden-section search in log w around the best grid point. Includes w = 0 and the
    high-frequency limit sigma_max(D)."""
    raise NotImplementedError


def hamiltonian_has_imaginary_eigs(A: np.ndarray, B: np.ndarray, C: np.ndarray, D: np.ndarray,
                                   gamma: float, tol: float = 1e-8) -> bool:
    """Return True if the Hamiltonian H(gamma) has an eigenvalue on the imaginary axis, i.e. if
    ||G||inf >= gamma for a stable G with sigma_max(D) < gamma.

    With R = gamma^2 I - D^T D, Rt = gamma^2 I - D D^T and F = A + B R^{-1} D^T C:
    H = [[F, B R^{-1} B^T], [-gamma^2 C^T Rt^{-1} C, -F^T]]
    (gamma^2 Rt^{-1} = I + D R^{-1} D^T; for D = 0 this is [[A, B B^T / gamma^2], [-C^T C, -A^T]]).
    """
    raise NotImplementedError


def hinf_norm_bisection(A: np.ndarray, B: np.ndarray, C: np.ndarray, D: np.ndarray,
                        rtol: float = 1e-8) -> float:
    """Return ||G||inf of a stable system by bisection on gamma using ``hamiltonian_has_imaginary_eigs``.

    Lower bound: max(sigma_max(D), sigma_max(G(0))); upper bound found by doubling. Raise ``ValueError``
    if A is not Hurwitz.
    """
    raise NotImplementedError


def tf_hinf_norm(num: np.ndarray, den: np.ndarray, method: str = "bisection") -> float:
    """Return ||num/den||inf of a stable proper transfer function, by ``method`` "bisection" or "sweep"."""
    A, B, C, D = tf_to_ss(num, den)
    if method == "bisection":
        return hinf_norm_bisection(A, B, C, D)
    return hinf_norm_sweep(A, B, C, D)[0]


# ---------------------------------------------------------------- margins

def disk_margin_bounds(Ms: float) -> tuple[float, float, float]:
    """Return (gm_low, gm_high, pm_deg) guaranteed by a sensitivity peak Ms > 1: the loop stays stable for any
    gain in (Ms/(Ms + 1), Ms/(Ms - 1)) and for any phase lag up to 2 arcsin(1/(2 Ms)).

    Reason: min_w |1 + L(jw)| = 1/Ms, so the Nyquist curve keeps out of the disk of radius 1/Ms about -1.
    """
    raise NotImplementedError


def margins(w: np.ndarray, L: np.ndarray) -> dict[str, float]:
    """Return measured margins of a loop sampled on an increasing grid.

    "pm" [deg] and "wgc": at the first fall of |L| through 1 (180 + angle L there, angle unwrapped);
    "gm_high": the smallest 1/|L| > 1 at a phase crossover (Im L changes sign with Re L < 0), inf if none;
    "gm_low": the largest 1/|L| < 1 at a phase crossover (gain-reduction margin), 0 if none.
    Crossings are located by linear interpolation between grid points; a loop whose lowest-frequency sample
    lies on the negative real axis (within 0.01 rad) counts as crossing there (L(0) real and negative, as for
    an unstable plant without an integrator).
    """
    raise NotImplementedError


# ---------------------------------------------------------------- uncertainty

def weight_covers(w: np.ndarray, delta: np.ndarray, Wnum: np.ndarray, Wden: np.ndarray) -> float:
    """Return max_w |delta(jw)| / |W(jw)|: the weight covers the relative error when this is <= 1."""
    return float(np.max(np.abs(delta) / np.abs(freqresp(Wnum, Wden, w))))


def robust_stability_index(nC: np.ndarray, dC: np.ndarray, nG: np.ndarray, dG: np.ndarray,
                           Wnum: np.ndarray, Wden: np.ndarray) -> float:
    """Return ||W T||inf (Hamiltonian bisection): the loop is stable for every G0 (1 + W Delta) with
    ||Delta||inf <= 1 if and only if this is < 1 (the small-gain theorem)."""
    raise NotImplementedError


def destabilising_perturbation(w_star: float, target: complex) -> tuple[np.ndarray, np.ndarray]:
    """Return a stable real-rational Delta(s) = c (a - s)/(a + s) (or a constant) with Delta(j w_star) = target.

    Used to show necessity in the small-gain theorem: with target = -1/(W T)(j w_star), the perturbed loop
    has closed-loop poles at +-j w_star. |Delta(jw)| = |target| at every w (all-pass), so ||Delta||inf =
    |target|. Raise ``ValueError`` for w_star <= 0 unless target is real.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- Bode sensitivity integral

def bode_sensitivity_integral(nS: np.ndarray, dS: np.ndarray, w_max: float = 1e6, decades: float = 12.0,
                              n: int = 240001) -> float:
    """Return int_0^inf ln|S(jw)| dw numerically.

    Trapezoidal rule on a grid uniform in log w from w_max 10^-decades to w_max (the integrand written as
    ln|S| w d(ln w)); the piece below the grid is neglected (it is O(w_lo ln w_lo)); the tail beyond w_max is
    added as w_max ln|S(j w_max)|, exact when ln|S| ~ c/w^2 there (a loop of relative degree two).
    """
    raise NotImplementedError


def unstable_airframe() -> tuple[np.ndarray, np.ndarray]:
    """Return G(s) = 20 / ((s - 2)(0.05 s + 1)): a statically unstable pitch mode (representative divergence
    rate p = 2 rad/s, time to double ln 2 / 2 = 0.35 s) with a fast actuator."""
    return np.array([20.0]), np.polymul([1.0, -2.0], [0.05, 1.0])


# ---------------------------------------------------------------- study

def design_table(designs: list[tuple[float, float, float]] | None = None) -> list[dict[str, float]]:
    """For each lead (Kc, z, p) on the nominal servo: Ms, Mt, the disk-margin bounds, ||W T||inf, the margins
    measured on the true servo (exact delay, dense grid), and whether the true servo loop is stable."""
    if designs is None:
        designs = [(0.5, 20.0, 80.0), (0.8, 20.0, 80.0), (1.2, 20.0, 80.0), (1.6, 25.0, 100.0), (2.5, 25.0, 150.0), (3.2, 25.0, 200.0)]
    nG, dG = servo_plant()
    nGt, dGt = true_servo_plant(pade_order=4)
    Wn, Wd = uncertainty_weight()
    w = np.logspace(-1, 4, 60000)
    rows = []
    for Kc, z, p in designs:
        nC, dC = lead(Kc, z, p)
        tfs = closed_loop_tfs(nC, dC, nG, dG)
        Ms = tf_hinf_norm(*tfs["S"])
        Mt = tf_hinf_norm(*tfs["T"])
        gl, gh, pmb = disk_margin_bounds(Ms)
        wt = robust_stability_index(nC, dC, nG, dG, Wn, Wd)
        Lnom = freqresp(nC, dC, w) * freqresp(nG, dG, w)
        Ltrue = Lnom * (1.0 + true_servo_relative_error(w))
        m_nom = margins(w, Lnom)
        m_true = margins(w, Ltrue)
        Ms_true = float(np.max(np.abs(1.0 / (1.0 + Ltrue))))
        rows.append({"Kc": Kc, "z": z, "p": p, "Ms": Ms, "Mt": Mt, "gm_low_bound": gl, "gm_high_bound": gh,
                     "pm_bound": pmb, "pm_nom": m_nom["pm"], "wgc": m_nom["wgc"], "WT": wt,
                     "pm_true": m_true["pm"], "gm_true": m_true["gm_high"], "Ms_true": Ms_true,
                     "true_stable": is_closed_loop_stable(nC, dC, nGt, dGt)})
    return rows


def main(out_dir: str | Path = ".") -> None:
    """Print the servo robustness study and save robustness.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    nG, dG = servo_plant()
    Wn, Wd = uncertainty_weight()
    w = np.logspace(-1, 4, 60000)
    delta = true_servo_relative_error(w)
    print(f"servo G0 = 30/(s(0.04 s + 1)); unmodelled mode 150 rad/s, zeta 0.3, delay 5 ms: peak |Delta| = {np.abs(delta).max():.3f}"
          f" at {w[np.argmax(np.abs(delta))]:.1f} rad/s")
    print(f"weight W = (s/40 + 0.05)/(s/160 + 1): max |Delta|/|W| = {weight_covers(w, delta, Wn, Wd):.3f} (<= 1: covered)")

    rows = design_table()
    print("lead designs: Kc, z, p | wgc, PM nominal | Ms, Mt | disk bounds GM in (lo, hi), PM >= | ||WT|| | true: stable, Ms, PM, GM, bounds from its Ms")
    for r in rows:
        line = (f"  {r['Kc']:.1f} {r['z']:4.0f} {r['p']:4.0f} | {r['wgc']:6.2f} {r['pm_nom']:5.1f} | {r['Ms']:.3f} {r['Mt']:.3f} |"
                f" ({r['gm_low_bound']:.3f}, {r['gm_high_bound']:.2f}) {r['pm_bound']:5.1f} | {r['WT']:.3f} | {r['true_stable']}")
        if r["true_stable"]:
            gl, gh, pb = disk_margin_bounds(r["Ms_true"])
            line += f" {r['Ms_true']:.3f} {r['pm_true']:5.1f} {r['gm_true']:5.2f} | GM >= {gh:.2f}, PM >= {pb:.1f}"
        else:
            line += f" (PM on the true loop {r['pm_true']:.1f})"
        print(line)

    # H-infinity of S for the chosen design, two ways
    nC, dC = lead(0.8, 20.0, 80.0)
    tfs = closed_loop_tfs(nC, dC, nG, dG)
    A, B, C, D = tf_to_ss(*tfs["S"])
    sw, wpk = hinf_norm_sweep(A, B, C, D)
    bi = hinf_norm_bisection(A, B, C, D)
    print(f"chosen design Kc = 0.8: Ms by sweep {sw:.8f} at {wpk:.2f} rad/s, by Hamiltonian bisection {bi:.8f}")
    ks = tf_hinf_norm(*tfs["KS"])
    print(f"   ||KS||inf = {ks:.3f} V/rad (peak control effort per unit reference), ||T||inf = {tf_hinf_norm(*tfs['T']):.4f}")

    # necessity: for the design with ||WT|| > 1, the destabilising Delta
    nC2, dC2 = lead(2.5, 25.0, 150.0)
    nT, dT = closed_loop_tfs(nC2, dC2, nG, dG)["T"]
    nWT, dWT = np.polymul(Wn, nT), np.polymul(Wd, dT)
    A2, B2, C2, D2 = tf_to_ss(nWT, dWT)
    peak, wst = hinf_norm_sweep(A2, B2, C2, D2)
    target = -1.0 / freqresp(nWT, dWT, [wst])[0]
    dn, dd = destabilising_perturbation(wst, target)
    nGp = np.polymul(nG, np.polyadd(np.polymul(Wd, dd), np.polymul(Wn, dn)))
    dGp = np.polymul(dG, np.polymul(Wd, dd))
    chi = np.polyadd(np.polymul(dC2, dGp), np.polymul(nC2, nGp))
    r = np.roots(chi)
    k = np.argmin(np.abs(r - 1j * wst))
    print(f"Kc = 2.5: ||WT|| = {peak:.3f} at {wst:.1f} rad/s; an all-pass Delta of size {abs(target):.3f} puts a closed-loop pole at {r[k]:.4f}")

    # Bode integral
    s_int = bode_sensitivity_integral(*tfs["S"])
    nU, dU = unstable_airframe()
    nCu, dCu = lead(0.6, 1.0, 40.0)
    print(f"unstable airframe with C = 0.6 (s + 1)/(s/40 + 1): closed loop stable {is_closed_loop_stable(nCu, dCu, nU, dU)}")
    tfu = closed_loop_tfs(nCu, dCu, nU, dU)
    u_int = bode_sensitivity_integral(*tfu["S"])
    Msu = tf_hinf_norm(*tfu["S"])
    print(f"Bode integral: servo {s_int:.2e} (theory 0); unstable airframe {u_int:.4f} (theory pi p = {np.pi * 2.0:.4f}); its Ms = {Msu:.3f}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for r_, col in zip(rows[:5], ("C0", "C1", "C2", "C3", "C4")):
        nC_, dC_ = lead(r_["Kc"], r_["z"], r_["p"])
        t_ = closed_loop_tfs(nC_, dC_, nG, dG)
        axes[0].loglog(w, np.abs(freqresp(*t_["T"], w)), col, lw=0.9, label=f"|T|, Kc = {r_['Kc']}")
    axes[0].loglog(w, 1.0 / np.abs(freqresp(Wn, Wd, w)), "k", lw=1.4, label="1/|W|")
    axes[0].loglog(w, np.abs(delta), "k:", lw=0.9, label="|Delta_true|")
    axes[0].set_ylim(1e-3, 10); axes[0].set_xlabel("w [rad/s]"); axes[0].legend(fontsize=7); axes[0].grid(True, which="both", lw=0.3)
    axes[0].set_title("small gain: |T| below 1/|W|", fontsize=9)
    wl = np.logspace(-2, 4, 4000)
    axes[1].semilogx(wl, np.log(np.abs(freqresp(*tfs["S"], wl))), label="servo, ln|S|")
    axes[1].semilogx(wl, np.log(np.abs(freqresp(*tfu["S"], wl))), label="unstable airframe, ln|S|")
    axes[1].axhline(0.0, color="k", lw=0.6)
    axes[1].set_xlabel("w [rad/s]"); axes[1].legend(fontsize=8); axes[1].grid(True, which="both", lw=0.3)
    axes[1].set_title("waterbed: area of ln|S| is fixed", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "robustness.png", dpi=150)
    plt.close(fig)
    print(f"saved {out_dir / 'robustness.png'}")


if __name__ == "__main__":
    main()
