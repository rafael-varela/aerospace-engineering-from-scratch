"""P09.17 Digital Control: starter code.

A flight computer samples the aircraft, computes, and holds its output constant until the next sample.
This module builds the three classic discretisations from scratch and uses them on a digital pitch
attitude loop:

    * zero-order hold (ZOH), exact for a held input, from the augmented matrix exponential
      exp([[A, B], [0, 0]] T) = [[Ad, Bd], [0, I]];
    * Tustin (bilinear) s = (2/T)(z - 1)/(z + 1), in state space and on polynomial transfer functions,
      with optional prewarping;
    * matched pole-zero, z = e^{sT} applied to poles and zeros;

plus the z-transform tools (difference-equation simulation, impulse response by long division,
frequency response on the unit circle, z-plane stability) and the ZOH discretisation with a
computational delay tau (0 <= tau <= T) that a real flight computer adds between sampling and output.
Tests check everything against ``scipy.signal.cont2discrete`` and ``scipy.linalg.expm``, and check step
invariance (ZOH reproduces the continuous step response exactly at the sample instants).

Conventions:
    * Polynomials are coefficient arrays in descending powers (NumPy ``polyval`` order), in s or in z.
    * Pitch state x = (alpha, q, theta, delta_e): angle of attack [rad], pitch rate [rad/s], pitch angle
      [rad], elevator [rad]. Input u = elevator command [rad]; positive elevator is nose-up (M_delta > 0).
    * The controller samples x at t_k = kT and its output u_k reaches the servo at t_k + tau.

Complete every function marked NotImplementedError.

Run:  python main.py          (prints the pitch-loop study, writes digital_pitch.png)
Test: pytest curriculum/phase-09/lesson-17-digital-control
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


# ---------------------------------------------------------------- given helpers

def expm(M: np.ndarray, n_terms: int = 20) -> np.ndarray:
    """Return e^M by scaling and squaring of a truncated Taylor series (as built in P02.05)."""
    M = np.asarray(M, dtype=float)
    norm = np.linalg.norm(M, 1)
    s = max(0, int(np.ceil(np.log2(norm / 0.5)))) if norm > 0.5 else 0
    Ms = M / 2.0 ** s
    term = np.eye(M.shape[0])
    E = np.eye(M.shape[0])
    for k in range(1, n_terms):
        term = term @ Ms / k
        E = E + term
    for _ in range(s):
        E = E @ E
    return E


def pitch_model() -> tuple[np.ndarray, np.ndarray]:
    """Return (A, B) of a representative small-UAV pitch model with a first-order elevator servo.

    alpha' = -L_alpha alpha + q;  q' = M_alpha alpha + M_q q + M_delta delta_e;  theta' = q;
    delta_e' = (delta_c - delta_e) / tau_s, with L_alpha = 1.2 1/s, M_alpha = -9 1/s^2, M_q = -2.5 1/s,
    M_delta = 12 1/s^2, tau_s = 0.05 s (short period: omega_n = 3.46 rad/s, zeta = 0.53). Representative
    numbers, not a specific aircraft's data.
    """
    La, Ma, Mq, Md, ts = 1.2, -9.0, -2.5, 12.0, 0.05
    A = np.array([[-La, 1.0, 0.0, 0.0],
                  [Ma, Mq, 0.0, Md],
                  [0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 0.0, -1.0 / ts]])
    B = np.array([[0.0], [0.0], [0.0], [1.0 / ts]])
    return A, B


def pitch_gains() -> np.ndarray:
    """Return the continuous-time pitch attitude gains K = [0, Kd, Kp, 0] (Kp = 2.0, Kd = 0.5 s) for
    delta_c = Kp (theta_ref - theta) - Kd q, i.e. u = -K x + Kp theta_ref."""
    return np.array([[0.0, 0.5, 2.0, 0.0]])


def freqresp_s(num: np.ndarray, den: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Return G(jw) = num(jw)/den(jw) for polynomial coefficient arrays in s."""
    s = 1j * np.asarray(w, dtype=float)
    return np.polyval(num, s) / np.polyval(den, s)


def step_metrics(t: np.ndarray, y: np.ndarray, y_final: float) -> dict[str, float]:
    """Return overshoot [fraction of y_final] and 2 % settling time of a sampled step response."""
    os_ = max(0.0, float(np.max(y) / y_final - 1.0))
    out = np.where(np.abs(y - y_final) > 0.02 * abs(y_final))[0]
    ts = float(t[out[-1] + 1]) if len(out) and out[-1] + 1 < len(t) else (float("nan") if len(out) else 0.0)
    return {"overshoot": os_, "settling": ts}


# ---------------------------------------------------------------- discretisation

def c2d_zoh(A: np.ndarray, B: np.ndarray, T: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (Ad, Bd) = (e^{AT}, int_0^T e^{As} ds B), the exact ZOH discretisation, from one exponential
    of the augmented matrix [[A, B], [0, 0]] T."""
    raise NotImplementedError


def c2d_tustin(A: np.ndarray, B: np.ndarray, C: np.ndarray, D: np.ndarray,
               T: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (Ad, Bd, Cd, Dd) of the Tustin (bilinear) discretisation in state space.

    With W = I - A T/2:  Ad = W^{-1}(I + A T/2),  Bd = W^{-1} B T,  Cd = C W^{-1},  Dd = D + C W^{-1} B T/2.
    (The same realisation ``scipy.signal.cont2discrete`` returns for method="bilinear".)
    """
    raise NotImplementedError


def c2d_zoh_delay(A: np.ndarray, B: np.ndarray, T: float, tau: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (Phi, Gam) of the ZOH discretisation with a computational delay 0 <= tau <= T.

    The input computed at t_k (from x_k) is applied at t_k + tau, so over [t_k, t_{k+1}] the plant sees
    u_{k-1} for tau seconds and u_k for T - tau seconds:
        x_{k+1} = e^{AT} x_k + G1 u_{k-1} + G0 u_k,
        G0 = int_0^{T-tau} e^{As} ds B,   G1 = e^{A(T-tau)} int_0^{tau} e^{As} ds B.
    With the augmented state z_k = (x_k, u_{k-1}): z_{k+1} = Phi z_k + Gam u_k, Phi = [[e^{AT}, G1], [0, 0]],
    Gam = [[G0], [I]]. tau = 0 reduces to plain ZOH (G1 = 0); tau = T to a pure one-sample delay (G0 = 0).
    Raise ``ValueError`` for tau outside [0, T].
    """
    raise NotImplementedError


def tustin_tf(num: np.ndarray, den: np.ndarray, T: float,
              w_prewarp: float | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Return (numz, denz) of G(s) with s = c (z - 1)/(z + 1), c = 2/T, or c = w0 / tan(w0 T / 2) if
    ``w_prewarp`` = w0 is given (then the discrete response equals G(j w0) exactly at w0).

    Multiply numerator and denominator by (z + 1)^n, n = max(deg num, deg den), so both are polynomials
    in z of degree n; normalise so denz[0] = 1.
    """
    raise NotImplementedError


def matched_tf(num: np.ndarray, den: np.ndarray, T: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (numz, denz) by matched pole-zero: every pole p -> e^{pT}, every finite zero q -> e^{qT}, every
    zero at infinity -> z = -1, and the gain chosen so the DC gains agree, Gd(1) = G(0).

    Raise ``ValueError`` if G has a pole at s = 0 (the DC gain is then infinite; match elsewhere).
    """
    raise NotImplementedError


# ---------------------------------------------------------------- z-domain tools

def dlsim(Ad: np.ndarray, Bd: np.ndarray, Cd: np.ndarray, Dd: np.ndarray, u: np.ndarray,
          x0: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Return (y, x) of x_{k+1} = Ad x_k + Bd u_k, y_k = Cd x_k + Dd u_k for k = 0..N-1.

    ``u`` has shape (N,) for one input or (N, m); y has shape (N, p) and x (N, n) (x[k] is x_k).
    """
    raise NotImplementedError


def impulse_by_long_division(numz: np.ndarray, denz: np.ndarray, n: int) -> np.ndarray:
    """Return h_0..h_{n-1}, the inverse z-transform of numz(z)/denz(z) as a power series in z^{-1}.

    Pad numz to the length of denz (a proper transfer function), then divide: with a = denz/denz[0] and
    b the padded numerator, h_k = b_k - sum_{i=1..min(k, len(a)-1)} a_i h_{k-i} (b_k = 0 beyond its length).
    Raise ``ValueError`` if the transfer function is improper.
    """
    raise NotImplementedError


def freqresp_z(numz: np.ndarray, denz: np.ndarray, w: np.ndarray, T: float) -> np.ndarray:
    """Return Gd(e^{jwT}) for polynomial coefficient arrays in z, at frequencies w [rad/s]."""
    raise NotImplementedError


def spectral_radius(M: np.ndarray) -> float:
    """Return max |eigenvalue| of M; the discrete system x_{k+1} = M x_k is asymptotically stable iff it is < 1."""
    raise NotImplementedError


def z_to_s(zp: np.ndarray, T: float) -> np.ndarray:
    """Return the equivalent s-plane poles s = ln(z)/T (principal branch) of z-plane poles."""
    raise NotImplementedError


# ---------------------------------------------------------------- the digital pitch loop

def discrete_pitch_loop(T: float, tau: float, K: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Return (Acl, Bcl) of the sampled pitch loop with computational delay tau, on the augmented state
    z_k = (x_k, u_{k-1}): z_{k+1} = Acl z_k + Bcl theta_ref, with u_k = -K x_k + Kp theta_ref."""
    raise NotImplementedError


def hybrid_simulation(T: float, tau: float, t_end: float, theta_ref: float = 1.0, n_sub: int = 20,
                      K: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Return (t, x) of the continuous pitch plant under the sampled controller, resolved between samples.

    The plant is advanced exactly (matrix exponential) on a fine grid of T/n_sub with the elevator command
    switching from u_{k-1} to u_k at t_k + tau. Independent of ``c2d_zoh_delay``: it never forms G0 or G1.
    Requires tau to be a multiple of T/n_sub.
    """
    A, B = pitch_model()
    K = pitch_gains() if K is None else np.atleast_2d(K)
    h = T / n_sub
    j_tau = int(round(tau / h))
    if abs(j_tau * h - tau) > 1e-12 * T:
        raise ValueError("tau must be a multiple of T / n_sub")
    Ah, Bh = c2d_zoh(A, B, h)
    N = int(round(t_end / T))
    x = np.zeros(4)
    u_prev = 0.0
    ts, xs = [0.0], [x.copy()]
    for k in range(N):
        u_k = float(-(K @ x)[0] + K[0, 2] * theta_ref)
        for j in range(n_sub):
            u = u_prev if j < j_tau else u_k
            x = Ah @ x + Bh[:, 0] * u
            ts.append(k * T + (j + 1) * h)
            xs.append(x.copy())
        u_prev = u_k
    return np.array(ts), np.array(xs)


def discrete_loop_margin(T: float, tau: float, K: np.ndarray | None = None) -> tuple[float, float]:
    """Return (phase margin [deg], crossover [rad/s]) of the sampled pitch loop broken at the controller output.

    L(e^{jwT}) = [K, 0] (e^{jwT} I - Phi)^{-1} Gam with (Phi, Gam) from ``c2d_zoh_delay``, evaluated on a
    log grid up to the Nyquist frequency pi/T; the crossover is the first fall of |L| through 1.
    """
    A, B = pitch_model()
    K = pitch_gains() if K is None else np.atleast_2d(K)
    Phi, Gam = c2d_zoh_delay(A, B, T, tau)
    Kaug = np.hstack([K, np.zeros((1, 1))])
    w = np.logspace(-2, np.log10(np.pi / T), 4000)
    n = Phi.shape[0]
    L = np.array([(Kaug @ np.linalg.solve(np.exp(1j * wi * T) * np.eye(n) - Phi, Gam))[0, 0] for wi in w])
    idx = np.where((np.abs(L[:-1]) >= 1.0) & (np.abs(L[1:]) < 1.0))[0]
    if len(idx) == 0:
        return float("inf"), float("nan")
    i = idx[0]
    ph = np.unwrap(np.angle(L))
    return float(180.0 + np.degrees(ph[i])), float(w[i])


def sample_rate_study(rates_hz: np.ndarray, delay_fracs: tuple[float, ...] = (0.0, 0.5, 1.0),
                      t_end: float = 6.0) -> list[dict[str, float]]:
    """For each sample rate and delay fraction tau/T, return the closed-loop spectral radius, the dominant
    equivalent damping ratio (of ln z / T for the least-damped complex pole), the theta step overshoot and
    the phase margin from ``discrete_loop_margin``."""
    rows = []
    for fs in rates_hz:
        T = 1.0 / fs
        for frac in delay_fracs:
            Acl, Bcl = discrete_pitch_loop(T, frac * T)
            lam = np.linalg.eigvals(Acl)
            rho = float(np.max(np.abs(lam)))
            s = z_to_s(lam[np.abs(lam) > 1e-9], T)
            cplx = s[np.abs(s.imag) > 1e-6]
            zeta = float(np.min(-cplx.real / np.abs(cplx))) if len(cplx) else 1.0
            N = int(round(t_end / T))
            _, zs = dlsim(Acl, Bcl, np.eye(5), np.zeros((5, 1)), np.ones(N))
            theta = zs[:, 2]
            os_ = float(max(0.0, theta.max() - 1.0)) if rho < 1 else float("inf")
            pm, _ = discrete_loop_margin(T, frac * T)
            rows.append({"fs": float(fs), "frac": float(frac), "rho": rho, "zeta": zeta, "overshoot": os_, "pm": pm})
    return rows


def main(out_dir: str | Path = ".") -> None:
    """Print the discretisation checks and the pitch-loop sample-rate study, and save digital_pitch.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    A, B = pitch_model()
    K = pitch_gains()
    lam_c = np.linalg.eigvals(A - B @ K)
    print("pitch model eigenvalues (open loop):", np.round(np.linalg.eigvals(A), 3))
    print("continuous closed loop (Kp = 2.0, Kd = 0.5 s):", np.round(lam_c, 3))
    w = np.logspace(-2, 3, 20000)
    L = np.array([(K @ np.linalg.solve(1j * wi * np.eye(4) - A, B))[0, 0] for wi in w])
    i = np.where(np.abs(L) < 1)[0][0]
    wc, pm = w[i], 180.0 + np.degrees(np.angle(L[i]))
    print(f"continuous loop at the servo input: crossover {wc:.2f} rad/s, phase margin {pm:.1f} deg")
    for fs in (50.0, 20.0, 10.0):
        T = 1.0 / fs
        print(f"  fs = {fs:4.0f} Hz: ZOH phase lag wc T/2 = {np.degrees(wc * T / 2):5.1f} deg, one-sample delay adds"
              f" {np.degrees(wc * T):5.1f} deg -> predicted margin {pm - np.degrees(1.5 * wc * T):5.1f} deg")

    # discretisation of the plant at 50 Hz
    T = 0.02
    Ad, Bd = c2d_zoh(A, B, T)
    print(f"\nZOH at 50 Hz: eig(Ad) = {np.round(np.linalg.eigvals(Ad), 4)} = e^(lambda T): "
          f"max error {np.max(np.abs(np.sort_complex(np.linalg.eigvals(Ad)) - np.sort_complex(np.exp(np.linalg.eigvals(A) * T)))):.1e}")

    # a lead compensator, Tustin and matched pole-zero, compared with the continuous one
    num, den = np.array([0.5, 2.0]) * 1.0, np.array([1.0 / 20.0, 1.0])   # C(s) = (0.5 s + 2)/(s/20 + 1)
    for fs in (50.0, 20.0):
        T = 1.0 / fs
        nt, dt_ = tustin_tf(num, den, T)
        npw, dpw = tustin_tf(num, den, T, w_prewarp=wc)
        nm, dm = matched_tf(num, den, T)
        Cc = freqresp_s(num, den, [wc])[0]
        for name, (nz, dz) in (("Tustin", (nt, dt_)), ("Tustin prewarped at wc", (npw, dpw)), ("matched", (nm, dm))):
            Cz = freqresp_z(nz, dz, [wc], T)[0]
            print(f"  lead at {fs:.0f} Hz, {name:22s}: |C| ratio {abs(Cz) / abs(Cc):.4f}, phase error {np.degrees(np.angle(Cz / Cc)):+6.2f} deg at wc")

    rates = np.array([100.0, 50.0, 20.0, 10.0, 7.0])
    rows = sample_rate_study(rates)
    print("\nsampled pitch loop: fs [Hz], tau/T, spectral radius, equivalent zeta, theta overshoot, PM measured, PM rule")
    for r in rows:
        T = 1.0 / r["fs"]
        rule = pm - np.degrees(wc * T / 2 + wc * r["frac"] * T)
        print(f"  {r['fs']:5.0f} {r['frac']:4.1f} {r['rho']:.4f} {r['zeta']:6.3f} {r['overshoot']:6.3f} {r['pm']:6.1f} {rule:6.1f}")
    # continuous-limit check
    Acl, _ = discrete_pitch_loop(1e-3, 0.0)
    lam_d = np.linalg.eigvals(Acl)
    s_eq = np.sort_complex(z_to_s(lam_d[np.abs(lam_d) > 1e-9], 1e-3))
    print(f"fs = 1 kHz, no delay: ln(z)/T = {np.round(s_eq, 2)} vs continuous {np.round(np.sort_complex(lam_c), 2)}")

    # hybrid simulation (between samples) against the augmented difference equation
    T, tau = 0.05, 0.025
    th, xh = hybrid_simulation(T, tau, 4.0)
    Acl, Bcl = discrete_pitch_loop(T, tau)
    _, zs = dlsim(Acl, Bcl, np.eye(5), np.zeros((5, 1)), np.ones(80))
    print(f"20 Hz, tau = T/2: hybrid simulation vs difference equation at the samples, max |dx| = {np.max(np.abs(xh[::20][:80] - zs[:, :4])):.1e}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    Acl_c = A - B @ K
    tc = np.linspace(0.0, 4.0, 801)
    yc = [np.linalg.solve(Acl_c, (expm(Acl_c * tt) - np.eye(4)) @ B * K[0, 2])[2, 0] for tt in tc]
    axes[0].plot(tc, yc, "k", lw=1.2, label="continuous")
    for fs, frac, col in ((50.0, 1.0, "C0"), (20.0, 1.0, "C1"), (10.0, 1.0, "C3")):
        T = 1.0 / fs
        th, xh = hybrid_simulation(T, frac * T, 4.0, n_sub=10)
        axes[0].plot(th, xh[:, 2], col, lw=0.9, label=f"{fs:.0f} Hz, delay {frac:.0f} T")
    axes[0].set_xlabel("t [s]"); axes[0].set_ylabel("theta / theta_ref"); axes[0].legend(fontsize=8); axes[0].grid(lw=0.4)
    circ = np.linspace(0, 2 * np.pi, 400)
    axes[1].plot(np.cos(circ), np.sin(circ), "k:", lw=0.8)
    for fs, col in ((50.0, "C0"), (20.0, "C1"), (10.0, "C3"), (7.0, "C4")):
        Acl, _ = discrete_pitch_loop(1.0 / fs, 1.0 / fs)
        lam = np.linalg.eigvals(Acl)
        axes[1].plot(lam.real, lam.imag, "x", color=col, label=f"{fs:.0f} Hz, delay T")
    axes[1].set_aspect("equal"); axes[1].legend(fontsize=8); axes[1].grid(lw=0.4)
    axes[1].set_title("closed-loop poles in the z-plane", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "digital_pitch.png", dpi=150)
    plt.close(fig)
    print(f"saved {out_dir / 'digital_pitch.png'}")


if __name__ == "__main__":
    main()
