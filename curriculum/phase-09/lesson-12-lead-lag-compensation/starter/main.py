"""P09.12 Lead-Lag Compensation: starter code.

Loop shaping a control-surface servo actuator, G(s) = K / (s (tau s + 1)): a lead compensator buys phase
margin at a chosen crossover frequency, a lag compensator buys low-frequency gain (the velocity error
constant Kv) without moving the crossover, and the design is checked by measuring the margins of the
compensated loop on its own frequency response and by its closed-loop step and ramp responses.

Conventions:
    * Polynomials are NumPy coefficient arrays, highest power first (as ``np.polyval``).
    * Lead: C_lead(s) = Kc (T s + 1) / (alpha T s + 1), 0 < alpha < 1 (zero below the pole).
    * Lag:  C_lag(s) = (s + z) / (s + z / beta), beta > 1: DC gain beta, high-frequency gain 1.
    * Angles in degrees at the interface (phase margin, maximum lead phase); radians inside.

Complete every function marked NotImplementedError.

Run:  python main.py              (prints the design, writes lead_lag.png)
Test: pytest curriculum/phase-09/lesson-12-lead-lag-compensation
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


# ---------------------------------------------------------------- given: the plant and small tools

def servo_plant() -> dict[str, float]:
    """Return a representative DC-motor control-surface servo: K = 20 1/s (surface rad/s per unit command), tau = 0.05 s."""
    return {"K": 20.0, "tau": 0.05}


def servo_tf(p: dict[str, float]) -> tuple[np.ndarray, np.ndarray]:
    """Return (num, den) of G(s) = K / (s (tau s + 1))."""
    return np.array([p["K"]]), np.array([p["tau"], 1.0, 0.0])


def series(n1: np.ndarray, d1: np.ndarray, n2: np.ndarray, d2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the product of two transfer functions (n1 n2, d1 d2)."""
    return np.polymul(n1, n2), np.polymul(d1, d2)


def unity_feedback(n: np.ndarray, d: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the closed loop L / (1 + L) of a loop transfer function L = n/d: (n, d + n)."""
    return np.asarray(n, dtype=float), np.polyadd(d, n)


def expm(M: np.ndarray, n_terms: int = 20) -> np.ndarray:
    """Return e^M by scaling and squaring a truncated Taylor series (the method of P02.05)."""
    M = np.asarray(M, dtype=float)
    norm = np.linalg.norm(M, 1)
    s = max(0, int(np.ceil(np.log2(norm / 0.5)))) if norm > 0.5 else 0
    X = M / 2 ** s
    term = np.eye(len(M))
    E = np.eye(len(M))
    for k in range(1, n_terms):
        term = term @ X / k
        E = E + term
    for _ in range(s):
        E = E @ E
    return E


def tf_response(num: np.ndarray, den: np.ndarray, t: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Return the response of the proper num/den to the input samples u held constant over each step of uniform t.

    Controllable canonical realisation, stepped exactly by the zero-order-hold matrix exponential.
    """
    num = np.atleast_1d(np.asarray(num, dtype=float))
    den = np.atleast_1d(np.asarray(den, dtype=float))
    num, den = num / den[0], den / den[0]
    n = len(den) - 1
    num = np.concatenate([np.zeros(n + 1 - len(num)), num])
    D = num[0]
    c = num[1:] - D * den[1:]
    A = np.zeros((n, n))
    A[0, :] = -den[1:]
    A[1:, :-1] = np.eye(n - 1)
    dt = t[1] - t[0]
    M = np.zeros((n + 1, n + 1))
    M[:n, :n] = A
    M[0, n] = 1.0
    E = expm(M * dt)
    Ad, Bd = E[:n, :n], E[:n, n]
    x = np.zeros(n)
    y = np.empty(len(t))
    for k in range(len(t)):
        y[k] = c @ x + D * u[k]
        x = Ad @ x + Bd * u[k]
    return y


def step_metrics(t: np.ndarray, y: np.ndarray, y_final: float = 1.0, band: float = 0.02) -> dict[str, float]:
    """Return overshoot [%], 10-90 % rise time [s] and settling time into +-band [s] of a step response."""
    i10 = int(np.argmax(y >= 0.1 * y_final))
    i90 = int(np.argmax(y >= 0.9 * y_final))
    outside = np.where(np.abs(y - y_final) > band * abs(y_final))[0]
    ts = t[min(outside[-1] + 1, len(t) - 1)] if len(outside) else t[0]
    return {"overshoot_pct": float(max(100.0 * (y.max() - y_final) / y_final, 0.0)),
            "rise_time": float(t[i90] - t[i10]), "settling_time": float(ts)}


# ---------------------------------------------------------------- frequency response and margins

def freq_response(num: np.ndarray, den: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Return G(jw) = num(jw) / den(jw) at the angular frequencies w [rad/s]."""
    raise NotImplementedError


def unwrapped_phase_deg(num: np.ndarray, den: np.ndarray, w: np.ndarray, w_start: float = 1e-3) -> np.ndarray:
    """Return the continuous phase of G(jw) in degrees, unwrapped along a log sweep from w_start up to each w."""
    raise NotImplementedError


def margins(num: np.ndarray, den: np.ndarray, w_lo: float = 1e-2, w_hi: float = 1e4) -> dict[str, float]:
    """Return gain crossover wc, phase margin pm [deg], phase crossover w180, gain margin gm (ratio) and delay margin [s].

    Sweep |L(jw)| and the unwrapped phase on a log grid, bracket every crossing of |L| = 1 and of
    phase = -180 deg (mod 360), and refine each by bisection. pm and gm are the smallest over the crossings;
    gm is inf when the phase never reaches -180 deg; the delay margin is pm (rad) / wc.
    """
    raise NotImplementedError


def velocity_constant(num: np.ndarray, den: np.ndarray) -> float:
    """Return Kv = lim_{s->0} s L(s) for a type-1 loop L = num/den (den has exactly one root at s = 0)."""
    raise NotImplementedError


# ---------------------------------------------------------------- lead and lag

def lead_max_phase(alpha: float) -> tuple[float, float]:
    """Return (phi_max [deg], |C| at w_max) of (T s + 1)/(alpha T s + 1): sin phi_max = (1 - alpha)/(1 + alpha), gain 1/sqrt(alpha).

    The maximum occurs at the geometric mean of the zero and pole, w_max = 1/(T sqrt(alpha)).
    """
    raise NotImplementedError


def lead_alpha(phi_max_deg: float) -> float:
    """Return alpha = (1 - sin phi)/(1 + sin phi), the pole-zero ratio that gives a maximum lead of phi_max_deg."""
    raise NotImplementedError


def lead_compensator(Kc: float, alpha: float, T: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (num, den) of Kc (T s + 1)/(alpha T s + 1)."""
    raise NotImplementedError


def lag_compensator(beta: float, z: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (num, den) of (s + z)/(s + z/beta): DC gain beta, unit gain well above z."""
    raise NotImplementedError


def design_lead(num: np.ndarray, den: np.ndarray, wc: float, pm_deg: float, max_lead_deg: float = 65.0) -> tuple[float, float, float]:
    """Return (Kc, alpha, T) of the lead that puts the gain crossover at wc with phase margin pm_deg.

    phi = pm - (180 + angle G(j wc)) is the phase the lead must add at wc; centring the lead's maximum there
    gives alpha = lead_alpha(phi), T = 1/(wc sqrt(alpha)), and Kc = 1/|C_lead,1(j wc) G(j wc)| makes |L| = 1 at wc.
    If no lead is needed (phi <= 0), alpha = 1 and the compensator is the gain Kc. Raise ``ValueError`` when
    phi exceeds max_lead_deg (one lead stage cannot do it sensibly; use two, or lower the crossover).
    """
    raise NotImplementedError


def design_lead_lag(num: np.ndarray, den: np.ndarray, wc: float, pm_deg: float, kv: float,
                    lag_ratio: float = 10.0) -> dict:
    """Design lead + lag for crossover wc, phase margin pm_deg and velocity constant at least kv.

    The lag zero sits at z = wc/lag_ratio. Because the lag changes the phase and gain at wc slightly, the lead is
    designed on G C_lag, and beta is iterated: beta <- beta kv / Kv(achieved) until it settles (beta = 1 if the
    lead alone already meets kv). Returns Kc, alpha, T, beta, z and the compensator and loop (num, den).
    """
    raise NotImplementedError


# ---------------------------------------------------------------- the example

def main(out_dir: str | Path = ".") -> None:
    """Shape the servo loop with lead and lag, measure the margins and the closed-loop step, save lead_lag.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    p = servo_plant()
    ng, dg = servo_tf(p)
    wc, pm, kv = 30.0, 55.0, 200.0
    print(f"representative servo: G(s) = {p['K']:.0f} / (s ({p['tau']} s + 1));  targets wc = {wc} rad/s, PM = {pm} deg, Kv >= {kv} 1/s")

    designs = {}
    Kp = kv / p["K"]
    designs["gain only (Kv met)"] = (np.array([Kp]), np.array([1.0]))
    m0 = margins(*series(np.array([Kp]), np.array([1.0]), ng, dg))
    print(f"gain only, Kc = {Kp:.1f} for Kv = {kv:.0f}: wc = {m0['wc']:.2f} rad/s, PM = {m0['pm']:.2f} deg")
    phase_g = unwrapped_phase_deg(ng, dg, np.array([wc]))[0]
    print(f"plant phase at wc = {phase_g:.2f} deg: a pure gain crossing at {wc} rad/s leaves PM = {180 + phase_g:.2f} deg")

    Kc, alpha, T = design_lead(ng, dg, wc, pm)
    phi, gmax = lead_max_phase(alpha)
    nl, dl = lead_compensator(Kc, alpha, T)
    designs["lead"] = (nl, dl)
    mL = margins(*series(nl, dl, ng, dg))
    kv_lead = velocity_constant(*series(nl, dl, ng, dg))
    print(f"lead: phi_max = {phi:.2f} deg, alpha = {alpha:.4f}, zero {1 / T:.2f} rad/s, pole {1 / (alpha * T):.2f} rad/s, Kc = {Kc:.4f};"
          f" measured wc = {mL['wc']:.4f}, PM = {mL['pm']:.4f} deg, GM = {mL['gm']}, Kv = {kv_lead:.2f} 1/s")

    d = design_lead_lag(ng, dg, wc, pm, kv)
    designs["lead-lag"] = d["comp"]
    mLL = margins(*d["loop"])
    print(f"lead-lag: beta = {d['beta']:.4f}, lag zero {d['z']:.2f}, lag pole {d['z'] / d['beta']:.4f} rad/s; lead alpha = {d['alpha']:.4f},"
          f" phi_max = {lead_max_phase(d['alpha'])[0]:.2f} deg, Kc = {d['Kc']:.4f}")
    print(f"          measured wc = {mLL['wc']:.4f} rad/s, PM = {mLL['pm']:.4f} deg, GM = {mLL['gm']}, delay margin {1e3 * mLL['delay_margin']:.2f} ms,"
          f" Kv = {velocity_constant(*d['loop']):.2f} 1/s")

    t = np.linspace(0.0, 1.5, 3001)
    steps = {}
    for name, (nc, dc) in designs.items():
        ncl, dcl = unity_feedback(*series(nc, dc, ng, dg))
        y = tf_response(ncl, dcl, t, np.ones_like(t))
        steps[name] = y
        m = step_metrics(t, y)
        nL, dL = series(nc, dc, ng, dg)
        tr = np.linspace(0.0, 10.0, 20001)
        # unit ramp command; held at each step's midpoint value so the staircase has no half-sample lag
        e_ramp = tr - tf_response(ncl, dcl, tr, tr + 0.5 * (tr[1] - tr[0]))
        print(f"  {name:20s}: overshoot {m['overshoot_pct']:5.1f} %, rise {1e3 * m['rise_time']:5.1f} ms, settle {1e3 * m['settling_time']:6.1f} ms,"
              f" ramp error at 10 s {1e3 * e_ramp[-1]:.3f} mrad per rad/s (1/Kv = {1e3 / velocity_constant(nL, dL):.3f})")

    print("\nphase-margin sweep (lead only, wc = 30 rad/s):")
    for pmi in (40.0, 50.0, 60.0, 70.0):
        try:
            Kc_i, a_i, T_i = design_lead(ng, dg, wc, pmi)
        except ValueError as exc:
            print(f"  PM {pmi:.0f}: {exc}")
            continue
        nli, dli = lead_compensator(Kc_i, a_i, T_i)
        ncl, dcl = unity_feedback(*series(nli, dli, ng, dg))
        m = step_metrics(t, tf_response(ncl, dcl, t, np.ones_like(t)))
        print(f"  PM {pmi:.0f} deg: alpha = {a_i:.4f}, overshoot {m['overshoot_pct']:.1f} %, settle {1e3 * m['settling_time']:.1f} ms")

    print("\nlag placement sweep (lead-lag, PM target 55 deg):")
    for ratio in (3.0, 5.0, 10.0, 20.0):
        di = design_lead_lag(ng, dg, wc, pm, kv, lag_ratio=ratio)
        ncl, dcl = unity_feedback(*di["loop"])
        tl = np.linspace(0.0, 4.0, 8001)
        m = step_metrics(tl, tf_response(ncl, dcl, tl, np.ones_like(tl)))
        lag_phase = np.degrees(np.angle(freq_response(*lag_compensator(di["beta"], di["z"]), np.array([wc]))[0]))
        print(f"  z = wc/{ratio:.0f}: lag phase at wc {lag_phase:.2f} deg, lead alpha {di['alpha']:.4f}, overshoot {m['overshoot_pct']:.1f} %,"
              f" settle {1e3 * m['settling_time']:.0f} ms")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    w = np.logspace(-1, 3, 800)
    for name, (nc, dc) in designs.items():
        nL, dL = series(nc, dc, ng, dg)
        axes[0].semilogx(w, 20 * np.log10(np.abs(freq_response(nL, dL, w))), label=name)
        axes[1].semilogx(w, unwrapped_phase_deg(nL, dL, w), label=name)
        axes[2].plot(t * 1e3, steps[name], label=name)
    axes[0].axhline(0.0, color="k", lw=0.6); axes[0].axvline(wc, color="k", lw=0.6, ls=":")
    axes[0].set_xlabel("w [rad/s]"); axes[0].set_ylabel("|L| [dB]"); axes[0].legend(fontsize=8); axes[0].grid(True, which="both", lw=0.3)
    axes[1].axhline(-180.0 + pm, color="k", lw=0.6, ls="--"); axes[1].axhline(-180.0, color="k", lw=0.6)
    axes[1].set_xlabel("w [rad/s]"); axes[1].set_ylabel("phase L [deg]"); axes[1].grid(True, which="both", lw=0.3)
    axes[2].set_xlabel("t [ms]"); axes[2].set_ylabel("surface angle / command"); axes[2].legend(fontsize=8); axes[2].grid(lw=0.4)
    axes[2].set_title("closed-loop step (simulation)", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "lead_lag.png", dpi=150)
    plt.close(fig)
    print(f"saved {out_dir / 'lead_lag.png'}")


if __name__ == "__main__":
    main()
