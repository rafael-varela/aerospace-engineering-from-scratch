"""P09.11 PID Design & Tuning: reference solution.

A pitch-attitude PID for a small fixed-wing aircraft, written as a discrete-time controller with a
filtered derivative, setpoint weighting, an elevator command limit and two anti-windup schemes
(clamping, also called conditional integration, and back-calculation). The plant is a representative
short-period model with an integrator to pitch angle and a first-order elevator servo; it is tuned by
Ziegler-Nichols (ultimate gain from the j-omega crossing of the root locus) and by assigning three
closed-loop poles exactly, and the simulated linear loop is checked against the closed-loop transfer
function.

Conventions:
    * Plant state x = (alpha, q, theta, delta): angle of attack [rad], pitch rate [rad/s], pitch angle
      [rad], elevator deflection [rad]. delta > 0 is trailing edge up (nose-up moment), so M_delta > 0.
    * Polynomials are NumPy coefficient arrays, highest power first (as ``np.polyval``).
    * gains = (Kp, Ki, Kd, Tf): C(s) = Kp + Ki/s + Kd s / (Tf s + 1).

Run:  python solution.py          (prints the design, writes pid_pitch.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-09/lesson-11-pid-design-tuning
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


# ---------------------------------------------------------------- given: the plant and small tools

def pitch_plant() -> dict[str, float]:
    """Return a representative (not aircraft-specific) short-period pitch model of a light aircraft / small UAV.

    Z_alpha [1/s], M_alpha [1/s^2], M_q [1/s], M_delta [1/s^2], servo bandwidth w_servo [rad/s] (time constant
    1/w_servo = 0.05 s) and the elevator command limit delta_max [rad] (15 degrees).
    """
    return {"Z_alpha": -1.2, "M_alpha": -6.0, "M_q": -2.5, "M_delta": 12.0, "w_servo": 20.0,
            "delta_max": float(np.radians(15.0))}


def pitch_state_space(p: dict[str, float]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (A, B, C) for x = (alpha, q, theta, delta), input the elevator command, output theta.

    alpha' = Z_alpha alpha + q;  q' = M_alpha alpha + M_q q + M_delta delta;  theta' = q;
    delta' = w_servo (delta_cmd - delta).
    """
    A = np.array([[p["Z_alpha"], 1.0, 0.0, 0.0],
                  [p["M_alpha"], p["M_q"], 0.0, p["M_delta"]],
                  [0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 0.0, -p["w_servo"]]])
    B = np.array([[0.0], [0.0], [0.0], [p["w_servo"]]])
    C = np.array([[0.0, 0.0, 1.0, 0.0]])
    return A, B, C


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


def discretize_zoh(A: np.ndarray, B: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (A_d, B_d) of the exact zero-order-hold discretisation from e^{[[A, B], [0, 0]] dt}."""
    n, m = B.shape
    M = np.zeros((n + m, n + m))
    M[:n, :n] = A
    M[:n, n:] = B
    E = expm(M * dt)
    return E[:n, :n], E[:n, n:]


def tf_step(num: np.ndarray, den: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Return the unit-step response of num/den (proper) at the uniformly spaced times t, exactly (ZOH of a constant).

    Builds the controllable canonical realisation and steps it with the matrix exponential.
    """
    num = np.atleast_1d(np.asarray(num, dtype=float))
    den = np.atleast_1d(np.asarray(den, dtype=float))
    num = num / den[0]
    den = den / den[0]
    n = len(den) - 1
    num = np.concatenate([np.zeros(n + 1 - len(num)), num])
    D = num[0]
    b = num[1:] - D * den[1:]
    A = np.zeros((n, n))
    A[0, :] = -den[1:]
    A[1:, :-1] = np.eye(n - 1)
    B = np.zeros((n, 1))
    B[0, 0] = 1.0
    C = b[None, :]
    dt = t[1] - t[0]
    Ad, Bd = discretize_zoh(A, B, dt)
    x = np.zeros(n)
    y = np.empty(len(t))
    for k in range(len(t)):
        y[k] = C[0] @ x + D
        x = Ad @ x + Bd[:, 0]
    return y


def step_metrics(t: np.ndarray, y: np.ndarray, y_final: float, band: float = 0.02) -> dict[str, float]:
    """Return overshoot [%], 10-90 % rise time [s], settling time into +-band of y_final [s] and the last value."""
    os_pct = 100.0 * (np.max(y) - y_final) / abs(y_final)
    i10 = int(np.argmax(y >= 0.1 * y_final))
    i90 = int(np.argmax(y >= 0.9 * y_final))
    outside = np.where(np.abs(y - y_final) > band * abs(y_final))[0]
    ts = t[min(outside[-1] + 1, len(t) - 1)] if len(outside) else t[0]
    return {"overshoot_pct": float(max(os_pct, 0.0)), "rise_time": float(t[i90] - t[i10]),
            "settling_time": float(ts), "final": float(y[-1])}


# ---------------------------------------------------------------- transfer functions

def pitch_plant_tf(p: dict[str, float]) -> tuple[np.ndarray, np.ndarray]:
    """Return (num, den) of theta/delta_cmd = M_delta w_s (s - Z_alpha) / (s (s^2 - (Z_alpha + M_q) s + Z_alpha M_q - M_alpha) (s + w_s)).

    Derived by eliminating alpha and q from the short-period equations (Z_delta = 0) and multiplying by the servo.
    """
    za, ma, mq, md, ws = p["Z_alpha"], p["M_alpha"], p["M_q"], p["M_delta"], p["w_servo"]
    num = md * ws * np.array([1.0, -za])
    short_period = np.array([1.0, -(za + mq), za * mq - ma])
    den = np.polymul(np.polymul([1.0, 0.0], short_period), [1.0, ws])
    return num, den


def pid_tf(Kp: float, Ki: float, Kd: float, Tf: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (num, den) of C(s) = Kp + Ki/s + Kd s/(Tf s + 1) = ((Kp Tf + Kd) s^2 + (Kp + Ki Tf) s + Ki) / (Tf s^2 + s)."""
    return np.array([Kp * Tf + Kd, Kp + Ki * Tf, Ki]), np.array([Tf, 1.0, 0.0])


def closed_loop_tf(num_p: np.ndarray, den_p: np.ndarray, gains: tuple[float, float, float, float],
                   b: float = 1.0, d_on_measurement: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """Return (num, den) of theta/theta_ref for the two-degree-of-freedom PID loop.

    u = Kp (b r - y) + Ki/s (r - y) + D(s) (x_d), with x_d = -y if ``d_on_measurement`` else r - y and
    D(s) = Kd s/(Tf s + 1). Then y = P (C_r r - C y) and T = P C_r / (1 + P C): the denominator (the
    characteristic polynomial) is the same for every choice of b and derivative input; only the zeros move.
    """
    Kp, Ki, Kd, Tf = gains
    nc, dc = pid_tf(Kp, Ki, Kd, Tf)
    if d_on_measurement:
        ncr = np.polymul([b * Kp, Ki], [Tf, 1.0])
    else:
        ncr = np.array([b * Kp * Tf + Kd, b * Kp + Ki * Tf, Ki])
    den = np.polyadd(np.polymul(den_p, dc), np.polymul(num_p, nc))
    num = np.polymul(num_p, ncr)
    return num, den


# ---------------------------------------------------------------- tuning

def ultimate_gain(num: np.ndarray, den: np.ndarray) -> tuple[float, float, float]:
    """Return (Ku, wu, Pu): the smallest positive proportional gain that puts closed-loop roots on the j-omega axis.

    On s = j w the characteristic equation den + K num = 0 needs K = -den(jw)/num(jw) real and positive, i.e.
    Im[den(jw) num(-jw)] = 0. That condition is a polynomial in w; its positive real roots are the candidate
    crossing frequencies. Pu = 2 pi / wu is the period of the sustained oscillation (Ziegler-Nichols test).
    """
    num = np.asarray(num, dtype=float)
    den = np.asarray(den, dtype=float)
    num_neg = num * (-1.0) ** np.arange(len(num) - 1, -1, -1)       # coefficients of num(-s)
    q = np.polymul(den, num_neg)                                     # q(s) = den(s) num(-s)
    deg = len(q) - 1
    # Im q(jw) = sum over odd k of c_k (-1)^((k-1)/2) w^k  -> polynomial in w (highest power first)
    poly_w = np.zeros(deg + 1)
    for i, c in enumerate(q):
        k = deg - i
        if k % 2 == 1:
            poly_w[i] = c * (-1.0) ** ((k - 1) // 2)
    best = None
    for w in np.roots(np.trim_zeros(poly_w, "f")):
        if abs(w.imag) > 1e-9 * max(1.0, abs(w)) or w.real <= 1e-12:
            continue
        w = w.real
        K = -np.polyval(den, 1j * w) / np.polyval(num, 1j * w)
        if abs(K.imag) < 1e-6 * max(1.0, abs(K)) and K.real > 0 and (best is None or K.real < best[0]):
            best = (float(K.real), float(w))
    if best is None:
        raise ValueError("no positive gain puts the closed loop on the imaginary axis")
    Ku, wu = best
    return Ku, wu, 2 * np.pi / wu


def ziegler_nichols(Ku: float, Pu: float) -> dict[str, tuple[float, float, float]]:
    """Return the classic Ziegler-Nichols ultimate-gain rules as (Kp, Ki, Kd) for 'P', 'PI' and 'PID'.

    P: Kp = 0.5 Ku. PI: Kp = 0.45 Ku, Ti = Pu/1.2. PID: Kp = 0.6 Ku, Ti = Pu/2, Td = Pu/8; Ki = Kp/Ti, Kd = Kp Td.
    """
    kp_pi = 0.45 * Ku
    kp = 0.6 * Ku
    return {"P": (0.5 * Ku, 0.0, 0.0),
            "PI": (kp_pi, kp_pi / (Pu / 1.2), 0.0),
            "PID": (kp, kp / (Pu / 2.0), kp * Pu / 8.0)}


def pid_place_three_poles(num_p: np.ndarray, den_p: np.ndarray, poles: list[complex], Tf: float) -> tuple[float, float, float]:
    """Return (Kp, Ki, Kd) that make three chosen closed-loop poles exact roots of den_p (Tf s^2 + s) + num_p N_c(s).

    N_c = (Kp Tf + Kd) s^2 + (Kp + Ki Tf) s + Ki is linear in the gains, so the characteristic polynomial
    evaluated at a chosen pole is an affine function of (Kp, Ki, Kd). One real pole gives one real equation;
    a complex pair (give both members or only one) gives two (real and imaginary parts). Three equations,
    three gains: a 3x3 linear solve. The remaining closed-loop poles fall where they fall and must be checked.
    """
    num_p = np.asarray(num_p, dtype=float)
    den_p = np.asarray(den_p, dtype=float)

    def char(Kp: float, Ki: float, Kd: float) -> np.ndarray:
        return np.polyadd(np.polymul(den_p, [Tf, 1.0, 0.0]), np.polymul(num_p, [Kp * Tf + Kd, Kp + Ki * Tf, Ki]))

    c0 = char(0.0, 0.0, 0.0)
    cols = [char(*e) - c0 for e in np.eye(3)]
    rows, rhs, used = [], [], []
    for pole in poles:
        pole = complex(pole)
        if any(abs(pole.conjugate() - u) < 1e-12 for u in used) and abs(pole.imag) > 0:
            continue                                                    # the conjugate gives no new equation
        used.append(pole)
        v0 = np.polyval(c0, pole)
        vs = np.array([np.polyval(c, pole) for c in cols])
        rows.append(vs.real)
        rhs.append(-v0.real)
        if abs(pole.imag) > 0:
            rows.append(vs.imag)
            rhs.append(-v0.imag)
    if len(rows) != 3:
        raise ValueError("give exactly three poles (a complex pair counts as two)")
    Kp, Ki, Kd = np.linalg.solve(np.array(rows), np.array(rhs))
    return float(Kp), float(Ki), float(Kd)


# ---------------------------------------------------------------- the discrete controller

def pid_step(state: np.ndarray, r: float, y: float, gains: tuple[float, float, float, float], dt: float,
             u_limit: float = np.inf, anti_windup: str = "none", Tt: float | None = None, b: float = 1.0,
             d_on_measurement: bool = True) -> tuple[float, float, np.ndarray]:
    """Advance the discrete PID one sample; return (u_applied, u_unsaturated, new_state).

    state = (I, D, x_prev): integrator value, filtered derivative, previous derivative input. Zero state means
    every signal was zero before the first sample. One sample:
        e = r - y;  x_d = -y if d_on_measurement else e
        D = (Tf D_prev + Kd (x_d - x_prev)) / (Tf + dt)          (backward-Euler derivative filter)
        v = Kp (b r - y) + I + D;  u = clip(v, -u_limit, u_limit)
        I_next = I + Ki dt e                                     ('none')
                 unchanged if saturated and e pushes further     ('clamp': conditional integration)
                 I + Ki dt e + (dt / Tt) (u - v)                 ('backcalc', Tt defaults to sqrt(Ti Td) or Ti)
    """
    Kp, Ki, Kd, Tf = gains
    I, D_prev, x_prev = state
    e = r - y
    xd = -y if d_on_measurement else e
    D = (Tf * D_prev + Kd * (xd - x_prev)) / (Tf + dt)
    v = Kp * (b * r - y) + I + D
    u = float(np.clip(v, -u_limit, u_limit))
    if anti_windup == "none":
        I_next = I + Ki * dt * e
    elif anti_windup == "clamp":
        saturated = u != v
        pushing = e * (v - u) > 0
        I_next = I if (saturated and pushing) else I + Ki * dt * e
    elif anti_windup == "backcalc":
        if Tt is None:
            Ti = Kp / Ki if Ki > 0 else np.inf
            Td = Kd / Kp if Kp > 0 else 0.0
            Tt = np.sqrt(Ti * Td) if Td > 0 else Ti
        I_next = I + Ki * dt * e + (dt / Tt) * (u - v)
    else:
        raise ValueError(f"unknown anti_windup {anti_windup!r}")
    return u, float(v), np.array([I_next, D, xd])


def simulate_pitch_loop(p: dict[str, float], gains: tuple[float, float, float, float], theta_ref: float,
                        t_end: float, dt: float, u_limit: float = np.inf, anti_windup: str = "none",
                        Tt: float | None = None, b: float = 1.0, d_on_measurement: bool = True,
                        disturbance: float = 0.0, t_disturbance: float = 0.0) -> dict[str, np.ndarray]:
    """Simulate the sampled PID on the continuous pitch plant; return t, theta, u (applied), v (unsaturated), I.

    The reference steps from 0 to theta_ref at t = 0. Each sample: measure theta, run ``pid_step``, hold the
    command for dt; the plant is advanced exactly by its zero-order-hold discretisation. ``disturbance`` is a
    constant pitching moment expressed as an equivalent elevator angle [rad], added at the plant input from
    t_disturbance on (a trim change, e.g. a CG shift or flap deployment).
    """
    A, B, C = pitch_state_space(p)
    # Disturbance enters as a moment M_delta * d on q, independent of the servo.
    Bw = np.array([[0.0], [p["M_delta"]], [0.0], [0.0]])
    Ad, Bd = discretize_zoh(A, np.hstack([B, Bw]), dt)
    n = int(round(t_end / dt))
    t = np.arange(n + 1) * dt
    x = np.zeros(4)
    state = np.zeros(3)
    theta = np.empty(n + 1)
    u_hist = np.empty(n + 1)
    v_hist = np.empty(n + 1)
    i_hist = np.empty(n + 1)
    for k in range(n + 1):
        y = C[0] @ x
        u, v, state_next = pid_step(state, theta_ref, y, gains, dt, u_limit, anti_windup, Tt, b, d_on_measurement)
        theta[k], u_hist[k], v_hist[k], i_hist[k] = y, u, v, state[0]
        state = state_next
        w = disturbance if t[k] >= t_disturbance else 0.0
        x = Ad @ x + Bd[:, 0] * u + Bd[:, 1] * w
    return {"t": t, "theta": theta, "u": u_hist, "v": v_hist, "I": i_hist}


# ---------------------------------------------------------------- the example

def main(out_dir: str | Path = ".") -> None:
    """Tune the pitch-attitude PID two ways, check the linear loop against its transfer function, show anti-windup."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    p = pitch_plant()
    num_p, den_p = pitch_plant_tf(p)
    print("representative pitch plant (not a specific aircraft): theta/delta_cmd =",
          f"{num_p[0]:.0f}(s + {num_p[1] / num_p[0]:.2f}) / [s (s^2 + {den_p[1] - p['w_servo']:.2f} s + {p['Z_alpha'] * p['M_q'] - p['M_alpha']:.2f})(s + {p['w_servo']:.0f})]")
    sp = np.roots([1.0, -(p["Z_alpha"] + p["M_q"]), p["Z_alpha"] * p["M_q"] - p["M_alpha"]])
    wn = np.sqrt(p["Z_alpha"] * p["M_q"] - p["M_alpha"])
    print(f"short-period poles {sp[0]:.3f}, {sp[1]:.3f}: w_n = {wn:.3f} rad/s, zeta = {-(p['Z_alpha'] + p['M_q']) / (2 * wn):.3f}")

    Ku, wu, Pu = ultimate_gain(num_p, den_p)
    print(f"\nZiegler-Nichols test: Ku = {Ku:.3f} rad/rad, wu = {wu:.3f} rad/s, Pu = {Pu:.3f} s")
    Tf = 0.02
    zn = ziegler_nichols(Ku, Pu)
    designs = {"ZN PID": (*zn["PID"], Tf)}
    poles = [-4 + 4j, -4 - 4j, -1.5]
    Kp, Ki, Kd = pid_place_three_poles(num_p, den_p, poles, Tf)
    designs["placed"] = (Kp, Ki, Kd, Tf)
    for name, g in designs.items():
        nc, dc = closed_loop_tf(num_p, den_p, g)
        print(f"{name:7s}: Kp = {g[0]:.3f}, Ki = {g[1]:.3f} 1/s, Kd = {g[2]:.3f} s, Tf = {g[3]} s; closed-loop poles",
              np.array2string(np.roots(dc), precision=2, suppress_small=True, max_line_width=200))

    print("\nlinear loop (no limit), 10 deg step, dt = 1 ms: simulation vs transfer function")
    dt, t_end, ref = 1e-3, 6.0, np.radians(10.0)
    results = {}
    for name, g in designs.items():
        for bw in (1.0, 0.0):
            sim = simulate_pitch_loop(p, g, ref, t_end, dt, b=bw)
            nc, dc = closed_loop_tf(num_p, den_p, g, b=bw)
            y_tf = ref * tf_step(nc, dc, sim["t"])
            m = step_metrics(sim["t"], sim["theta"], ref)
            results[(name, bw)] = sim
            print(f"  {name:7s} b = {bw:.0f}: overshoot {m['overshoot_pct']:5.1f} %, rise {m['rise_time']:.3f} s, settle {m['settling_time']:.2f} s,"
                  f" peak elevator {np.degrees(np.abs(sim['u']).max()):6.1f} deg; max |sim - TF| = {np.degrees(np.abs(sim['theta'] - y_tf).max()):.4f} deg")

    g = designs["placed"]
    nc, dc = closed_loop_tf(num_p, den_p, g, b=1.0)
    for dti in (4e-3, 2e-3, 1e-3):
        sim = simulate_pitch_loop(p, g, ref, t_end, dti)
        err = np.abs(sim["theta"] - ref * tf_step(nc, dc, sim["t"])).max()
        print(f"  convergence: dt = {dti * 1e3:.0f} ms, max |sim - TF| = {np.degrees(err):.4f} deg")

    lim = p["delta_max"]
    big = np.radians(20.0)
    print(f"\n20 deg pitch step with the elevator limited to +-{np.degrees(lim):.0f} deg (placed gains, b = 1):")
    aw = {}
    for scheme in ("none", "clamp", "backcalc"):
        sim = simulate_pitch_loop(p, g, big, 8.0, dt, u_limit=lim, anti_windup=scheme)
        m = step_metrics(sim["t"], sim["theta"], big)
        aw[scheme] = sim
        sat_time = np.sum(np.abs(sim["u"]) >= lim - 1e-12) * dt
        print(f"  {scheme:8s}: overshoot {m['overshoot_pct']:5.1f} %, settle {m['settling_time']:.2f} s, time saturated {sat_time:.2f} s,"
              f" peak integrator {np.degrees(np.abs(sim['I']).max()):.1f} deg")
    lin = simulate_pitch_loop(p, g, big, 8.0, dt)
    print(f"  no limit: overshoot {step_metrics(lin['t'], lin['theta'], big)['overshoot_pct']:.1f} %, peak elevator {np.degrees(np.abs(lin['u']).max()):.1f} deg")

    print("\ntrim disturbance of 2 deg equivalent elevator at t = 3 s, holding 5 deg pitch:")
    for name, gg in (("PD only", (g[0], 0.0, g[2], Tf)), ("PID", g)):
        sim = simulate_pitch_loop(p, gg, np.radians(5.0), 15.0, dt, disturbance=np.radians(-2.0), t_disturbance=3.0)
        print(f"  {name:7s}: pitch error at 15 s {np.degrees(np.radians(5.0) - sim['theta'][-1]):.3f} deg")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    for (name, bw), sim in results.items():
        axes[0].plot(sim["t"], np.degrees(sim["theta"]), label=f"{name}, b = {bw:.0f}")
    axes[0].axhline(10.0, color="k", lw=0.6)
    axes[0].set_xlabel("t [s]"); axes[0].set_ylabel("pitch [deg]"); axes[0].legend(fontsize=8); axes[0].grid(lw=0.4)
    axes[0].set_title("linear loop, 10 deg step (simulation)", fontsize=9)
    for scheme, sim in aw.items():
        axes[1].plot(sim["t"], np.degrees(sim["theta"]), label=scheme)
        axes[2].plot(sim["t"], np.degrees(sim["u"]), label=scheme)
    axes[1].plot(lin["t"], np.degrees(lin["theta"]), "k--", lw=0.8, label="no limit")
    axes[1].set_xlabel("t [s]"); axes[1].set_ylabel("pitch [deg]"); axes[1].legend(fontsize=8); axes[1].grid(lw=0.4)
    axes[1].set_title("20 deg step, elevator limit 15 deg", fontsize=9)
    axes[2].set_xlabel("t [s]"); axes[2].set_ylabel("elevator command [deg]"); axes[2].legend(fontsize=8); axes[2].grid(lw=0.4)
    fig.tight_layout()
    fig.savefig(out_dir / "pid_pitch.png", dpi=150)
    plt.close(fig)
    print(f"saved {out_dir / 'pid_pitch.png'}")


if __name__ == "__main__":
    main()
