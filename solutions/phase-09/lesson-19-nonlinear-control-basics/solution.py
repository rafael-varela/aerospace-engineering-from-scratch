"""P09.19 Nonlinear Control Basics: reference solution.

Three tools for loops that are not linear, each checked against simulation:

    * Lyapunov. The Lyapunov equation A^T P + P A = -Q for a linear(ised) loop (solved from scratch by the
      Kronecker form, checked against ``scipy.linalg.solve_continuous_lyapunov``); V = x^T P x as a
      Lyapunov function of the *nonlinear* loop; and a region-of-attraction estimate: the largest level set
      {V < c} on which dV/dt < 0, for an inverted pendulum stabilised by a torque-limited controller.
    * Describing functions. The first-harmonic gain N(A) of a saturation and of a relay (with hysteresis),
      and the harmonic-balance prediction of a limit cycle, G(jw) N(A) = -1, for a fin servo under relay
      feedback and under a saturating high-gain amplifier, compared with the amplitude and frequency
      measured in a simulation.
    * Feedback linearisation. A control law that cancels the pendulum's gravity torque so the closed loop is
      exactly linear, checked against the linear closed-form response at large angles.

Conventions:
    * Pendulum state x = (theta, omega): angle from upright [rad], rate [rad/s]; input u = torque [N m].
    * Describing functions are for a nonlinearity driven by e(t) = A sin(w t); N(A) is complex in general.
    * Transfer functions are (num, den) coefficient arrays in descending powers of s.

Run:  python solution.py          (prints the three studies, writes nonlinear.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-09/lesson-19-nonlinear-control-basics
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


# ---------------------------------------------------------------- given: models

def pendulum_params() -> dict[str, float]:
    """Return representative parameters of an inverted-pendulum test stand (a vehicle balanced on its pivot):
    m = 2 kg, l = 0.5 m (pivot to centre of mass), b = 0.1 N m s pivot damping, g = 9.81 m/s^2, torque limit
    6 N m (below m g l = 9.81 N m, so large angles cannot be held)."""
    return {"m": 2.0, "l": 0.5, "b": 0.1, "g": 9.81, "u_max": 6.0}


def pendulum_rhs(x: np.ndarray, u: float, p: dict[str, float]) -> np.ndarray:
    """Return dx/dt of m l^2 theta'' = m g l sin(theta) - b theta' + u (theta measured from upright)."""
    ml2 = p["m"] * p["l"] ** 2
    return np.array([x[1], (p["m"] * p["g"] * p["l"] * math.sin(x[0]) - p["b"] * x[1] + u) / ml2])


def pendulum_gains(p: dict[str, float], wn: float = 4.0, zeta: float = 0.7) -> np.ndarray:
    """Return K = [k1, k2] placing the linearised closed loop at s^2 + 2 zeta wn s + wn^2 for u = -K x:
    k1 = m l^2 (wn^2 + g/l), k2 = 2 zeta wn m l^2 - b."""
    ml2 = p["m"] * p["l"] ** 2
    return np.array([ml2 * (wn ** 2 + p["g"] / p["l"]), 2.0 * zeta * wn * ml2 - p["b"]])


def saturated_feedback_rhs(p: dict[str, float], K: np.ndarray):
    """Return f(x) for the pendulum with u = sat(-K x, u_max)."""
    k1, k2, um = float(K[0]), float(K[1]), p["u_max"]

    def f(x: np.ndarray) -> np.ndarray:
        u = min(um, max(-um, -(k1 * x[0] + k2 * x[1])))
        return pendulum_rhs(x, u, p)
    return f


def servo_loop_plant() -> tuple[np.ndarray, np.ndarray]:
    """Return G(s) = 30 / (s (0.04 s + 1)(0.01 s + 1)): the fin servo of P09.18 with a 10 ms amplifier lag,
    shaft angle [rad] per volt. Phase crossover at 1/sqrt(0.04 * 0.01) = 50 rad/s, |G| there 0.24 rad/V."""
    return np.array([30.0]), np.polymul(np.polymul([1.0, 0.0], [0.04, 1.0]), [0.01, 1.0])


def tf_to_ss(num: np.ndarray, den: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (A, B, C) of a strictly proper SISO transfer function in controllable canonical form."""
    den = np.asarray(den, dtype=float)
    num = np.asarray(num, dtype=float) / den[0]
    den = den / den[0]
    n = len(den) - 1
    A = np.zeros((n, n))
    A[0, :] = -den[1:]
    A[1:, :-1] = np.eye(n - 1)
    B = np.zeros(n)
    B[0] = 1.0
    C = np.concatenate([np.zeros(n - len(num)), num])
    return A, B, C


def simulate_relay_loop(num: np.ndarray, den: np.ndarray, M: float, hysteresis: float = 0.0,
                        t_end: float = 3.0, y0: float = 0.05,
                        max_switches: int = 5000) -> tuple[np.ndarray, np.ndarray]:
    """Simulate e = -y into a relay (output +-M, switching at e = +-hysteresis) driving G; return (t, e).

    The plant is linear between switches, so each segment is integrated by ``solve_ivp`` (RK45, tight
    tolerances) and stopped exactly at the switching event; nothing is smoothed. Raise ``RuntimeError`` after
    ``max_switches`` switches: an ideal relay can chatter (switch infinitely often in finite time).
    """
    A, B, C = tf_to_ss(num, den)
    x = np.zeros(len(B))
    x[-1] = y0 / C[-1]
    u = -M if C @ x > 0 else M
    ts, es = [0.0], [-(C @ x)]
    t = 0.0
    n_sw = 0
    while t < t_end:
        def ev(tt, xx, u=u):
            e = -(C @ xx)
            return e - hysteresis if u < 0 else e + hysteresis    # relay at -M flips when e rises to +h, etc.
        ev.terminal = True
        ev.direction = 1.0 if u < 0 else -1.0
        sol = solve_ivp(lambda tt, xx: A @ xx + B * u, (t, t_end), x, events=ev, rtol=1e-10, atol=1e-12,
                        max_step=2e-3, dense_output=False)
        ts.extend(sol.t[1:])
        es.extend(-(C @ sol.y[:, 1:]))
        t = sol.t[-1]
        x = sol.y[:, -1]
        if sol.status == 1:
            u = -u
            t += 1e-12
            n_sw += 1
            if n_sw > max_switches:
                raise RuntimeError(f"relay chattering: more than {max_switches} switches by t = {t:.4f} s")
    return np.array(ts), np.array(es)


def simulate_saturation_loop(num: np.ndarray, den: np.ndarray, gain: float, limit: float, t_end: float = 4.0,
                             y0: float = 0.01) -> tuple[np.ndarray, np.ndarray]:
    """Simulate e = -y, u = sat(gain e, limit) driving G with ``solve_ivp``; return (t, e) on a 0.5 ms grid."""
    A, B, C = tf_to_ss(num, den)
    x0 = np.zeros(len(B))
    x0[-1] = y0 / C[-1]

    def f(t, x):
        return A @ x + B * float(np.clip(-gain * (C @ x), -limit, limit))

    t_eval = np.arange(0.0, t_end, 5e-4)
    sol = solve_ivp(f, (0.0, t_end), x0, t_eval=t_eval, rtol=1e-9, atol=1e-12, max_step=1e-3)
    return sol.t, -(C @ sol.y)


# ---------------------------------------------------------------- Lyapunov

def lyapunov_solve(A: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Return P solving A^T P + P A + Q = 0 by the Kronecker form (I kron A^T + A^T kron I) vec P = -vec Q.

    Raise ``ValueError`` if the Kronecker matrix is singular (two eigenvalues of A sum to zero). The result
    is symmetrised when Q is symmetric. For Hurwitz A and Q > 0, P > 0 and V = x^T P x is a Lyapunov function.
    """
    A = np.asarray(A, dtype=float)
    Q = np.asarray(Q, dtype=float)
    n = A.shape[0]
    I = np.eye(n)
    Mk = np.kron(I, A.T) + np.kron(A.T, I)
    if np.linalg.cond(Mk) > 1e14:
        raise ValueError("singular Lyapunov equation")
    P = np.linalg.solve(Mk, -Q.reshape(-1, order="F")).reshape((n, n), order="F")
    return 0.5 * (P + P.T) if np.allclose(Q, Q.T) else P


def jacobian_fd(f, x0: np.ndarray, h: float = 1e-6) -> np.ndarray:
    """Return the Jacobian df/dx at x0 by central differences with step h (relative to max(1, |x0_i|))."""
    x0 = np.asarray(x0, dtype=float)
    f0 = np.asarray(f(x0))
    J = np.zeros((len(f0), len(x0)))
    for i in range(len(x0)):
        d = np.zeros_like(x0)
        d[i] = h * max(1.0, abs(x0[i]))
        J[:, i] = (np.asarray(f(x0 + d)) - np.asarray(f(x0 - d))) / (2 * d[i])
    return J


def lyapunov_derivative(f, P: np.ndarray, x: np.ndarray) -> float:
    """Return dV/dt = 2 x^T P f(x) for V = x^T P x along dx/dt = f(x)."""
    x = np.asarray(x, dtype=float)
    return float(2.0 * x @ P @ np.asarray(f(x)))


def roa_estimate(f, P: np.ndarray, n_dir: int = 360, n_r: int = 400, r_max: float = 3.0) -> float:
    """Return c* such that dV/dt < 0 on {0 < V(x) < c*}: an estimate of the region of attraction.

    Write x = r P^{-1/2} (cos phi, sin phi), so V = r^2 (2-state systems). For each of n_dir directions scan
    r on a uniform grid of n_r points in (0, r_max] and record the first r at which dV/dt >= 0; c* is the
    square of the smallest such r (the grid point before it, to stay inside), or r_max^2 if none is found.
    The set {V < c*} is then positively invariant and every trajectory starting in it converges.
    """
    P = np.asarray(P, dtype=float)
    lam, U = np.linalg.eigh(P)
    Pmh = U @ np.diag(1.0 / np.sqrt(lam)) @ U.T          # P^{-1/2}
    rs = np.linspace(r_max / n_r, r_max, n_r)
    r_first = r_max
    for phi in np.linspace(0.0, 2 * np.pi, n_dir, endpoint=False):
        d = Pmh @ np.array([np.cos(phi), np.sin(phi)])
        for j, r in enumerate(rs):
            if r >= r_first:
                break
            if lyapunov_derivative(f, P, r * d) >= 0.0:
                r_first = rs[j - 1] if j > 0 else 0.0
                break
    return float(r_first ** 2)


# ---------------------------------------------------------------- describing functions

def df_saturation(A: float | np.ndarray, delta: float, k: float = 1.0) -> np.ndarray:
    """Return the describing function of a saturation with slope k and linear range |e| <= delta.

    N(A) = k for A <= delta; N(A) = (2k/pi) [asin(delta/A) + (delta/A) sqrt(1 - (delta/A)^2)] for A > delta.
    Real (a memoryless odd nonlinearity has no phase shift).
    """
    A = np.atleast_1d(np.asarray(A, dtype=float))
    out = np.full_like(A, k, dtype=float)
    big = A > delta
    r = delta / A[big]
    out[big] = 2.0 * k / np.pi * (np.arcsin(r) + r * np.sqrt(1.0 - r * r))
    return out if out.size > 1 else out[:1]


def df_relay(A: float | np.ndarray, M: float, hysteresis: float = 0.0) -> np.ndarray:
    """Return the describing function of a relay +-M with hysteresis h (switches at e = +-h), A > h:

    N(A) = (4M/(pi A)) [sqrt(1 - (h/A)^2) - j h/A]; for h = 0, N = 4M/(pi A). Raise ``ValueError`` if A <= h.
    """
    A = np.atleast_1d(np.asarray(A, dtype=float))
    if np.any(A <= hysteresis):
        raise ValueError("amplitude must exceed the hysteresis")
    r = hysteresis / A
    return 4.0 * M / (np.pi * A) * (np.sqrt(1.0 - r * r) - 1j * r)


def first_harmonic_gain(nonlinearity, A: float, n: int = 4096) -> complex:
    """Return N(A) = (b1 + j a1)/A from the Fourier coefficients of y(t) = nonlinearity(A sin(theta)) over one
    period, computed numerically (rectangle rule over n points, which is spectrally accurate for periodic
    integrands). ``nonlinearity`` maps an array of theta values to outputs (to allow memory such as hysteresis)."""
    th = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    y = nonlinearity(th, A)
    b1 = 2.0 / n * np.sum(y * np.sin(th))
    a1 = 2.0 / n * np.sum(y * np.cos(th))
    return complex(b1, a1) / A


def phase_crossover(num: np.ndarray, den: np.ndarray, phase_deg: float = -180.0,
                    w_lo: float = 1e-3, w_hi: float = 1e5) -> float:
    """Return the frequency where the unwrapped phase of G(jw) equals ``phase_deg``, by bisection in log w
    after bracketing on a 4000-point log grid. Raise ``ValueError`` if there is no crossing."""
    w = np.logspace(np.log10(w_lo), np.log10(w_hi), 4000)
    ph = np.degrees(np.unwrap(np.angle(np.polyval(num, 1j * w) / np.polyval(den, 1j * w))))
    idx = np.where(np.sign(ph[:-1] - phase_deg) != np.sign(ph[1:] - phase_deg))[0]
    if len(idx) == 0:
        raise ValueError("no phase crossing")
    a, b = np.log(w[idx[0]]), np.log(w[idx[0] + 1])
    ref = ph[idx[0]]

    def phase(lw: float) -> float:
        g = np.polyval(num, 1j * np.exp(lw)) / np.polyval(den, 1j * np.exp(lw))
        p = np.degrees(np.angle(g))
        return p + 360.0 * np.round((ref - p) / 360.0)     # same branch as the grid

    fa = phase(a) - phase_deg
    for _ in range(100):
        m = 0.5 * (a + b)
        fm = phase(m) - phase_deg
        if np.sign(fm) == np.sign(fa):
            a, fa = m, fm
        else:
            b = m
    return float(np.exp(0.5 * (a + b)))


def predict_relay_limit_cycle(num: np.ndarray, den: np.ndarray, M: float) -> tuple[float, float]:
    """Return (A, w) of the limit cycle predicted for an ideal relay +-M around G: G(jw) = -1/N(A) is real and
    negative, so w is the phase crossover and A = 4 M |G(jw)| / pi."""
    w = phase_crossover(num, den)
    g = abs(np.polyval(num, 1j * w) / np.polyval(den, 1j * w))
    return 4.0 * M * g / np.pi, w


def predict_relay_hysteresis_limit_cycle(num: np.ndarray, den: np.ndarray, M: float, hysteresis: float,
                                        w_lo: float = 1e-2, w_hi: float = 1e4) -> tuple[float, float]:
    """Return (A, w) predicted for a relay +-M with hysteresis h around G.

    -1/N(A) = -(pi/(4M)) (sqrt(A^2 - h^2) + j h): a horizontal line at Im = -pi h/(4M). So w solves
    Im G(jw) = -pi h/(4M) with Re G(jw) < 0 (bisection in log w on the highest-frequency bracket of a 4000-point
    grid, where |G| is falling), and then A = sqrt(h^2 + (4M Re G(jw)/pi)^2). Raise ``ValueError`` if none.
    """
    target = -np.pi * hysteresis / (4.0 * M)
    w = np.logspace(np.log10(w_lo), np.log10(w_hi), 4000)
    G = np.polyval(num, 1j * w) / np.polyval(den, 1j * w)
    d = G.imag - target
    idx = [i for i in np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0] if G.real[i] < 0]
    if not idx:
        raise ValueError("no intersection with -1/N")
    a, b = np.log(w[idx[-1]]), np.log(w[idx[-1] + 1])

    def im(lw: float) -> float:
        return float((np.polyval(num, 1j * np.exp(lw)) / np.polyval(den, 1j * np.exp(lw))).imag - target)

    fa = im(a)
    for _ in range(100):
        m = 0.5 * (a + b)
        fm = im(m)
        if np.sign(fm) == np.sign(fa):
            a, fa = m, fm
        else:
            b = m
    wc = float(np.exp(0.5 * (a + b)))
    re = float((np.polyval(num, 1j * wc) / np.polyval(den, 1j * wc)).real)
    return float(np.sqrt(hysteresis ** 2 + (4.0 * M * re / np.pi) ** 2)), wc


def predict_saturation_limit_cycle(num: np.ndarray, den: np.ndarray, gain: float, limit: float) -> tuple[float, float]:
    """Return (A, w) of the limit cycle predicted for a saturating amplifier u = sat(gain e, limit) around G.

    The amplifier is a saturation of slope ``gain`` and linear range delta = limit/gain. At the phase
    crossover w, harmonic balance needs N(A) |G(jw)| = 1; N falls from ``gain`` to 0 as A grows, so a solution
    exists only if gain |G(jw)| > 1 (the linear loop is unstable); A is found by bisection.
    Raise ``ValueError`` otherwise.
    """
    w = phase_crossover(num, den)
    g = abs(np.polyval(num, 1j * w) / np.polyval(den, 1j * w))
    delta = limit / gain
    if gain * g <= 1.0:
        raise ValueError("linear loop is stable: no limit cycle predicted")
    lo, hi = delta, delta * 2.0
    while df_saturation(hi, delta, gain)[0] * g > 1.0:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if df_saturation(mid, delta, gain)[0] * g > 1.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), w


def oscillation_amplitude_frequency(t: np.ndarray, y: np.ndarray, t_settle: float) -> tuple[float, float]:
    """Return (amplitude, angular frequency) of a steady oscillation from samples after t_settle:
    amplitude = (max - min)/2; frequency = 2 pi / period, the period being the mean spacing of the upward
    zero crossings of y - mean(y) (located by linear interpolation)."""
    m = t >= t_settle
    tt, yy = t[m], y[m]
    amp = 0.5 * (yy.max() - yy.min())
    z = yy - 0.5 * (yy.max() + yy.min())
    idx = np.where((z[:-1] < 0) & (z[1:] >= 0))[0]
    tc = tt[idx] - z[idx] * (tt[idx + 1] - tt[idx]) / (z[idx + 1] - z[idx])
    if len(tc) < 2:
        raise ValueError("fewer than two crossings: not oscillating")
    return float(amp), float(2 * np.pi / np.mean(np.diff(tc)))


def converges(f, x0: np.ndarray, t_end: float = 8.0, tol: float = 1e-3) -> bool:
    """Return True if the trajectory of dx/dt = f(x) from x0 ends within ``tol`` of the origin at t_end."""
    sol = solve_ivp(lambda t, x: f(x), (0.0, t_end), np.asarray(x0, dtype=float), rtol=1e-8, atol=1e-10, max_step=0.05)
    return bool(np.linalg.norm(sol.y[:, -1]) < tol)


def largest_recoverable_angle(f, hi: float = np.radians(60.0), iters: int = 12) -> float:
    """Return (by bisection on theta0, starting at rest) the largest initial angle from which dx/dt = f(x)
    converges; assumes recovery is monotone in theta0 and fails at ``hi``."""
    lo = 0.0
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if converges(f, np.array([mid, 0.0]), t_end=10.0):
            lo = mid
        else:
            hi = mid
    return lo


# ---------------------------------------------------------------- feedback linearisation

def feedback_linearising_torque(x: np.ndarray, v: float, p: dict[str, float]) -> float:
    """Return u that makes theta'' = v exactly: u = m l^2 v - m g l sin(theta) + b theta'."""
    return float(p["m"] * p["l"] ** 2 * v - p["m"] * p["g"] * p["l"] * np.sin(x[0]) + p["b"] * x[1])


def simulate_pendulum(x0: np.ndarray, controller, p: dict[str, float], t_end: float,
                      n_out: int = 401) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (t, x, u) for the pendulum under u = controller(x) with ``solve_ivp`` (rtol 1e-10)."""
    t_eval = np.linspace(0.0, t_end, n_out)
    sol = solve_ivp(lambda t, x: pendulum_rhs(x, controller(x), p), (0.0, t_end), np.asarray(x0, dtype=float),
                    t_eval=t_eval, rtol=1e-10, atol=1e-12, max_step=0.01)
    u = np.array([controller(xx) for xx in sol.y.T])
    return sol.t, sol.y.T, u


# ---------------------------------------------------------------- the example

def main(out_dir: str | Path = ".") -> None:
    """Print the Lyapunov, describing-function and feedback-linearisation studies; save nonlinear.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    p = pendulum_params()
    K = pendulum_gains(p)
    ml2 = p["m"] * p["l"] ** 2
    print(f"pendulum: m = {p['m']} kg, l = {p['l']} m, b = {p['b']} N m s; m g l = {p['m'] * p['g'] * p['l']:.2f} N m,"
          f" torque limit {p['u_max']} N m -> largest holdable angle {np.degrees(np.arcsin(p['u_max'] / (p['m'] * p['g'] * p['l']))):.1f} deg")
    print(f"gains for wn = 4 rad/s, zeta = 0.7: k1 = {K[0]:.3f} N m/rad, k2 = {K[1]:.3f} N m s/rad")
    f_open = lambda x: pendulum_rhs(x, 0.0, p)
    J0 = jacobian_fd(f_open, np.zeros(2))
    print(f"open-loop Jacobian at upright {np.round(J0, 4).tolist()}, eigenvalues {np.round(np.linalg.eigvals(J0), 3)} (saddle: unstable)")
    f_cl = saturated_feedback_rhs(p, K)
    Acl = jacobian_fd(f_cl, np.zeros(2))
    print(f"closed-loop Jacobian eigenvalues {np.round(np.linalg.eigvals(Acl), 3)}")
    Q = np.diag([1.0, 0.1])
    P = lyapunov_solve(Acl, Q)
    print(f"Lyapunov P for Q = diag(1, 0.1): [[{P[0, 0]:.4f}, {P[0, 1]:.4f}], [{P[1, 0]:.4f}, {P[1, 1]:.4f}]],"
          f" residual {np.abs(Acl.T @ P + P @ Acl + Q).max():.1e}")
    c = roa_estimate(f_cl, P)
    th_edge = np.sqrt(c / (P[0, 0] - P[0, 1] ** 2 / P[1, 1]))   # widest theta on the ellipse
    print(f"region-of-attraction estimate: V < {c:.4f}, theta extent +-{np.degrees(th_edge):.1f} deg;"
          f" theta at rest inside it: +-{np.degrees(np.sqrt(c / P[0, 0])):.1f} deg")
    # check by simulation: samples inside converge; find how far beyond along omega = 0 still converges
    rng = np.random.default_rng(19)
    inside = []
    Lc = np.linalg.cholesky(np.linalg.inv(P))
    for _ in range(100):
        z = rng.standard_normal(2)
        z *= np.sqrt(c) * rng.uniform(0.0, 1.0) ** 0.5 / np.linalg.norm(z)
        inside.append(Lc @ z)
    conv = sum(converges(f_cl, x0) for x0 in inside)
    print(f"   {len(inside)} random starts inside the estimate (seed 19): {conv} converge")
    best = largest_recoverable_angle(f_cl)
    print(f"   from rest, simulation recovers up to {np.degrees(best):.2f} deg (bisection): the estimate is conservative")

    # describing functions
    nG, dG = servo_loop_plant()
    A_r, w_r = predict_relay_limit_cycle(nG, dG, 1.0)
    t, e = simulate_relay_loop(nG, dG, 1.0, t_end=3.0)
    A_s, w_s = oscillation_amplitude_frequency(t, e, 2.0)
    print(f"\nrelay +-1 V on G = 30/(s(0.04s+1)(0.01s+1)): predicted A = {A_r:.4f} rad at {w_r:.2f} rad/s;"
          f" simulated A = {A_s:.4f} rad at {w_s:.2f} rad/s (errors {100 * (A_r / A_s - 1):+.1f} %, {100 * (w_r / w_s - 1):+.1f} %)")
    h = 0.02
    t_h, e_h = simulate_relay_loop(nG, dG, 1.0, hysteresis=h, t_end=3.0)
    A_sh, w_sh = oscillation_amplitude_frequency(t_h, e_h, 2.0)
    A_hp, w_hp = predict_relay_hysteresis_limit_cycle(nG, dG, 1.0, h)
    print(f"relay with hysteresis {h} rad: predicted A = {A_hp:.4f} rad at {w_hp:.2f} rad/s; simulated A = {A_sh:.4f} at {w_sh:.2f} rad/s")
    k_crit = 1.0 / abs(np.polyval(nG, 1j * w_r) / np.polyval(dG, 1j * w_r))
    gain, limit = 2.0 * k_crit, 1.0
    A_p, w_p = predict_saturation_limit_cycle(nG, dG, gain, limit)
    t2, e2 = simulate_saturation_loop(nG, dG, gain, limit, t_end=4.0)
    A_m, w_m = oscillation_amplitude_frequency(t2, e2, 3.0)
    print(f"saturating amplifier, gain {gain:.3f} V/rad (= 2 x critical {k_crit:.3f}), limit {limit} V:"
          f" predicted A = {A_p:.4f} rad at {w_p:.2f} rad/s; simulated A = {A_m:.4f} rad at {w_m:.2f} rad/s")

    # feedback linearisation from 60 deg, without a torque limit, against the linear controller
    wn, zeta = 4.0, 0.7
    v_of = lambda x: -(wn ** 2) * x[0] - 2 * zeta * wn * x[1]
    fl = lambda x: feedback_linearising_torque(x, v_of(x), p)
    lin = lambda x: float(-K @ x)
    x0 = np.array([np.radians(60.0), 0.0])
    t3, x_fl, u_fl = simulate_pendulum(x0, fl, p, 4.0)
    t3, x_ln, u_ln = simulate_pendulum(x0, lin, p, 4.0)
    Alin = np.array([[0.0, 1.0], [-wn ** 2, -2 * zeta * wn]])
    lam, V = np.linalg.eig(Alin)
    x_exact = np.real(np.array([V @ np.diag(np.exp(lam * ti)) @ np.linalg.solve(V, x0) for ti in t3]))
    print(f"\nfeedback linearisation from 60 deg (no torque limit): max |theta - linear closed form| = {np.abs(x_fl[:, 0] - x_exact[:, 0]).max():.1e} rad;"
          f" peak torque {np.abs(u_fl).max():.2f} N m")
    print(f"   linear controller from 60 deg (no limit): overshoot {np.degrees(-x_ln[:, 0].min()):.2f} deg vs FL {np.degrees(-x_fl[:, 0].min()):.2f} deg;"
          f" peak torque {np.abs(u_ln).max():.2f} N m")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    th = np.linspace(0, 2 * np.pi, 400)
    ell = np.array([Lc @ (np.sqrt(c) * np.array([np.cos(a), np.sin(a)])) for a in th])
    TH, OM = np.meshgrid(np.linspace(-1.2, 1.2, 25), np.linspace(-6, 6, 25))
    dTH = np.zeros_like(TH); dOM = np.zeros_like(OM)
    for i in range(TH.shape[0]):
        for j in range(TH.shape[1]):
            d = f_cl(np.array([TH[i, j], OM[i, j]]))
            dTH[i, j], dOM[i, j] = d
    axes[0].streamplot(TH, OM, dTH, dOM, density=0.9, color="0.7", linewidth=0.6)
    axes[0].plot(ell[:, 0], ell[:, 1], "C3", lw=1.5, label="V < c* (estimate)")
    axes[0].set_xlabel("theta [rad]"); axes[0].set_ylabel("omega [rad/s]"); axes[0].legend(fontsize=8)
    axes[0].set_title("torque-limited pendulum", fontsize=9)
    axes[1].plot(t, e, lw=0.8, label="relay, simulated")
    axes[1].axhline(A_r, color="C3", ls="--", lw=0.8, label="DF prediction")
    axes[1].axhline(-A_r, color="C3", ls="--", lw=0.8)
    axes[1].plot(t2, e2, lw=0.8, label="saturation, simulated")
    axes[1].set_xlim(1.5, 3.0); axes[1].set_xlabel("t [s]"); axes[1].set_ylabel("e [rad]"); axes[1].legend(fontsize=7)
    axes[2].plot(t3, np.degrees(x_fl[:, 0]), label="feedback linearised")
    axes[2].plot(t3, np.degrees(x_exact[:, 0]), "k:", label="linear closed form")
    axes[2].plot(t3, np.degrees(x_ln[:, 0]), label="linear controller")
    axes[2].set_xlabel("t [s]"); axes[2].set_ylabel("theta [deg]"); axes[2].legend(fontsize=8)
    for ax in axes:
        ax.grid(lw=0.4)
    fig.tight_layout()
    fig.savefig(out_dir / "nonlinear.png", dpi=150)
    plt.close(fig)
    print(f"saved {out_dir / 'nonlinear.png'}")


if __name__ == "__main__":
    main()
