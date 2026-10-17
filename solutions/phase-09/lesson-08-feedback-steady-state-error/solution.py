"""P09.08 Feedback & Steady-State Error: reference solution.

An airspeed-hold (auto-throttle) loop on a light aircraft in level flight. The throttle changes the
propeller thrust; thrust minus drag accelerates the aircraft along its path. Linearised about a trim
airspeed V0, with altitude held by another loop and the thrust treated as independent of airspeed:

    m dv/dt = dT - c v + d,      c = dD/dV = 2 D0 / V0      (v, dT, d: deviations from trim)

with a first-order throttle servo (tau_a) and a first-order engine-and-propeller thrust lag (tau_e)
between the commanded thrust u and dT, and an input disturbance d (a drag change, for example from
icing or from a gear or flap extension, or the weight component m g sin(gamma) in a shallow climb):

    G(s) = v / u = 1 / ((tau_a s + 1)(tau_e s + 1)(m s + c)).

This module computes, from polynomial transfer functions:

* the loop L = C G, the closed loop T = L / (1 + L), the sensitivity S = 1 / (1 + L) (S + T = 1),
  and the disturbance path G / (1 + L);
* the system type and the error constants K_p = lim L, K_v = lim s L, K_a = lim s^2 L;
* the steady-state error to a step, a ramp and a parabola from them;
* the final-value theorem lim_{s->0} s Y(s), refused when s Y(s) has a pole in the closed right half
  plane (where the theorem does not apply);
* the largest stable proportional gain (by bisection on the closed-loop roots);
* a from-scratch time simulation (companion form, exact zero-order-hold steps) used as the
  independent check: the error of a long simulation must settle to the final-value prediction.

Numbers are representative of a light single-engine aircraft (stated as typical, not one type's data).
Units: SI. Airspeed in m/s; thrust and drag in N; gains of the controller in N per (m/s).

Run:  python solution.py [out_dir]     (prints the lesson's tables, writes cruise_control.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-09/lesson-08-feedback-steady-state-error
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# representative light aircraft at 55 m/s (about 107 kn) in level flight
MASS = 1100.0         # kg
V0 = 55.0             # m/s
D0 = 1080.0           # N, drag at trim (lift-to-drag ratio about 10)
TAU_A = 0.3           # s, throttle servo
TAU_E = 1.0           # s, engine and propeller thrust response


# ---------------------------------------------------------------- given helpers (polynomials highest power first)

def poly_add(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Sum of two polynomials of possibly different length (highest power first)."""
    p, q = np.atleast_1d(np.asarray(p, dtype=float)), np.atleast_1d(np.asarray(q, dtype=float))
    n = max(p.size, q.size)
    return np.concatenate([np.zeros(n - p.size), p]) + np.concatenate([np.zeros(n - q.size), q])


def poly_trim(p: np.ndarray) -> np.ndarray:
    """Drop leading (highest-power) zero coefficients, keeping at least one coefficient."""
    p = np.atleast_1d(np.asarray(p, dtype=float))
    nz = np.flatnonzero(p != 0.0)
    return p[nz[0]:] if nz.size else np.zeros(1)


def tf_eval(num: np.ndarray, den: np.ndarray, s: np.ndarray | complex) -> np.ndarray:
    """Evaluate num(s)/den(s) at complex points s."""
    return np.polyval(num, s) / np.polyval(den, s)


def cruise_plant(mass: float = MASS, v0: float = V0, d0: float = D0, tau_a: float = TAU_A,
                 tau_e: float = TAU_E) -> tuple[np.ndarray, np.ndarray]:
    """Return (num, den) of G(s) = v/u = 1 / ((tau_a s + 1)(tau_e s + 1)(m s + c)), c = 2 D0 / V0."""
    c = 2.0 * d0 / v0
    den = np.convolve(np.convolve([tau_a, 1.0], [tau_e, 1.0]), [mass, c])
    return np.array([1.0]), den


def tf_to_companion(b: np.ndarray, a: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """(A, B, C, D) of the controllable companion form of a proper b(s)/a(s)."""
    b, a = poly_trim(b), poly_trim(a)
    n = a.size - 1
    if b.size > a.size:
        raise ValueError("improper transfer function")
    b = np.concatenate([np.zeros(a.size - b.size), b]) / a[0]
    a = a / a[0]
    D = float(b[0])
    A = np.zeros((n, n))
    if n:
        A[0, :] = -a[1:]
        A[1:, :-1] = np.eye(n - 1)
    B = np.zeros((n, 1))
    if n:
        B[0, 0] = 1.0
    return A, B, (b[1:] - D * a[1:]).reshape(1, n), D


def expm_series(M: np.ndarray, n_terms: int = 20) -> np.ndarray:
    """Matrix exponential by scaling and squaring of a Taylor series (P02.05)."""
    norm = np.linalg.norm(M, 1)
    s = int(np.ceil(np.log2(norm / 0.5))) if norm > 0.5 else 0
    X = M / 2.0 ** s
    term, E = np.eye(len(M)), np.eye(len(M))
    for k in range(1, n_terms):
        term = term @ X / k
        E = E + term
    for _ in range(s):
        E = E @ E
    return E


def simulate_tf(num: np.ndarray, den: np.ndarray, t: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Response of num(s)/den(s) from rest to u on the uniform grid t (input held over each step; exact
    zero-order-hold steps from the augmented matrix exponential)."""
    A, B, C, D = tf_to_companion(num, den)
    n = A.shape[0]
    u = np.asarray(u, dtype=float)
    if n == 0:
        return D * u
    dt = t[1] - t[0]
    M = np.zeros((n + 1, n + 1))
    M[:n, :n], M[:n, n:] = A, B
    E = expm_series(M * dt)
    Ad, bd = E[:n, :n], E[:n, n]
    X = np.empty((t.size, n))
    x = np.zeros(n)
    for k in range(t.size):
        X[k] = x
        x = Ad @ x + bd * u[k]
    return X @ C[0] + D * u


# ---------------------------------------------------------------- loop algebra

def open_loop(num_c: np.ndarray, den_c: np.ndarray, num_p: np.ndarray, den_p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (num, den) of L = C P (controller in series with the plant)."""
    return np.convolve(num_c, num_p), np.convolve(den_c, den_p)


def closed_loop(num_l: np.ndarray, den_l: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (num, den) of T = L / (1 + L) = num_L / (den_L + num_L) (unity negative feedback)."""
    return poly_trim(num_l), poly_trim(poly_add(den_l, num_l))


def sensitivity(num_l: np.ndarray, den_l: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (num, den) of S = 1 / (1 + L) = den_L / (den_L + num_L): the map from reference to error."""
    return poly_trim(den_l), poly_trim(poly_add(den_l, num_l))


def disturbance_to_output(num_c: np.ndarray, den_c: np.ndarray, num_p: np.ndarray, den_p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (num, den) of Y/D = P / (1 + C P) for a disturbance d added at the plant input.

    = num_P den_C / (den_P den_C + num_P num_C).
    """
    num = np.convolve(num_p, den_c)
    den = poly_add(np.convolve(den_p, den_c), np.convolve(num_p, num_c))
    return poly_trim(num), poly_trim(den)


def is_stable(den: np.ndarray) -> bool:
    """True if every root of den(s) has a strictly negative real part."""
    den = poly_trim(den)
    if den.size == 1:
        return True
    return bool(np.all(np.roots(den).real < 0.0))


def system_type(num_l: np.ndarray, den_l: np.ndarray) -> int:
    """Number of free integrators in L: poles at s = 0 minus zeros at s = 0 (counted as trailing zero coefficients)."""
    def zeros_at_origin(p: np.ndarray) -> int:
        p = poly_trim(p)
        k = 0
        while k < p.size - 1 and p[p.size - 1 - k] == 0.0:
            k += 1
        return k

    return zeros_at_origin(den_l) - zeros_at_origin(num_l)


def error_constants(num_l: np.ndarray, den_l: np.ndarray) -> dict[str, float]:
    """Return {'Kp': lim L, 'Kv': lim s L, 'Ka': lim s^2 L} as s -> 0 (np.inf where the limit is infinite)."""
    n = system_type(num_l, den_l)
    num, den = poly_trim(num_l), poly_trim(den_l)
    # strip the common power of s so that L = s^-n * N0(s)/D0(s) with N0(0), D0(0) nonzero
    while num.size > 1 and num[-1] == 0.0 and den[-1] == 0.0:
        num, den = num[:-1], den[:-1]
    while den.size > 1 and den[-1] == 0.0:
        den = den[:-1]
    while num.size > 1 and num[-1] == 0.0:
        num = num[:-1]
    k0 = num[-1] / den[-1]                       # lim s^n L
    out = {}
    for k, key in enumerate(("Kp", "Kv", "Ka")):
        if k < n:
            out[key] = np.inf
        elif k == n:
            out[key] = float(k0)
        else:
            out[key] = 0.0
    return out


def steady_state_error(num_l: np.ndarray, den_l: np.ndarray, kind: str = "step", amplitude: float = 1.0) -> float:
    """Steady-state error of the unity-feedback loop L to r = A (step), A t (ramp) or A t^2 / 2 (parabola).

    step: A / (1 + Kp); ramp: A / Kv; parabola: A / Ka (np.inf when the constant is zero). Assumes the
    closed loop is stable; raises ValueError if it is not.
    """
    if not is_stable(closed_loop(num_l, den_l)[1]):
        raise ValueError("closed loop is unstable: there is no steady state")
    K = error_constants(num_l, den_l)
    if kind == "step":
        return 0.0 if np.isinf(K["Kp"]) else float(amplitude / (1.0 + K["Kp"]))
    key = {"ramp": "Kv", "parabola": "Ka"}[kind]
    if np.isinf(K[key]):
        return 0.0
    return np.inf if K[key] == 0.0 else float(amplitude / K[key])


def final_value(num: np.ndarray, den: np.ndarray) -> float:
    """Final-value theorem for Y(s) = num/den: lim_{t->inf} y(t) = lim_{s->0} s Y(s).

    Valid only if every pole of s Y(s) has a negative real part; raises ValueError otherwise (an
    unstable or oscillating y has no final value, whatever the formula would say).
    """
    num_s = np.convolve(poly_trim(num), [1.0, 0.0])
    den = poly_trim(den)
    while num_s.size > 1 and den.size > 1 and num_s[-1] == 0.0 and den[-1] == 0.0:
        num_s, den = num_s[:-1], den[:-1]
    if den.size > 1 and not np.all(np.roots(den).real < 0.0):
        raise ValueError("s Y(s) has a pole in the closed right half plane: the final-value theorem does not apply")
    if den[-1] == 0.0:
        raise ValueError("s Y(s) has a pole at the origin")
    return float(num_s[-1] / den[-1])


def critical_gain(num_p: np.ndarray, den_p: np.ndarray, k_hi: float = 1e9, rel_tol: float = 1e-10) -> float:
    """Largest K for which the loop K P is closed-loop stable (bisection on the roots of den + K num).

    Assumes the closed loop is stable for small K > 0. Returns np.inf if it is still stable at k_hi.
    """
    def stable(K: float) -> bool:
        return is_stable(poly_add(den_p, K * np.asarray(num_p, dtype=float)))

    if not stable(1e-12 * k_hi):
        raise ValueError("unstable even for small gain")
    if stable(k_hi):
        return np.inf
    lo, hi = 1e-12 * k_hi, k_hi
    while hi - lo > rel_tol * hi:
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if stable(mid) else (lo, mid)
    return lo


# ---------------------------------------------------------------- the lesson's numbers

def pi_controller(kp: float, ki: float) -> tuple[np.ndarray, np.ndarray]:
    """C(s) = kp + ki / s = (kp s + ki) / s."""
    return np.array([kp, ki]), np.array([1.0, 0.0])


def p_controller(kp: float) -> tuple[np.ndarray, np.ndarray]:
    """C(s) = kp."""
    return np.array([kp]), np.array([1.0])


KP_DESIGN = 400.0      # N per m/s
KI_DESIGN = 40.0       # N per m (per m/s per s)


def fvt_vs_simulation(num_c, den_c, kind: str, t_end: float = 600.0, dt: float = 0.05, amplitude: float = 1.0) -> tuple[float, float]:
    """Return (final-value prediction, simulated error at t_end) for the cruise loop with controller C.

    kind 'step' or 'ramp' (reference r of the given amplitude), or 'disturbance': a drag increase of
    `amplitude` newtons from t = 0 (d = -amplitude at the plant input) with r = 0, so the error r - v
    is the airspeed lost, positive."""
    num_p, den_p = cruise_plant()
    num_l, den_l = open_loop(num_c, den_c, num_p, den_p)
    t = np.arange(0.0, t_end + dt / 2, dt)
    if kind in ("step", "ramp"):
        ns, ds = sensitivity(num_l, den_l)
        if kind == "step":
            pred = final_value(amplitude * ns, np.convolve(ds, [1.0, 0.0]))
            r = amplitude * np.ones_like(t)
        else:
            pred = final_value(amplitude * ns, np.convolve(ds, [1.0, 0.0, 0.0]))
            r = amplitude * t
        nt, dtf = closed_loop(num_l, den_l)
        e = r - simulate_tf(nt, dtf, t, r)
    else:
        nd, dd = disturbance_to_output(num_c, den_c, num_p, den_p)
        d = -amplitude                                   # a drag increase decelerates the aircraft
        pred = -final_value(d * nd, np.convolve(dd, [1.0, 0.0]))
        e = -simulate_tf(nd, dd, t, d * np.ones_like(t))
    return pred, float(e[-1])


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's tables and save cruise_control.png in out_dir."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    num_p, den_p = cruise_plant()
    c = 2 * D0 / V0
    print(f"plant: m = {MASS:.0f} kg, V0 = {V0} m/s, D0 = {D0:.0f} N -> c = 2 D0/V0 = {c:.2f} N/(m/s), "
          f"m/c = {MASS / c:.1f} s, tau_a = {TAU_A} s, tau_e = {TAU_E} s")
    print("  den(s) =", np.array2string(den_p, precision=4), " poles", np.array2string(np.roots(den_p).real, precision=4))
    Kc = critical_gain(num_p, den_p)
    print(f"largest stable proportional gain: {Kc:.1f} N/(m/s) (L(0) = {Kc / c:.1f})")

    print("\nproportional control: type 0")
    for kp in (100.0, 400.0, 1000.0, 2000.0, 0.95 * Kc):
        nl, dl = open_loop(*p_controller(kp), num_p, den_p)
        K = error_constants(nl, dl)
        nt, dt_ = closed_loop(nl, dl)
        poles = np.roots(dt_)
        zeta = -poles.real / np.abs(poles)
        print(f"  Kp = {kp:7.1f}: L(0) = {K['Kp']:6.2f}, step error {100 * steady_state_error(nl, dl):5.2f} %, "
              f"ramp error inf, min damping ratio of closed-loop poles {zeta.min():.3f}")

    nc, dc = pi_controller(KP_DESIGN, KI_DESIGN)
    nl, dl = open_loop(nc, dc, num_p, den_p)
    K = error_constants(nl, dl)
    nt, dt_ = closed_loop(nl, dl)
    print(f"\nPI control, Kp = {KP_DESIGN:.0f} N/(m/s), Ki = {KI_DESIGN:.0f} N/m: type {system_type(nl, dl)}, "
          f"Kv = Ki/c = {K['Kv']:.4f} 1/s")
    print("  closed-loop poles", np.array2string(np.roots(dt_), precision=4))
    print(f"  step error 0; ramp of 0.5 m/s^2 (a 10 m/s speed change in 20 s): error {steady_state_error(nl, dl, 'ramp', 0.5):.3f} m/s")

    print("\nfinal-value theorem against a 600 s simulation (dt = 0.05 s):")
    cases = [("P 400, step 1 m/s", p_controller(400.0), "step", 1.0),
             ("P 400, 200 N drag step", p_controller(400.0), "disturbance", 200.0),
             ("PI, step 1 m/s", (nc, dc), "step", 1.0),
             ("PI, ramp 0.5 m/s^2", (nc, dc), "ramp", 0.5),
             ("PI, 200 N drag step", (nc, dc), "disturbance", 200.0)]
    for name, ctrl, kind, A in cases:
        pred, sim = fvt_vs_simulation(*ctrl, kind, amplitude=A)
        print(f"  {name:24s}: FVT {pred:+.5f} m/s, simulation at 600 s {sim:+.5f} m/s")
    print(f"  open loop, 200 N drag step: speed loss {200 / c:.3f} m/s")

    ns, ds = sensitivity(nl, dl)
    w = np.logspace(-3, 1, 400)
    S = tf_eval(ns, ds, 1j * w)
    T = tf_eval(nt, dt_, 1j * w)
    print(f"\nPI loop: |S| peak {np.abs(S).max():.3f} at {w[np.argmax(np.abs(S))]:.3f} rad/s; |S(j 0.01)| = "
          f"{abs(tf_eval(ns, ds, 0.01j)):.4f}; |T(j 5)| = {abs(tf_eval(nt, dt_, 5j)):.4f}; max |S + T - 1| = {np.abs(S + T - 1).max():.1e}")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    t = np.arange(0.0, 120.0, 0.05)
    for name, ctrl in (("P 400", p_controller(400.0)), ("P 1000", p_controller(1000.0)), ("PI 400/40", (nc, dc))):
        nl_, dl_ = open_loop(*ctrl, num_p, den_p)
        axes[0].plot(t, simulate_tf(*closed_loop(nl_, dl_), t, np.ones_like(t)), label=name)
        nd, dd = disturbance_to_output(*ctrl, num_p, den_p)
        axes[1].plot(t, simulate_tf(nd, dd, t, -200 * np.ones_like(t)), label=name)
    axes[0].axhline(1.0, color="k", lw=0.6)
    axes[0].set_xlabel("t [s]"); axes[0].set_ylabel("airspeed change [m/s]"); axes[0].set_title("1 m/s step command", fontsize=9)
    axes[0].legend(fontsize=7); axes[0].grid(lw=0.3)
    axes[1].set_xlabel("t [s]"); axes[1].set_ylabel("airspeed change [m/s]"); axes[1].set_title("200 N drag increase", fontsize=9)
    axes[1].legend(fontsize=7); axes[1].grid(lw=0.3)
    axes[2].loglog(w, np.abs(S), label="|S|")
    axes[2].loglog(w, np.abs(T), label="|T|")
    axes[2].set_xlabel("w [rad/s]"); axes[2].set_title("PI loop: S and T = 1 - S", fontsize=9); axes[2].legend(fontsize=7)
    axes[2].grid(True, which="both", lw=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "cruise_control.png", dpi=130)
    plt.close(fig)
    print(f"saved {out_dir / 'cruise_control.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
