"""P09.03 Transfer Functions & Block Diagrams: starter code.

A block-diagram reducer built on polynomial algebra. A transfer function is a pair (num, den) of
coefficient arrays, highest power of s first, normalised so that den is monic. Series, parallel and
feedback connections are polynomial products and sums; poles and zeros are polynomial roots; the DC gain
is G(0). Because the algebra never cancels anything by itself, ``minreal`` removes pole-zero pairs that
coincide. Every reduction is checked by an independent route: the unreduced diagram is simulated in the
time domain, block by block (each block in controllable canonical state-space form, the whole
interconnection integrated by fixed-step RK4), and the result is compared with the step response of the
reduced transfer function.

The example is a pitch-attitude autopilot for a light aircraft or small UAV: a short-period pitch-rate
model, a first-order elevator servo, an integrator from pitch rate to pitch angle, a rate-gyro filter, an
inner pitch-rate (damper) loop and an outer attitude loop. The numbers are a representative model, not
any specific aircraft's data.

Conventions:
    * TF = (num, den), numpy arrays, highest power first, den monic after ``tf``.
    * Elevator sign: positive delta_e gives a nose-up pitching moment (so M_delta > 0 here).
    * Angles in rad, rates in rad/s, time in s.

Complete every function marked NotImplementedError.

Run:  python main.py          (prints the reduction and its checks, writes block_diagram.png)
Test: pytest curriculum/phase-09/lesson-03-transfer-functions-block-diagrams
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np

TF = tuple[np.ndarray, np.ndarray]


# ---------------------------------------------------------------- polynomial helpers (given)

def trim(p: np.ndarray, tol: float = 0.0) -> np.ndarray:
    """Return p without leading coefficients whose magnitude is <= tol * max|p| (at least one coefficient kept)."""
    p = np.atleast_1d(np.asarray(p, dtype=complex if np.iscomplexobj(p) else float))
    if not np.any(p):
        return p[-1:] * 0
    scale = np.max(np.abs(p))
    i = 0
    while i < len(p) - 1 and abs(p[i]) <= tol * scale:
        i += 1
    return p[i:]


def poly_add(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Return the coefficient array of p(s) + q(s), aligning the constant terms (highest power first)."""
    p, q = np.atleast_1d(np.asarray(p, dtype=float)), np.atleast_1d(np.asarray(q, dtype=float))
    n = max(len(p), len(q))
    return np.concatenate([np.zeros(n - len(p)), p]) + np.concatenate([np.zeros(n - len(q)), q])


# ---------------------------------------------------------------- transfer-function algebra

def tf(num, den) -> TF:
    """Return the normalised TF (num, den): leading zeros trimmed, den divided by its leading coefficient.

    Raise ``ValueError`` if den is identically zero or the TF is improper (deg num > deg den).
    """
    raise NotImplementedError


def evaluate(G: TF, s) -> np.ndarray:
    """Return G(s) = num(s) / den(s) at the (complex) point(s) s."""
    raise NotImplementedError


def series(*Gs: TF) -> TF:
    """Return the cascade G1 G2 ... Gk: numerators and denominators multiply."""
    raise NotImplementedError


def parallel(G1: TF, G2: TF) -> TF:
    """Return G1 + G2 = (n1 d2 + n2 d1) / (d1 d2): two paths from one input summed at one output."""
    raise NotImplementedError


def feedback(G: TF, H: TF | None = None, sign: int = -1) -> TF:
    """Return the closed loop G / (1 - sign G H) = nG dH / (dG dH - sign nG nH).

    ``sign = -1`` (the default) is negative feedback; ``H = None`` means unity feedback H = 1.
    """
    raise NotImplementedError


def gain(k: float) -> TF:
    """Return the static gain k as a TF."""
    return tf([k], [1.0])


def poles(G: TF) -> np.ndarray:
    """Return the roots of den(s), sorted by real part (most negative first)."""
    raise NotImplementedError


def zeros(G: TF) -> np.ndarray:
    """Return the roots of num(s) (an empty array for a constant numerator)."""
    raise NotImplementedError


def dc_gain(G: TF) -> float:
    """Return G(0) = num(0)/den(0); ``inf`` (with the sign of num(0)) if there is a pole at s = 0."""
    raise NotImplementedError


def minreal(G: TF, tol: float = 1e-6) -> TF:
    """Return G with every zero that coincides with a pole (within tol, relative to max(1, |p|)) cancelled.

    The gain is fixed so that the reduced and the original TF agree at a test point away from all roots.
    """
    raise NotImplementedError


def tf_to_ss(G: TF) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (A, B, C, D) in controllable canonical form for a proper TF.

    With den = s^n + a_{n-1} s^{n-1} + ... + a_0 and D = num_n (the s^n coefficient of num), the strictly
    proper remainder num - D den = c_{n-1} s^{n-1} + ... + c_0 gives A = companion matrix (last row
    -a_0 ... -a_{n-1}), B = e_n, C = (c_0, ..., c_{n-1}).
    """
    raise NotImplementedError


def simulate_tf(G: TF, u: Callable[[float], float], t: np.ndarray, substeps: int = 4) -> np.ndarray:
    """Return y(t_k) for the TF G driven by the input function u(t), from rest, by fixed-step RK4 on (A, B, C, D)."""
    raise NotImplementedError


def simulate_diagram(blocks: dict[str, TF],
                     connect: Callable[[dict[str, float], float], dict[str, float]],
                     r: Callable[[float], float], t: np.ndarray, substeps: int = 4) -> dict[str, np.ndarray]:
    """Simulate an unreduced block diagram and return each block's output history {name: y(t_k)}.

    Every block is a strictly proper TF realised by ``tf_to_ss``, so its output C x depends on its state
    only and the diagram has no algebraic loop. ``connect(outputs, r_value)`` returns each block's input
    from the current block outputs and the reference (static gains and summing junctions live there).
    The joint state is integrated from rest by fixed-step RK4. Raise ``ValueError`` for a block with D != 0.
    """
    names = list(blocks)
    ss = {}
    for nm in names:
        A, B, C, D = tf_to_ss(blocks[nm])
        if abs(D[0, 0]) > 0.0:
            raise ValueError(f"block {nm} is not strictly proper")
        ss[nm] = (A, B[:, 0], C[0])
    sizes = [ss[nm][0].shape[0] for nm in names]
    offs = np.concatenate([[0], np.cumsum(sizes)])

    def outputs(x: np.ndarray) -> dict[str, float]:
        return {nm: float(ss[nm][2] @ x[offs[i]:offs[i + 1]]) for i, nm in enumerate(names)}

    def rhs(tt: float, x: np.ndarray) -> np.ndarray:
        u = connect(outputs(x), r(tt))
        dx = np.empty_like(x)
        for i, nm in enumerate(names):
            A, B, _ = ss[nm]
            xi = x[offs[i]:offs[i + 1]]
            dx[offs[i]:offs[i + 1]] = A @ xi + B * u[nm]
        return dx

    x = np.zeros(offs[-1])
    hist = {nm: np.empty(len(t)) for nm in names}
    out = outputs(x)
    for nm in names:
        hist[nm][0] = out[nm]
    for k in range(len(t) - 1):
        h = (t[k + 1] - t[k]) / substeps
        tt = t[k]
        for _ in range(substeps):
            k1 = rhs(tt, x)
            k2 = rhs(tt + h / 2, x + h / 2 * k1)
            k3 = rhs(tt + h / 2, x + h / 2 * k2)
            k4 = rhs(tt + h, x + h * k3)
            x = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            tt += h
        out = outputs(x)
        for nm in names:
            hist[nm][k + 1] = out[nm]
    return hist


# ---------------------------------------------------------------- the pitch-attitude loop

def pitch_model() -> dict[str, float]:
    """Return representative short-period, servo and sensor parameters (a light aircraft or small UAV, typical values).

    q/delta_e = M_delta (s + 1/T_theta2) / (s^2 + 2 zeta_sp w_sp s + w_sp^2); servo 1/(tau_s s + 1); rate-gyro
    filter 1/(tau_g s + 1); autopilot gains K_q (rate damper, s) and K_theta (attitude, rad/rad).
    """
    return {"w_sp": 4.0, "zeta_sp": 0.4, "inv_T_theta2": 1.5, "M_delta": 12.0, "tau_s": 0.05, "tau_g": 0.01,
            "K_q": 0.3, "K_theta": 3.0}


def pitch_blocks(p: dict[str, float]) -> dict[str, TF]:
    """Return the dynamic blocks of the pitch loop: servo, short_period (q/delta_e), integrator (theta/q), gyro (q filter)."""
    return {
        "servo": tf([1.0], [p["tau_s"], 1.0]),
        "short_period": tf([p["M_delta"], p["M_delta"] * p["inv_T_theta2"]], [1.0, 2 * p["zeta_sp"] * p["w_sp"], p["w_sp"] ** 2]),
        "integrator": tf([1.0], [1.0, 0.0]),
        "gyro": tf([1.0], [p["tau_g"], 1.0]),
    }


def pitch_connect(p: dict[str, float]) -> Callable[[dict[str, float], float], dict[str, float]]:
    """Return the wiring of the pitch loop for ``simulate_diagram``.

    servo command u = K_theta (theta_c - theta) - K_q q_measured; short period driven by delta_e; integrator
    by q; gyro filter by q (its output is the measured rate).
    """
    def connect(y: dict[str, float], theta_c: float) -> dict[str, float]:
        u = p["K_theta"] * (theta_c - y["integrator"]) - p["K_q"] * y["gyro"]
        return {"servo": u, "short_period": y["servo"], "integrator": y["short_period"], "gyro": y["short_period"]}
    return connect


def reduce_pitch_loop(p: dict[str, float]) -> dict[str, TF]:
    """Return the reduction step by step: forward (q/u open loop), inner (q/u with the damper), theta_u, closed (theta/theta_c).

    inner = feedback(servo * short_period, K_q * gyro); theta_u = inner * integrator;
    closed = feedback(K_theta * theta_u, 1).
    """
    raise NotImplementedError


def main(out_dir: str | Path = ".") -> None:
    """Print the pitch-loop reduction and its checks; save block_diagram.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    p = pitch_model()
    red = reduce_pitch_loop(p)
    np.set_printoptions(precision=4, suppress=True)
    for k, G in red.items():
        print(f"{k:8s}: num {G[0]}, den {G[1]}")
        print(f"          poles {np.round(poles(G), 3)}, zeros {np.round(zeros(G), 3)}, DC gain {dc_gain(G):.4g}")

    t = np.linspace(0.0, 6.0, 1201)
    step = lambda tt: 1.0
    theta_red = simulate_tf(red["closed"], step, t)
    hist = simulate_diagram(pitch_blocks(p), pitch_connect(p), step, t)
    theta_dia = hist["integrator"]
    print(f"\nunit step in theta_c: max |reduced TF - unreduced diagram| = {np.max(np.abs(theta_red - theta_dia)):.2e} rad")
    i_pk = int(np.argmax(theta_dia))
    out = np.where(np.abs(theta_dia - 1.0) > 0.02)[0]
    print(f"  peak {theta_dia[i_pk]:.4f} at t = {t[i_pk]:.2f} s (overshoot {100 * (theta_dia[i_pk] - 1):.1f} %),"
          f" 2 % settling {t[out[-1] + 1]:.2f} s, theta(6 s) = {theta_dia[-1]:.4f}")
    print(f"  peak elevator {np.max(np.abs(hist['servo'])):.3f} rad = {np.degrees(np.max(np.abs(hist['servo']))):.1f} deg per rad of command;"
          f" for a 5 deg step: {5 * np.max(np.abs(hist['servo'])):.1f} deg")

    print("\nclosed-loop poles against the damper gain K_q (K_theta = 3):")
    for Kq in (0.0, 0.1, 0.3, 0.5):
        q = dict(p, K_q=Kq)
        pl = poles(reduce_pitch_loop(q)["closed"])
        print(f"  K_q = {Kq:.1f}: max Re = {pl.real.max():+.3f}, poles {np.round(pl, 2)}")
    q = dict(p, K_q=0.0, K_theta=4.0)
    print(f"  K_q = 0, K_theta = 4: max Re = {poles(reduce_pitch_loop(q)['closed']).real.max():+.3f} (unstable without the damper)")

    cl = red["closed"]
    print(f"\nclosed-loop zero at {np.round(zeros(cl), 3)}: -1/T_theta2 = {-p['inv_T_theta2']} and the gyro-filter pole -1/tau_g = {-1 / p['tau_g']:.0f}")
    G1 = tf([1.0, 2.0], [1.0, 3.0, 2.0])
    print(f"(s + 2)/((s + 1)(s + 2)) before minreal: den {G1[1]}; after: {minreal(G1)[1]}")

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.8))
    axes[0].plot(t, theta_dia, label="unreduced diagram (RK4)")
    axes[0].plot(t, theta_red, "--", label="reduced TF")
    axes[0].axhline(1.0, color="k", lw=0.5)
    axes[0].set_xlabel("t [s]"); axes[0].set_ylabel("theta / theta_c"); axes[0].legend(fontsize=8); axes[0].grid(lw=0.4)
    axes[0].set_title("pitch-attitude step", fontsize=9)
    for Kq, mk in ((0.0, "x"), (0.1, "+"), (0.3, "o"), (0.5, "s")):
        pl = poles(reduce_pitch_loop(dict(p, K_q=Kq))["closed"])
        pl = pl[pl.real > -40]
        axes[1].plot(pl.real, pl.imag, mk, ls="none", label=f"K_q = {Kq}")
    axes[1].axvline(0, color="k", lw=0.5)
    axes[1].set_xlabel("Re s [1/s]"); axes[1].set_ylabel("Im s [rad/s]"); axes[1].legend(fontsize=8); axes[1].grid(lw=0.4)
    axes[1].set_title("closed-loop poles (those with Re > -40)", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "block_diagram.png", dpi=150)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'block_diagram.png'}")


if __name__ == "__main__":
    import sys

    main(sys.argv[1] if len(sys.argv) > 1 else ".")
