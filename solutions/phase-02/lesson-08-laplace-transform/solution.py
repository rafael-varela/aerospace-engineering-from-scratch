"""P02.08 Laplace Transform: reference solution.

A transfer-function solver for constant-coefficient ODEs: the Laplace transform evaluated by
numerical quadrature (to check the table), the transfer function of an actuator built from the ODE
coefficients, poles and residues by hand (partial fractions for simple poles), the inverse transform
as a sum of exponentials, step, impulse and delayed-step responses with initial conditions, and the
frequency response H(i omega) that links back to P02.07. The oracles are scipy.signal (lti step and
impulse), scipy.integrate.solve_ivp, and SymPy's symbolic transform.

Conventions:
    * F(s) = integral_0^inf f(t) e^{-s t} dt. Rational functions are given by coefficient arrays in
      *descending* powers, as numpy.polyval and scipy.signal use them: den = [a2, a1, a0] means a2 s^2 + a1 s + a0.
    * The actuator is a second-order system  y'' + 2 zeta wn y' + wn^2 y = wn^2 u(t), so G(s) = wn^2 / (s^2 + 2 zeta wn s + wn^2).
    * Time signals are zero for t < 0 (causal); u_c(t) is the unit step at t = c.

Run:  python solution.py          (prints the actuator example, writes laplace_transform.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-02/lesson-08-laplace-transform
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
from scipy.integrate import quad


def elevator_actuator() -> dict[str, float]:
    """Return the nominal electro-hydraulic elevator actuator: natural frequency 25 rad/s, damping ratio 0.6.

    The command u(t) is the demanded surface angle [deg]; y(t) is the achieved angle. The ODE is
    y'' + 2 zeta wn y' + wn^2 y = wn^2 u, so a constant command is reached with unit static gain.
    """
    return {"wn": 25.0, "zeta": 0.6}


# ---------------------------------------------------------------- the transform by quadrature

def laplace_numeric(f: Callable[[np.ndarray], np.ndarray], s: float, t_max: float | None = None) -> float:
    """Return integral_0^inf f(t) e^{-s t} dt by adaptive quadrature (scipy.integrate.quad), for real s > 0.

    The integral is split at t_max (default 60 / s) and the tail is integrated to infinity. Used to check
    the transform table and the derivative and shift rules numerically. Raise ``ValueError`` for s <= 0.
    """
    if s <= 0.0:
        raise ValueError("the numerical transform is written for real s > 0")
    if t_max is None:
        t_max = 60.0 / s
    g = lambda t: float(f(t)) * np.exp(-s * t)
    head, _ = quad(g, 0.0, t_max, limit=400)
    tail, _ = quad(g, t_max, np.inf, limit=200)
    return head + tail


def transform_table(name: str, s: complex, **p: float) -> complex:
    """Return the closed-form transform of a table entry at s.

    Names and parameters: 'step' (1/s), 'ramp' (1/s^2), 'exp' (a: 1/(s+a)), 'cos' (w: s/(s^2+w^2)),
    'sin' (w: w/(s^2+w^2)), 'damped_sin' (a, w: w/((s+a)^2+w^2)), 'damped_cos' (a, w: (s+a)/((s+a)^2+w^2)),
    'delayed_step' (c: e^{-cs}/s). Raise ``KeyError`` for an unknown name.
    """
    table = {
        "step": lambda: 1.0 / s,
        "ramp": lambda: 1.0 / s ** 2,
        "exp": lambda: 1.0 / (s + p["a"]),
        "cos": lambda: s / (s ** 2 + p["w"] ** 2),
        "sin": lambda: p["w"] / (s ** 2 + p["w"] ** 2),
        "damped_sin": lambda: p["w"] / ((s + p["a"]) ** 2 + p["w"] ** 2),
        "damped_cos": lambda: (s + p["a"]) / ((s + p["a"]) ** 2 + p["w"] ** 2),
        "delayed_step": lambda: np.exp(-p["c"] * s) / s,
    }
    if name not in table:
        raise KeyError(f"no table entry named {name!r}")
    return table[name]()


# ---------------------------------------------------------------- rational transfer functions

def transfer_function(wn: float, zeta: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (num, den) of G(s) = wn^2 / (s^2 + 2 zeta wn s + wn^2), descending powers."""
    if wn <= 0.0 or zeta < 0.0:
        raise ValueError("need wn > 0 and zeta >= 0")
    return np.array([wn * wn]), np.array([1.0, 2.0 * zeta * wn, wn * wn])


def evaluate(num: np.ndarray, den: np.ndarray, s: complex | np.ndarray) -> complex | np.ndarray:
    """Return N(s) / D(s) for coefficient arrays in descending powers (numpy.polyval convention)."""
    return np.polyval(num, s) / np.polyval(den, s)


def poles(den: np.ndarray) -> np.ndarray:
    """Return the roots of the denominator polynomial (numpy.roots), sorted by real part then imaginary part."""
    r = np.roots(np.asarray(den, dtype=float))
    return r[np.lexsort((r.imag, r.real))]


def is_stable(den: np.ndarray) -> bool:
    """Return True if every pole has strictly negative real part (every mode decays)."""
    return bool(np.all(poles(den).real < 0.0))


def residues(num: np.ndarray, den: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (poles, residues) of the strictly proper rational function N(s)/D(s) with *simple* poles.

    For a simple pole p the residue is N(p) / D'(p) (the cover-up rule in the limit). Raise ``ValueError``
    if deg N >= deg D or if two poles coincide to within 1e-8 relative (a repeated pole needs the
    derivative formula, left to SymPy in this lesson).
    """
    num = np.atleast_1d(np.asarray(num, dtype=float))
    den = np.atleast_1d(np.asarray(den, dtype=float))
    if num.size >= den.size:
        raise ValueError("residue expansion needs a strictly proper function (deg N < deg D)")
    p = poles(den)
    scale = max(1.0, np.abs(p).max())
    for i in range(p.size):
        for j in range(i + 1, p.size):
            if abs(p[i] - p[j]) < 1e-8 * scale:
                raise ValueError("repeated pole: the simple-pole residue formula does not apply")
    dden = np.polyder(den)
    r = np.polyval(num, p) / np.polyval(dden, p)
    return p, r


def inverse_laplace(num: np.ndarray, den: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Return f(t) = sum_k r_k e^{p_k t} for t >= 0 (zero for t < 0): the inverse transform of a strictly proper N/D with simple poles.

    The result is real for real coefficients (complex poles come in conjugate pairs); the imaginary part is
    dropped after checking that it is at rounding level.
    """
    t = np.asarray(t, dtype=float)
    p, r = residues(num, den)
    f = np.zeros(t.shape, dtype=complex)
    for pk, rk in zip(p, r):
        f += rk * np.exp(pk * t)
    f = np.where(t >= 0.0, f, 0.0)
    if np.max(np.abs(f.imag)) > 1e-8 * max(1.0, np.max(np.abs(f.real))):
        raise RuntimeError("inverse transform came out complex: coefficients were not real")
    return f.real


# ---------------------------------------------------------------- responses of the actuator

def response_with_initial_conditions(num: np.ndarray, den: np.ndarray, U_num: np.ndarray, U_den: np.ndarray,
                                      y0: float, ydot0: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (num, den) of Y(s) for the second-order ODE D(s) y = N u with y(0) = y0, y'(0) = ydot0 and input U(s) = U_num / U_den.

    Transforming y'' + a1 y' + a0 y = b0 u with the derivative rule gives
        (s^2 + a1 s + a0) Y = b0 U + (s + a1) y0 + ydot0,
    so Y = [N U + ((s + a1) y0 + ydot0)] / D. ``den`` must be monic and of degree 2.
    """
    num = np.atleast_1d(np.asarray(num, dtype=float))
    den = np.asarray(den, dtype=float)
    if den.size != 3 or den[0] != 1.0:
        raise ValueError("written for a monic second-order denominator [1, a1, a0]")
    ic_poly = np.array([y0, den[1] * y0 + ydot0])              # (s + a1) y0 + ydot0
    forced = np.polymul(num, np.atleast_1d(np.asarray(U_num, dtype=float)))
    free = np.polymul(ic_poly, np.atleast_1d(np.asarray(U_den, dtype=float)))
    Y_num = np.polyadd(forced, free)
    Y_den = np.polymul(den, np.atleast_1d(np.asarray(U_den, dtype=float)))
    return Y_num, Y_den


def step_response(wn: float, zeta: float, t: np.ndarray) -> np.ndarray:
    """Return the unit-step response of the actuator from rest, via Y(s) = G(s) / s and the inverse transform."""
    num, den = transfer_function(wn, zeta)
    Y_num, Y_den = response_with_initial_conditions(num, den, np.array([1.0]), np.array([1.0, 0.0]), 0.0, 0.0)
    return inverse_laplace(Y_num, Y_den, t)


def step_response_closed_form(wn: float, zeta: float, t: np.ndarray) -> np.ndarray:
    """Return the textbook underdamped step response 1 - e^{-zeta wn t} [cos(wd t) + (zeta / sqrt(1 - zeta^2)) sin(wd t)].

    Raise ``ValueError`` unless 0 <= zeta < 1 (the closed form used here is the underdamped one).
    """
    if not 0.0 <= zeta < 1.0:
        raise ValueError("closed form written for 0 <= zeta < 1")
    t = np.asarray(t, dtype=float)
    wd = wn * np.sqrt(1.0 - zeta ** 2)
    y = 1.0 - np.exp(-zeta * wn * t) * (np.cos(wd * t) + zeta / np.sqrt(1.0 - zeta ** 2) * np.sin(wd * t))
    return np.where(t >= 0.0, y, 0.0)


def impulse_response(wn: float, zeta: float, t: np.ndarray) -> np.ndarray:
    """Return the unit-impulse response g(t) = L^{-1}{G(s)}: the derivative of the step response."""
    num, den = transfer_function(wn, zeta)
    return inverse_laplace(num, den, t)


def delayed_step_response(wn: float, zeta: float, t: np.ndarray, delay: float) -> np.ndarray:
    """Return the response to the step u_c(t), c = delay >= 0: by the shift rule, y(t) = y_step(t - c) u_c(t)."""
    if delay < 0.0:
        raise ValueError("delay must be non-negative")
    return step_response(wn, zeta, np.asarray(t, dtype=float) - delay)


def actuator_rhs(t: float, y: np.ndarray, wn: float, zeta: float, u: Callable[[float], float]) -> np.ndarray:
    """Return dy/dt for the state (y, y') of y'' + 2 zeta wn y' + wn^2 y = wn^2 u(t): the solve_ivp check."""
    return np.array([y[1], wn * wn * u(t) - 2.0 * zeta * wn * y[1] - wn * wn * y[0]])


def convolve_causal(g: np.ndarray, u: np.ndarray, dt: float) -> np.ndarray:
    """Return the causal convolution y(t_n) = sum_{k<=n} g(t_k) u(t_n - t_k) dt on a uniform grid (trapezoid-free, first order).

    The time-domain meaning of Y(s) = G(s) U(s): the output is the impulse response weighted by the input's past.
    """
    g = np.asarray(g, dtype=float)
    u = np.asarray(u, dtype=float)
    return np.convolve(g, u)[: g.size] * dt


# ---------------------------------------------------------------- performance numbers and frequency response

def overshoot(zeta: float) -> float:
    """Return the fractional peak overshoot exp(-pi zeta / sqrt(1 - zeta^2)) of the underdamped step response (0 <= zeta < 1)."""
    if not 0.0 <= zeta < 1.0:
        raise ValueError("overshoot is defined for 0 <= zeta < 1")
    return float(np.exp(-np.pi * zeta / np.sqrt(1.0 - zeta ** 2)))


def peak_time(wn: float, zeta: float) -> float:
    """Return t_p = pi / (wn sqrt(1 - zeta^2)), the time of the first peak of the step response."""
    if not 0.0 <= zeta < 1.0:
        raise ValueError("peak time is defined for 0 <= zeta < 1")
    return float(np.pi / (wn * np.sqrt(1.0 - zeta ** 2)))


def settling_time(t: np.ndarray, y: np.ndarray, final: float, band: float = 0.02) -> float:
    """Return the first time after which |y - final| stays within band * |final| (2 % by default) on the sampled response.

    Raise ``ValueError`` if the response never settles on the grid.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    outside = np.abs(y - final) > band * abs(final)
    if outside[-1]:
        raise ValueError("response has not settled by the end of the grid")
    idx = np.where(outside)[0]
    return float(t[idx[-1] + 1]) if idx.size else float(t[0])


def frequency_response(num: np.ndarray, den: np.ndarray, omega: np.ndarray) -> np.ndarray:
    """Return H(i omega) = N(i omega) / D(i omega): the complex gain of P02.07 as the transfer function on the imaginary axis."""
    return evaluate(num, den, 1j * np.asarray(omega, dtype=float))


def laplace_symbolic_step_response(wn: float, zeta: float):
    """Return SymPy's y(t) for the actuator step response, via sympy.inverse_laplace_transform of G(s)/s.

    Used once as an independent check of the residue arithmetic; returns a SymPy expression in t.
    """
    import sympy as sp

    s, t = sp.symbols("s t", positive=True)
    wn_, z_ = sp.nsimplify(wn), sp.nsimplify(zeta)
    Y = wn_ ** 2 / (s * (s ** 2 + 2 * z_ * wn_ * s + wn_ ** 2))
    return sp.simplify(sp.inverse_laplace_transform(Y, s, t))


def main(out_dir: str | Path = ".") -> None:
    """Print the actuator analysis and save the step/impulse/delayed-step figure with the pole map."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.integrate import solve_ivp
    from scipy.signal import lti

    out_dir = Path(out_dir)
    act = elevator_actuator()
    wn, zeta = act["wn"], act["zeta"]
    print("transform table against quadrature at s = 2.0:")
    for name, f, p in (("step", lambda t: 1.0, {}), ("ramp", lambda t: t, {}), ("exp", lambda t: np.exp(-3.0 * t), {"a": 3.0}),
                       ("cos", lambda t: np.cos(5.0 * t), {"w": 5.0}),
                       ("damped_sin", lambda t: np.exp(-1.0 * t) * np.sin(4.0 * t), {"a": 1.0, "w": 4.0})):
        print(f"  {name:11s} table {transform_table(name, 2.0, **p):.10f}  quad {laplace_numeric(f, 2.0):.10f}")

    num, den = transfer_function(wn, zeta)
    p, r = residues(np.array([1.0]), np.polymul(den, [1.0, 0.0]))
    print(f"\nactuator G(s) = {wn ** 2:.0f} / (s^2 + {2 * zeta * wn:.0f} s + {wn ** 2:.0f}); poles {poles(den)}; stable: {is_stable(den)}")
    print(f"Y(s) = G/s has poles {p} with residues {wn ** 2 * r}")
    t = np.linspace(0.0, 0.5, 5001)
    y = step_response(wn, zeta, t)
    y_cf = step_response_closed_form(wn, zeta, t)
    _, y_lti = lti(num, den).step(T=t)
    sol = solve_ivp(actuator_rhs, (0.0, 0.5), [0.0, 0.0], t_eval=t, args=(wn, zeta, lambda tt: 1.0), rtol=1e-10, atol=1e-12)
    print(f"step response: max |residues - closed form| = {np.max(np.abs(y - y_cf)):.2e}, vs scipy.signal.lti {np.max(np.abs(y - y_lti)):.2e},"
          f" vs solve_ivp {np.max(np.abs(y - sol.y[0])):.2e}")
    print(f"  overshoot {100 * overshoot(zeta):.2f} % at t_p = {peak_time(wn, zeta) * 1e3:.1f} ms (grid peak {100 * (y.max() - 1):.2f} % at {t[np.argmax(y)] * 1e3:.1f} ms);"
          f" 2 % settling {settling_time(t, y, 1.0) * 1e3:.0f} ms (4 / (zeta wn) = {4e3 / (zeta * wn):.0f} ms)")
    g = impulse_response(wn, zeta, t)
    dt = t[1] - t[0]
    print(f"impulse response: max |g - dy_step/dt| = {np.max(np.abs(g[1:-1] - np.gradient(y, dt)[1:-1])):.2e};"
          f" integral of g = {np.trapezoid(g, t):.5f} (-> 1 = G(0))")
    y_conv = convolve_causal(g, np.ones_like(t), dt)
    print(f"convolution of g with a step vs step response: max error {np.max(np.abs(y_conv - y)):.2e} (first-order rule, dt = {dt * 1e3:.2f} ms)")
    print(f"symbolic (SymPy) step response: {laplace_symbolic_step_response(wn, zeta)}")
    omega = np.array([1.0, 10.0, 25.0, 50.0])
    H = frequency_response(num, den, omega)
    print("frequency response |H|, phase lag [deg]:", [(f"{abs(h):.3f}", f"{-np.degrees(np.angle(h)):.1f}") for h in H])

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    ax = axes[0]
    ax.plot(t, y, label="step (residues)")
    ax.plot(t, delayed_step_response(wn, zeta, t, 0.2), label="step delayed 0.2 s")
    ax.plot(t, g / wn, "--", label="impulse / wn")
    ax.axhline(1.0, color="0.6", lw=0.7)
    ax.set_xlabel("t [s]")
    ax.set_ylabel("y")
    ax.set_title(f"actuator wn = {wn} rad/s, zeta = {zeta}")
    ax.grid(True, lw=0.4)
    ax.legend(fontsize=8)
    ax = axes[1]
    for z in (0.2, 0.6, 1.0, 2.0):
        pz = poles(transfer_function(wn, z)[1])
        ax.plot(pz.real, pz.imag, "x", ms=8, mew=2, label=f"zeta = {z}")
    ax.axvline(0.0, color="k", lw=0.6)
    ax.axhline(0.0, color="k", lw=0.6)
    ax.set_xlabel("Re s [1/s]")
    ax.set_ylabel("Im s [rad/s]")
    ax.set_title("pole map: left half-plane = stable")
    ax.set_aspect("equal")
    ax.grid(True, lw=0.4)
    ax.legend(fontsize=8)
    ax = axes[2]
    w = np.logspace(0, 3, 400)
    for z in (0.2, 0.6, 1.0):
        Hz = frequency_response(*transfer_function(wn, z), w)
        ax.loglog(w, np.abs(Hz), label=f"zeta = {z}")
    ax.set_xlabel("omega [rad/s]")
    ax.set_ylabel("|G(i omega)|")
    ax.set_title("frequency response = G on the imaginary axis")
    ax.grid(True, which="both", lw=0.4)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "laplace_transform.png", dpi=150)
    print(f"saved {out_dir / 'laplace_transform.png'}")


if __name__ == "__main__":
    main()
