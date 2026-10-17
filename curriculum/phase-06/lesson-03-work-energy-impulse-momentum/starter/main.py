"""P06.03 Work-Energy & Impulse-Momentum: starter code.

A jet catches a cable on a carrier deck at 65 m/s. How far does it run out, how hard does it decelerate
at the worst moment, and how long does it take? The cable and its energy absorber are lumped into one
illustrative force law acting against the aircraft's motion,

    F(x, v) = F0 + k x + c v^2       (N; x = payout from engagement, v = speed),

a constant drag F0, a term k x that grows with payout, and a hydraulic-damper term c v^2. Three routes:

* work-energy: with u = v^2 as a function of x, d(m u / 2)/dx = -F gives a *linear* ODE
  u' + (2c/m) u = -(2/m)(F0 + k x), solved in closed form; the stopping distance is the root of u = 0,
  and the force along the stroke obeys F' = k - (2c/m) F, so it is monotonic and the peak is at
  engagement or at the stop;
* impulse-momentum: the force's impulse removes the momentum, J = integral F dt = m v0, and the
  stopping time is t = integral dx / v;
* time integration of m x'' = -F with `solve_ivp`, the independent check of both.

You write every function whose body is `raise NotImplementedError`. The helpers `_phi1` and `_phi2`
(for a closed form that stays accurate as c -> 0) and `main` are given.

The force law and all numbers are illustrative, not data for any real arresting gear. Units: SI.

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes arrestment.png)
Test: pytest curriculum/phase-06/lesson-03-work-energy-impulse-momentum
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

G0 = 9.80665
M_JET = 20000.0       # kg, illustrative carrier jet at landing weight
V_ENGAGE = 65.0       # m/s, illustrative engagement speed relative to the deck
X_RUNOUT = 100.0      # m, illustrative design runout


# ---------------------------------------------------------------- the force law

def arrestor_force(x, v, F0: float, k: float, c: float):
    """Retarding force F = F0 + k x + c v^2 (N) at payout x (m) and speed v >= 0 (m/s)."""
    raise NotImplementedError


# ---------------------------------------------------------------- work-energy

def _phi1(y):
    """(1 - e^-y) / y, = 1 at y = 0, accurate for small y."""
    y = np.asarray(y, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(y < 1e-8, 1.0 - y / 2.0, -np.expm1(-y) / np.where(y == 0, 1.0, y))


def _phi2(y):
    """(e^-y - 1 + y) / y^2, = 1/2 at y = 0; a series below y = 0.01 avoids the cancellation."""
    y = np.asarray(y, dtype=float)
    series = 0.5 - y / 6.0 + y ** 2 / 24.0 - y ** 3 / 120.0 + y ** 4 / 720.0 - y ** 5 / 5040.0
    with np.errstate(invalid="ignore", divide="ignore"):
        direct = (np.expm1(-y) + y) / np.where(y == 0, 1.0, y) ** 2
    return np.where(y < 1e-2, series, direct)


def v_squared(x, v0: float, m: float, F0: float, k: float, c: float):
    """Work-energy solution u(x) = v(x)^2. From d(m u/2)/dx = -(F0 + k x + c u), a linear first-order ODE
    in x with integrating factor e^(alpha x), alpha = 2c/m:

        u = (v0^2 - a) e^(-alpha x) + a + b x,   b = -k/c,  a = -F0/c + k m / (2 c^2)      (c > 0),
        u = v0^2 - (2/m)(F0 x + k x^2 / 2)                                                (c = 0).

    Computed in the equivalent form u = v0^2 e^(-y) - (2 F0/m) x phi1(y) - (2k/m) x^2 phi2(y), y = alpha x,
    with phi1 = (1 - e^-y)/y and phi2 = (e^-y - 1 + y)/y^2, which covers both cases and stays accurate
    when c is small (the first form then subtracts two huge numbers). Valid for 0 <= x <= stopping
    distance (beyond it u < 0: the aircraft has stopped)."""
    raise NotImplementedError


def force_along_stroke(x, v0: float, m: float, F0: float, k: float, c: float):
    """F(x) = F0 + k x + c u(x) along the stroke, from ``v_squared``."""
    raise NotImplementedError


def stopping_distance(v0: float, m: float, F0: float, k: float, c: float) -> float:
    """Distance at which u(x) = 0, by bracketing and `brentq` on the closed form. Returns inf when the
    force law cannot stop the aircraft (F0 = k = 0: a pure v^2 damper only slows it exponentially)."""
    raise NotImplementedError


def peak_force(v0: float, m: float, F0: float, k: float, c: float) -> tuple[float, float]:
    """Peak retarding force and where it acts, by work-energy alone. Differentiating F = F0 + k x + c u
    and using the energy ODE gives F'(x) = k - (2c/m) F: F relaxes monotonically toward k m / (2c), so
    the largest value is at engagement, F(0) = F0 + c v0^2, or at the stop, F(x_s) = F0 + k x_s.
    Returns (F_peak, x_at_peak)."""
    raise NotImplementedError


def constant_force_design(v0: float, m: float, x_stop: float, F0: float = 0.0) -> tuple[float, float]:
    """The (c, k) that make the force constant along the whole stroke: F' = k - (2c/m) F = 0 needs
    k = 2 c F / m, and a constant F stops the aircraft in x_stop when F = m v0^2 / (2 x_stop).
    With F = F0 + c v0^2 at engagement, c = (F - F0) / v0^2. Returns (c, k)."""
    raise NotImplementedError


# ---------------------------------------------------------------- impulse-momentum

def stopping_time(v0: float, m: float, F0: float, k: float, c: float) -> float:
    """Stopping time t = integral_0^{x_s} dx / v(x) with v = sqrt(u) from work-energy. The integrand has
    an integrable 1/sqrt singularity at x_s; the substitution x = x_s (1 - s^2) removes it:
    t = integral_0^1 2 x_s s / sqrt(u(x_s (1 - s^2))) ds."""
    raise NotImplementedError


def momentum_average_force(v0: float, m: float, t_stop: float) -> float:
    """Time-average force from impulse-momentum: integral F dt = m v0, so F_avg,t = m v0 / t_stop."""
    raise NotImplementedError


def energy_average_force(v0: float, m: float, x_stop: float) -> float:
    """Distance-average force from work-energy: integral F dx = m v0^2 / 2, so F_avg,x = m v0^2 / (2 x_stop)."""
    raise NotImplementedError


# ---------------------------------------------------------------- time integration (the independent check)

def simulate_arrestment(v0: float, m: float, F0: float, k: float, c: float, t_max: float = 60.0,
                        rtol: float = 1e-10, atol: float = 1e-10):
    """Integrate m x'' = -F(x, x') from x = 0, x' = v0 with `solve_ivp` (RK45, dense output) until the
    speed reaches zero (terminal event). Returns the result; ``sol.t_events[0]`` is the stopping time."""
    raise NotImplementedError


def energy_audit(sol, m: float, F0: float, k: float, c: float, n: int = 4001) -> dict:
    """Where the kinetic energy went, from a simulated run: the work of F0 (F0 x_s), the energy held by
    the k x term (k x_s^2 / 2: stored, and returned unless the gear locks), and the damper's dissipation
    integral c v^3 dt; plus the impulse integral F dt. Returns a dict of these and the initial KE."""
    raise NotImplementedError


def efficiency(v0: float, m: float, x_stop: float, F_peak: float) -> float:
    """Stroke efficiency: the distance-average force over the peak force. 1 for a constant force, which
    gives the smallest possible peak deceleration for a given runout."""
    raise NotImplementedError


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save arrestment.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    m, v0 = M_JET, V_ENGAGE
    print(f"illustrative jet {m:.0f} kg engaging at {v0} m/s: kinetic energy {0.5 * m * v0 ** 2 / 1e6:.2f} MJ, "
          f"momentum {m * v0 / 1e6:.3f} MN s")
    a_ideal = v0 ** 2 / (2 * X_RUNOUT)
    print(f"ideal constant force over {X_RUNOUT:.0f} m: {m * a_ideal / 1e3:.1f} kN, {a_ideal:.3f} m/s^2 = {a_ideal / G0:.3f} g, "
          f"time {2 * X_RUNOUT / v0:.4f} s")
    c0, k0 = constant_force_design(v0, m, X_RUNOUT)
    print(f"constant-force design (F0 = 0): c = {c0:.3f} N s^2/m^2, k = {k0:.3f} N/m")

    designs = {
        "constant force": (0.0, k0, c0),
        "damper-heavy": (0.0, 0.5 * k0, None),
        "spring-heavy": (0.0, 1.5 * k0, None),
        "spring only": (0.0, 2.0 * k0, 0.0),     # v0 sqrt(m / k) = 100 m with no damper at all
        "brake only": (m * a_ideal, 0.0, 0.0),    # constant force, all of it dissipated
    }
    print("\n design           F0 (kN)   k (N/m)    c (N s^2/m^2)  x_stop (m)  t_stop (s)  F_peak (kN)  at x (m)  peak g  efficiency")
    rows = {}
    for name, (F0, k, c) in designs.items():
        if c is None:   # choose c so that the runout is still X_RUNOUT
            c = brentq(lambda cc: stopping_distance(v0, m, F0, k, cc) - X_RUNOUT, 1e-6, 1e4, xtol=1e-12)
        xs = stopping_distance(v0, m, F0, k, c)
        ts = stopping_time(v0, m, F0, k, c)
        Fp, xp = peak_force(v0, m, F0, k, c)
        rows[name] = (F0, k, c, xs, ts, Fp)
        print(f" {name:15s}  {F0 / 1e3:7.1f}  {k:9.2f}  {c:12.3f}  {xs:10.3f}  {ts:10.4f}  {Fp / 1e3:10.1f}  "
              f"{xp:8.2f}  {Fp / (m * G0):6.3f}  {efficiency(v0, m, xs, Fp):8.3f}")

    print("\nindependent check, solve_ivp (rtol 1e-10) against work-energy and impulse-momentum:")
    print(" design           x_stop sim-WE (m)  t_stop sim-IM (s)  F_peak sim-WE (N)  impulse/(m v0)-1  energy balance/KE0-1")
    for name, (F0, k, c, xs, ts, Fp) in rows.items():
        sol = simulate_arrestment(v0, m, F0, k, c)
        au = energy_audit(sol, m, F0, k, c)
        bal = (au["W_F0"] + au["E_k"] + au["W_damper"]) / au["KE0"] - 1
        print(f" {name:15s}  {au['x_stop'] - xs:16.1e}  {au['t_stop'] - ts:17.1e}  {au['F_peak'] - Fp:17.1e}  "
              f"{au['impulse'] / (m * v0) - 1:16.1e}  {bal:19.1e}")
        print(f"   energy split: F0 {au['W_F0'] / au['KE0'] * 100:.1f} %, damper {au['W_damper'] / au['KE0'] * 100:.1f} %, "
              f"k-term stored {au['E_k'] / au['KE0'] * 100:.1f} %")
    F0, k, c, xs, ts, Fp = rows["damper-heavy"]
    print(f"\ndamper-heavy averages: by distance {energy_average_force(v0, m, xs) / 1e3:.1f} kN, "
          f"by time {momentum_average_force(v0, m, ts) / 1e3:.1f} kN, peak {Fp / 1e3:.1f} kN")

    print("\noff-design: constant-force gear set for 20 t at 65 m/s, other aircraft")
    print("   mass (t)  v0 (m/s)   x_stop (m)   peak g   peak at      | brake set for 20 t: x_stop (m)  g")
    for mm, vv in [(16000.0, 65.0), (20000.0, 65.0), (24000.0, 65.0), (20000.0, 60.0), (20000.0, 70.0), (24000.0, 70.0)]:
        xs = stopping_distance(vv, mm, 0.0, k0, c0)
        Fp, xp = peak_force(vv, mm, 0.0, k0, c0)
        Fb = m * a_ideal
        print(f"   {mm / 1e3:7.0f}  {vv:8.0f}   {xs:10.2f}   {Fp / (mm * G0):6.3f}   {'engagement' if xp == 0 else 'stop':10s}   | "
              f"{stopping_distance(vv, mm, Fb, 0.0, 0.0):27.2f}  {Fb / (mm * G0):.3f}")

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    for name, (F0, k, c, xs, ts, Fp) in rows.items():
        x = np.linspace(0.0, xs, 400)
        ax[0].plot(x, force_along_stroke(x, v0, m, F0, k, c) / (m * G0), label=name)
        sol = simulate_arrestment(v0, m, F0, k, c)
        tt = np.linspace(0.0, sol.t[-1], 400)
        xx, vv = sol.sol(tt)
        ax[1].plot(tt, arrestor_force(xx, np.clip(vv, 0, None), F0, k, c) / (m * G0), label=name)
        ax[2].plot(x, np.sqrt(np.clip(v_squared(x, v0, m, F0, k, c), 0, None)), label=f"{name}, work-energy")
        ax[2].plot(xx[::25], vv[::25], "k.", ms=3)
    ax[0].set_xlabel("payout x (m)")
    ax[0].set_ylabel("deceleration (g)")
    ax[0].set_title("work-energy: force along the stroke")
    ax[1].set_xlabel("t (s)")
    ax[1].set_ylabel("deceleration (g)")
    ax[1].set_title("solve_ivp: force against time")
    ax[2].set_xlabel("payout x (m)")
    ax[2].set_ylabel("speed (m/s)")
    ax[2].set_title("v(x): closed form (lines), solve_ivp (dots)")
    for a in ax:
        a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "arrestment.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'arrestment.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
