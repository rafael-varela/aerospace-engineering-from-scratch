"""P02.01 First-Order ODEs: starter code.

A rocket's mass, velocity and altitude during a burn and the coast that follows, modelled with
first-order ordinary differential equations: separable (mass depletion, the vacuum rocket equation,
quadratic-drag coast), linear with an integrating factor (linear-drag coast), and autonomous
(terminal velocity as a stable equilibrium). Closed forms are checked against Euler and RK4
steppers written from scratch and against ``scipy.integrate.solve_ivp`` as the oracle.

Conventions:
    * State vector y = (h, v, m): altitude [m, positive up], velocity [m/s, positive up], mass [kg].
    * Thrust T = mdot * ve acts while propellant remains (m > m_dry); gravity g is constant;
      drag D = 0.5 rho(h) C_D A v |v| with an exponential atmosphere rho = rho0 exp(-h / H).
    * The "vacuum" closed forms drop drag and keep gravity.

Complete every function marked NotImplementedError.

Run:  python main.py          (prints the sounding-rocket example, writes rocket_ascent.png)
Test: pytest curriculum/phase-02/lesson-01-first-order-odes
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

RHS = Callable[[float, np.ndarray], np.ndarray]

RHO0 = 1.225          # sea-level density [kg/m^3]
SCALE_HEIGHT = 8500.0  # exponential-atmosphere scale height [m]


def sounding_rocket() -> dict[str, float]:
    """Return the parameters of a single-stage sounding rocket (SI units).

    250 kg at lift-off, 100 kg dry, 7.5 kg/s mass flow at 2200 m/s exhaust velocity (16.5 kN of
    thrust, thrust-to-weight 6.7), 0.30 m diameter, drag coefficient 0.5 held constant.
    """
    return {
        "m0": 250.0, "m_dry": 100.0, "mdot": 7.5, "ve": 2200.0,
        "cd": 0.5, "area": float(np.pi * 0.15 ** 2), "g": 9.81,
    }


# ---------------------------------------------------------------- separable equations: mass and the vacuum burn

def burn_time(m0: float, m_dry: float, mdot: float) -> float:
    """Return t_b = (m0 - m_dry) / mdot [s], the time at which the propellant is gone.

    Raise ``ValueError`` if mdot <= 0 or m_dry >= m0.
    """
    raise NotImplementedError


def mass_at(t: np.ndarray | float, m0: float, m_dry: float, mdot: float) -> np.ndarray:
    """Return m(t) = m0 - mdot t during the burn and m_dry after it (solution of dm/dt = -mdot).

    Raise ``ValueError`` for negative times.
    """
    raise NotImplementedError


def tsiolkovsky_delta_v(ve: float, m0: float, mf: float) -> float:
    """Return the ideal velocity gain ve ln(m0 / mf) [m/s] of a burn from mass m0 to mass mf, no gravity, no drag."""
    raise NotImplementedError


def vacuum_velocity(t: np.ndarray | float, p: dict[str, float]) -> np.ndarray:
    """Return v(t) = ve ln(m0 / m(t)) - g t during the burn: the vertical rocket equation without drag.

    Valid for 0 <= t <= t_b; raise ``ValueError`` outside that interval (after burnout the thrust
    term is absent and the formula no longer applies; before lift-off there is no solution).
    """
    raise NotImplementedError


def vacuum_altitude(t: np.ndarray | float, p: dict[str, float]) -> np.ndarray:
    """Return h(t) = ve t - (ve m(t) / mdot) ln(m0 / m(t)) - g t^2 / 2 during the burn (integral of ``vacuum_velocity``)."""
    raise NotImplementedError


def apogee_from_burnout(h_b: float, v_b: float, g: float) -> float:
    """Return the vacuum apogee h_b + v_b^2 / (2 g) [m] reached after burnout at (h_b, v_b)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the full model: drag, atmosphere, right-hand side

def air_density(h: np.ndarray | float) -> np.ndarray:
    """Return rho(h) = rho0 exp(-h / H) [kg/m^3] with rho0 = 1.225 and H = 8500 m (exponential atmosphere)."""
    raise NotImplementedError


def rocket_rhs(t: float, y: np.ndarray, p: dict[str, float], burning: bool = True,
               drag: bool = True) -> np.ndarray:
    """Return dy/dt for y = (h, v, m): (v, (T - D) / m - g, -mdot) while ``burning``, (v, -D / m - g, 0) after.

    Thrust T = mdot ve. Drag D = 0.5 rho(h) cd area v |v| when ``drag`` is True, else 0. The burn/coast
    switch is a separate flag rather than a test on m so that each phase is a smooth right-hand side:
    a stepper that straddles the switch-off loses its order of accuracy.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- steppers from scratch

def euler_step(f: RHS, t: float, y: np.ndarray, h: float) -> np.ndarray:
    """Return y + h f(t, y): one explicit Euler step (first order)."""
    raise NotImplementedError


def rk4_step(f: RHS, t: float, y: np.ndarray, h: float) -> np.ndarray:
    """Return one classical Runge-Kutta step (fourth order): y + h (k1 + 2 k2 + 2 k3 + k4) / 6."""
    raise NotImplementedError


def integrate_fixed_step(f: RHS, t0: float, y0: np.ndarray, t_end: float, n_steps: int,
                         method: str = "rk4") -> tuple[np.ndarray, np.ndarray]:
    """Integrate dy/dt = f(t, y) from t0 to t_end with ``n_steps`` equal steps of ``euler_step`` or ``rk4_step``.

    Return (times of shape (n_steps + 1,), states of shape (n_steps + 1, len(y0))).
    Raise ``ValueError`` for an unknown method or n_steps < 1.
    """
    raise NotImplementedError


def simulate_rocket(p: dict[str, float], t_end: float, n_steps: int, method: str = "rk4",
                    drag: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """Integrate the rocket from rest on the pad, y0 = (0, 0, m0), to ``t_end`` with ``integrate_fixed_step``.

    The burn (0 to t_b) and the coast (t_b to t_end) are integrated as two phases with their own smooth
    right-hand sides, the step count split in proportion to their durations (at least one step each), so the
    thrust switch-off always falls on a grid point. If t_end <= t_b only the burn is integrated.
    """
    raise NotImplementedError


def convergence_study(method: str, ns: list[int], p: dict[str, float]) -> np.ndarray:
    """Return |v_numeric(t_b) - v_exact(t_b)| for each step count in ``ns``, integrating the drag-free burn to t_b."""
    raise NotImplementedError


# ---------------------------------------------------------------- coast phase: two more closed forms

def quadratic_drag_coast(v_b: float, k: float, g: float) -> tuple[float, float]:
    """Return (time to apex [s], height gained [m]) for an upward coast with dv/dt = -g - k v^2 (constant mass, constant density).

    Separable: v(t) = sqrt(g/k) tan(arctan(v_b sqrt(k/g)) - sqrt(g k) t); the apex is where v = 0 and the
    height gained is ln(1 + k v_b^2 / g) / (2 k). For k = 0 the vacuum values v_b / g and v_b^2 / (2 g) are returned.
    """
    raise NotImplementedError


def quadratic_drag_velocity(t: np.ndarray | float, v_b: float, k: float, g: float) -> np.ndarray:
    """Return v(t) = sqrt(g/k) tan(arctan(v_b sqrt(k/g)) - sqrt(g k) t) for the upward coast, valid until the apex."""
    raise NotImplementedError


def linear_drag_velocity(t: np.ndarray | float, v0: float, c: float, m: float, g: float) -> np.ndarray:
    """Return the solution of dv/dt + (c/m) v = -g, v(0) = v0 (constant mass, linear drag), by integrating factor.

    v(t) = (v0 + m g / c) exp(-c t / m) - m g / c. For c = 0 the free-fall line v0 - g t is returned.
    """
    raise NotImplementedError


def terminal_velocity(m: float, rho: float, cd: float, area: float, g: float) -> float:
    """Return sqrt(2 m g / (rho cd area)) [m/s]: the speed at which quadratic drag balances weight."""
    raise NotImplementedError


# ---------------------------------------------------------------- qualitative tools

def equilibria(f: Callable[[float], float], grid: np.ndarray) -> list[tuple[float, bool]]:
    """Return [(y*, stable), ...] for the autonomous equation dy/dt = f(y), one entry per sign change of f on ``grid``.

    Each root is located with ``scipy.optimize.brentq`` on the bracketing interval and is stable when the
    slope f'(y*) (central difference) is negative. Roots where the slope is zero are reported as unstable.
    """
    raise NotImplementedError


def direction_field(f: Callable[[float, float], float], t_grid: np.ndarray,
                    y_grid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return unit arrow components (U, V) of the slope field of dy/dt = f(t, y) on the mesh of ``t_grid`` x ``y_grid``.

    Each arrow points along (1, f(t, y)) and is normalised to unit length so the plot shows direction only.
    Shapes follow ``np.meshgrid(t_grid, y_grid)``: (len(y_grid), len(t_grid)).
    """
    raise NotImplementedError


# ---------------------------------------------------------------- the example

def main(out_dir: str | Path = ".") -> None:
    """Print the sounding-rocket burnout and apogee figures and save rocket_ascent.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.integrate import solve_ivp

    out_dir = Path(out_dir)
    p = sounding_rocket()
    tb = burn_time(p["m0"], p["m_dry"], p["mdot"])
    print(f"thrust {p['mdot'] * p['ve'] / 1e3:.1f} kN, T/W at lift-off {p['mdot'] * p['ve'] / (p['m0'] * p['g']):.2f},"
          f" burn time {tb:.1f} s")
    dv = tsiolkovsky_delta_v(p["ve"], p["m0"], p["m_dry"])
    vb, hb = float(vacuum_velocity(tb, p)), float(vacuum_altitude(tb, p))
    print(f"Tsiolkovsky delta-v {dv:.1f} m/s; gravity loss g t_b = {p['g'] * tb:.1f} m/s")
    print(f"vacuum burnout: v_b = {vb:.1f} m/s, h_b = {hb:.0f} m, apogee {apogee_from_burnout(hb, vb, p['g']) / 1e3:.2f} km")

    # Convergence of the two steppers on the drag-free burn.
    ns = [10, 20, 40, 80, 160]
    for method in ("euler", "rk4"):
        errs = convergence_study(method, ns, p)
        rates = np.log2(errs[:-1] / errs[1:])
        print(f"{method:5s} errors in v(t_b): {np.array2string(errs, precision=3)}  rates {np.round(rates, 2)}")

    # Full model with drag, RK4 against solve_ivp.
    t_end = 120.0
    n = 2400
    ts, ys = simulate_rocket(p, t_end, n, "rk4", drag=True)
    i_b = int(np.argmin(np.abs(ts - tb)))
    print(f"with drag (RK4, {n} steps): burnout v_b = {ys[i_b, 1]:.1f} m/s at h_b = {ys[i_b, 0]:.0f} m;"
          f" drag costs {vb - ys[i_b, 1]:.0f} m/s")
    i_apex = int(np.argmax(ys[:, 0]))
    print(f"apogee with drag {ys[i_apex, 0] / 1e3:.2f} km at t = {ts[i_apex]:.1f} s")
    burn = solve_ivp(lambda t, y: rocket_rhs(t, y, p, True, True), (0.0, tb), [0.0, 0.0, p["m0"]], rtol=1e-10, atol=1e-8)
    coast = solve_ivp(lambda t, y: rocket_rhs(t, y, p, False, True), (tb, t_end), burn.y[:, -1], rtol=1e-10, atol=1e-8,
                      dense_output=True)
    v_ivp = burn.y[1, -1]
    print(f"solve_ivp (RK45, rtol 1e-10): burnout v_b = {v_ivp:.1f} m/s, difference from RK4 {abs(v_ivp - ys[i_b, 1]):.2e} m/s")
    h_ivp = coast.sol(ts[ts >= tb])[0]
    print(f"max |h_RK4 - h_ivp| over the coast: {np.max(np.abs(ys[ts >= tb, 0] - h_ivp)):.2e} m")

    # Coast closed form with constant density (frozen at burnout altitude) against the variable-density model.
    rho_b = float(air_density(ys[i_b, 0]))
    k = rho_b * p["cd"] * p["area"] / (2.0 * p["m_dry"])
    t_apex, h_gain = quadratic_drag_coast(ys[i_b, 1], k, p["g"])
    print(f"coast with density frozen at {rho_b:.3f} kg/m^3: k = {k:.3e} 1/m, apex after {t_apex:.1f} s,"
          f" height gained {h_gain / 1e3:.2f} km (vacuum: {ys[i_b, 1] / p['g']:.1f} s, {ys[i_b, 1] ** 2 / 2 / p['g'] / 1e3:.2f} km)")
    print(f"full model gains {(ys[i_apex, 0] - ys[i_b, 0]) / 1e3:.2f} km in {ts[i_apex] - tb:.1f} s")
    vt = terminal_velocity(p["m_dry"], RHO0, p["cd"], p["area"], p["g"])
    print(f"terminal velocity of the empty rocket at sea level: {vt:.1f} m/s")
    f_auto = lambda v: -p["g"] - k * v * abs(v)
    print("equilibria of dv/dt = -g - k v|v| on [-800, 800]:", equilibria(f_auto, np.linspace(-800.0, 800.0, 161)))

    # Figure: trajectories, convergence, direction field of the vacuum rocket equation.
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))
    ts0, ys0 = simulate_rocket(p, t_end, n, "rk4", drag=False)
    axes[0].plot(ts, ys[:, 0] / 1e3, label="with drag (RK4)")
    axes[0].plot(ts0, ys0[:, 0] / 1e3, "--", label="no drag (RK4)")
    tt = np.linspace(0.0, tb, 50)
    axes[0].plot(tt, vacuum_altitude(tt, p) / 1e3, "k:", label="no drag (closed form)")
    axes[0].axvline(tb, color="grey", lw=0.6)
    axes[0].set_xlabel("t [s]"); axes[0].set_ylabel("altitude [km]"); axes[0].legend(fontsize=8); axes[0].grid(lw=0.4)
    for method, mk in (("euler", "o"), ("rk4", "s")):
        errs = convergence_study(method, ns, p)
        axes[1].loglog(ns, errs, mk + "-", label=method)
    axes[1].loglog(ns, errs[0] * (ns[0] / np.array(ns, float)) ** 4, "k:", lw=0.8, label="slope -4")
    axes[1].set_xlabel("steps over the burn"); axes[1].set_ylabel("|error in v(t_b)| [m/s]"); axes[1].legend(fontsize=8)
    axes[1].grid(True, which="both", lw=0.4)
    tg, vg = np.linspace(0.0, tb, 15), np.linspace(0.0, 2200.0, 15)
    U, V = direction_field(lambda t, v: (p["mdot"] * p["ve"] / (p["m0"] - p["mdot"] * t) - p["g"]) / 100.0, tg, vg)
    axes[2].quiver(tg, vg, U, V, angles="xy", pivot="mid", color="0.6", scale=25)
    axes[2].plot(tt, vacuum_velocity(tt, p), label="v(t), closed form")
    axes[2].set_xlabel("t [s]"); axes[2].set_ylabel("v [m/s]"); axes[2].legend(fontsize=8); axes[2].grid(lw=0.4)
    axes[2].set_title("direction field, vacuum burn", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "rocket_ascent.png", dpi=150)
    print(f"saved {out_dir / 'rocket_ascent.png'}")


if __name__ == "__main__":
    main()
