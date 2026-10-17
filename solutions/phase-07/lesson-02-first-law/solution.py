"""P07.02 The First Law: reference solution.

Blowdown of a pressurant tank (nitrogen or helium) through a choked orifice, written as an unsteady
control-volume energy balance and integrated in time:

    dm/dt = -mdot,     dU/dt = -mdot h + Qdot,

with the orifice mass flow choked (sonic throat) while p_tank / p_back exceeds ((gamma + 1)/2)^(gamma/(gamma - 1))
and subsonic after that. Three heat-transfer cases: adiabatic (Qdot = 0), isothermal (Qdot = mdot R T,
which holds T fixed) and a finite wall conductance Qdot = hA (T_wall - T) between the two.

Independent checks built into the module:
* the energy balance U(t) + integral(mdot h dt) - integral(Qdot dt) = U(0), carried as extra ODE states
  and closing to round-off (a linear invariant, which Runge-Kutta methods preserve exactly);
* the gas left in an adiabatic tank follows p ~ rho^gamma (it has expanded reversibly and adiabatically),
  and its specific internal-energy drop equals the closed-system boundary work integral(p dv);
* the closed-form choked adiabatic blowdown p/p0 = (1 + (gamma - 1) t / (2 tau))^(-2 gamma/(gamma - 1)),
  and the isothermal limit p/p0 = exp(-t / tau), with tau = V / (Cd A Gamma a0).

Model: ideal gas with constant cp and cv (gamma = 1.4 for N2, 5/3 for He; lesson P07.01 shows N2 has
cp/R = 3.503 at 300 K), uniform tank state, no kinetic energy inside the tank, an isentropic orifice with a
discharge coefficient, constant back pressure, constant wall temperature. Pressure vessels are NEVER built,
filled or vented at home in this course: this is a computational model only.

Units: SI (Pa, K, kg, m^3, J, s).

Run:  python solution.py          (prints the verification tables, writes tank_blowdown.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-07/lesson-02-first-law
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

R_U = 8.31446261815324            # J/(mol K), CODATA 2018

# Molar masses (kg/mol) from the NASA CEA thermo.inp database (P07.01); gamma from equipartition:
# 7/5 for N2 (vibration frozen near 300 K), 5/3 for monatomic He.
GASES = {
    "N2": {"M": 28.0134e-3, "gamma": 1.4},
    "He": {"M": 4.002602e-3, "gamma": 5.0 / 3.0},
}


@dataclass
class Tank:
    """A rigid tank: volume V [m^3], initial pressure p0 [Pa] and temperature T0 [K], gas name."""
    V: float
    p0: float
    T0: float
    gas: str = "N2"


@dataclass
class Orifice:
    """A sharp-edged orifice: discharge coefficient Cd [-] and throat area A [m^2]."""
    Cd: float
    A: float


# Typical, illustrative case: a 50 L nitrogen COPV at 30 MPa and 300 K discharging through a 2 mm
# orifice (Cd = 0.8) into a propellant-tank ullage held at 0.5 MPa; 3 MPa is taken as the lowest tank
# pressure the downstream regulator accepts.
CASE_TANK = Tank(V=0.050, p0=30.0e6, T0=300.0, gas="N2")
CASE_ORIFICE = Orifice(Cd=0.8, A=np.pi * 0.002 ** 2 / 4)
CASE_P_BACK = 0.5e6
CASE_P_STOP = 3.0e6


# ---------------------------------------------------------------- given helpers

def gas_properties(gas: str) -> tuple[float, float, float, float]:
    """(R, gamma, cv, cp) per unit mass for the named gas, constant-property model."""
    R = R_U / GASES[gas]["M"]
    g = GASES[gas]["gamma"]
    cv = R / (g - 1.0)
    return R, g, cv, cv + R


def initial_mass(tank: Tank) -> float:
    """m0 = p0 V / (R T0)."""
    R = gas_properties(tank.gas)[0]
    return tank.p0 * tank.V / (R * tank.T0)


# ---------------------------------------------------------------- orifice flow

def critical_pressure_ratio(gamma: float) -> float:
    """p*/p0 = (2/(gamma + 1))^(gamma/(gamma - 1)): the back-to-upstream pressure ratio below which the
    throat is sonic and the flow is choked (0.528 for gamma = 1.4, 0.487 for 5/3)."""
    return (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))


def choked_flow_function(gamma: float) -> float:
    """Gamma(gamma) = (2/(gamma + 1))^((gamma + 1)/(2 (gamma - 1))), so that mdot = Cd A Gamma p0 sqrt(gamma/(R T0))."""
    return (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))


def mass_flow(p0: float, T0: float, p_back: float, Cd: float, A: float, R: float, gamma: float) -> float:
    """Mass flow [kg/s] through an orifice from stagnation state (p0, T0) to back pressure p_back:
    choked, mdot = Cd A Gamma p0 sqrt(gamma / (R T0)), while p_back/p0 <= the critical ratio; otherwise
    the isentropic subsonic nozzle formula
    mdot = Cd A p0 sqrt(2 gamma / ((gamma - 1) R T0) [r^(2/gamma) - r^((gamma + 1)/gamma)]),  r = p_back/p0.
    Zero when p0 <= p_back (no reverse flow is modelled)."""
    if p0 <= p_back:
        return 0.0
    r = p_back / p0
    if r <= critical_pressure_ratio(gamma):
        return Cd * A * choked_flow_function(gamma) * p0 * np.sqrt(gamma / (R * T0))
    term = r ** (2.0 / gamma) - r ** ((gamma + 1.0) / gamma)
    return Cd * A * p0 * np.sqrt(2.0 * gamma / ((gamma - 1.0) * R * T0) * max(term, 0.0))


# ---------------------------------------------------------------- the control-volume model

def blowdown_rhs(t: float, y: np.ndarray, tank: Tank, orifice: Orifice, p_back: float,
                 mode: str, hA: float = 0.0, T_wall: float | None = None) -> list[float]:
    """Right-hand side of the tank ODEs with state y = [m, U, H_out, Q_in]:
    dm/dt = -mdot; dU/dt = -mdot h + Qdot; dH_out/dt = mdot h; dQ_in/dt = Qdot.
    mode: "adiabatic" (Qdot = 0), "isothermal" (Qdot = mdot (h - u) = mdot R T, which makes dU/dt = u dm/dt)
    or "wall" (Qdot = hA (T_wall - T), T_wall defaulting to T0). h = cp T and u = cv T (zero at 0 K)."""
    R, g, cv, cp = gas_properties(tank.gas)
    m, U = y[0], y[1]
    T = U / (m * cv)
    p = m * R * T / tank.V
    mdot = mass_flow(p, T, p_back, orifice.Cd, orifice.A, R, g)
    if mode == "adiabatic":
        Q = 0.0
    elif mode == "isothermal":
        Q = mdot * R * T
    elif mode == "wall":
        Q = hA * ((tank.T0 if T_wall is None else T_wall) - T)
    else:
        raise ValueError(f"unknown mode {mode!r}")
    h = cp * T
    return [-mdot, -mdot * h + Q, mdot * h, Q]


def simulate_blowdown(tank: Tank, orifice: Orifice, p_back: float, mode: str = "adiabatic",
                      t_end: float = 600.0, p_stop: float | None = None, hA: float = 0.0,
                      T_wall: float | None = None, rtol: float = 1e-10, n_out: int = 400) -> dict[str, np.ndarray]:
    """Integrate the blowdown with scipy's DOP853 from t = 0 to t_end, or until the tank pressure falls to
    p_stop (a terminal event). Returns a dict of arrays: t, m, U, H_out, Q_in, T, p, and the event time
    't_stop' (NaN if p_stop was not reached). Output is at the integrator's own steps plus n_out
    evenly spaced times (dense output), merged and sorted."""
    R, g, cv, cp = gas_properties(tank.gas)
    m0 = initial_mass(tank)
    y0 = [m0, m0 * cv * tank.T0, 0.0, 0.0]
    events = None
    if p_stop is not None:
        def hit(t, y, *args):
            return y[1] * (g - 1.0) / tank.V - p_stop          # p = (gamma - 1) U / V for constant cv
        hit.terminal, hit.direction = True, -1
        events = hit
    sol = solve_ivp(blowdown_rhs, (0.0, t_end), y0, method="DOP853", rtol=rtol,
                    atol=[1e-12 * m0, 1e-12 * y0[1], 1e-12 * y0[1], 1e-12 * y0[1]],
                    args=(tank, orifice, p_back, mode, hA, T_wall), events=events, dense_output=True)
    t_final = sol.t[-1]
    tt = np.unique(np.concatenate([sol.t, np.linspace(0.0, t_final, n_out)]))
    Y = sol.sol(tt)
    Y[:, -1] = sol.y[:, -1]
    Y[:, np.searchsorted(tt, sol.t)] = sol.y                  # use the integrator's own values at its steps
    m, U, H, Q = Y
    T = U / (m * cv)
    p = m * R * T / tank.V
    t_stop = float(sol.t_events[0][0]) if (p_stop is not None and len(sol.t_events[0])) else np.nan
    return {"t": tt, "m": m, "U": U, "H_out": H, "Q_in": Q, "T": T, "p": p, "t_stop": t_stop,
            "steps_t": sol.t, "steps_y": sol.y}


def energy_residual(res: dict[str, np.ndarray]) -> float:
    """Largest |U(t) + H_out(t) - Q_in(t) - U(0)| / U(0) at the integrator's steps: the first law for the
    tank as a control volume, closed over the whole run."""
    m, U, H, Q = res["steps_y"]
    return float(np.max(np.abs(U + H - Q - U[0])) / U[0])


# ---------------------------------------------------------------- closed forms and limits

def blowdown_time_constant(tank: Tank, orifice: Orifice) -> float:
    """tau = V / (Cd A Gamma a0), a0 = sqrt(gamma R T0): the e-folding time of the isothermal choked blowdown."""
    R, g, cv, cp = gas_properties(tank.gas)
    return tank.V / (orifice.Cd * orifice.A * choked_flow_function(g) * np.sqrt(g * R * tank.T0))


def adiabatic_choked_closed_form(t, tank: Tank, orifice: Orifice) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(p, T, m) of the adiabatic choked blowdown: rho/rho0 = (1 + (gamma - 1) t / (2 tau))^(-2/(gamma - 1)),
    p/p0 = (rho/rho0)^gamma, T/T0 = (rho/rho0)^(gamma - 1). Valid while the flow stays choked."""
    R, g, cv, cp = gas_properties(tank.gas)
    tau = blowdown_time_constant(tank, orifice)
    r = (1.0 + 0.5 * (g - 1.0) * np.asarray(t, float) / tau) ** (-2.0 / (g - 1.0))
    return tank.p0 * r ** g, tank.T0 * r ** (g - 1.0), initial_mass(tank) * r


def isothermal_choked_closed_form(t, tank: Tank, orifice: Orifice) -> tuple[np.ndarray, np.ndarray]:
    """(p, m) of the isothermal choked blowdown: p/p0 = m/m0 = exp(-t / tau)."""
    e = np.exp(-np.asarray(t, float) / blowdown_time_constant(tank, orifice))
    return tank.p0 * e, initial_mass(tank) * e


def adiabatic_time_to_pressure(p: float, tank: Tank, orifice: Orifice) -> float:
    """Time for the adiabatic choked blowdown to reach pressure p (inverting the closed form)."""
    g = gas_properties(tank.gas)[1]
    tau = blowdown_time_constant(tank, orifice)
    return 2.0 * tau / (g - 1.0) * ((p / tank.p0) ** (-(g - 1.0) / (2.0 * g)) - 1.0)


def remaining_mass_fraction(p_ratio: float, gamma: float, mode: str) -> float:
    """m/m0 left in the tank when its pressure has fallen to p_ratio p0: (p/p0)^(1/gamma) for an adiabatic
    blowdown (the remaining gas has expanded isentropically and is cold, hence dense), p/p0 for an
    isothermal one."""
    if mode == "adiabatic":
        return p_ratio ** (1.0 / gamma)
    if mode == "isothermal":
        return p_ratio
    raise ValueError(mode)


def polytropic_work(p1: float, v1: float, p2: float, v2: float, n: float) -> float:
    """Closed-system boundary work done BY a unit mass along p v^n = const, w = integral p dv [J/kg]:
    (p2 v2 - p1 v1)/(1 - n) for n != 1 and p1 v1 ln(v2/v1) for n = 1 (isothermal ideal gas)."""
    if abs(n - 1.0) < 1e-12:
        return p1 * v1 * np.log(v2 / v1)
    return (p2 * v2 - p1 * v1) / (1.0 - n)


def recovery_pressure(p_end: float, T_end: float, T_wall: float) -> float:
    """Pressure after the valve closes and the cold gas warms back to T_wall at constant volume (ideal gas)."""
    return p_end * T_wall / T_end


# ---------------------------------------------------------------- driver

def main(out_dir: str | Path = ".") -> None:
    """Print the verification tables and the engineering answer; write tank_blowdown.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    tank, orf, pb, pstop = CASE_TANK, CASE_ORIFICE, CASE_P_BACK, CASE_P_STOP
    R, g, cv, cp = gas_properties(tank.gas)
    m0 = initial_mass(tank)
    tau = blowdown_time_constant(tank, orf)
    mdot0 = mass_flow(tank.p0, tank.T0, pb, orf.Cd, orf.A, R, g)
    print("COMPUTATIONAL MODEL ONLY: no pressure vessel is built, filled or vented at home in this course.")
    print(f"N2 tank V = {tank.V*1e3:.0f} L, p0 = {tank.p0/1e6:.1f} MPa, T0 = {tank.T0:.0f} K: m0 = {m0:.4f} kg; "
          f"orifice d = 2 mm, Cd = {orf.Cd}; back pressure {pb/1e6:.1f} MPa")
    print(f"critical pressure ratio {critical_pressure_ratio(g):.4f}: choked while p > {pb/critical_pressure_ratio(g)/1e6:.3f} MPa; "
          f"initial mdot = {mdot0*1e3:.2f} g/s; tau = V/(Cd A Gamma a0) = {tau:.2f} s")

    out = {}
    for mode in ("adiabatic", "isothermal"):
        res = simulate_blowdown(tank, orf, pb, mode, t_end=5000.0, p_stop=pstop)
        out[mode] = res
        print(f"\n{mode}: reaches {pstop/1e6:.0f} MPa at t = {res['t_stop']:.2f} s, T = {res['T'][-1]:.2f} K, "
              f"m = {res['m'][-1]:.4f} kg ({100*res['m'][-1]/m0:.2f} % of m0 left), heat in {res['Q_in'][-1]/1e3:.2f} kJ")
        print(f"  energy balance residual max |U + int(mdot h) - int(Q) - U0| / U0 = {energy_residual(res):.1e}")
    ra = out["adiabatic"]
    pc, Tc, mc = adiabatic_choked_closed_form(ra["t"], tank, orf)
    print(f"\nadiabatic closed form: t to {pstop/1e6:.0f} MPa = {adiabatic_time_to_pressure(pstop, tank, orf):.2f} s; "
          f"max |p_num/p_closed - 1| = {np.max(np.abs(ra['p']/pc - 1)):.1e}")
    rho = ra["m"] / tank.V
    print(f"  p / rho^gamma varies by {np.ptp(ra['p']/rho**g)/(ra['p'][0]/rho[0]**g):.1e} (isentropic remaining gas)")
    v0, v = tank.V / m0, tank.V / ra["m"][-1]
    w = polytropic_work(tank.p0, v0, ra["p"][-1], v, g)
    print(f"  remaining gas, closed system: u0 - u = {cv*(tank.T0-ra['T'][-1])/1e3:.3f} kJ/kg, "
          f"boundary work int p dv = {w/1e3:.3f} kJ/kg")
    ri = out["isothermal"]
    pi, mi = isothermal_choked_closed_form(ri["t"], tank, orf)
    print(f"isothermal closed form: t to {pstop/1e6:.0f} MPa = tau ln(p0/p) = {tau*np.log(tank.p0/pstop):.2f} s; "
          f"max |p_num/p_closed - 1| = {np.max(np.abs(ri['p']/pi - 1)):.1e}; "
          f"heat = R T0 (m0 - m) = {R*tank.T0*(m0-ri['m'][-1])/1e3:.2f} kJ")
    print(f"\nstranded mass at {pstop/1e6:.0f} MPa: adiabatic {remaining_mass_fraction(pstop/tank.p0, g, 'adiabatic'):.4f} m0, "
          f"isothermal {remaining_mass_fraction(pstop/tank.p0, g, 'isothermal'):.4f} m0")
    print(f"after the valve closes the cold adiabatic gas warms back to {tank.T0:.0f} K: p rises from {pstop/1e6:.1f} to "
          f"{recovery_pressure(pstop, ra['T'][-1], tank.T0)/1e6:.2f} MPa")

    print("\nwall heat transfer between the limits (T_wall = 300 K):")
    print(f"{'hA (W/K)':>9} {'t to 3 MPa (s)':>15} {'T at 3 MPa (K)':>15} {'m left (kg)':>12} {'residual':>9}")
    for hA in (0.0, 1.0, 3.0, 10.0, 30.0, 100.0):
        r = simulate_blowdown(tank, orf, pb, "wall", t_end=5000.0, p_stop=pstop, hA=hA)
        print(f"{hA:>9.0f} {r['t_stop']:>15.2f} {r['T'][-1]:>15.2f} {r['m'][-1]:>12.4f} {energy_residual(r):>9.1e}")

    he = Tank(V=tank.V, p0=tank.p0, T0=tank.T0, gas="He")
    rh = simulate_blowdown(he, orf, pb, "adiabatic", t_end=5000.0, p_stop=pstop)
    print(f"\nhelium, same tank and orifice, adiabatic: m0 = {initial_mass(he):.4f} kg, tau = {blowdown_time_constant(he, orf):.2f} s, "
          f"reaches {pstop/1e6:.0f} MPa at {rh['t_stop']:.2f} s and {rh['T'][-1]:.2f} K")

    full = simulate_blowdown(tank, orf, pb, "adiabatic", t_end=600.0)
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
    tt = np.linspace(0.0, 120.0, 300)
    axs[0].plot(full["t"], full["p"] / 1e6, label="adiabatic, ODE")
    axs[0].plot(tt, adiabatic_choked_closed_form(tt, tank, orf)[0] / 1e6, "k:", label="adiabatic, closed form")
    rif = simulate_blowdown(tank, orf, pb, "isothermal", t_end=600.0)
    axs[0].plot(rif["t"], rif["p"] / 1e6, label="isothermal, ODE")
    axs[0].axhline(pb / critical_pressure_ratio(g) / 1e6, color="grey", ls="--", lw=0.8, label="unchoking pressure")
    axs[0].set_yscale("log")
    axs[0].set_xlim(0, 600)
    axs[0].set_xlabel("t (s)")
    axs[0].set_ylabel("tank pressure (MPa)")
    axs[0].set_title("N2 blowdown (computed)")
    axs[0].legend(fontsize=8)
    axs[1].plot(full["t"], full["T"], label="adiabatic")
    for hA in (3.0, 30.0):
        r = simulate_blowdown(tank, orf, pb, "wall", t_end=600.0, hA=hA)
        axs[1].plot(r["t"], r["T"], label=f"wall hA = {hA:.0f} W/K")
    axs[1].set_xlabel("t (s)")
    axs[1].set_ylabel("gas temperature (K)")
    axs[1].set_title("gas temperature")
    axs[1].legend(fontsize=8)
    m, U, H, Q = full["steps_y"]
    axs[2].semilogy(full["steps_t"][1:], np.abs(U + H - Q - U[0])[1:] / U[0] + 1e-18)
    axs[2].set_xlabel("t (s)")
    axs[2].set_ylabel("|U + int mdot h - U0| / U0")
    axs[2].set_title("first-law residual (adiabatic)")
    fig.tight_layout()
    fig.savefig(out_dir / "tank_blowdown.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
