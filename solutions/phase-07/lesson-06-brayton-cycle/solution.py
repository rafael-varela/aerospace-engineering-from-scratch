"""P07.06 Power Cycles: Brayton: reference solution.

What efficiency and specific work can a gas-turbine core with a pressure ratio of 12 and a turbine
inlet temperature of 1400 K deliver at best, what pressure ratio gives the most work, and how far is
it from Carnot? This module is an ideal (air-standard) Brayton cycle analyzer:

* the four-state table (T, p, v, h, s) for isentropic compression 1-2, constant-pressure heat
  addition 2-3, isentropic expansion 3-4 and constant-pressure heat rejection 4-1, for a
  constant-cp gas or the NASA variable-cp dry-air model;
* work, heat, thermal efficiency and back-work ratio, and their constant-cp closed forms
  eta = 1 - pi^(-(gamma-1)/gamma) and w_net = cp T1 [tau (1 - 1/x) - (x - 1)];
* the closed cycle path, discretised, and the net work as the enclosed area in the p-v and T-s
  planes by a numerical contour integral;
* the Carnot bound between the same extreme temperatures, and the pressure ratio that maximises
  specific work, (T3/T1)^(gamma/(2(gamma-1))) for constant cp, found also by numerical maximisation.

Property model: as in P07.05, NASA Glenn 9-coefficient polynomials for N2, O2, Ar and CO2 (NASA CEA
data/thermo.inp; McBride, Zehe & Gordon, NASA/TP-2002-211556) combined into the Phase 07 dry-air
model (mole fractions N2 0.7808, O2 0.2095, Ar 0.0093, CO2 0.0004); s° at p° = 1 bar; R = 8.31446261815324
J/(mol K) (CODATA 2018). The air-standard cycle replaces combustion by heat addition and the exhaust
by heat rejection to close the loop; the working fluid is air throughout.

Units: SI (K, Pa, m^3/kg, J/kg, J/(kg K)).

Run:  python solution.py [out_dir]      (prints the lesson's numbers, writes brayton_cycle.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-07/lesson-06-brayton-cycle
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.optimize import minimize_scalar

R_UNIVERSAL = 8.31446261815324          # J/(mol K), CODATA 2018 exact
P_STD = 1.0e5                           # Pa, reference pressure of the NASA s° data (1 bar)

# ---------------------------------------------------------------- property data (given)
# NASA Glenn 9-coefficient polynomials, copied unrounded from data/thermo.inp of NASA CEA
# (https://github.com/nasa/cea, Apache-2.0), documented in McBride, Zehe & Gordon, NASA/TP-2002-211556.
# Each interval: (T_low, T_high, [a1..a7], [b1, b2]).
#   cp/R  = a1 T^-2 + a2 T^-1 + a3 + a4 T + a5 T^2 + a6 T^3 + a7 T^4
#   h/RT  = -a1 T^-2 + a2 ln(T)/T + a3 + a4 T/2 + a5 T^2/3 + a6 T^3/4 + a7 T^4/5 + b1/T
#   s°/R  = -a1 T^-2/2 - a2 T^-1 + a3 ln T + a4 T + a5 T^2/2 + a6 T^3/3 + a7 T^4/4 + b2
NASA9 = {
    "N2": {"M": 28.0134e-3, "intervals": [
        (200.0, 1000.0, [22103.71497, -381.846182, 6.08273836, -0.00853091441, 1.384646189e-05, -9.62579362e-09, 2.519705809e-12], [710.846086, -10.76003744]),
        (1000.0, 6000.0, [587712.406, -2239.249073, 6.06694922, -0.00061396855, 1.491806679e-07, -1.923105485e-11, 1.061954386e-15], [12832.10415, -15.86640027])]},
    "O2": {"M": 31.9988e-3, "intervals": [
        (200.0, 1000.0, [-34255.6342, 484.700097, 1.119010961, 0.00429388924, -6.83630052e-07, -2.0233727e-09, 1.039040018e-12], [-3391.45487, 18.4969947]),
        (1000.0, 6000.0, [-1037939.022, 2344.830282, 1.819732036, 0.001267847582, -2.188067988e-07, 2.053719572e-11, -8.19346705e-16], [-16890.10929, 17.38716506])]},
    "Ar": {"M": 39.948e-3, "intervals": [
        (200.0, 1000.0, [0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0], [-745.375, 4.37967491]),
        (1000.0, 6000.0, [20.10538475, -0.0599266107, 2.500069401, -3.99214116e-08, 1.20527214e-11, -1.819015576e-15, 1.078576636e-19], [-744.993961, 4.37918011])]},
    "CO2": {"M": 44.0095e-3, "intervals": [
        (200.0, 1000.0, [49436.5054, -626.411601, 5.30172524, 0.002503813816, -2.127308728e-07, -7.68998878e-10, 2.849677801e-13], [-45281.9846, -7.04827944]),
        (1000.0, 6000.0, [117696.2419, -1788.791477, 8.29152319, -9.22315678e-05, 4.86367688e-09, -1.891053312e-12, 6.33003659e-16], [-39083.5059, -26.52669281])]},
}
AIR_X = {"N2": 0.7808, "O2": 0.2095, "Ar": 0.0093, "CO2": 0.0004}   # Phase 07 dry-air model (mole fractions)
M_AIR = sum(x * NASA9[k]["M"] for k, x in AIR_X.items())          # kg/mol, 28.965 g/mol
R_AIR = R_UNIVERSAL / M_AIR                                        # J/(kg K), 287.05


def _nasa_dimensionless(species: str, T: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """cp/R, h/(RT) and s°/R of one species at temperature(s) T (200-6000 K), picking the interval."""
    T = np.asarray(T, dtype=float)
    lo, hi = NASA9[species]["intervals"]
    out = []
    for (tl, th, a, b) in (lo, hi):
        a1, a2, a3, a4, a5, a6, a7 = a
        b1, b2 = b
        cp = a1 / T**2 + a2 / T + a3 + a4 * T + a5 * T**2 + a6 * T**3 + a7 * T**4
        h = -a1 / T**2 + a2 * np.log(T) / T + a3 + a4 * T / 2 + a5 * T**2 / 3 + a6 * T**3 / 4 + a7 * T**4 / 5 + b1 / T
        s = -a1 / T**2 / 2 - a2 / T + a3 * np.log(T) + a4 * T + a5 * T**2 / 2 + a6 * T**3 / 3 + a7 * T**4 / 4 + b2
        out.append((cp, h, s))
    use_lo = T <= 1000.0
    return tuple(np.where(use_lo, out[0][i], out[1][i]) for i in range(3))


def cp_air(T: np.ndarray) -> np.ndarray:
    """Specific heat at constant pressure of the dry-air model, J/(kg K)."""
    return sum(x * _nasa_dimensionless(k, T)[0] for k, x in AIR_X.items()) * R_AIR


def h_air(T: np.ndarray) -> np.ndarray:
    """Specific enthalpy of the dry-air model, J/kg (NASA reference: elements at 298.15 K; for air it is
    sensible enthalpy plus the small formation enthalpy of the CO2)."""
    T = np.asarray(T, dtype=float)
    return sum(x * _nasa_dimensionless(k, T)[1] for k, x in AIR_X.items()) * R_AIR * T


def s0_air(T: np.ndarray) -> np.ndarray:
    """Standard-state specific entropy s°(T) of the dry-air model at p° = 1 bar, J/(kg K), without the
    constant mixing term (it cancels in every difference this lesson takes)."""
    return sum(x * _nasa_dimensionless(k, T)[2] for k, x in AIR_X.items()) * R_AIR


@dataclass(frozen=True)
class Gas:
    """An ideal gas described by its gas constant and three property functions of T alone.
    s0 is the entropy at the reference pressure P_STD; s(T, p) = s0(T) - R ln(p / P_STD)."""
    name: str
    R: float
    cp: Callable[[np.ndarray], np.ndarray]
    h: Callable[[np.ndarray], np.ndarray]
    s0: Callable[[np.ndarray], np.ndarray]


AIR = Gas("dry air, NASA 9-coefficient", R_AIR, cp_air, h_air, s0_air)


def constant_cp_gas(gamma: float = 1.4, R: float = R_AIR) -> Gas:
    """A calorically perfect gas: cp = gamma R/(gamma - 1), h = cp (T - 298.15), s0 = cp ln(T/298.15)."""
    cp = gamma * R / (gamma - 1.0)
    return Gas(f"constant cp, gamma = {gamma}", R,
               lambda T: np.full_like(np.asarray(T, dtype=float), cp),
               lambda T: cp * (np.asarray(T, dtype=float) - 298.15),
               lambda T: cp * np.log(np.asarray(T, dtype=float) / 298.15))


# ---------------------------------------------------------------- inverting the property functions

def _newton_invert(f: Callable[[np.ndarray], np.ndarray], dfdT: Callable[[np.ndarray], np.ndarray],
                   target: np.ndarray, T_guess: float = 600.0) -> np.ndarray:
    """Solve f(T) = target for T (f monotonic increasing, dfdT > 0) by Newton's method on arrays,
    kept inside the 200-6000 K range of the data; converged when every step is below 1e-10 K."""
    target = np.asarray(target, dtype=float)
    T = np.full_like(target, T_guess)
    for _ in range(100):
        step = (f(T) - target) / dfdT(T)
        T = np.clip(T - step, 200.0, 6000.0)
        if np.all(np.abs(step) < 1e-10):
            return T
    raise RuntimeError("temperature inversion did not converge")


def temperature_from_s0(gas: Gas, s0_target: float | np.ndarray) -> float | np.ndarray:
    """The temperature at which gas.s0(T) = s0_target. s0 rises monotonically with T (ds0/dT = cp/T > 0),
    so Newton's method with that derivative converges from any start in 200-6000 K. Accepts arrays."""
    T = _newton_invert(gas.s0, lambda t: gas.cp(t) / t, s0_target)
    return float(T) if np.ndim(T) == 0 else T


def temperature_from_h(gas: Gas, h_target: float | np.ndarray) -> float | np.ndarray:
    """The temperature at which gas.h(T) = h_target (dh/dT = cp > 0), by Newton's method. Accepts arrays."""
    T = _newton_invert(gas.h, gas.cp, h_target)
    return float(T) if np.ndim(T) == 0 else T


# ---------------------------------------------------------------- the ideal Brayton cycle

def isentropic_temperature(gas: Gas, T1: float, pi: float) -> float:
    """Isentropic exit temperature for an ideal gas (given; built in P07.05): s0(T2) = s0(T1) + R ln(pi)."""
    return temperature_from_s0(gas, float(gas.s0(T1)) + gas.R * np.log(pi))


def specific_entropy(gas: Gas, T: np.ndarray, p: np.ndarray) -> np.ndarray:
    """Specific entropy s(T, p) = s0(T) - R ln(p / p°), J/(kg K), on the same reference as gas.s0."""
    return gas.s0(T) - gas.R * np.log(np.asarray(p, dtype=float) / P_STD)


def brayton_states(gas: Gas, T1: float, p1: float, pi: float, T3: float) -> dict[str, np.ndarray]:
    """State table of the ideal Brayton cycle, states 1-4 (index 0-3). 1: compressor inlet (T1, p1);
    2: after isentropic compression to pi p1; 3: after heat addition at constant pressure to T3;
    4: after isentropic expansion to p1. Returns arrays 'T' (K), 'p' (Pa), 'v' = RT/p (m^3/kg),
    'h' (J/kg) and 's' (J/(kg K))."""
    T2 = isentropic_temperature(gas, T1, pi)
    T4 = isentropic_temperature(gas, T3, 1.0 / pi)
    T = np.array([T1, T2, T3, T4], dtype=float)
    p = np.array([p1, pi * p1, pi * p1, p1], dtype=float)
    return {"T": T, "p": p, "v": gas.R * T / p, "h": gas.h(T), "s": specific_entropy(gas, T, p)}


def brayton_performance(gas: Gas, T1: float, p1: float, pi: float, T3: float) -> dict[str, float]:
    """Per unit mass of air: compressor work w_c = h2 - h1, turbine work w_t = h3 - h4, heat in
    q_in = h3 - h2, heat out q_out = h4 - h1, net work w_net = w_t - w_c, thermal efficiency
    eta = w_net / q_in, back-work ratio bwr = w_c / w_t, and the Carnot efficiency between T1 and T3."""
    h = brayton_states(gas, T1, p1, pi, T3)["h"]
    w_c, w_t = h[1] - h[0], h[2] - h[3]
    q_in, q_out = h[2] - h[1], h[3] - h[0]
    w_net = w_t - w_c
    return {"w_c": float(w_c), "w_t": float(w_t), "q_in": float(q_in), "q_out": float(q_out),
            "w_net": float(w_net), "eta": float(w_net / q_in), "bwr": float(w_c / w_t),
            "eta_carnot": carnot_efficiency(T1, T3)}


def efficiency_closed_form(pi: float | np.ndarray, gamma: float = 1.4) -> float | np.ndarray:
    """Constant-cp ideal Brayton efficiency eta = 1 - pi^(-(gamma-1)/gamma) = 1 - T1/T2."""
    return 1.0 - np.asarray(pi, dtype=float) ** (-(gamma - 1.0) / gamma)


def net_work_closed_form(T1: float, T3: float, pi: float | np.ndarray, gamma: float = 1.4,
                         R: float = R_AIR) -> float | np.ndarray:
    """Constant-cp ideal Brayton net work per unit mass, w_net = cp T1 [tau (1 - 1/x) - (x - 1)],
    with tau = T3/T1 and x = pi^((gamma-1)/gamma)."""
    cp = gamma * R / (gamma - 1.0)
    x = np.asarray(pi, dtype=float) ** ((gamma - 1.0) / gamma)
    return cp * T1 * ((T3 / T1) * (1.0 - 1.0 / x) - (x - 1.0))


def carnot_efficiency(T_cold: float, T_hot: float) -> float:
    """Carnot efficiency 1 - T_cold/T_hot between two reservoirs (temperatures in K)."""
    return 1.0 - T_cold / T_hot


def optimum_pressure_ratio_closed_form(T1: float, T3: float, gamma: float = 1.4) -> float:
    """Constant-cp pressure ratio of maximum net work, (T3/T1)^(gamma/(2(gamma-1))); at it T2 = T4 = sqrt(T1 T3)."""
    return (T3 / T1) ** (gamma / (2.0 * (gamma - 1.0)))


def optimum_pressure_ratio(gas: Gas, T1: float, p1: float, T3: float) -> float:
    """Pressure ratio of maximum net work found numerically: maximise brayton_performance(...)['w_net']
    over ln(pi) with scipy.optimize.minimize_scalar (bounded, between pi = 1.01 and the pi at which
    T2 reaches T3), to a tolerance of 1e-10 in ln(pi)."""
    ln_max = (float(gas.s0(T3)) - float(gas.s0(T1))) / gas.R          # pi at which T2s = T3
    res = minimize_scalar(lambda x: -brayton_performance(gas, T1, p1, float(np.exp(x)), T3)["w_net"],
                          bounds=(np.log(1.01), ln_max - 1e-6), method="bounded", options={"xatol": 1e-10})
    return float(np.exp(res.x))


def cycle_path(gas: Gas, T1: float, p1: float, pi: float, T3: float, n: int = 200) -> dict[str, np.ndarray]:
    """The closed cycle 1-2-3-4-1 as 4n + 1 points (the last repeats the first), n per process.
    Isentropic legs (1-2, 3-4): p geometric in n steps, T from s(T, p) = const. Isobaric legs (2-3, 4-1):
    T linear in n steps at constant p. Returns 'p', 'v', 'T', 's' arrays along the loop."""
    st = brayton_states(gas, T1, p1, pi, T3)
    T, p, s = st["T"], st["p"], st["s"]
    u = np.linspace(0.0, 1.0, n + 1)[:-1]
    legs_T, legs_p = [], []
    # 1-2 and 3-4: isentropic
    for (i, j) in ((0, 1), (2, 3)):
        pp = p[i] * (p[j] / p[i]) ** u
        TT = temperature_from_s0(gas, s[i] + gas.R * np.log(pp / P_STD))
        legs_p.append(pp)
        legs_T.append(TT)
        if i == 0:                                      # then 2-3 isobaric
            legs_p.append(np.full_like(u, p[1]))
            legs_T.append(T[1] + (T[2] - T[1]) * u)
    legs_p.append(np.full_like(u, p[3]))                # 4-1 isobaric
    legs_T.append(T[3] + (T[0] - T[3]) * u)
    P = np.concatenate(legs_p + [p[:1]])
    TT = np.concatenate(legs_T + [T[:1]])
    return {"p": P, "v": gas.R * TT / P, "T": TT, "s": specific_entropy(gas, TT, P)}


def enclosed_area(x: np.ndarray, y: np.ndarray) -> float:
    """Contour integral of y dx around a closed polygon (x[-1] = x[0]) by the trapezoid rule,
    sum of (y_k + y_{k+1})/2 (x_{k+1} - x_k). Positive for a clockwise loop in the (x, y) plane, which
    is a power cycle in both the p-v plane (integral of p dv) and the T-s plane (integral of T ds)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def observed_order(errors: np.ndarray, ratio: float = 2.0) -> np.ndarray:
    """Observed convergence orders log(e_k / e_{k+1}) / log(ratio) for refinements by ``ratio``."""
    e = np.asarray(errors, dtype=float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


# ---------------------------------------------------------------- the lesson's numbers

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save brayton_cycle.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    cold = constant_cp_gas(1.4)
    T1, p1, pi, T3 = 288.15, 101325.0, 12.0, 1400.0
    print(f"illustrative gas-turbine core (typical, illustrative): T1 = {T1} K, p1 = {p1:.0f} Pa, pi = {pi:.0f}, T3 = {T3:.0f} K")
    for gas in (cold, AIR):
        st = brayton_states(gas, T1, p1, pi, T3)
        pf = brayton_performance(gas, T1, p1, pi, T3)
        print(f"\n  {gas.name}")
        print("  state   T (K)     p (kPa)    v (m^3/kg)   h (kJ/kg)   s (J/(kg K))")
        for k in range(4):
            print(f"    {k + 1}   {st['T'][k]:8.2f}  {st['p'][k] / 1e3:9.2f}   {st['v'][k]:9.5f}   {st['h'][k] / 1e3:9.2f}   {st['s'][k]:9.2f}")
        print(f"  w_c = {pf['w_c'] / 1e3:.2f}, w_t = {pf['w_t'] / 1e3:.2f}, w_net = {pf['w_net'] / 1e3:.2f}, q_in = {pf['q_in'] / 1e3:.2f}, "
              f"q_out = {pf['q_out'] / 1e3:.2f} kJ/kg; eta = {pf['eta']:.4f}, back-work ratio = {pf['bwr']:.4f}, Carnot = {pf['eta_carnot']:.4f}")
        path = cycle_path(gas, T1, p1, pi, T3, 400)
        print(f"  enclosed area: p-v {enclosed_area(path['v'], path['p']) / 1e3:.4f} kJ/kg, T-s {enclosed_area(path['s'], path['T']) / 1e3:.4f} kJ/kg "
              f"(w_net {pf['w_net'] / 1e3:.4f}); first law q_in - q_out - w_net = {pf['q_in'] - pf['q_out'] - pf['w_net']:.2e} J/kg")
        pio = optimum_pressure_ratio(gas, T1, p1, T3)
        po = brayton_performance(gas, T1, p1, pio, T3)
        sto = brayton_states(gas, T1, p1, pio, T3)
        print(f"  optimum pi for work = {pio:.3f}: w_net = {po['w_net'] / 1e3:.2f} kJ/kg, eta = {po['eta']:.4f}, T2 = {sto['T'][1]:.2f} K, T4 = {sto['T'][3]:.2f} K")
    print(f"\n  closed forms, gamma 1.4: eta = {efficiency_closed_form(pi):.4f}, w_net = {net_work_closed_form(T1, T3, pi) / 1e3:.2f} kJ/kg, "
          f"pi_opt = {optimum_pressure_ratio_closed_form(T1, T3):.3f}, sqrt(T1 T3) = {np.sqrt(T1 * T3):.2f} K, "
          f"pi at zero work = {(T3 / T1) ** 3.5:.1f}")
    mdot = 5.0
    print(f"  at {mdot} kg/s (turboprop/APU scale, illustrative), ideal net power = {mdot * brayton_performance(AIR, T1, p1, pi, T3)['w_net'] / 1e6:.3f} MW")

    print("\ncontour-integral convergence (NASA air): n per process, |p-v area - w_net|, |T-s area - w_net| (J/kg)")
    wn = brayton_performance(AIR, T1, p1, pi, T3)["w_net"]
    e_pv, e_ts, ns = [], [], [25, 50, 100, 200, 400]
    for n in ns:
        path = cycle_path(AIR, T1, p1, pi, T3, n)
        e_pv.append(abs(enclosed_area(path["v"], path["p"]) - wn))
        e_ts.append(abs(enclosed_area(path["s"], path["T"]) - wn))
        print(f"  {n:4d}   {e_pv[-1]:.3e}   {e_ts[-1]:.3e}")
    print(f"  observed orders p-v {np.array2string(observed_order(e_pv), precision=2)}, T-s {np.array2string(observed_order(e_ts), precision=2)}")

    print("\nsweep, NASA air: T3 (K), pi, w_net (kJ/kg), eta, eta closed form (gamma 1.4), Carnot")
    for T3s in (1100.0, 1400.0, 1700.0):
        for p in (4.0, 8.0, 12.0, 20.0, 30.0, 40.0):
            pf = brayton_performance(AIR, T1, p1, p, T3s)
            print(f"  {T3s:6.0f}  {p:5.1f}   {pf['w_net'] / 1e3:7.2f}   {pf['eta']:.4f}   {float(efficiency_closed_form(p)):.4f}   {pf['eta_carnot']:.4f}")
        print(f"         optimum: closed form {optimum_pressure_ratio_closed_form(T1, T3s):.2f}, NASA air {optimum_pressure_ratio(AIR, T1, p1, T3s):.2f}")

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    path = cycle_path(AIR, T1, p1, pi, T3, 200)
    st = brayton_states(AIR, T1, p1, pi, T3)
    ax[0].plot(path["v"], path["p"] / 1e3)
    ax[0].plot(st["v"], st["p"] / 1e3, "ko")
    for k in range(4):
        ax[0].annotate(str(k + 1), (st["v"][k], st["p"][k] / 1e3), textcoords="offset points", xytext=(5, 5))
    ax[0].set_xlabel("v (m^3/kg)")
    ax[0].set_ylabel("p (kPa)")
    ax[0].set_title("ideal Brayton, p-v (NASA air)")
    ax[1].plot(path["s"], path["T"])
    ax[1].plot(st["s"], st["T"], "ko")
    for k in range(4):
        ax[1].annotate(str(k + 1), (st["s"][k], st["T"][k]), textcoords="offset points", xytext=(5, 5))
    ax[1].set_xlabel("s (J/(kg K))")
    ax[1].set_ylabel("T (K)")
    ax[1].set_title("ideal Brayton, T-s (NASA air)")
    pp = np.linspace(1.5, 60, 80)
    for T3s in (1100.0, 1400.0, 1700.0):
        ax[2].plot(pp, [brayton_performance(AIR, T1, p1, p, T3s)["w_net"] / 1e3 for p in pp], label=f"w_net, T3 = {T3s:.0f} K")
    ax2 = ax[2].twinx()
    ax2.plot(pp, efficiency_closed_form(pp), "k--", label="eta, gamma 1.4")
    ax2.set_ylabel("thermal efficiency")
    ax[2].set_xlabel("pressure ratio")
    ax[2].set_ylabel("net work (kJ/kg)")
    ax[2].set_title("work peaks, efficiency keeps rising")
    ax[2].legend(loc="lower right", fontsize=7)
    ax2.legend(loc="upper left", fontsize=7)
    fig.tight_layout()
    fig.savefig(out_dir / "brayton_cycle.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'brayton_cycle.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
