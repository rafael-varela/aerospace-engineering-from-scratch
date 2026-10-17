"""P07.07 Real Brayton Cycles: starter code.

A small gas-turbine core for a turboprop or an auxiliary power unit (APU) has a compressor and
turbine that are not isentropic, ducts, a combustor and heat exchangers that lose pressure, and
perhaps an intercooler between two compressors and a recuperator that preheats the compressed air
with the exhaust. Which of these features pays, by how much, and at what pressure ratio? This module
is an air-standard cycle model with

* two compressors (LP and HP) of isentropic efficiency eta_c, with an optional intercooler between
  them that cools the air to T_ic with a fractional pressure loss;
* an optional recuperator of effectiveness eps (enthalpy-based), with pressure losses on both sides;
* a combustor modelled as constant-pressure heat addition to the turbine inlet temperature, with a
  fractional pressure loss, and a turbine of isentropic efficiency eta_t;
* per-component first-law balances, the overall balance q_in - q_out = w_net, the Carnot bound,
  and the limits that recover the ideal cycle of P07.06 and the ideal recuperated-cycle closed form
  eta = 1 - (T1/T3) pi^((gamma-1)/gamma).

Station numbering: 1 LP compressor inlet, 2 LP exit, 3 intercooler exit (HP inlet), 4 HP exit,
5 recuperator cold-side exit (combustor inlet), 6 turbine inlet, 7 turbine exit (recuperator hot
inlet), 8 exhaust. Without an intercooler the LP compressor does the whole pressure ratio and
stations 3 and 4 equal 2; without a recuperator 5 equals 4 and 8 equals 7.

Property model: as in P07.05, the NASA Glenn 9-coefficient dry-air model (NASA CEA data/thermo.inp;
McBride, Zehe & Gordon, NASA/TP-2002-211556; N2 0.7808, O2 0.2095, Ar 0.0093, CO2 0.0004), or a
constant-cp gas for the closed-form checks. R = 8.31446261815324 J/(mol K) (CODATA 2018).

Units: SI (K, Pa, J/kg, W).

Run:  python main.py [out_dir]      (prints the lesson's numbers, writes real_brayton.png)
Test: pytest curriculum/phase-07/lesson-07-real-brayton-cycles
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, field
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


# ---------------------------------------------------------------- component models (given; built in P07.05)

def isentropic_temperature(gas: Gas, T1: float, pi: float) -> float:
    """Isentropic exit temperature for an ideal gas: s0(T2) = s0(T1) + R ln(pi)."""
    return temperature_from_s0(gas, float(gas.s0(T1)) + gas.R * np.log(pi))


def compressor_exit_temperature(gas: Gas, T1: float, pi: float, eta_s: float) -> float:
    """Compressor exit temperature for isentropic efficiency eta_s: h2 = h1 + (h2s - h1)/eta_s."""
    h1 = float(gas.h(T1))
    return temperature_from_h(gas, h1 + (float(gas.h(isentropic_temperature(gas, T1, pi))) - h1) / eta_s)


def turbine_exit_temperature(gas: Gas, T1: float, pi_t: float, eta_s: float) -> float:
    """Turbine exit temperature for expansion ratio pi_t = p_in/p_out > 1: h2 = h1 - eta_s (h1 - h2s)."""
    h1 = float(gas.h(T1))
    return temperature_from_h(gas, h1 - eta_s * (h1 - float(gas.h(isentropic_temperature(gas, T1, 1.0 / pi_t)))))


# ---------------------------------------------------------------- the cycle

@dataclass
class CycleSpec:
    """Inputs of the real cycle. Pressure losses are fractions of the inlet total pressure of each duct
    (dp = 0.03 means the exit pressure is 97 % of the inlet). pi is the overall compressor pressure
    ratio p4/p1; pi_lp is the LP compressor ratio p2/p1 (None: sqrt(pi), the optimum split)."""
    pi: float = 12.0
    T1: float = 288.15
    p1: float = 101325.0
    T_tit: float = 1400.0
    eta_c: float = 1.0
    eta_t: float = 1.0
    intercool: bool = False
    T_ic: float | None = None        # intercooler exit temperature, K (None: T1)
    pi_lp: float | None = None
    eps: float = 0.0                 # recuperator effectiveness, 0 = no recuperator
    dp_ic: float = 0.0               # intercooler
    dp_rc: float = 0.0               # recuperator cold side
    dp_cc: float = 0.0               # combustor
    dp_rh: float = 0.0               # recuperator hot side
    dp_ex: float = 0.0               # exhaust duct


def cycle_states(gas: Gas, spec: CycleSpec) -> dict[str, np.ndarray]:
    """Temperatures 'T', pressures 'p' and enthalpies 'h' at stations 1-8 (index 0-7).
    1 -> 2: LP compressor, ratio pi_lp (sqrt(pi) if None) when intercooled, else the whole pi.
    2 -> 3: intercooler to T_ic (T1 if None), p3 = p2 (1 - dp_ic); without one, 3 = 2.
    3 -> 4: HP compressor to p4 = pi p1 (ratio p4/p3); without an intercooler, 4 = 3.
    4 -> 5: recuperator cold side, h5 = h4 + eps (h(T7) - h4), p5 = p4 (1 - dp_rc).
    5 -> 6: combustor (heat addition) to T_tit, p6 = p5 (1 - dp_cc).
    6 -> 7: turbine, expanding to p7 = p1 / ((1 - dp_rh)(1 - dp_ex)) so that the exhaust leaves at p1.
    7 -> 8: recuperator hot side, h8 = h7 - (h5 - h4), p8 = p7 (1 - dp_rh).
    The turbine does not depend on the recuperator, so 7 is computed before 5."""
    raise NotImplementedError


def cycle_performance(gas: Gas, spec: CycleSpec) -> dict[str, float]:
    """Per unit mass of air: compressor work w_c = (h2 - h1) + (h4 - h3); turbine work w_t = h6 - h7;
    net work w_net = w_t - w_c; heat added q_in = h6 - h5; heat rejected in the intercooler
    q_ic = h2 - h3 and to the atmosphere q_ex = h8 - h1 (closing the air-standard loop); recuperated
    heat q_rec = h5 - h4; thermal efficiency eta = w_net / q_in; Carnot efficiency between T1 and T_tit;
    and 'residual' = q_in - q_ic - q_ex - w_net, the overall first-law closure (zero to round-off)."""
    raise NotImplementedError


def component_balances(gas: Gas, spec: CycleSpec) -> dict[str, float]:
    """Steady-flow energy balance residual (energy in - energy out, J/kg) of each component, computed from
    the station enthalpies and the work and heat of cycle_performance: LP and HP compressors
    (w + h_in - h_out), intercooler (h_in - q_ic - h_out), recuperator (hot-side loss minus cold-side gain),
    combustor (h5 + q_in - h6), turbine (h6 - w_t - h7), and the whole cycle. All should be zero."""
    raise NotImplementedError


def entropy_generation(gas: Gas, spec: CycleSpec, T_hot: float | None = None, T_cold: float | None = None) -> dict[str, float]:
    """Entropy generated per unit mass, J/(kg K), in each component and in total, with heat drawn from a
    reservoir at T_hot (default T_tit) and rejected to one at T_cold (default T1). Adiabatic components:
    s_out - s_in. Combustor: s6 - s5 - q_in/T_hot. Intercooler and exhaust: s_out - s_in + q/T_cold.
    For the complete cycle the sum equals q_out/T_cold - q_in/T_hot (the Clausius integral, negated)."""
    raise NotImplementedError


def ideal_recuperated_efficiency(pi: float, T1: float, T3: float, gamma: float = 1.4) -> float:
    """Closed form for the ideal recuperated cycle (constant cp, eps = 1, isentropic components,
    no losses): eta = 1 - (T1/T3) pi^((gamma-1)/gamma)."""
    raise NotImplementedError


def ideal_icr_efficiency(pi: float, T1: float, T3: float, gamma: float = 1.4) -> float:
    """Closed form for the ideal intercooled recuperated cycle with one intercooler back to T1 at the
    split sqrt(pi) (constant cp, eps = 1): w_t = cp T3 (1 - x^-1), w_c = 2 cp T1 (x^(1/2) - 1) with
    x = pi^((gamma-1)/gamma); q_in = w_t because the recuperated air reaches the turbine exit
    temperature; eta = 1 - 2 T1 (x^(1/2) - 1) / (T3 (1 - 1/x))."""
    raise NotImplementedError


def simple_cycle_efficiency(pi: float, gamma: float = 1.4) -> float:
    """The P07.06 closed form for the ideal simple cycle, eta = 1 - pi^(-(gamma-1)/gamma)."""
    raise NotImplementedError


def best_split(gas: Gas, spec: CycleSpec) -> float:
    """LP pressure ratio that minimises the total compressor work of an intercooled spec, by
    scipy.optimize.minimize_scalar (bounded) over ln(pi_lp) in (0, ln pi)."""
    raise NotImplementedError


def sweep(gas: Gas, spec: CycleSpec, pis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Thermal efficiency and net work of the spec over an array of overall pressure ratios."""
    eta, w = [], []
    for pi in pis:
        pf = cycle_performance(gas, CycleSpec(**{**spec.__dict__, "pi": float(pi)}))
        eta.append(pf["eta"])
        w.append(pf["w_net"])
    return np.array(eta), np.array(w)


# ---------------------------------------------------------------- the lesson's numbers

REAL = dict(eta_c=0.85, eta_t=0.88, dp_cc=0.04)                    # typical, illustrative
IC = dict(intercool=True, dp_ic=0.02)                               # typical, illustrative
REC = dict(eps=0.85, dp_rc=0.03, dp_rh=0.04)                        # typical, illustrative
CONFIGS = {"simple": {}, "intercooled": IC, "recuperated": REC, "intercooled + recuperated": {**IC, **REC}}


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save real_brayton.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    mdot = 5.0                                                       # kg/s, turboprop/APU scale, illustrative
    base = CycleSpec(pi=12.0, T_tit=1400.0)
    print(f"core (typical, illustrative): T1 = {base.T1} K, p1 = {base.p1:.0f} Pa, pi = {base.pi:.0f}, TIT = {base.T_tit:.0f} K, mdot = {mdot} kg/s")
    print(f"components: eta_c = {REAL['eta_c']} (each compressor), eta_t = {REAL['eta_t']}, combustor dp {REAL['dp_cc']:.0%}; "
          f"intercooler to T1 with dp {IC['dp_ic']:.0%}; recuperator eps = {REC['eps']}, dp {REC['dp_rc']:.0%} cold / {REC['dp_rh']:.0%} hot")
    print("\nNASA air, pi = 12:")
    print("  configuration                 eta     w_net (kJ/kg)  power (MW)  q_in (kJ/kg)  q_rec   T4 (K)   T7 (K)   T8 (K)  residual (J/kg)")
    for name, extra in CONFIGS.items():
        sp = CycleSpec(**{**base.__dict__, **REAL, **extra})
        pf = cycle_performance(AIR, sp)
        st = cycle_states(AIR, sp)
        print(f"  {name:28s}  {pf['eta']:.4f}   {pf['w_net'] / 1e3:8.2f}      {mdot * pf['w_net'] / 1e6:6.3f}     {pf['q_in'] / 1e3:8.2f}   "
              f"{pf['q_rec'] / 1e3:6.2f}  {st['T'][3]:7.2f}  {st['T'][6]:7.2f}  {st['T'][7]:7.2f}   {pf['residual']:.1e}")
    sp = CycleSpec(**{**base.__dict__, **REAL, **IC, **REC})
    st = cycle_states(AIR, sp)
    print("\n  intercooled + recuperated, station table:")
    for k in range(8):
        print(f"    {k + 1}  T = {st['T'][k]:8.2f} K   p = {st['p'][k] / 1e3:8.2f} kPa   h = {st['h'][k] / 1e3:8.2f} kJ/kg")
    bal = component_balances(AIR, sp)
    print("  component energy residuals (J/kg): " + ", ".join(f"{k} {v:.1e}" for k, v in bal.items()))
    sg = entropy_generation(AIR, sp)
    print("  entropy generated (J/(kg K)): " + ", ".join(f"{k} {v:.2f}" for k, v in sg.items()))
    pf = cycle_performance(AIR, sp)
    print(f"  Clausius check: q_out/T1 - q_in/TIT = {(pf['q_ic'] + pf['q_ex']) / sp.T1 - pf['q_in'] / sp.T_tit:.4f} J/(kg K); "
          f"Carnot {pf['eta_carnot']:.4f}")
    print(f"  best LP split by minimisation: {best_split(AIR, sp):.4f} (sqrt(pi) = {np.sqrt(sp.pi):.4f})")

    cold = constant_cp_gas(1.4)
    print("\nideal limits, constant cp (gamma 1.4), pi = 12, TIT = 1400 K:")
    pf0 = cycle_performance(cold, base)
    print(f"  eps = 0, ideal components: eta = {pf0['eta']:.6f}; P07.06 closed form {simple_cycle_efficiency(12.0):.6f}")
    pf1 = cycle_performance(cold, CycleSpec(**{**base.__dict__, 'eps': 1.0}))
    print(f"  eps = 1, ideal components: eta = {pf1['eta']:.6f}; 1 - (T1/T3) pi^m = {ideal_recuperated_efficiency(12.0, 288.15, 1400.0):.6f}")
    pf2 = cycle_performance(cold, CycleSpec(**{**base.__dict__, 'eps': 1.0, 'intercool': True}))
    print(f"  eps = 1 + intercooler, ideal: eta = {pf2['eta']:.6f}; closed form {ideal_icr_efficiency(12.0, 288.15, 1400.0):.6f}; "
          f"Carnot {1 - 288.15 / 1400.0:.6f}")
    print(f"  ideal recuperated cycle beats the simple cycle while pi < (T3/T1)^(gamma/(2(gamma-1))) = {(1400 / 288.15) ** 1.75:.2f}")

    print("\nsweep, NASA air, real components: pi, eta and w_net (kJ/kg) for each configuration")
    pis = np.array([3, 4, 5, 6, 8, 10, 12, 16, 20, 25, 30], dtype=float)
    fine = np.linspace(2.0, 40.0, 39)                                # step 1: contains every pi of the table
    peaks = {}
    for name, extra in CONFIGS.items():
        e, w = sweep(AIR, CycleSpec(**{**base.__dict__, **REAL, **extra}), fine)
        peaks[name] = (fine, e, w)
    rows = [int(np.argmin(np.abs(fine - p))) for p in pis]
    print("   pi   " + "   ".join(f"{n[:18]:>18s}" for n in CONFIGS))
    for i, p in zip(rows, pis):
        print(f"  {p:4.0f}  " + "   ".join(f"{peaks[n][1][i]:.4f} {peaks[n][2][i] / 1e3:7.1f}  " for n in CONFIGS))
    print("  peak efficiency and peak work (pi on a grid of step 1):")
    for name in CONFIGS:
        _, e, w = peaks[name]
        print(f"    {name:28s} eta_max {e.max():.4f} at pi {fine[e.argmax()]:.1f}; w_max {w.max() / 1e3:.1f} kJ/kg at pi {fine[w.argmax()]:.1f}")

    print("\nsensitivity, intercooled + recuperated at pi = 12 (NASA air): change in eta")
    ref = cycle_performance(AIR, sp)["eta"]
    for label, change in (("eta_c 0.85 -> 0.80", {"eta_c": 0.80}), ("eta_t 0.88 -> 0.83", {"eta_t": 0.83}),
                          ("eps 0.85 -> 0.75", {"eps": 0.75}), ("all pressure losses doubled", {"dp_ic": 0.04, "dp_rc": 0.06, "dp_cc": 0.08, "dp_rh": 0.08}),
                          ("TIT 1400 -> 1300 K", {"T_tit": 1300.0}), ("intercooler exit T1 + 20 K", {"T_ic": 308.15})):
        e = cycle_performance(AIR, CycleSpec(**{**sp.__dict__, **change}))["eta"]
        print(f"  {label:30s} eta {e:.4f} ({(e - ref) * 100:+.2f} points)")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for name in CONFIGS:
        f_, e, w = peaks[name]
        ax[0].plot(f_, e, label=name)
        ax[1].plot(f_, w / 1e3, label=name)
    ax[0].axhline(1 - 288.15 / 1400.0, color="k", ls=":", lw=0.8, label="Carnot, 288 K / 1400 K")
    ax[0].set_xlabel("overall pressure ratio")
    ax[0].set_ylabel("thermal efficiency")
    ax[0].set_title("real cycles, TIT 1400 K (NASA air)")
    ax[0].legend(fontsize=7)
    ax[1].set_xlabel("overall pressure ratio")
    ax[1].set_ylabel("net work (kJ/kg)")
    ax[1].set_title("specific work")
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out_dir / "real_brayton.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'real_brayton.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
