"""P07.05 Isentropic Processes & Efficiencies: starter code.

How hot does the air leave an aero-engine compressor, how much work does it take, and what do
"85 % efficient" and "90 % polytropic" actually mean? This module is a compressor (and turbine)
efficiency calculator for an ideal gas:

* the isentropic exit temperature, from s°(T2) - s°(T1) = R ln(p2/p1) with the NASA 9-coefficient
  dry-air model (variable cp), which reduces to T2 = T1 pi^((gamma-1)/gamma) for constant cp;
* the actual exit state for a given isentropic efficiency, through enthalpy, not temperature;
* the polytropic (small-stage) efficiency, s°(T2) - s°(T1) = (R/eta_p) ln pi, and its closed form
  for constant cp;
* a multi-stage compressor built from equal-pressure-ratio stages, whose overall efficiency falls
  below the stage efficiency and tends to the polytropic closed form as the stages get small;
* the entropy generated, the work, and two independent checks of the first law (the work as the
  integral of cp dT, and as the integral of v dp + T ds along the polytropic path).

Property model: NASA Glenn 9-coefficient polynomials for N2, O2, Ar and CO2 (NASA CEA data/thermo.inp,
McBride, Zehe & Gordon, NASA/TP-2002-211556), combined into the Phase 07 dry-air model with mole
fractions N2 0.7808, O2 0.2095, Ar 0.0093, CO2 0.0004. s° is at p° = 1 bar. R = 8.31446261815324 J/(mol K)
(CODATA 2018). The air is an ideal gas and its composition is frozen.

Units: SI (K, Pa, J/kg, J/(kg K)).

Run:  python main.py [out_dir]      (prints the lesson's numbers, writes compressor_efficiency.png)
Test: pytest curriculum/phase-07/lesson-05-isentropic-processes-efficiencies
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import quad

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
    raise NotImplementedError


def temperature_from_h(gas: Gas, h_target: float | np.ndarray) -> float | np.ndarray:
    """The temperature at which gas.h(T) = h_target (dh/dT = cp > 0), by Newton's method. Accepts arrays."""
    raise NotImplementedError


# ---------------------------------------------------------------- isentropic and polytropic relations

def isentropic_temperature_closed_form(T1: float, pi: float, gamma: float = 1.4) -> float:
    """Constant-cp isentropic exit temperature, T2s = T1 pi^((gamma-1)/gamma) (pi = p2/p1; pi < 1 expands)."""
    raise NotImplementedError


def isentropic_temperature(gas: Gas, T1: float, pi: float) -> float:
    """Isentropic exit temperature for any ideal gas: s0(T2s) = s0(T1) + R ln(pi)."""
    raise NotImplementedError


def compressor_exit_temperature(gas: Gas, T1: float, pi: float, eta_s: float) -> float:
    """Actual compressor exit temperature for isentropic efficiency eta_s = (h2s - h1)/(h2 - h1):
    h2 = h1 + (h2s - h1)/eta_s, then invert h."""
    raise NotImplementedError


def turbine_exit_temperature(gas: Gas, T1: float, pi_t: float, eta_s: float) -> float:
    """Actual turbine exit temperature for expansion ratio pi_t = p1/p2 > 1 and isentropic efficiency
    eta_s = (h1 - h2)/(h1 - h2s): h2 = h1 - eta_s (h1 - h2s)."""
    raise NotImplementedError


def polytropic_exit_temperature(gas: Gas, T1: float, pi: float, eta_p: float) -> float:
    """Exit temperature of a compression (pi > 1) or expansion (pi < 1) at constant polytropic efficiency.
    Compressor, dh = v dp / eta_p, so ds0 = (R / eta_p) d ln p; turbine, dh = eta_p v dp, so
    ds0 = eta_p R d ln p. Integrated: s0(T2) = s0(T1) + R ln(pi) / eta_p (pi > 1) or + eta_p R ln(pi) (pi < 1)."""
    raise NotImplementedError


def isentropic_efficiency(gas: Gas, T1: float, T2: float, pi: float) -> float:
    """Isentropic efficiency of a compressor (pi > 1: (h2s - h1)/(h2 - h1)) or turbine (pi < 1:
    (h1 - h2)/(h1 - h2s)) from its measured inlet and exit temperatures and pressure ratio pi = p2/p1."""
    raise NotImplementedError


def polytropic_efficiency(gas: Gas, T1: float, T2: float, pi: float) -> float:
    """Polytropic efficiency from inlet and exit temperatures: compressor R ln(pi)/(s0(T2) - s0(T1)),
    turbine (s0(T2) - s0(T1))/(R ln(pi))."""
    raise NotImplementedError


def isentropic_from_polytropic(pi: float, eta_p: float, gamma: float = 1.4) -> float:
    """Constant-cp closed form for a compressor of pressure ratio pi > 1:
    eta_s = (pi^m - 1)/(pi^(m/eta_p) - 1), m = (gamma-1)/gamma."""
    raise NotImplementedError


def isentropic_from_polytropic_turbine(pi_t: float, eta_p: float, gamma: float = 1.4) -> float:
    """Constant-cp closed form for a turbine of expansion ratio pi_t = p1/p2 > 1:
    eta_s = (1 - pi_t^(-eta_p m)) / (1 - pi_t^(-m)), m = (gamma-1)/gamma."""
    raise NotImplementedError


# ---------------------------------------------------------------- work, entropy, stages

def specific_work(gas: Gas, T1: float, T2: float) -> float:
    """Steady-flow work input per unit mass of an adiabatic compressor (first law, kinetic and potential
    energy changes neglected): w = h2 - h1, J/kg. Negative for a turbine (work out)."""
    raise NotImplementedError


def entropy_generation(gas: Gas, T1: float, T2: float, pi: float) -> float:
    """Entropy generated per unit mass in an adiabatic steady-flow device (second law, no heat):
    s_gen = s2 - s1 = s0(T2) - s0(T1) - R ln(pi), J/(kg K)."""
    raise NotImplementedError


def multistage_compressor(gas: Gas, T1: float, pi_total: float, n_stages: int, eta_stage: float) -> np.ndarray:
    """Exit temperatures of n equal-pressure-ratio stages, each of isentropic efficiency eta_stage
    (stage pressure ratio pi_total^(1/n)). Returns an array of n + 1 temperatures starting with T1."""
    raise NotImplementedError


def first_law_integrals(gas: Gas, T1: float, pi: float, eta_p: float, n: int = 2001) -> tuple[float, float, float, float]:
    """Two independent routes to the work of a polytropic compression from T1 through pi.
    Returns (h2 - h1, integral of cp dT by quad, integral of v dp along the path, integral of T ds along
    the path). Along the path ln p runs from 0 to ln pi; at each p the temperature satisfies
    s0(T) = s0(T1) + (R/eta_p) ln(p/p1); v dp = R T d ln p and T ds = T (ds0 - R d ln p).
    Gibbs (dh = T ds + v dp) says the last two sum to the first."""
    raise NotImplementedError


def observed_order(errors: np.ndarray, ratio: float = 2.0) -> np.ndarray:
    """Observed convergence orders log(e_k / e_{k+1}) / log(ratio) for refinements by ``ratio``."""
    e = np.asarray(errors, dtype=float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


# ---------------------------------------------------------------- the lesson's numbers

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save compressor_efficiency.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    cold = constant_cp_gas(1.4)
    T1, pi, eta_p, T0 = 288.15, 40.0, 0.90, 288.15
    print(f"dry-air model: M = {M_AIR * 1e3:.4f} g/mol, R = {R_AIR:.3f} J/(kg K), cp(300 K) = {float(cp_air(300.0)):.1f}, "
          f"cp(900 K) = {float(cp_air(900.0)):.1f} J/(kg K), gamma(300 K) = {float(cp_air(300.0) / (cp_air(300.0) - R_AIR)):.4f}, "
          f"gamma(900 K) = {float(cp_air(900.0) / (cp_air(900.0) - R_AIR)):.4f}")

    print(f"\naero-engine core compressor (typical, illustrative): T1 = {T1} K, overall pi = {pi:.0f}, eta_p = {eta_p}")
    for gas in (cold, AIR):
        T2s = isentropic_temperature(gas, T1, pi)
        T2 = polytropic_exit_temperature(gas, T1, pi, eta_p)
        eta_s = isentropic_efficiency(gas, T1, T2, pi)
        w, ws = specific_work(gas, T1, T2), specific_work(gas, T1, T2s)
        sg = entropy_generation(gas, T1, T2, pi)
        print(f"  {gas.name:28s}: T2s = {T2s:7.2f} K, T2 = {T2:7.2f} K, eta_s = {eta_s:.4f}, w = {w / 1e3:7.2f} kJ/kg "
              f"(isentropic {ws / 1e3:7.2f}), s_gen = {sg:6.2f} J/(kg K), T0 s_gen = {T0 * sg / 1e3:6.2f} kJ/kg")
    print(f"  closed form, gamma 1.4: T1 pi^0.2857 = {isentropic_temperature_closed_form(T1, pi):.2f} K, "
          f"eta_s(eta_p) = {isentropic_from_polytropic(pi, eta_p):.4f}")

    dh, cpdT, vdp, Tds = first_law_integrals(AIR, T1, pi, eta_p)
    print(f"  first law, NASA air: h2 - h1 = {dh / 1e3:.4f} kJ/kg; int cp dT = {cpdT / 1e3:.4f}; int v dp = {vdp / 1e3:.4f}; "
          f"int T ds = {Tds / 1e3:.4f}; v dp + T ds = {(vdp + Tds) / 1e3:.4f} kJ/kg")
    print(f"  T0 s_gen = {T0 * entropy_generation(AIR, T1, polytropic_exit_temperature(AIR, T1, pi, eta_p), pi) / 1e3:.2f} kJ/kg "
          f"is less than int T ds = {Tds / 1e3:.2f} kJ/kg: friction heat generated hot is partly recoverable")

    print("\nmulti-stage compressor, NASA air, pi = 40, stage isentropic efficiency 0.90:")
    print("   n    T_exit (K)   overall eta_s   overall eta_p   |T_exit - polytropic limit| (K)")
    T_lim = polytropic_exit_temperature(AIR, T1, pi, 0.90)
    errs = []
    ns = [1, 2, 4, 8, 16, 32, 64]
    for n in ns:
        T = multistage_compressor(AIR, T1, pi, n, 0.90)
        errs.append(abs(T[-1] - T_lim))
        print(f"  {n:3d}    {T[-1]:8.2f}      {isentropic_efficiency(AIR, T1, T[-1], pi):.4f}          "
              f"{polytropic_efficiency(AIR, T1, T[-1], pi):.4f}          {errs[-1]:.3f}")
    print(f"  polytropic limit {T_lim:.2f} K; observed orders {np.array2string(observed_order(errs[2:]), precision=2)}")

    print("\nsingle axial stage (typical, illustrative): pi = 1.35, eta_s = 0.90 at 288.15 K")
    Ts = compressor_exit_temperature(AIR, T1, 1.35, 0.90)
    print(f"  T2 = {Ts:.2f} K, eta_p = {polytropic_efficiency(AIR, T1, Ts, 1.35):.4f}, w = {specific_work(AIR, T1, Ts) / 1e3:.2f} kJ/kg")

    print("\nturbine (typical, illustrative): T1 = 1400 K, expansion ratio 12, eta_p = 0.90")
    Tt = polytropic_exit_temperature(AIR, 1400.0, 1 / 12, 0.90)
    print(f"  T2 = {Tt:.2f} K, eta_s = {isentropic_efficiency(AIR, 1400.0, Tt, 1 / 12):.4f} (closed form gamma 1.4: "
          f"{isentropic_from_polytropic_turbine(12.0, 0.90):.4f}), w = {specific_work(AIR, 1400.0, Tt) / 1e3:.2f} kJ/kg")

    print("\nsweep at eta_p = 0.90: pi, eta_s (NASA air), eta_s (closed form, gamma 1.4), T2 NASA - T2 constant cp (K), work error of constant cp")
    pis = np.array([1.5, 2, 5, 10, 20, 30, 40, 50, 60])
    eta_nasa, eta_cf = [], []
    for p in pis:
        T2n = polytropic_exit_temperature(AIR, T1, p, eta_p)
        T2c = polytropic_exit_temperature(cold, T1, p, eta_p)
        eta_nasa.append(isentropic_efficiency(AIR, T1, T2n, p))
        eta_cf.append(isentropic_from_polytropic(p, eta_p))
        wn, wc = specific_work(AIR, T1, T2n), specific_work(cold, T1, T2c)
        print(f"  {p:5.1f}   {eta_nasa[-1]:.4f}   {eta_cf[-1]:.4f}   {T2n - T2c:+7.2f}   {(wc - wn) / wn * 100:+6.2f} %")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    pp = np.linspace(1.05, 60, 200)
    for ep in (0.85, 0.90, 0.95):
        ax[0].plot(pp, [isentropic_from_polytropic(p, ep) for p in pp], label=f"compressor, eta_p = {ep}")
        ax[0].plot(pp, [isentropic_from_polytropic_turbine(p, ep) for p in pp], "--", label=f"turbine, eta_p = {ep}")
    ax[0].plot(pis, eta_nasa, "ko", ms=4, label="compressor, NASA air, eta_p = 0.90")
    ax[0].set_xlabel("pressure ratio")
    ax[0].set_ylabel("isentropic efficiency")
    ax[0].set_title("isentropic vs polytropic efficiency")
    ax[0].legend(fontsize=7)
    ax[1].loglog(ns[2:], errs[2:], "o-", label="NASA air, stage eta_s = 0.90")
    ax[1].loglog(ns[2:], errs[2] * ns[2] / np.array(ns[2:]), "k--", lw=0.8, label="slope -1")
    ax[1].set_xlabel("number of stages")
    ax[1].set_ylabel("|T_exit - polytropic limit| (K)")
    ax[1].set_title("many small stages -> polytropic closed form")
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(out_dir / "compressor_efficiency.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'compressor_efficiency.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
