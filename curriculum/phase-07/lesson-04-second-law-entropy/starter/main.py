"""P07.04 The Second Law & Entropy: starter code.

Complete every function marked NotImplementedError. When finished, running this file prints the
lesson's tables and writes its figure.

An entropy-generation tracker. Each irreversible process the lesson treats returns its entropy
generation S_gen (J/K, or W/K for a steady flow), computed from an entropy balance

    S_gen = (entropy out) - (entropy in) - sum(Q_k / T_k)  >= 0,

and an EntropyLedger collects them, refuses any negative entry (a second-law violation means a bug or an
impossible claim), ranks them and converts them to lost work T0 S_gen (the Gouy-Stodola preview of P07.08):

* heat crossing a finite temperature difference, and two bodies coming to thermal equilibrium;
* throttling and adiabatic duct friction of an ideal gas (h constant, so T constant);
* a two-stream heat exchanger (a bleed-air precooler);
* adiabatic mixing of streams of one gas at different temperatures, and mixing of different gases;
* friction in an incompressible liquid line (a fuel or propellant feed line);
* a cycle's Clausius integral, sum(Q_k / T_k) <= 0, and the Carnot bound on any cycle the tracker accepts.

Entropy is a property: ds_ideal_gas gives s2 - s1 from the end states, and ds_path_integral integrates
delta q_rev / T = cp dT / T - R dp / p along any path between them; the two must agree.

Air is dry air as fixed for Phase 07 (N2 0.7808, O2 0.2095, Ar 0.0093, CO2 0.0004), with variable cp from
the NASA 9-coefficient polynomials (lower interval, 200-1000 K), as built in P07.01.

Units: SI (K, Pa, J, kg, s). Specific entropy in J/(kg K).

Run:  python main.py              (prints the ledger and the checks, writes entropy_generation.png)
Test: pytest curriculum/phase-07/lesson-04-second-law-entropy
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.optimize import brentq

R_U = 8.31446261815324            # J/(mol K), CODATA 2018
P_REF = 1.0e5                     # Pa, NASA standard-state pressure
T0_DEAD = 288.15                  # K, dead-state temperature used for lost work T0 S_gen

# NASA Glenn 9-coefficient polynomials, lower interval (200-1000 K), from data/thermo.inp of NASA CEA
# (https://github.com/nasa/cea, Apache-2.0; McBride, Zehe & Gordon, NASA/TP-2002-211556): [a1..a7], [b1, b2], M in g/mol.
NASA9_LOW = {
    "N2": (28.0134, [22103.71497, -381.846182, 6.08273836, -0.00853091441, 1.384646189e-05, -9.62579362e-09, 2.519705809e-12], [710.846086, -10.76003744]),
    "O2": (31.9988, [-34255.6342, 484.700097, 1.119010961, 0.00429388924, -6.83630052e-07, -2.0233727e-09, 1.039040018e-12], [-3391.45487, 18.4969947]),
    "Ar": (39.948, [0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0], [-745.375, 4.37967491]),
    "CO2": (44.0095, [49436.5054, -626.411601, 5.30172524, 0.002503813816, -2.127308728e-07, -7.68998878e-10, 2.849677801e-13], [-45281.9846, -7.04827944]),
}
AIR = {"N2": 0.7808, "O2": 0.2095, "Ar": 0.0093, "CO2": 0.0004}
M_AIR = sum(x * NASA9_LOW[sp][0] for sp, x in AIR.items()) * 1e-3     # kg/mol
R_AIR = R_U / M_AIR                                                  # J/(kg K)


# ---------------------------------------------------------------- given helpers: dry air (from P07.01)

def _check_T(T):
    T = np.asarray(T, dtype=float)
    if np.any(T < 200.0) or np.any(T > 1000.0):
        raise ValueError("this lesson's air model covers 200-1000 K (the NASA lower interval)")
    return T


def air_cp(T):
    """cp of dry air [J/(kg K)], NASA polynomials, 200-1000 K."""
    T = _check_T(T)
    tot = 0.0
    for sp, x in AIR.items():
        M, a, b = NASA9_LOW[sp]
        tot = tot + x * (a[0] / T**2 + a[1] / T + a[2] + a[3] * T + a[4] * T**2 + a[5] * T**3 + a[6] * T**4)
    return R_U * tot / M_AIR


def air_h(T):
    """h of dry air [J/kg] (NASA formation-enthalpy datum), 200-1000 K."""
    T = _check_T(T)
    tot = 0.0
    for sp, x in AIR.items():
        M, a, b = NASA9_LOW[sp]
        tot = tot + x * (-a[0] / T**2 + a[1] * np.log(T) / T + a[2] + a[3] * T / 2 + a[4] * T**2 / 3
                         + a[5] * T**3 / 4 + a[6] * T**4 / 5 + b[0] / T) * T
    return R_U * tot / M_AIR


def air_s0(T):
    """Temperature part of the entropy of dry air at p deg = 1 bar [J/(kg K)], 200-1000 K (the constant
    mixing term cancels in every difference this lesson takes)."""
    T = _check_T(T)
    tot = 0.0
    for sp, x in AIR.items():
        M, a, b = NASA9_LOW[sp]
        tot = tot + x * (-a[0] / (2 * T**2) - a[1] / T + a[2] * np.log(T) + a[3] * T + a[4] * T**2 / 2
                         + a[5] * T**3 / 3 + a[6] * T**4 / 4 + b[1])
    return R_U * tot / M_AIR


def air_T_from_h(h: float) -> float:
    """Inverse of air_h on 200-1000 K."""
    return float(brentq(lambda T: air_h(T) - h, 200.0, 1000.0, xtol=1e-11))


# ---------------------------------------------------------------- entropy is a property

def ds_ideal_gas(T1: float, p1: float, T2: float, p2: float) -> float:
    """s2 - s1 of dry air [J/(kg K)] between two states: s deg(T2) - s deg(T1) - R ln(p2 / p1)."""
    raise NotImplementedError


def ds_path_integral(T_path, p_path, n: int = 2001) -> float:
    """Integral of delta q_rev / T = cp(T) dT / T - R dp / p along the path lambda -> (T(lambda), p(lambda)),
    lambda in [0, 1], by the trapezoid rule on n points. T_path and p_path are callables of lambda."""
    raise NotImplementedError


# ---------------------------------------------------------------- the irreversible processes

def heat_transfer_sgen(Q: float, T_hot: float, T_cold: float) -> float:
    """Entropy generated when heat Q [J] (or Qdot [W]) leaves a reservoir at T_hot and enters one at T_cold:
    S_gen = Q / T_cold - Q / T_hot (zero only when T_hot = T_cold)."""
    raise NotImplementedError


def bodies_equilibrate(m1: float, c1: float, T1: float, m2: float, c2: float, T2: float) -> tuple[float, float]:
    """Two incompressible bodies (mass m, specific heat c) in contact inside an insulated box reach a common
    temperature Tf = (m1 c1 T1 + m2 c2 T2) / (m1 c1 + m2 c2). Returns (Tf, S_gen), with
    S_gen = m1 c1 ln(Tf / T1) + m2 c2 ln(Tf / T2)."""
    raise NotImplementedError


def throttle_sgen(p1: float, p2: float) -> float:
    """Specific entropy generated by adiabatically throttling an ideal gas from p1 to p2 [J/(kg K)]:
    h2 = h1 means T2 = T1, so s_gen = s2 - s1 = R ln(p1 / p2). The same formula covers adiabatic duct
    friction of a gas at low Mach number."""
    raise NotImplementedError


def heat_exchanger(mdot_h: float, T_h_in: float, T_h_out: float, mdot_c: float, T_c_in: float) -> dict[str, float]:
    """Counterflow air-to-air heat exchanger, adiabatic outside, no pressure loss. The hot stream is cooled
    from T_h_in to T_h_out; the energy balance mdot_h (h_h,in - h_h,out) = mdot_c (h_c,out - h_c,in) gives
    T_c_out. Raises ValueError if the temperatures cross (T_c_out >= T_h_in or T_h_out <= T_c_in), which no
    counterflow exchanger can achieve. Returns Qdot [W], T_c_out [K] and Sdot_gen [W/K] =
    mdot_h (s_h,out - s_h,in) + mdot_c (s_c,out - s_c,in) (each stream isobaric)."""
    raise NotImplementedError


def mix_streams(mdots, Ts, ps, p_out: float) -> tuple[float, float]:
    """Adiabatic steady mixing of air streams (mass flows mdots, temperatures Ts, pressures ps) into one
    stream at p_out. Energy: sum mdot_i h_i = (sum mdot_i) h_out gives T_out. Entropy:
    Sdot_gen = sum mdot_i (s(T_out, p_out) - s(T_i, p_i)). Returns (T_out, Sdot_gen)."""
    raise NotImplementedError


def mix_gases_sgen(moles: dict[str, float]) -> float:
    """Entropy generated when different ideal gases, initially separate at the same T and p, mix at that T
    and p [J/K]: S_gen = -R_U sum n_i ln x_i (each gas expands into the whole volume)."""
    raise NotImplementedError


def liquid_line_friction(mdot: float, dp: float, rho: float, c: float, T_in: float) -> tuple[float, float]:
    """Adiabatic flow of an incompressible liquid through a line or filter with pressure drop dp [Pa].
    h = c T + p / rho is constant, so the friction heats the liquid: T_out = T_in + dp / (rho c).
    Returns (T_out, Sdot_gen) with Sdot_gen = mdot c ln(T_out / T_in) [W/K] (about mdot dp / (rho T))."""
    raise NotImplementedError


# ---------------------------------------------------------------- cycles

def clausius_integral(Q, T) -> float:
    """sum(Q_k / T_k) over a cycle's heat exchanges (Q into the cycle positive, T the boundary
    temperature where each crosses) [J/K]. Clausius: <= 0, with equality only for a reversible cycle."""
    raise NotImplementedError


def carnot_efficiency(T_hot: float, T_cold: float) -> float:
    """eta_Carnot = 1 - T_cold / T_hot."""
    raise NotImplementedError


def check_cycle(Q, T, tol: float = 1e-12) -> dict[str, float]:
    """Accept or reject a heat-engine cycle given its heat exchanges with reservoirs (Q_k > 0 into the
    cycle at T_k). The first law gives W = sum Q_k; the second law gives S_gen = -sum Q_k / T_k, which must
    be >= 0 (ValueError otherwise). Returns W, Q_in = sum of positive Q, eta = W / Q_in, S_gen and the
    Carnot efficiency between the hottest and coldest reservoirs used."""
    raise NotImplementedError


# ---------------------------------------------------------------- the tracker

class EntropyLedger:
    """Collects entropy-generation entries (name, S_gen) for a system, refusing negative ones."""

    def __init__(self, T0: float = T0_DEAD, tol: float = 1e-9):
        self.T0 = T0
        self.tol = tol
        self.entries: list[tuple[str, float]] = []

    def add(self, name: str, S_gen: float) -> float:
        """Record an entry; raise ValueError if S_gen < -tol (the second law forbids it). Returns S_gen."""
        raise NotImplementedError

    def total(self) -> float:
        """Sum of all entries."""
        raise NotImplementedError

    def lost_work(self) -> float:
        """T0 times the total entropy generation: the work the irreversibilities cost (P07.08)."""
        raise NotImplementedError

    def ranked(self) -> list[tuple[str, float, float]]:
        """Entries sorted from largest to smallest as (name, S_gen, fraction of the total)."""
        raise NotImplementedError


# Typical, illustrative bleed-air chain of an airliner environmental control system (no cited source;
# numbers chosen to be realistic in magnitude): 1 kg/s of bleed air at 1.0 MPa and 700 K.
BLEED = {"mdot": 1.0, "p_bleed": 1.0e6, "T_bleed": 700.0, "p_prv": 0.30e6, "T_precooled": 470.0,
         "mdot_fan": 1.2, "T_fan": 320.0, "p_duct_out": 0.28e6,
         "mdot_pack": 0.5, "T_pack": 270.0, "mdot_recirc": 0.5, "T_recirc": 297.0, "p_mix": 0.105e6}


def bleed_air_ledger(case: dict = BLEED, T0: float = T0_DEAD) -> EntropyLedger:
    """The entropy ledger of the bleed-air chain: pressure-regulating valve (throttle), precooler (finite
    temperature difference against fan air), duct friction, and the mix manifold (pack outlet air mixed
    with recirculated cabin air at the mix pressure)."""
    L = EntropyLedger(T0)
    c = case
    L.add("pressure-regulating valve (throttling)", c["mdot"] * throttle_sgen(c["p_bleed"], c["p_prv"]))
    L.add("precooler (heat across a finite dT)", heat_exchanger(c["mdot"], c["T_bleed"], c["T_precooled"], c["mdot_fan"], c["T_fan"])["S_gen"])
    L.add("duct friction", c["mdot"] * throttle_sgen(c["p_prv"], c["p_duct_out"]))
    L.add("mix manifold (mixing)", mix_streams([c["mdot_pack"], c["mdot_recirc"]], [c["T_pack"], c["T_recirc"]],
                                               [c["p_mix"], c["p_mix"]], c["p_mix"])[1])
    return L


# ---------------------------------------------------------------- driver

def main(out_dir: str | Path = ".") -> None:
    """Print the ledger, the checks and the cycle tests; write entropy_generation.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    print(f"dry air R = {R_AIR:.3f} J/(kg K); dead state T0 = {T0_DEAD} K")

    print("\nentropy is a property: air from (300 K, 100 kPa) to (800 K, 1 MPa), J/(kg K)")
    T1, p1, T2, p2 = 300.0, 1e5, 800.0, 1e6
    print(f"  end states:                      {ds_ideal_gas(T1, p1, T2, p2):.6f}")
    paths = {
        "isobaric then isothermal": (lambda l: np.where(l < 0.5, T1 + (T2 - T1) * 2 * l, T2),
                                     lambda l: np.where(l < 0.5, p1, p1 + (p2 - p1) * (2 * l - 1))),
        "isothermal then isobaric": (lambda l: np.where(l < 0.5, T1, T1 + (T2 - T1) * (2 * l - 1)),
                                     lambda l: np.where(l < 0.5, p1 + (p2 - p1) * 2 * l, p2)),
        "wiggly path": (lambda l: T1 + (T2 - T1) * l + 150 * np.sin(3 * np.pi * l),
                        lambda l: p1 * (p2 / p1) ** l * (1 + 0.5 * np.sin(2 * np.pi * l) ** 2)),
    }
    for name, (Tp, pp) in paths.items():
        print(f"  {name:<31}  {ds_path_integral(Tp, pp, 20001):.6f}")
    Tp, pp = paths["wiggly path"]
    exact = ds_ideal_gas(T1, p1, T2, p2)
    errs = [abs(ds_path_integral(Tp, pp, n) - exact) for n in (101, 201, 401, 801)]
    print("  wiggly-path trapezoid error, n = 101..801: " + ", ".join(f"{e:.2e}" for e in errs)
          + f"; observed orders {np.round(np.log2(np.array(errs[:-1]) / np.array(errs[1:])), 2)}")

    print("\nS_gen -> 0 as the temperature difference -> 0 (Q = 1 kJ into a reservoir at 300 K):")
    for dT in (100.0, 10.0, 1.0, 0.1):
        print(f"  dT = {dT:6.1f} K: heat transfer S_gen = {heat_transfer_sgen(1000.0, 300.0 + dT, 300.0):.3e} J/K; "
              f"two 1 kg aluminium blocks (c = 900) at 300 +/- dT/2: {bodies_equilibrate(1, 900, 300 + dT / 2, 1, 900, 300 - dT / 2)[1]:.3e} J/K")

    L = bleed_air_ledger()
    c = BLEED
    hx = heat_exchanger(c["mdot"], c["T_bleed"], c["T_precooled"], c["mdot_fan"], c["T_fan"])
    Tmix = mix_streams([c["mdot_pack"], c["mdot_recirc"]], [c["T_pack"], c["T_recirc"]], [c["p_mix"]] * 2, c["p_mix"])[0]
    print(f"\nbleed-air chain (typical, illustrative): precooler Q = {hx['Q']/1e3:.1f} kW, fan air leaves at {hx['T_c_out']:.2f} K; "
          f"mix manifold T_out = {Tmix:.2f} K")
    print(f"{'process':<42} {'S_gen (W/K)':>12} {'share':>7} {'T0 S_gen (kW)':>14}")
    for name, s, f in L.ranked():
        print(f"{name:<42} {s:>12.2f} {100*f:>6.1f}% {L.T0*s/1e3:>14.2f}")
    print(f"{'total':<42} {L.total():>12.2f} {'':>7} {L.lost_work()/1e3:>14.2f}")
    print(f"  for scale: the bleed air's enthalpy above T0 is {c['mdot']*(air_h(c['T_bleed'])-air_h(T0_DEAD))/1e3:.1f} kW")

    Tout, Sf = liquid_line_friction(2.0, 1.5e5, 800.0, 2000.0, 300.0)
    print(f"\nfuel line (typical, illustrative): 2 kg/s, dp = 150 kPa, rho = 800, c = 2000: T rises {Tout-300:.4f} K, "
          f"Sdot_gen = {Sf:.4f} W/K (mdot dp / (rho T) = {2.0*1.5e5/(800*300):.4f}); lost work {T0_DEAD*Sf:.1f} W "
          f"against pumping power mdot dp / rho = {2.0*1.5e5/800:.1f} W")
    print(f"mixing 0.7808 mol N2 with 0.2095 mol O2 and 0.0093 mol Ar at the same T and p: S_gen = "
          f"{mix_gases_sgen({'N2': 0.7808, 'O2': 0.2095, 'Ar': 0.0093}):.4f} J/K")

    print("\ncycle claims (typical, illustrative): heat in at 1050 K, rejected at 350 K, Carnot "
          f"{carnot_efficiency(1050.0, 350.0):.4f}")
    for W in (60.0, 100.0, 166.0, 180.0):
        Qin = 250.0
        try:
            r = check_cycle([Qin, -(Qin - W)], [1050.0, 350.0])
            print(f"  W = {W:5.1f} W from {Qin:.0f} W: accepted, eta = {r['eta']:.3f}, S_gen = {r['S_gen']:.4f} W/K")
        except ValueError as e:
            print(f"  W = {W:5.1f} W from {Qin:.0f} W: rejected ({e})")

    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
    names = [n.split(" (")[0] for n, _, _ in L.ranked()]
    axs[0].barh(names[::-1], [s for _, s, _ in L.ranked()][::-1])
    axs[0].set_xlabel("Sdot_gen (W/K)")
    axs[0].set_title("bleed-air chain (illustrative)")
    dTs = np.logspace(-2, 2, 30)
    axs[1].loglog(dTs, heat_transfer_sgen(1000.0, 300.0 + dTs, 300.0), label="Q = 1 kJ across dT")
    axs[1].loglog(dTs, [bodies_equilibrate(1, 900, 300 + d / 2, 1, 900, 300 - d / 2)[1] for d in dTs], label="two blocks, initial dT")
    axs[1].set_xlabel("dT (K)")
    axs[1].set_ylabel("S_gen (J/K)")
    axs[1].set_title("S_gen -> 0 as dT -> 0 (slopes 1 and 2)")
    axs[1].legend(fontsize=8)
    lam = np.linspace(0, 1, 400)
    for name, (Tp, pp) in paths.items():
        axs[2].plot(np.asarray(Tp(lam)) + 0 * lam, np.asarray(pp(lam)) / 1e3 + 0 * lam, label=name)
    axs[2].set_yscale("log")
    axs[2].set_xlabel("T (K)")
    axs[2].set_ylabel("p (kPa)")
    axs[2].set_title(f"three paths, one ds = {exact:.2f} J/(kg K)")
    axs[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "entropy_generation.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
