"""P08.04 Bernoulli & Energy Equation: reference solution.

The light aeroplane cruises at 3000 m at a true airspeed of 60 m/s. What do its pitot and static
ports sense, what does the airspeed indicator read, how big is the error from a static port that
sits in slightly fast air, and when does Bernoulli's equation stop being good enough? This module
applies conservation of momentum along a streamline (steady, inviscid, incompressible: Bernoulli)
and conservation of energy for a pipe or duct (with losses and shaft work):

    p + rho V^2 / 2 + rho g z = const along a streamline (no friction, no shaft work, rho const)
    p1/(rho g) + a1 V1^2/(2g) + z1 + h_pump = p2/(rho g) + a2 V2^2/(2g) + z2 + h_turbine + h_L

and builds the airspeed chain used in every cockpit:

    TAS (true airspeed) -> EAS = TAS sqrt(rho/rho0) (same dynamic pressure at sea-level density)
                        -> CAS (what a perfect instrument calibrated to the standard law reads)
                        -> IAS (what the dial shows, with position and instrument errors).

The standard calibration law is the compressible (isentropic, subsonic) pitot relation; it is
used here as given and derived in `P08.17`. The independent checks: every conversion is inverted
and the round trip closes to round-off; the compressible law tends to Bernoulli as M -> 0; the
Euler equation dp = -rho V dV integrated numerically through a wind-tunnel contraction against
Bernoulli; the venturi flow rate inverted; and an energy budget that balances.

Units: SI (m, s, kg, K, Pa, W). Knots appear only in printed output because pilots use them.

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes pitot_static.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-08/lesson-04-bernoulli-energy-equation
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable

import numpy as np

G0 = 9.80665
R_AIR = 287.05287
GAMMA = 1.4
T0, P0 = 288.15, 101325.0
RHO0 = P0 / (R_AIR * T0)                   # 1.225 kg/m^3
A0 = np.sqrt(GAMMA * R_AIR * T0)           # 340.29 m/s, sea-level speed of sound
KNOT = 1852.0 / 3600.0                     # m/s per knot (aviation unit)


# ---------------------------------------------------------------- atmosphere (troposphere only, given)

def isa_troposphere(h: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """ISA temperature (K), pressure (Pa) and density (kg/m^3) for 0 <= h <= 11 km (geopotential),
    the gradient-layer closed form of `P08.01`: T = T0 - 0.0065 h, p = p0 (T/T0)^(g0/(0.0065 R))."""
    h = np.asarray(h, dtype=float)
    if np.any(h < -500.0) or np.any(h > 11000.0):
        raise ValueError("troposphere only: -500 to 11 000 m")
    T = T0 - 0.0065 * h
    p = P0 * (T / T0) ** (G0 / (0.0065 * R_AIR))
    return T, p, p / (R_AIR * T)


# ---------------------------------------------------------------- Bernoulli

def dynamic_pressure(rho: float, V: np.ndarray) -> np.ndarray:
    """q = rho V^2 / 2 (Pa)."""
    return 0.5 * rho * np.asarray(V, dtype=float) ** 2


def bernoulli_airspeed(dp: np.ndarray, rho: float) -> np.ndarray:
    """Speed from a pitot-static pressure difference by incompressible Bernoulli, V = sqrt(2 dp / rho)."""
    dp = np.asarray(dp, dtype=float)
    if np.any(dp < 0):
        raise ValueError("a pitot-static difference cannot be negative in steady flow")
    return np.sqrt(2.0 * dp / rho)


def impact_pressure(V: np.ndarray, p: float, T: float) -> np.ndarray:
    """Pitot minus static pressure q_c (Pa) for subsonic isentropic deceleration (given; derived in
    `P08.17`): q_c = p [(1 + (gamma-1)/2 M^2)^(gamma/(gamma-1)) - 1], M = V / sqrt(gamma R T).
    For small M, q_c = q (1 + M^2/4 + ...), so Bernoulli underestimates it by about M^2/4."""
    M = np.asarray(V, dtype=float) / np.sqrt(GAMMA * R_AIR * T)
    if np.any(M >= 1.0):
        raise ValueError("subsonic only; above M = 1 a shock stands in front of the probe (P08.17)")
    return p * ((1.0 + 0.5 * (GAMMA - 1.0) * M ** 2) ** (GAMMA / (GAMMA - 1.0)) - 1.0)


def cas_from_impact_pressure(qc: np.ndarray) -> np.ndarray:
    """Calibrated airspeed (m/s): the speed that would give impact pressure q_c at sea-level standard
    conditions, the law an airspeed indicator is built to: CAS = a0 sqrt(5 [(q_c/p0 + 1)^(2/7) - 1])."""
    qc = np.asarray(qc, dtype=float)
    return A0 * np.sqrt(2.0 / (GAMMA - 1.0) * ((qc / P0 + 1.0) ** ((GAMMA - 1.0) / GAMMA) - 1.0))


def impact_pressure_from_cas(cas: np.ndarray) -> np.ndarray:
    """Inverse of ``cas_from_impact_pressure``: q_c = p0 [(1 + 0.2 (CAS/a0)^2)^3.5 - 1]."""
    return impact_pressure(cas, P0, T0)


def airspeeds_from_tas(tas: float, h: float, cp_static: float = 0.0, instrument_error: float = 0.0) -> dict[str, float]:
    """The airspeed chain at pressure altitude h (ISA day) for a true airspeed tas (m/s).

    The static port sits where the local pressure coefficient is cp_static, so it senses
    p_s = p + cp_static q (q = rho TAS^2 / 2); the pitot senses the total pressure p + q_c. The
    indicator works on the sensed difference q_c' = q_c - cp_static q and shows
    IAS = CAS(q_c') + instrument_error. Returns TAS, EAS, CAS (true q_c), IAS, the pressures and
    the Mach number, all SI."""
    T, p, rho = (float(v) for v in isa_troposphere(h))
    q = float(dynamic_pressure(rho, tas))
    qc = float(impact_pressure(tas, p, T))
    p_static_sensed = p + cp_static * q
    qc_sensed = qc - cp_static * q
    return {"TAS": float(tas), "EAS": float(tas * np.sqrt(rho / RHO0)), "CAS": float(cas_from_impact_pressure(qc)),
            "IAS": float(cas_from_impact_pressure(qc_sensed)) + instrument_error, "q": q, "qc": qc,
            "p": p, "p_static_sensed": p_static_sensed, "p_total": p + qc, "rho": rho, "T": T,
            "M": float(tas / np.sqrt(GAMMA * R_AIR * T))}


def tas_from_cas(cas: float, h: float) -> float:
    """True airspeed (m/s) from calibrated airspeed at pressure altitude h (ISA day): q_c from the
    calibration law, then the Mach number from q_c/p by inverting the isentropic relation, then
    TAS = M a. The inverse of the CAS step of ``airspeeds_from_tas``."""
    T, p, _ = (float(v) for v in isa_troposphere(h))
    qc = float(impact_pressure_from_cas(cas))
    M = np.sqrt(2.0 / (GAMMA - 1.0) * ((qc / p + 1.0) ** ((GAMMA - 1.0) / GAMMA) - 1.0))
    return float(M * np.sqrt(GAMMA * R_AIR * T))


def static_error_effects(cp_static: float, q: float, rho: float) -> dict[str, float]:
    """What a static-port pressure error dp = cp_static q does: the altimeter, which converts
    pressure to height with dp/dh = -rho g, reads high by dh = -dp/(rho g); the airspeed
    indicator sees q_c - dp. Returns {'dp', 'altitude_error'} in Pa and m."""
    dp = cp_static * q
    return {"dp": dp, "altitude_error": -dp / (rho * G0)}


# ---------------------------------------------------------------- Bernoulli in a duct: venturi and contraction

def venturi_flow_rate(dp: float, A1: float, A2: float, rho: float, Cd: float = 1.0) -> float:
    """Volume flow (m^3/s) through a venturi from the pressure drop between the inlet (A1) and the
    throat (A2): continuity V1 A1 = V2 A2 with Bernoulli gives Q = Cd A2 sqrt(2 dp / (rho (1 - (A2/A1)^2)))."""
    return Cd * A2 * np.sqrt(2.0 * dp / (rho * (1.0 - (A2 / A1) ** 2)))


def contraction_area(x: np.ndarray, A_in: float, A_out: float, Lc: float) -> np.ndarray:
    """Smooth wind-tunnel contraction area: a fifth-order polynomial blend from A_in at x = 0 to
    A_out at x = Lc with zero slope and curvature at both ends (given)."""
    s = np.clip(np.asarray(x, dtype=float) / Lc, 0.0, 1.0)
    blend = 10 * s ** 3 - 15 * s ** 4 + 6 * s ** 5
    return A_in + (A_out - A_in) * blend


def euler_pressure_along_duct(A_of_x: Callable[[np.ndarray], np.ndarray], x: np.ndarray, Q: float,
                              rho: float, p_in: float) -> np.ndarray:
    """Pressure along a duct by integrating the steady inviscid momentum (Euler) equation along the
    centre streamline, dp/dx = -rho V dV/dx with V = Q/A(x) from continuity, by the trapezoid rule
    on the grid x (V' by central differences). This does not use Bernoulli's algebra; it is the
    differential statement Bernoulli integrates."""
    V = Q / np.asarray(A_of_x(x), dtype=float)
    dVdx = np.gradient(V, x)
    f = -rho * V * dVdx
    return p_in + np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(x))])


def bernoulli_pressure_along_duct(A_of_x: Callable[[np.ndarray], np.ndarray], x: np.ndarray, Q: float,
                                  rho: float, p_in: float) -> np.ndarray:
    """Pressure along the duct from Bernoulli and continuity: p = p_in + rho (V_in^2 - V^2)/2."""
    V = Q / np.asarray(A_of_x(x), dtype=float)
    return p_in + 0.5 * rho * (V[0] ** 2 - V ** 2)


# ---------------------------------------------------------------- the energy equation

def shaft_head(p1: float, V1: float, z1: float, p2: float, V2: float, z2: float, rho: float,
               h_loss: float = 0.0, alpha1: float = 1.0, alpha2: float = 1.0, g: float = G0) -> float:
    """Head (m) a pump or fan must add between sections 1 and 2 (negative: a turbine extracts it),
    from the steady-flow energy equation with kinetic-energy correction factors alpha:
    h_shaft = (p2 - p1)/(rho g) + (alpha2 V2^2 - alpha1 V1^2)/(2 g) + (z2 - z1) + h_loss."""
    return (p2 - p1) / (rho * g) + (alpha2 * V2 ** 2 - alpha1 * V1 ** 2) / (2 * g) + (z2 - z1) + h_loss


def shaft_power(h_shaft: float, rho: float, Q: float, efficiency: float = 1.0, g: float = G0) -> float:
    """Shaft power (W) to supply a head h_shaft to a volume flow Q: rho g Q h / efficiency."""
    if not 0.0 < efficiency <= 1.0:
        raise ValueError("efficiency must be in (0, 1]")
    return rho * g * Q * h_shaft / efficiency


def open_tunnel_fan_power(V_test: float, A_test: float, rho: float, K_total: float, A_exit: float,
                          efficiency: float = 1.0) -> float:
    """Fan power (W) of an open-circuit wind tunnel drawing still air and discharging through an exit
    of area A_exit: the fan supplies the kinetic energy thrown away at the exit plus the losses,
    written as K_total test-section dynamic heads: P = rho Q [V_exit^2/2 + K_total V_test^2/2] / eta."""
    Q = V_test * A_test
    V_exit = Q / A_exit
    h = shaft_head(0.0, 0.0, 0.0, 0.0, V_exit, 0.0, rho, h_loss=K_total * V_test ** 2 / (2 * G0))
    return shaft_power(h, rho, Q, efficiency)


def observed_order(errors: np.ndarray, ratio: float = 2.0) -> np.ndarray:
    """Observed convergence orders log(e_k / e_{k+1}) / log(ratio) for grids refined by ``ratio``."""
    e = np.asarray(errors, dtype=float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


# ---------------------------------------------------------------- the lesson's run

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save pitot_static.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    h, tas = 3000.0, 60.0
    r = airspeeds_from_tas(tas, h)
    print(f"light aeroplane, {h:.0f} m ISA, TAS {tas} m/s ({tas / KNOT:.1f} kt), M = {r['M']:.4f}")
    print(f"  static p = {r['p']:.1f} Pa, rho = {r['rho']:.4f} kg/m^3; q = {r['q']:.2f} Pa; "
          f"impact pressure q_c = {r['qc']:.2f} Pa (q_c/q - 1 = {r['qc'] / r['q'] - 1:.5f}, M^2/4 = {r['M'] ** 2 / 4:.5f})")
    print(f"  EAS {r['EAS']:.3f} m/s ({r['EAS'] / KNOT:.1f} kt)   CAS {r['CAS']:.3f} m/s ({r['CAS'] / KNOT:.1f} kt)   "
          f"IAS (perfect ports) {r['IAS']:.3f} m/s")
    print(f"  inversion: tas_from_cas(CAS) = {tas_from_cas(r['CAS'], h):.12f} m/s")
    V_b = float(bernoulli_airspeed(r["qc"], r["rho"]))
    print(f"  Bernoulli with the true density on q_c: {V_b:.3f} m/s ({V_b / tas - 1:+.3%} vs TAS)")

    print("\nstatic-port error (port in air of pressure coefficient cp):")
    print("    cp      dp (Pa)   altimeter error (m)   IAS (m/s)   IAS - CAS (kt)")
    for cp in (-0.04, -0.02, 0.0, 0.02):
        rr = airspeeds_from_tas(tas, h, cp_static=cp)
        e = static_error_effects(cp, rr["q"], rr["rho"])
        print(f"   {cp:+.2f}   {e['dp']:+8.2f}       {e['altitude_error']:+7.2f}          {rr['IAS']:.3f}      "
              f"{(rr['IAS'] - rr['CAS']) / KNOT:+.2f}")

    print("\nBernoulli vs the compressible law: (Bernoulli V - TAS)/TAS from q_c at 3000 m")
    for V in (30.0, 60.0, 100.0, 150.0, 200.0):
        rr = airspeeds_from_tas(V, h)
        print(f"   TAS {V:5.0f} m/s  M = {rr['M']:.3f}  error {float(bernoulli_airspeed(rr['qc'], rr['rho'])) / V - 1:+.4%}")

    # wind-tunnel contraction: Euler integration vs Bernoulli
    A_in, A_out, Lc, rho = 4.0, 0.5, 2.0, 1.204
    Q = 30.0 * A_out
    print(f"\ncontraction {A_in} -> {A_out} m^2 over {Lc} m, 30 m/s in the test section, air 1.204 kg/m^3:")
    errs = []
    A_fun = lambda xx: contraction_area(xx, A_in, A_out, Lc)
    for n in (25, 50, 100, 200):
        x = np.linspace(0.0, Lc, n + 1)
        pe = euler_pressure_along_duct(A_fun, x, Q, rho, 0.0)
        pb = bernoulli_pressure_along_duct(A_fun, x, Q, rho, 0.0)
        errs.append(np.max(np.abs(pe - pb)))
        print(f"   n = {n:4d}: max |Euler - Bernoulli| = {errs[-1]:.3e} Pa (pressure drop {-pb[-1]:.2f} Pa)")
    print(f"   observed orders {np.array2string(observed_order(errs), precision=2)}")
    dp_c = -pb[-1]
    print(f"   venturi formula inverted from that drop: Q = {venturi_flow_rate(dp_c, A_in, A_out, rho):.6f} m^3/s (true {Q})")

    P = open_tunnel_fan_power(30.0, 0.5, 1.204, 0.6, 1.0, efficiency=0.7)
    print(f"\nopen-circuit tunnel, 30 m/s in 0.5 m^2, diffuser exit 1.0 m^2, losses 0.6 q_test, fan 70 %: shaft power {P:.0f} W")

    hs = np.linspace(0, 6000, 61)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for name in ("EAS", "CAS", "TAS"):
        ax[0].plot(hs, [airspeeds_from_tas(tas, hh)[name] if name != "TAS" else tas for hh in hs], label=name)
    ax[0].set_xlabel("pressure altitude (m)")
    ax[0].set_ylabel("airspeed (m/s)")
    ax[0].set_title("TAS 60 m/s: what the dial reads falls with height")
    ax[0].legend()
    x = np.linspace(0.0, Lc, 201)
    ax[1].plot(x, bernoulli_pressure_along_duct(A_fun, x, Q, rho, 0.0), label="Bernoulli")
    x25 = np.linspace(0.0, Lc, 26)
    ax[1].plot(x25, euler_pressure_along_duct(A_fun, x25, Q, rho, 0.0), "o", ms=3, mfc="none", label="Euler, n = 25")
    ax[1].set_xlabel("x along contraction (m)")
    ax[1].set_ylabel("p - p_in (Pa)")
    ax[1].set_title("wind-tunnel contraction")
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(out_dir / "pitot_static.png", dpi=110)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'pitot_static.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
