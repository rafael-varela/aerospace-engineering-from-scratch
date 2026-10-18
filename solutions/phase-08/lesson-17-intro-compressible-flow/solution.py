"""P08.17 Intro to Compressible Flow: reference solution.

How fast is the light aeroplane flying compared with the speed of sound, how wrong is the
incompressible Bernoulli airspeed of P08.04 there, at what Mach number does that error reach 1 %,
and what does a pitot probe read on a sounding rocket at Mach 2? This module answers with:

* a small International Standard Atmosphere to 47 km (the layers of P08.01, reimplemented here so the
  module stands alone) and the speed of sound a = sqrt(gamma R T);
* the isentropic total-to-static relations T0/T, p0/p, rho0/rho as functions of Mach number, their
  inversion, and the low-Mach limit in which p0 - p tends to the Bernoulli dynamic pressure;
* the compressible (Saint-Venant-Wantzel) pitot formula below M = 1 and the error of the
  incompressible formula;
* the normal shock relations, and the Rayleigh pitot formula for the total pressure behind the bow
  shock of a probe in supersonic flow, inverted for M by Brent's method.

The air is a calorically perfect gas with gamma = 1.4 and R = 287.05287 J/(kg K); the flow outside the
shock is adiabatic and reversible (isentropic); the shock is a thin adiabatic discontinuity.

Units: SI (K, Pa, kg/m^3, m/s).

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes compressible.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-08/lesson-17-intro-compressible-flow
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

GAMMA = 1.4
R_AIR = 287.05287          # J/(kg K)
G0 = 9.80665               # m/s^2
T_SL, P_SL = 288.15, 101325.0
CP_AIR = GAMMA * R_AIR / (GAMMA - 1.0)

# ISA layers: base geopotential altitude (m), lapse rate (K/m)
_ISA_LAYERS = [(0.0, -6.5e-3), (11000.0, 0.0), (20000.0, 1.0e-3), (32000.0, 2.8e-3), (47000.0, None)]


# ---------------------------------------------------------------- atmosphere and the speed of sound

def isa(h: float) -> tuple[float, float, float]:
    """Temperature (K), pressure (Pa) and density (kg/m^3) of the ISA at geopotential altitude
    0 <= h <= 47 000 m, layer by layer: T = Tb + L (h - hb); p = pb (T/Tb)^(-g0/(L R)) in a gradient
    layer and pb exp(-g0 (h - hb)/(R Tb)) in an isothermal one; rho = p/(R T)."""
    if not 0.0 <= h <= 47000.0:
        raise ValueError("this ISA covers 0 to 47 km geopotential altitude")
    Tb, pb = T_SL, P_SL
    for (hb, L), (htop, _) in zip(_ISA_LAYERS[:-1], _ISA_LAYERS[1:]):
        hh = min(h, htop) - hb
        if L == 0.0:
            T = Tb
            p = pb * np.exp(-G0 * hh / (R_AIR * Tb))
        else:
            T = Tb + L * hh
            p = pb * (T / Tb) ** (-G0 / (L * R_AIR))
        if h <= htop:
            return T, p, p / (R_AIR * T)
        Tb, pb = T, p
    raise AssertionError("unreachable")


def speed_of_sound(T: np.ndarray, gamma: float = GAMMA, R: float = R_AIR) -> np.ndarray:
    """Speed of sound of a perfect gas, a = sqrt(gamma R T) (m/s)."""
    return np.sqrt(gamma * R * np.asarray(T, dtype=float))


# ---------------------------------------------------------------- isentropic relations

def T0_over_T(M: np.ndarray, gamma: float = GAMMA) -> np.ndarray:
    """Total-to-static temperature ratio, 1 + (gamma - 1)/2 M^2 (adiabatic flow; no isentropy needed)."""
    M = np.asarray(M, dtype=float)
    return 1.0 + 0.5 * (gamma - 1.0) * M ** 2


def p0_over_p(M: np.ndarray, gamma: float = GAMMA) -> np.ndarray:
    """Isentropic total-to-static pressure ratio, (T0/T)^(gamma/(gamma - 1))."""
    return T0_over_T(M, gamma) ** (gamma / (gamma - 1.0))


def rho0_over_rho(M: np.ndarray, gamma: float = GAMMA) -> np.ndarray:
    """Isentropic total-to-static density ratio, (T0/T)^(1/(gamma - 1))."""
    return T0_over_T(M, gamma) ** (1.0 / (gamma - 1.0))


def mach_from_p0_over_p(ratio: np.ndarray, gamma: float = GAMMA) -> np.ndarray:
    """Inverse of ``p0_over_p``: M = sqrt(2/(gamma - 1) (ratio^((gamma - 1)/gamma) - 1))."""
    ratio = np.asarray(ratio, dtype=float)
    return np.sqrt(2.0 / (gamma - 1.0) * (ratio ** ((gamma - 1.0) / gamma) - 1.0))


def pressure_coefficient_ratio(M: np.ndarray, gamma: float = GAMMA) -> np.ndarray:
    """(p0 - p) / (rho V^2 / 2) for isentropic stagnation, = 2 (p0/p - 1) / (gamma M^2); it tends to 1
    (Bernoulli) as M -> 0, like 1 + M^2/4 + (2 - gamma) M^4 / 24. p0/p - 1 is evaluated as
    expm1(gamma/(gamma - 1) log1p((gamma - 1) M^2 / 2)), which keeps full precision at small M."""
    M = np.asarray(M, dtype=float)
    excess = np.expm1(gamma / (gamma - 1.0) * np.log1p(0.5 * (gamma - 1.0) * M ** 2))
    return 2.0 * excess / (gamma * M ** 2)


# ---------------------------------------------------------------- airspeed from a pitot-static probe

def bernoulli_airspeed(qc: float, rho: float) -> float:
    """Incompressible airspeed from the impact pressure qc = p0 - p: V = sqrt(2 qc / rho) (P08.04)."""
    return float(np.sqrt(2.0 * qc / rho))


def subsonic_pitot_airspeed(qc: float, p: float, T: float, gamma: float = GAMMA) -> tuple[float, float]:
    """True airspeed and Mach number from the impact pressure qc = p0 - p, static pressure p and static
    temperature T, for subsonic flow (isentropic compression into the probe):
    M = mach_from_p0_over_p(1 + qc/p), V = M a(T)."""
    M = float(mach_from_p0_over_p(1.0 + qc / p, gamma))
    return M * float(speed_of_sound(T, gamma)), M


def bernoulli_error(M: np.ndarray, gamma: float = GAMMA) -> np.ndarray:
    """Relative error V_Bernoulli / V - 1 of the incompressible formula, using the true impact pressure
    and the true static density: sqrt((p0 - p)/(rho V^2/2)) - 1, about M^2/8 at low Mach number."""
    return np.sqrt(pressure_coefficient_ratio(M, gamma)) - 1.0


def mach_for_bernoulli_error(err: float, gamma: float = GAMMA) -> float:
    """Mach number at which the incompressible airspeed is too high by the fraction ``err`` (subsonic)."""
    return brentq(lambda M: float(bernoulli_error(M, gamma)) - err, 1e-4, 1.0, xtol=1e-14)


def stagnation_temperature(T: float, V: float, cp: float = CP_AIR) -> float:
    """Total temperature from the steady adiabatic energy equation, T0 = T + V^2 / (2 cp)."""
    return T + V ** 2 / (2.0 * cp)


# ---------------------------------------------------------------- normal shock and the Rayleigh pitot formula

def normal_shock(M1: float, gamma: float = GAMMA) -> dict[str, float]:
    """Normal shock relations for upstream Mach number M1 >= 1 (from mass, momentum and energy across a
    thin adiabatic discontinuity): M2, p2/p1, rho2/rho1, T2/T1 and the total-pressure ratio p02/p01."""
    if M1 < 1.0:
        raise ValueError("a normal shock needs M1 >= 1")
    g = gamma
    M1s = M1 ** 2
    M2 = np.sqrt((1 + 0.5 * (g - 1) * M1s) / (g * M1s - 0.5 * (g - 1)))
    p21 = 1 + 2 * g / (g + 1) * (M1s - 1)
    r21 = (g + 1) * M1s / ((g - 1) * M1s + 2)
    T21 = p21 / r21
    p0201 = p21 * p0_over_p(M2, g) / p0_over_p(M1, g)
    return {"M2": float(M2), "p2_p1": float(p21), "rho2_rho1": float(r21), "T2_T1": float(T21),
            "p02_p01": float(p0201)}


def rayleigh_pitot(M1: np.ndarray, gamma: float = GAMMA) -> np.ndarray:
    """Rayleigh's supersonic pitot formula, the ratio of the pitot (total) pressure behind a normal
    shock to the free-stream static pressure:

        p02/p1 = [(gamma + 1)^2 M1^2 / (4 gamma M1^2 - 2 (gamma - 1))]^(gamma/(gamma - 1))
                 * (1 - gamma + 2 gamma M1^2) / (gamma + 1)."""
    M1 = np.asarray(M1, dtype=float)
    g = gamma
    a = ((g + 1) ** 2 * M1 ** 2 / (4 * g * M1 ** 2 - 2 * (g - 1))) ** (g / (g - 1))
    return a * (1 - g + 2 * g * M1 ** 2) / (g + 1)


def mach_from_pitot(p0_over_p_measured: float, gamma: float = GAMMA) -> float:
    """Mach number from the ratio of pitot pressure to static pressure: the isentropic inversion when
    the ratio is at most p0/p at M = 1 (1.8929 for gamma = 1.4), otherwise the root of the Rayleigh
    formula by Brent's method (it increases monotonically with M1)."""
    r = float(p0_over_p_measured)
    if r < 1.0:
        raise ValueError("pitot pressure cannot be below static pressure")
    if r <= float(p0_over_p(1.0, gamma)):
        return float(mach_from_p0_over_p(r, gamma))
    return brentq(lambda M: float(rayleigh_pitot(M, gamma)) - r, 1.0, 50.0, xtol=1e-14)


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save compressible.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    T, p, rho = isa(3000.0)
    a = float(speed_of_sound(T))
    V = 60.0
    M = V / a
    qc = p * (float(p0_over_p(M)) - 1.0)
    VB = bernoulli_airspeed(qc, rho)
    Vc, Mc = subsonic_pitot_airspeed(qc, p, T)
    print(f"ISA 3000 m: T = {T:.3f} K, p = {p:.1f} Pa, rho = {rho:.5f} kg/m^3, a = {a:.3f} m/s")
    print(f"light aeroplane at {V} m/s: M = {M:.4f}; impact pressure qc = {qc:.2f} Pa, rho V^2/2 = {0.5 * rho * V ** 2:.2f} Pa")
    print(f"  Bernoulli airspeed {VB:.4f} m/s (error {(VB / V - 1) * 100:.4f} %, M^2/8 = {M ** 2 / 8 * 100:.4f} %); "
          f"compressible formula {Vc:.6f} m/s")
    print(f"  total temperature rise V^2/(2 cp) = {stagnation_temperature(T, V) - T:.3f} K")
    for err in (0.001, 0.01, 0.05):
        Me = mach_for_bernoulli_error(err)
        print(f"  Bernoulli error reaches {err * 100:g} % at M = {Me:.4f} ({Me * a:.1f} m/s at 3000 m, "
              f"{Me * float(speed_of_sound(T_SL)):.1f} m/s at sea level)")
    T0, p0, rho0 = isa(0.0)
    print(f"  sea level: a = {float(speed_of_sound(T0)):.3f} m/s; p0/p at M = 1: {float(p0_over_p(1.0)):.5f}")

    # the sounding rocket
    Tr, pr, rhor = isa(11000.0)
    ar = float(speed_of_sound(Tr))
    M1 = 2.0
    V1 = M1 * ar
    ratio = float(rayleigh_pitot(M1))
    sh = normal_shock(M1)
    p02 = ratio * pr
    print(f"\nsounding rocket at M = {M1} at 11 km: T = {Tr:.2f} K, p = {pr:.1f} Pa, a = {ar:.2f} m/s, V = {V1:.2f} m/s")
    print(f"  normal shock: M2 = {sh['M2']:.5f}, p2/p1 = {sh['p2_p1']:.4f}, T2/T1 = {sh['T2_T1']:.5f}, "
          f"p02/p01 = {sh['p02_p01']:.5f}")
    print(f"  Rayleigh pitot p02/p1 = {ratio:.5f} -> pitot reads {p02:.1f} Pa; "
          f"isentropic (no shock) p01/p1 would be {float(p0_over_p(M1)):.5f}")
    print(f"  inverted: M = {mach_from_pitot(ratio):.10f}")
    M_wrong = float(mach_from_p0_over_p(ratio))
    VB_r = bernoulli_airspeed(p02 - pr, rhor)
    print(f"  subsonic formula on that reading: M = {M_wrong:.4f} ({(M_wrong / M1 - 1) * 100:.2f} %); "
          f"Bernoulli: V = {VB_r:.1f} m/s ({(VB_r / V1 - 1) * 100:.2f} %)")
    print(f"  total temperature {stagnation_temperature(Tr, V1):.2f} K")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    Ms = np.linspace(0.01, 1.0, 200)
    ax[0].plot(Ms, bernoulli_error(Ms) * 100, label="exact (isentropic)")
    ax[0].plot(Ms, Ms ** 2 / 8 * 100, "--", label="M^2 / 8")
    ax[0].axhline(1.0, color="k", lw=0.8, ls=":")
    ax[0].axvline(M, color="C2", lw=0.8, label="light aeroplane cruise")
    ax[0].set_xlabel("M")
    ax[0].set_ylabel("Bernoulli airspeed error (%)")
    ax[0].set_title("where incompressible Bernoulli fails")
    ax[0].legend(fontsize=8)
    Mt = np.linspace(0.01, 3.0, 300)
    sub = Mt <= 1
    ax[1].semilogy(Mt[sub], p0_over_p(Mt[sub]), "C0", label="isentropic p0/p (subsonic pitot)")
    ax[1].semilogy(Mt[~sub], p0_over_p(Mt[~sub]), "C0:", label="isentropic p0/p (no shock)")
    ax[1].semilogy(Mt[~sub], rayleigh_pitot(Mt[~sub]), "C1", label="Rayleigh pitot p02/p1")
    ax[1].set_xlabel("M")
    ax[1].set_ylabel("pitot / static pressure")
    ax[1].set_title("pitot reading")
    ax[1].legend(fontsize=8)
    hs = np.linspace(0, 47000, 200)
    ax[2].plot([float(speed_of_sound(isa(h)[0])) for h in hs], hs / 1e3)
    ax[2].set_xlabel("speed of sound (m/s)")
    ax[2].set_ylabel("geopotential altitude (km)")
    ax[2].set_title("ISA speed of sound")
    fig.tight_layout()
    fig.savefig(out_dir / "compressible.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'compressible.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
