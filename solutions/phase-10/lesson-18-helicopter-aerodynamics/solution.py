"""P10.18 Helicopter Aerodynamics: reference solution.

How much power does the reference light helicopter (1100 kg, R = 5 m, two blades of 0.25 m chord,
tip speed 210 m/s, NACA 0012) need to hover at sea level, at what airspeed is the power least and
the range greatest, and how fast does it come down in a power-off glide (autorotation)? This module
answers with a hierarchy of rotor models, each checked against the one beside it:

* momentum theory in hover: induced velocity v_h = sqrt(T / (2 rho A)), ideal power T v_h, and the
  figure of merit FM = P_ideal / P_actual;
* blade element theory with uniform inflow: C_T = (sigma a / 2)(theta_0/3 + theta_tw/4 - lambda/2),
  closed form and by strip integration, giving the collective pitch; the induced-power factor kappa
  and the profile power sigma C_d0 / 8;
* Glauert's forward-flight inflow, lambda = mu tan(alpha) + C_T / (2 sqrt(mu^2 + lambda^2)), solved
  by Newton iteration, and the power curve C_P = kappa C_T lambda_i + (sigma C_d0 / 8)(1 + K mu^2)
  + (f / 2A) mu^3, where K = 4.65 (the profile term with the drag of the radial flow; K = 3 without
  it, derived and checked numerically here);
* the minimum-power and maximum-range speeds, and the autorotative rate of descent from the energy
  balance W V_d = P_required.

Nondimensionalisation: C_T = T / (rho A (Omega R)^2), C_P = P / (rho A (Omega R)^3),
lambda = (V sin(alpha) + v) / (Omega R) (positive down through the disk), mu = V cos(alpha) / (Omega R),
alpha the tip-path-plane angle (positive nose down, tilting the thrust forward). Units SI.

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes helicopter_power.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-10/lesson-18-helicopter-aerodynamics
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar

RHO = 1.225          # kg/m^3, ISA sea level
G = 9.81             # m/s^2


@dataclass
class Rotor:
    """The reference rotor and the stated modelling assumptions.

    a: lift-curve slope of the NACA 0012 blade section (per rad), a thin-airfoil value reduced for
    viscosity. cd0: mean profile drag coefficient (assumed; NACA 0012 at Re ~ 3e6 has cd 0.006-0.008
    smooth and about 0.010 with leading-edge roughness, so 0.010 allows for real blades). kappa:
    induced-power factor for non-uniform inflow and tip loss (assumed 1.15). f: the airframe's
    equivalent flat-plate drag area (assumed 1.0 m^2 for a light helicopter). K: profile-power
    advance-ratio factor. theta_tw: linear blade twist (rad, tip minus root).
    """
    mass: float = 1100.0
    R: float = 5.0
    B: int = 2
    c: float = 0.25
    tip_speed: float = 210.0
    a: float = 5.73
    cd0: float = 0.010
    kappa: float = 1.15
    f: float = 1.0
    K: float = 4.65
    theta_tw: float = np.radians(-8.0)

    @property
    def A(self) -> float:
        """Disk area (m^2)."""
        return np.pi * self.R ** 2

    @property
    def sigma(self) -> float:
        """Solidity B c / (pi R)."""
        return self.B * self.c / (np.pi * self.R)

    @property
    def W(self) -> float:
        """Weight (N), equal to the thrust in hover and (to first order) in level flight."""
        return self.mass * G

    @property
    def Omega(self) -> float:
        """Rotor speed (rad/s)."""
        return self.tip_speed / self.R


def thrust_coefficient(T: float, rotor: Rotor, rho: float = RHO) -> float:
    """C_T = T / (rho A (Omega R)^2)."""
    return T / (rho * rotor.A * rotor.tip_speed ** 2)


def power_from_coefficient(CP: float | np.ndarray, rotor: Rotor, rho: float = RHO) -> np.ndarray:
    """P = C_P rho A (Omega R)^3 (W)."""
    return np.asarray(CP) * rho * rotor.A * rotor.tip_speed ** 3


# ---------------------------------------------------------------- hover: momentum theory

def hover_momentum(T: float, A: float, rho: float = RHO) -> dict:
    """Momentum theory in hover: v_h = sqrt(T / (2 rho A)), P_ideal = T v_h, far-wake velocity 2 v_h,
    wake contraction to half the disk area. Returns dict(v_h, P_ideal, v_wake, A_wake)."""
    v = np.sqrt(T / (2 * rho * A))
    return {"v_h": v, "P_ideal": T * v, "v_wake": 2 * v, "A_wake": A / 2}


def figure_of_merit(T: float, P: float, A: float, rho: float = RHO) -> float:
    """FM = ideal power / actual power = T^(3/2) / (sqrt(2 rho A) P)."""
    return T ** 1.5 / (np.sqrt(2 * rho * A) * P)


# ---------------------------------------------------------------- hover: blade element

def be_thrust_coefficient(theta0: float, lam: float, sigma: float, a: float, theta_tw: float = 0.0) -> float:
    """Blade-element C_T with uniform inflow and linear twist theta(r) = theta0 + theta_tw r:
    dC_T = (sigma a / 2)(theta r^2 - lambda r) dr integrated over 0 <= r <= 1 gives
    C_T = (sigma a / 2)(theta0 / 3 + theta_tw / 4 - lambda / 2)."""
    return 0.5 * sigma * a * (theta0 / 3 + theta_tw / 4 - lam / 2)


def be_thrust_strips(theta0: float, lam: float, sigma: float, a: float, theta_tw: float = 0.0,
                     n: int = 200) -> float:
    """The same C_T by strip integration (midpoint rule over n strips): the numerical cross-check."""
    r = (np.arange(n) + 0.5) / n
    theta = theta0 + theta_tw * r
    return float(np.sum(0.5 * sigma * a * (theta * r ** 2 - lam * r)) / n)


def hover_collective(CT: float, rotor: Rotor) -> tuple[float, float]:
    """Collective pitch theta0 (at the root, rad) that gives C_T in hover with the momentum inflow
    lambda_h = sqrt(C_T / 2): theta0 = 3 (2 C_T / (sigma a) - theta_tw / 4 + lambda / 2).
    Returns (theta0, theta_75) with theta_75 the pitch at 0.75 R."""
    lam = np.sqrt(CT / 2)
    th0 = 3 * (2 * CT / (rotor.sigma * rotor.a) - rotor.theta_tw / 4 + lam / 2)
    return th0, th0 + 0.75 * rotor.theta_tw


def hover_power(rotor: Rotor, T: float | None = None, rho: float = RHO) -> dict:
    """Hover power coefficient C_P = kappa C_T^(3/2) / sqrt(2) + sigma C_d0 / 8 and its parts (W)."""
    T = rotor.W if T is None else T
    CT = thrust_coefficient(T, rotor, rho)
    CPi = rotor.kappa * CT ** 1.5 / np.sqrt(2)
    CP0 = rotor.sigma * rotor.cd0 / 8
    P_i, P_0 = power_from_coefficient(CPi, rotor, rho), power_from_coefficient(CP0, rotor, rho)
    P = P_i + P_0
    return {"CT": CT, "CPi": CPi, "CP0": CP0, "CP": CPi + CP0, "P_induced": float(P_i),
            "P_profile": float(P_0), "P": float(P), "FM": figure_of_merit(T, float(P), rotor.A, rho),
            "CT_over_sigma": CT / rotor.sigma}


# ---------------------------------------------------------------- forward flight

def glauert_inflow(mu: float, CT: float, alpha_tpp: float = 0.0, tol: float = 1e-13,
                   max_iter: int = 100) -> tuple[float, float, int]:
    """Solve lambda = mu tan(alpha) + C_T / (2 sqrt(mu^2 + lambda^2)) by Newton's method from the
    hover value. Returns (lambda, lambda_i, iterations), lambda_i = C_T / (2 sqrt(mu^2 + lambda^2))."""
    lam = np.sqrt(CT / 2) + mu * np.tan(alpha_tpp)
    for it in range(1, max_iter + 1):
        s = np.sqrt(mu ** 2 + lam ** 2)
        f = lam - mu * np.tan(alpha_tpp) - CT / (2 * s)
        fp = 1 + CT * lam / (2 * s ** 3)
        step = f / fp
        lam -= step
        if abs(step) < tol:
            break
    lam_i = CT / (2 * np.sqrt(mu ** 2 + lam ** 2))
    return float(lam), float(lam_i), it


def profile_power_integral(mu: float, radial_flow: bool = True, n: int = 200) -> float:
    """Profile power of a rotor with constant chord and cd, normalised by its hover value:
    (1 / 2 pi) int int U^3 dr dpsi / (1/4), with U = |u_T| or sqrt(u_T^2 + u_R^2),
    u_T = r + mu sin(psi), u_R = mu cos(psi). Without radial flow and with u_T^3 taken signed the
    result is exactly 1 + 3 mu^2; with the radial flow it is close to 1 + 4.65 mu^2."""
    r = (np.arange(n) + 0.5) / n
    psi = (np.arange(2 * n) + 0.5) / (2 * n) * 2 * np.pi
    rr, pp = np.meshgrid(r, psi)
    ut = rr + mu * np.sin(pp)
    ur = mu * np.cos(pp)
    U3 = (ut ** 2 + ur ** 2) ** 1.5 if radial_flow else np.abs(ut) ** 3
    return float(np.mean(U3) / 0.25)


def forward_flight_power(V: float | np.ndarray, rotor: Rotor, rho: float = RHO, climb: float = 0.0) -> dict:
    """Power required in level flight (or climb at rate `climb`, m/s) at airspeed V.

    Trim: the rotor thrust balances weight and drag, T = sqrt(W^2 + D^2) with D = 0.5 rho V^2 f, so the
    tip-path plane tilts forward by alpha with tan(alpha) = D / W. Then
    C_P = kappa C_T lambda_i + (sigma C_d0 / 8)(1 + K mu^2) + D V / (rho A (Omega R)^3) + C_T lambda_c,
    with mu = V cos(alpha) / (Omega R), lambda from ``glauert_inflow`` and the parasite term
    D V / (rho A (Omega R)^3) = (f / 2A)(V / Omega R)^3 (about (f / 2A) mu^3).
    Returns arrays: mu, lambda, lambda_i, alpha, C_T, the C_P parts and the powers (W).
    """
    V = np.atleast_1d(np.asarray(V, dtype=float))
    keys = ("mu", "lam", "lam_i", "alpha", "CT", "CPi", "CP0", "CPp", "CPc", "CP")
    out = {k: np.zeros_like(V) for k in keys}
    for j, v in enumerate(V):
        D = 0.5 * rho * v ** 2 * rotor.f
        alpha = np.arctan2(D, rotor.W)
        CT = thrust_coefficient(np.hypot(rotor.W, D), rotor, rho)
        mu = v * np.cos(alpha) / rotor.tip_speed
        lam, lam_i, _ = glauert_inflow(mu, CT, alpha)
        CPi = rotor.kappa * CT * lam_i
        CP0 = rotor.sigma * rotor.cd0 / 8 * (1 + rotor.K * mu ** 2)
        CPp = 0.5 * rotor.f / rotor.A * (v / rotor.tip_speed) ** 3
        CPc = CT * climb / rotor.tip_speed
        for k, val in zip(keys, (mu, lam, lam_i, alpha, CT, CPi, CP0, CPp, CPc, CPi + CP0 + CPp + CPc)):
            out[k][j] = val
    scale = rho * rotor.A * rotor.tip_speed ** 3
    out.update({"V": V, "P_induced": out["CPi"] * scale, "P_profile": out["CP0"] * scale,
                "P_parasite": out["CPp"] * scale, "P": out["CP"] * scale})
    return out


def best_speeds(rotor: Rotor, rho: float = RHO, V_max: float = 80.0) -> dict:
    """Minimum-power speed (best endurance, least autorotative sink) and maximum-range speed
    (least power per unit speed, the tangent from the origin to the power curve)."""
    P = lambda v: float(forward_flight_power(v, rotor, rho)["P"][0])
    r1 = minimize_scalar(P, bounds=(1.0, V_max), method="bounded", options={"xatol": 1e-6})
    r2 = minimize_scalar(lambda v: P(v) / v, bounds=(5.0, V_max), method="bounded", options={"xatol": 1e-6})
    return {"V_min_power": r1.x, "P_min": r1.fun, "V_max_range": r2.x, "P_max_range": P(r2.x),
            "P_over_V_min": r2.fun}


def autorotation(rotor: Rotor, rho: float = RHO) -> dict:
    """Energy-balance estimate of steady autorotation: with no engine power, the descent supplies it,
    W V_d = P_required(V), so the least rate of descent is P_min / W at the minimum-power speed and
    the flattest glide is at the maximum-range speed (glide ratio V / V_d = W V / P). The descent
    changes the inflow slightly; that and the transmission losses are neglected."""
    bs = best_speeds(rotor, rho)
    Vd_min = bs["P_min"] / rotor.W
    Vd_range = bs["P_max_range"] / rotor.W
    return {"Vd_min": Vd_min, "V_at_Vd_min": bs["V_min_power"], "Vd_best_glide": Vd_range,
            "V_best_glide": bs["V_max_range"], "glide_ratio": bs["V_max_range"] / Vd_range}


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save helicopter_power.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    rot = Rotor()
    T = rot.W
    print(f"rotor: W = {T:.1f} N, A = {rot.A:.3f} m^2, disk loading {T / rot.A:.1f} N/m^2, sigma = {rot.sigma:.5f}, "
          f"Omega = {rot.Omega:.1f} rad/s ({rot.Omega * 60 / (2 * np.pi):.0f} rpm)")
    hm = hover_momentum(T, rot.A)
    CT = thrust_coefficient(T, rot)
    print(f"momentum hover: v_h = {hm['v_h']:.3f} m/s, lambda_h = {np.sqrt(CT / 2):.5f}, ideal power "
          f"{hm['P_ideal'] / 1e3:.2f} kW, C_T = {CT:.6f}, C_T/sigma = {CT / rot.sigma:.4f}")
    hp = hover_power(rot)
    print(f"hover power: induced {hp['P_induced'] / 1e3:.2f} kW (kappa = {rot.kappa}), profile "
          f"{hp['P_profile'] / 1e3:.2f} kW (sigma Cd0/8, Cd0 = {rot.cd0}), total {hp['P'] / 1e3:.2f} kW, FM = {hp['FM']:.3f}")
    th0, th75 = hover_collective(CT, rot)
    lam = np.sqrt(CT / 2)
    ct_cf = be_thrust_coefficient(th0, lam, rot.sigma, rot.a, rot.theta_tw)
    ct_st = be_thrust_strips(th0, lam, rot.sigma, rot.a, rot.theta_tw, 200)
    print(f"blade element, uniform inflow: theta0 = {np.degrees(th0):.3f} deg, theta_75 = {np.degrees(th75):.3f} deg "
          f"(twist {np.degrees(rot.theta_tw):.0f} deg); C_T closed form {ct_cf:.7f}, strips {ct_st:.7f}")
    # mean lift coefficient: C_T/sigma = C_L_bar / 6
    print(f"mean blade lift coefficient 6 C_T/sigma = {6 * CT / rot.sigma:.3f}")
    print("figure of merit vs thrust (same rotor):")
    for m in (700, 900, 1100, 1300):
        h = hover_power(rot, m * G)
        print(f"  {m:5d} kg: C_T/sigma = {h['CT_over_sigma']:.4f}, P = {h['P'] / 1e3:.1f} kW, FM = {h['FM']:.3f}")

    # forward flight
    print("Glauert inflow (alpha = 0): lambda_i / lambda_h and the high-speed limit C_T/(2 mu):")
    for mu in (0.0, 0.02, 0.05, 0.1, 0.2, 0.3):
        lam_, lam_i, it = glauert_inflow(mu, CT, 0.0)
        hi = CT / (2 * mu) if mu > 0 else float("nan")
        print(f"  mu = {mu:.2f}: lambda_i = {lam_i:.5f} ({lam_i / np.sqrt(CT / 2):.4f} lambda_h), "
              f"C_T/(2 mu) = {hi:.5f}, Newton iterations {it}")
    print("profile-power factor by integration over the disk (exact for the first column):")
    for mu in (0.1, 0.2, 0.3, 0.4):
        k3 = (profile_power_integral(mu, radial_flow=False) - 1) / mu ** 2
        kr = (profile_power_integral(mu, radial_flow=True) - 1) / mu ** 2
        print(f"  mu = {mu:.1f}: K without radial flow = {k3:.3f}, with radial flow = {kr:.3f}")
    V = np.linspace(0.0, 75.0, 151)
    ff = forward_flight_power(V, rot)
    bs = best_speeds(rot)
    print(f"minimum power {bs['P_min'] / 1e3:.2f} kW at V = {bs['V_min_power']:.2f} m/s "
          f"({bs['V_min_power'] / 0.5144:.0f} kn); maximum range at V = {bs['V_max_range']:.2f} m/s "
          f"({bs['V_max_range'] / 0.5144:.0f} kn), P = {bs['P_max_range'] / 1e3:.2f} kW")
    for v in (0.0, 20.0, 40.0, 60.0, 70.0):
        p = forward_flight_power(v, rot)
        print(f"  V = {v:4.0f} m/s: mu = {p['mu'][0]:.4f}, induced {p['P_induced'][0] / 1e3:6.2f} kW, profile "
              f"{p['P_profile'][0] / 1e3:6.2f} kW, parasite {p['P_parasite'][0] / 1e3:6.2f} kW, total {p['P'][0] / 1e3:6.2f} kW")
    pk3 = forward_flight_power(bs["V_max_range"], Rotor(K=3.0))
    print(f"with K = 3 instead of 4.65 at the max-range speed: P = {pk3['P'][0] / 1e3:.2f} kW")
    au = autorotation(rot)
    print(f"autorotation (energy balance): least descent rate {au['Vd_min']:.2f} m/s ({au['Vd_min'] * 196.85:.0f} ft/min) "
          f"at {au['V_at_Vd_min']:.1f} m/s; best glide {au['V_best_glide']:.1f} m/s, descent {au['Vd_best_glide']:.2f} m/s, "
          f"glide ratio {au['glide_ratio']:.2f}")
    # energy balance check at 50 m/s: C_T lambda equals induced (kappa = 1) + parasite
    p = forward_flight_power(50.0, Rotor(kappa=1.0))
    lhs = p["CT"][0] * p["lam"][0]
    rhs = p["CPi"][0] + p["CPp"][0]
    print(f"energy balance at 50 m/s (kappa = 1): C_T lambda = {lhs:.4e}, C_Pi + C_Pp = {rhs:.4e}")

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    ax[0].plot(V, ff["P_induced"] / 1e3, label="induced")
    ax[0].plot(V, ff["P_profile"] / 1e3, label="profile (K = 4.65)")
    ax[0].plot(V, ff["P_parasite"] / 1e3, label="parasite (f = 1 m^2)")
    ax[0].plot(V, ff["P"] / 1e3, "k", lw=2, label="total")
    ax[0].plot([0, bs["V_max_range"] * 1.6], [0, bs["P_over_V_min"] * bs["V_max_range"] * 1.6 / 1e3], "k:", lw=1,
               label="tangent: max range")
    ax[0].axvline(bs["V_min_power"], color="grey", ls="--", lw=1)
    ax[0].set_xlabel("airspeed V (m/s)"); ax[0].set_ylabel("power (kW)"); ax[0].set_ylim(0, None)
    ax[0].set_title("1100 kg light helicopter, sea level"); ax[0].legend(); ax[0].grid(alpha=.3)
    mus = np.linspace(0, 0.35, 141)
    li = [glauert_inflow(m, CT, 0.0)[1] for m in mus]
    ax[1].plot(mus, np.array(li) / np.sqrt(CT / 2), label="Glauert, alpha = 0")
    ax[1].plot(mus[10:], CT / (2 * mus[10:]) / np.sqrt(CT / 2), "--", label="C_T / (2 mu)")
    ax[1].set_xlabel("advance ratio mu"); ax[1].set_ylabel("lambda_i / lambda_h"); ax[1].legend(); ax[1].grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(out_dir / "helicopter_power.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
