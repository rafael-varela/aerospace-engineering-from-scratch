"""P10.17 Propeller & Rotor Theory: reference solution.

What rotational speed does the reference UAV's 10 x 4.5 in propeller need to give the cruise
thrust at 12 m/s, how efficient is it there, and how far below the momentum-theory ideal does it
sit? The module answers with two models, the second checked against the first:

* actuator-disk momentum theory: a thin disk that adds momentum to the stream it captures; it
  gives the induced velocity, the ideal power and the ideal (Froude) efficiency, and it is the
  upper bound no real propeller can beat at the same thrust and advance ratio;
* blade element momentum (BEM) theory: each annulus of the disk obeys momentum theory, while
  the blade element in it produces lift and drag from a section model; the two are made to agree
  by iterating on the axial and swirl induced velocities with relaxation, with Prandtl's
  tip-loss factor F.

Independent check: with zero profile drag, no tip loss and no swirl, a blade designed for uniform
inflow must reproduce the actuator-disk efficiency exactly at the same thrust.

Conventions: V the flight speed (m/s), n the rotational speed (rev/s), Omega = 2 pi n (rad/s),
D the diameter, J = V/(n D), C_T = T/(rho n^2 D^4), C_Q = Q/(rho n^2 D^5), C_P = P/(rho n^3 D^5)
= 2 pi C_Q, eta = J C_T / C_P. Induced velocities: u_a axial (the stream at the disk moves at
V + u_a), u_t swirl (the air at the disk turns with the blade at u_t, so the blade sees
Omega r - u_t). Units SI.

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes propeller_bem.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-10/lesson-17-propeller-rotor-theory
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.optimize import brentq

RHO = 1.225            # kg/m^3, ISA sea level
NU = 1.461e-5          # m^2/s
G = 9.81               # m/s^2

# reference small UAV (shared Phase 10 conventions)
UAV_MASS = 1.2         # kg
UAV_SPAN = 1.5         # m
UAV_CHORD = 0.20       # m
UAV_V = 12.0           # m/s
PROP_D = 0.254         # m (10 in)
PROP_PITCH = 0.114     # m (4.5 in)
PROP_B = 2


# ---------------------------------------------------------------- given helpers

@dataclass
class Section:
    """Stated section model for a low-Reynolds-number propeller section (Re ~ 5e4 - 1.5e5).

    cl rises with slope a0 from the zero-lift angle alpha0 and saturates smoothly at +-cl_max
    (a tanh cap); cd is a drag bucket cd0 + k (cl - cl_md)^2 plus a post-stall term that grows as
    the section passes its stall angle. These are modelling assumptions, not measured data: they
    are typical of thin cambered sections at these Reynolds numbers.
    """
    a0: float = 5.7                     # 1/rad, lift-curve slope (about 0.9 x 2 pi at low Re)
    alpha0: float = np.radians(-4.0)    # rad, zero-lift angle of a cambered section
    cl_max: float = 1.1
    cd0: float = 0.020
    k: float = 0.025
    cl_md: float = 0.3
    stall_drag: float = 1.8

    def cl(self, alpha: np.ndarray) -> np.ndarray:
        """Lift coefficient (exactly linear when cl_max is infinite)."""
        if not np.isfinite(self.cl_max):
            return self.a0 * (np.asarray(alpha) - self.alpha0)
        return self.cl_max * np.tanh(self.a0 * (np.asarray(alpha) - self.alpha0) / self.cl_max)

    def cd(self, alpha: np.ndarray) -> np.ndarray:
        """Drag coefficient."""
        alpha = np.asarray(alpha)
        cl = self.cl(alpha)
        if not np.isfinite(self.cl_max):
            return self.cd0 + self.k * (cl - self.cl_md) ** 2
        a_stall = self.cl_max / self.a0
        excess = np.clip(np.abs(alpha - self.alpha0) - a_stall, 0.0, None)
        return self.cd0 + self.k * (cl - self.cl_md) ** 2 + self.stall_drag * np.sin(excess) ** 2


@dataclass
class Blade:
    """Blade geometry: chord(r) and pitch angle beta(r) as callables of radius r (m)."""
    R: float
    B: int
    r_hub: float
    chord: Callable[[np.ndarray], np.ndarray]
    beta: Callable[[np.ndarray], np.ndarray]
    section: Section = field(default_factory=Section)


def reference_propeller(section: Section | None = None) -> Blade:
    """The reference UAV's 10 x 4.5 two-blade propeller with an assumed, stated geometry.

    Chord (assumed smooth): c/R = 0.045 + 0.55 x (1 - x)^1.3 with x = r/R, which peaks near
    c/R = 0.16 at mid-radius and ends at 0.045 R at the tip. Pitch: constant geometric pitch p =
    0.114 m, so beta(r) = atan(p / (2 pi r)). Hub radius 0.15 R.
    """
    R = PROP_D / 2

    def chord(r: np.ndarray) -> np.ndarray:
        x = np.asarray(r) / R
        return R * (0.045 + 0.55 * x * (1 - x) ** 1.3)

    def beta(r: np.ndarray) -> np.ndarray:
        return np.arctan(PROP_PITCH / (2 * np.pi * np.asarray(r)))

    return Blade(R=R, B=PROP_B, r_hub=0.15 * R, chord=chord, beta=beta, section=section or Section())


def annuli(R: float, r_hub: float, n: int = 60) -> tuple[np.ndarray, np.ndarray]:
    """Element midpoints and widths, clustered towards the tip with half-cosine spacing."""
    t = np.linspace(0.0, np.pi / 2, n + 1)
    edges = r_hub + (R - r_hub) * np.sin(t)
    return 0.5 * (edges[1:] + edges[:-1]), np.diff(edges)


def uav_cruise_thrust(cd0: float = 0.035, e: float = 0.8) -> float:
    """Thrust needed by the reference UAV in level cruise, from an assumed polar
    C_D = cd0 + C_L^2 / (pi e AR) (cd0 and e are stated assumptions)."""
    S = UAV_SPAN * UAV_CHORD
    AR = UAV_SPAN ** 2 / S
    q = 0.5 * RHO * UAV_V ** 2
    CL = UAV_MASS * G / (q * S)
    return q * S * (cd0 + CL ** 2 / (np.pi * e * AR))


# ---------------------------------------------------------------- momentum theory

def actuator_disk(T: float, V: float, A: float, rho: float = RHO) -> dict:
    """Actuator-disk momentum theory for thrust T at flight speed V through disk area A.

    T = 2 rho A (V + v) v gives the induced velocity v at the disk (the far wake moves at V + 2v);
    the ideal power is P = T (V + v) and the ideal (Froude) efficiency is eta = T V / P = V/(V + v).
    Returns dict(v, P, eta, mdot, v_wake).
    """
    v = -V / 2 + np.sqrt(V ** 2 / 4 + T / (2 * rho * A))
    P = T * (V + v)
    eta = V / (V + v) if V + v > 0 else 0.0
    return {"v": v, "P": P, "eta": eta, "mdot": rho * A * (V + v), "v_wake": V + 2 * v}


def ideal_efficiency(CT: float | np.ndarray, J: float | np.ndarray) -> np.ndarray:
    """Froude efficiency in propeller coefficients: eta_i = 2 / (1 + sqrt(1 + 8 C_T / (pi J^2)))."""
    return 2.0 / (1.0 + np.sqrt(1.0 + 8.0 * np.asarray(CT) / (np.pi * np.asarray(J) ** 2)))


def static_ideal_power(T: float, A: float, rho: float = RHO) -> float:
    """Ideal power to hold thrust T at zero speed: P = T^(3/2) / sqrt(2 rho A)."""
    return T ** 1.5 / np.sqrt(2 * rho * A)


# ---------------------------------------------------------------- blade element momentum

def prandtl_tip_loss(r: np.ndarray, R: float, B: int, phi: np.ndarray) -> np.ndarray:
    """Prandtl's tip-loss factor F = (2/pi) acos(exp(-B (R - r) / (2 r sin phi)))."""
    s = np.maximum(np.abs(np.sin(phi)), 1e-6)
    f = B * (R - r) / (2 * r * s)
    return (2 / np.pi) * np.arccos(np.clip(np.exp(-f), 0.0, 1.0))


def _loads(blade: Blade, r, c, beta, phi, W2, rho, drag):
    """Section forces per unit span resolved into thrust (normal) and torque (tangential) directions."""
    alpha = beta - phi
    cl = blade.section.cl(alpha)
    cd = blade.section.cd(alpha) if drag else np.zeros_like(r)
    Cn = cl * np.cos(phi) - cd * np.sin(phi)
    Ct = cl * np.sin(phi) + cd * np.cos(phi)
    return alpha, cl, cd, Cn, Ct


def _pack(blade, V, n_rps, rho, r, dr, c, beta, phi, ua, ut, F, drag, resid, iters) -> dict:
    """Integrate the annulus loads and form the coefficients."""
    Om = 2 * np.pi * n_rps
    W2 = (V + ua) ** 2 + (Om * r - ut) ** 2
    alpha, cl, cd, Cn, Ct = _loads(blade, r, c, beta, phi, W2, rho, drag)
    qW = 0.5 * rho * W2
    dT = blade.B * c * qW * Cn
    dQ = blade.B * c * qW * Ct * r
    T = float(np.sum(dT * dr))
    Q = float(np.sum(dQ * dr))
    P = Q * Om
    D = 2 * blade.R
    CT = T / (rho * n_rps ** 2 * D ** 4)
    CQ = Q / (rho * n_rps ** 2 * D ** 5)
    CP = P / (rho * n_rps ** 3 * D ** 5)
    J = V / (n_rps * D)
    eta = J * CT / CP if CP > 0 and T > 0 else 0.0
    return {"T": T, "Q": Q, "P": P, "CT": CT, "CQ": CQ, "CP": CP, "J": J, "eta": eta,
            "r": r, "dr": dr, "ua": ua, "ut": ut, "phi": phi, "alpha": alpha, "cl": cl, "cd": cd,
            "F": F, "dT": dT, "dQ": dQ, "residual": resid, "iterations": iters}


def element_residual(phi: np.ndarray, r: np.ndarray, c: np.ndarray, beta: np.ndarray, blade: Blade,
                     V: float, Omega: float, tip_loss: bool = True, swirl: bool = True,
                     drag: bool = True) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """BEM residual of each element as a function of its inflow angle phi alone.

    Equating blade-element and momentum thrust and torque on an annulus gives
      k  = sigma' Cn / (4 F sin^2 phi)        (u_a = V k / (1 - k)),
      k' = sigma' Ct / (4 F sin phi cos phi)  (u_t = Omega r k' / (1 + k')),
    with local solidity sigma' = B c / (2 pi r). The velocity triangle
    (V + u_a) cos phi = (Omega r - u_t) sin phi then becomes
      R(phi) = sin phi (1 - k) - (V / (Omega r)) cos phi (1 + k') = 0,
    which has no division by V, so the static case (V = 0, root at k = 1) needs no special code.
    Returns (R, k, k', F).
    """
    alpha, cl, cd, Cn, Ct = _loads(blade, r, c, beta, phi, None, None, drag)
    F = prandtl_tip_loss(r, blade.R, blade.B, phi) if tip_loss else np.ones_like(r)
    F = np.maximum(F, 1e-12)
    sig = blade.B * c / (2 * np.pi * r)
    sp, cp = np.sin(phi), np.cos(phi)
    k = sig * Cn / (4 * F * sp ** 2)
    kp = sig * Ct / (4 * F * sp * cp) if swirl else np.zeros_like(r)
    R = sp * (1 - k) - V / (Omega * r) * cp * (1 + kp)
    return R, k, kp, F


def bem(blade: Blade, V: float, n_rps: float, rho: float = RHO, n_el: int = 60,
        tip_loss: bool = True, swirl: bool = True, drag: bool = True, tol: float = 1e-12) -> dict:
    """Blade element momentum solution by bracketing the element residual R(phi) (robust).

    Every element's inflow angle is found by bisection of ``element_residual`` on
    (1e-6, pi/2], all elements at once (vectorised); R < 0 at small phi and R > 0 at phi = pi/2
    for a thrusting propeller element. Then u_t = Omega r k'/(1 + k') and the velocity triangle
    gives u_a = (Omega r - u_t) tan phi - V. Switches: tip_loss (F = 1 if off), swirl (u_t = 0 if
    off), drag (cd = 0 if off). Returns T, Q, P, CT, CQ, CP, eta, J, radial distributions and the
    largest |R| at the solution.
    """
    Om = 2 * np.pi * n_rps
    r, dr = annuli(blade.R, blade.r_hub, n_el)
    c = blade.chord(r)
    beta = blade.beta(r)
    lo = np.full_like(r, 1e-6)
    hi = np.full_like(r, np.pi / 2)
    R_lo = element_residual(lo, r, c, beta, blade, V, Om, tip_loss, swirl, drag)[0]
    it = 0
    while np.max(hi - lo) > tol and it < 200:
        mid = 0.5 * (lo + hi)
        R_mid = element_residual(mid, r, c, beta, blade, V, Om, tip_loss, swirl, drag)[0]
        same = np.sign(R_mid) == np.sign(R_lo)
        lo = np.where(same, mid, lo)
        R_lo = np.where(same, R_mid, R_lo)
        hi = np.where(same, hi, mid)
        it += 1
    phi = 0.5 * (lo + hi)
    R, k, kp, F = element_residual(phi, r, c, beta, blade, V, Om, tip_loss, swirl, drag)
    ut = Om * r * kp / (1 + kp)
    ua = (Om * r - ut) * np.tan(phi) - V
    return _pack(blade, V, n_rps, rho, r, dr, c, beta, phi, ua, ut, F, drag, float(np.max(np.abs(R))), it)


def bem_fixed_point(blade: Blade, V: float, n_rps: float, rho: float = RHO, n_el: int = 60,
                    tip_loss: bool = True, swirl: bool = True, drag: bool = True,
                    relax: float = 0.3, tol: float = 1e-10, max_iter: int = 2000) -> dict:
    """The same BEM equations solved by relaxed fixed-point iteration on (u_a, u_t).

    Each sweep: phi = atan2(V + u_a, Omega r - u_t); blade-element dT/dr and dQ/dr; then the
    momentum equations dT/dr = 4 pi r rho F (V + u_a) u_a and dQ/dr = 4 pi r^2 rho F (V + u_a) u_t
    solved for new (u_a, u_t); u <- u + relax (u_new - u). The classical method. It converges for
    moderate resolutions; elements very close to the tip (F -> 0) make it stiff, and it can then
    oscillate or diverge (the lesson's convergence study shows where). 'residual' is the final
    relative change of the induced velocities.
    """
    Om = 2 * np.pi * n_rps
    r, dr = annuli(blade.R, blade.r_hub, n_el)
    c = blade.chord(r)
    beta = blade.beta(r)
    ua = np.full_like(r, 0.1)
    ut = np.zeros_like(r)
    resid = np.inf
    it = 0
    with np.errstate(all="ignore"):
        for it in range(1, max_iter + 1):
            phi = np.arctan2(V + ua, Om * r - ut)
            W2 = (V + ua) ** 2 + (Om * r - ut) ** 2
            _, _, _, Cn, Ct = _loads(blade, r, c, beta, phi, W2, rho, drag)
            qW = 0.5 * rho * W2
            dT = blade.B * c * qW * Cn
            dQ = blade.B * c * qW * Ct * r
            F = prandtl_tip_loss(r, blade.R, blade.B, phi) if tip_loss else np.ones_like(r)
            F = np.maximum(F, 1e-12)
            g = dT / (4 * np.pi * r * rho * F)
            ua_new = -V / 2 + np.sqrt(np.maximum(V ** 2 / 4 + g, 0.0))
            ut_new = dQ / (4 * np.pi * r ** 2 * rho * F * np.maximum(V + ua_new, 1e-9)) if swirl else np.zeros_like(r)
            scale = max(np.max(np.abs(ua_new)), np.max(np.abs(ut_new)), 1e-12)
            resid = max(np.max(np.abs(ua_new - ua)), np.max(np.abs(ut_new - ut))) / scale
            ua = ua + relax * (ua_new - ua)
            ut = ut + relax * (ut_new - ut)
            if not np.isfinite(resid) or resid < tol:
                break
    phi = np.arctan2(V + ua, Om * r - ut)
    F = prandtl_tip_loss(r, blade.R, blade.B, phi) if tip_loss else np.ones_like(r)
    return _pack(blade, V, n_rps, rho, r, dr, c, beta, phi, ua, ut, F, drag, float(resid), it)


def design_uniform_inflow_blade(T: float, V: float, n_rps: float, R: float, B: int = 2,
                                r_hub_frac: float = 0.15, cl_design: float = 0.5,
                                rho: float = RHO) -> Blade:
    """A blade whose ideal BEM solution (no drag, no swirl, no tip loss) has uniform inflow.

    Momentum theory gives u_a for thrust T on the annulus area pi (R^2 - r_hub^2); every annulus
    must then carry dT/dr = 4 pi r rho (V + u_a) u_a. At each radius the inflow angle is
    phi = atan((V + u_a)/(Omega r)); the chord follows from the blade-element thrust at cl_design,
    c = dT/dr / (B q_W cl_design cos phi), and the pitch is beta = phi + alpha(cl_design) for a
    linear section cl = 2 pi alpha (no drag).
    """
    Om = 2 * np.pi * n_rps
    r_hub = r_hub_frac * R
    A = np.pi * (R ** 2 - r_hub ** 2)
    ua = actuator_disk(T, V, A, rho)["v"]
    sec = Section(a0=2 * np.pi, alpha0=0.0, cl_max=np.inf, cd0=0.0, k=0.0, cl_md=0.0, stall_drag=0.0)
    a_des = cl_design / (2 * np.pi)

    def phi(r):
        return np.arctan2(V + ua, Om * np.asarray(r))

    def chord(r):
        r = np.asarray(r)
        dTdr = 4 * np.pi * r * rho * (V + ua) * ua
        qW = 0.5 * rho * ((V + ua) ** 2 + (Om * r) ** 2)
        return dTdr / (B * qW * cl_design * np.cos(phi(r)))

    def beta(r):
        return phi(r) + a_des

    return Blade(R=R, B=B, r_hub=r_hub, chord=chord, beta=beta, section=sec)


def rpm_for_thrust(blade: Blade, T_req: float, V: float, n_lo: float = 20.0, n_hi: float = 400.0) -> float:
    """Rotational speed (rev/s) at which the BEM thrust equals T_req at flight speed V (brentq)."""
    return brentq(lambda n: bem(blade, V, n)["T"] - T_req, n_lo, n_hi, xtol=1e-8)


def sweep_advance_ratio(blade: Blade, n_rps: float, J: np.ndarray, **kw) -> dict:
    """CT, CP, CQ, eta over advance ratios J at fixed rotational speed (V = J n D)."""
    D = 2 * blade.R
    out = {k: [] for k in ("J", "CT", "CP", "CQ", "eta")}
    for j in J:
        res = bem(blade, j * n_rps * D, n_rps, **kw)
        for k in out:
            out[k].append(res[k])
    return {k: np.array(v) for k, v in out.items()}


def observed_order(errors: np.ndarray, ratio: float = 2.0) -> np.ndarray:
    """Observed order of convergence log(e_k / e_{k+1}) / log(ratio) for successive refinements."""
    e = np.asarray(errors, dtype=float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save propeller_bem.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    blade = reference_propeller()
    D, R = 2 * blade.R, blade.R
    A = np.pi * R ** 2
    T_req = uav_cruise_thrust()
    S = UAV_SPAN * UAV_CHORD
    CL = UAV_MASS * G / (0.5 * RHO * UAV_V ** 2 * S)
    print(f"UAV cruise: C_L = {CL:.3f}, thrust required = {T_req:.3f} N (assumed CD0 = 0.035, e = 0.8)")

    # 1. momentum theory at the cruise point
    ad = actuator_disk(T_req, UAV_V, A)
    print(f"actuator disk at 12 m/s: v = {ad['v']:.3f} m/s, wake {ad['v_wake']:.3f} m/s, "
          f"ideal power {ad['P']:.2f} W, ideal efficiency {ad['eta']:.4f}")

    # 2. BEM: rpm for the cruise thrust
    n = rpm_for_thrust(blade, T_req, UAV_V)
    res = bem(blade, UAV_V, n)
    eta_i = float(ideal_efficiency(res["CT"], res["J"]))
    print(f"BEM cruise: n = {n:.2f} rev/s = {60 * n:.0f} rpm, J = {res['J']:.4f}, CT = {res['CT']:.5f}, "
          f"CP = {res['CP']:.5f}, power {res['P']:.2f} W, eta = {res['eta']:.4f} "
          f"(momentum ideal at same CT, J: {eta_i:.4f}); residual {res['residual']:.1e} in {res['iterations']} it")
    tipRe = np.hypot(UAV_V, 2 * np.pi * n * 0.75 * R) * float(blade.chord(0.75 * R)) / NU
    print(f"  section Reynolds number at 0.75 R: {tipRe:.3g};  tip Mach {2*np.pi*n*R/340.3:.3f}")

    # 3. loss budget at the cruise rpm and speed: switch the losses on one at a time
    print("loss budget at the cruise point (same rpm, same V):")
    budget = []
    for label, kw in [("ideal: no drag, no swirl, no tip loss", dict(drag=False, swirl=False, tip_loss=False)),
                      ("+ swirl", dict(drag=False, swirl=True, tip_loss=False)),
                      ("+ tip loss", dict(drag=False, swirl=True, tip_loss=True)),
                      ("+ profile drag (full BEM)", dict(drag=True, swirl=True, tip_loss=True))]:
        rr = bem(blade, UAV_V, n, **kw)
        ei = float(ideal_efficiency(rr["CT"], rr["J"]))
        budget.append((label, rr["T"], rr["P"], rr["eta"], ei))
        print(f"  {label:40s} T = {rr['T']:.3f} N  P = {rr['P']:.2f} W  eta = {rr['eta']:.4f}  "
              f"(momentum ideal at this CT: {ei:.4f})")

    # 4. independent check: a uniform-inflow blade with no losses reproduces momentum theory
    nd = 100.0
    ub = design_uniform_inflow_blade(T_req, UAV_V, nd, R)
    ri = bem(ub, UAV_V, nd, drag=False, swirl=False, tip_loss=False)
    A_ann = np.pi * (R ** 2 - ub.r_hub ** 2)
    adu = actuator_disk(T_req, UAV_V, A_ann)
    print(f"uniform-inflow blade (n = {nd} rev/s), ideal BEM: T = {ri['T']:.5f} N (target {T_req:.5f}), "
          f"u_a spread {np.ptp(ri['ua']):.2e} m/s, eta = {ri['eta']:.6f}; actuator disk on the annulus: "
          f"eta = {adu['eta']:.6f}, v = {adu['v']:.5f} m/s; BEM u_a = {np.mean(ri['ua']):.5f} m/s")

    # 5. static thrust
    for rpm in (6000, 8000):
        st = bem(blade, 0.0, rpm / 60)
        FM = static_ideal_power(st["T"], A) / st["P"]
        print(f"static at {rpm} rpm: T = {st['T']:.3f} N, P = {st['P']:.2f} W, CT = {st['CT']:.4f}, "
              f"CP = {st['CP']:.4f}, figure of merit {FM:.3f}, ideal power {static_ideal_power(st['T'], A):.2f} W")

    # 6. sweep in advance ratio
    J = np.linspace(0.05, 0.60, 23)
    sw = sweep_advance_ratio(blade, n, J)
    i = int(np.argmax(sw["eta"]))
    j0 = brentq(lambda j: bem(blade, j * n * D, n)["T"], 0.3, 0.9)
    print(f"sweep at {60*n:.0f} rpm: eta_max = {sw['eta'][i]:.4f} at J = {sw['J'][i]:.3f}; zero thrust at J = {j0:.4f} "
          f"(geometric pitch/D = {PROP_PITCH / D:.4f})")
    for jj in (0.1, 0.2, 0.3, 0.4, 0.5):
        rr = bem(blade, jj * n * D, n)
        print(f"  J = {jj:.1f}: CT = {rr['CT']:.5f}, CP = {rr['CP']:.5f}, eta = {rr['eta']:.4f}, "
              f"eta_ideal = {float(ideal_efficiency(rr['CT'], jj)):.4f}")

    # 7. convergence in the number of elements (reference: n_el = 1280)
    print("element convergence of cruise thrust (error against n_el = 1280):")
    T_ref = bem(blade, UAV_V, n, n_el=1280)["T"]
    nes = (10, 20, 40, 80, 160)
    errs = np.array([abs(bem(blade, UAV_V, n, n_el=ne)["T"] - T_ref) for ne in nes])
    for ne, e in zip(nes, errs):
        print(f"  n_el = {ne:4d}: |T - T_ref| = {e:.3e} N")
    print(f"  observed order {np.array2string(observed_order(errs), precision=2)}")

    # 8. the classical relaxed fixed-point iteration against the bracketing solver
    print("relaxed fixed-point iteration (cruise point):")
    for ne in (20, 60, 160, 320):
        ref = bem(blade, UAV_V, n, n_el=ne)["T"]
        for rl in (0.3, 0.1):
            fp = bem_fixed_point(blade, UAV_V, n, n_el=ne, relax=rl)
            ok = np.isfinite(fp["residual"]) and fp["residual"] < 1e-10
            dT = abs(fp["T"] - ref) if ok else float("nan")
            print(f"  n_el = {ne:3d}, relax = {rl}: {'converged' if ok else 'NOT converged'} after {fp['iterations']} "
                  f"iterations, |T_fp - T_bracket| = {dT:.1e} N")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    ax[0].plot(sw["J"], sw["CT"], label="C_T")
    ax[0].plot(sw["J"], sw["CP"], label="C_P")
    ax[0].set_xlabel("advance ratio J"); ax[0].set_title("10 x 4.5 propeller, BEM"); ax[0].legend(); ax[0].grid(alpha=.3)
    ax[1].plot(sw["J"], sw["eta"], label="BEM")
    mask = sw["CT"] > 0
    ax[1].plot(sw["J"][mask], ideal_efficiency(sw["CT"][mask], sw["J"][mask]), "--", label="momentum ideal, same C_T")
    ax[1].set_ylim(0, 1); ax[1].set_xlabel("J"); ax[1].set_ylabel("efficiency"); ax[1].legend(); ax[1].grid(alpha=.3)
    ax[2].plot(res["r"] / R, res["dT"], label="dT/dr (N/m)")
    ax[2].plot(res["r"] / R, res["ua"] * 10, label="u_a x 10 (m/s)")
    ax[2].plot(res["r"] / R, res["F"] * 10, ":", label="F x 10")
    ax[2].set_xlabel("r/R"); ax[2].set_title("cruise radial distributions"); ax[2].legend(); ax[2].grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(out_dir / "propeller_bem.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
