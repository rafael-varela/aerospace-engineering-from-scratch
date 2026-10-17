"""P10.16 Fuselage & Body Aerodynamics: reference solution.

How much does the reference trainer's 7 m fuselage destabilise it in pitch, and how far forward
does it move the neutral point? This module answers with slender-body theory, checked against the
exact potential flow about a prolate spheroid:

* slender-body (Munk) theory: in a cross-flow V sin(alpha) each section of a slender body carries
  apparent mass rho S(x) per unit length; its rate of change along the body gives the local normal
  force dN/dx = 2 q alpha dS/dx. Integrating gives N = 2 q alpha S_base (C_N_alpha = 2 on the base
  area; zero for a closed body) and, for a closed body, a pure nose-up couple M = 2 q alpha Vol
  (Munk's moment, dM/dalpha = rho V^2 Vol), whatever the moment reference;
* the exact potential solution for a prolate spheroid (Lamb's added-mass coefficients k1, k2):
  M = q Vol (k2 - k1) sin(2 alpha), which tends to the slender-body moment as the fineness ratio grows;
* the Allen-Perkins viscous cross-flow correction (empirical) for the normal force at larger alpha;
* body drag from a turbulent flat-plate skin friction, the wetted area and an empirical form
  factor of the fineness ratio; the fineness ratio of least drag for a given volume;
* the fuselage's Cm_alpha on the trainer and the neutral-point shift.

Coordinates: x from the nose, aft positive; normal force N positive up; pitching moment positive
nose-up about x_ref: an upward force at x contributes N (x_ref - x). Units SI, alpha in radians.

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes fuselage_slender_body.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-10/lesson-16-fuselage-body-aerodynamics
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

RHO = 1.225
NU = 1.461e-5

# reference trainer (shared Phase 10 conventions)
S_REF = 12.0         # m^2 wing area
B_SPAN = 10.0        # m
C_REF = 1.2          # m chord
V_CRUISE = 55.0      # m/s
L_FUS = 7.0          # m fuselage length
D_FUS = 1.1          # m maximum diameter
S_TAIL = 2.0         # m^2 horizontal tail area
L_TAIL = 4.5         # m tail arm behind the wing quarter chord
X_WING_LE = 2.0      # m from the nose (assumed layout)


# ---------------------------------------------------------------- given helpers: geometry

def ellipsoid_radius(x: np.ndarray, L: float, R: float) -> np.ndarray:
    """Prolate spheroid (ellipsoid of revolution) of length L and maximum radius R."""
    xi = 2 * np.asarray(x) / L - 1
    return R * np.sqrt(np.clip(1 - xi ** 2, 0.0, None))


def paraboloid_radius(x: np.ndarray, L: float, R: float) -> np.ndarray:
    """Paraboloid of revolution r = R sqrt(x / L), 0 <= x <= L: a nose ending in a full base of radius R."""
    return R * np.sqrt(np.clip(np.asarray(x) / L, 0.0, None))


def trainer_fuselage_radius(x: np.ndarray) -> np.ndarray:
    """Equivalent-radius profile of the reference trainer's fuselage (assumed shape, stated).

    Nose 0-1.8 m: elliptic, from 0 to the maximum radius 0.55 m; cabin 1.8-3.2 m: constant; tail cone
    3.2-7.0 m: r = 0.55 - 0.47 ((x - 3.2)/3.8)^1.3, ending in a small base of radius 0.08 m. The real
    fuselage is not round; its cross-section area is represented by a circle of the same area.
    """
    x = np.asarray(x, dtype=float)
    Rm = D_FUS / 2
    r = np.empty_like(x)
    nose = x < 1.8
    cabin = (x >= 1.8) & (x <= 3.2)
    tail = x > 3.2
    r[nose] = Rm * np.sqrt(np.clip(1 - (1 - x[nose] / 1.8) ** 2, 0.0, None))
    r[cabin] = Rm
    r[tail] = Rm - (Rm - 0.08) * ((x[tail] - 3.2) / 3.8) ** 1.3
    return r


def turbulent_cf(Re: float) -> float:
    """Prandtl-Schlichting turbulent flat-plate mean skin friction 0.455 / (log10 Re)^2.58."""
    return 0.455 / np.log10(Re) ** 2.58


def helmbold_lift_slope(AR: float, a0: float = 2 * np.pi) -> float:
    """Lift-curve slope of an unswept wing (Helmbold's formula, a lifting-line-based approximation)."""
    k = a0 / (np.pi * AR)
    return a0 / (np.sqrt(1 + k ** 2) + k)


# ---------------------------------------------------------------- slender-body theory

def section_area(r: np.ndarray) -> np.ndarray:
    """Cross-section area S = pi r^2."""
    return np.pi * np.asarray(r) ** 2


def body_volume(x: np.ndarray, r: np.ndarray) -> float:
    """Volume = integral of S dx (trapezoid rule)."""
    return float(np.trapezoid(section_area(r), x))


def wetted_area(x: np.ndarray, r: np.ndarray) -> float:
    """Wetted area = integral of 2 pi r sqrt(1 + r'^2) dx (trapezoid rule on the polyline: exact for the
    chain of frustums through the points), plus a flat base disk is NOT included."""
    dx = np.diff(x)
    dr = np.diff(r)
    return float(np.sum(np.pi * (r[1:] + r[:-1]) * np.sqrt(dx ** 2 + dr ** 2)))


def slender_body_loading(x: np.ndarray, r: np.ndarray, alpha: float, q: float) -> np.ndarray:
    """Local normal force per unit length dN/dx = 2 q alpha dS/dx (N/m), dS/dx by central differences.

    Derivation: a section of area S in the cross-flow w = V alpha carries vertical momentum
    rho S w per unit length (its apparent mass is rho S). Air passing the station at speed V carries it
    aft, so the body must supply d/dx(rho S w) V = rho V^2 alpha dS/dx = 2 q alpha dS/dx.
    """
    S = section_area(r)
    return 2 * q * alpha * np.gradient(S, x)


def normal_force_and_moment(x: np.ndarray, r: np.ndarray, alpha: float, q: float, x_ref: float) -> tuple[float, float]:
    """Integrate the slender-body loading: N = int dN/dx dx, M = int dN/dx (x_ref - x) dx (nose-up +)."""
    dN = slender_body_loading(x, r, alpha, q)
    N = float(np.trapezoid(dN, x))
    M = float(np.trapezoid(dN * (x_ref - x), x))
    return N, M


def slender_body_closed_form(S_base: float, volume: float, L: float, alpha: float, q: float,
                             x_ref: float) -> tuple[float, float]:
    """Closed forms from integrating 2 q alpha S'(x) by parts over 0 <= x <= L:
    N = 2 q alpha S_base, M = 2 q alpha (Vol - S_base (L - x_ref))."""
    N = 2 * q * alpha * S_base
    M = 2 * q * alpha * (volume - S_base * (L - x_ref))
    return N, M


def lamb_coefficients(fineness: float) -> tuple[float, float]:
    """Lamb's added-mass coefficients (k1 axial, k2 transverse) of a prolate spheroid of fineness
    ratio L/D, from the exact potential solution:
      e = sqrt(1 - (D/L)^2), a0 = (2 (1 - e^2) / e^3)(artanh e - e),
      b0 = 1/e^2 - ((1 - e^2) / e^3) artanh e,  k1 = a0 / (2 - a0),  k2 = b0 / (2 - b0).
    A sphere (L/D = 1) has k1 = k2 = 1/2; L/D -> infinity gives k1 -> 0, k2 -> 1."""
    if abs(fineness - 1.0) < 1e-9:
        return 0.5, 0.5
    e = np.sqrt(1 - 1 / fineness ** 2)
    at = np.arctanh(e)
    a0 = 2 * (1 - e ** 2) / e ** 3 * (at - e)
    b0 = 1 / e ** 2 - (1 - e ** 2) / e ** 3 * at
    return float(a0 / (2 - a0)), float(b0 / (2 - b0))


def spheroid_munk_moment(fineness: float, volume: float, alpha: float, q: float) -> float:
    """Exact potential-flow moment on a prolate spheroid: M = q Vol (k2 - k1) sin(2 alpha)."""
    k1, k2 = lamb_coefficients(fineness)
    return q * volume * (k2 - k1) * np.sin(2 * alpha)


def allen_perkins_cn(alpha: float, S_base: float, A_plan: float, S_ref: float,
                     eta: float = 0.65, cdc: float = 1.2) -> float:
    """Allen-Perkins normal-force coefficient (empirical): potential part plus viscous cross-flow,
    C_N = (S_base / S_ref) sin(2 alpha) cos(alpha / 2) + eta c_dc (A_plan / S_ref) sin^2(alpha),
    A_plan the planform (side-view) area, c_dc ~ 1.2 the drag coefficient of a circular cylinder in
    cross-flow and eta ~ 0.6-0.7 a finite-length factor (both empirical, stated)."""
    return (S_base / S_ref) * np.sin(2 * alpha) * np.cos(alpha / 2) + eta * cdc * (A_plan / S_ref) * np.sin(alpha) ** 2


# ---------------------------------------------------------------- drag and fineness

def body_form_factor(fineness: float) -> float:
    """Empirical body form factor FF = 1 + 60 / f^3 + f / 400 (a design-textbook correlation of body
    drag data: the 60/f^3 term is the pressure drag of the bluntness, f/400 the extra wetted length)."""
    return 1 + 60 / fineness ** 3 + fineness / 400


def spheroid_drag_area(fineness: float, volume: float, V: float, nu: float = NU, n: int = 801) -> dict:
    """Drag area D/q = C_f FF S_wet of a spheroid of given volume and fineness, fully turbulent."""
    # Vol = (pi/6) L D^2 and L = f D  ->  D = (6 Vol / (pi f))^(1/3)
    D = (6 * volume / (np.pi * fineness)) ** (1 / 3)
    L = fineness * D
    t = np.linspace(0, np.pi, n)
    x = L * (1 - np.cos(t)) / 2
    r = ellipsoid_radius(x, L, D / 2)
    Sw = wetted_area(x, r)
    Cf = turbulent_cf(V * L / nu)
    FF = body_form_factor(fineness)
    return {"L": L, "D": D, "S_wet": Sw, "Cf": Cf, "FF": FF, "D_over_q": Cf * FF * Sw}


def best_fineness(volume: float, V: float) -> float:
    """Fineness ratio that minimises the drag area of a spheroid of fixed volume."""
    res = minimize_scalar(lambda f: spheroid_drag_area(f, volume, V)["D_over_q"], bounds=(1.5, 20.0),
                          method="bounded", options={"xatol": 1e-6})
    return float(res.x)


# ---------------------------------------------------------------- the trainer

def fuselage_pitch_derivatives(x: np.ndarray, r: np.ndarray, x_ref: float, viscous_cut: bool = False) -> dict:
    """Fuselage C_N_alpha and C_m_alpha (per rad, on the wing area S and chord c) about x_ref.

    viscous_cut=False: the full slender-body (potential) loading. viscous_cut=True: the loading is kept
    only ahead of the station of maximum area and set to zero behind it, a common engineering
    approximation for the separated, thick boundary layer of a tail cone, which cannot carry the
    potential-flow negative lift."""
    q, alpha = 1.0, 1.0       # linear: per unit q and per radian
    dN = slender_body_loading(x, r, alpha, q)
    if viscous_cut:
        dN = np.where(x <= x[np.argmax(r)], dN, 0.0)
    N = float(np.trapezoid(dN, x))
    M = float(np.trapezoid(dN * (x_ref - x), x))
    return {"CN_alpha": N / S_REF, "Cm_alpha": M / (S_REF * C_REF)}


def trainer_neutral_point(Cm_alpha_fus: float = 0.0, CN_alpha_fus: float = 0.0, x_fus_ref: float | None = None,
                          AR_tail: float = 4.0, eta_t: float = 0.9) -> dict:
    """Stick-fixed neutral point of the trainer, h_n in chords aft of the wing leading edge.

    Wing: Helmbold lift slope at AR = b^2/S, aerodynamic centre 0.25 c. Tail (assumed AR 4, efficiency
    eta_t): Helmbold slope, aerodynamic centre L_TAIL behind the wing quarter chord, downwash gradient
    d eps / d alpha = 2 a_w / (pi AR) (elliptic-wing estimate). Fuselage: a moment derivative about the
    wing quarter chord (x_fus_ref) plus, optionally, a normal force acting there. h_n is where
    dCm/dalpha = 0: h_n = (sum a_i h_i - Cm_alpha_fus) / sum a_i, with a_i the lift slopes on S.
    """
    AR = B_SPAN ** 2 / S_REF
    a_w = helmbold_lift_slope(AR)
    a_t = helmbold_lift_slope(AR_tail)
    deps = 2 * a_w / (np.pi * AR)
    a_tail_eff = eta_t * S_TAIL / S_REF * a_t * (1 - deps)
    h_w = 0.25
    h_t = 0.25 + L_TAIL / C_REF
    num = a_w * h_w + a_tail_eff * h_t + CN_alpha_fus * h_w - Cm_alpha_fus
    den = a_w + a_tail_eff + CN_alpha_fus
    return {"h_n": num / den, "a_w": a_w, "a_t": a_t, "deps_dalpha": deps, "CL_alpha": den,
            "V_H": S_TAIL * L_TAIL / (S_REF * C_REF)}


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save fuselage_slender_body.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    q = 0.5 * RHO * V_CRUISE ** 2
    alpha = np.radians(4.0)
    print(f"cruise q = {q:.1f} Pa, alpha = 4 deg")

    # 1. closed bodies and the paraboloid: numerical loading vs closed forms
    L, R = 7.0, 0.55
    t = np.linspace(0, np.pi, 4001)
    x = L * (1 - np.cos(t)) / 2
    r = ellipsoid_radius(x, L, R)
    vol_e = 4 / 3 * np.pi * (L / 2) * R ** 2
    for xr in (0.0, 3.5, 7.0):
        N, M = normal_force_and_moment(x, r, alpha, q, xr)
        print(f"spheroid L = {L}, D = {2 * R}: x_ref = {xr}: N = {N:.3e} N, M = {M:.2f} N m; closed form 2 q alpha Vol "
              f"= {2 * q * alpha * vol_e:.2f} N m")
    xp = np.linspace(0, 2.0, 4001)
    rp = paraboloid_radius(xp, 2.0, R)
    Np, Mp = normal_force_and_moment(xp, rp, alpha, q, 0.0)
    Nc, Mc = slender_body_closed_form(np.pi * R ** 2, body_volume(xp, rp), 2.0, alpha, q, 0.0)
    print(f"paraboloid nose 2.0 m: N = {Np:.2f} N (closed {Nc:.2f}), centre of pressure {-Mp / Np:.4f} m "
          f"(closed {-Mc / Nc:.4f} m = L/2); C_N_alpha on base area = {Np / (q * alpha * np.pi * R ** 2):.4f}")

    # 2. exact spheroid vs slender body
    print("prolate spheroid: exact potential moment / slender-body moment = k2 - k1:")
    for f in (1.0, 2.0, 4.0, 6.36, 10.0, 20.0, 50.0):
        k1, k2 = lamb_coefficients(f)
        print(f"  L/D = {f:6.2f}: k1 = {k1:.4f}, k2 = {k2:.4f}, k2 - k1 = {k2 - k1:.4f}")

    # 3. the trainer fuselage
    xf = np.linspace(0, L_FUS, 7001)
    rf = trainer_fuselage_radius(xf)
    vol = body_volume(xf, rf)
    Sw = wetted_area(xf, rf)
    fin = L_FUS / D_FUS
    print(f"trainer fuselage: volume {vol:.3f} m^3, wetted area {Sw:.2f} m^2, fineness {fin:.2f}, base area "
          f"{np.pi * 0.08 ** 2:.4f} m^2")
    x_qc = X_WING_LE + 0.25 * C_REF
    k1, k2 = lamb_coefficients(fin)
    full = fuselage_pitch_derivatives(xf, rf, x_qc)
    cut = fuselage_pitch_derivatives(xf, rf, x_qc, viscous_cut=True)
    print(f"slender body, about the wing quarter chord (x = {x_qc} m): C_N_alpha = {full['CN_alpha']:.4f}, "
          f"C_m_alpha = {full['Cm_alpha']:.4f} /rad; 2 Vol/(S c) = {2 * vol / (S_REF * C_REF):.4f}")
    print(f"  with the spheroid correction k2 - k1 = {k2 - k1:.4f}: C_m_alpha = {(k2 - k1) * full['Cm_alpha']:.4f} /rad")
    print(f"  loading cut behind the maximum section: C_N_alpha = {cut['CN_alpha']:.4f}, C_m_alpha = {cut['Cm_alpha']:.4f} /rad")
    np0 = trainer_neutral_point()
    np1 = trainer_neutral_point((k2 - k1) * full["Cm_alpha"], (k2 - k1) * full["CN_alpha"])
    np2 = trainer_neutral_point(cut["Cm_alpha"], cut["CN_alpha"])
    print(f"wing a_w = {np0['a_w']:.3f} /rad, tail a_t = {np0['a_t']:.3f} /rad, d eps/d alpha = {np0['deps_dalpha']:.3f}, "
          f"V_H = {np0['V_H']:.3f}")
    print(f"neutral point (chords aft of wing LE): wing + tail {np0['h_n']:.4f}; + fuselage (spheroid-corrected "
          f"slender body) {np1['h_n']:.4f} (shift {np1['h_n'] - np0['h_n']:+.4f} c = {(np1['h_n'] - np0['h_n']) * C_REF * 100:+.1f} cm); "
          f"+ fuselage (cut loading) {np2['h_n']:.4f} (shift {np2['h_n'] - np0['h_n']:+.4f} c)")

    # 4. drag and fineness
    Re = V_CRUISE * L_FUS / NU
    Cf = turbulent_cf(Re)
    FF = body_form_factor(fin)
    print(f"fuselage drag: Re = {Re:.3g}, Cf = {Cf:.5f}, FF = {FF:.4f}, D/q = {Cf * FF * Sw:.4f} m^2, "
          f"Delta C_D = {Cf * FF * Sw / S_REF:.5f}")
    fb = best_fineness(vol, V_CRUISE)
    print(f"spheroid of the same volume: least drag at L/D = {fb:.2f}")
    for f in (2.0, 3.0, fb, 6.36, 10.0):
        d = spheroid_drag_area(f, vol, V_CRUISE)
        print(f"  L/D = {f:5.2f}: L = {d['L']:.2f} m, D = {d['D']:.3f} m, S_wet = {d['S_wet']:.2f} m^2, FF = {d['FF']:.3f}, "
              f"Cf = {d['Cf']:.5f}, D/q = {d['D_over_q']:.4f} m^2")

    # 5. Allen-Perkins
    A_plan = float(np.trapezoid(2 * rf, xf))
    S_b = np.pi * R ** 2
    print(f"Allen-Perkins on the spheroid-like body (S_base = max section {S_b:.3f} m^2, planform {A_plan:.3f} m^2, on S_max):")
    for a in (2, 5, 10, 20, 30):
        ar = np.radians(a)
        pot = (S_b / S_b) * np.sin(2 * ar) * np.cos(ar / 2)
        print(f"  alpha = {a:2d} deg: C_N potential {pot:.4f}, total {allen_perkins_cn(ar, S_b, A_plan, S_b):.4f}")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    ax[0].plot(xf, rf, "k"); ax[0].plot(xf, -rf, "k"); ax[0].set_aspect("equal")
    ax[0].set_title("trainer fuselage (equivalent radius)"); ax[0].set_xlabel("x (m)")
    dN = slender_body_loading(xf, rf, alpha, q)
    ax[1].plot(xf, dN, label="slender body dN/dx")
    ax[1].axhline(0, color="grey", lw=.5); ax[1].axvline(x_qc, color="grey", ls="--", lw=1, label="wing c/4")
    ax[1].set_xlabel("x (m)"); ax[1].set_ylabel("N/m at 4 deg, 55 m/s"); ax[1].legend(); ax[1].grid(alpha=.3)
    fs = np.linspace(1.0, 30.0, 100)
    ax[2].plot(fs, [np.subtract(*lamb_coefficients(f)[::-1]) for f in fs], label="exact k2 - k1")
    ax[2].axhline(1.0, color="k", ls="--", label="slender body")
    ax[2].set_xlabel("fineness L/D"); ax[2].set_ylabel("moment / slender-body moment"); ax[2].legend(); ax[2].grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(out_dir / "fuselage_slender_body.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
