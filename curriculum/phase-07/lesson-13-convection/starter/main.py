"""P07.13 Heat Transfer: Convection: starter code.

How hot does an avionics box get when it is cooled only by the still air around it, how much air
speed or coolant flow does it take to hold it at a given temperature, and how far can the
correlations that answer these questions be trusted? This module answers with the standard
convection correlations, each written as a Nusselt number Nu = h L / k in terms of the Reynolds,
Prandtl and Rayleigh numbers:

* air and water properties interpolated linearly in temperature from Lienhard & Lienhard,
  A Heat Transfer Textbook, 6th ed. (v6.00, 2024), Tables A.6 (air at 1 atm) and A.3 (saturated
  liquid water); every node of the copied table is reproduced exactly;
* a flat plate in laminar and turbulent flow (local and average Nu, skin friction), with the
  mixed laminar-turbulent average for a sharp transition, checked against the Reynolds-Colburn
  analogy St Pr^(2/3) = c_f/2 and against numerical integration of the local coefficient;
* a cylinder in cross-flow (Churchill-Bernstein); pipe flow (fully developed laminar,
  Gnielinski with Filonenko's friction factor, Dittus-Boelter for comparison);
* natural convection on a vertical plate and a horizontal cylinder (Churchill-Chu) and on the
  upper and lower faces of a horizontal plate (Raithby-Hollands, Fujii-Imura);
* three applications: an avionics box in still air (natural convection plus a linearised
  radiation coefficient), the same box in a forced air stream, and a water-cooled cold plate.

Correlation numbers in comments are Lienhard's equation numbers. Properties are evaluated at the
film temperature (Tw + Tinf)/2 unless a function says otherwise. Units: SI (K, W, m, kg, s).

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes convection.png)
Test: pytest curriculum/phase-07/lesson-13-convection
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

G0 = 9.80665                 # m/s^2, standard gravity
SIGMA = 5.670374419e-8       # W/(m^2 K^4), Stefan-Boltzmann constant (CODATA 2018)

# ---------------------------------------------------------------- property data (given)
# Air at 1 atm (101 325 Pa), Lienhard & Lienhard, A Heat Transfer Textbook 6th ed., Table A.6.
# Columns: T [K], rho [kg/m^3], cp [J/(kg K)], mu [1e-5 kg/(m s)], nu [1e-5 m^2/s], k [W/(m K)],
#          alpha [1e-5 m^2/s], Pr [-]
AIR_TABLE_A6 = np.array([
    [200.0, 1.769, 1007.0, 1.333, 0.7537, 0.0185, 1.04, 0.726],
    [250.0, 1.413, 1006.0, 1.604, 1.135, 0.0226, 1.59, 0.715],
    [260.0, 1.359, 1006.0, 1.655, 1.218, 0.0233, 1.71, 0.713],
    [270.0, 1.308, 1006.0, 1.706, 1.304, 0.0241, 1.83, 0.711],
    [280.0, 1.261, 1006.0, 1.756, 1.392, 0.0249, 1.96, 0.710],
    [290.0, 1.218, 1006.0, 1.805, 1.482, 0.0256, 2.09, 0.708],
    [300.0, 1.177, 1006.0, 1.854, 1.575, 0.0264, 2.23, 0.707],
    [310.0, 1.139, 1007.0, 1.902, 1.670, 0.0271, 2.37, 0.706],
    [320.0, 1.103, 1007.0, 1.949, 1.766, 0.0279, 2.51, 0.705],
    [330.0, 1.070, 1008.0, 1.995, 1.865, 0.0286, 2.65, 0.704],
    [340.0, 1.038, 1009.0, 2.041, 1.966, 0.0293, 2.80, 0.703],
    [350.0, 1.009, 1009.0, 2.087, 2.069, 0.0300, 2.95, 0.702],
    [400.0, 0.8823, 1014.0, 2.306, 2.613, 0.0335, 3.74, 0.699],
    [450.0, 0.7842, 1021.0, 2.512, 3.204, 0.0368, 4.59, 0.698],
    [500.0, 0.7058, 1030.0, 2.709, 3.839, 0.0399, 5.50, 0.698],
    [550.0, 0.6416, 1040.0, 2.897, 4.515, 0.0430, 6.45, 0.700],
    [600.0, 0.5881, 1051.0, 3.077, 5.232, 0.0460, 7.44, 0.703],
])
# Saturated liquid water, Lienhard & Lienhard, Table A.3.
# Columns: T [K], rho [kg/m^3], cp [J/(kg K)], k [W/(m K)], alpha [1e-7 m^2/s], nu [1e-7 m^2/s],
#          Pr [-], beta [1/K]
WATER_TABLE_A3 = np.array([
    [275.0, 999.9, 4214.0, 0.5645, 1.340, 16.82, 12.55, -3.55e-5],
    [280.0, 999.9, 4201.0, 0.5740, 1.366, 14.34, 10.63, 4.36e-5],
    [285.0, 999.5, 4193.0, 0.5835, 1.392, 12.40, 8.91, 0.000112],
    [290.0, 998.8, 4187.0, 0.5927, 1.417, 10.85, 7.66, 0.000172],
    [295.0, 997.8, 4183.0, 0.6017, 1.442, 9.600, 6.66, 0.000226],
    [300.0, 996.5, 4181.0, 0.6103, 1.465, 8.568, 5.85, 0.000275],
    [305.0, 995.0, 4180.0, 0.6184, 1.487, 7.708, 5.18, 0.000319],
    [310.0, 993.3, 4179.0, 0.6260, 1.508, 6.982, 4.63, 0.000361],
    [320.0, 989.3, 4181.0, 0.6396, 1.546, 5.832, 3.77, 0.000436],
    [340.0, 979.5, 4189.0, 0.6605, 1.610, 4.308, 2.68, 0.000565],
    [360.0, 967.4, 4202.0, 0.6737, 1.657, 3.371, 2.03, 0.000679],
    [373.15, 958.3, 4216.0, 0.6791, 1.681, 2.940, 1.75, 0.000751],
])
# Total emittance of black gloss paint near 40 C, Lienhard Table 10.1: 0.90.
EPS_BLACK_GLOSS = 0.90

RE_CRIT_PLATE = 5.0e5        # conventional transition Reynolds number for a flat plate (sharp transition model)


# ---------------------------------------------------------------- properties

def air_properties(T: float | np.ndarray) -> dict[str, np.ndarray]:
    """Air at 1 atm from Table A.6, linear interpolation in T (200-600 K; ValueError outside).
    Returns rho, cp, mu, nu, k, alpha, Pr (SI) and beta = 1/T (ideal gas)."""
    raise NotImplementedError


def water_properties(T: float | np.ndarray) -> dict[str, np.ndarray]:
    """Saturated liquid water from Table A.3, linear interpolation in T (275-373.15 K; ValueError
    outside). Returns rho, cp, k, alpha, nu, mu (= nu rho), Pr, beta (SI)."""
    raise NotImplementedError


def reynolds(U: float, L: float, nu: float) -> float:
    """Re = U L / nu."""
    raise NotImplementedError


def rayleigh(beta: float, dT: float, L: float, nu: float, alpha: float, g: float = G0) -> float:
    """Ra = g beta |dT| L^3 / (nu alpha)."""
    raise NotImplementedError


def h_from_nu(Nu: float, k: float, L: float) -> float:
    """h = Nu k / L."""
    raise NotImplementedError


def stanton(Nu: float, Re: float, Pr: float) -> float:
    """St = Nu / (Re Pr) = h / (rho cp U)."""
    raise NotImplementedError


# ---------------------------------------------------------------- flat plate (forced)

def plate_laminar_local_nu(Re_x, Pr):
    """Laminar flat plate, uniform wall temperature, local: Nu_x = 0.332 Re_x^(1/2) Pr^(1/3), Pr >= 0.6 (6.58)."""
    raise NotImplementedError


def plate_laminar_average_nu(Re_L, Pr):
    """Laminar flat plate, uniform wall temperature, average over 0..L: Nu_L = 0.664 Re_L^(1/2) Pr^(1/3) (6.68)."""
    raise NotImplementedError


def plate_laminar_cf_local(Re_x):
    """Blasius local skin-friction coefficient c_f = tau_w / (rho U^2 / 2) = 0.664 Re_x^(-1/2)."""
    raise NotImplementedError


def plate_turbulent_local_nu(Re_x, Pr):
    """Turbulent flat plate (gases), local, from the leading edge: Nu_x = 0.0296 Re_x^0.8 Pr^0.6 (6.112)."""
    raise NotImplementedError


def plate_turbulent_cf_local(Re_x):
    """One-fifth power law for the turbulent local skin friction, c_f = 0.0592 Re_x^(-1/5): the friction
    that the Reynolds-Colburn analogy pairs with 0.0296 Re_x^0.8 at Pr = 1."""
    raise NotImplementedError


def plate_local_nu(Re_x, Pr, Re_crit: float = RE_CRIT_PLATE):
    """Local Nu on a plate with a sharp transition at Re_crit: laminar (6.58) below, turbulent (6.112) above."""
    raise NotImplementedError


def plate_average_nu(Re_L, Pr, Re_crit: float = RE_CRIT_PLATE):
    """Average Nu over a plate of length L with a sharp transition at Re_crit (Lienhard 6.120 with
    Re_l = Re_u = Re_crit): laminar 0.664 Re_L^(1/2) Pr^(1/3) for Re_L <= Re_crit, otherwise
    0.037 Pr^0.6 (Re_L^0.8 - Re_crit^0.8) + 0.664 Re_crit^(1/2) Pr^(1/3)."""
    raise NotImplementedError


def plate_average_nu_numerical(Re_L: float, Pr: float, Re_crit: float = RE_CRIT_PLATE) -> float:
    """Independent route: average the local coefficient over the plate numerically,
    Nu_L = (L/k) (1/L) int_0^L h(x) dx = int_0^1 Nu_x(Re_L s) / s ds (s = x/L), by adaptive quadrature
    with the transition point as a breakpoint."""
    raise NotImplementedError


# ---------------------------------------------------------------- cylinder and pipe (forced)

def cylinder_crossflow_nu(Re_D, Pr):
    """Churchill-Bernstein, cylinder in cross-flow, all Re with Re Pr > 0.2 (7.65):
    Nu = 0.3 + 0.62 Re^(1/2) Pr^(1/3) / [1 + (0.4/Pr)^(2/3)]^(1/4) * [1 + (Re/282000)^(5/8)]^(4/5)."""
    raise NotImplementedError


def filonenko_f(Re_D):
    """Darcy friction factor for a smooth pipe, f = (1.82 log10 Re - 1.64)^(-2) (7.42)."""
    raise NotImplementedError


def gnielinski_nu(Re_D, Pr):
    """Gnielinski, turbulent pipe flow, 2300 <= Re <= 5e6, 0.6 <= Pr <= 1e5 (7.41):
    Nu = (f/8)(Re - 1000) Pr / (1 + 12.7 sqrt(f/8) (Pr^(2/3) - 1)). ValueError outside the Re range."""
    raise NotImplementedError


def dittus_boelter_nu(Re_D, Pr, heating: bool = True):
    """Dittus-Boelter, fully turbulent pipe flow: Nu = 0.023 Re^0.8 Pr^n, n = 0.4 heating the fluid, 0.3 cooling."""
    raise NotImplementedError


def pipe_nu(Re_D: float, Pr: float, wall: str = "q") -> float:
    """Fully developed pipe flow: laminar (Re < 2300) 4.364 for uniform heat flux ('q') or 3.657 for
    uniform wall temperature ('T') (7.22, 7.23); Gnielinski above. The jump at 2300 is real: the
    transition range is not predictable to better than a factor of two, and designs avoid it."""
    raise NotImplementedError


# ---------------------------------------------------------------- natural convection

def vertical_plate_nu(Ra_L, Pr):
    """Churchill-Chu, vertical isothermal plate of height L: laminar form (8.13a) for Ra <= 1e9,
    Nu = 0.68 + 0.670 Ra^(1/4) / [1 + (0.492/Pr)^(9/16)]^(4/9); the all-Ra form (8.13b) above,
    Nu = {0.825 + 0.387 Ra^(1/6) / [1 + (0.492/Pr)^(9/16)]^(8/27)}^2."""
    raise NotImplementedError


def vertical_plate_nu_allra(Ra_L, Pr):
    """Churchill-Chu all-Ra form (8.13b) on its own, for comparison with (8.13a)."""
    raise NotImplementedError


def horizontal_cylinder_nu(Ra_D, Pr):
    """Churchill-Chu, horizontal isothermal cylinder (8.28), 1e-6 <= Ra_D <= 1e9:
    Nu = 0.36 + 0.518 Ra^(1/4) / [1 + (0.559/Pr)^(9/16)]^(4/9)."""
    raise NotImplementedError


def horizontal_plate_up_nu(Ra_Ls, Pr):
    """Upper face of a hot plate (or lower face of a cold one), length scale L* = A/P.
    Raithby-Hollands (8.37a) Nu = 0.560 Ra^(1/4) / [1 + (0.492/Pr)^(9/16)]^(4/9), with the thick
    boundary-layer correction (8.37b) Nu = 1.4 / ln(1 + 1.4/Nu) when Nu < 10, for 1 < Ra* < 1e7;
    (8.35) Nu = 0.14 Ra^(1/3) (1 + 0.0107 Pr)/(1 + 0.01 Pr) above 1e7."""
    raise NotImplementedError


def horizontal_plate_down_nu(Ra_L):
    """Lower face of a hot plate (stable), plate width L: Fujii-Imura Nu = 0.58 Ra^(1/5) (8.38),
    stated for 1e6 <= Ra <= 1e11."""
    raise NotImplementedError


def h_natural(geometry: str, Tw: float, Tinf: float, L: float) -> float:
    """Natural-convection h for 'vertical' (L = height), 'cylinder' (L = diameter), 'up' (upper face of
    a hot plate, L = A/P) or 'down' (lower face of a hot plate, L = width), air properties at the film
    temperature, beta = 1/T_film. A cold surface (Tw < Tinf) swaps 'up' and 'down'."""
    raise NotImplementedError


def h_radiation_linear(eps: float, Tw: float, Tsur: float) -> float:
    """Exact linearised radiation coefficient for a small gray body in large surroundings (lesson P07.14):
    q = eps sigma (Tw^4 - Tsur^4) = h_rad (Tw - Tsur), h_rad = eps sigma (Tw^2 + Tsur^2)(Tw + Tsur)."""
    raise NotImplementedError


# ---------------------------------------------------------------- applications (given geometry)

def box_faces(L: float, W: float, H: float) -> list[tuple[str, float, float]]:
    """Faces of a box of depth L (front to back), width W and height H standing on a rack tray whose
    bottom is treated as adiabatic: (geometry, area, length scale) for the natural-convection
    correlations. Two side faces H x L and two end faces H x W (vertical, height H) and the top
    L x W (upper face of a hot plate, L* = A/P)."""
    top = L * W
    return [("vertical", 2 * H * L, H), ("vertical", 2 * H * W, H), ("up", top, top / (2 * (L + W)))]


def box_heat_loss_natural(Tw: float, Tinf: float, L: float, W: float, H: float,
                          eps: float = 0.0, Tsur: float | None = None) -> dict[str, float]:
    """Heat lost by an isothermal box at Tw in still air at Tinf: natural convection from each face plus
    gray radiation to large surroundings at Tsur (default Tinf). Returns Q_conv, Q_rad, Q, A, h_conv
    (area-weighted mean) and h_rad."""
    raise NotImplementedError


def box_surface_temperature(Q: float, Tinf: float, L: float, W: float, H: float, eps: float = 0.0) -> float:
    """Surface temperature at which the box rejects Q watts by natural convection and radiation:
    the root of box_heat_loss_natural(Tw) = Q (brentq on Tinf .. 590 K)."""
    raise NotImplementedError


def box_heat_loss_forced(U: float, Tw: float, Tinf: float, L: float, W: float, H: float,
                         eps: float = 0.0, Re_crit: float = RE_CRIT_PLATE) -> dict[str, float]:
    """Heat lost by the box in an air stream of speed U along its depth L: the two sides and the top
    are flat plates of length L (average Nu with a sharp transition); the ends and the tray are left
    out (a conservative simplification). Radiation as in box_heat_loss_natural."""
    raise NotImplementedError


def box_required_speed(Q: float, Tw: float, Tinf: float, L: float, W: float, H: float, eps: float = 0.0) -> float:
    """Air speed that holds the box at Tw while it dissipates Q (brentq on 1e-4 .. 100 m/s). Returns 0.0
    if natural convection and radiation alone already reject Q at Tw."""
    raise NotImplementedError


def cold_plate(Q: float, mdot: float, D: float, length: float, T_in: float) -> dict[str, float]:
    """A cold plate: Q watts into water flowing at mdot through a circular channel of diameter D and
    the given length, uniform heat flux along the channel. Energy balance T_out = T_in + Q/(mdot cp);
    fully developed h from pipe_nu (uniform flux) at the mean bulk temperature; the hottest wall is at
    the outlet, T_w,out = T_out + q''/h. Entry-length enhancement is ignored (conservative)."""
    raise NotImplementedError


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save convection.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    p300 = air_properties(300.0)
    print(f"air at 300 K (Table A.6): rho {float(p300['rho']):.4f} kg/m^3, k {float(p300['k']):.4f} W/m K, "
          f"nu {float(p300['nu']):.4e} m^2/s, Pr {float(p300['Pr']):.3f}")
    p315 = air_properties(315.0)
    print(f"air at 315 K (interpolated): nu {float(p315['nu']):.4e} m^2/s, k {float(p315['k']):.5f} W/m K, Pr {float(p315['Pr']):.4f}")

    # flat plate checks
    Re, Pr = 1.0e5, 0.707
    print(f"\nflat plate, Re_x = {Re:.0e}, Pr = {Pr}: laminar Nu_x {float(plate_laminar_local_nu(Re, Pr)):.2f}, "
          f"average Nu_L {float(plate_laminar_average_nu(Re, Pr)):.2f} (= 2 x local)")
    lhs = stanton(plate_laminar_local_nu(Re, Pr), Re, Pr) * Pr ** (2 / 3)
    print(f"Reynolds-Colburn, laminar: St Pr^(2/3) = {float(lhs):.6e}, c_f/2 = {float(plate_laminar_cf_local(Re)) / 2:.6e}")
    for Ret in (1e6, 1e7):
        lhs_t = stanton(plate_turbulent_local_nu(Ret, Pr), Ret, Pr) * Pr ** (2 / 3)
        print(f"Reynolds-Colburn, turbulent Re_x = {Ret:.0e}: St Pr^(2/3) / (c_f/2) = "
              f"{float(lhs_t / (plate_turbulent_cf_local(Ret) / 2)):.4f} (= Pr^(0.6 - 1/3) = {Pr ** (0.6 - 1/3):.4f})")
    for ReL in (2e5, 5e5, 2e6):
        print(f"plate average Re_L = {ReL:.0e}: closed form {float(plate_average_nu(ReL, Pr)):.2f}, "
              f"numerical average of local {plate_average_nu_numerical(ReL, Pr):.2f}")

    # avionics box
    L, W, H, Tinf, Q = 0.30, 0.12, 0.20, 298.15, 60.0
    print(f"\navionics box {L} x {W} x {H} m (bottom on an adiabatic tray), {Q:.0f} W, still air {Tinf - 273.15:.0f} C")
    Tw_c = box_surface_temperature(Q, Tinf, L, W, H, eps=0.0)
    Tw_r = box_surface_temperature(Q, Tinf, L, W, H, eps=EPS_BLACK_GLOSS)
    bc = box_heat_loss_natural(Tw_c, Tinf, L, W, H)
    br = box_heat_loss_natural(Tw_r, Tinf, L, W, H, EPS_BLACK_GLOSS)
    print(f"  area {bc['A']:.4f} m^2; convection only: Tw = {Tw_c - 273.15:.1f} C (h_conv {bc['h_conv']:.2f} W/m^2 K)")
    print(f"  black paint eps = {EPS_BLACK_GLOSS}: Tw = {Tw_r - 273.15:.1f} C, h_conv {br['h_conv']:.2f}, h_rad {br['h_rad']:.2f} W/m^2 K, "
          f"radiation carries {100 * br['Q_rad'] / br['Q']:.0f} %")
    for geom, area, Ls in box_faces(L, W, H):
        Tf = 0.5 * (Tw_r + Tinf)
        pp = air_properties(Tf)
        Ra = rayleigh(float(pp["beta"]), Tw_r - Tinf, Ls, float(pp["nu"]), float(pp["alpha"]))
        print(f"    face {geom:8s} area {area:.4f} m^2, L = {Ls:.4f} m, Ra = {Ra:.3e}, h = {h_natural(geom, Tw_r, Tinf, Ls):.2f} W/m^2 K")
    T_lim = 273.15 + 50.0
    U_req = box_required_speed(Q, T_lim, Tinf, L, W, H, EPS_BLACK_GLOSS)
    fr = box_heat_loss_forced(U_req, T_lim, Tinf, L, W, H, EPS_BLACK_GLOSS)
    print(f"  to hold Tw = 50 C: air speed {U_req:.3f} m/s (Re_L {fr['Re_L']:.0f}, h {fr['h']:.2f} W/m^2 K)")
    Tw_film = box_surface_temperature(Q, Tinf, L, W, H, eps=EPS_BLACK_GLOSS)
    print("  power sweep (black paint): Q [W] -> Tw [C]: " + ", ".join(
        f"{q:.0f} -> {box_surface_temperature(q, Tinf, L, W, H, EPS_BLACK_GLOSS) - 273.15:.1f}" for q in (20, 40, 60, 80, 100)))

    # cold plate
    Qc, Dc, Lc, Tin = 400.0, 0.004, 1.0, 293.15
    print(f"\ncold plate {Qc:.0f} W, water channel D = {Dc * 1e3:.0f} mm, length {Lc} m, inlet {Tin - 273.15:.0f} C")
    for mdot in (0.005, 0.010, 0.015, 0.020, 0.030):
        c = cold_plate(Qc, mdot, Dc, Lc, Tin)
        print(f"  mdot {mdot * 1e3:5.1f} g/s: Re {c['Re']:7.0f}, Nu {c['Nu']:6.2f}, h {c['h']:7.0f} W/m^2 K, "
              f"T_out {c['T_out'] - 273.15:5.2f} C, hottest wall {c['T_wall_out'] - 273.15:6.2f} C")

    # correlation comparison
    print("\nGnielinski against Dittus-Boelter (heating):")
    for Pr_ in (0.707, 5.85):
        row = []
        for Re_ in (1e4, 3e4, 1e5, 1e6):
            row.append(f"Re {Re_:.0e}: {float(dittus_boelter_nu(Re_, Pr_) / gnielinski_nu(Re_, Pr_)):.3f}")
        print(f"  Pr {Pr_}: DB/Gn " + ", ".join(row))
    print(f"Churchill-Chu vertical plate at Ra = 1e9, Pr = 0.71: (8.13a) {float(vertical_plate_nu(1e9, 0.71)):.1f}, "
          f"(8.13b) {float(vertical_plate_nu_allra(1e9, 0.71)):.1f}")

    # figure
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    ReL = np.logspace(3, 7, 200)
    ax[0].loglog(ReL, plate_laminar_average_nu(ReL, 0.707), "--", label="laminar average (6.68)")
    ax[0].loglog(ReL, plate_average_nu(ReL, 0.707), label="mixed average, transition 5e5")
    ax[0].loglog(ReL, plate_local_nu(ReL, 0.707), ":", label="local, sharp transition")
    ax[0].set_xlabel("Re"); ax[0].set_ylabel("Nu"); ax[0].set_title("flat plate, air"); ax[0].legend(fontsize=8)
    Qs = np.linspace(5, 120, 40)
    ax[1].plot(Qs, [box_surface_temperature(q, Tinf, L, W, H, 0.0) - 273.15 for q in Qs], label="convection only")
    ax[1].plot(Qs, [box_surface_temperature(q, Tinf, L, W, H, EPS_BLACK_GLOSS) - 273.15 for q in Qs], label="convection + radiation, eps 0.90")
    ax[1].axhline(50, color="k", lw=0.8, ls=":")
    ax[1].set_xlabel("dissipation Q (W)"); ax[1].set_ylabel("surface temperature (C)")
    ax[1].set_title("avionics box in still 25 C air"); ax[1].legend(fontsize=8)
    Res = np.logspace(np.log10(2300), 6, 100)
    for Pr_ in (0.707, 5.85):
        ax[2].loglog(Res, gnielinski_nu(Res, Pr_), label=f"Gnielinski Pr {Pr_}")
        ax[2].loglog(Res, dittus_boelter_nu(Res, Pr_), "--", label=f"Dittus-Boelter Pr {Pr_}")
    ax[2].set_xlabel("Re_D"); ax[2].set_ylabel("Nu_D"); ax[2].set_title("turbulent pipe flow"); ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "convection.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'convection.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
