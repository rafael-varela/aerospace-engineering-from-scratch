"""P08.16 Flow Separation & Drag: starter code.

How much drag does a bluff body make, why, and how much does streamlining save? This module builds the
tools for bluff-body drag at low speed:

* pressure drag from a surface pressure distribution, by integrating c_p cos(theta) round a circular
  cylinder: zero for the attached potential flow (d'Alembert's paradox, lesson P08.08) and close to the
  measured value once a separated wake with a low base pressure is modelled;
* the drag coefficient of a sphere as a function of Reynolds number from three CORRELATIONS (curve fits
  to experiments, not measurements): Stokes's 24/Re (exact as Re -> 0), Schiller-Naumann (Re < 1000,
  then 0.44) and Morrison (2013, up to Re = 1e6, including the drag crisis), and White's fit for a
  cylinder;
* vortex shedding: the Strouhal number of a circular cylinder from Roshko's correlations, and the
  shedding frequency f = St V/d;
* terminal velocity by root finding on the force balance, checked against Stokes's closed form in
  the creeping limit and against a time integration of the equation of motion;
* streamlining: a round strut against a faired one of the same thickness, with the fairing's drag from
  a turbulent flat-plate skin friction and a classical form-factor estimate (a correlation).

Units: SI. theta is measured from the front stagnation point.

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes drag.png)
Test: pytest curriculum/phase-08/lesson-16-flow-separation-drag
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

G0 = 9.80665
RHO_AIR, MU_AIR = 1.204, 1.81e-5            # air at 20 C, 101 325 Pa
RHO_WATER, MU_WATER = 998.2, 1.002e-3       # water at 20 C


# ---------------------------------------------------------------- given helpers

def isa_troposphere(h: float) -> tuple[float, float, float, float]:
    """ISA troposphere (geopotential h <= 11 km): T, p, rho and Sutherland viscosity mu (lesson P08.01)."""
    T = 288.15 - 0.0065 * h
    p = 101325.0 * (T / 288.15) ** (G0 / (287.05287 * 0.0065))
    rho = p / (287.05287 * T)
    mu = 1.458e-6 * T ** 1.5 / (T + 110.4)
    return T, p, rho, mu


def reynolds(rho: float, V: float | np.ndarray, d: float, mu: float) -> float | np.ndarray:
    """Reynolds number rho V d / mu."""
    return rho * np.asarray(V, dtype=float) * d / mu


def flat_plate_cf_turbulent(re_x: float | np.ndarray) -> float | np.ndarray:
    """Mean skin-friction coefficient of a fully turbulent flat plate, both sides counted separately:
    C_f = 0.074 Re^(-1/5) per side (Prandtl's 1/7-power-law estimate, lesson P08.14; a correlation)."""
    return 0.074 * np.asarray(re_x, dtype=float) ** -0.2


# ---------------------------------------------------------------- pressure drag

def cp_potential_cylinder(theta: np.ndarray) -> np.ndarray:
    """Surface pressure coefficient of attached potential flow round a circular cylinder,
    c_p = 1 - 4 sin^2(theta) (lesson P08.08)."""
    raise NotImplementedError


def cp_separated_model(theta: np.ndarray, theta_sep: float, cp_base: float) -> np.ndarray:
    """A simple separated-flow model: potential-flow c_p ahead of the separation angle theta_sep on each
    side (|theta| < theta_sep, theta folded into [0, pi]), and a constant base pressure cp_base in the
    wake behind it."""
    raise NotImplementedError


def pressure_drag_coefficient(theta: np.ndarray, cp: np.ndarray) -> float:
    """Pressure-drag coefficient per unit span of a circular cylinder (reference area d x span) from a
    surface c_p(theta) sampled on a uniform periodic grid theta_k = 2 pi k/n, k = 0..n-1:
    c_d = (1/2) integral_0^{2 pi} c_p cos(theta) dtheta (the trapezoid rule, spectrally accurate for a
    smooth periodic integrand)."""
    raise NotImplementedError


def pressure_lift_coefficient(theta: np.ndarray, cp: np.ndarray) -> float:
    """Side-force coefficient from the same integral with -sin(theta): zero for any c_p symmetric about
    the stagnation line."""
    raise NotImplementedError


# ---------------------------------------------------------------- drag correlations

def cd_sphere_stokes(re: float | np.ndarray) -> np.ndarray:
    """Stokes drag of a sphere, c_d = 24/Re: the exact creeping-flow result (Re << 1)."""
    raise NotImplementedError


def cd_sphere_schiller_naumann(re: float | np.ndarray) -> np.ndarray:
    """Schiller-Naumann CORRELATION: c_d = (24/Re)(1 + 0.15 Re^0.687) for Re <= 1000, and the
    Newton-regime constant 0.44 above (it knows nothing about the drag crisis)."""
    raise NotImplementedError


def cd_sphere_morrison(re: float | np.ndarray) -> np.ndarray:
    """Morrison's (2013) CORRELATION for a smooth sphere, fitted for Re up to 1e6, including the drag
    crisis near Re = 3e5:
    c_d = 24/Re + 2.6 (Re/5)/(1 + (Re/5)^1.52) + 0.411 (Re/2.63e5)^-7.94/(1 + (Re/2.63e5)^-8)
          + 0.25 (Re/1e6)/(1 + Re/1e6)."""
    raise NotImplementedError


def cd_cylinder_white(re: float | np.ndarray) -> np.ndarray:
    """White's CORRELATION for a smooth circular cylinder in cross-flow, c_d = 1 + 10 Re^(-2/3),
    for 1 < Re < 2e5 (subcritical; it does not contain the cylinder's drag crisis)."""
    raise NotImplementedError


def drag_crisis(cd: Callable = cd_sphere_morrison, lo: float = 1e5, hi: float = 1e6) -> tuple[float, float]:
    """The Reynolds number and value of the minimum of c_d(Re) between lo and hi (the bottom of the drag
    crisis), by a bounded scalar minimisation in log10(Re)."""
    raise NotImplementedError


# ---------------------------------------------------------------- vortex shedding

def strouhal_cylinder(re: float | np.ndarray) -> np.ndarray:
    """Strouhal number St = f d/V of a circular cylinder from Roshko's (1954) CORRELATIONS:
    0.212 (1 - 21.2/Re) for 50 <= Re <= 150 (laminar shedding) and 0.212 (1 - 12.7/Re) for
    300 <= Re (fitted to 2000, used here as the typical subcritical value up to 2e5, where St is
    0.20-0.21). Linear interpolation across the irregular range 150 < Re < 300; NaN below Re = 47,
    where the wake is steady and nothing is shed."""
    raise NotImplementedError


def shedding_frequency(V: float, d: float, rho: float, mu: float) -> float:
    """Vortex-shedding frequency f = St(Re) V/d (Hz) of a circular cylinder of diameter d."""
    raise NotImplementedError


# ---------------------------------------------------------------- terminal velocity

def stokes_terminal_velocity(d: float, rho_p: float, rho_f: float, mu: float) -> float:
    """Creeping-flow terminal velocity of a sphere, V = (rho_p - rho_f) g d^2/(18 mu) (from 24/Re)."""
    raise NotImplementedError


def terminal_velocity(d: float, rho_p: float, rho_f: float, mu: float,
                      cd: Callable = cd_sphere_morrison) -> tuple[float, float, float]:
    """Terminal velocity of a sphere: the root of (rho_p - rho_f) g (pi d^3/6) = (1/2) rho_f V^2 c_d(Re)
    (pi d^2/4), Re = rho_f V d/mu, found by brentq between the Stokes velocity's 1e-6 and 1e3 times
    (the residual changes sign once). Returns (V, Re, c_d)."""
    raise NotImplementedError


def simulate_fall(d: float, rho_p: float, rho_f: float, mu: float, t_end: float,
                  cd: Callable = cd_sphere_morrison) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Integrate the equation of motion of a sphere released from rest, m dV/dt = (m - m_f) g - D,
    D = (1/2) rho_f V^2 c_d(Re) A (added mass and history force neglected), with dz/dt = V (z down).
    Returns t, V, z on 400 output times (solve_ivp, RK45, rtol 1e-9)."""
    raise NotImplementedError


# ---------------------------------------------------------------- streamlining

def fairing_cd(t_over_c: float, re_c: float) -> float:
    """Drag coefficient of a 2-D streamlined (faired) section based on its THICKNESS t, from a classical
    form-factor CORRELATION: c_d,chord = 2 C_f (1 + 2 t/c + 60 (t/c)^4), with C_f the turbulent flat-plate
    coefficient at the chord Reynolds number re_c (both surfaces: factor 2), then c_d,t = c_d,chord c/t.
    Valid for attached flow, t/c roughly 0.1-0.5."""
    raise NotImplementedError


def strut_drag(d: float, length: float, V: float, rho: float, mu: float, t_over_c: float = 0.25) -> dict:
    """A round strut (diameter d) against a faired strut of the same thickness and chord d/t_over_c.
    Returns Re_d, c_d (cylinder, White), drag (N) and power (W) for each, the drag ratio, the shedding
    frequency of the round strut, and the diameter of a round wire with the same drag as the fairing."""
    raise NotImplementedError


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save drag.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    T, p, rho, mu = isa_troposphere(3000.0)
    print(f"ISA 3000 m: T = {T:.2f} K, p = {p:.0f} Pa, rho = {rho:.4f} kg/m^3, mu = {mu:.4e} Pa s")
    s = strut_drag(0.025, 1.0, 60.0, rho, mu)
    print("1 m of 25 mm round strut on the light aeroplane at 60 m/s TAS, 3000 m:")
    print(f"  Re_d = {s['re_d']:.0f}, c_d = {s['cd_round']:.3f} (White correlation), drag {s['drag_round']:.1f} N, "
          f"power {s['power_round'] / 1e3:.2f} kW")
    print(f"  shedding: St = {float(strouhal_cylinder(s['re_d'])):.3f} (Roshko correlation), f = {s['f_shed']:.0f} Hz")
    print(f"  faired, t/c = 0.25 (chord 100 mm, Re_c = {s['re_chord']:.0f}): c_d (on thickness) = {s['cd_faired']:.4f}, "
          f"drag {s['drag_faired']:.2f} N, power {s['power_faired']:.0f} W")
    print(f"  fairing cuts the drag by a factor {s['ratio']:.1f}; same drag as a round wire of "
          f"{s['equivalent_wire_d'] * 1e3:.2f} mm")

    th = 2 * np.pi * np.arange(2048) / 2048
    print("\npressure drag of a cylinder from c_p(theta):")
    print(f"  attached potential flow: c_d = {pressure_drag_coefficient(th, cp_potential_cylinder(th)):.1e} "
          "(d'Alembert)")
    for ts, cpb in ((np.radians(82), -1.0), (np.radians(82), -1.2), (np.radians(120), -0.2)):
        print(f"  separation at {np.degrees(ts):.0f} deg, base c_p {cpb}: c_d = "
              f"{pressure_drag_coefficient(th, cp_separated_model(th, ts, cpb)):.3f}")

    print("\nsphere c_d correlations:")
    print("        Re      Stokes   Schiller-N   Morrison")
    for re in (0.1, 1.0, 10.0, 100.0, 1e3, 1e4, 1e5, 2.63e5, 4e5, 1e6):
        print(f"  {re:9.3g}   {float(cd_sphere_stokes(re)):8.3f}   {float(cd_sphere_schiller_naumann(re)):8.3f}"
              f"   {float(cd_sphere_morrison(re)):8.3f}")
    rc, cmin = drag_crisis()
    print(f"  Morrison drag-crisis minimum: c_d = {cmin:.4f} at Re = {rc:.3g}; "
          f"c_d at 1e5 is {float(cd_sphere_morrison(1e5)) / cmin:.1f} times larger")

    print("\nterminal velocities (Morrison correlation):")
    cases = [("ping-pong ball, 40 mm, 2.7 g, air", 0.040, 2.7e-3 / (np.pi * 0.040 ** 3 / 6), RHO_AIR, MU_AIR),
             ("glass bead 0.1 mm in water", 1e-4, 2500.0, RHO_WATER, MU_WATER),
             ("steel ball 10 mm in water", 0.010, 7850.0, RHO_WATER, MU_WATER),
             ("water drop 1 mm in air at 3000 m", 0.001, RHO_WATER, rho, mu)]
    for name, d, rp, rf, m in cases:
        V, re, cd = terminal_velocity(d, rp, rf, m)
        Vsn, _, _ = terminal_velocity(d, rp, rf, m, cd=cd_sphere_schiller_naumann)
        t, Vt, z = simulate_fall(d, rp, rf, m, 15.0 * V / (G0 * (1.0 - rf / rp)))
        print(f"  {name}: V_t = {V:.4g} m/s (Re = {re:.3g}, c_d = {cd:.3f}); Schiller-Naumann {Vsn:.4g} m/s; "
              f"Stokes {stokes_terminal_velocity(d, rp, rf, m):.4g} m/s; ODE at t = {t[-1]:.3g} s: {Vt[-1]:.4g} m/s")

    re = np.logspace(-1, 6, 400)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].loglog(re, cd_sphere_morrison(re), label="sphere, Morrison 2013 (correlation)")
    ax[0].loglog(re, cd_sphere_schiller_naumann(re), "--", label="sphere, Schiller-Naumann (correlation)")
    ax[0].loglog(re[re < 10], cd_sphere_stokes(re[re < 10]), "k:", label="Stokes 24/Re (exact, Re << 1)")
    rc_ = re[(re > 1) & (re < 2e5)]
    ax[0].loglog(rc_, cd_cylinder_white(rc_), "-.", label="cylinder, White (correlation)")
    ax[0].set_xlabel("Re")
    ax[0].set_ylabel("c_d")
    ax[0].set_title("drag coefficient (correlations, not data)")
    ax[0].legend(fontsize=8)
    deg = np.degrees(th)
    ax[1].plot(deg, cp_potential_cylinder(th), "k--", label="attached potential flow, c_d = 0")
    ax[1].plot(deg, cp_separated_model(th, np.radians(82), -1.2), label="separation at 82 deg, c_pb = -1.2 (model)")
    ax[1].set_xlim(0, 180)
    ax[1].set_xlabel("theta from the front stagnation point (deg)")
    ax[1].set_ylabel("c_p")
    ax[1].set_title("cylinder surface pressure")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "drag.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'drag.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
