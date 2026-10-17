"""P04.02 Work, Energy & Power: starter code.

Energy bookkeeping for an aircraft treated as a point mass. The total mechanical energy per unit
weight is the energy height h_E = h + V^2 / (2 g); the rate at which the engine-minus-drag adds to it
is the specific excess power P_s = (T - D) V / W = dh_E/dt. From these two definitions follow the
rate of climb, the ceilings, the zoom climb, the glide, the phugoid, and the minimum-time climb
path (Rutowski's energy method). Every result is checked against a closed form or against a
`scipy.integrate.solve_ivp` integration of the longitudinal point-mass equations of motion.

Model: ISA troposphere and lower stratosphere (to 20 km), a parabolic drag polar
C_D = C_D0 + k C_L^2, thrust T = T0 (rho/rho0)^n independent of speed, no compressibility (the
speeds used stay below Mach 0.75), flat non-rotating Earth, constant g and constant mass.

Run:  python main.py [out_dir]     (prints the tables, writes energy_height.png)
Test: pytest curriculum/phase-04/lesson-02-work-energy-power
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

G0 = 9.80665              # m/s^2
R_AIR = 287.05287         # J/(kg K)
RHO0 = 1.225              # kg/m^3, ISA sea level


# ---------------------------------------------------------------- the atmosphere and the aircraft

def isa(h: float) -> tuple[float, float, float]:
    """ISA temperature [K], pressure [Pa] and density [kg/m^3] at geopotential altitude h [m],
    0 <= h <= 20 km: troposphere with lapse rate 6.5 K/km to 11 km, isothermal 216.65 K above."""
    if h < 0.0 or h > 20_000.0:
        raise ValueError("isa() covers 0-20 km")
    if h <= 11_000.0:
        T = 288.15 - 0.0065 * h
        p = 101325.0 * (T / 288.15) ** (G0 / (0.0065 * R_AIR))
    else:
        T = 216.65
        p11 = 101325.0 * (216.65 / 288.15) ** (G0 / (0.0065 * R_AIR))
        p = p11 * math.exp(-G0 * (h - 11_000.0) / (R_AIR * T))
    return T, p, p / (R_AIR * T)


def speed_of_sound(h: float) -> float:
    """a = sqrt(gamma R T) with gamma = 1.4 [m/s]."""
    return math.sqrt(1.4 * R_AIR * isa(h)[0])


@dataclass(frozen=True)
class Aircraft:
    """A twin-engine jet trainer, point-mass performance model."""
    mass: float = 5000.0        # kg
    S: float = 20.0             # m^2 wing area
    CD0: float = 0.022          # zero-lift drag coefficient
    k: float = 0.060            # induced-drag factor
    T0: float = 16_000.0        # N, total sea-level static thrust
    thrust_exp: float = 1.0     # T = T0 (rho/rho0)^n
    CLmax: float = 1.4
    M_max: float = 0.75         # speed limit as a Mach number (no compressibility model beyond it)

    @property
    def W(self) -> float:
        return self.mass * G0


# ---------------------------------------------------------------- work and energy

def work_along_path(force: Callable[[np.ndarray, np.ndarray], np.ndarray], path: np.ndarray) -> float:
    """Work W = integral of F . dr along a sampled path (shape (n, d)), by the trapezoid rule on the
    segments: sum over i of (F(r_i, dr_i) + F(r_{i+1}, dr_i)) / 2 . dr_i with dr_i = r_{i+1} - r_i.

    `force(r, dr)` returns the force at position r; it receives the segment direction dr so that a
    velocity-dependent force such as drag (opposite to the motion) can be expressed.
    """
    raise NotImplementedError


def kinetic_energy(mass: float, V: np.ndarray | float) -> np.ndarray | float:
    """K = m V^2 / 2 [J]."""
    raise NotImplementedError


def energy_height(h: np.ndarray | float, V: np.ndarray | float, g: float = G0) -> np.ndarray | float:
    """Specific energy (energy per unit weight) h_E = h + V^2 / (2 g) [m]."""
    raise NotImplementedError


def zoom_altitude(h: float, V: float, V_end: float, g: float = G0) -> float:
    """Altitude reached in an ideal zoom (no thrust, no drag) from (h, V) to speed V_end:
    h_end = h_E - V_end^2 / (2 g). The energy height is conserved."""
    raise NotImplementedError


# ---------------------------------------------------------------- forces and power

def drag(ac: Aircraft, rho: float, V: np.ndarray | float, n: float = 1.0) -> np.ndarray | float:
    """D = q S (C_D0 + k C_L^2) with C_L = n W / (q S), q = rho V^2 / 2 [N]. n is the load factor L/W."""
    raise NotImplementedError


def thrust(ac: Aircraft, rho: float) -> float:
    """T = T0 (rho / rho0)^n [N], independent of speed."""
    raise NotImplementedError


def specific_excess_power(ac: Aircraft, h: np.ndarray | float, V: np.ndarray | float, n: float = 1.0) -> np.ndarray:
    """P_s = (T - D) V / W [m/s] at altitude h and true airspeed V; vectorised over h and V
    (broadcast together). The rate of change of the energy height in flight at load factor n."""
    raise NotImplementedError


def power_required_available(ac: Aircraft, h: float, V: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Power required P_r = D V and available P_a = T V [W] in level flight at altitude h."""
    raise NotImplementedError


def min_drag_speed(ac: Aircraft, h: float) -> float:
    """Speed of minimum drag in level flight, where C_D0 = k C_L^2: V_md = sqrt(2W/(rho S)) (k/C_D0)^(1/4)."""
    raise NotImplementedError


def min_power_speed(ac: Aircraft, h: float) -> float:
    """Speed of minimum power required, where C_D0 = k C_L^2 / 3: V_mp = V_md / 3^(1/4)."""
    raise NotImplementedError


def max_lift_to_drag(ac: Aircraft) -> float:
    """(L/D)_max = 1 / (2 sqrt(k C_D0)): the best glide ratio."""
    raise NotImplementedError


def stall_speed(ac: Aircraft, h: float) -> float:
    """V_s = sqrt(2 W / (rho S C_Lmax))."""
    return math.sqrt(2.0 * ac.W / (isa(h)[2] * ac.S * ac.CLmax))


def best_climb(ac: Aircraft, h: float) -> tuple[float, float]:
    """Speed and value of the maximum specific excess power at altitude h, over stall speed to the
    speed limit V_max = M_max a(h) (bounded scalar search, xatol 1e-4 m/s). In a steady climb at constant
    speed this P_s is the rate of climb."""
    raise NotImplementedError


def ceiling(ac: Aircraft, ps_target: float = 0.0) -> float:
    """Altitude at which the best P_s equals `ps_target` (0: absolute ceiling; 0.508 m/s = 100 ft/min:
    service ceiling), by `brentq` on 0-20 km. Raises ValueError if the aircraft cannot reach the
    target even at sea level, or exceeds it at 20 km."""
    raise NotImplementedError


# ---------------------------------------------------------------- the equations of motion

def point_mass_rhs(t: float, y: np.ndarray, ac: Aircraft, CL: float, thrust_on: bool = True,
                   drag_on: bool = True, const_rho: float | None = None) -> np.ndarray:
    """Longitudinal point-mass equations for y = [x, h, V, gamma] at constant lift coefficient:

        dx/dt = V cos(gamma),   dh/dt = V sin(gamma),
        dV/dt = (T - D)/m - g sin(gamma),   dgamma/dt = (L - W cos(gamma)) / (m V),

    with L = q S CL, D = q S (CD0 + k CL^2), T = thrust(rho(h)) (thrust along the flight path).
    `thrust_on=False` sets T = 0; `drag_on=False` sets D = 0 (an idealised glider); `const_rho`, if given,
    replaces the ISA density by a constant (the assumption behind Lanchester's phugoid).
    """
    raise NotImplementedError


@dataclass
class Flight:
    t: np.ndarray
    x: np.ndarray
    h: np.ndarray
    V: np.ndarray
    gamma: np.ndarray


def simulate_flight(ac: Aircraft, h0: float, V0: float, gamma0: float, CL: float, t_end: float,
                    thrust_on: bool = True, drag_on: bool = True, rtol: float = 1e-10, n_samples: int = 2001,
                    const_rho: float | None = None) -> Flight:
    """Integrate `point_mass_rhs` from (x, h, V, gamma) = (0, h0, V0, gamma0) to t_end with
    `solve_ivp` (DOP853, rtol given, atol 1e-9), sampled at n_samples equally spaced times."""
    raise NotImplementedError


def phugoid_period_lanchester(V: float, g: float = G0) -> float:
    """Lanchester's phugoid period for a drag-free aircraft at constant C_L: T = sqrt(2) pi V / g."""
    raise NotImplementedError


def density_scale_gradient(h: float, dh: float = 1.0) -> float:
    """beta = -(1/rho) d rho/dh [1/m] by a central difference of the ISA density."""
    raise NotImplementedError


def phugoid_period_with_gradient(V: float, beta: float, g: float = G0) -> float:
    """Small-amplitude, drag-free phugoid period at constant C_L when the density falls with height:
    omega^2 = 2 g^2 / V^2 + g beta, T = 2 pi / omega. beta = 0 recovers Lanchester."""
    raise NotImplementedError


def oscillation_period(t: np.ndarray, signal: np.ndarray) -> float:
    """Mean period of an oscillating signal from its upward mean-crossings (linear interpolation
    between samples). Needs at least two upward crossings; raises ValueError otherwise."""
    raise NotImplementedError


# ---------------------------------------------------------------- minimum-time climb

def rutowski_path(ac: Aircraft, hE_grid: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """For each energy height h_E in the grid, the speed V* that maximises P_s on the constant-h_E
    curve h = h_E - V^2/(2g) (subject to V >= stall speed at that altitude, V <= M_max a(h) and h >= 0),
    and the altitude and P_s there. Returns (V*, h*, P_s*). This is the minimum-time climb schedule of the
    energy method."""
    raise NotImplementedError


def time_to_energy(hE_grid: np.ndarray, ps: np.ndarray) -> float:
    """Time to climb through the energy heights in the grid at the given P_s at each: the integral of
    dh_E / P_s by the trapezoid rule. Raises ValueError if any P_s <= 0 (the energy cannot be reached)."""
    raise NotImplementedError


def constant_speed_climb_ps(ac: Aircraft, V: float, hE_grid: np.ndarray) -> np.ndarray:
    """P_s along a climb at constant true airspeed V through the given energy heights (h = h_E - V^2/2g)."""
    raise NotImplementedError


def climb_time_accelerate_then_constant(ac: Aircraft, V_start: float, V_climb: float, hE_end: float, n: int = 161) -> float:
    """Time for a two-segment climb: change speed in level flight at sea level from V_start to V_climb
    (if V_climb < V_start the energy given up is regained later, so this segment takes the time
    to decelerate as zero and starts the climb from the lower energy), then climb at constant V_climb
    to the energy height hE_end. Uses `time_to_energy` on both segments."""
    hE0 = float(energy_height(0.0, V_start))
    hE1 = float(energy_height(0.0, V_climb))
    t1 = 0.0
    if V_climb > V_start:
        Vs = np.linspace(V_start, V_climb, n)
        t1 = time_to_energy(energy_height(0.0, Vs), specific_excess_power(ac, 0.0, Vs))
    grid = np.linspace(hE1, hE_end, n)
    return t1 + time_to_energy(grid, constant_speed_climb_ps(ac, V_climb, grid))


def climb_rate_with_acceleration(ps: float, V: float, dV_dh: float, g: float = G0) -> float:
    """dh/dt when the speed changes with altitude along the climb: from dh_E/dt = P_s and
    h_E = h + V^2/2g, dh/dt = P_s / (1 + (V/g) dV/dh). The bracket is the acceleration factor."""
    raise NotImplementedError


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save energy_height.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    ac = Aircraft()
    print(f"jet trainer: m = {ac.mass:.0f} kg, W = {ac.W / 1e3:.2f} kN, S = {ac.S} m^2, CD0 = {ac.CD0}, k = {ac.k}, "
          f"T0 = {ac.T0 / 1e3:.0f} kN (T/W = {ac.T0 / ac.W:.3f}), (L/D)max = {max_lift_to_drag(ac):.2f}")
    h, V = 0.0, 150.0
    rho = isa(h)[2]
    T, D = thrust(ac, rho), float(drag(ac, rho, V))
    print(f"at sea level, 150 m/s: T = {T / 1e3:.2f} kN, D = {D / 1e3:.3f} kN, h_E = {float(energy_height(h, V)):.1f} m, "
          f"P_s = {float(specific_excess_power(ac, h, V)):.2f} m/s; KE = {float(kinetic_energy(ac.mass, V)) / 1e6:.2f} MJ")
    print(f"ideal zoom from (0 m, 150 m/s) to 80 m/s: {zoom_altitude(0.0, 150.0, 80.0):.0f} m")
    print(f"V_md(0) = {min_drag_speed(ac, 0.0):.1f} m/s, V_mp(0) = {min_power_speed(ac, 0.0):.1f} m/s, V_stall(0) = {stall_speed(ac, 0.0):.1f} m/s")
    print("\nbest climb against altitude:")
    for hh in (0.0, 3000.0, 6000.0, 9000.0, 11000.0, 12000.0):
        vb, pb = best_climb(ac, hh)
        print(f"  h = {hh / 1e3:4.1f} km: V* = {vb:6.1f} m/s, P_s,max = {pb:6.2f} m/s, T = {thrust(ac, isa(hh)[2]) / 1e3:5.2f} kN")
    h_abs, h_srv = ceiling(ac), ceiling(ac, 0.508)
    print(f"absolute ceiling {h_abs:.0f} m, service ceiling (100 ft/min) {h_srv:.0f} m")

    # the work integral: gravity is conservative, drag is not
    m = ac.mass
    grav = lambda r, dr: np.array([0.0, -m * G0])
    straight = np.column_stack([np.linspace(0, 1000, 201), np.linspace(0, 500, 201)])
    curved = np.column_stack([np.linspace(0, 1000, 201), 500 * np.sin(np.linspace(0, np.pi / 2, 201)) ** 2 + 300 * np.sin(np.linspace(0, 3 * np.pi, 201))])
    print(f"\nwork of gravity, straight vs wavy path to the same point: {work_along_path(grav, straight) / 1e6:.4f} MJ, "
          f"{work_along_path(grav, curved) / 1e6:.4f} MJ (-m g dz = {-m * G0 * 500 / 1e6:.4f} MJ)")
    Dc = 5000.0
    dragf = lambda r, dr: -Dc * dr / np.linalg.norm(dr)
    print(f"work of a constant 5 kN drag: straight {work_along_path(dragf, straight) / 1e6:.3f} MJ, wavy {work_along_path(dragf, curved) / 1e6:.3f} MJ")

    # phugoid: drag-free and thrust-free, constant CL at 150 m/s level trim
    V0 = 150.0; h0 = 3000.0
    rho0 = isa(h0)[2]
    CL0 = ac.W / (0.5 * rho0 * V0 ** 2 * ac.S)
    beta = density_scale_gradient(h0)
    flc = simulate_flight(ac, h0, V0 * 1.01, 0.0, CL0, 300.0, thrust_on=False, drag_on=False, const_rho=rho0)
    fl = simulate_flight(ac, h0, V0 * 1.01, 0.0, CL0, 300.0, thrust_on=False, drag_on=False)
    hE = energy_height(fl.h, fl.V)
    print(f"\nphugoid at 150 m/s, 3 km, CL = {CL0:.4f}, started 1 % fast, no thrust and no drag: h_E drift {np.ptp(hE):.1e} m over 300 s; "
          f"altitude swing {np.ptp(fl.h):.1f} m, speed swing {np.ptp(fl.V):.2f} m/s")
    print(f"  period with constant density {oscillation_period(flc.t, flc.h):.2f} s (Lanchester sqrt(2) pi V/g = {phugoid_period_lanchester(V0):.2f} s); "
          f"with the ISA density {oscillation_period(fl.t, fl.h):.2f} s (beta = {beta * 1e3:.4f} /km, formula {phugoid_period_with_gradient(V0, beta):.2f} s)")
    fl5 = simulate_flight(ac, h0, V0 * 1.10, 0.0, CL0, 300.0, thrust_on=False, drag_on=False, const_rho=rho0)
    print(f"  started 10 % fast (constant density): period {oscillation_period(fl5.t, fl5.h):.2f} s, altitude swing {np.ptp(fl5.h):.0f} m")
    fl2 = simulate_flight(ac, h0, V0 * 1.05, 0.0, CL0, 200.0, thrust_on=True, drag_on=True, n_samples=4001)
    hE2 = energy_height(fl2.h, fl2.V)
    dhE = np.gradient(hE2, fl2.t)
    ps_path = np.array([(thrust(ac, isa(hh)[2]) - 0.5 * isa(hh)[2] * v * v * ac.S * (ac.CD0 + ac.k * CL0 ** 2)) * v / ac.W for hh, v in zip(fl2.h, fl2.V)])
    print(f"  with thrust and drag: max |dh_E/dt - (T - D)V/W| = {np.max(np.abs(dhE[2:-2] - ps_path[2:-2])):.1e} m/s (P_s from {ps_path.min():.1f} to {ps_path.max():.1f} m/s)")

    # glide
    rho3 = isa(3000.0)[2]
    CLg = math.sqrt(ac.CD0 / ac.k)
    gam_g = -math.atan(1.0 / max_lift_to_drag(ac))
    Vg = math.sqrt(2.0 * ac.W * math.cos(gam_g) / (rho3 * ac.S * CLg))
    flg = simulate_flight(ac, 3000.0, Vg, gam_g, CLg, 60.0, thrust_on=False, const_rho=rho3)
    flgi = simulate_flight(ac, 3000.0, Vg, gam_g, CLg, 300.0, thrust_on=False)
    print(f"glide at (L/D)max = {max_lift_to_drag(ac):.2f}, CL = {CLg:.3f}, V = {Vg:.1f} m/s at 3 km: glide angle {math.degrees(gam_g):.2f} deg; "
          f"distance/height lost over 60 s at constant density {flg.x[-1] / (3000.0 - flg.h[-1]):.4f}; in the ISA over 300 s "
          f"{flgi.x[-1] / (3000.0 - flgi.h[-1]):.2f} (height lost {3000.0 - flgi.h[-1]:.0f} m, speed {flgi.V[-1]:.1f} m/s at the end)")

    # minimum-time climb to h_E = 10 km, against accelerate-then-constant-speed climbs
    hE_start, hE_end = float(energy_height(0.0, 150.0)), 10_000.0
    hE_grid = np.linspace(hE_start, hE_end, 161)
    Vr, hr, pr = rutowski_path(ac, hE_grid)
    t_r = time_to_energy(hE_grid, pr)
    print(f"\nminimum-time (energy) climb from (0 m, 150 m/s) to h_E = 10 km: {t_r:.0f} s; schedule:")
    for e in (2000.0, 4000.0, 6000.0, 8000.0, 10000.0):
        i = int(np.argmin(np.abs(hE_grid - e)))
        print(f"  h_E = {hE_grid[i]:6.0f} m: V* = {Vr[i]:6.1f} m/s, h = {hr[i]:6.0f} m, P_s = {pr[i]:5.2f} m/s")
    for Vc in (120.0, 150.0, 170.0):
        t_tot = climb_time_accelerate_then_constant(ac, 150.0, Vc, hE_end)
        print(f"  level speed change at sea level to {Vc:.0f} m/s, then constant {Vc:.0f} m/s: {t_tot:.0f} s ({100 * (t_tot / t_r - 1):.0f} % longer)")
    print(f"acceleration factor: climbing at P_s = 20 m/s while accelerating 5 m/s per km at 150 m/s gives dh/dt = "
          f"{climb_rate_with_acceleration(20.0, 150.0, 5e-3):.2f} m/s")

    # sensitivity: ceiling against weight and thrust
    print("\nsensitivity of the service ceiling:")
    for f_m, f_t in ((1.0, 1.0), (1.1, 1.0), (0.9, 1.0), (1.0, 1.1), (1.0, 0.9)):
        a2 = replace(ac, mass=ac.mass * f_m, T0=ac.T0 * f_t)
        print(f"  mass x{f_m}, thrust x{f_t}: {ceiling(a2, 0.508):.0f} m")

    # figure
    Vg_ = np.linspace(50.0, 240.0, 77)
    hg_ = np.linspace(0.0, 13_000.0, 66)
    VV, HH = np.meshgrid(Vg_, hg_)
    PS = specific_excess_power(ac, HH, VV)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    ax = axes[0]
    cs = ax.contour(VV, HH / 1e3, PS, levels=[0, 5, 10, 20, 30, 40, 50], colors="C0")
    ax.clabel(cs, fmt="%g m/s")
    ce = ax.contour(VV, HH / 1e3, energy_height(HH, VV), levels=np.arange(2000, 14001, 2000), colors="0.6", linestyles=":")
    ax.clabel(ce, fmt="%g m")
    ax.plot(Vr, hr / 1e3, "C3", lw=2, label="minimum-time path")
    vs = [stall_speed(ac, hh) for hh in hg_]
    ax.plot(vs, hg_ / 1e3, "k--", label="stall")
    ax.plot([ac.M_max * speed_of_sound(hh) for hh in hg_], hg_ / 1e3, "k-.", label=f"M = {ac.M_max}")
    ax.set_xlabel("true airspeed (m/s)"); ax.set_ylabel("altitude (km)"); ax.set_title("P_s contours (blue), h_E (dotted)")
    ax.legend(loc="upper left")
    ax = axes[1]
    ax.plot(fl.t, fl.h - h0, label="h - h0 (m), no thrust/drag, ISA")
    ax.plot(fl.t, (hE - hE[0]) * 1e6, label="h_E drift (micrometres)")
    ax.plot(fl2.t, fl2.h - h0, label="h - h0 (m), with thrust and drag")
    ax.set_xlabel("t (s)"); ax.set_title("phugoid: energy exchange"); ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "energy_height.png", dpi=120)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'energy_height.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
