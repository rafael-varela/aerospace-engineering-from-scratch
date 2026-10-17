"""P08.03 Control Volume: Mass & Momentum: starter code.

How much thrust does a small turbojet give on the test stand and in cruise, how much of it is
pressure thrust, and does the answer depend on where the control volume is drawn? This module
applies the integral (control-volume) forms of conservation of mass and of linear momentum, from
the Reynolds transport theorem, to steady flows through a fixed or uniformly moving control volume:

    mass:      d/dt int rho dV + sum_faces rho u_n A = 0
    momentum:  sum F_x (pressure on the surface + forces through cut solids) = d/dt int rho u_x dV
               + sum_faces rho u_n u_x A

with u_n the velocity component along the outward normal, relative to the control surface. The
stress model is the simplest one: on every face that the fluid crosses, only the pressure acts
(viscous normal stresses neglected, faces drawn where the flow is uniform or its profile known).

The tool is a control volume written as a list of flat faces, each carrying uniform properties
(``Face``), with functions that sum the mass flux and the x-momentum flux plus pressure force.
The independent checks: a box drawn around a turbojet, of any size, whose mass balance closes to
round-off and whose momentum balance returns the same thrust as the textbook formula
F = mdot_e V_e - mdot_a V_0 + (p_e - p_a) A_e; momentum-flux correction factors of non-uniform
jet profiles by quadrature against their closed forms; the rocket equation integrated by RK4
against Tsiolkovsky with gravity loss; and a moving vane whose power peaks at one third of the jet
speed.

Units: SI (kg, m, s, N, Pa, W).

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes jet_thrust.png)
Test: pytest curriculum/phase-08/lesson-03-control-volume-mass-momentum
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

G0 = 9.80665
R_AIR = 287.05287


@dataclass(frozen=True)
class Face:
    """A flat piece of control surface with uniform properties.

    area (m^2); nx: x-component of the outward unit normal; rho (kg/m^3); un: velocity component
    along the outward normal, relative to the face (m/s, positive out); ux: x-velocity of the fluid
    (m/s, in the frame of the momentum balance); p: pressure (Pa)."""
    area: float
    nx: float
    rho: float
    un: float
    ux: float
    p: float


# ---------------------------------------------------------------- the control-volume sums

def mass_outflow(faces: list[Face]) -> float:
    """Net mass flow out of the control surface, sum of rho u_n A (kg/s). Zero for a steady flow."""
    raise NotImplementedError


def momentum_outflow_x(faces: list[Face]) -> float:
    """Net x-momentum flow out of the control surface, sum of rho u_n u_x A (N)."""
    raise NotImplementedError


def pressure_force_x(faces: list[Face]) -> float:
    """x-force of the surrounding fluid's pressure on the control-volume contents, -sum p n_x A (N)."""
    raise NotImplementedError


def solid_force_x(faces: list[Face]) -> float:
    """x-force that solids cut by the control surface (a pylon, a test-stand load cell) must apply
    to the contents to keep a steady flow steady: momentum outflow minus the pressure force,
    F = sum rho u_n u_x A + sum p n_x A."""
    raise NotImplementedError


# ---------------------------------------------------------------- the turbojet

def exit_area(mdot_e: float, Ve: float, pe: float, Te: float, R: float = R_AIR) -> float:
    """Nozzle exit area (m^2) from continuity, A_e = mdot_e / (rho_e V_e), with rho_e = p_e/(R T_e)."""
    raise NotImplementedError


def turbojet_thrust(mdot_a: float, f: float, V0: float, Ve: float, pe: float, pa: float, Ae: float) -> dict[str, float]:
    """Net thrust of a jet engine from the momentum balance (textbook formula):

        F = (1 + f) mdot_a V_e - mdot_a V_0 + (p_e - p_a) A_e.

    Returns the net thrust 'F', the 'gross' momentum thrust (1+f) mdot_a V_e, the 'ram_drag'
    mdot_a V_0 and the 'pressure' thrust (p_e - p_a) A_e, all in N."""
    raise NotImplementedError


def turbojet_control_volume(mdot_a: float, f: float, V0: float, Ve: float, pe: float, Te: float,
                            pa: float, rho_a: float, A_box: float) -> list[Face]:
    """A cylindrical control volume of cross-section A_box around the engine, in the engine's frame:
    air approaches at V0 along +x (the aeroplane flies towards -x). Faces:

    * inlet disc: area A_box, n_x = -1, ambient air at V0, p_a;
    * exit disc, jet core: area A_e from continuity, n_x = +1, exhaust at V_e, p_e;
    * exit disc, outer part: area A_box - A_e, n_x = +1, ambient air at V0, p_a;
    * side wall: n_x = 0, carrying whatever mass the other faces leave unbalanced (the stream
      tube around the engine spills sideways), with x-velocity V0;
    * fuel line: n_x = 0, fuel entering at zero x-velocity (the tanks are outside the box).

    The side face's flow follows from conservation of mass; it is the face that makes the box
    independent of its size."""
    raise NotImplementedError


# ---------------------------------------------------------------- non-uniform profiles

def profile_fluxes(u_of_r: Callable[[np.ndarray], np.ndarray], R: float, rho: float, n: int = 200) -> dict[str, float]:
    """Mass flow, momentum flow and the momentum-flux correction factor of an axisymmetric velocity
    profile u(r) on a circular face of radius R, by Gauss-Legendre quadrature in r on n points:

        mdot = int rho u 2 pi r dr,  J = int rho u^2 2 pi r dr,  beta = J A / (rho Q^2)

    with Q the volume flow and A = pi R^2. beta = 1 for a uniform profile and > 1 otherwise: using
    mdot * V_mean underestimates the momentum a non-uniform jet carries."""
    raise NotImplementedError


def beta_laminar() -> float:
    """Momentum-flux correction factor of the laminar pipe profile u = u_max (1 - r^2/R^2): 4/3."""
    raise NotImplementedError


def beta_power_law(n: float) -> float:
    """Momentum-flux correction factor of the power-law profile u = u_max (1 - r/R)^(1/n):
    beta = (n + 1)(2n + 1)^2 / (4 n^2 (n + 2)); 1.020 for the 1/7th law."""
    raise NotImplementedError


# ---------------------------------------------------------------- rockets and moving control volumes

def rocket_thrust(mdot: float, Ve: float, pe: float, pa: float, Ae: float) -> float:
    """Rocket thrust from the momentum balance of a control volume moving with the vehicle (no air
    enters): F = mdot V_e + (p_e - p_a) A_e."""
    raise NotImplementedError


def rocket_burn(m0: float, mdot: float, F: float, t_burn: float, n: int, g: float = G0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Vertical flight with constant thrust F and mass flow mdot, no drag, uniform g, integrated by
    RK4 in n steps: dv/dt = F/m - g, dh/dt = v, m = m0 - mdot t. Returns (t, v, h)."""
    raise NotImplementedError


def tsiolkovsky_delta_v(c: float, m0: float, mf: float, t_burn: float = 0.0, g: float = G0) -> float:
    """Velocity gained in a vertical burn with effective exhaust velocity c = F/mdot:
    dv = c ln(m0/mf) - g t_burn (the second term is the gravity loss)."""
    raise NotImplementedError


def vane_force(rho: float, A: float, V: float, u_vane: float, theta_deg: float) -> tuple[float, float]:
    """Force (Fx, Fy) in N that a jet of area A and speed V exerts on a single vane moving
    downstream at u_vane and turning the jet through theta, from a control volume moving with the
    vane (steady in that frame; relative speed V - u_vane, unchanged in magnitude across the vane):
    Fx = rho A (V - u)^2 (1 - cos theta), Fy = -rho A (V - u)^2 sin theta."""
    raise NotImplementedError


def thrust_reverser_fraction(theta_deg: float) -> float:
    """Fraction of the gross jet thrust that a reverser turning the jet through theta delivers as
    net reverse force (static, same exit speed): -cos(theta). 135 deg gives 0.707."""
    raise NotImplementedError


def observed_order(errors: np.ndarray, ratio: float = 2.0) -> np.ndarray:
    """Observed convergence orders log(e_k / e_{k+1}) / log(ratio) for steps refined by ``ratio``."""
    e = np.asarray(errors, dtype=float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


# ---------------------------------------------------------------- the lesson's engine (illustrative numbers)

STATIC = dict(mdot_a=20.0, f=0.02, V0=0.0, Ve=600.0, pe=101325.0, Te=800.0, pa=101325.0, rho_a=1.225)
CRUISE = dict(mdot_a=10.0, f=0.02, V0=236.0, Ve=650.0, pe=40000.0, Te=900.0, pa=22632.0, rho_a=0.36392)


def engine_case(case: dict[str, float]) -> dict[str, float]:
    """Thrust breakdown of an engine case dict (keys as STATIC / CRUISE) with A_e from continuity."""
    Ae = exit_area((1 + case["f"]) * case["mdot_a"], case["Ve"], case["pe"], case["Te"])
    out = turbojet_thrust(case["mdot_a"], case["f"], case["V0"], case["Ve"], case["pe"], case["pa"], Ae)
    out["Ae"] = Ae
    return out


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save jet_thrust.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    print("illustrative turbojet (numbers chosen for the lesson, not a specific engine)")
    for name, case in (("static, sea level", STATIC), ("cruise, 11 km, 236 m/s", CRUISE)):
        r = engine_case(case)
        print(f"  {name:24s} A_e = {r['Ae']:.4f} m^2  gross {r['gross']:8.1f} N  ram drag {r['ram_drag']:7.1f} N  "
              f"pressure {r['pressure']:7.1f} N  net F = {r['F']:8.1f} N")
    c = CRUISE
    print("\ncontrol-volume box around the cruise engine, any size:")
    print("   A_box (m^2)   mass residual (kg/s)   F from box (N)   F - formula (N)")
    F_formula = engine_case(c)["F"]
    for A_box in (0.5, 2.0, 10.0, 100.0):
        faces = turbojet_control_volume(c["mdot_a"], c["f"], c["V0"], c["Ve"], c["pe"], c["Te"], c["pa"], c["rho_a"], A_box)
        Fb = solid_force_x(faces)
        print(f"   {A_box:9.1f}     {mass_outflow(faces):+.2e}           {Fb:9.3f}      {Fb - F_formula:+.1e}")

    print("\nmomentum-flux correction factors (quadrature vs closed form):")
    lam = profile_fluxes(lambda r: 2.0 * (1 - r ** 2), 1.0, 1.0)
    print(f"  laminar parabola: beta = {lam['beta']:.6f} (4/3 = {beta_laminar():.6f})")
    for nn in (6.0, 7.0, 10.0):
        pw = profile_fluxes(lambda r: (1 - r) ** (1 / nn), 1.0, 1.0, n=200)
        print(f"  power law 1/{nn:.0f}: beta = {pw['beta']:.6f} (closed form {beta_power_law(nn):.6f}; "
              f"the profile's infinite slope at the wall limits Gauss quadrature)")

    print("\nsounding rocket (illustrative): m0 = 100 kg, 60 kg propellant, 8 s burn, c = 2000 m/s")
    m0, mp, tb, cc = 100.0, 60.0, 8.0, 2000.0
    md = mp / tb
    F = md * cc
    dv = tsiolkovsky_delta_v(cc, m0, m0 - mp, tb)
    print(f"  thrust {F:.0f} N, initial thrust/weight {F / (m0 * G0):.2f}")
    errs = []
    for n in (10, 20, 40, 80):
        t, v, h = rocket_burn(m0, md, F, tb, n)
        errs.append(abs(v[-1] - dv))
    print(f"  burnout speed: Tsiolkovsky - g t_b = {dv:.3f} m/s (gravity loss {G0 * tb:.2f} m/s); RK4 n = 80: {v[-1]:.6f} m/s")
    print(f"  RK4 errors {np.array2string(np.array(errs), precision=2)}, orders {np.array2string(observed_order(errs), precision=2)}")
    print(f"  burnout altitude {h[-1]:.1f} m")

    rho_w, Aj, Vj = 998.2, 1e-4, 20.0
    us = np.linspace(0, Vj, 201)
    P = np.array([vane_force(rho_w, Aj, Vj, u, 180.0)[0] * u for u in us])
    print(f"\nvane turning a water jet (A = 1 cm^2, V = 20 m/s) by 180 deg: power peaks at u = {us[np.argmax(P)]:.2f} m/s "
          f"(V/3 = {Vj / 3:.2f}), {P.max():.1f} W")
    print(f"thrust reverser turning the jet by 135 deg: {thrust_reverser_fraction(135.0):.3f} of gross thrust reversed")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    V0s = np.linspace(0, 300, 61)
    Fs = [turbojet_thrust(20.0, 0.02, V, 600.0, 101325.0, 101325.0, 0.0)["F"] for V in V0s]
    ax[0].plot(V0s, np.array(Fs) / 1e3, label="fixed mdot and V_e")
    ax[0].plot(V0s, 20.0 * V0s / 1e3, "--", label="ram drag mdot V0")
    ax[0].set_xlabel("flight speed V0 (m/s)")
    ax[0].set_ylabel("force (kN)")
    ax[0].set_title("net thrust falls with ram drag")
    ax[0].legend(fontsize=8)
    t, v, h = rocket_burn(m0, md, F, tb, 200)
    ax[1].plot(t, v, label="RK4")
    tt = np.linspace(0.01, tb, 50)
    ax[1].plot(tt, [tsiolkovsky_delta_v(cc, m0, m0 - md * x, x) for x in tt], "o", ms=3, mfc="none", label="Tsiolkovsky - g t")
    ax[1].set_xlabel("t (s)")
    ax[1].set_ylabel("v (m/s)")
    ax[1].set_title("sounding rocket burn")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "jet_thrust.png", dpi=110)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'jet_thrust.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
