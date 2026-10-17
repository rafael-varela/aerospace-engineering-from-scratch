"""P08.07 Pipe Flow & Losses: reference solution.

What pressure must a pressure-fed rocket's fuel tank hold to push the propellant through its feed
line and injector into the chamber, and how large must the line be? This module builds the answer
from the energy equation for a pipe with losses:

* the Darcy-Weisbach friction loss dp = f (L/D) rho V^2 / 2, with the laminar friction factor
  f = 64/Re (from the exact Hagen-Poiseuille solution of P08.06);
* the Colebrook equation for turbulent flow, solved by Newton's method on x = 1/sqrt(f), checked
  against the explicit Haaland formula, the smooth-pipe and fully rough limits and the Blasius
  correlation;
* minor losses dp = K rho V^2 / 2 for fittings (typical loss coefficients, labelled typical);
* an injector modelled as an orifice, dp = (mdot / (Cd A))^2 / (2 rho);
* a tank-to-injector pressure budget for a line made of several segments, with the hydrostatic head
  of the propellant column under the vehicle's acceleration;
* a root-find for the smallest line diameter that meets a loss limit.

The Colebrook and Haaland formulas are correlations fitted to measurements (Moody's chart), not
exact results; in the transition range 2300 < Re < 4000 the friction factor is genuinely uncertain
and is only interpolated here.

Units: SI (kg, m, s, Pa).

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes feed_line.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-08/lesson-07-pipe-flow-losses
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

G0 = 9.80665
RE_LAMINAR = 2300.0           # below this, laminar (64/Re)
RE_TURBULENT = 4000.0         # above this, Colebrook

# Typical propellant properties (order-of-magnitude design values; check a data sheet for design):
RHO_RP1, MU_RP1 = 810.0, 2.0e-3          # RP-1 kerosene near 20 degC: kg/m^3, Pa s (typical)
RHO_LOX, MU_LOX = 1141.0, 1.95e-4        # liquid oxygen near 90 K: kg/m^3, Pa s (typical)
ROUGHNESS_DRAWN = 1.5e-6                 # m, drawn stainless tubing (typical)

# Typical minor-loss coefficients K (handbook values; they vary with the fitting's make and size)
K_TYPICAL = {
    "tank outlet, sharp-edged": 0.5,
    "tank outlet, rounded": 0.05,
    "90 deg bend, r/D = 1.5": 0.3,
    "ball valve, fully open": 0.05,
    "swing check valve": 2.0,
    "filter element": 1.5,
    "exit into manifold": 1.0,
}


# ---------------------------------------------------------------- friction factor

def reynolds(rho: float, V: float, D: float, mu: float) -> float:
    """Pipe Reynolds number Re = rho V D / mu (dimensionless)."""
    return rho * V * D / mu


def friction_laminar(Re: float) -> float:
    """Darcy friction factor of fully developed laminar flow, f = 64 / Re (Hagen-Poiseuille)."""
    return 64.0 / Re


def friction_haaland(Re: float, rel_rough: float) -> float:
    """Haaland's explicit approximation to Colebrook:
    1/sqrt(f) = -1.8 log10[(rel_rough/3.7)^1.11 + 6.9/Re]."""
    return (-1.8 * np.log10((rel_rough / 3.7) ** 1.11 + 6.9 / Re)) ** -2


def colebrook_residual(f: float, Re: float, rel_rough: float) -> float:
    """Residual of the Colebrook equation, 1/sqrt(f) + 2 log10(rel_rough/3.7 + 2.51/(Re sqrt(f)))."""
    return 1.0 / np.sqrt(f) + 2.0 * np.log10(rel_rough / 3.7 + 2.51 / (Re * np.sqrt(f)))


def friction_colebrook(Re: float, rel_rough: float, x0: float | None = None, tol: float = 1e-13,
                       maxiter: int = 50) -> tuple[float, np.ndarray]:
    """Solve Colebrook for the Darcy friction factor by Newton's method on x = 1/sqrt(f):

        g(x) = x + 2 log10(a + b x) = 0,  a = rel_rough/3.7,  b = 2.51/Re,
        g'(x) = 1 + 2 b / ((a + b x) ln 10),

    starting from Haaland (or from ``x0``), until the step is below tol * x. g is increasing and
    concave in x, so Newton converges from any x > 0. Returns (f, iterates of x including x0)."""
    a, b = rel_rough / 3.7, 2.51 / Re
    x = 1.0 / np.sqrt(friction_haaland(Re, rel_rough)) if x0 is None else float(x0)
    hist = [x]
    for _ in range(maxiter):
        g = x + 2.0 * np.log10(a + b * x)
        dg = 1.0 + 2.0 * b / ((a + b * x) * np.log(10.0))
        step = g / dg
        x = x - step
        hist.append(x)
        if abs(step) <= tol * abs(x):
            break
    else:
        raise RuntimeError("Colebrook Newton iteration did not converge")
    return 1.0 / x ** 2, np.array(hist)


def friction_blasius(Re: float) -> float:
    """Blasius's smooth-pipe correlation f = 0.3164 Re^-0.25, for roughly 4000 < Re < 1e5."""
    return 0.3164 * Re ** -0.25


def friction_factor(Re: float, rel_rough: float) -> float:
    """Darcy friction factor for any Re: 64/Re below 2300, Colebrook above 4000, and a linear
    interpolation in Re between the two in the transition range (an engineering patch, not physics:
    the real value there depends on disturbances and may be either)."""
    if Re <= 0:
        raise ValueError("Re must be positive")
    if Re <= RE_LAMINAR:
        return friction_laminar(Re)
    if Re >= RE_TURBULENT:
        return friction_colebrook(Re, rel_rough)[0]
    fl = friction_laminar(RE_LAMINAR)
    ft = friction_colebrook(RE_TURBULENT, rel_rough)[0]
    s = (Re - RE_LAMINAR) / (RE_TURBULENT - RE_LAMINAR)
    return (1 - s) * fl + s * ft


# ---------------------------------------------------------------- losses

def darcy_weisbach_dp(f: float, L: float, D: float, rho: float, V: float) -> float:
    """Friction pressure drop of a straight pipe, dp = f (L/D) rho V^2 / 2 (Pa)."""
    return f * L / D * 0.5 * rho * V ** 2


def minor_loss_dp(K: float, rho: float, V: float) -> float:
    """Pressure drop of fittings with total loss coefficient K, dp = K rho V^2 / 2 (Pa)."""
    return K * 0.5 * rho * V ** 2


def segment_loss(mdot: float, rho: float, mu: float, D: float, L: float, roughness: float,
                 K: float = 0.0) -> dict[str, float]:
    """Velocity, Reynolds number, friction factor, friction loss and minor loss of one line segment
    of diameter D and length L carrying the mass flow mdot."""
    A = np.pi * D ** 2 / 4.0
    V = mdot / (rho * A)
    Re = reynolds(rho, V, D, mu)
    f = friction_factor(Re, roughness / D)
    return {"V": V, "Re": Re, "f": f, "dp_friction": darcy_weisbach_dp(f, L, D, rho, V),
            "dp_minor": minor_loss_dp(K, rho, V)}


def injector_dp(mdot: float, rho: float, Cd: float, A_inj: float) -> float:
    """Pressure drop across an injector treated as an orifice: dp = (mdot / (Cd A_inj))^2 / (2 rho)."""
    return (mdot / (Cd * A_inj)) ** 2 / (2.0 * rho)


def injector_area_for_dp(mdot: float, rho: float, Cd: float, dp: float) -> float:
    """Total injector orifice area giving the design pressure drop dp at mdot."""
    return mdot / (Cd * np.sqrt(2.0 * rho * dp))


def feed_line_budget(mdot: float, rho: float, mu: float, segments: list[dict], p_chamber: float,
                     Cd: float, A_inj: float, head: float = 0.0, accel_g: float = 1.0) -> dict:
    """Tank pressure needed to deliver mdot into a chamber at p_chamber:

        p_tank = p_chamber + dp_injector + sum(dp_friction + dp_minor) - rho * accel_g * g0 * head,

    where each segment is a dict with keys D, L, roughness and K (total minor-loss coefficient), and
    ``head`` is the height of the propellant column above the injector, under an axial acceleration of
    ``accel_g`` g (the column pushes, so it reduces the tank pressure needed). The tank's ullage is
    taken at rest, so the dynamic pressure at the line exit appears through the exit loss (K = 1) in
    the last segment. Returns the terms, the per-segment results and p_tank."""
    segs = [segment_loss(mdot, rho, mu, s["D"], s["L"], s["roughness"], s.get("K", 0.0)) for s in segments]
    dp_f = sum(s["dp_friction"] for s in segs)
    dp_m = sum(s["dp_minor"] for s in segs)
    dp_i = injector_dp(mdot, rho, Cd, A_inj)
    dp_h = rho * accel_g * G0 * head
    return {"segments": segs, "dp_friction": dp_f, "dp_minor": dp_m, "dp_injector": dp_i,
            "dp_head": dp_h, "p_chamber": p_chamber,
            "p_tank": p_chamber + dp_i + dp_f + dp_m - dp_h}


def line_diameter_for_loss(mdot: float, rho: float, mu: float, L: float, roughness: float, K: float,
                           dp_max: float, D_lo: float = 1e-3, D_hi: float = 0.5) -> float:
    """Smallest diameter of a single line (length L, roughness, total K) whose friction plus minor
    loss at mdot equals dp_max, by Brent's method (the loss falls monotonically with D)."""
    def excess(D: float) -> float:
        s = segment_loss(mdot, rho, mu, D, L, roughness, K)
        return s["dp_friction"] + s["dp_minor"] - dp_max
    return brentq(excess, D_lo, D_hi, xtol=1e-12, rtol=1e-12)


# ---------------------------------------------------------------- the stage

def fuel_line(D: float = 0.020) -> list[dict]:
    """The lesson's fuel line: a sharp-edged tank outlet, 2.5 m of drawn tubing of diameter D with
    three bends, a ball valve, a filter and a check valve, ending in the injector manifold."""
    K = (K_TYPICAL["tank outlet, sharp-edged"] + 3 * K_TYPICAL["90 deg bend, r/D = 1.5"]
         + K_TYPICAL["ball valve, fully open"] + K_TYPICAL["filter element"]
         + K_TYPICAL["swing check valve"] + K_TYPICAL["exit into manifold"])
    return [{"D": D, "L": 2.5, "roughness": ROUGHNESS_DRAWN, "K": K}]


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save feed_line.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    mdot, pc, dp_inj_design, Cd = 1.6, 2.0e6, 0.4e6, 0.7
    head, accel = 1.0, 3.0
    A_inj = injector_area_for_dp(mdot, RHO_RP1, Cd, dp_inj_design)
    line = fuel_line()
    b = feed_line_budget(mdot, RHO_RP1, MU_RP1, line, pc, Cd, A_inj, head, accel)
    s = b["segments"][0]
    print(f"RP-1 (typical): rho = {RHO_RP1} kg/m^3, mu = {MU_RP1} Pa s; mdot = {mdot} kg/s, chamber {pc / 1e6} MPa")
    print(f"line D = {line[0]['D'] * 1e3:.0f} mm, L = {line[0]['L']} m, roughness {ROUGHNESS_DRAWN * 1e6} um, "
          f"total K = {line[0]['K']:.2f}")
    print(f"  V = {s['V']:.3f} m/s, Re = {s['Re']:.0f}, f (Colebrook) = {s['f']:.5f}, "
          f"Haaland {friction_haaland(s['Re'], ROUGHNESS_DRAWN / line[0]['D']):.5f}, "
          f"Blasius {friction_blasius(s['Re']):.5f}")
    f_n, hist = friction_colebrook(s["Re"], ROUGHNESS_DRAWN / line[0]["D"])
    print(f"  Newton from Haaland: {hist.size - 1} iterations, residual {colebrook_residual(f_n, s['Re'], ROUGHNESS_DRAWN / line[0]['D']):.1e}")
    print(f"  dynamic pressure rho V^2/2 = {0.5 * RHO_RP1 * s['V'] ** 2 / 1e3:.2f} kPa")
    print(f"budget (kPa): chamber {pc / 1e3:.1f} + injector {b['dp_injector'] / 1e3:.1f} (Cd {Cd}, "
          f"A = {A_inj * 1e6:.2f} mm^2) + friction {b['dp_friction'] / 1e3:.2f} + minor {b['dp_minor'] / 1e3:.2f} "
          f"- head {b['dp_head'] / 1e3:.2f} ({head} m at {accel} g) = tank {b['p_tank'] / 1e3:.1f} kPa")
    line_loss = b["dp_friction"] + b["dp_minor"]
    print(f"  line loss {line_loss / 1e3:.2f} kPa = {line_loss / pc * 100:.2f} % of chamber pressure")
    D_min = line_diameter_for_loss(mdot, RHO_RP1, MU_RP1, 2.5, ROUGHNESS_DRAWN, line[0]["K"], 0.1 * pc)
    print(f"  smallest line for a loss of 10 % of p_c: D = {D_min * 1e3:.2f} mm")
    print("\ndiameter sweep: D (mm), V (m/s), Re, f, friction (kPa), minor (kPa), tank (MPa)")
    Ds = np.array([10, 12, 16, 20, 25, 32]) * 1e-3
    tanks, losses = [], []
    for D in Ds:
        bb = feed_line_budget(mdot, RHO_RP1, MU_RP1, fuel_line(D), pc, Cd, A_inj, head, accel)
        ss = bb["segments"][0]
        tanks.append(bb["p_tank"])
        losses.append(bb["dp_friction"] + bb["dp_minor"])
        print(f"  {D * 1e3:5.0f}  {ss['V']:7.3f}  {ss['Re']:8.0f}  {ss['f']:.5f}  {bb['dp_friction'] / 1e3:9.2f}  "
              f"{bb['dp_minor'] / 1e3:9.2f}  {bb['p_tank'] / 1e6:.4f}")
    slope = np.polyfit(np.log(Ds[2:]), np.log(np.array(losses)[2:]), 1)[0]
    print(f"  log-log slope of line loss against D (16-32 mm): {slope:.2f} (V^2 alone gives -4; f(Re) adds the rest)")

    # Colebrook against Haaland across the Moody chart
    Res = np.logspace(np.log10(4e3), 8, 60)
    rels = [0.0, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2]
    worst = max(abs(friction_haaland(Re, e) / friction_colebrook(Re, e)[0] - 1) for Re in Res for e in rels)
    iters = max(friction_colebrook(Re, e)[1].size - 1 for Re in Res for e in rels)
    print(f"\nHaaland vs Colebrook over 4e3 < Re < 1e8, 0 <= eps/D <= 0.05: max |difference| {worst * 100:.2f} %; "
          f"Newton iterations from Haaland <= {iters}")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    Re_l = np.logspace(np.log10(600), np.log10(RE_LAMINAR), 30)
    ax[0].loglog(Re_l, friction_laminar(Re_l), "k-", label="laminar 64/Re")
    for e in rels[1:]:
        ax[0].loglog(Res, [friction_colebrook(Re, e)[0] for Re in Res], label=f"Colebrook eps/D = {e:g}")
    ax[0].loglog(Res, [friction_colebrook(Re, 0.0)[0] for Re in Res], "k--", label="smooth")
    ax[0].loglog(Res[Res < 1e5], friction_blasius(Res[Res < 1e5]), "k:", label="Blasius (smooth)")
    ax[0].set_xlabel("Re")
    ax[0].set_ylabel("Darcy friction factor f")
    ax[0].set_title("Moody chart from the correlations")
    ax[0].legend(fontsize=7)
    ax[1].plot(Ds * 1e3, np.array(tanks) / 1e6, "o-", label="tank pressure needed")
    ax[1].axhline((pc + dp_inj_design - RHO_RP1 * accel * G0 * head) / 1e6, color="k", lw=0.8, ls="--",
                  label="chamber + injector - head (no line loss)")
    ax[1].set_xlabel("line diameter (mm)")
    ax[1].set_ylabel("p_tank (MPa)")
    ax[1].set_title(f"RP-1 feed line, {mdot} kg/s")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "feed_line.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'feed_line.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
