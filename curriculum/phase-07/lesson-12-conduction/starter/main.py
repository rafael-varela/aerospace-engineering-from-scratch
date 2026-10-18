"""P07.12 Heat Transfer: Conduction: starter code.

Two aerospace conduction questions, answered with a finite-volume solver checked against closed forms:

* steady, two-dimensional: how hot does a processor get on a conduction-cooled avionics circuit card
  clamped along two edges to a cold chassis wall, and does the card need an aluminium heat frame?
  (Laplace/Poisson equation for the in-plane temperature of a thin sheet, k t (T_xx + T_yy) + q'' = 0,
  with optional heat loss from the faces, h_f (T - T_inf), which turns the sheet into a 2D fin);
* transient, one-dimensional: how many seconds can an uncooled ("heat-sink") rocket test chamber wall
  fire before its hot-gas surface reaches its temperature limit, and when does the semi-infinite
  (erfc) solution stop describing a wall of finite thickness?

The numerics are vertex-centred finite volumes: every node owns a control volume (half at an edge,
a quarter at a corner), neighbours exchange heat through conductances G = k t (face length)/(spacing),
and each node's energy balance is one row of K T = f. Steady: solve K T = f with fixed-temperature
rows replaced. Transient: C dT/dt = f - K T, marched with the theta method (0 explicit, 1/2
Crank-Nicolson, 1 implicit). Because the rows are energy balances, the heat leaving through every
boundary adds up to the heat generated to round-off.

Checks: the Fourier-series solution of a rectangle with one heated edge; the 1D wall and fin closed
forms; second-order convergence; the global energy balance; the semi-infinite erfc solutions with a
fixed surface temperature and with convection (Lienhard eq. 5.53).

Units: SI (m, s, K, W). Material data: Lienhard & Lienhard, A Heat Transfer Textbook, 6th ed.
(v6.00, 2024), Appendix A, Tables A.1 (metals) and A.2 (non-metals), properties at 20 C.

Run:  python main.py [out_dir]   (prints the lesson's numbers, writes conduction.png)
Test: pytest curriculum/phase-07/lesson-12-conduction
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.special import erf, erfc

# Lienhard & Lienhard, A Heat Transfer Textbook, 6th ed. (2024), Appendix A, properties at 20 C:
# density rho (kg/m^3), specific heat cp (J/(kg K)), conductivity k (W/(m K)).
MATERIALS = {
    "copper":         {"rho": 8954.0, "cp": 384.0, "k": 398.0},    # Table A.1, pure copper
    "aluminium 6061": {"rho": 2700.0, "cp": 896.0, "k": 167.0},    # Table A.1, alloy 6061-T6
    "Inconel X-750":  {"rho": 8510.0, "cp": 442.0, "k": 11.6},     # Table A.1, nickel alloy
    "Ti-6Al-4V":      {"rho": 4430.0, "cp": 580.0, "k": 7.1},      # Table A.1, titanium alloy
    "stainless 304":  {"rho": 8000.0, "cp": 400.0, "k": 13.8},     # Table A.1
    "FR4":            {"rho": 1800.0, "cp": 1600.0, "k": 0.29},    # Table A.2, epoxy/glass-cloth laminate (G-10, FR4); cp approximate
    "fused silica":   {"rho": 2200.0, "cp": 745.0, "k": 1.38},     # Table A.2, silicon dioxide, fused silica glass (27 C)
}


def diffusivity(material: str) -> float:
    """Thermal diffusivity alpha = k / (rho cp) [m^2/s]."""
    m = MATERIALS[material]
    return m["k"] / (m["rho"] * m["cp"])


# ---------------------------------------------------------------- closed forms

def wall_heat_flux(layers: list[tuple[float, float]], T_hot: float, T_cold: float,
                   h_hot: float = np.inf, h_cold: float = np.inf) -> tuple[float, list[float]]:
    """Steady 1D heat flux [W/m^2] through a plane wall of layers [(thickness m, k W/(m K)), ...] in
    series, with optional convection h [W/(m^2 K)] on each side (np.inf = surface at the fluid
    temperature). Thermal resistances per unit area add: R'' = 1/h_hot + sum L/k + 1/h_cold, and
    q'' = (T_hot - T_cold) / R''. Returns (q'', interface temperatures from the hot surface to the cold)."""
    raise NotImplementedError


def fin_closed_form(x: np.ndarray, L: float, m: float, theta_b: float) -> np.ndarray:
    """Temperature excess theta = T - T_inf of a straight fin with an adiabatic tip,
    theta'' = m^2 theta, theta(0) = theta_b, theta'(L) = 0:  theta = theta_b cosh(m (L - x)) / cosh(m L).
    For a sheet of thickness t losing heat from both faces, m^2 = 2 h / (k t)."""
    raise NotImplementedError


def fin_heat_rate(kA: float, L: float, m: float, theta_b: float) -> float:
    """Heat entering the base of an adiabatic-tip fin [W]: k A_c m theta_b tanh(m L), where kA = k A_c."""
    raise NotImplementedError


def rectangle_fourier(x: np.ndarray, y: np.ndarray, W: float, H: float, T_top: float = 1.0,
                      n_terms: int = 401) -> np.ndarray:
    """Laplace's equation on 0 < x < W, 0 < y < H with T = 0 on the left, right and bottom edges and
    T = T_top on the top edge (y = H), by separation of variables:
    T = (4 T_top / pi) sum_{n odd} sin(n pi x / W) sinh(n pi y / W) / (n sinh(n pi H / W)).
    The ratio of sinh is evaluated as exp(-n pi (H - y)/W)(1 - e^{-2 n pi y/W})/(1 - e^{-2 n pi H/W})
    so that large n do not overflow. Converges slowly near the two top corners (a jump there)."""
    raise NotImplementedError


def semi_infinite_step(x: np.ndarray, t: float, alpha: float, T_i: float, T_s: float) -> np.ndarray:
    """Semi-infinite solid at T_i whose surface is held at T_s from t = 0:
    (T - T_s) / (T_i - T_s) = erf(x / (2 sqrt(alpha t)))."""
    raise NotImplementedError


def semi_infinite_step_energy(t: float, k: float, alpha: float, T_i: float, T_s: float) -> float:
    """Heat absorbed per unit area by that solid up to time t [J/m^2]: the integral of the surface
    flux k (T_s - T_i) / sqrt(pi alpha t), i.e. 2 k (T_s - T_i) sqrt(t / (pi alpha))."""
    raise NotImplementedError


def semi_infinite_convection(x: np.ndarray, t: float, alpha: float, k: float, h: float, T_i: float,
                             T_inf: float) -> np.ndarray:
    """Semi-infinite solid at T_i, its surface exposed from t = 0 to a fluid at T_inf with coefficient h
    (Lienhard eq. 5.53): Theta = (T - T_inf)/(T_i - T_inf) = erf(zeta/2) + exp(beta zeta + beta^2)
    erfc(zeta/2 + beta), zeta = x / sqrt(alpha t), beta = h sqrt(alpha t) / k. The product exp * erfc
    is evaluated with the scaled function erfcx to avoid overflow."""
    raise NotImplementedError


# ---------------------------------------------------------------- the finite-volume system

def conduction_system_2d(W: float, H: float, nx: int, ny: int, kt: float | Callable,
                         q_gen: float | Callable = 0.0, bcs: dict | None = None,
                         h_face: float = 0.0, T_face: float = 0.0, rho_c_t: float = 1.0) -> dict:
    """Assemble the vertex-centred finite-volume system for the in-plane temperature of a sheet,
    div(k t grad T) + q'' - h_face (T - T_face) = rho c t dT/dt, on [0, W] x [0, H] with nx x ny nodes
    (spacing dx = W/(nx-1), dy = H/(ny-1)).

    kt: sheet conductance k t [W/K] (a number, or a function kt(x, y) evaluated at face midpoints);
    q_gen: heat generated per unit sheet area [W/m^2] (a number, a function of (x, y) evaluated at the
           nodes, or an (ny, nx) array of control-volume averages);
    h_face: total face coefficient [W/(m^2 K)] to T_face (both faces together), the 2D-fin term;
    rho_c_t: heat capacity per unit area [J/(m^2 K)] (for transients).
    bcs: for each edge 'left' (x = 0), 'right' (x = W), 'bottom' (y = 0), 'top' (y = H) one of
         ('adiabatic',), ('T', value or function of the coordinate along the edge),
         ('flux', q into the sheet per unit edge length [W/m]), ('h', conductance per unit edge length
         [W/(m K)] (h times the sheet thickness), T_fluid). Missing edges are adiabatic.
         Fixed-temperature edges win at corners.

    Node i has the control volume V_i = (dx or dx/2)(dy or dy/2). Row i of K T = f is its energy
    balance: sum over neighbours G_ij (T_i - T_j) + h_face V_i T_i (+ edge conductances) = q'' V_i +
    h_face V_i T_face (+ edge inflows). K is symmetric. Returns {'K' (csr), 'f', 'C' (= rho_c_t V),
    'fixed' (bool mask), 'T_fixed', 'x', 'y', 'V', 'dx', 'dy', 'nx', 'ny', 'gen' (q'' V), 'h_face', 'T_face',
    'edge_q', 'edge_G', 'edge_Tf' (edge inflows and conductances per node, for the energy balance)}."""
    raise NotImplementedError


def conduction_system_1d(L: float, n: int, k: float | np.ndarray, rho_c: float | np.ndarray,
                         left: tuple, right: tuple) -> dict:
    """The same finite-volume system for a 1D slab 0 <= x <= L, per unit area (n nodes, dx = L/(n-1)):
    conductance k/dx between neighbours, capacity rho c times the control-volume width (half at the
    ends). k may be an array of n - 1 values, one per interval (a layered wall: put the interfaces on
    nodes); rho_c may be an array of n node values. Boundary conditions as in conduction_system_2d:
    ('adiabatic',), ('T', value), ('flux', q'' into the slab [W/m^2]), ('h', h [W/(m^2 K)], T_fluid).
    Returns the same keys (with 'x', and 'V' the control-volume widths)."""
    raise NotImplementedError


def solve_steady(sys_: dict) -> np.ndarray:
    """Steady temperatures: solve K T = f on the free nodes with the fixed nodes moved to the right-hand
    side (sparse direct solve). Returns T as a 2D array (ny, nx) for a 2D system, 1D otherwise."""
    raise NotImplementedError


def heat_balance(sys_: dict, T: np.ndarray) -> dict[str, float]:
    """Global energy balance of a steady solution [W] (per unit depth or area for 1D):
    'generated' = sum q'' V; 'in_edges' = flux edges plus convective edges G_e (T_fluid - T);
    'out_faces' = sum h_face V (T - T_face); 'out_fixed' = heat drawn out through the fixed-temperature
    nodes, minus the sum of their residuals (K T - f) (the residual is what each fixed node's energy
    balance lacks, which the reservoir holding it at T supplies). Because the free rows are satisfied
    and the conduction terms cancel in pairs (K symmetric with zero row sums), generated + in_edges =
    out_fixed + out_faces exactly: 'imbalance' is round-off. 'residual_free' is the largest free-row residual."""
    raise NotImplementedError


# ---------------------------------------------------------------- transients

def march(sys_: dict, T0: float | np.ndarray, dt: float, n_steps: int, theta: float = 1.0,
          stop: Callable[[float, np.ndarray], bool] | None = None) -> dict:
    """March C dT/dt = f - K T with the theta method on the free nodes (fixed nodes held at T_fixed
    from t = 0+):  (C/dt + theta K) T^{n+1} = (C/dt - (1 - theta) K) T^n + f  on the free rows.
    theta = 0 explicit (stable only for dt <= min C_i / K_ii, Fourier number alpha dt / dx^2 <= 1/2 in
    1D), 1/2 Crank-Nicolson (second order in dt), 1 implicit Euler (first order, unconditionally
    stable and free of oscillations). Energy bookkeeping, exact for the scheme: the energy stored in
    the free nodes changes by dt x (generation + edge inflows - face losses + heat from the fixed nodes),
    all evaluated at T^theta = theta T^{n+1} + (1 - theta) T^n.
    stop(t, T), if given, is called after every step and ends the march early when it returns True.
    Returns {'T' (final, flat), 't', 'E_stored' (J, free nodes, relative to t = 0), 'Q_in' (J)}."""
    raise NotImplementedError


def time_to_surface_limit(material: str, thickness: float, h: float, T_gas: float, T_init: float,
                          T_limit: float, n: int = 101, dt: float = 0.01, t_max: float = 120.0) -> tuple[float, dict]:
    """Uncooled ("heat-sink") wall: a slab of the given material and thickness, initially at T_init,
    hot-gas side (x = 0) convecting with coefficient h to the gas at T_gas, outside (x = thickness)
    adiabatic. Marches with Crank-Nicolson until the hot-side surface reaches T_limit and returns the
    time [s] (linear interpolation within the last step; np.inf if not reached by t_max) and the march
    result (with 'x' added)."""
    raise NotImplementedError


def semi_infinite_time_to_limit(material: str, h: float, T_gas: float, T_init: float, T_limit: float) -> float:
    """The same time [s] from the semi-infinite convective solution (x = 0 of semi_infinite_convection),
    by brentq: valid only while the heat has not reached the back of the wall (sqrt(alpha t) << thickness)."""
    raise NotImplementedError


def observed_order(errors: np.ndarray, ratio: float = 2.0) -> np.ndarray:
    """Observed convergence orders log(e_k / e_{k+1}) / log(ratio) for grids refined by `ratio`."""
    e = np.asarray(errors, dtype=float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


# ---------------------------------------------------------------- the avionics card

def board_sheet_conductance(t_board: float = 1.6e-3, n_planes: int = 4, t_cu: float = 35e-6,
                            t_frame: float = 0.0) -> float:
    """In-plane sheet conductance k t [W/K] of a circuit card: copper planes, FR4 laminate and an
    optional aluminium 6061 heat frame conduct side by side (parallel), so their k t add."""
    raise NotImplementedError


def card_heat_load(x: np.ndarray, y: np.ndarray, P_chip: float = 10.0, P_spread: float = 5.0,
                   chip: float = 0.030) -> np.ndarray:
    """Heat generated per unit card area [W/m^2], averaged over each node's control volume, on the grid
    x (nx), y (ny) of a card W = x[-1] by H = y[-1]: P_spread spread evenly over the card plus a
    processor dissipating P_chip on a chip x chip square at the centre. The chip's share of each
    control volume is its exact overlap area, so the total is P_chip + P_spread on any grid.
    (Typical, illustrative loads.) Returns an (ny, nx) array."""
    W, H = x[-1], y[-1]

    def overlap(c: np.ndarray, L: float, centre: float) -> np.ndarray:
        d = c[1] - c[0]
        lo, hi = np.maximum(c - d / 2, 0.0), np.minimum(c + d / 2, L)
        o = np.clip(np.minimum(hi, centre + chip / 2) - np.maximum(lo, centre - chip / 2), 0.0, None)
        return o / (hi - lo)

    frac = np.outer(overlap(y, H, H / 2), overlap(x, W, W / 2))
    return P_spread / (W * H) + P_chip / chip ** 2 * frac


def solve_card(kt: float, n: int = 81, T_wall: float = 343.15, P_chip: float = 10.0, P_spread: float = 5.0,
               W: float = 0.160, H: float = 0.100) -> tuple[np.ndarray, dict, dict]:
    """Conduction-cooled card W x H clamped (wedge locks) along its two long edges y = 0 and y = H to a
    chassis wall at T_wall; short edges and faces adiabatic (sealed box, no airflow). Returns
    (T (ny, nx), system, heat balance). The grid has n nodes along x and a matching spacing along y."""
    ny = int(round((n - 1) * H / W)) + 1
    q = card_heat_load(np.linspace(0.0, W, n), np.linspace(0.0, H, ny), P_chip, P_spread)
    s = conduction_system_2d(W, H, n, ny, kt, q, {"bottom": ("T", T_wall), "top": ("T", T_wall)})
    T = solve_steady(s)
    return T, s, heat_balance(s, T)


# ---------------------------------------------------------------- the lesson's numbers

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save conduction.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    print("materials (Lienhard Tables A.1/A.2, 20 C): k (W/m K), alpha (mm^2/s), effusivity sqrt(k rho cp) (J/m^2 K s^0.5)")
    for name, m in MATERIALS.items():
        print(f"  {name:15s} k = {m['k']:7.2f}   alpha = {diffusivity(name) * 1e6:7.3f}   e = {np.sqrt(m['k'] * m['rho'] * m['cp']):8.0f}")

    # 1. the 1D wall: a cooled chamber wall, closed form against finite volumes
    layers = [(0.001, MATERIALS["Inconel X-750"]["k"]), (0.004, MATERIALS["copper"]["k"])]
    hg, Tg, hc, Tc = 3000.0, 3000.0, 2.0e4, 300.0
    q, Ts = wall_heat_flux(layers, Tg, Tc, hg, hc)
    n = 51
    kk = np.where(np.linspace(0, 0.005, n)[:-1] + 0.0001 / 2 < 0.001, layers[0][1], layers[1][1])
    s1 = conduction_system_1d(0.005, n, kk, 1.0, ("h", hg, Tg), ("h", hc, Tc))
    T1 = solve_steady(s1)
    print(f"\n1D cooled wall (1 mm Inconel X-750 on 4 mm copper; gas h = {hg:.0f}, {Tg:.0f} K; coolant h = {hc:.0f}, {Tc:.0f} K; illustrative):")
    print(f"  closed form q'' = {q / 1e6:.4f} MW/m^2, hot face {Ts[0]:.2f} K, interface {Ts[1]:.2f} K, cold face {Ts[2]:.2f} K")
    print(f"  finite volumes: hot face {T1[0]:.2f} K, interface {T1[10]:.2f} K, cold face {T1[-1]:.2f} K; "
          f"max |diff| {max(abs(T1[0] - Ts[0]), abs(T1[10] - Ts[1]), abs(T1[-1] - Ts[2])):.1e} K")

    # 2. the fin: a 2D sheet losing heat from its faces against the 1D fin closed form
    tf, Lf, Hf, hf = 1.0e-3, 0.040, 0.020, 50.0
    kt_f = MATERIALS["aluminium 6061"]["k"] * tf
    mfin = np.sqrt(hf / kt_f)                       # h_face is the total of both faces: m^2 = h_face / (k t)
    print(f"\nfin: 6061 sheet {tf * 1e3:.0f} mm thick, {Lf * 1e3:.0f} mm long, faces h = {hf:.0f} W/m^2K in total, "
          f"base 40 K above the air (illustrative): m L = {mfin * Lf:.4f}")
    errs_fin = []
    for nn in (11, 21, 41, 81):
        sf = conduction_system_2d(Lf, Hf, nn, 5, kt_f, 0.0, {"left": ("T", 40.0)}, h_face=hf, T_face=0.0)
        Tf = solve_steady(sf)
        errs_fin.append(np.max(np.abs(Tf[2] - fin_closed_form(sf["x"], Lf, mfin, 40.0))))
        hb = heat_balance(sf, Tf)
    Q_fin = fin_heat_rate(kt_f * Hf, Lf, mfin, 40.0)
    print(f"  max |FV - closed form| for n = 11, 21, 41, 81: " + ", ".join(f"{e:.2e}" for e in errs_fin)
          + f" K; orders {np.array2string(observed_order(errs_fin), precision=2)}")
    print(f"  heat rate: closed form {Q_fin:.4f} W, FV {hb['out_faces']:.4f} W (energy imbalance {hb['imbalance']:.1e} W); "
          f"efficiency {Q_fin / (hf * Lf * Hf * 40.0):.3f}")

    # 3. the rectangle with one heated edge against its Fourier series
    print("\nunit square, T = 1 on the top edge, 0 on the others: FV against the Fourier series (y <= 0.75)")
    for nn in (21, 41, 81):
        s = conduction_system_2d(1.0, 1.0, nn, nn, 1.0, 0.0, {"left": ("T", 0.0), "right": ("T", 0.0),
                                                            "bottom": ("T", 0.0), "top": ("T", 1.0)})
        T = solve_steady(s)
        X, Y = np.meshgrid(s["x"], s["y"])
        e = np.max(np.abs(T - rectangle_fourier(X, Y, 1.0, 1.0))[Y <= 0.75])
        print(f"  n = {nn:3d}: max error {e:.3e}; centre T = {T[nn // 2, nn // 2]:.6f} (series {float(rectangle_fourier(0.5, 0.5, 1, 1)):.6f})")
    errs = []
    for nn in (11, 21, 41, 81):
        s = conduction_system_2d(2.0, 1.0, 2 * nn - 1, nn, 1.0, 0.0, {"left": ("T", 0.0), "right": ("T", 0.0),
                                 "bottom": ("T", 0.0), "top": ("T", lambda xx: np.sin(np.pi * xx / 2))})
        T = solve_steady(s)
        X, Y = np.meshgrid(s["x"], s["y"])
        errs.append(np.max(np.abs(T - np.sin(np.pi * X / 2) * np.sinh(np.pi * Y / 2) / np.sinh(np.pi / 2))))
    print(f"  smooth edge sin(pi x / W): errors " + ", ".join(f"{e:.2e}" for e in errs)
          + f"; observed orders {np.array2string(observed_order(errs), precision=3)}")

    # 4. the avionics card
    T_wall = 343.15
    kt0 = board_sheet_conductance()
    print(f"\nconduction-cooled card 160 x 100 mm, clamped on both long edges at {T_wall - 273.15:.0f} C; 10 W processor "
          f"(30 x 30 mm) + 5 W spread (illustrative)")
    print(f"  sheet conductance k t: 1.6 mm FR4 alone {MATERIALS['FR4']['k'] * 1.6e-3:.5f} W/K; with four 35 um copper planes "
          f"{kt0:.5f} W/K (effective in-plane k = {kt0 / 1.6e-3:.1f} W/m K)")
    for nn in (41, 81, 161):
        T, s, hb = solve_card(kt0, nn, T_wall)
        print(f"  n = {nn:3d}: hottest point {T.max() - 273.15:.2f} C; heat out through the clamps {hb['out_fixed']:.6f} W "
              f"of {hb['generated']:.6f} W generated (imbalance {hb['imbalance']:.1e} W)")
    print("  aluminium 6061 heat frame thickness (mm), k t (W/K), hottest point (C):")
    for tfr in (0.0, 0.5, 1.0, 1.5, 2.0, 3.0):
        kt = board_sheet_conductance(t_frame=tfr * 1e-3)
        print(f"    {tfr:.1f}   {kt:.4f}   {solve_card(kt, 81, T_wall)[0].max() - 273.15:.2f}")
    from scipy.optimize import brentq
    t_need = brentq(lambda tt: solve_card(board_sheet_conductance(t_frame=tt), 81, T_wall)[0].max() - (373.15), 0.0, 3e-3, xtol=1e-9)
    print(f"  frame needed for a 100 C board limit under the processor: {t_need * 1e3:.3f} mm")
    Tc_card, sc, _ = solve_card(kt0, 81, T_wall)

    # 5. transients: the semi-infinite checks
    cu = MATERIALS["copper"]
    a_cu = diffusivity("copper")
    print(f"\ncopper slab 50 mm, at 300 K, surface stepped to 1000 K; t = 0.5 s (sqrt(alpha t) = {np.sqrt(a_cu * 0.5) * 1e3:.2f} mm), Fo = 1/2 steps")
    print("  n     implicit: max error (K), energy error      Crank-Nicolson: max error (K), energy error")
    errs_i, errs_c = [], []
    for nn in (26, 51, 101, 201):
        s = conduction_system_1d(0.05, nn, cu["k"], cu["rho"] * cu["cp"], ("T", 1000.0), ("adiabatic",))
        dxx = 0.05 / (nn - 1)
        steps = int(round(0.5 / (0.5 * dxx ** 2 / a_cu)))
        row = []
        for th, store in ((1.0, errs_i), (0.5, errs_c)):
            out = march(s, 300.0, 0.5 / steps, steps, theta=th)
            ex = semi_infinite_step(s["x"], 0.5, a_cu, 300.0, 1000.0)
            E = out["E_stored"] + s["C"][0] * 700.0
            store.append(np.max(np.abs(out["T"] - ex)))
            row.append(f"{store[-1]:.3e}, {E / semi_infinite_step_energy(0.5, cu['k'], a_cu, 300.0, 1000.0) - 1:+.2e}")
        print(f"  {nn:3d}   {row[0]}                 {row[1]}")
    print(f"  orders: implicit {np.array2string(observed_order(errs_i), precision=2)}, CN {np.array2string(observed_order(errs_c), precision=2)}")
    s = conduction_system_1d(0.05, 101, cu["k"], cu["rho"] * cu["cp"], ("h", 3000.0, 3000.0), ("adiabatic",))
    out = march(s, 300.0, 1e-3, 1000, theta=0.5)
    print(f"  convection (h = 3000, 3000 K) at t = 1 s: max |FV - eq. 5.53| = "
          f"{np.max(np.abs(out['T'] - semi_infinite_convection(s['x'], 1.0, a_cu, cu['k'], 3000.0, 300.0, 3000.0))):.2e} K; "
          f"stored {out['E_stored']:.6e} J/m^2 vs heat in {out['Q_in']:.6e} J/m^2")
    for Fo in (0.45, 0.55):
        s = conduction_system_1d(0.05, 51, cu["k"], cu["rho"] * cu["cp"], ("T", 1000.0), ("adiabatic",))
        dt = Fo * (0.05 / 50) ** 2 / a_cu
        out = march(s, 300.0, dt, 400, theta=0.0)
        print(f"  explicit, Fo = {Fo}: max T after 400 steps {np.max(out['T']):.4g} K")

    # 6. the heat-sink chamber wall
    h, Tgas, Tlim = 3000.0, 3000.0, 800.0
    print(f"\nheat-sink chamber wall, gas h = {h:.0f} W/m^2K at {Tgas:.0f} K, initial 300 K, hot-face limit {Tlim:.0f} K (illustrative):")
    # the low-conductivity walls are heated through under 1 mm, so they need a much finer grid than copper:
    # at n = 101 (0.25 mm cells) their times come out about 1 % long, above the semi-infinite bound
    for mat, thick, n, dt in (("copper", 0.025, 101, 0.01), ("Inconel X-750", 0.025, 1601, 1e-4),
                              ("stainless 304", 0.025, 1601, 1e-4)):
        t_fd, res = time_to_surface_limit(mat, thick, h, Tgas, 300.0, Tlim, n=n, dt=dt)
        t_si = semi_infinite_time_to_limit(mat, h, Tgas, 300.0, Tlim)
        m = MATERIALS[mat]
        t_lump = m["rho"] * m["cp"] * thick / h * np.log((Tgas - 300.0) / (Tgas - Tlim))
        print(f"  {mat:14s} {thick * 1e3:.0f} mm: finite wall {t_fd:.3f} s, semi-infinite {t_si:.3f} s, lumped {t_lump:.2f} s; "
              f"sqrt(alpha t) = {np.sqrt(diffusivity(mat) * t_fd) * 1e3:.2f} mm, Bi = {h * thick / m['k']:.3f}; "
              f"energy check {res['E_stored'] - res['Q_in']:.1e} J/m^2")
    print("  copper wall thickness (mm): time to the limit (s)")
    for thick in (0.005, 0.010, 0.020, 0.025, 0.040, 0.060):
        print(f"    {thick * 1e3:4.0f}: {time_to_surface_limit('copper', thick, h, Tgas, 300.0, Tlim)[0]:.3f}")

    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2))
    im = ax[0].contourf(sc["x"] * 1e3, sc["y"] * 1e3, Tc_card - 273.15, 20, cmap="inferno")
    fig.colorbar(im, ax=ax[0], label="T (C)")
    ax[0].set_xlabel("x (mm)")
    ax[0].set_ylabel("y (mm), clamped at y = 0 and 100")
    ax[0].set_title("card with four copper planes, no heat frame")
    ax[0].set_aspect("equal")
    s = conduction_system_1d(0.025, 101, cu["k"], cu["rho"] * cu["cp"], ("h", h, Tgas), ("adiabatic",))
    for tt in (0.5, 1.0, 2.0, 4.0):
        out = march(s, 300.0, 0.01, int(round(tt / 0.01)), theta=0.5)
        line, = ax[1].plot(s["x"] * 1e3, out["T"], label=f"t = {tt:g} s, 25 mm wall")
        ax[1].plot(s["x"] * 1e3, semi_infinite_convection(s["x"], tt, a_cu, cu["k"], h, 300.0, Tgas), "--", color=line.get_color(), lw=0.9)
    ax[1].axhline(Tlim, color="k", lw=0.6)
    ax[1].set_xlabel("depth from the hot face (mm)")
    ax[1].set_ylabel("T (K)")
    ax[1].set_title("copper heat-sink wall (dashed: semi-infinite)")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "conduction.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'conduction.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
