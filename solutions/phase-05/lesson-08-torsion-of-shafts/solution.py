"""P05.08 Torsion of Shafts: reference solution.

A shaft carries power as torque times speed. Its cross-sections rotate rigidly about the axis, the
shear strain grows linearly with radius, and Hooke's law in shear turns that into
tau = T r / J and a twist rate T / (G J). Thin-walled closed sections (a wing box, a thin drive
tube) carry torque as a constant shear flow q = T / (2 A_m) around the cell (Bredt-Batho);
cutting the same wall open loses almost all of its torsional stiffness. This module builds:

* power, speed and torque; the polar second moment of solid and hollow circles; stress and twist;
  stepped shafts and a statically indeterminate shaft between two walls;
* sizing of solid and hollow shafts to an allowable shear stress, and their mass;
* Bredt-Batho for any single-cell polygonal thin-walled section, and the open-section result
  sum(b t^3 / 3) for comparison;
* an independent check: Prandtl's stress function solved by finite differences on a solid
  rectangle (against the exact series) and on a hollow rectangle (against Bredt-Batho).

Units: mm, N, MPa (N/mm^2), torque in N mm, angles in radians unless named _deg, power in W,
speed in rev/min. Material values are typical, not design allowables.

Run:  python solution.py [out_dir]     (prints the lesson's tables, writes torsion.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-05/lesson-08-torsion-of-shafts
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

# Typical properties: E, G (MPa; G = E / 2(1 + nu)), Fty, Ftu (MPa), density (kg/m^3).
MATERIALS: dict[str, dict[str, float]] = {
    "2024-T3":   {"E": 73.1e3, "G": 73.1e3 / 2.66, "Fty": 345.0, "Ftu": 483.0, "rho": 2780.0},
    "7075-T6":   {"E": 71.7e3, "G": 71.7e3 / 2.66, "Fty": 503.0, "Ftu": 572.0, "rho": 2810.0},
    "Ti-6Al-4V": {"E": 113.8e3, "G": 113.8e3 / 2.684, "Fty": 880.0, "Ftu": 950.0, "rho": 4430.0},
    "4130":      {"E": 205.0e3, "G": 205.0e3 / 2.58, "Fty": 435.0, "Ftu": 670.0, "rho": 7850.0},
}


def shear_allowables(material: str) -> tuple[float, float]:
    """(yield, ultimate) shear strengths used in this lesson: Fsy = Fty / sqrt(3) (the von Mises
    estimate, derived in P05.14) and Fsu = 0.6 Ftu (a common approximation for ductile metals)."""
    m = MATERIALS[material]
    return m["Fty"] / np.sqrt(3.0), 0.6 * m["Ftu"]


# ---------------------------------------------------------------- circular shafts

def shaft_torque(power_W: float, rpm: float) -> float:
    """Torque [N mm] transmitted by a shaft carrying power_W [W] at rpm [rev/min]: T = P / omega."""
    return power_W / (rpm * 2.0 * np.pi / 60.0) * 1e3


def polar_moment(d_o: float, d_i: float = 0.0) -> float:
    """Polar second moment of area J = pi (d_o^4 - d_i^4) / 32 [mm^4] of a solid or hollow circle."""
    return np.pi * (d_o ** 4 - d_i ** 4) / 32.0


def torsion_shear_stress(T: float, r: float, J: float) -> float:
    """Shear stress at radius r in a circular shaft: tau = T r / J [MPa]."""
    return T * r / J


def twist_angle(T: float, L: float, G: float, J: float) -> float:
    """Angle of twist [rad] of a uniform shaft of length L: phi = T L / (G J)."""
    return T * L / (G * J)


def stepped_shaft_twist(T: np.ndarray, L: np.ndarray, G: np.ndarray, J: np.ndarray) -> float:
    """Total twist of segments in series, each with its own internal torque: sum T_i L_i / (G_i J_i)."""
    T, L, G, J = (np.asarray(v, dtype=float) for v in (T, L, G, J))
    return float(np.sum(T * L / (G * J)))


def fixed_fixed_shaft(T0: float, L1: float, GJ1: float, L2: float, GJ2: float) -> tuple[float, float, float]:
    """A shaft of two segments between two rigid walls with torque T0 applied at the junction.
    Compatibility (zero total twist) shares T0 in proportion to the segments' torsional
    stiffnesses k_i = GJ_i / L_i: T1 = T0 k1 / (k1 + k2) is carried by segment 1 (to wall A) and
    T2 = -T0 k2 / (k1 + k2) by segment 2. Returns (phi_junction, T1, T2), internal torques with
    the sign convention of the right-hand rule along +x."""
    k1, k2 = GJ1 / L1, GJ2 / L2
    phi = T0 / (k1 + k2)
    return phi, k1 * phi, -k2 * phi


def size_solid_shaft(T: float, tau_allow: float) -> float:
    """Smallest solid diameter whose maximum shear stress is tau_allow: d = (16 T / (pi tau))^(1/3)."""
    return (16.0 * T / (np.pi * tau_allow)) ** (1.0 / 3.0)


def size_hollow_shaft(T: float, tau_allow: float, ratio: float) -> float:
    """Smallest outer diameter of a hollow shaft with d_i = ratio d_o:
    d_o = (16 T / (pi tau (1 - ratio^4)))^(1/3)."""
    return (16.0 * T / (np.pi * tau_allow * (1.0 - ratio ** 4))) ** (1.0 / 3.0)


def size_for_twist(T: float, L: float, G: float, phi_allow: float, ratio: float = 0.0) -> float:
    """Smallest outer diameter for which the twist over L does not exceed phi_allow [rad]."""
    return (32.0 * T * L / (np.pi * G * phi_allow * (1.0 - ratio ** 4))) ** 0.25


def shaft_mass_per_length(d_o: float, d_i: float, rho: float) -> float:
    """Mass per unit length [kg/m] of a tube with diameters in mm and density in kg/m^3."""
    return rho * np.pi * (d_o ** 2 - d_i ** 2) / 4.0 * 1e-6


# ---------------------------------------------------------------- thin-walled sections

def polygon_area(vertices: np.ndarray) -> float:
    """Area enclosed by a closed polygon (shoelace formula), positive for anticlockwise order."""
    v = np.asarray(vertices, dtype=float)
    x, y = v[:, 0], v[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def bredt_batho(T: float, vertices: np.ndarray, t: np.ndarray, G: float) -> dict[str, np.ndarray | float]:
    """Single-cell closed thin-walled section under torque T.

    vertices: (n, 2) corners of the wall MID-LINE, in order; wall k runs from vertex k to k+1 (the
    last closes the cell) with thickness t[k]. Equilibrium of a wall element shows the shear flow
    q = tau t is the same everywhere; its moment about any point is q 2 A_m, so

        q = T / (2 A_m),   tau_k = q / t_k,   dtheta/dx = T / (4 A_m^2 G) sum(s_k / t_k),

    and the torsion constant is J = 4 A_m^2 / sum(s_k / t_k). Returns a dict with A_m, q, tau
    (per wall), twist_rate [rad/mm] and J.
    """
    v = np.asarray(vertices, dtype=float)
    t = np.broadcast_to(np.asarray(t, dtype=float), (len(v),))
    s = np.linalg.norm(np.roll(v, -1, axis=0) - v, axis=1)
    Am = abs(polygon_area(v))
    line_integral = float(np.sum(s / t))
    J = 4.0 * Am ** 2 / line_integral
    q = T / (2.0 * Am)
    return {"A_m": Am, "q": q, "tau": q / t, "twist_rate": T / (G * J), "J": J, "lengths": s}


def thin_tube_J(r_m: float, t: float) -> float:
    """Bredt-Batho torsion constant of a thin circular tube of mean radius r_m: 2 pi r_m^3 t."""
    return 2.0 * np.pi * r_m ** 3 * t


def open_section_J(b: np.ndarray, t: np.ndarray) -> float:
    """Torsion constant of an open thin-walled section made of thin strips: sum b t^3 / 3."""
    return float(np.sum(np.asarray(b, dtype=float) * np.asarray(t, dtype=float) ** 3) / 3.0)


def rectangle_torsion_constant(a: float, b: float, nterms: int = 50) -> float:
    """Exact (series) torsion constant of a solid a x b rectangle, a >= b:
    J = (a b^3 / 3) [1 - (192 / pi^5)(b / a) sum_{n odd} tanh(n pi a / (2 b)) / n^5]."""
    a, b = max(a, b), min(a, b)
    n = np.arange(1, 2 * nterms, 2)
    return a * b ** 3 / 3.0 * (1.0 - 192.0 / np.pi ** 5 * (b / a) * np.sum(np.tanh(n * np.pi * a / (2 * b)) / n ** 5))


def prandtl_fd(a: float, b: float, h: float, hole: tuple[float, float] | None = None) -> tuple[float, np.ndarray]:
    """Torsion constant of an a x b rectangle (optionally with a centred rectangular hole of size
    hole = (ai, bi)) from Prandtl's stress function by finite differences.

    The stress function phi (scaled by G theta) satisfies laplacian(phi) = -2 in the material,
    phi = 0 on the outer boundary and phi = const on the hole boundary; J = 2 integral(phi dA) over
    the full outline (hole included). Discretised as the minimum of the energy
    sum_edges (phi_i - phi_j)^2 / 2 - 2 h^2 sum_i phi_i on a square grid of spacing h, in which
    every node of the hole and its boundary is one unknown (this enforces the constant value and
    the circulation condition automatically). For a solid section this is the usual 5-point scheme.
    The dimensions must be integer multiples of h. Returns (J, phi on the grid).
    """
    nx, ny = int(round(a / h)), int(round(b / h))
    if not (np.isclose(nx * h, a) and np.isclose(ny * h, b)):
        raise ValueError("a and b must be integer multiples of h")
    ix, iy = np.meshgrid(np.arange(nx + 1), np.arange(ny + 1), indexing="ij")
    outer = (ix == 0) | (ix == nx) | (iy == 0) | (iy == ny)
    in_hole = np.zeros_like(outer)
    if hole is not None:
        hx, hy = int(round(hole[0] / h)), int(round(hole[1] / h))
        if (nx - hx) % 2 or (ny - hy) % 2:
            raise ValueError("the hole must be centred on grid lines")
        x0, y0 = (nx - hx) // 2, (ny - hy) // 2
        in_hole = (ix >= x0) & (ix <= x0 + hx) & (iy >= y0) & (iy <= y0 + hy)
    dof = -np.ones(ix.shape, int)
    free = ~outer & ~in_hole
    dof[free] = np.arange(free.sum())
    n_dof = int(free.sum())
    if in_hole.any():
        dof[in_hole] = n_dof
        n_dof += 1
    rows, cols, vals = [], [], []
    for sx, sy in ((1, 0), (0, 1)):                            # every horizontal and vertical edge once
        p = dof[: nx + 1 - sx, : ny + 1 - sy].ravel()
        q = dof[sx:, sy:].ravel()
        for i, j in ((p, q), (q, p)):
            m = i >= 0
            rows.append(i[m]); cols.append(i[m]); vals.append(np.ones(m.sum()))
            mm = m & (j >= 0) & (i != j)
            rows.append(i[mm]); cols.append(j[mm]); vals.append(-np.ones(mm.sum()))
        same = (p >= 0) & (p == q)                             # edges inside the hole: no energy
        rows.append(p[same]); cols.append(p[same]); vals.append(-2.0 * np.ones(same.sum()))
    K = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n_dof, n_dof))
    f = np.full(n_dof, 2.0 * h * h)
    if in_hole.any():
        f[-1] = 2.0 * h * h * in_hole.sum()
    phi_dof = spla.spsolve(K.tocsc(), f)
    phi = np.zeros(ix.shape)
    phi[dof >= 0] = phi_dof[dof[dof >= 0]]
    return float(2.0 * h * h * phi.sum()), phi


# ---------------------------------------------------------------- the lesson's design case

POWER_W, RPM = 1.4e6, 20_900          # a medium turboshaft output (typical class of engine)
LIMIT_FACTOR = 1.25                    # transient over-torque assumed for limit load
FS_ULT = 1.5
SHAFT_L = 600.0                        # mm, length between couplings
TWIST_LIMIT_DEG = 5.0                  # assumed drive-train stiffness requirement: twist at limit torque
WING_BOX = np.array([[0.0, 0.0], [400.0, 0.0], [400.0, 120.0], [0.0, 120.0]])   # mid-line, mm
WING_BOX_T = np.array([1.6, 3.0, 1.6, 3.0])                                      # skin, rear spar, skin, front spar


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save torsion.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)

    T_rated = shaft_torque(POWER_W, RPM)
    T_lim = LIMIT_FACTOR * T_rated
    T_ult = FS_ULT * T_lim
    print(f"turboshaft output {POWER_W / 1e3:.0f} kW at {RPM} rpm: rated torque {T_rated / 1e3:.1f} N m, "
          f"limit {T_lim / 1e3:.1f} N m, ultimate {T_ult / 1e3:.1f} N m")
    print("material    ratio  Fsy   Fsu   d_o (limit/yield)  d_o (ult)  d_o (twist 5 deg)  chosen d_o  d_i   "
          "kg/m   tau_lim  twist (deg)")
    rows = {}
    for mat in ("Ti-6Al-4V", "4130"):
        Fsy, Fsu = shear_allowables(mat)
        m = MATERIALS[mat]
        for ratio in (0.0, 0.6, 0.8, 0.9):
            d_y = size_hollow_shaft(T_lim, Fsy, ratio)
            d_u = size_hollow_shaft(T_ult, Fsu, ratio)
            d_t = size_for_twist(T_lim, SHAFT_L, m["G"], np.radians(TWIST_LIMIT_DEG), ratio)
            d = max(d_y, d_u, d_t)
            J = polar_moment(d, ratio * d)
            rows[(mat, ratio)] = d
            print(f"{mat:10s}  {ratio:.1f}  {Fsy:4.0f}  {Fsu:4.0f}   {d_y:8.2f}          {d_u:7.2f}     {d_t:8.2f}"
                  f"          {d:7.2f}  {ratio * d:5.2f}  {shaft_mass_per_length(d, ratio * d, m['rho']):5.3f}"
                  f"   {torsion_shear_stress(T_lim, d / 2, J):5.0f}   {np.degrees(twist_angle(T_lim, SHAFT_L, m['G'], J)):5.3f}")
    d0, d8 = rows[("Ti-6Al-4V", 0.0)], rows[("Ti-6Al-4V", 0.8)]
    rho = MATERIALS["Ti-6Al-4V"]["rho"]
    print(f"Ti-6Al-4V hollow (0.8) / solid mass: {shaft_mass_per_length(d8, 0.8 * d8, rho) / shaft_mass_per_length(d0, 0, rho):.3f}")

    # ---- indeterminate shaft
    G = MATERIALS["Ti-6Al-4V"]["G"]
    phi, T1, T2 = fixed_fixed_shaft(T_lim, 200.0, G * polar_moment(24.0), 400.0, G * polar_moment(30.0, 24.0))
    print(f"\nshaft fixed at both ends, limit torque at the junction: T1 = {T1 / 1e3:.1f} N m, T2 = {T2 / 1e3:.1f} N m, "
          f"junction rotation {np.degrees(phi):.3f} deg")

    # ---- Bredt-Batho: thin tube against exact, wing box, open section
    print("\nthin circular tube: Bredt J = 2 pi r_m^3 t against exact pi (r_o^4 - r_i^4) / 2")
    for tr in (0.02, 0.05, 0.1, 0.2, 0.5):
        rm, t = 1.0, tr
        exact = polar_moment(2 * (rm + t / 2), 2 * (rm - t / 2))
        print(f"  t / r_m = {tr:4.2f}: ratio Bredt / exact = {thin_tube_J(rm, t) / exact:.5f}")
    al = MATERIALS["2024-T3"]
    T_box = 10e6
    res = bredt_batho(T_box, WING_BOX, WING_BOX_T, al["G"])
    print(f"wing box 400 x 120 mm, skins 1.6 mm, spars 3.0 mm, 2024-T3, T = {T_box / 1e6:.0f} kN m:")
    print(f"  A_m = {res['A_m']:.0f} mm^2, q = {res['q']:.1f} N/mm, tau skin {res['tau'][0]:.1f} MPa, "
          f"tau spar {res['tau'][1]:.1f} MPa, J = {res['J']:.4e} mm^4, twist {np.degrees(res['twist_rate']) * 1e3:.3f} deg/m")
    J_open = open_section_J(res["lengths"], WING_BOX_T)
    tau_open = T_box * WING_BOX_T.max() / J_open
    print(f"  cut open (same walls): J = {J_open:.0f} mm^4, stiffness ratio closed/open = {res['J'] / J_open:.0f}, "
          f"open-section stress at the same torque {tau_open:.0f} MPa")

    # ---- independent check: Prandtl stress function by finite differences
    print("\nPrandtl FD, solid 2 x 1 rectangle: series J = %.6f" % rectangle_torsion_constant(2.0, 1.0))
    errs = []
    for h in (0.1, 0.05, 0.025, 0.0125):
        J, _ = prandtl_fd(2.0, 1.0, h)
        errs.append(abs(J - rectangle_torsion_constant(2.0, 1.0)))
        print(f"  h = {h:.4f}: J = {J:.6f}, error {errs[-1]:.2e}")
    print(f"  observed orders {np.array2string(np.log2(np.array(errs[:-1]) / errs[1:]), precision=2)}")
    print("square tube, outer side 1, wall t: FD (h = t/8) against Bredt (a - t)^3 t")
    for t in (0.05, 0.1, 0.2):
        J, _ = prandtl_fd(1.0, 1.0, t / 8, hole=(1 - 2 * t, 1 - 2 * t))
        jb = bredt_batho(1.0, np.array([[0, 0], [1 - t, 0], [1 - t, 1 - t], [0, 1 - t]]), t, 1.0)["J"]
        print(f"  t = {t:.2f}: FD J = {J:.5f}, Bredt J = {jb:.5f}, Bredt / FD = {jb / J:.4f}")

    # ---- figure
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    for mat in ("Ti-6Al-4V", "4130"):
        ratios = np.linspace(0, 0.95, 60)
        Fsy, Fsu = shear_allowables(mat)
        m = MATERIALS[mat]
        mass = []
        for r in ratios:
            d = max(size_hollow_shaft(T_lim, Fsy, r), size_hollow_shaft(T_ult, Fsu, r),
                    size_for_twist(T_lim, SHAFT_L, m["G"], np.radians(TWIST_LIMIT_DEG), r))
            mass.append(shaft_mass_per_length(d, r * d, m["rho"]))
        ax[0].plot(ratios, mass, label=mat)
    ax[0].set_xlabel("d_i / d_o"); ax[0].set_ylabel("mass per length (kg/m)"); ax[0].legend()
    ax[0].set_title("shaft mass for the same torque")
    r = np.linspace(0, d8 / 2, 50)
    ax[1].plot(r, torsion_shear_stress(T_lim, r, polar_moment(d0)), label=f"solid d = {d0:.1f} mm")
    rr = np.linspace(0.4 * d8, d8 / 2, 20)
    ax[1].plot(rr, torsion_shear_stress(T_lim, rr, polar_moment(d8, 0.8 * d8)), label=f"hollow d_o = {d8:.1f} mm")
    ax[1].set_xlabel("radius (mm)"); ax[1].set_ylabel("shear stress at limit (MPa)"); ax[1].legend()
    ax[1].set_title("tau = T r / J")
    _, phi = prandtl_fd(1.0, 1.0, 0.0125, hole=(0.6, 0.6))
    ax[2].contourf(phi.T, levels=20, extent=(0, 1, 0, 1))
    ax[2].set_aspect("equal"); ax[2].set_title("Prandtl stress function, square tube")
    fig.tight_layout()
    fig.savefig(out_dir / "torsion.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'torsion.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
