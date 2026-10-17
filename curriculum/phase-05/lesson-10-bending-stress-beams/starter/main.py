"""P05.10 Bending Stress in Beams: starter code.

The main spar of a wing is a cantilever built into the fuselage and loaded upward by lift. This
module computes the bending stress in it, from the section up:

* section properties of a section built from rectangles (area, centroid, I_xx, I_yy, I_xy about
  the centroid), with ready-made I-section and Z-section builders;
* the flexure formula sigma = -M y / I for symmetric bending, and the general (unsymmetric) formula
  for a section with a product of inertia or a moment about two axes, with the neutral-axis angle;
* an independent "fibre" model of the section: chop it into small cells, impose a plane strain field
  eps = a + b x + c y, and solve the three equilibrium equations N = 0, M, M_h for a, b, c by plain
  sums over the cells, with no section-property formula anywhere;
* the spanwise load of an elliptically loaded wing, its shear force and bending moment in closed
  form and by numerical integration;
* an independent 1-D finite-element cantilever (Hermite cubic Euler-Bernoulli elements, Phase 03):
  its curvature times E c gives the bending stress without using the flexure formula;
* spar-cap sizing with the two-flange idealisation and aerospace margins of safety.

Sign convention (stated once, used everywhere): x is horizontal in the section (chordwise), y is
vertical (up), the beam axis runs along the span. A bending moment M is positive when it compresses
the fibres at positive y (the top): sigma = -M y / I. Lift on a cantilever wing is a positive M.
M_h is the moment about the vertical axis, positive when it compresses the fibres at positive x.
The stress resultants are therefore M = -integral(sigma y dA) and M_h = -integral(sigma x dA).

Units: SI (N, m, Pa). Material properties are typical room-temperature values, not design allowables.

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes wing_spar_stress.png)
Test: pytest curriculum/phase-05/lesson-10-bending-stress-beams
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

G0 = 9.80665                          # m/s^2, standard gravity


@dataclass(frozen=True)
class Material:
    """Typical room-temperature properties (publicly viewable ASM Aerospace Specification Metals data sheets).
    E and G in Pa, F_ty (0.2 % offset yield) and F_tu (ultimate) in Pa, density in kg/m^3."""
    name: str
    E: float
    G: float
    F_ty: float
    F_tu: float
    rho: float


AL7075_T6 = Material("7075-T6", 71.7e9, 26.9e9, 503e6, 572e6, 2810.0)
AL2024_T3 = Material("2024-T3", 73.1e9, 28.0e9, 345e6, 483e6, 2780.0)
TI6AL4V = Material("Ti-6Al-4V", 113.8e9, 44.0e9, 880e6, 950e6, 4430.0)
STEEL4130 = Material("4130 normalized", 205e9, 80e9, 435e6, 670e6, 7850.0)


# ---------------------------------------------------------------- section properties

def section_properties(rects: np.ndarray) -> dict:
    """Properties of a section made of rectangles.

    ``rects`` has shape (n, 4); each row is (b, h, xc, yc): width along x, height along y and the
    centre of one rectangle. Returns a dict with A (m^2), the centroid xc, yc (m) and the second
    moments about centroidal axes parallel to x and y:
    Ixx = integral(y^2 dA), Iyy = integral(x^2 dA), Ixy = integral(x y dA)  (m^4),
    each rectangle's own b h^3 / 12 (or h b^3 / 12, or 0) moved by the parallel-axis theorem.
    """
    raise NotImplementedError


def i_section(h: float, b_f: float, t_f: float, t_w: float) -> np.ndarray:
    """Rectangles of a doubly symmetric I-section of overall depth h, flange width b_f, flange
    thickness t_f and web thickness t_w, centred on the origin (rows as in ``section_properties``)."""
    hw = h - 2 * t_f
    return np.array([[b_f, t_f, 0.0, (h - t_f) / 2],
                     [t_w, hw, 0.0, 0.0],
                     [b_f, t_f, 0.0, -(h - t_f) / 2]])


def z_section(h: float, b_f: float, t: float) -> np.ndarray:
    """Rectangles of a Z-section of overall depth h, flange width b_f (overall, including the web
    thickness) and uniform thickness t. The web is centred on x = 0; the top flange runs from
    x = -t/2 to b_f - t/2 (pointing to +x), the bottom flange from -(b_f - t/2) to t/2. The section is
    antisymmetric, so its centroid is at the origin, and it has a nonzero product of inertia."""
    hw = h - 2 * t
    xf = (b_f - t) / 2
    return np.array([[b_f, t, xf, (h - t) / 2],
                     [t, hw, 0.0, 0.0],
                     [b_f, t, -xf, -(h - t) / 2]])


# ---------------------------------------------------------------- bending stress

def flexure_stress(M: float, y: np.ndarray, I: float) -> np.ndarray:
    """Flexure formula for symmetric bending: sigma = -M y / I, with y measured from the neutral
    (centroidal) axis and M positive when it compresses the top (y > 0)."""
    raise NotImplementedError


def bending_stress_general(M: float, M_h: float, x: np.ndarray, y: np.ndarray,
                           Ixx: float, Iyy: float, Ixy: float) -> np.ndarray:
    """Bending stress for unsymmetric bending, x and y measured from the centroid:

        sigma = -[ (M Iyy - M_h Ixy) y + (M_h Ixx - M Ixy) x ] / (Ixx Iyy - Ixy^2).

    It reduces to -M y / Ixx when Ixy = 0 and M_h = 0."""
    raise NotImplementedError


def neutral_axis_angle(M: float, M_h: float, Ixx: float, Iyy: float, Ixy: float) -> float:
    """Angle (rad, from the +x axis, in (-pi/2, pi/2]) of the neutral axis, the line through the
    centroid on which ``bending_stress_general`` is zero."""
    raise NotImplementedError


def fiber_mesh(rects: np.ndarray, n: int = 40) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Split every rectangle into n x n equal cells; return the cell centres x, y and areas dA."""
    xs, ys, das = [], [], []
    for b, h, xc, yc in np.atleast_2d(np.asarray(rects, dtype=float)):
        u = (np.arange(n) + 0.5) / n - 0.5
        X, Y = np.meshgrid(xc + u * b, yc + u * h)
        xs.append(X.ravel())
        ys.append(Y.ravel())
        das.append(np.full(n * n, b * h / n ** 2))
    return np.concatenate(xs), np.concatenate(ys), np.concatenate(das)


def fiber_strain_plane(rects: np.ndarray, M: float, M_h: float, E: float, N: float = 0.0,
                       n: int = 40) -> np.ndarray:
    """Independent check of the flexure formulas. Assume plane sections stay plane, so the axial
    strain is eps = a + b x + c y in the section's own (absolute) coordinates, and sigma = E eps.
    Impose the three equilibrium equations with sums over the fibre cells

        sum(sigma dA) = N,   -sum(sigma y dA) = M,   -sum(sigma x dA) = M_h

    and solve the 3 x 3 linear system for (a, b, c). No centroid or I formula is used."""
    raise NotImplementedError


def fiber_stress(rects: np.ndarray, M: float, M_h: float, x: np.ndarray, y: np.ndarray,
                 E: float = 70e9, n: int = 40) -> np.ndarray:
    """Stress at points (x, y) (absolute section coordinates) from ``fiber_strain_plane``."""
    raise NotImplementedError


# ---------------------------------------------------------------- the wing load

def elliptic_lift(y: np.ndarray, L_semi: float, s: float) -> np.ndarray:
    """Lift per unit span (N/m) of an elliptically loaded semi-span s carrying total lift L_semi:
    w(y) = w0 sqrt(1 - (y/s)^2), w0 = 4 L_semi / (pi s), y from the root (0) to the tip (s)."""
    raise NotImplementedError


def elliptic_shear_moment(y: np.ndarray, L_semi: float, s: float) -> tuple[np.ndarray, np.ndarray]:
    """Closed-form shear force V(y) = integral_y^s w and bending moment M(y) = integral_y^s w (eta - y)
    of the elliptic load (the load outboard of station y), with u = y/s and w0 = 4 L_semi/(pi s):

        V = w0 s / 2 [acos u - u sqrt(1 - u^2)],
        M = w0 s^2 [ (1 - u^2)^{3/2} / 3 - (u/2)(acos u - u sqrt(1 - u^2)) ].

    At the root V = L_semi and M = L_semi * 4 s / (3 pi)."""
    raise NotImplementedError


def shear_moment_numeric(y: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Shear force and bending moment of a cantilever loaded by w(y) (y ascending, free tip at y[-1])
    by integrating dV/dy = -w and dM/dy = -V inward from V = M = 0 at the tip (trapezoid rule)."""
    raise NotImplementedError


# ---------------------------------------------------------------- independent check: 1-D FEM

def cantilever_fem(L: float, EI: float, n_el: int, load: Callable[[np.ndarray], np.ndarray]) -> dict:
    """Euler-Bernoulli cantilever (clamped at x = 0, free at x = L) by n_el Hermite cubic elements.

    ``load(x)`` is the upward distributed load (N/m); consistent nodal loads use 5-point Gauss
    quadrature on each element. Returns a dict with the node positions ``x``, deflection ``w`` (up),
    slope ``theta``, the bending moment at the nodes recovered from the element end forces
    ``M_nodes`` (sagging positive), and the curvature at element midpoints ``kappa_mid`` with their
    positions ``x_mid``. The curvature is the FE route to the stress: sigma = -E kappa y."""
    raise NotImplementedError


# ---------------------------------------------------------------- spar caps and margins

def two_flange_stress(M: float, A_cap: float, h_c: float) -> float:
    """Cap stress of the idealised two-flange spar: the caps (area A_cap each, centroids h_c apart)
    carry all the bending as a couple, sigma = M / (A_cap h_c). The web is assumed to carry shear only."""
    raise NotImplementedError


def cap_area_required(M: float, h_c: float, allowable: float, factor: float = 1.5) -> float:
    """Smallest cap area for which factor * M / (A h_c) equals the allowable stress."""
    raise NotImplementedError


def margin_of_safety(allowable: float, applied: float, factor: float = 1.0) -> float:
    """Aerospace margin of safety MS = allowable / (factor * applied) - 1 (MS >= 0 passes)."""
    raise NotImplementedError


@dataclass
class Spar:
    """Root section of the lesson's wing spar: a 7075-T6 I-section."""
    h: float = 0.180          # m, overall depth
    b_f: float = 0.060        # m, cap (flange) width
    t_f: float = 0.010        # m, cap thickness
    t_w: float = 0.003        # m, web thickness
    material: Material = AL7075_T6

    @property
    def rects(self) -> np.ndarray:
        return i_section(self.h, self.b_f, self.t_f, self.t_w)

    @property
    def I(self) -> float:
        return section_properties(self.rects)["Ixx"]

    @property
    def c(self) -> float:
        return self.h / 2

    @property
    def A_cap(self) -> float:
        return self.b_f * self.t_f

    @property
    def h_c(self) -> float:
        return self.h - self.t_f


def wing_semispan_lift(mass: float = 750.0, n_limit: float = 3.8) -> float:
    """Limit lift on one wing of an aircraft of ``mass`` kg at load factor n: n m g / 2 (tail load
    and wing weight relief neglected)."""
    return n_limit * mass * G0 / 2.0


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save wing_spar_stress.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    mass, n_lim, s = 750.0, 3.8, 5.0
    L = wing_semispan_lift(mass, n_lim)
    V0, M0 = elliptic_shear_moment(0.0, L, s)
    print(f"aircraft {mass:.0f} kg at n = {n_lim}: lift per wing {L:.1f} N, semi-span {s} m")
    print(f"elliptic load: w0 = {4 * L / (np.pi * s):.1f} N/m, root shear {float(V0):.1f} N, "
          f"root moment {float(M0):.1f} N m (centre of lift at {4 * s / (3 * np.pi):.4f} m)")
    y = np.linspace(0.0, s, 2001)
    Vn, Mn = shear_moment_numeric(y, elliptic_lift(y, L, s))
    print(f"numerical integration (2001 points): root moment {Mn[0]:.1f} N m "
          f"(rel. error {abs(Mn[0] - M0) / M0:.1e})")

    sp = Spar()
    p = section_properties(sp.rects)
    mat = sp.material
    sig = M0 * sp.c / sp.I
    print(f"\nroot spar {mat.name} I-section h = {sp.h * 1e3:.0f} mm, caps {sp.b_f * 1e3:.0f} x "
          f"{sp.t_f * 1e3:.0f} mm, web {sp.t_w * 1e3:.0f} mm: A = {p['A'] * 1e6:.0f} mm^2, "
          f"I = {sp.I * 1e12:.4g} mm^4, S = I/c = {sp.I / sp.c * 1e9:.1f} mm^3")
    print(f"flexure formula at the root: top {flexure_stress(M0, sp.c, sp.I) / 1e6:.1f} MPa, "
          f"bottom {flexure_stress(M0, -sp.c, sp.I) / 1e6:.1f} MPa (limit load)")
    sf = fiber_stress(sp.rects, float(M0), 0.0, np.array([0.0, 0.0]), np.array([sp.c, -sp.c]), mat.E, 40)
    print(f"fibre model (40 x 40 cells per rectangle): top {sf[0] / 1e6:.2f}, bottom {sf[1] / 1e6:.2f} MPa")
    s2 = two_flange_stress(M0, sp.A_cap, sp.h_c)
    print(f"two-flange idealisation: cap stress {s2 / 1e6:.1f} MPa (flexure at cap centroid "
          f"{M0 * (sp.h_c / 2) / sp.I / 1e6:.1f} MPa)")
    print(f"margins: yield (limit, factor 1.0) MS = {margin_of_safety(mat.F_ty, sig, 1.0):.3f}; "
          f"ultimate (factor 1.5) MS = {margin_of_safety(mat.F_tu, sig, 1.5):.3f}")
    A_req = cap_area_required(M0, sp.h_c, mat.F_tu, 1.5)
    A_req_y = cap_area_required(M0, sp.h_c, mat.F_ty, 1.0)
    print(f"two-flange cap area for MS = 0: ultimate {A_req * 1e6:.0f} mm^2, yield {A_req_y * 1e6:.0f} mm^2 "
          f"(fitted {sp.A_cap * 1e6:.0f} mm^2)")

    # independent check: FEM curvature -> stress
    EI = mat.E * sp.I
    print("\nFEM cantilever under the elliptic load, root stress from curvature at the first element midpoint")
    print("  n_el   sigma_FEM(x_mid) MPa   flexure M(x_mid) c/I MPa   rel. diff     M_nodes(root) rel. diff")
    for n_el in (5, 10, 20, 40, 80):
        r = cantilever_fem(s, EI, n_el, lambda x: elliptic_lift(x, L, s))
        sf = -mat.E * r["kappa_mid"][0] * (-sp.c)        # bottom fibre, tension
        sm = float(elliptic_shear_moment(r["x_mid"][0], L, s)[1]) * sp.c / sp.I
        print(f"  {n_el:4d}   {sf / 1e6:10.3f}             {sm / 1e6:10.3f}             "
              f"{abs(sf - sm) / sm:.2e}     {abs(r['M_nodes'][0] - M0) / M0:.1e}")

    # unsymmetric preview: Z stringer
    zr = z_section(0.040, 0.020, 0.002)
    zp = section_properties(zr)
    Mz = 100.0
    xt, yt = 0.019, 0.019                                # top-flange tip, mid-thickness
    s_gen = bending_stress_general(Mz, 0.0, xt, yt, zp["Ixx"], zp["Iyy"], zp["Ixy"])
    s_naive = flexure_stress(Mz, yt, zp["Ixx"])
    ang = np.degrees(neutral_axis_angle(Mz, 0.0, zp["Ixx"], zp["Iyy"], zp["Ixy"]))
    print(f"\nZ stringer 40 x 20 x 2 mm: Ixx = {zp['Ixx'] * 1e12:.0f}, Iyy = {zp['Iyy'] * 1e12:.0f}, "
          f"Ixy = {zp['Ixy'] * 1e12:.0f} mm^4")
    print(f"M = {Mz:.0f} N m about x: stress at top-flange tip (19, 19) mm: general {s_gen / 1e6:.1f} MPa, "
          f"naive -My/Ixx {s_naive / 1e6:.1f} MPa; neutral axis at {ang:.1f} deg to x")
    s_fib = fiber_stress(zr, Mz, 0.0, np.array([xt]), np.array([yt]), 70e9, 40)[0]
    print(f"fibre model at the same point: {s_fib / 1e6:.2f} MPa")

    # figure
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    Vc, Mc = elliptic_shear_moment(y, L, s)
    ax[0].plot(y, elliptic_lift(y, L, s) / 1e3, label="w (kN/m)")
    ax[0].plot(y, Vc / 1e3, label="V (kN)")
    ax[0].plot(y, Mc / 1e3, label="M (kN m)")
    ax[0].set_xlabel("span station y (m)")
    ax[0].set_title("elliptic lift: load, shear, moment")
    ax[0].legend()
    r = cantilever_fem(s, EI, 20, lambda x: elliptic_lift(x, L, s))
    ax[1].plot(y, Mc * sp.c / sp.I / 1e6, label="flexure formula M c / I")
    ax[1].plot(r["x_mid"], mat.E * r["kappa_mid"] * sp.c / 1e6, "o", ms=4, label="FEM E kappa c (20 el.)")
    ax[1].axhline(mat.F_tu / 1.5 / 1e6, ls="--", c="k", lw=0.8, label="F_tu / 1.5")
    ax[1].set_xlabel("span station y (m)")
    ax[1].set_ylabel("bottom-cap stress (MPa)")
    ax[1].set_title("spar stress at limit load")
    ax[1].legend()
    yy = np.linspace(-sp.c, sp.c, 50)
    ax[2].plot(flexure_stress(M0, yy, sp.I) / 1e6, yy * 1e3)
    ax[2].axvline(0, c="k", lw=0.5)
    ax[2].set_xlabel("sigma (MPa)")
    ax[2].set_ylabel("y (mm)")
    ax[2].set_title("root section: linear stress, zero at centroid")
    fig.tight_layout()
    fig.savefig(out_dir / "wing_spar_stress.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'wing_spar_stress.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
