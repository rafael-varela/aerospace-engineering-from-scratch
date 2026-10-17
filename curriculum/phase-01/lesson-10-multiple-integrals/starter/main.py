"""P01.10 Multiple Integrals: starter code.

Mass properties of a tapered, swept wing box modelled as a solid of uniform effective density
rho occupying, for spanwise station 0 <= y <= L,

    x from x_le(y) + f_0 c(y) to x_le(y) + f_1 c(y),     (box between two chord fractions)
    z from -1/2 tau c(y) to +1/2 tau c(y),               (box depth tau * chord)
    c(y) = c_r + (c_t - c_r) y / L,  x_le(y) = y tan(sweep).

Mass, centre of mass and the inertia tensor are triple integrals. They are evaluated three ways:
closed forms (polynomials in y, integrated by SymPy), a composite midpoint rule on the unit cube
after the change of variables (xi, eta, zeta) -> (x, y, z) with Jacobian c(y)**2 (f_1 - f_0) tau L,
and scipy.integrate.tplquad as the library oracle. A polar double integral gives the inertia of a
cylindrical fuselage frame, and the parallel-axis theorem moves tensors between reference points.

Complete every function marked NotImplementedError.

Run:  python main.py          (prints the mass-property table, writes wing_box_mass.png)
Test: pytest curriculum/phase-01/lesson-10-multiple-integrals
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import sympy as sp


@dataclass(frozen=True)
class WingBox:
    """Wing-box geometry and effective density in SI units.

    L : half-span (root to tip), m
    c_root, c_tip : root and tip chords, m
    f0, f1 : chord fractions of the box's front and rear walls (0 <= f0 < f1 <= 1)
    tau : box depth as a fraction of the local chord
    sweep : leading-edge sweep angle, rad (x_le = y tan(sweep))
    rho : effective (smeared) density of the box, kg/m**3
    """

    L: float = 5.0
    c_root: float = 1.6
    c_tip: float = 0.8
    f0: float = 0.15
    f1: float = 0.65
    tau: float = 0.10
    sweep: float = np.radians(10.0)
    rho: float = 180.0


DEFAULT_BOX = WingBox()


def _check(box: WingBox) -> None:
    if box.L <= 0 or box.c_root <= 0 or box.c_tip <= 0 or box.tau <= 0 or box.rho <= 0:
        raise ValueError("lengths, depth fraction and density must be positive")
    if not (0.0 <= box.f0 < box.f1 <= 1.0):
        raise ValueError("need 0 <= f0 < f1 <= 1")


def chord(y: np.ndarray | float, box: WingBox = DEFAULT_BOX) -> np.ndarray | float:
    """Return the local chord c(y) = c_root + (c_tip - c_root) y / L for 0 <= y <= L. Raise ``ValueError`` outside."""
    raise NotImplementedError


def map_to_box(xi: np.ndarray, eta: np.ndarray, zeta: np.ndarray, box: WingBox = DEFAULT_BOX) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Map unit-cube coordinates (xi, eta, zeta) in [0, 1]**3 to (x, y, z) in the wing box; also return the Jacobian.

    y = eta L; x = y tan(sweep) + (f0 + (f1 - f0) xi) c(y); z = (zeta - 1/2) tau c(y).
    The Jacobian determinant is (f1 - f0) c(y) * L * tau c(y) = (f1 - f0) tau L c(y)**2.
    """
    raise NotImplementedError


def midpoint_3d(f, n: int, box: WingBox = DEFAULT_BOX) -> float:
    """Integrate f(x, y, z) over the wing box with an n x n x n composite midpoint rule on the unit cube.

    ``f`` receives broadcastable arrays and must return an array of the same shape. The result is
    sum f(x_ijk, y_ijk, z_ijk) * J_ijk / n**3. Raise ``ValueError`` for n < 1.
    """
    raise NotImplementedError


def mass_properties_midpoint(n: int, box: WingBox = DEFAULT_BOX) -> tuple[float, np.ndarray, np.ndarray]:
    """Return (mass, centre of mass, inertia tensor about the origin) by the midpoint rule.

    Inertia tensor I_ij = rho * integral (|r|**2 delta_ij - r_i r_j) dV, in kg m**2, about the
    origin (root leading edge), axes x aft, y spanwise, z up.
    """
    raise NotImplementedError


def mass_properties_analytic(box: WingBox = DEFAULT_BOX) -> tuple[float, np.ndarray, np.ndarray]:
    """Return (mass, centre of mass, inertia tensor about the origin) in closed form via SymPy.

    The integrand is a polynomial in (xi, eta, zeta) after the change of variables, so every
    integral is exact. Numbers are returned as floats.
    """
    raise NotImplementedError


def mass_closed_form(box: WingBox = DEFAULT_BOX) -> float:
    """Return the mass from the one-dimensional closed form rho (f1 - f0) tau L (c_r**2 + c_r c_t + c_t**2) / 3.

    The cross-section area is (f1 - f0) tau c(y)**2 and integral_0^L c(y)**2 dy for a linear chord
    is L (c_r**2 + c_r c_t + c_t**2) / 3.
    """
    raise NotImplementedError


def spanwise_cg_closed_form(box: WingBox = DEFAULT_BOX) -> float:
    """Return y_cg = L (c_r**2 + 2 c_r c_t + 3 c_t**2) / (4 (c_r**2 + c_r c_t + c_t**2)) for a linear chord."""
    raise NotImplementedError


def parallel_axis(I_cg: np.ndarray, m: float, d: np.ndarray) -> np.ndarray:
    """Return the inertia tensor about a point displaced by ``d`` from the centre of mass: I_cg + m (|d|**2 E - d d^T)."""
    raise NotImplementedError


def inertia_about_cg(I_origin: np.ndarray, m: float, cg: np.ndarray) -> np.ndarray:
    """Return the inertia tensor about the centre of mass from the tensor about the origin (inverse parallel-axis shift)."""
    raise NotImplementedError


def ring_inertia_polar(r_in: float, r_out: float, length: float, rho: float, n: int) -> tuple[float, float]:
    """Return (mass, I_axis) of a hollow cylinder (fuselage frame) by an n x n midpoint rule in polar coordinates.

    m = rho * length * integral integral r dr dtheta,  I_axis = rho * length * integral integral r**2 * r dr dtheta,
    over r in [r_in, r_out], theta in [0, 2 pi]. The factor r is the polar Jacobian. Raise ``ValueError``
    unless 0 <= r_in < r_out, length > 0, rho > 0, n >= 1.
    """
    raise NotImplementedError


def ring_inertia_closed_form(r_in: float, r_out: float, length: float, rho: float) -> tuple[float, float]:
    """Return (mass, I_axis) = (rho L pi (r_out**2 - r_in**2), 1/2 m (r_in**2 + r_out**2)) for a hollow cylinder."""
    raise NotImplementedError


def convergence_study(ns: list[int], box: WingBox = DEFAULT_BOX) -> np.ndarray:
    """Return the relative error of the midpoint mass against the closed form for each n in ``ns``.

    The mapped integrand is the polynomial c(eta)**2, so the midpoint error is a clean second order:
    doubling n divides the error by four.
    """
    raise NotImplementedError


def tplquad_mass_and_ixx(box: WingBox = DEFAULT_BOX) -> tuple[float, float]:
    """Return (mass, I_xx about the origin) from scipy.integrate.tplquad in physical coordinates (the library oracle).

    Integration order: z innermost, then x, then y, with limits that depend on y (a general region).
    """
    raise NotImplementedError


def main(out_dir: str | Path = ".") -> None:
    """Print the mass-property table and save wing_box_mass.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    box = DEFAULT_BOX
    m_a, cg_a, I_a = mass_properties_analytic(box)
    print(f"Wing box: L = {box.L} m, c_root = {box.c_root} m, c_tip = {box.c_tip} m, walls at {box.f0:.2f}c and {box.f1:.2f}c, depth {box.tau:.2f}c, sweep {np.degrees(box.sweep):.0f} deg, rho_eff = {box.rho} kg/m^3")
    print(f"  mass: closed form {mass_closed_form(box):.4f} kg, SymPy {m_a:.4f} kg")
    print(f"  centre of mass (m): {np.round(cg_a, 5)}; closed-form y_cg = {spanwise_cg_closed_form(box):.5f} m ({spanwise_cg_closed_form(box)/box.L:.4f} L)")
    print("  inertia tensor about the root leading edge (kg m^2):")
    print(np.array2string(I_a, precision=4, suppress_small=True))
    I_cg = inertia_about_cg(I_a, m_a, cg_a)
    print("  inertia tensor about the centre of mass (kg m^2):")
    print(np.array2string(I_cg, precision=4, suppress_small=True))
    print(f"  principal moments about the CG: {np.round(np.linalg.eigvalsh(I_cg), 4)}")
    for n in (2, 4, 8, 16, 32):
        m_n, cg_n, I_n = mass_properties_midpoint(n, box)
        print(f"  midpoint n = {n:2d}: mass {m_n:.6f} kg (err {m_n - m_a:+.2e}), y_cg {cg_n[1]:.6f}, I_xx {I_n[0,0]:.5f} (rel err {abs(I_n[0,0]-I_a[0,0])/I_a[0,0]:.2e})")
    errs = convergence_study([2, 4, 8, 16, 32], box)
    print("  observed orders for the mass:", np.round(np.log2(errs[:-1] / errs[1:]), 2))
    m_q, ixx_q = tplquad_mass_and_ixx(box)
    print(f"  scipy tplquad: mass {m_q:.6f} kg, I_xx {ixx_q:.6f} kg m^2 (SymPy {I_a[0,0]:.6f})")
    m_r, I_r = ring_inertia_polar(0.95, 1.0, 0.05, 2700.0, 64)
    m_rc, I_rc = ring_inertia_closed_form(0.95, 1.0, 0.05, 2700.0)
    print(f"  fuselage frame (r 0.95-1.00 m, 50 mm wide, aluminium): polar midpoint m = {m_r:.4f} kg, I = {I_r:.5f} kg m^2; closed form {m_rc:.4f} kg, {I_rc:.5f} kg m^2")

    fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))
    y = np.linspace(0.0, box.L, 200)
    c = chord(y, box)
    xle = y * np.tan(box.sweep)
    axs[0].fill_between(y, xle + box.f0 * c, xle + box.f1 * c, color="0.8", label="wing box (plan view)")
    axs[0].plot(y, xle, "k-", lw=0.8, label="leading edge")
    axs[0].plot(y, xle + c, "k--", lw=0.8, label="trailing edge")
    axs[0].plot(cg_a[1], cg_a[0], "r*", ms=12, label=f"CG ({cg_a[0]:.3f}, {cg_a[1]:.3f}) m")
    axs[0].set_xlabel("spanwise y [m]")
    axs[0].set_ylabel("chordwise x (aft) [m]")
    axs[0].invert_yaxis()
    axs[0].set_aspect("equal")
    axs[0].legend(fontsize=7)
    axs[0].set_title(f"Planform; mass {m_a:.1f} kg")
    ns = [2, 4, 8, 16, 32]
    axs[1].loglog(ns, errs, "o-", label="midpoint mass relative error")
    axs[1].loglog(ns, errs[0] * (ns[0] / np.array(ns)) ** 2, "k--", lw=0.8, label="n^-2 reference")
    axs[1].set_xlabel("points per axis n")
    axs[1].set_ylabel("relative error")
    axs[1].grid(True, which="both")
    axs[1].legend(fontsize=8)
    axs[1].set_title("Midpoint-rule convergence")
    fig.tight_layout()
    fig.savefig(out_dir / "wing_box_mass.png", dpi=150)
    print(f"saved {out_dir / 'wing_box_mass.png'}")


if __name__ == "__main__":
    main()
