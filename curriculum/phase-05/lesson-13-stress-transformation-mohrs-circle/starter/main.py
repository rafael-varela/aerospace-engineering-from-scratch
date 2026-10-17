"""P05.13 Stress Transformation & Mohr's Circle: starter code.

A pressurised fuselage skin carries hoop stress, axial stress and shear at the same point. Which way
does the skin want to crack, how large is the largest normal stress, and how large the largest shear?
The answers come from rotating the stress tensor. This module provides:

* the 2-D transformation equations and their geometric form, Mohr's circle: principal stresses and
  their angle, the maximum in-plane shear, and the absolute maximum shear once the third (out-of-plane)
  principal stress is remembered;
* 3-D stress: tensor rotation sigma' = R sigma R^T, the three invariants, the principal stresses from
  the characteristic cubic solved in closed (trigonometric) form, principal directions by cross
  products, the traction on any plane and the three Mohr circles that bound every (sigma_n, tau_n);
* the strain-gauge rosette: three normal strains -> the strain state -> plane-stress Hooke's law.

The independent check is numpy.linalg.eigh: the closed-form cubic and the Mohr-circle formulas must
agree with a numerical symmetric eigen-solver on any tensor.

Sign convention: tension positive. tau_xy is the shear on the face whose normal is +x, acting in +y.
theta is measured anticlockwise from the x axis to the rotated x' axis. On Mohr's circle we plot
(sigma, tau) with tau positive upward for the x face, so a rotation of theta in the body is a rotation
of 2 theta on the circle, in the *opposite* sense (clockwise) with this plotting choice.

Units: Pa (the printed numbers are MPa).

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes mohr_circle.png)
Test: pytest curriculum/phase-05/lesson-13-stress-transformation-mohrs-circle
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

E_2024 = 73.1e9            # Pa, typical Young's modulus of 2024-T3 (ASM Aerospace Specification Metals data sheet)
NU_2024 = 0.33             # typical Poisson's ratio of 2024-T3


# ---------------------------------------------------------------- 2-D

def transform_stress_2d(sx: float, sy: float, txy: float, theta: np.ndarray) -> tuple:
    """Stresses on axes x', y' rotated anticlockwise by theta (rad) from x, y:

        sx' = (sx + sy)/2 + (sx - sy)/2 cos 2t + txy sin 2t
        sy' = (sx + sy)/2 - (sx - sy)/2 cos 2t - txy sin 2t
        tx'y' = -(sx - sy)/2 sin 2t + txy cos 2t

    Returns (sx', sy', tx'y'); theta may be an array."""
    raise NotImplementedError


def mohr_circle_2d(sx: float, sy: float, txy: float) -> tuple[float, float]:
    """Centre C = (sx + sy)/2 and radius R = sqrt(((sx - sy)/2)^2 + txy^2) of Mohr's circle."""
    raise NotImplementedError


def principal_stresses_2d(sx: float, sy: float, txy: float) -> tuple[float, float, float]:
    """In-plane principal stresses s1 >= s2 and the angle theta_p (rad, in (-pi/2, pi/2]) from x
    to the direction of s1: theta_p = atan2(2 txy, sx - sy) / 2."""
    raise NotImplementedError


def max_inplane_shear(sx: float, sy: float, txy: float) -> tuple[float, float]:
    """Maximum in-plane shear stress R and the angle theta_s = theta_p - pi/4 of a plane carrying it
    (on that plane the normal stress is the circle's centre)."""
    raise NotImplementedError


def plane_stress_abs_max_shear(sx: float, sy: float, txy: float) -> float:
    """Absolute maximum shear of a plane-stress state, remembering sigma_z = 0 is a principal stress:
    (max(s1, s2, 0) - min(s1, s2, 0)) / 2. When s1 and s2 have the same sign it exceeds the in-plane R."""
    raise NotImplementedError


# ---------------------------------------------------------------- 3-D

def stress_tensor(sx: float, sy: float, sz: float, txy: float, tyz: float, txz: float) -> np.ndarray:
    """Symmetric 3 x 3 Cauchy stress tensor from its six components."""
    return np.array([[sx, txy, txz], [txy, sy, tyz], [txz, tyz, sz]], dtype=float)


def rotation_z(theta: float) -> np.ndarray:
    """Matrix whose rows are the new axes x', y', z' after rotating x, y anticlockwise by theta about z."""
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])


def transform_tensor(sigma: np.ndarray, R: np.ndarray) -> np.ndarray:
    """Components of the tensor on new axes: sigma' = R sigma R^T, the rows of R being the new axes."""
    raise NotImplementedError


def invariants(sigma: np.ndarray) -> tuple[float, float, float]:
    """I1 = tr sigma, I2 = sum of principal 2 x 2 minors, I3 = det sigma; the characteristic
    equation is s^3 - I1 s^2 + I2 s - I3 = 0, unchanged by any rotation."""
    raise NotImplementedError


def principal_stresses_3d(sigma: np.ndarray) -> np.ndarray:
    """Principal stresses s1 >= s2 >= s3 from the characteristic cubic in closed form.

    Remove the mean p = I1/3: the deviator has invariants J2 = (I1^2 - 3 I2)/3 and
    J3 = det(sigma - p E). With r = 2 sqrt(J2/3) and cos 3a = (J3/2)(3/J2)^{3/2} (clipped to [-1, 1]),
    the roots are p + r cos(a - 2 pi k / 3), k = 0, 1, 2. No iterative solver is used."""
    raise NotImplementedError


def principal_direction(sigma: np.ndarray, s: float) -> np.ndarray:
    """Unit vector n with (sigma - s E) n = 0 for a simple (non-repeated) principal stress s: the
    largest cross product of two rows of sigma - s E (each row is perpendicular to n)."""
    raise NotImplementedError


def absolute_max_shear(principal: np.ndarray) -> float:
    """Absolute maximum shear (s_max - s_min) / 2, on planes at 45 degrees to the s1 and s3 directions."""
    raise NotImplementedError


def traction_on_plane(sigma: np.ndarray, n: np.ndarray) -> tuple[float, float]:
    """Normal and shear stress on the plane with unit normal n: t = sigma n, sigma_n = n . t,
    tau_n = sqrt(|t|^2 - sigma_n^2)."""
    raise NotImplementedError


def mohr_circles_3d(principal: np.ndarray) -> list[tuple[float, float]]:
    """(centre, radius) of the three Mohr circles through s1, s2, s3 (sorted descending):
    circles on [s2, s1], [s3, s2] and the big one on [s3, s1]."""
    raise NotImplementedError


def in_mohr_region(sn: float, tn: float, principal: np.ndarray, tol: float = 1e-9) -> bool:
    """True when (sn, tn) lies inside the big circle and outside both small circles (the admissible
    region for tractions on all planes through a point)."""
    raise NotImplementedError


# ---------------------------------------------------------------- strain rosette

def rosette_strains(eps: np.ndarray, angles_deg: np.ndarray = (0.0, 45.0, 90.0)) -> tuple[float, float, float]:
    """Strain state (eps_x, eps_y, gamma_xy) from three gauge readings at the given angles, solving
    eps_theta = eps_x cos^2 t + eps_y sin^2 t + gamma_xy sin t cos t for the three unknowns."""
    raise NotImplementedError


def plane_stress_from_strain(ex: float, ey: float, gxy: float, E: float, nu: float) -> tuple[float, float, float]:
    """Plane-stress Hooke's law: sx = E/(1 - nu^2)(ex + nu ey), sy = E/(1 - nu^2)(ey + nu ex),
    txy = G gxy with G = E / (2 (1 + nu))."""
    raise NotImplementedError


# ---------------------------------------------------------------- the fuselage skin

def fuselage_skin_stress(dp: float, r: float, t: float, sigma_bending: float = 0.0,
                         tau: float = 0.0) -> tuple[float, float, float]:
    """Membrane stresses in a thin cylindrical fuselage skin (x axial, y hoop): axial dp r / (2t)
    plus a fuselage-bending stress, hoop dp r / t, and a shear from torsion or transverse shear."""
    raise NotImplementedError


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save mohr_circle.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    MPa = 1e6
    dp, r, t = 60e3, 2.0, 1.6e-3
    print(f"fuselage r = {r} m, skin t = {t * 1e3} mm, dp = {dp / 1e3:.0f} kPa: hoop {dp * r / t / MPa:.1f} MPa, "
          f"axial {dp * r / (2 * t) / MPa:.1f} MPa")
    cases = {"side panel": (0.0, 30e6), "crown panel": (40e6, 8e6)}
    for name, (sb, tau) in cases.items():
        sx, sy, txy = fuselage_skin_stress(dp, r, t, sb, tau)
        s1, s2, tp = principal_stresses_2d(sx, sy, txy)
        C, R = mohr_circle_2d(sx, sy, txy)
        print(f"\n{name}: sx = {sx / MPa:.1f}, sy = {sy / MPa:.1f}, txy = {txy / MPa:.1f} MPa")
        print(f"  Mohr circle centre {C / MPa:.2f}, radius {R / MPa:.2f} MPa; s1 = {s1 / MPa:.2f}, s2 = {s2 / MPa:.2f} MPa "
              f"at theta_p = {np.degrees(tp):.2f} deg")
        print(f"  max in-plane shear {R / MPa:.2f} MPa; absolute max shear (with sz = 0) "
              f"{plane_stress_abs_max_shear(sx, sy, txy) / MPa:.2f} MPa")
        sig = stress_tensor(sx, sy, 0.0, txy, 0.0, 0.0)
        ev = np.linalg.eigvalsh(sig)[::-1]
        print(f"  eigh check: {np.array2string(ev / MPa, precision=4)} MPa")

    # 3-D: an illustrative stress state at a fitting (synthetic numbers, not a measurement)
    sig = stress_tensor(180e6, -40e6, 90e6, 60e6, 30e6, -25e6)
    I1, I2, I3 = invariants(sig)
    p = principal_stresses_3d(sig)
    ev, vec = np.linalg.eigh(sig)
    print(f"\n3-D fitting point (synthetic): I1 = {I1 / MPa:.1f} MPa, I2 = {I2 / MPa ** 2:.1f} MPa^2, "
          f"I3 = {I3 / MPa ** 3:.1f} MPa^3")
    print(f"  closed-form cubic: {np.array2string(p / MPa, precision=6)} MPa")
    print(f"  numpy eigh:        {np.array2string(ev[::-1] / MPa, precision=6)} MPa; "
          f"max difference {np.max(np.abs(p - ev[::-1])):.1e} Pa")
    n1 = principal_direction(sig, p[0])
    print(f"  s1 direction by cross products {np.array2string(n1, precision=4)}, "
          f"|n1 . v1(eigh)| = {abs(n1 @ vec[:, -1]):.12f}")
    print(f"  absolute max shear {absolute_max_shear(p) / MPa:.2f} MPa")
    R = rotation_z(np.radians(30.0))
    print(f"  invariants after a 30 deg rotation: I1 {invariants(transform_tensor(sig, R))[0] / MPa:.1f}, "
          f"I2 {invariants(transform_tensor(sig, R))[1] / MPa ** 2:.1f}")

    # rosette (synthetic readings)
    eps = np.array([1200e-6, 900e-6, -300e-6])
    ex, ey, g = rosette_strains(eps)
    sx, sy, txy = plane_stress_from_strain(ex, ey, g, E_2024, NU_2024)
    s1, s2, tp = principal_stresses_2d(sx, sy, txy)
    print(f"\nrectangular rosette 0/45/90 (synthetic readings {', '.join(f'{e * 1e6:.0f}' for e in eps)} microstrain): "
          f"ex = {ex * 1e6:.0f}, ey = {ey * 1e6:.0f}, gxy = {g * 1e6:.0f} microstrain")
    print(f"  2024-T3: sx = {sx / MPa:.1f}, sy = {sy / MPa:.1f}, txy = {txy / MPa:.1f} MPa; "
          f"s1 = {s1 / MPa:.1f}, s2 = {s2 / MPa:.1f} MPa at {np.degrees(tp):.1f} deg")

    # figure: Mohr circle of the side panel, and 3-D circles with random planes
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    sx, sy, txy = fuselage_skin_stress(dp, r, t, 0.0, 30e6)
    C, Rc = mohr_circle_2d(sx, sy, txy)
    th = np.linspace(0, 2 * np.pi, 361)
    ax[0].plot((C + Rc * np.cos(th)) / MPa, Rc * np.sin(th) / MPa)
    ax[0].plot([sx / MPa, sy / MPa], [txy / MPa, -txy / MPa], "o-", label="x face, y face")
    ang = np.linspace(0, np.pi, 13)
    a, b_, c = transform_stress_2d(sx, sy, txy, ang)
    ax[0].plot(a / MPa, c / MPa, "k.", label="rotated x' face, 15 deg steps")
    ax[0].axhline(0, c="k", lw=0.5)
    ax[0].set_aspect("equal")
    ax[0].set_xlabel("sigma (MPa)")
    ax[0].set_ylabel("tau (MPa)")
    ax[0].set_title("fuselage side panel")
    ax[0].legend(fontsize=8)
    rng = np.random.default_rng(13)
    sig = stress_tensor(180e6, -40e6, 90e6, 60e6, 30e6, -25e6)
    p = principal_stresses_3d(sig)
    nn = rng.normal(size=(1500, 3))
    pts = np.array([traction_on_plane(sig, v) for v in nn])
    ax[1].plot(pts[:, 0] / MPa, pts[:, 1] / MPa, ".", ms=2, alpha=0.5, label="random planes")
    for cc, rr in mohr_circles_3d(p):
        ax[1].plot((cc + rr * np.cos(th[:181])) / MPa, rr * np.sin(th[:181]) / MPa, "k")
    ax[1].set_aspect("equal")
    ax[1].set_xlabel("sigma_n (MPa)")
    ax[1].set_ylabel("tau_n (MPa)")
    ax[1].set_title("3-D Mohr circles (synthetic fitting stress)")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "mohr_circle.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'mohr_circle.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
