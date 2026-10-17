"""P02.14 Tensor Notation & Index Gymnastics: starter code.

Index notation made executable: the Kronecker delta and Levi-Civita symbol as arrays, every contraction
written as an ``numpy.einsum`` string that mirrors the index equation, and the transformation laws of
vectors, second-order tensors (stress, strain) and the fourth-order stiffness tensor under a rotation of
axes. Applied to the skin of a pressurised fuselage under a fin-induced torque: stresses on an inclined
riveted seam, principal stresses, von Mises stress and margin against yield.

Conventions:
    * Summation convention: a repeated index in a term is summed over 1..3 (0..2 in code).
    * Q is the direction-cosine matrix of P01.15 (passive): Q[i, j] = e'_i . e_j, rows are the NEW basis
      vectors in OLD components. Components transform as v'_i = Q_ij v_j, sigma'_ij = Q_ip Q_jq sigma_pq,
      C'_ijkl = Q_ia Q_jb Q_kc Q_ld C_abcd.
    * Stress sigma_ij: force in direction i on the face whose outward normal is j (symmetric, so the order
      does not matter). Cauchy's formula: traction t_i = sigma_ij n_j. Tension positive. Units: Pa.
    * Fuselage axes: x_1 along the fuselage (longitudinal), x_2 around it (hoop), x_3 through the skin
      (radial); thin-walled, so sigma_33 = sigma_13 = sigma_23 = 0 (plane stress).

Run:  python main.py          (prints the fuselage-seam analysis and identity checks, writes stress_rotation.png)
Test: pytest curriculum/phase-02/lesson-14-tensor-notation
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


def fuselage_case() -> dict[str, float]:
    """Return the lesson's fuselage skin case (plausible transport-aircraft numbers, not a specific type).

    Cabin differential pressure p = 0.060 MPa, radius r = 2.0 m, skin thickness t = 1.6 mm, a torque
    T = 0.80 MN m from a fin side load, 2024-T3 sheet with yield stress 325 MPa, and a riveted lap seam
    inclined 30 degrees to the longitudinal axis.
    """
    return {"p": 0.060e6, "r": 2.0, "t": 1.6e-3, "torque": 0.80e6, "sigma_yield": 325e6, "seam_angle_deg": 30.0}


# ---------------------------------------------------------------- the two symbols

def kronecker_delta() -> np.ndarray:
    """Return delta_ij as a 3x3 array (1 if i == j else 0)."""
    raise NotImplementedError


def levi_civita() -> np.ndarray:
    """Return epsilon_ijk as a 3x3x3 array: +1 for even permutations of (0,1,2), -1 for odd, 0 otherwise."""
    raise NotImplementedError


def epsilon_delta_identity_residual() -> float:
    """Return max |eps_ijk eps_imn - (delta_jm delta_kn - delta_jn delta_km)| over all j, k, m, n (zero if the identity holds)."""
    raise NotImplementedError


# ---------------------------------------------------------------- contractions as einsum strings

def dot(a: np.ndarray, b: np.ndarray) -> float:
    """Return a_i b_i."""
    raise NotImplementedError


def cross(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Return c_i = eps_ijk a_j b_k."""
    raise NotImplementedError


def determinant(A: np.ndarray) -> float:
    """Return det A = eps_ijk A_1i A_2j A_3k (rows 0, 1, 2 in code)."""
    raise NotImplementedError


def trace(A: np.ndarray) -> float:
    """Return A_ii."""
    raise NotImplementedError


def double_dot(A: np.ndarray, B: np.ndarray) -> float:
    """Return A_ij B_ij (the Frobenius inner product)."""
    raise NotImplementedError


def matrix_product(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Return (AB)_ik = A_ij B_jk."""
    raise NotImplementedError


def symmetric_part(A: np.ndarray) -> np.ndarray:
    """Return A_(ij) = (A_ij + A_ji) / 2."""
    raise NotImplementedError


def antisymmetric_part(A: np.ndarray) -> np.ndarray:
    """Return A_[ij] = (A_ij - A_ji) / 2."""
    raise NotImplementedError


def axial_vector(W: np.ndarray) -> np.ndarray:
    """Return w_i = -(1/2) eps_ijk W_jk, the vector dual to an antisymmetric W (so that W_jk v_k = (w x v)_j)."""
    raise NotImplementedError


# ---------------------------------------------------------------- rotations and transformation laws

def rotation_matrix(axis: np.ndarray, angle: float) -> np.ndarray:
    """Return the passive DCM Q for a rotation of the AXES by ``angle`` about the unit vector ``axis``.

    Rodrigues: Q_ij = cos(a) delta_ij + (1 - cos a) n_i n_j + sin(a) eps_ijk n_k. (Rotating the axes by +a
    is the same as rotating a vector's components by -a.) Raise ``ValueError`` for a zero axis.
    """
    raise NotImplementedError


def is_orthogonal(Q: np.ndarray, tol: float = 1e-12) -> bool:
    """Return True if Q_ik Q_jk = delta_ij and det Q = +1 (a proper rotation)."""
    raise NotImplementedError


def transform_vector(v: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Return v'_i = Q_ij v_j: the components of the same vector in the rotated axes."""
    raise NotImplementedError


def transform_tensor2(S: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Return S'_ij = Q_ip Q_jq S_pq: a second-order tensor in the rotated axes (equals Q S Q^T)."""
    raise NotImplementedError


def transform_tensor4(C: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Return C'_ijkl = Q_ia Q_jb Q_kc Q_ld C_abcd: a fourth-order tensor in the rotated axes."""
    raise NotImplementedError


# ---------------------------------------------------------------- the stress tensor

def traction(sigma: np.ndarray, n: np.ndarray) -> np.ndarray:
    """Return Cauchy's traction t_i = sigma_ij n_j on the plane with unit normal n. Raise ``ValueError`` if sigma is not symmetric."""
    raise NotImplementedError


def normal_and_shear_on_plane(sigma: np.ndarray, n: np.ndarray) -> tuple[float, float]:
    """Return (sigma_n, tau): normal stress t_i n_i and shear magnitude |t - sigma_n n| on the plane with normal n."""
    raise NotImplementedError


def stress_invariants(sigma: np.ndarray) -> tuple[float, float, float]:
    """Return (I1, I2, I3) = (sigma_ii, (sigma_ii sigma_jj - sigma_ij sigma_ji) / 2, det sigma): unchanged by any rotation."""
    raise NotImplementedError


def principal_stresses(sigma: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (values, directions): principal stresses in descending order and the unit normals as columns.

    sigma_ij n_j = s n_i: the eigenproblem of the symmetric stress tensor (``numpy.linalg.eigh``). The
    direction triad is made right-handed (det = +1).
    """
    raise NotImplementedError


def hydrostatic_and_deviatoric(sigma: np.ndarray) -> tuple[float, np.ndarray]:
    """Return (p, s): the mean stress p = sigma_kk / 3 and the deviator s_ij = sigma_ij - p delta_ij (traceless)."""
    raise NotImplementedError


def von_mises(sigma: np.ndarray) -> float:
    """Return sigma_vm = sqrt(3/2 s_ij s_ij) with s the deviator; equals sqrt(((s1-s2)^2 + (s2-s3)^2 + (s3-s1)^2)/2)."""
    raise NotImplementedError


def max_shear(sigma: np.ndarray) -> float:
    """Return tau_max = (sigma_1 - sigma_3) / 2, the largest shear stress on any plane (Tresca)."""
    raise NotImplementedError


def mohr_2d(sx: float, sy: float, txy: float, theta: float) -> tuple[float, float, float]:
    """Return (sx', sy', tx'y') on axes rotated by theta (counter-clockwise, about x3) for plane stress.

    sx' = (sx + sy)/2 + (sx - sy)/2 cos 2theta + txy sin 2theta,
    sy' = (sx + sy)/2 - (sx - sy)/2 cos 2theta - txy sin 2theta,
    tx'y' = -(sx - sy)/2 sin 2theta + txy cos 2theta.
    These are the i, j = 1, 2 entries of Q_ip Q_jq sigma_pq with Q the rotation about x3; the lesson
    derives them from the tensor law and the test checks that they agree.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- strain and Hooke's law

def strain_from_displacement_gradient(H: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (epsilon, omega): the small strain eps_ij = (u_i,j + u_j,i)/2 and rotation omega_ij = (u_i,j - u_j,i)/2 for H_ij = u_i,j."""
    raise NotImplementedError


def isotropic_stiffness(E: float, nu: float) -> np.ndarray:
    """Return C_ijkl = lambda delta_ij delta_kl + mu (delta_ik delta_jl + delta_il delta_jk) for an isotropic solid.

    lambda = E nu / ((1 + nu)(1 - 2 nu)), mu = E / (2 (1 + nu)). Raise ``ValueError`` unless E > 0 and -1 < nu < 1/2.
    """
    raise NotImplementedError


def hooke_stress(C: np.ndarray, eps: np.ndarray) -> np.ndarray:
    """Return sigma_ij = C_ijkl eps_kl."""
    raise NotImplementedError


def isotropic_strain(sigma: np.ndarray, E: float, nu: float) -> np.ndarray:
    """Return eps_ij = ((1 + nu) sigma_ij - nu sigma_kk delta_ij) / E, the inverse of Hooke's law for an isotropic solid."""
    raise NotImplementedError


# ---------------------------------------------------------------- the fuselage skin

def fuselage_skin_stress(p: float, r: float, t: float, torque: float) -> np.ndarray:
    """Return the plane-stress tensor [Pa] in the skin of a thin pressurised cylinder carrying a torque.

    Axes: 1 longitudinal, 2 hoop, 3 radial. sigma_11 = p r / (2 t), sigma_22 = p r / t (membrane equilibrium),
    sigma_12 = T / (2 pi r^2 t) (Bredt's formula for a closed thin tube), all others zero.
    Raise ``ValueError`` unless r > 0 and 0 < t < 0.1 r (thin-walled).
    """
    raise NotImplementedError


def seam_stresses(sigma: np.ndarray, seam_angle: float) -> tuple[float, float]:
    """Return (normal, shear) stress [Pa] across a seam whose line makes ``seam_angle`` [rad] with axis 1, in the 1-2 plane.

    The seam's normal is n = (-sin a, cos a, 0); the normal stress opens the seam (tension positive), the
    shear slides it. Rivets along the seam carry both.
    """
    raise NotImplementedError


def yield_margin(sigma: np.ndarray, sigma_yield: float) -> float:
    """Return the margin of safety sigma_yield / sigma_vm - 1 under the von Mises criterion (negative means yield)."""
    raise NotImplementedError


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Check the identities, analyse the fuselage seam, and plot stresses against seam angle."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    rng = np.random.default_rng(14)
    np.set_printoptions(precision=4, suppress=True)
    print(f"epsilon-delta identity residual: {epsilon_delta_identity_residual():.1e}")
    a, b = rng.normal(size=3), rng.normal(size=3)
    print(f"cross via eps_ijk a_j b_k vs numpy: {np.max(np.abs(cross(a, b) - np.cross(a, b))):.1e};"
          f" a . (a x b) = {dot(a, cross(a, b)):.1e}")
    A = rng.normal(size=(3, 3))
    print(f"det via eps_ijk A_1i A_2j A_3k vs numpy: {abs(determinant(A) - np.linalg.det(A)):.1e}")
    Q = rotation_matrix(rng.normal(size=3), 0.9)
    print(f"Rodrigues DCM is a proper rotation: {is_orthogonal(Q)}; Q_ik Q_jk - delta_ij max {np.max(np.abs(Q @ Q.T - np.eye(3))):.1e}")

    case = fuselage_case()
    sig = fuselage_skin_stress(case["p"], case["r"], case["t"], case["torque"])
    MPa = 1e-6
    print(f"\nFUSELAGE SKIN (p = {case['p'] * MPa:.3f} MPa, r = {case['r']} m, t = {case['t'] * 1e3} mm, T = {case['torque'] * MPa:.2f} MN m)")
    print("stress tensor [MPa], axes 1 longitudinal, 2 hoop, 3 radial:\n", sig * MPa)
    ang = np.deg2rad(case["seam_angle_deg"])
    sn, tau = seam_stresses(sig, ang)
    print(f"seam at {case['seam_angle_deg']:.0f} deg to the axis: normal (opening) stress {sn * MPa:.2f} MPa, shear along the seam {tau * MPa:.2f} MPa")
    Qz = rotation_matrix([0, 0, 1], ang)
    sig_r = transform_tensor2(sig, Qz)
    print(f"same by rotating the tensor (Q_ip Q_jq sigma_pq) to seam axes:\n", sig_r * MPa)
    m2 = mohr_2d(sig[0, 0], sig[1, 1], sig[0, 1], ang)
    print(f"same by the 2-D Mohr formulas: sx' = {m2[0] * MPa:.2f}, sy' = {m2[1] * MPa:.2f}, tx'y' = {m2[2] * MPa:.2f} MPa")
    I = stress_invariants(sig)
    Ir = stress_invariants(sig_r)
    print(f"invariants I1, I2, I3 [MPa, MPa^2, MPa^3]: {I[0] * MPa:.3f}, {I[1] * MPa ** 2:.3f}, {I[2] * MPa ** 3:.3f}; after rotation {Ir[0] * MPa:.3f}, {Ir[1] * MPa ** 2:.3f}, {Ir[2] * MPa ** 3:.3f}")
    vals, dirs = principal_stresses(sig)
    print(f"principal stresses [MPa]: {vals * MPa}; sigma_1 direction {dirs[:, 0]} ({np.degrees(np.arctan2(dirs[1, 0], dirs[0, 0])):.1f} deg from the axis)")
    print(f"max shear {max_shear(sig) * MPa:.2f} MPa on planes at 45 deg to the principal axes")
    p_mean, s = hydrostatic_and_deviatoric(sig)
    print(f"mean stress {p_mean * MPa:.2f} MPa; von Mises {von_mises(sig) * MPa:.2f} MPa; yield {case['sigma_yield'] * MPa:.0f} MPa;"
          f" margin of safety {yield_margin(sig, case['sigma_yield']):.2f}")
    sig_p = fuselage_skin_stress(case["p"], case["r"], case["t"], 0.0)
    print(f"pressure alone: von Mises {von_mises(sig_p) * MPa:.2f} MPa (hoop {sig_p[1, 1] * MPa:.1f}, long {sig_p[0, 0] * MPa:.1f}), margin {yield_margin(sig_p, case['sigma_yield']):.2f}")

    E, nu = 73e9, 0.33
    C = isotropic_stiffness(E, nu)
    eps = isotropic_strain(sig, E, nu)
    print(f"\n2024-T3 (E = {E * 1e-9:.0f} GPa, nu = {nu}): strains from the inverse law, eps_11 = {eps[0, 0] * 1e6:.0f}, eps_22 = {eps[1, 1] * 1e6:.0f},"
          f" eps_33 = {eps[2, 2] * 1e6:.0f}, eps_12 = {eps[0, 1] * 1e6:.0f} microstrain; C_ijkl eps_kl recovers sigma to {np.max(np.abs(hooke_stress(C, eps) - sig)) * MPa:.1e} MPa")
    Cr = transform_tensor4(C, Q)
    print(f"isotropic stiffness under a random rotation: max |C' - C| / max |C| = {np.max(np.abs(Cr - C)) / np.max(np.abs(C)):.1e}")
    C_ortho = C.copy()
    C_ortho[0, 0, 0, 0] *= 2.0                                  # stiffen one direction: no longer isotropic
    print(f"after stiffening C_1111 only: max |C' - C| / max |C| = {np.max(np.abs(transform_tensor4(C_ortho, Q) - C_ortho)) / np.max(np.abs(C_ortho)):.2f}")

    # Figure: seam normal and shear stress against seam angle, and Mohr's circle.
    angles = np.linspace(0.0, np.pi, 361)
    sn_a = np.array([seam_stresses(sig, a_)[0] for a_ in angles]) * MPa
    tau_a = np.array([seam_stresses(sig, a_)[1] for a_ in angles]) * MPa
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    ax1.plot(np.degrees(angles), sn_a, label="normal stress across the seam")
    ax1.plot(np.degrees(angles), tau_a, label="shear along the seam")
    ax1.axvline(case["seam_angle_deg"], color="k", lw=0.8, ls="--")
    ax1.set_xlabel("seam angle from the longitudinal axis [deg]"); ax1.set_ylabel("stress [MPa]")
    ax1.set_title("stresses on a riveted seam"); ax1.legend(fontsize=8); ax1.grid(True, lw=0.4)
    centre, radius = 0.5 * (sig[0, 0] + sig[1, 1]) * MPa, np.hypot(0.5 * (sig[0, 0] - sig[1, 1]), sig[0, 1]) * MPa
    th = np.linspace(0, 2 * np.pi, 361)
    ax2.plot(centre + radius * np.cos(th), radius * np.sin(th))
    ax2.plot([sig[0, 0] * MPa, sig[1, 1] * MPa], [sig[0, 1] * MPa, -sig[0, 1] * MPa], "o-", label="axes 1-2")
    ax2.plot([m2[0] * MPa, m2[1] * MPa], [m2[2] * MPa, -m2[2] * MPa], "s-", label=f"seam axes ({case['seam_angle_deg']:.0f} deg)")
    ax2.plot(vals[:2] * MPa, [0, 0], "k^", label="principal")
    ax2.set_aspect("equal"); ax2.set_xlabel("normal stress [MPa]"); ax2.set_ylabel("shear stress [MPa]")
    ax2.set_title("Mohr's circle (plane stress)"); ax2.legend(fontsize=8); ax2.grid(True, lw=0.4)
    fig.suptitle("P02.14  Stress tensor transformations: fuselage skin")
    fig.tight_layout()
    fig.savefig(out_dir / "stress_rotation.png", dpi=150)
    print(f"\nsaved {out_dir / 'stress_rotation.png'}")


if __name__ == "__main__":
    main()
