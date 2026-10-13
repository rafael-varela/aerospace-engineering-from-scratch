"""P01.16 Eigenvalues & Eigenvectors: starter code.

Principal axes of a spacecraft inertia tensor, found with eigenvalue algorithms written from
scratch in NumPy (characteristic cubic, power iteration, shifted inverse iteration, cyclic Jacobi
rotations) and checked against ``numpy.linalg.eigh`` / ``numpy.linalg.eigvals`` as the oracle.

Conventions:
    * Inertia tensor about a point O for point masses m_k at r_k = (x, y, z) relative to O:
          I = sum_k m_k [ (r_k . r_k) I_3 - r_k r_k^T ]        [kg m**2]
      so the diagonal entries are moments of inertia and the off-diagonal entries are minus the
      products of inertia (I_xy = -sum m x y).
    * ``principal_axes`` returns the DCM C_pb from body to principal axes in the passive convention
      of P01.15: rows are the principal axes in body coordinates, det C = +1, and C I C^T is diagonal.

Complete every function marked NotImplementedError.

Run:  python main.py          (prints the satellite example, writes eigen_convergence.png)
Test: pytest curriculum/phase-01/lesson-16-eigenvalues-eigenvectors
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np


def satellite_point_masses() -> tuple[np.ndarray, np.ndarray]:
    """Return (masses [kg], positions [m]) of a small-satellite point-mass model in body axes.

    Body axes: x along the main antenna, y along the solar-array axis, z completing the triad.
    The bus is four 22 kg corner masses, the arrays two 7 kg masses on the y axis, plus a 4 kg
    antenna boom tip and a 16 kg propellant tank that sit off every symmetry plane, so the
    products of inertia are not zero and the principal axes are tilted from the body axes.
    """
    masses = np.array([22.0, 22.0, 22.0, 22.0, 7.0, 7.0, 4.0, 16.0])
    positions = np.array([
        [0.30, 0.20, 0.15], [0.30, -0.20, -0.15], [-0.30, 0.20, -0.15], [-0.30, -0.20, 0.15],
        [0.00, 1.60, 0.10], [0.00, -1.60, 0.10],
        [0.90, 0.25, -0.55],
        [-0.20, 0.10, 0.28],
    ])
    return masses, positions


# ---------------------------------------------------------------- the inertia tensor

def center_of_mass(masses: np.ndarray, positions: np.ndarray) -> np.ndarray:
    """Return the mass-weighted mean position, shape (3,)."""
    raise NotImplementedError


def inertia_tensor(masses: np.ndarray, positions: np.ndarray, about: np.ndarray | None = None) -> np.ndarray:
    """Return the 3x3 inertia tensor of point masses about the point ``about`` (default: the centre of mass).

    I = sum m [ (r.r) I_3 - r r^T ] with r measured from ``about``. Raise ``ValueError`` if any mass is
    not positive or the shapes do not match.
    """
    raise NotImplementedError


def moment_about_axis(masses: np.ndarray, positions: np.ndarray, axis: np.ndarray,
                      about: np.ndarray | None = None) -> float:
    """Return sum m |r_perp|**2, the moment of inertia about the line through ``about`` along ``axis``.

    Computed directly from the geometry, |r_perp|**2 = |r|**2 - (r.n)**2, so it is an independent check
    of n^T I n.
    """
    raise NotImplementedError


def rotate_tensor(I: np.ndarray, C: np.ndarray) -> np.ndarray:
    """Return C I C^T: the tensor expressed in the frame that ``C`` maps components into."""
    raise NotImplementedError


# ---------------------------------------------------------------- the characteristic polynomial

def characteristic_coefficients(A: np.ndarray) -> np.ndarray:
    """Return [1, -I1, I2, -I3], the coefficients of det(lambda I - A) for a 3x3 matrix.

    I1 = trace, I2 = sum of the three principal 2x2 minors, I3 = determinant (computed by cofactor
    expansion, not by numpy). These are the three rotation invariants of the tensor.
    """
    raise NotImplementedError


def symmetric_eigenvalues_3x3(A: np.ndarray) -> np.ndarray:
    """Return the three real eigenvalues of a symmetric 3x3 matrix in ascending order, in closed form.

    Shift by the mean to remove the quadratic term, then use the trigonometric solution of the depressed
    cubic: for B = A - (tr A / 3) I, lambda_k = (tr A)/3 + 2 sqrt(p) cos((phi + 2 pi k)/3) with
    p = tr(B^2)/6 and cos(3 phi) ... = det(B) / (2 p^{3/2}). Raise ``ValueError`` if A is not symmetric.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- iterative algorithms

def power_iteration(A: np.ndarray, x0: np.ndarray | None = None, tol: float = 1e-12,
                    max_iter: int = 10000) -> tuple[float, np.ndarray, int]:
    """Return (eigenvalue, unit eigenvector, iterations) for the eigenvalue of largest magnitude.

    Repeat x <- A x / |A x| until the change in the unit vector is below ``tol``; the eigenvalue is the
    Rayleigh quotient x^T A x. Starting vector defaults to ones. Raise ``RuntimeError`` if ``max_iter``
    is reached. Converges linearly with ratio |lambda_2 / lambda_1|.
    """
    raise NotImplementedError


def shifted_inverse_iteration(A: np.ndarray, shift: float, x0: np.ndarray | None = None,
                              tol: float = 1e-12, max_iter: int = 10000) -> tuple[float, np.ndarray, int]:
    """Return (eigenvalue, unit eigenvector, iterations) for the eigenvalue of ``A`` closest to ``shift``.

    Power iteration applied to (A - shift I)^{-1}, done with ``np.linalg.solve`` each step rather than an
    explicit inverse. Raise ``RuntimeError`` if ``max_iter`` is reached or the shift is exactly an
    eigenvalue (singular system).
    """
    raise NotImplementedError


def off_diagonal_norm(A: np.ndarray) -> float:
    """Return sqrt(sum of squares of the off-diagonal entries): Jacobi's measure of 'how far from diagonal'."""
    raise NotImplementedError


def jacobi_eigen(A: np.ndarray, tol: float = 1e-12, max_sweeps: int = 50,
                 history: list[float] | None = None) -> tuple[np.ndarray, np.ndarray, int]:
    """Return (eigenvalues ascending, V with eigenvectors as columns, sweeps) of a symmetric matrix by cyclic Jacobi.

    Each rotation zeroes one off-diagonal pair (p, q) exactly: with tau = (a_qq - a_pp) / (2 a_pq),
    t = sign(tau) / (|tau| + sqrt(1 + tau**2)), c = 1 / sqrt(1 + t**2), s = t c, the update is
    A <- J^T A J and V <- V J. A sweep visits every pair once; stop when ``off_diagonal_norm`` is below
    ``tol`` times the Frobenius norm (rounding limits ``tol`` to about n * 1e-16). If ``history`` is a list, the off-diagonal norm after each sweep is
    appended to it. Raise ``ValueError`` for a non-symmetric input, ``RuntimeError`` if ``max_sweeps`` is hit.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- principal axes and their meaning

def rayleigh_quotient(A: np.ndarray, x: np.ndarray) -> float:
    """Return x^T A x / x^T x. For an inertia tensor and a unit axis this is the moment of inertia about that axis."""
    raise NotImplementedError


def principal_axes(I: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (principal moments ascending, C_pb) for a symmetric inertia tensor ``I`` in body axes.

    C_pb is a proper rotation (det +1): row k is principal axis k in body coordinates, and
    rotate_tensor(I, C_pb) is diagonal with the principal moments on the diagonal. Built with ``jacobi_eigen``;
    if the eigenvector triad comes out left-handed, the sign of the third axis is flipped.
    """
    raise NotImplementedError


def triangle_inequalities_hold(moments: np.ndarray, rtol: float = 1e-12) -> bool:
    """Return True if every principal moment is at most the sum of the other two (and none is negative).

    Any distribution of positive mass satisfies I_1 + I_2 >= I_3 etc., with equality only for a flat
    (planar) body; a tensor that violates this cannot belong to a real body.
    """
    raise NotImplementedError


def spin_axis_stability(moments: np.ndarray, axis_index: int) -> dict[str, bool]:
    """Classify a pure spin about principal axis ``axis_index`` (0, 1 or 2 in the ordering of ``moments``).

    Returns {'rigid': ..., 'with_dissipation': ...}. For a torque-free rigid body, spin about the axis of
    largest or smallest moment is stable and spin about the intermediate axis is unstable. If any energy
    is dissipated (fuel slosh, flexing), only spin about the axis of largest moment is stable.
    Raise ``ValueError`` for a bad index or for a tensor with two equal moments touching that axis.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- the spectral theorem at work

def matrix_function_symmetric(A: np.ndarray, f: Callable[[np.ndarray], np.ndarray]) -> np.ndarray:
    """Return V f(Lambda) V^T for symmetric ``A`` = V Lambda V^T, with ``f`` applied entrywise to the eigenvalues.

    With f = np.exp this is the matrix exponential; with f = np.sqrt the symmetric square root.
    """
    raise NotImplementedError


def linear_system_decay_rates(A: np.ndarray) -> np.ndarray:
    """Return the real parts of the eigenvalues of a general (possibly non-symmetric) state matrix.

    For du/dt = A u every solution is a sum of exp(lambda t) terms, so the system is asymptotically stable
    exactly when every real part is negative. Non-symmetric eigenvalues are left to ``numpy.linalg.eigvals``:
    the lesson explains why.
    """
    raise NotImplementedError


def main(out_dir: str | Path = ".") -> None:
    """Print the satellite's principal axes and plot the convergence of power iteration and Jacobi sweeps."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    np.set_printoptions(precision=4, suppress=True)
    masses, positions = satellite_point_masses()
    I = inertia_tensor(masses, positions)
    print(f"total mass {masses.sum():.1f} kg, centre of mass {center_of_mass(masses, positions)} m")
    print("inertia tensor about the CM [kg m^2]:\n", I)
    coeffs = characteristic_coefficients(I)
    print("characteristic polynomial coefficients [1, -I1, I2, -I3]:", coeffs)
    print("numpy.poly agrees to", np.max(np.abs(coeffs - np.poly(I)) / np.abs(coeffs)))
    print("closed-form eigenvalues   :", symmetric_eigenvalues_3x3(I))
    moments, C_pb = principal_axes(I)
    print("Jacobi principal moments  :", moments)
    print("numpy.linalg.eigh         :", np.linalg.eigvalsh(I))
    print("C_pb (rows = principal axes in body coordinates):\n", C_pb)
    print("det C_pb =", np.linalg.det(C_pb))
    I_p = rotate_tensor(I, C_pb)
    print("C I C^T =\n", I_p, "\n   off-diagonal norm", off_diagonal_norm(I_p))
    deg = 180.0 / np.pi
    for k in range(3):
        nearest = int(np.argmax(np.abs(C_pb[k])))
        tilt = np.arccos(abs(C_pb[k, nearest])) * deg
        print(f"principal axis {k}: moment {moments[k]:9.3f} kg m^2, nearest body axis {'xyz'[nearest]},"
              f" tilted {tilt:6.2f} deg from it; check by geometry {moment_about_axis(masses, positions, C_pb[k]):9.3f}")
    print("triangle inequalities hold:", triangle_inequalities_hold(moments))
    for k, name in enumerate(["minor (smallest moment)", "intermediate", "major (largest moment)"]):
        print(f"spin about the {name} axis:", spin_axis_stability(moments, k))

    lam1, v1, n1 = power_iteration(I)
    lam_min, v_min, n_min = shifted_inverse_iteration(I, 0.0)
    print(f"\npower iteration: lambda_max = {lam1:.6f} in {n1} iterations;"
          f" inverse iteration: lambda_min = {lam_min:.6f} in {n_min} iterations")
    ratio = moments[1] / moments[2]
    print(f"predicted iterations for 1e-12 from the ratio lambda_2/lambda_1 = {ratio:.4f}:"
          f" {np.log(1e-12) / np.log(ratio):.0f}")

    # Convergence histories for the figure.
    x = np.ones(3) / np.sqrt(3.0)
    errs_power = []
    for _ in range(60):
        x = I @ x
        x /= np.linalg.norm(x)
        errs_power.append(abs(rayleigh_quotient(I, x) - moments[2]) / moments[2])
    hist: list[float] = []
    jacobi_eigen(I, tol=1e-15, history=hist)
    print(f"Jacobi on the 3x3 tensor: off-diagonal norm per sweep {np.array(hist)}")
    # A larger symmetric matrix (a 10-mass spring chain stiffness pattern) shows the quadratic convergence.
    n = 10
    K = 2.0 * np.eye(n) - np.eye(n, k=1) - np.eye(n, k=-1) + np.diag(np.linspace(0.0, 1.0, n))
    hist10: list[float] = []
    jacobi_eigen(K, tol=1e-13, history=hist10)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.8))
    ax1.semilogy(np.arange(1, 61), errs_power, "o-", ms=3, label="power iteration")
    ax1.semilogy(np.arange(1, 61), ratio ** (2 * np.arange(1, 61)), "--", label=r"$(\lambda_2/\lambda_1)^{2k}$")
    ax1.set_xlabel("iteration k")
    ax1.set_ylabel("relative error in $\\lambda_{max}$")
    ax1.set_ylim(1e-17, 1)
    ax1.grid(True, which="both", lw=0.4)
    ax1.legend()
    ax2.semilogy(np.arange(1, len(hist10) + 1), np.array(hist10) / np.linalg.norm(K), "s-")
    ax2.set_xlabel("Jacobi sweep (10x10 stiffness-pattern matrix)")
    ax2.set_ylabel("off-diagonal norm / |K|")
    ax2.grid(True, which="both", lw=0.4)
    fig.suptitle("Eigenvalue algorithms on the satellite inertia tensor")
    fig.tight_layout()
    fig.savefig(out_dir / "eigen_convergence.png", dpi=150)
    print(f"Jacobi on the 10x10 matrix: off-diagonal norm per sweep {np.array(hist10)}")
    print(f"saved {out_dir / 'eigen_convergence.png'}")

    # A non-symmetric example: the damped oscillator x'' + 2 zeta w x' + w^2 x = 0 as du/dt = A u.
    w, zeta = 2.0, 0.1
    A = np.array([[0.0, 1.0], [-w * w, -2.0 * zeta * w]])
    print("damped-oscillator decay rates:", linear_system_decay_rates(A), " (expected", -zeta * w, ")")


if __name__ == "__main__":
    main()
