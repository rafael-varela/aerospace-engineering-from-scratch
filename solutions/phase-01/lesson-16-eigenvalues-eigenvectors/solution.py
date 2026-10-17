"""P01.16 Eigenvalues & Eigenvectors: reference solution.

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

Run:  python solution.py          (prints the satellite example, writes eigen_convergence.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-01/lesson-16-eigenvalues-eigenvectors
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
    m = np.asarray(masses, dtype=float)
    r = np.asarray(positions, dtype=float)
    return (m[:, None] * r).sum(axis=0) / m.sum()


def inertia_tensor(masses: np.ndarray, positions: np.ndarray, about: np.ndarray | None = None) -> np.ndarray:
    """Return the 3x3 inertia tensor of point masses about the point ``about`` (default: the centre of mass).

    I = sum m [ (r.r) I_3 - r r^T ] with r measured from ``about``. Raise ``ValueError`` if any mass is
    not positive or the shapes do not match.
    """
    m = np.asarray(masses, dtype=float)
    r = np.asarray(positions, dtype=float)
    if r.ndim != 2 or r.shape[1] != 3 or r.shape[0] != m.shape[0]:
        raise ValueError("positions must have shape (N, 3) matching the N masses")
    if np.any(m <= 0.0):
        raise ValueError("every mass must be positive")
    origin = center_of_mass(m, r) if about is None else np.asarray(about, dtype=float)
    d = r - origin
    I = np.zeros((3, 3))
    for mk, rk in zip(m, d):
        I += mk * (np.dot(rk, rk) * np.eye(3) - np.outer(rk, rk))
    return I


def moment_about_axis(masses: np.ndarray, positions: np.ndarray, axis: np.ndarray,
                      about: np.ndarray | None = None) -> float:
    """Return sum m |r_perp|**2, the moment of inertia about the line through ``about`` along ``axis``.

    Computed directly from the geometry, |r_perp|**2 = |r|**2 - (r.n)**2, so it is an independent check
    of n^T I n.
    """
    m = np.asarray(masses, dtype=float)
    r = np.asarray(positions, dtype=float)
    n = np.asarray(axis, dtype=float)
    n = n / np.linalg.norm(n)
    origin = center_of_mass(m, r) if about is None else np.asarray(about, dtype=float)
    d = r - origin
    return float(np.sum(m * (np.sum(d * d, axis=1) - (d @ n) ** 2)))


def rotate_tensor(I: np.ndarray, C: np.ndarray) -> np.ndarray:
    """Return C I C^T: the tensor expressed in the frame that ``C`` maps components into."""
    C = np.asarray(C, dtype=float)
    return C @ np.asarray(I, dtype=float) @ C.T


# ---------------------------------------------------------------- the characteristic polynomial

def characteristic_coefficients(A: np.ndarray) -> np.ndarray:
    """Return [1, -I1, I2, -I3], the coefficients of det(lambda I - A) for a 3x3 matrix.

    I1 = trace, I2 = sum of the three principal 2x2 minors, I3 = determinant (computed by cofactor
    expansion, not by numpy). These are the three rotation invariants of the tensor.
    """
    A = np.asarray(A, dtype=float)
    if A.shape != (3, 3):
        raise ValueError("characteristic_coefficients is written for 3x3 matrices")
    i1 = A[0, 0] + A[1, 1] + A[2, 2]
    i2 = ((A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]) + (A[0, 0] * A[2, 2] - A[0, 2] * A[2, 0])
          + (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]))
    i3 = (A[0, 0] * (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1])
          - A[0, 1] * (A[1, 0] * A[2, 2] - A[1, 2] * A[2, 0])
          + A[0, 2] * (A[1, 0] * A[2, 1] - A[1, 1] * A[2, 0]))
    return np.array([1.0, -i1, i2, -i3])


def symmetric_eigenvalues_3x3(A: np.ndarray) -> np.ndarray:
    """Return the three real eigenvalues of a symmetric 3x3 matrix in ascending order, in closed form.

    Shift by the mean to remove the quadratic term, then use the trigonometric solution of the depressed
    cubic: for B = A - (tr A / 3) I, lambda_k = (tr A)/3 + 2 sqrt(p) cos((phi + 2 pi k)/3) with
    p = tr(B^2)/6 and cos(3 phi) ... = det(B) / (2 p^{3/2}). Raise ``ValueError`` if A is not symmetric.
    """
    A = np.asarray(A, dtype=float)
    if A.shape != (3, 3) or not np.allclose(A, A.T, atol=1e-12 * max(1.0, np.abs(A).max())):
        raise ValueError("expected a symmetric 3x3 matrix")
    q = np.trace(A) / 3.0
    B = A - q * np.eye(3)
    p = np.trace(B @ B) / 6.0
    if p < 1e-300:                                   # A is already a multiple of the identity
        return np.array([q, q, q])
    r = np.linalg.det(B) / (2.0 * p ** 1.5)
    r = np.clip(r, -1.0, 1.0)
    phi = np.arccos(r) / 3.0
    lam = q + 2.0 * np.sqrt(p) * np.cos(phi + 2.0 * np.pi * np.arange(3) / 3.0)
    return np.sort(lam)


# ---------------------------------------------------------------- iterative algorithms

def power_iteration(A: np.ndarray, x0: np.ndarray | None = None, tol: float = 1e-12,
                    max_iter: int = 10000) -> tuple[float, np.ndarray, int]:
    """Return (eigenvalue, unit eigenvector, iterations) for the eigenvalue of largest magnitude.

    Repeat x <- A x / |A x| until the change in the unit vector is below ``tol``; the eigenvalue is the
    Rayleigh quotient x^T A x. Starting vector defaults to ones. Raise ``RuntimeError`` if ``max_iter``
    is reached. Converges linearly with ratio |lambda_2 / lambda_1|.
    """
    A = np.asarray(A, dtype=float)
    n = A.shape[0]
    x = np.ones(n) if x0 is None else np.asarray(x0, dtype=float).copy()
    x /= np.linalg.norm(x)
    for k in range(1, max_iter + 1):
        y = A @ x
        ny = np.linalg.norm(y)
        if ny == 0.0:
            return 0.0, x, k
        y /= ny
        if y @ x < 0.0:                               # keep the sign steady for a negative eigenvalue
            y = -y
        delta = np.linalg.norm(y - x)
        x = y
        if delta < tol:
            return float(x @ A @ x), x, k
    raise RuntimeError(f"power iteration did not converge in {max_iter} iterations")


def shifted_inverse_iteration(A: np.ndarray, shift: float, x0: np.ndarray | None = None,
                              tol: float = 1e-12, max_iter: int = 10000) -> tuple[float, np.ndarray, int]:
    """Return (eigenvalue, unit eigenvector, iterations) for the eigenvalue of ``A`` closest to ``shift``.

    Power iteration applied to (A - shift I)^{-1}, done with ``np.linalg.solve`` each step rather than an
    explicit inverse. Raise ``RuntimeError`` if ``max_iter`` is reached or the shift is exactly an
    eigenvalue (singular system).
    """
    A = np.asarray(A, dtype=float)
    n = A.shape[0]
    M = A - shift * np.eye(n)
    x = np.ones(n) if x0 is None else np.asarray(x0, dtype=float).copy()
    x /= np.linalg.norm(x)
    for k in range(1, max_iter + 1):
        try:
            y = np.linalg.solve(M, x)
        except np.linalg.LinAlgError as exc:
            raise RuntimeError("shift is an eigenvalue: the shifted matrix is singular") from exc
        y /= np.linalg.norm(y)
        if y @ x < 0.0:
            y = -y
        delta = np.linalg.norm(y - x)
        x = y
        if delta < tol:
            return float(x @ A @ x), x, k
    raise RuntimeError(f"inverse iteration did not converge in {max_iter} iterations")


def off_diagonal_norm(A: np.ndarray) -> float:
    """Return sqrt(sum of squares of the off-diagonal entries): Jacobi's measure of 'how far from diagonal'."""
    A = np.asarray(A, dtype=float)
    off = A - np.diag(np.diag(A))
    return float(np.sqrt(np.sum(off ** 2)))


def jacobi_eigen(A: np.ndarray, tol: float = 1e-12, max_sweeps: int = 50,
                 history: list[float] | None = None) -> tuple[np.ndarray, np.ndarray, int]:
    """Return (eigenvalues ascending, V with eigenvectors as columns, sweeps) of a symmetric matrix by cyclic Jacobi.

    Each rotation zeroes one off-diagonal pair (p, q) exactly: with tau = (a_qq - a_pp) / (2 a_pq),
    t = sign(tau) / (|tau| + sqrt(1 + tau**2)), c = 1 / sqrt(1 + t**2), s = t c, the update is
    A <- J^T A J and V <- V J. A sweep visits every pair once; stop when ``off_diagonal_norm`` is below
    ``tol`` times the Frobenius norm (rounding limits ``tol`` to about n * 1e-16). If ``history`` is a list, the off-diagonal norm after each sweep is
    appended to it. Raise ``ValueError`` for a non-symmetric input, ``RuntimeError`` if ``max_sweeps`` is hit.
    """
    A = np.array(A, dtype=float)
    n = A.shape[0]
    if A.shape != (n, n) or not np.allclose(A, A.T, atol=1e-12 * max(1.0, np.abs(A).max())):
        raise ValueError("jacobi_eigen needs a symmetric matrix")
    V = np.eye(n)
    scale = np.linalg.norm(A)
    if scale == 0.0:
        return np.zeros(n), V, 0
    for sweep in range(1, max_sweeps + 1):
        for p in range(n - 1):
            for q in range(p + 1, n):
                if abs(A[p, q]) < 1e-300 * scale:    # already zero to within underflow: nothing to rotate
                    A[p, q] = A[q, p] = 0.0
                    continue
                tau = (A[q, q] - A[p, p]) / (2.0 * A[p, q])
                if abs(tau) > 1e150:                  # rotation angle below rounding: t = 1 / (2 tau)
                    t = 0.5 / tau
                else:
                    t = np.sign(tau) / (abs(tau) + np.sqrt(1.0 + tau * tau)) if tau != 0.0 else 1.0
                c = 1.0 / np.sqrt(1.0 + t * t)
                s = t * c
                J = np.eye(n)
                J[p, p] = J[q, q] = c
                J[p, q] = s
                J[q, p] = -s
                A = J.T @ A @ J
                A[p, q] = A[q, p] = 0.0               # exact by construction; remove rounding dust
                V = V @ J
        off = off_diagonal_norm(A)
        if history is not None:
            history.append(off)
        if off < tol * scale:
            order = np.argsort(np.diag(A))
            return np.diag(A)[order], V[:, order], sweep
    raise RuntimeError(f"Jacobi did not converge in {max_sweeps} sweeps")


# ---------------------------------------------------------------- principal axes and their meaning

def rayleigh_quotient(A: np.ndarray, x: np.ndarray) -> float:
    """Return x^T A x / x^T x. For an inertia tensor and a unit axis this is the moment of inertia about that axis."""
    A = np.asarray(A, dtype=float)
    x = np.asarray(x, dtype=float)
    return float(x @ A @ x / (x @ x))


def principal_axes(I: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (principal moments ascending, C_pb) for a symmetric inertia tensor ``I`` in body axes.

    C_pb is a proper rotation (det +1): row k is principal axis k in body coordinates, and
    rotate_tensor(I, C_pb) is diagonal with the principal moments on the diagonal. Built with ``jacobi_eigen``;
    if the eigenvector triad comes out left-handed, the sign of the third axis is flipped.
    """
    moments, V, _ = jacobi_eigen(I)
    if np.linalg.det(V) < 0.0:
        V[:, 2] *= -1.0
    return moments, V.T


def triangle_inequalities_hold(moments: np.ndarray, rtol: float = 1e-12) -> bool:
    """Return True if every principal moment is at most the sum of the other two (and none is negative).

    Any distribution of positive mass satisfies I_1 + I_2 >= I_3 etc., with equality only for a flat
    (planar) body; a tensor that violates this cannot belong to a real body.
    """
    a, b, c = np.sort(np.asarray(moments, dtype=float))
    if a < -rtol * c:
        return False
    return bool(a + b >= c * (1.0 - rtol))


def spin_axis_stability(moments: np.ndarray, axis_index: int) -> dict[str, bool]:
    """Classify a pure spin about principal axis ``axis_index`` (0, 1 or 2 in the ordering of ``moments``).

    Returns {'rigid': ..., 'with_dissipation': ...}. For a torque-free rigid body, spin about the axis of
    largest or smallest moment is stable and spin about the intermediate axis is unstable. If any energy
    is dissipated (fuel slosh, flexing), only spin about the axis of largest moment is stable.
    Raise ``ValueError`` for a bad index or for a tensor with two equal moments touching that axis.
    """
    m = np.asarray(moments, dtype=float)
    if axis_index not in (0, 1, 2):
        raise ValueError("axis_index must be 0, 1 or 2")
    others = np.delete(m, axis_index)
    me = m[axis_index]
    if np.any(np.isclose(others, me, rtol=1e-9)):
        raise ValueError("two principal moments are equal: the axis is not isolated and the rule does not apply")
    is_max = bool(me > others.max())
    is_min = bool(me < others.min())
    return {"rigid": is_max or is_min, "with_dissipation": is_max}


# ---------------------------------------------------------------- the spectral theorem at work

def matrix_function_symmetric(A: np.ndarray, f: Callable[[np.ndarray], np.ndarray]) -> np.ndarray:
    """Return V f(Lambda) V^T for symmetric ``A`` = V Lambda V^T, with ``f`` applied entrywise to the eigenvalues.

    With f = np.exp this is the matrix exponential; with f = np.sqrt the symmetric square root.
    """
    lam, V, _ = jacobi_eigen(A)
    return V @ np.diag(f(lam)) @ V.T


def linear_system_decay_rates(A: np.ndarray) -> np.ndarray:
    """Return the real parts of the eigenvalues of a general (possibly non-symmetric) state matrix.

    For du/dt = A u every solution is a sum of exp(lambda t) terms, so the system is asymptotically stable
    exactly when every real part is negative. Non-symmetric eigenvalues are left to ``numpy.linalg.eigvals``:
    the lesson explains why.
    """
    return np.sort(np.real(np.linalg.eigvals(np.asarray(A, dtype=float))))


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
