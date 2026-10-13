"""P01.18 SVD & Matrix Norms: reference solution.

Conditioning of a structural stiffness matrix: vector and matrix norms, the singular value decomposition
by one-sided Jacobi rotations, the condition number, the pseudoinverse and low-rank approximation, all
written from scratch in NumPy and checked against ``numpy.linalg.svd`` / ``norm`` / ``cond`` / ``pinv``.

Structural model: an axially loaded bar (a spar cap, a truss member, a landing-gear strut) discretised into
n elements in series, node 0 clamped, nodes 1..n free. Element i (between nodes i-1 and i) has stiffness
k_i = E A_i / L_i [N/m]. Nodal loads f [N] and displacements u [m] satisfy K u = f with the tridiagonal
    K[i, i] = k_i + k_{i+1} (k_{n+1} = 0),   K[i, i+1] = K[i+1, i] = -k_{i+1}      (0-based node i+1).

Run:  python solution.py          (prints the bar example, writes stiffness_conditioning.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-01/lesson-18-svd-and-norms
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

# Aluminium alloy bar, 10 mm x 10 mm, 1 m long in 10 elements (illustrative sizing).
E_ALUMINIUM = 70.0e9      # Pa
BAR_AREA = 1.0e-4         # m**2
BAR_LENGTH = 1.0          # m
N_ELEMENTS = 10


# ---------------------------------------------------------------- norms (from scratch)

def vector_norm(x: np.ndarray, p: float = 2) -> float:
    """Return the p-norm of a vector for p = 1, 2 or inf (``np.inf``). Raise ``ValueError`` for any other p.

    |x|_1 = sum |x_i|,  |x|_2 = sqrt(sum x_i**2),  |x|_inf = max |x_i|.
    """
    x = np.asarray(x, dtype=float).ravel()
    if p == 1:
        return float(np.sum(np.abs(x)))
    if p == 2:
        return float(np.sqrt(np.sum(x * x)))
    if p == np.inf:
        return float(np.max(np.abs(x))) if x.size else 0.0
    raise ValueError("p must be 1, 2 or inf")


def matrix_norm_1(A: np.ndarray) -> float:
    """Return the induced 1-norm: the largest column sum of absolute values."""
    A = np.asarray(A, dtype=float)
    return float(np.max(np.sum(np.abs(A), axis=0)))


def matrix_norm_inf(A: np.ndarray) -> float:
    """Return the induced infinity-norm: the largest row sum of absolute values."""
    A = np.asarray(A, dtype=float)
    return float(np.max(np.sum(np.abs(A), axis=1)))


def frobenius_norm(A: np.ndarray) -> float:
    """Return sqrt(sum of all squared entries) = sqrt(sum of squared singular values)."""
    A = np.asarray(A, dtype=float)
    return float(np.sqrt(np.sum(A * A)))


def matrix_norm_2_power(A: np.ndarray, tol: float = 1e-13, max_iter: int = 20000) -> float:
    """Return the induced 2-norm, sigma_max(A), as the square root of the dominant eigenvalue of A^T A by power iteration.

    Starting vector: ones. Stops when the eigenvalue estimate changes by less than ``tol`` relative. Raise
    ``RuntimeError`` at ``max_iter``. The 2-norm is the only one here that needs an iteration: it is the
    largest stretch max |A x| / |x| over all directions, and that maximum is not read off the entries.
    """
    A = np.asarray(A, dtype=float)
    M = A.T @ A
    x = np.ones(M.shape[0])
    x /= np.sqrt(x @ x)
    lam_old = 0.0
    for _ in range(max_iter):
        y = M @ x
        lam = x @ y
        ny = np.sqrt(y @ y)
        if ny == 0.0:
            return 0.0
        x = y / ny
        if abs(lam - lam_old) <= tol * max(lam, 1e-300):
            return float(np.sqrt(max(lam, 0.0)))
        lam_old = lam
    raise RuntimeError("power iteration for the 2-norm did not converge")


# ---------------------------------------------------------------- the SVD (from scratch)

def svd_jacobi(A: np.ndarray, tol: float = 1e-14, max_sweeps: int = 60) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (U, s, Vt) with A = U diag(s) Vt, s descending, by one-sided (Hestenes) Jacobi rotations.

    For an m x n matrix with m >= n: rotate pairs of columns of a working copy W = A V until every pair is
    orthogonal (|w_i . w_j| <= tol |w_i| |w_j|). Then s_i = |w_i|, u_i = w_i / s_i, and V accumulates the rotations.
    U is m x n with orthonormal columns (a column for a zero singular value is completed by Gram-Schmidt against
    the others), Vt is n x n orthogonal. Raise ``ValueError`` if m < n, ``RuntimeError`` if ``max_sweeps`` is hit.
    """
    W = np.array(A, dtype=float)
    m, n = W.shape
    if m < n:
        raise ValueError("svd_jacobi expects m >= n (transpose a wide matrix first)")
    V = np.eye(n)
    for sweep in range(1, max_sweeps + 1):
        rotated = False
        for p in range(n - 1):
            for q in range(p + 1, n):
                alpha = W[:, p] @ W[:, p]
                beta = W[:, q] @ W[:, q]
                gamma = W[:, p] @ W[:, q]
                if abs(gamma) <= tol * np.sqrt(alpha * beta) or gamma == 0.0:
                    continue
                rotated = True
                zeta = (beta - alpha) / (2.0 * gamma)
                t = np.sign(zeta) / (abs(zeta) + np.sqrt(1.0 + zeta * zeta)) if zeta != 0.0 else 1.0
                c = 1.0 / np.sqrt(1.0 + t * t)
                s = c * t
                wp, wq = W[:, p].copy(), W[:, q].copy()
                W[:, p], W[:, q] = c * wp - s * wq, s * wp + c * wq
                vp, vq = V[:, p].copy(), V[:, q].copy()
                V[:, p], V[:, q] = c * vp - s * vq, s * vp + c * vq
        if not rotated:
            break
    else:
        raise RuntimeError(f"one-sided Jacobi did not converge in {max_sweeps} sweeps")
    sigma = np.sqrt(np.sum(W * W, axis=0))
    order = np.argsort(-sigma)
    sigma, W, V = sigma[order], W[:, order], V[:, order]
    U = np.zeros((m, n))
    smax = sigma[0] if sigma.size else 0.0
    for j in range(n):
        if sigma[j] > 1e-300 and sigma[j] > 1e-15 * smax:
            U[:, j] = W[:, j] / sigma[j]
        else:
            sigma[j] = 0.0
    # Complete U for zero singular values: Gram-Schmidt a unit vector against the columns found so far.
    for j in range(n):
        if np.all(U[:, j] == 0.0):
            for e in np.eye(m):
                v = e - U @ (U.T @ e)
                v -= U @ (U.T @ v)
                if np.linalg.norm(v) > 1e-8:
                    U[:, j] = v / np.linalg.norm(v)
                    break
    return U, sigma, V.T


def condition_number(A: np.ndarray) -> float:
    """Return sigma_max / sigma_min from ``svd_jacobi``; ``inf`` if sigma_min is zero (singular matrix)."""
    _, s, _ = svd_jacobi(A)
    return float("inf") if s[-1] == 0.0 else float(s[0] / s[-1])


def pseudoinverse(A: np.ndarray, rcond: float = 1e-12) -> np.ndarray:
    """Return A^+ = V diag(1/s) U^T with singular values below ``rcond * s_max`` treated as zero.

    For full-column-rank A this is (A^T A)^{-1} A^T (least squares); for full-row-rank A it is A^T (A A^T)^{-1}
    (minimum-norm solution); for a singular square A it inverts on the range and ignores the null space.
    """
    A = np.asarray(A, dtype=float)
    transpose = A.shape[0] < A.shape[1]
    U, s, Vt = svd_jacobi(A.T if transpose else A)
    s_inv = np.array([1.0 / si if si > rcond * s[0] else 0.0 for si in s])
    P = Vt.T @ np.diag(s_inv) @ U.T
    return P.T if transpose else P


def low_rank_approximation(A: np.ndarray, k: int) -> np.ndarray:
    """Return the best rank-k approximation sum_{i<k} s_i u_i v_i^T (Eckart-Young). Raise ``ValueError`` for bad k."""
    A = np.asarray(A, dtype=float)
    if k < 1 or k > min(A.shape):
        raise ValueError("k must be between 1 and min(m, n)")
    U, s, Vt = svd_jacobi(A)
    return (U[:, :k] * s[:k]) @ Vt[:k, :]


def null_space(A: np.ndarray, rtol: float = 1e-12) -> np.ndarray:
    """Return an orthonormal basis (columns, shape (n, n - r)) of the null space of an m x n matrix of rank r.

    The row space is spanned by the right singular vectors with s > rtol * s_max (from ``svd_jacobi`` of A, or of
    A^T for a wide matrix); the null space is its orthogonal complement, completed by Gram-Schmidt against the
    coordinate vectors so that a wide matrix gets all n - r directions, not only those the reduced SVD lists.
    """
    A = np.asarray(A, dtype=float)
    m, n = A.shape
    if m >= n:
        _, s, Vt = svd_jacobi(A)
        rows = Vt
    else:
        U, s, _ = svd_jacobi(A.T)          # A^T = U s V^T  =>  A = V s U^T: row space of A is spanned by U
        rows = U.T
    smax = s[0] if s.size else 0.0
    basis = rows[s > rtol * smax] if smax > 0.0 else np.zeros((0, n))
    R = basis.T.copy()                      # orthonormal basis of the row space, as columns
    Z = []
    for e in np.eye(n):
        v = e - R @ (R.T @ e)
        for z in Z:
            v -= (z @ v) * z
        v -= R @ (R.T @ v)
        nv = np.linalg.norm(v)
        if nv > 1e-8:
            Z.append(v / nv)
        if len(Z) == n - R.shape[1]:
            break
    return np.array(Z).T if Z else np.zeros((n, 0))


# ---------------------------------------------------------------- the structural model

def element_stiffness(E: float, area: float, length: float) -> float:
    """Return k = E A / L in N/m for a uniform axial element. Raise ``ValueError`` unless all inputs are positive."""
    if E <= 0.0 or area <= 0.0 or length <= 0.0:
        raise ValueError("E, area and length must be positive")
    return E * area / length


def chain_stiffness_matrix(k: np.ndarray) -> np.ndarray:
    """Assemble the n x n stiffness matrix of n elements in series with node 0 clamped, from element stiffnesses ``k``.

    Non-negative stiffnesses are allowed (a zero models a failed element) so that a mechanism can be studied;
    a negative value raises ``ValueError``.
    """
    k = np.asarray(k, dtype=float)
    if k.ndim != 1 or k.size == 0 or np.any(k < 0.0):
        raise ValueError("k must be a 1-D array of non-negative element stiffnesses")
    n = k.size
    K = np.zeros((n, n))
    for i in range(n):
        K[i, i] += k[i]
        if i + 1 < n:
            K[i, i] += k[i + 1]
            K[i, i + 1] = K[i + 1, i] = -k[i + 1]
    return K


def uniform_chain_condition_analytic(n: int) -> float:
    """Return cond_2(K) for n identical elements, from the exact eigenvalues 2k(1 - cos((2j-1) pi / (2n+1))).

    The ratio does not depend on k. For large n it grows as 16 n**2 / pi**2.
    """
    if n < 1:
        raise ValueError("n must be at least 1")
    j = np.arange(1, n + 1)
    lam = 1.0 - np.cos((2 * j - 1) * np.pi / (2 * n + 1))
    return float(lam.max() / lam.min())


def tip_deflection_analytic(force: float, E: float, area: float, length: float) -> float:
    """Return the end displacement F L / (E A) of a clamped-free uniform bar under an axial tip load."""
    return force * length / (E * area)


def error_amplification(K: np.ndarray, f: np.ndarray, df: np.ndarray) -> tuple[float, float]:
    """Return (actual relative error in u, the bound cond(K) |df| / |f|) for K u = f perturbed to K (u + du) = f + df.

    Relative errors in the 2-norm. The bound is attained when f lies along the first left singular vector
    of K and df along the last.
    """
    K = np.asarray(K, dtype=float)
    f = np.asarray(f, dtype=float)
    df = np.asarray(df, dtype=float)
    u = np.linalg.solve(K, f)
    du = np.linalg.solve(K, df)
    actual = vector_norm(du) / vector_norm(u)
    bound = condition_number(K) * vector_norm(df) / vector_norm(f)
    return float(actual), float(bound)


def worst_case_loads(K: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (f, df) unit vectors that attain the error-amplification bound: u_1 and u_n of the SVD of K."""
    U, _, _ = svd_jacobi(K)
    return U[:, 0].copy(), U[:, -1].copy()


def mechanism_modes(K: np.ndarray, rtol: float = 1e-10) -> np.ndarray:
    """Return the null-space basis of a stiffness matrix: displacement patterns that need no force (mechanisms).

    An empty (n x 0) array means the structure is properly constrained.
    """
    return null_space(K, rtol)


def main(out_dir: str | Path = ".") -> None:
    """Analyse the 10-element aluminium bar and plot the condition number against n and against a soft element."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    np.set_printoptions(precision=4, suppress=True)
    k = element_stiffness(E_ALUMINIUM, BAR_AREA, BAR_LENGTH / N_ELEMENTS)
    K = chain_stiffness_matrix(np.full(N_ELEMENTS, k))
    print(f"element stiffness k = E A / L_e = {k:.3e} N/m; K is {N_ELEMENTS}x{N_ELEMENTS}, tridiagonal")
    print(f"norms: |K|_1 = {matrix_norm_1(K):.3e}, |K|_inf = {matrix_norm_inf(K):.3e}, |K|_F = {frobenius_norm(K):.3e},"
          f" |K|_2 (power) = {matrix_norm_2_power(K):.3e}, numpy |K|_2 = {np.linalg.norm(K, 2):.3e}")
    U, s, Vt = svd_jacobi(K)
    print(f"singular values [N/m]: {s}")
    print(f"max |K - U S V^T| = {np.max(np.abs(K - (U * s) @ Vt)):.1e}, max |U^T U - I| = {np.max(np.abs(U.T @ U - np.eye(N_ELEMENTS))):.1e}")
    print(f"numpy.linalg.svd singular values: {np.linalg.svd(K, compute_uv=False)}")
    cond = condition_number(K)
    print(f"cond_2(K) = {cond:.3f};  numpy.linalg.cond = {np.linalg.cond(K):.3f};  analytic = {uniform_chain_condition_analytic(N_ELEMENTS):.3f};"
          f"  16 n^2 / pi^2 = {16 * N_ELEMENTS ** 2 / np.pi ** 2:.1f}")

    f = np.zeros(N_ELEMENTS)
    f[-1] = 1000.0
    u = np.linalg.solve(K, f)
    print(f"\ntip load 1 kN: tip displacement {u[-1] * 1e3:.4f} mm (analytic F L / E A = {tip_deflection_analytic(1000.0, E_ALUMINIUM, BAR_AREA, BAR_LENGTH) * 1e3:.4f} mm)")
    rng = np.random.default_rng(18)
    df = 0.01 * vector_norm(f) * rng.standard_normal(N_ELEMENTS)
    df *= 0.01 * vector_norm(f) / vector_norm(df)
    actual, bound = error_amplification(K, f, df)
    print(f"1 % random load error: relative displacement error {actual * 100:.2f} % (bound {bound * 100:.1f} %)")
    f_w, df_w = worst_case_loads(K)
    actual_w, bound_w = error_amplification(K, f_w, 0.01 * df_w)
    print(f"1 % worst-case load error: relative displacement error {actual_w * 100:.1f} % (bound {bound_w * 100:.1f} %)")
    print(f"worst-case load pattern f (first left singular vector): {f_w}")
    kk = np.full(N_ELEMENTS, k)
    kk[4] *= 1e-3
    K_damaged = chain_stiffness_matrix(kk)
    actual_d, bound_d = error_amplification(K_damaged, f, df)
    f_wd, df_wd = worst_case_loads(K_damaged)
    actual_wd, bound_wd = error_amplification(K_damaged, f_wd, 0.01 * df_wd)
    print(f"element 5 at 0.1 % stiffness, same 1 % random load error: {actual_d * 100:.3f} % (bound {bound_d * 100:.3e} %);"
          f" worst case {actual_wd * 100:.3e} % (bound {bound_wd * 100:.3e} %)")

    print("\nelement 5 damaged: stiffness ratio r = k_5 / k")
    ratios = np.array([1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 0.0])
    conds = []
    for r in ratios:
        kk = np.full(N_ELEMENTS, k)
        kk[4] = r * k
        c = condition_number(chain_stiffness_matrix(kk))
        conds.append(c)
        print(f"   r = {r:7.0e}: cond_2(K) = {c:.3e}")
    modes = mechanism_modes(chain_stiffness_matrix(np.where(np.arange(N_ELEMENTS) == 4, 0.0, k)))
    print(f"with k_5 = 0 the null space has dimension {modes.shape[1]}; mode shape (nodes 1..10): {modes[:, 0] if modes.size else '-'}")

    ns = np.array([2, 4, 8, 16, 32, 64, 128])
    cond_n = [condition_number(chain_stiffness_matrix(np.full(n, k))) for n in ns]
    print("\ncondition number against element count:")
    for n, c in zip(ns, cond_n):
        print(f"   n = {n:3d}: cond_2 = {c:9.1f}   analytic {uniform_chain_condition_analytic(int(n)):9.1f}   16n^2/pi^2 = {16 * n ** 2 / np.pi ** 2:9.1f}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 3.8))
    ax1.loglog(ns, cond_n, "o", label="one-sided Jacobi SVD")
    ax1.loglog(ns, [uniform_chain_condition_analytic(int(n)) for n in ns], "-", label="exact eigenvalues")
    ax1.loglog(ns, 16 * ns ** 2 / np.pi ** 2, "--", label=r"$16 n^2 / \pi^2$")
    ax1.set_xlabel("number of elements n")
    ax1.set_ylabel(r"cond$_2(K)$")
    ax1.grid(True, which="both", lw=0.4)
    ax1.legend(fontsize=8)
    ax2.loglog(ratios[1:-1], conds[1:-1], "s-")
    ax2.set_xlabel("stiffness ratio of damaged element")
    ax2.set_ylabel(r"cond$_2(K)$, n = 10")
    ax2.grid(True, which="both", lw=0.4)
    fig.suptitle("Conditioning of the bar stiffness matrix")
    fig.tight_layout()
    fig.savefig(out_dir / "stiffness_conditioning.png", dpi=150)
    print(f"saved {out_dir / 'stiffness_conditioning.png'}")


if __name__ == "__main__":
    main()
