"""P03.09 Iterative & Sparse Solvers: reference solution.

Steady heat conduction in a thin square plate (a 2-D Poisson problem) discretised by the five-point
stencil into a sparse SPD system, solved by iterative methods written from scratch:

* Jacobi and Gauss-Seidel (stationary splittings), steepest descent, and the conjugate-gradient
  method (``conjugate_gradient``), with and without preconditioning (Jacobi, incomplete Cholesky);
* the matrix is stored in ``scipy.sparse`` CSR format (storage only: every solver here is ours), and
  ``scipy.sparse.linalg.cg`` / ``spsolve`` are the oracles.

Unknowns are the n x n interior temperatures, ordered row by row (index i + n j); the boundary is at
the wall temperature and enters the right-hand side. The variable-conductivity version uses harmonic
means of k at the cell faces, which keeps the matrix symmetric.

Run:  python solution.py          (prints the verification tables, writes iterative_solvers.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-03/lesson-09-iterative-sparse-solvers
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np
import scipy.sparse as sp

Precond = Callable[[np.ndarray], np.ndarray]


# ---------------------------------------------------------------- the discrete Poisson problem

def grid(n: int, L: float) -> tuple[np.ndarray, np.ndarray, float]:
    """Interior grid of n x n points on (0, L)^2 with spacing h = L/(n+1). Returns (X, Y, h), X[j, i] = x_i."""
    h = L / (n + 1)
    x = h * np.arange(1, n + 1)
    X, Y = np.meshgrid(x, x)
    return X, Y, h


def poisson_matrix(n: int, h: float, k: np.ndarray | float = 1.0) -> sp.csr_matrix:
    """Sparse matrix A of -div(k grad T) on n x n interior points (five-point stencil, Dirichlet edges).

    With constant k the stencil is (k/h^2) [-1, -1, 4, -1, -1]. With a conductivity field k given
    on the (n+2) x (n+2) grid including the boundary ring, the face conductivities are harmonic
    means of the two neighbouring cells, so A stays symmetric positive definite. Row-major ordering:
    unknown (i, j) (x index i, y index j) has index i + n j.
    """
    N = n * n
    if np.isscalar(k):
        kf = np.full((n + 2, n + 2), float(k))
    else:
        kf = np.asarray(k, dtype=float)
        if kf.shape != (n + 2, n + 2):
            raise ValueError("k must be given on the (n+2) x (n+2) grid including the boundary")
    rows, cols, vals = [], [], []
    for j in range(n):
        for i in range(n):
            p = i + n * j
            I, J = i + 1, j + 1                       # indices in the (n+2) grid
            diag = 0.0
            for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ka, kb = kf[J, I], kf[J + dj, I + di]
                kface = 2.0 * ka * kb / (ka + kb)
                diag += kface
                ii, jj = i + di, j + dj
                if 0 <= ii < n and 0 <= jj < n:
                    rows.append(p); cols.append(ii + n * jj); vals.append(-kface / h ** 2)
            rows.append(p); cols.append(p); vals.append(diag / h ** 2)
    return sp.csr_matrix((vals, (rows, cols)), shape=(N, N))


def poisson_rhs(n: int, h: float, source: np.ndarray, T_wall: float, k: np.ndarray | float = 1.0) -> np.ndarray:
    """Right-hand side: the source q''' (W/m^3, on the interior grid, shape (n, n)) plus the wall
    temperature carried over from the boundary neighbours (k_face T_wall / h^2)."""
    b = np.asarray(source, dtype=float).reshape(n * n).copy()
    kf = np.full((n + 2, n + 2), float(k)) if np.isscalar(k) else np.asarray(k, dtype=float)
    for j in range(n):
        for i in range(n):
            I, J = i + 1, j + 1
            for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ii, jj = i + di, j + dj
                if not (0 <= ii < n and 0 <= jj < n):
                    ka, kb = kf[J, I], kf[J + dj, I + di]
                    b[i + n * j] += (2.0 * ka * kb / (ka + kb)) * T_wall / h ** 2
    return b


def manufactured_solution(X: np.ndarray, Y: np.ndarray, L: float, k: float) -> tuple[np.ndarray, np.ndarray]:
    """T = (1 + x/L) sin(pi x/L) sin(pi y/L) (zero on the edges, not a discrete eigenvector) and the
    source f = -k lap T = k [2 (pi/L)^2 (1 + x/L) sin(pi x/L) - (2 pi/L^2) cos(pi x/L)] sin(pi y/L)."""
    s, c, S = np.sin(np.pi * X / L), np.cos(np.pi * X / L), np.sin(np.pi * Y / L)
    g = 1.0 + X / L
    T = g * s * S
    f = k * (2.0 * (np.pi / L) ** 2 * g * s - 2.0 * np.pi / L ** 2 * c) * S
    return T, f


def poisson_eigen_bounds(n: int, h: float, k: float = 1.0) -> tuple[float, float]:
    """Smallest and largest eigenvalues of the constant-k five-point matrix:
    lambda_{p,q} = (4k/h^2)[sin^2(p pi/(2(n+1))) + sin^2(q pi/(2(n+1)))], p, q = 1..n."""
    s_min = np.sin(np.pi / (2 * (n + 1))) ** 2
    s_max = np.sin(n * np.pi / (2 * (n + 1))) ** 2
    return 8.0 * k / h ** 2 * s_min, 8.0 * k / h ** 2 * s_max


def jacobi_spectral_radius(n: int) -> float:
    """rho(I - D^{-1} A) = cos(pi/(n+1)) for the five-point Laplacian: the per-iteration error factor of Jacobi."""
    return float(np.cos(np.pi / (n + 1)))


# ---------------------------------------------------------------- stationary iterations

def jacobi(A: sp.csr_matrix, b: np.ndarray, x0: np.ndarray | None = None, tol: float = 1e-8,
           maxiter: int = 20000) -> tuple[np.ndarray, np.ndarray]:
    """x <- x + D^{-1}(b - A x). Stops when ||r|| <= tol ||b||. Returns (x, residual norms per iteration)."""
    d = A.diagonal()
    x = np.zeros(A.shape[0]) if x0 is None else np.array(x0, dtype=float)
    bnorm = np.linalg.norm(b)
    hist = []
    for _ in range(maxiter):
        r = b - A @ x
        hist.append(np.linalg.norm(r))
        if hist[-1] <= tol * bnorm:
            break
        x = x + r / d
    return x, np.array(hist)


def gauss_seidel(A: sp.csr_matrix, b: np.ndarray, x0: np.ndarray | None = None, tol: float = 1e-8,
                 maxiter: int = 20000) -> tuple[np.ndarray, np.ndarray]:
    """Row-by-row sweep using the newest values: x_i <- (b_i - sum_{j != i} a_ij x_j) / a_ii.
    Same stopping rule as ``jacobi``. Written with the CSR arrays directly (slow in Python, but honest)."""
    A = sp.csr_matrix(A)
    indptr, indices, data = A.indptr, A.indices, A.data
    n = A.shape[0]
    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float)
    bnorm = np.linalg.norm(b)
    hist = []
    for _ in range(maxiter):
        r = b - A @ x
        hist.append(np.linalg.norm(r))
        if hist[-1] <= tol * bnorm:
            break
        for i in range(n):
            s, aii = b[i], 0.0
            for idx in range(indptr[i], indptr[i + 1]):
                j = indices[idx]
                if j == i:
                    aii = data[idx]
                else:
                    s -= data[idx] * x[j]
            x[i] = s / aii
    return x, np.array(hist)


# ---------------------------------------------------------------- Krylov methods

def steepest_descent(A: sp.csr_matrix, b: np.ndarray, x0: np.ndarray | None = None, tol: float = 1e-8,
                     maxiter: int = 20000) -> tuple[np.ndarray, np.ndarray]:
    """Minimise phi(x) = x^T A x / 2 - b^T x along the residual: alpha = r.r / r.Ar. One matvec per step."""
    x = np.zeros(A.shape[0]) if x0 is None else np.array(x0, dtype=float)
    r = b - A @ x
    bnorm = np.linalg.norm(b)
    hist = [np.linalg.norm(r)]
    for _ in range(maxiter):
        if hist[-1] <= tol * bnorm:
            break
        Ar = A @ r
        alpha = (r @ r) / (r @ Ar)
        x = x + alpha * r
        r = r - alpha * Ar
        hist.append(np.linalg.norm(r))
    return x, np.array(hist)


def conjugate_gradient(A: sp.csr_matrix, b: np.ndarray, x0: np.ndarray | None = None, tol: float = 1e-8,
                       maxiter: int | None = None, M: Precond | None = None) -> tuple[np.ndarray, np.ndarray]:
    """(Preconditioned) conjugate gradients for SPD A. ``M(r)`` applies the preconditioner (returns
    an approximation to A^{-1} r); None means no preconditioning.

        r0 = b - A x0, z0 = M r0, p0 = z0
        alpha = (r.z) / (p.Ap);  x += alpha p;  r_new = r - alpha A p
        z_new = M r_new;  beta = (r_new.z_new) / (r.z);  p = z_new + beta p

    One matvec and one preconditioner application per iteration. Stops when ||r|| <= tol ||b||.
    Returns (x, residual norms per iteration, the first being ||r0||).
    """
    n = A.shape[0]
    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float)
    maxiter = 10 * n if maxiter is None else maxiter
    r = b - A @ x
    z = r if M is None else M(r)
    p = z.copy()
    rz = r @ z
    bnorm = np.linalg.norm(b)
    hist = [np.linalg.norm(r)]
    for _ in range(maxiter):
        if hist[-1] <= tol * bnorm:
            break
        Ap = A @ p
        alpha = rz / (p @ Ap)
        x = x + alpha * p
        r = r - alpha * Ap
        hist.append(np.linalg.norm(r))
        z = r if M is None else M(r)
        rz_new = r @ z
        beta = rz_new / rz
        p = z + beta * p
        rz = rz_new
    return x, np.array(hist)


def cg_error_bound(kappa: float, k: np.ndarray | int) -> np.ndarray:
    """||e_k||_A / ||e_0||_A <= 2 ((sqrt(kappa) - 1)/(sqrt(kappa) + 1))^k: the classical CG bound."""
    s = np.sqrt(kappa)
    return 2.0 * ((s - 1.0) / (s + 1.0)) ** np.asarray(k, dtype=float)


def a_norm(A: sp.csr_matrix, e: np.ndarray) -> float:
    """The energy norm sqrt(e^T A e), in which CG is optimal at every step."""
    return float(np.sqrt(e @ (A @ e)))


# ---------------------------------------------------------------- preconditioners

def jacobi_preconditioner(A: sp.csr_matrix) -> Precond:
    """M r = D^{-1} r: scales each row by its diagonal (useful when the diagonal varies, as with variable k)."""
    d = A.diagonal()
    return lambda r: r / d


def incomplete_cholesky(A: sp.csr_matrix) -> sp.csr_matrix:
    """IC(0): the Cholesky recursion restricted to the sparsity pattern of A (no fill-in is kept).

    L_jj = sqrt(a_jj - sum_k L_jk^2), L_ij = (a_ij - sum_k L_ik L_jk)/L_jj for (i, j) in the pattern
    of A only. L L^T approximates A with the same nonzero cost as A itself. Returns L as CSR.
    Raises ValueError if a pivot is not positive (can happen for IC(0) even on SPD matrices).
    """
    A = sp.csc_matrix(A)
    n = A.shape[0]
    Ld = A.toarray()                                    # dense work array; fine for the sizes here
    pattern = Ld != 0.0
    L = np.zeros_like(Ld)
    for j in range(n):
        d = Ld[j, j] - L[j, :j] @ L[j, :j]
        if d <= 0.0:
            raise ValueError(f"incomplete Cholesky broke down at column {j}")
        L[j, j] = np.sqrt(d)
        rows = np.nonzero(pattern[j + 1:, j])[0] + j + 1
        for i in rows:
            L[i, j] = (Ld[i, j] - L[i, :j] @ L[j, :j]) / L[j, j]
    return sp.csr_matrix(np.tril(L))


def sparse_triangular_solve(L: sp.csr_matrix, b: np.ndarray, lower: bool = True) -> np.ndarray:
    """Forward (lower=True) or backward (lower=False) substitution on a sparse triangular CSR matrix."""
    L = sp.csr_matrix(L)
    indptr, indices, data = L.indptr, L.indices, L.data
    n = L.shape[0]
    x = np.zeros(n)
    order = range(n) if lower else range(n - 1, -1, -1)
    for i in order:
        s, dii = b[i], 0.0
        for idx in range(indptr[i], indptr[i + 1]):
            j = indices[idx]
            if j == i:
                dii = data[idx]
            else:
                s -= data[idx] * x[j]
        x[i] = s / dii
    return x


def ic_preconditioner(A: sp.csr_matrix) -> Precond:
    """M r = (L L^T)^{-1} r with L from ``incomplete_cholesky``: two sparse triangular solves per application."""
    L = incomplete_cholesky(A)
    Lt = sp.csr_matrix(L.T)
    return lambda r: sparse_triangular_solve(Lt, sparse_triangular_solve(L, r, lower=True), lower=False)


# ---------------------------------------------------------------- the plate and the studies

def plate_conductivity(n: int, L: float, k_plate: float, k_insert: float, insert: tuple[float, float, float, float]
                       ) -> np.ndarray:
    """Conductivity on the (n+2) x (n+2) grid: k_plate everywhere, k_insert inside the rectangle
    insert = (x0, x1, y0, y1) in metres."""
    h = L / (n + 1)
    x = h * np.arange(n + 2)
    X, Y = np.meshgrid(x, x)
    k = np.full((n + 2, n + 2), float(k_plate))
    x0, x1, y0, y1 = insert
    k[(X >= x0) & (X <= x1) & (Y >= y0) & (Y <= y1)] = k_insert
    return k


def solve_plate(n: int, L: float, kt: float | np.ndarray, q_flux: float, T_wall: float, tol: float = 1e-12
                ) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Steady temperature of a thin plate with uniform absorbed flux q_flux [W/m^2], edges at T_wall,
    in-plane conductance k t [W/K] (thickness times conductivity): -div(k t grad T) = q_flux.
    Returns (X, Y, T (n, n), CG iterations)."""
    X, Y, h = grid(n, L)
    A = poisson_matrix(n, h, kt)
    b = poisson_rhs(n, h, np.full((n, n), q_flux), T_wall, kt)
    T, hist = conjugate_gradient(A, b, tol=tol)
    return X, Y, T.reshape(n, n), hist.size - 1


def discretisation_errors(ns: list[int], L: float = 0.2, k: float = 1.0) -> np.ndarray:
    """Max-norm error against the manufactured solution for each n (CG solved to 1e-12): order study."""
    errs = []
    for n in ns:
        X, Y, h = grid(n, L)
        T_ex, f = manufactured_solution(X, Y, L, k)
        A = poisson_matrix(n, h, k)
        b = poisson_rhs(n, h, f, 0.0, k)
        T, _ = conjugate_gradient(A, b, tol=1e-12)
        errs.append(np.max(np.abs(T.reshape(n, n) - T_ex)))
    return np.array(errs)


def iterations_vs_n(method: Callable, ns: list[int], tol: float = 1e-6, L: float = 0.2) -> np.ndarray:
    """Iterations to reduce the residual by ``tol`` on the manufactured problem, for each n."""
    its = []
    for n in ns:
        X, Y, h = grid(n, L)
        _, f = manufactured_solution(X, Y, L, 1.0)
        A = poisson_matrix(n, h, 1.0)
        b = poisson_rhs(n, h, f, 0.0, 1.0)
        _, hist = method(A, b, tol=tol)
        its.append(hist.size - 1)
    return np.array(its)


def scaling_exponent(ns: np.ndarray, its: np.ndarray) -> float:
    """Least-squares slope of log(iterations) against log(n)."""
    return float(np.polyfit(np.log(np.asarray(ns, dtype=float)), np.log(np.asarray(its, dtype=float)), 1)[0])


def convergence_order(hs: np.ndarray, errs: np.ndarray) -> float:
    """Least-squares slope of log(err) against log(h)."""
    hs, errs = np.asarray(hs, dtype=float), np.asarray(errs, dtype=float)
    if np.any(errs <= 0):
        raise ValueError("errors must be positive to take logarithms")
    return float(np.polyfit(np.log(hs), np.log(errs), 1)[0])


def main(out_dir: str | Path = ".") -> None:
    """Run the lesson's experiments, print the tables and save ``iterative_solvers.png``."""
    import time
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.sparse.linalg import cg as scipy_cg, spsolve, eigsh

    out_dir = Path(out_dir)
    L_plate, kt, q_flux, T_wall = 0.2, 170.0 * 2e-3, 2000.0, 300.0

    # ---- the plate
    for n in (16, 32, 64):
        X, Y, T, its = solve_plate(n, L_plate, kt, q_flux, T_wall)
        print(f"plate n = {n:3d} ({n * n:5d} unknowns): centre rise {T[n // 2, n // 2] - T_wall:.3f} K "
              f"(grid point nearest the centre), max rise {T.max() - T_wall:.3f} K, CG iterations {its}")
    print(f"series solution for a square plate: centre rise 0.07367 q L^2 / (k t) = {0.07367 * q_flux * L_plate ** 2 / kt:.3f} K")

    # ---- matrix facts and the solver race at n = 24
    n = 24
    X, Y, h = grid(n, L_plate)
    A = poisson_matrix(n, h, 1.0)
    lam_min, lam_max = poisson_eigen_bounds(n, h)
    ev_min = eigsh(A, k=1, which="SA", return_eigenvectors=False)[0]
    ev_max = eigsh(A, k=1, which="LA", return_eigenvectors=False)[0]
    print(f"\nn = {n}: nnz = {A.nnz} of {n ** 4} entries ({A.nnz / n ** 4 * 100:.2f} %); eigenvalues {lam_min:.4e} .. {lam_max:.4e} "
          f"(formula) vs {ev_min:.4e} .. {ev_max:.4e} (eigsh); kappa = {lam_max / lam_min:.1f} ~ (2(n+1)/pi)^2 = {(2 * (n + 1) / np.pi) ** 2:.1f}")
    T_ex, f = manufactured_solution(X, Y, L_plate, 1.0)
    b = poisson_rhs(n, h, f, 0.0, 1.0)
    x_direct = spsolve(sp.csc_matrix(A), b)
    t0 = time.perf_counter(); x_j, h_j = jacobi(A, b, tol=1e-6); t_j = time.perf_counter() - t0
    t0 = time.perf_counter(); x_gs, h_gs = gauss_seidel(A, b, tol=1e-6, maxiter=3000); t_gs = time.perf_counter() - t0
    t0 = time.perf_counter(); x_sd, h_sd = steepest_descent(A, b, tol=1e-6); t_sd = time.perf_counter() - t0
    t0 = time.perf_counter(); x_cg, h_cg = conjugate_gradient(A, b, tol=1e-6); t_cg = time.perf_counter() - t0
    t0 = time.perf_counter(); x_sp, _ = scipy_cg(A, b, rtol=1e-6, atol=0.0); t_sp = time.perf_counter() - t0
    rho_j = jacobi_spectral_radius(n)
    print(f"residual reduction 1e-6 at n = {n}: Jacobi {h_j.size - 1} iterations ({t_j * 1e3:.0f} ms; theory ~ ln(1e-6)/ln(rho_J) = "
          f"{np.log(1e-6) / np.log(rho_j):.0f} with rho_J = {rho_j:.5f}), Gauss-Seidel {h_gs.size - 1} ({t_gs * 1e3:.0f} ms; rho_GS = rho_J^2 -> "
          f"{np.log(1e-6) / np.log(rho_j ** 2):.0f}), steepest descent {h_sd.size - 1} ({t_sd * 1e3:.0f} ms), CG {h_cg.size - 1} ({t_cg * 1e3:.1f} ms), "
          f"scipy cg {t_sp * 1e3:.1f} ms")
    print(f"max |x - spsolve|: Jacobi {np.max(np.abs(x_j - x_direct)):.1e}, GS {np.max(np.abs(x_gs - x_direct)):.1e}, "
          f"SD {np.max(np.abs(x_sd - x_direct)):.1e}, CG {np.max(np.abs(x_cg - x_direct)):.1e}, scipy {np.max(np.abs(x_sp - x_direct)):.1e}")
    # observed Jacobi contraction
    ratios = h_j[-50:] / h_j[-51:-1]
    print(f"observed Jacobi residual ratio over the last 50 iterations {ratios.mean():.5f} vs rho_J = {rho_j:.5f}")
    # CG bound in the A-norm
    kappa = lam_max / lam_min
    e0 = a_norm(A, x_direct)
    xk, hk = conjugate_gradient(A, b, tol=1e-14, maxiter=40)
    print(f"CG A-norm error after 40 iterations: {a_norm(A, xk - x_direct) / e0:.2e}; bound 2((sqrt(k)-1)/(sqrt(k)+1))^40 = {cg_error_bound(kappa, 40):.2e}")

    # ---- scaling with n and discretisation order
    ns = [8, 16, 32, 64]
    its_cg = iterations_vs_n(conjugate_gradient, ns)
    its_j = iterations_vs_n(jacobi, ns[:3])
    its_sd = iterations_vs_n(steepest_descent, ns[:3])
    print(f"\niterations for 1e-6 at n = {ns}: CG {its_cg.tolist()} (exponent {scaling_exponent(ns, its_cg):.2f}, theory 1: ~sqrt(kappa) ~ n); "
          f"Jacobi {its_j.tolist()} (exponent {scaling_exponent(ns[:3], its_j):.2f}, theory 2); steepest descent {its_sd.tolist()} "
          f"(exponent {scaling_exponent(ns[:3], its_sd):.2f}, theory 2)")
    errs = discretisation_errors(ns)
    hs = L_plate / (np.array(ns) + 1)
    print(f"discretisation error vs manufactured solution: {np.array2string(errs, precision=3, formatter={'float_kind': lambda v: f'{v:.3e}'})}; observed order {convergence_order(hs, errs):.2f}")

    # ---- preconditioning on the plate with an insert
    n = 32
    X, Y, h = grid(n, L_plate)
    for ratio in (1.0, 10.0, 100.0):
        kfield = plate_conductivity(n, L_plate, kt, kt / ratio, (0.06, 0.14, 0.06, 0.14))
        A = poisson_matrix(n, h, kfield)
        b = poisson_rhs(n, h, np.full((n, n), q_flux), T_wall, kfield)
        x_ref = spsolve(sp.csc_matrix(A), b)
        _, h0 = conjugate_gradient(A, b, tol=1e-8)
        _, h1 = conjugate_gradient(A, b, tol=1e-8, M=jacobi_preconditioner(A))
        t0 = time.perf_counter(); x2, h2 = conjugate_gradient(A, b, tol=1e-8, M=ic_preconditioner(A)); t2 = time.perf_counter() - t0
        lmin = eigsh(A, k=1, which="SA", return_eigenvectors=False)[0]; lmax = eigsh(A, k=1, which="LA", return_eigenvectors=False)[0]
        print(f"insert k ratio {ratio:5.0f}: kappa = {lmax / lmin:.2e}; CG {h0.size - 1} iterations, Jacobi-PCG {h1.size - 1}, IC(0)-PCG {h2.size - 1} "
              f"({t2 * 1e3:.0f} ms incl. factorisation); max rise {x_ref.max() - T_wall:.2f} K; PCG error {np.max(np.abs(x2 - x_ref)):.1e} K")

    # ---- figure
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    ax = axes[0, 0]
    X, Y, T, _ = solve_plate(64, L_plate, kt, q_flux, T_wall)
    cs = ax.contourf(X, Y, T - T_wall, 20)
    fig.colorbar(cs, ax=ax, label="temperature rise [K]")
    ax.set_aspect("equal"); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.set_title("plate, uniform 2 kW/m^2, edges at 300 K")
    ax = axes[0, 1]
    for name, hh in (("Jacobi", h_j), ("Gauss-Seidel", h_gs), ("steepest descent", h_sd), ("CG", h_cg)):
        ax.semilogy(hh / hh[0], label=name)
    ax.set_xlim(0, 1000); ax.set_xlabel("iteration"); ax.set_ylabel("||r|| / ||r0||"); ax.set_title("n = 24: residual histories"); ax.legend(fontsize=8)
    ax = axes[1, 0]
    ax.loglog(ns, its_cg, "o-", label=f"CG (slope {scaling_exponent(ns, its_cg):.2f})")
    ax.loglog(ns[:3], its_j, "s-", label=f"Jacobi (slope {scaling_exponent(ns[:3], its_j):.2f})")
    ax.loglog(ns[:3], its_sd, "^-", label=f"steepest descent (slope {scaling_exponent(ns[:3], its_sd):.2f})")
    ax.set_xlabel("n (grid points per side)"); ax.set_ylabel("iterations for 1e-6"); ax.set_title("work against resolution"); ax.legend(fontsize=8)
    ax = axes[1, 1]
    ax.loglog(hs, errs, "o-", label=f"order {convergence_order(hs, errs):.2f}")
    ax.loglog(hs, errs[0] * (hs / hs[0]) ** 2, "k--", lw=0.8, label="h^2")
    ax.set_xlabel("h [m]"); ax.set_ylabel("max |T - exact|"); ax.set_title("manufactured solution: discretisation error"); ax.legend(fontsize=8)
    fig.suptitle("P03.09  Iterative solvers on a 2-D Poisson problem")
    fig.tight_layout()
    fig.savefig(out_dir / "iterative_solvers.png", dpi=150)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'iterative_solvers.png'}")


if __name__ == "__main__":
    main()
