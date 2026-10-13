"""P01.13 Systems of Linear Equations: starter code.

Gaussian elimination with partial pivoting, LU factorisation, determinant and inverse, all written
from scratch in NumPy and checked against ``numpy.linalg``; applied to the member forces of a
statically determinate planar truss (the method of joints as one linear system).

Conventions:
    * ``A`` is an (n, n) float array, ``b`` an (n,) array; every routine copies its inputs.
    * Pivoting is partial (row swaps on the largest |entry| in the column); ``P A = L U`` with
      ``P`` a permutation matrix, ``L`` unit lower triangular, ``U`` upper triangular.
    * Truss: nodes (j, 2) in metres; members as (i, k) index pairs; supports as (node, dof) with
      dof 0 = x, 1 = y; loads as {node: (Fx, Fy)} in newtons. Unknowns are the m member forces
      (tension positive) followed by the r reactions; the system is A x = -P.

Complete every function marked NotImplementedError.

Run:  python main.py          (prints the truss solution and the verification table, writes truss_forces.png)
Test: pytest curriculum/phase-01/lesson-13-linear-systems
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

SINGULAR_TOL = 1e-14      # a pivot below this fraction of max|A| is treated as zero


# ----------------------------------------------------------------------------- input checks

def _as_system(A: np.ndarray, b: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray | None]:
    """Copy A (and b) as float arrays and check that A is square and b matches it."""
    A = np.array(A, dtype=float, copy=True)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError(f"A must be square, got shape {A.shape}")
    if b is None:
        return A, None
    b = np.array(b, dtype=float, copy=True)
    if b.shape != (A.shape[0],):
        raise ValueError(f"b must have shape ({A.shape[0]},), got {b.shape}")
    return A, b


# ----------------------------------------------------------------------------- triangular solves

def back_substitution(U: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Solve U x = c for upper-triangular U, last row first. Raise ``ValueError`` on a zero diagonal."""
    raise NotImplementedError


def forward_substitution(L: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Solve L y = c for lower-triangular L, first row first. Raise ``ValueError`` on a zero diagonal."""
    raise NotImplementedError


# ----------------------------------------------------------------------------- elimination

def forward_elimination(A: np.ndarray, b: np.ndarray, pivot: bool = True) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reduce A x = b to U x = c by row operations; return (U, c, perm).

    ``perm`` records the row order: row i of U came from row perm[i] of A. With ``pivot=True`` each
    column's largest remaining |entry| is swapped into the pivot position first (partial pivoting).
    Raise ``ValueError`` if a pivot smaller than ``SINGULAR_TOL`` times max|A| is met (singular matrix).
    """
    raise NotImplementedError


def gaussian_elimination(A: np.ndarray, b: np.ndarray, pivot: bool = True) -> np.ndarray:
    """Solve A x = b by forward elimination (partial pivoting by default) and back substitution."""
    raise NotImplementedError


def lu_decompose(A: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Factor P A = L U with partial pivoting; return (P, L, U). Raise ``ValueError`` if A is singular."""
    raise NotImplementedError


def lu_solve(P: np.ndarray, L: np.ndarray, U: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Solve A x = b given P A = L U: forward-substitute L y = P b, then back-substitute U x = y."""
    raise NotImplementedError


def determinant(A: np.ndarray) -> float:
    """det(A) from the LU factors: the product of U's diagonal times the sign of the permutation."""
    raise NotImplementedError


def inverse_gauss_jordan(A: np.ndarray) -> np.ndarray:
    """Invert A by reducing [A | I] to [I | A^-1] with partial pivoting. Raise ``ValueError`` if singular."""
    raise NotImplementedError


def residual_norm(A: np.ndarray, x: np.ndarray, b: np.ndarray) -> float:
    """Relative residual ||b - A x|| / ||b|| (2-norm); for b = 0 return the absolute residual."""
    raise NotImplementedError


def iterative_refinement(A: np.ndarray, b: np.ndarray, steps: int = 3, pivot: bool = True) -> np.ndarray:
    """Solve A x = b, then ``steps`` times solve A d = r for the residual r = b - A x and update x += d.

    Each correction reuses the same (possibly poor) elimination; the residual is computed with the
    full matrix, so a backward-unstable solve is repaired but an ill-conditioned problem is not.
    """
    raise NotImplementedError


def elimination_flops(n: int) -> float:
    """Floating-point operations of forward elimination plus back substitution on one right-hand side.

    Elimination: sum over k of 2 (n - k)(n - k + 1) ~ 2 n**3 / 3; back substitution: n**2.
    """
    raise NotImplementedError


def hilbert(n: int) -> np.ndarray:
    """The n x n Hilbert matrix H[i, j] = 1 / (i + j + 1): symmetric, positive definite and terribly conditioned."""
    i = np.arange(n)
    return 1.0 / (i[:, None] + i[None, :] + 1.0)


# ----------------------------------------------------------------------------- the truss

def two_bay_truss() -> tuple[np.ndarray, list[tuple[int, int]], list[tuple[int, int]]]:
    """A two-bay planar truss, 4 m long and 1.5 m tall: 6 nodes, 9 members, pin at node 0, roller at node 2.

    Nodes 0, 1, 2 along the bottom chord at x = 0, 2, 4 m; nodes 3, 4, 5 above them at y = 1.5 m.
    Members: bottom chord, top chord, three verticals, two diagonals (0-4 and 1-5).
    """
    nodes = np.array([[0.0, 0.0], [2.0, 0.0], [4.0, 0.0], [0.0, 1.5], [2.0, 1.5], [4.0, 1.5]])
    members = [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5), (0, 4), (1, 5)]
    supports = [(0, 0), (0, 1), (2, 1)]
    return nodes, members, supports


def truss_equilibrium_matrix(nodes: np.ndarray, members: list[tuple[int, int]],
                             supports: list[tuple[int, int]]) -> np.ndarray:
    """Equilibrium matrix A, shape (2 j, m + r): row 2 i + d is the d-component of force balance at node i.

    Column of member (i, k): +e at node i and -e at node k, with e the unit vector from i to k, so a positive
    (tension) force pulls each end toward the other. Column of reaction (node, dof): 1 in that row.
    """
    raise NotImplementedError


def truss_load_vector(n_nodes: int, loads: dict[int, tuple[float, float]]) -> np.ndarray:
    """Right-hand side -P of the joint equations: external loads enter with a minus sign."""
    raise NotImplementedError


def solve_truss(nodes: np.ndarray, members: list[tuple[int, int]], supports: list[tuple[int, int]],
                loads: dict[int, tuple[float, float]]) -> tuple[np.ndarray, np.ndarray]:
    """Member forces (tension positive, N) and reactions (N) of a statically determinate truss.

    Raise ``ValueError`` if the equilibrium matrix is not square (m + r != 2 j) or is singular (a mechanism).
    """
    raise NotImplementedError


def main(out_dir: str | Path = ".") -> None:
    """Solve the two-bay truss, verify every routine against numpy.linalg and save truss_forces.png."""
    import time

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    nodes, members, supports = two_bay_truss()
    loads = {1: (0.0, -15000.0)}          # a 1.5 t engine hangs from the lower chord's middle joint
    A = truss_equilibrium_matrix(nodes, members, supports)
    b = truss_load_vector(len(nodes), loads)
    forces, reactions = solve_truss(nodes, members, supports, loads)
    x_np = np.linalg.solve(A, b)
    print("Two-bay truss, 15 kN down at node 1 (tension positive):")
    for (i, k), f in zip(members, forces):
        print(f"  member {i}-{k}: {f / 1e3:+8.3f} kN")
    for (node, dof), rr in zip(supports, reactions):
        print(f"  reaction node {node} {'xy'[dof]}: {rr / 1e3:+8.3f} kN")
    print(f"  max |ours - numpy.linalg.solve| = {np.max(np.abs(np.concatenate([forces, reactions]) - x_np)):.2e} N")
    print(f"  relative residual ||b - A x|| / ||b|| = {residual_norm(A, np.concatenate([forces, reactions]), b):.2e}")
    print(f"  det A = {determinant(A):.6f}   (numpy: {np.linalg.det(A):.6f})")

    # Mechanism: remove a diagonal and the matrix becomes singular.
    try:
        solve_truss(nodes, members[:-1] + [(3, 5)], supports, loads)   # replace diagonal 1-5 by a duplicate top chord
    except ValueError as err:
        print(f"  without diagonal 1-5: {err}")

    # Verification against numpy on random systems.
    rng = np.random.default_rng(13)
    print("\nRandom systems, ours against numpy.linalg:")
    for n in (3, 10, 50):
        M = rng.standard_normal((n, n))
        v = rng.standard_normal(n)
        x = gaussian_elimination(M, v)
        P, L, U = lu_decompose(M)
        print(f"  n = {n:3d}: |x - numpy| = {np.max(np.abs(x - np.linalg.solve(M, v))):.1e}, "
              f"|PA - LU| = {np.max(np.abs(P @ M - L @ U)):.1e}, |det - numpy| / |det| = "
              f"{abs(determinant(M) - np.linalg.det(M)) / abs(np.linalg.det(M)):.1e}, "
              f"|A A^-1 - I| = {np.max(np.abs(M @ inverse_gauss_jordan(M) - np.eye(n))):.1e}")

    # Why pivot: a tiny pivot without swapping destroys the answer.
    T = np.array([[1e-13, 1.0], [1.0, 1.0]])
    t = np.array([1.0, 2.0])
    print(f"\nTiny pivot: without pivoting x = {gaussian_elimination(T, t, pivot=False)}, with pivoting x = {gaussian_elimination(T, t)}, exact ~ [1, 1]")
    print(f"  no pivoting but 2 steps of iterative refinement: x = {iterative_refinement(T, t, steps=2, pivot=False)}")

    # Conditioning: the Hilbert matrix.
    print("\nHilbert matrix H x = H 1 (exact solution: all ones):")
    for n in (4, 8, 10):
        H = hilbert(n)
        ones = np.ones(n)
        x = gaussian_elimination(H, H @ ones)
        x_ref = iterative_refinement(H, H @ ones, steps=3)
        print(f"  n = {n:2d}: cond = {np.linalg.cond(H):.1e}, residual = {residual_norm(H, x, H @ ones):.1e}, "
              f"max error = {np.max(np.abs(x - ones)):.1e}, after refinement = {np.max(np.abs(x_ref - ones)):.1e}")

    # Cost: operation count and wall time against LAPACK.
    print("\nCost: 2 n^3 / 3 operations, and a Python loop against LAPACK:")
    for n in (100, 200, 400):
        M = rng.standard_normal((n, n))
        v = rng.standard_normal(n)
        t0 = time.perf_counter()
        gaussian_elimination(M, v)
        t1 = time.perf_counter()
        np.linalg.solve(M, v)
        t2 = time.perf_counter()
        print(f"  n = {n}: {elimination_flops(n):.2e} flops, ours {1e3 * (t1 - t0):7.1f} ms, numpy {1e3 * (t2 - t1):6.2f} ms")

    # Figure: the truss with member forces coloured by sign.
    fig, ax = plt.subplots(figsize=(7, 3.5))
    fmax = np.max(np.abs(forces))
    for (i, k), f in zip(members, forces):
        xs, ys = nodes[[i, k], 0], nodes[[i, k], 1]
        ax.plot(xs, ys, color="C3" if f > 0 else "C0", lw=1.0 + 4.0 * abs(f) / fmax)
        ax.text(xs.mean(), ys.mean(), f"{f / 1e3:+.1f}", fontsize=8, ha="center", va="bottom")
    ax.plot(nodes[:, 0], nodes[:, 1], "ko")
    ax.annotate("", xy=(2.0, -0.6), xytext=(2.0, 0.0), arrowprops=dict(arrowstyle="->", lw=2))
    ax.text(2.1, -0.5, "15 kN")
    ax.set_aspect("equal")
    ax.set_xlim(-0.5, 4.5)
    ax.set_ylim(-0.9, 2.0)
    ax.set_title("Member forces [kN], red = tension, blue = compression")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_dir / "truss_forces.png", dpi=150)
    print(f"\nsaved {out_dir / 'truss_forces.png'}")


if __name__ == "__main__":
    main()
