"""P01.14 Vector Spaces, Basis & Rank: reference solution.

Reduced row echelon form, rank, null space and column space written from scratch in NumPy and
checked against ``numpy.linalg.matrix_rank`` and an SVD null space; applied to the equilibrium
matrix of a planar truss to detect mechanisms (left null space) and states of self-stress
(null space), which Maxwell's counting rule alone cannot see.

Conventions:
    * Matrices are float arrays (m, n); bases are returned as matrices whose *columns* are the basis vectors.
    * ``tol`` is a relative tolerance: an entry below tol * max|A| is treated as zero during elimination.
    * Truss: nodes (j, 2) in metres; members as (i, k) index pairs; supports as (node, dof), dof 0 = x, 1 = y.
      Equilibrium matrix A is (2 j) x (m + r) with member columns +e at node i, -e at node k (tension positive)
      and a 1 for each reaction. Mechanisms are vectors d with A^T d = 0; self-stress states are x with A x = 0.

Run:  python solution.py          (prints the classification of four trusses, writes truss_mechanism.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-01/lesson-14-vector-spaces-basis-rank
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

DEFAULT_TOL = 1e-10


# ----------------------------------------------------------------------------- row reduction

def rref(A: np.ndarray, tol: float = DEFAULT_TOL) -> tuple[np.ndarray, list[int]]:
    """Reduced row echelon form R of A and the list of pivot columns.

    Gauss-Jordan with partial pivoting: in each column the largest remaining |entry| is the pivot; a column whose
    largest remaining entry is below ``tol * max|A|`` has no pivot (it is a free column). Pivot rows are scaled
    to 1 and the pivot column is cleared above and below. Returns a new array; A is not modified.
    """
    R = np.array(A, dtype=float, copy=True)
    if R.ndim != 2:
        raise ValueError("A must be a 2-D array")
    m, n = R.shape
    scale = float(np.max(np.abs(R))) if R.size else 0.0
    if scale == 0.0:
        return R, []
    pivots: list[int] = []
    row = 0
    for col in range(n):
        if row >= m:
            break
        p = row + int(np.argmax(np.abs(R[row:, col])))
        if abs(R[p, col]) <= tol * scale:
            R[row:, col] = 0.0                   # numerically zero column below the current row
            continue
        if p != row:
            R[[row, p]] = R[[p, row]]
        R[row] /= R[row, col]
        for i in range(m):
            if i != row and R[i, col] != 0.0:
                R[i] -= R[i, col] * R[row]
        pivots.append(col)
        row += 1
    R[np.abs(R) <= tol * scale] = 0.0
    return R, pivots


def rank(A: np.ndarray, tol: float = DEFAULT_TOL) -> int:
    """Rank of A: the number of pivots in its reduced row echelon form."""
    return len(rref(A, tol)[1])


def nullspace(A: np.ndarray, tol: float = DEFAULT_TOL) -> np.ndarray:
    """Basis of the null space {x : A x = 0} as the columns of an (n, n - rank) array (the special solutions).

    One special solution per free column: set that free variable to 1, the other free variables to 0, and read
    the pivot variables off the reduced rows. Returns shape (n, 0) when A has full column rank.
    """
    R, pivots = rref(A, tol)
    n = R.shape[1]
    free = [c for c in range(n) if c not in pivots]
    N = np.zeros((n, len(free)))
    for j, f in enumerate(free):
        N[f, j] = 1.0
        for i, p in enumerate(pivots):
            N[p, j] = -R[i, f]
    return N


def left_nullspace(A: np.ndarray, tol: float = DEFAULT_TOL) -> np.ndarray:
    """Basis of {y : A^T y = 0} (equivalently y^T A = 0) as columns of an (m, m - rank) array."""
    return nullspace(np.asarray(A, dtype=float).T, tol)


def column_space_basis(A: np.ndarray, tol: float = DEFAULT_TOL) -> np.ndarray:
    """The pivot columns of A, as an (m, rank) array: a basis of the column space made of columns of A itself."""
    A = np.asarray(A, dtype=float)
    return A[:, rref(A, tol)[1]]


def row_space_basis(A: np.ndarray, tol: float = DEFAULT_TOL) -> np.ndarray:
    """The non-zero rows of rref(A), returned as the columns of an (n, rank) array: a basis of the row space."""
    R, pivots = rref(A, tol)
    return R[:len(pivots)].T


def is_independent(vectors: np.ndarray, tol: float = DEFAULT_TOL) -> bool:
    """True if the columns of ``vectors`` are linearly independent (rank equals the number of columns)."""
    V = np.asarray(vectors, dtype=float)
    if V.ndim != 2:
        raise ValueError("pass the vectors as the columns of a 2-D array")
    return rank(V, tol) == V.shape[1]


def four_subspace_dimensions(A: np.ndarray, tol: float = DEFAULT_TOL) -> dict[str, int]:
    """Dimensions of the row space, column space, null space and left null space of an (m, n) matrix.

    Keys: 'm', 'n', 'rank', 'nullity' (= n - rank), 'left_nullity' (= m - rank).
    """
    A = np.asarray(A, dtype=float)
    r = rank(A, tol)
    return {"m": A.shape[0], "n": A.shape[1], "rank": r, "nullity": A.shape[1] - r, "left_nullity": A.shape[0] - r}


def is_consistent(A: np.ndarray, b: np.ndarray, tol: float = DEFAULT_TOL) -> bool:
    """True if b lies in the column space of A, i.e. rank([A | b]) == rank(A)."""
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float).reshape(-1, 1)
    return rank(np.hstack([A, b]), tol) == rank(A, tol)


def solve_general(A: np.ndarray, b: np.ndarray, tol: float = DEFAULT_TOL) -> tuple[np.ndarray, np.ndarray]:
    """All solutions of A x = b as (x_particular, N): x = x_particular + N c for any c.

    x_particular has the free variables set to zero. Raise ``ValueError`` if the system is inconsistent.
    """
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float)
    if b.shape != (A.shape[0],):
        raise ValueError(f"b must have shape ({A.shape[0]},), got {b.shape}")
    aug = np.hstack([A, b.reshape(-1, 1)])
    R, pivots = rref(aug, tol)
    n = A.shape[1]
    if n in pivots:
        raise ValueError("inconsistent system: b is not in the column space of A")
    x = np.zeros(n)
    for i, p in enumerate(pivots):
        x[p] = R[i, n]
    return x, nullspace(A, tol)


# ----------------------------------------------------------------------------- the truss

def truss_equilibrium_matrix(nodes: np.ndarray, members: list[tuple[int, int]],
                             supports: list[tuple[int, int]]) -> np.ndarray:
    """Equilibrium matrix A, shape (2 j, m + r): row 2 i + d is the d-component of force balance at node i.

    Member column: +e in the two rows of node i and -e in those of node k (e the unit vector from i to k);
    reaction column: a single 1 in the row of its node and direction.
    """
    nodes = np.asarray(nodes, dtype=float)
    j = nodes.shape[0]
    m, r = len(members), len(supports)
    A = np.zeros((2 * j, m + r))
    for col, (i, k) in enumerate(members):
        d = nodes[k] - nodes[i]
        length = np.linalg.norm(d)
        if length == 0.0:
            raise ValueError(f"member {col} joins coincident nodes {i} and {k}")
        e = d / length
        A[2 * i:2 * i + 2, col] += e
        A[2 * k:2 * k + 2, col] -= e
    for col, (node, dof) in enumerate(supports):
        if dof not in (0, 1):
            raise ValueError("support dof must be 0 (x) or 1 (y)")
        A[2 * node + dof, m + col] = 1.0
    return A


def truss_variant(name: str) -> tuple[np.ndarray, list[tuple[int, int]], list[tuple[int, int]]]:
    """Four two-bay trusses on the same six nodes (bottom 0, 1, 2 at y = 0; top 3, 4, 5 at y = 1.5 m; 2 m bays).

    'determinate'      : chords, verticals, diagonals 0-4 and 1-5 (9 members)        -> rigid, no redundancy.
    'missing-diagonal' : the same without 1-5 (8 members)                              -> the right bay is a mechanism.
    'extra-diagonal'   : 'determinate' plus 1-3 (10 members)                           -> the left bay is redundant.
    'maxwell-trap'     : both diagonals 0-4 and 1-3 in the left bay, none in the right -> 9 members, still a mechanism.
    Supports in every case: pin at node 0, roller (y) at node 2.
    """
    nodes = np.array([[0.0, 0.0], [2.0, 0.0], [4.0, 0.0], [0.0, 1.5], [2.0, 1.5], [4.0, 1.5]])
    frame = [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)]
    variants = {
        "determinate": frame + [(0, 4), (1, 5)],
        "missing-diagonal": frame + [(0, 4)],
        "extra-diagonal": frame + [(0, 4), (1, 5), (1, 3)],
        "maxwell-trap": frame + [(0, 4), (1, 3)],
    }
    if name not in variants:
        raise ValueError(f"unknown variant {name!r}; choose from {sorted(variants)}")
    return nodes, variants[name], [(0, 0), (0, 1), (2, 1)]


def maxwell_count(m: int, r: int, j: int) -> int:
    """Maxwell's count m + r - 2 j for a planar pin-jointed truss: equals s - k (self-stress states minus mechanisms)."""
    return m + r - 2 * j


def truss_mechanisms(A: np.ndarray, tol: float = DEFAULT_TOL) -> np.ndarray:
    """Basis of the mechanisms: joint displacement vectors d (2 j, k) with A^T d = 0.

    A^T d = 0 says every member has zero elongation and every restrained degree of freedom is zero, yet d is
    not zero: the truss can move without stretching anything, at least to first order.
    """
    return left_nullspace(A, tol)


def truss_self_stress(A: np.ndarray, tol: float = DEFAULT_TOL) -> np.ndarray:
    """Basis of the states of self-stress: force vectors x (m + r, s) with A x = 0 (equilibrium with no load)."""
    return nullspace(A, tol)


def member_elongations(nodes: np.ndarray, members: list[tuple[int, int]], d: np.ndarray) -> np.ndarray:
    """First-order elongation of every member under the joint displacement vector d (2 j,): e . (d_k - d_i)."""
    nodes = np.asarray(nodes, dtype=float)
    d = np.asarray(d, dtype=float).reshape(-1, 2)
    out = np.empty(len(members))
    for col, (i, k) in enumerate(members):
        e = nodes[k] - nodes[i]
        e = e / np.linalg.norm(e)
        out[col] = e @ (d[k] - d[i])
    return out


def classify_truss(nodes: np.ndarray, members: list[tuple[int, int]], supports: list[tuple[int, int]],
                   tol: float = DEFAULT_TOL) -> dict[str, object]:
    """Rank analysis of a truss: counts, rank, number of mechanisms k and self-stress states s, and a verdict.

    Keys: 'm', 'r', 'j', 'maxwell' (m + r - 2 j), 'rank', 'k', 's', 'status' in
    {'determinate', 'mechanism', 'indeterminate', 'mechanism and indeterminate'}.
    """
    A = truss_equilibrium_matrix(nodes, members, supports)
    j, m, r = len(nodes), len(members), len(supports)
    rk = rank(A, tol)
    k, s = 2 * j - rk, m + r - rk
    if k == 0 and s == 0:
        status = "determinate"
    elif k > 0 and s == 0:
        status = "mechanism"
    elif k == 0 and s > 0:
        status = "indeterminate"
    else:
        status = "mechanism and indeterminate"
    return {"m": m, "r": r, "j": j, "maxwell": maxwell_count(m, r, j), "rank": rk, "k": k, "s": s, "status": status}


def main(out_dir: str | Path = ".") -> None:
    """Classify the four trusses, verify rank and null spaces against numpy, and draw the mechanism."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    print("Rank analysis of four two-bay trusses (j = 6 joints, r = 3 reactions):")
    print(f"  {'variant':18s} {'m':>2s} {'m+r-2j':>7s} {'rank':>4s} {'k':>2s} {'s':>2s}  verdict          numpy rank")
    results = {}
    for name in ("determinate", "missing-diagonal", "extra-diagonal", "maxwell-trap"):
        nodes, members, supports = truss_variant(name)
        info = classify_truss(nodes, members, supports)
        A = truss_equilibrium_matrix(nodes, members, supports)
        results[name] = (nodes, members, supports, A, info)
        print(f"  {name:18s} {info['m']:2d} {info['maxwell']:7d} {info['rank']:4d} {info['k']:2d} {info['s']:2d}  "
              f"{info['status']:28s} {np.linalg.matrix_rank(A)}")

    # The mechanism of the Maxwell trap: zero elongation everywhere, zero at the supports, yet a real motion.
    nodes, members, supports, A, info = results["maxwell-trap"]
    D = truss_mechanisms(A)
    d = D[:, 0] / np.max(np.abs(D[:, 0]))
    print("\nMaxwell trap, mechanism mode (joint displacements, normalised):")
    for i in range(len(nodes)):
        print(f"  node {i}: ({d[2 * i]:+.3f}, {d[2 * i + 1]:+.3f})")
    print(f"  max |A^T d| = {np.max(np.abs(A.T @ d)):.1e}, max member elongation = {np.max(np.abs(member_elongations(nodes, members, d))):.1e}")
    S = truss_self_stress(A)
    x = S[:, 0] / np.max(np.abs(S[:, 0]))
    print("Maxwell trap, state of self-stress (member forces, normalised, tension positive):")
    for (i, k), f in zip(members, x[:len(members)]):
        if abs(f) > 1e-12:
            print(f"  member {i}-{k}: {f:+.3f}")
    print(f"  max |A x| = {np.max(np.abs(A @ x)):.1e}")

    # Compare the null space with an SVD null space (numpy): same subspace, so the projectors agree.
    _, sv, vt = np.linalg.svd(A)
    N_svd = vt[np.sum(sv > 1e-10 * sv[0]):].T
    P_ours = S @ np.linalg.pinv(S)
    P_svd = N_svd @ N_svd.T
    print(f"  null-space projector, ours vs SVD: max difference {np.max(np.abs(P_ours - P_svd)):.1e}")

    # A random rank-deficient matrix, as a generic check.
    rng = np.random.default_rng(14)
    B = rng.standard_normal((7, 4)) @ rng.standard_normal((4, 9))
    R, piv = rref(B)
    print(f"\nRandom 7 x 9 matrix of rank 4: rref pivots {piv}, our rank {rank(B)}, numpy rank {np.linalg.matrix_rank(B)}, "
          f"nullity {nullspace(B).shape[1]}, left nullity {left_nullspace(B).shape[1]}, |B N| = {np.max(np.abs(B @ nullspace(B))):.1e}")

    # Figure: the determinate truss and the trap's mechanism, exaggerated.
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.2))
    for ax, name, title in ((axes[0], "determinate", "determinate: rank 12, no mechanism"),
                            (axes[1], "maxwell-trap", "Maxwell trap: m + r = 2j but rank 11")):
        nodes, members, supports, A, info = results[name]
        for i, k in members:
            ax.plot(nodes[[i, k], 0], nodes[[i, k], 1], color="0.6", lw=1.5)
        if name == "maxwell-trap":
            moved = nodes + 0.3 * d.reshape(-1, 2)
            for i, k in members:
                ax.plot(moved[[i, k], 0], moved[[i, k], 1], color="C3", lw=1.5)
            ax.plot(moved[:, 0], moved[:, 1], "o", color="C3", ms=4)
        ax.plot(nodes[:, 0], nodes[:, 1], "ko", ms=4)
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=10)
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_dir / "truss_mechanism.png", dpi=150)
    print(f"\nsaved {out_dir / 'truss_mechanism.png'}")


if __name__ == "__main__":
    main()
