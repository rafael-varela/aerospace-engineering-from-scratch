# P01.14  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) Row 2 ← row 2 − 2 row 1 (becomes zero); swap the zero row to the bottom; pivot in column 2 of the remaining row; clear column 2 from row 1. Pivot columns {1, 2}, free {3, 4}, rank 2. (b) Null space: two special solutions (free x₃ = 1 and free x₄ = 1). Column space: columns 1 and 2 of A. Row space: the two non-zero rows of R. Left null space: the combination of rows that gave the zero row (row 2 − 2 row 1). (c) `rref(A)`, `nullspace(A)`, etc.; to show one of your vectors is in the code's span, solve `N c = v` by least squares and check the residual, or check `rank([N, v]) == rank(N)`. (d) b must be orthogonal to the left null vector: b₂ − 2b₁ = 0.

## 2.
(a) Independence: restrict any combination Σ c_j **s**_j to the free rows; the identity block gives c_j = 0. Spanning: for **x** ∈ N(A), set **z** = **x** − Σ x_{f_j} **s**_j; **z** is in N(A) with all free components zero; the reduced equations then force every pivot component to zero, so **z** = 0. (b) A x = 0 says row_i · x = 0 for every row; a vector in both subspaces is orthogonal to itself. (c) C(Aᵀ) = row space of A, dimension r by rref (the non-zero rows of R are independent and span it); so rank Aᵀ = r. (d) A is (2j) × (m + r): dim N(A) = m + r − rank = s, dim N(Aᵀ) = 2j − rank = k; subtract. N(A): member forces (and reactions) in equilibrium with zero load. N(Aᵀ): joint displacements with zero first-order elongation in every member and zero at every support.

## 3.
(a) `rank(B)`; full rank 3 ⇒ C(B) = R³. (b) `rank`, `nullspace(B)`, `left_nullspace(B)`; the left null vector is the unreachable torque direction; the null vector gives wheel torques that cancel (wheel 3 is a combination of 1 and 2). (c) `solve_general(B, [0.1, 0.2, 0.3])` returns x_p and N; the null direction changes wheel speeds without torquing the body, which lets you redistribute stored momentum. (d) Delete each column in turn and compute the rank.

## 4.
(a) Loop over y₄; build `nodes` with node 4 moved; `truss_equilibrium_matrix`; `rank(A, tol)` for the three tolerances; `np.linalg.matrix_rank`; `np.linalg.svd(A, compute_uv=False)`. (b) `solve_general(A, -P)` with P the load vector (or the P01.13 solver); `np.max(np.abs(x))`. (c) Compare the pivot sizes (before scaling) with σ_min; relate the force growth to 1/σ_min and 1/y₄. (d) Read the row where the tol = 10⁻³ rank first becomes 11; compare y₄/(diagonal length) with 10⁻³.

## 5.
(a) `members = list(members); del members[s % 9]`; `classify_truss`; `truss_mechanisms(A)`; divide by its max absolute entry; print per node. (b) `member_elongations(nodes, members, d)`; `d[2*node + dof]` for each support. (c) `itertools.combinations(range(6), 2)` minus existing members; `classify_truss(nodes, members + [cand], supports)` for each. (d) Compare the candidates' lengths and whether they would be in compression under the design load.

## 6.
For each statement name which of the three questions (exact rank, distance to singularity, load path and stability) it answers and which it silently assumes; name the number to compute (σ_min, κ, a buckling load, a second-order analysis) or the test to run.

## 7.
(a) `rank(H)`; `left_nullspace(H)` gives **y**⊥; expand **y**⊥ᵀ(H f + e). (b) `y = H @ f`; residual `Y.T @ y`; add 0.5 to the chosen component and recompute. (c) `Y.T @ np.eye(4)`; detectable if the signature is non-zero; isolable if no other gauge has a parallel signature. (d) Append a row, recompute `left_nullspace` (shape (5, 2)), form the 2 × 5 signature matrix, check every 2 × 2 determinant of pairs of columns is non-zero.
