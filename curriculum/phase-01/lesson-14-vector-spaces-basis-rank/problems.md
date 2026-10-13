# P01.14  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. Bases are not unique: when you compare yours with the code's, compare the *subspaces* (does each of your vectors satisfy the defining equation, and do you have the right number?), not the numbers.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: the four subspaces by hand)

A = [[1, 2, 3, s], [2, 4, 6, 2s], [0, 1, 1, 1]].

(a) Row-reduce to R by hand. List the pivot columns, the free columns and the rank.
(b) Write a basis of each of the four fundamental subspaces, with its dimension. For the left null space, find the row combination that produced the zero row.
(c) Check with `rref`, `nullspace`, `column_space_basis`, `left_nullspace` and `four_subspace_dimensions`. If your null-space basis differs from the code's, show that each of your vectors is in the code's span.
(d) For which right-hand sides b is A x = b solvable? Give the condition on (b₁, b₂, b₃) and verify it with `is_consistent` on one solvable and one unsolvable b.

Answer form: a 3 × 4 matrix; four bases with dimensions; one condition; two booleans.

## 2. (Derivation: rank–nullity and Maxwell's rule)

(a) Prove that the special solutions of A x = 0 built from rref(A) are linearly independent and span the null space. Conclude dim N(A) = n − r.
(b) Prove that every vector in N(A) is orthogonal to every vector in the row space, and that N(A) ∩ C(Aᵀ) = {0}.
(c) Show that rank(Aᵀ) = rank(A). (Hint: the row space of A has dimension r and is the column space of Aᵀ.)
(d) For a planar truss with j joints, m members and r reactions, derive m + r − 2j = s − k where s = dim N(A) and k = dim N(Aᵀ), and explain in one sentence each what a vector of N(A) and a vector of N(Aᵀ) mean physically.

Answer form: four short proofs with each step named.

## 3. (Unfamiliar scenario: reaction wheels on a satellite)

A satellite controls its attitude with reaction wheels. Wheel i, spun up with torque τ_i about its axis **a**_i, applies torque −τ_i **a**_i to the satellite body. The achievable body torques are **T** = B **τ** with B = [**a**₁ **a**₂ ⋯] (3 × n_wheels, columns of unit length), up to the sign.

(a) Three wheels with axes (1, 0, 0), (0, 1, 0), (1, 1, 1)/√3. Find the rank of B and state whether every body torque is achievable.
(b) Three wheels with axes (1, 0, 0), (0, 1, 0), (1, 1, 0)/√2. Find the rank, a basis of N(B) and of N(Bᵀ). Which torque direction cannot be produced, and what does the null-space vector say about wheel 3?
(c) Four wheels: (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)/√3. Find N(B). Use `solve_general` to write every wheel-torque vector that produces the body torque (0.1, 0.2, 0.3) N m. Why is a non-trivial null space *useful* here? (Think about what happens to the wheels' stored momentum over months.)
(d) Any one of the four wheels in (c) fails. Which single failures leave the satellite fully controllable?

Answer form: a rank and a yes/no; a rank, two vectors and two sentences; a general solution and a paragraph; a list.

## 4. (Numerical, uses your code: what numerical rank means)

Take the determinate two-bay truss and move node 4 from y = 1.5 m to 0.15, 0.015, 1.5 × 10⁻³, 1.5 × 10⁻⁴, 1.5 × 10⁻⁵ and 1.5 × 10⁻⁶ m.

(a) For each geometry tabulate `rank(A)` at tol = 10⁻¹⁰, 10⁻⁶ and 10⁻³, `np.linalg.matrix_rank(A)`, the smallest singular value σ_min from `np.linalg.svd(A, compute_uv=False)`, and κ = σ_max/σ_min.
(b) Solve for the member forces under 15 kN at node 1 (use `solve_general` or the solver from `P01.13`) and add the largest force to the table.
(c) Explain why the rref rank at tol = 10⁻¹⁰ never drops although σ_min reaches 3.6 × 10⁻⁷, and why `matrix_rank` does not drop either. Which quantity in the table tracks the physical danger, and how does it scale with y₄?
(d) At which y₄ does the rank at tol = 10⁻³ drop, and what does that threshold mean geometrically?

Answer form: a table of seven rows; one paragraph; one number with its meaning.

## 5. (Repair a mechanism)

Remove member number (s mod 9) from the determinate truss's member list (0-based; the list order is chords 0-1, 1-2, 3-4, 4-5, verticals 0-3, 1-4, 2-5, diagonals 0-4, 1-5).

(a) Classify the result with `classify_truss`. If it is a mechanism, print the mode with `truss_mechanisms`, normalise it, and describe the motion in words.
(b) Verify with `member_elongations` that no member stretches and that the supported degrees of freedom are zero.
(c) Try every member you could add between two currently unconnected joints, one at a time, and classify each. List the ones that restore a determinate truss and the ones that do not, and explain the pattern geometrically.
(d) Among the repairs, which would you choose, and why? (Length, buckling of a compression member, symmetry.)

Answer form: a classification and a description; two residuals; two lists with an explanation; one judgment.

## 6. (Judgment: what rank can and cannot tell you)

For each statement, say whether it is justified, give the measurable quantity you would compute or measure to decide, and name the assumption or tool involved:

(a) "The rank test says k = 0, so the structure is rigid and safe."
(b) "The rank test says k = 1, so the structure will collapse."
(c) "`matrix_rank` returns full rank, so the member forces from the solver are reliable."
(d) "s = 2 means the structure is over-designed; remove two members."

Answer form: four short paragraphs, each with a number or a measurable quantity.

## 7. (Hard, uses your code: a parity check for a faulty sensor)

A wind-tunnel balance measures three loads (normal force N, axial force A, pitching moment M) with four strain gauges whose readings are **y** = H **f** + **e**, with

H = [[1, 0.5, 0], [0, 1, 0.5], [0.5, 0, 1], [1, 1, 1]] (gauge sensitivities, dimensionless after calibration),

**f** = (N, A, M) and **e** the gauge errors.

(a) Show that rank H = 3 and find a basis vector **y**⊥ of the left null space. Explain why **y**⊥ᵀ **y** does not depend on the loads at all: it is a **parity** residual that is zero for a healthy balance whatever the loads.
(b) Take **f** = (2, −1, 3) and compute the clean readings and the parity residual. Then add a fault of +0.5 to gauge 2 and recompute. Then to gauge 4.
(c) The fault "signature" of gauge i is **y**⊥ᵀ **e**_i. Compute all four. Which faults can be *detected*, and which can be *isolated* (told apart from the others)? Why does one parity relation never suffice to isolate among three gauges?
(d) Propose a fifth gauge row that makes every single-gauge fault isolable, compute the new left null space (now two-dimensional) and the 2-vector signatures, and verify that no two signatures are parallel.

Answer form: a rank and a vector with an argument; three residuals; four signatures and a paragraph; a 5 × 3 matrix with five 2-vector signatures.
