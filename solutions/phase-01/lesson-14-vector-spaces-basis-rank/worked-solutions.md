# P01.14  Worked solutions

Every number below is produced by `solution.py` (the matrices are written out explicitly in the text, so you can reproduce them without a seed). Where a seed is needed the solutions use **s = 15**; substitute your own. Bases are not unique: compare subspaces, not numbers.

## 1. The four subspaces by hand (s = 15)

**Assumptions.** A = [[1, 2, 3, 15], [2, 4, 6, 30], [0, 1, 1, 1]].

**Derivation.**
(a) Row 2 − 2 row 1 → (0, 0, 0, 0). Swap it to the bottom; the second row is now (0, 1, 1, 1), already a pivot row with pivot in column 2. Row 1 − 2 row 2 → (1, 0, 1, 13).

R = [[1, 0, 1, 13], [0, 1, 1, 1], [0, 0, 0, 0]], pivot columns {1, 2} (0-based [0, 1]), free columns {3, 4}, rank 2. In general the (1, 4) entry is s − 2.

(b) Null space (dimension n − r = 2). Free x₃ = 1, x₄ = 0: pivot rows give x₁ = −1, x₂ = −1 → (−1, −1, 1, 0). Free x₄ = 1, x₃ = 0: x₁ = −13, x₂ = −1 → (−13, −1, 0, 1) (general: (−(s − 2), −1, 0, 1)).
Column space (dimension 2): columns 1 and 2 of **A**, (1, 2, 0) and (2, 4, 1).
Row space (dimension 2): (1, 0, 1, 13) and (0, 1, 1, 1).
Left null space (dimension m − r = 1): the zero row came from row 2 − 2 row 1, so **y** = (−2, 1, 0); check **y**ᵀA = (−2 + 2, −4 + 4, −6 + 6, −30 + 30) = 0.

(c) `rref(A)` returns exactly the R above with pivots `[0, 1]`; `nullspace(A)` returns the two special solutions as columns; `column_space_basis(A)` returns columns (1, 2, 0), (2, 4, 1); `left_nullspace(A)` returns (−2, 1, 0); `four_subspace_dimensions(A)` → `{'m': 3, 'n': 4, 'rank': 2, 'nullity': 2, 'left_nullity': 1}`. If your null vectors are, say, (1, 1, −1, 0) and (−12, 0, −1, 1), check `rank(np.column_stack([N, v])) == 2` for each.

(d) **b** must be orthogonal to the left null space: −2 b₁ + b₂ = 0, i.e. b₂ = 2 b₁ (b₃ is free because the third row is independent). `is_consistent(A, [1, 2, 5])` → `True`; `is_consistent(A, [1, 3, 5])` → `False`.

**Result.** rank 2; dims (row, column, null, left-null) = (2, 2, 2, 1); condition b₂ = 2 b₁.

**Check.** 2 + 2 = 4 = n and 2 + 1 = 3 = m. Each null vector dotted with each row of R is zero: (1, 0, 1, 13)·(−13, −1, 0, 1) = −13 + 13 = 0.

**Common mistakes.** Writing the null vector with +13 (forgetting the sign flip of the rref entry); giving the rref rows as a basis of the *column* space (they span the row space; the column space of R is not the column space of A — compare column 3 of A, (3, 6, 1), with column 3 of R, (1, 1, 0)).

## 2. Rank–nullity and Maxwell's rule

**(a) Special solutions are a basis of N(A).** Let the free columns be f₁ < … < f_{n−r} and **s**_j the special solution with (s_j)_{f_j} = 1, (s_j)_{f_i} = 0 for i ≠ j, and pivot components read from −R. *Independence:* suppose Σ c_j **s**_j = 0. Look at component f_i: only **s**_i has a non-zero entry there (it is 1), so c_i = 0 for every i. *Spanning:* take any **x** ∈ N(A) and form **z** = **x** − Σ x_{f_j} **s**_j. **z** ∈ N(A) (a combination of null vectors) and every free component of **z** is zero by construction. Row i of R **z** = 0 reads z_{p_i} + Σ_free R_{i f} z_f = z_{p_i} = 0, so every pivot component vanishes too: **z** = 0 and **x** = Σ x_{f_j} **s**_j. Hence dim N(A) = number of free columns = n − r. (R **x** = 0 ⇔ A **x** = 0 because R = E A with E invertible.)

**(b) N(A) ⊥ C(Aᵀ).** A **x** = 0 says every row **a**_iᵀ **x** = 0, so **x** is orthogonal to each row, hence to every combination of rows, i.e. to the row space C(Aᵀ). If **v** is in both, **v**ᵀ**v** = 0 (take the row-space vector to be **v** itself), so **v** = 0.

**(c) rank Aᵀ = rank A.** The row space of A equals the row space of R (row operations are invertible combinations of rows). The r non-zero rows of R are independent (each has a 1 in its pivot column where the others have 0), so dim C(Aᵀ) = r. But C(Aᵀ) is, by definition, the column space of Aᵀ, whose dimension is rank Aᵀ. So rank Aᵀ = r = rank A.

**(d) Maxwell.** A is (2j) × (m + r): one row per joint equilibrium equation, one column per unknown force. rank–nullity in the columns: s = dim N(A) = (m + r) − rank. In the rows (applied to Aᵀ, using (c)): k = dim N(Aᵀ) = 2j − rank. Subtract: s − k = m + r − 2j. ∎
A vector of N(A) is a set of member forces and reactions in equilibrium with *no* applied load: a self-stress state (the structure can be pre-tensioned). A vector of N(Aᵀ) is a set of joint displacements that stretch no member to first order and move no support: a mechanism.

**Common mistakes.** Proving spanning by "there are n − r of them and the dimension is n − r" (circular — the dimension is what is being established); forgetting that rank Aᵀ = rank A is needed to write k = 2j − rank.

## 3. Reaction wheels on a satellite

**(a)** B₁ = [[1, 0, 0.577], [0, 1, 0.577], [0, 0, 0.577]]. `rank(B1)` = 3: C(B₁) = R³, every body torque is achievable (the third wheel supplies the only z-component).

**(b)** B₂ = [[1, 0, 0.707], [0, 1, 0.707], [0, 0, 0]]. `rank(B2)` = 2. N(B₂) is spanned by (−0.707, −0.707, 1) (from `nullspace`): spinning wheel 3 up with torque τ while braking wheels 1 and 2 by 0.707 τ each produces no body torque — wheel 3 is a combination of wheels 1 and 2. N(B₂ᵀ) is spanned by (0, 0, 1): the z-row of B₂ is zero, so no torque about z can ever be produced; the satellite cannot be controlled in yaw about that axis.

**(c)** B₃ = [[1, 0, 0, 0.577], [0, 1, 0, 0.577], [0, 0, 1, 0.577]], rank 3, N(B₃) spanned by **n** = (−0.577, −0.577, −0.577, 1). `solve_general(B3, [0.1, 0.2, 0.3])` → x_p = (0.1, 0.2, 0.3, 0) N m and the general solution **τ** = (0.1, 0.2, 0.3, 0) + c (−0.577, −0.577, −0.577, 1). Over months, disturbance torques (gravity gradient, solar pressure) are absorbed as wheel momentum; a wheel that saturates can no longer deliver torque. The null direction lets the controller shift momentum between wheels (e.g. unload wheel 4 while loading 1–3) *without* disturbing the attitude — so a 4-wheel set can balance its momentum distribution. Without a null space (three wheels) every momentum change is also an attitude torque.

**(d)** Deleting any single column of B₃ leaves rank 3 (`rank(np.delete(B3, i, axis=1))` = 3 for i = 0…3): any three of the four axes span R³, because (1, 1, 1)/√3 has a non-zero component along every coordinate axis. All four single failures keep full controllability — this is why the skewed fourth wheel is the standard redundant layout.

**Common mistakes.** Treating the sign of **T** = −B **τ** as changing the subspaces (it does not: a sign flips the vectors, not the span); concluding from (b) that "wheel 3 is useless" — it is useful for redundancy in the x–y plane, it just adds nothing about z.

## 4. What numerical rank means

**Assumptions.** Determinate two-bay truss (nodes at (0, 0), (2, 0), (4, 0), (0, 1.5), (2, 1.5), (4, 1.5); members chords 0-1, 1-2, 3-4, 4-5, verticals 0-3, 1-4, 2-5, diagonals 0-4, 1-5; supports node 0 x, node 0 y, node 2 y), 15 kN downward at node 1; node 4 lowered to y₄.

**Result (a–b).**

| y₄ (m) | rank tol 10⁻¹⁰ | tol 10⁻⁶ | tol 10⁻³ | `matrix_rank` | σ_min | κ | max force (kN) |
|---|---|---|---|---|---|---|---|
| 1.5 | 12 | 12 | 12 | 12 | 0.288 | 7.4 | 12.5 |
| 0.15 | 12 | 12 | 12 | 12 | 0.0370 | 58.6 | 125 |
| 0.015 | 12 | 12 | 12 | 12 | 3.57 × 10⁻³ | 606 | 1 250 |
| 1.5 × 10⁻³ | 12 | 12 | 12 | 12 | 3.56 × 10⁻⁴ | 6 080 | 12 500 |
| 1.5 × 10⁻⁴ | 12 | 12 | 11 | 12 | 3.56 × 10⁻⁵ | 6.08 × 10⁴ | 125 000 |
| 1.5 × 10⁻⁵ | 12 | 12 | 11 | 12 | 3.56 × 10⁻⁶ | 6.08 × 10⁵ | 1.25 × 10⁶ |
| 1.5 × 10⁻⁶ | 12 | 12 | 11 | 12 | 3.56 × 10⁻⁷ | 6.08 × 10⁶ | 1.25 × 10⁷ |

The maximum force sits in the two diagonals (and the top chord), and scales exactly as 1/y₄: at y₄ = 1.5 m the 12.5 kN diagonal carries the 7.5 kN shear with a sine of 0.6; at y₄ → 0 the diagonal becomes horizontal and must carry the same vertical component with a vanishing sine.

A practical note: `solve_general(A, -P)` with P in newtons raises "inconsistent" on the last row. The reason is instructive: `rref` scales its tolerance by the largest entry of the *augmented* matrix, which is now 15 000, so the tolerance is 1.5 × 10⁻⁶ and the 6 × 10⁻⁷ geometry entries are zeroed. Give the load in kilonewtons (`-P/1000`), or use the P01.13 solver, and the row is produced. Equilibrate your units before you trust a tolerance.

**(c)** rref divides each pivot row by its pivot. A pivot of 6 × 10⁻⁷ (the vertical component of a diagonal's unit vector, y₄/2.5) is still 6 × 10⁻⁷ ≫ 10⁻¹⁰ · max|A| = 10⁻¹⁰, so it is accepted and the rank stays 12; with tol = 10⁻⁶ the comparison narrowly passes even at y₄ = 1.5 × 10⁻⁶: by the time that column is pivoted, earlier eliminations have combined the two diagonals' rows and the smallest pivot is 1.2 × 10⁻⁶, just above the threshold 10⁻⁶. `matrix_rank` uses σ_min against ε · σ_max · max(m, n) ≈ 2.2 × 10⁻¹⁶ · 2.4 · 12 ≈ 6 × 10⁻¹⁵, and σ_min = 3.6 × 10⁻⁷ is far above that. Both are *correct*: the matrix is not singular. What tracks the danger is σ_min (equivalently κ): it falls as y₄ and the forces grow as 1/σ_min ∝ 1/y₄. Exact rank is a yes/no question; σ_min answers "how close to no".

**(d)** The rank at tol = 10⁻³ drops to 11 between y₄ = 1.5 × 10⁻³ m and 1.5 × 10⁻⁴ m. Geometrically the pivot being discarded is the sine of the diagonal's angle, y₄/√(4 + y₄²) ≈ y₄/2: tol = 10⁻³ declares the diagonal "horizontal" when its slope is below about one part in a thousand, i.e. node 4 less than about 2 mm above the bottom chord. The threshold is a statement about your tolerance, not about the truss: at y₄ = 1.5 mm the structure is still a (terrible) structure.

**Common mistakes.** Reading κ as a force amplification (it bounds the *relative error* amplification; the force growth here comes from geometry, and happens to scale the same way); concluding from `matrix_rank = 12` that the design is fine.

## 5. Repair a mechanism (s = 15)

**Assumptions.** s mod 9 = 6 → remove member index 6, the vertical 2-5. Remaining 8 members, 3 reactions, 6 joints: Maxwell count 8 + 3 − 12 = −1.

**(a)** `classify_truss` → rank 11, k = 1, s = 0, status `mechanism`. The normalised mode (`truss_mechanisms(A)[:, 0]`, divided by its largest entry):

| node | u_x | u_y |
|---|---|---|
| 0 | 0 | 0 |
| 1 | 0 | 0.5 |
| 2 | 0 | 0 |
| 3 | −0.375 | 0 |
| 4 | −0.375 | 0.5 |
| 5 | −0.375 | 1 |

The whole top chord (3, 4, 5) translates left by 0.375 while rotating: node 5, no longer tied to node 2, rises by 1; node 4 rises by 0.5 and node 1 with it (vertical 1-4 keeps their vertical separation); node 3 stays at height 0 (the vertical 0-3 is pinned at node 0). It is the right bay folding: diagonal 1-5 and chord 4-5 swing about node 4, with the left bay rotating in sympathy about the support at node 0.

**(b)** `member_elongations(nodes, members, d)` returns eight zeros (max |e| = 0.0 to machine precision) and the supported components d[0], d[1], d[5] are all 0.0. Both residuals vanish because the mode is, by construction, in N(Aᵀ).

**(c)** Seven candidate members are not yet present. `classify_truss` on each:

| candidate | length (m) | status |
|---|---|---|
| 2-3 | 4.27 | determinate |
| 2-4 | 2.50 | determinate |
| 2-5 | 1.50 | determinate (the member you removed) |
| 0-2 | 4.00 | mechanism and indeterminate (k = 1, s = 1) |
| 0-5 | 4.27 | mechanism and indeterminate |
| 1-3 | 2.50 | mechanism and indeterminate |
| 3-5 | 4.00 | mechanism and indeterminate |

Pattern: a repair must stretch under the mode, i.e. join two joints whose *relative* displacement has a component along the member. Node 2 does not move; nodes 3, 4, 5 all move. 2-3: node 3 moves (−0.375, 0), member direction (−4, 1.5)/4.27 → elongation ≠ 0 ✓. 2-4 and 2-5 likewise ✓. The failures: 0-2 joins two fixed joints; 3-5 joins two joints that move by the same horizontal amount, so a horizontal member sees no change in length; 0-5 is the direction (4, 1.5) while node 5 moves (−0.375, 1) — the dot product is −1.5 + 1.5 = 0, the member is tangent to the motion; 1-3: node 1 moves (0, 0.5), node 3 moves (−0.375, 0), relative (−0.375, −0.5), member direction (−2, 1.5)/2.5 → 0.75 − 0.75 = 0, again tangent. Each failed candidate adds a column in the column space of the existing A (hence s rises by 1) without raising the rank (k stays 1).

**(d)** Among 2-3, 2-4 and 2-5: 2-5 is shortest (1.5 m), restores the symmetric layout, and under the design load carries 0 (it carried 0 before removal) — the obvious choice. 2-4 (2.5 m) works but is a second diagonal in the right bay crossing 1-5, inconvenient to fabricate and in compression if the load moves to node 2. 2-3 (4.27 m) spans the whole top, is the longest and would be the first to buckle in compression (Euler load ∝ 1/L²: a quarter of that of a 2.1 m member for the same section). Choose 2-5.

**Common mistakes.** Testing candidates by Maxwell's count alone (every candidate gives m + r − 2j = 0, and four of them are still mechanisms); describing the mode as "node 5 moves up" without noticing the whole top chord shifts.

## 6. What rank can and cannot tell you

**(a)** Not justified. k = 0 is a *first-order* statement about an idealised pin-jointed geometry. It says nothing about distance to singularity (compute σ_min and κ of A — in P4 a k = 0 truss had κ = 6 × 10⁶), about member stresses (solve for the forces and compare with yield), or about buckling (Euler load π²EI/L² of each compression member). "Safe" needs all four numbers; rank provides none of them.

**(b)** Not justified either. k = 1 means a displacement mode with zero *first-order* elongation. Whether it collapses depends on second-order behaviour: a pre-tensioned cable net (s > 0 and k > 0, like the Maxwell trap with the diagonals in tension) can be stiff — it is a tensegrity-like "infinitesimal mechanism". The measurable quantity is the second-order stiffness along the mode (the geometric stiffness from pre-stress, or a test: load along the mode and measure the displacement). Also check whether the mode is actually excited by any realistic load (a mode purely horizontal under purely vertical loads).

**(c)** Not justified. `matrix_rank` tells you the solve will not divide by zero. Reliability of the forces is a conditioning question: forward error ≲ κ · (ε_machine + relative error in the loads and geometry). With κ = 6 × 10⁴ and a 1 % geometric tolerance on y₄, the forces are uncertain by a factor of order 1 (P4: forces ∝ 1/y₄). Compute κ and propagate the input uncertainties, or compare the residual with the forces' sensitivity.

**(d)** Not justified. s = 2 means two independent self-stress states, i.e. two redundant load paths. Redundancy is not over-design: it is what lets a structure survive one member failure (the "fail-safe" requirement in aircraft structures) and lets pre-stress stiffen it. Removing two members may also create mechanisms (k can rise as s falls: Maxwell's count is a difference). The measurable quantities are the member utilisations (force/allowable) under every load case and the rank *after* each hypothetical removal. Remove a member only if the remaining structure stays k = 0 under every load case with acceptable utilisations.

## 7. A parity check for a faulty sensor

**Assumptions.** H = [[1, 0.5, 0], [0, 1, 0.5], [0.5, 0, 1], [1, 1, 1]]; **f** = (N, A, M); **y** = H **f** + **e**.

**(a)** `rank(H)` = 3 (the top 3 × 3 block has determinant 1 + 0.125 = 1.125 ≠ 0). `left_nullspace(H)` returns **y**⊥ ∝ (−2/3, −2/3, −2/3, 1): check **y**⊥ᵀH = (−2/3 − 1/3 + 1, −1/3 − 2/3 + 1, −1/3 − 2/3 + 1) = 0 — row 4 is 2/3 of the sum of rows 1–3. Then **y**⊥ᵀ**y** = **y**⊥ᵀH **f** + **y**⊥ᵀ**e** = 0 · **f** + **y**⊥ᵀ**e**: the loads drop out identically, so the residual depends only on the errors. For a healthy balance (**e** = 0) it is zero whatever the loads — a parity relation.

**(b)** **f** = (2, −1, 3) → **y** = H**f** = (1.5, 0.5, 4, 4) and residual **y**⊥ᵀ**y** = −(1.5 + 0.5 + 4)·2/3 + 4 = −4 + 4 = 0. Fault +0.5 on gauge 2: residual = (−2/3)(0.5) = −0.333. Fault +0.5 on gauge 4: residual = (1)(0.5) = +0.5.

**(c)** Signatures **y**⊥ᵀ**e**_i = (−2/3, −2/3, −2/3, 1) for gauges 1–4. All four are non-zero, so every single-gauge fault is *detectable* (the residual leaves zero). Isolation requires telling the signatures apart up to scale: with a one-dimensional left null space every signature is a scalar, and any two non-zero scalars differ only by a factor, which is indistinguishable from a different fault size. Worse, gauges 1–3 share exactly the same signature, so even the sign gives no help; gauge 4 is distinguishable from the others only by sign, and only if the fault sign is known. One parity relation can never isolate among three (or more) candidates: you need the signatures to be vectors in a space of dimension ≥ 2.

**(d)** Append a fifth row **h**₅ = (0, 0, 1) (a gauge reading the moment alone): H₅ is 5 × 3, rank 3, left null space of dimension 2. `left_nullspace(H5)` gives (as rows of the 2 × 5 signature matrix)

| | g1 | g2 | g3 | g4 | g5 |
|---|---|---|---|---|---|
| relation 1 | −0.667 | −0.667 | −0.667 | 1 | 0 |
| relation 2 | 0.444 | −0.222 | −0.889 | 0 | 1 |

The five signatures are the columns. All ten 2 × 2 determinants of column pairs are non-zero (values ±0.444, 0.889, −0.667, 0.222, …, 1.0), so no two signatures are parallel: a fault's direction in the residual plane identifies the gauge, and its length gives the size. The choice matters: appending (0.5, 1, 0) instead gives signatures with det(g1, g2) = 0 — gauges 1 and 2 still not isolable — because that row is a combination that leaves rows 1 and 2 symmetric.

**Common mistakes.** Normalising **y**⊥ differently and reporting "different" residuals — only ratios are meaningful; claiming a fault is isolated because its residual is non-zero (that is detection).

## Check-yourself answers

1. For a 3 × 5 matrix of rank 3: column space dimension 3 (all of R³), row space 3, null space 5 − 3 = 2, left null space 3 − 3 = 0.
2. The special solutions have an identity block in the free rows: the j-th has a 1 in free column f_j and 0 in the others, so a combination summing to zero must have every coefficient zero.
3. Row operations preserve the row space and the null space but change the column space: read the pivot *column indices* from R, then take those columns from the original A, never from R.
4. Maxwell's count gives s − k = m + r − 2j: a count of zero only says s = k. The "maxwell-trap" variant (chords, verticals, diagonals 0-4 and 1-3) has m + r − 2j = 0 but rank 11, so k = s = 1: a mechanism *and* a self-stress state.
5. Check that (i) `member_elongations` is zero for every member and (ii) the displacement component at every supported degree of freedom is zero. Both follow from Aᵀ**d** = 0; the elongation rows are the member columns of A, the support rows are the reaction columns.
6. A diagonal of a 2 m × 1.5 m bay has length 2.5 m and unit vector (0.8, 0.6): its column holds ±0.8 in the x-rows and ±0.6 in the y-rows of its two nodes.
7. Not safe to conclude anything from full rank alone: quote σ_min and κ = σ_max/σ_min (the truss with y₄ = 1.5 µm has full rank with κ = 6 × 10⁶ and forces of 1.25 × 10⁷ kN from a 15 kN load).
8. With axes (1, 0, 0), (0, 1, 0), (1, 1, 0)/√2 no torque about z can be produced; that direction is the left null space of B, N(Bᵀ) = span{(0, 0, 1)}.
