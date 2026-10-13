# P01.13  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) Augmented array; clear column 1 with multipliers −3/2 and −1 (row 2 ← row 2 + 1.5 row 1, row 3 ← row 3 + row 1); clear column 2 with the multiplier (entry)/(pivot 2); back-substitute from z. Pivots are the diagonal of the final array. (b) Swap rows 1 and 2 first (largest |entry| in column 1), clear, then check whether column 2 needs a swap; the row order is the permutation; product of pivots times (−1)^swaps. (c) L is unit lower with the multipliers of (a); U is the final array without b. `lu_decompose` returns the *pivoted* factors; verify `P @ A ≈ L @ U`.

## 2.
(a) Step k: n − k rows, each n − k + 1 multiply-subtract pairs (columns k..n plus b) = 2(n − k)(n − k + 1) flops. Put j = n − k, Σ_{j=0}^{n−1} 2j(j + 1) = 2[(n−1)n(2n−1)/6 + (n−1)n/2] = 2n³/3 + O(n²). Back substitution: Σ_i 2(n − i) ≈ n². Gauss–Jordan on [A | I]: eliminate above and below, n rows × 2n columns per step, ≈ 2n³, i.e. three times 2n³/3. (b) (I − ℓ **e**_i **e**_kᵀ)(I + ℓ **e**_i **e**_kᵀ) = I − ℓ² **e**_i (**e**_kᵀ **e**_i) **e**_kᵀ = I because **e**_kᵀ **e**_i = 0 for i ≠ k. The product of two such inverses with different (i, k) in the elimination order has no cross term for the same reason, so entries just accumulate below the diagonal. (c) Row operations of the form "add a multiple of a row" leave det unchanged; a swap multiplies by −1; det of a triangular matrix is the product of its diagonal; therefore det A = (−1)^swaps Π pivots.

## 3.
(a) Node 1: G₁₂(T₂ − T₁) + Q₁ = 0; node 2: G₁₂(T₁ − T₂) + G₂₃(T₃ − T₂) = 0; node 3: G₂₃(T₂ − T₃) + G₃ₛ(0 − T₃) = 0. Move to G **T** = **Q** with positive diagonals. (b) Eliminate column 1 with multiplier −1, then column 2; back-substitute; check G₃ₛ T₃ = Q₁. (c) Rebuild G with G₂₃ = 3 and solve; node 3's temperature is fixed by Q₁/G₃ₛ regardless. (d) G⁻¹[i, j] = ∂T_i/∂Q_j: the temperature rise at node i per watt injected at node j. Positivity: G is an M-matrix (diagonally dominant with non-positive off-diagonals), whose inverse is entrywise non-negative.

## 4.
(a) `min(timeit.repeat(lambda: f(A, b), number=1, repeat=3))` for each n; `np.polyfit(np.log(ns), np.log(ts), 1)[0]`; flop/s = `elimination_flops(n) / t`. (b) Loop over ε; `exact = np.linalg.solve`; errors in x₁; `iterative_refinement(T, t, steps=k, pivot=False)`. (c) Loop over n; `H @ np.ones(n)` as the right-hand side; `residual_norm`; `np.max(np.abs(x - 1))`.

## 5.
(a) ΣF_x: R0x = 0. ΣM about node 0: R2y × 4 = 15 × 2 + s × 2. ΣF_y: R0y = 15 + s − R2y. (b) `solve_truss(nodes, members, supports, {1: (0, -15000), 4: (0, -s*1000)})`. Zero-force: node 3 (two members, no load); node 2's x equation once the roller and vertical are known. (c) Loads `{1: (3000, -15000)}`; the pin's x reaction equals −3 kN; only the members in the load path from node 1 to node 0 along the bottom chord change. (d) `members[:-1] + [(3, 5)]`: the right bay 1-2-5-4 has no diagonal; elimination finds an exact zero pivot.

## 6.
For each: identify whether the evidence concerns the *solve* (residual) or the *problem* (condition number, physical plausibility), then name the cause (near-mechanism geometry; wrong assembly or unconverged solver; ill-conditioning; a sign error).

## 7.
(a) Tension in segment i (counted from the top, i = 1..n) is the weight below it, (n − i + 1) m g; extension of segment i is that over k; tip displacement is the sum, m g n(n+1)/(2k); top tension n m g. (b) Forward sweep: c′_i = c_i/(b_i − a_i c′_{i−1}), d′_i = (d_i − a_i d′_{i−1})/(b_i − a_i c′_{i−1}); back sweep x_i = d′_i − c′_i x_{i+1}. About 8n flops. (c) `time.perf_counter` around each solver on the full dense K (build K with `np.diag`); compare tip displacements. (d) |2k| > |−k| + |−k| is only equality, but after the first elimination step the diagonal entries become strictly larger than the remaining off-diagonal, and the largest |entry| in each column is always on the diagonal, so `argmax` returns the current row.
