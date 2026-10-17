# P01.18  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. For a 2×2 symmetric matrix the eigenvalues are (tr ± √(tr² − 4 det))/2; the singular values are their square roots. The inverse of an upper-triangular 2×2 is easy: [[a, b], [0, d]]⁻¹ = [[1/a, −b/(ad)], [0, 1/d]].
2. Every inequality in (a)–(c) comes from either "a maximum is at least any particular value" or "orthogonal matrices do not change the 2-norm or the Frobenius norm". For (d), write δ**u** = K⁻¹δ**f** and **f** = K**u** and bound each with consistency.
3. The pseudoinverse of a wide matrix gives the shortest solution, which lives in the row space; the true load has a component outside it that no two gauges can see. For (d), think about which singular direction of C is weakest and add a row along it.
4. Everything is a loop over r with `chain_stiffness_matrix`, `svd_jacobi`, `condition_number` and `error_amplification`. For (d), D⁻¹KD⁻¹ has unit diagonal; compute its κ with `condition_number`.
5. Separate three things: the bound (a theorem), the typical case (depends on the load direction), and the rounding loss (log₁₀ κ digits). Each sentence in the problem mixes two of them.
6. The noise floor of a 50 × 30 Gaussian matrix with entry standard deviation σ is roughly σ(√50 + √30) for its largest singular value. The two built-in terms are each rank one: one depends on α, the other does not.
