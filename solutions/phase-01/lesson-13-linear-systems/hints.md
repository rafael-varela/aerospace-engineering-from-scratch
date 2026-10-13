# P01.13  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Keep the multipliers you use (−3/2, −1, then one more) as you go: they are the entries of L. With pivoting, the largest |entry| in column 1 is −3, so rows 1 and 2 swap before anything else. `lu_decompose` returns P as well; compare P A with L U, not A with L U.
2. For (a), substitute j = n − k and use Σ j² ≈ n³/3. For (b), check E E⁻¹ = I using **e**_kᵀ **e**_i = 0 for i ≠ k, and multiply two inverses to see that the multipliers do not interact. For (c), the determinant is unchanged by adding a multiple of one row to another, changes sign under a swap, and is the product of the diagonal for a triangular matrix.
3. Each diagonal entry is the sum of the conductances at that node; each off-diagonal is minus the conductance between the two nodes. The radiator's row has an extra G₃ₛ on the diagonal because space is a node held at the reference. Power leaving through the radiator is G₃ₛ T₃.
4. Use `time.perf_counter` and `min` over three repeats. The unpivoted error is the rounding of (1 − y) divided by ε. For Hilbert, compare the error with κ × 10⁻¹⁶.
5. Reactions first: the pin takes the horizontal load; moments about node 0 give the roller reaction. Zero-force members come from unloaded joints with only two non-collinear members (node 3) or three members two of which are collinear with no load in the third direction.
6. Residual tiny and forces huge: look at the condition number or at the determinant relative to the column norms. Residual large: the solver is not converged or the matrix is not what you think. Condition number 10¹⁴ means at most two trustworthy digits. A negative temperature rise with positive sources violates the maximum principle: a sign is wrong in G.
7. Equilibrium of the whole tether below node i gives the tension in segment i; the extension of each segment is tension/k; sum. The Thomas algorithm is elimination where each row has only one entry to clear and each back-substitution step has only one term.
