# P03.08  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The counts are in 2.1–2.3; the time is flops/rate; the number of load cases in the same time is (factorisation flops)/(substitution flops). The second reason is accuracy: A⁻¹b has a larger backward error than two substitutions.
2. (a) det of a product of triangular matrices is the product of diagonals; the leading block of LU is the product of the leading blocks. (b) The Schur complement S = A₂₂ − A₂₁A₁₁⁻¹A₁₂ satisfies xᵀAx = (x₁ + A₁₁⁻¹A₁₂x₂)ᵀA₁₁(…) + x₂ᵀSx₂. (c) Divide A(x̂ − x) = −r by ‖x‖ after taking norms and inserting ‖A‖/‖A‖.
3. Row i of G has G_ii = Σ_j G_ij (including the radiator conductance) and G_ij = −G_ij; the radiator conductances appear on the diagonal but their off-diagonal partners are gone (the radiator is not an unknown), which is what makes G strictly diagonally dominant in those rows. Cholesky by hand: four square roots and six divisions.
4. The chord force is a moment divided by the depth; the deflection of a truss has a shear part (the diagonals) that does not scale as 1/d², which is why the pure bending scaling fails. κ₂ comes from `np.linalg.eigvalsh(Kf)`.
5. κ(VᵀV) = κ(V)² because the singular values of VᵀV are the squares of V's; the normal equations lose twice as many digits. When κ(V)² exceeds 1/ε, VᵀV is numerically singular and Cholesky meets a negative pivot.
6. Match the structure of the matrix (dense, SPD, banded, sparse, rectangular) and the number of right-hand sides to the cost table; for 10⁵–10⁶ unknowns, dense storage is out before cost is even considered.
