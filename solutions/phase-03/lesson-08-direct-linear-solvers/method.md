# P03.08  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. n = 100 s; LU 2n³/3; substitutions 2n²; Cholesky n³/3; inverse 2n³; Thomas 8n. 2. t = flops/10⁹; load cases in that time = (2n³/3)/(2n²) = n/3. 3. Inverse plan: 2n³ + m·2n² against 2n³/3 + m·2n²; the inverse costs three factorisations and its product A⁻¹b is not backward stable.

## 2.
(a) A_k = L_k U_k (leading blocks of a triangular product are products of leading blocks); det A_k = Π_{i≤k} u_ii; divide consecutive determinants. For [[0, 1], [1, 0]], det A₁ = 0 so u₁₁ = 0 and the second pivot is 0/0. Positive leading minors give positive ratios. (b) Partition A = [[A₁₁, A₁₂], [A₂₁, A₂₂]] after j − 1 columns; the j-th pivot is the (1,1) entry of S = A₂₂ − A₂₁A₁₁⁻¹A₁₂; with x = (−A₁₁⁻¹A₁₂ e₁, e₁), xᵀAx = e₁ᵀSe₁ = pivot ≤ 0. (c) x̂ − x = −A⁻¹r; ‖x̂ − x‖ ≤ ‖A⁻¹‖‖r‖ = κ‖r‖/‖A‖; divide by ‖x‖. For n = 4: κ₁ × η with η ≈ 10⁻¹⁶ (relative to ‖A‖‖x‖).

## 3.
(a) G_ii = sum of all conductances touching node i (radiator ones included); G_ij = −G_ij for unknown neighbours; rhs_i = Q_i + G_ir T_r. Rows 2 and 4 are strictly dominant; without the radiator every row sum is zero and G·(1, 1, 1, 1) = 0. (b) Column 1: L₁₁ = √G₁₁, L_i1 = G_i1/L₁₁; column 2: L₂₂ = √(G₂₂ − L₂₁²), …; then forward and back substitution. (c) Radiator heat = G₂r(T₂ − T_r) + G₄r(T₄ − T_r); compare with 35 W. (d) G is a stiffness (conductance) matrix: symmetric, SPD once a temperature is fixed; sparse Cholesky or CG for 5 000 nodes.

## 4.
(a) N_upper = 4 kN × 2 m/d; N_diag = −4 kN × √(0.25 + d²)/d. (b) `solve_truss(*spar_truss(depth=d), 70e9, 2e-4, [0, 1, 10, 11], f)`, `member_forces`, `np.linalg.eigvalsh(Kf)`. (c) Ratio of deflections against (0.3/d)²; the shear (diagonal and vertical) contribution scales differently. (d) Loop over d; count eigenvalues below 10⁻⁶ λ_max of the unsupported K.

## 5.
(a) `AtA = V.T @ V`, `c = cholesky_solve(cholesky(AtA), V.T @ y)`; `Q, R = np.linalg.qr(V)`, `c = back_substitution(R, Q.T @ y)`; `np.linalg.lstsq`. (b) `np.linalg.cond`; errors against κε. (c) Loop p; catch `ValueError` from `cholesky`. (d) Squaring the condition number: never form AᵀA for an ill-conditioned A; use QR or SVD.

## 6.
(a) Cholesky once, 200 substitutions. (b) Thomas (8n flops, O(n) memory). (c) Sparse-direct (nested dissection) or CG with a preconditioner (`P03.09`); dense is 10¹² entries. (d) Dense LU (or Cholesky if symmetric): 6×6 is 144 flops, overhead dominates. (e) Sparse-direct Cholesky (one factorisation, bandwidth-reduced) or PCG. (f) QR of the 10⁴ × 8 matrix (never the normal equations unless the columns are well conditioned).
