# P01.17  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. A = [[1, 0], [1, 4], [1, 8]], **b** = the three C_L values. 2. AᵀA = [[3, 12], [12, 80]], Aᵀ**b** = (ΣC_L, 4C_L,2 + 8C_L,3). 3. Solve by elimination. 4. **r** = **b** − A**x̂**; check Σr and Σαr. 5. `least_squares_qr(A, b)`; multiply c₁ by 180/π; compare with 6.28.

## 2.
(a) f(**x**) = **x**ᵀAᵀA**x** − 2**b**ᵀA**x** + **b**ᵀ**b**; ∇f = 2AᵀA**x** − 2Aᵀ**b**. Geometric: Aᵀ(**b** − A**x**) = 0 directly. (b) Transpose P using (AᵀA)⁻ᵀ = (AᵀA)⁻¹; P² collapses the middle (AᵀA)⁻¹AᵀA; PA = A(AᵀA)⁻¹(AᵀA) = A; (I − P)² = I − 2P + P² = I − P and (I − P)P = 0. (c) Substitute A = QR into AᵀA**x** = Aᵀ**b** and cancel Rᵀ. (d) cond₂(A) = σ_max/σ_min; eigenvalues of AᵀA are σ²; cond₂(AᵀA) = σ_max²/σ_min². For Läuchli: AᵀA = J + ε²I exactly; in floating point 1 + ε² = 1, so AᵀA = J, rank 1.

## 3.
(a) Generate C_L, C_D as specified; `least_squares_qr(design_matrix(cl), cd)` and the same with `"general"`; rms of each residual vector. (b) Complete the square; standard errors from `parameter_standard_errors(Ag, rg)`; first-order propagation: δC_L0 ≈ |∂/∂c₁|δc₁ + |∂/∂c₂|δc₂ with the partials of −c₁/(2c₂), or run 200 Monte Carlo trials. (c) `np.sign(residuals)` ordered by C_L; count the runs. (d) Substitute the general polar into the parabolic design: the parabolic fit absorbs the c₁C_L term into c₀ and c₂C_L² as best it can; the linear term with negative c₁ (C_L0 > 0) reduces the apparent curvature.

## 4.
(a) `fit_drag_polar`, `oswald_efficiency`, SE(e) = e·k_se/k; repeat for `wind_tunnel_data(n_points=...)`; `np.polyfit(log m, log SE, 1)`. (b) `monte_carlo_k_scatter(400, seed=s, noise_cl=0.0)` and `noise_cl=0.005`; ratio sd/se. (c) 2K C_L σ_CL at C_L = 0.8 and at C_L,max = 0.095 × 14; add in quadrature to 0.0008. (d) m_needed = 17 (SE_17/0.005)²; assumption: same α range, same noise, fixed-regressor formula.

## 5.
(a) Stalled points have C_D far above the parabola at high C_L; the fit raises K and lowers C_D0; residuals become negative at low C_L, positive at high; a runs test or a plot of residual against C_L shows the trend. (b) Equal variance breaks; weighted least squares with weights 1/σ_i², i.e. scale each row of A and **b** by 1/σ_i. (c) `condition_squaring(np.vander(alpha, 11, increasing=True))`; digits lost ≈ log₁₀ cond(AᵀA) for the normal equations, log₁₀ cond(A) for QR; rescale α to [−1, 1] and use QR or `lstsq`. (d) Anything about the *form* of the polar (its parabolic shape, the absence of a linear term) or about the noise model: the simulation assumed them.

## 6.
(a) Two calls to `classical_gram_schmidt`; check `np.max(np.abs(A - Q2 @ (R2 @ R1)))`. (b) Loop over ε with `laeuchli_matrix`; `orthogonality_loss` of each Q; `2.2e-16 * np.linalg.cond(A)`. (c) MGS: each projection coefficient is computed from an already-reduced vector, so the error is first order in ε_mach cond(A); CGS: the subtraction of the large common component happens after all coefficients are computed from the original column, which contains it, so errors compound to second order. (d) `back_substitution(R2 @ R1, Q2.T @ b)`; compare the error against x = ones with the single-pass MGS explicit route (2.0) and the augmented route (10⁻¹⁶); conclude which factor in the lesson 2.5 argument was responsible.
