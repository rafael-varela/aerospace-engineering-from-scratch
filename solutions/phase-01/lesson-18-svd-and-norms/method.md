# P01.18  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. Column sums and row sums of |A|; Frobenius from the squares. 2. AᵀA = [[s², s], [s, 5]]; eigenvalues from trace s² + 5 and determinant 4s²; σ = √λ. 3. κ₂ = σ₁/σ₂; A⁻¹ = [[1/s, −1/(2s)], [0, 1/2]]; κ₁ = ‖A‖₁‖A⁻¹‖₁, κ∞ likewise. 4. `matrix_norm_*`, `condition_number`, `np.linalg.cond(A, p)`.

## 2.
(a) ‖A**x**‖/‖**x**‖ ≤ max over all **x**; then ‖AB**x**‖ ≤ ‖A‖‖B**x**‖ ≤ ‖A‖‖B‖‖**x**‖ and divide by ‖**x**‖. (b) Bound Σ_i|Σ_j a_ij x_j| by Σ_j|x_j|Σ_i|a_ij| ≤ (max column sum)‖**x**‖₁; attain with **e**_j. (c) ‖A**x**‖₂ = ‖ΣVᵀ**x**‖₂ with **y** = Vᵀ**x**, ‖**y**‖ = ‖**x**‖; ‖Σ**y**‖² = Σσ_i²y_i² ≤ σ₁²‖**y**‖², equality at **y** = **e**₁; ‖A‖_F² = tr(AᵀA) = tr(VΣ²Vᵀ) = Σσ_i²; then σ₁² ≤ Σσ_i² ≤ nσ₁². (d) ‖δ**u**‖ ≤ ‖K⁻¹‖‖δ**f**‖ and ‖**f**‖ ≤ ‖K‖‖**u**‖; multiply. Attainment: K⁻¹**u**_n = **v**_n/σ_n and K⁻¹**u**₁ = **v**₁/σ₁. Reverse pairing: ratio 1/κ. (e) σ₁ ≥ σ_n; Q has all σ = 1; cA has σ scaled by |c|.

## 3.
(a) `g = C @ F + 2 * rng.standard_normal(4)`; `pseudoinverse(C) @ g`; `condition_number(C)`; residual **g** − C**F̂**; s² = **r**ᵀ**r**/(4 − 3); standard errors √diag(s²(CᵀC)⁻¹). (b) `np.linalg.solve(C[:3], g[:3])`; κ₂ of the 3×3; the square system has no redundancy, so no averaging of noise and no residual to estimate the noise from. (c) `pseudoinverse(C[:2]) @ g[:2]`; check C₂**F̂** = **g**₂; `null_space(C[:2])`; project the true **F** onto the null vector; compare norms. (d) `svd_jacobi(C)`: take **v**₃ (the row direction C senses least); append it as a row; recompute κ₂; compare with appending **v**₁ or a random unit row.

## 4.
(a) Loop; `svd_jacobi` gives σ₁ and σ_n; `np.polyfit(log r, log κ, 1)`; σ₁ barely moves while σ_n ∝ r. (b) `error_amplification` for ten `rng.standard_normal` perturbations scaled to 1 %, then with `worst_case_loads`. (c) `np.linalg.solve(K, f)[-1]` against F(9/k + 1/(rk)); the soft element's share is (1/(rk))/(9/k + 1/(rk)). (d) `D = np.sqrt(np.diag(K))`; `Ks = K / np.outer(D, D)`; `condition_number(Ks)`; compare the r-dependence of κ(K) and κ(Kₛ).

## 5.
(a) κ bounds the worst case; the 1.36 % is the typical case for a smooth load; "fine up to 10⁴" confuses a bound with a prediction and ignores the digits lost to rounding. (b) log₁₀(1.6 × 10⁶) ≈ 6.2 digits of 16; compare 10 remaining digits with 3 required. (c) `rcond * σ₁` against σ_n ≈ 22.4⁻¹ r σ₁ ... compute σ_n/σ₁ for r = 10⁻⁵ and compare with 10⁻⁴; the threshold declares the element failed. (d) Think of a panel under acoustic loading or a checkerboard pressure pattern on a fine mesh; the mesh must resolve the pattern and the conditioning then matters.

## 6.
(a) Build the matrix; `svd_jacobi(P.T)`; count σ_i above σ_noise(√50 + √30). (b) `low_rank_approximation` for k = 1, 2, 3; `frobenius_norm` ratios; `matrix_norm_2_power(P.T - Ak)` against σ_{k+1}. (c) Plot **u**₁, **u**₂ against x and compare with √((1 − x)/x) and √(x(1 − x)) (they will be mixtures, because the SVD orthogonalises them); **v**_i against α. (d) Σ_{i≤2}σ_i²/Σσ_i²; storage 2 × 50 + 2 × 30 numbers against 1500; a flow with separation, shock motion or non-linear α dependence breaks the two-mode structure.
