# P01.18  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown.

## 1. Four norms and a condition number by hand

**Derivation** ($s = 15$: A = [[15, 1], [0, 2]]). (a) Column sums of |A|: 15 and 3, so ‖A‖₁ = 15. Row sums: 16 and 2, so ‖A‖∞ = 16. ‖A‖_F = √(225 + 1 + 4) = √230 = 15.166.

(b) AᵀA = [[225, 15], [15, 5]]: trace 230, determinant 1125 − 225 = 900. λ = (230 ± √(52 900 − 3600))/2 = (230 ± 222.04)/2 = 226.02 and 3.982. σ₁ = √226.02 = 15.034, σ₂ = √3.982 = 1.995. ‖A‖₂ = 15.034.

(c) κ₂ = 15.034/1.995 = 7.534. A⁻¹ = [[1/15, −1/30], [0, 1/2]]: ‖A⁻¹‖₁ = max(1/15, 1/30 + 1/2) = 0.5333, so κ₁ = 15 × 0.5333 = 8.0; ‖A⁻¹‖∞ = max(1/15 + 1/30, 1/2) = 0.5, so κ∞ = 16 × 0.5 = 8.0.

(d) `matrix_norm_2_power(A)` = 15.0339, `condition_number(A)` = 7.5339; `np.linalg.cond(A)` = 7.5339, `np.linalg.cond(A, 1)` = 8.0, `np.linalg.cond(A, np.inf)` = 8.0.

**Check.** ‖A‖₂ ≤ ‖A‖_F: 15.034 ≤ 15.166 ✓; ‖A‖_F² = σ₁² + σ₂² = 226.02 + 3.98 = 230 ✓; det A = 30 = σ₁σ₂ = 15.034 × 1.995 = 29.99 ✓. The three condition numbers differ by less than the factor n = 2 allowed.

**Common mistakes.** Taking the 2-norm as the largest entry or the largest eigenvalue of A (the eigenvalues of this triangular A are 15 and 2, and neither is a singular value); forgetting the square root.

## 2. The inequalities

**Derivation.** (a) ‖A‖ is the maximum of ‖A**x**‖/‖**x**‖ over **x** ≠ 0, so for any particular **x**, ‖A**x**‖/‖**x**‖ ≤ ‖A‖, i.e. ‖A**x**‖ ≤ ‖A‖‖**x**‖. Then ‖AB**x**‖ = ‖A(B**x**)‖ ≤ ‖A‖‖B**x**‖ ≤ ‖A‖‖B‖‖**x**‖; dividing by ‖**x**‖ and maximising gives ‖AB‖ ≤ ‖A‖‖B‖.

(b) ‖A**x**‖₁ = Σ_i |Σ_j a_ij x_j| ≤ Σ_i Σ_j |a_ij||x_j| = Σ_j |x_j| (Σ_i |a_ij|) ≤ (max_j Σ_i |a_ij|) Σ_j |x_j| = (max column sum) ‖**x**‖₁. With **x** = **e**_J for the maximising column J, A**x** is column J and ‖A**e**_J‖₁ = Σ_i |a_iJ| equals the bound, so the maximum is attained.

(c) With A = UΣVᵀ and **y** = Vᵀ**x** (‖**y**‖₂ = ‖**x**‖₂ because V is orthogonal): ‖A**x**‖₂ = ‖UΣ**y**‖₂ = ‖Σ**y**‖₂ (U has orthonormal columns), and ‖Σ**y**‖₂² = Σ σ_i² y_i² ≤ σ₁² Σ y_i² = σ₁²‖**y**‖₂², with equality at **y** = **e**₁, i.e. **x** = **v**₁. So ‖A‖₂ = σ₁. Frobenius: ‖A‖_F² = tr(AᵀA) = tr(VΣ²Vᵀ) = tr(Σ²VᵀV) = Σ σ_i². Then σ₁² ≤ Σ σ_i² ≤ n σ₁² gives ‖A‖₂ ≤ ‖A‖_F ≤ √n ‖A‖₂.

(d) δ**u** = K⁻¹δ**f** so ‖δ**u**‖ ≤ ‖K⁻¹‖‖δ**f**‖ = ‖δ**f**‖/σ_n; **f** = K**u** so ‖**f**‖ ≤ σ₁‖**u**‖, i.e. 1/‖**u**‖ ≤ σ₁/‖**f**‖. Multiply: ‖δ**u**‖/‖**u**‖ ≤ (σ₁/σ_n)‖δ**f**‖/‖**f**‖. Equality needs both inequalities tight: ‖K**u**‖ = σ₁‖**u**‖ requires **u** ∝ **v**₁, i.e. **f** ∝ **u**₁; ‖K⁻¹δ**f**‖ = ‖δ**f**‖/σ_n requires δ**f** ∝ **u**_n. The reverse pairing (**f** ∝ **u**_n so **u** = **f**/σ_n is as *long* as possible; δ**f** ∝ **u**₁ so δ**u** = δ**f**/σ₁ is as *short* as possible) gives ‖δ**u**‖/‖**u**‖ = (σ_n/σ₁)‖δ**f**‖/‖**f**‖: the amplification is 1/κ, the best case.

(e) σ₁ ≥ σ_n by ordering, so κ ≥ 1. An orthogonal Q has QᵀQ = I, eigenvalues all 1, so every σ_i = 1 and κ = 1. cA has singular values |c|σ_i, and the ratio is unchanged.

**Check.** The bar: `error_amplification(K, u_1, 0.01 u_n)` returns (1.751, 1.751) and the reverse pairing returns 0.01/175.09 = 5.71 × 10⁻⁵. ✓

**Common mistakes.** Proving (a) for the Frobenius norm by the same argument (it is not an induced norm; its submultiplicativity needs Cauchy–Schwarz); in (d), using **v**₁ where **u**₁ is needed: the load lives in the output space of K.

## 3. A strain-gauge balance with a dead gauge ($s = 15$)

**Assumptions.** Linear gauges, calibration matrix exact, independent Gaussian noise of 2 units.

**Result.** The simulated readings are **g** = (931.07, −341.06, 449.64, 443.74).

(a) `pseudoinverse(C) @ g` = (1001.17, −501.26, 299.17) N, identical to `np.linalg.pinv(C) @ g` to 10⁻¹⁰. κ₂(C) = 1.889 (singular values 1.614, 0.859, 0.854). Residual **r** = (0.24, 0.25, 0.36, −0.68): the fourth gauge disagrees most. With m − n = 1 degree of freedom, s² = **r**ᵀ**r**/1 = 0.711 and the standard errors √diag(s²(CᵀC)⁻¹) are (0.87, 0.87, 0.83) N. All three loads are within 1.5 standard errors of the truth; note that a single degree of freedom makes s² itself uncertain by a factor of about 2.

(b) Gauges 1–3: `np.linalg.solve(C[:3], g[:3])` = (1001.34, −501.09, 299.49) N, κ₂ = 1.522. The condition number *fell* (the fourth row, nearly parallel to the sum of the others, was the one that stretched the singular-value range), yet the recovery is no better and there is no longer any residual: with m = n every reading is reproduced exactly, so the noise goes straight into the loads with nothing averaged out and nothing left over to estimate σ from. Redundancy, not conditioning, is what a fourth gauge buys.

(c) Gauges 1–2 only: `pseudoinverse(C[:2]) @ g[:2]` = (1021.06, −438.66, −22.53) N. It solves the two equations exactly (residual (0, 0) to rounding) and has norm 1111.5 N against ‖**F**‖ = 1157.6 N for the truth: it is the shortest solution. The null-space direction is **z** = (0.06, 0.19, −0.98) (unit), along which the two gauges read nothing (C₂**z** = 0 to 10⁻¹⁶). The true load's component along it is **z**·**F** = −329 N, and that is exactly what is missing: projecting the truth onto the row space gives (1019.74, −437.50, −22.39), which is the pseudoinverse answer up to the 2-unit noise. Two gauges cannot see F_z, and the minimum-norm answer quietly sets the invisible part to zero.

(d) The weakest sensed direction is **v**₃ = (0.707, −0.707, 0)/‖·‖ with σ₃ = 0.854, but σ₂ = 0.859 is almost as small. Adding a row along **v**₃ lifts that direction (new singular values 1.614, 1.315, 0.859) and κ₂ falls only from 1.889 to 1.879, because the former σ₂ direction now sets the floor. Adding **v**₂ instead gives the same κ (1.889 → 1.889, the roles swap); adding **v**₁ makes things worse (κ = 2.222); a random unit row gives 2.17; the row (1, 0, 0) gives 2.03. The reasoned choice is any unit row in the plane of **v**₂ and **v**₃, and the honest conclusion is that one more gauge cannot lift *two* nearly degenerate weak directions: the balance needs two more, one along each.

**Check.** κ₂ with the **v**₃ row cannot be below σ₁/σ₂(old) = 1.614/0.859 = 1.879 whatever row is added, and that is what the code returns.

**Common mistakes.** Reading the pseudoinverse answer in (c) as "the loads"; expecting a square system to be better than an overdetermined one because it is "exactly solvable".

## 4. Your own damaged element ($s = 15$: element e = 6, index 5)

**Result.** (a)

| r | κ₂ | σ₁ [N/m] | σ_n [N/m] |
|---|---|---|---|
| 10⁻¹ | 3.13 × 10² | 2.580 × 10⁸ | 8.24 × 10⁵ |
| 10⁻² | 1.96 × 10³ | 2.578 × 10⁸ | 1.32 × 10⁵ |
| 10⁻³ | 1.85 × 10⁴ | 2.578 × 10⁸ | 1.39 × 10⁴ |
| 10⁻⁴ | 1.84 × 10⁵ | 2.578 × 10⁸ | 1.40 × 10³ |

Fitted exponent of κ against r over the last three rows: −0.99. σ₁ does not move (the stiff part of the bar sets it) while σ_n ∝ r: the soft element in series with the five nodes beyond it behaves as one spring of stiffness rk loaded by a rigid block, so σ_n ≈ 0.2 rk and κ ≈ 18.5/r. (With element 5 damaged, as in the lesson, the constant is 22.4; it depends on how many nodes sit beyond the break.)

(b) r = 10⁻³, 1 kN tip load, ten random 1 % perturbations (`np.random.default_rng(1815)`): actual relative displacement errors between 0.09 % and 1.14 %, median 0.5 %, against a bound of 185 % (1.85 × 10⁴ × 0.01). The scatter among draws is itself informative: the error depends on how much of the random perturbation happens to land on the free block beyond the break. Worst-case pair: 185.27 % actual against 185.27 % bound, equal to ten figures. The worst-case load **f**_w = (−0.33, 0.55, −0.60, 0.46, −0.17, 0, 0, 0, 0, 0) alternates on the *stiff* side of the break and the worst-case perturbation δ**f**_w = (0, 0, 0, 0, 0, 0.45, 0.45, 0.45, 0.45, 0.45) is a uniform push on the free block: the one load pattern a damaged bar amplifies is exactly the one that moves the mechanism.

(c) Tip deflection for r = 10⁻³: 14.414 mm, against 0.1429 mm intact; F(9/k + 1/(rk)) = 1000 × (9/7 × 10⁷ + 1/(7 × 10⁴)) = 0.1286 + 14.286 mm = 14.414 mm ✓. The damaged element carries 99.1 % of the elongation, 111 times the other nine together.

(d) Diagonal scaling: κ(D⁻¹KD⁻¹) = 275, 1.70 × 10³, 1.61 × 10⁴, 1.60 × 10⁵ for the four r values, against 313, 1.96 × 10³, 1.85 × 10⁴, 1.84 × 10⁵ unscaled: a 12–13 % reduction at every r. Almost all of the conditioning is geometry (a near-mechanism), not stiffness contrast; the scaling fixes the entries, and the soft mode is still there.

**Check.** The exponent −1 and the 1/r scaling of σ_n are consistent with the lesson's element-5 case; the tip-deflection identity holds to all printed figures.

**Common mistakes.** Fitting the exponent including r = 1 (the curve is not yet a power law there); reporting the bound as the expected error.

## 5. What a condition number is and is not

(a) Right: κ = 175 is modest, and the typical amplification of a smooth load error is close to 1 (1.36 % from 1 %). Wrong: "fine up to 10⁴" treats κ as a predictor when it is a bound, and ignores two separate facts. For a load with a component along the high mode **u**₁ the amplification *is* κ (the worst-case pair gave 175 %); and the rounding error of the solve itself is ε_mach κ ≈ 2 × 10⁻¹⁴ at κ = 175 and 2 × 10⁻¹² at κ = 10⁴, which is still negligible, so the *rounding* part of the sentence is right by accident. The honest statement: "for tip and distributed loads the displacement error is about the load error; κ tells us the worst load pattern could be amplified 175×, and we have checked that no test load resembles it."

(b) log₁₀(1.62 × 10⁶) = 6.2 digits of the 15.9 available, leaving about 9.7 correct digits in the worst case. A displacement to three figures is safe by six orders of magnitude. (Single precision, with 7.2 digits, would be marginal: one good digit left.)

(c) For r = 10⁻⁵ (element 5), σ_n/σ₁ = 4.5 × 10⁻⁷ < 10⁻⁴, so `pseudoinverse(K, rcond=1e-4)` sets 1/σ_n to zero: rank 9, and the solve drops any load component along the soft mode. Physically that declares element 5 *failed* and the six nodes beyond it a free body, which is false: the element is damaged but intact. But look at what the intact model predicts for it under the 1 kN tip load: an elongation of 1000/(10⁻⁵ × 7 × 10⁷) = 1.43 m in an element 0.1 m long. At r = 10⁻⁵ the model's own answer says linear elasticity has already failed, so the threshold question has become a modelling question, not a numerical one. Set `rcond` from physics (the smallest stiffness you believe can still be a linear spring), not from habit.

(d) Acoustic or buffet loading on a fuselage panel at a frequency whose half-wavelength matches the element size: the pressure alternates in sign from node to node, exactly the high-mode pattern. The mesh must then be fine enough that the alternating load is resolved (several elements per half-wavelength), the stiffness matrix is correspondingly larger and worse conditioned (κ ∝ n²), and the error amplification of the high mode is no longer unphysical. That is the case where iterative solvers need a preconditioner and the condition number is the first number to look at.

**Common mistakes.** Treating κ as an error estimate; choosing `rcond` to make a warning go away.

## 6. Modes of a pressure field ($s = 15$: noise 0.002)

**Result.** (a) Singular values of the 50 × 30 matrix: 52.447, 10.731, 0.0225, 0.0219, 0.0209, … Two are far above the rest; the expected noise floor for a 50 × 30 Gaussian matrix with entry standard deviation 0.002 is about 0.002(√50 + √30) = 0.025, and the third singular value, 0.0225, sits right at it. The noise-free matrix has rank exactly 2 (`np.linalg.matrix_rank` = 2).

(b)

| k | relative Frobenius error | 2-norm error | σ_{k+1} |
|---|---|---|---|
| 1 | 2.00 × 10⁻¹ | 10.731 | 10.731 |
| 2 | 1.38 × 10⁻³ | 0.02251 | 0.02251 |
| 3 | 1.31 × 10⁻³ | 0.02195 | 0.02195 |

The 2-norm error equals σ_{k+1} to all printed figures (Eckart–Young); the Frobenius error from rank 2 to rank 3 barely moves because everything beyond rank 2 is noise of nearly equal size in every direction.

(c) The matrix is **a**(α) ⊗ f₁(x) + **1** ⊗ f₂(x) with f₁ = √((1 − x)/x) (the flat-plate loading, singular at the leading edge) and f₂ = 2√(x(1 − x)) (a camber-like loading, symmetric about mid-chord). The SVD orthogonalises: **u**₁ is dominated by f₁ with an admixture of f₂ (the two functions are not orthogonal), and **u**₂ is the part of f₂ orthogonal to **u**₁. Correspondingly **v**₁ is nearly linear in α (it carries the 4α factor plus a constant offset) and **v**₂ is nearly constant with a linear correction. The SVD recovers the *subspace* of the two built-in terms exactly; the individual terms only up to a rotation within it.

(d) (σ₁² + σ₂²)/Σσ_i² = 0.999998: 99.9998 % of the energy in two modes. Storing **u**₁, **u**₂ (2 × 50 numbers) and the coefficients σ_i v_ij (2 × 30) costs 160 numbers instead of 1500, a factor of 9, with a relative error of 10⁻³ (the noise). It breaks when the physics stops being a sum of a few fixed shapes with varying amplitudes: separation that moves with α, a shock whose position changes, or any non-linear dependence spreads the energy over many singular values and the rank-2 picture fails. A singular-value plot that does not fall off a cliff is the diagnostic.

**Check.** ‖P‖_F² = 52.447² + 10.731² + (noise) = 2750.7 + 115.2 + 0.006 ≈ 2865.9 and the noise-free part has ‖·‖_F = 53.530, 53.530² = 2865.5 ✓.

**Common mistakes.** Expecting **u**₁ to equal f₁ exactly; counting the noise singular values as "modes".

## Check-yourself answers

1. ‖**x**‖ ≥ 0 with equality only at **0**; ‖c**x**‖ = |c|‖**x**‖; ‖**x** + **y**‖ ≤ ‖**x**‖ + ‖**y**‖. For (3, −4): 7, 5, 4.
2. ‖A‖ = max_{**x**≠0} ‖A**x**‖/‖**x**‖. ‖AB**x**‖ ≤ ‖A‖‖B**x**‖ ≤ ‖A‖‖B‖‖**x**‖; divide by ‖**x**‖ and take the maximum.
3. ‖A**x**‖₂² = **x**ᵀAᵀA**x** is a Rayleigh quotient of AᵀA, maximised by its largest eigenvalue σ₁²; ‖A‖_F² = tr(AᵀA) = sum of its eigenvalues = Σσ_i².
4. ‖δ**u**‖/‖**u**‖ ≤ κ₂(K)‖δ**f**‖/‖**f**‖ with κ₂ = σ₁/σ_n; attained for **f** ∝ **u**₁ and δ**f** ∝ **u**_n; the pairing **f** ∝ **u**_n, δ**f** ∝ **u**₁ gives the minimum amplification 1/κ.
5. κ₂ = 4.0 × 10⁴ (`condition_number` gives 40 002): up to log₁₀(4 × 10⁴) = 4.6 of the 16 digits can be lost.
6. κ₂ ≈ 16n²/π², quadratic in the element count. A finer mesh has stiffer local behaviour relative to its softest global mode, so the ratio of extreme stiffnesses grows; iterative solvers need more iterations and direct solvers lose more digits.
7. The set of displacement patterns that need no force: the nodes beyond the failed element translating together while those before it stay put. The right singular vector with σ = 0 (the last column of V) is that pattern.
8. Among all matrices of rank k it has the smallest error in both the 2-norm and the Frobenius norm (Eckart–Young); its 2-norm error is σ_{k+1}.
