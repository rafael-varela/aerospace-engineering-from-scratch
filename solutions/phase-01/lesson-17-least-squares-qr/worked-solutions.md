# P01.17  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown. All aerodynamic coefficients are dimensionless.

## 1. A lift curve by hand

**Derivation** ($s = 15$: C_L = 0.215, 0.585, 0.995 at α = 0°, 4°, 8°). A = [[1, 0], [1, 4], [1, 8]]. AᵀA = [[3, 12], [12, 80]]; Aᵀ**b** = (0.215 + 0.585 + 0.995, 4 × 0.585 + 8 × 0.995) = (1.795, 10.30). Eliminate: from the first row c₀ = (1.795 − 12c₁)/3; substitute into 12c₀ + 80c₁ = 10.30: 4(1.795 − 12c₁) + 80c₁ = 10.30 → 7.18 + 32c₁ = 10.30 → c₁ = 0.0975 per degree, c₀ = 0.2083.

(b) Fitted values 0.2083, 0.5983, 0.9883; residuals **r** = (0.0067, −0.0133, 0.0067). Σr = 0 ✓ (orthogonal to the ones column); Σαr = 0 × 0.0067 + 4 × (−0.0133) + 8 × 0.0067 = −0.0533 + 0.0533 = 0 ✓.

(c) `least_squares_qr(A, b)` returns (0.20833, 0.09750); `np.linalg.lstsq` the same. c₁ = 0.0975 × 180/π = 5.59 per radian, 11 % below 2π = 6.28, which is typical of a finite wing (the lift-curve slope falls with aspect ratio: 2π/(1 + 2/AR) = 5.03 for AR = 8, so this wing's slope sits between the finite-wing estimate and the 2D value).

**Result.** C_L = 0.2083 + 0.0975 α[deg]; residual rms 0.0094.

**Check.** The middle residual is minus twice each end residual: for three equally spaced points a straight-line fit always gives residuals in the ratio (1, −2, 1) times the deviation from linearity, which is a limiting-case identity worth remembering.

**Common mistakes.** Forgetting the column of ones (forces the line through the origin); mixing degrees and radians in the slope.

## 2. Projections, normal equations, conditioning

**Derivation.** (a) f(**x**) = (A**x** − **b**)ᵀ(A**x** − **b**) = **x**ᵀAᵀA**x** − 2**b**ᵀA**x** + **b**ᵀ**b**. With M = AᵀA symmetric, ∇(**x**ᵀM**x**) = 2M**x** and ∇(**b**ᵀA**x**) = Aᵀ**b**, so ∇f = 2AᵀA**x** − 2Aᵀ**b** = 0 gives AᵀA**x** = Aᵀ**b**. Geometrically: the residual **b** − A**x** must be orthogonal to the column space, i.e. to each column: Aᵀ(**b** − A**x**) = 0, the same equation in one line.

(b) Pᵀ = (Aᵀ)ᵀ((AᵀA)⁻¹)ᵀAᵀ = A(AᵀA)⁻¹Aᵀ = P, since (AᵀA)⁻¹ is symmetric. P² = A(AᵀA)⁻¹(AᵀA)(AᵀA)⁻¹Aᵀ = A(AᵀA)⁻¹Aᵀ = P. PA = A(AᵀA)⁻¹(AᵀA) = A. For I − P: (I − P)² = I − 2P + P² = I − P and (I − P)ᵀ = I − P, so it is a projector; (I − P)A = 0, so its range is orthogonal to the column space, and P + (I − P) = I splits every vector into the two parts.

(c) AᵀA = RᵀQᵀQR = RᵀR and Aᵀ**b** = RᵀQᵀ**b**. The normal equations read RᵀR**x** = RᵀQᵀ**b**; Rᵀ is invertible (R has non-zero diagonal), so R**x** = Qᵀ**b** and **x** = R⁻¹Qᵀ**b**.

(d) cond₂(A) = σ_max/σ_min. AᵀA is symmetric with eigenvalues σ_i², so cond₂(AᵀA) = λ_max/λ_min = σ_max²/σ_min² = cond₂(A)². Läuchli, ε = 10⁻⁸: AᵀA = J + ε²I exactly, where J is the all-ones matrix, with eigenvalues 3 + ε² (once) and ε² (twice): invertible, with cond = (3 + 10⁻¹⁶)/10⁻¹⁶ ≈ 3 × 10¹⁶. In double precision 1 + ε² = 1 + 10⁻¹⁶ rounds to 1 (the spacing of doubles near 1 is 2.2 × 10⁻¹⁶), so the stored AᵀA is exactly J, of rank 1: `numpy.linalg.solve` finds a zero pivot and raises `LinAlgError`. The information that distinguished the three columns was in ε, and squaring sent it below the rounding floor.

**Check.** `condition_squaring(design_matrix(cl))` on the polar returns (2.43, 5.92) = (2.43, 2.43²). ✓

**Common mistakes.** Writing ∇(**x**ᵀM**x**) = M**x** (loses the factor 2, harmless here but wrong); claiming P⁻¹ exists (a projector onto a proper subspace is singular).

## 3. A cambered wing ($s = 15$: C_L0 = 0.15)

**Assumptions.** Same α grid and noise as the lesson; the general polar is the true model.

**Result.** (a) Parabolic fit: C_D0 = 0.02102, K = 0.03592, residual rms 2.30 × 10⁻³. General fit: c₀ = 0.02304, c₁ = −0.01335, c₂ = 0.04627, residual rms 0.90 × 10⁻³. The general fit's rms is at the noise level (0.8 × 10⁻³ in C_D plus the C_L contribution); the parabolic fit's is 2.6 times larger.

(b) Completing the square: c₀ + c₁C_L + c₂C_L² = c₂(C_L + c₁/2c₂)² + c₀ − c₁²/(4c₂), so C_L0 = −c₁/(2c₂) = 0.1443, C_D,min = c₀ − c₁²/(4c₂) = 0.02207, K = c₂ = 0.04627. Standard errors of (c₀, c₁, c₂) from `parameter_standard_errors`: (0.00040, 0.00151, 0.00125). First-order propagation for C_L0: ∂/∂c₁ = −1/(2c₂) = −10.8, ∂/∂c₂ = c₁/(2c₂²) = −3.12, so δC_L0 ≈ √((10.8 × 0.00151)² + (3.12 × 0.00125)²) = 0.017 (ignoring the c₁–c₂ correlation, which tightens it somewhat). So C_L0 = 0.144 ± 0.017 (true 0.15), C_D,min = 0.0221 ± 0.0004 (true 0.022), K = 0.0463 ± 0.0012 (true 0.0468). All within one standard error.

(c) Parabolic residual signs in C_L order: + + + − + − − − − − − − − − + + +. Three positive at the bottom, a long negative run in the middle, three positive at the top: 4 runs in 17 points, where random noise would give about 9. The rms says the fit is poor; the sign pattern says *how*: the parabola is too flat at the ends and too high in the middle, the signature of a missing linear term.

(d) Biased low. The true C_D = (C_D,min + K C_L0²) − 2K C_L0 C_L + K C_L². The parabolic model has no C_L term, so the −2K C_L0 C_L piece (negative for C_L0 > 0) must be absorbed by the other two: over a range where C_L is mostly positive, a negative linear trend is mimicked by a *smaller* curvature and a smaller constant. Hence K_fit = 0.0359 < 0.0468 (23 % low) and C_D0,fit = 0.0210 < 0.0221. The Oswald efficiency from the parabolic fit would be 1.11, an impossible value above 1, which is itself a warning sign.

**Check.** The general fit reproduces the rms expected from the noise model; the parabolic fit's K bias has the sign predicted by the algebra.

**Common mistakes.** Reporting K from the two-parameter fit for a cambered wing; propagating errors without noticing that an e > 1 result is unphysical.

## 4. How many points?

**Result.** (a) Default data: e = 0.857, SE(e) = e × 0.00054/0.04642 = 0.0100. With 34, 68, 136 points: SE(e) = 0.0051, 0.0034, 0.0024. Fitting log SE against log m gives an exponent of −0.67 over this particular set of seeds; the expected −0.5 is recovered on average over seeds (each single realisation also carries its own s², which scatters by about 1/√(2(m − 2))). (b) Seed 15, 400 trials: C_L noise 0: sd = 0.000352, mean SE = 0.000337, ratio 1.04. C_L noise 0.005: sd = 0.000406, mean SE = 0.000361, ratio 1.12. (c) 2K C_L σ_CL = 2 × 0.0468 × 0.8 × 0.005 = 3.7 × 10⁻⁴ at C_L = 0.8 and 6.2 × 10⁻⁴ at C_L,max = 1.33. Adding the C_L = 0.8 value in quadrature to 8 × 10⁻⁴ gives 8.8 × 10⁻⁴; the seed-17 fit's rms of 1.19 × 10⁻³ is at the high end because that realisation happens to be noisy: the Monte Carlo mean rms is 0.79 × 10⁻³, below the total noise because the fit absorbs 2 of 17 degrees of freedom (E[rms²] = σ²(m − n)/m gives 0.94 × 8.5 × 10⁻⁴ ≈ 0.80 × 10⁻³). (d) SE ∝ 1/√m: m = 17 × (0.0100/0.005)² = 68 points, assuming the same α range, the same noise and the fixed-regressor formula; with the 12 % extra scatter from C_L noise, about 85.

**Check.** The 68-point run gave SE(e) = 0.0034, better than the 0.005 target: the single-seed s² happened to be low; across seeds the formula is right.

**Common mistakes.** Fitting the exponent to one seed and reporting −0.68 as "the" law; forgetting that the target standard error must be stated with a confidence level.

## 5. When the fit lies

(a) Stalled points at α = 14°–18° have C_D far above the parabola (drag can double past stall) at the highest C_L (which itself stops growing). The fit pivots: K rises to meet them, C_D0 falls to compensate at the low end, and the residuals show a trend: negative at low C_L, positive at high C_L, with perhaps 3–4 sign runs over 20 points where noise gives 10. A runs count or a plot of residual against C_L detects it without knowing the truth.

(b) Equal variance (homoscedasticity) breaks: the high-α points carry 2.5 times the noise. Ordinary least squares then over-weights them. Weighted least squares divides each row of A and **b** by its σ_i (so each equation has unit noise) and the standard-error formula becomes s² (AᵀWA)⁻¹ with W = diag(1/σ_i²). The weights change the fit by about one standard error here; they change the error bars more.

(c) `condition_squaring(np.vander(alpha, 11, increasing=True))` gives cond(A) = 5.8 × 10¹⁰ and cond(AᵀA) = 9.6 × 10²¹ (the second exceeds 10¹⁶, so AᵀA is singular to rounding). The normal equations can lose up to 22 digits: none of the 8 reported figures is guaranteed, and the result depends on the rounding of the solver. QR or `lstsq` loses up to 11 digits: about 5 trustworthy figures. Rescaling α to t = (α − 4)/8 ∈ [−1, 1] brings cond(A) to 3.9 × 10³ (cond(AᵀA) = 1.5 × 10⁷): QR then keeps about 12 figures and even the normal equations keep about 9, so all 8 reported figures are back.

(d) The parabolic *form* of the polar, the absence of a linear term, the Gaussian and independent noise, and the value of AR: all were put in. The simulation validates the *solver and the error bars given the model*; it cannot validate the model. Only tunnel data can.

**Common mistakes.** Treating a small rms as proof of the model; trusting coefficient digits beyond what the condition number allows.

## 6. Repairing Gram–Schmidt

**Result.** (a) `cgs2` returns Q₂ and R = R₂R₁ with max|A − Q₂R| = 1.1 × 10⁻¹⁶ at ε = 10⁻⁸ (and smaller at larger ε).

(b)

| ε | cond(A) | CGS | MGS | CGS2 | ε_mach cond(A) |
|---|---|---|---|---|---|
| 10⁻⁵ | 1.7 × 10⁵ | 4.8 × 10⁻⁸ | 5.9 × 10⁻¹³ | 1.1 × 10⁻¹⁶ | 3.8 × 10⁻¹¹ |
| 10⁻⁷ | 1.7 × 10⁷ | 1.3 × 10⁻² | 1.6 × 10⁻⁹ | 2.2 × 10⁻¹⁶ | 3.8 × 10⁻⁹ |
| 10⁻⁸ | 1.7 × 10⁸ | 5.0 × 10⁻¹ | 7.1 × 10⁻⁹ | 2.2 × 10⁻¹⁶ | 3.8 × 10⁻⁸ |

(c) MGS tracks ε_mach cond(A) (within a factor of 10): each coefficient is computed from a vector from which the earlier components have already been removed, so the error made in computing it is relative to what is *left*, and the bound is first order in the conditioning. CGS computes every coefficient from the original column, which still contains the large component along q₀; the rounding error in that large component, of relative size ε_mach, is then re-subtracted with the wrong coefficient and grows as cond(A)²: 10⁻⁸ at cond 10⁵, 10⁻² at cond 10⁷, order one at cond 10⁸, matching the squares. A second pass removes what the first left, because the leftover is now small and the second pass's error is ε_mach times a small quantity: "twice is enough".

(d) CGS2 with an explicit Qᵀ**b** at ε = 10⁻⁸ gives an error of 4.4 × 10⁻¹⁶ in **x**: the re-orthogonalisation rescues the explicit route completely, matching the augmented-column route (2 × 10⁻¹⁶) and Householder. This confirms the lesson 2.5 diagnosis: the single-pass MGS explicit route failed (error 2.0) because its Q had lost orthogonality by about ε_mach cond(A) ≈ 7 × 10⁻⁹, and that error, divided by R₃₃ ≈ ε, became order one. With Q orthogonal to 10⁻¹⁶ the small components of Qᵀ**b** are computed to relative precision (each is an accurately orthogonal q_i dotted with **b**), and dividing them by a small R₃₃ costs nothing. The augmented-column trick is the cheap fix when you only have one MGS pass; a second pass or Householder is the robust one.

**Check.** `np.linalg.qr(A)` (Householder) gives orthogonality loss 2.2 × 10⁻¹⁶ at all three ε, same as CGS2.

**Common mistakes.** Assuming that because CGS fails, CGS2 must also fail (it does not: the second pass works on a nearly orthogonal input and its error is ε_mach times a small quantity); conflating loss of orthogonality with loss of the factorisation A = QR, which all three methods keep to 10⁻¹⁶ or better.

## Check-yourself answers

1. AᵀA**x̂** = Aᵀ**b**. They say the residual **b** − A**x̂** is orthogonal to every column of A, i.e. A**x̂** is the projection of **b** onto the column space.
2. x = **a**ᵀ**b**/**a**ᵀ**a** = 3/2; **p** = (1.5, 1.5, 0); **e** = (−0.5, 0.5, 2); **e**·**a** = −0.5 + 0.5 = 0.
3. Pᵀ = P, P² = P, PA = A (also trace P = n, eigenvalues 0 or 1). Proof of P² = P: A(AᵀA)⁻¹AᵀA(AᵀA)⁻¹Aᵀ = A(AᵀA)⁻¹Aᵀ.
4. R is the upper-triangular matrix of Gram–Schmidt coefficients (R_ij = q_iᵀa_j). R**x** = Qᵀ**b** is solved by back substitution in n² operations with no AᵀA formed, so the condition number of the problem is cond(A), not cond(A)².
5. Classical computes each projection coefficient from the original column; modified computes it from the remainder after the earlier projections are subtracted. In floating point the modified version keeps Q orthogonal to about ε_mach cond(A), the classical to ε_mach cond(A)².
6. The eigenvalues of AᵀA are the squared singular values of A, so the ratio of extreme eigenvalues is the squared ratio of extreme singular values. Hopeless when cond(A)² ≳ 10¹⁶, i.e. cond(A) ≳ 10⁸.
7. d/dC_L [C_L/(C_D0 + K C_L²)] = 0 ⇒ C_D0 + K C_L² − 2K C_L² = 0 ⇒ C_L* = √(C_D0/K); then L/D = C_L*/(2C_D0) = 1/(2√(C_D0 K)).
8. e = 1/(π × 8 × 0.04642) = 0.857; δe = e × 0.00054/0.04642 = 0.010; e = 0.857 ± 0.010.
