# P00.14  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. Mean of the eight readings. Deviations from the mean, squared, summed, divided by n − 1, square-rooted: s. Divide by √n: u_A.
2. u_B(resolution) = 0.1 / √12. u_B(tolerance) = 0.05 / √3.
3. u_c = √(u_A² + u_B,res² + u_B,tol²). U = 2 u_c.
4. Round U to two significant figures, the mean to the same decimal place.
5. Offset = mean − 20.00 g; compare |offset| with u_c.

## 2.
1. (a) ∂y/∂a = 1, ∂y/∂b = ±1; square and sum.
2. (b) For y = ab: ∂y/∂a = b, ∂y/∂b = a; substitute into the law, divide by y² = a²b². For a/b the same with ∂y/∂b = −a/b².
3. (c) ∂y/∂a = p aᵖ⁻¹; divide by y = aᵖ.
4. (d) ρ ∝ p¹ T⁻¹: relative uncertainties of p and T add in quadrature (R exact). Multiply by ρ; double for k = 2; compare the two squared relative terms.

## 3.
1. Write C_D = D · (½ρ)⁻¹ · v⁻² · S⁻¹ and apply the product/power rule: relative uncertainties with exponents 1, 1, 2, 1.
2. Square each term; the largest square is the dominant input.
3. Multiply the relative u by C_D, then by k = 2; round.
4. Halve the velocity relative uncertainty, recompute the root-sum-square, compare.

## 4.
1. (a) Call `calibration_fit`; compute |slope − 1| / u_slope and |offset| / u_offset.
2. (b) `combine(type_b_resolution(0.1), fit.s_residual)`; errorbar plot of `fit.residuals` against reference; look for a bow, a step, or a trend.
3. (c) `corrected_value(23.4, 0.05, fit)`; report `.value`, `expanded(.u)`, and `.contributions` in order reading, slope, offset.
4. (d) `np.polyfit(reference, reading, 1, cov=True)` returns coefficients highest power first and a 2×2 covariance; square-root the diagonal. If they differ from yours by a constant factor, check the `cov` argument's documented scaling.

## 5.
1. List mechanisms: elastic non-linearity, creep or hysteresis, drift with temperature or battery, auto-zero behaviour, overload, friction in the mechanism.
2. For each, say what shape it leaves: a bow in residuals vs reference; loading/unloading curves that differ; a trend vs time or order; a jump.
3. Match the response: non-linearity → higher-order model or restricted range; hysteresis/drift → procedural fix (warm-up, same loading direction, re-zero); jump → find and remove the cause, re-take.

## 6.
1. (a)(i) Treat each reading's offset as its own input: four inputs, all independent; root-sum-square. (ii) Make the offset one input entering both readings: y = (r₁ + z) − (r₂ + z); the z terms cancel analytically. In code: inputs (r₁, r₂, z₁, z₂) with corr[2,3] = corr[3,2] = 1.
2. (b) Define uncertainty (a spread about a believed value) versus a blunder (a wrong value with no spread). Note that the thruster data's own uncertainty was small and internally consistent in either unit. The check that catches it is dimensional: carry units in code (pint) so lbf·s cannot be added to N·s.

## 7.
1. (a) The mean of n normal readings divided by its estimated standard error follows Student's t with n − 1 degrees of freedom; for ν = 2 the 97.5 % quantile is far above 2.
2. (b) Tabulate t₀.₉₇₅(ν) for ν = 2 … 20 and find where it is within 10 % of 2.
3. (c) `rng.normal(0.30, 0.10, 100_000)`; discard non-positive samples (and say how many); `t = sqrt(2 h / g)`; compare `t.std()` with the linear law's ½ · t · u_h / h; comment on the skew.
