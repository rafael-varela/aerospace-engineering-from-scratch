# P01.03  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; exact results are marked, numerical ones are given to four significant figures and come from `solution.py`.

## 1. Your own section

**Assumptions.** NACA four-digit thickness (1), open trailing edge, solid symmetric section.

**Derivation** (s = 15: t = 0.08 + 15/400 = 0.1175, c = 1.75 m). (a) A = 0.685083 × 0.1175 × 1.75² = 0.246523 m² = 0.6851 t c². (b) x̄ = 0.42044 × 1.75 = 0.7358 m = 0.4204 c. (c) `section_area(0.1175, 1.75, n=1000)` = 0.246525 m², 2.07 × 10⁻⁶ m² above the exact value. (d) Predicted: 0.022 × (0.1175/0.12) × 1.75² × 1000^(−1.5) = 2.08 × 10⁻⁶ m². The law reproduces the error to 1 %, which confirms that the leading edge, not the smooth part, controls the error. (e) `section_area_substituted(0.1175, 1.75, n=1000)` is 1.4 × 10⁻⁷ m² off: fifteen times better at the same n.

**Result.** A = 0.2465 m² (0.6851 t c²), x̄ = 0.7358 m (0.4204 c); midpoint error 2.07 × 10⁻⁶ m² (predicted 2.08 × 10⁻⁶); substituted error 1.4 × 10⁻⁷ m².

**Check.** Dimensions: t c² is an area, c a length. Limiting case: t → 0 gives A → 0 with x̄ unchanged, as the closed forms show. The NACA 0012 at c = 1 reproduces 0.082210 m².

**Common mistakes.** Forgetting the factor 2 for both surfaces (A = 0.3425 t c²); using c instead of c² in A; scaling the error estimate by c instead of c².

## 2. The closed forms and the error law

**Derivation.** (a) A = 2∫₀ᶜ y_t dx. With ξ = x/c, dx = c dξ, the limits become 0 and 1: A = 2 · 5tc · c ∫₀¹ Σ aₖ ξ^(eₖ) dξ = 10tc² Σ aₖ ∫₀¹ ξ^(eₖ) dξ. By FTC part 2 with the antiderivative G(ξ) = ξ^(e+1)/(e + 1) (power rule, valid for e = ½ as for integers), ∫₀¹ ξ^e dξ = G(1) − G(0) = 1/(e + 1). Hence A = 10tc² Σ aₖ/(eₖ + 1) = 10tc² × 0.0685083 (exact: the sum of five rationals). For x̄: ∫₀ᶜ x · 2y_t dx = 10tc³ Σ aₖ ∫₀¹ ξ^(eₖ+1) dξ = 10tc³ Σ aₖ/(eₖ + 2); divide by A: x̄ = c Σ aₖ/(eₖ + 2)/Σ aₖ/(eₖ + 1) = c × 0.0288033/0.0685083 = 0.42044 c.

(b) Linear: on one panel ∫(αx + β) dx = αh x_m + βh (since ∫(x − x_m) dx = 0 over a symmetric interval), which equals h f(x_m). Error: f(x) = f(x_m) + f′(x_m)(x − x_m) + ½f″(x_m)(x − x_m)² + …; integrating over [x_m − h/2, x_m + h/2] gives h f(x_m) + 0 + ½f″(x_m) · (2/3)(h/2)³ + … = h f(x_m) + f″(x_m)h³/24 + O(h⁵). Summing: ∫ₐᵇ f − M_n = (h²/24) Σ f″(x_m) h + … ≈ (h²/24) ∫ₐᵇ f″ dx = (h²/24)[f′(b) − f′(a)], using the midpoint sum of f″ as an approximation to its integral (accurate to O(h²) itself). Assumed: f has a continuous second derivative on [a, b] and a finite f′ at both ends.

(c) x = cu², dx = 2cu du, (x/c)^e = u^(2e): A = 2∫₀¹ 5tc Σ aₖ u^(2eₖ) · 2cu du = 20tc² ∫₀¹ Σ aₖ u^(2eₖ+1) du. For t = 0.12, c = 1 the prefactor is 2.4 and g(u) = 2.4 Σ aₖ u^(2eₖ+1) = 0.71256u² − 0.3024u³ − 0.84384u⁵ + 0.68232u⁷ − 0.2436u⁹. g′(u) = 1.42512u − 0.9072u² − 4.2192u⁴ + 4.77624u⁶ − 2.1924u⁸, so g′(0) = 0 and g′(1) = −1.11744. Predicted error at n = 1280: (1/1280²)/24 × 1.11744 = 2.84 × 10⁻⁸ in magnitude; since ∫g − M_n = (h²/24)(g′(1) − g′(0)) < 0, the midpoint sum lies *above* the integral. The code measures 2.84 × 10⁻⁸.

**Check.** Integrating g(u) term by term gives 2.4 × (0.2969/3 − 0.1260/4 − 0.3516/6 + 0.2843/8 − 0.1015/10) = 2.4 × 0.0342542 = 0.082210: the same area, as a substitution must give. ✓

**Common mistakes.** Dropping the 2cu Jacobian; using ∫₀¹ ξ^(1/2) dξ = 1/2 instead of 2/3; claiming the midpoint rule is exact for quadratics (it is not: f″ ≠ 0).

## 3. Total impulse of a model-rocket motor

**Assumptions.** Thrust exactly piecewise linear between the five points; g = 9.80665 m/s²; vertical flight, no drag.

**Derivation** (s = 15: F_peak = 25 N). (a) Trapezoids: ½(0 + 25)(0.1) + ½(25 + 18)(0.2) + ½(18 + 12)(0.9) + ½(12 + 0)(0.3) = 1.25 + 4.30 + 13.50 + 1.80 = 20.85 N s (exact). Average thrust 20.85/1.5 = 13.90 N.

(b) All nine sums return 20.85 N s. With n = 15 the panel width is 0.1 s, so every breakpoint (0.1, 0.3, 1.2, 1.5 s) is a panel edge; on each panel F is a straight line, and the midpoint rule is exact for straight lines, so M₁₅ is exact. The left and right sums differ by h[F(1.5) − F(0)], which is zero because the thrust starts and ends at zero, and their average is the trapezoid rule, which is exact for piecewise-linear data whose breakpoints are panel edges; so L_n = R_n = T_n = exact. Refining to 150 and 1500 keeps the breakpoints on edges and changes nothing. (With n = 14 the breakpoints would fall inside panels and the sums would differ at the fourth figure.)

(c) Constant mass: Δv = 20.85/0.25 − 9.80665 × 1.5 = 83.40 − 14.71 = 68.69 m/s. With m(t) = 0.25 − 0.02t kg: ∫₀^1.5 F/m dt = 87.86 m/s by a 3000-panel midpoint sum, so Δv = 73.15 m/s. With the average mass 0.235 kg: 20.85/0.235 − 14.71 = 74.01 m/s. The constant-mass estimate is 6 % low, the average-mass estimate 1.2 % high; the integral is the right one because F/m is not linear in t.

(d) Quote 73 m/s with the mass-change integral. All three are wrong by more than 10 % if drag is included (at 70 m/s a 0.25 kg model rocket with C_D A ≈ 10⁻³ m² has D ≈ 3 N, a quarter of the average thrust) or if the published thrust curve is a nominal one (motor-to-motor impulse variation of ±5–10 % is normal).

**Check.** Units: N s / kg = m/s. Limiting case: m constant and F constant gives Δv = Ft/m − gt, the textbook result. The total impulse 20.85 N s falls in the 20–40 N s "E" class of model-rocket motors, so the number is realistic.

**Common mistakes.** Forgetting gravity during the burn; using n = 14 and reporting the sums as "slightly different because of numerical error" instead of asking why; dividing the impulse by the burn time to get a "speed".

## 4. Convergence and the Fundamental Theorem

**Result** (NACA 0012, c = 1; errors in m²).

(a) Orders on n ≥ 40: left 1.40, right 1.73, midpoint 1.51, substituted 2.00. On n ≥ 320 only: left 1.35, right 1.92 (the partial cancellation described in lesson §5 is passing through its best point), midpoint 1.50, substituted 2.00. The left and right orders are the ones that move, because their error is a sum of two power laws (h and h^(3/2)) whose balance shifts with n.

(b) Signed errors (sum − A): at n = 40, left −3.31 × 10⁻⁴, right −2.68 × 10⁻⁴; at n = 80, left −1.21 × 10⁻⁴, right −8.94 × 10⁻⁵. Both negative: the singular term dominates and has the same sign for both. Their difference, L_n − R_n = −h[f(c) − f(0)] = −0.00252h, is the first-order part: −6.3 × 10⁻⁵ at n = 40, −3.2 × 10⁻⁵ at n = 80, halving with h as it should. The average (L_n + R_n)/2 − A is −3.00 × 10⁻⁴, −1.05 × 10⁻⁴, … and converges with order 1.5, because the first-order terms cancel exactly and only the singular term is left; the average *is* the trapezoid rule, and it is still only order 1.5 on this integrand.

(c) FTC residuals |dF/dx − 2y_t| with n = 2000: at h = 10⁻³ they are 1.4 × 10⁻⁶ at 0.05 c, 1.9 × 10⁻⁷ at 0.5 c and 3.1 × 10⁻⁷ at 0.95 c; at h = 10⁻⁴ the leading-edge residual falls to 6.6 × 10⁻⁸ while the other two stay at 2.6 × 10⁻⁷ and 3.6 × 10⁻⁷. Near the leading edge the residual is central-difference truncation, (h²/6)f‴ with f‴ ∝ x^(−5/2) largest there, and it shrinks a hundredfold with h as it should; at mid-chord and aft the residual is already at the floor set by the 2000-panel sum's own h^(3/2) error, which the finite difference cannot remove however small h is. The error is largest at the leading edge for the coarse h and is limited by the integration, not the differentiation, elsewhere.

**Check.** The midpoint error ratio between n = 640 and 1280 is 1.35 × 10⁻⁶/4.77 × 10⁻⁷ = 2.83 = 2^(1.50). ✓

**Common mistakes.** Fitting the left and right orders through the whole range and reporting 1.4 as "the order of the left rule"; taking absolute values before averaging L and R; using h = 10⁻⁵ in (c) and blaming the Fundamental Theorem for the noise.

## 5. Centroid of a cambered section

**Derivation.** (a) ȳ = 0.015430 c (2412), 0.030859 c (4412), 0.046289 c (6412): exactly linear, ȳ/m = 0.7715.

(b) ȳ = ∫ y_c · 2y_t dx / ∫ 2y_t dx. The denominator is A, independent of m; the camber line is y_c = m × (shape depending on p only), so the numerator is m times an integral that does not depend on m. The ratio ȳ/m is a weighted mean of the shape function with weights 2y_t > 0, and a weighted mean of a function that rises from 0 to 1 (at x = pc) and falls back to 0 is strictly less than 1; here 0.77 because the weights favour the front of the section where the camber line is near its peak.

(c) The camber slope is largest at the leading edge: dy_c/dx = 2m/p = 0.1 for NACA 2412, θ_max = arctan(0.1) = 5.7°. Under the perpendicular construction the upper surface point moves back by y_t sin θ and the lower one forward by the same amount; to first order in θ the area gained on one surface is lost on the other, so the area change is second order in θ: of order (y_t θ)·(dθ/dx)·c integrated, roughly θ_max² × (a few per cent) ≈ 10⁻³ of A at most. The same cancellation applies to x̄: the upper surface moves aft and the lower forward by equal amounts, so their mean position is unchanged to first order. The first-order effect of the perpendicular construction is on the surface *shape* (the nose becomes blunter on the upper side), not on A or x̄.

(d) ★ Shoelace on 2000 points per surface: the perpendicular construction gives an area 0.087 % larger than the vertical-thickness value for NACA 2412, 0.35 % for 4412 and 0.78 % for 6412. The ratios 1 : 4 : 9 show it is second order in the camber, as argued; the sign is positive because the perpendicular offset lengthens the surface around the convex nose.

**Check.** m = 0 reproduces ȳ = 0 to 10⁻¹² (symmetric). The 4412 value is exactly twice the 2412 value in the code (ratio 2.0000).

**Common mistakes.** Reporting ȳ = m because "the centroid is on the camber line"; assuming the perpendicular construction changes the area at first order in θ.

## 6. When a Riemann sum misleads

(a) **x^(−1/2).** Integrable (value 2) but unbounded: the first panel's midpoint value √(2/h) is finite, so the sum converges, but only as h^(1/2) (order 0.5). Observable: the error falls by only √2 per doubling. Do: substitute x = u², which gives ∫₀¹ 2 du, exact.

(b) **1/x.** Not integrable: ∫ₕ¹ dx/x = −ln h → ∞. The midpoint sum grows like ln n without bound. Observable: the sum increases by a constant (about 0.69) each time n doubles and never settles. Do: nothing numerical can help; the integral does not exist, and the model that produced it needs a cutoff with a physical justification.

(c) **Step at 0.37.** Bounded and integrable (value 0.63), but the panel containing the jump contributes an error of up to h: order 1 regardless of rule. Observable: the error does not decrease smoothly; it jumps as the step moves across panel boundaries, and can be exactly zero by luck. Do: split the interval at the discontinuity and integrate each piece (exact for this f).

(d) **sin(200x).** Smooth, but the period is 0.0314; with n < 100 panels the midpoint sum samples fewer than four points per period and returns nonsense (aliasing: the error oscillates with n and can exceed the integral itself, 0.0026). Once h ≪ period the order is 2 with f′(1) − f′(0) = 200(cos 200 − 1). Observable: non-monotone error against n for small n. Do: use at least 10 panels per period, or integrate analytically.

(e) **11 noisy points.** The sum is limited by the data, not the rule: with ±1 % noise on each of 11 values the midpoint sum has about 1 %/√11 ≈ 0.3 % random error, and refining n is impossible (there is no f between the points). Observable: a different noise realisation changes the answer at the third figure. Do: fit a smooth model, integrate the fit, and quote the uncertainty from the noise, not from the rule.

**Common mistakes.** Treating (b) as a convergence problem; refining (c) blindly and finding the error has not improved; reporting the integral of (e) to five figures.

## 7. Second moment of area of the section

**Assumptions.** Solid symmetric section; bending about the chord line (which passes through the centroid by symmetry).

**Derivation.** (a) y_t³ = 125 t³ c³ (Σ aₖ ξ^(eₖ))³, a polynomial in ξ with half-integer powers from 1.5 to 12. SymPy expands and integrates it exactly: I = (2/3)∫₀ᶜ y_t³ dx = κ t³ c⁴ with κ = 203 861 397 623 761/5 173 168 000 000 000 = **0.039407**. For the NACA 0012 at c = 1 m: I = 0.039407 × 0.12³ = 6.810 × 10⁻⁵ m⁴.

(b) Midpoint errors for n = 40 … 1280: 6.1 × 10⁻⁹, 1.1 × 10⁻⁹, 1.9 × 10⁻¹⁰, 3.4 × 10⁻¹¹, 6.0 × 10⁻¹², 1.1 × 10⁻¹²: order **2.50**. Two facts combine. Near the leading edge f = (2/3)y_t³ ∝ x^(3/2), whose midpoint error on the first panel scales as h^(5/2), so the singular contribution is order 2.5. The smooth term (h²/24)[f′(c) − f′(0)] would be order 2 and would normally win, but f′(0) = 0 (the derivative of x^(3/2) vanishes at 0) and f′(c) = 2y_t(c)² y_t′(c) = 2 × (0.00126)² × (−0.140) = −4.5 × 10⁻⁷, so the second-order coefficient is 1.9 × 10⁻⁸ h², which at n = 40 is 1.2 × 10⁻¹¹, far below the measured 6 × 10⁻⁹. With the h² term almost switched off by the nearly closed trailing edge, the h^(5/2) term is all that is visible.

(c) √(I/A) = √(6.810 × 10⁻⁵/0.082210) = 0.02878 c. A rectangle c × tc has I = c(tc)³/12 = 1.440 × 10⁻⁴ m⁴ for t = 0.12, c = 1; the section's I is 0.473 of that, while its area is 0.685 of the rectangle's. The second moment weights material by the square of its distance from the chord line; the section's material is concentrated where y_t is less than its maximum (everywhere except x = 0.3c), and y³ punishes that far more than y does.

(d) I ∝ t³ c⁴. NACA 0012 → 0015: (0.15/0.12)³ = 1.953, mass × 1.25. Chord 0.25 → 0.30 m: (0.30/0.25)⁴ = 2.074, mass × 1.44. The thickness change gives almost the same stiffness gain for a 25 % mass increase against 44 %: thickening is cheaper in mass (which is why structural designers like thick sections and aerodynamicists push back).

**Check.** Dimensions: t³c⁴ has units m⁴ with t dimensionless. ✓ The rectangle comparison: I_section/I_rect must be less than A_section/A_rect because the section is "thinner on average", and 0.473 < 0.685. ✓

**Common mistakes.** Using (2y_t)³/12 without the factor 2/3 bookkeeping (I = ∫(2y_t)³/12 dx = ∫(8/12)y_t³ dx = ∫(2/3)y_t³ dx: the same thing, but easy to do twice); expecting order 1.5 because "it is the same square root"; comparing I values of different chords in different units.

## Check-yourself answers

1. Part 1: the accumulation function F(x) = ∫ₐˣ f has derivative F′(x) = f(x). Part 2: for any antiderivative G of f, ∫ₐᵇ f = G(b) − G(a).
2. On a panel the midpoint value times h equals the exact area under a straight line, because the linear term's deviation from the midpoint value integrates to zero by symmetry; the left value misses the (h/2) × slope × h triangle.
3. Error ≈ (h²/24)[f′(1) − f′(0)] = (0.0625/24) × 2 = 0.005208; and 1/3 − 0.328125 = 0.005208.
4. Because the integrand's derivative is infinite at the leading edge, the error law's f′(a) is undefined and the first panel alone contributes an error ∝ h^(3/2); the term responsible is 0.2969 √(x/c).
5. x = c u², which gives dx = 2cu du; the √(x/c) becomes u and every other term a power of u, so the integrand is a polynomial with the extra factor 2cu.
6. A ∝ t: NACA 0018 has 1.5 × 0.0822 = 0.123 c². A ∝ c²: doubling the chord gives 4 × 0.0822 c₀² = 0.329 c₀².
7. The section is thicker ahead of mid-chord (maximum at 0.30 c), so the first moment is weighted forward. With vertical thickness the thickness distribution, and therefore x̄, is unchanged by camber; only ȳ moves.
8. −96 m is the net change in height (the altimeter reads it); 165 m is the total vertical distance travelled, climbing and sinking counted positively (no instrument reads it directly; it is what the variometer's needle swept).
