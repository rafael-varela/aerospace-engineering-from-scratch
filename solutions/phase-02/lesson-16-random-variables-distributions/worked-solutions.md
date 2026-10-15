# P02.16  Worked solutions

Each solution states: assumptions, units, the derivation, the result (exact, numerical with tolerance, or simulated with standard error), at least one independent check, and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py`. Every simulated dataset is synthetic.

## 1. Two discrete variables by hand ($s = 15$: p = 0.915, λ = 2.5 per month)

**Assumptions.** Thrusters independent and identical; alerts a Poisson process at constant rate.

**Derivation.** (a) pmf C(4, k)(0.915)ᵏ(0.085)⁴⁻ᵏ: k = 0..4 → 5.22 × 10⁻⁵, 0.002248, 0.03629, 0.26046, 0.70095 (sum 1.0000). E[X] = Σ k P(k) = **3.660** = 4p ✓. Var X = Σ (k − 3.66)² P(k) = **0.3111** = 4p(1 − p) = 4 × 0.915 × 0.085 ✓. P(X ≥ 3) = 0.26046 + 0.70095 = **0.9614**.
(b) P(0) = e^{−2.5} = **0.0821**; P(≥ 3) = 1 − e^{−2.5}(1 + 2.5 + 3.125) = **0.4562**; a quarter has rate 7.5, so P(0 in a quarter) = e^{−7.5} = **5.53 × 10⁻⁴**, which equals (e^{−2.5})³: the three months are independent and the multiplication rule gives the same answer.
(c) `binomial_pmf(np.arange(5), 4, 0.915)`, `discrete_mean`, `discrete_variance` and `poisson_pmf` reproduce every figure.

**Check.** Units: dimensionless counts. Limiting case p → 1: pmf → (0, 0, 0, 0, 1), Var → 0.

**Common mistakes.** Using 1 − P(0) for P(≥ 3); forgetting that the Poisson rate scales with the interval length.

## 2. Moments and the ellipse

**Derivation.** (a) With μ = E[X], E[aX + b] = aμ + b, so Var(aX + b) = E[(aX + b − aμ − b)²] = a² E[(X − μ)²] = a² Var X. For a sum, Var(X + Y) = E[((X − μ_X) + (Y − μ_Y))²] = Var X + Var Y + 2 E[(X − μ_X)(Y − μ_Y)] = Var X + Var Y + 2 Cov(X, Y).
(b) E[X] = ∫ₐᵇ x dx/(b − a) = (b² − a²)/(2(b − a)) = (a + b)/2. Var X = ∫ₐᵇ (x − μ)² dx/(b − a); with u = x − μ running over [−w/2, w/2], w = b − a, this is (1/w)[u³/3] = (1/w)(2(w/2)³/3) = w²/12.
(c) E[T] = ∫₀^∞ t λe^{−λt} dt; by parts with u = t, dv = λe^{−λt} dt: [−t e^{−λt}]₀^∞ + ∫₀^∞ e^{−λt} dt = 0 + 1/λ. E[T²] = ∫₀^∞ t² λe^{−λt} dt = [−t² e^{−λt}]₀^∞ + 2∫₀^∞ t e^{−λt} dt = 2/λ², so Var T = 2/λ² − 1/λ² = 1/λ². Memoryless: P(T > s + t | T > s) = P(T > s + t)/P(T > s) = e^{−λ(s+t)}/e^{−λs} = e^{−λt} = P(T > t).
(d) Cov = VΛVᵀ with V orthogonal. **u** = Vᵀ**y** has covariance VᵀCov V = Λ, diagonal, and a linear map of a Gaussian vector is Gaussian, so u₁, u₂ are independent (for Gaussians, zero covariance means independence) with variances λ₁, λ₂. Then zᵢ = uᵢ/√λᵢ are independent N(0, 1) with joint density e^{−(z₁² + z₂²)/2}/(2π). The ellipse **y**ᵀCov⁻¹**y** ≤ k² is z₁² + z₂² ≤ k². In polar coordinates (z₁, z₂) = ρ(cos θ, sin θ), dz₁ dz₂ = ρ dρ dθ: P = ∫₀^{2π} dθ ∫₀ᵏ e^{−ρ²/2} ρ dρ/(2π) = [−e^{−ρ²/2}]₀ᵏ = 1 − e^{−k²/2}.

**Check.** (d) at k = 1 gives 0.3935, matching `scipy.stats.chi2.cdf(1, 2)`. (b) with [−½, ½] gives 1/12 ✓ (lesson). Dimensions: Var has units of X².

**Common mistakes.** Writing Var(aX + b) = a Var X + b; in (d), forgetting the Jacobian ρ, which gives a wrong 1 − e^{−k} form.

## 3. Micrometeoroid impacts ($s = 15$: A = 35 m²)

**Assumptions.** Constant flux, independent impacts, all of the array exposed; the flux figure is synthetic.

**Result.** (a) λ = 2 × 10⁻⁵ × 35 × 1095 = **0.7665** impacts; the count is Poisson(0.7665). (b) P(0) = e^{−0.7665} = **0.4646**; P(≥ 3) = 1 − e^{−λ}(1 + λ + λ²/2) = **0.0427**; the most probable count is **0** (mode = ⌊λ⌋ for non-integer λ). (c) Rate per day r = 2 × 10⁻⁵ × 35 = 7 × 10⁻⁴ d⁻¹; mean waiting time 1/r = **1429 days** (longer than the mission); P(first impact within 365 d) = 1 − e^{−0.2555} = **0.2255**. (d) z = (9 − 0.7665)/√0.7665 = **9.40**: nine impacts would be wildly inconsistent with the assumed flux (the exact Poisson tail P(≥ 9) is below 10⁻⁷); either the flux is an order of magnitude higher or the counting includes smaller particles. (e) At 30 (z = 33) the first assumption to question is *independence and constant rate*: a single fragmentation event or a meteor shower produces a cluster, and a clustered count has variance far above its mean, which is the signature to look for in the timing of the impacts.

**Check.** P(0) from the Poisson count, e^{−λ}, must equal the exponential waiting time's P(T > 1095 d) = e^{−1095 r} = e^{−0.7665} ✓.

**Common mistakes.** Using the Poisson rate per day in the Poisson formula for the whole mission; quoting the most probable count as λ.

## 4. A windier site ($s = 15$: σ_w = 4.5 m/s)

**Result.** (a) Cov = [[473 625, 121 500], [121 500, 48 025]] m²; σ_x = **688.2 m**, σ_y = **219.1 m**, ρ = **0.806** (up from 0.664: the wind now dominates both axes). (b) 99 % ellipse: k = 3.035, semi-axes **2159 m** and **381 m**, orientation **14.86°**, area **2.585 km²**. CEP = **497 m**; 99 % circle radius **1837 m**, area 10.60 km², **4.10** times the ellipse (up from 2.81: the ellipse got more elongated, aspect ratio 5.66). (c) 20 000 synthetic landings with seed 15: fraction inside the 99 % ellipse **0.9908 ± 0.0007**, 1.1 SE from 0.99. (d) The major axis rotated by only **+0.09°** (14.86° against 14.78°). The off-diagonal covariance comes entirely from the wind column of J, whose direction is (150, 40), at atan(40/150) = 14.93° from down-range; as the wind dominates, the ellipse aligns ever more closely with that direction, and it was already within 0.15° of it.

**Check.** Only Σ₁₁ changed, so Cov_xy = 150 × 40 × 4.5² = 121 500 ✓ and Cov_xx − 220 500 = 150²(4.5² − 3²) = 253 125 ✓. The semi-minor axis is unchanged to 0.1 % because the minor direction is almost perpendicular to the wind's.

**Common mistakes.** Scaling the whole Σ by (4.5/3)²; comparing the simulated fraction to 0.99 without a standard error.

## 5. How fast is the central limit theorem? ($s = 15$)

**Result.** (a) Uniform sums, 200 000 draws: gaps 0.0584, 0.0167, 0.0086, 0.0047, 0.0031, 0.0025 for n = 1, 2, 4, 8, 16, 32; power-law fit over n ≥ 2 gives α = **0.69**. (b) The 200 000-draw empirical cdf has a sampling scatter of about 0.002 in its largest deviation (the Kolmogorov statistic scales as 1/√N); the n = 16 and n = 32 points (0.0031, 0.0025) are within a factor 1.5 of that floor, which is why the fitted α is below the asymptotic 1 that symmetric distributions show. To measure further, raise `n_trials` to 10⁷ or more, or compare at a few fixed quantiles where the sampling error is smaller. (c) Exponential sums: gaps 0.159, 0.0947, 0.0679, 0.0485, 0.0341, 0.0239; α = **0.50**. The approach is slower and the exponent is exactly the Berry–Esseen 1/√n: the exponential is skewed (skewness 2), the standardised sum's skewness decays as 2/√n, and that skewness term dominates the cdf error. The uniform has zero skewness, so its leading correction is the kurtosis term, which decays as 1/n.

**Check.** At n = 1 the uniform gap 0.0584 is the exact maximum distance between the standardised uniform cdf and Φ (attained at z = ±√3... computed analytically, max|F_U − Φ| ≈ 0.058) and the exponential gap 0.159 is max|1 − e^{−(1+z)} − Φ(z)| ≈ 0.159, both independent of the sample.

**Common mistakes.** Including n = 1 in the power-law fit (the distribution at n = 1 is not in the asymptotic regime); reading a slope from points that sit on the sampling floor.

## 6. What the ellipse does not know

(a) z = 1200/469.6 = 2.555, so P(x > 1200 m) = 1 − Φ(2.555) = **0.0053** under the Gaussian model. With a t₄ wind error of the same standard deviation the probability is **higher**: a t distribution with ν = 4 has tails falling as |x|⁻⁵ rather than e^{−x²/2}, so beyond about 2.5σ it carries several times the Gaussian probability. A 2 × 10⁶-draw simulation (seed 16) with the wind drawn from t₄ scaled to 3 m/s gives 0.0106, a factor **2.0** higher. The 99.5 % statement became a 98.9 % one.

(b) The new element is Σ₁₄ = Σ₄₁ = 0.5 × 3.0 × 0.5 = 0.75 (m/s·deg). The off-diagonal of Cov gains J₁₁J₂₄Σ₁₄ = 150 × 250 × 0.75 = 28 125 and the cross-range variance gains 2 J₂₁J₂₄Σ₁₄ = 2 × 40 × 250 × 0.75 = 15 000: Cov = [[220 500, 82 125], [82 125, 45 025]]. `dispersion_ellipse` gives a major axis at **21.6°**, rotated *towards* cross-range from 14.8°, because a positive wind error now also drags the landing to the right through the heading.

(c) The 99.9 % ellipse (k = 3.717) is less trustworthy than the 99 % one for three reasons. First, it is a statement about 1 landing in 1000, and the whole figure rests on the Gaussian shape of the input tails at 3.7σ, which nothing in the analysis tests and which gusty winds are known to violate (part (a) showed a factor 2 at 2.5σ). Second, even the synthetic check is weak: of 20 000 simulated landings only 13 fall outside that ellipse, so the simulation confirms the number to ±30 % at best, and a real campaign of a few dozen flights cannot confirm it at all. Third, the inputs' standard deviations are themselves estimates with perhaps 20 % uncertainty, and at k = 3.7 a 20 % error in σ changes the exceedance probability by a factor of about 5. Evidence to support a 99.9 % claim would be measured wind-error distributions out to 4σ (hundreds of soundings) and a nonlinear trajectory Monte Carlo with those distributions as inputs; without them, publish the 99 % ellipse and state the assumption.

## Check-yourself answers

1. pmf: P(X = x) for a discrete variable; pdf: density f(x) whose integrals are probabilities, for a continuous one; cdf: F(x) = P(X ≤ x), defined for both.
2. E[Y] = 2 × 3 − 5 = 1, Var Y = 4 × 4 = 16. Mean of 16 copies: E = 3, Var = 4/16 = 0.25 (σ = 1).
3. With p = λ/n, np = λ exactly and np(1 − p) = λ(1 − λ/n) → λ. On count data, a sample variance much larger than the sample mean ("over-dispersion") says the events are clustered or the rate varies, so the Poisson model is wrong.
4. The probability of surviving a further time t does not depend on the age already reached. It rules out wear-out (a hazard rate that rises with age) and infant mortality (one that falls).
5. The standardised sum (or mean) of n independent identically distributed variables with finite variance tends to a standard normal as n → ∞; the mean of n measurements has standard deviation σ/√n.
6. In two dimensions the probability inside the k-ellipse is 1 − e^{−k²/2}, which at k = 1 is 0.3935: the radial density ρ e^{−ρ²/2} has little weight near the centre, so more probability sits beyond one standard deviation than in one dimension.
7. Cov(**y**) = J Cov(**x**) Jᵀ; the derivation uses linearity of expectation, E[J **x** **x**ᵀ Jᵀ] = J E[**x** **x**ᵀ] Jᵀ.
8. CEP = 1.1774 × 100 = 117.7 m; 99 % circle radius 100√(−2 ln 0.01) = 303.5 m. For an elongated covariance with the same σ_x but a smaller σ_y, the CEP is *smaller* (87.0 m for σ_y = 50 m): less total variance means a tighter circle.
