# P02.16  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) p = 0.90 + s/1000; pmf C(4, k)pᵏ(1 − p)⁴⁻ᵏ for k = 0..4; E = Σ k P(k) and 4p; Var = Σ (k − E)² P(k) and 4p(1 − p); P(X ≥ 3) = P(3) + P(4). (b) λ = 1 + s/10; P(0) = e^{−λ}; P(≥ 3) = 1 − P(0) − P(1) − P(2); quarter: e^{−3λ}. (c) `binomial_pmf(np.arange(5), 4, p)`, `poisson_pmf(k, lam)`, `discrete_mean`, `discrete_variance`.

## 2.
(a) Let Y = aX + b, E[Y] = aμ + b; Var Y = E[(a(X − μ))²] = a² Var X. For the sum, expand E[((X − μ_X) + (Y − μ_Y))²]. (b) E = ∫ₐᵇ x/(b − a) dx; Var = ∫ₐᵇ (x − μ)²/(b − a) dx; substitute u = x − μ over [−w/2, w/2]. (c) ∫₀^∞ t λe^{−λt} dt by parts with u = t; E[T²] = 2/λ² by parts twice; Var = 2/λ² − 1/λ². Memoryless: P(T > s + t)/P(T > s) with the survival function e^{−λt}. (d) Cov = VΛVᵀ; u = Vᵀy has independent components with variances λᵢ; zᵢ = uᵢ/√λᵢ are independent N(0, 1); the joint density is e^{−ρ²/2}/(2π); integrate over the disk ρ ≤ k in polar coordinates.

## 3.
(a) λ = 2 × 10⁻⁵ × A × 1095; Poisson(λ). (b) e^{−λ}; 1 − Σₖ₌₀² P(k); mode = ⌊λ⌋. (c) Rate per day r = 2 × 10⁻⁵ A; mean 1/r days; P(T ≤ 365) = 1 − e^{−365r}. (d) z = (9 − λ)/√λ; |z| < 2 is unremarkable. (e) Constant flux (debris environment changes; a shower), or independence (a fragmentation event produces a cluster).

## 4.
(a) `Sigma = np.diag([(3 + s/10)**2, 0.05**2, 0.3**2, 0.5**2])`; `cov = propagate_covariance(J, Sigma)`; σ = √diag; ρ = cov₀₁/(σ_xσ_y). (b) `dispersion_ellipse(cov, 0.99)`; `radius_for_probability(cov, 0.5)` and `(cov, 0.99)`; ratio πR²/area. (c) `simulate_landings(20_000, J, Sigma, seed=s)`; `fraction_inside_ellipse(..., k_for_probability(0.99))`; SE = √(0.99 × 0.01/20 000). (d) New angle minus 14.8°; the wind column of J has direction (150, 40), at atan(40/150) = 14.9° from down-range, so a wind-dominated covariance aligns with it.

## 5.
(a) Loop over n, collect gaps, `np.polyfit(np.log(n[1:]), np.log(gap[1:]), 1)`. (b) Points below ≈ 0.003 are at the sampling floor; increase `n_trials` or use a smoother statistic (e.g. the gap at fixed quantiles). (c) `(rng.exponential(size=(n_trials, n)).sum(axis=1) - n) / np.sqrt(n)`; the exponential has skewness 2, which decays as 2/√n and controls the gap; the uniform has skewness 0 and its first non-Gaussian correction is the kurtosis term, decaying as 1/n.

## 6.
(a) `1 - normal_cdf(1200, 0, 469.6)`; z = 2.56; the t₄ tail beyond 2.5σ is a few times heavier than the Gaussian tail, so the probability rises by a factor of about 2–4. (b) `Sigma[0, 3] = Sigma[3, 0] = 0.5 * 3.0 * 0.5`; the extra term in Cov_xy is J₀₀J₁₃Σ₀₃ + J₁₀J₁₃Σ₃₃... write it out and see which components grow; recompute the angle. (c) Out of 20 000 synthetic landings only about 20 fall outside the 99.9 % ellipse, and no real campaign has thousands of flights; the figure rests entirely on the Gaussian tail assumption for the inputs.
