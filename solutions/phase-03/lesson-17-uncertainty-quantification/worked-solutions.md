# P03.17  Worked solutions

Each solution states: assumptions, units, the derivation, the result with its uncertainty, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown. Every input distribution is an assumption of the problem statement.

## 1. Linear propagation by hand ($s = 15$: I_sp = 315 s, MR = 3.5)

**Derivation.** Δv₀ = 9.80665 × 315 × ln 3.5 = 9.80665 × 315 × 1.2528 = **3869.9 m/s**. In logarithms, ln Δv = ln g₀ + ln I_sp + ln(ln MR), so the relative sensitivity to I_sp is 1 and to MR is 1/ln MR = 0.798. u_rel² = (0.01)² + (0.798 × 0.02)² = 1.0 × 10⁻⁴ + 2.549 × 10⁻⁴ = 3.549 × 10⁻⁴, u_rel = 1.884 %, **u = 72.9 m/s**. Shares: I_sp **0.282**, mass ratio **0.718**.

(b) z = (3700 − 3869.9)/72.9 = −2.33, **P(Δv < 3700) = Φ(−2.33) = 0.0099 ≈ 1 %**.

(c) The mass ratio owns 72 % of the variance; halving its uncertainty (to 1 %) gives u_rel² = 1.0 × 10⁻⁴ + 0.637 × 10⁻⁴, u_rel = 1.28 %, u = 49.5 m/s: a **32 % reduction**. Halving the I_sp uncertainty instead would give u = 64.6 m/s, only 11 %.

**Check.** The code's `linear_propagation` on f(ξ) = g₀ I_sp(1 + 0.01ξ₁) ln(MR(1 + 0.02ξ₂)) returns 72.90 m/s and shares (0.282, 0.718) ✓; a degree-4 PCE gives mean 3869.3 m/s (the logarithm is concave, so the mean shifts *down* by 0.6 m/s) and std 72.93 m/s; Monte Carlo with 10⁵ samples gives P(Δv < 3700) = 0.0107 ± 0.0003, 8 % above the Gaussian 0.0099: the concavity thickens the lower tail slightly.

**Common mistakes.** Using the relative uncertainty of MR directly instead of dividing by ln MR; adding standard deviations instead of variances.

## 2. Errors, orthogonality, coefficients

(a) E[μ̂] = (1/N) Σ E[y_n] = μ: unbiased. Var μ̂ = (1/N²) Σ Var y_n = σ²/N by independence (covariances vanish). For p̂ = (1/N) Σ 1[y_n < y_crit], each indicator is Bernoulli(p) with variance p(1 − p), so Var p̂ = p(1 − p)/N. Relative accuracy ε means √(p(1 − p)/N) = εp, i.e. **N = (1 − p)/(p ε²)** ≈ 1/(p ε²) for small p: 10⁴/p at ε = 1 %, 100/p at ε = 10 %.

(b) Let G(t) = e^{ξt − t²/2} = Σ_n He_n(ξ) tⁿ/n!. Then G(s)G(t) = e^{ξ(s + t) − (s² + t²)/2}, and since E[e^{ξu}] = e^{u²/2} for a standard normal, E[G(s)G(t)] = e^{(s + t)²/2 − (s² + t²)/2} = e^{st} = Σ_n sⁿtⁿ/n!. Expanding the left side, Σ_{m,n} E[He_m He_n] sᵐtⁿ/(m! n!), and matching the coefficient of sᵐtⁿ: E[He_m He_n]/(m! n!) = δ_mn/n!, so **E[He_m He_n] = n! δ_mn**.

(c) e^{aξ} = e^{a²/2} · e^{aξ − a²/2} = e^{a²/2} G(a) = e^{a²/2} Σ He_n(ξ) aⁿ/n!, so **c_n = e^{a²/2} aⁿ/n!**. Then Σ c_n² n! = e^{a²} Σ a^{2n}/n! = e^{a²} e^{a²} = e^{2a²} = E[e^{2aξ}] ✓ (Parseval for the expansion). Numerically for a = 0.7: the code's coefficients match to 10⁻⁸ and the variance e^{a²}(e^{a²} − 1) = 1.0276 to 10⁻⁵ at degree 6 (the truncation drops terms of relative size a¹⁴/7! ≈ 10⁻⁶).

(d) E[f] = Σ c_α E[He_α] = c₀, because E[He_α] = E[He_α He₀] = 0 for α ≠ 0. E[f²] = Σ_{α,β} c_α c_β E[He_α He_β] = Σ c_α² α!, so Var f = Σ_{α ≠ 0} c_α² α!. For an additive model every α ≠ 0 in the expansion has exactly one nonzero component, so every variance term is counted in exactly one S_i and Σ S_i = Σ_{α ≠ 0} c_α² α!/Var f = 1.

**Common mistakes.** In (b), forgetting that E[e^{ξu}] = e^{u²/2} needs the −t²/2 in the generating function; in (d), dividing by E[f²] instead of the variance.

## 3. Landing distance ($s = 15$: V = 63 m/s, a = 2.5 m/s², runway 900 m)

(a) s₀ = 63²/5 = **793.8 m**. Relative sensitivities 2 × 0.03 = 0.06 (speed) and −0.08 (deceleration): u_rel = √(0.0036 + 0.0064) = **0.100**, u = **79.4 m**, shares **0.36 / 0.64**. Second derivatives in the standardised variables: ∂²s/∂ξ_V² = 2s₀ s_V² = 1.43 m, ∂²s/∂ξ_a² = 2s₀ s_a² = 10.16 m (from 1/(1 + s_aξ) ≈ 1 − s_aξ + s_a²ξ²), mean shift ½(1.43 + 10.16) = **+5.8 m** (799.6 m): the inverse deceleration is convex, so the mean roll exceeds the nominal.

(b)

| method | runs | mean [m] | std [m] |
|---|---|---|---|
| linear (first / second order) | 5 | 793.8 / 799.6 | 79.4 |
| PCE degree 1 | 4 | 799.63 | 79.94 |
| PCE degree 2 | 9 | 799.70 | 81.04 |
| PCE degree 3 | 16 | 799.70 | 81.08 |
| PCE degree 4 | 25 | 799.70 | 81.08 |
| Monte Carlo, 2 × 10⁵ | 2 × 10⁵ | 799.60 ± 0.18 | 81.05 |

The PCE standard deviation is 2 % above the linear one (the 1/a curvature adds variance); degree 2 is converged to 0.05 %.

(c) P(s_L > 900 m): Monte Carlo **0.1107 ± 0.0007**; degree-4 surrogate (4 × 10⁵ samples) **0.1116**; Gaussian approximation **0.0905**. The samples have skewness **+0.43**: the 1/a factor stretches the long-roll tail, and the Gaussian, which has none, underestimates the overrun probability by 18 %. For a runway requirement this is the difference between "1 in 11" and "1 in 9".

(d) Sobol: S = (0.350, 0.648), S_T = (0.352, 0.650); interactions 0.2 %. The deceleration (brakes, runway friction, reversers) dominates, and that is what the designer must allow for with runway margin; the touchdown speed, 35 % of the variance, is the pilot's: a speed discipline of ±1.5 % instead of ±3 % would cut its share to a tenth.

**Check.** Mean shift by hand 799.59 m vs PCE 799.70 m: the 0.1 m difference is the third- and fourth-order terms of 1/(1 + 0.08ξ). Units: m²/s² ÷ m/s² = m ✓.

**Common mistakes.** Using the relative sensitivity 1 instead of 2 for V²; comparing the Gaussian tail with the Monte Carlo tail and blaming the sampling.

## 4. Your own drag-polar uncertainty ($s = 15$: L/D ± 3.5 %, required 1400 km)

(a) Linear: mean 1502.22 km (second order 1503.53), std **74.69 km**, shares (0.364, 0.496, 0.140). PCE degree 3 (64 runs): mean **1503.53 km**, std **74.86 km**, Sobol S = (0.365, 0.494, 0.140), S_T = (0.366, 0.495, 0.140). Monte Carlo (seed 17):

| N | mean [km] | std [km] | P(R < 1400 km) |
|---|---|---|---|
| 10³ | 1503.4 ± 2.4 | 75.5 | 0.087 ± 0.009 |
| 10⁴ | 1503.5 ± 0.8 | 75.4 | 0.085 ± 0.003 |
| 10⁵ | 1503.6 ± 0.2 | 75.0 | 0.0810 ± 0.0009 |

(b) Surrogate (4 × 10⁵ samples): **0.0797**; direct Monte Carlo 0.0810 ± 0.0009; Gaussian 0.0856 (6 % high: with a smaller L/D spread the skewness of the curved inputs matters relatively more). With the drag polar known to 3.5 % instead of 5 % the required range moved from 1350 to 1400 km in the problem, so the probabilities are not comparable to the lesson's; what is comparable is the L/D share, which fell from 0.67 to 0.49.

(c) N = (1 − p)/(p × 0.01) = 0.92/(0.08 × 0.01) ≈ **1150 samples** for 10 % relative on p; the PCE used **64 model runs** (27 would have done: degree 2 is converged).

(d) "R = 1504 ± 75 km (one standard deviation), P(R < 1400 km) = 8.0 % ± 0.1 %, under assumed independent Gaussian inputs of 3 % (c_T), 3.5 % (L/D) and 2 % (W_fuel); L/D contributes 49 % of the variance, c_T 37 %; 64 model runs (PCE), checked by 10⁵ Monte Carlo runs."

**Common mistakes.** Quoting the N = 10³ probability without its ±0.009; forgetting to change the required range.

## 5. Variance reduction and the bootstrap (lesson inputs, seed 17)

(a)

| N | se plain [km] | se antithetic [km] | variance ratio | equivalent plain N |
|---|---|---|---|---|
| 10³ | 2.92 | 0.160 | 335 | 3.4 × 10⁵ |
| 10⁴ | 0.926 | 0.049 | 356 | 3.6 × 10⁶ |

(b) 5th percentile from 10⁴ samples: **1350.4 km**, bootstrap 95 % interval **[1346.8, 1354.0] km**. Degree-3 surrogate (4 × 10⁵ samples): **1355.1 km**; Gaussian approximation 1503.5 − 1.645 × 91.95 = **1351.0 km**. The interval does *not* cover the surrogate's value (it misses by 1.1 km), and 10⁵ direct samples give 1354.9 km, agreeing with the surrogate: the 10⁴-sample estimate happens to sit low. A 95 % interval fails one time in twenty by construction; this is such a time, and the only defence is more data (the 10⁵ interval is about ±1.2 km and does cover 1355).

(c) The pair average (f(ξ) + f(−ξ))/2 is the even part of f, whose distribution is far narrower than f's: its 5th percentile is near 1503 km, nowhere near the range's 1355 km. Antithetic sampling estimates the *mean* with the even part's small variance; it does not produce samples of the range distribution at all, so quantiles, tail probabilities and histograms must come from plain (or surrogate) sampling. The 350-fold gain is real and specific to the mean of a nearly linear model.

**Common mistakes.** Treating the pair averages as range samples; concluding from one missed interval that the bootstrap is wrong.

## 6. Certifying a probability

(a) **Methods, in order.** Linear propagation first (13 runs for six inputs) to get the scale and the variance shares; then a PCE by regression (degree 2 in six variables is 28 terms, so about 60–100 runs; degree 3 is 84 terms, 200 runs) or a sparse-grid projection, giving moments, Sobol indices and a surrogate; then Monte Carlo on the surrogate with at least 10⁶ samples for a 10⁻⁴ level (expected 100 exceedances, ±10 %); then 100–300 *direct* model runs placed in the tail region the surrogate identifies, to check that the polynomial is right where it matters. Confidence: observing k exceedances in N surrogate samples bounds p; 90 % confidence that p < 10⁻⁴ with, say, 60 observed exceedances in 10⁶ samples follows from the Poisson interval (upper 90 % bound ≈ 70/10⁶), so the surrogate must be sampled at least that much and the direct runs must confirm its tail.

(b) **Correlation and bounds.** Correlated inputs break the product structure of the Hermite basis: transform to independent standard normals with a Cholesky factor of the covariance (ξ = L⁻¹(x − μ)/σ) before building the expansion; the Sobol indices then refer to the transformed, independent variables and no longer attribute variance to the original two inputs separately. A physically bounded input (a scale height cannot be negative or exceed a known maximum) must get a bounded distribution; pushing a truncated Gaussian through Hermite polynomials converges slowly because the truncation is a discontinuity in ξ. Use the orthogonal family of the bounded distribution (Legendre for uniform, Jacobi for beta) or Monte Carlo for that input.

(c) **Tail disagreement.** Believe the 200 direct runs: they are the model, while the surrogate is a polynomial whose quadrature or regression points lie mostly within ±2σ, so its behaviour at ±4σ is extrapolation. A factor 2 at 10⁻⁴ is typical of a low-degree PCE in a tail. Remedies: raise the degree and watch the tail converge, fit a local surrogate in the tail region, or use the surrogate only to find the tail region and then run the model directly there (importance sampling). Do not average the two.

(d) **Statement.** "Under the expert-elicited input distributions (six inputs, independent after a Cholesky decorrelation of the two correlated ones, with the density scale height bounded by a beta distribution), the probability that the peak heat flux exceeds the allowable is 6 × 10⁻⁵, with a statistical uncertainty of ±1 × 10⁻⁵ from 10⁶ samples of a degree-3 surrogate, confirmed by 250 direct model runs in the tail that agree with the surrogate to 20 %. The 90 % upper confidence bound is 8 × 10⁻⁵, below the 10⁻⁴ requirement. The result is conditional on the elicited distributions, in particular on the ±6 % entry-angle spread, which carries 55 % of the variance (Sobol index); a tracking-based entry-angle dispersion from the previous three flights, or an arc-jet test of the allowable, would change this number more than any further computation."

**Common mistakes.** Reporting the surrogate's probability without direct runs in the tail; treating the elicited distributions as data; reporting a probability without its confidence bound.

## Check-yourself answers

1. μ̂ = (1/N)Σy_n, se = s/√N with s the sample standard deviation; halving the error needs four times the samples.
2. N ≈ (1 − p)/(pε²) = 0.999/(10⁻³ × 0.01) ≈ 10⁵ direct samples; the alternative is a PCE surrogate built from tens of runs and sampled 10⁵ times, checked by a few hundred direct runs in the tail.
3. X = F⁻¹(U) with U uniform on (0, 1) has distribution F. The normal distribution function has no closed-form inverse, so Box–Muller uses two uniforms and the polar form (radius² = −2 ln U₁, angle 2πU₂) to produce two independent normals exactly.
4. u² = Σ(∂f/∂ξ_i)² for independent unit-variance inputs; E[f] ≈ f(0) + ½Σ∂²f/∂ξ_i². For the range the correction is positive (+1.3 km) because 1/c_T is convex and dominates the small concave logarithm.
5. E[He_m He_n] = n! δ_mn; He₂ = ξ² − 1, He₃ = ξ³ − 3ξ.
6. p + 1 nodes per dimension (the integrand f He_α has degree ≤ 2p ≤ 2(p + 1) − 1), that is (p + 1)^d model runs: 27 for p = 2, d = 3.
7. Mean = c₀; Var = Σ_{α≠0} c_α² α!; S_i = (Σ over α with only α_i ≠ 0 of c_α² α!)/Var.
8. The variance carried by interaction terms (mixed multi-indices such as (1, 1, 0)): 0.09 % of the total, the degree to which the range is not an additive function of its three inputs.
