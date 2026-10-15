# P02.20  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. Prior vector (π, 1 − π), likelihood vector (0.95, 0.05 + s/1000). 2. `bayes_discrete` or by hand: posterior = 0.95π/[0.95π + (0.05 + s/1000)(1 − π)]. 3. Odds: π/(1 − π) × 0.95/(0.05 + s/1000); convert back with p = O/(1 + O). 4. Second alarm: feed the posterior back as the prior (or square the Bayes factor). 5. The explanation turns on the base rate: with a rare fault the false alarms outnumber the true ones even for a good detector; a second independent alarm squares the Bayes factor.

## 2.
(a) Multiply (4) by (5), collect λ^{α+k−1}e^{−(β+T)λ}, and integrate with the Gamma-function identity; the constants cancel in the ratio. (b) Write (α + k)/(β + T) = (β/(β + T))(α/β) + (T/(β + T))(k/T). (c) ∫e^{−λt}β^αλ^{α−1}e^{−βλ}/Γ(α) dλ = β^α/(β + t)^α; Jensen: for convex f, E[f(λ)] ≥ f(E[λ]). (d) log p = k ln λ − λT + const; ∂²/∂λ² = −k/λ²; I = E[k]/λ² = T/λ; Jeffreys ∝ √I ∝ λ^{−1/2}. For θ = 1/λ: log p = −k ln θ − T/θ; I(θ) = T/θ³; prior ∝ θ^{−3/2}; transform with |dθ/dλ| = λ^{−2}: θ^{−3/2}λ^{−2} = λ^{3/2}λ^{−2} = λ^{−1/2}. (e) Gamma(k + 1, T) density T^{k+1}λ^k e^{−Tλ}/k!; substitute λ = x/(2T) and compare with the χ²_{2k+2} density x^k e^{−x/2}/(2^{k+1}k!).

## 3.
(a) log posterior = −(b − μ₀)²/(2σ₀²) − (z − b)²/(2σ_z²) + const; expand, collect b²: −b²(1/σ₀² + 1/σ_z²)/2; collect b: b(μ₀/σ₀² + z/σ_z²); complete the square. (b) Numbers from the formulas. (c) Grid on [−3, 5] with 80 001 points; `grid_posterior`, `grid_mean`, SD by `np.trapezoid((g − mean)² pdf)`, `credible_interval(…, 0.9)`; compare with `norm.ppf`. (d) K → 1 (trust the measurement entirely) and K → 0 (ignore it).

## 4.
(a) `gamma_poisson_update(2.5, 1.2e6, 1, T)`; mean α/β; SD √α/β; `gamma.ppf([0.05, 0.95], α, scale=1/β)`; `prob_rate_below`, `posterior_predictive_reliability`, `hours_to_demonstrate(…, 0.9)`. (b) `grid_posterior` with `log_gamma_pdf` and `log_poisson_likelihood`; `grid_mean`, `credible_interval`, `hpd_interval`. (c) `metropolis` with a log target that returns −∞ for λ ≤ 0; `autocorrelation_time`; n_eff = n/τ; scatter of the mean σ/√n_eff, of the SD σ/√(2n_eff), of a probability p: √(p(1 − p)/n_eff). (d) `frequentist_upper_bound(1, T, 0.95)` vs `gamma.ppf(0.95, α, scale=1/β)`; the gap is the heritage hours.

## 5.
(a) Same failure mechanisms with the same stresses; e.g. a software change outside the bearing/lubricant path preserves, a new bearing supplier destroys; the heritage speaks for the mechanisms it exercised, not for new ones. (b) `gamma_poisson_update(0.5, 0, 0, 1.2e6)` then update with (1, 6e4); `prob_rate_below`, `posterior_predictive_reliability`. (c) P(λ ≤ λ_req) is the probability that the *rate* meets the requirement given everything known; it is neither a pass/fail verdict nor the reliability itself; R(3 yr) is the probability the wheel survives, which is what the mission cares about and is 0.93 even though the rate is "probably" above the requirement. (d) 1 − (β/(β + 6 × 10⁴))^α with the heritage prior; compare with a conventional surprise threshold; what would be surprising is several failures.

## 6.
(a) Three calls to `credible_interval_coverage`; scatter √(level(1 − level)/4000). (b) Copy the function body, replace the draw of λ by `rng.gamma(2.5, 5 / 1.2e6)` (NumPy's gamma takes shape and *scale*), keep the heritage prior in the update. (c) With T ≪ β the posterior is the prior and the intervals sit near 2 per 10⁶ h while the true rates average 10; coverage recovers only when T is large enough for the data to outweigh β = 1.2 × 10⁶ h. (d) The prior predictive probability of k failures or more in T hours: for k = 1 it is 1 − (β/(β + T))^α; evaluate for the lesson fleet, and note what value of k would have given a probability below, say, 0.05.
