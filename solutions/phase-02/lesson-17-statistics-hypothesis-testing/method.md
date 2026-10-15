# P02.17  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) x̄ = Σx/6; s² = Σ(x − x̄)²/5; SE = s/√6. (b) x̄ ± 2.571 SE. (c) t = (x̄ − 440)/SE; p = P(T₅ > t) = `stats.t.sf(t, 5)`; compare with 0.05. (d) two-sided p = 2 × one-sided p. (e) `confidence_interval_mean(x)`, `one_sample_t_test(x, 440, "greater")`.

## 2.
(a) Σ(xᵢ − x̄)² = Σ(xᵢ − μ)² − n(x̄ − μ)² (the cross term is −2(x̄ − μ)Σ(xᵢ − μ) = −2n(x̄ − μ)²). Take expectations: nσ² − n·σ²/n = (n − 1)σ². (b) Var(Σxᵢ/n) = (1/n²)Σ Var xᵢ because covariances vanish for independent variables. (c) Reject ⇔ t > t_{1−α} ⇔ x̄ − μ₀ > t_{1−α}s/√n ⇔ x̄ − t_{1−α}s/√n > μ₀. (d) Var(x̄₁ − x̄₂) = Var x̄₁ + Var x̄₂ = σ₁²/n₁ + σ₂²/n₂ by independence; estimate each σ by its sample s. The pooled test sets σ₁ = σ₂ and uses a weighted average of s₁², s₂².

## 3.
(a) `sample_mean`, `sample_variance` on each lot. (b) `welch_t_test(a, b)`; `permutation_test(a, b, 5000, seed=s)`. (c) With n = 6 normality cannot be checked; the permutation test needs only exchangeability under H0, but has resolution 1/5001 and tests "same distribution", not "same mean". (d) Both compare means; the requirement is on every motor, i.e. a lower tail quantile: `lower_tolerance_bound` with p = 0.99 (or whatever fraction the requirement specifies).

## 4.
(a) `power_analytic(8.0, sigma, n, 0.05, "greater")` for the four n. (b) `required_sample_size(8.0, sigma, 0.05, 0.9, "greater")`. (c) `type_i_error_rate_mc(440, sigma, 8, 0.05, 4000, seed=s, alternative="greater")`; `power_monte_carlo(8.0, sigma, 8, 0.05, 4000, seed=s, alternative="greater")`; SE = √(p(1 − p)/4000). (d) Loop: draw `rng.normal(452, sigma, 8)`, call `confidence_interval_mean`, count hits; SE = √(0.95 × 0.05/2000).

## 5.
(a) `sample_mean`, `sample_variance`, `one_sided_tolerance_factor(12, 0.90, 0.95)`, `lower_tolerance_bound(batch)`. (b) Naive: x̄ − 1.282 s; CI lower end from `confidence_interval_mean`. The three are: a 95 %-confident bound on the 10th percentile; a point estimate of the 10th percentile; a 97.5 %-confident bound on the mean. (c) Call the function; copy its body with k = 1.282. (d) Evaluate `one_sided_tolerance_factor(n)` for increasing n until k − 1.282 < 0.1.

## 6.
(a) p = P(data at least this extreme | μ = 440), not P(μ ≤ 440 | data). (b) Each look is a 5 % chance of a false acceptance; looks accumulate; fix n from `required_sample_size` and test once. (c) Solve 6.12/(11.38/√n) > t_{0.95, n−1} by trying n; it is a what-if on the power, not a stopping rule, because the next coupons will not reproduce x̄ and s. (d) n p₀ ≥ 10 (and n(1 − p₀) ≥ 10): n ≥ 10/p₀.
