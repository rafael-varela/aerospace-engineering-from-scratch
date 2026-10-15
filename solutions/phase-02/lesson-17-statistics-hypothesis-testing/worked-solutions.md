# P02.17  Worked solutions

Each solution states: assumptions, units, the derivation, the result (with p-values, confidence levels or Monte Carlo standard errors), at least one independent check, and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py`. Every dataset is synthetic.

## 1. Six coupons by hand ($s = 15$: sixth coupon 442 MPa)

**Assumptions.** Independent coupons from one normal population.

**Derivation.** Data 448, 455, 441, 462, 450, 442. x̄ = 2698/6 = **449.67 MPa**. Deviations −1.67, 5.33, −8.67, 12.33, 0.33, −7.67; squares sum to 317.33; s² = 317.33/5 = 63.47, s = **7.97 MPa**; SE = 7.97/√6 = **3.25 MPa**. (b) 95 % CI: 449.67 ± 2.571 × 3.25 = **(441.3, 458.0) MPa**. (c) H0: μ = 440, H1: μ > 440; t = (449.67 − 440)/3.25 = **2.972** with 5 degrees of freedom; p = P(T₅ > 2.972) = **0.0155** < 0.05: reject H0, the batch mean is demonstrated above 440 MPa. (d) Two-sided p = 2 × 0.0155 = 0.031, also below 0.05: same decision here, but with a smaller margin (a one-sided test has more power for the direction it looks in).

**Check.** Duality: the one-sided 95 % lower bound 449.67 − 2.015 × 3.25 = 443.1 exceeds 440 ✓. `one_sample_t_test(x, 440, "greater")` returns t = 2.9722, p = 0.015539.

**Common mistakes.** Dividing by 6; using the normal 1.96 instead of t₅ = 2.571; reporting the two-sided p for a one-sided question.

## 2. The machinery

**Derivation.** (a) xᵢ − x̄ = (xᵢ − μ) − (x̄ − μ). Squaring and summing: Σ(xᵢ − x̄)² = Σ(xᵢ − μ)² − 2(x̄ − μ)Σ(xᵢ − μ) + n(x̄ − μ)² = Σ(xᵢ − μ)² − n(x̄ − μ)², since Σ(xᵢ − μ) = n(x̄ − μ). Expectations: E[Σ(xᵢ − μ)²] = nσ², E[n(x̄ − μ)²] = n Var x̄ = σ². So E[Σ(xᵢ − x̄)²] = (n − 1)σ² and E[s²] = σ².
(b) Var(Σxᵢ/n) = (1/n²) Var(Σxᵢ) = (1/n²)[Σ Var xᵢ + Σᵢ≠ⱼ Cov(xᵢ, xⱼ)]; independence kills every covariance, leaving nσ²/n² = σ²/n.
(c) Reject ⇔ t = (x̄ − μ₀)/(s/√n) > t_{1−α} ⇔ x̄ − μ₀ > t_{1−α} s/√n ⇔ x̄ − t_{1−α} s/√n > μ₀. The last expression is the one-sided lower (1 − α) confidence bound.
(d) x̄₁ and x̄₂ are independent with variances σ₁²/n₁ and σ₂²/n₂ (by (b)), so Var(x̄₁ − x̄₂) = σ₁²/n₁ + σ₂²/n₂; replacing each σ² by s² gives the Welch standard error. The pooled test assumes σ₁ = σ₂ = σ and estimates σ² by the weighted average [(n₁ − 1)s₁² + (n₂ − 1)s₂²]/(n₁ + n₂ − 2); when s₁ = s₂ the two standard errors coincide.

**Check.** For the lesson's batch, the one-sided bound 438.5 < 440 and p = 0.086 > 0.05 agree. Units: variances in MPa² throughout.

**Common mistakes.** In (a), treating x̄ as a constant; in (c), using the two-sided t*.

## 3. Two propellant lots ($s = 15$)

**Result.** (a) Lot A (seed 15): 2364.1, 2393.8, 2473.6, 2418.6, 2481.5, 2498.4 N·s; mean **2438.4**, s **54.0 N·s**. Lot B (seed 115): 2461.6, 2282.2, 2356.5, 2455.5, 2400.7, 2346.1; mean **2383.8**, s **69.2 N·s**. (b) Welch: difference 54.6 ± 35.9 N·s, t = **1.523**, df = 9.44, two-sided p = **0.161**; permutation (5000 shuffles, seed 15): p = **0.149**. Both fail to reject H0 at α = 0.05: the lots are not shown to differ. (c) With six motors per lot a normality check is impossible, so the t test rests on an untestable assumption; the permutation test needs only that under H0 the twelve labels are exchangeable, which is what "same lot behaviour" means. Trust the permutation p more, with two caveats: its resolution is 1/5001, and it tests "same distribution", not "same mean", so a difference in scatter alone could reject it. (d) Both tests compare *means*. The requirement is that *every* motor exceed 2300 N·s, a statement about the lower tail of the population: the lower tolerance bound (p = 0.99 or higher, 95 % confidence) is the quantity to compute, and with n = 6 the factor k is large (k = 5.06 for 99 %/95 %), so the bound for lot B would be far below 2300 N·s and the lot could not be cleared from six firings.

**Check.** The 95 % CI for the difference, 54.6 ± 2.25 × 35.9 = (−26, 135) N·s, contains zero, consistent with p = 0.16.

**Common mistakes.** Concluding "the lots are the same" from p = 0.16; using the pooled test without justification; confusing a mean comparison with a per-unit requirement.

## 4. Planning the acceptance test ($s = 15$: σ = 13 MPa)

**Result.** (a) Power to demonstrate an 8 MPa margin, one-sided α = 0.05: n = 4: **0.249**; 8: **0.470**; 16: **0.759**; 32: **0.961**. (b) 90 % power needs **n = 25**. (c) 4000 synthetic batches, seed 15: Type I rate **0.0523 ± 0.0034** (0.7 SE from 0.05); power at n = 8 **0.475 ± 0.008** against the analytic 0.470 (0.6 SE). (d) Coverage of the 95 % interval over 2000 batches: **0.948 ± 0.005**, 0.4 SE from 0.95.

**Check.** Power depends on δ√n/σ: the lesson's (δ, σ, n) = (10, 12, 8) gives 2.36 and power 0.683; here (8, 13, 8) gives 1.74 and power 0.470, lower as it must be. Doubling n from 8 to 16 raises the noncentrality by √2.

**Common mistakes.** Quoting a Monte Carlo rate without its standard error; using the two-sided alternative in the power calculation for a one-sided test.

## 5. The allowable ($s = 15$: σ = 11.5 MPa, 12 coupons, seed 15)

**Result.** (a) x̄ = **455.09 MPa**, s = **12.48 MPa**; k(12; 0.90, 0.95) = **2.210**; lower tolerance bound 455.09 − 2.210 × 12.48 = **427.5 MPa**. (b) Naive x̄ − 1.282 s = **439.1 MPa**; lower end of the 95 % CI for the mean = **447.2 MPa**. The tolerance bound is a value that, with 95 % confidence, at least 90 % of the population exceeds; the naive number is a point estimate of the 10th percentile with no confidence attached; the CI end is a 97.5 %-confident lower bound on the *mean*, which says nothing about weak coupons. (c) Coverage of the proper bound over 4000 synthetic samples: **0.9485** (target 0.95, SE 0.0034); coverage of the naive bound: **0.46**: it lies above the true 10th percentile more than half the time. (d) k falls to 1.382 (within 0.1 of 1.282) at **n ≈ 530** coupons (k = 1.385 at n = 500, 1.376 at 600).

**Check.** k(12) = 2.210 > z₀.₉₀ = 1.282 ✓, and the bound is below the CI's lower end ✓ (a tail bound must sit below a mean bound).

**Common mistakes.** Using the two-sided t* in place of k; interpreting the tolerance bound as a bound on the mean; expecting the naive bound to be "a bit optimistic" rather than useless as a bound.

## 6. What the p-value does and does not say

(a) The p-value 0.086 is the probability of observing a sample mean at least 6.1 MPa above 440 *if the true batch mean were exactly 440 MPa*; it conditions on the hypothesis, not on the data. The probability that the batch is bad given the data needs a prior probability for batches being bad and Bayes' rule (Lesson 20), and is not 8.6 %.

(b) Each look at the data with a 5 % threshold is a fresh 5 % chance of a false acceptance when the batch is at the limit; with three looks (8, 12, 16 coupons) the overall Type I rate rises to roughly 11 %, and with unlimited looks it tends to 1: a batch at exactly 440 MPa will eventually be "accepted". The correct procedure fixes n *before* testing from the power calculation, here **n = 14** for 90 % power at a 10 MPa margin, tests once, and reports the result whatever it is.

(c) With x̄ = 446.12 and s = 11.38 held fixed, t = 6.12√n/11.38 exceeds the one-sided 5 % critical value at **n = 12** (t = 1.86, p = 0.045; at n = 11, p = 0.052). It is not a plan because the next coupons will not reproduce this x̄ and s: it is the same optional stopping as (b), dressed up. It is useful as a what-if: it says that *if* the batch really is 6 MPa above the limit with this scatter, about 12 coupons are what a demonstration costs, which matches the power calculation's 11–14.

(d) Use the normal approximation only when n p₀ ≥ 10 and n(1 − p₀) ≥ 10 (some texts say 5); otherwise use the exact binomial test. At p₀ = 0.02 that needs **n ≥ 500**. With n = 120 the two p-values were 0.045 and 0.094, and the exact one is the one to report.

## Check-yourself answers

1. Because x̄ is closer to the data than μ is, Σ(xᵢ − x̄)² under-estimates Σ(xᵢ − μ)² by exactly σ² on average; dividing by n − 1 makes E[s²] = σ² (unbiasedness).
2. Only the third. The first confuses the interval for the mean with the spread of individual coupons; the second treats the fixed μ as random.
3. H0: μ = 440 MPa (or μ ≤ 440), H1: μ > 440 MPa. The batch must *disprove* the null of non-compliance, so a Type I error (accepting a bad batch) is the error whose rate is controlled at α.
4. Type I: accepting a batch whose mean is at (or below) 440; Type II: holding a batch whose mean is above 440. The supplier cares about Type II, which is the risk of a good batch being held.
5. It is the probability, if the batch mean were exactly 440 MPa, of a sample mean at least as far above 440 as the one observed. It is not the probability that the batch is bad.
6. About 14 (`required_sample_size` gives 14 for 90 %). Power depends on δ√n/σ, so the noncentrality grows as √n: going from 0.68 to 0.90 needs the noncentrality to rise from 2.36 to about 3.1, i.e. n up by (3.1/2.36)² ≈ 1.75.
7. It estimates each population variance separately and adds s₁²/n₁ + s₂²/n₂, so unequal variances are allowed; the Welch–Satterthwaite df lies between min(n₁, n₂) − 1 and n₁ + n₂ − 2.
8. With n → ∞ the mean and standard deviation are known exactly and the 10th percentile is μ − 1.282σ; with 8 coupons both are uncertain, and k must be inflated to 2.58 so that the bound still lies below the true percentile in 95 % of samples.
