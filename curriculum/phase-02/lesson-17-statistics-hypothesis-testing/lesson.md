# P02.17  Statistics & Hypothesis Testing

> Eight coupons cannot tell you the strength of a batch. They can tell you how wrong you might be, and whether that is small enough to sign.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the sampling distribution of the mean, the unbiasedness of the n − 1 sample variance, the t-based confidence interval and its equivalence to a one-sided test, the Welch standard error, and the noncentral-t expression for power. 2. Implement and test, from scratch, sample statistics, confidence intervals, one-sample and Welch t-tests, a proportion test, Monte Carlo estimates of the Type I error rate and power, a sample-size planner, a percentile bootstrap, a permutation test and a one-sided normal tolerance bound, against `scipy.stats.ttest_1samp`, `ttest_ind`, `binomtest` and `nct` as oracles. 3. Measure, by seeded simulation, the coverage of a 95 % interval, the Type I error rate of a 5 % test, the power of an 8-coupon acceptance test and the coverage of a tolerance bound. 4. Decide whether a synthetic batch of aluminium sheet meets a 440 MPa requirement, how many coupons a supplier must test to demonstrate a 10 MPa margin with 90 % power, and when the normal approximation to a proportion test gives the wrong verdict. |
| **Prerequisites** | `P02.15` (probability rules, Monte Carlo error), `P02.16` (normal distribution, central limit theorem, sample mean and variance), `P00.14` (measurement uncertainty), `P00.03`, `P00.06`, `P00.07`. Skill check: you can compute a mean and standard deviation by hand for six numbers, read a z-score, and state what σ/√n is. |
| **Engineering question** | Eight tensile coupons from a batch of aluminium sheet average 446 MPa against a design requirement that the batch mean exceed 440 MPa and a supplier certificate claiming 470 MPa: can the batch be accepted, is the certificate credible, and how many coupons would settle the question? |
| **Assumptions and limits** | Coupons are independent random draws from one normally distributed population (the t procedures are exact only then, and approximately right for moderate n by the CLT; the bootstrap and permutation tests relax normality but not independence). Historical scatter σ = 12 MPa is treated as known only for power planning, never in the test itself. A hypothesis test controls the probability of a *wrong rejection*; it does not measure the probability that a hypothesis is true (Lesson 20 does). All coupon data are *synthetic*, drawn from stated normal distributions with stated seeds; the 440 MPa requirement and the 12 MPa scatter are illustrative, not values from any material specification. Tolerance bounds assume normality in the tail, where it is least checkable. |
| **Analytical work** | Sampling distribution of x̄; E[s²] = σ²; Student's t and why it replaces z; confidence interval and its interpretation; null and alternative hypotheses; test statistic, p-value, significance level; Type I and II errors; the test–interval duality; power as a noncentral-t tail; sample-size from the power equation; Welch–Satterthwaite degrees of freedom; the normal approximation to a proportion and its validity condition; the one-sided tolerance factor from the noncentral t. |
| **Computational deliverable** | `starter/main.py`: `sample_mean`, `sample_variance`, `standard_error`, `synthetic_batch`, `t_critical`, `confidence_interval_mean`, `p_value_from_t`, `one_sample_t_test`, `welch_t_test`, `proportion_z_test`, `decide`, `type_i_error_rate_mc`, `power_monte_carlo`, `power_analytic`, `required_sample_size`, `bootstrap_ci_mean`, `permutation_test`, `one_sided_tolerance_factor`, `lower_tolerance_bound`, `tolerance_bound_coverage_mc`. Tests prove: statistics match NumPy; interval coverage is 95 % within Monte Carlo error over 2000 synthetic batches; t and p match SciPy to 10⁻¹⁰ for every alternative; two-sided p is twice the smaller one-sided p; Welch df lies between its bounds and matches SciPy; the permutation p agrees with Welch on normal data; the Type I error rate is α within error; Monte Carlo power matches the noncentral t and is monotone in effect and n; the sample-size planner is minimal; the tolerance factor exceeds z_p, decreases with n and has the stated coverage. |
| **Practical activity** | **computational experiment**: coverage of confidence intervals and tolerance bounds, Type I error and power against the analytic predictions, and the sample-size curve; all on synthetic coupon data. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `acceptance_testing.png`, the hand derivations, the four simulation tables, and a one-page acceptance memo stating the decision, its p-value, the power of the test and the recommended coupon count. |
| **Assessment** | `assessments/phase-02/grading-rubric.md`; competencies 1–4. Pass: derivations complete, including the test–interval duality and E[s²] = σ²; all public tests pass; the hand interval and t statistic for the 8-coupon batch agree with the code to three figures; the memo names H0 and H1, gives the p-value and the power at a 10 MPa margin, and distinguishes "not demonstrated" from "not in spec". |
| **References** | `@openstax-introductory-statistics-2e` §8.1 "A Single Population Mean using the Normal Distribution", §8.2 "A Single Population Mean using the Student t-Distribution", §9.1 "Null and Alternative Hypotheses", §9.2 "Outcomes and the Type I and Type II Errors", §9.3 "Probability Distribution Needed for Hypothesis Testing", §9.4 "Rare Events, the Sample, Decision and Conclusion", §9.5 "Additional Information and Full Hypothesis Test Examples", §9.6 "Hypothesis Testing of a Single Mean and Single Proportion" (all required-free); Chapter 2 "Descriptive Statistics" (supplementary-free, for mean and standard deviation). Power, the Welch test, the bootstrap, permutation tests and tolerance bounds are not covered there; parts 2.5–2.8 are original instruction. |
| **Extensions** | ★ Implement the paired t-test for coupons cut from the same sheet before and after a heat treatment, and show when pairing beats the two-sample test. ★ Replace the normal population by a lognormal with the same mean and variance and measure how the coverage of the t interval and of the tolerance bound degrade at n = 8. ★ Implement the two-sided tolerance interval and compare its factor with the one-sided one. |
| **Time** | unmeasured |

## 1. The engineering problem

A batch of aluminium alloy sheet arrives for a wing-rib production run. The design requires that the batch's mean ultimate tensile strength exceed 440 MPa; the supplier's certificate states 470 MPa; the quality engineer cuts eight coupons and pulls them. The eight results (synthetic, generated from N(452, 12²) with seed 17 to stand in for a real batch) average 446.12 MPa with a sample standard deviation of 11.38 MPa. Three decisions hang on them. Can the batch be *accepted*, meaning the evidence shows the mean exceeds 440 MPa? Is the certificate's 470 MPa *credible*? And if neither can be settled with eight coupons, how many would it take?

Each of these is a question about a population (the whole batch) asked of a sample (eight coupons), and the honest answer has two parts: a number, and a statement of how far the number could be from the truth. Lesson 16 gave the second part its foundation: the mean of n samples from a population with standard deviation σ has standard deviation σ/√n. This lesson turns that fact into **confidence intervals** (a range that captures the population mean with stated probability), **hypothesis tests** (a rule for deciding, with a stated rate of wrong decisions), **power** (the probability that the test catches what it is meant to catch) and **tolerance bounds** (the strength that most of the population exceeds, which is what a structural allowable actually is).

Success criteria, written down first: the from-scratch t statistics and p-values match `scipy.stats` to 10⁻¹⁰; 2000 synthetic 95 % intervals cover the true mean 95 % ± 2 % of the time; the simulated Type I error rate of the 5 % test is 0.050 ± 0.006; the Monte Carlo power at a 10 MPa margin agrees with the noncentral-t value 0.683 within four standard errors (0.013); and the one-sided tolerance bound covers its target fraction with the stated 95 % confidence to within 1 %.

## 2. Understand and derive

### 2.1 From a sample to a population

A **population** is everything we could measure (every possible coupon from the batch); a **sample** is the n we did. The population has a mean μ and standard deviation σ we do not know; the sample gives estimates, the **sample mean** x̄ = Σxᵢ/n and the **sample variance**

**s² = Σ(xᵢ − x̄)²/(n − 1)**.   (1)

Why n − 1? Because x̄ is closer to the data than μ is: Σ(xᵢ − x̄)² ≤ Σ(xᵢ − μ)² always, by the minimising property of the mean. Taking expectations, E[Σ(xᵢ − x̄)²] = E[Σ(xᵢ − μ)²] − n E[(x̄ − μ)²] = nσ² − n(σ²/n) = (n − 1)σ², so dividing by n − 1 makes E[s²] = σ²: the estimate is **unbiased**. For the eight coupons, x̄ = 446.12 MPa and s = 11.38 MPa.

By the central limit theorem (P02.16), x̄ is approximately N(μ, σ²/n); for a normal population it is exactly so. The **standard error** of the mean is σ/√n, estimated by s/√n = 11.38/√8 = 4.02 MPa. Eight coupons locate the batch mean to about ±4 MPa, not to ±11.

### 2.2 Confidence intervals and Student's t

If σ were known, (x̄ − μ)/(σ/√n) would be standard normal, and P(−1.96 < z < 1.96) = 0.95 rearranges to the 95 % **confidence interval** x̄ ± 1.96σ/√n (`@openstax-introductory-statistics-2e` §8.1). We do not know σ; replacing it by s adds uncertainty, and the ratio

**t = (x̄ − μ)/(s/√n)**   (2)

follows **Student's t distribution** with n − 1 **degrees of freedom** (§8.2): bell-shaped, symmetric, but with heavier tails than the normal because a small sample can under-estimate σ. Its two-sided 95 % point is 2.365 for 7 degrees of freedom (against 1.960 for the normal), 2.045 for 29, and it tends to 1.960 as n → ∞. The interval is

**x̄ ± t*_{n−1} s/√n**,   (3)

for the batch 446.12 ± 2.365 × 4.02 = (436.6, 455.6) MPa. The interval's meaning is a statement about the *procedure*: if we drew many samples of 8 and formed (3) each time, 95 % of the intervals would contain μ. It is not "μ has a 95 % probability of lying in (436.6, 455.6)": μ is a fixed number, and this interval either contains it or does not. The third panel of `acceptance_testing.png` draws 40 synthetic intervals; about two miss.

A **bootstrap** interval (2.7) makes the same statement without assuming a normal population: (439.0, 453.7) MPa here, slightly narrower because the percentile bootstrap at n = 8 under-covers a little.

### 2.3 Hypothesis tests: the vocabulary

A **hypothesis test** (§9.1) sets up a **null hypothesis** H0, the claim that the data must *disprove*, against an **alternative** H1. For acceptance the burden of proof lies on the batch: H0: μ = 440 MPa (the batch is only just at the limit, or below), H1: μ > 440 MPa. This is a **one-sided** test; a question like "does the batch match the certificate?" is **two-sided**, H0: μ = 470, H1: μ ≠ 470.

The **test statistic** is (2) with μ = μ₀; under H0 it follows the t distribution of §9.3. The **p-value** (§9.4) is the probability, *if H0 were true*, of a statistic at least as extreme as the one observed: P(T ≥ t) for H1: μ > μ₀, P(T ≤ t) for μ < μ₀, and 2P(T ≥ |t|) for the two-sided case. The **significance level** α (usually 0.05) is the largest p-value we will treat as "too rare to believe H0": p < α rejects H0; otherwise we **fail to reject** it, which is not the same as accepting it.

For the batch against 440 MPa: t = (446.12 − 440)/4.02 = 1.521, df = 7, p = 0.086. At α = 0.05 we fail to reject H0: eight coupons averaging 6 MPa above the limit, with a standard error of 4 MPa, have *not demonstrated* that the batch mean exceeds 440. Against the certificate's 470 MPa: t = −5.935, two-sided p = 0.0006: the certificate is not credible for this batch.

**Test–interval duality.** Rejecting H0: μ = μ₀ at level α in a two-sided test is exactly the statement that μ₀ lies outside the (1 − α) confidence interval (3), because both say |x̄ − μ₀| > t* s/√n. The one-sided test at α corresponds to the one-sided (1 − α) lower bound x̄ − t_{1−α} s/√n = 446.12 − 1.895 × 4.02 = 438.5 MPa, which is below 440: the same verdict, read off an interval. Report the interval; it carries the decision and the size of the effect at once.

### 2.4 Two kinds of error, and power

The decision can be wrong in two ways (§9.2):

| | H0 true (batch at 440) | H1 true (batch above 440) |
|---|---|---|
| reject H0 (accept the batch) | **Type I error**, probability α | correct, probability = **power** 1 − β |
| fail to reject (hold the batch) | correct | **Type II error**, probability β |

α is set by the analyst; β depends on how far the truth is from H0, on σ and on n. The **power** is the probability of rejecting H0 when the true mean is μ₀ + δ. Under H1 the statistic (2) with μ₀ in the numerator follows a **noncentral t** distribution with noncentrality δ√n/σ, so

**power(δ) = P(T′_{n−1, δ√n/σ} > t_{1−α, n−1})**.   (4)

`power_analytic` evaluates (4) with `scipy.stats.nct`; `power_monte_carlo` simulates it by drawing synthetic batches from N(μ₀ + δ, σ²) and counting rejections. For n = 8, σ = 12 MPa, α = 0.05 one-sided:

| true margin δ above 440 MPa | power, noncentral t (4) | power, 20 000 synthetic batches |
|---|---|---|
| 0 (H0 true) | 0.050 | 0.0496 ± 0.0015 |
| 5 | 0.281 | 0.281 |
| 10 | 0.683 | 0.686 |
| 15 | 0.936 | 0.938 |
| 20 | 0.995 | 0.995 |

A batch whose true mean is 450 MPa, comfortably in spec, passes this eight-coupon test only 68 % of the time; one at 445 MPa passes 28 % of the time. The test protects the designer (Type I at 5 %) and exposes the supplier to a large Type II risk. The remedy is n: solving (4) for the smallest n with power ≥ 0.9 at δ = 10 MPa gives **n = 14** (`required_sample_size`); for 0.8, 11; for 0.95, 18. Power roughly follows δ√n/σ, so halving the detectable margin costs four times the coupons.

### 2.5 Two samples: Welch's t

To compare the batch with ten reference coupons from a certified lot (synthetic, N(470, 12²), seed 1701, mean 470.19 MPa) we test H0: μ₁ = μ₂. The difference of means has variance σ₁²/n₁ + σ₂²/n₂ (variances of independent quantities add, P02.16), estimated by s₁²/n₁ + s₂²/n₂, so

**t = (x̄₁ − x̄₂)/√(s₁²/n₁ + s₂²/n₂)**,   df = (s₁²/n₁ + s₂²/n₂)² / [(s₁²/n₁)²/(n₁ − 1) + (s₂²/n₂)²/(n₂ − 1)].   (5)

The **Welch–Satterthwaite** degrees of freedom in (5) lie between min(n₁, n₂) − 1 and n₁ + n₂ − 2 and do not assume equal variances; the classical "pooled" test does, and there is rarely a reason to make that assumption. Batch against reference: t = −4.147, df = 15.88, p = 0.0008. The batch is weaker than the reference lot by 24.1 ± 5.8 MPa (the standard error from (5)).

### 2.6 A proportion: rivets, and a warning

Of 120 rivets inspected, 5 are defective; the process is supposed to run at 2 % or better. H0: p = 0.02, H1: p > 0.02. The count is binomial (P02.15), and for large n the sample proportion p̂ = k/n is approximately normal with standard deviation √(p₀(1 − p₀)/n) (§9.6):

**z = (p̂ − p₀)/√(p₀(1 − p₀)/n)** = (0.0417 − 0.02)/0.01278 = 1.695,   p = 0.045.   (6)

Reject at 5 %? Not so fast. The normal approximation needs n p₀ and n(1 − p₀) both at least about 10, and n p₀ = 2.4. The exact binomial probability of 5 or more defects out of 120 at p = 0.02, `scipy.stats.binomtest`, is 0.094: *fail to reject*. The approximation flips the verdict. When the count is small, use the exact test; `proportion_z_test` is kept in the module because it is what the textbook teaches and what the test compares against the exact result.

### 2.7 Resampling: the bootstrap and the permutation test

Both t procedures assume a normal population. Two **resampling** methods make weaker assumptions and are worth having from scratch.

The **percentile bootstrap** treats the sample as the best available picture of the population: draw n values from it *with replacement*, compute the mean, repeat B times, and take the 2.5 % and 97.5 % quantiles of the B means as the interval. It is exact as n → ∞ and slightly too narrow at n = 8 (it cannot see tails the sample did not contain), which is why its interval (439.0, 453.7) is 4.3 MPa narrower than the t interval.

The **permutation test** asks: if the batch and reference coupons came from the *same* distribution (the sharp null), how often would a random split of the 18 pooled values into groups of 8 and 10 produce a mean difference as large as the one observed? Shuffle, split, count. With 5000 shuffles the p-value is 0.0022: the same verdict as Welch (p = 0.0008), with no normality assumption. The permutation p-value cannot fall below 1/(B + 1), so B sets the resolution.

### 2.8 Tolerance bounds: what an allowable really is

A confidence interval locates the *mean*. A structural allowable must bound the *weak coupons*: the number that 90 % (or 99 %) of the population exceeds. Were μ and σ known, that would be μ − z₀.₉₀ σ = μ − 1.282σ. With estimates from n coupons it must be

**lower tolerance bound = x̄ − k s**,   k = t′_{n−1, z_p√n}(conf)/√n,   (7)

where k comes from the noncentral t quantile and is chosen so that, with probability conf over repeated samples, the bound lies below the true p-quantile. For p = 0.90 and conf = 0.95 (the definition used for "B-basis"-style allowables), k = 2.582 for n = 8, 1.777 for n = 30 and 1.354 for n = 1000, tending to z₀.₉₀ = 1.282 as the estimates sharpen. For the batch: 446.12 − 2.582 × 11.38 = **416.7 MPa**, against 431.5 MPa from the naive x̄ − 1.282 s. The 15 MPa gap is the price of estimating μ and σ from eight coupons; a naive bound lies below the true 10th percentile only about 46 % of the time instead of 95 %, which the simulation in part 4 shows: it is a coin toss, not a bound. A tolerance bound is a statement about the tail of a *normal* population, and the tail is exactly where normality is least verifiable: that is the assumption to flag in the memo.

### 2.9 Where the model stops being valid

- **Independence.** Eight coupons cut from the same corner of one sheet are not eight independent draws from the batch; the effective n is smaller and every interval here is too narrow. Sample across sheets and positions.
- **Normality.** The t interval is robust for the mean at n ≥ 8 unless the population is strongly skewed; the tolerance bound is not robust at all (extension 2).
- **Known σ for planning.** Power and sample size used σ = 12 MPa from history; if the new batch scatters more, the power is lower than planned.
- **Normal approximation to a proportion.** Fails for n p₀ < 10, as 2.6 showed.
- **What a p-value is not.** p = 0.086 is not "an 8.6 % probability that the batch is bad". It is the probability of data this extreme *if* the batch were exactly at the limit. The probability that H0 is true needs a prior and Bayes' rule (Lesson 20).

## 3. Build the computational model

Open `starter/main.py`. Student-t and noncentral-t quantiles come from `scipy.stats`; everything else is yours.

### Step 1: statistics and intervals

```python
def sample_variance(x):
    m = sample_mean(x)
    return float(np.sum((x - m) ** 2) / (x.size - 1))                 # (1), two-pass

def confidence_interval_mean(x, conf=0.95):
    m, se = sample_mean(x), standard_error(x)
    half = t_critical(conf, x.size - 1) * se                           # (3)
    return m - half, m + half
```

`t_critical(conf, df)` is `stats.t.ppf(0.5 + 0.5 * conf, df)`. Raise `ValueError` for fewer than two observations.

### Step 2: the tests

```python
def p_value_from_t(t, df, alternative):
    if alternative == "two-sided": return float(2.0 * stats.t.sf(abs(t), df))
    if alternative == "greater":   return float(stats.t.sf(t, df))
    return float(stats.t.cdf(t, df))

def one_sample_t_test(x, mu0, alternative="two-sided"):
    m, se = sample_mean(x), standard_error(x)
    t = (m - mu0) / se                                                 # (2)
    return {"t": t, "df": float(x.size - 1), "p": p_value_from_t(t, x.size - 1, alternative), "mean": m, "se": se}
```

`welch_t_test` implements (5) and returns `diff` and `se` as well; `proportion_z_test` implements (6) with `stats.norm`; `decide(p, alpha)` returns the two strings. Validate `alternative` against the three allowed values.

### Step 3: errors and power

```python
def power_monte_carlo(delta, sigma, n, alpha=0.05, n_trials=4000, seed=0, alternative="greater", mu0=0.0):
    rng = np.random.default_rng(seed)
    samples = rng.normal(mu0 + delta, sigma, size=(n_trials, n))     # synthetic batches
    t = (samples.mean(axis=1) - mu0) / (samples.std(axis=1, ddof=1) / math.sqrt(n))
    p = stats.t.sf(t, n - 1) if alternative == "greater" else ...
    return float(np.mean(p < alpha))
```

Vectorise over trials: 20 000 batches of 8 take milliseconds. `type_i_error_rate_mc` is the same call with δ = 0. `power_analytic` is (4) with `stats.nct.sf`, handling the two-sided case as the sum of both tails; `required_sample_size` loops n upward until it reaches the target.

### Step 4: resampling

`bootstrap_ci_mean` draws an `(n_boot, n)` array of indices with `rng.integers`, takes row means and `np.quantile`. `permutation_test` pools, permutes with `rng.permutation`, splits at `x.size`, and counts |difference| ≥ observed; return (count + 1)/(n_perm + 1).

### Step 5: tolerance bounds

```python
def one_sided_tolerance_factor(n, p=0.90, conf=0.95):
    zp = stats.norm.ppf(p)
    return float(stats.nct.ppf(conf, n - 1, zp * math.sqrt(n)) / math.sqrt(n))      # (7)
```

`lower_tolerance_bound` is x̄ − k s; `tolerance_bound_coverage_mc` draws N(0, 1) samples and counts how often the bound lies below the true p-quantile −z_p.

Running `main()` prints the tables of 2.2–2.8 and writes `acceptance_testing.png` (the t densities under H0 and H1 with the critical and observed values, the power curves for n = 4, 8, 16, 32 with the Monte Carlo points, and 40 synthetic confidence intervals).

### Verify

Run `python -m pytest curriculum/phase-02/lesson-17-statistics-hypothesis-testing`. The tests:

1. **Analytic cases.** Mean, variance and standard error match NumPy; `t_critical(0.95, 7)` = 2.365 and → 1.960 for large df; t and p match `scipy.stats.ttest_1samp` for all three alternatives to 10⁻¹⁰; the two-sided p is twice the smaller one-sided p; Welch matches `ttest_ind(equal_var=False)`; the proportion z matches the normal approximation and its p is within 0.06 of the exact binomial for n p₀ ≥ 10.
2. **Limiting cases.** μ₀ = x̄ gives t = 0 and two-sided p = 1; the statistic is invariant to adding a constant to both data and μ₀ and to scaling both; Welch df lies in [min(n) − 1, n₁ + n₂ − 2]; power at δ = 0 equals α; power → 1 for large δ; the tolerance factor exceeds z_p and tends to it as n grows.
3. **Coverage and error rates (seeded, 4 SE tolerances stated in comments).** 2000 synthetic 95 % intervals cover the true mean 95 % ± 2 %; the Type I error rate at n = 8 is 0.05 ± 0.014 over 4000 batches; Monte Carlo power matches (4) within 4 SE at four margins; the tolerance bound covers its quantile 95 % ± 1.4 % of the time.
4. **Duality and monotonicity.** The one-sided test rejects exactly when the one-sided lower bound exceeds μ₀; power increases with δ and with n; `required_sample_size` returns an n whose power reaches the target while n − 1 does not; interval half-width scales as 1/√n for fixed s.
5. **Resampling.** The bootstrap interval of a 200-point normal sample agrees with the t interval to 10 %; the permutation p-value agrees with the Welch p within a factor 3 on normal data and is ≥ 1/(B + 1).
6. `main()` writes `acceptance_testing.png` to the directory it is given.

### Compare

`scipy.stats.ttest_1samp`, `ttest_ind(equal_var=False)`, `binomtest` and `scipy.stats.bootstrap` do all of this and more, with the numerical care (log-space tails, exact binomial sums) that matters for p-values near 10⁻¹⁰. Our versions agree to 10⁻¹⁰ where both exist. We wrote them to make three things visible that the library hides: that a p-value is a tail area of a specific distribution under a specific null; that power is the same tail area under a shifted distribution, which is why it depends on n through δ√n/σ; and that "confidence" is a property of a procedure measured over repetitions, which is why we simulated 2000 of them. In practice: call the library for the test, and keep a seeded Monte Carlo to check that your α and power are what you think they are.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All data are synthetic, drawn from stated distributions with stated seeds.

1. **Coverage.** Draw 2000 batches of n = 8 from N(452, 12²) with seed 173 and count how many 95 % intervals (3) contain 452. Predict 0.95 ± 0.005 (one SE). Then repeat with a *lognormal* population of the same mean and variance and record the coverage.
2. **Type I error.** `type_i_error_rate_mc(440, 12, 8, 0.05, 20000, seed=170, alternative="greater")`: 0.0496 ± 0.0015. Repeat at n = 3 and n = 50 to see that the t test holds α at every n (the normal-z test, by contrast, over-rejects at small n: try replacing `stats.t` by `stats.norm` and watch n = 3 give about 0.11).
3. **Power.** Tabulate `power_analytic` and `power_monte_carlo` at δ = 0, 5, 10, 15, 20 MPa for n = 8 (table in 2.4) and plot the power curves for n = 4, 8, 16, 32. Read off the margin that each n detects with 90 % power: about 21, 14, 9.5, 6.5 MPa.
4. **Tolerance-bound coverage.** `tolerance_bound_coverage_mc(8, 0.90, 0.95, 20000, seed=172)` gives 0.9501. Replace k by the naive 1.282 and measure the coverage: about 0.46.

Raw data format: CSV tables `population, coverage`; `n, type_I_rate, se`; `delta, power_analytic, power_mc`; `method, coverage`.

What you should find: coverage is nominal for the normal population and a little low for the lognormal at n = 8 (the mean of a skewed sample is itself skewed); the Type I rate is α at every n for the t test; Monte Carlo power sits on the noncentral-t curve within its standard error; and the naive tolerance bound sits above the true 10th percentile more than half the time, because x̄ − 1.282 s is an unbiased *estimate* of that percentile, not a bound on it.

## 5. Validate and improve

| quantity | hand (§2) | code / library / simulation | difference | cause |
|---|---|---|---|---|
| s of the batch | 11.38 MPa | 11.38 (`sample_variance`) = `np.std(ddof=1)` | 0 | same formula |
| 95 % CI | 446.12 ± 2.365 × 4.02 = (436.6, 455.6) | (436.61, 455.63) | rounding | — |
| t, p against 440 (one-sided) | 1.52, 0.086 | 1.521, 0.08598 = `ttest_1samp` | < 10⁻¹⁰ | same distribution |
| Welch t, df, p | — | −4.147, 15.88, 0.00077 = `ttest_ind` | < 10⁻¹⁰ | — |
| permutation p | — | 0.0022 (5000 shuffles) | factor 2.9 above Welch | resolution 1/5001 and no normality assumed |
| proportion test p | 0.045 (normal) | 0.094 (exact binomial) | verdict flips | n p₀ = 2.4 < 10 |
| Type I rate, n = 8 | 0.050 | 0.0496 ± 0.0015 | 0.3 SE | sampling |
| power at δ = 10 | 0.683 (noncentral t) | 0.686 ± 0.003 | 0.8 SE | sampling |
| tolerance factor k, n = 8 | 2.582 | coverage 0.9501 ± 0.0015 | 0.1 SE | sampling |

Every discrepancy is either Monte Carlo sampling error of the predicted size or a named failure of an approximation (the normal approximation to the binomial). The permutation and Welch p-values differ, as they should: they test different nulls (identical distributions against equal means) with different assumptions, and both are far below 0.05.

What the agreement does not show: that real coupons behave like independent normal draws. The simulation draws from the model the analysis assumes, so this is verification of the code, not validation of the procedure on real material. Real acceptance testing adds things this lesson omits: coupons from several sheets and positions to make independence plausible, a normality check (a probability plot) before quoting a tolerance bound, and an allowable derived from pooled historical data rather than from one batch of eight.

Improvements worth making now: implement the paired test for before/after comparisons on the same coupons (extension 1); run the coverage experiment on a skewed population to see the tolerance bound fail (extension 2); replace the normal approximation in `proportion_z_test` by an exact binomial sum. Improvement worth *not* making: lowering α to 0.01 to "be safe": it raises β, so good batches are held more often, and it does nothing about the independence assumption, which is the real weakness.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivations of E[s²] = σ², the test–interval duality and (5); `acceptance_testing.png`; the four tables from part 4; and a one-page acceptance memo (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating H0 and H1, the t statistic and p-value, the decision at α = 0.05, the power of the eight-coupon test at a 10 MPa margin, the number of coupons recommended for 90 % power, the tolerance bound with its assumption, and the difference between "not demonstrated" and "not in spec".

### Adversarial self-checks

- **What if the number of coupons doubles?** The standard error falls by √2, the interval half-width from 9.5 to about 6.1 MPa (t* also drops, 2.365 → 2.131), and the same 6.12 MPa margin would give t = 2.15, p = 0.024: the batch would pass. Did the batch get better? No. The *evidence* got better. That is the whole point of power.
- **What if σ goes to zero? To infinity?** σ → 0: any margin is detected with one coupon; σ → ∞: no finite n detects anything. Power depends on δ/σ, not on δ.
- **Which assumption fails first in the real world?** Independence of the coupons (same sheet, same position, same operator). Second, normality of the tail for the tolerance bound.
- **What measurement could falsify the model?** Coupons from several sheets whose between-sheet variance exceeds the within-sheet variance: then the eight coupons from one sheet were never a sample of the batch, and the intervals here are too narrow by a factor you can compute from a one-way analysis of variance (Phase 09 revisits this).
- **Does the result survive a change of method?** t, Welch and permutation agree on the certificate question; t and bootstrap interval endpoints agree within 2.5 MPa; the proportion question does *not* survive the change from the normal approximation to the exact test, and the exact test is the one to trust.
- **Which decision would change if the uncertainty doubled?** With s = 22.8 MPa the one-sided p would rise to 0.24 and the lower tolerance bound would fall to 387 MPa; the batch would be held, and the supplier would need 46 coupons rather than 14 to demonstrate a 10 MPa margin with 90 % power. The design would also be asked whether a 440 MPa *mean* requirement, rather than a tolerance-bound allowable, was ever the right specification.

## Check yourself

1. Why is the sample variance divided by n − 1, and what property does that buy?
2. A 95 % confidence interval for a batch mean is (436.6, 455.6) MPa. Which of these is correct: "95 % of coupons lie in this range", "the batch mean lies in this range with probability 0.95", or "95 % of intervals built this way contain the batch mean"?
3. Write H0 and H1 for an acceptance test that must *demonstrate* a batch mean above 440 MPa, and say why the burden of proof is placed that way.
4. Define Type I and Type II errors for that test, and say which one the supplier cares about.
5. The one-sided p-value is 0.086. State in one sentence what that number is the probability of, and in another what it is not.
6. Power at a 10 MPa margin is 0.68 with 8 coupons. Roughly how many coupons give 0.9, and why does the answer scale the way it does?
7. Why does Welch's test not assume equal variances, and between which two values do its degrees of freedom lie?
8. A tolerance bound for 90 % of the population at 95 % confidence uses k = 2.58 for n = 8 but only 1.28 for n → ∞. Explain the difference in one sentence.
