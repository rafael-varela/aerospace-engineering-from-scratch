# P02.17  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. Every test must state H0, H1, the statistic, its distribution under H0, the p-value and the decision at the stated α; every interval must state its confidence level. Label every simulated dataset as synthetic.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: six coupons by hand)

Six tensile coupons (synthetic) give, in MPa: 448, 455, 441, 462, 450 and 439 + s/5.

(a) Compute x̄, s and the standard error by hand. (b) Form the 95 % confidence interval for the batch mean (t* for 5 degrees of freedom is 2.571). (c) Test H0: μ = 440 against H1: μ > 440 at α = 0.05: statistic, p-value (use `scipy.stats.t.sf` or a table), decision. (d) Would the two-sided test against 440 reach the same decision? (e) Confirm every number with `confidence_interval_mean` and `one_sample_t_test`.

Answer form: three statistics, an interval, two test results with p-values.

## 2. (Derivation: the machinery)

(a) Prove that E[s²] = σ² for s² = Σ(xᵢ − x̄)²/(n − 1), using Var(x̄) = σ²/n.
(b) Prove Var(x̄) = σ²/n for independent observations, stating where independence is used.
(c) Show that the one-sided test of H0: μ = μ₀ against μ > μ₀ at level α rejects exactly when the one-sided (1 − α) lower confidence bound x̄ − t_{1−α} s/√n exceeds μ₀.
(d) Derive the Welch standard error √(s₁²/n₁ + s₂²/n₂) for the difference of two independent sample means, and explain in one sentence why the pooled-variance alternative is a special case.

Answer form: four derivations.

## 3. (Unfamiliar scenario: two propellant lots)

Two lots of a solid rocket motor are static-fired, six motors each, and the total impulse is recorded in N·s. Generate the (synthetic) data with `synthetic_batch(6, 2450.0, 60.0, seed=s)` for lot A and `synthetic_batch(6, 2400.0, 60.0, seed=s + 100)` for lot B.

(a) Give both means and standard deviations. (b) Test whether the lots differ with `welch_t_test` (two-sided) and with `permutation_test` (5000 shuffles, seed s); report both p-values and the decision at α = 0.05. (c) Which of the two tests would you trust more with n = 6 per lot, and why? (d) The flight requirement is a total impulse of at least 2300 N·s for *every* motor. Explain why neither test answers that question, and which quantity from the lesson does.

Answer form: four statistics, two p-values with a decision, two short arguments.

## 4. (Numerical, uses your code: planning the acceptance test)

A new alloy scatters more than the lesson's: σ = 10 + s/5 MPa. The acceptance test is one-sided at α = 0.05 against 440 MPa.

(a) With `power_analytic`, tabulate the power to demonstrate an 8 MPa margin for n = 4, 8, 16, 32. (b) With `required_sample_size`, find the n for 90 % power at 8 MPa. (c) Check the Type I error rate and the power at n = 8 by `type_i_error_rate_mc` and `power_monte_carlo` with 4000 synthetic batches and seed s; give each with its standard error and compare with α and with (a). (d) Draw 2000 synthetic batches of 8 from N(452, σ²) with seed s and measure the coverage of the 95 % interval; compare with 0.95 and its standard error.

Answer form: a four-entry table, one n, two Monte Carlo rates with standard errors, one coverage with standard error.

## 5. (Numerical, uses your code: the allowable)

Twelve synthetic coupons: `synthetic_batch(12, 455.0, 10 + s/10, seed=s)`.

(a) Give x̄ and s, the tolerance factor k for 90 % of the population at 95 % confidence, and the lower tolerance bound. (b) Compare it with the naive bound x̄ − 1.282 s and with the lower end of the 95 % confidence interval for the mean, and explain in one sentence what each of the three numbers is a statement about. (c) With `tolerance_bound_coverage_mc(12, 0.90, 0.95, 4000, seed=s)` measure the coverage of the bound; then measure the coverage of the naive bound by adapting the function (replace k by 1.282). (d) How many coupons would bring k within 0.1 of 1.282?

Answer form: four numbers, three one-sentence interpretations, two coverages, one n.

## 6. (Judgment: what the p-value does and does not say)

(a) The lesson's batch gave p = 0.086 against 440 MPa and was not accepted. A colleague says "so there is an 8.6 % chance the batch is bad". Correct the statement in two sentences.
(b) The same colleague proposes to pull four more coupons and, if the combined p falls below 0.05, accept the batch; otherwise pull four more, and so on. Explain what this does to the Type I error rate, and propose a correct procedure (with a number of coupons from the lesson).
(c) With the observed x̄ = 446.12 and s = 11.38 MPa held fixed, how many coupons would make the one-sided p-value fall below 0.05? Why is this calculation *not* a valid plan (compare with (b)), and what is it good for?
(d) The rivet example gave p = 0.045 by the normal approximation and 0.094 exactly. State the rule that tells you which to use, and give the smallest sample size at p₀ = 0.02 for which the approximation is acceptable.

Answer form: two sentences, a short paragraph with a number, a number with two sentences, a rule with a number.
