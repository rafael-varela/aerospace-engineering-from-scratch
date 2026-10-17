# P02.20  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer, say which results are exact and which are numerical or statistical, and name the prior behind every posterior.

Several problems use your personal seed **s**: the day of the month you were born (1–31). Every fleet or record in this set is synthetic, drawn with your own seeded code; say so wherever you report one.

## 1. (Warm-up: a telemetry alarm)

A fault-detection monitor on a satellite raises an alarm when it believes a thruster valve has stuck. From the fleet's history the prior probability that the valve is stuck at any given check is s/1000. The monitor's true-positive rate is 0.95 and its false-positive rate is 0.05 + s/1000.

(a) The alarm fires once. Compute P(stuck | alarm) with Bayes' theorem and with the odds form, and show they agree. (b) A second, independent check also alarms. Update again. (c) Explain in two sentences why a 95 % detector gives so little confidence after one alarm and so much after two.

Answer form: three probabilities, one Bayes factor, a short explanation.

## 2. (Derivation: the Gamma–Poisson machinery)

(a) Starting from the Poisson likelihood (4) and the Gamma prior (5), derive the posterior (6) and show the normalising integral is Γ(α + k)/(β + T)^{α+k}.
(b) Show that the posterior mean is the weighted average (7) and state the weights for the lesson's heritage prior and fleet data.
(c) Derive the posterior predictive (9) by evaluating the integral, and prove R(t) ≥ e^{−(α/β)t} from the convexity of e^{−λt}.
(d) Compute the Fisher information of the Poisson likelihood for λ and show that Jeffreys' prior is proportional to λ^{−1/2}. Show that it is invariant: applying the same rule to the MTBF θ = 1/λ gives a prior on θ that transforms to λ^{−1/2} on λ.
(e) Show that if λ ~ Gamma(k + 1, T) then 2Tλ ~ χ²_{2k+2}, and hence that the flat-prior 95 % upper credible bound is (8).

Answer form: five short derivations.

## 3. (Unfamiliar scenario: the normal–normal update is a Kalman gain)

A gyro's bias is believed, before launch, to be b ~ N(μ₀, σ₀²) with μ₀ = 0.4 °/h and σ₀ = 0.5 + s/100 °/h. On orbit a star-tracker comparison measures the bias as z = 1.2 °/h with Gaussian error of SD σ_z = 0.3 °/h.

(a) Show that the posterior is N(μ₁, σ₁²) with 1/σ₁² = 1/σ₀² + 1/σ_z² (precisions add) and μ₁ = μ₀ + K(z − μ₀) with K = σ₀²/(σ₀² + σ_z²). (This K is the scalar Kalman gain; Phase 07 builds the vector version.)
(b) Compute K, μ₁ and σ₁.
(c) Check with the code: put the normal log-densities on a grid, run `grid_posterior`, and compare `grid_mean`, the grid SD and the 90 % `credible_interval` with the closed form.
(d) What happens to K as σ_z → 0 and as σ_z → ∞, and what does each limit say about the measurement?

Answer form: a derivation, three numbers with units, a code check, two limits.

## 4. (Numerical, uses your code: your own fleet)

The new wheel's fleet has accumulated T = 40 000 + 2000s unit-hours with one failure. Use the lesson's heritage prior Gamma(2.5, 1.2 × 10⁶ h), the mission of 26 298 h and the requirement R ≥ 0.95.

(a) Give the posterior (α, β), its mean and SD, its 90 % equal-tailed interval (from `scipy.stats.gamma.ppf`), P(λ ≤ λ_req), R(3 yr) and the additional failure-free hours needed to reach P = 0.9.
(b) Compute the grid posterior on [0, 15 × 10⁻⁶] with 6001 points and report its mean, 90 % equal-tailed and HPD intervals.
(c) Run `metropolis` with 40 000 draws, step 1.5 × posterior SD and seed s; report the acceptance rate, the autocorrelation time, the effective sample size, the mean, SD and P(λ ≤ λ_req) with the scatter you expect from the effective sample size.
(d) Compare the frequentist 95 % upper bound from the fleet data alone with the heritage posterior's 95 % quantile, and explain the gap in one sentence.

Answer form: a table with scatter where appropriate, one sentence.

## 5. (Judgment: priors and decisions)

(a) Under what conditions is it legitimate to use the predecessor's 1.2 million hours as a prior for the new wheel? Name one change to the design that would preserve the legitimacy and one that would destroy it, and say which parts of the failure-rate the heritage can and cannot speak for.
(b) Suppose the heritage had been 0 failures in 1.2 million hours instead of 2. Compute the heritage prior, the posterior after the lesson's fleet data, P(λ ≤ λ_req) and R(3 yr). Does the requirement become demonstrated?
(c) The heritage posterior gives P(λ ≤ λ_req) = 0.33. A manager reads this as "the wheel fails the requirement". Explain what the number does and does not say, and what the posterior predictive R(3 yr) = 0.93 adds.
(d) A posterior predictive check: under the heritage prior, what was the probability of seeing at least one failure in 60 000 fleet hours? Was the observed failure evidence against the heritage? What observation *would* be?

Answer form: four short paragraphs, each with a number or a named condition.

## 6. (Hard, uses your code: calibration under a wrong prior)

(a) Run `credible_interval_coverage` with the heritage prior Gamma(2.5, 1.2 × 10⁶) for T = 6 × 10⁴ h and T = 6 × 10⁵ h at level 0.9, and for T = 6 × 10⁴ at level 0.5, with 4000 trials and seed s. Compare each with its level and the binomial scatter.
(b) Write your own loop in which the *true* rates are drawn from a Gamma with the same shape 2.5 but five times the heritage mean, while the update still uses the heritage prior. Measure the 90 % coverage for T = 6 × 10⁴ and T = 6 × 10⁵ h with seed s.
(c) Explain the two coverages in (b): why is it far below 0.9, and why does it recover only slowly with T?
(d) Propose a diagnostic that an engineer could run on the *real* fleet data to detect this situation without knowing the true rate, and compute it for the lesson's fleet (one failure in 60 000 h) under the heritage prior.

Answer form: five coverages with scatter, an explanation, a diagnostic with a number.
