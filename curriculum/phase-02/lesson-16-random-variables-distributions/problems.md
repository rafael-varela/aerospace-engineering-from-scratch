# P02.16  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. Give units on every quantity, say whether a number is exact, numerical (with its tolerance) or simulated (with its standard error), and label every simulated dataset as synthetic.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: two discrete variables by hand)

(a) A spacecraft has four identical thrusters, each working on command with probability p = 0.90 + s/1000. Write the pmf of the number X that work, and compute E[X] and Var X two ways: from the pmf by (1) and from the binomial formulas (2). Give P(X ≥ 3), the probability that the attitude manoeuvre (which needs three) succeeds.
(b) Conjunction alerts arrive as a Poisson process at λ = 1 + s/10 per month. Give P(no alert in a month), P(at least 3 in a month), and the probability of no alert in a *quarter*.
(c) Confirm every number with `binomial_pmf`, `poisson_pmf`, `discrete_mean` and `discrete_variance`.

Answer form: a pmf table, two means, two variances, four probabilities.

## 2. (Derivation: moments and the ellipse)

(a) Prove Var(aX + b) = a² Var X from the definition, and Var(X + Y) = Var X + Var Y + 2 Cov(X, Y).
(b) Derive E[X] and Var X for the uniform distribution on [a, b] by integration.
(c) Derive E[T] = 1/λ and Var T = 1/λ² for the exponential distribution (integration by parts), and prove the memoryless property.
(d) Derive P(inside the k-sigma ellipse) = 1 − e^{−k²/2} for a two-dimensional Gaussian, stating where you use independence of the principal components and where the polar-coordinate change of variables.

Answer form: four derivations.

## 3. (Unfamiliar scenario: micrometeoroid impacts on a solar array)

A solar array of area A = 20 + s m² flies for three years (1095 days) through a flux of 2 × 10⁻⁵ impacts per m² per day by particles large enough to matter (synthetic figure). Impacts are independent and the flux is constant.

(a) What is the expected number of impacts over the mission, and what distribution does the count follow? (b) Give P(no impact), P(at least 3 impacts) and the most probable count. (c) The waiting time to the first impact is exponential: give its mean in days and the probability that the first impact occurs within the first year. (d) The engineers observe 9 impacts in the three years. Using Var X = E[X] for a Poisson count, is 9 surprising? Give a z-score and a one-line judgement. (e) Which assumption would you question first if the count were 30?

Answer form: one expectation with its distribution, three probabilities, one mean and one probability, a z-score with a sentence, one assumption.

## 4. (Numerical, uses your code: a windier site)

The payload of the lesson is flown from a site where the wind-speed uncertainty is σ_w = 3 + s/10 m/s instead of 3 m/s; the other inputs are unchanged.

(a) Build the new Σ, propagate it with `propagate_covariance` and give σ_x, σ_y and ρ. (b) Give the 99 % ellipse (semi-axes, orientation, area) with `dispersion_ellipse`, and the CEP and 99 % circle radius with `radius_for_probability`. By what factor is the circle's area larger than the ellipse's? (c) Simulate 20 000 landings with seed s and report the fraction inside the 99 % ellipse with its standard error. (d) By how many degrees did the major axis rotate relative to the lesson's 14.8°, and why does more wind rotate it in that direction?

Answer form: three statistics, an ellipse, two radii, an area ratio, a simulated fraction with standard error, an angle with a one-sentence explanation.

## 5. (Numerical, uses your code: how fast is the central limit theorem?)

(a) With `clt_distance_to_normal(n, 200_000, seed=s)`, tabulate the maximum cdf gap for n = 1, 2, 4, 8, 16, 32 uniform terms and fit a power law gap ∝ n^(−α) to the points with n ≥ 2. Report α. (b) The empirical cdf of 200 000 draws has a sampling scatter of about 1/√200 000 = 0.0022 in its largest deviation; which of your points are at that floor, and what would you change to measure the rate further? (c) Repeat (a) for the sum of n *exponential* variables (write a four-line variant of `clt_sum_of_uniforms` with `rng.exponential`, standardised by mean n and standard deviation √n). Is the approach to normality faster or slower than for uniforms, and why?

Answer form: two tables, two exponents, a one-paragraph comparison.

## 6. (Judgment: what the ellipse does not know)

(a) The range-safety limit is a road 1200 m down-range of the aim point (cross-range irrelevant). Using the lesson's σ_x = 469.6 m, give the probability of landing beyond the road under the Gaussian model. Then suppose the wind error is actually Student-t with 4 degrees of freedom scaled to the same 3 m/s standard deviation; without computing, say whether the true probability is higher or lower and roughly by how much (a factor), and justify with one property of the t distribution.
(b) The model assumes the heading error and the wind error are independent. If the heading is estimated from the same wind forecast, with a correlation of +0.5 between the two errors, which element of Σ changes, and does the ellipse's major axis rotate towards or away from cross-range? Verify numerically with `propagate_covariance` and `dispersion_ellipse`.
(c) A reviewer proposes to publish the 99.9 % ellipse (k = 3.717) as the "guaranteed" zone. Write three sentences on why the 99.9 % figure is less trustworthy than the 99 % one, and what evidence would be needed to support it.

Answer form: one probability and one factor with justification, one element and one direction with a numerical check, three sentences.
