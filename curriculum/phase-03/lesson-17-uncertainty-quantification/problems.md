# P03.17  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer, give every statistical result with its standard error or interval, and say which input distributions you assumed.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: linear propagation by hand)

A rocket stage has specific impulse I_sp = (300 + s) s known to ±1 % and mass ratio m₀/m_f = 3.5 known to ±2 % (one standard deviation each, independent, Gaussian). The ideal velocity increment is Δv = g₀ I_sp ln(m₀/m_f).

(a) Compute the nominal Δv and, by first-order propagation, its standard deviation and the share of variance from each input.
(b) Using a Gaussian approximation, estimate the probability that Δv falls below 3700 m/s.
(c) Which of the two inputs is worth measuring better, and by how much would halving its uncertainty reduce the standard deviation of Δv?

Answer form: Δv ± u in m/s, two variance shares, one probability, one sentence with a number.

## 2. (Derivation: errors, orthogonality, coefficients)

(a) Show that the Monte Carlo mean estimator is unbiased and has variance σ²/N, and that the estimator of a probability p has variance p(1 − p)/N. Hence derive the sample size needed to know a small p to a relative accuracy ε.
(b) Using the generating function e^{ξt − t²/2} = Σ_n He_n(ξ) tⁿ/n!, show that E[He_m(ξ) He_n(ξ)] = n! δ_mn. (Hint: E[e^{ξ(s + t)}] = e^{(s + t)²/2}.)
(c) From the same generating function, derive the PCE coefficients of f(ξ) = e^{aξ}: c_n = e^{a²/2} aⁿ/n!. Check that Σ c_n² n! equals E[e^{2aξ}] = e^{2a²}.
(d) Starting from (5), derive the PCE mean and variance formulas and show that for an additive model (every α with at most one nonzero entry) the first-order Sobol indices sum to exactly 1.

Answer form: four short derivations; (c) ends with an identity you verify numerically.

## 3. (Unfamiliar scenario: landing distance)

An aircraft's ground roll is s_L = V²/(2a), with the touchdown speed V = 60 + s/5 m/s uncertain by ±3 % and the mean deceleration a = 2.5 m/s² uncertain by ±8 % (Gaussian, independent). The available runway is 900 m.

(a) Compute the nominal ground roll and the linear standard deviation and variance shares by hand, and the second-order mean shift (both terms).
(b) Build f(ξ) in the lesson's standard-normal form and run `linear_propagation`, `pce_projection` at degrees 1–4 and `monte_carlo` with 2 × 10⁵ samples. Report mean and standard deviation from each.
(c) Compute P(s_L > 900 m) from the Monte Carlo samples (with its standard error), from the degree-4 surrogate, and from the Gaussian approximation. Explain the disagreement with the skewness of the samples.
(d) Give the Sobol indices and say which input a pilot can influence and which the runway designer must allow for.

Answer form: a table (method, mean, std), three probabilities with errors, two Sobol indices, a short paragraph.

## 4. (Numerical, uses your code: your own drag-polar uncertainty)

Repeat the lesson's range propagation with the lift-to-drag relative standard deviation set to 0.02 + s/1000 and the required range to 1400 km, everything else unchanged.

(a) Report the linear mean, std and variance shares; the degree-3 PCE mean, std and Sobol indices; and the Monte Carlo mean, std and shortfall probability at N = 10³, 10⁴, 10⁵ with standard errors.
(b) Compare the shortfall probability from the surrogate, from direct Monte Carlo and from the Gaussian approximation.
(c) How many direct Monte Carlo samples would you need to know the shortfall probability to 10 % relative? How many model runs did the PCE need?
(d) Write the two-line uncertainty statement of part 5 of the lesson for your case.

Answer form: three tables, one count, two lines.

## 5. (Numerical, uses your code: variance reduction and the bootstrap)

For the lesson's inputs:

(a) Compare the standard error of the mean from `monte_carlo` and `antithetic_monte_carlo` at N = 10³ and 10⁴, and compute the equivalent plain sample size that antithetic sampling buys.
(b) Estimate the 5th percentile of the range from 10⁴ Monte Carlo samples with a bootstrap 95 % interval (2000 resamples), and compare with the 5th percentile of the degree-3 surrogate sampled 4 × 10⁵ times and with the Gaussian approximation. Does the interval cover the surrogate's value? Should it always?
(c) A colleague says "antithetic sampling made the mean 350 times more accurate, so I will use it for the 5th percentile too". Test the claim and explain the result.

Answer form: a table of standard errors, three percentiles with an interval, one short paragraph.

## 6. (Judgment: certifying a probability)

A regulator asks for the probability that a re-entry capsule's peak heat flux exceeds the heat-shield allowable, computed from a model with six uncertain inputs (entry angle, velocity, atmospheric density scale height, ballistic coefficient, two material properties) whose distributions were elicited from experts. The requirement is a probability below 10⁻⁴ with 90 % confidence.

(a) Which propagation method(s) would you use, and in what order? Give the model-run counts.
(b) Two of the inputs are known to be correlated, and one (the density scale height) is bounded physically. What does each do to the lesson's machinery?
(c) The PCE surrogate's tail disagrees with 200 direct runs placed in the tail by a factor 2. Which do you believe, and what do you do?
(d) The computed probability is 6 × 10⁻⁵. Write the one-paragraph statement you would give the regulator, including what it is conditional on and what evidence would change it.

Answer form: four short paragraphs, each with at least one number.
