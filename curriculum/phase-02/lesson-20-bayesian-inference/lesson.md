# P02.20  Bayesian Inference

> A prior is what you knew, a likelihood is what the data say, and the posterior is the only honest answer to "how reliable is it now?"

| | |
|---|---|
| **Learning outcomes** | 1. Derive the Gamma–Poisson and Beta–binomial conjugate updates, the posterior predictive reliability (β/(β + t))^α, Jeffreys' prior for a rate, and the equivalence of the flat-prior credible bound with the frequentist χ² bound. 2. Implement and test, from scratch, a grid posterior and a random-walk Metropolis sampler, with log-space likelihoods that survive exposures of 10⁶ hours, against `scipy.stats` as the oracle. 3. Measure the calibration of 90 % credible intervals on seeded synthetic fleets and the effective sample size of a Metropolis chain. 4. Decide whether a new reaction wheel is shown to meet a 95 % three-year reliability requirement, how many more failure-free hours would demonstrate it, and how much the decision depends on the prior. |
| **Prerequisites** | `P02.15` (Bayes' theorem, conditional probability), `P02.16` (Poisson, exponential, binomial and Gamma distributions, expectation), `P02.17` (likelihood, confidence intervals, the χ² distribution), `P02.18` (least squares as maximum likelihood), `P01.04` (improper integrals, the Gamma function), `P00.03`, `P00.07`. Skill check: you can evaluate ∫₀^∞ x^{a−1}e^{−bx}dx = Γ(a)/b^a and explain why 0.95^n → 0. |
| **Engineering question** | A redesigned reaction wheel has run 60 000 unit-hours with one failure, and its predecessor ran 1.2 million unit-hours with two: what is the failure rate now, is the mission requirement P(no failure in 3 years) ≥ 0.95 demonstrated with 90 % confidence, and if not, how many more failure-free hours would do it? |
| **Assumptions and limits** | Failures are a Poisson process with a constant rate λ (no infant mortality, no wear-out, no common-cause failures, units independent); the heritage data describe the same failure mechanisms as the new design (the prior's legitimacy is a judgment, examined in part 5); the Gamma family is rich enough to hold what is known; every count is a lesson parameter or a seeded synthetic draw, not a measurement. The model says nothing about *why* the wheel failed, and a Bayesian interval is calibrated only on average over the prior: if the prior is wrong, so is the coverage, and part 4 measures by how much. |
| **Analytical work** | Bayes' theorem in probability and odds forms; a discrete two-hypothesis update; the Poisson likelihood and its log; conjugacy of the Gamma prior; the posterior mean as a precision-weighted average; the posterior predictive by the moment-generating function and Jensen's inequality; Jeffreys' prior from the Fisher information; the χ²_{2k+2} bound; the zero-failure demonstration 1 − (1 − p)^{n+1}; the normal–normal update as the scalar Kalman gain (problem 3). |
| **Computational deliverable** | `starter/main.py`: `bayes_discrete`, `odds_update`, `log_binomial_likelihood`, `log_poisson_likelihood`, `log_gamma_pdf`, `log_beta_pdf`, `gamma_poisson_update`, `beta_binomial_update`, `sequential_update`, `grid_posterior`, `grid_mean`, `grid_cdf`, `grid_quantile`, `credible_interval`, `hpd_interval`, `metropolis`, `autocorrelation_time`, `prob_rate_below`, `posterior_predictive_reliability`, `hours_to_demonstrate`, `frequentist_upper_bound`, `zero_failure_tests_required`, `beta_cdf_grid`, `simulate_fleet_failures`, `credible_interval_coverage`. Tests prove: the discrete posterior matches a hand calculation and composes over independent data; log-likelihoods and log-densities match `scipy.stats` to 10⁻¹² including the edges and integrate to one; the grid posterior reproduces the conjugate Gamma pointwise, in mean and in quantiles; the posterior mean is the precision-weighted average and the SD falls as 1/√T; sequential updates are order-independent; the HPD interval is shorter than the equal-tailed one and holds its mass; Metropolis reproduces the mean, SD and a tail probability within tolerances derived from its effective sample size; the predictive matches Monte Carlo and Jensen; `hours_to_demonstrate` closes the loop; the flat-prior bound equals χ²_{2k+2}/(2T); the zero-failure count is 58; and 90 % intervals cover the true rate in 90 % ± 2 % of 4000 synthetic fleets. |
| **Practical activity** | **computational experiment**: calibration of credible intervals on seeded synthetic fleets drawn from the prior (and, in problem 6, from a wrong prior), and a Metropolis chain's effective sample size against its step. All data are synthetic. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `failure_rate_posterior.png`, the hand derivations of the conjugate update and the predictive, the prior-sensitivity table, and a short note with the reliability decision, the hours needed, and the assumption you trust least. |
| **Assessment** | `assessments/phase-02/grading-rubric.md`; competencies 1–4. Pass: the conjugate and predictive derivations are complete; all public tests pass; the hand values for the heritage posterior agree with the code to three figures; the note states a decision with a probability, the hours that would change it, and the sensitivity to the prior. |
| **References** | `@mit-ocw-18-05-introduction-to-probability-and-statistics` Reading 3 "Conditional Probability, Independence and Bayes' Theorem", Reading 11 "Bayesian Updating with Discrete Priors", Reading 12a "Bayesian Updating: Probabilistic Prediction", Reading 12b "Bayesian Updating: Odds", Reading 13a "Bayesian Updating with Continuous Priors", Reading 14 "Continuous Data with Continuous Priors", Reading 15 "Conjugate priors: Beta and normal", Reading 16a "Choosing priors", Reading 16b "Probability intervals", Reading 20 "Comparison of frequentist and Bayesian inference" (required-free). Reading 10a "Introduction to Statistics", Reading 10b "Maximum Likelihood Estimates", Reading 13b "Notational conventions" (supplementary-free). The Gamma–Poisson pair, Jeffreys' prior and the reliability framing are not in these readings; part 2 derives them. |
| **Extensions** | ★ Replace the constant rate by a Weibull hazard with shape β and put a prior on both parameters; sample the joint posterior with Metropolis and show how a wear-out shape changes the 3-year prediction. ★ Treat the heritage prior's weight as uncertain (a hierarchical model: the new rate is the old rate times an unknown factor with its own prior) and show how the answer moves between the Jeffreys and heritage columns. ★ Implement the sequential probability ratio test and compare the hours it needs with `hours_to_demonstrate`. |
| **Time** | unmeasured |

## 1. The engineering problem

A spacecraft needs three reaction wheels working for three years. The wheel is a redesign: same bearings and lubricant as the predecessor, a new motor controller. The predecessor flew 1.2 million unit-hours and failed twice. The new wheel has accumulated 60 000 unit-hours across ten units in life test and on orbit, and one of them has failed. The customer's requirement is that the probability of a wheel surviving the mission, 26 298 h, is at least 0.95, and they want that demonstrated with 90 % confidence.

Count the ways a frequentist analysis struggles with this. The new data alone give a maximum-likelihood rate of 1/60 000 = 16.7 per 10⁶ h, eight times too high for the requirement (which needs λ ≤ −ln 0.95/26 298 = 1.95 per 10⁶ h), but one failure is one failure: the 95 % upper confidence bound is 79 per 10⁶ h and the lower bound is essentially zero. The heritage data are far more informative, 1.2 million hours, but they are not the new wheel's data, and a frequentist test has no principled way to weigh them. And the requirement is a statement about a *probability of a probability*, P(reliability ≥ 0.95) ≥ 0.9, which a confidence interval does not express.

Bayesian inference does all three things in one formula. The failure rate λ is treated as uncertain and given a **prior** distribution that encodes the heritage; the new data enter through the **likelihood**; Bayes' theorem combines them into the **posterior**, a distribution over λ from which "P(λ ≤ requirement)", "P(no failure in 3 years)" and "how many more hours" are all direct computations. The price is that the prior is a choice and the answer depends on it: this lesson computes the answer three times (a non-informative prior, a flat prior, and the heritage prior) and treats the spread as part of the result.

Success criteria, written first: the from-scratch grid posterior agrees with the conjugate closed form to 10⁻⁴ in its mean and 2 × 10⁻³ in its quantiles; the Metropolis chain reproduces the posterior mean within 5 % and P(λ ≤ requirement) within 0.03 using its measured effective sample size; the flat-prior upper credible bound equals the frequentist χ² bound to 10⁻¹⁰; and 90 % credible intervals cover the true rate in 90 % ± 2 % of 4000 synthetic fleets drawn from the prior.

## 2. Understand and derive

### 2.1 Bayes' theorem, in probabilities and in odds

For hypotheses H_i that are mutually exclusive and exhaustive, and data D (`@mit-ocw-18-05-introduction-to-probability-and-statistics` Reading 3, Reading 11),

**P(H_i | D) = P(D | H_i) P(H_i) / P(D)**,   P(D) = Σ_j P(D | H_j) P(H_j).   (1)

P(H_i) is the prior, P(D | H_i) the likelihood of the data under H_i, P(D) the **evidence** (the probability of the data averaged over the hypotheses), and P(H_i | D) the posterior. A worked table: a lot of pyrotechnic initiators is either "good" (per-unit failure probability p = 0.01) or "bad" (p = 0.10), with prior probabilities 0.9 and 0.1 from the supplier's history; 20 units are fired and one fails. The binomial likelihoods are 20 × 0.01 × 0.99¹⁹ = 0.1652 and 20 × 0.10 × 0.90¹⁹ = 0.2702. The evidence is 0.9 × 0.1652 + 0.1 × 0.2702 = 0.1757 and the posterior probability that the lot is good is 0.9 × 0.1652/0.1757 = **0.846**. One failure in twenty moved the belief from 90 % to 85 %: the data were more likely under "bad", but not by enough to overturn a 9:1 prior.

The **odds form** (Reading 12b) makes that arithmetic transparent. With odds O = P(H₁)/P(H₂),

**posterior odds = prior odds × Bayes factor**,   Bayes factor = P(D | H₁)/P(D | H₂).   (2)

Here 9 × (0.1652/0.2702) = 9 × 0.612 = 5.50, and 5.50/(1 + 5.50) = 0.846 again. The Bayes factor is the data's whole contribution; two independent data sets multiply their Bayes factors, which is why `test_discrete_update_matches_hand_calculation_and_sums_to_one` checks that updating twice equals updating once with the product of likelihoods.

### 2.2 A continuous unknown: the failure rate

When the unknown is a continuous quantity such as λ, the sum in (1) becomes an integral and the probabilities become densities (Reading 13a, 14):

**p(λ | D) = p(D | λ) p(λ) / ∫ p(D | λ′) p(λ′) dλ′**.   (3)

For a constant failure rate λ (units: 1/h) the number of failures in exposure T hours is Poisson with mean λT (P02.16), so the likelihood of k failures is

p(k | λ, T) = (λT)^k e^{−λT}/k!,   log p = k ln(λT) − λT − ln k!.   (4)

Two things about (4) matter in practice. It is a function of λ for fixed data, not a probability distribution over λ; its integral over λ is not 1. And its size is absurd in floating point: for T = 1.2 × 10⁶ h and λ = 2 × 10⁻⁶ the factor e^{−λT} is e^{−2.4}, harmless, but for λ = 10⁻⁴ it is e^{−120} ≈ 10⁻⁵², and the normalising integral in (3) sums terms over 50 decades. Every density in the code is therefore computed as a log and exponentiated only after subtracting the maximum (`grid_posterior`); the tests put λ = 10⁻⁹ and T = 10⁶ into the same function to make sure.

### 2.3 The conjugate prior and the update rule

Choose the prior from the **Gamma** family in its rate parametrisation,

p(λ) = Gamma(λ; α, β) = β^α λ^{α−1} e^{−βλ}/Γ(α),   mean α/β, variance α/β², mode (α − 1)/β.   (5)

Multiply by (4): p(λ | k, T) ∝ λ^{α−1}e^{−βλ} · λ^k e^{−λT} = λ^{α+k−1} e^{−(β+T)λ}. That is again a Gamma density, so the normalising integral is known without integrating anything:

**p(λ | k, T) = Gamma(λ; α + k, β + T)**.   (6)

A family that reproduces itself under the likelihood is called **conjugate** (Reading 15 does this for the Beta–binomial and normal–normal pairs; the Gamma–Poisson pair works the same way and is derived here because reliability needs it). The parameters have a plain reading: α is a count of prior failures and β a count of prior hours, so the update (6) says *add the failures to α and the hours to β*. The Beta–binomial pair for a per-demand failure probability p reads the same way: Beta(a, b) with a prior failures and b prior successes becomes Beta(a + k, b + n − k) after k failures in n demands.

Three consequences, each a test. **Order does not matter**: updating with (k₁, T₁) then (k₂, T₂) gives α + k₁ + k₂, β + T₁ + T₂, whichever comes first (`sequential_update`). **The posterior mean is a weighted average**:

(α + k)/(β + T) = [β/(β + T)] · (α/β) + [T/(β + T)] · (k/T),   (7)

the prior mean α/β and the maximum-likelihood estimate k/T, weighted by their hours. For the heritage prior β = 1.2 × 10⁶ and the fleet T = 6 × 10⁴, the prior carries 95 % of the weight: that is the arithmetic behind the engineering question. **Data eventually win**: as T → ∞ with k ≈ λT the mean tends to k/T and the SD √(α + k)/(β + T) falls as 1/√T; the prior's pull on the mean is O(α/T).

### 2.4 Which prior?

A prior is a modelling choice and must be stated (Reading 16a). Three are used here.

- **Heritage.** Start from nothing (the Jeffreys prior below) and update with the predecessor's 2 failures in 1.2 × 10⁶ h: Gamma(2.5, 1.2 × 10⁶), mean 2.08 per 10⁶ h, which says "the new wheel is the old wheel until shown otherwise". Legitimate if the failure mechanisms are shared; the new controller is the reason to doubt it.
- **Flat.** p(λ) = constant, Gamma(1, 0): improper (it does not integrate), but the posterior Gamma(1 + k, T) is proper as soon as there is exposure. It is not as innocent as it looks: a flat prior on λ puts most of its mass at large rates, and its posterior mean (k + 1)/T is twice the ML estimate when k = 1.
- **Jeffreys.** The prior proportional to the square root of the Fisher information, which for the Poisson likelihood is I(λ) = E[−∂² log p/∂λ²] = E[k]/λ² = T/λ, so p(λ) ∝ λ^{−1/2}, Gamma(1/2, 0). It is invariant to reparametrisation (the same rule applied to the MTBF 1/λ gives the same posterior for λ), and its posterior mean (k + ½)/T is the ML estimate plus half a failure.

With the fleet data (k = 1, T = 6 × 10⁴ h) and λ_req = 1.95 per 10⁶ h, the three posteriors are:

| prior | posterior (α, β) | mean (10⁻⁶/h) | 90 % interval | P(λ ≤ λ_req) | R(3 yr) | hours to 90 % |
|---|---|---|---|---|---|---|
| Jeffreys (½, 0) | (1.5, 6.0 × 10⁴) | 25.0 | [2.93, 65.1] | 0.028 | 0.580 | 1.54 × 10⁶ |
| flat (1, 0) | (2, 6.0 × 10⁴) | 33.3 | [5.92, 79.1] | 0.006 | 0.483 | 1.93 × 10⁶ |
| heritage (2.5, 1.2 × 10⁶) | (3.5, 1.26 × 10⁶) | 2.78 | [0.860, 5.58] | 0.330 | 0.930 | 1.82 × 10⁶ |

The heritage posterior has mean 2.78 and SD √3.5/1.26 × 10⁶ = 1.49 per 10⁶ h. Under no prior is the requirement demonstrated: even granting the heritage, the probability that λ is below the requirement is 0.33, not 0.9. Note that the heritage column needs *more* additional hours than Jeffreys to reach 90 %, because its posterior already sits near the requirement with a narrow spread, while the Jeffreys posterior, though centred far too high, is so wide that failure-free hours move it quickly. The spread between columns is the sensitivity of the decision to the prior, and it is the first thing a reviewer will ask for.

### 2.5 Credible intervals and the frequentist comparison

A **credible interval** (Reading 16b) is a range holding a stated posterior probability. The **equal-tailed** 90 % interval runs from the 5 % to the 95 % quantile; the **highest-posterior-density** (HPD) interval is the shortest range holding 90 %, found by lowering a density threshold until the mass above it is 0.9. For a right-skewed Gamma the HPD sits to the left: heritage posterior ETI [0.860, 5.58], HPD [0.538, 4.94] per 10⁶ h, 7 % shorter. Both are probability statements about λ given the data, which is what engineers usually *think* a confidence interval is.

What a confidence interval actually is (P02.17, Reading 20): a procedure that, over repeated samples at fixed λ, traps λ 95 % of the time. For k Poisson events in exposure T the one-sided 95 % upper bound is

**λ_U = χ²_{2k+2}(0.95)/(2T)**,   (8)

because 2λT has a χ²_{2k+2} distribution at the bound. Now compare with the flat-prior posterior Gamma(k + 1, T): a Gamma(k + 1, T) variable times 2T is χ²_{2k+2}, so its 95 % quantile is exactly (8). The frequentist bound *is* the flat-prior Bayesian bound: 79.06 per 10⁶ h for the fleet data, by both routes, to 10⁻¹⁰ in the test. The two schools differ not in arithmetic but in what the number means and in whether other information (the heritage) may be used. When there is no other information, a flat or Jeffreys prior reproduces the classical answer; when there is, the Bayesian route is the only one that can use it, and it must say so.

### 2.6 What the mission needs: the posterior predictive

The requirement is not about λ but about the next 26 298 h. The probability of no failure, averaged over what we believe about λ, is the **posterior predictive** (Reading 12a):

R(t) = ∫₀^∞ e^{−λt} Gamma(λ; α, β) dλ = (β/(β + t))^α.   (9)

The integral is the Gamma moment-generating function at −t: β^α/Γ(α) ∫λ^{α−1}e^{−(β+t)λ}dλ = β^α/(β + t)^α. Because e^{−λt} is convex, Jensen's inequality gives R(t) ≥ e^{−(α/β)t}: the reliability at the *mean* rate understates the predictive, slightly (0.9296 against 0.9303 here), because the chance that the rate is low helps survival more than the chance that it is high hurts it. For the heritage posterior R(3 yr) = (1.26/1.286)^3.5 = **0.930**, below the required 0.95; for the Jeffreys posterior 0.58. The test checks (9) against a Monte Carlo average of e^{−λt} over 200 000 seeded Gamma draws.

**Test planning.** How many additional failure-free hours T′ make P(λ ≤ λ_req) = 0.9? The posterior after T′ is Gamma(α, β + T′) and P(λ ≤ λ_req) = P(α, (β + T′)λ_req), the regularised incomplete gamma function, which rises with T′. Inverting, (β + T′)λ_req = P⁻¹(α, 0.9), so T′ = P⁻¹(α, 0.9)/λ_req − β: 1.82 × 10⁶ h for the heritage posterior, which is 208 wheel-years. The requirement cannot be demonstrated by testing in any practical time; it can only be *inherited* from heritage, which makes the legitimacy of the prior the whole decision. That is not a weakness of the method: it is the method making the real issue visible.

**Zero-failure demonstration.** For a per-demand device, separation nuts fired n times with no failure, the flat prior Beta(1, 1) becomes Beta(1, n + 1) and P(p ≤ p_max) = 1 − (1 − p_max)^{n+1}. For p_max = 0.05 at 95 %: n + 1 ≥ ln 0.05/ln 0.95 = 58.4, so n = 58 firings (the classical rule gives 59, because it uses (1 − p_max)^n ≤ 0.05: the flat prior counts as one free success). Forty firings give 0.878. Jeffreys' Beta(½, ½) gives Beta(½, 40.5) and 0.958: the two non-informative priors differ by 8 % in the probability of meeting the requirement, which is the kind of sensitivity that decides whether a lot ships.

### 2.7 Where the model stops being valid

A constant rate means no wear-out and no infant mortality; a wheel's bearings wear out, so a three-year prediction from a 6000-hour-per-unit test assumes the hazard is flat beyond anything observed. Independence between units fails under a common-cause defect (a bad lubricant batch), in which case the fleet's ten units are closer to one unit. The heritage prior assumes shared mechanisms; a new motor controller adds mechanisms the heritage knows nothing about, and no choice of α and β can express "the old data might be irrelevant" (extension 2 can). Finally, calibration: Bayesian 90 % intervals cover the true value 90 % of the time *when the true values are drawn from the prior*. Part 4 shows this on synthetic fleets and problem 6 shows what happens when the prior is wrong by a factor of five: coverage falls to about 20 %. The method's honesty depends entirely on the honesty of the prior.

## 3. Build the computational model

Open `starter/main.py`. The dictionaries `HERITAGE`, `FLEET`, `PYRO` and the mission constants hold the lesson values; rates are stored in 1/h and printed per 10⁶ h.

### Step 1: Bayes' theorem on a table, and the odds form

```python
def bayes_discrete(prior, likelihood):
    joint = prior * likelihood
    evidence = joint.sum()
    return joint / evidence, evidence            # (1)

def odds_update(prior_odds, bayes_factor):
    return prior_odds * bayes_factor             # (2)
```

Validate: the prior sums to 1, nothing is negative, the evidence is not zero (data impossible under every hypothesis is a modelling error, not a result).

### Step 2: likelihoods and densities in log space

```python
def log_poisson_likelihood(k, T, lam):
    lam = np.asarray(lam, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ll = np.where(lam > 0.0, k * np.log(lam * T), 0.0 if k == 0 else -np.inf)
    return ll - lam * T - lgamma(k + 1)          # (4)

def log_gamma_pdf(lam, alpha, beta):
    ll = (alpha - 1) * log(lam) - beta * lam      # -inf for lam <= 0
    if beta > 0: ll += alpha * np.log(beta) - lgamma(alpha)   # beta = 0: improper, unnormalised
    return ll
```

`math.lgamma` is the log of the Gamma function; `log_beta_pdf` and `log_binomial_likelihood` follow the same pattern, with the edge cases 0 log 0 = 0 handled explicitly so that p = 0 and p = 1 return the right limits (the tests compare with `scipy.stats.binom.logpmf` at both). The conjugate updates are two-line functions returning the new (α, β) or (a, b).

### Step 3: the grid posterior

```python
def grid_posterior(grid, log_prior, log_likelihood):
    lp = log_prior + log_likelihood
    dens = np.where(np.isfinite(lp), np.exp(lp - lp[np.isfinite(lp)].max()), 0.0)
    return dens / np.trapezoid(dens, grid)        # (3), numerically

def grid_cdf(grid, pdf):
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(grid))])
    return cdf / cdf[-1]
```

`grid_quantile` interpolates the tabulated CDF; `credible_interval` takes the two tail quantiles; `hpd_interval` sorts the density values and lowers the threshold until the mass between the first and last grid point above it reaches the level. The grid is the general-purpose tool: it needs no conjugacy, only a prior and a likelihood you can evaluate, and for one parameter 6001 points cost nothing. Its limits are the grid's: it must cover the posterior (15 per 10⁶ h here, where the heritage posterior has 10⁻⁹ of its mass beyond), and it fails above three or four parameters, which is what Metropolis is for.

### Step 4: Metropolis

```python
def metropolis(log_target, x0, n_samples, step, rng, burn_in=1000):
    x, lx = x0, log_target(x0)
    innov = rng.normal(0.0, step, size=n_samples + burn_in)
    uni = np.log(rng.uniform(size=n_samples + burn_in))
    for i in range(n_samples + burn_in):
        xp = x + innov[i]
        lxp = log_target(xp)
        if np.isfinite(lxp) and uni[i] < lxp - lx:
            x, lx = xp, lxp
        if i >= burn_in: out[i - burn_in] = x
    return out, acceptance_rate
```

Propose a Gaussian step, accept with probability min(1, p(x′)/p(x)), otherwise stay (and *record the stay*: a common bug is to record only moves). Only ratios of the target appear, so the normalising integral in (3) is never needed; a proposal outside the support (λ < 0, log density −∞) is simply rejected. The chain's samples are correlated; `autocorrelation_time` returns 1 + 2Σρ_k and the **effective sample size** n/τ is what sets the statistical tolerance. With step 1.5σ the heritage chain accepts 55 % of proposals, has τ = 6.5, and 40 000 draws are worth 6100 independent ones: mean 2.770 (closed form 2.778), SD 1.464 (1.485), P(λ ≤ λ_req) 0.330 (0.330). The test also runs the chain with steps 100 times smaller and larger to show acceptance near 1 and near 0, both of which mix badly.

### Step 5: the mission quantities

`prob_rate_below` is `scipy.special.gammainc(alpha, beta * lam)`; `posterior_predictive_reliability` is (9); `hours_to_demonstrate` inverts it with `gammaincinv`; `frequentist_upper_bound` is (8) with `scipy.stats.chi2.ppf`; `zero_failure_tests_required` is ⌈ln(1 − c)/ln(1 − p_max) − 1⌉; `beta_cdf_grid` integrates `log_beta_pdf` on a grid with the substitution p = sin²(πu/2) so that Jeffreys' p^{−1/2} singularity is tamed. `simulate_fleet_failures` is one seeded Poisson draw and `credible_interval_coverage` runs the calibration experiment of part 4.

### Verify

Run `python -m pytest curriculum/phase-02/lesson-20-bayesian-inference`. The 16 tests:

1. **Analytic cases.** The discrete posterior equals the hand calculation; equal likelihoods leave the prior alone; the odds form agrees; log-likelihoods and log-densities match `scipy.stats` to 10⁻¹² at p = 0, 1 and λ = 0, 10⁻⁹; each density integrates to one; the Jeffreys and flat improper priors have the right shape.
2. **Grid against the closed form.** The grid posterior equals Gamma(α + k, β + T) pointwise to 2 × 10⁻⁵; its mean to 10⁻⁴; quantiles and the equal-tailed interval to 2 × 10⁻³; the CDF to 2 × 10⁻⁵; improper priors give proper posteriors with mean (α₀ + k)/T.
3. **Properties.** (7) holds; the SD falls by √0.1 per decade of T once T ≫ β; sequential equals batch in any order; the HPD holds 90 % ± 0.2 %, is shorter than the ETI and contains the mode; P(λ ≤ x) is monotone with the right limits; the predictive is decreasing, tends to e^{−λt} for a known rate, and exceeds e^{−mean t}; `hours_to_demonstrate` returns 0 when the requirement is met and makes P(λ ≤ λ_req) = 0.9 otherwise; (8) equals the Gamma quantile; the zero-failure count is 58 and 298 for p_max = 0.01.
4. **Statistics** (tolerances derived in comments from the effective sample size): Metropolis mean within 5 % (7σ), SD within 10 % (11σ), tail probability within 0.03 (5σ), reproducible under the same seed; the predictive against 200 000 Monte Carlo draws within 10⁻³ (8σ); Poisson draws with mean 2 ± 0.1 (4.5σ); coverage 0.90 ± 0.02 and 0.50 ± 0.035 (4σ).
5. **Errors.** `ValueError` on a prior that does not sum to 1, impossible data, k > n, T = 0, α = 0, level = 1, a negative rate.

### Compare

`scipy.stats.gamma` is the oracle for everything the conjugate form gives in closed form: pdf, cdf, ppf and random draws, all of which the grid and the chain reproduce. Why write a grid and a sampler at all? Because conjugacy is a lucky accident of the Poisson–Gamma and binomial–Beta pairs. The moment the model changes, a Weibull hazard, a rate that differs between units, a prior that is a mixture of "heritage applies" and "heritage does not", there is no closed form, and the grid (for one or two parameters) or Metropolis (for more) is the only route. `scipy.special.gammainc` and `gammaincinv` are used directly for P(λ ≤ x) and its inverse because they are the same function the grid CDF approximates, computed to machine precision; the test confirms the grid CDF agrees with them to 10⁻⁴. In practice one reaches for a library sampler (an adaptive or Hamiltonian method) over a hand-written random walk once the parameters exceed a handful, for the same reason one reaches for `solve_ivp` over hand-written RK4: the algorithm is the same idea with twenty years of step-size wisdom attached.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Every fleet is synthetic: a rate drawn from a stated distribution and a seeded Poisson count. Nothing here is a flight record.

1. **Calibration under the prior.** For 4000 trials: draw λ from the heritage prior Gamma(2.5, 1.2 × 10⁶), draw k ~ Poisson(λ × 60 000), update, compute the 90 % equal-tailed interval, and record whether it contains λ. Prediction to write down first: the coverage is 0.90 within the binomial scatter √(0.9 × 0.1/4000) = 0.5 %. Repeat with level 0.5 and with T = 10⁶ h.
2. **Metropolis step size.** For steps of 0.1σ, 0.5σ, 1.5σ, 5σ and 20σ (σ the posterior SD), record the acceptance rate and the autocorrelation time of 40 000-draw chains. Prediction: acceptance falls monotonically with the step; τ is large at both ends and smallest near 1–5σ.
3. **Prior sensitivity.** Reproduce the table of 2.4 and add a fourth prior: the heritage prior with its weight halved (α = 1.25, β = 0.6 × 10⁶, same mean, twice the variance). Where does its column land?
4. **A wrong prior** (problem 6). Draw λ from a Gamma with five times the heritage mean, update with the heritage prior, and measure the coverage.

Raw data format: CSV `trial, lambda_true, k, lo, hi, covered` for item 1 and `step, acceptance, tau` for item 2; the figure `failure_rate_posterior.png` (prior, likelihood, grid posterior and Metropolis histogram; posterior CDFs for the three priors against the requirement; the zero-failure curve; the test-planning curves).

What you should find: coverage 0.902 for 4000 fleets at the 90 % level and 0.492 at the 50 % level (both within 1σ of the binomial scatter); acceptance 0.55 and τ = 6.5 at step 1.5σ, against τ ≈ 950 at 0.1σ (acceptance 0.96, the chain creeps) and τ ≈ 34 at 20σ (acceptance 0.06, the chain sticks); the half-weight heritage column lands between the heritage and Jeffreys columns (mean 3.41, P(λ ≤ λ_req) = 0.30, and *fewer* hours to demonstrate, 1.52 × 10⁶, because a wider posterior is moved more by each failure-free hour); and coverage near 0.2 when the prior is wrong by a factor of five, the number that belongs in the honesty box of every Bayesian report.

## 5. Validate and improve

| quantity (heritage prior, fleet data) | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| λ_req (per 10⁶ h) | −ln 0.95/26 298 = 1.950 | 1.950 | 0 | exact |
| P(lot good | 1 of 20) | 0.846 | 0.8463 | 0 | exact |
| posterior (α, β) | (3.5, 1.26 × 10⁶) | same | 0 | (6) |
| posterior mean, SD (per 10⁶ h) | 2.778, 1.485 | grid 2.7777; chain 2.770, 1.464 | 3 × 10⁻⁵; 0.3 %, 1.4 % | trapezoid; 6100 effective draws (0.7 %, 0.9 % scatter) |
| 90 % interval | [0.860, 5.58] (Gamma quantiles) | grid [0.860, 5.582] | < 10⁻³ | interpolated CDF |
| P(λ ≤ λ_req) | 0.330 | grid 0.330, chain 0.330 | 0 | — |
| R(3 yr) | (1.26/1.286)^3.5 = 0.930 | 0.9303; Monte Carlo 0.9303 | 0 | (9) |
| e^{−mean t} | 0.9296 | 0.9296 | below R by 0.07 % | Jensen |
| frequentist 95 % bound (fleet only) | χ²₄(0.95)/(1.2 × 10⁵) = 79.06 | 79.064 = flat-prior quantile | 0 | (8) |
| zero-failure firings for 0.05 at 95 % | 58 | 58 | 0 | 1 − 0.95^{59} = 0.951 |
| coverage of 90 % intervals, 4000 fleets | 0.90 ± 0.005 | 0.902 | 0.4σ | binomial scatter |

Every difference has a cause and a size. The chain's 1.4 % error in the SD is a 1.5σ draw for 6100 effective samples, and is the reason the test tolerances are quoted in multiples of the effective-sample scatter rather than as round numbers: a chain that reported 40 000 samples' worth of precision would be lying by a factor √6.5.

What the agreement does not show: that the posterior is right about the wheel. It shows that three implementations of the same model agree. The model's assumptions, a constant rate, independence, and above all the heritage prior's claim that 1.2 million hours of the old design speak for the new one, are not tested by any of this. The sensitivity table is the only evidence on that point, and it says the decision flips from "0.33" to "0.03" depending on the prior: the honest report states both and explains what the new controller changes.

Improvements worth making: a mixture prior, p(λ) = w × heritage + (1 − w) × Jeffreys, with w the engineering judgment of how much the designs share, which puts the sensitivity inside the model instead of beside it; a Weibull hazard (extension 1) to test the constant-rate assumption against the heritage failure times, if they are known; a posterior predictive check: under the heritage prior the probability of seeing at least one failure in the fleet's 60 000 h was 1 − (1.2/1.26)^2.5 = 0.115, so the observed failure was unlucky but not surprising, and the model is not contradicted. Improvement worth *not* making: a finer grid or a longer chain, which refine the third decimal of a number whose first decimal depends on the prior.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivations of (6), (7), (9) and the Jeffreys prior; `failure_rate_posterior.png`; the calibration CSV; and a one-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the posterior rate with its 90 % interval under each prior, R(3 yr), the hours needed to demonstrate the requirement, the calibration result, and the assumption you trust least.

### Adversarial self-checks

- **What if the exposure doubles with no new failure?** β → 1.32 × 10⁶, mean 2.65, P(λ ≤ λ_req) rises from 0.33 to 0.36 only: 60 000 hours are 5 % of the prior's weight. What if the fleet had seen two failures? α → 4.5, mean 3.57, P → 0.16. Check both with the code.
- **What if T → 0? T → ∞?** The posterior is the prior; the posterior collapses on k/T and the prior is forgotten, at a rate √(α + k)/(β + T).
- **Which assumption fails first in the real world?** The transferability of the heritage: the new controller is not in the 1.2 million hours. Second, the constant rate: bearings wear.
- **What measurement could falsify the model?** Failure *times*, not just counts. A constant rate predicts exponential inter-failure times; a cluster of failures late in life (or early) rejects it, and the Weibull extension quantifies by how much.
- **Does the result survive a change in grid, chain length, seed or prior?** Grid and chain agree with the closed form within their stated scatter for any seed; the *prior* changes the answer by a factor of ten in P(λ ≤ λ_req). State which of those is a numerical issue and which is the engineering.
- **Which decision would change if the uncertainty doubled?** If the heritage were worth half its hours (α = 1.25, β = 0.6 × 10⁶), the probability of meeting the requirement would fall to 0.30 and the hours needed to 1.52 × 10⁶, still 174 wheel-years; the decision "accept on heritage" would not survive, and the honest alternative is a redundancy or a design change rather than more testing, because 208 wheel-years of testing is not available at any uncertainty.

## Check yourself

1. State Bayes' theorem for a continuous parameter and name the four pieces; which one is a density over the parameter and which is not?
2. What is a conjugate prior, and what do α and β of a Gamma prior count?
3. Write the Gamma–Poisson update and the posterior mean as a weighted average; what are the weights for β = 1.2 × 10⁶ h and T = 6 × 10⁴ h?
4. Why does the code evaluate every likelihood as a logarithm?
5. Give the posterior predictive probability of no failure in t hours and say which direction Jensen's inequality points.
6. In what sense is the frequentist 95 % upper bound on a Poisson rate a Bayesian answer, and with which prior?
7. A Metropolis chain of 40 000 draws has autocorrelation time 6.5. How many independent samples is it worth, and what relative scatter should you attach to its estimate of a probability near 0.33?
8. What does "calibrated" mean for a credible interval, and under what condition does the calibration hold?
