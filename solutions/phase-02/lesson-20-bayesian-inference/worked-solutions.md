# P02.20  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown. Every fleet is synthetic.

## 1. A telemetry alarm ($s = 15$: prior 0.015, false-positive rate 0.065)

**Derivation.** (a) Evidence P(alarm) = 0.95 × 0.015 + 0.065 × 0.985 = 0.01425 + 0.06403 = 0.07828. Posterior P(stuck | alarm) = 0.01425/0.07828 = **0.182**. Odds form: prior odds 0.015/0.985 = 0.01523; Bayes factor 0.95/0.065 = **14.6**; posterior odds 0.2226; probability 0.2226/1.2226 = 0.182 ✓.

(b) Second independent alarm: prior odds now 0.2226, times 14.6 again gives 3.254, probability 3.254/4.254 = **0.765**. (Equivalently prior odds × 14.6².)

(c) With a fault present 1.5 % of the time, 98.5 % of checks are on a healthy valve, and 6.5 % of those alarm: the false alarms (0.064 per check) outnumber the true ones (0.014) four to one, so one alarm is more likely false than true however good the detector is at catching real faults. The second independent alarm multiplies the odds by the Bayes factor again, and 14.6² = 214 overwhelms a 1:66 prior.

**Check.** Limiting case: a perfect detector (false-positive rate 0) gives P = 1 after one alarm; a useless one (false-positive rate 0.95) gives a Bayes factor of 1 and the prior unchanged. Code: `bayes_discrete([0.015, 0.985], [0.95, 0.065])` returns 0.1821.

**Common mistakes.** Confusing P(alarm | stuck) = 0.95 with P(stuck | alarm); forgetting the base rate; treating the second alarm as a fresh update from the original prior (which gives 0.182 again instead of 0.765).

## 2. The Gamma–Poisson machinery

**Derivation.** (a) p(λ | k, T) ∝ [λ^k e^{−λT}T^k/k!] × [β^αλ^{α−1}e^{−βλ}/Γ(α)] ∝ λ^{α+k−1}e^{−(β+T)λ}. The normalising integral is ∫₀^∞λ^{α+k−1}e^{−(β+T)λ}dλ = Γ(α + k)/(β + T)^{α+k} by the Gamma-function identity (substitute u = (β + T)λ), so the posterior is (β + T)^{α+k}λ^{α+k−1}e^{−(β+T)λ}/Γ(α + k) = Gamma(α + k, β + T) ∎.

(b) (α + k)/(β + T) = α/(β + T) + k/(β + T) = [β/(β + T)](α/β) + [T/(β + T)](k/T) ∎. Weights for β = 1.2 × 10⁶, T = 6 × 10⁴: 0.952 on the prior mean and 0.048 on k/T.

(c) R(t) = ∫₀^∞e^{−λt}β^αλ^{α−1}e^{−βλ}/Γ(α) dλ = β^α/Γ(α) × Γ(α)/(β + t)^α = (β/(β + t))^α ∎. Jensen: f(λ) = e^{−λt} has f″ = t²e^{−λt} > 0, so for any distribution E[f(λ)] ≥ f(E[λ]) = e^{−(α/β)t} ∎.

(d) log p(k | λ) = k ln λ − λT + k ln T − ln k!. ∂/∂λ = k/λ − T; ∂²/∂λ² = −k/λ². Fisher information I(λ) = E[k/λ²] = λT/λ² = T/λ. Jeffreys' prior ∝ √I ∝ λ^{−1/2} (the constant T is irrelevant) ∎. Invariance: in θ = 1/λ the log-likelihood is −k ln θ − T/θ + const; ∂²/∂θ² = k/θ² − 2T/θ³; I(θ) = −E[…] = −(λT/θ² − 2T/θ³) = −T/θ³ + 2T/θ³ = T/θ³, so the Jeffreys prior on θ is ∝ θ^{−3/2}. Transforming a density: p_λ(λ) = p_θ(θ)|dθ/dλ| = θ^{−3/2}λ^{−2} = λ^{3/2}λ^{−2} = λ^{−1/2}, the same prior ∎. (A flat prior on θ would transform to λ^{−2} on λ: not the same as flat on λ, which is why "flat" is not a neutral choice.)

(e) Gamma(k + 1, T): p(λ) = T^{k+1}λ^k e^{−Tλ}/k!. Let x = 2Tλ, dλ = dx/(2T): p(x) = T^{k+1}(x/2T)^k e^{−x/2}/(k! × 2T) = x^k e^{−x/2}/(2^{k+1}k!). The χ²_ν density is x^{ν/2−1}e^{−x/2}/(2^{ν/2}Γ(ν/2)); with ν = 2k + 2, ν/2 − 1 = k, 2^{ν/2} = 2^{k+1}, Γ(k + 1) = k!: identical ∎. So P(λ ≤ λ_U) = 0.95 ⇔ P(χ²_{2k+2} ≤ 2Tλ_U) = 0.95 ⇔ λ_U = χ²_{2k+2}(0.95)/(2T), which is (8).

**Check.** Dimensions in (c): β and t are both hours, the ratio is dimensionless ✓. Limit of (c): α, β → ∞ at fixed α/β = λ₀ gives (1 + t/β)^{−α} → e^{−λ₀t} ✓ (the code's test). Code: `frequentist_upper_bound(1, 6e4, 0.95)` = 79.06 per 10⁶ h = `gamma.ppf(0.95, 2, scale=1/6e4)` ✓.

**Common mistakes.** Using the scale parametrisation of the Gamma (mean αθ) where the rate one is meant; forgetting that the likelihood's k!/T^k constants cancel but the prior's β^α does not when comparing *evidences*; in (d), differentiating with respect to θ but taking the expectation in λ.

## 3. The normal–normal update ($s = 15$: σ₀ = 0.65 °/h)

**Assumptions.** Gaussian prior and Gaussian measurement error, independent; a scalar bias constant over the measurement.

**Derivation.** (a) log p(b | z) = −(b − μ₀)²/(2σ₀²) − (z − b)²/(2σ_z²) + const. The b² terms: −b²/2 × (1/σ₀² + 1/σ_z²), so the posterior precision is 1/σ₁² = 1/σ₀² + 1/σ_z². The b terms: b(μ₀/σ₀² + z/σ_z²), so μ₁ = σ₁²(μ₀/σ₀² + z/σ_z²). Rearranging, μ₁ = μ₀ + K(z − μ₀) with K = σ₁²/σ_z² = σ₀²/(σ₀² + σ_z²), and σ₁² = (1 − K)σ₀² ∎.

(b) σ₀² = 0.4225, σ_z² = 0.09: K = 0.4225/0.5125 = **0.824**; μ₁ = 0.4 + 0.824 × 0.8 = **1.060 °/h**; 1/σ₁² = 2.367 + 11.11 = 13.48 (°/h)⁻², σ₁ = **0.272 °/h**.

(c) Grid on [−3, 5] with 80 001 points, log prior `norm.logpdf(g, 0.4, 0.65)`, log likelihood `norm.logpdf(1.2, g, 0.3)`: `grid_mean` = 1.05951, grid SD = 0.27239, 90 % interval [0.6115, 1.5076] against `norm.ppf([0.05, 0.95], 1.0595, 0.2724)` = [0.6115, 1.5076]: agreement to five figures.

(d) σ_z → 0: K → 1, μ₁ → z, σ₁ → 0: a perfect measurement replaces the prior. σ_z → ∞: K → 0, μ₁ → μ₀, σ₁ → σ₀: a worthless measurement leaves the prior alone. The gain is the fraction of the innovation z − μ₀ that is believed; it is also the weight of the measurement in the precision-weighted average, which is the normal–normal version of (7).

**Check.** Precision adds: 13.48 = 2.37 + 11.11 ✓. The posterior SD is below both input SDs (0.272 < 0.3 < 0.65), as it must be: combining two independent estimates never makes you less certain.

**Common mistakes.** Adding variances instead of precisions; writing K = σ_z²/(σ₀² + σ_z²) (the weight of the prior, not of the measurement); on the grid, using `norm.logpdf(grid, z, sz)` for the likelihood, which happens to be the same number here because the Gaussian is symmetric, but is the wrong habit.

## 4. Your own fleet ($s = 15$: T = 70 000 h, one failure)

**Result.** (a) Posterior Gamma(**3.5, 1.27 × 10⁶**): mean 2.756 per 10⁶ h, SD 1.473, 90 % interval [0.853, 5.538]; P(λ ≤ λ_req) = **0.334**; R(3 yr) = **0.9308**; hours to reach P = 0.9: **1.81 × 10⁶** additional failure-free hours. The fleet's own ML rate is 14.3 per 10⁶ h; the posterior mean is pulled to 2.76 by the 1.2 million heritage hours (weight 0.945).

(b) Grid (6001 points on [0, 15 × 10⁻⁶]): mean 2.7559 (closed form 2.7559), 90 % equal-tailed [0.853, 5.538], HPD [0.533, 4.90] (11 % shorter).

(c) Metropolis, seed 15, step 1.5 × 1.473 × 10⁻⁶: acceptance 0.55, autocorrelation time 6.7, effective samples 5934. Mean 2.780 ± 0.019 (scatter σ/√n_eff = 1.473/77 = 0.019; closed form 2.756, a 1.2σ draw), SD 1.479 ± 0.014 (σ/√(2n_eff); closed form 1.473, 0.4σ), P(λ ≤ λ_req) = 0.328 ± 0.006 (√(0.33 × 0.67/5934); closed form 0.334, 1σ).

(d) Frequentist 95 % upper bound from the fleet alone: χ²₄(0.95)/(2 × 70 000) = 9.488/140 000 = **67.8 per 10⁶ h**; the heritage posterior's 95 % quantile is 5.54. The factor of twelve between them is the 1.2 million heritage hours that the frequentist bound does not use.

**Check.** Weighted average: 0.945 × 2.083 + 0.055 × 14.29 = 1.968 + 0.786 = 2.754 ≈ 2.756 ✓ (rounding). Jensen: e^{−2.756 × 10⁻⁶ × 26 298} = 0.9301 < 0.9308 ✓.

**Common mistakes.** Quoting the chain's numbers to four figures without the scatter; computing hours-to-demonstrate from the Jeffreys posterior and attributing it to the heritage one; using `scipy.stats.gamma(a, scale=b)` with the rate in place of the scale (a factor β² error in the mean).

## 5. Priors and decisions

(a) The heritage prior is legitimate when the new wheel's failures arise from the same mechanisms under the same stresses as the old one's: same bearings, lubricant, preload, speed range, thermal environment and duty cycle. A change confined to the motor-controller firmware, with the same current and speed profile, preserves it (the controller is not what failed in the heritage, if its failures were bearing failures); a new bearing supplier or a higher operating speed destroys it, because the 1.2 million hours exercised neither. The heritage can speak for the mechanisms it exercised, and only for those; it cannot speak for a mechanism it never contained, and no value of (α, β) can say "I may be irrelevant". Mixture priors (extension 2) can.

(b) Heritage of 0 failures in 1.2 × 10⁶ h from the Jeffreys start: prior Gamma(0.5, 1.2 × 10⁶), mean 0.42 per 10⁶ h. After the fleet's 1 failure in 60 000 h: Gamma(1.5, 1.26 × 10⁶), mean **1.19 per 10⁶ h**, P(λ ≤ λ_req) = **0.822**, R(3 yr) = **0.969**. The predictive reliability now exceeds 0.95, but the rate requirement is demonstrated at 82 %, not 90 %: still short, and the single fleet failure now carries two thirds of the posterior's "failure count".

(c) P(λ ≤ λ_req) = 0.33 says: given the heritage, the fleet data and the model, there is a one-in-three chance that the wheel's constant failure rate meets the requirement, and a two-in-three chance that it does not. It is not a verdict that the wheel *fails*, and it is not the probability that the wheel fails in flight. That second quantity is 1 − R(3 yr) = 0.07: the wheel most probably survives the mission, because even a rate somewhat above the requirement gives a survival probability near 0.9. The requirement "R ≥ 0.95 at 90 % confidence" and the prediction "R = 0.93" are different statements; the first is about demonstrating a bound, the second about expecting an outcome, and a manager needs both.

(d) Under the heritage prior Gamma(2.5, 1.2 × 10⁶), the probability of no failure in 60 000 h is (1.2/1.26)^2.5 = 0.885, so P(at least one failure) = **0.115**. The observed failure had an 11 % chance under the model: unlucky, not surprising, and not evidence against the heritage by any usual threshold. Three or more failures in 60 000 h (prior predictive probability 0.07 %, from the negative-binomial prior predictive) would be; so would a failure of the new controller, whatever the count, because the heritage prior says nothing about that mechanism.

**Common mistakes.** Reading P(λ ≤ λ_req) as a reliability; treating "not demonstrated" as "fails"; using the posterior rather than the prior to judge whether the data were surprising (the posterior has already absorbed them).

## 6. Calibration under a wrong prior ($s = 15$)

**Result.** (a) Coverage of 90 % intervals, true rates from the heritage prior: **0.892** at T = 6 × 10⁴ h and **0.896** at T = 6 × 10⁵ h (binomial scatter √(0.9 × 0.1/4000) = 0.005: 1.6σ and 0.9σ); coverage of 50 % intervals at T = 6 × 10⁴: **0.492** (scatter 0.008, 1σ). Calibrated, as the theorem says they must be when the truth is drawn from the prior.

(b) True rates from Gamma(2.5, scale 5/1.2 × 10⁶), mean 10.4 per 10⁶ h, updates still with the heritage prior: coverage **0.198** at T = 6 × 10⁴ h and **0.289** at T = 6 × 10⁵ h.

(c) At T = 6 × 10⁴ the data carry 5 % of the weight; the posterior is essentially the prior, centred at 2 per 10⁶ h with a 90 % interval reaching to about 5.6, while the true rates average 10: only the lucky low draws from the true distribution are covered. At T = 6 × 10⁵ the data weight is 33 % and the intervals move toward the truth, but the prior still pulls the centre down by two thirds of the gap, so most intervals still miss. Coverage recovers as T/β → ∞, i.e. once the fleet has accumulated several million hours, which is the same statement as (7): a wrong prior is forgotten at the rate the data outweigh its hours.

(d) The prior predictive check of problem 5(d): compute, before trusting the posterior, the probability under the prior of a count at least as large as the one observed. For the lesson fleet, P(k ≥ 1 | heritage, 60 000 h) = 0.115: not alarming. Had the true rate been 10 per 10⁶ h the expected count would be 0.6 and a value of k ≥ 2 (prior predictive probability 0.9 % under the heritage, but 13 % at the true rate) would have appeared in one fleet in eight; k ≥ 2 is the signal to stop trusting the heritage. The engineer does not need the true rate to run this check, only the prior and the data, which is why it belongs in every Bayesian report.

**Check.** Limit: with T → ∞ in (b) the coverage must return to 0.9, since the posterior then ignores the prior; the trend 0.198 → 0.289 from 6 × 10⁴ to 6 × 10⁵ h is in that direction. Scatter: 4000 trials give ±0.006 on a coverage near 0.2.

**Common mistakes.** Passing the rate instead of the scale to `rng.gamma`; drawing the true rate from the heritage prior in (b) (which just repeats (a)); concluding from (a) that Bayesian intervals "are calibrated" without the qualifier *under the prior*.

## Check-yourself answers

1. p(λ | D) = p(D | λ)p(λ)/∫p(D | λ′)p(λ′)dλ′: prior p(λ) (a density over λ), likelihood p(D | λ) (a function of λ but a probability of the data, not a density over λ), evidence (a number), posterior (a density over λ).
2. A prior family that the likelihood maps into itself, so the posterior has the same form with updated parameters. For Gamma(α, β) with a Poisson likelihood, α counts prior failures and β prior hours.
3. (α, β) → (α + k, β + T); mean (α + k)/(β + T) = [β/(β + T)](α/β) + [T/(β + T)](k/T); weights 0.952 and 0.048.
4. Because e^{−λT} with T ~ 10⁶ h spans fifty decades across the λ range of interest and products of such factors underflow to zero; logs turn the products into sums and the maximum is subtracted before exponentiating.
5. R(t) = (β/(β + t))^α; Jensen says R(t) ≥ e^{−(α/β)t}, the reliability at the mean rate understates the predictive.
6. It is the 95 % upper credible bound of the flat-prior posterior Gamma(k + 1, T), because 2Tλ is then χ²_{2k+2}; the prior is Gamma(1, 0), flat on λ.
7. 40 000/6.5 ≈ 6150 independent samples; a probability near 0.33 has scatter √(0.33 × 0.67/6150) ≈ 0.006.
8. A 90 % credible interval is calibrated if, over repeated problems, it contains the true value 90 % of the time; this holds when the true values are distributed as the prior says, and fails (problem 6: 20 % coverage) when the prior is wrong.
