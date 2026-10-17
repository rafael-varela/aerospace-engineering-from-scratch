# P02.16  Random Variables & Distributions

> A random variable is a number you cannot predict but whose distribution you can; the distribution, not the number, is what you design against.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the mean and variance of the binomial, Poisson, uniform and exponential distributions, the memoryless property, the linear-propagation law Cov(Jx) = J Σ Jᵀ, and the probability 1 − e^{−k²/2} inside the k-sigma ellipse of a two-dimensional Gaussian. 2. Implement and test, from scratch, pmf/pdf/cdf routines, Simpson-rule moments, a bisection normal quantile, covariance propagation, the closed-form 2×2 eigenproblem, dispersion ellipses and circular containment by quadrature, against `scipy.stats` as the oracle. 3. Measure, by seeded simulation, how fast a sum of uniform variables approaches the normal distribution and how closely 20 000 synthetic landings fill the predicted ellipses. 4. Decide the size and orientation of a 99 % recovery zone for a parachute payload from four input uncertainties, compare an ellipse with a circle of the same probability, and name the input whose improvement shrinks the zone most. |
| **Prerequisites** | `P02.15` (sample spaces, binomial counting, Monte Carlo error), `P01.04` (integration by parts, for the exponential moments), `P01.16` (eigenvalues of a symmetric matrix), `P00.03`, `P00.06`, `P00.07`, `P00.14`. Skill check: you can integrate x e^{−x} by parts, find the eigenvalues of a 2×2 symmetric matrix, and explain what a standard deviation is. |
| **Engineering question** | A parachute-recovered sounding-rocket payload has four uncertain inputs (wind speed ±3 m/s, drag coefficient ±5 %, entry flight-path angle ±0.3°, heading ±0.5°, all 1σ): how large, how elongated and how oriented is the zone that contains 99 % of landings, how much larger must a *circular* zone be, and which input should the team spend money on? |
| **Assumptions and limits** | Inputs are independent, Gaussian and small enough that the landing offset is *linear* in them (a fixed sensitivity matrix J); then the landing offset is exactly bivariate Gaussian and every ellipse result is exact. The model fails when the sensitivities depend on the inputs (strong winds, large angle errors), when inputs are correlated (a wind forecast error that also biases the heading), or when an input has heavier tails than a Gaussian (gusts). The sensitivities and standard deviations are illustrative, and every landing in the lesson is *synthetic*, drawn from the model with a stated seed. The Poisson model for alerts assumes independent events at a constant rate; the exponential lifetime assumes no wear-out. |
| **Analytical work** | Definitions of pmf, pdf, cdf; expectation and variance and their linearity; binomial mean np and variance np(1 − p); Poisson as the limit of the binomial and its equal mean and variance; uniform and exponential moments by integration; the memoryless property; z-scores and the empirical rule; the central limit theorem (statement, and the 1/√n scaling of the mean's standard deviation); Cov(Jx) = J Σ Jᵀ; the eigen-decomposition of a 2×2 covariance; the chi-square-with-two-degrees-of-freedom argument that gives 1 − e^{−k²/2}; the CEP of a circular Gaussian, 1.1774σ. |
| **Computational deliverable** | `starter/main.py`: `simpson_integral`, `discrete_mean`, `discrete_variance`, `binomial_pmf`, `poisson_pmf`, `binomial_poisson_gap`, `uniform_pdf`, `uniform_cdf`, `exponential_pdf`, `exponential_cdf`, `memoryless_gap`, `normal_pdf`, `normal_cdf`, `normal_quantile`, `z_score`, `normal_interval`, `continuous_mean`, `continuous_variance`, `sample_mean_sd`, `empirical_rule_fractions`, `clt_sum_of_uniforms`, `clt_distance_to_normal`, `sample_covariance`, `propagate_covariance`, `landing_sensitivities`, `simulate_landings`, `eig2x2_symmetric`, `ellipse_containment_probability`, `k_for_probability`, `dispersion_ellipse`, `mahalanobis_squared`, `fraction_inside_ellipse`, `circle_containment_probability`, `radius_for_probability`. Tests prove: every pmf sums to one and every pdf integrates to one; moments match the closed forms and `scipy.stats`; the cdf is the integral of the pdf; the quantile inverts the cdf; the binomial–Poisson gap falls as 1/n; the CLT distance falls with n; propagated covariance matches the sample covariance of simulated landings; the 2×2 eigensolver matches `numpy.linalg.eigh`; ellipse containment matches the chi-square cdf and the simulated fraction inside; circular containment reduces to the closed form. |
| **Practical activity** | **computational experiment**: a CLT convergence study (distance to the normal cdf against the number of summed terms), a 20 000-landing Monte Carlo check of the ellipse probabilities, and an uncertainty budget that halves each input in turn. All data are synthetic. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `landing_dispersion.png`, the hand derivations of the four moments and of the ellipse probability, the CLT table, the ellipse-against-simulation table, the uncertainty budget, and a short note recommending a recovery zone and one input to improve. |
| **Assessment** | `assessments/phase-02/grading-rubric.md`; competencies 1–4. Pass: the derivations are complete including the ellipse probability; all public tests pass; the hand covariance for the landing model agrees with the code to four figures; the note gives the 99 % zone's dimensions and orientation, the circle-to-ellipse area ratio, and the input to improve with the area reduction it buys. |
| **References** | `@openstax-introductory-statistics-2e` §4.1 "Probability Distribution Function (PDF) for a Discrete Random Variable", §4.2 "Mean or Expected Value and Standard Deviation", §4.3 "Binomial Distribution", §4.6 "Poisson Distribution", §5.1 "Continuous Probability Functions", §5.2 "The Uniform Distribution", §5.3 "The Exponential Distribution", §6.1 "The Standard Normal Distribution", §6.2 "Using the Normal Distribution", §7.1 "The Central Limit Theorem for Sample Means (Averages)", §7.2 "The Central Limit Theorem for Sums" (all required-free). Multivariate Gaussians, covariance propagation and dispersion ellipses are not covered by Introductory Statistics; part 2.7 is original instruction building on `P01.16`. |
| **Extensions** | ★ Replace the Gaussian wind error by a Student-t with 4 degrees of freedom of the same standard deviation and measure how much the simulated 99 % containment of the Gaussian ellipse drops. ★ Add a correlation of 0.5 between the wind-speed and heading errors to Σ and show, by hand and by code, how the ellipse rotates. ★ Derive the Rayleigh distribution of the miss distance for a circular Gaussian and fit it to the simulated radii. |
| **Time** | unmeasured |

## 1. The engineering problem

Lesson 15 ended with a payload drifting through coin-flip winds and a scatter of landing points that looked like a smeared circle. Real dispersion analyses start from a different place: a list of the things the flight does not know exactly, each with a standard deviation, and a sensitivity of the landing point to each one. For the parachute-recovered payload in this lesson the list is: the wind speed along the descent (±3 m/s), the drag coefficient of the parachute (±5 %), the flight-path angle at the start of the descent (±0.3°) and the heading (±0.5°). A trajectory code has been run with each input perturbed and gives, in metres of landing offset per unit of input, the **sensitivity matrix** J. All of these numbers are illustrative; the point is the method.

The recovery team wants a zone that contains 99 % of landings. The range-safety officer wants to know the probability of landing beyond a road 1.2 km down-range. The project manager wants to know which of the four uncertainties to spend money on. All three questions need the **distribution** of the landing point, a random quantity in two dimensions. This lesson supplies the vocabulary (random variable, pmf, pdf, cdf, expectation, variance), the five distributions every engineer must know, the theorem that explains why so many measured errors are Gaussian, and the two-dimensional result that turns a covariance matrix into an ellipse on a map.

Success criteria, written down first: every from-scratch distribution matches `scipy.stats` to 10⁻¹⁰; Simpson-rule moments of the uniform and exponential pdf agree with the closed forms to 10⁻⁶; the covariance propagated through J agrees with the sample covariance of 20 000 synthetic landings to within 1 %; the fraction of simulated landings inside the 1σ, 2σ and 99 % ellipses agrees with 0.3935, 0.8647 and 0.99 to within four standard errors (0.0035, 0.0024 and 0.0007); and the circular-error-probable routine reproduces 1.1774σ for a circular covariance to five figures.

## 2. Understand and derive

### 2.1 Random variables and their distributions

A **random variable** X is a number assigned to each outcome of an experiment: the number of working IMUs, the time to a component failure, the down-range landing offset. If X takes a finite or countable set of values it is **discrete** and is described by its **probability mass function** (OpenStax says "probability distribution function") P(X = x), with every value non-negative and Σ P(x) = 1 (`@openstax-introductory-statistics-2e` §4.1). If X takes values in a continuum it is **continuous** and P(X = x) = 0 for every single x; it is described instead by a **probability density function** f(x) ≥ 0 with ∫ f dx = 1, and probabilities are areas: P(a ≤ X ≤ b) = ∫ₐᵇ f(x) dx (§5.1). Either way the **cumulative distribution function** F(x) = P(X ≤ x) exists, is non-decreasing from 0 to 1, and for a continuous variable F′ = f.

The **expectation** (mean) and **variance** are

**E[X] = Σ x P(x)** or **∫ x f(x) dx**,   **Var X = E[(X − μ)²] = E[X²] − μ²**,   σ = √Var X.   (1)

Expectation is linear, E[aX + bY] = a E[X] + b E[Y], always; variance is not: Var(aX + b) = a² Var X, and Var(X + Y) = Var X + Var Y only when X and Y are **independent** (in general the covariance 2 Cov(X, Y) is added). These two facts carry the whole lesson. Units: E[X] and σ have the units of X; Var X has them squared.

### 2.2 Binomial and Poisson

The **binomial** variable (§4.3) counts successes in n independent trials of probability p; Lesson 15 derived its pmf C(n, k) pᵏ(1 − p)ⁿ⁻ᵏ by counting. Its mean and variance follow from linearity: X is a sum of n indicator variables Iᵢ ∈ {0, 1} with E[Iᵢ] = p and Var Iᵢ = p − p² = p(1 − p), so

**E[X] = np**,   **Var X = np(1 − p)**.   (2)

For the three IMUs at p = 0.98: pmf (0.000008, 0.001176, 0.057624, 0.941192), mean 2.94, variance 0.0588; the code's `discrete_mean` and `discrete_variance` reproduce (2) from the pmf.

The **Poisson** variable (§4.6) counts events that occur independently at a constant average rate: conjunction alerts per month, micrometeoroid impacts on an array, photons on a detector pixel. Take the binomial with n trials and p = λ/n and let n → ∞ with λ fixed: C(n, k)(λ/n)ᵏ(1 − λ/n)ⁿ⁻ᵏ → λᵏ e^{−λ}/k!, because C(n, k)/nᵏ → 1/k! and (1 − λ/n)ⁿ → e^{−λ} (P01.05). So

**P(X = k) = λᵏ e^{−λ}/k!**,   **E[X] = Var X = λ**.   (3)

The mean and variance follow from (2) in the limit: np = λ and np(1 − p) → λ. A satellite operator receiving 2.4 conjunction alerts a month has P(0) = e^{−2.4} = 0.0907 and P(≥ 5) = 0.0959; a month with five alerts is not evidence that anything changed. The code evaluates (3) through logarithms, k ln λ − λ − ln Γ(k + 1), because λᵏ and k! overflow separately long before their ratio does. `binomial_poisson_gap` measures the approach to the limit: the largest pmf difference between Bin(n, 2.4/n) and Poisson(2.4) is 0.0339, 0.0028 and 0.00028 for n = 10, 100, 1000, a clean 1/n.

### 2.3 Uniform and exponential

The **uniform** density on [a, b] (§5.2) is 1/(b − a). By (1), E[X] = (a + b)/2 and Var X = (b − a)²/12: integrate (x − μ)² over a symmetric interval of width w = b − a to get w³/12 and divide by w. The quantisation error of an analogue-to-digital converter is uniform over ±½ least-significant bit, so its standard deviation is 1/√12 = 0.2887 LSB; `continuous_variance` by Simpson's rule gives 0.083333 for the variance.

The **exponential** density (§5.3) f(t) = λ e^{−λt} for t ≥ 0 describes the waiting time to the next Poisson event, or the lifetime of a component that fails at a constant **hazard rate** λ (no infant mortality, no wear-out). Integrating by parts (P01.04), E[T] = 1/λ (the mean time between failures, MTBF) and Var T = 1/λ². Its cdf is F(t) = 1 − e^{−λt}. A unit with MTBF 5000 h fails within 1000 h with probability 1 − e^{−0.2} = 0.1813 and survives 10 000 h with probability e^{−2} = 0.1353: the exponential has a long tail, and "MTBF" is not "the time it will last".

Its defining feature is **memorylessness**: P(T > s + t | T > s) = e^{−λ(s+t)}/e^{−λs} = e^{−λt} = P(T > t). A unit that has run 3000 h is exactly as likely to survive the next 1000 h as a new one; the code's `memoryless_gap` returns 0 to rounding. Real hardware wears out (the hazard rises with age), and that is exactly what the exponential cannot model; Lesson 20 fits a failure rate from data.

### 2.4 The normal distribution

The **normal** (Gaussian) density with mean μ and standard deviation σ is

**f(x) = exp(−(x − μ)²/(2σ²)) / (σ√(2π))**,   **F(x) = ½[1 + erf((x − μ)/(σ√2))]**.   (4)

The cdf has no elementary closed form; the **error function** erf is the one primitive this lesson borrows from `scipy.special`, and a test checks that F equals the Simpson integral of f. The **z-score** z = (x − μ)/σ (§6.1) measures distance from the mean in standard deviations and reduces every normal to the **standard normal** N(0, 1). The **empirical rule** (§6.2): P(|z| ≤ 1, 2, 3) = 0.6827, 0.9545, 0.9973. The two-sided 95 % and 99 % points are z = 1.9600 and 2.5758; `normal_quantile` finds them by bisection on F and the test compares with `scipy.stats.norm.ppf`. An accelerometer whose noise is N(0, 0.02 m/s²) produces 95 % of its samples within ±0.0392 m/s².

### 2.5 Why so much is Gaussian: the central limit theorem

The **central limit theorem** (§7.1, §7.2): if X₁, …, Xₙ are independent with the same distribution, mean μ and variance σ², then as n grows the sum S = Σ Xᵢ is approximately N(nμ, nσ²) and the mean X̄ = S/n approximately N(μ, σ²/n), *whatever the distribution of the Xᵢ* (as long as its variance is finite). The standard deviation of a mean of n samples is σ/√n: averaging 100 accelerometer samples cuts the noise ten-fold. Measurement errors are sums of many small independent contributions (thermal noise, vibration, quantisation, reading error), which is why they are so often Gaussian; the landing offset of 2.7 is Gaussian for the same reason, as a sum of four Gaussian contributions.

How fast? `clt_sum_of_uniforms` standardises S for n uniform terms (mean n/2, variance n/12) and `clt_distance_to_normal` measures the largest gap between its empirical cdf and Φ, from 200 000 synthetic draws with seed 16:

| n | max cdf gap | fraction within 1σ (Gaussian 0.6827) |
|---|---|---|
| 1 | 0.0580 | 0.5782 |
| 2 | 0.0162 | 0.6520 |
| 5 | 0.0084 | 0.6725 |
| 12 | 0.0044 | 0.6777 |
| 30 | 0.0026 | 0.6793 |

Two uniforms already give a triangle that is within 0.016 of Gaussian; twelve are within 0.004. (The sum of twelve uniforms minus six was once a standard way to generate Gaussian random numbers.) The remaining gap at n = 30 is partly the finite sample: 200 000 draws give an empirical cdf with its own scatter of about 1/√200 000 = 0.002.

### 2.6 Samples against populations

Everything above is about the **population**: the distribution the variable *has*. A **sample** of n draws gives estimates: the sample mean x̄ and the sample standard deviation s with the n − 1 divisor (`sample_mean_sd`), and for several variables at once the **sample covariance matrix** S = Σ (xᵢ − x̄)(xᵢ − x̄)ᵀ/(n − 1). Lesson 17 asks how far these estimates can be from the truth; here we only need to know that for n = 20 000 the sample covariance is within about √(2/n) = 1 % of the population covariance, which is what the test asserts.

### 2.7 Two dimensions: covariance propagation and the dispersion ellipse

Let the vector of input errors **x** (wind speed, drag coefficient, flight-path angle, heading) have zero mean and covariance matrix Σ, here diagonal because the inputs are independent: Σ = diag(3², 0.05², 0.3², 0.5²). If the landing offset is linear in the inputs, **y** = J **x**, then by linearity of expectation E[**y**] = 0 and

**Cov(y) = E[J x xᵀ Jᵀ] = J Σ Jᵀ**.   (5)

This is the matrix form of Var(aX) = a² Var X, and it is exact for any input distribution; if the inputs are Gaussian, **y** is Gaussian too. With the model's

J = [[150, −1200, −400, 0], [40, 0, 0, 250]] (m per m/s, per unit C_D, per degree, per degree),

each input's 1σ contribution to (x, y) is J_{·i} σᵢ: wind (450, 120) m, drag (−60, 0), flight-path angle (−120, 0), heading (0, 125). Summing squares down each column gives the diagonal of (5), and the wind alone gives the off-diagonal:

Cov = [[220 500, 54 000], [54 000, 30 025]] m²,   σ_x = 469.6 m, σ_y = 173.3 m, ρ = 54 000/(469.6 × 173.3) = 0.664.

The correlation is there because the wind pushes both down-range and cross-range: a landing that is long is also likely to be right of track.

```
      cross-range y
          ▲          ╱ major axis, 14.8° from down-range
          │     ╱ ╱ ╱
          │  ╱  ╱  ╱   ← 1σ ellipse: semi-axes 485 m and 126 m
  ────────┼──────────────▶ down-range x
       ╱  │  ╱
    ╱  ╱  │
```

The curves of constant density of a two-dimensional Gaussian are the ellipses **y**ᵀ Cov⁻¹ **y** = k² (the left side is the squared **Mahalanobis distance**). Write Cov = V Λ Vᵀ with V the eigenvectors and Λ = diag(λ₁, λ₂) (P01.16). In the rotated coordinates **u** = Vᵀ**y** the components are independent with variances λ₁ and λ₂, so the ellipse has semi-axes k√λ₂ ≥ k√λ₁ along the eigenvectors. For a 2×2 symmetric matrix the eigenvalues are m ± r with m = (a + d)/2, r = √(((a − d)/2)² + b²); here λ = 15 781 and 234 744 m², semi-axes 126 and 485 m at k = 1, major axis tilted atan2(v_y, v_x) = 14.8° from down-range (`eig2x2_symmetric`, `dispersion_ellipse`).

**How much probability is inside the k-sigma ellipse?** In the rotated coordinates u₁/√λ₁ and u₂/√λ₂ are independent standard normals, and the ellipse is the disk z₁² + z₂² ≤ k². The density e^{−(z₁² + z₂²)/2}/(2π) depends only on the radius ρ, so integrating in polar coordinates,

**P(inside k-sigma ellipse) = ∫₀ᵏ e^{−ρ²/2} ρ dρ = 1 − e^{−k²/2}**.   (6)

(Statisticians say z₁² + z₂² is chi-square with two degrees of freedom; `scipy.stats.chi2.cdf(k², 2)` is the oracle in the tests.) The 1σ ellipse contains only **39.35 %** (not 68 %: in two dimensions there are more ways to be far away), the 2σ ellipse 86.47 %, and 95 % and 99 % need k = √(−2 ln 0.05) = 2.448 and √(−2 ln 0.01) = 3.035. The 99 % recovery zone is therefore an ellipse with semi-axes 1470 m and 381 m, tilted 14.8° right of down-range, area 1.761 km².

**Why not a circle?** The smallest circle containing 99 % has radius 1255 m (found by `radius_for_probability`, which integrates the Gaussian over a disk by Simpson's rule and bisects on the radius) and area 4.945 km²: 2.81 times the ellipse. For a *circular* covariance σ²I the disk integral has the closed form 1 − e^{−R²/(2σ²)} and the 50 % circle, the **circular error probable** (CEP), has radius σ√(2 ln 2) = 1.1774σ; for the model's elongated covariance the CEP is 354 m. A recovery team searching a circle around the aim point is searching the wrong shape.

### 2.8 Where the model stops being valid

- **Linearity.** J is a first derivative. A 3 m/s wind error moves the landing 450 m; a 9 m/s error does not move it 1350 m if the parachute's descent time changes with the wind. Check the trajectory code at ±3σ before trusting (5).
- **Independent inputs.** If the wind forecast error also biases the heading estimate, Σ gains an off-diagonal term and the ellipse rotates (extension 2).
- **Gaussian inputs.** Gust statistics have heavier tails than a Gaussian; the 1σ and 2σ ellipses barely change but the 99 % ellipse under-covers (extension 1). The 99 % zone is a *tail* statement and tails are where Gaussian assumptions fail first.
- **Constant rate (Poisson, exponential).** Alert rates change with solar activity, and hardware wears out. Both distributions describe the no-memory case only.
- **CLT with few terms.** The sum of two uniforms is a triangle, not a Gaussian; a lognormal or heavy-tailed term needs far more than 12 terms. The theorem is about the limit; the table in 2.5 is about the rate.

## 3. Build the computational model

Open `starter/main.py`. Build in the order below; each step's tests are independent.

### Step 1: Simpson's rule and the discrete moments

```python
def simpson_integral(f, a, b, n=2000):
    if n < 2 or n % 2: raise ValueError("n must be an even integer >= 2")
    x = np.linspace(a, b, n + 1); y = f(x); h = (b - a) / n
    return float(h / 3.0 * (y[0] + y[-1] + 4.0 * y[1:-1:2].sum() + 2.0 * y[2:-1:2].sum()))

def discrete_mean(values, probs):       # check probs >= 0 and sum to 1, then
    return float(np.dot(values, probs))
```

`discrete_variance` is Σ (x − μ)² P(x). `binomial_pmf` uses `math.comb`; `poisson_pmf` uses `math.exp(k * math.log(lam) - lam - math.lgamma(k + 1))`, with the λ = 0 case handled separately.

### Step 2: the continuous distributions

Write `uniform_pdf/cdf`, `exponential_pdf/cdf`, `normal_pdf` and `normal_cdf` from (4) as vectorised NumPy functions that return zero density outside the support (use `np.where`, and clip the argument of `np.exp` so that x < 0 does not produce an overflow warning in the exponential). `memoryless_gap(lam, s, t)` returns S(s + t)/S(s) − S(t) with S = 1 − F. `continuous_mean` and `continuous_variance` hand x f(x) and (x − μ)² f(x) to `simpson_integral` over a finite interval the caller chooses to hold the density.

### Step 3: the normal quantile by bisection

```python
def normal_quantile(p, tol=1e-12):
    lo, hi = -40.0, 40.0
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if normal_cdf(mid) < p: lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)
```

Fifty halvings take [−40, 40] to 10⁻¹²; `normal_interval(mu, sigma, prob)` uses it with p = ½ + prob/2.

### Step 4: samples and the CLT

`sample_mean_sd` is the two-pass formula (mean first, then Σ(x − x̄)² /(n − 1)); `empirical_rule_fractions` counts |z| ≤ 1, 2, 3. `clt_sum_of_uniforms(n_terms, n_trials, seed)` draws an `(n_trials, n_terms)` uniform array, sums the rows and standardises by n/2 and √(n/12); `clt_distance_to_normal` sorts the result, forms the empirical cdf k/N and returns max |ecdf − Φ|.

### Step 5: covariance propagation and the ellipse

```python
def propagate_covariance(J, Sigma):     return J @ Sigma @ J.T                     # (5)

def eig2x2_symmetric(C):
    a, b, d = C[0, 0], C[0, 1], C[1, 1]
    m, r = 0.5 * (a + d), math.sqrt((0.5 * (a - d)) ** 2 + b * b)
    lam = np.array([m - r, m + r])
    v2 = np.array([b, lam[1] - a]); v2 /= np.linalg.norm(v2)   # eigenvector of the larger eigenvalue (b != 0)
    V = np.column_stack([[-v2[1], v2[0]], v2])
    return lam, V

def dispersion_ellipse(cov, prob):
    lam, V = eig2x2_symmetric(cov)
    k = k_for_probability(prob)                                   # sqrt(-2 ln(1 - prob)), from (6)
    a, b = k * math.sqrt(lam[1]), k * math.sqrt(lam[0])
    angle = math.degrees(math.atan2(V[1, 1], V[0, 1]))            # fold into (-90, 90]
    return {"k": k, "a": a, "b": b, "angle_deg": angle, "area": math.pi * a * b}
```

Handle b = 0 (already diagonal) separately in the eigensolver, and raise `ValueError` for a non-symmetric input. `mahalanobis_squared` is **d**ᵀ Cov⁻¹ **d** per row (`np.einsum("ij,jk,ik->i", d, Ci, d)`); `fraction_inside_ellipse` counts rows with it ≤ k². `simulate_landings` draws inputs with `rng.multivariate_normal(0, Sigma)` and multiplies by Jᵀ; `landing_sensitivities` returns the model's J, Σ and input names.

### Step 6: circles

`circle_containment_probability(cov, R)` integrates, in principal axes, f₁(u) [Φ(h/s₂) − Φ(−h/s₂)] with h = √(R² − u²) over u ∈ [−R, R] by Simpson; `radius_for_probability` bisects on R until the containment equals the target. With prob = 0.5 this is the CEP.

Running `main()` prints the tables of 2.2–2.7 and, for the landing model:

```
sigma_x = 469.6 m, sigma_y = 173.3 m, correlation 0.664
  0.3935 ellipse: k = 1.000, a = 485 m, b = 126 m, angle 14.8 deg, area 0.191 km^2; simulated fraction inside 0.3926
  0.8647 ellipse: k = 2.000, a = 969 m, b = 251 m, angle 14.8 deg, area 0.765 km^2; simulated fraction inside 0.8649
  0.9900 ellipse: k = 3.035, a = 1470 m, b = 381 m, angle 14.8 deg, area 1.761 km^2; simulated fraction inside 0.9890
CEP (50 % circle) = 354 m; 99 % circle radius = 1255 m, area 4.945 km^2 vs 99 % ellipse area 1.761 km^2 (ratio 2.81)
```

### Verify

Run `python -m pytest curriculum/phase-02/lesson-16-random-variables-distributions`. The tests:

1. **Analytic cases.** Simpson integrates a cubic exactly; the binomial pmf matches `scipy.stats.binom` and has mean np and variance np(1 − p); the Poisson pmf matches `scipy.stats.poisson` and has mean = variance = λ; uniform and exponential moments by Simpson match (a + b)/2, (b − a)²/12, 1/λ, 1/λ² to 10⁻⁶; the normal cdf equals the Simpson integral of the pdf to 10⁻⁹; the quantile inverts the cdf and matches `norm.ppf` to 10⁻⁹.
2. **Limiting cases.** The binomial–Poisson gap falls by at least a factor 5 per decade of n; the memoryless gap is zero for the exponential and non-zero for a uniform lifetime; the CLT distance falls monotonically with n; the 2×2 eigensolver on a diagonal matrix returns the diagonal and the identity.
3. **Conservation.** Every pmf sums to 1, every pdf integrates to 1; trace(J Σ Jᵀ) equals the sum of the squared 1σ contributions; the eigenvalues of the covariance sum to its trace and the eigenvectors are orthonormal and match `numpy.linalg.eigh` up to sign.
4. **Agreement with simulation.** The sample covariance of 20 000 landings is within 3 % of (5); the fraction inside the k = 1, 2, 2.448, 3.035 ellipses is within four standard errors of (6); the fraction within the CEP is within four standard errors of 0.5.
5. **Circles.** For a circular covariance, `circle_containment_probability` equals 1 − e^{−R²/(2σ²)} to 10⁻⁵ (the square-root end-point singularity of the disk integrand limits Simpson's rule to about 10⁻⁶) and `radius_for_probability(…, 0.5)` equals 1.1774σ to 10⁻⁴; for the elongated covariance the 99 % circle contains more area than the 99 % ellipse; containment is monotone in R.
6. `main()` writes `landing_dispersion.png` to the directory it is given.

### Compare

`scipy.stats` provides every distribution here as an object with `.pmf`, `.pdf`, `.cdf`, `.ppf`, `.mean()`, `.var()` and `.rvs()`, and `scipy.stats.multivariate_normal` handles the two-dimensional Gaussian; `scipy.stats.chi2.cdf(k**2, 2)` is (6). Our versions match them to 10⁻¹⁰ wherever both exist. The library is the right tool in practice because it handles the numerical corners (extreme tails, large parameters, the inverse cdf of every distribution) that a bisection on an erf cannot, and it is tested on millions of cases. We wrote our own to know that a cdf is an integral, a quantile is an inverse, a covariance is J Σ Jᵀ and an ellipse is an eigenproblem, so that when the library returns a 99 % zone we can say what it assumed.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All data are synthetic, generated by the stated models with stated seeds.

1. **CLT rate.** Run `clt_distance_to_normal(n, 200_000, seed=16)` for n = 1, 2, 3, 5, 8, 12, 20, 30 and plot the distance against n on log-log axes. Before running, predict the slope (the Berry–Esseen bound gives a 1/√n *upper* bound; for symmetric distributions the actual rate is closer to 1/n). Measured at n = 1, 2, 5, 12, 30: 0.0580, 0.0162, 0.0084, 0.0044, 0.0026.
2. **Ellipse fill.** Simulate 20 000 landings with seed 2016 and tabulate the fraction inside the k = 0.5, 1, 1.5, 2, 2.5, 3 ellipses against 1 − e^{−k²/2}, with the standard error √(p(1 − p)/N) for each. Measured at k = 1, 2, 3.035: 0.3926, 0.8649, 0.9890 against 0.3935, 0.8647, 0.9900.
3. **Uncertainty budget.** Halve each input standard deviation in turn (wind to 1.5 m/s, drag to 2.5 %, angle to 0.15°, heading to 0.25°), recompute the 99 % ellipse area, and rank the four. Predict first from the 1σ contributions in 2.7 (wind 450 m down-range dominates). The budget decides where the money goes.
4. **Heavy tails (extension 1 made concrete).** Replace the wind draw by `rng.standard_t(4) * 3.0 / sqrt(2)` (a t₄ scaled to σ = 3 m/s) and re-measure the fraction inside the Gaussian 99 % ellipse.

Raw data format: CSV tables `n, cdf_gap`; `k, fraction, predicted, se`; `input_halved, area_km2`; and a one-line result for item 4.

What you should find: the CLT distance falls faster than 1/√n for the uniform (symmetric) case; every ellipse fraction lies within two standard errors of (6); halving the wind error shrinks the 99 % area by far more than halving any other input, and halving the drag error does almost nothing; the heavy-tailed wind leaves the 1σ fraction nearly unchanged but drops the 99 % containment to about 98 %: the zone that was "99 %" is no longer.

## 5. Validate and improve

| quantity | hand (§2) | code / simulation | difference | cause |
|---|---|---|---|---|
| Var U(−½, ½) | 1/12 = 0.083333 | 0.083333 (Simpson, 4000 panels) | < 10⁻⁹ | Simpson exact for quadratics |
| P(T < 1000 h), MTBF 5000 h | 1 − e^{−0.2} = 0.1813 | 0.1813 | 0 | closed form |
| z for two-sided 95 % | 1.96 (table) | 1.9600 (bisection) | < 10⁻⁹ | both exact |
| σ_x, σ_y, ρ of landing | 469.6 m, 173.3 m, 0.664 | 471.4 m, 173.9 m, 0.672 (20 000 landings) | 0.4 %, 0.4 %, 1.2 % | sampling: SE of an SD ≈ σ/√(2N) = 0.5 % |
| fraction in 1σ ellipse | 0.3935 | 0.3926 ± 0.0035 | 0.3 SE | sampling |
| fraction in 99 % ellipse | 0.9900 | 0.9890 ± 0.0007 | 1.4 SE | sampling |
| CEP, circular σ = 469.6 m | 1.1774σ = 552.88 m | 552.88 m (quadrature + bisection) | < 0.01 m | quadrature tolerance |
| 99 % circle area / ellipse area | — | 4.945 / 1.761 = 2.81 | — | elongation √(λ₂/λ₁) = 3.86 |

Every difference is within the sampling error of 20 000 draws or the stated numerical tolerance. The area ratio deserves a comment: the ellipse has aspect ratio 3.86, and a circle must reach the far end of the major axis, so a crude estimate of the ratio is πR²/(πab) with R ≈ a × (something below 1, because the circle need not reach the 99 % *ellipse* tip to hold 99 % *probability*): the quadrature gives R = 1255 m against a = 1470 m, so the ratio 2.81 is less than the aspect ratio 3.86. Searching a circle wastes 64 % of the effort.

What the agreement does not show: that the inputs are Gaussian, independent and small. The simulation draws from exactly the model the analysis assumes, so agreement proves the code, not the model; it is a *verification*, in the language of P00.12, not a validation. Validation would need measured landings (a few dozen flights give a covariance with 15–20 % uncertainty, enough to catch a factor-2 error in J but not a 20 % one) or, failing that, a comparison of the linear model with the full trajectory code at ±3σ inputs.

Improvements worth making now: carry the full trajectory code's nonlinearity by sampling the inputs and running the code (the Monte Carlo route; the ellipse then becomes a fitted, not derived, quantity); add the wind–heading correlation of extension 2; replace the Gaussian wind by the measured gust distribution once Lesson 19 provides one. Improvement worth *not* making: quoting the 99.9 % ellipse. At k = 3.717 the Gaussian tail assumption is the whole answer, and nothing in this lesson has tested it.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivations of (2), (3), the uniform and exponential moments, the memoryless property and (6); `landing_dispersion.png`; the four tables from part 4; and a one-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) that specifies the 99 % recovery zone (semi-axes, orientation, area), states the penalty of a circular zone, recommends one input to improve with the area reduction it buys, and names the assumption that would fail first.

### Adversarial self-checks

- **What if every input standard deviation doubles?** Σ scales by 4, Cov by 4, every semi-axis by 2, every area by 4; the orientation and the probabilities inside each k-ellipse are unchanged. What if only the wind doubles? The ellipse grows mostly along its major axis *and rotates*: compute the new angle.
- **What if the number of landings simulated goes to infinity? To ten?** The sample covariance converges to (5) as 1/√N; with ten landings the eigenvectors of the sample covariance are unreliable (the 1σ ellipse's orientation can be off by tens of degrees). Try it with `simulate_landings(10, …)`.
- **Which assumption fails first in the real world?** Gaussian tails on the wind, at the 99 % level. Second, linearity of J for the angle inputs. Third, independence of wind and heading.
- **What measurement could falsify the model?** The down-range–cross-range correlation of measured landings: the model predicts ρ = 0.66 from the wind sensitivities alone. A measured ρ near zero would mean the wind contributes less than assumed; a negative ρ would mean a sign error in J.
- **Does the result survive a change of method?** Ellipse semi-axes from `eig2x2_symmetric` and from `numpy.linalg.eigh` agree to 10⁻¹²; the containment probability from (6), from `scipy.stats.chi2.cdf` and from the simulated fraction agree within sampling error; CEP by quadrature agrees with the closed form in the circular case. If the Simpson panel count is halved the CEP moves by less than 0.01 m.
- **Which decision would change if the uncertainty doubled?** If the wind standard deviation were 6 m/s instead of 3, the 99 % ellipse's semi-major axis would grow from 1470 m to about 2800 m, and a road 1.2 km down-range would move from well outside the 95 % ellipse to inside it: the launch would need a different aim point or a wind limit. The *decision* (go, or wait for calmer winds) depends on the tail of the wind distribution, which is exactly the part the model knows least.

## Check yourself

1. Define pmf, pdf and cdf, and say which of the three exists for both discrete and continuous random variables.
2. A variable has E[X] = 3 and Var X = 4. What are the mean and variance of Y = 2X − 5? Of the mean of 16 independent copies of X?
3. Why do the binomial mean np and variance np(1 − p) both tend to λ in the Poisson limit, and what does "Var X = E[X]" let you check on a dataset of counts?
4. What does "memoryless" mean for an exponential lifetime, and what real feature of hardware does it rule out?
5. State the central limit theorem in one sentence and give the standard deviation of the mean of n measurements with standard deviation σ.
6. Why does the 1σ ellipse of a two-dimensional Gaussian contain 39 % of the probability rather than 68 %?
7. Write Cov(y) for y = J x in terms of Cov(x), and say which property of expectation the derivation uses.
8. For a circular Gaussian with σ = 100 m, give the CEP and the radius of the 99 % circle; for an elongated one with the same σ_x, is the CEP larger or smaller?
