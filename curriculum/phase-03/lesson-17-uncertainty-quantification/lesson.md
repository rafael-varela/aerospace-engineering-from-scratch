# P03.17  Uncertainty Quantification

> A prediction without an uncertainty is a guess with a decimal point. Push the input distributions through the model, by sampling or by polynomials, and report what comes out: a mean, a spread, a tail probability, and which input is to blame.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the Monte Carlo estimators of a mean and of a tail probability with their 1/√N standard errors, the inversion and Box–Muller samplers, the percentile bootstrap, and the antithetic-variate reduction. 2. Derive first-order (linear) propagation with its second-order mean correction, and state when it fails. 3. Derive the Hermite polynomial-chaos expansion from orthogonality: the basis, Gauss–Hermite quadrature by Golub–Welsch, non-intrusive projection, moments and Sobol indices from the coefficients, and its spectral convergence. 4. Implement and test all of it from scratch with `scipy.stats` and `numpy.polynomial.hermite_e` as oracles, measure the observed N^{−1/2} order of Monte Carlo, and decide which route to use for a given model, budget and question. |
| **Prerequisites** | `P02.16` (random variables, normal distribution, moments, the central limit theorem), `P02.15` (probability, independence), `P01.16` (symmetric eigenvalue problems, used by Golub–Welsch), `P01.17` (least squares, for the regression PCE), `P03.04` (Gauss quadrature), `P00.14` (first-order uncertainty propagation). Skill check: you can write the variance of a sum of independent variables and the Taylor expansion of f(x₀ + δ) to second order. |
| **Engineering question** | The light twin of `P03.16` has a nominal Breguet range of 1502 km, but its specific fuel consumption is known to ±3 %, its cruise lift-to-drag ratio to ±5 % and its usable fuel to ±2 % (one standard deviation each, assumed Gaussian and independent). What are the mean and standard deviation of the range, what is the probability that it falls short of the 1350 km the mission with reserves requires, which input dominates that probability, and how many model runs does each answer cost? |
| **Assumptions and limits** | The input distributions are *assumed*: independent Gaussians with stated relative standard deviations, chosen as plausible engineering judgements, not measured statistics; nothing here is a flight measurement. The model is the Breguet cruise equation with fixed speed and zero-fuel weight. Monte Carlo assumes only that the model can be evaluated; the linear method assumes the model is nearly linear over the input spread; polynomial chaos assumes the model is smooth in the inputs (it fails on discontinuities) and, with tensor quadrature, that the number of inputs is small (cost n^d). |
| **Analytical work** | The Monte Carlo mean estimator, its variance s²/N and the binomial error of a probability; inversion sampling and the Box–Muller transform; the bootstrap principle; why antithetic pairs cancel the odd part of a model; the linear propagation u² = Σ(∂f/∂ξ_i)² and the mean shift ½ tr H; the Hermite recurrence, orthogonality E[He_m He_n] = n!δ_mn, the Jacobi matrix whose eigenvalues are the Gauss–Hermite nodes, the projection formula c_α = E[f He_α]/α!, mean = c₀, variance = Σ c_α² α!, Sobol indices as partial sums, the closed-form expansion of e^{aξ}; the hand sensitivities of the range. |
| **Computational deliverable** | `starter/main.py`: `standard_normals_box_muller`, `inverse_cdf_sample`, `monte_carlo`, `mc_probability`, `bootstrap_ci`, `antithetic_monte_carlo`, `mc_error_study`, `linear_propagation`, `hermite_e`, `gauss_hermite_e`, `pce_basis`, `pce_norms`, `pce_projection`, `pce_regression`, `pce_evaluate`, `pce_moments`, `sobol_indices` (the dataclass, `range_model`, `total_degree_multi_index` and `main` are given). Tests prove: Box–Muller samples pass moment and Kolmogorov–Smirnov checks at tolerances set by their own standard errors; inversion reproduces an exponential and a normal; the MC mean lies within 4 standard errors of the quadrature-exact mean and its tail probability within 0.01 of the normal approximation; the bootstrap interval covers the reference and matches ±1.96 se for the mean; antithetic sampling cuts the standard error of this nearly linear model by more than 5× and is exact for odd functions; the RMS error of the MC mean over replications has observed order −0.5 ± 0.12 in N; the linear sensitivities match the hand formulas to 10⁻⁶; He_n matches `numpy.polynomial.hermite_e`, is orthogonal with norm n!, and the from-scratch Gauss–Hermite nodes and weights match `hermegauss` to 10⁻¹⁰ while integrating E[ξ^{2k}] = (2k − 1)!! exactly up to degree 2n − 1 and not beyond; the PCE reproduces a quadratic exactly, gives the closed-form coefficients e^{a²/2}aⁿ/n! of e^{aξ}, converges spectrally in degree, agrees with the regression fit, and its Sobol indices match the linear variance shares and the hand values of two test functions; the surrogate's tail probability agrees with Monte Carlo within 4 standard errors. |
| **Practical activity** | **computational experiment**: the Monte Carlo convergence study over replications; the cost-accuracy table of Monte Carlo against PCE against linear propagation; the shortfall probability by three methods and its sensitivity to the assumed input spreads. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `range_uq.png`, the hand derivations of the MC error, of the Hermite orthogonality and of the projection formula, the hand sensitivities, the three experiment tables, and a half-page uncertainty statement for the range in the form "R = … ± … km (k = 1), P(R < 1350 km) = …, dominated by …, inputs assumed, N model runs". |
| **Assessment** | `assessments/phase-03/grading-rubric.md`; competencies 1–4. Pass: the MC error law and the PCE moment formulas are derived; all public tests pass; the observed MC order is measured; the uncertainty statement names the dominant input with its Sobol index and says which numbers are assumptions. |
| **References** | `@ocw-16-90` Unit 3 "Probabilistic Methods and Optimization": "3.2 Review of Probability and Statistics", "3.3 Monte Carlo Methods" (3.3.4 Inversion Method for Sampling), "3.4 Error Estimates for the Monte Carlo Method" (3.4.4 Bootstrapping), "3.5 Variance Reduction Techniques for the Monte Carlo Method" (required-free). `@jcgm-101-2008` "Propagation of distributions using a Monte Carlo method" (required-free: the metrology standard for what this lesson's Monte Carlo does, and the coverage-interval language of the uncertainty statement). `@chaospy-docs` "Main usage" > "Pseudo-spectral projection" ("Generate quadrature nodes and weights", "Expansion of orthogonal polynomials", "Fourier coefficients"), "Point collocation" ("Select polynomial expansion", "Solve the linear regression problem"); "Advanced topics" > "Generalized polynomial chaos" (supplementary-free: a software view of the two PCE constructions built here from scratch; read it after part 3 to see the same steps as library calls). |
| **Extensions** | ★ Replace the tensor Gauss–Hermite grid by a Smolyak sparse grid and count model runs at d = 6. ★ Add a correlation between c_T and L/D through a Cholesky factor of the input covariance and show how the Sobol indices lose their meaning. ★ Use the PCE surrogate inside the `P03.16` optimiser to minimise the 95th percentile of fuel instead of its mean (robust design). |
| **Time** | unmeasured |

## 1. The engineering problem

`P03.16` optimised a wing for a block fuel of 460 kg and a range of 1500 km with every input known exactly. None of them is. The engine's specific fuel consumption varies from unit to unit and with age by a few per cent; the cruise lift-to-drag ratio comes from a drag polar whose error is several per cent; the fuel actually loaded is known to the gauging accuracy. The operator's question is not "what is the range?" but **"what is the probability that this aircraft, on this mission, cannot reach its alternate with reserves?"**, and the answer must come with a statement of what was assumed and of how many model evaluations it cost, because for a real aircraft each evaluation might be a CFD run.

This lesson propagates three uncertain inputs through the Breguet range R = (V/(g c_T))(L/D) ln((W_zf + W_fuel)/W_zf) by the three methods every engineer should own. **Linear propagation** (`P00.14`) costs a handful of runs and gives a mean and a standard deviation, but assumes the model is straight over the input spread. **Monte Carlo** (16.90 Unit 3, JCGM 101) costs thousands of runs, gives the whole distribution and every tail probability, and its error falls as 1/√N whatever the number of inputs. **Polynomial chaos** costs tens of runs for a smooth model with few inputs, gives moments and sensitivities exactly from a handful of coefficients, and provides a surrogate that Monte Carlo can then sample for free. All three are built from scratch; `scipy.stats` and `numpy.polynomial.hermite_e` check them.

The input distributions are **assumed** throughout: c_T = 1.4 × 10⁻⁵ kg/(N s) ± 3 %, L/D = 17.9 ± 5 %, W_fuel = 4512 N ± 2 %, each Gaussian, independent, one standard deviation. They are plausible numbers for a light twin, not measurements, and the lesson's results are only as good as that choice.

Success criteria, written first: the three methods agree on the mean to 0.1 % and on the standard deviation to 1 %; Monte Carlo's observed error order is −0.5 ± 0.1; the from-scratch Gauss–Hermite rule matches NumPy's to 10⁻¹⁰; the PCE reproduces the closed-form coefficients of e^{aξ} to 10⁻⁸; the shortfall probability from the PCE surrogate agrees with Monte Carlo within its standard error; and the Sobol indices name the dominant input.

## 2. Understand and derive

### 2.1 Standardised inputs

Write each input as X_i = μ_i(1 + s_i ξ_i) with ξ_i a standard normal: ξ = 0 is the nominal design, ξ_i = 1 is one standard deviation above. The model becomes a function f(ξ) of a standard normal vector, which is the natural setting for Hermite polynomials and lets every method share the same sampler. For the range,

f(ξ) = (V/(g c_T,0 (1 + s₁ξ₁))) · (L/D)₀(1 + s₂ξ₂) · ln(1 + (W_fuel,0/W_zf)(1 + s₃ξ₃)),   (1)

linear in ξ₂, convex in ξ₁ (the 1/x) and concave in ξ₃ (the logarithm). f(0) = 1502.2 km.

### 2.2 Monte Carlo and its error

Draw N independent ξ⁽ⁿ⁾, evaluate y_n = f(ξ⁽ⁿ⁾), and estimate

**μ̂ = (1/N) Σ y_n, s² = (1/(N − 1)) Σ (y_n − μ̂)².**   (2)

By independence Var μ̂ = σ²/N, so the **standard error** of the mean is s/√N, and by the central limit theorem (`P02.16`) μ̂ is normal around the true mean for N beyond a few tens. The error falls as N^{−1/2} *regardless of the number of inputs*: ten times the accuracy costs a hundred times the runs, in three dimensions or three hundred. A probability p = P(y < y_crit) is the mean of the indicator 1[y < y_crit], a Bernoulli variable with variance p(1 − p), so **p̂ = (number below)/N with standard error √(p̂(1 − p̂)/N)**. For a small p this is ≈ √(p/N): to know a 5 % probability to 10 % relative needs N ≈ (1 − p)/(p · 0.01) ≈ 1900 samples, and a 10⁻⁴ probability to the same relative accuracy needs 10⁶. The run (seed 17) gives, for N = 10², 10³, 10⁴, 10⁵: means 1495.7 ± 9.4, 1503.0 ± 2.9, 1503.3 ± 0.9, 1503.7 ± 0.3 km and shortfall probabilities 0.010 ± 0.010, 0.043 ± 0.006, 0.049 ± 0.002, 0.0443 ± 0.0007. The N = 100 run says "1 %"; it is wrong by a factor 4 and its own error bar says so.

**Measuring the order.** Repeat the N-sample estimate 200 times with fresh seeds and take the root-mean-square error against the exact mean (which the PCE of 2.6 provides): 8.64, 4.35, 2.22, 1.21 km for N = 100, 400, 1600, 6400, a log-log slope of −0.474 against the theoretical −0.5. With 200 replications the RMS itself is known to about 5 %, so a slope within ±0.05 of −0.5 is what one expects, and that is what the test demands (±0.12 with 100 replications).

**Sampling.** A uniform U on (0, 1) pushed through the inverse of a distribution function, X = F⁻¹(U), has distribution F (16.90 3.3.4, the **inversion method**): P(X ≤ x) = P(U ≤ F(x)) = F(x). For the normal, F⁻¹ has no closed form, and the **Box–Muller** transform does the job with two uniforms: Z₁ = √(−2 ln U₁) cos 2πU₂, Z₂ = √(−2 ln U₁) sin 2πU₂ are independent standard normals (the radius² = −2 ln U₁ is exponential with mean 2, which is the chi-square distribution with two degrees of freedom of Z₁² + Z₂², and the angle is uniform). `standard_normals_box_muller` uses NumPy only for the uniforms; the test checks its first four moments and a Kolmogorov–Smirnov test against `scipy.stats.norm` with tolerances of five standard errors each.

**Bootstrap** (16.90 3.4.4). The formula s/√N is for the mean. For a median, a quantile or a tail probability the standard error has no simple formula, and the **percentile bootstrap** estimates it from the data alone: resample the N values with replacement B times, compute the statistic each time, and take the 2.5 % and 97.5 % quantiles of the B results as a 95 % interval. For the mean it reproduces μ̂ ± 1.96 s/√N to 10 % (the test checks); for the shortfall probability at N = 10⁴ it gives [0.0449, 0.0534], consistent with the binomial 0.049 ± 0.0022.

**Antithetic variates** (16.90 3.5). Evaluate f at ξ and at −ξ and average the pair. Write f = f_even + f_odd; the pair average is f_even(ξ) exactly, so the odd part, which for a nearly linear model is almost all of the variance, is removed, and only the variance of the even part remains. For the range at N = 10⁴ the standard error of the mean drops from 0.93 km to 0.049 km: the same accuracy as 3.6 × 10⁶ plain samples, for the price of 10⁴. The trick does nothing for the distribution's shape or for tail probabilities, which depend on the odd part; it is a mean-estimation device.

### 2.3 Linear propagation and its limit

Expand f about ξ = 0 (`P00.14`): f ≈ f(0) + Σ g_i ξ_i with g_i = ∂f/∂ξ_i. The ξ_i are independent with unit variance, so

**u² = Σ g_i², mean ≈ f(0), share of input i = g_i²/u².**   (3)

For the range the hand sensitivities are g₁ = −s₁ f(0) (from 1/x), g₂ = +s₂ f(0) (linear), g₃ = s₃ f(0) · φ/((1 + φ) ln(1 + φ)) with φ = W_fuel/W_zf = 0.1549 (differentiate the logarithm): −45.07, +75.11, +27.98 km, u = 91.95 km, shares 0.240, 0.667, 0.093. The lift-to-drag ratio owns two thirds of the variance. Keep the second-order terms and the **mean shifts**: E[f] ≈ f(0) + ½ Σ ∂²f/∂ξ_i² (the diagonal of the Hessian, because E[ξ_i ξ_j] = δ_ij and E[ξ_i²] = 1), here +1.3 km: the convexity of 1/c_T raises the mean range slightly above the nominal. The linear method's limits are visible in the run: it predicts a Gaussian range and a shortfall probability of 0.0489; the true distribution has skewness +0.10 and a shortfall probability of 0.044. Ten per cent on a probability from a nearly linear model; on a strongly nonlinear one (problem 3, a landing distance) the linear tail probability is wrong by 20 %.

### 2.4 Hermite polynomials and Gauss–Hermite quadrature

The **probabilists' Hermite polynomials** He_n are defined by He₀ = 1, He₁ = ξ, and

**He_{n+1}(ξ) = ξ He_n(ξ) − n He_{n−1}(ξ),**   (4)

so He₂ = ξ² − 1, He₃ = ξ³ − 3ξ, He₄ = ξ⁴ − 6ξ² + 3. They are orthogonal under the standard normal density: **E[He_m(ξ) He_n(ξ)] = n! δ_mn**. (The clean proof uses the generating function e^{ξt − t²/2} = Σ He_n(ξ) tⁿ/n!: the expectation of the product of two generating functions is e^{st}, and matching powers of s and t gives the result; problem 2 walks through it.) Orthogonality is the whole point: a function of ξ expanded in He_n has coefficients that are *independent* of each other, and the mean and variance drop out term by term.

The expectations will be computed by **Gauss–Hermite quadrature**, E[g(ξ)] ≈ Σ_k w_k g(ξ_k), the Gaussian rule of `P03.04` for the normal weight: n nodes integrate every polynomial of degree ≤ 2n − 1 exactly. The nodes are the roots of He_n, and the **Golub–Welsch** construction (`P01.16`) finds them without root-finding: the recurrence (4), normalised, is the eigenvalue equation of the symmetric tridiagonal **Jacobi matrix** with zeros on the diagonal and √k (k = 1 … n − 1) on the off-diagonals; its eigenvalues are the nodes and the squared first components of its unit eigenvectors are the weights (which sum to 1 with this normalisation). `gauss_hermite_e` is six lines of `np.linalg.eigh`; the tests check it against `numpy.polynomial.hermite_e.hermegauss` (whose weights carry the factor √(2π) of the unnormalised weight e^{−x²/2}) to 10⁻¹⁰, check E[ξ^{2k}] = (2k − 1)!! exactly for 2k ≤ 2n − 1, and check that the rule is *not* exact one degree higher, which is how you know the exactness claim is sharp.

### 2.5 The polynomial chaos expansion

Expand the model in the orthogonal basis. For d inputs the basis functions are products, He_α(ξ) = Π_i He_{α_i}(ξ_i), indexed by multi-indices α ∈ ℕ^d, with E[He_α He_β] = α! δ_αβ where α! = Π α_i!. A **total-degree** truncation keeps |α| = Σ α_i ≤ p, which is (d + p)!/(d! p!) terms: 10 for d = 3, p = 2; 35 for p = 4. The expansion

**f(ξ) ≈ Σ_{|α| ≤ p} c_α He_α(ξ), c_α = E[f He_α]/α!,**   (5)

is the least-squares projection of f onto the truncated basis in the norm E[·²] (Fourier analysis, `P02.09`, with the normal density as the inner product). **Non-intrusive spectral projection** computes each c_α by Gauss–Hermite quadrature on a tensor grid of n^d points: with n = p + 1 points per axis the integrand f He_α for a polynomial f of degree ≤ p has degree ≤ 2p ≤ 2n − 1, so the expansion is *exact* for polynomial models, as the test on a quadratic in two variables checks to 10⁻¹². For a smooth non-polynomial model the coefficients decay faster than any power of |α| (spectral convergence): the exponential e^{aξ} has the closed form c_n = e^{a²/2} aⁿ/n! (from the generating function with t = a), which the test reproduces to 10⁻⁸, and its coefficients fall factorially.

Once the coefficients are known everything else is free:

- **Mean** = c₀ (every other basis function has zero mean, by orthogonality with He₀ = 1).
- **Variance** = Σ_{α ≠ 0} c_α² α! (orthogonality kills the cross terms).
- **Sobol sensitivity indices.** The variance of the terms that involve *only* ξ_i, divided by the total, is the first-order index S_i: the fraction of variance explained by input i alone. The sum over every term that involves ξ_i at all is the total index S_Ti, which includes interactions. For an additive model (no mixed terms) S_i = S_Ti and Σ S_i = 1; the shortfall of Σ S_i below 1 measures interaction.
- **A surrogate.** (5) can be evaluated at a million sample points in a second, so tail probabilities, quantiles and histograms come from Monte Carlo *on the polynomial*, with the polynomial's truncation error as the only price.

For the range at degree 2 (27 model runs) the coefficients in kilometres are c₀ = 1503.535, c_{(1,0,0)} = −45.187, c_{(0,1,0)} = 75.177, c_{(0,0,1)} = 28.006 (the linear sensitivities of 2.3 to three figures, now exact rather than first-order), c_{(2,0,0)} = 1.356, c_{(0,0,2)} = −0.038, c_{(0,2,0)} = 0 (f is linear in L/D), and mixed terms c_{(1,1,0)} = −2.259, c_{(1,0,1)} = −0.842, c_{(0,1,1)} = 1.400. Two checks you can do by hand: the mean shift c₀ − f(0) = 1.318 km equals c_{(2,0,0)} + c_{(0,0,2)}, because the coefficient of He₂ is half the second derivative, which is exactly the ½ tr H of 2.3; and the variance Σ c_α² α! = 8489 km², standard deviation 92.14 km, against the linear 91.95 km: the curvature and the interactions add 0.2 % to the spread. Degrees 3 and 4 change the standard deviation by less than 10⁻⁶ relative: the model is nearly quadratic, and the expansion knows it. Sobol indices: S = (0.241, 0.666, 0.092), S_T = (0.242, 0.667, 0.093), Σ S_i = 0.9991: an almost additive model dominated by the lift-to-drag ratio. Sampling the degree-4 surrogate 2 × 10⁵ times gives P(R < 1350 km) = 0.0440, against 0.0443 ± 0.0007 from 10⁵ direct Monte Carlo runs and 0.0489 from the Gaussian approximation.

**Regression instead of projection.** If the model can only be evaluated at points you do not choose (archived runs, a sampled database), fit the coefficients by least squares (`P01.17`) to N ≥ 2 × (number of terms) random evaluations: `pce_regression`. With 400 samples at degree 2 the coefficients agree with projection to 2 × 10⁻⁵ of c₀. Point collocation is the method of choice for high-dimensional or noisy models; projection is exact and cheap for low-dimensional smooth ones.

### 2.6 Choosing a method

| method | runs for this problem | gives | fails when |
|---|---|---|---|
| linear propagation | 2d + 1 = 7 | mean, std, shares | model curved over the input spread |
| Monte Carlo | N; error ∝ N^{−1/2} | everything, any d | budget: 10⁴ runs for 1 % on the std, 10⁶ for a 10⁻⁴ probability |
| PCE projection, degree p | (p + 1)^d = 27 at p = 2 | mean, variance, Sobol, surrogate | d large (cost n^d), model non-smooth |
| PCE regression | ≥ 2 (d + p)!/(d! p!) | the same | ill-conditioned for large p |

The question decides. A mean and a standard deviation for a design review: linear, checked once against PCE. A certification probability at the 10⁻⁵ level: Monte Carlo on a PCE surrogate, with the surrogate's accuracy in the tail checked against a few hundred direct runs there. Twenty uncertain inputs: Monte Carlo, or sparse-grid or regression PCE with care.

### 2.7 Where the model stops being valid

- **The inputs are assumed.** Change the lift-to-drag spread from 5 % to 3.5 % and the shortfall probability doubles the other way (problem 4). The honest output of this lesson is conditional on the assumed distributions, and the uncertainty statement must say so.
- **Independence.** c_T and L/D may be correlated through the engine's installation; correlation changes the variance (and makes Sobol indices ambiguous), and needs a covariance, which a single engineer's judgement rarely supplies.
- **Gaussian tails.** A Gaussian fuel load admits negative fuel at −50σ; harmless here, fatal for a quantity with a hard physical bound and a wide spread.
- **Smoothness.** A model with a threshold (a stall margin that switches a penalty on) has a discontinuity in ξ; its PCE converges algebraically and rings, and Monte Carlo is the only safe route.
- **Simulation is not measurement.** Every number here is a propagation of assumptions through a model; a flight-test range with its own measurement uncertainty is what would validate it.

## 3. Build the computational model

Open `starter/main.py`. `RangeInputs`, `range_model`, `total_degree_multi_index` and `main` are given.

### Step 1: samplers and Monte Carlo

```python
def monte_carlo(f, d, n, seed):
    rng = np.random.default_rng(seed)
    y = np.asarray(f(rng.standard_normal((n, d))), dtype=float)
    s = float(np.std(y, ddof=1))
    return {"mean": float(np.mean(y)), "std": s, "se_mean": s / math.sqrt(n), "samples": y}
```

`standard_normals_box_muller` draws ⌈n/2⌉ pairs of uniforms with `rng.random` and applies the transform (guard U₁ = 0 with the smallest positive float); `inverse_cdf_sample` applies a given F⁻¹ and raises for u outside (0, 1); `mc_probability` counts and returns the binomial error; `bootstrap_ci` resamples with `rng.integers`; `antithetic_monte_carlo` averages f(ξ) and f(−ξ) over n/2 pairs and takes the standard error over the pair averages; `mc_error_study` loops replications.

### Step 2: linear propagation

Central differences at ξ = 0 with h = 10⁻⁴ for the gradient and the diagonal second derivatives; evaluate f on single-row (1, d) arrays so that the same vectorised model serves every method. Return the first-order mean, the second-order mean, the standard deviation, the gradient and the variance shares.

### Step 3: Hermite polynomials and the quadrature rule

```python
def gauss_hermite_e(n):
    J = np.diag(np.sqrt(np.arange(1, n)), 1)
    J = J + J.T                                     # Jacobi matrix of the He recurrence
    lam, V = np.linalg.eigh(J)
    return lam, V[0, :] ** 2                        # nodes, weights (sum to 1)
```

`hermite_e(n, x)` is the recurrence (4) with two running values.

### Step 4: the expansion

```python
def pce_projection(f, d, p, n_quad):
    nodes, weights = gauss_hermite_e(n_quad)
    grids = np.meshgrid(*([nodes] * d), indexing="ij")
    xi = np.stack([g.ravel() for g in grids], axis=1)                       # n_quad^d points
    w = np.prod(np.stack(np.meshgrid(*([weights] * d), indexing="ij")), axis=0).ravel()
    index = total_degree_multi_index(d, p)
    Psi = pce_basis(index, xi)                                              # He_alpha at every point
    coeffs = (Psi.T @ (w * np.asarray(f(xi)))) / pce_norms(index)          # (5)
    return coeffs, index
```

`pce_basis` multiplies the one-dimensional polynomials across the columns of ξ; `pce_norms` returns Π α_i!; `pce_moments` reads c₀ and Σ c_α² α!; `sobol_indices` sums the variance terms by which inputs they involve; `pce_regression` solves `np.linalg.lstsq(Psi, f(xi))` on seeded samples; `pce_evaluate` is `Psi @ coeffs`.

Running `main()` prints

```
linear (GUM) propagation: mean 1502.22 km (second order 1503.53), std 91.95 km, variance shares [0.24  0.667 0.093]
Monte Carlo (seed 17):
  N =     100: mean 1495.71 +- 9.42 km, std 94.16 km, P(R < required) = 0.0100 +- 0.0099
  N =  100000: mean 1503.67 +- 0.29 km, std 92.26 km, P(R < required) = 0.0443 +- 0.0007
  antithetic N = 10^4: mean 1503.65 +- 0.049 km (plain +- 0.93)
  RMS error of the MC mean over 200 replications: [8.64 4.35 2.22 1.21] km for N = [100, 400, 1600, 6400]; slope -0.474 (theory -0.5)
polynomial chaos (Hermite, tensor Gauss-Hermite with p + 1 points per dimension):
  degree 2:  10 terms,   27 model runs, mean 1503.535 km, std 92.137 km
  Sobol first-order (c_T, L/D, W_fuel): [0.241  0.6657 0.0924], total: [0.2417 0.6666 0.0927], sum of first-order 0.9991
  shortfall probability from 2 x 10^5 samples of the degree-4 surrogate: 0.0440 (normal approximation with the linear std: 0.0489)
```

and the figure: the histogram against the Gaussian and the surrogate, the N^{−1/2} line, and the coefficient decay by degree.

### Verify

- **Analytic cases:** the Breguet formula and its monotonicities; He₃ = ξ³ − 3ξ; the exponential's coefficients; a quadratic's exact expansion with mean 1.5 and variance 14.5 by hand; Sobol indices 1/19 and 18/19 for ξ₁ + 3ξ₂², and (0, 0)/(1, 1) for the pure interaction ξ₁ξ₂.
- **Limiting and failure cases:** an odd function has an antithetic mean of exactly zero with zero spread; the two-node rule is not exact for ξ⁴; too few regression samples raise; u outside (0, 1) raises.
- **Statistical checks with stated tolerances:** Box–Muller moments at 5 standard errors and a KS p-value above 10⁻³; MC mean within 4 standard errors of the PCE mean (false-alarm rate 6 × 10⁻⁵); sample std within 5 % (its own sd is 0.5 %); the bootstrap half-width within 10 % of 1.96 se; the surrogate's tail probability within 4 binomial standard errors of direct Monte Carlo.
- **Convergence:** MC RMS error slope −0.5 ± 0.12 and level σ/√N to 20 %; PCE standard deviation differences shrinking with degree and below 10⁻⁶ relative by degree 4.
- **Agreement with the oracles:** `hermeval`, `hermegauss`, `scipy.stats.norm`, `scipy.stats.expon`, `kstest`.

### Compare

`scipy.stats` provides every one-dimensional distribution with `rvs`, `ppf` and `cdf`, and `numpy.random.default_rng` draws normals by the ziggurat method, faster and with better tail behaviour than Box–Muller; use them. `numpy.polynomial.hermite_e.hermegauss` returns the same nodes as our Golub–Welsch to 10⁻¹⁰ (it uses the same idea, followed by a Newton polish). Libraries such as the one documented in `@chaospy-docs` wrap the whole of 2.5 (basis generation, quadrature, projection, regression, Sobol indices) for arbitrary input distributions through the generalised polynomial chaos (Legendre for uniform, Laguerre for gamma, and so on). What writing it bought is the knowledge that a PCE "mean" is a quadrature of the model against a polynomial, exact only to the degree the rule supports, and that a "99.99 % reliable" number from a surrogate is a statement about the polynomial's tail, not the model's, until checked.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment** on the assumed input distributions; write each prediction first.

1. **Monte Carlo convergence.** Run `mc_error_study` with N = 100, 400, 1600, 6400 and 200 replications; fit the slope; compare the RMS errors with σ/√N using the PCE standard deviation. Predict −0.5 and the four levels 9.2, 4.6, 2.3, 1.15 km.
2. **Cost and accuracy.** Tabulate mean, std and P(R < 1350 km) against model runs for: linear (7 runs), PCE degrees 1–4 (8, 27, 64, 125 runs) with the surrogate sampled 2 × 10⁵ times, and direct Monte Carlo at N = 27, 125, 10³, 10⁴, 10⁵. Predict that PCE at 27 runs beats direct Monte Carlo at 10⁴ on every quantity.
3. **Sensitivity of the answer to the assumptions.** Halve and double each input's relative standard deviation in turn and record P(R < 1350 km) from the degree-3 surrogate. Predict, from the Sobol indices, that the L/D spread moves the probability most.
4. **Variance reduction.** Compare the standard error of the mean for plain and antithetic sampling at N = 10³ and 10⁴, and the bootstrap interval of the 5th percentile of the range at N = 10⁴ against the surrogate's and the Gaussian approximation's 5th percentile. Predict that antithetic sampling helps the mean enormously and the percentile not at all.

Raw data format: CSV tables (`N, rms_error_km, sigma_over_sqrtN`; `method, runs, mean_km, std_km, p_short`; `input, factor, p_short`; `N, se_plain, se_antithetic`), each with a header naming the assumed distributions; `range_uq.png`.

What you should find: slope −0.47 and levels within 10 % of σ/√N; the degree-2 PCE (27 runs) giving mean 1503.5, std 92.1 km and p = 0.044 while 125 direct runs give the mean to ±8 km and the probability as "somewhere between 0.01 and 0.08"; the L/D factor moving p from 0.044 to about 0.013 (half) and 0.12 (double) while the fuel factor hardly moves it; antithetic standard errors 0.16 and 0.049 km against 2.9 and 0.93 km (a variance ratio near 350), and a 5th percentile of 1350 ± 4 km from the bootstrap against 1355 km from the surrogate and 1351 km from the Gaussian.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| nominal range | 1502.2 km | 1502.22 | 0 | — |
| linear sensitivities (km) | −45.07, 75.11, 27.98 | same to 10⁻⁶ | central differences | — |
| linear std | 91.95 km | 91.95 | 0 | — |
| mean shift | ½ tr H = c₍₂₀₀₎ + c₍₀₀₂₎ = 1.318 km | 1503.535 − 1502.217 = 1.318 | 0 | — |
| PCE std vs linear | slightly above | 92.14 vs 91.95 | 0.2 % | curvature and interactions |
| MC mean, N = 10⁵ | 1503.53 (PCE) | 1503.67 ± 0.29 | 0.5 se | sampling |
| MC std, N = 10⁵ | 92.14 | 92.26 | 0.1 % (its sd is 0.2 %) | sampling |
| MC order | −0.5 | −0.474 | within the replication noise | 200 replications |
| Gauss–Hermite vs NumPy | — | 10⁻¹⁰ | rounding | Newton polish in NumPy |
| P(R < 1350 km): MC / surrogate / Gaussian | — | 0.0443 ± 0.0007 / 0.0440 / 0.0489 | Gaussian 10 % high | skewness +0.10 |
| Σ first-order Sobol | ≤ 1 | 0.9991 | — | small interactions |

Three methods, three independent codes, one answer: mean 1503.5 km, standard deviation 92 km, shortfall probability 4.4 %, two thirds of the variance from the lift-to-drag ratio. That agreement validates the *arithmetic*. It says nothing about whether c_T really scatters by 3 % or L/D by 5 %, and the experiment in part 4 shows that the shortfall probability is three times more sensitive to the L/D assumption than to anything the propagation method could get wrong. The uncertainty statement must therefore read: "R = 1504 ± 92 km (one standard deviation), P(R < 1350 km) = 4.4 % ± 0.1 %, under assumed independent Gaussian inputs of 3 %, 5 % and 2 %; dominated by L/D (S = 0.67); 27 model runs for the PCE, 10⁵ for the Monte Carlo check". The thing to improve is not the method but the inputs: a drag polar from a wind-tunnel or flight test with its own uncertainty would replace the 5 % guess with a number, and that is where `P03.18` picks up.

Improvements worth making: a sparse grid for more inputs (extension 1); a correlation model (extension 2); robust optimisation on the surrogate (extension 3). Improvement worth *not* making: degree 6 PCE on this model. Degree 2 is converged to 10⁻⁶; the remaining error is in the assumptions, by a factor of a thousand.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of the Monte Carlo error, the Box–Muller transform, Hermite orthogonality, the projection and moment formulas and the hand sensitivities; `range_uq.png`; the four experiment tables; and a half-page uncertainty statement (template in `projects/_template/report-template.md`, sections 1, 5, 6 and 7) in the form of part 5, naming every assumption and the number of model runs.

### Adversarial self-checks

- **What if the spreads double?** All three relative standard deviations doubled: the linear std doubles to 184 km, the mean shift quadruples (it is quadratic in the spreads) to 5.3 km, and the shortfall probability rises from 4.4 % to about 20 %; the PCE needs degree 3 for 10⁻⁶ instead of degree 2 because the curvature now matters over the wider range. Check all three with the code.
- **What if a spread goes to zero?** With s₂ = 0 the model loses its dominant and only exactly linear input: the variance falls to a third, the skewness grows (the remaining inputs are the curved ones), and the Gaussian approximation gets worse, not better, although the spread shrank.
- **Which assumption fails first in the real world?** Independence of c_T and L/D (both depend on the installation and the cruise altitude), then the Gaussian shape of the fuel load (gauging errors are bounded, not Gaussian), then the fixed zero-fuel weight.
- **What measurement could falsify the model?** A fleet's recorded ranges on the same mission: a scatter far from 92 km falsifies the assumed input spreads; a mean far from 1504 km falsifies the model or its nominal inputs; both are routine operational data, and neither is produced by this lesson.
- **Does the result survive a change of method, seed or degree?** Linear, Monte Carlo (four seeds, five sample sizes) and PCE (degrees 1–4, projection and regression) agree to the tolerances stated in part 5. A disagreement beyond them would be a bug, with one exception: the Gaussian tail probability is systematically 10 % high, and that is the method, not a bug.
- **Which decision changes if the uncertainty doubles?** At 4.4 % shortfall an operator might accept the mission with a fuel stop as a contingency; at 20 % (doubled spreads) the mission would not be planned without more fuel, which `P03.16` says costs about 10 kg per 1 % of range. The decision hinges on the L/D spread, so the money goes into measuring the drag polar, not into more Monte Carlo samples.

## Check yourself

1. Write the Monte Carlo estimator of a mean and its standard error. How many samples are needed to halve the error?
2. A tail probability of 10⁻³ must be known to 10 % relative. How many direct Monte Carlo samples does that take, and what alternative does this lesson offer?
3. State the inversion method and explain why the Box–Muller transform is needed for the normal distribution.
4. Write the first-order propagation formula for independent inputs and the second-order correction to the mean. Which sign does the correction have for the range, and why?
5. State the orthogonality relation of the probabilists' Hermite polynomials and give He₂ and He₃.
6. How many Gauss–Hermite nodes are needed per dimension to compute the degree-p projection coefficients exactly for a polynomial model of degree p, and how many model runs is that in d dimensions?
7. Write the PCE mean and variance in terms of the coefficients, and define the first-order Sobol index of input i.
8. The lesson's Sobol indices sum to 0.9991. What does the missing 0.0009 represent?
