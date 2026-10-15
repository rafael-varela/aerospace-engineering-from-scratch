# P02.18  Regression & Curve Fitting

> A fitted curve is a claim about the mean. The band around it says how far the mean might be; the wider band says how far the next measurement might be; and neither says anything outside the data.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the least-squares line from the normal equations, the orthogonality of its residuals, the variance of the slope and of a predicted mean, the confidence and prediction bands, the t-test for the correlation coefficient, and the delta-method variance of a function of fitted parameters. 2. Implement and test, from scratch, simple linear regression with bands and outlier flags, general linear least squares by QR with the parameter covariance, leave-one-out cross-validation for choosing a polynomial degree, a damped Gauss–Newton solver with delta-method bands, and the log-transform exponential fit, against `scipy.stats.linregress`, `numpy.linalg.lstsq` and `scipy.optimize.curve_fit` as oracles. 3. Measure, by seeded simulation, the coverage of a slope interval and of confidence and prediction bands, and how leave-one-out error separates fitting from over-fitting. 4. Decide the thrust coefficient and total impulse of a synthetic motor with their uncertainties, whether a flagged point is a fault, which model to report the peak thrust from, and what the fit does not promise beyond the data. |
| **Prerequisites** | `P01.17` (normal equations, QR, projection, parameter covariance), `P02.16` (variance of linear combinations, normal distribution), `P02.17` (t intervals, standard errors, hypothesis tests), `P00.03`, `P00.06`, `P00.07`, `P00.14`. Skill check: you can solve a 2×2 linear system, differentiate a sum of squares, and state what a standard error is. |
| **Engineering question** | From twelve synthetic static fires at different chamber pressures, what is the motor's thrust coefficient times throat area and how well is it known; from sixty noisy thrust samples, what are the peak thrust and total impulse with their uncertainties, and which fitted curve should the data sheet carry? |
| **Assumptions and limits** | Ordinary least squares: x exact, y with independent errors of equal variance (and normal, for the t-based bands). Fails when x carries noise (slope biased toward zero), when the scatter grows with y (weighted least squares is then needed), when errors are correlated in time (Lesson 19), and whenever the fitted form is extrapolated. The three-parameter thrust shape a t^b e^{−ct} is an empirical convenience, not combustion physics, and the pressure–thrust line is exact only while C_F is constant. All data are *synthetic*: generated from the stated models plus Gaussian noise with stated seeds. Gauss–Newton needs a start inside the basin of the minimum. |
| **Analytical work** | The normal equations and their solution; residual orthogonality and r² = 1 − SS_res/SS_tot; Var(b) = σ²/S_xx; confidence and prediction bands; the correlation t-test; the Vandermonde design matrix and β = R⁻¹Qᵀy; Cov(β) = s_e²(AᵀA)⁻¹; the leave-one-out residual e/(1 − h); the Gauss–Newton step from the linearised model; the delta method; the closed-form impulse a Γ(b + 1)/c^{b+1} and its gradient; the attenuation bias of x-noise. |
| **Computational deliverable** | `starter/main.py`: `synthetic_pressure_thrust`, `thrust_model`, `thrust_model_jacobian`, `synthetic_thrust_curve`, `linear_regression`, `predict_line`, `residuals`, `correlation_coefficient`, `correlation_t_test`, `slope_confidence_interval`, `confidence_band`, `prediction_band`, `outlier_flags`, `design_matrix_poly`, `linear_least_squares`, `general_band`, `polynomial_fit`, `loo_cv_error`, `choose_degree`, `gauss_newton`, `nonlinear_band`, `total_impulse_trapezoid`, `total_impulse_model`, `impulse_gradient`, `impulse_uncertainty_delta`, `impulse_uncertainty_mc`, `exponential_fit_by_log`. Tests prove: the line matches `linregress` to 10⁻¹²; residuals sum to zero and are orthogonal to x; an exact line gives r = 1; the slope interval and both bands cover at their nominal rate over hundreds of synthetic datasets; the prediction band exceeds the confidence band everywhere and both are narrowest at x̄; the injected fault is the only flagged point; QR least squares matches `lstsq` and its covariance matches s_e²(AᵀA)⁻¹; LOO picks the true degree of a synthetic cubic and never falls below the in-sample RMS; the Jacobian matches finite differences; Gauss–Newton matches `curve_fit` and solves a linear model in one step; the impulse closed form matches quadrature and its delta-method uncertainty matches Monte Carlo. |
| **Practical activity** | **computational experiment**: coverage of intervals and bands over synthetic replicates, the residual-against-LOO curve over polynomial degree, the leverage of an injected fault, and the attenuation of the slope by noise in x. All data are synthetic. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `thrust_curve_fit.png`, the hand derivations, the four tables from part 4, and a one-page motor data-sheet note giving C_F A_t with its interval, the peak thrust and total impulse with uncertainties, the fitted form used and its range of validity. |
| **Assessment** | `assessments/phase-02/grading-rubric.md`; competencies 1–4. Pass: derivations complete including the bands and the delta method; all public tests pass; the hand fit of problem 1 agrees with the code to three figures; the data-sheet note states which band is which, flags the extrapolation, and gives the impulse uncertainty by two methods. |
| **References** | `@openstax-introductory-statistics-2e` §12.1 "Linear Equations", §12.2 "Scatter Plots", §12.3 "The Regression Equation", §12.4 "Testing the Significance of the Correlation Coefficient", §12.5 "Prediction", §12.6 "Outliers" (all required-free). `@mit-ocw-18-06-linear-algebra` Lecture 16 "Projection matrices and least squares" (required-free), Lecture 17 "Orthogonal matrices and Gram-Schmidt" (supplementary-free; QR for the general fit). `@openstax-precalculus-2e` §4.8 "Fitting Exponential Models to Data" (supplementary-free; the log-transform fit). Confidence and prediction bands for general linear models, cross-validation, Gauss–Newton and the delta method are not covered by these sources; parts 2.4–2.7 are original instruction. |
| **Extensions** | ★ Implement weighted least squares for a thrust transducer whose noise grows with thrust (σ ∝ F) and show how the bands change. ★ Fit the thrust curve with a Levenberg–Marquardt damping term and compare the basin of convergence with plain Gauss–Newton from twenty random starts. ★ Derive the Working–Hotelling simultaneous band and compare its width with the pointwise band at the ends of the data. |
| **Time** | unmeasured |

## 1. The engineering problem

A small solid motor is static-fired twelve times at chamber pressures from 2 to 6 MPa, and the thrust is logged. Rocket theory says F = C_F A_t p_c: thrust is proportional to chamber pressure, with the thrust coefficient C_F and throat area A_t as the constant of proportionality. The data (synthetic: 290 N/MPa plus 25 N of Gaussian noise, seed 18) do not lie on a line; they scatter around one. The questions are: what is the slope, how well is it known, what thrust should we expect at a pressure we have not fired, and how sure are we of that? One of the twelve points sits well above the others; is it a logging fault or a real event?

The same motor's thrust–time trace (synthetic: a three-parameter curve with 10 N of noise, 60 samples over 4 s) has no straight line in it. The data sheet needs the peak thrust and the total impulse ∫F dt, with uncertainties, and a fitted curve that a simulation can use. There are two routes: a polynomial of some degree, which is linear in its coefficients and fits by the machinery of P01.17, and a physically shaped model F = a t^b e^{−ct}, which is not, and needs an iterative solver. Both routes raise the same question, which this lesson answers with a number: how do you know you are fitting the motor and not the noise?

Success criteria, written down first: the from-scratch line matches `scipy.stats.linregress` to 10⁻¹²; over 500 synthetic datasets the 95 % slope interval contains the true slope 95 % ± 4 % of the time; the injected fault is flagged and nothing else is; Gauss–Newton matches `scipy.optimize.curve_fit` to 10⁻⁶ in the parameters; the delta-method impulse uncertainty matches a 20 000-draw Monte Carlo within 15 %; and leave-one-out cross-validation picks degree 3 for a synthetic cubic.

## 2. Understand and derive

### 2.1 The least-squares line

For n pairs (xᵢ, yᵢ) the **regression line** ŷ = a + b x minimises S(a, b) = Σ(yᵢ − a − bxᵢ)² (`@openstax-introductory-statistics-2e` §12.3; `@mit-ocw-18-06-linear-algebra` Lecture 16 gives the same result as a projection). Setting ∂S/∂a = 0 and ∂S/∂b = 0 gives the two **normal equations**, Σeᵢ = 0 and Σxᵢeᵢ = 0 with eᵢ = yᵢ − ŷᵢ the **residuals**. Their solution, with S_xx = Σ(xᵢ − x̄)², S_xy = Σ(xᵢ − x̄)(yᵢ − ȳ):

**b = S_xy/S_xx**,   **a = ȳ − b x̄**.   (1)

The line passes through (x̄, ȳ). The residuals sum to zero and are orthogonal to x, so the total scatter splits exactly: Σ(yᵢ − ȳ)² = Σ(ŷᵢ − ȳ)² + Σeᵢ², and the **coefficient of determination**

**r² = 1 − SS_res/SS_tot**   (2)

is the fraction of the variance of y that the line accounts for. Its square root, with the sign of b, is the **correlation coefficient** r = S_xy/√(S_xx S_yy). The residual standard deviation is s_e = √(SS_res/(n − 2)): two parameters were fitted, so two degrees of freedom are gone.

For the twelve static fires: b = **299.96 N/MPa**, a = −27.9 N, r = 0.99719, r² = 0.99438, s_e = 31.0 N (the noise was 25 N; the estimate scatters). Dividing the slope by 10⁶ Pa/MPa gives C_F A_t = 3.000 × 10⁻⁴ m² against the synthetic truth 2.900 × 10⁻⁴.

### 2.2 How well is the slope known?

Treat the xᵢ as exact and the yᵢ as independent with variance σ². The slope b = Σ(xᵢ − x̄)yᵢ/S_xx is a linear combination of the yᵢ, so

**Var b = σ² Σ(xᵢ − x̄)²/S_xx² = σ²/S_xx**,   **se_b = s_e/√S_xx**.   (3)

For the fires se_b = 7.13 N/MPa, and the 95 % interval is b ± t*₁₀ se_b = (284.1, 315.9) N/MPa: the true 290 is inside. The spread of the x values is in the denominator: firing at pressures far apart buys precision in the slope. The intercept's standard error is s_e√(1/n + x̄²/S_xx) = 29.9 N, and the interval for it, −27.9 ± 66.6 N, contains the theoretical zero.

**Is there a relationship at all?** Under H0: the population correlation is zero, t = r√(n − 2)/√(1 − r²) follows t_{n−2} (§12.4). Here t = 42.1, p = 1.4 × 10⁻¹²: the question is not in doubt for this dataset, and the test matters only for weak correlations.

### 2.3 Two bands, two questions

The fitted mean at x₀ is ŷ(x₀) = ȳ + b(x₀ − x̄). Because Σ(xᵢ − x̄) = 0, ȳ and b are uncorrelated, and

**Var ŷ(x₀) = σ²[1/n + (x₀ − x̄)²/S_xx]**.   (4)

The **confidence band** for the *mean response* is ŷ ± t* s_e √(1/n + (x₀ − x̄)²/S_xx); it is narrowest at x̄ (where it equals t* s_e/√n) and widens as a hyperbola away from the data. A *new measurement* at x₀ scatters about the mean by a further σ², so the **prediction band** is ŷ ± t* s_e √(1 + 1/n + (x₀ − x̄)²/S_xx) (§12.5). For the fires:

| p_c (MPa) | ŷ (N) | 95 % confidence ± | 95 % prediction ± |
|---|---|---|---|
| 2 | 572 | 37.5 | 78.6 |
| 4 (= x̄) | 1172 | 19.9 | 71.9 |
| 6 | 1772 | 37.5 | 78.6 |
| 8 | 2372 | 66.6 | 96.0 |

The row at 8 MPa is an **extrapolation**: the arithmetic is the same, the meaning is not. The band says where the *line* would be; it does not know whether the motor's C_F stays constant at a pressure never fired (erosive burning, nozzle erosion and transducer range all say it might not). Report it with the words "beyond the data".

Both bands are **pointwise**: at each x₀, 95 % of repeated experiments give a band containing the true mean *there*. The probability that the whole true line lies inside the whole band is lower; a simultaneous band is wider (extension 3).

### 2.4 Outliers and leverage

A point whose residual exceeds about two residual standard deviations is suspect (§12.6 uses 1.9 s; we use 2 s_e). With a +120 N logging fault injected at the ninth fire, the slope moves to 305.7 N/MPa, s_e inflates to 51.8 N, and the fault's residual is 124.0 N = 2.40 s_e: flagged, and nothing else is. Refitting without it gives 298.9 N/MPa and s_e = 31.9 N, back where the clean data were. Two cautions. First, a point at the end of the x range has high **leverage** (hᵢᵢ = 1/n + (xᵢ − x̄)²/S_xx): it pulls the line toward itself and its own residual comes out *smaller* than the fault that caused it, so faults at the ends are the hardest to flag (problem 5). Second, with twelve clean normal points the expected number beyond 2 s_e by chance is 12 × 0.046 = 0.55: flags are prompts to look at the logging record, not deletion commands.

### 2.5 General linear least squares: polynomials by QR

Any model that is linear in its parameters, y ≈ Σ βⱼ φⱼ(x), is a least-squares problem A β ≈ y with design matrix Aᵢⱼ = φⱼ(xᵢ); for a polynomial of degree d, A is the Vandermonde matrix [1, x, x², …, x^d]. P01.17 derived the solution two ways; the stable one is A = QR, β = R⁻¹Qᵀy, which avoids squaring the condition number (`@mit-ocw-18-06-linear-algebra` Lecture 17). The parameter covariance is

**Cov β = s_e² (AᵀA)⁻¹ = s_e² R⁻¹R⁻ᵀ**,   s_e² = SS_res/(n − p),   (5)

and the bands generalise by the same argument as (4): for a new row **a** of the design matrix, Var(**a**ᵀβ) = **a**ᵀ Cov β **a**, plus s_e² for a prediction. `general_band` does this for any linear model and the test checks that it reproduces 2.3 for d = 1.

**Which degree?** Fitting the 60-sample thrust trace with polynomials of degree 1 to 10:

| degree | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 |
|---|---|---|---|---|---|---|---|---|
| residual s_e (N) | 59.6 | 58.4 | 45.7 | 26.3 | 15.8 | 9.60 | 9.29 | 9.40 |
| leave-one-out RMS (N) | 61.5 | 63.1 | 51.4 | 30.4 | 19.1 | 10.46 | 10.08 | 10.39 |

The residual s_e always falls (or nearly) with degree: more parameters can only fit the data better, and past degree 6 it sits at the 10 N noise floor, fitting noise. **Leave-one-out cross-validation** asks the right question: fit without point i, predict yᵢ, repeat for every i, and take the RMS of those *prediction* errors. The LOO residual is eᵢ/(1 − hᵢᵢ), never smaller than the ordinary one, and it turns up again when the degree is too high (10.39 at degree 10 against 10.08 at degree 8). `choose_degree` picks 8 here, with 6 within 4 %; the pick between 6, 8 and 10 is inside the noise, and a polynomial of degree 8 is a poor thing to put on a data sheet: it has no physical parameters, it oscillates between samples, and beyond 4 s it diverges. The polynomial route tells us the noise level (about 10 N) and that six or more coefficients are needed; the physical route below needs three.

### 2.6 Nonlinear least squares: Gauss–Newton

The thrust shape F(t; a, b, c) = a t^b e^{−ct} rises as a power, peaks at t = b/c and tails off exponentially. It is linear in a but not in b or c, so (1)–(5) do not apply directly. **Gauss–Newton** linearises: near the current parameters **p**, F(t; **p** + δ) ≈ F(t; **p**) + J δ with J the Jacobian ∂F/∂**p** (an n × 3 matrix, here analytic: columns F/a, F ln t, −tF). The step δ that minimises |**r** − Jδ|² with **r** = **y** − F(**p**) is a *linear* least-squares problem, solved by QR as in 2.5, and the process repeats. Two safeguards make it robust: halve the step until the sum of squares actually decreases (damping), and stop when the relative step falls below a tolerance. At the solution the same linearisation gives the parameter covariance s_e²(JᵀJ)⁻¹ and, through `general_band` on J, the confidence and prediction bands (the **delta method**: Var g(**p̂**) ≈ ∇gᵀ Cov **p̂** ∇g for any smooth g).

From the start (1000, 0.5, 1.0), seven iterations give a = 1653 ± 54 N s^{−0.833}, b = 0.833 ± 0.022, c = 1.646 ± 0.030 s⁻¹, s_e = 9.77 N; `scipy.optimize.curve_fit` returns the same parameters to six figures and the same standard errors to three. The truth (1550, 0.8, 1.6) lies within 1.9, 1.5 and 1.5 standard errors; the three estimates are strongly correlated (a larger a is compensated by a larger c), which is why individual parameters look further off than the curve does. Peak thrust: t = b/c = 0.506 s, F = 407.4 N.

Gauss–Newton is fast (quadratic convergence near the minimum for small residuals) and not safe: from a start outside the basin it diverges, and a model with a flat direction gives a singular J. `curve_fit` uses Levenberg–Marquardt, which adds a damping term to JᵀJ and trades speed for a wider basin (extension 2); it is the right tool in practice, and we wrote Gauss–Newton to see what "the covariance of the fit" means operationally: the inverse of JᵀJ at the solution, scaled by the residual variance.

### 2.7 Total impulse and its uncertainty

The data sheet wants I = ∫₀^∞ F dt. Three estimates:

1. **Trapezoid of the raw samples**, 614.8 N·s: honest about the data but misses the first 0.05 s and the tail beyond 4 s, and carries the noise.
2. **Closed form of the fitted model**: ∫₀^∞ a t^b e^{−ct} dt = a Γ(b + 1)/c^{b+1} = **623.8 N·s**, which uses the model to fill the gaps.
3. The truth from the generating parameters: 619.5 N·s.

Uncertainty of (2) by the **delta method**: ∇I = (I/a, I[ψ(b + 1) − ln c], −I(b + 1)/c) with ψ the digamma function, and σ_I = √(∇Iᵀ Cov ∇I) = **4.9 N·s**. By Monte Carlo (20 000 parameter draws from N(**p̂**, Cov), seed 18): 4.9 N·s. The two agree because the model is nearly linear in its parameters over the width of the covariance; when they disagree the Monte Carlo is right and the fit is poorly conditioned. The model's 623.8 is 0.9σ_I from the truth.

### 2.8 The log-transform fit, and its trap

A pure exponential tail F = F₀e^{−t/τ} becomes a straight line in ln F = ln F₀ − t/τ, which `linear_regression` fits directly (`@openstax-precalculus-2e` §4.8). On the window 1.0 ≤ t ≤ 2.5 s (thrust 60–310 N) it gives τ = 0.893 s with r² = 0.991 in log space, bracketed by the fitted model's local decay time 1/(c − b/t), 1.23 s at 1.0 s and 0.76 s at 2.5 s. On t > 2.5 s, where the thrust is below 60 N and the noise is 10 N, the same fit gives τ = 0.479 s with r² = 0.575: the logarithm turns 10 N of noise on a 15 N reading into a factor-two error, and the fit is wrecked. Minimising squared errors in ln F weights the smallest values most heavily, which is the opposite of what a thrust transducer deserves. Use the transform to find a starting point, then fit in the original units.

### 2.9 Where the model stops being valid

- **x exact.** Transducer noise in the chamber pressure biases the slope toward zero by the factor S_xx/(S_xx + nσ_x²) (regression dilution); for 0.1 MPa noise on this dataset that is 0.994, negligible, but for noisy x it is not (problem 6).
- **Equal variances.** If the thrust noise grows with thrust, the high points deserve less weight: weighted least squares (extension 1).
- **Independent errors.** A thrust trace sampled at 1 kHz has correlated noise; sixty samples at 15 Hz are nearly independent here by construction. Lesson 19 treats correlated noise.
- **The fitted form.** A line in p_c assumes constant C_F; the three-parameter shape assumes one rise and one decay. A residual plot with a *pattern* (problem 3) says the form is wrong, and no band can repair that.
- **Inside the data.** Every band is conditional on the model and valid between the smallest and largest x.

## 3. Build the computational model

Open `starter/main.py`. Build in the order below.

### Step 1: the line

```python
def linear_regression(x, y):
    xm, ym = x.mean(), y.mean()
    dx, dy = x - xm, y - ym
    sxx, syy, sxy = dx @ dx, dy @ dy, dx @ dy
    slope = sxy / sxx; intercept = ym - slope * xm                                  # (1)
    resid = y - (intercept + slope * x); ss_res = resid @ resid
    s_e = math.sqrt(ss_res / (n - 2))
    return {"slope": slope, "intercept": intercept, "r": sxy / math.sqrt(sxx * syy),
            "r2": 1.0 - ss_res / syy, "s_e": s_e, "se_slope": s_e / math.sqrt(sxx),
            "se_intercept": s_e * math.sqrt(1.0 / n + xm ** 2 / sxx), "n": n, "df": n - 2, "x_mean": xm, "sxx": sxx}
```

Raise `ValueError` for fewer than three points or no spread in x. `confidence_band` and `prediction_band` return half-widths from (4); `outlier_flags` compares |residual| with k s_e; `correlation_t_test` is 2.2.

### Step 2: general linear least squares

```python
def linear_least_squares(A, y):
    Q, R = np.linalg.qr(A)
    beta = np.linalg.solve(R, Q.T @ y)
    resid = y - A @ beta
    s_e = math.sqrt(resid @ resid / (n - p))
    Rinv = np.linalg.solve(R, np.eye(p))
    return {"beta": beta, "cov": s_e ** 2 * Rinv @ Rinv.T, "s_e": s_e, "df": n - p, "residuals": resid}   # (5)
```

Check for n ≤ p and for a tiny diagonal entry of R (rank deficiency). `general_band` is t* √(**a**ᵀ Cov **a** [+ s_e²]) per row, with `np.einsum("ij,jk,ik->i", A_new, cov, A_new)`. `loo_cv_error` loops over i with a boolean mask; `choose_degree` takes the argmin.

### Step 3: Gauss–Newton

```python
def gauss_newton(model, jac, x, y, p0, tol=1e-10, max_iter=100):
    p = np.array(p0, float); r = y - model(x, p); ss = r @ r
    for it in range(1, max_iter + 1):
        Q, R = np.linalg.qr(jac(x, p)); dp = np.linalg.solve(R, Q.T @ r)
        step = 1.0
        while True:                                     # damping: halve until the sum of squares falls
            p_new = p + step * dp; r_new = y - model(x, p_new); ss_new = r_new @ r_new
            if ss_new <= ss or step < 1e-8: break
            step *= 0.5
        converged = np.linalg.norm(step * dp) <= tol * max(1.0, np.linalg.norm(p_new))
        p, r, ss = p_new, r_new, ss_new
        if converged: ...                               # covariance s_e^2 (J^T J)^{-1} from a final QR of J
    raise RuntimeError("Gauss-Newton did not converge")
```

Return the parameters, covariance, s_e, df, iteration count and SS_res. `nonlinear_band` hands the Jacobian rows at new t to `general_band`. `thrust_model_jacobian` is analytic; the test compares it with central differences.

### Step 4: impulse and the exponential tail

`total_impulse_trapezoid` is Σ ½(Fᵢ + Fᵢ₊₁)Δtᵢ written out; `total_impulse_model` is a Γ(b + 1)/c^{b+1} with `math.gamma`; `impulse_gradient` uses `scipy.special.digamma`; `impulse_uncertainty_delta` is √(gᵀ Cov g); `impulse_uncertainty_mc` draws from `rng.multivariate_normal` and takes the standard deviation of the closed form. `exponential_fit_by_log` fits `linear_regression(t, np.log(F))` and returns F₀ = e^a, τ = −1/b, raising for non-positive F or a non-decaying slope.

Running `main()` prints the tables of 2.1–2.8 and writes `thrust_curve_fit.png`: the pressure–thrust line with both bands and the 6 MPa data edge marked, the thrust curve with its Gauss–Newton fit, bands and the true curve, and the residual-against-LOO plot over degree.

### Verify

Run `python -m pytest curriculum/phase-02/lesson-18-regression-curve-fitting`. The tests:

1. **Analytic cases.** Slope, intercept, r, standard errors and p match `linregress` to 10⁻¹²; an exact line gives r = ±1 and s_e = 0; an exact cubic is recovered by `polynomial_fit` to 10⁻¹⁰; a noiseless exponential is recovered by the log fit to 10⁻¹⁰; the impulse closed form matches the trapezoid of the exact curve to 10⁻⁴ and reduces to a/c for b = 0; the analytic Jacobian and impulse gradient match central differences.
2. **Conservation and structure.** Residuals sum to zero and are orthogonal to x (and to every column of A); r² = 1 − SS_res/SS_tot = r²; the slope is invariant to shifting x and scales correctly with y and with the units of x; degree-1 `polynomial_fit` equals `linear_regression`; `general_band` on the degree-1 design equals the two line bands.
3. **Oracles.** `linear_least_squares` matches `lstsq` and its covariance matches s_e²(AᵀA)⁻¹; Gauss–Newton matches `curve_fit` to 10⁻⁶ in the parameters and 10⁻³ in the covariance, and solves a linear model in at most two iterations.
4. **Coverage (seeded, 4 SE tolerances in comments).** The 95 % slope interval covers 290 N/MPa 95 % ± 4 % of the time over 500 synthetic datasets; the confidence band at 5 MPa covers the true mean, and the prediction band a fresh observation, 95 % ± 4.4 % over 400 datasets.
5. **Model selection and faults.** LOO picks degree 3 for a synthetic cubic and never falls below the in-sample RMS; the injected +120 N fault is the only point flagged at k = 2 and nothing is flagged at k = 3; the fault inflates s_e.
6. **Uncertainty propagation.** Delta-method and Monte Carlo impulse uncertainties agree within 15 %; the fitted impulse is within 3σ of the truth; the prediction band exceeds the confidence band everywhere.
7. `main()` writes `thrust_curve_fit.png` to the directory it is given.

### Compare

`scipy.stats.linregress` gives the line with its standard errors and p-value in one call; `numpy.polyfit(cov=True)` gives polynomial coefficients with covariance; `scipy.optimize.curve_fit` gives nonlinear fits with covariance by Levenberg–Marquardt, with bounds, with weights (`sigma=`), and with a far wider basin than our Gauss–Newton. All three agree with our code on these data to the stated tolerances. We wrote our own to see three things: that every band is **a**ᵀ Cov **a** for some row **a**, whether the model is a line, a polynomial or a linearised nonlinear curve; that the covariance of a nonlinear fit is a linearisation, trustworthy only when the delta method and the Monte Carlo agree; and that `curve_fit` returning without complaint is not evidence that it found the right minimum, which only a residual plot and a second start can provide.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All data are synthetic, generated by the stated models with stated seeds.

1. **Coverage.** Over 500 datasets `synthetic_pressure_thrust(seed=1000 + k)`, count how often the 95 % slope interval contains 290 N/MPa and how often the confidence band at 5 MPa contains 1450 N. Predict 0.95 ± 0.01 for both.
2. **Residual against LOO.** Tabulate s_e and the LOO RMS for degrees 1–10 on `synthetic_thrust_curve(seed=18)` (table in 2.5) and plot both against degree on log axes. Mark the 10 N noise floor. Where does s_e cross it, and where does LOO turn up?
3. **Leverage.** Inject a +120 N fault at index 8 (p_c = 4.9 MPa) and then at index 0 (2.0 MPa) and index 11 (6.0 MPa); record the slope shift and the fault's residual in units of s_e each time. Predict first which position is flagged and which moves the slope most.
4. **Attenuation.** Add N(0, σ_x²) noise to the pressures with σ_x = 0.1, 0.3, 0.6 MPa and refit 2000 times each (seed 18); record the mean fitted slope against the prediction 290 × S_xx/(S_xx + nσ_x²).

Raw data format: CSV tables `quantity, coverage`; `degree, s_e, loo_rms`; `index, slope_shift, residual_over_se, flagged`; `sigma_x, mean_slope, predicted_slope`.

What you should find: coverage within one or two standard errors of 0.95; s_e reaches the noise floor at degree 6 while LOO bottoms at 6–8 and rises after; the end-point faults move the slope most and are flagged least (their residual is diluted by leverage); and the fitted slope falls below 290 by the predicted factor, 0.994 at 0.1 MPa but 0.81 at 0.6 MPa.

## 5. Validate and improve

| quantity | hand / theory (§2) | code / library / simulation | difference | cause |
|---|---|---|---|---|
| slope (N/MPa) | 290 (generating value) | 299.96 ± 7.13 = `linregress` | 1.4 se | sampling (12 points, 25 N noise) |
| intercept (N) | 0 (generating value) | −27.9 ± 29.9 | 0.9 se | sampling |
| band at x̄ | t* s_e/√n = 2.228 × 31.0/√12 = 19.9 N | 19.9 N | 0 | same formula |
| slope-interval coverage | 0.95 | 0.95 ± 0.01 (500 datasets) | sampling | — |
| fault residual | +120 N injected | 124.0 N = 2.40 s_e, flagged | leverage 0.03 | the fit absorbs a little |
| Gauss–Newton vs `curve_fit` | — | parameters to 10⁻⁶, sd to 10⁻³ | < stated tolerance | different stopping rules |
| impulse (N·s) | 619.5 (generating) | 623.8 ± 4.9 (delta) ± 4.9 (MC) | 0.9σ | parameter sampling |
| impulse (N·s) | 619.5 | 614.8 (trapezoid of data) | −0.8 % | missing first 0.05 s and tail beyond 4 s |
| tail decay τ (s) | 0.76–1.23 (local, model) | 0.893 (log fit, 1–2.5 s) | inside the range | decay time varies along the window |

Every difference is within sampling error or has a named cause. The parameters a, b, c individually sit 1.5–1.9 standard errors from the generating values while the curve itself, the peak and the impulse sit within one: the three estimates are correlated (the off-diagonal entries of Cov are not small), and a quantity that depends on all three can be better known than any one of them. The delta method handles that automatically through the full covariance; quoting three independent ± values would not.

What the agreement does not show: that a real motor follows either model. The pressure–thrust line and the a t^b e^{−ct} shape were *the generating models*, so recovering them proves the solver, not the physics. Real thrust traces have ignition transients, chuffing and a tail-off whose shape depends on the grain geometry, and a residual plot is the first thing to look at before any band is drawn.

Improvements worth making now: weighted least squares for a thrust-proportional noise model (extension 1); a Levenberg–Marquardt damping term and a basin-of-convergence map (extension 2); a residual-pattern test (the runs test, or a fit of the residuals against a quadratic) to catch a wrong model form before fitting bands to it. Improvement worth *not* making: reporting the degree-8 polynomial because it has the lowest LOO error. It predicts the samples marginally better and the physics not at all.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivations of (1)–(4) and of the delta-method gradient; `thrust_curve_fit.png`; the four tables from part 4; and a one-page motor data-sheet note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving C_F A_t with its 95 % interval, the peak thrust and total impulse with uncertainties by two methods, the fitted form and its range of validity, a statement of which band is which, and the one extrapolation you refused to make.

### Adversarial self-checks

- **What if the noise doubles?** s_e doubles, every band and standard error doubles, r² falls (from 0.994 to 0.979 for this x range), the slope point estimate changes only by sampling; the number of fires needed for the same slope precision quadruples. Check with `synthetic_pressure_thrust(noise=50.0)`.
- **What if the pressure range shrinks to 3.5–4.5 MPa?** S_xx falls 16-fold and se_b rises 4-fold: se_b becomes 28.5 N/MPa instead of 7.1 and the intercept is almost unconstrained. Precision in a slope comes from spread in x, not from more points at the same x.
- **Which assumption fails first in the real world?** Constant C_F (the fitted form) at the ends of the pressure range, and equal variances for the thrust trace (transducer noise often scales with signal). Independence of errors fails next for densely sampled traces.
- **What measurement could falsify the model?** A residual plot against pressure with curvature, or a fire at 7 MPa landing outside the prediction band: the band *is* the falsifiable statement. For the thrust curve, a second motor from the same lot whose trace lies outside the prediction band in the tail-off.
- **Does the result survive a change of method?** Line: `linregress` and QR agree to 10⁻¹². Curve: Gauss–Newton and `curve_fit` agree from the same start; from the start (3000, 2.0, 4.0) Gauss–Newton still converges here but a start with b < 0 fails (the model is not defined at t = 0 for b < 0): say which starts you allow. Impulse: closed form, trapezoid of the fitted curve on a fine grid and Monte Carlo agree to 0.1 %.
- **Which decision would change if the uncertainty doubled?** With σ_I = 9.8 N·s instead of 4.9 the impulse interval 624 ± 20 N·s would straddle a motor-class boundary set at, say, 640 N·s, and the data sheet could no longer place the motor in a class without more fires. The *decision* (how many fires to buy) follows from the band width, not from the point estimate.

## Check yourself

1. Write the two normal equations for the line and say what each one states about the residuals.
2. r² = 0.994 for the pressure–thrust fit. What does that number mean, and why does it not by itself tell you the line is the right model?
3. Give Var(b) in terms of σ and the x values, and say how to halve the standard error of the slope without reducing the noise.
4. State the difference between the confidence band and the prediction band, and say which one a range-safety engineer needs for a single future firing.
5. Why is leave-one-out error a better guide to polynomial degree than the residual standard deviation?
6. What does one Gauss–Newton step solve, and why can the solution's covariance be written s_e²(JᵀJ)⁻¹?
7. State the delta method and apply it to g = 1/c for a fitted parameter c with standard error σ_c.
8. The log-transform fit of a decaying tail minimises squared errors in ln F. Which data points does that weight most heavily, and when does that go wrong?
