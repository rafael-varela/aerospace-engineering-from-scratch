# P02.18  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. Every fitted parameter carries a standard error; every band states its confidence level and whether it is for the mean response or a new observation. Label every simulated dataset as synthetic.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a four-point line by hand)

Four static fires (synthetic) at chamber pressures 2, 3, 4, 5 MPa give thrusts 590 + s, 885, 1175 and 1480 − s N.

(a) Compute x̄, ȳ, S_xx, S_xy, S_yy, the slope and the intercept. (b) Compute the four residuals, check that they sum to zero, and compute s_e and the standard error of the slope. (c) Compute r and r², and test the significance of r at α = 0.05 (df = 2). (d) Give the 95 % confidence interval for the slope and state what physical quantity the slope estimates (in m²) if F = C_F A_t p_c. (e) Confirm with `linear_regression`.

Answer form: seven numbers, four residuals, two statistics, a test, an interval with a physical reading.

## 2. (Derivation: least squares and its bands)

(a) Minimise S(a, b) = Σ(yᵢ − a − b xᵢ)² over a and b; show that the two normal equations give b = S_xy/S_xx and a = ȳ − b x̄.
(b) Show that the residuals sum to zero and are orthogonal to x, and explain in one sentence why this means r² = 1 − SS_res/SS_tot.
(c) Treating the xᵢ as fixed and the yᵢ as independent with variance σ², derive Var(b) = σ²/S_xx and Var(ŷ(x₀)) = σ²[1/n + (x₀ − x̄)²/S_xx]. Hence write the confidence band for the mean response and the prediction band for a new observation, and say where each is narrowest.
(d) State the delta method: for a smooth function g of fitted parameters with covariance C, Var g ≈ ∇gᵀ C ∇g. Apply it to the total impulse I(a, b, c) = a Γ(b + 1)/c^(b+1) to obtain the gradient used by `impulse_gradient`.

Answer form: four derivations.

## 3. (Unfamiliar scenario: a lift curve that bends)

A wind-tunnel run (synthetic) measures the lift coefficient at angles of attack α = 0, 2, …, 16°. Generate it with `rng = np.random.default_rng(s)` and C_L = 0.10 + 0.100 α − 0.0012 max(0, α − 8)² + rng.normal(0, 0.02, size=9), with α in degrees.

(a) Fit a line to all nine points with `linear_regression`; report the slope (per degree) with its standard error and r². (b) Plot the residuals against α. What pattern do you see, and what does it mean? (c) Refit using only α ≤ 8°; compare the slope and s_e with (a). (d) Which slope would you report as the lift-curve slope of the linear range, and how would you report the range of validity? (e) Would `outlier_flags` with k = 2 have caught the problem? Explain why a residual *pattern* is a different failure from an outlier.

Answer form: two fits with standard errors, a residual plot with a one-sentence reading, a recommendation, a short explanation.

## 4. (Numerical, uses your code: a noisier thrust curve)

Generate a thrust curve with `synthetic_thrust_curve(seed=s, n=60, noise=10 + s/2)`.

(a) Fit it with `gauss_newton` from p0 = (1000, 0.5, 1.0); report the three parameters with standard errors and the iteration count, and compare with `scipy.optimize.curve_fit`. (b) Give the total impulse from the fitted model with its delta-method and Monte Carlo (20 000 draws, seed s) uncertainties, and compare both with the trapezoid integral of the raw data and with the true value 619.5 N·s. (c) Run `choose_degree` over degrees 1–10 and report the chosen degree, the LOO error at that degree, and the residual s_e of the Gauss–Newton fit. Which model would you use to report the peak thrust and why? (d) Plot the fit with its 95 % confidence and prediction bands (`nonlinear_band`) and the true curve; is the true curve inside the confidence band everywhere?

Answer form: three parameters with standard errors, an impulse with two uncertainties and two comparisons, a degree with two errors and a choice, one plot with a yes/no.

## 5. (Numerical, uses your code: outliers and leverage)

Start from `synthetic_pressure_thrust(seed=s)`.

(a) Add 10 s N to the 9th point (index 8) and fit. Is it flagged by `outlier_flags` with k = 2? With k = 3? Report its residual in units of s_e. (b) Move the same fault to the first point (index 0, the lowest pressure) and repeat. Compare the change in slope in the two cases and explain the difference in terms of leverage (distance of x from x̄). (c) Remove the flagged point (if any) and refit; report the slope with its 95 % interval before and after. (d) The textbook rule flags residuals beyond about 2 s_e. With 12 clean normally distributed points, how many would you expect to be flagged by chance? What does that say about deleting flagged points automatically?

Answer form: two flag verdicts with residuals, two slope changes with an explanation, two intervals, a number with a one-sentence judgement.

## 6. (Judgment: what the fit cannot tell you)

(a) The lesson's pressure–thrust fit was made between 2 and 6 MPa. A colleague reads off 2372 ± 67 N at 8 MPa from the confidence band. Explain in three sentences what that interval does and does not promise, and name one physical reason the thrust might not follow the line at 8 MPa.
(b) The fit treats chamber pressure as exact, but the transducer has a noise of 0.1 MPa. In which direction does noise in x bias the estimated slope, and roughly by what factor here (compare the x-noise variance with S_xx/n)? Check by simulation with seed s: add N(0, 0.1²) to p before fitting, 2000 times, and report the mean fitted slope.
(c) A 95 % confidence band is drawn on the thrust curve. A reviewer says "so the true curve lies entirely inside the band with probability 0.95". Correct the statement.

Answer form: three sentences, a direction with a factor and a simulated mean slope, a corrected statement.
