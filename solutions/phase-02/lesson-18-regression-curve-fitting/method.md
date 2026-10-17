# P02.18  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. x̄ = 3.5, ȳ = Σy/4, S_xx = Σ(x − 3.5)² = 5, S_xy = Σ(x − 3.5)(y − ȳ), S_yy likewise. 2. b = S_xy/S_xx, a = ȳ − b x̄. 3. Residuals eᵢ = yᵢ − a − b xᵢ; Σe = 0; s_e = √(Σe²/2); se_b = s_e/√S_xx. 4. r = S_xy/√(S_xx S_yy); t = r√2/√(1 − r²) against t₂. 5. CI: b ± 4.303 se_b; C_F A_t = b × 10⁻⁶ m². 6. `linear_regression(x, y)`.

## 2.
(a) ∂S/∂a = −2Σ(yᵢ − a − bxᵢ) = 0 and ∂S/∂b = −2Σxᵢ(yᵢ − a − bxᵢ) = 0; subtract x̄ times the first from the second. (b) The two normal equations *are* Σe = 0 and Σxe = 0; then Σ(y − ȳ)² = Σ(ŷ − ȳ)² + Σe² because the cross term Σ(ŷ − ȳ)e = bΣ(x − x̄)e = 0. (c) b = Σ(xᵢ − x̄)yᵢ/S_xx is a linear combination of the yᵢ with weights summing to zero, so Var b = σ²Σ(xᵢ − x̄)²/S_xx² = σ²/S_xx; ŷ(x₀) = ȳ + b(x₀ − x̄) with Cov(ȳ, b) = 0; add variances; replace σ by s_e and multiply by t*; for a new observation add σ². (d) First-order Taylor expansion g(θ̂) ≈ g(θ) + ∇gᵀ(θ̂ − θ); take the variance. Gradient components: I/a, I[ψ(b + 1) − ln c], −I(b + 1)/c.

## 3.
(a) Build α and C_L with the given formula and seed; `linear_regression`. (b) `residuals`; plot against α; look for curvature. (c) Mask α ≤ 8; refit; compare slope, se_slope, s_e. (d) Report the linear-range slope with its standard error and the range (α ≤ 8°) over which residuals show no pattern. (e) Apply `outlier_flags`; a smooth pattern spreads the misfit over many points so none exceeds 2 s_e.

## 4.
(a) `synthetic_thrust_curve(seed=s, noise=10 + s/2)`; `gauss_newton(thrust_model, thrust_model_jacobian, t, F, [1000, 0.5, 1.0])`; `np.sqrt(np.diag(cov))`; `curve_fit` with the same p0. (b) `total_impulse_model`, `impulse_uncertainty_delta`, `impulse_uncertainty_mc(..., seed=s)`, `total_impulse_trapezoid(t, F)`, `total_impulse_model(TRUE_THRUST_PARAMS)`. (c) `choose_degree(t, F, list(range(1, 11)))`; compare the LOO error with `gn["s_e"]`; the physical model has 3 parameters and a closed-form integral. (d) `nonlinear_band(thrust_model_jacobian, tt, params, cov, s_e, df)` with and without `prediction=True`; overlay `thrust_model(tt, TRUE_THRUST_PARAMS)`.

## 5.
(a) `F[8] += 10 * s`; `linear_regression`; `outlier_flags(..., k=2)` and `k=3`; residual/s_e. (b) Same with `F[0]`; compare |Δslope|; leverage hᵢᵢ = 1/n + (xᵢ − x̄)²/S_xx is largest at the ends. (c) `slope_confidence_interval` before and after removing the flagged point. (d) Expected count = 12 × P(|z| > 2) under normality; deleting automatically throws away real data about half the time you do it.

## 6.
(a) Confidence band for the mean response of the *fitted linear model*; outside 2–6 MPa it says nothing about the real motor (erosive burning, nozzle erosion, flow separation, transducer range). (b) Errors in x bias the slope *toward zero* by λ = S_xx/(S_xx + n σ_x²) (regression dilution); S_xx/n for 12 points on 2–6 MPa is 1.59 MPa², σ_x² = 0.01; simulate: 2000 fits with `p + rng.normal(0, 0.1, size=12)`. (c) Pointwise: at each x, 95 % of repeated experiments give a band containing the true mean there; a *simultaneous* band (e.g. Working–Hotelling) is wider.
