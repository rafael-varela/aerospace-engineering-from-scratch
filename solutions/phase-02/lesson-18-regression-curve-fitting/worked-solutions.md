# P02.18  Worked solutions

Numbers below use **s = 15**. Yours differ if your seed differs; the structure of the reasoning should not. Every number here was produced by `solution.py` or by the short scripts quoted in each problem. All data are synthetic.

## 1. Four-point line by hand

Thrusts: 605, 885, 1175, 1465 N at 2, 3, 4, 5 MPa.

(a) x̄ = 3.5, ȳ = 1032.5. S_xx = 2.25 + 0.25 + 0.25 + 2.25 = 5. S_xy = (−1.5)(−427.5) + (−0.5)(−147.5) + (0.5)(142.5) + (1.5)(432.5) = 1435. S_yy = 411 875. **b = 1435/5 = 287.0 N/MPa**, **a = 1032.5 − 287 × 3.5 = 28.0 N**.
(b) Residuals: 605 − 602 = 3, 885 − 889 = −4, 1175 − 1176 = −1, 1465 − 1463 = 2; sum 0. s_e = √(30/2) = 3.873 N; se_b = 3.873/√5 = 1.732 N/MPa.
(c) r = 1435/√(5 × 411 875) = 0.999964, r² = 0.999927. t = r√2/√(1 − r²) = 165.7 with 2 df, p = 3.6 × 10⁻⁵: significant at 0.05.
(d) t*₂ = 4.303; 95 % CI = 287.0 ± 4.303 × 1.732 = **(279.5, 294.5) N/MPa**. With F = C_F A_t p_c, the slope in N/Pa is C_F A_t = 2.87 × 10⁻⁴ m² (for a C_F of 1.45, a throat area of 1.98 cm²).
(e) `linear_regression(np.array([2.,3.,4.,5.]), np.array([605.,885.,1175.,1465.]))` returns the same slope, intercept, r and standard errors to machine precision.

## 2. Derivations

(a) S = Σ(yᵢ − a − bxᵢ)². ∂S/∂a = −2Σ(yᵢ − a − bxᵢ) = 0 gives ȳ = a + b x̄. ∂S/∂b = −2Σxᵢ(yᵢ − a − bxᵢ) = 0 gives Σxᵢyᵢ = aΣxᵢ + bΣxᵢ². Substitute a = ȳ − b x̄: Σxᵢyᵢ − n x̄ ȳ = b(Σxᵢ² − n x̄²), i.e. S_xy = b S_xx.
(b) The two normal equations read Σeᵢ = 0 and Σxᵢeᵢ = 0 directly (eᵢ = yᵢ − a − bxᵢ). Then yᵢ − ȳ = (ŷᵢ − ȳ) + eᵢ, and the cross term Σ(ŷᵢ − ȳ)eᵢ = bΣ(xᵢ − x̄)eᵢ = b(Σxᵢeᵢ − x̄Σeᵢ) = 0, so SS_tot = SS_reg + SS_res and r² = SS_reg/SS_tot = 1 − SS_res/SS_tot.
(c) b = Σ(xᵢ − x̄)(yᵢ − ȳ)/S_xx = Σ(xᵢ − x̄)yᵢ/S_xx because Σ(xᵢ − x̄)ȳ = 0. The yᵢ are independent with variance σ², so Var b = σ²Σ(xᵢ − x̄)²/S_xx² = σ²/S_xx. ŷ(x₀) = ȳ + b(x₀ − x̄); Cov(ȳ, b) = (σ²/n)Σ(xᵢ − x̄)/S_xx = 0, so Var ŷ(x₀) = σ²/n + (x₀ − x̄)²σ²/S_xx. Confidence band: ŷ ± t*_{n−2} s_e√(1/n + (x₀ − x̄)²/S_xx), narrowest at x̄. Prediction band: add σ² for the new observation's own scatter, ŷ ± t* s_e√(1 + 1/n + (x₀ − x̄)²/S_xx), also narrowest at x̄ but never narrower than t* s_e.
(d) g(θ̂) ≈ g(θ) + ∇gᵀ(θ̂ − θ), so Var g(θ̂) ≈ ∇gᵀ C ∇g. For I = aΓ(b + 1)c^{−(b+1)}: ∂I/∂a = I/a; ∂I/∂b = I[Γ'(b + 1)/Γ(b + 1) − ln c] = I[ψ(b + 1) − ln c]; ∂I/∂c = −(b + 1)I/c. These are the three components `impulse_gradient` returns.

## 3. A lift curve that bends

```python
rng = np.random.default_rng(15); al = np.arange(0, 17, 2.0)
cl = 0.10 + 0.100*al - 0.0012*np.maximum(0, al - 8)**2 + rng.normal(0, 0.02, size=9)
```

(a) All nine points: slope **0.0971 ± 0.0020 per degree**, r² = 0.9972, s_e = 0.0303.
(b) Residuals in order of α: −0.031, −0.015, +0.017, +0.005, +0.031, +0.038, −0.016, +0.015, −0.045. Negative at both ends, positive in the middle: the line is too low in the middle and too high at the ends, which is the signature of a curve that bends downward at high α (the onset of stall) being fitted by a straight line. r² = 0.997 looked excellent and said nothing about this.
(c) α ≤ 8° only (five points): slope **0.1043 ± 0.0018 per degree**, s_e = 0.0114. The slope rises by 7 % and s_e falls by a factor 2.6: the five-point fit is the one consistent with 0.02 noise.
(d) Report 0.104 ± 0.002 per degree (5.98 per radian) as the lift-curve slope *for 0 ≤ α ≤ 8°*, and state that the residuals above 8° show systematic curvature. The nine-point slope is a compromise between two regimes and belongs to neither.
(e) No. The largest |residual| in the nine-point fit is 1.49 s_e (the 16° point); no point reaches 2 s_e. A wrong model form spreads its misfit over all points and inflates s_e, which lowers every standardised residual. An outlier is one point that disagrees with a right model; a pattern is every point disagreeing a little with a wrong one. `outlier_flags` tests the first and cannot see the second; a residual plot sees both.

## 4. A noisier thrust curve

`t, F = synthetic_thrust_curve(seed=15, n=60, noise=17.5)`.

(a) Gauss–Newton from (1000, 0.5, 1.0): **8 iterations**; a = 1734 ± 104, b = 0.872 ± 0.040, c = 1.691 ± 0.054 s⁻¹; s_e = 17.36 N (noise 17.5). `curve_fit` with the same p0: a = 1733.84, b = 0.87244, c = 1.69101 with standard errors 104.31, 0.0401, 0.0543: identical to six figures in the parameters. The generating values (1550, 0.8, 1.6) sit 1.8, 1.8 and 1.7 standard errors away, all in the same direction because the three estimates are strongly correlated.
(b) I_model = 1734 × Γ(1.872)/1.691^{1.872} = **617.6 N·s**; delta method ± 8.72, Monte Carlo (20 000 draws, seed 15) ± 8.72. Trapezoid of the raw samples 613.6 N·s (misses the first 0.05 s and the tail beyond 4 s). True 619.5 N·s: the model value is 0.2σ from the truth; the trapezoid is 1 % low for structural reasons, not noise.
(c) `choose_degree` picks **degree 9**, LOO RMS 18.9 N, against the Gauss–Newton residual s_e of 17.4 N. Degrees 6–10 all have LOO between 18.9 and 19.8 N; the pick is inside the noise. Report the peak from the physical model: t_peak = b/c = 0.516 s, F_peak = 406.8 N (true 0.500 s, 405.8 N); it has three parameters with standard errors and a closed-form impulse, whereas the degree-9 polynomial has ten coefficients, no physical reading, and oscillates between samples.
(d) The true curve lies inside the 95 % confidence band at every one of 40 test points on 0.1–4 s. That is a little better than the nominal pointwise 95 % promises, which is expected: the band's errors at neighbouring t are strongly correlated, so a band tends to contain the whole truth or miss it over a stretch.

## 5. Outliers and leverage

`p, F = synthetic_pressure_thrust(seed=15)`: clean slope 298.5 N/MPa, s_e 25.9 N, 95 % CI (285.3, 311.8).

(a) F[8] += 150 (p = 4.91 MPa). Slope 305.8 (+7.2); s_e 44.2; residual of point 8 = 108.1 N = **2.45 s_e: flagged at k = 2, not at k = 3**. Of the 150 N injected, 108 N shows in the residual; the fit absorbed the rest (leverage h₈ = 0.127).
(b) F[0] += 150 (p = 2.0 MPa, the low end). Slope 282.7 (**−15.9**, more than twice the shift in (a)); s_e 41.1; residual 86.9 N = 2.11 s_e: flagged at k = 2 only just, not at k = 3. Leverage h₀ = 1/12 + (2 − 4)²/18.9 = 0.295: the end point pulls the line toward itself, which lowers its own residual (87 N of 150) and tilts the slope hard. A fault far from x̄ does the most damage and is the hardest to see.
(c) Case (a): before (283.1, 328.4), after removing point 8: 299.8, (285.9, 313.7). Case (b): before (261.6, 303.8), after removing point 0: 295.7, (280.1, 311.3). Both refits land back on the clean interval.
(d) Expected number beyond 2 s_e with 12 clean normal points: 12 × 0.0455 = 0.55. About every other clean dataset flags a genuine point. Deleting flagged points automatically would, roughly half the time, throw away a true measurement and shrink s_e below the real noise; a flag is a reason to check the logging record, not a deletion rule.

## 6. What the fit cannot tell you

(a) The 2372 ± 67 N interval says: *if the thrust stays a straight line in chamber pressure with the same slope and the same scatter, then the mean thrust at 8 MPa lies here with 95 % confidence*. It does not promise that the line continues past 6 MPa, because no data exist there to test it, and it is not a bound on any single firing (that is the prediction band, ± 96 N). Physical reasons it might fail: erosive burning raises the burn rate and the pressure–thrust relation at high mass flux; nozzle throat erosion changes A_t during a hotter firing; the thrust stand or transducer may be out of its calibrated range at 2.4 kN.
(b) Noise in x biases the slope **toward zero** (regression dilution), by λ = S_xx/(S_xx + nσ_x²). Here S_xx/n = 1.58 MPa² and σ_x² = 0.01 MPa², so λ = 0.9937: a 0.6 % attenuation, from 300.0 to a predicted 298.1 N/MPa. Simulation (seed 15, 2000 fits with `p + rng.normal(0, 0.1, size=12)`): mean fitted slope **298.5 N/MPa**, 0.5 % below the noiseless value and within 0.2 % of the prediction. For σ_x = 0.6 MPa the factor is 0.81 and the simulated mean slope 255: not negligible.
(c) The band is **pointwise**: at each fixed t, the interval contains the true mean thrust at that t in 95 % of repeated experiments. The probability that the whole true curve lies inside the whole band is smaller (how much smaller depends on how correlated the band's errors are along t). A statement about the entire curve needs a simultaneous band, which is wider.

## Check-yourself answers

1. Σeᵢ = 0 (the residuals have zero mean; the line passes through (x̄, ȳ)) and Σxᵢeᵢ = 0 (the residuals are uncorrelated with x; no linear trend remains).
2. 99.4 % of the variance of the thrust values is accounted for by the line. It does not say the line is the right form: a gently curved relation also gives r² near 1, and only the residual plot shows the curvature.
3. Var b = σ²/S_xx. Double the spread of x: S_xx grows 4-fold and se_b halves. Quadrupling n at the same x values does the same.
4. Confidence band: where the *mean* response lies. Prediction band: where a *single new observation* lies; it adds the measurement's own σ². Range safety needs the prediction band for a single future firing.
5. The residual s_e can only fall as parameters are added, because the fit is evaluated on the same points it was tuned to. LOO evaluates each point by a fit that never saw it, so it rises when added parameters fit noise.
6. One step solves the linear least-squares problem min |r − Jδ|², the fit of the residuals by the linearised model. At the solution the model is locally y ≈ F(p̂) + J(p − p̂), a linear model with design matrix J, whose parameter covariance is s_e²(JᵀJ)⁻¹ by (5).
7. Var g ≈ (dg/dc)² σ_c²; for g = 1/c, dg/dc = −1/c², so σ_g = σ_c/c². For c = 1.646 ± 0.030 s⁻¹: τ = 0.6075 ± 0.011 s.
8. The smallest F values: a fixed noise σ_F is an error of σ_F/F in ln F, largest where F is smallest. It goes wrong when the noise is comparable with the signal, as in the thrust tail below 60 N, where the log fit's τ fell from 0.89 s to 0.48 s and r² from 0.99 to 0.58.
