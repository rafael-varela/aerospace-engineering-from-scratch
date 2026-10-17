# P10.22  Worked solutions

Each solution states assumptions, units, the derivation, the result, an independent check, and the common mistakes. Seeded problems are worked for s = 15. Computed numbers come from `solution.py` and a probe script that calls it: a flat vortex lattice on the half-wing with a symmetry image, uniform spacing, lattices (chordwise × spanwise per half) 2 × 6, 4 × 12, 8 × 24 (written 8 × 48 for the full span). Measured values are from NACA TN-1669, TN-1671 and TN-1422 as transcribed in `finite_wing_cla.csv`, most digitized from figures.

## 1. Which kind of evidence

(a) **Code verification**: the symmetric problem must give a symmetric answer exactly; it catches indexing and sign errors, not modelling errors. (b) **Solution verification**: the change with lattice size estimates the discretisation error of this answer; it says nothing about whether the model is right. (c) **Validation**: a prediction against a measurement; it can show the model is adequate for this wing within the uncertainty, and only for wings like it. (d) **Code verification**: a known limit (2π for a two-dimensional flat plate; 6.263 is 0.3 % below); it cannot test three-dimensional effects. (e) **Solution verification**: a banded estimate of the numerical error; not the model error. (f) **Validation data quality**: the scatter between two reports on one planform (3 %) is part of the experimental uncertainty; any validation finer than that is not supported.

## 2. Solution verification and the section correction

(a) With f_i = f* + C h_iᵖ, h₂ = r h₁, h₃ = r² h₁: f₂ − f₁ = C h₁ᵖ (rᵖ − 1), f₃ − f₂ = C h₁ᵖ rᵖ (rᵖ − 1). Their ratio is rᵖ, so **p = ln((f₃ − f₂)/(f₂ − f₁))/ln r**. Then f* = f₁ − C h₁ᵖ = **f₁ + (f₁ − f₂)/(rᵖ − 1)**. The three values must converge monotonically ((f₃ − f₂)(f₂ − f₁) > 0) and be in the asymptotic range (p roughly constant under further refinement); oscillating convergence makes p meaningless.

(b) The GCI is meant to bound |f₁ − f*|/|f₁| with roughly 95 % confidence (the factor Fs = 1.25 is a safety factor for a three-grid study with a measured p). Treated as an expanded uncertainty with coverage factor 2, the standard uncertainty is about GCI/2, which is what is combined with the other standard uncertainties.

(c) Γ(y) = ½V c(y) a₀[α − α_i(y)] with a₀ = κ2π equals ½V (κc(y)) 2π [α − α_i(y)]: the same equation for a wing of chord κc and slope 2π. The induced angle α_i depends only on Γ(y) and the span, which are unchanged, so Γ(y) is the same and so is the lift L = ρV∫Γ dy. The scaled wing's area is κS, so its CL′ = L/(qκS) = CL/κ, and its aspect ratio is b²/(κS) = AR/κ. Hence CL_α(AR, a₀) = κ CL_α,2π(AR/κ). In a lattice the chord also sets the chordwise distribution of vorticity and the distance from bound vortices to control points, so the two problems are not identical; the transformation is an approximation that is exact in the lifting-line limit.

## 3. A prediction for the trainer

(a) κ = 5.9329/2π = 0.9443. CL_α = 4.6359, 4.5454, **4.4942** per radian on the three lattices; **p = 0.82**; Richardson **4.428**; GCI 1.85 %, **u_num = 0.042** per radian.

(b) ∂CL_α/∂κ = 3.617 per radian (central difference); u_κ = 0.0303, so **u_input = 0.110** per radian.

(c) u = √(0.042² + 0.110²) = 0.118; **CL_α = 4.49 ± 0.24 per radian (k = 2)**, 0.0784 ± 0.0041 per degree. The section slope dominates. Measure (or compute with coupling) a₀ of the NACA 2412 at the trainer's Re (about 4.5 × 10⁶ at cruise); the lattice error is second.

(d) Across the 15 rows the mean signed error is +0.2 % and the mean absolute error 4.3 %, with failures concentrated at AR 1.34 and the tapered TN-1671 wings. For an unswept AR 8 wing the TN-1422 rows (−3 to −4 %) and wing 7 of TN-1669 (−1 %) are the relevant evidence: both small and of the same sign as the lattice's convergence (the prediction falls with refinement). A model-error term of about ±4 % (one standard uncertainty) is defensible and would make the prediction 4.49 ± 0.43 per radian (k = 2); omitting it claims more than the validation shows.

**Common mistakes.** Quoting the coarse-lattice value; forgetting that the trainer's section is not the proxy section.

## 4. How the verdict depends on what you assumed

(a) Fraction within k = 2: **0.73** (1.5 % scatter), **0.80** (3 %), **0.93** (6 %). For s = 15, (3 + 1.5) % = 4.5 % lies between; run it for your own value.

(b) κ = 1: mean absolute error 3.9 %, mean signed **+3.7 %**, 73 % within; with the correction: 4.3 % absolute, **+0.2 %** signed, 80 % within. The correction removes the bias; the absolute error hardly changes because the scatter of the residuals (AR 1.34, taper) dominates it.

(c) Honest if the 6 % scatter is justified from the data (it is double the documented TN-1669/TN-1671 difference) and the report says the fraction depends on it, giving the 3 % result too. Misleading if the uncertainty was enlarged after seeing the residuals, so that "within uncertainty" was guaranteed: a validation metric must be computed with an uncertainty fixed before the comparison.

## 5. Does the GCI cover the error?

| planform | 2 × 12 | 4 × 24 | 8 × 48 | 16 × 96 | p (first three) | GCI | Richardson |
|---|---|---|---|---|---|---|---|
| AR 1.34, unswept | 0.03407 | 0.03307 | 0.03253 | 0.03224 | 0.88 | 2.5 % | 0.03188 |
| AR 9, taper 0.4 | 0.08399 | 0.08295 | 0.08232 | 0.08197 | 0.71 | 1.5 % | 0.08133 |
| AR 2.61, 60° | 0.03728 | 0.03623 | 0.03572 | 0.03548 | 1.06 | 1.6 % | 0.03525 |

(b) Actual error of the 8 × 48 value against 16 × 96: 0.88 %, 0.43 %, 0.68 %; all are inside the GCI (2.5, 1.5, 1.6 %). The 16 × 96 value is itself not converged, but the Richardson values lie beyond it in the same direction, so the GCI plausibly covers the true error too.

(c) The span load falls to zero at the tips like √(distance from the tip); a uniform lattice resolves that square-root behaviour with an error that falls roughly as the panel width, so p ≈ 1 (cosine spacing would do better).

## 6. Where the linear model fails

(a) At a low aspect ratio, vorticity rolls up along the side edges (and, for swept wings, the leading edge) into vortices that sit above the wing and induce extra suction: vortex lift, growing faster than linearly with α. A flat-wake linear lattice keeps the wake in the wing plane and cannot produce it. Because the measured CL(α) curves upward, the slope fitted to the lowest CL is the smallest; the report does not state the range used, so the 0.0280 value is at the low end and the VLM's linear 0.0325 (+16 %) compares a tangent at zero with a chord-like slope over an unknown range. The data README also notes the effect makes the 45° wing read higher than the unswept one at AR 1.34.

(b) Model: the chord-scaling correction and the half-span-uniform lattice treat the narrow tip of a taper-0.25 wing (with a large local chord-to-span gradient) poorly; tip stall or local nonlinearity at the highly loaded tip is also outside the linear model. Data: digitized from a figure with ±0.0003 per degree reading uncertainty, over C_L < 0.6, with only an approximate jet-boundary correction based on unswept-wing theory, and the documented 3 % scatter against TN-1669.

(c) It is an input uncertainty: u_κ (the lesson uses the 2 → 6 × 10⁶ difference, 3 %) should be widened by the expected difference between NACA 0012 and 64-/65-210 slopes. NACA TR-824, *Summary of Airfoil Data*, "Lift Characteristics of Smooth Airfoils", gives section lift data for the 6-series sections and is free (`@ntrs-naca-tr-824-summary-of-airfoil-data`).

(d) "For the twelve low-speed planforms tested in NACA TN-1669, TN-1671 and TN-1422, the vortex lattice with a section-slope correction predicts CL_α with a mean absolute error of 4.3 % and a largest error of 16 % (AR 1.34, where the measured lift is nonlinear); 80 % of the 15 comparisons agree within the k = 2 combined uncertainty. The validation supports its use for attached-flow lift slopes of wings with AR above about 2.5 and sweep up to 60°, not for low-aspect-ratio wings at moderate CL."

## Check-yourself answers

1. Code verification: is the code solving its equations correctly (exact properties, limits, observed order)? Solution verification: how large is the numerical error of a particular answer (grid study, GCI)? Validation: does the model represent reality for this purpose (comparison with measurements, with uncertainty)?
2. From three grids refined by r: p = ln((f₃ − f₂)/(f₂ − f₁))/ln r; the Richardson value f₁ + (f₁ − f₂)/(rᵖ − 1); the GCI Fs|f₂ − f₁|/|f₁|/(rᵖ − 1). About first order (0.7–1.1) here.
3. Multiply every chord by κ = a₀/(2π) about the quarter-chord line, run the 2π lattice, multiply the result by κ.
4. Numerical (GCI/2), input (section slope), and experimental (figure-reading uncertainty and run-to-run scatter, each as a uniform band).
5. Mean |E/D| = 4.3 %, max 16.2 %, 80 % within k = 2.
6. The AR 1.34 unswept wing: nonlinear vortex lift that a linear lattice with a flat wake cannot produce.
7. Because its value depends on the uncertainty assumed: with a 1.5 % scatter it is 73 %, with 6 % it is 93 %. It must be reported with that assumption and with the error metrics themselves.
