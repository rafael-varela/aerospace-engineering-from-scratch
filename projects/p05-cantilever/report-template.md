# Cantilever strip: Technical report

**Author:** · **Date:** · **Lesson:** P05.18 · **Equipment tier used:** · **Mode:** physical experiment (or: computational, SYNTHETIC dataset seed ___, if no measurement was made)

## Summary (5 sentences)
Question, method (strip, clamp, reading method, loads and marks), headline result (compliance and E ± U with k = 2, E_n against the prediction; θ₀ detected or not, with its uncertainty), whether the model held, main limitation.

## 1. Objective and requirements
Test δ = Pa³/(3EI) and the deflected shape w(x) of a clamped strip with a typical modulus; success criterion: E_n < 2 for the compliance against the prediction, χ² consistent with the reading uncertainty, and the clamp-rotation detectability predicted correctly.

## 2. Theory and assumptions
E I w″ = M integrated with the clamp conditions; the straight tail beyond the load; self-weight and why differencing removes it; I = bh³/12; the clamp stress and the safe load; the uncertainty budget; the clamp-rotation term θ₀x; the large-deflection limit. State where each assumption (rigid clamp, small deflection, uniform section, vertical load) is used.

## 3. Prediction (recorded before loading)
Table: measured b, h, a, L with uncertainties; I, EI; safe load and chosen top load; δ ± u at each planned load with the budget by term; predicted compliance; predicted profile at the top-but-one load; predicted u(θ₀) and the detectable / not-detectable verdict.

## 4. Apparatus, calibration, and safety
Strip (material, source, dimensions), clamp and clamping arrangement, reading method and its calibrated uncertainty, mass calibration, the completed safety box from `safety.md`, deviations from the procedure.

## 5. Results
Load series: differences with uncertainties; the fitted slope, intercept, χ²/dof; E ± U (k = 2) with its budget. Profile: EI, θ₀ ± U, χ²/dof. Load–deflection plot with the prediction band; profile plot with both curves.

## 6. Discussion
E_n values; the intercept; χ² against the dof; the clamp result against the detectability prediction; ranked causes for any disagreement; the follow-up test and its result; whether permanent set or creep was seen.

## 7. Conclusions and redesign
What the evidence supports about the beam model, the modulus and the clamp; what you would change (reading method, thickness measurement, clamp, load range).

## 8. Reproducibility
Code commit, raw data files, photo file names, calibration file, BOM tier. If any part used SYNTHETIC data, say which and give the seed.

## Appendix: raw tables and uncertainty budget
