# Airfoil bench rig: Technical report

**Author:** · **Date:** · **Lesson:** P10.23 · **Equipment tier used:** · **Mode:** physical experiment (or: computational, SYNTHETIC dataset seed ___, if no measurement was made)

## Summary (5 sentences)
Question, rig (fan, straightener, model, balances, airspeed method), headline result (CL slope, CL_max and CD at the Reynolds number reached, with expanded uncertainties k = 2; E_n against the prediction), how the results compare with thin-airfoil theory, the panel code and Ladson's NACA data, and the main limitation.

## 1. Objective and requirements
Measure CL(α) and CD(α) of a NACA 0012 section at Re ≈ 3 × 10⁴; success criteria: a dated prediction with uncertainty, CL to ±10 % or better at α ≥ 4°, CD with a stated uncertainty, every disagreement with theory and published data explained.

## 2. Theory and assumptions
Thin-airfoil cl = 2πα; the panel code (inviscid); the finite-span correction and the effective aspect ratio you assumed for the end-plated model; the low-Re drag estimate; Thwaites and Michel (laminar separation before transition); the manometer relation q = (ρ_w − ρ)g Δl sin θ; the balance equations with tares and the arm ratio; the blockage estimate. Say where each assumption is used.

## 3. Prediction (recorded before measuring)
Table from `data/prediction.md`: model dimensions, room air, predicted V and q, Re, and for each α the predicted CL band, CD band and scale changes in grams, with the points flagged as resolution-limited; the expected direction of the disagreement and why.

## 4. Apparatus, calibration, and safety
Fan (model, rated power from its label, guard condition), straightener and settling length, model construction and measured dimensions, balance geometry, calibration results (lift-scale slope, arm ratio R ± u with up/down friction, manometer incline ± u, airspeed scatter), the completed safety box from `safety.md`, deviations from the procedure.

## 5. Results
Tares (sting lift and drag, zero-lift angle); the CL and CD table against α with uncertainties and the number of repeats; the lift slope, flow angle, CL_max and stall angle; plots of CL and CD against α with k = 2 bars, the prediction bands, the panel curve and Ladson's data.

## 6. Discussion
E_n values; the Reynolds-number gap (Re here against Ladson's), with the Thwaites/Michel argument and the published comparison of slopes, CL_max and CD; ranked alternative explanations; the follow-up test and its result; what the fan's unsteadiness and non-uniformity did.

## 7. Conclusions and redesign
What the evidence supports about low-Reynolds-number aerodynamics of this section; what you would change (scales, airspeed measurement, model, straightener, end plates).

## 8. Reproducibility
Code version, raw data files, photo file names, calibration files, BOM tier. If any part used SYNTHETIC data, say which and give the seed.

## Appendix: raw tables and uncertainty budget
