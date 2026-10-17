# Cooling curve: Technical report

**Author:** · **Date:** · **Lesson:** P07.15 · **Equipment tier used:** · **Mode:** physical experiment (or: computational, SYNTHETIC dataset seed ___, if no measurement was made)

## Summary (5 sentences)
Question, method (bar, warm-up water temperature as measured, sling, thermometer and reading interval), headline result (τ and h ± U with k = 2, E_n against the prediction; χ²/dof of one exponential; Bi), whether the lumped model held, main limitation.

## 1. Objective and requirements
Test τ = ρcV/(hA) with h from natural convection (Churchill–Chu) plus linearised radiation, and the lumped-capacitance model; success criterion: E_n < 2 for τ against the dated prediction, Bi < 0.1 confirmed, χ² explained.

## 2. Theory and assumptions
Lumped energy balance and solution; the Biot number; h_conv for the horizontal cylinder and the end faces; h_rad and its linearisation; the probe-lag model; the variable-h model; the series solutions and the product rule. State where each assumption (still air, uniform h, known ε, constant T∞, negligible probe) is used.

## 3. Prediction (recorded before measuring)
Table: measured D, L, mass with uncertainties; V, A, L_c, ρcV; room temperature; h_conv, h_rad, h; Bi; τ ± u with the budget by term; τ predicted for each third of the run.

## 4. Apparatus, calibration, and safety
Bar (material, source, dimensions, mass), sling and stand, thermometer or thermistor and how it was attached, room thermometer, calibration results (offset, reading uncertainty, probe response), the measured water temperature (≤ 60 °C), the completed safety box from `safety.md`, deviations from the procedure.

## 5. Results
Raw-data plot; semi-log plot with the fits and the prediction; the table of fits (log-linear, nonlinear with T∞ fixed and free, lagged, variable-h) with τ, uncertainties and χ²/dof; windowed τ against mean excess; h ± U (k = 2); E_n; Bi with the measured h; the distributed-model comparison.

## 6. Discussion
E_n values; what χ² says about a single exponential; the windowed trend against the correlation's θ^(1/4); the free-T∞ fit and why it differs; f_conv and its trade-off with ε; ranked causes for any disagreement; the follow-up test and its result.

## 7. Conclusions and redesign
What the evidence supports about the lumped model, the correlation in your room and the emittance; what you would change (draught shield, logging, a second bar, a better room thermometer).

## 8. Reproducibility
Code commit, raw data files, calibration file, BOM tier. If any part used SYNTHETIC data, say which and give the seed.

## Appendix: raw tables and uncertainty budget
