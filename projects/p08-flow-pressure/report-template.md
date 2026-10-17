# Flow and pressure measurement: Technical report

**Author:** · **Date:** · **Lesson:** P08.19 · **Equipment tier used:** · **Mode:** physical experiment (or: computational, SYNTHETIC dataset seed ___, if no measurement was made)

## Summary (5 sentences)
Question, method (bottle and hole, timing and level reading; manometer, probe, sources, upright and tilted), headline results (C_d ± U with k = 2 and E_n against the typical value and against the prediction; airspeeds ± U for each source, or "below resolution"), whether the models held, main limitation.

## 1. Objective and requirements
Test Torricelli's law with a discharge coefficient, and measure airspeed with a pitot probe and a water U-tube; success criteria: E_n < 2 for the draining time against the frozen prediction, χ² consistent with the reading uncertainties, airspeed uncertainties stated, and the resolution prediction for the fan confirmed or refuted.

## 2. Theory and assumptions
Bernoulli from the free surface to the vena contracta and mass conservation for the bottle (quasi-steady, (a/A)² ≪ 1); C_d as the product of the contraction and velocity coefficients; √h linear in t. Pitot total pressure, free-jet static pressure equal to the room's, Δp = (ρ_w − ρ_a) g Δh, V = √(2Δp/ρ_a); the tilted tube; incompressibility at these speeds; Reynolds numbers of the jet, the probe and the source. State where each assumption is used.

## 3. Prediction (recorded before measuring)
Table: D, d with uncertainties; C_d assumed; the 200 → 30 mm time ± u with the budget by term; jet Reynolds numbers; air density; the manometer's resolution limit upright and tilted; the predicted Δh for each source; the decision on tilting.

## 4. Apparatus, calibration, and safety
Bottle (type, wall), hole (how made, how gauged), ruler zero, camera and timing check; tube bore, water fill, probe, board and tilt measurement; air sources (rating plates); the calibration results (u_t, u_h, reading uncertainty, tilt); the completed safety box from `safety.md`; deviations from the procedure.

## 5. Results
Drain: C_d ± U (k = 2) for each run, budget, χ²/dof, mean and scatter; the measured time against the prediction. Manometer: table of source, setting, tilt, L, Δh, Δp, V ± U, Re of probe and source; the anemometer readings (enhanced). The √h–t plot and the airspeed plot.

## 6. Discussion
E_n values; χ² against the dof; why repeat runs agree far better than the stated uncertainty (a shared systematic); the fan against the resolution prediction; ranked causes for any disagreement; the follow-up test and its result; the effect of viscosity (Reynolds numbers) on both experiments.

## 7. Conclusions and redesign
What the evidence supports about Torricelli's law and the orifice, and about the pitot–manometer airspeed; what you would change (hole gauging, a bigger bottle, a longer tilted tube, a different liquid).

## 8. Reproducibility
Code commit, raw data files, photo and video file names, calibration file, BOM tier. If any part used SYNTHETIC data, say which and give the seed.

## Appendix: raw tables and uncertainty budgets
