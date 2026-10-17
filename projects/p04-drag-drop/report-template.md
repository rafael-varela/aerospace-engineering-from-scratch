# Projectile with drag: Technical report

**Author:** · **Date:** · **Lesson:** P04.01 · **Equipment tier used:** · **Mode:** physical experiment

## Summary (5 sentences)
Question, method (two objects, frame rate, number of drops and tosses), headline result (C_D ± u from the drop fit; the toss range against prediction), whether the drag model held within uncertainty and whether the vacuum model was rejected, main limitation.

## 1. Objective and requirements
Test the quadratic-drag fall-time law for a ping-pong ball over 0.6–2.0 m against a dense control object; measure C_D with a stated uncertainty; compare a gentle toss with the simulated range. Success criteria set before the work: E_n < 2 against the drag model at every height, or an identified cause; the vacuum model rejected at the heights the prediction said it would be; u_CD/C_D ≤ 20 %.

## 2. Theory and assumptions
Newton's second law with weight and quadratic drag; terminal speed; the closed forms v(t), s(t), t(h) and the first-order excess gh/(6v_t²); the coupled two-dimensional equations solved numerically; the assumptions (constant C_D, no spin, still air, constant ρ, release from rest) and where each is expected to fail.

## 3. Prediction (recorded before measurement)
The prediction table: height, vacuum time, drag time for each object with its uncertainty budget (h, g, C_D, ρ, m, d); the expected E_n against the vacuum model at your frame rate; the chosen frame rate and number of drops with the reason; the toss range band and the expected effect of spin.

## 4. Apparatus, calibration, and safety
Objects (mass, diameter, how measured), air density and its source, phone and measured frame rate, frame rule, release method, height marks and their check, the completed safety box, deviations from the procedure.

## 5. Results
Raw frame tables (appendix); per-height statistics for both objects; E_n tables against both models; the fit of v_t with residuals and C_D ± u; the toss table (launch speed and angle, measured and predicted range, E_n, spin note); the plot with k = 2 error bars. If any SYNTHETIC data were used (no kit), label every table and plot that contains them SYNTHETIC and state that no measurement was made.

## 6. Discussion
Agreement and disagreement; the control object's role; ranked causes for any E_n > 2 and the follow-up test; the toss against prediction and what spin did; how the uncertainty budget held up (which term dominated, whether the reported u_CD was conservative).

## 7. Conclusions and redesign
What the evidence supports about the drag model and about C_D; what you would change (frame rate, release, heights). A mechanical launcher is not permitted in this project: say what permitted change would reduce spin instead.

## 8. Reproducibility
Code commit, raw data files, video file names, calibration file, BOM tier.

## Appendix: frame tables, fit residuals and uncertainty budget
