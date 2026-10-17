# String pendulum: Technical report

**Author:** · **Date:** · **Lesson:** P04.06 · **Equipment tier used:** · **Mode:** physical experiment (or: computational, SYNTHETIC dataset seed ___, if no measurement was made)

## Summary (5 sentences)
Question, method (video frame rate, swings per timing, lengths and amplitudes), headline result (g ± U with k = 2, E_n against the reference; the amplitude effect detected or not at each amplitude), whether the model held, main limitation.

## 1. Objective and requirements
Test T₀ = 2π√(L/g) over your range of lengths and the exact amplitude dependence at L ≈ 1 m; success criterion: E_n < 2 at every length after the amplitude correction, and predicted detections confirmed.

## 2. Theory and assumptions
The pendulum equation and its small-angle limit; the exact period T = 4√(L/g)K(k²) and its series; the frame-quantisation model; g from L and T and its uncertainty. State where each assumption (point bob, rigid support, inextensible string, small amplitude) is used.

## 3. Prediction (recorded before measurement)
Table: target length, T₀ ± u, frame term per period, and at L ≈ 1 m the predicted excess at 10°, 20°, 30° with a detectable / not-detectable verdict.

## 4. Apparatus, calibration, and safety
Support, string, bob (mass and size), phone and calibrated frame rate, length reference and its calibration, the completed safety box from `safety.md`, deviations from the procedure.

## 5. Results
Per-length T ± u_T and g ± u_g (corrected); weighted mean g ± U (k = 2); fitted exponent p and g from the fit; T against L plot; T/T₀ against amplitude plot with the exact curve.

## 6. Discussion
E_n values and their pattern across lengths; the amplitude comparison against the prediction; ranked causes for any disagreement; the follow-up test and its result; how the uncertainty budget held up (was the Type A scatter consistent with the frame term?).

## 7. Conclusions and redesign
What the evidence supports about the model and your measurement; what you would change (frame rate, swings per timing, length reference, support).

## 8. Reproducibility
Code commit, raw data files, video file names, calibration file, BOM tier. If any part used SYNTHETIC data, say which and give the seed.

## Appendix: raw timing tables and uncertainty budget
