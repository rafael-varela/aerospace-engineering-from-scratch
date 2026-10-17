# Pendulum system identification: Technical report

**Author:** · **Date:** · **Lesson:** P06.19 · **Equipment tier used:** · **Mode:** physical experiment (or: computational, SYNTHETIC dataset seed ___, if no measurement was made)

## Summary (5 sentences)
Question, method (logging tier, sample rate, record length, number of runs, release amplitudes), headline results (f₀ ± U and ζ ± U with k = 2, with the uncertainty method named), agreement with the prediction (E_n), the damping mechanism and the amplitude effect, main limitation.

## 1. Objective and requirements
Identify f₀ and ζ with defensible uncertainties and diagnose the damping. Success criteria set before the work: E_n < 2 for f₀ against the prediction; ζ within the predicted range or a tested cause; the uncertainty of each parameter checked by at least two methods.

## 2. Theory and assumptions
The pendulum equation with viscous, quadratic and Coulomb damping; the exact amplitude correction; the damped-cosine model and what it leaves out (constant frequency and damping ratio); the log decrement; the envelope laws for each damping type. State where each assumption is used.

## 3. Prediction (recorded before measurement)
f₀ ± u from L and g; f at the release amplitude and at half of it; b₂ and ζ_eq from air drag; the predicted amplitude after 120 s. Copy the dated `data/prediction.md`.

## 4. Apparatus, calibration, and safety
Support, string, bob or phone (mass with pouch), L and its readings, logging method and its time-base check, angle or rate noise and gyro bias, the completed safety box from `safety.md`, deviations from the procedure.

## 5. Results
Per run: damped frequency (crossings, FFT), ζ (log decrement over n cycles, ln-peak line), NLS parameters with covariance uncertainties, residual RMS against the noise, bootstrap sd (plain and block), envelope exponent p ± u, b₁ and b₂ with their correlation. Across runs: mean and scatter of f₀ and ζ against the per-run uncertainties. Plots: data with fit, decay envelope, frequency against amplitude with the exact law. Mark every SYNTHETIC item.

## 6. Discussion
E_n for f₀ and the comparison of ζ with the predicted range; which uncertainty estimate you report and why (covariance, bootstrap, repeat-run scatter); the damping diagnosis and how confident it is; the amplitude effect detected or not, against the prediction; ranked causes of any disagreement and the follow-up test.

## 7. Conclusions and redesign
What the evidence supports; what you would change (record length, amplitude range, logging tier, bob shape, pivot).

## 8. Reproducibility
Code commit, raw data files, video or app export files, calibration file, BOM tier. If any part used SYNTHETIC data, say which and give the seed.

## Appendix: run table and uncertainty budget
