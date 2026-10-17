# Servo feedback control: Technical report

**Author:** · **Date:** · **Lesson:** P09.21 · **Option (A: feedback servo / B: gear motor):** · **Equipment tier used:** · **Board and code path (MicroPython / Arduino-style):** · **Mode:** physical experiment / analysis of the SYNTHETIC dataset (say which; if synthetic, write SYNTHETIC in the title)

## Summary (5 sentences)
Question; plant and sensor; identified model with uncertainties; designed gain and its predicted rise time, overshoot and settling time (with the date of the prediction); measured values and the main discrepancy with its cause.

## 1. Objective and requirements
The design target (overshoot or rise time), the success criteria set before measuring (for example: measured overshoot within 3 percentage points of the prediction, rise and settling within 15 %, the duty or command limit never exceeded), and the safety envelope.

## 2. Theory and assumptions
The plant model (first order with dead time, or first order plus integrator with a dead zone), the identification method (log-linear, tail line, Levenberg–Marquardt), the closed-loop transfer function, the second-order or first-order metric formulas, the sampled-data effects (one-sample delay, hold, quantisation, saturation) and which are in which model.

## 3. Prediction (recorded before closing the loop)
Table: run, step, gain(s), predicted rise time, overshoot, settling time, steady-state error, peak duty or command; the date and time; the identified parameters used. Then the second, improved prediction for the redesign, also dated.

## 4. Apparatus, calibration, and safety
Board, actuator, driver, sensor and its calibration (degrees per count, noise, resolution), the pack voltage, the logged results of the four pre-run tests (stop, limit, sign, travel), the completed safety box, deviations from the procedure.

## 5. Results
The open-loop step data with the fitted model and residuals; the identified parameters with standard errors and their spread over repeats; the closed-loop runs with the prediction overlaid; the comparison table (predicted, measured mean ± standard deviation over repeats, difference). Label every synthetic figure SYNTHETIC.

## 6. Discussion
Each discrepancy larger than the repeat spread; ranked causes; the improved model (delay as a lag, sampled simulation) and how much of the discrepancy it removes; what remains and why.

## 7. Conclusions and redesign
What the evidence supports; the redesign and its measured result; what you would change in the hardware, the firmware (loop rate, limit, compensation) or the identification.

## 8. Reproducibility
Code commit, firmware listing as run, raw data files, calibration file, prediction files with dates, BOM tier.

## Appendices
Derivations, full tables, the firmware, the safety-test log.
