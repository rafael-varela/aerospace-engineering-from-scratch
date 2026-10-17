# Electronics bench lab: Technical report

**Author:** · **Date:** · **Lesson:** P04.15 · **Equipment tier used:** · **Board and code path (MicroPython / Arduino-style):** · **Mode:** physical experiment / analysis of the SYNTHETIC dataset (say which; if synthetic, write SYNTHETIC in the title)

## Summary (5 sentences)
Question; circuits and instruments; headline results (divider within band or not, meter input resistance with uncertainty, a temperature as value ± U with k = 2); agreement with the predictions; the main limitation.

## 1. Objective and requirements
The predictions to be tested and the success criteria set before measuring (E_n < 2, calibration residual below 2 true LSB, fit residuals within the reference thermometer's uncertainty).

## 2. Theory and assumptions
Divider and loaded divider, LED fixed-V_f model, voltmeter and ammeter models, ADC quantization, ratiometric thermistor divider, beta and Steinhart–Hart models; each assumption and where it is used.

## 3. Prediction (recorded before measurement)
Table: quantity, predicted value, uncertainty, source of each input (datasheet, manual, measured resistor).

## 4. Apparatus, calibration, and safety
Board, meter (with its accuracy specification), resistor measurements, ADC calibration (gain, offset, residuals), the completed safety box, deviations from the procedure.

## 5. Results
Raw-data summary; the divider, loading and LED tables with E_n; the calibration plot with residuals; the thermistor fits with residuals; the uncertainty budget for one temperature point. Label every synthetic figure SYNTHETIC.

## 6. Discussion
Agreement or disagreement beyond uncertainty; ranked causes; the follow-up test and what it showed; which input dominates each budget.

## 7. Conclusions and redesign
What the evidence supports; what you would change in the circuit, the instruments or the procedure (for example a 0.1 % fixed resistor, a better reference thermometer, a different divider resistance).

## 8. Reproducibility
Code commit, raw data files, calibration file, board firmware or sketch, BOM tier.

## Appendices
Derivations, full data tables, the board program listing, uncertainty budgets.
