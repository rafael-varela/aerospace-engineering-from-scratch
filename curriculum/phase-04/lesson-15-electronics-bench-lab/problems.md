# P04.15  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Problems that describe a circuit you may build stay inside the limits of `projects/p04-electronics-bench/safety.md`; Problem 3 is a paper exercise only.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a battery monitor and an indicator)

A 9 V battery is to be read by a 3.3 V board's 12-bit ADC, so the divider output must not exceed 3.0 V at 9.0 V.

(a) Design the divider with a total resistance of 10(1 + s/10) kΩ. Choose E12 values that keep the output at or below 3.0 V (round the top resistor up and the bottom one down). Report the output at 9.0 V, the current drawn from the battery and the power.
(b) What battery voltage change does one ADC LSB correspond to?
(c) A green LED (V_f = 2.2 V) on the board's 3.3 V rail is to run at (2 + s/5) mA. Choose the resistor and give the actual current.

Answer form: two resistors, three numbers with units; one voltage; one resistor and one current.

## 2. (Derivation: four results the lab relies on)

(a) Show that the rms quantization error of an ideal rounding quantizer, for a signal spread over many codes, is LSB/√12.
(b) Derive the voltmeter loading error (2) of the lesson from Thevenin's theorem, and show that it reduces to −R/(R + 2R_m) for an equal divider of R + R.
(c) For the thermistor-on-top divider read ratiometrically, show that V_ref cancels and derive (6). Which component's error passes straight into R_th?
(d) Show that the fixed resistor maximising the change of the divider ratio between R_lo = R(T_hi) and R_hi = R(T_lo) is √(R_lo R_hi).

Answer form: four derivations.

## 3. (Unfamiliar scenario, paper only: a spacecraft bus-voltage monitor)

A satellite's 28 V nominal power bus is monitored by a 12-bit, 3.3 V ADC. The monitor must read up to (30 + s/5) V with the ADC input at or below 3.0 V. **Do not build this: it exceeds the lab's 9 V limit.**

(a) Design the divider with about 100 kΩ total, using E12 values (nearest in log ratio). What bus voltage gives full scale, and what bus voltage is one LSB?
(b) With 1 % resistors, what is the worst-case and the RSS (rectangular, k = 1) relative error of the divider ratio?
(c) The ADC's sampling input needs a source resistance well below 10 kΩ. Is the divider's Thevenin resistance acceptable?
(d) Is this measurement ratiometric? What, therefore, must the telemetry calibration include that the thermistor channel does not need?

Answer form: two resistors, two voltages; two percentages; one resistance and a yes/no; two sentences.

## 4. (Numerical, uses your code: the thermistor analysis)

Generate the SYNTHETIC dataset with `make_synthetic_bench_data(seed=s)` (or use your own measured sweep, and say which).

(a) Fit the beta and Steinhart–Hart models; report R₂₅, β, a, b, c and the residuals at each point.
(b) Build the uncertainty budget for the 25 °C point: Type A from the 16 samples, quantization of a 12-bit converter on the 16-bit scale, R_f = 10 kΩ ± 1 % rectangular. Report T, u(T), the contribution of each input and the dominant one.
(c) Calibrate the ADC against the multimeter column; report gain, offset and residual standard deviation, and compare the residual with one true LSB.
(d) Label every number and the figure as synthetic if they came from the generator.

Answer form: a fit table, a three-line budget, three numbers, and the label.

## 5. (Numerical with uncertainty: measuring your meter's input resistance)

With two 1.00 MΩ resistors (assume exact for (a)) as a divider from a supply you measure as 5.00 V, your meter reads V_r = 2.30 + s/100 V at the midpoint.

(a) Derive R_m from V_r, V_in and the resistors, and compute it.
(b) Propagate a meter specification of ±(0.5 % + 2 digits) at 1 mV resolution (rectangular) on both V_r and V_in into u(R_m). Then add 1 % (rectangular) tolerances on both resistors. Which term dominates?
(c) How would you redesign the experiment to measure R_m to ±10 %? Consider larger resistors and measuring the resistors first.

Answer form: one resistance, two uncertainties and a sentence, a short redesign.

## 6. (Judgment: trusting the measurement chain)

For each case, say which instrument, model or component you suspect and the one measurement that would decide it: (a) your 8.2 k/12 k divider reads 3.05 V from 5.00 V; (b) your LED draws 7.5 mA instead of the predicted 9.1 mA, with the meter in series; (c) the ADC calibration residuals form a smile, largest at the two ends; (d) the thermistor reads 0.4 °C higher in still air than in stirred water at the same reference temperature; (e) a friend says their board has a 16-bit ADC because `read_u16()` returns 0–65 535; (f) a temperature telemetry channel on a satellite drifts 1.5 °C over a year while the bus voltage drifts 2 %.

Answer form: six labelled lines.
