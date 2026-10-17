# P04.10  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Any simulated signal you generate is SYNTHETIC: label it so in every table and plot.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a load on the 400 Hz bus)

A galley load on a 115 V rms, 400 Hz bus is R = 8 + s/5 Ω in series with L = 1.5 mH.

(a) Compute X_L, |Z|, the phase, the RMS current, P, Q and the power factor. Check with `ac_power`.
(b) What capacitance in parallel with the load makes the power factor 1? What is the line current then, and what is the load's real power?
(c) Why do aircraft electrical specifications limit power factor? One sentence.

Answer form: seven numbers with units, a capacitance and a current, one sentence.

## 2. (Derivation: the RC filter, the loaded ladder and the op-amp)

(a) From the voltage-divider rule with impedances, derive (2). Show that |H| = 1/√2 and the phase is −45° at f_c = 1/(2πRC), that the gain falls at 20 dB per decade far above f_c, and that the low-frequency group delay is RC.
(b) Write the nodal equations of two RC sections connected directly (source → R → node 1 with C to ground → R → node 2 with C to ground) and derive (4). At what frequency does the unbuffered pair reach −3 dB, in terms of 1/(RC)?
(c) Derive the non-inverting gain (5) from the two ideal-op-amp rules. Then keep A finite and show G = A/(1 + Aβ); expand for large A to show the error is about G²/A.

Answer form: three derivations, one frequency in closed form.

## 3. (Unfamiliar scenario: an LC filter on a 28 V DC bus)

A satellite's 28 V DC bus feeds a noisy payload through an input filter: a series inductor L = 100 μH (with its winding resistance R = 0.5 Ω) and a shunt capacitor C = (10 + s) μF across the payload's input. The filter's transfer function from bus noise to payload voltage is H = Z_C/(R + Z_L + Z_C).

(a) Find the resonant frequency, Q, and |H| at f₀. Explain why the filter *amplifies* noise near f₀.
(b) Sweep |H| from 100 Hz to 100 kHz with your impedance functions; find the peak and its frequency, and |H| at 100 kHz in dB. What is the roll-off rate above f₀?
(c) Repeat for R = 0.05 Ω and R = 2 Ω. What value of R gives Q = 1, and why do power-supply designers add damping resistance even though it costs power?

Answer form: three numbers, a sweep summary, a table of two peaks, a resistance and a sentence.

## 4. (Numerical, uses your code: your own anti-aliasing filter)

A strain-gauge channel is sampled at f_s = 200 Hz by a 10-bit ADC with v_ref = 3.3 V. The signal band is 0–10 Hz; a vibration couples into the wiring at f_vib = (610 + 20s) Hz with amplitude 20 mV.

(a) Where does the vibration appear after sampling, and by how much must it be attenuated to stay below 1 LSB?
(b) With `required_stages_and_corner`, find the number of buffered RC sections and the largest corner for a 0.5 dB droop limit at 10 Hz. Tabulate the corner and the droop for one, two and three sections.
(c) Choose R from the E12 series and C ∈ {0.1, 0.22, 0.47, 1} μF; report f_c, the residual in LSB, the droop and the delay.
(d) Check the design with `synthetic_sensor_signal` (2 Hz, 0.5 V, plus your vibration) and `filter_and_sample`; report the alias tone in the counts with and without the filter, labelled SYNTHETIC. Your corner may lie above the Nyquist frequency: what does that leave unprotected?

Answer form: a frequency and a factor, a 3-row table, four numbers, two numbers and a sentence.

## 5. (Numerical: what averaging an ADC can and cannot do)

A constant voltage v₀ = 1.000 V + s × 0.1 mV is read by the 12-bit, 3.3 V ADC of the lesson, and the converted values (counts + ½) LSB are averaged over N readings. Before conversion, Gaussian noise of standard deviation σ is added (seed `numpy.random.default_rng(s)`, 400 repeats per case).

(a) For σ = 0, what count does the ADC return, and what is the error of the average for any N? Why does averaging not help?
(b) For σ = 0.3, 1 and 2 LSB and N = 1, 16, 256, tabulate the bias and the RMS error of the average (in LSB) over the 400 repeats.
(c) Explain the results: when does noise *improve* the averaged reading, how does the RMS error scale with N, and what is the price?

Answer form: a count and an error, a 3 × 3 table of (bias, RMS), three sentences.

## 6. (Judgment: choosing the signal chain)

For each, say what you would do and why in two or three sentences, with a number where one is natural.

(a) A colleague proposes to remove the 30 Hz alias in the lesson's channel with a digital notch filter after the ADC.
(b) A 16-bit ADC is offered as an "upgrade" for the lesson's channel, with the same 3.3 V reference and the same RC design.
(c) The two RC sections of the design are to be wired without the op-amp to save a part.
(d) A thermocouple-like sensor produces 0–5 mV and is wired directly to the 12-bit, 3.3 V ADC.
(e) The attitude-control loop that uses the filtered sensor has a crossover frequency of 8 Hz.

Answer form: five short paragraphs.
