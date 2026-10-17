# P04.15  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` (and the probe scripts used to write this file) and given to the figures shown. Results computed from `make_synthetic_bench_data` are SYNTHETIC: they show the analysis, not a measurement.

## 1. A battery monitor and an indicator ($s = 15$)

**Assumptions.** Ideal resistors at nominal values; the ADC draws no current; fixed-V_f LED model.

(a) Total 25 kΩ: exact R_t = 16.67 kΩ, R_b = 8.33 kΩ. Rounding the top up and the bottom down: **R_t = 18 kΩ, R_b = 8.2 kΩ**. At 9.0 V the output is **2.817 V** (below 3.0 V ✓), the current **0.344 mA**, the power **3.09 mW**. Over a 500 mAh battery's life the monitor alone would take 1450 h to drain it: put it behind a switch or a GPIO-controlled transistor if the battery matters.

(b) One LSB at the ADC is 0.806 mV; on the battery it is 0.806 mV × 26.2/8.2 = **2.57 mV**.

(c) (3.3 − 2.2)/5 mA = 220 Ω exact; **220 Ω** is an E12 value, so the current is **5.00 mA**.

**Check.** Rounding the top up and the bottom down both lower the ratio, so the 3.0 V limit is met with any rounding error of the same sign.

**Common mistakes.** Rounding both to nearest (8.2 kΩ/18 kΩ happens to work here, but 15 kΩ/8.2 kΩ would give 3.18 V); forgetting that the divider's own current is a permanent load on the battery.

## 2. Four results the lab relies on

(a) The error e = Q(v) − v of a rounding quantizer lies in [−Δ/2, Δ/2], Δ = LSB. For a signal spread over many codes, e is uniform on that interval: E[e] = 0, E[e²] = (1/Δ)∫_{−Δ/2}^{Δ/2} e² de = Δ²/12, so **σ_q = Δ/√12**. ∎ (Not true for a constant input, where the error is a fixed bias; that is why it is a noise *model*.)

(b) Thevenin: the divider's open-circuit output is V_th = V_inR_b/(R_t + R_b) and its output resistance (sources shorted) is R_th = R_t ∥ R_b. The meter R_m forms a divider with R_th: reading = V_thR_m/(R_th + R_m), so **reading/V_th − 1 = −R_th/(R_th + R_m)**. ∎ For R + R, R_th = R/2: −(R/2)/(R/2 + R_m) = **−R/(R + 2R_m)**; at R = 1 MΩ, R_m = 10 MΩ: −1/21 = −4.76 %.

(c) V_out = V_refR_f/(R_th + R_f) and the ADC reports raw = FS × V_out/V_ref (both referred to the same V_ref), so raw/FS = R_f/(R_th + R_f): **V_ref has cancelled**. Solving: **R_th = R_f(FS/raw − 1)**. The fixed resistor's error passes straight into R_th: a 1 % error in R_f is a 1 % error in R_th.

(d) Δ(R) = R/(R_lo + R) − R/(R_hi + R). dΔ/dR = R_lo/(R_lo + R)² − R_hi/(R_hi + R)² = 0 ⇒ R_lo(R_hi + R)² = R_hi(R_lo + R)² ⇒ expanding, R_loR_hi² + R_loR² = R_hiR_lo² + R_hiR² ⇒ (R_hi − R_lo)R² = R_loR_hi(R_hi − R_lo) ⇒ **R = √(R_loR_hi)**. ∎ For a β = 3950 K, 10 kΩ part over 0–50 °C: √(33.6 kΩ × 3.59 kΩ) ≈ 11.0 kΩ (the code's 11 011 Ω from the synthetic fit).

## 3. A spacecraft bus-voltage monitor ($s = 15$: 33 V full range; paper only)

(a) k = 33/3.0 = 11: exact R_b = 9.09 kΩ, R_t = 90.9 kΩ; nearest E12 in log ratio: **R_t = 100 kΩ, R_b = 10 kΩ**, k = 11.0 exactly. Full scale (3.3 V at the ADC) is **36.3 V** on the bus; one LSB is **8.86 mV** of bus voltage.

(b) Worst case (top +1 %, bottom −1 %, and the reverse): **+1.84 % / −1.80 %**. RSS with rectangular 1 % tolerances: **0.74 %**. At 33 V the worst case puts 3.055 V on the ADC: still below 3.3 V, but the 3.0 V headroom is mostly used.

(c) R_th = 100 kΩ ∥ 10 kΩ = **9.09 kΩ**: **not** well below 10 kΩ. Either halve both resistors (and accept 0.3 mA from the bus), or add a capacitor at the ADC pin to supply the sampling charge, or buffer with an op-amp follower (`P04.10`).

(d) No: the bus is not the ADC's reference, so V_ref's own error and drift enter the reading directly. The telemetry calibration must include the reference voltage (measured, or a precision reference monitored by a second channel); the thermistor channel needs only R_f.

**Common mistakes.** Using the 28 V nominal as the full scale (no headroom for a charging bus); ignoring the sampling input's source-resistance requirement.

## 4. The thermistor analysis ($s = 15$: SYNTHETIC data, seed 15)

(a) **SYNTHETIC.** Beta fit: **R₂₅ = 10 066 Ω, β = 3937 K**, residuals (°C) +0.032, −0.076, +0.121, −0.160, +0.084. Steinhart–Hart: **a = 1.0593 × 10⁻³, b = 2.4645 × 10⁻⁴, c = 2.918 × 10⁻⁸**, residuals +0.008, −0.059, +0.149, −0.148, +0.050. Both are within the ±0.25 °C (one standard deviation) scatter the generator put on the reference temperatures; with two degrees of freedom left, Steinhart–Hart's smaller residuals do not prove it is the better model.

(b) **SYNTHETIC.** 25 °C point: raw 32 747.25, Type A 8.00 counts, quantization 4.62 counts, combined **u_raw = √(8.00² + 4.62²) = 9.24 counts**; R_f = 10 kΩ ± 57.7 Ω. **T = 25.15 °C, u(T) = 0.131 °C**; contributions: counts **0.011 °C**, fixed resistor **0.130 °C**. Dominant: the fixed resistor's 1 % tolerance.

(c) **SYNTHETIC.** Gain **49.92 µV/count** (ideal 50.35), offset **−1.76 mV**, residual standard deviation **1.82 mV** = 2.3 LSB of a 12-bit converter: the calibration is limited by the multimeter readings' scatter, not by the ADC's resolution.

(d) Every number above is from simulated data and is labelled so.

**Common mistakes.** Treating the 16-bit scale as the resolution (quantization of 1/√12 count instead of 16/√12); fitting Steinhart–Hart in °C instead of kelvin; reporting residuals smaller than the reference's uncertainty as evidence of accuracy.

## 5. Measuring your meter's input resistance ($s = 15$: V_r = 2.45 V)

(a) The bottom resistor in parallel with the meter is R_b′ = V_rR_t/(V_in − V_r) = 2.45/2.55 MΩ = 0.9608 MΩ. Then R_m = R_b′R_b/(R_b − R_b′) = **24.5 MΩ**.

(b) Meter spec ±(0.5 % + 2 digits) at 1 mV: u(V_r) = 8.23 mV, u(V_in) = 15.6 mV (rectangular). Propagated: **u(R_m) = 5.6 MΩ** from the meter readings alone; adding 1 % rectangular tolerances on both resistors (5.0 MΩ from them alone): **7.5 MΩ** in total, 31 %. The readings dominate slightly. The reason is structural: R_m is computed from the small difference R_b − R_b′ (4 % of R_b), which amplifies every input error by about 25.

(c) Make the meter's effect larger and the inputs better known: use two 10 MΩ resistors (the reading drops to 2.076 V, a 17 % effect instead of 2 %), and measure each resistor with the ohms range before building (say to 0.5 %). With those changes the propagated u(R_m) is **3.0 %**, inside the ±10 % target; with 4.7 MΩ resistors it is 5.8 %.

**Common mistakes.** Forgetting that V_in must also be measured (it carries the same meter's error); computing R_m from the 2.500 V nominal rather than the measured supply.

## 6. Trusting the measurement chain

(a) 3.05 V against 2.970 V predicted is +2.7 %: two 5 % resistors at the ends of their bands can give up to 3.09 V, two 1 % resistors at most 2.994 V. Suspect the resistor values (5 % parts, or a wrongly read colour band). Decide by measuring both resistors unpowered.
(b) −18 %: far more than a shunt (−3 %); suspect the LED's forward voltage (a red LED at 2.3 V instead of 2.0 V gives 8.2 mA; a different colour more). Decide by measuring the voltage across the 330 Ω resistor without the ammeter: I = V_R/330.
(c) A smile is nonlinearity near the rails; decide by recalibrating over 10–90 % of the range only and checking that the residuals become scatter.
(d) Self-heating: still air removes heat more slowly than stirred water. Decide by repeating in still air with a 47 kΩ fixed resistor (much lower current): the offset should shrink.
(e) `read_u16()` scales to 16 bits regardless of the hardware; look at the low four bits of many readings (always zero or repeating the top bits) or the board's documentation; the resolution is the hardware's.
(f) A 2 % reference drift would not move a ratiometric thermistor channel at all, so either the channel is not ratiometric (its divider is powered from a different supply than the ADC reference) or the fixed resistor or the thermistor itself has aged. Decide by checking the schematic, then correlating the temperature drift with the reference telemetry over the year.

## Check-yourself answers

1. Voltmeter in parallel across the part; ammeter in series, with the circuit broken; resistance unpowered because the meter supplies its own test current and any other current corrupts the reading (and can damage the meter).
2. (5 − 2)/0.010 = 300 Ω; choose 330 Ω (9.09 mA). Rounding up guarantees the current stays below the target.
3. 2.381 V (−4.76 %); R_th = 500 kΩ.
4. 0.806 mV; 0.233 mV rms. No: most ports scale a 10- or 12-bit conversion to the 16-bit range.
5. The divider output and the ADC full scale are both proportional to V_ref, so their ratio is independent of it. The fixed resistor's value (and the ADC's own offset, gain and nonlinearity) do not cancel.
6. R = R₀exp(β(1/T − 1/T₀)), linear in ln R₀ and β with x = 1/T − 1/T₀; 1/T = a + b ln R + c(ln R)³, linear in a, b, c.
7. The fixed resistor's 1 % tolerance (0.130 of 0.131 °C); measure R_f with the multimeter (or buy a 0.1 % resistor).
8. Battery or USB only, never mains (electric shock); ≤ 9 V DC and ≤ 500 mA (overheating, burns, fire from a short); every LED through a resistor (LED destruction and overheated supply); no lithium-cell charging circuits (thermal runaway, fire); capacitors ≤ 25 V and ≤ 1000 µF, discharged before touching (stored-energy discharge, burst capacitor).
