# P04.15  Electronics Bench Lab 🔬

> Every number a spacecraft sends home passed through a resistor divider, an ADC and a calibration someone wrote down. Build one, measure it with a meter you have learned to distrust correctly, and you will read telemetry differently for the rest of your career.

| | |
|---|---|
| **Learning outcomes** | 1. Design a voltage divider and an LED series resistor from Ohm's and Kirchhoff's laws, choose standard (E12) values, and check the power in each part. 2. Measure voltage (in parallel), current (in series) and resistance (unpowered) with a multimeter, and predict and measure the error a voltmeter's input resistance and an ammeter's shunt introduce. 3. Read a voltage with a microcontroller's ADC, convert counts to volts, state its resolution (LSB) and quantization noise (LSB/√12), and calibrate it against the multimeter. 4. Read a thermistor through a ratiometric divider, fit the beta and Steinhart–Hart models, and report a temperature with a propagated uncertainty whose dominant term you can name. 5. Run the whole experiment inside a written safety envelope, with a prediction recorded before each measurement. |
| **Prerequisites** | `P04.08` (Ohm's law, Kirchhoff's rules, power, measuring instruments), `P04.10` (the microcontroller as CPU + GPIO + ADC; diodes), `P04.12` (sensors; ADC sampling and resolution), `P00.13` (safety box), `P00.14` (uncertainty: Type A and B, propagation, calibration fit), `P02.16` (the uniform distribution and its variance). Skill check: you can compute the current through 330 Ω with 3 V across it and say which way a voltmeter is connected. |
| **Engineering question** | With a breadboard, a multimeter and any 3.3 V or 5 V microcontroller board, can you predict, then measure, a divider's output, an LED's current and a thermistor's temperature to within your stated uncertainties, and when they disagree, can you name which instrument, model or component is responsible? |
| **Assumptions and limits** | Ideal resistors at their measured values (tolerance 1 % or 5 % until measured); an LED modelled as a fixed forward voltage (wrong near its threshold: the real I–V curve is exponential); a voltmeter modelled as a resistance of about 10 MΩ and an ammeter as a small series resistance (both values from your meter's manual, or measured); an ideal ADC with a linear transfer and quantization (real ADCs have offset, gain error and nonlinearity, largest near the rails, which the calibration corrects to first order); a thermistor at the temperature of the water around it (self-heating and lag are the first failures); the beta model over a narrow range, Steinhart–Hart over 0–50 °C. Never mains, never above 9 V DC, never more than 500 mA. |
| **Analytical work** | Divider and loaded divider; Thevenin resistance; LED resistor and power; the voltmeter loading error −R_th/(R_th + R_m) and the ammeter burden; LSB and quantization noise from the uniform distribution; the ratiometric divider and why V_ref cancels; the beta and Steinhart–Hart models as linear least-squares problems; the best fixed resistor (geometric mean); a full uncertainty budget for one temperature point. |
| **Computational deliverable** | `starter/main.py`: `nearest_e_series`, `divider_output`, `divider_design`, `led_resistor`, `led_current`, `resistor_power`, `loading_error`, `thevenin_resistance`, `lsb`, `quantization_noise_rms`, `volts_to_counts`, `counts_to_volts`, `u16_to_volts`, `adc_calibration_fit`, `thermistor_beta_r`, `thermistor_beta_t`, `fit_beta`, `steinhart_hart_t`, `steinhart_hart_r`, `fit_steinhart_hart`, `thermistor_r_from_ratio`, `temperature_from_u16`, `best_fixed_resistor`, `dmm_standard_uncertainty`, `propagate`, `analyse_thermistor`; `make_synthetic_bench_data`, `write_synthetic_csv`, `read_bench_csv` are given. Tests prove: divider limits, complementarity and scale invariance; design round trips and E12 rounding; the LED current never exceeds the target after rounding up; the loading error equals −R_th/(R_th + R_m) to 10⁻¹² and −1/21 for 1 MΩ/1 MΩ; the ammeter burden −10/340; LSB scaling and the LSB/√12 noise to 1 % on a 200 001-point ramp; clipping; a known calibration line recovered exactly; beta and Steinhart–Hart round trips to 10⁻¹² and their fits recovering a known truth; beta as Steinhart–Hart with c = 0; the ratiometric inverse; the geometric-mean resistor maximising the swing; temperature from counts to 10⁻⁹; propagation against the product and power rules; the synthetic data reproducible, physical, analysable to within 1 °C, labelled synthetic and readable back; the figure. |
| **Practical activity** | **physical experiment** (optional: a learner without the kit does the whole analysis with the supplied SYNTHETIC dataset): build a divider, an LED circuit, a potentiometer and a thermistor divider on a breadboard; measure them with a multimeter and the board's ADC; calibrate the ADC; sweep the thermistor through 5–45 °C in water. See `projects/p04-electronics-bench/`. |
| **Equipment** | Level 1–2: solderless breadboard, jumper wires, assorted 5 % or 1 % resistors, two 1 MΩ resistors, red LEDs, a 10 kΩ potentiometer, a 10 kΩ NTC thermistor, a digital multimeter, any 3.3 V or 5 V microcontroller board with an ADC programmable over USB, a USB cable; a kitchen thermometer, mugs of cold and warm water, a small zip bag. Tiers and substitutes in `projects/p04-electronics-bench/bom.yaml`. |
| **Safety** | See `projects/p04-electronics-bench/safety.md` (authoritative). Hard limits: battery or USB power only, never mains; supply ≤ 9 V DC; current ≤ 500 mA; every LED through a resistor; no lithium-cell charging circuits; no soldering required; no capacitors above 25 V or 1000 µF, discharged before touching; electronics kept dry. |
| **Evidence** | `data/prediction.md` (dated before measuring), the raw CSV files, the ADC calibration table and residuals, the thermistor fit with residuals, the uncertainty budget for one temperature, the loading-error measurement, `bench_lab.png` from your data (or from the synthetic data, labelled as such), and a report following `projects/p04-electronics-bench/report-template.md`. |
| **Assessment** | `assessments/phase-04/grading-rubric.md`; all four competencies. Threshold: a prediction with uncertainty recorded before each measurement; the loading error predicted and measured; a temperature reported as value ± U (k = 2) with a budget naming the dominant term; every disagreement beyond 2 standard uncertainties explained by a named, tested cause; the completed safety box. A report built on the synthetic dataset must say so in its title and summary. |
| **References** | `@openstax-university-physics-vol2` Ch. 10 "Direct-Current Circuits": 10.4 Electrical Measuring Instruments (required-free): how voltmeters and ammeters work and why they load a circuit; 10.6 Household Wiring and Electrical Safety (required-free): why this lab stays at battery and USB voltages. `@sparkfun-learn-tutorials` "How to Use a Breadboard": "Anatomy of a Breadboard" (Terminal Strips; Power Rails; Rows and Columns), "Providing Power to a Breadboard", "Building Your First Breadboard Circuit"; "How to Use a Multimeter": "Measuring Voltage", "Measuring Resistance", "Measuring Current", "Continuity", "Changing the Fuse"; "Analog to Digital Conversion": "Relating ADC Value to Voltage" (required-free; ignore the purchase sections). `@micropython-docs` "Quick reference for the RP2": "Pins and GPIO", "ADC (analog to digital conversion)"; "class ADC - analog to digital conversion": Constructors (`machine.ADC(id, *, sample_ns, atten)`), Methods (`ADC.read_u16()`) (required-free; the lesson does not require that board). |
| **Extensions** | ★ Measure the LED's I–V curve with a potentiometer and fit the Shockley diode equation; compare with the fixed-V_f model. ★ Measure the thermistor's time constant by moving it from cold to warm water and fitting an exponential (an RC-like first-order system, `P04.10`). ★ Estimate the self-heating by reading the thermistor with 10 kΩ and with 1 kΩ fixed resistors in still air. |
| **Time** | unmeasured |

## 1. The engineering problem

A satellite's housekeeping telemetry is a few hundred numbers sent every few seconds: bus voltages, currents, temperatures, pressures. Almost every one of them is a sensor, a resistor network that brings its signal into an ADC's range, a conversion to counts, and a calibration that turns counts into engineering units. When a temperature reads 2 °C high, the operators must decide whether the box is hot, the thermistor is wrong, the divider resistor drifted, or the calibration was never right. That decision needs exactly the skills of this lab: knowing what a measurement chain does to the number, and having measured it.

The lab is the smallest real version of that chain. You will design and build a voltage divider and an LED circuit and check them with a multimeter, then discover that the multimeter itself changes what it measures (it draws current), predict by how much, and confirm it. You will read a potentiometer with a microcontroller's ADC and calibrate the ADC against the meter. And you will read a thermistor, a resistor whose value falls steeply with temperature, through a divider, fit its two standard models, and report a temperature with an uncertainty you can defend term by term.

Success criteria, written first: a dated prediction for every measured quantity; the divider within its predicted tolerance band; the loading error at 1 MΩ predicted and measured within 1 percentage point; the ADC calibrated with residuals below 2 LSB of its true resolution; the thermistor fit with residuals within the reference thermometer's uncertainty; and a temperature as value ± U (k = 2) with its dominant term named. A learner without the kit meets the same criteria on the supplied SYNTHETIC dataset, says so on every page, and does not claim a measurement.

## 2. Understand and derive

### 2.1 The breadboard and the meter

A solderless breadboard is rows of spring clips. In the central area each short row of five holes is one electrical node; the long rails along the edges are the power and ground buses (SparkFun, "Anatomy of a Breadboard"). A circuit is a set of nodes; a wrong row is a wrong node, and **continuity mode** (the meter beeps below a few tens of ohms) is how you check a node without power.

A multimeter measures **voltage across** two nodes (leads in parallel with the part), **current through** a branch (the circuit broken and the meter inserted in series, red lead moved to the mA jack), and **resistance** of a part that is not powered (the meter supplies its own small current, which a powered circuit would corrupt). The single most common way to blow a meter's fuse is to leave the leads in the current jack and touch them across a supply: the ammeter is nearly a short circuit (SparkFun, "Measuring Current", "Changing the Fuse").

### 2.2 The divider and the LED

Two resistors in series across V_in carry the same current V_in/(R_t + R_b), so the voltage across the bottom one is

**V_out = V_in R_b/(R_t + R_b)**.   (1)

It depends only on the ratio: 8.2 kΩ/12 kΩ and 8.2 Ω/12 Ω give the same 2.970 V from 5 V. The design rule is the inverse: for V_out from V_in with a total R, R_b = R V_out/V_in, R_t = R − R_b. A 5 V → 3.0 V divider of 20 kΩ needs 8 kΩ over 12 kΩ; the nearest E12 values (the 12 standard values per decade 1.0, 1.2, 1.5, … 8.2) give 8.2 kΩ/12 kΩ and **2.970 V**. Choosing R is a trade: smaller R draws more current (5 V/20 kΩ = 0.25 mA here); larger R makes the output easier to disturb (2.3).

An LED is a diode: almost no current below its forward voltage V_f (about 2.0 V for red, 2.2–3.2 V for other colours), and a steeply rising current above it. A voltage source connected directly to an LED therefore drives a current limited only by the source and the wiring, and destroys the LED; **every LED gets a series resistor**. With the fixed-V_f model, R = (V_s − V_f)/I. For 10 mA from 5 V: 300 Ω exact; round **up** to the E12 value 330 Ω, which gives **9.09 mA** and dissipates I²R = **27.3 mW** in the resistor (a ¼ W part is ample). From 3.3 V: 130 Ω exact, 150 Ω chosen, **8.67 mA**, 11.3 mW. Rounding up never overdrives the LED; rounding down might. Many microcontroller pins are rated for only a few milliamps; if you drive the LED from a pin rather than from a supply rail, design for the pin's rated current from your board's documentation, not for 10 mA.

### 2.3 What the meter does to the circuit

A digital voltmeter behaves as a resistance R_m, typically 10 MΩ, between its leads (OpenStax 10.4). Across the bottom resistor of a divider it is in parallel with R_b, so (1) becomes V_in R_b′/(R_t + R_b′) with R_b′ = R_b ∥ R_m. The cleanest way to see the error is Thevenin's theorem: the divider looks, from its output, like a source of V_out behind a resistance R_th = R_t ∥ R_b. The meter and R_th form a second divider, so

**reading/true − 1 = −R_th/(R_th + R_m)**.   (2)

For an equal divider of R + R, R_th = R/2:

| each resistor | R_th | error with a 10 MΩ meter |
|---|---|---|
| 1 kΩ | 500 Ω | −0.005 % |
| 10 kΩ | 5 kΩ | −0.05 % |
| 100 kΩ | 50 kΩ | −0.50 % |
| 1 MΩ | 500 kΩ | **−4.76 %** |
| 10 MΩ | 5 MΩ | −33 % |

From 5 V, a 1 MΩ/1 MΩ divider should give 2.500 V and a 10 MΩ meter should read **2.381 V**. That is a measurement you can make, and it is a prediction you can write down before you make it. Turned around, the reading measures R_m: if it reads V_r, then R_b′ = V_rR_t/(V_in − V_r) and R_m = R_b′R_b/(R_b − R_b′).

An ammeter is a small resistance (its shunt plus the leads), so inserting it adds resistance in series and lowers the current it measures. With a 10 Ω shunt (an illustrative value; your meter's manual or a measurement gives yours), the 330 Ω LED circuit's current falls from 9.091 to **8.824 mA, −2.94 %**: the fractional error is −R_shunt/(R + R_shunt) in the fixed-V_f model. The fix in both cases is the same as in every measurement: know the instrument's effect, correct for it or make it negligible.

### 2.4 The ADC

An n-bit ADC with reference V_ref maps 0…V_ref onto 2ⁿ codes; one code (the **LSB**) is V_ref/2ⁿ: 3.22 mV at 10 bits, **0.806 mV at 12 bits** and 3.3 V. Converting counts back to volts is counts × LSB (some documentation divides by 2ⁿ − 1 instead; the difference is under one LSB). A rounding quantizer makes an error between −½ and +½ LSB; for a signal that wanders over many codes the error is uniformly distributed, and the standard deviation of a uniform distribution of width LSB is

**σ_q = LSB/√12**   (3)

(`P02.16`: variance = width²/12), 0.233 mV at 12 bits. The tests verify (3) on a 200 001-point ramp.

MicroPython gives every port the same call (MicroPython docs, "class ADC"): `adc = machine.ADC(pin)` then `adc.read_u16()`, which returns 0–65 535 whatever the hardware's resolution; most boards scale a 10- or 12-bit conversion up to 16 bits, so **16 bits of range is not 16 bits of resolution**. Volts = raw/65 535 × V_ref. On an Arduino-style board the equivalent is the analog-read function, which returns counts at the hardware's own resolution (often 10 bits, 0–1023, on classic 5 V boards); volts = counts × V_ref/2ⁿ. The analysis code accepts either through `u16_to_volts` or `counts_to_volts`. Whatever the board, check two things in its documentation before connecting anything: the maximum voltage an analog pin tolerates (never more than the board's supply, often 3.3 V), and what it uses as V_ref.

Real ADCs have an offset, a gain error and nonlinearity. The first two are removed by a **calibration** against the multimeter: read several potentiometer settings with both, fit volts = g × raw + o by least squares (`P00.14`), and look at the residuals. Nonlinearity shows as structure in the residuals, often near the rails; stay in the middle 10–90 % of the range or calibrate there.

### 2.5 The thermistor and the ratiometric divider

An NTC thermistor's resistance falls roughly exponentially with temperature. The **beta model**

**R(T) = R₀ exp(β(1/T − 1/T₀))**,   (4)

T in kelvin, R₀ at T₀ = 298.15 K, is ln R linear in 1/T; a 10 kΩ part with β = 3950 K goes from 33.6 kΩ at 0 °C to 3.59 kΩ at 50 °C. Over a wider range ln R is not quite linear in 1/T, and the **Steinhart–Hart equation**

**1/T = a + b ln R + c (ln R)³**   (5)

absorbs the curvature. Both are linear in their parameters (ln R₀ and β in (4) against 1/T; a, b, c in (5)), so both fit by linear least squares in one line; (4) is (5) with c = 0 (the tests check the identity).

Put the thermistor on top and a fixed resistor R_f to ground, power the divider from the ADC's own reference, and read the junction. Then raw/full-scale = R_f/(R_th + R_f), and

**R_th = R_f (full-scale/raw − 1)**,   (6)

with **no V_ref in it**: this is a **ratiometric** measurement, and it is why a thermistor read this way does not care whether the 3.3 V rail is really 3.28 V. What it does care about is R_f; its tolerance passes straight into R_th. The best R_f for a range T_lo…T_hi maximises the change in ratio and is the geometric mean √(R(T_lo)R(T_hi)) (Problem 2): 11.0 kΩ for 0–50 °C, so a 10 kΩ part is close.

Two physical limits. **Self-heating**: the divider current heats the thermistor; at 25 °C with R_f = 10 kΩ from 3.3 V the thermistor dissipates 0.27 mW, which, with a dissipation constant of order 1.5 mW/°C in still air (illustrative; your part's datasheet gives it), raises it by about 0.18 °C; in stirred water much less. **Lag**: a thermistor in a bag in water needs tens of seconds to settle; read until the value stops drifting.

### 2.6 The uncertainty budget of one temperature

Take the 25 °C point of the SYNTHETIC dataset: mean raw 32 351 counts from 16 samples with a standard deviation of 31.5 counts. Inputs and standard uncertainties:

- Type A: 31.5/√16 = 7.87 counts;
- quantization: the hardware is 12-bit, so one true LSB is 16 counts on the 16-bit scale, and (3) gives 16/√12 = 4.62 counts (a conservative choice, since averaging 16 noisy samples reduces it);
- combined raw: √(7.87² + 4.62²) = **9.1 counts**;
- R_f = 10 kΩ ± 1 % (rectangular): 57.7 Ω.

Propagating through (6) and the fitted (5) with the finite-difference Jacobian of `P00.14`: T = **24.63 °C**, u(T) = **0.131 °C**, of which the counts contribute 0.013 °C and the resistor tolerance 0.130 °C. One true LSB is worth 0.022 °C here. The budget names the culprit: the fixed resistor, not the ADC. Measuring R_f with the multimeter (±(0.5 % + 2 digits) typical for a hobby meter on its resistance range) would roughly halve it; buying a 0.1 % resistor would make the ADC noise and the reference thermometer the limit. (The synthetic dataset's R_f is really 9.95 kΩ, 0.5 % low: with the true value the point reads 24.74 °C. That 0.11 °C shift is exactly the kind of systematic effect the budget warned about.)

## 3. Build the computational model

Open `starter/main.py`. The synthetic-data generator, its writer and the CSV reader are given; everything else is yours.

### Step 1: design

`nearest_e_series(value, series, direction)`: build candidates from the E12 list over three decades around the value, return the nearest in log ratio, or the smallest not below ('up') or the largest not above ('down'). `divider_output(vin, r_top, r_bottom, r_load)` from (1) with the load in parallel; `divider_design`; `led_resistor` (raise for V_s ≤ V_f), `led_current` (zero below V_f, with an optional extra series resistance for the ammeter), `resistor_power`.

### Step 2: instruments

`thevenin_resistance` and `loading_error` from (2), computed through `divider_output` so that the test against the Thevenin formula is an independent check.

### Step 3: the ADC

`lsb`, `quantization_noise_rms` (3), `volts_to_counts` (round and clip), `counts_to_volts`, `u16_to_volts`, `adc_calibration_fit` (least squares with residual standard deviation on n − 2 degrees of freedom).

### Step 4: the thermistor

`thermistor_beta_r`, `thermistor_beta_t`, `fit_beta`; `steinhart_hart_t`, `steinhart_hart_r` (Cardano's formula for the depressed cubic c L³ + b L + (a − 1/T) = 0, which has one real root because b, c > 0), `fit_steinhart_hart`; `thermistor_r_from_ratio` (6), `temperature_from_u16`, `best_fixed_resistor`, `analyse_thermistor` (resistances, both fits, residuals in °C).

### Step 5: uncertainty

`dmm_standard_uncertainty(reading, pct, digits, resolution)` = (pct % of reading + digits × resolution)/√3, a Type B term from the accuracy specification treated as a rectangular half-width; `propagate(f, x, u)` with a central-difference gradient.

Running `main()` prints

```
DESIGN
  divider 5 V -> 3.0 V, 20 kOhm total: exact 8000 / 12000 ohm; E12 8200 / 12000 -> 2.970 V
  red LED (Vf 2.0 V, 10 mA) from 5.0 V: exact 300 ohm, E12 up 330 ohm -> 9.09 mA, 27.3 mW in the resistor
  red LED (Vf 2.0 V, 10 mA) from 3.3 V: exact 130 ohm, E12 up 150 ohm -> 8.67 mA, 11.3 mW in the resistor
METER LOADING (10 MOhm voltmeter across the bottom resistor of an equal divider)
       1000 ohm +      1000 ohm: R_th       500 ohm, error -0.0050 %
      10000 ohm +     10000 ohm: R_th      5000 ohm, error -0.0500 %
     100000 ohm +    100000 ohm: R_th     50000 ohm, error -0.4975 %
    1000000 ohm +   1000000 ohm: R_th    500000 ohm, error -4.7619 %
   10000000 ohm +  10000000 ohm: R_th   5000000 ohm, error -33.3333 %
AMMETER BURDEN: 330 ohm LED circuit at 5 V: 9.091 mA without, 8.824 mA with a 10 ohm shunt (-2.94 %)
ADC
  10 bit at 3.3 V: LSB  3.2227 mV, quantization noise  0.9303 mV rms
  12 bit at 3.3 V: LSB  0.8057 mV, quantization noise  0.2326 mV rms
  16 bit at 3.3 V: LSB  0.0504 mV, quantization noise  0.0145 mV rms
SYNTHETIC DATASET (simulation, not a measurement)
  ADC calibration vs multimeter: gain 50.0149 uV/count (ideal 50.3548), offset -8.26 mV, residual sd 1.20 mV
  beta fit: R25 = 10030 ohm, beta = 3931 K; residuals [ 0.073  0.135 -0.507  0.226  0.076] C
  Steinhart-Hart fit: a = 1.2382e-03, b = 2.1734e-04, c = 1.430e-07; residuals [-0.044  0.219 -0.368  0.278 -0.087] C
  25 C point: raw 32351.0 +- 9.1 counts, R_fixed 10 kOhm +-1 % (rectangular): T = 24.63 C, u(T) = 0.131 C
  best fixed resistor for 0-50 C (beta model): 11011 ohm
```

and writes `bench_lab.png`: the thermistor's resistance against temperature with both fits (the points labelled SYNTHETIC data), the fit residuals, and the loading error against divider resistance.

### Verify

Run `python -m pytest curriculum/phase-04/lesson-15-electronics-bench-lab`. The tests check the divider's limits (0 Ω top, 0 Ω bottom), complementarity (swapping the resistors gives V_in − V_out) and scale invariance; the design round trip and the E12 rounding in all three directions; the LED current and power, and that rounding up never overdrives; the loading error against the Thevenin formula for four dividers to 10⁻¹², −1/21 at 1 MΩ, zero for an ideal meter, and monotone in R; the ammeter burden; the LSB and its factor-4 per two bits, the quantization noise on a ramp within 1 % and the error bounded by ½ LSB; clipping and the u16 scale; the calibration fit on an exact line; the beta model's round trip, reference point, negative slope and exact fit; the Steinhart–Hart inverse and fit recovering a known truth to 10⁻⁶; beta as Steinhart–Hart with c = 0; the ratiometric inverse, its midpoint and its domain; the geometric-mean resistor beating ±20 % alternatives; temperature from u16 to 10⁻⁹ °C; the multimeter Type B term; propagation against the product rule and a power law and zero-in, zero-out; the synthetic data reproducible by seed, different between seeds, monotone and in range; its analysis within 1 °C and its β within 3 % of the truth; the CSV files named and headed SYNTHETIC and read back; the figure.

### Compare

Three independent checks. Formula (2) against the circuit computed by `divider_output` with the meter in parallel: two different calculations of the same loading. The Steinhart–Hart fit against the beta fit on the same data: at the 25 °C point they differ by 0.14 °C (−0.507 against −0.368 °C residual), a measure of model error comparable to the reference thermometer's uncertainty, so neither is wrong. And the ADC against the multimeter: the calibration's residual standard deviation (1.20 mV on the synthetic data) is the honest resolution of your ADC chain, whatever its datasheet says.

## 4. Experiment at home, or investigate real data

**Mode: physical experiment** (optional). The bill of materials, the authoritative safety box and the full procedure are in `projects/p04-electronics-bench/` (`bom.yaml`, `safety.md`, `procedure.md`). Read the safety box first, every session. **Without the kit:** run the analysis on `projects/p04-electronics-bench/data/synthetic/` (SYNTHETIC data, generated by `make_synthetic_bench_data(seed=2026)`), or generate your own with your seed; label every table, plot and report "synthetic", and do not describe any of it as measured.

> **Safety (summary; the authoritative box is `projects/p04-electronics-bench/safety.md`).** Hazards: short-circuited batteries and USB ports (heat, burns, fire), an LED or resistor overheating, a meter fuse blown by a current-range mistake, water near electronics, sharp wire ends. Energy limits: **battery or USB power only, never mains; supply ≤ 9 V DC; current ≤ 500 mA (what a USB port supplies); no lithium-cell charging circuits; every LED through a resistor; no capacitors above 25 V or 1000 µF, and discharge any capacitor through a resistor before touching it; no soldering required** (if you choose to solder, the soldering-iron controls in the safety box apply). Controls: build and change circuits with the power disconnected; check with continuity before powering; feel for warm parts in the first minute; thermistor in a sealed bag, the board and breadboard on a dry bench away from the water; meter leads in the voltage jack except when deliberately measuring current. Supervision: none needed for an adult; children with an adult. Stop criteria: any part hot to the touch, any smell of burning, a USB port that disconnects, water on the bench, a reading outside the envelope. Not permitted: mains, wall adapters you have opened, lithium cells outside their device, high-voltage modules, anything above 9 V.

The bench layout (power off while building):

```
 USB 5 V or 3.3 V rail ──┬───────────────┬──────────────────┬───────────── (+) rail
                         │               │                  │
                       [R_t]          [330 Ω]            [NTC 10k]  ← in a zip bag in the mug
                         │               │                  │
                         ├── meter V     ▼ LED              ├────── ADC pin (A)
                         │               │                  │
                       [R_b]             │               [R_f 10k]
                         │               │                  │
 GND ────────────────────┴───────────────┴──────────────────┴───────────── (−) rail
           divider            LED (ammeter in series            thermistor divider
                              here for current)                 (powered from the ADC's reference)
```

Procedure in outline (full version in `procedure.md`):

1. **Predict first** in `data/prediction.md`: the divider output with its tolerance band; the 1 MΩ/1 MΩ reading with your meter's R_m (2.381 V for 10 MΩ from 5.00 V); the LED current with and without the ammeter's shunt; the ADC's LSB in volts; the thermistor ratio at each target temperature from its datasheet β.
2. **Resistances unpowered**: measure every resistor with the meter; record value and the meter's specification.
3. **Divider and loading**: build the 8.2 k/12 k divider, measure V_in and V_out; then the 1 MΩ/1 MΩ divider; compute R_m from (2).
4. **LED**: measure the supply, the resistor voltage, then the current with the meter in series (power off to rewire).
5. **ADC calibration**: potentiometer between the board's 3.3 V (or 5 V on a 5 V board) and GND, wiper to the ADC pin; at five settings record the meter's volts and the mean of 16 ADC readings.
6. **Thermistor sweep**: ratiometric divider; thermistor in a sealed bag in water at about 5, 15, 25, 35 and 45 °C (ice water and warm tap water mixed; never boiling or kettle-hot water); kitchen thermometer as reference; wait until the reading is steady; record 16 readings at each.

Raw data go in `data/raw/` in the formats of `data/README.md`, never edited afterwards.

## 5. Validate and improve

Put the three comparisons side by side.

**Prediction against measurement.** For each quantity, E_n = |pred − meas|/√(u_pred² + u_meas²) as in `P00.15`. A divider outside its tolerance band points at a wrong resistor (read its bands again, or measure it). A 1 MΩ divider reading near 2.38 V confirms (2) and your meter's 10 MΩ; a reading near 2.45 V means R_m ≈ 24.5 MΩ, which some meters have, and the formula still holds. An LED current lower than predicted by a few percent, with the meter in series, is the shunt; by much more, the LED's forward voltage is higher than the 2.0 V of the model (measure it).

**ADC against meter.** On the synthetic data the calibration gives a gain 0.68 % below ideal and an offset of −8.3 mV, with residuals of 1.2 mV, about 1.5 true LSB. On real data, residual structure (a curve rather than scatter) is nonlinearity: restrict the calibrated range or fit a quadratic and say so.

**Model against model.** On the synthetic data both thermistor fits have residuals within ±0.51 °C, consistent with a reference thermometer of ±0.5 °C. The fits cannot be better than the reference; if your residuals are much smaller than your thermometer's uncertainty, your thermometer is better than you said, or you have fitted too many parameters to too few points (three parameters to five points leaves two degrees of freedom).

| quantity | prediction (§2) | synthetic analysis | what limits it |
|---|---|---|---|
| divider output, 8.2 k/12 k at 5 V | 2.970 V | — (measure it) | resistor tolerance |
| 1 MΩ/1 MΩ reading, 10 MΩ meter | 2.381 V (−4.76 %) | — (measure it) | the meter's R_m |
| LED current, 330 Ω, 5 V | 9.09 mA; 8.82 mA with 10 Ω shunt | — (measure it) | V_f model, shunt |
| ADC gain | 50.35 µV/count (ideal) | 50.01 µV/count | the board's reference and ADC gain |
| ADC calibration residual | ≥ quantization 0.23 mV | 1.20 mV | noise, nonlinearity, meter |
| thermistor fit residuals | ≤ reference ±0.5 °C | ≤ 0.51 °C | the reference thermometer |
| u(T) at 25 °C | — | 0.131 °C | the fixed resistor's 1 % |

Improve and re-run: measure R_f and use the measured value (the budget's dominant term); average more samples only if the Type A term dominates (here it does not); recalibrate the ADC only over the range the thermistor uses. Then re-run `analyse_thermistor` and put the two budgets side by side.

## 6. Defend your work

Submit: `data/prediction.md`, the raw CSV files, your analysis script (using this lesson's functions), the calibration table with residuals, the two thermistor fits with residuals, the uncertainty budget for the 25 °C point, the loading measurement with the inferred R_m and its uncertainty, `bench_lab.png`, the completed safety box, and a report following `projects/p04-electronics-bench/report-template.md`. If you used the synthetic dataset, the title, the summary and every figure say SYNTHETIC.

### Adversarial self-checks

- **What if the divider resistors are ten times larger?** The output is unchanged in the model, but the loading error with a 10 MΩ meter grows tenfold (0.05 % → 0.5 % at 100 kΩ), and the ADC's own sampling input starts to load the divider too. Larger is not free.
- **What if the supply goes to zero?** Every divider output goes to zero, the LED current goes to zero below V_f (not linearly: the fixed-V_f model's switch-off), and the ratiometric thermistor reading becomes 0/0: the code must refuse a ratio of 0 or 1, and it does.
- **Which assumption fails first in the real world?** The LED's fixed forward voltage (V_f rises with current, so the real current is a little lower than the model's), then the thermistor's thermal equilibrium (lag and self-heating), then the ADC's linearity near the rails.
- **What measurement could falsify the loading model?** Reading the 1 MΩ/1 MΩ divider with two different meters: (2) predicts a different reading for each that is consistent with each meter's R_m; a reading above 2.500 V would falsify it (a passive meter cannot raise the voltage).
- **Does the result survive a change in model or instrument?** The temperature at 25 °C changes by 0.14 °C between beta and Steinhart–Hart and by 0.11 °C between the nominal and the true R_f; a change of ADC board changes the calibration, not the temperature, because the measurement is ratiometric.
- **Which decision changes if the uncertainty doubles?** At u(T) = 0.26 °C the thermistor still resolves a 1 °C change in a satellite box with margin, but a claim that Steinhart–Hart is better than beta on these data is no longer supported (their difference is within the uncertainty).

## Check yourself

1. Which way is a voltmeter connected, which way an ammeter, and why must resistance be measured unpowered?
2. Design an LED resistor for a red LED (V_f = 2.0 V) at 10 mA from 5 V. Which E12 value do you choose, and why round up?
3. A 10 MΩ voltmeter reads the midpoint of a 1 MΩ/1 MΩ divider from 5 V. What does it read, and what is R_th?
4. What is the LSB of a 12-bit ADC at 3.3 V, and the rms quantization noise? Does `read_u16()` returning values up to 65 535 mean 16-bit resolution?
5. Why does V_ref cancel in the ratiometric thermistor divider, and what does not cancel?
6. Write the beta model and the Steinhart–Hart equation. Which is linear in which parameters?
7. In the 25 °C budget, which input dominates u(T), and what single action would reduce it most?
8. Name three hard electrical limits of this lab and the hazard each prevents.
