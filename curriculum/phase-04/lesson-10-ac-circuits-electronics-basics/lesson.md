# P04.10  AC Circuits & Electronics Basics

> A capacitor is a resistor whose value depends on frequency. Put one after a resistor and you have the filter that decides whether a sensor's readings mean anything after a computer has sampled them.

| | |
|---|---|
| **Learning outcomes** | 1. Represent R, L and C as complex impedances, analyse a series RLC circuit with phasors, and compute its resonant frequency, quality factor and bandwidth and the real power, reactive power and power factor of a load on a 115 V, 400 Hz aircraft bus. 2. Derive the transfer function, corner frequency, phase and group delay of an RC low-pass filter, show why a second section must be buffered, and verify the response against `scipy.signal`. 3. Apply the two ideal-op-amp rules to the buffer, the non-inverting and the inverting amplifier, and show the closed-loop gain A/(1 + Aβ) tends to the ideal result as the open-loop gain grows; describe a diode and a transistor switch at the conceptual level. 4. Explain what a microcontroller ADC returns (counts, LSB, quantisation error, aliasing), and design an anti-aliasing filter for a sensor sampled at 100 Hz by a 12-bit ADC that keeps a 430 Hz vibration below one LSB, checked in the time domain on SYNTHETIC signals. |
| **Prerequisites** | `P04.08` (Ohm's law, Kirchhoff's rules, the RC time constant), `P02.07` (complex numbers and phasors), `P02.08` (the Laplace transform and transfer functions), `P02.03` (forced response and resonance), `P02.10` (the FFT, the sampling theorem), `P00.14` (measurement uncertainty, needed for the optional measurement). Skill check: you can multiply and divide complex numbers in polar form and solve a 2 × 2 linear system. |
| **Engineering question** | A pressure sensor on a propeller aircraft carries a slow signal (0–5 Hz) and a 50 mV vibration pickup at 430 Hz; it is sampled at 100 Hz by a 12-bit ADC with a 3.3 V reference. Unfiltered, the vibration aliases to a false 30 Hz tone of 62 counts. What RC filter in front of the ADC brings that alias below one count (0.81 mV) while costing less than 0.5 dB at 5 Hz, and what delay does it add? Answer: one RC section cannot (its corner would have to be 6.9 Hz, losing 1.8 dB at 5 Hz); two buffered sections of 3.3 kΩ and 1 μF (corner 48.2 Hz) give 0.77 counts of alias, 0.09 dB of droop and 6.6 ms of delay. |
| **Assumptions and limits** | Linear, time-invariant circuits with ideal lumped R, L, C (no parasitic resistance of the inductor, no leakage of the capacitor, no stray capacitance), sinusoidal steady state for phasors. The ideal op-amp: infinite gain, infinite input impedance, zero output impedance, unlimited bandwidth and swing; real parts fail at high frequency, near the supply rails and above their output current. The ADC is ideal: uniform steps, no noise, instantaneous sampling. The vibration is a single sinusoid; real vibration is broadband and the design must then cover a band. Diodes and transistors are described, not modelled. |
| **Analytical work** | Impedances of R, L and C from their differential laws; series RLC: resonance, Q, bandwidth; RMS and the power triangle; the RC divider as H(jω) = 1/(1 + jωRC): corner, phase, roll-off, group delay; the loaded two-section ladder 1/(1 + 3x + x²) by nodal analysis; the attenuation-and-droop design rule for n buffered sections; the virtual short and the non-inverting gain 1 + R_f/R_g; A/(1 + Aβ); aliasing |f − f_s round(f/f_s)|; LSB and the LSB/√12 quantisation noise. |
| **Computational deliverable** | `starter/main.py`: impedances, `series`, `parallel`, `rlc_series_current`, `resonant_frequency`, `quality_factor`, `rms`, `ac_power`, `rc_corner_frequency`, `rc_lowpass`, `rc_highpass`, `rc2_unbuffered`, `rc_loaded`, `gain_db`, `group_delay_dc`, `rc_step_response`, `rc_lti`, `required_stages_and_corner`, the op-amp gains and `closed_loop_gain`, `thermal_voltage`, `diode_current`, `bjt_base_resistor`, `gate`, `half_adder`, `aliased_frequency`, `adc_lsb`, `adc_convert`, `counts_to_volts`, `quantisation_noise_rms`, `filter_and_sample` (`synthetic_sensor_signal`, `tone_amplitude` and `main` are given). Tests prove: impedance phases ±90° and the 1/f and f laws; the series RLC current is V/R at f₀ with −3 dB width f₀/Q; RMS of a sine is A/√2; P = I²R and pf = cos arg Z; the RC response is −3.01 dB and −45° at 1/(2πRC), rolls off at 20n dB/decade and matches `scipy.signal.freqs` to 10⁻¹⁰; the step response matches `scipy.signal.step`; the group delay is nRC; the unbuffered ladder equals nodal analysis; the design rule meets attenuation exactly and the droop bound; the closed-loop gain error falls as G²/A; 60 mV per decade in the diode; truth tables, NAND universality, De Morgan and the half adder; aliases fold into [0, f_s/2]; quantisation error within ±½ LSB with RMS LSB/√12; the 30 Hz alias in SYNTHETIC sampled data falls from 62 to under 1 LSB with two stages. |
| **Practical activity** | **computational experiment**: an RLC resonance sweep, the 400 Hz bus power study, the anti-aliasing design for the 100 Hz, 12-bit channel with a time-domain check on SYNTHETIC signals, and a buffered-versus-loaded comparison. An optional measurement of an RC step response with a microcontroller ADC uses the `projects/p04-electronics-bench` kit and is clearly marked optional. |
| **Equipment** | None required: a computer with the environment from `P00.01`. Optional (Level 1): the `projects/p04-electronics-bench` kit (breadboard, a 10 kΩ resistor, a 10 μF capacitor rated ≥ 16 V, jumper wires, a 3.3 V or 5 V microcontroller board with an ADC, USB cable). |
| **Safety** | None beyond normal computer use for the required work. The optional measurement follows `projects/p04-electronics-bench/safety.md`: USB power only, never mains; ≤ 9 V DC; no capacitor above 25 V or 1000 μF; observe electrolytic polarity; discharge capacitors before handling. |
| **Evidence** | Completed `main.py`, passing tests, `rc_filter.png`, the derivations of the RC transfer function and its corner, of the loaded ladder and of the non-inverting gain, the four experiment tables, and a half-page design note for the anti-aliasing filter with component values, attenuation, droop and delay; if the optional measurement was done, the raw CSV and the fitted time constant with its uncertainty. |
| **Assessment** | `assessments/phase-04/grading-rubric.md`; competencies 1–4. Pass: the RC response and the ideal-op-amp gains are derived; all public tests pass; the design note shows that one stage fails and two pass, with numbers; the SYNTHETIC check is labelled as such and its agreement with the prediction is explained. |
| **References** | `@openstax-university-physics-vol2` Ch. 14 "Inductance": 14.5 Oscillations in an LC Circuit, 14.6 RLC Series Circuits; Ch. 15 "Alternating-Current Circuits": 15.1 AC Sources, 15.2 Simple AC Circuits, 15.3 RLC Series Circuits with AC, 15.4 Power in an AC Circuit, 15.5 Resonance in an AC Circuit (required-free); 15.6 Transformers (supplementary-free). `@openstax-university-physics-vol3` Ch. 9 "Condensed Matter Physics": 9.6 Semiconductors and Doping, 9.7 Semiconductor Devices (required-free: the p-n junction, the diode and the transistor); 9.5 Band Theory of Solids (supplementary-free). `@fiore-opamps-linear-ics` Ch. 4 "Basic Op Amp Circuits": 4.2 Inverting and Noninverting Amplifiers (required-free); Ch. 12 "Analog-to-Digital-to-Analog Conversion": 12.2 The Sampling Theorem, 12.3 Resolution and Sampling Rate (required-free). `@kahn-digital-circuit-projects` Ch. 4 "Gates", Ch. 6 "Adders" (required-free); Ch. 9 "Memory Basics - Flip-Flops and Latches" (supplementary-free). `@sparkfun-learn-tutorials` "Digital Logic": "Combinational Logic", "Sequential Logic" (D-type Flip-Flop); "What is an Arduino?": "What's on the board?"; "Analog to Digital Conversion": "What is the ADC?", "Relating ADC Value to Voltage" (required-free). `@micropython-docs` "Quick reference for the RP2": "ADC (analog to digital conversion)" (supplementary-free; for the optional measurement). |
| **Extensions** | ★ Replace the two RC sections by a unity-gain Sallen–Key second-order section and show that, for the same 430 Hz attenuation, a Butterworth response (Q = 0.707) loses less at 5 Hz than two equal RC sections (Q = 0.5). ★ Oversample at 1 kHz and decimate by 10 with a digital FIR filter (`scipy.signal.decimate`): how much analog filtering is still needed? ★ Add white noise from 0 to 10 kHz to the SYNTHETIC signal and compute how much aliased noise power each design admits. |
| **Time** | unmeasured |

## 1. The engineering problem

Every number a flight computer knows about the physical world arrives as a voltage, is filtered, is sampled by an analog-to-digital converter (ADC), and becomes an integer. The pitot-static pressure that becomes the airspeed, the strain gauge on a wing spar, the accelerometers of an inertial unit, the thermistors on a battery pack: all follow that chain. Most of its failures are not in the sensor but in the step between the wire and the integer. An ADC that samples at 100 Hz cannot tell a 430 Hz vibration from a 30 Hz one: the samples are identical. Once sampled, the false 30 Hz tone is in the data and no software can remove it, because the information that would distinguish it was thrown away at the sampling instant. The only remedy is to remove the high frequencies before the ADC, with an analog low-pass filter, and the cheapest analog filter is a resistor and a capacitor.

To design one you need the language of AC circuits: sinusoidal signals as phasors, capacitors and inductors as frequency-dependent impedances, resonance, RMS power. The same language describes the 115 V, 400 Hz AC bus of a transport aircraft, so the lesson starts there. Then it builds the RC filter, finds why two filter sections need an amplifier between them (the op-amp buffer), meets the diode and the transistor that make digital electronics possible, the logic gates built from them, and finally the ADC and what it returns.

Success criteria, written first: the RC response matches `scipy.signal.freqs` to 10⁻¹⁰; the series RLC peaks at 1/(2π√LC) with bandwidth f₀/Q; the design rule picks the smallest number of sections that meets both the 430 Hz attenuation (62 to 1) and the 0.5 dB droop limit at 5 Hz; in a time-domain simulation with SYNTHETIC signals the 30 Hz alias in the ADC counts falls from about 62 counts to below 1.

## 2. Understand and derive

### 2.1 Phasors and impedance

A sinusoidal voltage v(t) = V_p cos(ωt + φ), ω = 2πf, is the real part of V e^{jωt} with the complex amplitude, or **phasor**, V = V_p e^{jφ} (`P02.07`). In a linear circuit driven at one frequency, every voltage and current is a sinusoid at that frequency, and d/dt becomes multiplication by jω. The element laws become Ohm's law with complex resistances, the **impedances** (OpenStax 15.2):

- resistor, v = Ri: **Z_R = R**, voltage and current in phase;
- inductor, v = L di/dt: **Z_L = jωL**, voltage leads current by 90°, |Z| rises with f;
- capacitor, i = C dv/dt: **Z_C = 1/(jωC)**, voltage lags current by 90°, |Z| falls with f, infinite at DC.

Impedances combine like resistances: Z = Z₁ + Z₂ in series and 1/Z = 1/Z₁ + 1/Z₂ in parallel, and Kirchhoff's laws of `P04.08` hold for phasors. Everything in this lesson is that sentence applied.

### 2.2 The series RLC circuit and resonance

A source V drives R, L and C in series: I = V/(R + j(ωL − 1/(ωC))). The reactances cancel at the **resonant frequency**

**f₀ = 1/(2π√(LC)),**   (1)

where the circuit is purely resistive and |I| = |V|/R is largest (15.3, 15.5). Below f₀ the capacitor dominates and the current leads; above, the inductor dominates and it lags. The sharpness is the **quality factor** Q = ω₀L/R = (1/R)√(L/C), and the band where |I| exceeds its peak over √2 has width **f₀/Q**. With R = 10 Ω, L = 10 mH, C = 1 μF: f₀ = 1591.5 Hz, Q = 10.00, f₀/Q = 159.2 Hz; a numerical sweep on a 10 mHz grid gives the peak at 1591.5 Hz and a −3 dB width of 159.1 Hz. (The band is not symmetric about f₀, its edges lie at a geometric rather than arithmetic mean, but for the series RLC its width is exactly R/(2πL) = f₀/Q at any Q; the sweep differs only by its grid.) Energy sloshes between the inductor's field and the capacitor's field (14.5) and R removes a fraction 2π/Q of it per cycle. Resonance is wanted in a radio tuner and unwanted in a power-line filter, where an undamped LC section amplifies whatever arrives at f₀ (Problem 3).

### 2.3 RMS and AC power on a 400 Hz bus

A periodic voltage heats a resistor as much as a DC voltage equal to its root-mean-square value, V_rms = √(mean v²); for a sine, V_rms = V_p/√2 (15.4). AC power ratings are RMS: transport-aircraft AC buses run at 115 V rms, 400 Hz (the high frequency makes transformers and generators lighter, since core size falls with frequency for the same power). With phasors in RMS units, the complex power drawn by a load Z is S = V I* = |V|²/Z*: the real part P (watts) is dissipated, the imaginary part Q (volt-ampere reactive, var) sloshes back and forth, and the **power factor** P/|S| = cos(arg Z).

A load of R = 10 Ω in series with L = 2 mH at 400 Hz: X_L = ωL = 5.027 Ω, |Z| = 11.192 Ω, phase 26.7°, I = 10.275 A rms, **P = 1056 W, Q = 531 var, pf = 0.893**. The wiring must carry 10.3 A although only 92 % of what a 1056 W resistive load would need is useful, because the current is out of phase: that is why aircraft electrical loads have power-factor limits and why a capacitor can correct them (Problem 1).

### 2.4 The RC low-pass filter

A resistor R in series and a capacitor C to ground form a voltage divider with the output across C:

**H(jω) = V_out/V_in = Z_C/(R + Z_C) = 1/(1 + jωRC).**   (2)

At low frequency |H| → 1; at high frequency |H| → 1/(ωRC), falling by a factor 10 per decade of frequency (**−20 dB/decade**). The two regimes meet at the **corner frequency**

**f_c = 1/(2πRC),**   (3)

where |H| = 1/√2 (−3.01 dB) and the phase is −45°. In Laplace form (`P02.08`) H(s) = 1/(RCs + 1), and the step response is the RC charging curve of `P04.08`, 1 − e^{−t/RC}. The phase −arctan(ωRC) is close to −ωRC at low frequency, so slow signals come out delayed by RC: the **group delay** −dφ/dω = RC. For a control loop, that delay is the price of the filter. The swap of R and C makes the high-pass filter jωRC/(1 + jωRC), the complement of (2).

**Two sections need a buffer.** Cascading two RC sections should give H² = 1/(1 + jωRC)², −40 dB/decade. Connected directly, it does not: the second section draws current from the first and changes its output. Nodal analysis with x = jωRC on the two capacitor nodes gives

**H_unbuffered = 1/(1 + 3x + x²),**   (4)

against 1/(1 + 2x + x²) for the buffered pair. With R = 3.3 kΩ and C = 1 μF the loaded pair loses 0.316 dB at 5 Hz instead of 0.093 dB and has its −3 dB point at 18.0 Hz instead of 31.0 Hz; the high-frequency attenuation is about the same. The fix is an amplifier with very high input impedance and very low output impedance between the sections: the op-amp buffer of 2.5.

### 2.5 The ideal op-amp

An operational amplifier is a differential amplifier with a very large gain A (10⁵–10⁶ at low frequency): v_out = A(v₊ − v₋). Used with **negative feedback** (part of the output returned to the inverting input), its behaviour is set by the feedback network, not by A, and two rules describe it (Fiore 4.2):

1. **No current flows into either input** (input impedance infinite).
2. **With negative feedback, the output makes v₋ = v₊** (the virtual short).

**Buffer**: output tied to v₋, input to v₊: v_out = v_in, gain 1, but the input draws no current and the output can drive a load. **Non-inverting amplifier**: R_f from output to v₋, R_g from v₋ to ground; rule 2 gives v₋ = v_in, rule 1 makes R_f and R_g a divider, v_in = v_out R_g/(R_f + R_g):

**G = 1 + R_f/R_g.**   (5)

**Inverting amplifier**: input through R_in to v₋, R_f to the output, v₊ grounded: v₋ is a virtual ground and G = −R_f/R_in. To see why the rules work, keep A finite: with feedback fraction β = R_g/(R_f + R_g), v_out = A(v_in − βv_out), so G = A/(1 + Aβ) → 1/β as A → ∞, with error ≈ G²/A. For R_f = 27 kΩ, R_g = 10 kΩ: ideal G = 3.70, with A = 10⁵ G = 3.69986, with A = 100 G = 3.5680. A sensor whose output spans 0–0.8 V would be amplified by this stage to span most of a 3.3 V ADC range, using more of its counts.

### 2.6 Diodes, transistors, and logic

Doping silicon with a little phosphorus gives free electrons (n-type); with boron, free holes (p-type) (OpenStax Vol. 3, 9.6). At a **p-n junction** the carriers diffuse across and leave a depleted region with a built-in field; a forward voltage lowers the barrier and the current grows exponentially, a reverse voltage raises it and only a tiny saturation current flows (9.7). The Shockley equation I = I_s(e^{V/(nV_T)} − 1), with V_T = kT/q = 25.85 mV at 300 K, puts numbers on it: with I_s = 10⁻¹² A, 0.60 V gives 12.0 mA and 0.66 V gives 122 mA: ten times the current for 60 mV more. A diode is a one-way valve; across a coil it is the flyback path that absorbs the L di/dt spike of `P04.09`.

A **bipolar transistor** is two junctions back to back; a small base current controls a collector current β times larger (β ≈ 100). Driven hard, it is a switch: to switch the 0.12 A torquer coil of `P04.09` from a 3.3 V microcontroller pin with β = 100 and a fivefold overdrive, the base resistor is (3.3 − 0.7)/(5 × 0.12/100) = **433 Ω**. A field-effect transistor (MOSFET) does the same with a voltage on its gate and is the switch in every logic chip.

Switches make **logic gates**: AND, OR, NOT and their combinations NAND, NOR, XOR, each defined by a truth table (Kahn, Ch. 4 "Gates"; SparkFun "Combinational Logic"). NAND alone is universal (NOT a = a NAND a, and every other gate follows). Two gates make a **half adder**, sum = a XOR b, carry = a AND b (Kahn, Ch. 6 "Adders"), and chains of adders make arithmetic. Feed a gate's output back to its input and the circuit remembers: the **D flip-flop** stores one bit at each clock edge (Ch. 9; SparkFun "Sequential Logic"), and registers, counters and memory are flip-flops in arrays. A microcontroller is these pieces around a processor: digital input/output pins, timers, serial buses, and an ADC (SparkFun "What's on the board?").

### 2.7 What an ADC returns

An n-bit ADC with reference v_ref divides [0, v_ref) into 2ⁿ steps of one **LSB** (least significant bit) = v_ref/2ⁿ and returns the integer step number (Fiore 12.3; SparkFun "Relating ADC Value to Voltage"). A 12-bit ADC at 3.3 V has LSB = 0.8057 mV and returns 0–4095; 1.000 V gives 1241. Converting back with (counts + ½) LSB leaves an error within ±½ LSB, uniformly spread for a busy signal, with RMS **LSB/√12** = 0.2326 mV. Some software interfaces rescale the counts to a fixed range (MicroPython's `ADC.read_u16()` returns 0–65535 whatever the converter's resolution); the resolution is still the hardware's.

The ADC samples at f_s. By the sampling theorem (Fiore 12.2; `P02.10`), only frequencies below the **Nyquist frequency** f_s/2 are represented unambiguously; a sinusoid at f produces exactly the same samples as one at

**f_alias = |f − f_s round(f/f_s)|,**   (6)

which lies in [0, f_s/2]. At f_s = 100 Hz, 430 Hz becomes 30 Hz, 70 Hz becomes 30 Hz, 400 Hz becomes 0 Hz (a constant offset). The ADC cannot know.

### 2.8 Designing the anti-aliasing filter

The channel: signal 0–5 Hz, f_s = 100 Hz, 12 bits, 3.3 V. The vibration pickup, 50 mV at 430 Hz, is 62.1 LSB. Requirement: the residual at 430 Hz below 1 LSB, so |H(430 Hz)| ≤ 1/62.1 (−35.9 dB); and at most 0.5 dB loss at 5 Hz. For n identical buffered RC sections |H(f)| = (1 + (f/f_c)²)^{−n/2}; the attenuation requirement gives the largest allowed corner

**f_c = f_reject/√(a^{2/n} − 1),**   (7)

with a = 62.1, and the droop at 5 Hz follows:

| sections | largest f_c | droop at 5 Hz |
|---|---|---|
| 1 | 6.93 Hz | 1.820 dB |
| 2 | 55.03 Hz | 0.071 dB |
| 3 | 112.25 Hz | 0.026 dB |

One section fails the droop requirement; **two sections** meet both. Choosing standard values below the 55 Hz limit, R = 3.3 kΩ and C = 1 μF give f_c = 48.23 Hz, |H(430 Hz)| = 0.01242 (−38.1 dB), a residual of 0.621 mV = **0.77 LSB**, a droop of **0.093 dB** at 5 Hz and a low-frequency delay of 2RC = **6.6 ms**. Note that the corner lies close to the Nyquist frequency: the filter is sized for the 430 Hz tone, and broadband noise between 50 and 100 Hz would still alias; a real design specifies a band to reject, not one tone.

### 2.9 Where the models stop

Real capacitors are ±10–20 % (ceramic, electrolytic) and change with temperature and voltage; the corner moves accordingly, and the design above has a factor 55/48 of margin. The ideal op-amp fails above its gain-bandwidth product (a 1 MHz part has A ≈ 2300 at 430 Hz: still ideal here) and near the rails (a buffer on a 3.3 V supply needs a rail-to-rail part). An ADC's sample-and-hold capacitor draws a pulse of current at every conversion; the buffer drives it, an RC alone with 3.3 kΩ might not settle. The single-tone vibration is the weakest assumption.

## 3. Build the computational model

Open `starter/main.py`. `synthetic_sensor_signal`, `tone_amplitude` and `main` are given.

### Step 1: impedances and the circuit algebra

`impedance_R`, `impedance_L`, `impedance_C` (return `inf` at f = 0 for the capacitor, without a division warning), `series`, `parallel`. `rlc_series_current`, `resonant_frequency`, `quality_factor` follow from 2.2; `rms` and `ac_power` from 2.3.

### Step 2: the filter

`rc_corner_frequency` is (3); `rc_lowpass(f, R, C, stages)` is (2) raised to the number of buffered stages; `rc_highpass`; `rc2_unbuffered` is (4); `rc_loaded` is the divider with C ∥ R_load; `gain_db`, `group_delay_dc`, `rc_step_response`. `rc_lti` builds the `scipy.signal.lti` object with denominator (RCs + 1)ⁿ by repeated `np.polymul`. `required_stages_and_corner` loops n = 1, 2, … applying (7) and the droop check.

### Step 3: op-amp, diode, transistor, gates

One line each from 2.5 and 2.6. `gate(name, a, b)` on 0/1 integer arrays with `&`, `|`, `^`; `half_adder`.

### Step 4: the ADC and the time-domain check

`aliased_frequency` is (6); `adc_lsb`, `adc_convert` (floor and clip), `counts_to_volts`, `quantisation_noise_rms`. `filter_and_sample(t, v, R, C, stages, fs)` runs `scipy.signal.lsim` on the finely sampled SYNTHETIC waveform (20 kHz), starting the filter in its DC steady state, then keeps every (20 000/f_s)-th point and quantises it:

```python
ss = rc_lti(R, C, stages).to_ss()
x0 = np.linalg.solve(ss.A, -ss.B[:, 0] * v[0])
_, vf, _ = signal.lsim(ss, v, t, X0=x0)
```

Running `main()` prints

```
400 Hz bus, 115 V rms, load R = 10.0 ohm + L = 2.0 mH: X_L = 5.027 ohm, |Z| = 11.192 ohm, phase 26.7 deg, I = 10.275 A rms, P = 1056 W, Q = 531 var, pf = 0.893
rms of a 115 V rms sine sampled over one period: 115.0000 V

series RLC R = 10.0 ohm, L = 10.0 mH, C = 1.0 uF: f0 = 1591.5 Hz, Q = 10.00, f0/Q = 159.2 Hz; swept peak at 1591.5 Hz, -3 dB width 159.1 Hz

ADC: 12 bits, 3.3 V: LSB = 0.8057 mV, quantisation noise 0.2326 mV rms
vibration 430 Hz sampled at 100 Hz appears at 30 Hz; 50 mV = 62.1 LSB, so the filter must attenuate by 62.1 (35.9 dB)
  1 stage(s): largest corner   6.93 Hz, droop at 5 Hz 1.820 dB
  2 stage(s): largest corner  55.03 Hz, droop at 5 Hz 0.071 dB
  3 stage(s): largest corner 112.25 Hz, droop at 5 Hz 0.026 dB
design: 2 buffered stages, corner <= 55.03 Hz; chosen R = 3.3 kohm, C = 1 uF: f_c = 48.23 Hz, |H(430)| = 0.01242 (-38.1 dB), residual 0.621 mV = 0.77 LSB, droop at 5 Hz 0.093 dB, delay 6.6 ms
scipy.signal.freqs check: |H| = 0.989366, 0.012424; ours 0.989366, 0.012424
the same two sections without a buffer: |H(430)| = 0.01206, droop at 5 Hz 0.316 dB, -3 dB at 18.01 Hz
```

then the op-amp, diode and transistor lines of 2.5–2.6, and the SYNTHETIC check:

```
SYNTHETIC test signal (not a measurement): 1.65 V + 0.5 V at 2 Hz + 50 mV at 430 Hz, sampled at 100 Hz, 12 bits
  0 stage(s): 30 Hz alias in the ADC counts   62.22 LSB (predicted  62.06 LSB)
  1 stage(s): 30 Hz alias in the ADC counts    7.05 LSB (predicted   6.92 LSB)
  2 stage(s): 30 Hz alias in the ADC counts    0.69 LSB (predicted   0.77 LSB)
```

and writes `rc_filter.png` (gain curves of one, two and three buffered stages and the unbuffered pair with the 430 Hz requirement marked; the RLC resonance; the ADC counts with and without the filter).

### Verify

Run `python -m pytest curriculum/phase-04/lesson-10-ac-circuits-electronics-basics`. The tests:

1. **Impedance**: phases 0, +90°, −90°; |Z_L| ∝ f, |Z_C| ∝ 1/f; series and parallel.
2. **RLC**: I = V/R exactly at f₀; Q = ω₀L/R; the −3 dB width f₀/Q within 0.2 % on an 8 mHz grid; leading below f₀, lagging above.
3. **Power**: RMS of a sine A/√2; P = I²R, Q = I²X, pf = cos arg Z; a capacitive load has Q < 0.
4. **RC**: −3.0103 dB and −45° at f_c; unit DC gain; low-pass + high-pass = 1.
5. **Roll-off**: 20n dB per decade for n = 1, 2, 3.
6. **Oracle**: `rc_lowpass` against `scipy.signal.freqs` for n = 1–3 to 10⁻¹⁰, and against an independently built (RCs + 1)².
7. **Step and delay**: the closed form against `scipy.signal.step`; −dφ/dω = nRC.
8. **Loading**: (4) against a 2 × 2 nodal solve at three frequencies; the loaded pair below the buffered pair at f_c; the loaded single stage's DC divider and its infinite-load limit.
9. **Design rule**: two stages for the lesson's numbers, attenuation met exactly at the returned corner, droop below 0.5 dB; one stage for a lax droop spec; an impossible spec raises.
10. **Op-amp**: 3.7, 1 and −10; the closed-loop error falling monotonically to < 10⁻⁷ and ≈ G²/A.
11. **Diode and transistor**: V_T, ten times the current per V_T ln 10, −I_s in reverse; the base current.
12. **Logic**: the five truth tables, NAND universality, De Morgan, the half adder's 2c + s = a + b.
13. **Aliasing**: 430 → 30, 70 → 30, 400 → 0, all aliases in [0, f_s/2]; identical samples of 430 and 30 Hz cosines.
14. **ADC**: clipping; error within ½ LSB; RMS error LSB/√12 within 3 % (seeded uniform voltages).
15. **SYNTHETIC end to end**: the measured alias tone within 0.6 LSB of the prediction for 0, 1, 2 stages; above 50 LSB without filter and below 1 with two stages; the 2 Hz signal passed with the predicted gain.
16. **main** writes the figure.

### Compare

`scipy.signal` is the oracle. `freqs(num, den, worN)` evaluates an analog transfer function from its polynomial coefficients and agrees with (2)ⁿ to 10⁻¹⁰; `step` and `lsim` integrate the state-space form of the filter, and `lsim` is what turns the frequency-domain design into a time-domain check: on the SYNTHETIC signal, the alias tone measured in the ADC counts (62.22, 7.05, 0.69 LSB) agrees with the prediction from |H(430 Hz)| (62.06, 6.92, 0.77 LSB) to better than one LSB, the residual being quantisation itself (the alias is now smaller than the step it is measured in). What the hand-built functions add is the design rule (7) and the explanation of each number; what `scipy.signal` adds is the guarantee that the algebra was done right and, for any filter more complicated than RC sections (`butter`, `cheby1`, `bessel`), the design itself.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment**, with an optional measurement. Write each prediction first.

1. **Resonance sweep.** For the series RLC with R = 10 Ω, L = 10 mH, C = 1 μF and then R = 1 Ω and 100 Ω, sweep 100 Hz–10 kHz, record f₀, the peak current for 1 V, the −3 dB width and the phase at f₀ ± bandwidth/2. Predict f₀ unchanged, width ∝ R, and phases of ±45°.
2. **Power factor.** On the 400 Hz bus, sweep the load inductance 0–5 mH with R = 10 Ω; record P, Q, pf and line current. Predict the line current needed for a fixed real power to grow as 1/pf.
3. **Anti-aliasing design.** Repeat 2.8 for vibration at 230, 430 and 830 Hz and droop limits 0.1, 0.5 and 1 dB; record the number of sections, the corner, the chosen R and C (from the E12 series, C ∈ {0.1, 0.22, 0.47, 1} μF), the residual in LSB, and the delay. Predict which cases need three sections. Check each design with `filter_and_sample` on the SYNTHETIC signal.
4. **Buffered versus loaded.** For the chosen design, the gain at 5 Hz, 48 Hz and 430 Hz buffered and unbuffered, and a loaded single stage feeding a 10 kΩ input. Predict where loading matters (in the passband and near the corner, not far above it).

**Optional measurement (physical; uses the `projects/p04-electronics-bench` kit; skip it without loss).** Safety summary first: USB power from a computer or a USB power bank only, never mains; nothing above 9 V DC; a 10 μF capacitor rated at least 16 V, inserted with its marked negative lead to ground; touch the circuit only with the USB cable unplugged; the full safety box is `projects/p04-electronics-bench/safety.md`. Build R = 10 kΩ in series from a digital output pin to a 10 μF capacitor to ground, and connect the capacitor's top to an ADC pin. **Prediction, written before measuring:** τ = RC = 0.100 s ± 20 % (the capacitor's tolerance dominates), corner 1.59 Hz. Set the pin high and read the ADC as fast as the board allows for 0.6 s, then set it low and repeat. Fit v(t) = V(1 − e^{−t/τ}) with `scipy.optimize.curve_fit`, and report τ with its uncertainty from the fit and from the component tolerances (`P00.14`). Raw data format: CSV with header `t_s, counts`, plus a header line stating the board, its ADC bits and v_ref, R and C as marked, and the date.

Raw data format for the computational parts: CSV tables `R_ohm, f0_Hz, I_peak_A, width_Hz` (1); `L_mH, P_W, Q_var, pf, I_A` (2); `f_vib_Hz, droop_limit_dB, stages, fc_Hz, R_ohm, C_F, residual_LSB, delay_ms, synthetic_alias_LSB` (3); `f_Hz, H_buffered_dB, H_unbuffered_dB, H_loaded_dB` (4); plus `rc_filter.png`. Every table from the SYNTHETIC signal carries the word SYNTHETIC in its header.

What you should find: widths 15.9, 159.1 and 1591 Hz for R = 1, 10, 100 Ω (at R = 100 Ω, Q = 1 and the "band" is barely a peak), phases ±45.0° at the band edges; pf falling from 1 to 0.97, 0.893, 0.798, 0.705 and 0.623 as L goes from 0 to 5 mH, with P falling from 1322 to 513 W; three sections needed only at 230 Hz with the 0.1 dB limit (corner ≤ 60.0 Hz, for example 27 kΩ and 0.1 μF, 8.1 ms delay); at 430 Hz two sections for every droop limit (for example 6.8 kΩ and 0.47 μF, 49.8 Hz, 0.82 LSB, 6.4 ms); at 830 Hz two sections for 0.1 and 0.5 dB but one section (120 kΩ, 0.1 μF, 13.3 Hz, 12 ms delay) passes the 1 dB limit. The SYNTHETIC check reproduces each predicted residual within about 0.15 LSB, which means a design predicted at 0.95–1.0 LSB can read 1.05–1.1 LSB: leave margin below 1. Loading is visible at 5 and 48 Hz (0.316 against 0.093 dB at 5 Hz) and negligible at 430 Hz (0.01206 against 0.01242).

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| RLC f₀ | 1/(2π√LC) = 1591.5 Hz | swept peak 1591.5 Hz | < 0.1 Hz | sweep grid |
| RLC −3 dB width | f₀/Q = 159.15 Hz | 159.1 Hz | < 0.05 Hz | sweep grid |
| RMS of a sine | 115 V | 115.0000 V | 0 | — |
| bus power factor | cos 26.7° = 0.893 | 0.893 | 0 | — |
| RC at f_c | −3.01 dB, −45° | same to 10⁻¹² | 0 | — |
| two stages at 430 Hz | 0.012424 | `freqs` 0.012424 | 10⁻¹⁰ | — |
| unbuffered droop at 5 Hz | (4): 0.316 dB | nodal solve: same | 10⁻¹² | — |
| alias, no filter | 62.06 LSB | SYNTHETIC 62.22 LSB | 0.16 LSB | quantisation |
| alias, one stage | 6.92 LSB | SYNTHETIC 7.05 LSB | 0.13 LSB | quantisation |
| alias, two stages | 0.77 LSB | SYNTHETIC 0.69 LSB | 0.08 LSB | the tone is below one step |
| non-inverting gain, A = 10⁵ | 3.70 | 3.69986 | G²/A = 1.4 × 10⁻⁴ | finite A |

Everything agrees to the precision the theory claims, and the SYNTHETIC time-domain check shows the design does what the frequency response promised, inside the resolution of the ADC itself. These comparisons are between a derivation and a simulation, not a measurement: the SYNTHETIC signal tests the analysis pipeline. The engineering answer is **two buffered RC sections of 3.3 kΩ and 1 μF, corner 48 Hz: alias 0.77 LSB, droop 0.09 dB, delay 6.6 ms**, with a 14 % margin on the corner for component tolerance.

What to improve: a second-order active filter (Extension 1) for a sharper knee; oversampling with digital decimation, which moves most of the filtering into software where it is exact; and a band specification for the vibration rather than one tone. What not to improve: the ADC resolution. A 16-bit converter has 16 times smaller LSBs, so the same 50 mV pickup is 993 LSB and the filter must attenuate 16 times more; more bits make aliasing worse, not better.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2)–(5) and (7); `rc_filter.png`; the four experiment tables, each SYNTHETIC one labelled; and a half-page design note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) for an anti-aliasing filter on your own channel: signal band, sample rate, ADC, interference, the chosen circuit and values, attenuation, droop, delay and the margin for component tolerance. If you did the optional measurement, add the raw CSV and the fitted τ against the prediction.

### Adversarial self-checks

- **What if the vibration amplitude doubles?** The alias doubles to 1.5 LSB with the chosen design; (7) with a = 124 gives a largest corner of 38.9 Hz for two sections, so C = 1.5 μF (f_c = 32.2 Hz) restores the margin (0.69 LSB) at a cost of 9.9 ms delay instead of 6.6.
- **What if the vibration frequency falls to 105 Hz?** Two sections at 48 Hz attenuate it only by 5.7 (|H| = 0.174), leaving 10.8 LSB. Rule (7) says three sections at a 27.4 Hz corner reach 62 with 0.43 dB of droop, but the delay grows to 3RC = 17.4 ms, and as the interference approaches the signal band the delay and the sensitivity to component values grow without bound. The better remedy is a faster ADC (sample at 1 kHz and decimate digitally) or a sharper active filter.
- **Which assumption fails first in the real world?** The single tone: an engine or propeller produces harmonics and a broadband floor; the next is the capacitor value (±20 % for many parts), then the op-amp's output swing near a 3.3 V rail.
- **What measurement could falsify the design?** Inject a known 430 Hz, 50 mV sine at the filter input on the bench and look for a 30 Hz tone in the logged counts; more than 1 LSB means the filter or its wiring is wrong; or record the counts with the aircraft's engine running and look at the spectrum.
- **Does the result survive a change of tool?** `freqs`, `step`, `lsim` and the closed forms agree to 10⁻⁹ or better; the time-domain SYNTHETIC check agrees with the frequency-domain prediction to within one LSB.
- **Which decision changes if the component tolerance doubles?** With ±40 % on C, the corner can rise to 80 Hz and the 430 Hz attenuation falls to 1/29: 2.1 LSB, failing. A tighter capacitor (C0G/NP0 ceramic, ±5 %) or a third section restores the margin.

## Check yourself

1. Write the impedances of R, L and C. Which one is an open circuit at DC and which a short circuit?
2. A series RLC has L = 10 mH and C = 1 μF. What is f₀, and what does R change?
3. A 400 Hz load draws 10 A at 115 V rms with power factor 0.9. What real power does it consume, and why does the wiring see 10 A?
4. State the corner frequency of an RC low-pass, its gain and phase there, and the roll-off above it.
5. Why does connecting a second RC section directly to the first change the response, and what fixes it?
6. Give the two ideal-op-amp rules and derive the gain of the non-inverting amplifier from them.
7. A 430 Hz tone is sampled at 100 Hz. At what frequency does it appear, and why can software not remove it afterwards?
8. A 12-bit ADC has v_ref = 3.3 V. What is one LSB, what count does 1.000 V produce, and what is the RMS quantisation error?
