# P04.10  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown. Every simulated signal here is SYNTHETIC.

## 1. A load on the 400 Hz bus (s = 15: R = 11 Ω)

**Assumptions.** Sinusoidal steady state, ideal components, RMS phasors.

(a) X_L = 2π × 400 × 1.5 × 10⁻³ = **3.770 Ω**; |Z| = √(11² + 3.770²) = **11.628 Ω**; phase **18.9°**; I = 115/11.628 = **9.890 A rms**; **P = I²R = 1076 W**, **Q = I²X = 369 var**, **pf = cos 18.9° = 0.946**.

(b) A parallel capacitor draws −ωCV² vars; cancelling Q needs **C = Q/(ωV²) = 369/(2π·400·115²) = 11.09 μF**. The combined impedance is 12.29 Ω, purely resistive; the line current falls to **9.356 A** while the load still takes **1076 W** (the capacitor dissipates nothing).

(c) The generator, the wires and the circuit breakers are sized by the current, i.e. by |S| = VI, not by P; a poor power factor makes them carry current that does no work.

**Check.** P = V I pf = 115 × 9.890 × 0.946 = 1076 W.

**Common mistakes.** Using peak instead of RMS; adding the capacitor in series; expecting the correction to change P.

## 2. The RC filter, the loaded ladder and the op-amp

(a) H = Z_C/(R + Z_C) = 1/(1 + R/Z_C) = **1/(1 + jωRC)**. |H|² = 1/(1 + (ωRC)²) equals ½ at ωRC = 1, i.e. **f_c = 1/(2πRC)**; there arg H = −arctan 1 = **−45°**. For ωRC ≫ 1, |H| ≈ 1/(ωRC): a factor 10 in ω is a factor 10 in |H|, **−20 dB/decade**. φ(ω) = −arctan(ωRC) ≈ −ωRC for small ω, so **−dφ/dω = RC**.

(b) With admittances: node 1, (v₁ − v_in)/R + jωCv₁ + (v₁ − v₂)/R = 0; node 2, (v₂ − v₁)/R + jωCv₂ = 0, so v₁ = v₂(1 + x), x = jωRC. Substituting in node 1 (times R): v₁(2 + x) − v₂ = v_in → v₂[(1 + x)(2 + x) − 1] = v_in → **H = 1/(1 + 3x + x²)**. For −3 dB, x = jy: |1 − y² + 3jy|² = 2 → y⁴ + 7y² − 1 = 0 → y² = (√53 − 7)/2, **y = 0.3742**, f = 0.3742/(2πRC), against y = √(√2 − 1) = 0.6436 for the buffered pair. For R = 3.3 kΩ, C = 1 μF: 18.05 and 31.04 Hz, as the code's sweep finds.

(c) Rule 1: no current into v₋, so R_f and R_g carry the same current and v₋ = v_out R_g/(R_f + R_g) = βv_out. Rule 2: v₋ = v₊ = v_in. Hence **G = 1/β = 1 + R_f/R_g**. With finite A: v_out = A(v_in − βv_out) → **G = A/(1 + Aβ)**. Then 1/β − G = 1/(β(1 + Aβ)) ≈ 1/(Aβ²) = **G²/A**: for G = 3.7 and A = 10⁵, 1.37 × 10⁻⁴, the code's 3.70 − 3.69986.

**Check.** (4) at x → 0 is 1 and at large x is 1/x², the same −40 dB/decade as the buffered pair: loading changes the knee, not the asymptote.

**Common mistakes.** Treating the second section as if the first were unloaded; taking the −3 dB point of a pair at f_c (it is at 0.64 f_c buffered).

## 3. An LC filter on a 28 V bus (s = 15: C = 25 μF)

(a) **f₀ = 1/(2π√(100 μH × 25 μF)) = 3183 Hz**; √(L/C) = 2.0 Ω, so **Q = 2.0/0.5 = 4.0**; at f₀ the reactances cancel, the current is V_in/R and the voltage across C is |Z_C|/R = **Q = 4.0** times the input. Noise at 3.2 kHz arrives at the payload four times larger: an undamped filter is a resonator.

(b) Sweep: **peak |H| = 4.03 at 3133 Hz** (slightly below f₀ because the output is taken across C); **|H(100 kHz)| = 1.01 × 10⁻³, −59.9 dB**; **−40 dB/decade** above f₀ (two reactive elements).

(c) R = 0.05 Ω: Q = 40, peak 40.0; R = 2 Ω: Q = 1, peak 1.155. **Q = 1 at R = √(L/C) = 2 Ω.** Designers add a damping resistor (usually in series with a large capacitor in parallel with C, so it dissipates nothing at DC) because an undamped input filter peaks, and interacting with the negative input resistance of a switching regulator it can oscillate.

**Check.** At f₀, |Z_C| = 1/(2π × 3183 × 25 × 10⁻⁶) = 2.0 Ω = √(L/C).

**Common mistakes.** Forgetting the winding resistance (Q becomes infinite); reading the peak as attenuation.

## 4. Your own anti-aliasing filter (s = 15: f_vib = 910 Hz)

(a) 910 Hz sampled at 200 Hz appears at **90 Hz**. The 10-bit LSB is 3.2227 mV, so 20 mV is **6.21 LSB**: attenuation ≥ 6.21 (15.9 dB).

(b) Largest corners and droop at 10 Hz: **1 section 148.57 Hz, 0.020 dB; 2 sections 398.83 Hz, 0.005 dB; 3 sections 590.2 Hz, 0.004 dB.** One section suffices.

(c) R = 1.5 kΩ, C = 1 μF: **f_c = 106.1 Hz**, |H(910)| = 0.1158, residual **0.72 LSB**, droop **0.038 dB**, delay **1.5 ms**.

(d) SYNTHETIC check (2 Hz, 0.5 V signal plus 20 mV at 910 Hz, 40 kHz simulation grid): alias tone in the counts **6.24 LSB without the filter, 0.78 LSB with it**, against 6.21 and 0.72 predicted. The corner is above the 100 Hz Nyquist frequency: anything between 100 Hz and a few hundred hertz passes almost unattenuated and aliases into 0–100 Hz. The design handles the one tone it was given and nothing else; a band specification (or a corner below Nyquist) is needed if broadband noise is present.

**Common mistakes.** Using the 12-bit LSB; designing for the alias frequency (90 Hz) instead of the real frequency (910 Hz) the analog filter sees.

## 5. What averaging can and cannot do (s = 15: v₀ = 1.0015 V)

(a) v₀ = 1243.07 LSB, so the ADC returns **1243** every time; the converted value 1243.5 LSB is **+0.426 LSB** high, and the average of N identical readings has the same error for every N: averaging cannot recover information that quantisation removed.

(b) Bias and RMS error of the average, in LSB, 400 SYNTHETIC repeats, seed 15:

| σ | N = 1 | N = 16 | N = 256 |
|---|---|---|---|
| 0.3 LSB | +0.041, 0.493 | +0.031, 0.124 | +0.025, 0.040 |
| 1 LSB | +0.024, 1.032 | +0.009, 0.262 | −0.001, 0.064 |
| 2 LSB | −0.024, 1.960 | +0.011, 0.519 | +0.012, 0.123 |

(c) Noise of about one LSB or more **dithers** the input across neighbouring codes, so the fraction of readings above and below encodes where v₀ lies between them and the average converges on v₀: the bias falls from 0.43 to a few hundredths of an LSB. The RMS error then falls as 1/√N (1.03 → 0.26 → 0.064 at σ = 1, factors of 4 for 16 times more samples), from √(σ² + LSB²/12). At σ = 0.3 LSB the dither is too weak to remove the bias fully (+0.025 at N = 256). The price is time: 256 samples at 100 Hz take 2.56 s, which is a 256-fold reduction in bandwidth.

**Check.** σ = 1, N = 256: √(1 + 1/12)/16 = 0.065 LSB, against 0.064.

**Common mistakes.** Expecting averaging alone to add resolution; forgetting that averaging is a low-pass filter with its own delay.

## 6. Choosing the signal chain

(a) No. After sampling, the 430 Hz vibration *is* a 30 Hz tone in the data; a notch would also remove any true 30 Hz content and would miss the alias as soon as the propeller speed changes and the alias moves. Aliases must be removed before the ADC.

(b) More bits make the problem worse with the same filter: at 16 bits one LSB is 50 μV, the 50 mV pickup is 993 LSB, and the two-section design leaves **12.3 LSB** of alias. Either the filter is redesigned for an attenuation of about 1000, or the extra bits are noise.

(c) The droop at 5 Hz grows from 0.093 to 0.316 dB (still within 0.5 dB) and the −3 dB point drops from 31 to 18 Hz; the 430 Hz attenuation is nearly unchanged. Acceptable for this specification if the extra phase lag is acceptable, but the ADC's sample-and-hold now loads a 3.3 kΩ source directly, which may not settle in the acquisition time; check the datasheet before deleting the buffer.

(d) 5 mV is **6.2 LSB** of the 12-bit, 3.3 V converter: about 2.6 bits of resolution. Add a non-inverting amplifier of gain about 600 (for example R_f = 600 kΩ, R_g = 1 kΩ: G = 601, 5 mV → 3.005 V), with a low-offset op-amp since its input offset is amplified 601 times too.

(e) Two sections at 48.2 Hz lag by 2 arctan(8/48.2) = **18.8°** at 8 Hz, about 6.5 ms; that eats roughly a third of a typical 45–60° phase margin. Either the filter corner is raised (a sharper active filter with less passband phase) or the controller is designed with the filter in the loop.

## Check-yourself answers

1. Z_R = R, Z_L = jωL, Z_C = 1/(jωC). The capacitor is an open circuit at DC, the inductor a short circuit.
2. f₀ = 1/(2π√(10 mH × 1 μF)) = 1591.5 Hz; R sets the peak current V/R and the bandwidth f₀/Q = R/(2πL), not f₀.
3. P = 115 × 10 × 0.9 = 1035 W; the wiring carries the full 10 A because the current includes a reactive part that flows back and forth without delivering energy.
4. f_c = 1/(2πRC); −3.01 dB and −45°; −20 dB per decade per section.
5. The second section draws current from the first, so the first's output is no longer its unloaded value (H becomes 1/(1 + 3x + x²)); a unity-gain buffer between them, which draws no current and drives the second section from a low impedance, restores 1/(1 + x)².
6. No current into the inputs; with negative feedback v₋ = v₊. Then v_in = v₋ = v_out R_g/(R_f + R_g), G = 1 + R_f/R_g.
7. At |430 − 4 × 100| = 30 Hz. The samples of the 430 Hz and 30 Hz sinusoids are identical, so nothing in the data distinguishes them.
8. 3.3/4096 = 0.8057 mV; count 1241; LSB/√12 = 0.2326 mV.
