# P04.10  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. R = 8 + s/5; X_L = 2π·400·1.5 × 10⁻³; Z = R + jX_L; |Z|, arg Z. 2. I = 115/|Z|; `ac_power(115, Z)` → P, Q, pf. 3. C = Q/(ωV²); new total Z_t = Z ∥ Z_C (`parallel`); line current 115/|Z_t|; P unchanged (the capacitor dissipates nothing). 4. Line current, hence wire mass and generator rating, scales with |S|, not P.

## 2.
(a) H = 1/(1 + R/Z_C) = 1/(1 + jωRC); |H|² = 1/(1 + (ωRC)²); phase −arctan ωRC; for ωRC ≫ 1, |H| ≈ 1/(ωRC); −dφ/dω at ω → 0 equals RC. (b) Node 1: (v₁ − v_in)/R + jωCv₁ + (v₁ − v₂)/R = 0; node 2: (v₂ − v₁)/R + jωCv₂ = 0 → v₁ = v₂(1 + x); substitute; collect. −3 dB: with x = jy, |1 − y² + 3jy|² = 2 → y⁴ + 7y² − 1 = 0. (c) v₋ = v_out R_g/(R_f + R_g) = v_in → G; with v_out = A(v_in − βv_out) solve for G; 1/β − A/(1 + Aβ) = 1/(β(1 + Aβ)) ≈ 1/(Aβ²).

## 3.
f₀ = 1/(2π√LC); Q = √(L/C)/R; |H(f₀)| = |Z_C|/R = Q. Sweep `ff = np.logspace(2, 5, 200001)`, H = Z_C/(R + Z_L + Z_C); `argmax`; `gain_db(H(1e5))`. Repeat for R = 0.05 and 2 Ω. Q = 1 at R = √(L/C).

## 4.
f_vib = 610 + 20s; alias = `aliased_frequency(f_vib, 200)`; a = 0.02/`adc_lsb(3.3, 10)`. `required_stages_and_corner(10, f_vib, a, 0.5)`; for n = 1, 2, 3 evaluate (7) and the droop. Pick E12 R and C with f_c ≤ the limit; residual 0.02|H(f_vib)|/LSB; `group_delay_dc`. `t = np.arange(0, 3, 1/40000)`; `v = synthetic_sensor_signal(t, 2.0, 0.5, f_vib=f_vib, a_vib=0.02)`; `filter_and_sample(t, v, R, C, n, 200, 3.3, 10)` for n = 0 and your n; `tone_amplitude(counts[ts >= 1], 200, alias)`.

## 5.
v₀ = 1 + s × 10⁻⁴ V; `adc_convert(v0)`; error (counts + ½)LSB − v₀. Loop σ ∈ {0.3, 1, 2} LSB, N ∈ {1, 16, 256}: 400 repeats of the mean of `counts_to_volts(adc_convert(v0 + rng.normal(0, σ·LSB, N)))` − v₀; bias = mean, RMS = √(mean square), both in LSB. Compare the RMS with σ_q/√N where σ_q² ≈ σ² + LSB²/12.

## 6.
(a) Aliasing is irreversible: a notch at 30 Hz also removes any true 30 Hz content, and the alias moves if the engine speed changes. (b) Count the interference in LSB at 16 bits. (c) Compare 0.316 and 0.093 dB droop and the corner. (d) 5 mV in LSB of the 12-bit converter; add a non-inverting gain stage (5) to fill the range. (e) Compare the filter's phase lag at 8 Hz with a typical phase margin.
