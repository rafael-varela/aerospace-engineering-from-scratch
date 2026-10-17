# P09.06  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown. Every gyro and accelerometer record here is SYNTHETIC.

## 1. Where the propeller lines go (s = 15: f₁ = 385 Hz)

**Assumptions.** Ideal sampling at 400 Hz; the lines are pure sinusoids.

(a) 385 Hz is 15 Hz below 400: **alias 15 Hz**. 770 Hz is 30 Hz below 800: **30 Hz**. 1155 Hz is 45 Hz below 1200: **45 Hz**.

(b) The fundamental (15 Hz) and the second harmonic (30 Hz, right at the band edge) fall in the 0–30 Hz band; the third harmonic lands just outside, at 45 Hz.

(c) `dominant_frequency` returns **15.00, 30.00 and 45.00 Hz**.

**Check.** The samples of cos(2π·385n/400) and cos(2π·15n/400) are identical because 385 = 400 − 15 and cos is even.

**Common mistakes.** Subtracting f_s only once (1155 − 400 = 755 Hz is still above Nyquist); forgetting that harmonics alias too.

## 2. The Butterworth prototype and its delay

(a) With s = jω, ω² = −s², so H(s)H(−s) = 1/(1 + (−s²/ω_c²)ⁿ). Its poles solve (−1)ⁿ(s/ω_c)²ⁿ = −1, i.e. (s/ω_c)²ⁿ = e^{jπ(n+1)}. Hence \|s\| = ω_c and s = ω_c e^{jθ_m}, θ_m = π(n + 1 + 2m)/(2n), m = 0, …, 2n − 1: **2n points spaced π/n apart on the circle**. Those with cos θ < 0 lie in the left half plane: the angles strictly between π/2 and 3π/2. These are m = 0, …, n − 1, and writing k = m + 1 gives θ_k = π(2k + n − 1)/(2n), k = 1..n (for k = 1, θ = π(n + 1)/(2n), just past π/2; for k = n, θ = π(3n − 1)/(2n), just short of 3π/2). H(s) takes these n so that it is stable; H(−s) the mirror images.

(b) From (4): 2n log(f_p/f_c) ≤ log ε_p gives log f_c ≥ log f_p − log(ε_p)/(2n); 2n log(f_st/f_c) ≥ log ε_s gives log f_c ≤ log f_st − log(ε_s)/(2n). A corner exists when log f_p − log(ε_p)/(2n) ≤ log f_st − log(ε_s)/(2n), i.e. **n ≥ log(ε_s/ε_p)/(2 log(f_st/f_p))**, and then **f_p ε_p^{−1/2n} ≤ f_c ≤ f_st ε_s^{−1/2n}**.

(c) N(jω) = b₀(1 + jω b₁/b₀ + O(ω²)) and similarly for D, so log H(jω) = log(b₀/a₀) + jω(b₁/b₀ − a₁/a₀) + O(ω²). The phase is the imaginary part, so φ ≈ −ω(a₁/a₀ − b₁/b₀) and **τ(0) = −dφ/dω = a₁/a₀ − b₁/b₀**. For the Butterworth, b is a constant (b₁ = 0) and a(s) = Π(s − p_k): a₀ = Π(−p_k), a₁ = Σ_j Π_{k≠j}(−p_k) = a₀ Σ_k 1/(−p_k). So τ(0) = Σ_k(−1/p_k) = Σ_k −cos θ_k/ω_c (the imaginary parts cancel in conjugate pairs).
- n = 1: p = −ω_c, **τ = 1/ω_c**.
- n = 2: p = ω_c e^{±j3π/4}, −cos(3π/4) = 1/√2 twice, **τ = √2/ω_c**.
- n = 3: p = −ω_c and ω_c e^{±j2π/3}: 1 + ½ + ½, **τ = 2/ω_c**.

These equal 1/(ω_c sin(π/2n)) = 1/ω_c, 1/(ω_c sin 45°) = √2/ω_c, 1/(ω_c sin 30°) = 2/ω_c. The code gives τω_c = 1.0, 1.4142, 2.0, 2.6131, 3.2361 for n = 1–5, matching 1/sin(π/2n).

**Check.** For the lesson's design (n = 3, ω_c = 314.3 rad/s): 2/314.3 = 6.36 ms, as `group_delay_dc` prints.

**Common mistakes.** Taking all 2n poles; forgetting that the phase of a constant numerator is zero; confusing group delay with phase delay (they coincide only at DC).

## 3. A launch-vehicle accelerometer

**Assumptions.** Butterworth prototype; the vibration can be anywhere above the band; the 60 dB requirement applies at f_s − f_band.

(a) f_s = 1000 Hz: stop edge **900 Hz**; (5) gives **n = 4**; corner **118.4–160.0 Hz**; at the geometric mean 137.7 Hz the DC delay is **3.02 ms** and the phase at 20 Hz is **−21.82°**.

(b) f_s = 4000 Hz: stop edge **3900 Hz**; **n = 3**; corner **125.3–390.0 Hz**; at the geometric mean 221.0 Hz the delay is **1.44 ms** and the phase at 20 Hz **−10.38°**.

(c) The 4 kHz rate halves the phase error at 20 Hz with a lower order and a wider corner range (more tolerance to component values). The cost is four times the data rate and storage, and, if the analysts want 1 kHz data, a digital decimation filter (`P09.07`), whose delay must also be counted. If phase accuracy matters, recommend 4 kHz.

**Check.** Ratio f_st/f_p = 9 in (a) and 39 in (b); log₁₀(ε_s/ε_p) = log₁₀((10⁶ − 1)/(10^{0.1} − 1)) = 6.587; 6.587/(2 log₁₀ 9) = 3.45 → 4, and 6.587/(2 log₁₀ 39) = 2.07 → 3.

**Common mistakes.** Putting the stop edge at f_s/2 (500 Hz), which over-designs; quoting phase without saying it depends on the corner chosen.

## 4. Your SYNTHETIC gyro record (s = 15)

The design range from `antialias_spec`: order 3, f_c ∈ [42.60, 58.73] Hz, geometric mean 50.02 Hz.

(a) SYNTHETIC record, seed 15 (not a measurement):

| f_c | predicted 10 Hz alias | measured | DC delay | phase at 5 Hz |
|---|---|---|---|---|
| no filter | 5.0000 deg/s | 5.0032 | 0 | 0° |
| 42.60 Hz | 0.0065 | 0.0065 | 7.47 ms | −13.48° |
| 50.02 Hz | 0.0105 | 0.0106 | 6.36 ms | −11.47° |
| 58.73 Hz | 0.0171 | 0.0171 | 5.42 ms | −9.77° |

All three corners meet the 0.02 deg/s requirement; the trade is 2 ms of delay for a factor 2.6 in alias.

(b) The 780 Hz harmonic (2 deg/s) appears at 20 Hz: **2.0085 deg/s** without a filter. With the filter the predictions are 0.00033, 0.00053 and 0.00085 deg/s; the measurements are 0.00060, 0.00047 and 0.00042 deg/s.

(c) (a) agrees within 3 %. (b) does not, and should not: the predicted harmonic is at or below the white-noise level in one DFT bin. Running the filtered noise alone through the same sampler for 20 seeds puts 0.00062 deg/s (mean; up to 0.0017) into the 20 Hz bin. The measurement is noise plus a tiny line, so it does not decrease with f_c as the prediction does. The right statement is "below the noise floor", not a number.

**Common mistakes.** Predicting with \|H\| at the alias (10 Hz) instead of the true frequency (390 Hz); reading a noise-floor measurement as confirmation.

## 5. Reconstruction from a finite record

Maximum error over the middle half of the record (f_s = 50 Hz):

| N | 5 Hz | 12.5 Hz | 22.5 Hz |
|---|---|---|---|
| 50 | 0.0183 | 0.0117 | 0.0399 |
| 200 | 0.0042 | 0.0051 | 0.0123 |
| 800 | 0.0011 | 0.0013 | 0.0029 |

(b) Formula (7) assumes infinitely many samples. The samples outside the record are treated as zero, and each one's sinc tail decays only as 1/(distance in samples); the error at a point is roughly the sum of the missing tails, which falls about as 1/(distance to the record's ends), i.e. as 1/N in the middle. At 22.5 Hz = 0.45 f_s the samples alternate in sign almost every step; the missing tails, which alternate too, add rather than cancel, and the error is largest. Near Nyquist, reconstruction needs long records.

(c) Upsampling by sinc interpolation needs samples from the future (the sinc is non-causal) and a long window for accuracy near 200 Hz; in a real-time loop it adds delay and cannot recover anything the 400 Hz sampler did not capture. It is useful for display or analysis, not for a faster inner loop.

**Check.** At every sample instant the reconstruction equals the sample to 10⁻¹², as the test asserts.

**Common mistakes.** Measuring the error at the ends of the record; expecting reconstruction to recover frequencies above f_s/2.

## 6. Is the design good enough?

(a) Over 351–429 Hz the alias runs from 49 Hz (351 Hz) down to 0 Hz (400 Hz, which appears as a constant bias) and back up to 29 Hz (429 Hz). The order-3, 50.02 Hz filter leaves **0.0145 deg/s at 351 Hz**, 0.0124 at 370 Hz, 0.0105 at 390 Hz, 0.0091 at 410 Hz and 0.0079 at 429 Hz: every case is below 0.02 deg/s. The design is robust because it is specified at the stop edge 370 Hz and the filter is monotonic; only the 351–370 Hz part of the range sits below the stop edge, and its aliases (30–49 Hz) fall outside the band anyway.

(b) First, after sampling, the alias and a real 10 Hz motion have identical samples, so a notch at 10 Hz also removes real 10 Hz motion (a structural mode, a gust response). Second, the alias moves with rotor speed (part (a)), including to 0 Hz, where a notch would remove the true low-frequency rate, and a fixed notch misses it elsewhere.

(c) The phase is almost exactly a delay: −2.29°, −4.58°, −11.47° and −23.07° at 1, 2, 5 and 10 Hz are 6.36, 6.37, 6.37 and 6.41 ms. The delay is 2/ω_c, so it is set by the corner choice and the order together. With only 35° of phase margin, 11.5° is a third of it; move the corner to the top of the range (58.7 Hz: −9.77°), or oversample (part d), or include the filter in the loop design and add phase lead (`P09.12`).

(d) At 2 kHz the stop edge is 1970 Hz, **order 2** suffices, the corner range is 50.8–124.6 Hz, and at the geometric mean 79.5 Hz the delay is **2.83 ms** (−5.10° at 5 Hz). The new work: five times the samples, and a digital decimation filter to bring the rate down to what the controller runs at, whose own delay (`P09.07`) must be counted.

**Common mistakes.** Designing at the nominal line only; treating filter phase as unrelated to delay.

## Check-yourself answers

1. At \|1190 − 3 × 400\| = **10 Hz**. Its samples are identical to those of a real 10 Hz signal, so any filter after the sampler that removes it removes real 10 Hz content too.
2. Frequencies between f_s/2 and f_s − f_band alias to between f_band and f_s/2: outside the band, harmless. The first frequency that folds *into* the band is f_s − f_band; that is where attenuation is needed.
3. \|H\|² = 1/(1 + (ω/ω_c)²ⁿ); −3.01 dB at ω_c for every n; −20n dB per decade far above.
4. On the unit circle at angles 120°, 180° and 240°: −1 and −½ ± j(√3/2).
5. Its delay, about 1/(ω_c sin(π/2n)), grows with n at a fixed corner, and that delay is phase lag inside the loop at the crossover frequency.
6. x(t) = Σ x[n] sinc(f_s t − n). It returns a 10 Hz cosine, the band-limited signal with those samples.
7. 80 dB: 20n log₁₀ 10 = 80 for n = 4 (the code gives −80.00 dB at 10ω_c).
