# P09.06  Sampling & Aliasing

> A sampler cannot tell f from f − k f_s. Whatever arrives above the Nyquist frequency is folded into the band you care about, and nothing after the sampler can unfold it. The filter has to come first.

| | |
|---|---|
| **Learning outcomes** | 1. Compute where any frequency lands after sampling, f_a = \|f − f_s round(f/f_s)\|, and confirm it from the FFT of the samples. 2. Derive the analog Butterworth low-pass from \|H\|² = 1/(1 + (ω/ω_c)²ⁿ): its poles on a circle, its −3 dB corner and its 20n dB/decade roll-off; implement it from its poles and verify it against `scipy.signal.butter(analog=True)`. 3. Choose the order and the corner of an anti-aliasing filter from a passband droop and a stopband attenuation, with the stopband edge at f_s − f_band, and state the delay it costs the control loop (τ(0) = a₁/a₀). 4. Simulate the analog filter on a SYNTHETIC gyro record, sample it, and show that the measured alias equals the predicted a_vib\|H(f_vib)\|. 5. Reconstruct a band-limited signal from its samples by sinc interpolation and explain why an undersampled tone comes back as its alias. |
| **Prerequisites** | `P09.05` (frequency response, Bode magnitude and phase), `P09.02` (impulse response), `P09.01` (signals), `P04.10` (RC filters, the ADC, the first look at aliasing), `P02.10` (the Fourier transform, the FFT), `P02.07` (complex numbers and phasors), `P02.05` (the matrix exponential and the zero-order-hold step). Skill check: you can find the roots of s² + √2 s + 1 by hand and read a magnitude in dB. |
| **Engineering question** | A MEMS rate gyro on a small UAV is sampled at f_s = 400 Hz; the attitude loops need 0–30 Hz with at most 0.5 dB of loss, and the airframe shakes at a 390 Hz propeller blade-pass line of 5 deg/s. Where does that line appear in the samples, what analog Butterworth filter in front of the sampler keeps every alias below 0.02 deg/s, and what does that filter cost the attitude loop in delay? Answer: at 10 Hz, inside the band; order 3 with a corner between 42.6 and 58.7 Hz (50.0 Hz chosen); 6.4 ms of delay, 11.5° of phase at 5 Hz. |
| **Assumptions and limits** | Ideal, instantaneous sampling at a constant rate; linear time-invariant analog filter with exact component values; vibration as discrete sinusoidal lines plus white noise (real vibration also has broadband content and lines that move with rotor speed); the filter is simulated on a 20 kHz grid with the input held over each 50 μs step, which lowers the 390 Hz line by 0.06 % in amplitude (the hold's sinc(f Δt) factor). Quantisation, gyro bandwidth limits and the gyro's own internal filtering are ignored. Every gyro record is **SYNTHETIC** (seeded, formula stated in `synthetic_imu_record`). |
| **Analytical work** | The folding rule from the identity cos(2π f n/f_s) = cos(2π(f − k f_s) n/f_s); the Nyquist frequency and the stopband edge f_s − f_band; Butterworth poles from 1 + (s/jω_c)²ⁿ = 0; the order inequality and the allowed corner range; the DC group delay τ(0) = a₁/a₀ − b₁/b₀ and its Butterworth value 1/(ω_c sin(π/2n)); the sinc interpolation formula. |
| **Computational deliverable** | `starter/main.py`: `aliased_frequency`, `dominant_frequency`, `sinc_reconstruct`, `butterworth_poles`, `butterworth_analog`, `freq_response`, `butterworth_order`, `antialias_spec`, `group_delay_dc`, `simulate_analog` (`poly_from_roots`, `tf_to_companion`, `expm_series`, `discretize_zoh`, `tone_amplitude`, `synthetic_imu_record`, `sample_through_filter`, `imu_experiment` and `main` are given). Tests prove: aliases fold into [0, f_s/2] and give identical samples; the FFT of a sampled tone peaks at the predicted alias; the prototype equals `scipy.signal.butter(analog=True)` for orders 1–6; \|H\|² = 1/(1 + (ω/ω_c)²ⁿ), −3.01 dB at ω_c, 20n dB/decade; the order rule meets both specifications and n − 1 cannot; the delay formula equals the phase slope; the simulator equals `scipy.signal.lsim` and is linear; sinc reconstruction interpolates the samples and returns the alias of an undersampled tone; in the SYNTHETIC record the measured alias equals the prediction within 3 %. |
| **Practical activity** | **computational experiment** on a SYNTHETIC gyro record: run the record through Butterworth filters of order 0–4, sample at 400 Hz, measure the 10 Hz alias and the delay, and sweep the corner across the allowed range. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `sampling_aliasing.png`, the derivations of the folding rule, the Butterworth poles and the order rule, the experiment table (labelled SYNTHETIC), and a half-page design note giving the order, corner, attenuation, droop and delay of your anti-aliasing filter. |
| **Assessment** | `assessments/phase-09/grading-rubric.md`; competencies 1–4. Pass: the folding rule and the pole formula are derived; all public tests pass; the note shows that order 2 fails and order 3 passes, with numbers; the SYNTHETIC check is labelled as such and its agreement with the prediction explained. |
| **References** | `@mit-ocw-6-003-signals-and-systems` Lecture 21 "Sampling" and Lecture 22 "Sampling and Quantization" (required-free: the sampling theorem, aliasing, reconstruction); Lecture 17 "Discrete-Time (DT) Frequency Representations" (supplementary-free). `@scipy-userguide` "Signal Processing (scipy.signal)" (supplementary-free: `butter`, `freqs`, `lsim`, used here only as oracles). The Butterworth design rules and the anti-aliasing order selection are original instruction (no free source mapped). Map `references/phase-09.md` not yet available when written; only ids already in `sources.yaml` are cited, with lecture titles as recorded in `references/phase-02.md`. |
| **Extensions** | ★ Replace the Butterworth with a Bessel prototype (maximally flat delay) of the same order and compare the attenuation at 370 Hz and the phase at 5 Hz. ★ Let the blade-pass line sweep from 351 to 429 Hz (± 10 % rotor speed) and plot where its alias moves and how large it is. ★ Add quantisation (a 16-bit gyro output) and find the order at which the alias falls below one LSB. |
| **Time** | unmeasured |

## 1. The engineering problem

A flight computer does not see the gyro's voltage or its digital output continuously; it reads one number every 1/f_s seconds. For a small UAV with f_s = 400 Hz that is every 2.5 ms, and it is all the attitude controller ever knows about the motion. The question this lesson answers is what happens to the parts of the signal the sampler is too slow for. The answer is the uncomfortable one already met in `P04.10`: they do not disappear. A 390 Hz vibration line from the propellers produces exactly the same 400 samples per second as a 10 Hz rotation would. The controller sees a 10 Hz oscillation of 5 deg/s, well inside the 0–30 Hz band it acts on, and drives the motors to fight a motion that is not there. No digital filter after the sampler can help, because the samples of the alias and of a real 10 Hz motion are identical.

The only defence is an analog low-pass filter before the sampler, the **anti-aliasing filter**, and it is a trade. A sharp filter removes the vibration but delays the signal the controller needs, and every millisecond of delay costs phase margin (`P09.10`). A gentle filter keeps the delay small and lets the alias through. This lesson turns the trade into numbers: it finds where every frequency lands after sampling, builds the Butterworth filter (the standard "maximally flat" analog prototype) from its poles, picks the smallest order that meets a passband and a stopband specification, and checks the design by running a SYNTHETIC gyro record through a simulated filter and a simulated sampler and measuring what comes out.

Success criteria, written first: the FFT of a sampled 390 Hz tone peaks at the predicted 10 Hz; the Butterworth prototype built from its poles matches `scipy.signal.butter(analog=True)` to 10⁻¹⁰; the design keeps the passband loss at 30 Hz below 0.5 dB and the attenuation at 370 Hz above 47.96 dB; in the SYNTHETIC record the 10 Hz alias falls from about 5 deg/s to below 0.02 deg/s and agrees with the prediction a_vib\|H(390 Hz)\| within 3 %.

## 2. Understand and derive

### 2.1 What a sampler cannot see

Sample x(t) = cos(2πft + φ) at t = n/f_s: x[n] = cos(2π(f/f_s)n + φ). Replace f by f − k f_s for any integer k. The argument changes by −2πkn, a whole number of turns, so **the samples are identical**. Replace f by k f_s − f and use cos(−θ) = cos θ: the samples are those of a cosine at f with phase −φ, again identical in frequency. Every frequency therefore has a twin in the first **Nyquist zone** [0, f_s/2], the **alias**

**f_a = \|f − f_s round(f/f_s)\|.**   (1)

The **Nyquist frequency** f_s/2 is the highest frequency a sampler represents without ambiguity (the sampling theorem of `P02.10`). For the gyro, f_s = 400 Hz and f_s/2 = 200 Hz:

| true frequency | 10 | 190 | 210 | 390 | 410 | 780 | 1190 Hz |
|---|---|---|---|---|---|---|---|
| alias at 400 Hz | 10 | 190 | 190 | **10** | **10** | **20** | **10** Hz |

The 390 Hz line and its second harmonic at 780 Hz both land inside 0–30 Hz. The table is reproduced by `main()`; the FFT of 400 samples of a 390 Hz cosine peaks at 10.00 Hz, which is the independent check: the spectrum of the samples does not know about formula (1).

**Which frequencies threaten the band?** A frequency f aliases into [0, f_band] exactly when it lies within f_band of a multiple of f_s. The lowest such frequency above the band is **f_s − f_band** = 370 Hz. That is the **stopband edge** of the anti-aliasing filter: everything at or above 370 Hz must be attenuated enough, and a monotonic filter that is attenuated enough at 370 Hz is attenuated at least as much everywhere above it. Frequencies between f_band and f_s − f_band (here 30–370 Hz) alias to between 30 and 200 Hz: out of the band, harmless to the controller if it ignores that range.

### 2.2 The Butterworth prototype

Ask for a low-pass whose magnitude is as flat as possible at ω = 0 and falls monotonically:

**\|H(jω)\|² = 1 / (1 + (ω/ω_c)²ⁿ).**   (2)

At the **corner** ω_c, \|H\|² = ½: −3.01 dB for every order n. Far above it \|H\| ≈ (ω_c/ω)ⁿ, a slope of **−20n dB per decade**. All 2n − 1 derivatives of \|H\|² vanish at ω = 0, which is what "maximally flat" means.

To find H(s), write \|H(jω)\|² = H(s)H(−s) at s = jω, so ω² = −s² and H(s)H(−s) = 1/(1 + (−s²/ω_c²)ⁿ). Its poles solve (s/(jω_c))²ⁿ = −1: 2n points equally spaced on the circle \|s\| = ω_c. H(s) takes the n in the left half plane (so the filter is stable), H(−s) the mirror images:

**p_k = ω_c exp(jπ(2k + n − 1)/(2n)), k = 1, …, n,  H(s) = ω_cⁿ / Π(s − p_k).**   (3)

For n = 3 and ω_c = 2π·50.02 Hz = 314.3 rad/s, the poles are −314.3 and −157.1 ± 272.2j rad/s, and the denominator is s³ + 628.6 s² + 1.975 × 10⁵ s + 3.104 × 10⁷ (`main()` prints it). Odd orders have one real pole at −ω_c; conjugate pairs make every coefficient real. `butterworth_analog` builds (3) by multiplying out (s − p_k) one root at a time, and the test compares it with `scipy.signal.butter(n, ω_c, analog=True)` for n = 1–6.

### 2.3 Choosing the order and the corner

Two requirements: a loss of at most A_p dB at the band edge f_p (here 0.5 dB at 30 Hz) and an attenuation of at least A_s dB at the stopband edge f_st (here 370 Hz). The attenuation needed is set by the vibration and by what we accept: 20 log₁₀(5/0.02) = **47.96 dB**. With ε_p = 10^{A_p/10} − 1 and ε_s = 10^{A_s/10} − 1, (2) turns the two requirements into

**(f_p/f_c)²ⁿ ≤ ε_p  and  (f_st/f_c)²ⁿ ≥ ε_s,**   (4)

that is f_p ε_p^{−1/2n} ≤ f_c ≤ f_st ε_s^{−1/2n}. The range is non-empty when

**n ≥ log₁₀(ε_s/ε_p) / (2 log₁₀(f_st/f_p)).**   (5)

For the gyro (5) gives n ≥ 2.62, so **n = 3**:

| order | passband needs f_c ≥ | stopband needs f_c ≤ | |
|---|---|---|---|
| 1 | 85.88 Hz | 1.48 Hz | impossible |
| 2 | 50.76 Hz | 23.40 Hz | impossible |
| 3 | 42.60 Hz | 58.73 Hz | OK |
| 4 | 39.02 Hz | 93.05 Hz | OK |

Any corner from 42.60 to 58.73 Hz works. At the low end the passband is met exactly and the stopband has margin; at the high end the reverse. `antialias_spec` takes the geometric mean, **f_c = 50.02 Hz**, giving \|H(30 Hz)\| = 0.9775 (−0.20 dB), \|H(370 Hz)\| = 0.00247 (−52.14 dB) and \|H(390 Hz)\| = 0.00211 (−53.52 dB): the 5 deg/s line becomes 0.0105 deg/s.

### 2.4 What the filter costs: delay

For slow signals a filter acts like a pure delay. Write H(jω) = N(jω)/D(jω) with N(s) = … + b₁s + b₀ and D(s) = … + a₁s + a₀. Near ω = 0, log H(jω) ≈ log(b₀/a₀) + jω(b₁/b₀ − a₁/a₀), so the phase is ≈ −ω(a₁/a₀ − b₁/b₀) and the **group delay** −dφ/dω at DC is

**τ(0) = a₁/a₀ − b₁/b₀.**   (6)

For a first-order 1/(τs + 1) that is τ, the RC result of `P04.10`. For the Butterworth, a₁/a₀ = Σ(−1/p_k), which works out to **τ(0) = 1/(ω_c sin(π/2n))**: 1/ω_c, √2/ω_c, 2/ω_c, 2.613/ω_c for n = 1–4. The design's 6.36 ms means **−11.47° at 5 Hz** and −74.1° at the 30 Hz band edge. That lag sits inside the attitude loop. A 5 Hz crossover with 45° of phase margin loses a quarter of it to the filter. Higher order buys attenuation with delay: at the same 50 Hz corner, order 1 costs 3.18 ms and order 4 costs 8.31 ms.

### 2.5 Reconstruction

If x(t) contains nothing at or above f_s/2, the samples determine it completely, and the **Whittaker–Shannon formula** rebuilds it:

**x(t) = Σₙ x[n] sinc(f_s t − n),  sinc u = sin(πu)/(πu).**   (7)

Each sample contributes a sinc pulse that is 1 at its own instant and 0 at every other sample instant, so (7) passes through every sample; between samples the pulses add up to the unique band-limited signal with those values. Two consequences matter here. First, a real record is finite, and the missing sinc tails make (7) inexact near the ends: for 3 Hz + 11 Hz sampled at 50 Hz for 4 s, the reconstruction in the middle 2 s is within 0.0013 of the truth. Second, (7) always returns a signal below f_s/2. Feed it the samples of a 40 Hz cosine taken at 50 Hz and it returns a 10 Hz cosine (within 0.0043), because that *is* the band-limited signal with those samples. Reconstruction cannot undo aliasing; it confirms it.

### 2.6 Where the model stops

The design assumes the vibration is a line at 390 Hz. Real blade-pass frequency follows rotor speed: a 10 % change moves the line anywhere from 351 to 429 Hz and its alias anywhere from 0 Hz (a line at exactly 400 Hz looks like a bias) to 49 Hz (Problem 6). Because the filter is designed at the stopband edge and is monotonic, the residual stays below 0.0145 deg/s over that whole range, which is the reason to design at f_s − f_band rather than at the line. Real capacitors and resistors move the corner by their tolerance (the range 42.6–58.7 Hz is the margin for that). Real gyros often contain their own digital low-pass and decimation, which changes the analysis: the "analog" filter of this lesson is then the gyro's internal one, and its delay is in the data sheet. Broadband vibration above 370 Hz is attenuated by at least 52 dB here, but its folded power adds up over all the Nyquist zones it covers.

## 3. Build the computational model

Open `starter/main.py`. The polynomial, state-space and simulation helpers and `main()` are given.

### Step 1: aliasing

`aliased_frequency` is (1), vectorised with `np.round`. `dominant_frequency` removes the mean, applies a Hann window, zero-pads to 8N points and returns the frequency of the largest `rfft` bin. `sinc_reconstruct` is (7) as one matrix product: `np.sinc(fs * (t[:, None] - t0) - n[None, :]) @ samples`.

### Step 2: the prototype

```python
def butterworth_poles(n, wc):
    k = np.arange(1, n + 1)
    return wc * np.exp(1j * np.pi * (2 * k + n - 1) / (2 * n))

def butterworth_analog(n, wc):
    a = poly_from_roots(butterworth_poles(n, wc))      # given: multiplies out (s - p_k)
    return np.array([wc ** n]), np.real(a)
```

`freq_response(b, a, w)` evaluates `np.polyval(b, 1j*w) / np.polyval(a, 1j*w)`. `butterworth_order` is (4)–(5) and returns (n, f_c,min, f_c,max); subtract 10⁻¹² before `ceil` so an exact integer is not rounded up. `antialias_spec` sets f_st = f_s − f_band and A_s = 20 log₁₀(a_vib/a_allow) and raises `ValueError` if f_s ≤ 2 f_band. `group_delay_dc` is (6).

### Step 3: simulate, then sample

`simulate_analog(b, a, t, u)` converts b/a to the companion form with `tf_to_companion`, takes the exact zero-order-hold step (A_d, B_d) from `discretize_zoh` (the augmented matrix exponential of `P02.05`), and loops x ← A_d x + B_d u_k, y_k = Cx_k + Du_k. Starting from rest, it is the response a real filter gives to a signal switched on at t = 0, so the experiment discards the first 0.25 s. `sample_through_filter` keeps every 50th point of the 20 kHz grid: the sampler.

Running `main()` prints the alias table and FFT check of 2.1, the order table of 2.3, then

```
design: order 3, fc in [42.60, 58.73] Hz, chosen fc = 50.02 Hz
  |H(30 Hz)| = 0.97751 (  -0.20 dB), phase  -74.08 deg
  |H(370 Hz)| = 0.00247 ( -52.14 dB), phase  105.54 deg
  |H(390 Hz)| = 0.00211 ( -53.52 dB), phase  104.74 deg
  group delay at DC 6.36 ms; phase at 5 Hz -11.47 deg
```

followed by the SYNTHETIC experiment of part 4 and the reconstruction errors of 2.5, and writes `sampling_aliasing.png` (the Butterworth family with the specification marked, the sampled SYNTHETIC record with and without the filter, and the spectrum of the samples).

### Verify

Run `python -m pytest curriculum/phase-09/lesson-06-sampling-aliasing`. The 29 tests check:

1. **Folding**: 390 and 1190 → 10 Hz, 150 Hz unchanged, 800 → 0 Hz, every alias in [0, f_s/2]; the samples of a tone and of its alias agree to 10⁻⁹ (four frequencies).
2. **FFT check**: for 37, 390, 557 and 1210 Hz the spectral peak of 400 samples is at the predicted alias within a quarter bin.
3. **Oracle**: `butterworth_analog` against `scipy.signal.butter(analog=True)` (coefficients and `freqs`) for n = 1–6; `freq_response` against `freqs` for a general rational function.
4. **Prototype properties**: n poles on \|s\| = ω_c, all in the left half plane, in conjugate pairs; (2) to 10⁻⁹; −3.0103 dB at ω_c; 20n dB/decade.
5. **Order rule**: for three specifications, the returned n meets both at f_c,min, the mean and f_c,max, and n − 1 has an empty corner range; `antialias_spec` keeps every line above 370 Hz under 0.02 deg/s and refuses f_s < 2 f_band.
6. **Delay**: (6) against a central difference of the phase; 1/(τs + 1) gives τ.
7. **Simulator**: equal to `scipy.signal.lsim(..., interp=False)` to 10⁻⁹; unit DC gain; zero in, zero out; superposition.
8. **Reconstruction**: interpolates the samples exactly; recovers 3 + 11 Hz within 0.01; a 40 Hz tone sampled at 50 Hz returns as 10 Hz.
9. **SYNTHETIC end to end**: for orders 0, 1, 3 the measured alias equals a_vib\|H(390 Hz)\| within 3 % (plus 0.006 deg/s for the noise in one DFT bin); above 4.5 deg/s unfiltered, below 0.02 deg/s with order 3, delay increasing with order.

### Compare

`scipy.signal.butter(n, Wn, analog=True)` is the production design; it uses the same pole formula and agrees to rounding. `scipy.signal.lsim` agrees with our simulator to 10⁻⁹ when told to hold the input (`interp=False`); with its default linear interpolation of the input it differs by the difference between the two input models: 0.06 % in the filtered 390 Hz line on a 20 kHz grid (0.010542 against 0.010536 deg/s). What the hand-built version adds is (4)–(6): `butter` designs a filter when you already know n and ω_c, and this lesson is about finding them and pricing them in delay.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment on a SYNTHETIC gyro record.** The record is generated by `synthetic_imu_record` (seed 9): a true rate 2 sin(2π·1 t) + sin(2π·3 t) deg/s, a 5 deg/s line at 390 Hz, a 2 deg/s line at 780 Hz and white noise of 0.05 deg/s per 50 μs grid point. It is not a measurement, and every table and plot from it says so. Write each prediction down before running.

1. **Order sweep.** At f_c = 50.02 Hz, run orders 0–4 (`imu_experiment`), and record the predicted alias 5\|H(390 Hz)\|, the alias measured in the samples, the DC delay and the phase at 5 Hz. Predict: each order divides the alias by roughly (390/50)= 7.8; the delay follows 1/(ω_c sin(π/2n)).
2. **Corner sweep.** At order 3, repeat with f_c = 42.60, 50.02 and 58.73 Hz. Predict the trade: lower corner, smaller alias, more delay.
3. **The second harmonic.** Measure the 20 Hz tone (the 780 Hz line's alias) without a filter and with the order-3 design. Predict 2 deg/s and 2\|H(780 Hz)\| = 0.0005 deg/s; then explain why the measured value sits near 0.0005 deg/s but does not track the prediction closely.
4. **Oversampling.** Redo the design of 2.3 for f_s = 2000 Hz (stop edge 1970 Hz) and compare the order, the corner range and the delay.

Raw data format: CSV `order, fc_Hz, alias_predicted_deg_s, alias_measured_deg_s, delay_ms, phase_5Hz_deg`, with the header line "SYNTHETIC gyro record, seed 9, not a measurement".

What you should find (from `main()` and the probe of Problem 4):

| order | predicted alias | measured alias (SYNTHETIC) | DC delay | phase at 5 Hz |
|---|---|---|---|---|
| 0 | 5.0000 deg/s | 4.9969 deg/s | 0 | 0° |
| 1 | 0.6361 | 0.6361 | 3.18 ms | −5.71° |
| 2 | 0.0822 | 0.0820 | 4.50 ms | −8.13° |
| 3 | 0.0105 | 0.0109 | 6.36 ms | −11.47° |
| 4 | 0.0014 | 0.0015 | 8.31 ms | −14.99° |

Across the corner range at order 3 (seed 15 for Problem 4): 0.0065, 0.0105, 0.0171 deg/s of alias for 7.47, 6.36, 5.42 ms of delay. The 780 Hz harmonic falls from 2.01 deg/s to the noise floor of a single DFT bin, about 0.0005 deg/s. Oversampling at 2 kHz needs only order 2 (corner 50.8–124.6 Hz, delay 2.83 ms, −5.10° at 5 Hz): a faster sampler buys a gentler filter and less delay, at the price of more samples to process and a digital filter afterwards (`P09.07`).

## 5. Validate and improve

| quantity | analysis (§2) | code | independent check | difference |
|---|---|---|---|---|
| alias of 390 Hz at 400 Hz | (1): 10 Hz | 10 Hz | FFT peak 10.00 Hz | 0 |
| Butterworth n = 3 denominator | (3) | from poles | `butter(analog=True)` | < 10⁻¹⁰ relative |
| minimum order | (5): 2.62 → 3 | 3 | n = 2 corner range empty | — |
| \|H(390 Hz)\|, order 3, 50.02 Hz | (2): 0.00211 | 0.00211 | `freqs` | < 10⁻⁹ |
| DC delay, order 3 | 2/ω_c = 6.36 ms | 6.36 ms | phase slope | < 10⁻⁵ relative |
| alias after the filter | 0.0105 deg/s | — | SYNTHETIC record 0.0109 deg/s | 0.0004 deg/s, noise in one bin |
| reconstruction, 3 + 11 Hz | exact (infinite record) | 0.0013 max error | — | finite record |

The frequency-domain design and the time-domain simulation agree to within the noise in one DFT bin. That is a check of the *analysis pipeline*, not of an aircraft: the record is SYNTHETIC, and a real gyro adds things this model leaves out (moving lines, broadband vibration, its own internal filter). The engineering answer is **a third-order Butterworth at 50 Hz in front of the 400 Hz sampler: alias 0.0105 deg/s, passband loss 0.20 dB, delay 6.4 ms**.

What to improve: design against a band of vibration, not a line (the stopband-edge design already does this for anything above 370 Hz); consider oversampling plus a digital filter (experiment 4 and `P09.07`), which moves most of the filtering into software where it is exact and can be changed; and put the filter's phase into the loop design (`P09.10`, `P09.12`) instead of treating it as free.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1), (3), (5) and (6); `sampling_aliasing.png`; the experiment tables, each labelled SYNTHETIC; and a half-page design note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) for an anti-aliasing filter on a channel of your choice, giving the sample rate, band, stop edge, attenuation, order, corner, droop and delay.

### Adversarial self-checks

- **What if the vibration doubles to 10 deg/s?** The attenuation needed rises by 6.02 dB to 53.98 dB; (5) still gives n = 3, but the allowed corner range narrows (the upper limit falls by 2^{1/3}, from 58.7 to 46.6 Hz) and the chosen 50.0 Hz no longer qualifies; the new geometric-mean corner is 44.6 Hz. Recompute rather than assume.
- **What if the sample rate falls to 200 Hz?** The stop edge drops to 170 Hz, only 5.7 times the band edge: order 4 with a corner of 39.0–42.8 Hz, 10.2 ms of delay and −18.4° at 5 Hz. A slower computer makes the analog filter harder, not easier.
- **Which assumption fails first?** The single, fixed line: rotor speed changes move it, and harmonics and broadband content add aliases at other frequencies. The stopband-edge design covers all of them above 370 Hz; nothing covers content between 30 and 370 Hz except that it aliases out of the band.
- **What measurement could falsify the design?** Log the gyro on the vehicle with the motors running and the airframe still: any spectral line in 0–30 Hz that follows rotor speed is an alias that got through.
- **Does the result survive a change of time step?** Halving the simulation grid step (40 kHz) changes the measured alias from 0.0109 to 0.0104 deg/s against the unchanged prediction of 0.0105: both within the noise in one DFT bin (the noise realisation changes with the grid).
- **Which decision changes if the noise allowance halves?** At 0.01 deg/s the attenuation needed is 53.98 dB, the upper corner limit drops to 46.6 Hz, and the geometric-mean corner moves from 50.0 to 44.6 Hz, raising the delay from 6.36 to 7.14 ms; the order stays 3.

## Check yourself

1. A 1190 Hz vibration is sampled at 400 Hz. At what frequency does it appear, and why can a digital filter after the sampler not remove it?
2. Why is the stopband edge of an anti-aliasing filter f_s − f_band and not f_s/2?
3. State \|H(jω)\|² for an n-th order Butterworth low-pass. What is its gain at ω_c in dB, and how fast does it fall far above ω_c?
4. Where are the poles of a third-order Butterworth with ω_c = 1 rad/s?
5. Why does a higher-order anti-aliasing filter hurt the attitude loop even though it passes the 0–30 Hz band almost unchanged?
6. Write the sinc reconstruction formula. What does it return when fed the samples of a 40 Hz tone taken at 50 Hz?
7. How many dB down is a fourth-order Butterworth at ten times its corner?
