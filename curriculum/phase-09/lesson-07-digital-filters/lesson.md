# P09.07  Digital Filters

> A digital filter is a difference equation run once per sample. Its coefficients fix its frequency response exactly, and its delay is paid inside every loop that uses it.

| | |
|---|---|
| **Learning outcomes** | 1. Evaluate the frequency response H(e^{jω}) of a digital filter from its coefficients, implement the direct-form difference equation, and verify both against `scipy.signal.freqz` and against the DFT of the filter's own impulse response. 2. Design FIR filters (the moving average and the windowed-sinc low-pass with rectangular, Hann, Hamming and Blackman windows), prove they have exactly linear phase, and compute their delay of (N − 1)/2 samples. 3. Derive the bilinear (Tustin) transform s = K(1 − z⁻¹)/(1 + z⁻¹), its frequency warping and the prewarped K that maps one chosen frequency exactly; design a digital Butterworth and check it against `scipy.signal.butter(fs=…)` and `cont2discrete(method='bilinear')`. 4. Compare filters on a SYNTHETIC gyro record by attenuation, output noise (checked by Parseval), delay and tracking error, and choose one for an attitude loop. |
| **Prerequisites** | `P09.06` (sampling, aliasing, the Butterworth prototype), `P09.05` (frequency response), `P09.02` (convolution and the impulse response), `P02.10` (the DFT and the FFT), `P02.07` (complex exponentials), `P02.08` (the Laplace transform). Skill check: you can evaluate Σ_{k=0}^{N−1} e^{−jωk} as a geometric series and find the roots of a quadratic. |
| **Engineering question** | A gyro is sampled at 1000 Hz behind its anti-aliasing filter; the attitude loops need 0–10 Hz, and a 3 deg/s vibration line at 180 Hz plus 0.3 deg/s of white noise ride on top. Which digital filter, a 50-point moving average, a 101-tap windowed-sinc FIR or a second- or fourth-order Butterworth IIR at 40 Hz, removes the line and the noise at the smallest cost in delay, and how big is that delay? Answer: the moving average nulls the line exactly but costs 24.5 ms and 36 % of the signal at 10 Hz; the FIR is nearly perfect in magnitude but costs 50 ms; the fourth-order IIR leaves 0.0047 deg/s of the line for 10.4 ms; the second-order IIR leaves 0.12 deg/s for 5.6 ms and tracks the true rate best. |
| **Assumptions and limits** | Linear time-invariant filters in exact floating-point arithmetic (fixed-point coefficients and rounding are treated only in Problem 6); a constant sample rate with no jitter; filters start from rest, so the first 0.5 s of every record is discarded as start-up transient; the input is already free of aliases (`P09.06`). The gyro record is **SYNTHETIC** (seeded; formula in `synthetic_gyro_record`): a smooth manoeuvre, one vibration line and white noise. Real vibration lines move with rotor speed, and real gyros add bias, drift and scale-factor errors. |
| **Analytical work** | H(e^{jω}) = B(e^{jω})/A(e^{jω}) from the z-transform of the difference equation; the moving average as a geometric series (Dirichlet kernel) and its zeros at k f_s/N; linear phase of symmetric FIR taps; the ideal low-pass impulse response 2f_c/f_s sinc(2f_c(k − M)/f_s) and why it must be windowed; the bilinear map of the left half plane into the unit circle, the warping law f_d = (f_s/π) arctan(π f_a/f_s) and the prewarped K = 2πf_p/tan(πf_p/f_s); group delay from the coefficients; Parseval's relation for filtered white noise. |
| **Computational deliverable** | `starter/main.py`: `digital_freq_response`, `difference_equation`, `impulse_response`, `group_delay`, `is_stable_digital`, `moving_average`, `window`, `windowed_sinc_lowpass`, `bilinear`, `butterworth_digital` (`poly_from_roots`, `butterworth_analog`, `analog_freq_response`, `tone_amplitude`, `synthetic_gyro_record`, `lesson_filters`, `gyro_filter_comparison` and `main` are given). Tests prove: the response equals `freqz` for FIR and IIR filters, including a[0] ≠ 1; the impulse response is the taps of an FIR and aⁿ for a first-order IIR, and its DFT equals `freqz`; the difference equation is linear, time-invariant and causal; the group-delay formula equals the phase slope; the moving average nulls k f_s/N; the windows match NumPy's; the windowed sinc is symmetric, −6 dB at f_c, delayed by (N − 1)/2 at every passband frequency and below its window's stopband floor; the bilinear transform equals `cont2discrete`, preserves stability, and with prewarping matches the analog response at f_p; the digital Butterworth equals `butter(fs=…)` for n = 1–4; and on the SYNTHETIC record the filtered noise obeys Parseval within 10 %. |
| **Practical activity** | **computational experiment** on a SYNTHETIC gyro record: run four filters (and a notch, Problem 3), measure the residual vibration, the output noise, the delay and the tracking error, and map attenuation against delay. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `digital_filters.png`, the derivations of the moving-average response, the linear-phase property and the bilinear warping law, the comparison table (labelled SYNTHETIC), and a half-page recommendation of a filter for the attitude loop with its attenuation, delay and phase at the crossover frequency. |
| **Assessment** | `assessments/phase-09/grading-rubric.md`; competencies 1–4. Pass: the three derivations are correct; all public tests pass; the recommendation weighs delay against attenuation with numbers; the SYNTHETIC record is labelled as such. |
| **References** | `@mit-ocw-6-003-signals-and-systems` Lecture 17 "Discrete-Time (DT) Frequency Representations", Lecture 18 "Discrete-Time (DT) Fourier Representations" and Lecture 20 "Applications of Fourier Transforms" (required-free: frequency response of difference equations, filtering); Lecture 21 "Sampling" (supplementary-free). `@scipy-userguide` "Signal Processing (scipy.signal)" (supplementary-free: `freqz`, `butter`, `cont2discrete`, used only as oracles). Windowed-sinc design and the prewarped bilinear transform are original instruction (no free source mapped). Map `references/phase-09.md` not yet available when written; only ids already in `sources.yaml` are cited, with lecture titles as recorded in `references/phase-02.md`. |
| **Extensions** | ★ Implement the transposed direct form II and show it gives the same output as direct form I with half the memory. ★ Run the fourth-order Butterworth at 5 Hz with its coefficients rounded to 16-bit fixed point, then as two second-order sections, and compare pole positions. ★ Filter a logged record offline forwards and then backwards (zero-phase filtering) and show the delay vanishes and the magnitude is squared; explain why a flight controller cannot do this. |
| **Time** | unmeasured |

## 1. The engineering problem

`P09.06` put an analog filter in front of the sampler so that nothing aliases into the band. What arrives in the flight computer is clean of aliases but not clean: a 180 Hz vibration line, below the 500 Hz Nyquist frequency of a 1 kHz sampler, is represented faithfully, and so is broadband sensor noise. The attitude loops need 0–10 Hz. The rest must be removed in software, by a **digital filter**: each new output is a weighted sum of recent inputs and, for a recursive filter, of recent outputs. The flight software runs that sum once per sample, 1000 times a second.

A digital filter is exact in a way an analog one never is: its coefficients are numbers, not components with tolerances, and its frequency response follows from them by a formula. It can do things an analog filter cannot, such as exactly linear phase or an exact zero at 180 Hz. But it obeys the same trade as its analog cousin. Every filter that removes high frequencies delays the low ones, and in a control loop that delay is phase lag at the crossover frequency. A 50 ms delay is 18° at a 1 Hz crossover; for the 5 Hz crossover of a small multirotor's attitude loop it is 90°, which no loop survives.

This lesson builds the two families, FIR (finite impulse response: no feedback, output a weighted sum of inputs) and IIR (infinite impulse response: recursive), from the difference equation up, checks every one against `scipy.signal`, and runs them on a SYNTHETIC gyro record. Success criteria, written first: every frequency response matches `freqz` to 10⁻¹⁰; the DFT of the impulse response computed by our own difference equation matches `freqz` to 10⁻⁹; the bilinear transform matches `cont2discrete` and the digital Butterworth matches `butter(fs=…)` to 10⁻⁸; the measured output noise on the SYNTHETIC record matches σ√(Σh²) within 10 %.

## 2. Understand and derive

### 2.1 The difference equation and its frequency response

A linear time-invariant digital filter is

**a₀y[n] + a₁y[n−1] + … + a_M y[n−M] = b₀x[n] + b₁x[n−1] + … + b_L x[n−L].**   (1)

Apply it to x[n] = e^{jωn}, ω = 2πf/f_s in radians per sample. Linearity and time invariance make the output H e^{jωn} with a single complex number H. Substituting and dividing by e^{jωn}:

**H(e^{jω}) = (Σ_k b_k e^{−jωk}) / (Σ_k a_k e^{−jωk}).**   (2)

This is the z-transform transfer function H(z) = B(z⁻¹)/A(z⁻¹) on the unit circle z = e^{jω}. It is periodic in f with period f_s and, for real coefficients, \|H\| is even, so 0 ≤ f ≤ f_s/2 tells everything. A filter with all a_k = 0 for k ≥ 1 is **FIR**: its **impulse response** h[n] (the output for x = δ[n]) is the list of taps b_n, and it is always stable. A filter with feedback is **IIR**: its impulse response never ends, and it is stable when every root of A(z) (every **pole**) lies inside the unit circle \|z\| < 1. The pole at z_p contributes a term z_pⁿ to h[n], which decays only if \|z_p\| < 1.

Two independent checks follow from (2). The DFT of the impulse response, h[n] computed by running (1) on δ[n], equals H at the DFT frequencies (to the extent the impulse response has died out within the record). And the **group delay** −dφ/dω can be read from the coefficients: for C(ω) = Σc_k e^{−jωk}, −d arg C/dω = Re[Σ k c_k e^{−jωk} / C], and the delay of B/A is that of B minus that of A.

### 2.2 FIR filters

**Moving average.** b_k = 1/N for k = 0, …, N − 1. The geometric series gives

**H(e^{jω}) = (1/N) e^{−jω(N−1)/2} sin(Nω/2)/sin(ω/2),**   (3)

a real function (the Dirichlet kernel) times a pure delay of (N − 1)/2 samples. Its zeros are at ω = 2πk/N, that is **f = k f_s/N**. With N = 50 at 1 kHz they are at 20, 40, …, 180 Hz: the 180 Hz line falls exactly on a zero, \|H(180 Hz)\| = 3 × 10⁻¹⁷. The price is in the band: \|H(10 Hz)\| = 0.637, a 36 % loss of a 10 Hz motion, and 24.5 ms of delay.

**Linear phase.** If the taps are symmetric, h[k] = h[N−1−k], pair the terms k and N−1−k in (2): each pair is e^{−jω(N−1)/2}·2h[k]cos(ω((N−1)/2 − k)). So H = e^{−jω(N−1)/2}·(real function): **the phase is exactly −ω(N−1)/2 wherever the real factor is positive**, a constant delay of (N − 1)/2 samples at every frequency. Every frequency component is delayed equally and the waveform is not distorted. No causal IIR filter has this property.

**Windowed sinc.** The ideal low-pass with cutoff f_c has impulse response 2(f_c/f_s) sinc(2(f_c/f_s)k), infinitely long in both directions. Truncate it to N = 2M + 1 taps, shift it by M samples to make it causal, multiply by a **window** w[k] that tapers to the ends, and normalise to unit DC gain:

**h[k] = 2(f_c/f_s) sinc(2(f_c/f_s)(k − M)) w[k] / Σ(…).**   (4)

Abrupt truncation (the rectangular window) causes ripple whose worst sidelobe stays at about −35 dB however many taps are used (the Gibbs phenomenon of `P02.09`). Tapering trades a wider transition band for lower sidelobes: Hann 0.5 − 0.5 cos(2πk/(N−1)), Hamming 0.54 − 0.46 cos(…), Blackman 0.42 − 0.5 cos(…) + 0.08 cos(2…). For 101 taps at f_c = 40 Hz, `main()` measures:

| window | \|H(40 Hz)\| | worst gain above 80 Hz | \|H(180 Hz)\| |
|---|---|---|---|
| rectangular | 0.513 | −35.0 dB | −50.3 dB |
| Hann | 0.499 | −62.4 dB | −97.3 dB |
| Hamming | 0.500 | −63.6 dB | −73.1 dB |
| Blackman | 0.500 | −75.8 dB | −106.4 dB |

All are −6 dB (gain ½) at f_c, not −3 dB: the windowed sinc's transition band is centred on the cutoff. Every one of them delays the signal by exactly 50 samples, **50 ms**.

### 2.3 IIR filters by the bilinear transform

An analog prototype H_a(s) already has the sharpness we want for few coefficients (`P09.06`). To turn it into a difference equation, approximate the integrator 1/s by the trapezoidal rule. Over one step, y[n] = y[n−1] + (T/2)(x[n] + x[n−1]), so Y/X = (T/2)(1 + z⁻¹)/(1 − z⁻¹), and 1/s ↔ that, i.e.

**s = K (1 − z⁻¹)/(1 + z⁻¹),  K = 2f_s.**   (5)

Substituting (5) into H_a(s) gives the **bilinear** (Tustin) digital filter. Three properties:

- **Stability is preserved.** Solving (5) for z gives z = (K + s)/(K − s). For Re s < 0, \|K + s\| < \|K − s\|, so \|z\| < 1: the left half plane maps inside the unit circle and the imaginary axis onto it.
- **Frequencies are warped.** On the circle, z = e^{jω} gives s = jK tan(ω/2). The whole analog axis 0 ≤ Ω < ∞ is squeezed into 0 ≤ ω < π, so the digital response at f_d equals the analog response at f_a = (K/2π) tan(πf_d/f_s). With K = 2f_s, an analog corner at f_a lands at **f_d = (f_s/π) arctan(πf_a/f_s)**, always lower: a 250 Hz Butterworth at 1 kHz lands at 211.9 Hz and has \|H(250 Hz)\| = 0.525 instead of 0.707.
- **Prewarping fixes one frequency.** Choose K = 2πf_p/tan(πf_p/f_s); then tan(πf_p/f_s)K/(2π) = f_p and the digital response at f_p equals the analog response at f_p *exactly*. The digital Butterworth with its −3 dB point exactly at f_c is the analog prototype at f_c mapped with K prewarped at f_c, which is what `scipy.signal.butter(n, f_c, fs=f_s)` returns.

For n = 2, f_c = 40 Hz, f_s = 1 kHz: b = [0.013359, 0.026718, 0.013359], a = [1, −1.64746, 0.700897], poles at \|z\| = 0.8372. The analog prototype's n zeros at s = ∞ all map to z = −1: the response is exactly zero at f_s/2.

### 2.4 What each filter costs

The delay of an IIR filter is not constant with frequency, but at low frequency it matches the analog prototype: 1/(ω_c sin(π/2n)) from `P09.06`, √2/(2π·40) = 5.63 ms for n = 2 against 5.60 ms computed from (2.1) on the digital filter. At 5 Hz, where an attitude loop might cross over:

| filter | phase at 5 Hz | delay at 2 Hz |
|---|---|---|
| MA 50 | −44.1° | 24.5 ms |
| FIR 101 Hamming | −90.0° | 50.0 ms |
| Butterworth 2, 40 Hz | −10.1° | 5.6 ms |
| Butterworth 4, 40 Hz | −18.7° | 10.4 ms |

### 2.5 Noise through a filter

White noise of standard deviation σ (independent samples) through a filter with impulse response h has output variance σ²Σh[n]², because the output is Σh[k]x[n−k] and the cross terms average to zero. By Parseval this equals σ² times the average of \|H\|² over the band. The moving average has Σh² = 1/N, so σ_out = σ/√N = 0.3/√50 = 0.0424 deg/s. This is the third independent check of the build: the measured noise on the SYNTHETIC record against σ√(Σh²).

### 2.6 Where the model stops

Floating-point arithmetic hides the main practical danger of IIR filters: **coefficient sensitivity**. When f_c ≪ f_s the poles crowd towards z = 1, and a high-order direct-form denominator becomes extremely sensitive to rounding. A fourth-order Butterworth at 5 Hz and 1 kHz has poles at \|z\| = 0.988; rounding its coefficients to four decimals moves a pole to \|z\| = 1.051 and the filter becomes unstable (Problem 6). Real implementations use cascaded second-order sections. Filters also start from rest: the first few time constants of the output are a transient, which a flight controller sees at every reset.

## 3. Build the computational model

Open `starter/main.py`. The analog prototype, the SYNTHETIC record and the comparison are given.

### Step 1: evaluation and the difference equation

`digital_freq_response(b, a, f, fs)` is (2): with z⁻¹ = e^{−j2πf/f_s}, Σc_k z^{−k} is `np.polyval(c[::-1], zinv)`. `difference_equation` loops over n, keeping zero-padded copies of x and y so that every index is valid:

```python
for n in range(x.size):
    acc = b_rev @ xp[n:n + nb]                 # sum b_k x[n-k]
    if na > 1:
        acc -= a_rev @ yp[n:n + na - 1]        # minus sum_{k>=1} a_k y[n-k]
    yp[n + na - 1] = acc / a[0]
```

`impulse_response` runs it on δ[n]. `group_delay` is the formula of 2.1 applied to b and a. `is_stable_digital` checks `np.abs(np.roots(a)) < 1`.

### Step 2: FIR design

`moving_average(N)` is `np.full(N, 1/N)`. `window(name, n)` evaluates the four formulas with x = 2πk/(n − 1). `windowed_sinc_lowpass` is (4); require an odd number of taps (a whole-sample delay) and 0 < f_c < f_s/2.

### Step 3: bilinear transform

`bilinear(b, a, fs, f_prewarp=None)` substitutes (5) term by term. For an analog polynomial of degree ≤ N, the term c_k s^k becomes c_k K^k (1 − z⁻¹)^k (1 + z⁻¹)^{N−k} after multiplying numerator and denominator by (1 + z⁻¹)^N; build each with repeated `np.convolve` and add. Normalise so a[0] = 1. `butterworth_digital(n, fc, fs)` is `bilinear(*butterworth_analog(n, 2*np.pi*fc), fs, f_prewarp=fc)`.

Running `main()` prints the window table of 2.2, the warping lines of 2.3, the Butterworth coefficients, the SYNTHETIC comparison of part 4, and writes `digital_filters.png` (magnitude of the four filters with the 180 Hz line marked, their group delay over 0–100 Hz, and a 0.6 s window of the SYNTHETIC record with the FIR and IIR outputs).

### Verify

Run `python -m pytest curriculum/phase-09/lesson-07-digital-filters`. The 31 tests check:

1. **Oracle**: `digital_freq_response` against `scipy.signal.freqz` for an FIR, a first-order IIR, a `butter`-designed IIR and a filter with a[0] = 2, at 77 frequencies.
2. **Impulse response**: the taps of an FIR; 0.8ⁿ and 1.5·0.5ⁿ for first-order IIRs; and, the independent check, its 4096-point DFT equal to `freqz` at the DFT bins to 10⁻⁹ for all four filters.
3. **LTI**: superposition, time shift by 25 samples, causality.
4. **Group delay** against a central difference of the phase; **stability** of poles inside, outside and on the circle.
5. **Moving average**: unit DC gain, zeros at 20, 40, …, 480 Hz, delay 24.5 samples.
6. **Windows** equal to `np.hanning`, `np.hamming`, `np.blackman`.
7. **Windowed sinc** (three windows): symmetric taps, unit DC gain, gain ½ at f_c, delay (N − 1)/2 at 20 passband frequencies, stopband below −50, −70 and −40 dB beyond f_c + 6f_s/N.
8. **Bilinear** against `cont2discrete(method='bilinear')` for three analog filters; stable poles inside, unstable outside; prewarped response equal to the analog one at f_p (40 and 250 Hz) while plain Tustin falls below −3 dB.
9. **Digital Butterworth** against `butter(n, 40, fs=1000)` for n = 1–4, with \|H(f_s/2)\| = 0.
10. **SYNTHETIC comparison**: output noise within 10 % of σ√(Σh²) for all four filters; the moving average nulls the line; the FIR delay is 50.0 ms; IIR 2 < IIR 4 < FIR in delay; the residual line equals 3\|H(180 Hz)\|.

### Compare

`scipy.signal.freqz`, `butter(fs=…)` and `cont2discrete(method='bilinear')` agree with the hand-built versions to rounding: they implement the same formulas, so the agreement shows our algebra is right, not that the formulas are. The independent physics checks are the impulse response (our own difference equation, transformed, reproduces `freqz`) and Parseval on the SYNTHETIC noise. SciPy would normally run the filter with `lfilter` or, better for high orders, `sosfilt` on second-order sections; our direct-form loop is slower and, for poles near z = 1, numerically worse (2.6).

## 4. Experiment at home, or investigate real data

**Mode: computational experiment on a SYNTHETIC gyro record.** `synthetic_gyro_record(seed=7)`: 4 s at 1 kHz of a true rate 5 sin(2π·0.5t) + 2 sin(2π·2t) deg/s, a 3 deg/s line at 180 Hz and 0.3 deg/s white noise. It is not a measurement, and every table and plot from it says so. Write each prediction first.

1. **Four filters.** Run `gyro_filter_comparison()`. Predict for each: the residual line from 3\|H(180 Hz)\|; the output noise from 0.3√(Σh²); the delay at 2 Hz; and which tracks the true rate best (the root-mean-square error against the true rate, which counts delay as error).
2. **Taps against delay.** For the Hamming FIR at 40 Hz, use 51, 101 and 201 taps. Predict the delay (25, 50, 100 ms) and roughly how the 180 Hz attenuation changes.
3. **A notch instead** (Problem 3): an IIR notch at 180 Hz from the bilinear transform of (s² + ω₀²)/(s² + (ω₀/Q)s + ω₀²) with prewarping. Predict its delay at 2 Hz and its effect on the noise.

Raw data format: CSV `filter, atten_180_dB, residual_deg_s, noise_out_deg_s, noise_parseval_deg_s, delay_2Hz_ms, rms_error_deg_s`, with the header line "SYNTHETIC gyro record, seed 7, not a measurement".

What you should find (seed 7, from `main()`):

| filter | \|H(180 Hz)\| | residual line | noise out (Parseval) | delay at 2 Hz | RMS error vs true |
|---|---|---|---|---|---|
| MA 50 | zero | 0.0000 deg/s | 0.0427 (0.0424) | 24.5 ms | 0.513 deg/s |
| FIR 101 Hamming, 40 Hz | −73.1 dB | 0.0007 | 0.0851 (0.0808) | 50.0 ms | 1.038 |
| Butterworth 2, 40 Hz | −28.0 dB | 0.1188 | 0.0921 (0.0891) | 5.6 ms | 0.172 |
| Butterworth 4, 40 Hz | −56.1 dB | 0.0047 | 0.0895 (0.0859) | 10.4 ms | 0.236 |

The FIR is the best filter by magnitude and the worst by tracking error: its 50 ms delay applied to a 5 deg/s, 0.5 Hz motion is an error of up to 2π·0.5·5·0.05 = 0.79 deg/s. The moving average has the least noise (its 50 taps average the most) but loses 36 % at 10 Hz. The second-order Butterworth leaves the most vibration and still has the smallest error against the true rate, because delay costs more than 0.12 deg/s of residual line. For taps 51, 101, 201: −68.8, −73.1, −78.8 dB at 180 Hz for 25, 50 and 100 ms; the attenuation grows slowly, the delay in proportion. The notch (Q = 2) removes the line completely for 0.39 ms of delay at 2 Hz but leaves the noise untouched (0.2675 deg/s out for 0.3 in).

## 5. Validate and improve

| quantity | analysis (§2) | code | independent check | difference |
|---|---|---|---|---|
| MA zeros | k f_s/N = 20 k Hz | \|H\| < 10⁻¹² at 20…480 Hz | — | rounding |
| FIR delay | (N − 1)/2 = 50 samples | group delay 50.000 | phase slope | < 10⁻⁸ |
| Butterworth 2, 40 Hz coefficients | (5) prewarped | b, a above | `butter(fs=1000)` | < 10⁻⁸ relative |
| response of all filters | (2) | `digital_freq_response` | `freqz`; DFT of impulse response | < 10⁻⁹ |
| Tustin corner, 250 Hz | (f_s/π) arctan(π·250/1000) = 211.9 Hz | \|H(250)\| = 0.525 | prewarped: 0.7071 | warping |
| output noise, MA 50 | 0.3/√50 = 0.0424 | — | SYNTHETIC 0.0427 | 0.7 %, finite record |
| output noise, FIR 101 | 0.0808 | — | SYNTHETIC 0.0851 | 5 %, finite record |

Every response agrees with two independent routes. The noise agreement is statistical: 3500 output samples of a correlated process estimate a standard deviation to a few per cent, and the FIR, whose output is correlated over 101 samples, has the fewest independent samples and the largest deviation. These checks test the code and the theory on a SYNTHETIC record; a real gyro adds bias and moving lines.

The engineering answer for an attitude loop with a crossover near 5 Hz: **the second-order Butterworth at 40 Hz (5.6 ms, −10.1° at 5 Hz), followed by a notch at the vibration line if the residual 0.12 deg/s matters**. That combination costs about 6 ms of delay. The FIR and the moving average are good for logging and for navigation filters that run slower and can tolerate (or compensate) a known delay; they are poor choices inside a fast loop.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2), (3), linear phase and the warping law; `digital_filters.png`; the comparison tables labelled SYNTHETIC; and a half-page recommendation (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) of a filter for a loop of your choice, with its attenuation at the line, its delay and its phase at the crossover frequency.

### Adversarial self-checks

- **What if the vibration line moves to 170 Hz?** The moving average's zero is at 160 or 180 Hz, and \|H(170 Hz)\| of the 50-point average is no longer zero; the notch with Q = 10 passes 0.82 of the line at 170 Hz; the Butterworth filters change little. Exact nulls are fragile; monotonic roll-off is robust.
- **What if the sample rate doubles to 2 kHz with the same coefficients?** Every frequency in Hz doubles: the Butterworth corner moves to 80 Hz and the moving-average zeros to multiples of 40 Hz. Digital filters are designed in cycles per sample, not in hertz.
- **Which assumption fails first?** Exact arithmetic: a low-corner, high-order IIR in direct form on a fixed-point processor (Problem 6). Then the single, fixed vibration line.
- **What measurement could falsify the filter model?** Feed a recorded step or chirp through the flight code and compare with `difference_equation` on the same input; any difference is an implementation error (a coefficient typo, a wrong sign convention for a[1:], a reset).
- **Does the result survive a change of record length?** The residual line and the delay do not depend on it; the measured noise approaches the Parseval value as the record grows.
- **Which decision changes if the delay budget halves to 3 ms?** No filter in the table meets it except the notch; the second-order Butterworth would need its corner raised to about 75 Hz, and the residual line grows accordingly.

## Check yourself

1. Write H(e^{jω}) for a filter with coefficients b and a. What makes a filter FIR, and why is an FIR always stable?
2. A 10-tap moving average runs at 500 Hz. Where are its first three zeros?
3. Why does a symmetric FIR have exactly linear phase, and what is its delay in samples?
4. Why must the ideal sinc impulse response be windowed, and what does a Blackman window buy compared with a rectangular one?
5. State the bilinear transform. Where does the left half of the s-plane go?
6. A 250 Hz analog corner is mapped by plain Tustin at 1 kHz. Where does it land, and how does prewarping fix it?
7. White noise of σ = 0.3 deg/s passes a 50-point moving average. What is the output standard deviation?
