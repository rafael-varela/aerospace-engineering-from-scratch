# P09.06  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Every gyro or accelerometer record you generate is SYNTHETIC: label it so in every table and plot.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: where do the propeller lines go?)

A gyro is sampled at f_s = 400 Hz. A propeller produces a blade-pass line at f₁ = 355 + 2s Hz and harmonics at 2f₁ and 3f₁.

(a) Compute the alias of each line with (1) of the lesson.
(b) Which of them fall inside the 0–30 Hz control band?
(c) Confirm each alias with `dominant_frequency` on 400 samples of a cosine at the true frequency.

Answer form: three frequencies, a list, three FFT peaks.

## 2. (Derivation: the Butterworth prototype and its delay)

(a) Starting from \|H(jω)\|² = 1/(1 + (ω/ω_c)²ⁿ) and H(s)H(−s) = \|H(jω)\|² at s = jω, show that the poles of H(s)H(−s) are 2n points equally spaced on \|s\| = ω_c, and that the left-half-plane ones are p_k = ω_c exp(jπ(2k + n − 1)/(2n)), k = 1, …, n.
(b) From the two inequalities (4), derive the minimum-order formula (5) and the allowed corner range.
(c) Show that the DC group delay of any rational H = N/D is a₁/a₀ − b₁/b₀, and that for the Butterworth it equals Σ_k(−1/p_k). Evaluate it for n = 1, 2, 3 and compare with 1/(ω_c sin(π/2n)).

Answer form: three derivations, three closed-form delays.

## 3. (Unfamiliar scenario: a launch-vehicle accelerometer)

A structural accelerometer on a launch vehicle records 0–100 Hz with at most 1 dB of loss. Engine and turbopump vibration reach 1 g at high frequency, and the alias of anything must stay below 1 mg (60 dB of attenuation).

(a) With f_s = 1000 Hz, find the stopband edge, the Butterworth order and the allowed corner range. What is the DC delay at the geometric-mean corner, and the phase at 20 Hz?
(b) Repeat for f_s = 4000 Hz.
(c) The structural-dynamics team wants the smallest possible phase error at 20 Hz. Which sample rate do you recommend, and what does it cost?

Answer form: two designs (order, range, delay, phase), one recommendation with a reason.

## 4. (Numerical, uses your code: your own SYNTHETIC gyro record)

Run `imu_experiment(orders=(3,), fc=..., seed=s)` with the lesson's channel for f_c = f_c,min, the geometric mean and f_c,max, plus `orders=(0,)`.

(a) Tabulate the predicted and measured 10 Hz alias, the DC delay and the phase at 5 Hz for each corner.
(b) Measure the 20 Hz alias of the 780 Hz harmonic (use `sample_through_filter` and `tone_amplitude`) with and without the filter, and compare with 2\|H(780 Hz)\|.
(c) Explain any disagreement larger than 3 % in (a) or (b).

Answer form: a table labelled SYNTHETIC, two numbers per corner for (b), a short explanation.

## 5. (Reconstruction from a finite record)

Sample cos(2πft + 0.3) at f_s = 50 Hz for N = 50, 200 and 800 samples, with f = 5, 12.5 and 22.5 Hz. Reconstruct with `sinc_reconstruct` on the middle half of each record and record the maximum error.

(a) Tabulate the nine errors.
(b) Why does the error fall as N grows, and why is it largest for 22.5 Hz?
(c) What does this imply for a flight computer that wants to "upsample" a 400 Hz gyro to 1600 Hz for a fast inner loop?

Answer form: a 3 × 3 table, two short arguments, one sentence.

## 6. (Judgment: is the design good enough?)

(a) Rotor speed varies by ±10 % about the nominal, so the blade-pass line moves between 351 and 429 Hz. Where can its alias lie, and what is the worst residual with the order-3, 50.02 Hz design? Is the stopband-edge design robust to this?
(b) A colleague proposes to drop the analog filter and remove the 10 Hz alias with a digital notch after the sampler. Give two reasons this fails.
(c) The order-3 filter lags 11.5° at 5 Hz. Is that a property of the filter or of the corner choice? What would you change if the attitude loop's phase margin were only 35°?
(d) What does oversampling at 2 kHz change (order, corner, delay)? What new work does it create?

Answer form: a frequency range and a residual; two reasons; two sentences; numbers and one sentence.
