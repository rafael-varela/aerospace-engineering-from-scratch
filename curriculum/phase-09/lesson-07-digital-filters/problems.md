# P09.07  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Every gyro record you generate is SYNTHETIC: label it so in every table and plot.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a moving average on a 1 kHz gyro)

A flight controller averages the last N = 20 + s gyro samples at f_s = 1000 Hz.

(a) Where are its first seven zeros? Does it null a 180 Hz line?
(b) Give \|H\| at 2 Hz, 10 Hz and 180 Hz, and the delay in samples and milliseconds.
(c) Check (a) and (b) with `moving_average`, `digital_freq_response` and `group_delay`.

Answer form: seven frequencies and a yes/no; three gains and a delay; a confirmation.

## 2. (Derivation: three properties)

(a) Sum the geometric series to derive (3), the moving-average response, and its zeros.
(b) Prove that symmetric taps h[k] = h[N−1−k] give a phase of exactly −ω(N − 1)/2 wherever the real factor is positive. What happens to the phase where it is negative?
(c) Derive the bilinear transform from the trapezoidal rule, show that it maps Re s < 0 into \|z\| < 1 and the imaginary axis onto the unit circle, and derive the warping law f_d = (f_s/π) arctan(π f_a/f_s) for K = 2f_s and the prewarped K = 2πf_p/tan(πf_p/f_s).

Answer form: three derivations.

## 3. (Unfamiliar scenario: a notch for the vibration line)

Instead of a low-pass, the vibration line at 180 Hz is removed with a second-order notch, the bilinear transform (prewarped at 180 Hz) of the analog notch H_a(s) = (s² + ω₀²)/(s² + (ω₀/Q)s + ω₀²), ω₀ = 2π·180 rad/s.

(a) For Q = 2 and Q = 10 compute the digital coefficients, \|H(180 Hz)\|, \|H(2 Hz)\|, the −3 dB band of the notch, \|H(170 Hz)\| and \|H(190 Hz)\|, the delay at 2 Hz and the pole radius.
(b) Run each notch on `synthetic_gyro_record(seed=s)` and report the residual line, the output noise and the RMS error against the true rate (after 0.5 s).
(c) When is a notch the right tool, and when is it not?

Answer form: a table, three numbers per Q, two sentences.

## 4. (Numerical, uses your code: your own SYNTHETIC comparison)

(a) Run `gyro_filter_comparison(seed=s)` and tabulate the residual line, the output noise against Parseval, the delay at 2 Hz and the RMS error for the four filters.
(b) For the Hamming FIR at 40 Hz with 51, 101 and 201 taps, give \|H(180 Hz)\| in dB and the delay.
(c) Rank the four filters for an attitude loop with a 5 Hz crossover, and say which number decided the ranking.

Answer form: a table labelled SYNTHETIC, three pairs, a ranking with a reason.

## 5. (Frequency warping)

Map a second-order analog Butterworth with corner f_c = 50, 150, 300 and 450 Hz to f_s = 1000 Hz with plain Tustin and with prewarping.

(a) For each, find \|H(f_c)\| of the Tustin filter and the frequency where it is actually −3 dB, and compare with the warping law.
(b) For which corners is prewarping essential, and for which is it a refinement?

Answer form: a table of four rows, one sentence.

## 6. (Judgment: arithmetic and delay)

(a) Design a fourth-order Butterworth at f_c = 40, 5 and 1 Hz (f_s = 1000 Hz). Give the largest pole radius of each, then round the denominator coefficients to four decimal places and give the new largest pole radius and DC gain. Explain what happens and what a flight-software engineer does about it.
(b) Give the phase at 5 Hz of each of the four lesson filters. Which are acceptable inside a loop with a 5 Hz crossover and 45° of phase margin?
(c) An analyst wants to remove the 180 Hz line from a logged flight record with no delay at all. How can she, and why can the flight controller not do the same?

Answer form: a table and a paragraph; four phases and a verdict; two sentences.
