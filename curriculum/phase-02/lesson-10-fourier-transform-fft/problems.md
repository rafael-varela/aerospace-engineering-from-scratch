# P02.10  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: the three numbers of a record)

A strain-gauge channel is logged for N = 4096 samples at f_s = 2000 + 10s Hz.

(a) Give the record length, the frequency resolution Δf, the Nyquist frequency and the number of bins in the single-sided amplitude spectrum.
(b) A tone of 3.0 με amplitude sits exactly on bin k = 100. What are |X[100]| and |X[N − 100]|, and what does `amplitude_spectrum` return at that bin?
(c) How long a record would resolve two tones 0.1 Hz apart with a rectangular window? With a Hann window (main lobe twice as wide)?

Answer form: six numbers with units; one sentence on (c).

## 2. (Derivation: orthogonality, Parseval, leakage, and a DFT by hand)

(a) Prove that Σ_{k=0}^{N−1} e^{2πik(m−n)/N} = N δ_{mn} using the finite geometric series, and hence derive the inverse DFT (2) from the forward one (1).
(b) Prove Parseval's theorem (3) for the DFT.
(c) A complex tone e^{2πi(k₀+δ)n/N}, 0 ≤ δ < 1, is transformed. Show that |X[k₀]| = |sin πδ / sin(πδ/N)| and that for large N this is N sin(πδ)/(πδ). Evaluate the factor at δ = 0, ¼, ½.
(d) Compute by hand, from (1), the DFT of x = [1, 0, −1, 0]. Interpret the result as the spectrum of a sampled cosine: which bins, which amplitude, and why two bins?

Answer form: three short proofs, one 4-element vector with interpretation.

## 3. (Unfamiliar scenario: a helicopter's vibration logged too slowly)

A helicopter's main rotor turns at 258 rpm with four blades; its tail rotor turns at 1190 rpm with two blades. The dominant vibration tones are the blade-passing frequencies (rotor speed × blade count) and the once-per-revolution tones. A maintenance logger is set to f_s = 40 + s Hz.

(a) List the four tones in hertz.
(b) For each, state whether it is below the Nyquist frequency and, if not, where it appears in the logged spectrum (equation (6)). Check with `alias_frequency`.
(c) The maintenance manual says "a peak near 15 Hz indicates a loose tail-rotor gearbox". Explain how an aliased blade-passing tone could produce a false alarm, and name the single hardware change that prevents it.
(d) What is the lowest sampling rate that puts every tone below Nyquist with a 20 % margin?

Answer form: a table of four tones with alias frequencies, two short paragraphs, one number.

## 4. (Numerical, uses your code: the same tone, three window lengths)

A pure tone of amplitude 2.0 m/s² at f₀ = 80 + s/7 Hz (no noise) is sampled at 1024 Hz for 1 s, 2 s and 4 s.

(a) For each record length give Δf, the nearest bin to f₀, and the offset δ in bins.
(b) Predict the rectangular-window reading of the peak bin from 2 sin(πδ)/(πδ), then run `amplitude_spectrum` and compare.
(c) Repeat with `window="hann"`. Which window's reading varies less with record length?
(d) Explain in two sentences why the "measured" amplitude changed while the signal did not, and what you would tell a colleague who reports a vibration amplitude without stating the window and record length.

Answer form: a table with three rows (T, Δf, δ, predicted, rectangular, Hann) and a short paragraph.

## 5. (Numerical, uses your code: filtering through the convolution theorem)

A first-order low-pass filter with corner f_c = 40 Hz, sampled at f_s = 1024 Hz, has impulse response h[n] = (1 − α) αⁿ with α = e^{−2π f_c/f_s}, truncated after 200 samples.

(a) Show that Σ h[n] → 1 (unit DC gain) and that the filter's frequency response is H(f) = (1 − α)/(1 − α e^{−2πif/f_s}). Give |H| at 25, 50 and 83.3 Hz.
(b) Apply the filter to the 2 s synthetic record with `convolve_fft`, discard the first 256 samples of the output (why?), and read the three tone amplitudes before and after with a Hann window.
(c) Compare each measured attenuation ratio with |H(f₀)| from (a), and with the continuous first-order lag 1/√(1 + (f/f_c)²).
(d) Verify that `convolve_fft` and `numpy.convolve` agree to 10⁻⁹ and compare their run times for this case.

Answer form: a table of three tones (before, after, ratio, |H|, continuous), two numbers for timing.

## 6. (Judgment: specifying the data acquisition)

You must specify the accelerometer channel for a panel-vibration survey. The tones of interest lie between 20 and 400 Hz; the two closest tones you must separate are 0.5 Hz apart; the amplitude must be reported to ±5 %; the engine speed wanders by ±1 % during a run.

(a) Choose the sampling rate and the anti-aliasing filter corner, with reasons.
(b) Choose the record length and the window, with reasons, and state the amplitude bias your choice carries.
(c) Which of the four requirements conflicts with the engine-speed wander, and what would you do about it?
(d) A colleague proposes zero-padding each record to ten times its length "to get 0.05 Hz resolution". What does zero-padding actually do, and does it meet the 0.5 Hz separation requirement?

Answer form: four short paragraphs, each with at least one number.
