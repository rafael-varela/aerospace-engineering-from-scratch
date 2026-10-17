# P02.10  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown.

## 1. The three numbers of a record

**Derivation** ($s = 15$: f_s = 2150 Hz, N = 4096). (a) T = N/f_s = 1.905 s; Δf = f_s/N = **0.5249 Hz**; f_Nyquist = **1075 Hz**; bins = N/2 + 1 = **2049** (0 to 1075 Hz inclusive).

(b) 3.0 cos(2πf₁₀₀t + φ) = 1.5 e^{iφ} e^{2πi·100n/N} + 1.5 e^{−iφ} e^{−2πi·100n/N}. By orthogonality |X[100]| = |X[N − 100]| = 1.5N = **6144 με** (the units of the DFT are the units of x times a count), and `amplitude_spectrum` returns 2 × 6144/4096 = **3.0 με**, exactly.

(c) Rectangular window: main lobe 2Δf wide, two tones separable when Δf ≤ 0.1 Hz, so **T ≥ 10 s**. Hann: main lobe 4Δf, so **T ≥ 20 s**. Longer records buy resolution; the sample rate buys bandwidth. They are independent.

**Check.** Δf × N = f_s ✓; 2049 bins × 0.5249 Hz = 1075.5 Hz, the Nyquist bin at k = 2048 ✓.

**Common mistakes.** Taking Δf = f_s/2N; forgetting the factor 2 in (5) and reporting 1.5 με.

## 2. Orthogonality, Parseval, leakage, a DFT by hand

(a) Let r = e^{2πi(m−n)/N}. If m = n, every term is 1 and the sum is N. Otherwise r ≠ 1 (because |m − n| < N), rᴺ = e^{2πi(m−n)} = 1, and Σ_{k=0}^{N−1} rᵏ = (1 − rᴺ)/(1 − r) = 0. Now multiply (1) by e^{2πikm/N} and sum over k: Σ_k X[k] e^{2πikm/N} = Σ_n x[n] Σ_k e^{2πik(m−n)/N} = Σ_n x[n] N δ_{mn} = N x[m], which is (2).

(b) Σ_n |x[n]|² = Σ_n x[n] x̄[n] = (1/N²) Σ_n Σ_k Σ_l X[k] X̄[l] e^{2πi(k−l)n/N} = (1/N²) Σ_k Σ_l X[k] X̄[l] N δ_{kl} = (1/N) Σ_k |X[k]|².

(c) X[k₀] = Σ_{n=0}^{N−1} e^{2πi(k₀+δ)n/N} e^{−2πik₀n/N} = Σ_n e^{2πiδn/N} = (1 − e^{2πiδ})/(1 − e^{2πiδ/N}). With |1 − e^{iθ}| = 2|sin(θ/2)|, |X[k₀]| = |sin πδ|/|sin(πδ/N)|. For N ≫ 1, sin(πδ/N) ≈ πδ/N, so |X[k₀]| ≈ N sin(πδ)/(πδ). The factor is **1** at δ = 0, **0.9003** at δ = ¼, **0.6366 = 2/π** at δ = ½.

(d) w = e^{−2πi/4} = −i. X[k] = x[0] + x[2] w^{2k} = 1 − (−i)^{2k} = 1 − (−1)^k, so X = **[0, 2, 0, 2]**. The sequence is cos(2πn/4) = cos(πn/2) sampled at four points per cycle: a cosine of amplitude 1 on bin k = 1; the two bins (1 and N − 1 = 3) are the e^{+iθ} and e^{−iθ} halves of the cosine, each of magnitude N/2 = 2. `dft([1, 0, −1, 0])` returns [0, 2, 0, 2] ✓, and `amplitude_spectrum` would read 2 × 2/4 = 1.

**Common mistakes.** In (a), forgetting that the identity needs |m − n| < N; in (c), using sin(πδ/N) ≈ πδ/N at the numerator too; in (d), using w = +i (the inverse-transform sign).

## 3. Helicopter logged too slowly ($s = 15$: f_s = 55 Hz, Nyquist 27.5 Hz)

**Derivation.** (a) Main rotor 1P = 258/60 = **4.30 Hz**; blade passing 4 × 4.30 = **17.2 Hz**; tail rotor 1P = 1190/60 = **19.83 Hz**; tail blade passing 2 × 19.83 = **39.67 Hz**.

(b) 4.30, 17.2 and 19.83 Hz are below 27.5 Hz and appear where they are. 39.67 Hz is above Nyquist: r = 39.67 mod 55 = 39.67, alias = min(39.67, 55 − 39.67) = **15.33 Hz**. `alias_frequency(39.667, 55)` returns 15.333 ✓.

(c) The tail-rotor blade-passing tone, folded to 15.3 Hz, lands in the band the manual associates with a loose tail-rotor gearbox. The logged spectrum cannot tell a genuine 15.3 Hz tone from the alias: both are the same set of numbers. A maintainer would chase a fault that does not exist. The fix is an analogue anti-aliasing low-pass filter ahead of the sampler, with its corner below 27.5 Hz; software after the fact cannot separate the two.

(d) f_s ≥ 2 × 1.2 × 39.67 = **95.2 Hz**; the next standard rate (100 Hz) would do, with a filter corner near 48 Hz.

**Check.** At f_s = 100 Hz the 39.67 Hz tone is below Nyquist and `alias_frequency(39.667, 100)` returns 39.667 ✓.

**Common mistakes.** Folding about f_s instead of f_s/2 (giving 55 − 39.67 = 15.33 by luck here, but 2f_s − f for tones between f_s and 1.5f_s); forgetting the 1P tones.

## 4. The same tone, three window lengths ($s = 15$: f₀ = 82.1429 Hz, amplitude 2.0 m/s²)

| T [s] | Δf [Hz] | nearest bin [Hz] | δ | predicted 2 sin(πδ)/(πδ) | rectangular | Hann |
|---|---|---|---|---|---|---|
| 1 | 1.00 | 82.00 | 0.143 | 1.934 | 1.935 | 1.974 |
| 2 | 0.50 | 82.00 | 0.286 | 1.742 | 1.743 | 1.897 |
| 4 | 0.25 | 82.25 | 0.429 | 1.448 | 1.449 | 1.774 |

(a)–(c) as tabulated (from `amplitude_spectrum` on 2 cos(2πf₀t) with no noise). The prediction holds to 0.1 %. The Hann readings spread 1.77–1.97 (10 %) while the rectangular ones spread 1.45–1.94 (25 %): the Hann window varies less with record length.

(d) The DFT analyses the periodic extension of the record, and the jump at the window edge depends on how many cycles fit in T: the bin grid moved, the tone did not. An amplitude from a DFT is meaningless without the window and the record length (equivalently δ), and a report should give both, or give the Hann reading with its 0–15 % bias.

**Check.** The longer record here gives the *worse* rectangular reading (δ grew from 0.14 to 0.43): resolution and amplitude accuracy are different things.

**Common mistakes.** Reporting δ from the bin below instead of the nearest bin; applying the Hann window without the gain correction (readings halve).

## 5. Filtering through the convolution theorem (f_c = 40 Hz, f_s = 1024 Hz, α = 0.78236)

**Derivation.** (a) Σ_{n≥0}(1 − α)αⁿ = (1 − α)/(1 − α) = 1; the truncation at 200 samples loses α²⁰⁰ ≈ 10⁻²¹, so `h.sum()` = 1.0000. The frequency response is the same geometric series with ratio α e^{−2πif/f_s}: H(f) = (1 − α)/(1 − α e^{−2πif/f_s}). Magnitudes: **0.8488** at 25 Hz, **0.6271** at 50 Hz, **0.4376** at 83.3 Hz.

(b) The first 256 output samples contain the start-up transient (time constant 1/(1 − α) ≈ 4.1 samples; by 256 samples α²⁵⁶ ≈ 5 × 10⁻²⁸): the filter output has not yet "forgotten" that the input began abruptly, and a window over a transient leaks. Hann readings of `a[256:]` and `y[256:]`:

| tone [Hz] | before [m/s²] | after [m/s²] | ratio | \|H(f₀)\| | 1/√(1 + (f/f_c)²) |
|---|---|---|---|---|---|
| 25 | 3.825 | 3.246 | 0.8488 | 0.8488 | 0.848 |
| 50 | 1.256 | 0.787 | 0.6271 | 0.6271 | 0.625 |
| 83.3 | 1.922 | 0.841 | 0.4376 | 0.4376 | 0.433 |

(c) The measured ratios equal |H(f₀)| to four figures (the same leakage affects numerator and denominator). The continuous lag agrees to 1 %, the difference being the discretisation of the pole (α = e^{−2πf_c/f_s} is the exact sampling of a continuous first-order lag at zero-order hold, which matches at DC and departs slowly above f_c).

(d) max |`convolve_fft` − `np.convolve`| = 2 × 10⁻¹⁵. Timing on the author's machine: `convolve_fft` 70 ms (three Python-recursion FFTs of length 4096), `np.convolve` 0.09 ms (C, direct, 2048 × 200 = 4 × 10⁵ multiply-adds). The FFT route wins only when N_x N_h is far larger than N log N *and* the FFT is a compiled one; with `scipy.signal.fftconvolve` the crossover is near N_h ≈ 50–100 taps.

**Common mistakes.** Not discarding the transient (the 25 Hz reading drops by several per cent); comparing a windowed "after" with an unwindowed "before".

## 6. Specifying the data acquisition

(a) **Sampling rate and filter.** Tones to 400 Hz: Nyquist must exceed 400 Hz with margin, so f_s ≥ 1000 Hz (Nyquist 500 Hz); pick 1024 Hz or 2048 Hz. The anti-aliasing corner goes at about 400–450 Hz with enough roll-off that whatever is at 600 Hz (which would fold to 424 Hz at f_s = 1024 Hz) is attenuated below the ±5 % amplitude requirement; a 4-pole filter at 420 Hz gives 24 dB per octave, which is marginal at 1024 Hz and comfortable at 2048 Hz. Choose 2048 Hz.

(b) **Record length and window.** 0.5 Hz separation with a Hann window (main lobe 4Δf) needs Δf ≤ 0.125 Hz, so T ≥ 8 s; with a rectangular window T ≥ 4 s. Amplitude to ±5 % rules out the rectangular window (36 % worst case) and, strictly, the Hann (15 % worst case). Use Hann for separation and report the bias, or use a flat-top window for amplitude at the cost of resolution. One record length cannot optimise both; two analyses of the same record can.

(c) **Drift.** ±1 % of a 400 Hz tone is ±4 Hz, which over an 8 s record (Δf = 0.125 Hz) smears the tone across 60 bins: the peak reading collapses and the 0.5 Hz separation is lost. The conflict is between resolution and stationarity. Remedies: order tracking (resample the record in engine-angle rather than time, so the tone is fixed in "orders"), or short records (1 s, Δf = 1 Hz) accepting the loss of the 0.5 Hz requirement, or stabilising the engine speed for the run.

(d) **Zero-padding.** Appending zeros to a record of length T evaluates the same windowed spectrum at a finer frequency grid: the main lobe is still 2/T (rectangular) or 4/T (Hann) wide in hertz, so two tones closer than that still merge. Zero-padding helps *locate* a single peak between bins (interpolation) and makes plots smoother; it does not meet a separation requirement. Ten-fold padding gives 0.05 Hz *grid spacing*, not 0.05 Hz *resolution*.

## Check-yourself answers

1. X[k] = Σ_{n=0}^{N−1} x[n] e^{−2πikn/N}; x[n] = (1/N) Σ_k X[k] e^{+2πikn/N}. The 1/N is on the inverse (the `scipy.fft` convention).
2. Δf = 1/T = 0.5 Hz; f_Nyquist = f_s/2 = 512 Hz; N/2 + 1 = 1025 bins.
3. Because x̄[n] = x[n], so X[N − k] = X̄[k]: the tone's energy is split equally between bins k and N − k, each with |X| = Na/2, which is why the single-sided amplitude doubles |X[k]|/N.
4. Rectangular: 3 × 2/π = 1.91 m/s². Hann: about 15 % low, roughly 2.55 m/s².
5. A circular shift by m multiplies X[k] by e^{−2πikm/N}, a pure phase of unit magnitude; |X[k]| is unchanged, so the magnitude spectrum does not depend on the start time (for a stationary signal).
6. r = 410 mod 300 = 110; min(110, 190) = 110 Hz. An anti-aliasing low-pass filter with its corner below 150 Hz, before the sampler.
7. The rectangle sum has error ½Δt·x(0⁺) from the jump at t = 0, a constant of order Δt. Halving the first (and last) sample, the trapezoidal rule, cancels it and leaves an O(Δt²) error.
8. Σ_n |x[n]|² = (1/N)Σ_k |X[k]|²; dividing by N and substituting (5): RMS = √(A[0]² + ½Σ_{k>0} A[k]²), with the Nyquist bin counted once (A²) rather than halved.
