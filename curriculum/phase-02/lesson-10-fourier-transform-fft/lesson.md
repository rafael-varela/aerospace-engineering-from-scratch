# P02.10  Fourier Transform & the FFT

> A record of N samples is N numbers in time or N numbers in frequency; the DFT is the change of basis between them, and the FFT is the only reason it is affordable.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the discrete Fourier transform from the Fourier series of a sampled, periodically extended record, state its inverse, and prove Parseval's theorem and the shift theorem from the definition. 2. Implement and test, from scratch, the O(N²) DFT as a matrix product and the radix-2 Cooley–Tukey FFT, against `scipy.fft` as the oracle. 3. Convert a DFT into a single-sided amplitude spectrum in physical units, measure how spectral leakage distorts an off-bin tone, and decide when a Hann window is worth its cost. 4. Predict the alias frequency of an undersampled tone, and estimate a continuous Fourier transform from samples with a known, measured convergence order. |
| **Prerequisites** | `P02.09` (Fourier series, orthogonality of complex exponentials), `P02.08` (the Laplace transform; the Fourier transform is its value on the imaginary axis), `P02.07` (complex exponentials, phasors), `P01.05` (geometric series), `P00.03`, `P00.07`. Skill check: you can sum a finite geometric series and evaluate e^{−2πi k} for integer k. |
| **Engineering question** | A (synthetic) accelerometer record from a light aircraft's skin panel is 2 s long at 1024 samples per second: which tones does it contain, at what amplitude, how badly does a 2 s window mis-read a tone that does not sit on a frequency bin, and what would the same record say if the data logger had been set to 120 Hz? |
| **Assumptions and limits** | The record is stationary (the tones do not change over the 2 s); sampling is uniform and the clock is exact; the DFT treats the record as one period of a periodic signal, so a tone that does not complete a whole number of cycles in the window leaks into neighbouring bins; any content above half the sampling rate is folded, not removed; the from-scratch FFT handles only power-of-two lengths. The record itself is synthetic: three tones plus seeded Gaussian noise, not a flight measurement. |
| **Analytical work** | DFT from the Fourier series of the periodic extension; inverse by orthogonality of the sampled exponentials (finite geometric series); linearity, conjugate symmetry, shift theorem, Parseval; the radix-2 split and its operation count; the leakage factor for a tone half-way between bins (2/π); the alias formula; the Riemann-sum estimate of the continuous transform and its first- and second-order error; the convolution theorem. |
| **Computational deliverable** | `starter/main.py`: `dft_matrix`, `dft`, `idft`, `fft_radix2`, `fft_any`, `fft_frequencies`, `hann_window`, `amplitude_spectrum`, `parseval`, `ft_exp_decay`, `ft_from_samples`, `alias_frequency`, `synthetic_vibration`, `dominant_peaks`, `rms_from_spectrum`, `convolve_fft`. Tests prove: the DFT matrix satisfies FᴴF = NI; DFT and FFT match `scipy.fft` to 10⁻⁹ for every tested length; the inverse round-trips; linearity, the shift theorem, conjugate symmetry and Parseval hold; an on-bin sinusoid is recovered to 10⁻¹⁰ in amplitude; a mid-bin tone reads 2/π low with a rectangular window and 10–20 % low with a Hann window; the undersampled record peaks where the alias formula says; the sampled transform of e^{−at} converges to 1/(a + 2πif) at first order (rectangle) and second order (trapezoid); FFT convolution equals `numpy.convolve`. |
| **Practical activity** | **data-based investigation** on a labelled synthetic record: identify its tones, compare rectangular and Hann windows, undersample it and verify the alias prediction, and measure the convergence of the sampled transform. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `vibration_fft.png`, the spectrum table (true tones against rectangular and Hann readings), the alias table, the convergence table for e^{−at}, and a short note on which window you would use for the panel-mode amplitude and why. |
| **Assessment** | `assessments/phase-02/grading-rubric.md`; competencies 1–4. Pass: the derivation of the inverse DFT and of Parseval's theorem is complete; all public tests pass; the hand values for resolution, Nyquist frequency, leakage factor and alias frequencies agree with the code; the window note states a number for the amplitude error of each window. |
| **References** | `@mit-ocw-6-003-signals-and-systems` Lecture 16 "Fourier Transform", Lecture 19 "Relations Among Fourier Representations", Lecture 21 "Sampling" (required-free); `@mit-ocw-18-06-linear-algebra` Lecture 26 "Complex matrices; fast fourier transform" (required-free: the DFT matrix and the FFT factorisation); `@mit-ocw-6-003-signals-and-systems` Lecture 17 "Discrete-Time (DT) Frequency Representations", Lecture 18 "Discrete-Time (DT) Fourier Representations", Lecture 20 "Applications of Fourier Transforms", Lecture 22 "Sampling and Quantization" (supplementary-free); `@mit-ocw-18-03-differential-equations` notes 20 "Fourier series", 21 "Operations on Fourier series" (prerequisite, see `P02.09`). |
| **Extensions** | ★ Implement an iterative (bit-reversal) radix-2 FFT and compare its speed with the recursive one. ★ Add zero-padding and parabolic peak interpolation to `dominant_peaks` and show the 83.3 Hz tone can be located to 0.05 Hz from the same 2 s record. ★ Estimate a power spectral density by averaging Hann-windowed segments (Welch's method) and compare with `scipy.signal.welch`; this is the tool `P02.19` needs. |
| **Time** | unmeasured |

## 1. The engineering problem

An accelerometer bonded to a skin panel behind the propeller of a light aircraft records acceleration for 2 s at 1024 samples per second. The flight-test engineer wants three things from the record: the frequency of each tone in it, the amplitude of each tone in m/s², and a statement of what the record cannot tell her. The propeller turns at 1500 rpm, so a once-per-revolution (**1P**) tone at 25 Hz and its first harmonic (**2P**) at 50 Hz are expected; the panel has a first mode near 83 Hz, which may be excited by the engine, by the slipstream, or by nothing at all. If the 83 Hz amplitude is large the panel is being worked near resonance and will crack at the rivet line; if the instrument *reads* it large because of a processing artefact, a sound panel gets an unnecessary stiffener. Telling those two cases apart is what this lesson is about.

The record in this lesson is **synthetic**: `synthetic_vibration` adds three cosines of known frequency, amplitude and phase to seeded Gaussian noise, so that every reading can be compared with the truth. Nothing in it was measured on an aircraft. The numbers (25 Hz at 4.0 m/s², 50 Hz at 1.5 m/s², 83.3 Hz at 2.0 m/s², noise 0.5 m/s² RMS) are chosen to be plausible, and the 83.3 Hz tone is placed deliberately *between* the frequency bins the analysis can resolve, because that is the case where the naive answer is wrong by a quarter.

Success criteria, written down first: the from-scratch DFT and FFT agree with `scipy.fft` to 10⁻⁹; the on-bin tones are read to within the noise (a few per cent); the off-bin tone is read 20–30 % low with a rectangular window and about 10 % low with a Hann window, and both errors are predicted before the run; sampling the same signal at 120 Hz folds the 83.3 Hz tone to the frequency the alias formula predicts, to within one bin; and the estimate of the continuous transform of e^{−at} converges at the order the derivation predicts.

## 2. Understand and derive

### 2.1 From the Fourier series to the DFT

`P02.09` showed that a periodic signal with period T is a sum of complex exponentials, x(t) = Σ_k c_k e^{2πikt/T}, with c_k = (1/T)∫₀ᵀ x(t) e^{−2πikt/T} dt. Now sample x at N equally spaced times t_n = n Δt, n = 0 … N−1, with N Δt = T, and replace the integral by a sum over the samples (the rectangle rule):

c_k ≈ (1/T) Σ_n x[n] e^{−2πik n Δt/T} Δt = (1/N) Σ_{n=0}^{N−1} x[n] e^{−2πi k n/N}.

The sum, without the 1/N, is the **discrete Fourier transform**:

**X[k] = Σ_{n=0}^{N−1} x[n] e^{−2πi k n/N}, k = 0 … N−1.**   (1)

Two things happened in that step. First, the exponential e^{−2πi k n/N} is unchanged when k is replaced by k + N, so only N distinct values of k exist: the DFT of N samples is N numbers, not infinitely many. Second, the DFT knows nothing about x(t) between the samples or outside the window: it is the Fourier series of the **periodic extension** of the sampled record, the signal you get by repeating the N samples forever. Every artefact in this lesson (leakage, aliasing) is a consequence of that fact, not a bug in the algorithm (`@mit-ocw-6-003-signals-and-systems` Lecture 19).

### 2.2 The inverse, by orthogonality

Multiply (1) by e^{+2πi k m/N} and sum over k:

Σ_k X[k] e^{2πikm/N} = Σ_n x[n] Σ_{k=0}^{N−1} e^{2πi k(m−n)/N}.

The inner sum is a finite geometric series with ratio r = e^{2πi(m−n)/N} (`P01.05`). If m = n every term is 1 and the sum is N. Otherwise r ≠ 1 but r^N = e^{2πi(m−n)} = 1, so the sum is (1 − r^N)/(1 − r) = 0. Hence Σ_k e^{2πik(m−n)/N} = N δ_{mn}, the sampled exponentials are **orthogonal**, and

**x[n] = (1/N) Σ_{k=0}^{N−1} X[k] e^{+2πi k n/N}.**   (2)

In matrix language (`@mit-ocw-18-06-linear-algebra` Lecture 26): with w = e^{−2πi/N}, the **DFT matrix** F has entries F_{kn} = w^{kn}. It is symmetric, and the orthogonality just proved says FᴴF = N I, so F⁻¹ = F̄/N: the inverse transform is the conjugate matrix divided by N. F/√N is unitary, which is why the transform preserves lengths (2.3) and why rounding error does not grow through it. `dft` is literally `F @ x` and `idft` is `conj(F) @ X / N`; the tests check FᴴF = NI directly.

### 2.3 Four properties you will use every time

**Linearity.** DFT(ax + by) = a DFT(x) + b DFT(y), immediate from (1).

**Conjugate symmetry.** If x is real, X[N−k] = Σ x[n] e^{−2πi(N−k)n/N} = Σ x[n] e^{+2πikn/N} = X̄[k]. Half the bins repeat the other half, which is why a real record of N samples gives a one-sided spectrum of N/2 + 1 bins (2.5).

**Shift theorem.** A circular shift of the samples by m, x[n − m], multiplies X[k] by e^{−2πikm/N}: substitute n' = n − m in (1) and use the periodicity of the exponential. Magnitudes are unchanged; only phases move. (A shift in *time* is a linear phase in *frequency*; the mirror statement, that a modulation e^{2πik₀n/N} in time is a shift in frequency, is why a tone appears as a spike at k₀.)

**Parseval's theorem.** Using (2) and orthogonality, Σ_n |x[n]|² = (1/N) Σ_n Σ_k Σ_l X[k] X̄[l] e^{2πi(k−l)n/N} = (1/N) Σ_k |X[k]|²:

**Σ_n |x[n]|² = (1/N) Σ_k |X[k]|².**   (3)

Energy is the same in both domains. For the vibration record Σ|x|² = 23 151.61 and (1/N)Σ|X|² = 23 151.61; `parseval` returns both and the test demands equality to 10⁻¹². The practical use is a check that your amplitude scaling (2.5) is right: the RMS of the record must equal the RMS reconstructed from the spectrum, 3.3622 m/s² both ways here.

### 2.4 The fast Fourier transform

Equation (1) costs N complex multiplications per bin, N² in all: for N = 2048 that is 4.2 million, for a 2²⁰-sample record (a 17-minute flight-test channel at 1 kHz) 1.1 × 10¹². The **Cooley–Tukey** observation is that for even N the sum splits into even and odd samples:

X[k] = Σ_{m=0}^{N/2−1} x[2m] e^{−2πi k m/(N/2)} + e^{−2πik/N} Σ_{m=0}^{N/2−1} x[2m+1] e^{−2πi k m/(N/2)} = E[k] + w^k O[k],

where E and O are DFTs of length N/2 and w^k = e^{−2πik/N} is the **twiddle factor**. Because E and O are periodic with period N/2, the bins above N/2 come free:

**X[k] = E[k] + w^k O[k], X[k + N/2] = E[k] − w^k O[k], k = 0 … N/2 − 1.**   (4)

Applied recursively (N must be a power of two for the pure radix-2 form), the cost obeys C(N) = 2C(N/2) + N/2 multiplications, whose solution is C(N) = (N/2) log₂ N. For N = 2048 that is 11 264 against 4 194 304: 372 times fewer; for N = 2²⁰ it is 10⁵ times fewer. Lecture 26 of `@mit-ocw-18-06-linear-algebra` writes the same thing as a matrix factorisation, F_N = [I D; I −D] [F_{N/2} 0; 0 F_{N/2}] P with D = diag(w^k) and P the even–odd permutation. The algorithm changes nothing about the answer: `fft_radix2` and `dft` agree to 1.2 × 10⁻¹² and 9.8 × 10⁻¹⁰ respectively against `scipy.fft` on the vibration record (the O(N²) form is *less* accurate because it sums 2048 rounded products per bin). On the author's machine the DFT took about 0.6 s, the Python recursive FFT about 15 ms and `scipy.fft` about 80 µs: the library is written in C and uses a mixed-radix plan, which is why part 3 calls it in practice once the from-scratch version has proved what it computes.

### 2.5 Frequencies, resolution, and amplitude in physical units

Bin k of a record of N samples spaced Δt corresponds to frequency f_k = k/(N Δt) = k/T. Three numbers follow (`@mit-ocw-6-003-signals-and-systems` Lecture 21):

- **Resolution** Δf = 1/T = 1/(N Δt): 0.5 Hz for our 2 s record. Twice the record length gives twice the resolution; the sampling rate does not enter.
- **Nyquist frequency** f_s/2 = 1/(2 Δt): 512 Hz. Bins above N/2 are the negative frequencies (conjugate symmetry), which `fft_frequencies` lists in the `scipy.fft.fftfreq` order: 0, Δf, …, (N/2 − 1)Δf, −N/2 Δf, …, −Δf.
- **Amplitude.** A cosine a cos(2πf_k t + φ) that sits exactly on bin k is ½a e^{iφ} e^{2πikn/N} + ½a e^{−iφ} e^{−2πikn/N}; by orthogonality X[k] = (N/2) a e^{iφ} and X[N−k] is its conjugate. So the **single-sided amplitude spectrum** is

**A[k] = 2|X[k]|/N for 0 < k < N/2, A[0] = |X[0]|/N,**   (5)

and the Nyquist bin, which has no mirror image, is |X[N/2]|/N. With (5), a 4.0 m/s² tone reads 4.0 m/s², and the mean of the record reads as A[0]. The RMS of the record is √(A[0]² + ½Σ_{k>0} A[k]²) (plus ½A[N/2]² for the unpaired Nyquist bin), which is Parseval again in engineering clothes.

### 2.6 Leakage: what happens between the bins

The 83.3 Hz tone completes 166.6 cycles in 2 s. Its periodic extension therefore has a jump at the window edge, and a jump has a spectrum that falls off only as 1/f. Concretely, for a tone at k₀ + δ bins (0 ≤ δ < 1) the DFT of the complex exponential is a geometric series again, Σ_n e^{2πi(k₀+δ−k)n/N} = (1 − e^{2πiδ})/(1 − e^{2πi(k₀+δ−k)/N}), with magnitude ≈ N |sin πδ|/(π|k₀ + δ − k|) for bins far from the tone: the **sidelobes** decay as 1/|k − k₀|. The peak bin reads |sin πδ|/(π δ) times the true amplitude: 1 at δ = 0, **2/π = 0.637** at δ = ½. The 83.3 Hz tone is at δ = 0.6 from bin 83.0 (0.4 from 83.5), so the predicted reading at 83.5 Hz is 2.0 × sin(0.4π)/(0.4π) = 1.514 m/s². The run gives **1.523 m/s²**: a 24 % under-read, and the 0.01 difference is noise.

A **window** tapers the record to zero at both ends so that the periodic extension has no jump. The Hann window w[n] = ½(1 − cos 2πn/N) has mean (coherent gain) exactly ½, so the windowed spectrum is divided by ½ to restore on-bin amplitudes; its sidelobes fall as 1/f³ instead of 1/f, and the worst-case loss at δ = ½ is about 15 % instead of 36 %. The price is a main lobe twice as wide (two tones closer than 2Δf merge) and a noise bandwidth 1.5 times larger. With the Hann window the 83.3 Hz tone reads **1.798 m/s²**, 10 % low; the on-bin tones read 3.977 and 1.480 against 4.0 and 1.5, the deficit being noise and the leakage of the 83.3 Hz tone's own sidelobes. For a mid-bin test tone the code gives 1.275 (rectangular, predicted 2 × 0.637 = 1.273) and 1.698 (Hann) for a true 2.0.

Which window? For *finding* a tone's frequency, either; for *reading its amplitude*, the Hann window unless you know the tone is on a bin; for *separating two close tones*, the rectangular window. The engineering decision in part 1 (is the panel near resonance?) needs the amplitude, so the Hann reading of 1.8 m/s², with its known 0–15 % bias, is the one to report, together with the bias.

### 2.7 Aliasing: what sampling throws away

The sampled exponential e^{2πi f n Δt} is identical to e^{2πi (f + m f_s) n Δt} for every integer m, because the difference in the exponent is 2πi m n, a whole number of turns. A sampled record therefore cannot distinguish f from f ± f_s, f ± 2f_s, …, and for a real signal cannot distinguish f from −f either. Every frequency folds into [0, f_s/2]:

**f_alias = min(f mod f_s, f_s − f mod f_s).**   (6)

This is the **sampling theorem** read backwards (`@mit-ocw-6-003-signals-and-systems` Lecture 21): a signal is recoverable from its samples only if it contains nothing above f_s/2. Set the logger to 120 Hz and the 83.3 Hz panel mode appears at |83.3 − 120| = 36.7 Hz, indistinguishable from a real 36.7 Hz tone; the 25 and 50 Hz tones, below the 60 Hz Nyquist frequency, stay put. The run confirms it: peaks at 25.0, 36.5 and 50.0 Hz with amplitudes 4.002, 1.506 and 1.477 m/s² (the 36.5 is the nearest 0.5 Hz bin to 36.7, and the 1.5 is the same leakage loss as before). No processing after the fact can undo this: the only defence is an analogue low-pass (**anti-aliasing**) filter before the sampler, with its corner below f_s/2, and that is why a data logger's bandwidth specification is always lower than half its sample rate.

### 2.8 The continuous transform from samples

For an aperiodic signal, let the period T → ∞ in the Fourier series: the harmonics k/T crowd together into a continuous variable f, the coefficients c_k T tend to a function, and the series becomes an integral (`@mit-ocw-6-003-signals-and-systems` Lecture 16):

**X(f) = ∫ x(t) e^{−2πift} dt, x(t) = ∫ X(f) e^{2πift} df.**   (7)

Compare with the Laplace transform of `P02.08`, X(s) = ∫₀^∞ x(t) e^{−st} dt: for a signal that is zero before t = 0 and decays, the Fourier transform is the Laplace transform evaluated on the imaginary axis, s = 2πif. So the causal exponential x(t) = e^{−at} has X(f) = 1/(a + 2πif): the frequency response of a first-order lag, with magnitude 1/a at f = 0 and 1/(a√2) (−3 dB) at f = a/2π.

Sampling turns (7) into a rectangle-rule sum, X(f_k) ≈ Δt Σ_n x[n] e^{−2πi f_k n Δt} = Δt · DFT[k], which is `ft_from_samples`. Two errors enter. **Truncation**: the record must be long enough that x has decayed (e^{−100} at t = 2 s for a = 50 s⁻¹: negligible). **Quadrature**: the rectangle rule has error of order Δt, and for a signal with a jump at t = 0 that error is ½Δt·x(0), a constant offset that becomes a large *relative* error where |X| is small. The trapezoidal rule (half weight on the first and last samples) removes it and leaves an error of order Δt². The run, for a = 50 s⁻¹ and f ≤ 100 Hz, gives maximum relative errors of 0.312, 0.155 and 0.077 for the rectangle sum at Δt = 1/1024, 1/2048 and 1/4096 s (ratio 2.0: first order) and 0.0318, 0.0079 and 0.0020 for the trapezoid (ratio 4.0: second order). Measuring the order is the test that the error is understood rather than merely small.

### 2.9 The convolution theorem

The response of a linear time-invariant system is the convolution y = x ∗ h of the input with the impulse response (`P02.08`). In the frequency domain convolution is multiplication: Y(f) = X(f) H(f) (`@mit-ocw-6-003-signals-and-systems` Lecture 20). For sampled records the DFT version holds for *circular* convolution; to get the ordinary (linear) convolution of lengths N_x and N_h, zero-pad both to at least N_x + N_h − 1 so that the wrap-around lands in the padding. `convolve_fft` does exactly that and matches `numpy.convolve` to rounding, at O(N log N) instead of O(N_x N_h): for a 10⁶-sample record and a 10⁴-tap filter that is the difference between a second and an hour.

### 2.10 Where the model stops being valid

- **Stationarity.** A tone whose frequency drifts during the window (an engine spooling up) smears across bins; the remedy is short windows and a spectrogram, at the cost of resolution.
- **Uniform sampling.** Logger jitter adds phase noise; dropped samples invalidate (1) entirely.
- **Periodic extension.** Leakage (2.6) is not noise and does not average away.
- **Nothing above f_s/2.** Aliasing (2.7) is irreversible. The anti-aliasing filter is part of the instrument, not of the analysis.
- **Power of two.** `fft_radix2` refuses other lengths; `fft_any` falls back to the O(N²) DFT, and `scipy.fft` handles any length with mixed radices.

## 3. Build the computational model

Open `starter/main.py`. `vibration_components` and `main` are given; everything else is yours.

### Step 1: the DFT matrix and the O(N²) transform

```python
def dft_matrix(n):
    k = np.arange(n)
    return np.exp(-2j * np.pi * np.outer(k, k) / n)                       # F[k, m] = w^{km}

def dft(x):
    x = np.asarray(x, dtype=complex)
    return dft_matrix(x.size) @ x                                         # (1)

def idft(X):
    X = np.asarray(X, dtype=complex)
    return np.conj(dft_matrix(X.size)) @ X / X.size                       # (2)
```

Reject n < 1 and non-1D input with `ValueError`. The matrix form costs N² memory as well as time; for N = 2048 that is 64 MB of complex numbers, which is the last size at which it is sensible.

### Step 2: the radix-2 FFT

```python
def fft_radix2(x):
    x = np.asarray(x, dtype=complex)
    n = x.size
    if n < 1 or (n & (n - 1)) != 0:
        raise ValueError("fft_radix2 needs a length that is a power of two")
    if n == 1:
        return x.copy()
    even, odd = fft_radix2(x[0::2]), fft_radix2(x[1::2])
    twiddle = np.exp(-2j * np.pi * np.arange(n // 2) / n) * odd
    return np.concatenate([even + twiddle, even - twiddle])              # (4)
```

`n & (n - 1) == 0` is the power-of-two test. `fft_any` dispatches to this when it can and to `dft` otherwise; every later function calls `fft_any`.

### Step 3: frequencies, window, amplitude spectrum

`fft_frequencies(n, dt)` returns k/(n dt) with k = 0 … ⌈n/2⌉ − 1 followed by the negative half, matching `scipy.fft.fftfreq`. `hann_window(n)` is ½(1 − cos 2πm/n). `amplitude_spectrum(x, dt, window)` multiplies by the window, transforms with `fft_any`, keeps bins 0 … n/2, applies (5), and divides by the coherent gain (1 or ½); the Nyquist bin, present only for even n, is *not* doubled. `parseval` returns the two sides of (3); `rms_from_spectrum(A, n)` inverts (5) into an RMS.

### Step 4: the sampled continuous transform and the alias formula

`ft_exp_decay(f, a)` is 1/(a + 2πif) and must raise for a ≤ 0. `ft_from_samples(x, dt, trapezoid)` halves the first and last samples when asked and returns Δt times the transform, with frequencies in `fftfreq` order. `alias_frequency(f, fs)` is (6): `r = f mod fs; min(r, fs − r)`.

### Step 5: the record, its peaks, and convolution

`synthetic_vibration(fs, duration, seed, noise_std)` builds t = n/f_s, sums the tones of `vibration_components`, and adds `noise_std * rng.standard_normal(n)` from `default_rng(seed)`. `dominant_peaks(f, A, n_peaks, min_separation)` finds local maxima (greater than the left neighbour, not less than the right), sorts them by amplitude, and skips any within `min_separation` of a stronger peak already chosen: that is what keeps leakage sidelobes out of the tone list. `convolve_fft` zero-pads to the next power of two ≥ N_x + N_h − 1, multiplies the transforms, inverts (you can invert with the forward FFT: ifft(Y) = conj(fft(conj(Y)))/N) and trims.

Running `main()` prints

```
SYNTHETIC record: fs = 1024 Hz, N = 2048, duration 2.0 s, resolution 0.500 Hz, Nyquist 512 Hz
max |DFT - scipy| = 9.80e-10, max |FFT - scipy| = 1.22e-12
Parseval: sum|x|^2 = 23151.6076, (1/N) sum|X|^2 = 23151.6076
true tones (frequency Hz, amplitude m/s^2): [(25.0, 4.0), (50.0, 1.5), (83.3, 2.0)]
rectangular window peaks: [(25.0, 3.983), (83.5, 1.523), (50.0, 1.482)]
Hann window peaks       : [(25.0, 3.977), (83.5, 1.798), (50.0, 1.48)]
RMS from time record 3.3622 m/s^2, from rectangular spectrum 3.3622 m/s^2
sampled at 120 Hz (Nyquist 60 Hz): peaks [(25.0, 4.002), (36.5, 1.506), (50.0, 1.477)]
predicted alias of each tone: [(25.0, 25.0), (50.0, 50.0), (83.3, 36.7)]
```

followed by the convergence table of 2.8 and the figure.

### Verify

- **Analytic cases:** FᴴF = NI for N = 1, 4, 7, 16; an on-bin cosine of amplitude a and phase φ returns A[k] = a to 10⁻¹⁰ with every other bin below 10⁻⁹; the mean of the record appears in A[0]; a rectangular pulse of width T transforms to T sinc(fT) with its first zero at 1/T; the transform of e^{−at} at f = 0 is 1/a and its magnitude at f = a/2π is 1/(a√2).
- **Limiting and failure cases:** lengths 1 and 2 transform correctly; non-power-of-two lengths make `fft_radix2` raise and `fft_any` fall back; a tone at exactly f_s stays at 0 Hz and one at f_s/2 stays at f_s/2; an unknown window name raises.
- **Conservation:** Parseval to 10⁻¹² for both DFT and FFT; the RMS from the spectrum equals the time-domain RMS to 10⁻¹⁰.
- **Symmetry and structure:** linearity, conjugate symmetry of a real record, the shift theorem, convolution with a delayed impulse is a delay.
- **Convergence:** the sampled transform error halves per halving of Δt for the rectangle sum and quarters for the trapezoid, to within 5 % and 10 % respectively.
- **Agreement with the oracle:** `scipy.fft.fft`, `fftfreq` and `numpy.convolve` to 10⁻⁹.

### Compare

`scipy.fft.fft` uses the pocketfft library: mixed-radix (2, 3, 5, …) butterflies chosen by a plan, with an O(N log N) Bluestein fallback for prime lengths, all in C. It gives the same numbers as (1) to rounding, in about a thousandth of the time of the Python recursion and a ten-thousandth of the matrix product. Call it. What writing our own bought is the knowledge of what the numbers *are*: a change of basis to sampled exponentials with the conventions of 2.1 and 2.5, so that a factor N, a factor 2, a missing window gain or a folded frequency is recognised as an error in the analysis rather than accepted as a property of the panel.

## 4. Experiment at home, or investigate real data

**Mode: data-based investigation** on the labelled synthetic record. Four studies; write the prediction before each run.

1. **Spectrum table.** For the 2 s record at 1024 Hz, read the three tones with the rectangular and the Hann window. Predict the rectangular reading of the 83.3 Hz tone from δ = 0.4: 2.0 × sin(0.4π)/(0.4π). Record `tone, true_amp, rect_amp, hann_amp`.
2. **Window length.** Repeat with 1 s and 4 s records (`synthetic_vibration(1024, 1.0)`, `(1024, 4.0)`). Predict first which readings change: the resolution becomes 1 Hz and 0.25 Hz, so the 83.3 Hz tone moves to δ = 0.3 and δ = 0.2 bins respectively; the on-bin tones do not move.
3. **Aliasing.** Sample at 120, 150 and 200 Hz (noise off) and tabulate the three peak frequencies against (6). At which of these rates does the 50 Hz tone alias? (At 100 Hz it would sit *on* the Nyquist frequency; avoid that rate and say why.)
4. **Convergence.** For e^{−50t}, record the maximum relative error up to 100 Hz against Δt = 1/1024, 1/2048, 1/4096 for both sums and fit the slopes on log-log axes.

Raw data format: four CSV tables with a one-line header naming the fixed parameters; `vibration_fft.png` with the four panels of `main()`.

What you should find: the rectangular reading of the off-bin tone depends on δ, not on the tone, so changing the window length changes the error even though nothing physical changed (1.51 at δ = 0.4, 1.72 at δ = 0.3, 1.87 at δ = 0.2, from 2 sin(πδ)/(πδ); the runs give 1.52, 1.75 and 1.88); the Hann reading is less sensitive to δ; every alias lands within one bin of (6); the convergence slopes are 1.0 and 2.0.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| resolution, Nyquist | 0.5 Hz, 512 Hz | 0.500 Hz, 512 Hz | 0 | definitions |
| FᴴF − NI | 0 | < 10⁻¹⁰ | rounding | — |
| 25 Hz amplitude, rectangular | 4.0 | 3.983 | 0.4 % | noise (0.5 m/s² RMS spread over 1025 bins gives about 0.02 per bin) |
| 83.3 Hz amplitude, rectangular | 1.514 (δ = 0.4) | 1.523 | 0.6 % | noise |
| 83.3 Hz amplitude, Hann | 1.7–2.0 (0–15 % low) | 1.798 | within the predicted band | Hann main-lobe shape at δ = 0.4 |
| mid-bin test tone, rectangular | 2 × 2/π = 1.273 | 1.275 | 0.2 % | finite N in the sinc approximation |
| RMS | 3.3622 (time) | 3.3622 (spectrum) | 10⁻¹⁰ | Parseval |
| alias of 83.3 Hz at f_s = 120 Hz | 36.7 Hz | 36.5 Hz | 0.2 Hz | nearest 0.5 Hz bin |
| convergence order, rectangle / trapezoid | 1 / 2 | 2.0 / 4.0 error ratios | 0 | — |

Every discrepancy is named, and all are below the uncertainty a real accelerometer would add (a typical piezoelectric sensor is specified to ±5 % in amplitude across its band). The one number that would mislead an engineer is the rectangular-window 1.52 m/s² for the panel mode: it is 24 % low, the error is systematic, and no amount of averaging removes it. A Hann window with its known bias, or a longer record that moves the tone onto a bin, or the interpolation of the extension, are the fixes; the first is the one to make by default.

What the agreement does not show: that a real panel record has three tones. Real spectra have broad humps (random excitation of lightly damped modes, the subject of `P02.19`), harmonics of the engine firing frequency, and electrical pickup at the mains frequency. The tools are the same; the interpretation is harder, and the first question to ask of any peak is "is this on a bin, and what window was used?".

Improvements worth making: Welch averaging of Hann-windowed segments for a noise floor that falls as 1/√(number of segments); parabolic interpolation of the peak bin for frequency to 0.05 Hz; an iterative FFT with bit reversal to remove the recursion overhead. Improvement worth *not* making: zero-padding to "increase resolution". Zero-padding interpolates the spectrum between bins and makes plots smoother; it does not separate two tones that the record length cannot separate.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivation of (2) and (3); `vibration_fft.png`; the four tables of part 4; and a one-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating the three tones with amplitudes and the bias of each reading, the sampling rate you would set the logger to and the anti-aliasing corner you would ask for, and the assumption you trust least.

### Adversarial self-checks

- **What if the record is twice as long?** Resolution halves to 0.25 Hz, the noise per bin falls by √2, and the 83.3 Hz tone sits at δ = 0.2 instead of 0.4, so the rectangular reading rises from 1.51 to 1.87 m/s² though the panel did nothing different. Which reading would you have believed if you had not known the truth?
- **What if the amplitude goes to zero?** The noise floor remains: 0.5 m/s² RMS over 1025 bins gives about 0.02 m/s² per bin, so tones below about 0.05 m/s² are not detectable from 2 s of data. State the detection threshold in your note.
- **Which assumption fails first in the real world?** Stationarity: the engine speed wanders by a per cent or two, which at 50 Hz is a drift across two bins during a 2 s window. The spectrum smears, and the reading drops even for an "on-bin" tone.
- **What measurement could falsify the model?** Re-sample the same signal at a higher rate: if a peak moves, it was an alias. Change the window length: if a peak's amplitude moves, it was leakage. Both are routine checks in flight-test practice.
- **Does the result survive a change in solver or input?** DFT, radix-2 FFT and `scipy.fft` agree to 10⁻⁹. Change the seed: the tone readings move by the noise, about ±0.02 m/s², and nothing else.
- **Which decision would change if the uncertainty doubled?** If the sensor were ±10 % instead of ±5 %, the 83 Hz reading would be 1.8 ± 0.3 m/s² with a Hann window, and the question "is the panel above its 2.0 m/s² endurance limit?" could not be answered either way from this record; you would ask for a longer record and an on-bin frequency grid before adding a stiffener.

## Check yourself

1. Write the DFT and its inverse for N samples, and say where the 1/N goes in the convention used here.
2. A 2 s record at 1024 Hz: what are its frequency resolution, its Nyquist frequency, and the number of bins in the single-sided spectrum?
3. Why does a real record have a conjugate-symmetric spectrum, and what does that mean for the amplitude of a tone split between bins k and N − k?
4. A tone of amplitude 3 m/s² sits exactly half-way between two bins. What does the peak bin read with a rectangular window? With a Hann window, roughly?
5. State the shift theorem and explain why the magnitude spectrum of a record does not depend on when the recording started.
6. A 410 Hz tone is sampled at 300 Hz. At what frequency does it appear, and what would prevent it?
7. Why is the sampled estimate Δt·DFT of a continuous transform first-order accurate for e^{−at}, and what one change makes it second order?
8. State Parseval's theorem for the DFT and use it to write the RMS of a record in terms of its single-sided amplitude spectrum.
