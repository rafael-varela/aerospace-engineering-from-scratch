# P08.13  Turbulence Fundamentals

> Turbulence is random in detail and orderly in its statistics: a mean, a variance, a correlation time and a spectrum, each with an error bar set by how long you listened.

| | |
|---|---|
| **Learning outcomes** | 1. Apply the Reynolds decomposition to a velocity record and to the Navier–Stokes equations, and explain the Reynolds stress and the closure problem. 2. Implement and test the statistics of a turbulent signal from scratch: moments, the autocorrelation and integral scale, Taylor's frozen-turbulence hypothesis, the Welch power spectrum, the −5/3 inertial range, the dissipation rate and the Kolmogorov scales. 3. Show that the analysis recovers the known variance, integral scale and dissipation rate of a seeded synthetic signal with a stated model spectrum, within the statistical scatter that a finite record predicts, and derive that scatter. 4. Decide what record length, sampling rate and probe size a measurement needs, and what a synthetic Gaussian signal can and cannot test. |
| **Prerequisites** | `P08.12` (instability and transition), `P08.10` (boundary layers), `P08.06` (Navier–Stokes), `P02.19` (random processes, autocorrelation, power spectra), `P02.10` (FFT), `P00.14` (measurement uncertainty), `P02.18` (regression on a log–log plot). Skill check: you can compute an FFT with NumPy and say what its frequency bins are. |
| **Engineering question** | Behind a turbulence grid in a wind tunnel running at 10 m/s, a hot wire sampled at 20 kHz records u(t); here the record is a **seeded synthetic signal** with a von Kármán model spectrum (σ = 0.3 m/s, L = 50 mm). What are the turbulence intensity, the integral scale, the dissipation rate and the Kolmogorov scale, is a 1 mm wire small enough, and how long must the record be? Answer: Tu = 2.97 %, L = 54 mm from the autocorrelation (model 50 mm, inside the estimator's ±5 % scatter), ε = 0.385 m²/s³ (model 0.3855), η = 0.31 mm; a 1 mm wire spans 3.3η and under-reads ε by 11 %; 4 s of record gives the variance to ±5 % (one standard deviation), 52 s to ±1.3 %. |
| **Assumptions and limits** | Statistically stationary, single-point record; Taylor's hypothesis (u′/U ≪ 1; 3 % here); Kolmogorov's inertial-range law with C1 = 0.49. **The data are synthetic**: a Gaussian random process with the von Kármán spectrum, built in the frequency domain; it has the right spectrum and two-point statistics but no intermittency, no derivative skewness and no dissipation range, and Taylor's hypothesis is built in, so it cannot be tested with it. |
| **Analytical work** | Reynolds averaging of the momentum equation; the integral scale and S(0) = 4σ²T; the variance of the sample variance and the record length T = 2∫ρ²dτ/ε_r²; Kolmogorov's scales by dimensional analysis; the von Kármán spectrum's variance, integral scale and implied ε. |
| **Computational deliverable** | `starter/main.py`: `reynolds_decomposition`, `moments`, `autocorrelation`, `integral_time_scale`, `welch_psd`, `integral_scale_from_spectrum`, `fit_inertial_range`, `kolmogorov_scales`, `record_length_for` (`von_karman_e11`, `model_epsilon`, `synthetic_signal`, `band_variance` and `main` are given). Tests prove: the decomposition is exact; the moments of a square wave, a sinusoid and a spike train are right; the autocorrelation of a cosine and of white noise match theory and the estimator is linear, not circular; the integral scale of an exponential is exact; the Welch spectrum equals `scipy.signal.welch` and integrates to the variance; and on the synthetic record the variance, integral scale (two ways), −5/3 slope and ε are recovered within about four standard deviations of their estimators; Kolmogorov scales have unit Reynolds number; results are unit-invariant. |
| **Practical activity** | **data-based investigation** of a **seeded synthetic** turbulence record (labelled synthetic throughout): error bars from an ensemble of records, two estimators of the integral scale, the inertial range and ε, and what a finite probe and a low sampling rate do to the spectrum; optionally, the same analysis on a line of the Johns Hopkins Turbulence Database (simulation data, ODC-By) that the learner fetches. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `turbulence_synthetic.png`, the derivations of the record-length formula and the Kolmogorov scales, the experiment tables, and a half-page note specifying a hot-wire measurement (wire length, sampling rate, filter, record length) for the tunnel, with the error bar of each statistic. |
| **Assessment** | `assessments/phase-08/grading-rubric.md`; competencies 1–4. Pass: the Reynolds decomposition and the record-length formula are derived; all public tests pass; every statistic is reported with an uncertainty that the ensemble experiment confirms; the note justifies the wire length and sampling rate from η and f_η; synthetic data are labelled synthetic everywhere. |
| **References** | `@leishman-intro-aerospace-flight-vehicles` "Turbulent Flows" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/turbulent-flows/): Visualizing Turbulent Flows; Reynolds Decomposition; Turbulence Scales & Turbulent Kinetic Energy (required-free; link only). `@mit-ocw-16-100-aerodynamics` Lecture 26 "Effect of Turbulent Fluctuations on Mean Flow: Reynolds-Averaging" (required-free). `@mit-ocw-2-29-numerical-fluid-mechanics` Lecture 24 (Properties of Turbulent Flows; Energy Cascade and Scales) and Lecture 25 (Kolmogorov microscale; Direct Numerical Simulations (DNS)) (required-free). `@libretexts-liburdy-intermediate-fluid-mechanics` Ch. 10 "Introduction to Turbulence Effects": Reynolds Averaging (required-free). `@mit-ocw-2-27-turbulent-flow-transport` "Concepts in Turbulence": 3.2 Features of turbulent flows; 3.3 The energy cascade and the smallest scales of turbulent motion; 3.5 Statistical averages of random quantities in turbulent flow (supplementary-free). `@libretexts-smyth-all-things-flow` 5.4 "Pancakes and noodles- the geometry of turbulence" (supplementary-free). `@jhtdb-johns-hopkins-turbulence-database` forced isotropic turbulence, https://turbulence.idies.jhu.edu/home (supplementary-free; ODC-By; simulation data, fetched by the learner, never by the tests). Optional textbooks for deeper study: Pope, *Turbulent Flows* (paid; library copy; see `references/phase-08.md`). |
| **Extensions** | ★ Add a dissipation range to the model spectrum (multiply by exp(−β(kη)) or Pope's f_η) and recover η from the roll-off rather than from ε. ★ Compute the second-order structure function ⟨(u(t + τ) − u(t))²⟩ and show its τ^(2/3) range gives the same ε. ★ Fetch one line of the JHTDB isotropic field and compare the flatness of velocity increments at small separations with the synthetic signal's 3. |
| **Time** | unmeasured |

## 1. The engineering problem

`P08.12` ended where the boundary layer turns turbulent. From there on, the velocity at a point is never steady: it fluctuates randomly around a mean over a wide range of frequencies. Engineers cannot predict the individual fluctuations and do not need to. They need their **statistics**: the mean (the drag of `P08.14` comes from the mean wall shear), the variance (the turbulence intensity of a wind tunnel, the gust loads on an aircraft), how long the signal stays correlated (the size of the energetic eddies, which sets mixing), and how the energy is spread over scales (which tells you how small a probe or a CFD cell must be).

Every one of those comes from a finite record, and a finite record of a random signal gives a random answer. A tunnel engineer who reports "Tu = 2.97 %" without an error bar has not finished; one who reports an integral scale from a 3-second record to three figures is fooling himself. This lesson builds the standard analysis of a single-point velocity record and, more importantly, checks it.

The check needs a signal whose true statistics are known. Real data do not come with their answers, so the lesson uses a **synthetic** signal: a seeded Gaussian random process whose spectrum is the von Kármán model with a chosen variance σ² = 0.09 m²/s², integral scale L = 50 mm and mean speed U = 10 m/s, sampled at 20 kHz like a hot wire behind a turbulence grid. **It is synthetic everywhere it appears**: in its generator's docstring, in the figure titles, in the experiment tables and in the worked solutions. Success criteria, written first: the analysis recovers the variance within the scatter √(2∫ρ²dτ/T_rec) it predicts (1.3 % for a 52 s record), the integral scale within its estimator's scatter, the −5/3 slope within 0.05 and the model's ε within 3 %, and an ensemble of 40 records confirms the predicted error bars.

## 2. Understand and derive

### 2.1 Reynolds decomposition and the closure problem

Split the velocity into a mean and a fluctuation, **u = Ū + u′**, where Ū is the time average over a record (for a stationary signal) and ⟨u′⟩ = 0 by definition. Averaging the momentum equation (`P08.06`) with this split (Problem 2a) gives the **Reynolds-averaged** equations:

ρ(∂Ū_i/∂t + Ū_j ∂Ū_i/∂x_j) = −∂P̄/∂x_i + ∂/∂x_j(μ ∂Ū_i/∂x_j − ρ⟨u′_i u′_j⟩).   (1)

The mean flow obeys Navier–Stokes with an extra stress, the **Reynolds stress** −ρ⟨u′_i u′_j⟩. Mass conservation holds for the mean and for the fluctuation separately. The averaged equations contain six new unknowns and no equations for them: the **closure problem**, which turbulence models (`P08.15`) address. In a boundary layer −ρ⟨u′v′⟩ is many times the viscous stress away from the wall; it is why a turbulent layer is so much fuller and draggier than a laminar one (`P08.14`).

### 2.2 Moments, correlation and the integral scale

From a record of n samples: the mean, the variance σ² = ⟨u′²⟩, the **turbulence intensity** Tu = σ/Ū, the **skewness** ⟨u′³⟩/σ³ (asymmetry) and the **flatness** ⟨u′⁴⟩/σ⁴ (tails; 3 for a Gaussian). The **autocorrelation**

ρ(τ) = ⟨u′(t)u′(t + τ)⟩/σ²

tells how long the signal remembers itself. Its integral is the **integral time scale** T_int = ∫₀^∞ρ dτ. The estimate of ρ from a finite record is noisy at large lags, so the integral is taken to the first zero crossing. The biased estimator R(k) = (1/n)Σu′_j u′_(j+k) is computed by FFT with zero padding to 2n, so that it is the linear, not the circular, correlation.

**Taylor's hypothesis.** If u′ ≪ Ū, the eddies are swept past the probe faster than they change, so a time record is a spatial cut: x = Ūt. Then the **integral length scale** is L = ŪT_int, wavenumber k = 2πf/Ū, and the wavenumber spectrum E11(k) = ŪS(f)/(2π). In the tunnel u′/Ū = 3 %.

### 2.3 Spectrum, cascade and Kolmogorov scales

The one-sided power spectral density S(f) spreads the variance over frequency: ∫₀^∞S df = σ². It is the Fourier transform of the autocorrelation (Wiener–Khinchin, `P02.19`), so S(0) = 4σ²T_int: the low-frequency plateau is a second route to the integral scale. **Welch's method** estimates S by cutting the record into segments with 50 % overlap, removing each segment's mean, applying a Hann window, and averaging |FFT|², scaled so that it integrates to the variance.

Energy enters at the large scales (L), is passed to smaller and smaller eddies (the **cascade**) and is dissipated by viscosity at the smallest. In between, at scales much smaller than L and much larger than the dissipative ones, the spectrum can depend only on the energy flux ε (W/kg = m²/s³), so by dimensional analysis

**E11(k) = C1 ε^(2/3) k^(−5/3),  C1 ≈ 0.49:**   (2)

Kolmogorov's **−5/3 law**. Fitting the spectral level in that range gives ε. The dissipative scale can depend only on ν and ε (Problem 2c):

**η = (ν³/ε)^(1/4),  τ_η = (ν/ε)^(1/2),  v_η = (νε)^(1/4),**   (3)

with η v_η/ν = 1: the smallest eddies have a Reynolds number of one.

### 2.4 The model spectrum and what it implies

The synthetic signal uses the von Kármán longitudinal spectrum

E11(k) = (2σ²L/π)/(1 + (1.339kL)²)^(5/6),

which integrates to σ² and has integral scale exactly L (Problem 2d). At large k it becomes (2σ²L/π)(1.339kL)^(−5/3), which matches (2) when **ε = [(2σ²L/π)(1.339L)^(−5/3)/C1]^(3/2) = 0.3855 m²/s³** for σ = 0.3 m/s and L = 0.05 m. So the model fixes the three answers the analysis must recover: σ² (less the 1.3 % above the 10 kHz Nyquist frequency, 0.08884 m²/s² in a 52 s record), L = 50 mm and ε = 0.3855 m²/s³. With ν = 1.5033 × 10⁻⁵ m²/s, (3) gives η = 0.307 mm and a Kolmogorov frequency Ū/(2πη) = 5.2 kHz.

### 2.5 How long to listen

The sample variance of a record of length T is itself a random variable. For a Gaussian signal (Problem 2b)

**Var(σ̂²)/σ⁴ = (2/T)∫_{−∞}^{∞}ρ²dτ,  so  T = 2∫ρ²dτ/ε_r²**   (4)

for a relative standard error ε_r. With an exponential ρ, ∫ρ²dτ = T_int and T = 2T_int/ε_r²; for the von Kármán model ∫ρ²dτ = 0.866T_int (computed from ½∫S²df/σ⁴). At T_int = L/U = 5 ms: ±10 % needs 1 s, ±5 % needs 4 s, ±2 % needs 25 s. The requirement is in integral scales, not samples: sampling faster does not shorten it.

### 2.6 Where the analysis stops being valid

- **Stationarity**: a tunnel that warms up, or the atmosphere over a long flight leg, drifts; (4) assumes it does not.
- **Taylor's hypothesis** fails when u′/Ū is not small (jets, wakes, separated flow) or when eddies evolve as fast as they pass.
- **Probe and sampling**: a wire longer than η averages the small eddies; a sampling rate below twice the highest frequency present aliases energy into the band (part 4).
- **Synthetic data**: Gaussian by construction, so skewness 0 and flatness 3; real turbulence is intermittent, its derivative is negatively skewed, and its spectrum rolls off below η. This signal tests the *analysis*, not turbulence.

## 3. Build the computational model

Open `starter/main.py`. `von_karman_e11`, `model_epsilon`, `synthetic_signal` (the seeded synthetic generator), `band_variance` and `main` are given.

### Step 1: decomposition and moments

`reynolds_decomposition(u)` returns (Ū, u − Ū). `moments(u)` returns the mean, variance, standard deviation, skewness and flatness (population moments, divided by n).

### Step 2: autocorrelation and integral scale

```python
nfft = 1 << int(np.ceil(np.log2(2 * n)))
F = np.fft.rfft(up, nfft)
R = np.fft.irfft(F * np.conj(F), nfft)[: max_lag + 1] / n
rho = R / R[0]
```

`integral_time_scale(rho, dt)` integrates with the trapezoid rule up to and including the first non-positive value.

### Step 3: the spectrum

`welch_psd(u, fs, nperseg)` loops over 50 %-overlapping segments, removes each mean, applies the periodic Hann window `np.hanning(nperseg + 1)[:-1]`, accumulates |rfft|², divides by the number of segments and by f_s Σw², and doubles every bin except f = 0 and the Nyquist bin. It must agree with `scipy.signal.welch(u, fs, window="hann", nperseg=..., noverlap=nperseg//2)` to round-off; that is the first test.

### Step 4: from the spectrum to L, ε and η

`integral_scale_from_spectrum` uses S(0) = 4σ²T on the mean of the first three nonzero bins. `fit_inertial_range(f, S, U, f_lo, f_hi)` converts to (k, E11) with Taylor's hypothesis, fits a straight line to log E11 against log k, and averages (E11k^(5/3)/C1)^(3/2) over the band. `kolmogorov_scales` is (3); `record_length_for` is (4) with the exponential ρ.

Running `main()` prints

```
SYNTHETIC grid-turbulence record (von Karman spectrum; U = 10.0 m/s, sigma = 0.3 m/s, L = 0.05 m),
  1048576 samples at 20000 Hz = 52.4 s, seed 2026
mean 10.0000 m/s; variance 0.08844 m^2/s^2 (model in band 0.08884, total 0.09000); Tu = 2.97 %
skewness 0.019, flatness 3.014 (Gaussian: 0 and 3)
integral time scale 5.402 ms -> L = U T = 54.02 mm (model 50.0 mm, T = 5.0 ms)
Welch PSD (nperseg 2^14): integral of S = 0.08778 m^2/s^2; L from S(0) = 56.54 mm
inertial range 300-3000 Hz: slope -1.658 (Kolmogorov -5/3 = -1.667); epsilon 0.3849 m^2/s^3 (model 0.3855)
Kolmogorov scales (nu = 1.5033e-05): eta = 0.307 mm, tau = 6.249 ms, v = 0.0490 m/s; f_eta = U/(2 pi eta) = 5192 Hz
ratio L/eta = 163; Re_L = sigma L / nu = 998; a 1 mm wire spans 3.3 eta
record length for 10 % (1 sigma) on the variance: 1.00 s
record length for 5 % (1 sigma) on the variance: 4.00 s
record length for 2 % (1 sigma) on the variance: 25.00 s
```

and writes `turbulence_synthetic.png` (every panel titled synthetic): the first 0.2 s of the signal, the autocorrelation with T marked, and the Welch spectrum on the model spectrum with its −5/3 range.

### Verify

Run `python -m pytest curriculum/phase-08/lesson-13-turbulence-fundamentals`. The tests:

1. **Decomposition and moments**: Ū + u′ = u and ⟨u′⟩ = 0 to 10⁻¹²; a square wave has variance a², skewness 0, flatness 1; a sinusoid flatness 1.5; a spike train is strongly skewed.
2. **Autocorrelation**: (1 − k/n)cos(2πk/50) for a cosine; |ρ(k > 0)| < 4/√n for white noise; equal to the direct linear sum for a ramp (no wrap-around); T_int of an exponential to 10⁻³, and the tail after the first zero ignored.
3. **Spectrum**: equal to `scipy.signal.welch` to 10⁻⁸; integrates to the variance within 3 % (the segment means and window remove a little at the lowest bins); a pure tone appears at its frequency with its variance.
4. **The synthetic record** (seeded): reproducible from its seed with the requested mean; variance within 6 % of the band variance (over four standard deviations); skewness and flatness Gaussian; L from the autocorrelation and from S(0) within 20 % of 50 mm (about four standard deviations of each estimator, measured over 40 seeds); slope −5/3 ± 0.05 and ε within 3 % of the model.
5. **Scales and units**: η v_η/ν = 1, v_η²/τ_η = ε, 16ε halves η; record length 4 s for ±5 % at 5 ms and four times as long for half the error; the same record in cm/s has 10⁴ times the variance and the same ρ.
6. `main` writes the figure.

### Compare

`scipy.signal.welch` is the established tool for step 3 and agrees to round-off. The real comparison is with the model: every statistic lands within its predicted scatter. The variance is 0.45 % below the band variance (0.35 standard deviations of 1.29 %). L from the autocorrelation is 54.0 mm, +8 %, within two standard deviations of that estimator (4.7 % over 40 seeds); from the spectral plateau 56.5 mm. ε is within 0.2 % of the model, because the inertial range contains hundreds of nearly independent spectral estimates. The integral scale is the hardest statistic to measure; ε and the mean the easiest.

## 4. Experiment at home, or investigate real data

**Mode: data-based investigation** of a **seeded synthetic** signal (and, optionally, of simulation data the learner fetches from the JHTDB). Write each prediction before running it.

1. **Are the error bars right?** 40 synthetic records of 2¹⁶ samples (3.3 s) and 40 of 2²⁰ samples (52 s), seeds 0–39. Predict the scatter of the variance from (4). *Finding:* 2¹⁶: relative scatter 4.40 %, predicted 5.14 %; 2²⁰: 1.28 %, predicted 1.29 %. The largest deviation among the 40 long records was 3.2 %, about 2.5 standard deviations.
2. **Two estimators of L.** The same ensembles. Predict which is biased. *Finding:* 2¹⁶ records: autocorrelation 51.9 ± 9.9 mm, spectral plateau (nperseg 2¹²) 40.1 ± 3.9 mm, biased 20 % low because its first bins (4.9–14.6 Hz) sit on the spectrum's roll-off (corner at 23.8 Hz); 2²⁰ records: 51.4 ± 2.4 mm and 47.6 ± 2.8 mm (nperseg 2¹⁴). The autocorrelation estimator is unbiased within its scatter and twice as noisy; the spectral one trades bias for variance through the segment length.
3. **The inertial range.** Slope and ε for 300–3000 Hz over the ensembles. *Finding:* 2²⁰ records: slope −1.6641 ± 0.0040, ε = 0.3860 ± 0.0013 m²/s³ (model 0.3855); 2¹⁶ records: −1.669 ± 0.018 and 0.3896 ± 0.0067. ε is the best-determined quantity: 0.35 % scatter from a 52 s record.
4. **What the instrument does.** The seed-2026 record filtered by a hot wire of ℓ = 1 mm and 3 mm (a moving average over ℓ/U), and decimated to 5 kHz without an anti-aliasing filter. Predict the sign of each error in ε. *Finding:* 1 mm: variance −1.1 %, slope −1.76, ε −11.5 %; 3 mm: variance −4.1 %, slope −3.27, ε −64 %; decimated: slope −1.46, ε +25 %, and S at 2 kHz 46 % above the full-rate value. A probe too long under-reads ε; a sampling rate too low over-reads it.
5. **Optional, real simulation data.** Fetch one line of the forced isotropic turbulence field from the Johns Hopkins Turbulence Database (ODC-By; cite it as its Citing page asks; do not commit the file) and run the same analysis on it, with x in place of Ūt. Compare the slope and the flatness of velocity increments at small separations with the synthetic record's. Label the result DNS (simulation) data, not measurement.

Raw data format: CSV tables `seed, n, var, L_acf_m, L_spec_m, slope, eps` (1–3); `case, window_samples, var_ratio, slope, eps` (4); each file's header says "SYNTHETIC: von Karman spectrum, U, sigma, L, fs, seed" (or the JHTDB dataset, time and line for 5).

## 5. Validate and improve

| quantity | model (§2.4) | analysis of the seed-2026 record | difference | expected scatter |
|---|---|---|---|---|
| variance (m²/s²) | 0.08884 in band | 0.08844 | −0.45 % | 1.29 % (4) |
| L (mm), autocorrelation | 50.0 | 54.0 | +8 % | 4.7 % (ensemble) |
| L (mm), S(0), nperseg 2¹⁴ | 50.0 | 56.5 | +13 % | 5.9 % (ensemble), bias −5 % |
| inertial slope | −1.667 | −1.658 | 0.009 | 0.004 (ensemble) |
| ε (m²/s³) | 0.3855 | 0.3849 | −0.2 % | 0.35 % (ensemble) |
| skewness, flatness | 0, 3 | 0.019, 3.014 | — | Gaussian by construction |

The spectral L of this record lies 3.2 ensemble standard deviations above the ensemble mean (47.6 mm) and the slope 1.5: this seed happens to carry a little extra low-frequency energy, and the 40-record ensemble (part 4) shows that the estimators themselves are sound. One record is one sample of a random variable. **The engineering answer, for this synthetic tunnel record: Tu = 2.97 % (±0.6 % of itself), L = 54 ± 2.5 mm, ε = 0.385 ± 0.0013 m²/s³, η = 0.31 mm and f_η = 5.2 kHz. A 1 mm wire spans 3.3η and reads ε 11 % low: use a wire of about η, 0.3 mm, sampled at 20 kHz with a filter below 10 kHz. Four seconds pins the variance to ±5 %; the integral scale needs the full minute.**

What to improve: a dissipation range in the model (Extension 1) so the probe-length experiment is realistic; real data (experiment 5) for anything the Gaussian signal cannot show. What not to improve: the number of samples per integral scale (sampling faster does not shorten the record).

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (4) and (3); `turbulence_synthetic.png`; the four experiment tables (five with the optional JHTDB analysis); and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) specifying the hot-wire measurement for the tunnel and the error bar on each statistic, with every synthetic number labelled synthetic.

### Adversarial self-checks

- **What if the tunnel speed doubles?** With the same grid, σ doubles (Tu fixed), ε rises 8× (ε ∝ σ³/L), η falls by 8^(1/4) = 1.7 to 0.18 mm, and f_η rises 3.4× to about 17.5 kHz, beyond the 10 kHz Nyquist frequency: the probe and the sampling rate both stop being adequate.
- **What if the record is short?** At 3.3 s the variance is uncertain by ±5 % and the integral scale by ±19 % (9.8 mm in 51); a single short record can be "wrong" by 30 % with nothing wrong in the code.
- **Which assumption fails first in the real world?** Gaussianity: real turbulence is intermittent, so the record length needed for higher moments (flatness) is far longer than (4) suggests; then Taylor's hypothesis in any flow with u′/Ū above about 10 %.
- **What measurement could falsify the analysis?** Two hot wires a known distance apart along the stream: the time delay of their cross-correlation peak should equal the separation over Ū if Taylor's hypothesis holds.
- **Does the result survive a change in method?** The integral scale from ρ and from S(0) agree within their scatter on long records; ε from the spectrum could be checked by the structure function (Extension 2).
- **Which decision changes if the uncertainty doubles?** Tunnel acceptance tests that require Tu below a limit: a ±10 % error bar on Tu forces a record four times as long, not a better instrument.

## Check yourself

1. Define the Reynolds decomposition. Which new term appears when the momentum equation is averaged, and why does it create a closure problem?
2. What is the integral time scale, how is it computed from a record, and why is the integral stopped at the first zero crossing?
3. State Taylor's hypothesis and the condition for it to hold. What is u′/Ū in the tunnel of this lesson?
4. What is the −5/3 law, over which range does it hold, and how does it give ε?
5. Write the Kolmogorov length, time and velocity scales. How does L/η scale with the Reynolds number?
6. Why does the synthetic signal have skewness 0 and flatness 3, and what does that mean for what it can test?
7. How long a record gives the variance to ±5 % if T_int = 5 ms?
8. What do a probe longer than η and a sampling rate below twice the highest frequency each do to the measured ε?
