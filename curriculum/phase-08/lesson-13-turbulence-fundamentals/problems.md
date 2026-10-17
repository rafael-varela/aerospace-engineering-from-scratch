# P08.13  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer, give an uncertainty for every statistic you estimate from a record, and label every signal you analyse as **synthetic** or as **measured** (with its source and license).

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Unless stated otherwise the synthetic grid-turbulence record of the lesson is used: U = 10 m/s, σ = 0.3 m/s, L = 0.05 m, sampled at 20 kHz, air at 20 °C (ν = 1.5033 × 10⁻⁵ m²/s).

## 1. (Warm-up: statistics of one short record)

Generate the **synthetic** record `synthetic_signal(2**16, 20000.0, 10.0, 0.3, 0.05, seed=s)` (3.28 s).

(a) Compute the mean, variance, turbulence intensity, skewness and flatness. Compare the variance with `band_variance` for the same record length, and express the difference in units of the expected standard deviation √(2∫ρ²dτ/T_rec), with ∫ρ²dτ = 0.866 L/U for this model.
(b) Compute the integral time scale and length scale from the autocorrelation, and L from the spectral plateau with `welch_psd(u, fs, 2**12)`. Which estimate is biased, and why at this segment length?
(c) Fit the inertial range from 300 to 3000 Hz. Report the slope and ε, and compare ε with `model_epsilon`.

Answer form: five numbers, a comparison in standard deviations; three lengths and a sentence; a slope, two dissipation rates and a percentage.

## 2. (Derivation: what averaging does, and what a record can tell you)

(a) Substitute u_i = Ū_i + u′_i into the incompressible x-momentum equation and average. Show that the mean flow obeys the Navier–Stokes equations with an extra stress −ρ⟨u′_i u′_j⟩, the Reynolds stress, and that continuity holds separately for the mean and the fluctuation. Why is this a closure problem?
(b) For a stationary Gaussian signal with autocorrelation ρ(τ), the variance of the sample variance over a record of length T ≫ T_int is (2σ⁴/T)∫_{−∞}^{∞}ρ²dτ. Use it to derive the record length for a relative standard error ε_r in the variance. Show that for an exponential ρ = e^(−|τ|/T_int) this gives T = 2T_int/ε_r², the formula in `record_length_for`.
(c) Find the only combinations of ν (m²/s) and ε (m²/s³) with units of length, time and velocity, and show that the Reynolds number built from them is 1.
(d) Show that the one-sided frequency spectrum of a signal with integral time scale T_int satisfies S(0) = 4σ²T_int, and that the von Kármán spectrum in the lesson integrates to σ² and has integral length scale L. (∫₀^∞(1 + s²)^(−5/6)ds = √π Γ(1/3)/(2Γ(5/6)).)

Answer form: four derivations; the final formulas boxed.

## 3. (Unfamiliar scenario: the light aeroplane in rough air)

The light aeroplane cruises at 60 m/s at 3000 m (ν = 1.8630 × 10⁻⁵ m²/s) through atmospheric turbulence that a gust probe on its nose records. Model the longitudinal gusts with the same von Kármán spectrum, σ = 0.5 + s/30 m/s and L = 200 + 10s m.

(a) Compute the turbulence intensity seen by the aircraft, the integral time scale in the probe's record, and the frequency U/L of the energy-containing eddies.
(b) Compute ε from `model_epsilon` and the Kolmogorov length, time and the frequency U/(2πη) at which the smallest eddies pass the probe. What sampling rate would the probe need to resolve them?
(c) How long (in seconds and in kilometres flown) must a straight flight leg be to estimate σ² to ±10 % and to ±5 % (one standard deviation)? Use `record_length_for` and state the correction for the von Kármán ρ.
(d) Is Taylor's hypothesis safer here or in the wind tunnel of the lesson? Give the number that decides it.

Answer form: three numbers; four numbers and a rate; four durations and two distances with a correction factor; a sentence and a number.

## 4. (Numerical, uses your code: is my error bar right?)

Generate 40 **synthetic** records of 2¹⁶ samples with seeds s, s + 1, …, s + 39 (U, σ, L, fs as in the lesson).

(a) For each, compute the variance, the autocorrelation length scale, the spectral-plateau length scale (nperseg 2¹²) and ε (300–3000 Hz).
(b) Compare the scatter of the variance (standard deviation over the 40, divided by the band variance) with the prediction √(2∫ρ²dτ/T_rec). Do they agree within the uncertainty of a standard deviation estimated from 40 samples (about 11 %)?
(c) Report the mean and standard deviation of the two length-scale estimates. Which estimator is biased at this record length, by how much, and why?
(d) Repeat (b) with 2²⁰-sample records (you may use fewer seeds). Does the scatter fall as 1/√T_rec?

Answer form: four columns of 40 numbers summarised by mean and standard deviation; a comparison; two means with standard deviations and a sentence; a ratio.

## 5. (Numerical, uses your code: what the instrument does to the spectrum)

Use the lesson's 2²⁰-sample **synthetic** record (seed 2026).

(a) A hot wire of length ℓ averages the velocity along its length. Model it, with Taylor's hypothesis, as a moving average over ℓ/U in time (round to a whole number of samples). For ℓ = 1 mm and 3 mm, recompute the variance, the inertial-range slope and ε (300–3000 Hz). How large is the error in ε, and what does it do to your estimate of η?
(b) Now sample the record at 5 kHz by keeping every fourth sample, **without** a low-pass filter. Recompute the spectrum (nperseg 2¹²) and fit 300–1500 Hz. Explain the change in slope and ε. Compare S at 2000 Hz with the full-rate spectrum and the model.
(c) A colleague proposes a 1 mm wire sampled at 5 kHz for the lesson's tunnel. Give your verdict with numbers.

Answer form: a table for two wire lengths; a slope, an ε and a paragraph; a verdict.

## 6. (Judgment: what a synthetic signal cannot tell you)

(a) The synthetic signal has skewness 0 and flatness 3 by construction. Name two statistics of real turbulence that it cannot reproduce, and why each matters for an engineer (for example, for fatigue under gust loads or for mixing in a combustor).
(b) The Johns Hopkins Turbulence Database (`@jhtdb-johns-hopkins-turbulence-database`, ODC-By) offers a forced isotropic turbulence simulation. Describe, without downloading anything into this repository, how you would use one line of its velocity field to test the lesson's analysis: what you would extract, what you would compare with the synthetic results, and how you would credit the source. Is it a measurement?
(c) Taylor's hypothesis is tested in a jet shear layer where the local turbulence intensity is 20 %. Would you trust an integral scale computed with it there? Give the reason in one sentence.

Answer form: two statistics with a reason each; a short plan; a sentence.
