# P08.13  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. `synthetic_signal`, `moments`, `band_variance`; relative difference divided by √(2 × 0.866(L/U)/T_rec). 2. `reynolds_decomposition`, `autocorrelation(u′, 4000)`, `integral_time_scale`; L = UT. 3. `welch_psd(u, fs, 2**12)`, `integral_scale_from_spectrum`. 4. `fit_inertial_range(f, S, U, 300, 3000)`; `model_epsilon(0.3, 0.05)`.

## 2.
(a) Average each term; ⟨Ū u′⟩ = 0; ⟨u′_j ∂u′_i/∂x_j⟩ = ∂⟨u′_i u′_j⟩/∂x_j using continuity of u′; move to the stress side. (b) Var(σ̂²) = (1/T²)∫∫2R(t − t′)²dt dt′ ≈ (2/T)∫R²dτ; set the relative standard error to ε_r and solve for T; ∫e^(−2|τ|/T)dτ = T. (c) [ν^a ε^b] = m^(2a+2b) s^(−a−3b). (d) S(f) = 4∫₀^∞R(τ)cos(2πfτ)dτ; f = 0; substitute s = 1.339kL.

## 3.
1. σ/U; L/U; U/L. 2. `model_epsilon(σ, L)`; `kolmogorov_scales(ε, ν)`; f = U/(2πη). 3. `record_length_for(L/U, 0.1 and 0.05)`, times 0.866 for the von Kármán ρ; distance = U T. 4. Compare σ/U in the two cases.

## 4.
1. Loop over 40 seeds; store the four statistics. 2. Std of the variances (ddof = 1) over band variance against the prediction. 3. Means and stds of both L estimates against 50 mm. 4. Repeat with 2²⁰ samples; ratio of scatters against √16.

## 5.
1. w = round(ℓ/U × f_s); `np.convolve(u, ones(w)/w, "same")`; variance ratio, `welch_psd`, `fit_inertial_range`; η ∝ ε^(−1/4). 2. `u[::4]`, f_s/4, `welch_psd(..., 2**12)`, fit 300–1500 Hz; compare S at 2 kHz. 3. Compare ℓ with η and the Nyquist frequency with f_η.

## 6.
(a) Intermittency (flatness of velocity increments growing at small scales) and the negative skewness of du/dx (the signature of the cascade); fat-tailed gust distributions for fatigue; fine-scale mixing. (b) Use the database's cutout service or client to extract one line of u, treat the line as a "frozen" record, run the lesson's analysis, compare E11(k) with −5/3 and with the synthetic results; cite per the database's Citing page; label it simulation (DNS) data, not measurement. (c) Taylor needs u′/U ≪ 1.
