# P06.16  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. n = (slope in dB/oct)/3.0103; G(20) = G_p(20/50)², G(2000) = G_p(800/2000)². 2. Segment areas G₁f₁[(f₂/f₁)^{n+1} − 1]/(n + 1) for n = 2, 0, −2; sum; square root. 3. `psd_area`, `grms` on the four breakpoints. 4. `np.trapezoid(G_bp, f_bp)`, square root, ratio to the exact Grms; explain with the convexity of a power law on linear axes.

## 2.
(a) Integrate G₁f₁⁻ⁿfⁿ; the antiderivative is fⁿ⁺¹/(n + 1) except at n = −1, where it is ln f; n = −1 is −3.01 dB/oct. (b) Substitute r = 1/u in I₂; show the integrand becomes 1/D(u); add I₀ + 4ζ²I₂ = (1 + 4ζ²)I₀; multiply by fₙ and G; drop 4ζ², replace G by G(fₙ). (c) Parseval with h(t) = e^{−σt} sin(ω_d t)/ω_d, σ = ζωₙ; integrate e^{−2σt} sin² with the half-angle identity; convert ∫dω over the whole axis to ∫dr from 0 to ∞. (d) |Z/Ÿ|² = 1/[(ωₙ² − ω²)² + (2ζωₙω)²]; substitute ω = rωₙ, df = ωₙdr/(2π); use I₀; compare with Miles divided by ωₙ².

## 3.
1. fₙ, G(fₙ) from `psd_loglog`; `miles_grms`, `response_grms`; ratio. 2. Integrate |H|²G below 0.7fₙ and over the whole band on a fine log grid; ratio. 3. 3σ = 3 Grms; force = 3 Grms m g₀; `upcrossing_rate`, `expected_peak_factor(nu, 120)`. 4. `relative_displacement_rms`; 3σ; compare with 0.10 mm.

## 4.
(a) For each T: `synthesize_random(..., T + 0.5, seed=s)`, `simulate_sdof_base`, discard 0.5 s, RMS / `response_grms`; u = `variance_estimate_relative_sd(fn, zeta, T)/2`; (ratio − 1)/u. (b) `csd(a, y, fs, nperseg=8192)` and `welch(a, ...)`; H = P_xy/P_xx; `curve_fit` of √`transmissibility_sq` against |H| on 60–400 Hz; uncertainties from the covariance. Half-power: find the bins where |H|² crosses half its peak on each side, interpolate, ζ ≈ Δf/(2f_peak). (c) Compare the frequency resolution with the half-power bandwidth; say what the covariance does not include.

## 5.
(a) 10 log₁₀((7.2/9.66)²) dB; compare with ±1.5 dB. (b) Solve `response_grms(…, 150, 1/(2Q)) = 7.2` for Q with `brentq`; then solve for fₙ at Q = 10. (c) Candidates: Q lower at full level (amplitude-dependent damping), fₙ lower than predicted, accelerometer not at the mass centre or on a node of a higher mode, input at the low edge of tolerance. Tests: a low-level and a full-level sine survey; check the accelerometer location; read the control PSD. (d) Weigh test-correlated values against flight-to-test differences and the scatter between units.

## 6.
(a) Synthesize, simulate, discard 0.5 s; σ from the record; find i with y[i] > y[i−1], y[i] > y[i+1], y[i] > 0; fractions above kσ; compare with exp(−k²/2). (b) m₀, m₂ from `spectral_moments`; m₄ by trapezoid of f⁴|H|²G on a fine grid; α and √(m₄/m₂); count maxima per second in the record. (c) ν₀⁺T e^{−9/2}; solve ν₀⁺T e^{−k²/2} = 1 for k; compare with `expected_peak_factor`.
