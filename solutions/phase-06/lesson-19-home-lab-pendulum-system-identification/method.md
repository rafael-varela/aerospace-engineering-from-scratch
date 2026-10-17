# P06.19  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) `predict_frequency(L, 0.003, 9.81, 0.01)`; compare u_L/L with u_g/g. (b) f₀ × `amplitude_frequency_ratio` at 10°, 20°, 40°; differences relative to f₀ against 7 × 10⁻⁵. (c) `quadratic_drag_coefficient(1.2, 0.47, 0.04, 0.10, L)`; `equivalent_zeta(20°, ω₀, b2=...)`; with β = 4b₂ω/(3π), A(t) = A₀/(1 + βA₀t).

## 2.
(a) Integrate the damping power b₁θ̇², b₂|θ̇|³, b₀|θ̇| over a period with θ̇ = −Aω sin ωt; divide by the period and by ω²A to get dA/dt; ζ_eq = −(dA/dt)/(ωA). (b) Separate variables: A^(−p) dA = −κ dt. (c) u = 1/A: du/dt = αu + β; solve and invert. (d) δ = (ln A₀ − ln Aₙ)/n; propagate two independent relative errors s; expand δ/√(4π² + δ²).

## 3.
(a) Build t, y with the seed; `fit_damped_cosine`; f_d from the result; RMS against 0.02 g. (b) `find_peaks_interp`; `log_decrement_fit`; `log_decrement(ap, 20)`. (c) `bootstrap_fit(t, y, fit, 100, seed=s)`; compare sd with the covariance; compare the RMS with the noise. (d) Cycles = f × T; envelope at 10 s; think about how much information the late, noise-dominated cycles carry.

## 4.
(a) `read_swing_csv`; mask t ≤ T; `fit_damped_cosine`; `bootstrap_fit(..., block=60, seed=2)`; record RMS. (b) Compare RMS with 0.30° for each length; relate the growing RMS to the amplitude-dependent frequency and ζ. (c) Weigh a longer record against a better model (windows, or the full ODE fit of Extension 1).

## 5.
(a) `make_synthetic_swing(s, theta0_deg=30, duration=60, noise_deg=0.1, offset_deg=0, ...)` for each case; `find_peaks_interp`; `fit_envelope`. (b) Same with the mixed coefficients at 40° for 150 s; `fit_envelope_vq(tp, ap, √(g/L))` and `fit_envelope`. (c) Repeat at 20° for 120 s; compare uncertainties and the correlation; relate to the amplitude range.

## 6.
(a) Compare RMS with noise; decide which estimate includes the model error. (b) For i = 0..3, mask 30i ≤ t < 30(i + 1), shift t, fit, find peaks, compute the mean ratio and f₀; compute ζ_eq from the generator's b₁, b₂ at the window's peaks. (c) Weighted mean and χ² = Σ((f₀ − mean)/u)²; compare χ² with 3 degrees of freedom; if inconsistent, use the scatter of the four values (standard error) instead.
