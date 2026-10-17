# P02.19  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. b = V/L; f_b = b/2π; Φ_w(0) = σ_w²L/(πV). 2. Let g(x) = (1 + 3x²)/(1 + x²)²; g′(x) = [6x(1 + x²)² − (1 + 3x²)·2(1 + x²)·2x]/(1 + x²)⁴ = 2x(1 − 3x²)/(1 + x²)³, zero at x = 1/√3; g(1/√3) = 9/8. 3. T_int = ∫₀^∞ R_w dτ/σ² = L/(2V); zero of R_w at aτ/2 = 1, τ = 2L/V. 4. Evaluate `dryden_psd` at b/√3 and its neighbours; evaluate `dryden_autocorrelation` at 2L/V.

## 2.
(a) Substitute ω = bx; with x = tan θ, dx = sec²θ dθ, 1 + x² = sec²θ: ∫dx/(1 + x²)² = ∫cos²θ dθ = π/4 over [0, π/2], and ∫x²dx/(1 + x²)² = ∫sin²θ dθ = π/4. (b) F{e^{−a|τ|}} = 2a/(a² + ω²); differentiate with respect to a: F{−|τ|e^{−a|τ|}} = 2(a² + ω²)⁻¹ − 4a²(a² + ω²)⁻² = 2(ω² − a²)/(a² + ω²)². Combine for σ²[(1 − a|τ|/2)e^{−a|τ|}] and divide by π. (c) ∫₀^∞(1 − aτ + a²τ²/4)e^{−2aτ}dτ = 1/(2a) − a/(4a²) + (a²/4)(2/(8a³)); then Var(σ̂²) = (2/T)·2·σ⁴·(that), relative SD = √Var/σ². Assumption: T ≫ L/V so the double integral reduces to a single one. (d) ω²Φ_w → 3σ²V/(πL) as ω → ∞, so the integrand does not decay.

## 3.
(a) Output spectrum Φ_θ = Φ_n/ω² = G₀/ω²; ∫₀^∞ dω/ω² diverges at 0. For fixed t, θ(t) = Σ increments, Var = Σ Var(increments) = S₀t with S₀ = G₀/2 (two-sided = one-sided/2 in this lesson's convention, since the one-sided density is S/π per rad/s, i.e. 2S per Hz). (b) σ_θ = √(S₀t). (c) Loop over seeds; `np.cumsum(e) * dt`; SD across seeds; expected relative scatter of an SD from n samples ≈ 1/√(2n). (d) `welch_psd(e, dt, 4096)`; `band_average` over 1–40 Hz.

## 4.
(a) τ = 2m/(ρVSC_Lα); (11) with b = V/L, c = 1/τ; σ_n = σ_a/g; compare with σ_w/(gτ). (b) `generate_gust(2000, 0.02, 1.5, 300, 100, "w", seed=s)`, `simulate_heave`; scatter: variance √(5L/(4VT)), RMS half of that. (c) `gaussian_exceedance_fraction(0.25 g, σ_a)` vs `np.mean(np.abs(a) > 0.25 g)`. (d) `welch_psd(a/g, 0.02, 8192)`; model |H_a|²·`dryden_psd_hz`/g²; `band_average` in each decade; the peak: for b ≪ ω ≪ c, |H_a|² ≈ ω² and Φ_w ≈ 3σ²V/(πLω²) would give a flat product, so the peak sits where the two corner frequencies pull it, near the geometric mean √(bc).

## 5.
(a) Write each quantity as an integral over ω and identify which frequency range dominates it; for σ_a note that the integrand of |H_a|²Φ_w decays as ω⁻² beyond c, so it converges, but more than half of it comes from above c. (b) Compare the sample variance (or the Welch level) of successive 5-minute blocks; the record you trust is the longest run of blocks whose variances agree within the scatter (13) predicts. (c) σ_ȧ² = ∫ω²|H_a|²Φ_w dω with an integrand tending to a constant; name the missing roll-off (unsteady lift, structural filtering). (d) σ_a² ≈ (3/2)σ_w²bc with b, c ∝ V, so σ_n ∝ V; V_new = V·(0.10/σ_n); name the assumption (quasi-steady lift at the new speed, or the plunge-only model).

## 6.
(a) T = 5L/(4V·rel_sd²). (b) Segment length = 1/Δf seconds = that × 50 samples; K = floor((N − nperseg/2)/(nperseg/2)); per-bin relative SD ≈ 1/√K. (c) `generate_gust(T, 0.02, ...)`, `welch_psd(w, 0.02, nperseg)`, ratios against `dryden_psd_hz`; per-bin SD/mean of the ratio. (d) Both requirements scale with T but in opposite directions for a fixed T: more segments means shorter segments. Report the variance with (13) and either the fine-resolution spectrum with its 1/√K bin error or a band-averaged spectrum with the averaged error.
