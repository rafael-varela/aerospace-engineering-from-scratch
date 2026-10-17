# P05.06  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. d = 6 + s/10 mm; A = πd²/4. 2. σ = F/A (MPa); ε = σ/E with E of 7075-T6; ΔL = εL; Δd = −νεd. 3. `true_stress_strain(σ, ε)`; ratio σ_t/σ = 1 + ε. 4. F_y = Fty × A.

## 2.
(a) AL = A₀L₀ ⇒ A = A₀/(1 + ε); σ_t = F/A; ε_t = ∫dL/L from L₀ to L. Uniform deformation fails once a neck forms. (b) dF = A dσ_t + σ_t dA = 0 with dA = −A dε_t. Power law: mKε_t^(m−1) = Kε_tᵐ ⇒ ε_t = m. Evaluate ln(1 + ε_u) from `considere_point(73.1e3, 345.0, 9.0)` and compare with 1/9; the Ramberg–Osgood plastic part behaves like a power law with exponent 1/n, but the elastic part and the 0.002 coefficient shift the point. (c) Put σ = σ_y in the plastic term; for offset o, 0.002(σ/σ_y)ⁿ = o ⇒ σ = σ_y(o/0.002)^(1/n). (d) Diagonal strain ε_d = (τ + ντ)/E; the angle change of the square's corner is γ = 2ε_d (small-strain geometry); τ/γ = E/2(1 + ν).

## 3.
(a) Hemisphere FBD: p πR² = σ 2πRt. (b) `hooke_strain(np.diag([σ, σ, 0]), E, ν)`; ΔR = ε₁₁R; Δt = ε₃₃t; V = 4πR³/3, ΔV = 3ε₁₁V. (c) Uniaxial: ΔR = σR/E; ratio 1/(1 − ν). (d) σ/Fty; real tanks carry factors for proof and burst pressure, fatigue from fill cycles, welds and temperature.

## 4.
(a) Generate, `analyze_tensile(x, y, window=21)`. (b) Table of inputs (E, σ_y, Considère ε_u and σ_u, ε_f − σ_f/E) against results. (c) Loop over (lo, hi). Plastic strain at the top of each window from (6) with n = 20. (d) Loop over seeds 0–99; `np.mean`, `np.std`; |your − mean| < 2 sd?

## 5.
(a) `analyze_tensile(*synthetic_tensile_test(**COUPON, compliance=c))`; compare E with any aluminium alloy. (b) c_est = 1/E_app − 1/72.49 GPa. (c) ε_rec = σ/E_app + ε_p; offset line σ = E_app(ε_rec − ε₀ − 0.002); solve for the crossing. Run with two values of c. (d) `offset_yield(x, y, 73.1e3, r.eps0)`; difference from the fitted-modulus result.

## 6.
(a) Strain source, modulus used for the offset line. (b) U ≈ (Fty + Ftu)/2 × elongation (MJ/m³); per mass U/ρ; m = KE/(U/ρ) with KE = ½ × 90 × 4². Resilience σ_y²/2E: what an elastic spring returns. (c) Allowables as statistical lower bounds. (d) Forming, stretch-forming, drawing, crash: strains to tens of per cent, where the engineering stress no longer refers to the current section.
