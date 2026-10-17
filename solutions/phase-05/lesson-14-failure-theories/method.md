# P05.14  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. σ_x = (120 + 5s) MPa. Centre c and radius R of Mohr's circle; σ_a, σ_b = c ± R. 2. Sort (σ_a, 0, σ_b) to get σ₁ ≥ σ₂ ≥ σ₃. 3. τ_max = (σ₁ − σ₃)/2; σ_T = σ₁ − σ₃; σ_vM from (5) or √(σ_x² − σ_xσ_y + σ_y² + 3τ_xy²). 4. `margins(plane_stress_tensor(...), MATERIALS["7075-T6"])` and check by hand with (9).

## 2.
(a) ε = ((1 + ν)σ − ν tr(σ)I)/E; u = ½ σ : ε = ((1 + ν)σ : σ − ν(tr σ)²)/(2E). Substitute σ = pI + s with p = tr(σ)/3: σ : σ = 3p² + s : s. Collect the p² terms into (1 − 2ν)(tr σ)²/(6E) and the s : s term into (1 + ν)s : s/(2E). Express s : s in principal stresses: s : s = ⅓[(σ₁ − σ₂)² + (σ₂ − σ₃)² + (σ₃ − σ₁)²]. Tensile test: σ₁ = F_ty, others zero. (b) With σ₂ = σ₃ + f(σ₁ − σ₃): σ_vM² = (σ₁ − σ₃)²(1 − f + f²). The bracket ranges between 3/4 (f = ½) and 1 (f = 0 or 1). Take the square root and the reciprocal. (c) Insert σ₁ = τ, σ₃ = −τ in (7) with n = 1; solve for τ. Set S_t = S_c = S and compare with (2).

## 3.
(a) `principal_stresses(sig)`, `von_mises`, `tresca`, `margins(sig, MATERIALS["Ti-6Al-4V"])`. (b) Recompute with sig − 100e6·I. Compare principal stresses (all shift by −100 MPa) and the two equivalent stresses (unchanged). (c) `strain_energy_density`, `volumetric_energy_density`, `distortion_energy_density` with E = 113.8 GPa, ν = 0.342; fraction u_d/u.

## 4.
1. Scale the loads by (1 + s/100). 2. g(t) = F_ty/σ_vM,max(t) − 1 − 0.15 with σ_vM,max from `critical_point`; `brentq(g, 0.8e-3, 10e-3)`. 3. The same with the Tresca stress at θ = π (the compression side is critical for compressive P; confirm it). 4. MS_ult = F_tu/(1.5σ_vM) − 1 at the chosen t. 5. Mass per metre ρ·A(t) for each t.

## 5.
(a) Principal stresses of (σ_x, 0, τ) including the zero. n_CM = 1/max(σ₁/S_t − σ₃/S_c, σ₁/S_t); n_R = 1/max(σ₁/S_t, −σ₃/S_c). (b) Same with σ_x positive. (c) Compare the size of the tensile principal stress in the two cases; Tresca assumes equal tension and compression strengths, which cast iron does not have.

## 6.
Points to cover: typical versus statistically based allowables and the size of the difference; the 1.5σ_eq stress at ultimate as a linear extrapolation (conservative for a ductile part once it yields, because plastic redistribution caps the peak, but the extrapolation is not what is being disputed); the criterion choice changing the sign of the margin is itself a warning; what the nominal analysis omits (weld or fitting, K_t, fatigue under landing cycles, corrosion of 2024 in a gear bay); MS = 0 means the design has no room for any of these. Recommend, and name the evidence (allowables from MMPDS for the actual product form, a detail stress analysis of the joint, or a test).
