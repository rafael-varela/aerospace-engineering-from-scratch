# P05.13  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. C = (σ_x + σ_y)/2, R = √(((σ_x − σ_y)/2)² + τ_xy²); σ₁,₂ = C ± R. 2. θ_p = ½ atan2(2τ_xy, σ_x − σ_y); note the sign: a negative τ_xy turns σ₁ clockwise from x. 3. Maximum in-plane shear R at θ_p − 45°, normal stress C. 4. σ₂ < 0 < σ₁, so the third principal stress (0) lies between them and the absolute maximum equals R. 5. (1) at θ = 30°; check σ_x′ + σ_y′ = σ_x + σ_y. 6. `principal_stresses_2d`, `transform_stress_2d`, `eigvalsh` of the 2 × 2 (or 3 × 3 with zeros) matrix.

## 2.
(a) Wedge: forces along n: σ_x′ dA = σ_x cos θ (dA cos θ) + τ_xy sin θ (dA cos θ) + σ_y sin θ (dA sin θ) + τ_xy cos θ (dA sin θ); along the face similarly; double-angle identities. Equilibrium only. (b) Square and add; differentiate σ_x′ with respect to θ and set to zero; tan 2θ has period 90° in θ; the shear expression's derivative gives the 45° planes. (c) Expand det(σ − sI); identify I₁, I₂, I₃. Deviator: p = I₁/3, J₂ = ½ s:s, J₃ = det s; x = r cos α, r² = 4J₂/3 makes the coefficients match 4cos³ − 3cos; cos 3α = 4J₃/r³. Real roots of a symmetric matrix ⇒ |cos 3α| ≤ 1. (d) Principal stresses σ₁, σ₂, 0; largest circle on [0, σ₁]; planes at 45° between the σ₁ direction and the normal to the skin.

## 3.
1. Axial pr/(2t), hoop pr/t; r/t and p/σ_hoop. 2. Weld normal at β + 90°; `transform_stress_2d(σ_a, σ_h, τ, β + 90°)` with and without τ. 3. `principal_stresses_2d`; zero-shear weld angles are θ_p and θ_p + 90°. 4. Outer: `plane_stress_abs_max_shear`; inner: (σ₁ − (−p))/2. 5. Compare the normal stress across the weld with that across a longitudinal seam (σ_hoop).

## 4.
1. `stress_tensor`, `invariants`, `principal_stresses_3d`, `eigh` and the difference. 2. `principal_direction` for each root; N Nᵀ = I; |n_k · v_k| = 1. 3. `absolute_max_shear`; normal stress (σ₁ + σ₃)/2. 4. `traction_on_plane`, `in_mohr_region`; compare σ_n with I₁/3. 5. `transform_tensor(σ, rotation_z(25°))`, `invariants` again.

## 5.
1. `rosette_strains(eps, (0, 60, 120))`; derive the closed forms from (8). 2. `plane_stress_from_strain`, `principal_stresses_2d`. 3. Principal strains C_ε ± R_ε with R_ε = √(((ε_x − ε_y)/2)² + (γ/2)²); σ₁ = E/(1 − ν²)(ε₁ + νε₂); ε_z = −ν(ε_x + ε_y)/(1 − ν). 4. Seeded loop with noise; generate the 0/45/90 readings from the same (ε_x, ε_y, γ) with (8); same loop.

## 6.
1. Pressure-only stresses; θ_p; crack direction θ_p ± 90°. 2. Required θ_p = 25° + 90° (equivalently −65°); τ = (σ_x − σ_y)/2 · tan 2θ_p; torsional shear T/(2πr²t). 3. When σ_x > σ_y, θ_p from (4) lies within ±45° of x; compare with the required angle. 4. List alternatives (local stress concentrations, residual stresses from riveting, a repair or doubler edge, multiple-site damage linking, a different load case); name the evidence (fractography, strain survey, structural drawings, repair history).
