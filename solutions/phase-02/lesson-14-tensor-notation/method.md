# P02.14  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) σ₁₁ = pr/(2t), σ₂₂ = pr/t, σ₁₂ = T/(2πr²t); zeros elsewhere. (b) c = (σ_x + σ_y)/2, d = (σ_x − σ_y)/2; σ_n = c − d cos 2a − τ sin 2a; τ_seam = |−d sin 2a + τ cos 2a|; `seam_stresses(S, a)`, `transform_tensor2(S, rotation_matrix([0,0,1], a))[1,1]` and `[0,1]`, `mohr_2d(...)[1:]`. (c) I₁ = σ_x + σ_y, I₂ = σ_xσ_y − τ², I₃ = 0; R = √(d² + τ²); σ₁,₂ = c ± R; tan 2θ_p = τ/d, pick the quadrant with `atan2(2τ, σ_x − σ_y)`. (d) τ_max = σ₁/2 (σ₃ = 0); σ_vm = √(σ₁² − σ₁σ₂ + σ₂²); margin = σ_y/σ_vm − 1.

## 2.
(a) Both sides vanish unless j ≠ k, m ≠ n, and {j,k} = {m,n} = complement of i; then both are +1 for (j,k) = (m,n) and −1 for (j,k) = (n,m). Contract (a × b)_i (c × d)_i = ε_ijk ε_imn a_j b_k c_m d_n and apply (2); set c = a, d = b for the second. (b) New basis e′_i = cos a e_i + (1 − cos a)(e_i · n)n + sin a (n × e_i); take the e_j component: Q_ij = e′_i · e_j; the cross term is ε_jkl n_k (e_i)_l = ε_jki n_k = ε_ijk n_k. (c) t′_i = Q_ip t_p = Q_ip σ_pq n_q; invert (4) by multiplying by Q_jq and using Q_iq Q_jq = δ_ij: n_q = Q_jq n′_j; then t′_i = (Q_ip Q_jq σ_pq) n′_j for all n′ and Cauchy's formula in the new frame identifies the bracket as σ′_ij. (d) Contract i = j in (5) and use Q_ip Q_iq = δ_pq; for σ′_ij σ′_ji use it twice; det(QσQᵀ) = det σ. For (9), Q = [[c, s, 0], [−s, c, 0], [0, 0, 1]]; σ′₁₁ = Q₁p Q₁q σ_pq = c²σ_x + 2csτ + s²σ_y, etc.; convert with c² = (1 + cos 2θ)/2, cs = sin 2θ/2.

## 3.
(a) ε_b = (ε_a + ε_c)/2 + ε_xy. (b) Mohr's circle on strain with ε_xy (not γ/2 twice); `principal_stresses(eps)` works for any symmetric tensor. (c) From the inverse law with σ₃₃ = 0 and the plane-stress formulas, ε₃₃ = −ν(ε_x + ε_y)/(1 − ν). Build eps as 3 × 3 with ε₃₃, then `hooke_stress(isotropic_stiffness(E, nu), eps)`; the 3,3 component must come out zero. (d) `hooke_stress` with ε₃₃ = 0 gives λ(ε_x + ε_y) extra on every diagonal term and a non-zero σ₃₃.

## 4.
(a) `seam_stresses` over `np.arange(0, 180, 0.1)`; argmax/argmin; compare with θ_p and θ_p + 90°; shear maximum at θ_p ± 45° (seam-angle 45° from the principal seams). (b) `fuselage_skin_stress(p, r, t, -T)`; compare `seam_stresses(S, -a)` and `seam_stresses(S, π − a)`. (c) `brentq(lambda T: von_mises(fuselage_skin_stress(p, r, t, T)) - sy/1.5, 0, 50e6)`; by hand τ* = √[((σ_y/1.5)² − (σ_x² − σ_xσ_y + σ_y²))/3], T* = 2πr²t τ*. (d) Loop; `principal_stresses` gives all three; the third goes negative once σ_xσ_y < τ² (I₂ < 0); the mode is shear buckling of the panel between frames and stringers.

## 5.
(a) Loop over the four index orders of each symmetric pair; `np.allclose(C, C.transpose(1,0,2,3))` etc.; `hooke_stress(C, diag(e,0,0))`. (b) `transform_tensor4(C, rotation_matrix([0,0,1], -theta))`; read the five entries; the formulas use c = cos θ, s = sin θ. (c) Q̄ = [[C̄₁₁₁₁, C̄₁₁₂₂, 2C̄₁₁₁₂], [C̄₂₂₁₁, C̄₂₂₂₂, 2C̄₂₂₁₂], [C̄₁₂₁₁, C̄₁₂₂₂, 2C̄₁₂₁₂]]; S̄ = Q̄⁻¹; E_x = 1/S̄₁₁, E_y = 1/S̄₂₂, G_xy = 1/S̄₃₃, ν_xy = −S̄₁₂/S̄₁₁; strains = S̄ [100 MPa, 0, 0]; γ_xy ≠ 0 because S̄₁₃ ≠ 0 (shear–extension coupling). (d) Loop over θ; interpolate for E_x = 14 GPa.

## 6.
(a) `yield_margin` and σ_y/(2 `max_shear`) − 1; the ratio σ_vm/(2τ_max) lies in [√3/2, 1]. (b) K_t at a loaded hole ≈ 3 on the seam normal stress; lap-joint bending ≈ 1.5 on one surface; frames share load (reduces); thin wall is exact to t/r. (c) Fatigue: the stress *range* at the hole edge per pressurisation cycle, times the number of cycles. (d) Invariant: I₁, I₂, I₃, principal values, τ_max, σ_vm, the margins. Frame dependent: the components σ_ij, the principal *directions' components*, the seam stresses (which belong to a specific plane, not a frame).
