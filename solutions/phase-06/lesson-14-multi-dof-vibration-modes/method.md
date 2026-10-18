# P06.14  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. M = diag(m₁, m₂), K = [[k₁ + k₂, −k₂], [−k₂, k₂]]. 2. m₁m₂ω⁴ − (m₁k₂ + m₂(k₁ + k₂))ω² + k₁k₂ = 0; quadratic formula; f = ω/2π. 3. x₂/x₁ = (k₁ + k₂ − ω²m₁)/k₂; normalise by √(φᵀMφ). 4. Γ = ΦᵀM(1, 1); squares; sum. 5. `chain_matrices`, `modes`, `participation`.

## 2.
(a) φⱼᵀKφᵢ = ωᵢ²φⱼᵀMφᵢ and φᵢᵀKφⱼ = ωⱼ²φᵢᵀMφⱼ; subtract using symmetry. Repeated frequencies: any combination is a mode; choose an M-orthogonal basis. (b) rᵀMr = rᵀMΦΦᵀMr = ΓᵀΓ. (c) η̂ᵢ = −Γᵢ A/(ωᵢ² − ω² + 2iζᵢωᵢω); ẍ = A r − ω²Φη̂; F = rᵀMẍ; use rᵀMΦ = Γᵀ and simplify each term. (d) ΦᵀCΦ = aI + b diag(ω²); divide by 2ωᵢ.

## 3.
(a) K from k₂, k₃ only; K(1, 1, 1) = 0. (b) `eigh(Kf, M)`; compare with `modes` of the fixed stack. (c) mᵀφ for each flexible mode; orthogonality with the rigid mode (1, 1, 1). (d) Project the pulse on each mode.

## 4.
(a) Each spring carries the inertia force of everything above it: (Σ masses above) × a₀. (b) `direct_integration(M, C, K, lambda t: -M @ ones * a0, 0, 0, t)` and `modal_integration` with n_modes = 1, 2, 3; spring forces −kᵢ(zᵢ − zᵢ₋₁); peaks; divide by static. (c) Sum of modal contributions each up to twice their static share. (d) a₀ + z̈₃; differentiate z numerically or use the force on the payload spring divided by m₃.

## 5.
(a) C with c_d at the (2, 3) block; `modal_damping`; eigenvalues of the 6 × 6 matrix; ζ = −Re λ/|λ|, f_d = Im λ/2π. (b) Scan ω near ω₁ with both FRFs; peaks; ratio. (c) `base_sine_accelerations(w1, g, M, C, K)`. (d) Compare φ₃ − φ₂ across modes.

## 6.
(a) Weigh model detail against what the notch depends on. (b) Damping, interface flexibility (testable); mass properties or test hardware differences (not testable by the model alone). (c) What a coupled-loads analysis needs from each side. (d) The role of a low-level sweep.
