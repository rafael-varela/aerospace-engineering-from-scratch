# P09.18  Method

The route to each answer, without the final numbers. Pair with `hints.md`; `worked-solutions.md` has the numbers for s = 15.

## 1. Peaks and margins

(a) Evaluate (Ms/(Ms + 1), Ms/(Ms − 1)) and 2 arcsin(1/(2Ms)). (b) |3/(jω + 2)| is largest at ω = 0. For the resonance, set d|G|²/dω = 0: ω_peak = ω_n√(1 − 2ζ²), peak 1/(2ζ√(1 − ζ²)). (c) `tf_hinf_norm(num, den, "bisection")` and `"sweep"`. (d) The bound is a lower bound on both the upward and the downward gain margin, and a separate lower bound on the phase margin (the two hold one at a time; a combined change is covered only inside the disk |1 − e^{jφ}/k| < 1/Ms); convert both ends of the gain range to dB.

## 2. Derivations

(a) y = Gu + d, u = C(r − y − n) ⇒ (1 + L)y = L(r − n) + d ⇒ y = T(r − n) + Sd; r − y = Sr − Sd + Tn (using 1 − T = S); u = C(r − n − y) = KS(r − n − d). S + T = (1 + L)/(1 + L) = 1. (b) min|1 + L| = 1/max|S|. A real-axis crossing at −x is outside the disk iff |x − 1| ≥ 1/Ms; scaling by k reaches −1 at k = 1/x, so the allowed k range follows from the extreme admissible x. At |L| = 1, the triangle 0, −1, L has sides 1, 1 and |1 + L| ≥ 1/Ms; the angle at 0 is 2 arcsin(|1 + L|/2). (c) 1 + L₀(1 + WΔ) = (1 + L₀)(1 + L₀WΔ/(1 + L₀)) = (1 + L₀)(1 + TWΔ). The second factor's Nyquist curve stays inside the disk of radius |TWΔ| < 1 about 1, so it never encircles 0; the zeros of 1 + CG are those of 1 + L₀. Necessity: at ω* with |WT| ≥ 1 choose Δ with Δ(jω*) = −1/(WT)(jω*), realised as an all-pass. (d) H v = jωv with v = (p, q): (jωI − A)p = BBᵀq/γ², (jωI + Aᵀ)q = −CᵀCp. Substitute G(jω) = C(jωI − A)⁻¹B and its conjugate transpose to obtain G(jω)ᴴG(jω)w = γ²w with w = Bᵀq/γ (or the reverse); conversely build v from a singular vector.

## 3. Waterbed

πp = ∫₀^{ω₁} + ∫_{ω₁}^{ω₂} + ∫_{ω₂}^∞ ≤ ω₁ ln ε + (ω₂ − ω₁) ln Ms + 0; rearrange for ln Ms. Repeat with the smaller ω₂. For (c), choose stabilising leads (check with `is_closed_loop_stable`) and call `bode_sensitivity_integral` on S. (d) The bound contains only p, ε, ω₁, ω₂: change p (more static stability) or ω₂ (a faster actuator).

## 4. Your servo

K_c = 0.4 + s/50. Build the loop, compute the three norms and `disk_margin_bounds(Ms)`; both norm routes on `tf_to_ss(*S)`; `robust_stability_index` with `uncertainty_weight()`; for the true loop compute Ms_true = max|1/(1 + L_true)| on the grid, its bounds, and `margins` of L_true.

## 5. Judgment

(a) Separate "stable for every plant in the set" from "stable for this plant". (b) Coverage ratio ≤ 1 means covered. (c) A complex Δ of size r includes rotations of the gain; the true real gain error is a single direction. Compare 1/Mt (complex tolerance) with the range of real k that keeps the loop stable. (d) Weigh the 0.974 against how well the weight was measured; name a frequency-response measurement of the actuator.

## 6. Robust versus true limit

Two bisections on K_c; then at 1.05 × K_c,robust compute ‖WT‖∞ and its frequency with `hinf_norm_sweep`, the target −1/(WT), the all-pass with `destabilising_perturbation`, and the roots of the perturbed characteristic polynomial. Compare Δ_true/W at ω* with the target.
