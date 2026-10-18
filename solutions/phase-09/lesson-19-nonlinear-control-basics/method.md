# P09.19  Method

The route to each answer, without the final numbers. Pair with `hints.md`; `worked-solutions.md` has the numbers for s = 15.

## 1. Lyapunov by hand

(a) AᵀP + PA + I = 0 gives (1,1): −2kp₁₂ + 1 = 0; (2,2): 2p₁₂ − 2cp₂₂ + 1 = 0; (1,2): p₁₁ − cp₁₂ − kp₂₂ = 0. So p₁₂ = 1/(2k), p₂₂ = (1 + 2p₁₂)/(2c), p₁₁ = cp₁₂ + kp₂₂. (b) Eigenvalues of P (`np.linalg.eigvalsh`); V̇ = xᵀ(AᵀP + PA)x = −xᵀx. (c) `lyapunov_solve(A, I)` and `solve_continuous_lyapunov(A.T, -I)`. (d) p₂₂ and p₁₁ grow like 1/c: an undamped oscillator conserves energy, so no V can decrease at rate |x|².

## 2. Derivations

(a) d/dt(e^{Aᵀt}Qe^{At}) = Aᵀ(·) + (·)A; integrate from 0 to ∞ (the integrand decays because A is Hurwitz): [·]₀^∞ = −Q = AᵀP + PA. P > 0 because xᵀPx = ∫|Q^{1/2}e^{At}x|²dt > 0. Converse: if Av = λv with Re λ ≥ 0, then v*(AᵀP + PA)v = 2 Re λ v*Pv ≥ 0 but must equal −v*Qv < 0. (b) With β = arcsin(δ/A), b₁ = (4/π)∫₀^{π/2} y sin θ dθ = (4/π)[∫₀^β kA sin²θ dθ + ∫_β^{π/2} kδ sin θ dθ]; evaluate and divide by A. (c) Output +M for θ ∈ [β, π + β), −M otherwise, β = arcsin(h/A): b₁ = (4M/π)cos β, a₁ = −(4M/π)sin β; then −1/N = −(πA/(4M))/(cos β − j sin β) = −(πA/(4M))(cos β + j sin β). (d) u = ml²v − mgl sin θ + bθ̇; it requires m, l, b, g, and the measured θ and θ̇, exactly, and unlimited torque.

## 3. Thruster limit cycle

(a) arg G = −90° − arctan(ωτ) > −180° for all ω, so −1/N (the negative real axis) is never reached; simulate with h = 0 and count switches (chattering). (b) `predict_relay_hysteresis_limit_cycle([K], [τ, 1, 0], 1, 0.01)`. (c) `simulate_relay_loop(..., hysteresis=0.01, t_end=6)` and `oscillation_amplitude_frequency(t, e, 4)`. (d) Reduce h: smaller amplitude, higher frequency, so more firings per second (two per cycle: ω/π); compute both for two other h.

## 4. Your loops

(a) `predict_saturation_limit_cycle(nG, dG, m/0.24, 1.0)`, `simulate_saturation_loop(..., t_end=4)`, measure after 3 s; percentage errors. (b) For each Q: `lyapunov_solve(Acl, Q)`, `roa_estimate`, at-rest angle √(c*/P₁₁), and 100 seeded random starts inside the ellipse through `converges`.

## 5. Judgment

(a) Plants that do not filter harmonics (a resonance near 3ω, a pure integrator chain with weak roll-off), asymmetric nonlinearities, several nonlinearities, or loops whose −1/N locus touches G at a grazing angle. (b) Not wrong: a sufficient condition with an ellipse-shaped set; less conservative with non-quadratic V or by enlarging Q in the direction of the true region. (c) The command is clipped, the cancellation is incomplete, the closed loop is nonlinear again; recovery is bounded by the static limit. (d) Linearise the mismatched loop at 0 and apply Routh to s² + a₁s + a₀.

## 6. Mismatch and torque limit

(a) Two simulations from 60° with the plant mass scaled; subtract the linear closed form; max absolute difference. (b) Bisection on θ₀ ∈ [0, 80°] with `simulate_pendulum` over 10 s, convergence if ‖x(10)‖ < 10⁻³. (c) At rest the controller can apply at most 6 N·m, gravity applies mgl sin θ₀.
