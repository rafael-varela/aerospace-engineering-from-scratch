# P09.17  Method

The route to each answer, without the final numbers. Pair with `hints.md`; `worked-solutions.md` has the numbers for s = 15.

## 1. One servo pole

(a) A = −1/τ, B = 1/τ: A_d = e^{−T/τ}, B_d = 1 − e^{−T/τ}, C = 1: G_d = B_d/(z − A_d). (b) G(s) = 1/(τs + 1) with s = (2/T)(z − 1)/(z + 1): G_d = (z + 1)/((2τ/T + 1)z + (1 − 2τ/T)), pole (1 − T/(2τ))/(1 + T/(2τ)). Matched: pole e^{−T/τ}, the zero at infinity to z = −1, gain chosen so G_d(1) = 1. (c) `c2d_zoh` on the 1 × 1 matrices, `tustin_tf([1], [τ, 1], T)`, `matched_tf([1], [τ, 1], T)`; DC gain = G_d(1). (d) T/(2τ) > 1 makes the Tustin pole negative: the response contains (negative)^k, an oscillation at the Nyquist frequency that the continuous servo does not have.

## 2. Derivations

(a) Mᵏ = [[Aᵏ, A^{k−1}B], [0, 0]] for k ≥ 1; sum Mᵏ Tᵏ/k!: upper-left e^{AT}, upper-right Σ A^{k−1}TᵏB/k! = ∫₀ᵀ e^{As}ds B, lower-right I. (b) x(t_k + τ) = e^{Aτ}x_k + ∫₀^τ e^{As}ds B u_{k−1}; then x_{k+1} = e^{A(T−τ)}x(t_k + τ) + ∫₀^{T−τ} e^{As}ds B u_k. Combine. For G₀ + G₁, split ∫₀ᵀ e^{As}ds at T − τ and substitute s = T − τ + σ in the upper piece. (c) z = e^{jθ}: (z − 1)/(z + 1) = j tan(θ/2), purely imaginary; |z| < 1 ⇔ Re s < 0 (write z = (1 + sT/2)/(1 − sT/2) and compare distances to ±2/T). With s = jω_a and θ = ω_dT, ω_a = (2/T) tan(ω_dT/2). (d) h(t) = 1 for 0 ≤ t < T: H(s) = (1 − e^{−sT})/s = e^{−sT/2}(e^{sT/2} − e^{−sT/2})/s; at s = jω this is e^{−jωT/2}·2 sin(ωT/2)/ω = T e^{−jωT/2} sinc(ωT/2).

## 3. Notch

(a) `tustin_tf(num, den, 0.01)`, `tustin_tf(num, den, 0.01, w_prewarp=ω_b)`, `matched_tf(num, den, 0.01)`. (b) On f from 1 to 49.9 Hz (a fine grid), |`freqresp_z`| has its minimum at the notch; evaluate 20 log₁₀|N_d| at f_b and compare with the continuous value. (c) Prewarped Tustin and matched pole-zero both place the notch at f_b; prefer prewarped Tustin (it also keeps the gain shape symmetric and maps the whole axis to the circle). Plain Tustin compresses frequencies by the arctangent law. (d) Solve (1/(πT)) arctan(πf_bT) = 0.99 f_b for 1/T, or scan f_s.

## 4. Your flight computer

(a) `sample_rate_study` with your rate returns spectral radius, ζ_eq, overshoot and PM for τ/T = 0, 0.5, 1. (b) Rule: 65.75° − (180/π)·7.12·T·(0.5 + τ/T). (c) Loop f upward from 10 Hz and stop at the first f with `discrete_loop_margin(1/f, 1/f)[0] ≥ 45`. (d) `hybrid_simulation(T, T/2, 3.0)` against `dlsim(*discrete_pitch_loop(T, T/2), ...)` at the samples.

## 5. Judgment

(a) Loss ≈ 1.5ω_cT rad = 540°/(ω_s/ω_c); set it to 10° and compare with the measured margins at a few rates. (b) Anti-aliasing filter, sensor dynamics (gyro bandwidth), actuator electronics, digital filtering of the measurement, network latency; evaluate one at ω_c. (c) Design for τ = T, the worst case, because the margin is monotone in τ and jitter is rarely bounded below T in practice. (d) When ω_s/ω_c is small (below about 10), when the plant has dynamics near the Nyquist frequency, or when deadbeat or minimum-time behaviour is wanted.

## 6. Digital gain margin

(a) For each (f_s, τ): bracket k by doubling, then bisect 60 times on `spectral_radius(discrete_pitch_loop(T, τ, kK)[0]) < 1`. (b) At the critical k take the eigenvalue of largest modulus and compute Im(ln z)/T. (c) Compare with π/T; the extra delay adds phase ωτ, so the phase reaches −180° at a lower frequency, where the plant's gain is larger. (d) Scan k on a grid to 100 and record max Re eig(A − kBK); look at the root locus asymptotes for large k.
