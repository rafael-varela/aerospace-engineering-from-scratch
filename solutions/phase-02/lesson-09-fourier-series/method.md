# P02.09  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) a₀ = (2/T)∫₀^{dT} dt = 2d; a_n = (2/T)∫₀^{dT} cos(nω₀t)dt = sin(2πnd)/(πn); b_n = (1 − cos(2πnd))/(πn); amplitude² = [sin² + (1 − cos)²]/(πn)² = 2(1 − cos 2πnd)/(πn)² = 4 sin²(πnd)/(πn)². (b) `fourier_coefficients(lambda t: np.where(np.mod(t, 1.0) < d, 1.0, 0.0), 1.0, 6, n_samples=8192)`. (c) `parseval_power` with 6 and with 200 harmonics against d. (d) Harmonic n vanishes when nd is an integer; the second vanishes for d = 0.5.

## 2.
(a) Product-to-sum; integrate; the cos((m − n)ω₀t) term gives T/2 only when m = n. (b) Multiply (1) by cos(nω₀t) or sin(nω₀t), integrate over a period, use (a). (c) Split the integral at T/2; each piece gives (A/(nπ))(1 − cos nπ) with the same sign. (d) Square (1), integrate, keep only the diagonal terms; for the square wave the left side is A², the right side (8A²/π²)Σ_odd 1/n²; then Σ_all = Σ_odd + ¼Σ_all.

## 3.
(a) Period 1/(2f₀); v(−t) = v(t); with ω₀ = 4πf₀ and the substitution θ = 2πf₀t the coefficient is (4/π)∫₀^{π/2}... carry the integral ∫₀^π sin θ cos(2nθ) dθ = 2/(1 − 4n²). (b) `fourier_coefficients(lambda t: np.abs(np.sin(2*np.pi*400*t)), 1/800, 6)`. (c) Lowest ripple harmonic is n = 1 at 800 Hz; |a₁|/(a₀/2) = 2/3; rms ripple² = ½Σ a_n² = mean square − (a₀/2)² = ½ − 4/π². (d) Attenuation needed = 0.01/(2/3); solve 1/√(1 + (ωτ)²) = that at ω = 2π × 800; higher harmonics are both smaller and more attenuated.

## 4.
(a) 3/T, 5/T, 7/T Hz; `response_coefficients(*square_wave_coefficients(51), T, 2*pi*5, 0.01)`; `dominant_harmonic`; `harmonic_amplitudes(An, Bn) * wn**2`; peak via `partial_sum` on a fine grid over one period. (b) Loop over `np.arange(0.90, 1.301, 0.01)`. (c) `solve_ivp(forced_mode_rhs, (0, 40), [0, 0], args=(f, wn, zeta), max_step=2e-3, rtol=1e-9)`; compare on t ≥ 35. (d) Resonant T for harmonic n: n/f with f in [4.8, 5.2]: T ∈ [n/5.2, n/4.8]; your T is safe if outside all three bands.

## 5.
(a) Pointwise convergence away from the jump; the overshoot's width → 0; the torque is the square wave itself. (b) Coefficients exact; f_mode ± few %; ζ ± factor 2 (response ∝ 1/ζ); rise time reduces high harmonics (5th barely). (c) n ≈ T/(π t_r) ≈ 8; the 5th harmonic is below the knee, so the response changes by a few per cent at most. (d) Compare 20 s with the rise time of the resonant response 1/(ζω_n) = 3.2 s: the amplitude reaches 1 − e^{−20/3.2} of the steady value.

## 6.
(a) Σ_{k odd ≤ N} e^{ikx} = e^{ix}(1 − e^{i(N+1)x})/(1 − e^{2ix}); multiply numerator and denominator to extract sines; take the real part. (b) S_N′ = 0 first when (N + 1)ω₀t = π. (c) S_N(t*) = ∫₀^{t*} S_N′ dt; substitute; the integrand becomes (2/π) sin u/u du. (d) Tail Σ_{odd n>N} 1/n² ≈ ∫_N^∞ dn/(2n²) = 1/(2N); rms = √(½ × 16/π² × 1/(2N)); the overshoot region has width ∝ 1/N and bounded height, so its contribution to the mean square is O(1/N), while... compare its order with the tail's and note the mean-square error already includes it.
