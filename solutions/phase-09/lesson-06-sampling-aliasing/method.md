# P09.06  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. f₁ = 355 + 2s; compute k f₁ for k = 1, 2, 3. 2. `aliased_frequency(k*f1, 400)`. 3. Compare with 30 Hz. 4. `n = np.arange(400); dominant_frequency(np.cos(2*np.pi*f*n/400), 400)` for each f.

## 2.
(a) \|H(jω)\|² = 1/(1 + (ω²/ω_c²)ⁿ); with ω² = −s² the denominator of H(s)H(−s) is 1 + (−1)ⁿ(s/ω_c)²ⁿ; set it to zero: (s/ω_c)²ⁿ = (−1)^{n+1}, so s = ω_c e^{jθ} with 2nθ = (n + 1)π + 2πm. List θ, keep the n with cos θ < 0, re-index to k = 1..n. (b) From (4): 2n log(f_c/f_p) ≥ log ε_p⁻¹... rearrange each to a bound on log f_c; the lower bound must not exceed the upper bound; solve for n. (c) log H(jω) = log N(jω) − log D(jω); N(jω) = b₀(1 + jωb₁/b₀ + O(ω²)), log(1 + x) ≈ x; the imaginary part is the phase. For a(s) = Π(s − p_k): a₀ = Π(−p_k), a₁ = a₀ Σ_k 1/(−p_k). Evaluate with the pole lists for n = 1, 2, 3 (cos 60° = ½ for n = 3).

## 3.
`antialias_spec(100, fs, 1.0, 1e-3, 1.0)` for fs = 1000 and 4000 (the amplitudes enter only through their ratio). Build `butterworth_analog(n, 2*np.pi*fc)`; `group_delay_dc`; phase `np.angle(freq_response(b, a, 2*np.pi*20))`. Compare and weigh delay against data rate.

## 4.
`antialias_spec(30, 400, 5, 0.02, 0.5)` for the corner range. For each corner, `imu_experiment(orders=(3,), fc=fc, seed=s)`; then `t = np.arange(25000)/20000`, `_, v = synthetic_imu_record(t, seed=s)`, `ts, ys = sample_through_filter(t, v, 400, b, a)`, `tone_amplitude(ys[ts >= 0.25], 400, 20)`; same without b, a. Prediction 2·\|H(2π·780)\|. Estimate the single-bin noise: the filtered noise is white noise of σ = 0.05 deg/s through a filter of noise bandwidth about (π/2)f_c, sampled; one bin of 400 samples carries about 2σ_f/√400.

## 5.
Loop N ∈ {50, 200, 800}, f ∈ {5, 12.5, 22.5}: `n = np.arange(N)`, `x = np.cos(2*np.pi*f*n/50 + 0.3)`, `tt = np.linspace(0.25*N/50, 0.75*N/50, 501)`, maximum of \|`sinc_reconstruct(x, 50, tt)` − cos(2πf tt + 0.3)\|.

## 6.
(a) `aliased_frequency` over 351–429 Hz (note 400 Hz itself); 5\|H(2πf)\| at the ends and the middle; the minimum attenuation over the range is at its lowest frequency. (b) After sampling, alias and true motion have identical samples; the line moves. (c) `freq_response` phase at 1, 2, 5, 10 Hz divided by ω: nearly constant, the DC delay 2/ω_c; delay scales as 1/ω_c and with n. (d) `antialias_spec(30, 2000, 5, 0.02, 0.5)` and its delay and 5 Hz phase.
