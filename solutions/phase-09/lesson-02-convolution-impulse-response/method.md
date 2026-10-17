# P09.02  Method

The route to each answer, without the final numbers. `worked-solutions.md` has the numbers for s = 15.

## 1. Convolution by hand

(a) y[n] = Σ_k x[k]h[n − k] for n = 0..4: y[0] = x₀h₀; y[1] = x₀h₁ + x₁h₀; y[2] = x₀h₂ + x₁h₁ + x₂h₀; y[3] = x₁h₂ + x₂h₁; y[4] = x₂h₂. Check Σy against (Σx)(Σh). (b) c_L[n] = Σ_m y[n + mL] for n < L: for L = 3 add y[3], y[4] onto y[0], y[1]; for L ≥ 5 nothing folds and the extra samples are zero. (c) `convolve_direct(x, h)`, `convolve_fft(x, h)`, `circular_convolve(x, h, L)`.

## 2. Derivations

(a) Staircase x ≈ Σ_k x(τ_k) p_Δ(t − τ_k) with p_Δ a unit-area pulse of width Δ; time invariance makes each response x(τ_k)Δ h_Δ(t − τ_k); linearity adds them; Δ → 0 turns the sum into the integral and h_Δ into h. Causality restricts τ to [0, t]. (b) ∫₀ᵗ e^{−aτ}e^{−b(t − τ)}dτ = e^{−bt}∫₀ᵗ e^{(b − a)τ}dτ; for b = a the integrand is e^{−at}. Laplace: 1/((s + a)(s + b)) = [1/(s + a) − 1/(s + b)]/(b − a). (c) Impulse: τẏ + y = δ ⇒ y(0⁺) = 1/τ, then free decay. Step: particular solution 1, homogeneous e^{−t/τ}, y(0) = 0. Centroid: ∫₀^∞ t e^{−t/τ}dt/τ = τ, ∫₀^∞ e^{−t/τ}dt/τ = 1. (d) Trapezoid = Σ with weights ½, 1, …, 1, ½ and error O(dt²) for smooth integrands; the rectangle sum puts weight 1 on both end points, so it equals trapezoid + ½dt(f(0)g(t) + f(t)g(0)).

## 3. Pitot-static lag

(a) y = h * v with h = e^{−t/τ}/τ and v = 2t: y(t) = 2[t − τ(1 − e^{−t/τ})]. (b) For t ≫ τ, v − y → 2τ. (c) `convolve_continuous(v, h, 1e-3, "trap")` on t = arange(0, 12, 1e-3) with v = 2 min(t, 10); compare with (a) for t ≤ 10. (d) After the ramp stops the lag is 2τ(1 − e^{−10/τ}) e^{−(t − 10)/τ}; evaluate at t − 10 = 3τ. (e) Compare the lag with the tolerance of a speed callout.

## 4. FFT crossover

(a) Loop over M; `operation_counts(N, M)`; first M with fft < direct; confirm it stays cheaper for larger M in the range. (b) L = `next_power_of_two(N + M − 1)` is the same power of two for all M in the range while N + M − 1 stays below it. (c) `default_rng(s)`; `circular_convolve(x, h, 1000)` against `convolve_direct(x, h)[:1000]`; errors in the first M − 1 samples only. (d) `time.perf_counter()` around each call, several repeats; compare the time ratio with the count ratio and say why they differ (compiled FFT, memory traffic, Python loop overhead over M).

## 5. Judgment

(a) Largest rate of the linear model: A/τ at the first step, 2A/τ at the reversal; set 2A/τ equal to the limit. Above it the actuator ramps at its limit and lags the linear prediction. (b) dt/τ against 1 − e^{−dt/τ} at dt = 20 ms; the exact discretisation keeps the pole at e^{−dt/τ} and is exact for piecewise-constant inputs. (c) DC gain of the truncated record = 1 − e^{−T/τ}. (d) Oil or grease viscosity, motor winding resistance, battery voltage: τ and K change. Fit τ from the servo's response to commands in flight data and compare.

## 6. Step-invariant convolution

(a) On [t_j, t_{j+1}) the integrand is x_j h(t_k − τ); ∫ h(t_k − τ)dτ over that strip is s(t_k − t_j) − s(t_k − t_{j+1}); set m = k − j. (b) d = s on the grid, differenced, with d₀ = s(0) = 0; `convolve_direct(cmd, d)[:N]`; compare with the closed form and with both Riemann rules at three dt. (c) The doublet is piecewise constant on the grid (its jumps fall on grid points), so the assumption holds exactly; 2.5 is about rules that sample x and h at grid points, not about integrating h exactly. (d) Smooth inputs (a sinusoid): the formula holds x constant over each strip, a first-order error, while the trapezoid rule is second order. Compute both for a 2 Hz sine.
