# P09.02  Worked solutions

Worked for **s = 15**. Every number below is produced by `solution.py` (its `main()` or the functions named) or by a probe script that calls them; the "independent check" line says how it was verified by a different route. Every simulated record is SYNTHETIC.

## 1. Convolution by hand (x = [1, 2, 4], h = [3, −1, 1])

**Result.** y[0] = 1·3 = 3; y[1] = 1·(−1) + 2·3 = 5; y[2] = 1·1 + 2·(−1) + 4·3 = 11; y[3] = 2·1 + 4·(−1) = −2; y[4] = 4·1 = 4. **y = [3, 5, 11, −2, 4]**, Σy = 21 = (Σx)(Σh) = 7 × 3. `convolve_direct`, `convolve_fft` and `numpy.convolve` all return [3, 5, 11, −2, 4].

(b) **L = 3**: [y₀ + y₃, y₁ + y₄, y₂] = **[1, 9, 11]**: the last two linear samples have wrapped onto the first two. **L = 5** = N + M − 1: nothing folds, **[3, 5, 11, −2, 4]**. **L = 8**: **[3, 5, 11, −2, 4, 0, 0, 0]**, the padded linear result. `circular_convolve` gives the same three sequences.

**Independent check.** The sum of every circular result is also 21: folding moves samples, it does not create or destroy area.

**Common mistakes.** Not reversing h; giving a 3-sample answer (the result has N + M − 1 = 5); believing the FFT computes the linear convolution without padding.

## 2. Derivations

**(a)** Approximate x by Σ_k x(τ_k) p_Δ(t − τ_k), where p_Δ is a pulse of width Δ and height 1/Δ. Time invariance: the response to p_Δ(t − τ_k) is h_Δ(t − τ_k), where h_Δ is the response to p_Δ. Linearity: the response to the staircase is Σ_k x(τ_k)Δ h_Δ(t − τ_k). As Δ → 0, p_Δ → δ, h_Δ → h and the sum becomes ∫x(τ)h(t − τ)dτ. For a causal system h(t − τ) = 0 for τ > t, and for a causal input x(τ) = 0 for τ < 0, so the limits are 0 and t.
**(b)** ∫₀ᵗ e^{−aτ}e^{−b(t−τ)}dτ = e^{−bt}[e^{(b−a)τ}/(b − a)]₀ᵗ = (e^{−at} − e^{−bt})/(b − a). For b = a the integrand is e^{−at} for every τ, so the integral is t e^{−at}; the formula's limit (by l'Hôpital in b) is the same. Laplace: the transforms are 1/(s + a) and 1/(s + b), the product 1/((s + a)(s + b)) = [1/(s + a) − 1/(s + b)]/(b − a), whose inverse is (3).
**(c)** τẏ + y = δ(t): integrating across t = 0 gives τ y(0⁺) = 1, then ẏ = −y/τ, so h = (1/τ)e^{−t/τ}u(t). Step: y = 1 − e^{−t/τ}. Centroid: ∫₀^∞ (t/τ)e^{−t/τ}dt = τ (by parts) and ∫₀^∞ (1/τ)e^{−t/τ}dt = 1, ratio **τ**.
**(d)** The trapezoid rule on [0, t_k] with k strips is dt[½f₀g_k + f₁g_{k−1} + … + ½f_kg₀] and has error O(dt²) for smooth integrands. The rectangle sum dt Σ_{j=0}^{k} f_j g_{k−j} equals the trapezoid sum plus ½dt(f₀g_k + f_kg₀). Hence the rectangle error is ½dt[f(0)g(t_k) + f(t_k)g(0)] + O(dt²), first order, and subtracting it, (6), leaves O(dt²). At t = 0 the rectangle sum is dt f₀g₀ against the exact 0: the code's maximum rectangle error for (e^{−t}u)*(e^{−3t}u) is exactly dt (0.1000 at dt = 0.1).

## 3. Pitot-static lag (τ = 0.35 s)

**Assumptions.** First-order pneumatic lag, linear; airspeed v(t) = 2t m/s for 0 ≤ t ≤ 10 s, then 20 m/s.

**(a)** y(t) = ∫₀ᵗ 2τ′ (1/τ)e^{−(t−τ′)/τ}dτ′ = **2[t − τ(1 − e^{−t/τ})]** m/s. **(b)** For t ≫ τ, v − y → 2τ = **0.70 m/s**: a first-order lag behind a ramp is a constant delay τ, so the indicated speed is the true speed of τ seconds ago. **(c)** `convolve_continuous` with the trapezoid rule at dt = 1 ms differs from (a) by at most **1.4 × 10⁻⁵ m/s**; at t = 10 s the lag is 0.69999 m/s. **(d)** After the ramp stops, the lag decays as 2τ(1 − e^{−10/τ}) e^{−(t − 10)/τ}; at 3τ it is 0.70 e^{−3} = **0.035 m/s** (the code gives 0.0348).
**(e)** A callout based on the indicated speed comes about τ = 0.35 s (0.7 m/s) late during the acceleration: small against the tolerance of a hobby UAV's take-off logic but not negligible for a high-acceleration launch, where the lag is rate × τ.

**Independent check.** The mean-delay result of 2.3: the centroid of h is τ, and for a ramp the whole output is delayed by exactly that centroid.

**Common mistakes.** Computing the lag as τ in m/s (it is rate × τ); forgetting the transient term τe^{−t/τ} near t = 0.

## 4. When does the FFT pay? (N = 5500)

| M | 8 | 16 | 32 | 64 | 128 | 256 | 512 | 1024 |
|---|---|---|---|---|---|---|---|---|
| direct (N·M) | 4.40 × 10⁴ | 8.80 × 10⁴ | 1.76 × 10⁵ | 3.52 × 10⁵ | 7.04 × 10⁵ | 1.41 × 10⁶ | 2.82 × 10⁶ | 5.63 × 10⁶ |
| FFT | 6.72 × 10⁵ | 6.72 × 10⁵ | 6.72 × 10⁵ | 6.72 × 10⁵ | 6.72 × 10⁵ | 6.72 × 10⁵ | 6.72 × 10⁵ | 6.72 × 10⁵ |

**(a)** The FFT route is cheaper from **M = 123** on (and stays cheaper for every larger M up to 2000). **(b)** For every M in the table N + M − 1 lies between 4097 and 8192, so L = 8192 and the count 4[3(L/2)log₂L + L] = 671 744 is fixed; the direct count grows linearly in M. **(c)** With M = 65 and L = 1000 (no padding), the largest error is **12.1**, in the first M − 1 = 64 samples (the wrap-around region); elsewhere the error is 1.2 × 10⁻¹⁴. **(d)** Timings are machine-specific and are not reproduced here; expect the compiled FFT to win by more than the counts suggest at (10⁵, 2000) and the direct loop to win at (10⁵, 10), where it is ten vectorised array additions.

**Common mistakes.** Treating the operation count as a timing; forgetting that the FFT length jumps by a factor of two at a power-of-two boundary.

## 5. Judgment

**(a)** The linear response to a doublet of amplitude A has its fastest rate at the reversal, where the target jumps by 2A from a deflection close to +A: 2A/τ. For the lesson's ±2 deg doublet that is 80 deg/s (the code's closed-form response peaks at 79.9 deg/s). With a 60 deg/s rate limit the prediction is valid only up to **A = 60 × 0.05/2 = 1.5 deg**; the lesson's 2 deg doublet already exceeds it at the reversal, where the real surface ramps at 60 deg/s and lags the linear prediction. Above the limit the system is not linear, and convolution with h does not apply.
**(b)** At dt = 20 ms: dt/τ = **0.400**, against 1 − e^{−0.4} = **0.330**. The colleague's recursion is forward Euler; its pole 1 − dt/τ = 0.600 differs from e^{−dt/τ} = 0.670, so the simulated servo is faster than the real one, and at dt/τ > 2 it would be unstable. Fly the exact recursion with a = e^{−dt/τ}: it is exact for a command that is held constant between updates, which is what a flight computer does.
**(c)** 1 − e^{−2/0.5} = 0.9817: the truncated record has a **1.8 %** DC-gain error, which appears as a scale-factor error on every steady response predicted from it.
**(d)** Viscosity of grease or hydraulic fluid falls with temperature, the motor winding resistance rises, the battery voltage changes: τ and K drift by tens of percent. Fit τ and K to the commanded and measured surface positions in flight (the data of `P09.21` is the home-lab version) and compare with the ground value; a residual that grows with temperature is the signature.

## 6. Exact convolution for piecewise-constant commands

**(a)** On [t_j, t_{j+1}) the input is x_j, so ∫ x(τ)h(t_k − τ)dτ over that strip = x_j ∫_{t_j}^{t_{j+1}} h(t_k − τ)dτ = x_j[s(t_k − t_j) − s(t_k − t_{j+1})], using s′ = h. Summing over the strips and setting m = k − j, y(t_k) = Σ_j x_j d_{k−j} with **d_m = s(m dt) − s((m − 1)dt)** (d₀ = s(0) − s(−dt) = 0): a discrete convolution of x with the step-response increments.

**(b)** Maximum error against the closed form (deg):

| dt | rectangle | trapezoid | step-invariant |
|---|---|---|---|
| 2 ms | 0.120 | 0.080 | 2.2 × 10⁻¹⁵ |
| 1 ms | 0.060 | 0.040 | 4.4 × 10⁻¹⁵ |
| 0.5 ms | 0.030 | 0.020 | 6.4 × 10⁻¹⁵ |

**(c)** The doublet's jumps fall on grid points, so the input really is constant on every strip and the formula integrates h exactly: no approximation is made. This does not contradict 2.5, which concerns rules that *sample* x and h; here h is integrated analytically, which is the time-domain form of the zero-order-hold discretisation of `P02.05` and `P09.17`.
**(d)** For smooth inputs the formula's assumption (x constant over each strip) is itself a first-order error. For a 2 Hz sine through the same servo at dt = 1 ms, the step-invariant result is off by **5.4 × 10⁻³** while the trapezoid rule is off by **3.9 × 10⁻⁵**. Use the step-invariant form for commands that are genuinely held (digital commands, doublets on the grid) and the trapezoid rule for smooth signals.

**Common mistakes.** Off-by-one in the increments (shifting d by one sample doubles the error instead of removing it); applying the formula to a sampled smooth signal and expecting exactness.

## Check-yourself answers (lesson.md)

1. [1, 3, 3, 2]: 3 + 2 − 1 = 4 samples, sum 9 = 3 × 3.
2. Any input is a superposition of shifted, scaled impulses; time invariance makes each response a shifted copy of h, and linearity lets the responses add with the same scale factors.
3. h(t) = 500 e^{−500t} u(t) [1/s], s(t) = 1 − e^{−500t}, mean delay τ = 2 ms.
4. The first 63 output samples are corrupted by wrap-around of the last 63 linear samples (and the last 63 linear samples are missing). Pad both to at least 512 + 64 − 1 = 575 (in practice 1024).
5. The ratio dt/τ: the rectangle rule's DC gain is (dt/τ)/(1 − e^{−dt/τ}) ≈ 1 + dt/(2τ), here dt/(2τ) = 0.0156. The trapezoid rule removes it (to O((dt/τ)²)).
6. The integrand is discontinuous inside the strips that contain a jump of the input, and neither rule knows where in the strip the jump lies; the local error there is O(dt) and does not average out.
7. a = e^{−0.02/0.05} = 0.670; DC gain (1 − a)/(1 − a) = 1.
