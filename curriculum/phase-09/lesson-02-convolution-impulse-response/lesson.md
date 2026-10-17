# P09.02  Convolution & Impulse Response

> Hit a linear, time-invariant system with an impulse and write down what comes out. That one curve predicts its response to every other input: slide it, scale it, add it up.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the convolution sum and integral from superposition and time invariance, prove its algebraic properties, and compute short convolutions by hand. 2. Implement and test discrete convolution from the defining sum and by FFT multiplication, explain the wrap-around of circular convolution, and show that zero padding to N + M − 1 removes it. 3. Compute continuous convolutions by Riemann sums, measure their first- and second-order convergence against closed forms, and explain why a jump in the input keeps the error first order. 4. Derive the impulse and step responses of a first-order system and use them to predict the output of a gyro's RC front end and of a first-order elevator servo, including the mean delay τ. |
| **Prerequisites** | `P09.01` (LTI systems, impulse response, BIBO), `P02.08` (Laplace transform, step and impulse responses, the convolution formula), `P02.10` (DFT, FFT, the convolution theorem), `P04.10` (the RC low-pass), `P03.02` or `P00.07` (Riemann sums, trapezoid rule). Skill check: you can integrate e^{−at} over [0, t] and multiply two polynomials by hand. |
| **Engineering question** | A MEMS gyro feeds its rate signal through an RC low-pass with R = 1.6 kΩ and C = 1 μF, and an elevator servo behaves as K/(τs + 1) with τ = 50 ms. Given only their impulse responses, what do they do to a 10 deg/s pitch-rate step buried in a 1200 Hz vibration line (a SYNTHETIC record) and to a ±2 deg flight-test doublet? How should the convolution be computed: directly or by FFT, with what padding, with which integration rule, and how accurate is it? |
| **Assumptions and limits** | Linear, time-invariant systems starting from rest; causal signals (zero for t < 0); the RC filter is ideal (no source or load resistance, an ideal capacitor) and buffered; the servo is linear (no rate limit, no saturation, no backlash), which a real control-surface actuator violates at large commands. Continuous convolution is approximated on a uniform grid; its accuracy depends on the rule and on the smoothness of the signals. The gyro record is SYNTHETIC (seeded), not a measurement. |
| **Analytical work** | The convolution sum from x = Σ x[k] δ[n − k]; commutativity, associativity, distributivity, identity and delay by δ; the convolution integral as the limit of a staircase; (e^{−at}u) * (e^{−bt}u) and its a → b limit; rect * rect; the RC filter's ODE RC ẏ + y = x, its impulse and step responses; the mean delay ∫t h/∫h = τ; circular versus linear convolution and the wrap-around formula; operation counts; the leading error of the rectangle rule and why the trapezoid rule removes it; the sampled first-order recursion. |
| **Computational deliverable** | `starter/main.py`: `convolve_direct`, `convolve_fft`, `circular_convolve`, `causal_lti_output`, `convolve_continuous`, `exp_conv_exp`, `rect_conv_rect`, `first_order_impulse_response`, `first_order_step_response`, `step_from_impulse`, `impulse_from_step`, `impulse_response_centroid`, `sampled_first_order` (`next_power_of_two`, `operation_counts`, `rc_time_constant`, `doublet`, `doublet_response_closed_form`, `synthetic_gyro_record` and `main` are given). Tests prove: the direct sum equals FFT multiplication with zero padding and `numpy.convolve` for four length pairs; your FFT convolution equals the direct sum; unpadded circular convolution wraps the last M − 1 samples onto the first and padding removes it; the convolution algebra; the causal output equals the recursion; the Riemann sum converges to the closed form at order 1 (rectangle) and 2 (trapezoid); the a = b case is exact by trapezoid; rect * rect is a triangle; the step response by convolution matches the closed form and `scipy.signal.step`; step and impulse responses are integral and derivative of each other; the mean delay is τ; the sampled recursion has unit DC gain and its pole at e^{−dt/τ}; the doublet response by convolution matches superposed shifted steps and `scipy.signal.lsim` within the first-order bound 2A dt/τ. |
| **Practical activity** | **computational experiment**: direct against FFT convolution (accuracy, wrap-around, operation counts); convergence of the rectangle and trapezoid rules; the RC front end applied by convolution to a SYNTHETIC gyro record; the servo's doublet response by convolution against superposed steps. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `convolution.png`, the hand convolution of Problem 1, the derivations of the closed forms and of the wrap-around formula, the convergence table, the SYNTHETIC gyro check (labelled), and a short note recommending how to compute the servo's response in a test-planning tool. |
| **Assessment** | `assessments/phase-09/grading-rubric.md`; competencies 1–4. Pass: the convolution sum and the first-order responses are derived; all public tests pass; the convergence orders are measured and explained, including the first-order error at a jump; the SYNTHETIC record is labelled as such. |
| **References** | `@mit-ocw-18-03-differential-equations` L21 "Convolution Formula", L23 "Use with Impulse Inputs", notes 24 "Step response, impulse response" and 25 "Convolution" (required-free); `@mit-ocw-6-003-signals-and-systems` Lecture 20 "Applications of Fourier Transforms" (supplementary-free: convolution and filtering in the frequency domain); `@mit-ocw-18-06-linear-algebra` Lecture 26 "Complex matrices; fast fourier transform" (supplementary-free). GAP: circular convolution, zero padding and Riemann-sum error analysis are original instruction here. The Phase 09 reference map (`references/phase-09.md`) was not available when this lesson was written; the topics marked GAP are covered by original instruction in this lesson. |
| **Extensions** | ★ Implement overlap-add block convolution and show that it beats both the direct sum and a single FFT for a long record and a 65-tap filter. ★ Estimate an unknown impulse response from a SYNTHETIC noisy input-output record by least squares (deconvolution) and show how noise limits it. ★ Make `convolve_continuous` step-invariant for piecewise-constant inputs (Problem 6) and use it for the servo. |
| **Time** | unmeasured |

## 1. The engineering problem

`P09.01` showed that a linear, time-invariant (LTI) system is completely described by one signal, its **impulse response** h. That is a strong statement, and this lesson cashes it in. Two systems recur throughout this phase. The first is the analog anti-aliasing front end of a MEMS rate gyro: an RC low-pass with R = 1.6 kΩ and C = 1 μF, time constant 1.6 ms, corner 99.5 Hz, placed before the converter so that an engine or rotor vibration line at 1200 Hz does not alias into the pitch-rate signal. The second is the electro-mechanical servo that moves an elevator: for small commands it behaves as a first-order system K/(τs + 1) with τ = 50 ms. A flight-test engineer planning a doublet manoeuvre (+2 deg for 0.5 s, then −2 deg for 0.5 s) wants the surface deflection before the flight; a sensor engineer wants to know what the RC filter does to a pitch-rate step and to the vibration.

Both answers are a **convolution**: y = h * x. The computation looks trivial, and three things go wrong in practice. The direct sum costs N·M multiplications, which for a long record and a long impulse response is slow, so engineers use the FFT, and an FFT computes a *circular* convolution whose end wraps around onto its beginning unless the sequences are padded. Continuous convolution must be discretised, and the obvious rectangle rule is only first-order accurate; it also overestimates the DC gain by dt/(2τ), 1.6 % for the RC filter at 20 kHz, which is a calibration error, not a rounding error. And a sampled first-order system is not the continuous one sampled: its first sample is (1 − e^{−dt/τ}), not dt/τ.

Success criteria, written first: the direct sum equals FFT multiplication with zero padding to 10⁻¹¹; the Riemann sums converge to the closed forms at measured orders 1.0 and 2.0; the RC step response computed by convolution matches 1 − e^{−t/τ} to 10⁻⁴; the servo's doublet response computed by convolution matches the closed form within the first-order error bound; and on the SYNTHETIC gyro record the vibration amplitude is reduced by the predicted factor |H(1200 Hz)| = 0.083.

## 2. Understand and derive

### 2.1 The convolution sum

Any sequence is a sum of scaled, shifted impulses: x[n] = Σ_k x[k] δ[n − k]. Feed this to an LTI system S with h = S(δ). Time invariance gives S(δ[n − k]) = h[n − k]; linearity lets the sum and the scale factors through:

**y[n] = Σ_k x[k] h[n − k] = (x * h)[n].**   (1)

Each input sample launches a copy of h, scaled by x[k] and delayed by k, and the output is the sum of the copies. For sequences of lengths N and M starting at 0 the result has **N + M − 1** samples. Substituting m = n − k shows **x * h = h * x** (commutativity); grouping the triple sum shows associativity, so a cascade of LTI blocks has impulse response h₁ * h₂ whatever the order; the sum distributes over addition; **δ is the identity** and δ[n − d] delays by d. Summing (1) over n gives Σy = (Σx)(Σh): areas multiply, which is the statement that DC gains multiply in a cascade.

Two small hand examples fix the mechanics. [1, 2, 4] * [3, −1, 1]: slide, multiply, add, to get [3, 5, 11, −2, 4] (Problem 1). A moving average of M samples has h = [1/M, …, 1/M]; convolving it with a step gives a ramp that reaches 1 after M samples.

### 2.2 The convolution integral

Approximate a continuous input by a staircase of width Δ: the strip at τ_k is a pulse of area x(τ_k)Δ, which for small Δ acts like x(τ_k)Δ δ(t − τ_k) and produces x(τ_k)Δ h(t − τ_k). Summing and letting Δ → 0:

**y(t) = ∫ x(τ) h(t − τ) dτ,**   and for causal signals   **y(t) = ∫₀ᵗ x(τ) h(t − τ) dτ.**   (2)

The properties of 2.1 carry over. Two closed forms are the test oracles. For a, b > 0,

**(e^{−at}u) * (e^{−bt}u) = ∫₀ᵗ e^{−aτ} e^{−b(t−τ)} dτ = (e^{−at} − e^{−bt})/(b − a),**   (3)

whose limit as b → a is t e^{−at} (the integrand becomes e^{−at}, constant in τ). Two unit pulses of width w convolve to a **triangle** of height w on [0, 2w], area w².

### 2.3 The first-order system: an RC front end and a servo

The RC low-pass: the current through R charges C, (x − y)/R = C dy/dt, so

**τ dy/dt + y = x,  τ = RC,  H(s) = 1/(τs + 1).**   (4)

Its impulse response is the solution with y(0⁺) = 1/τ (an impulse of input charges the capacitor instantly to 1/τ): **h(t) = (1/τ) e^{−t/τ} u(t)**. Convolving h with a unit step, (2) gives the **step response s(t) = ∫₀ᵗ h = 1 − e^{−t/τ}**: the step response is the running integral of the impulse response, and h = ds/dt. A servo whose velocity loop makes the surface behave as K/(τs + 1) has h = (K/τ)e^{−t/τ} and s = K(1 − e^{−t/τ}).

Two numbers describe a low-pass impulse response. Its **area** ∫h = H(0) is the DC gain (1 for the RC, K for the servo). Its **centroid** ∫t h dt / ∫h dt = τ is the **mean delay**: slow signals come out delayed by τ, the group delay of `P04.10` seen in the time domain. For the gyro front end, τ = 1.6 ms of delay is the price of attenuating the 1200 Hz line to |H| = 1/√(1 + (2π·1200·0.0016)²) = 0.083.

A first-order system run in a computer at step dt is the recursion y[n] = a y[n−1] + (1 − a)K x[n] with **a = e^{−dt/τ}**, so h[n] = (1 − a)K aⁿ, unit DC gain times K, pole at a. For dt ≪ τ, 1 − a ≈ dt/τ and h[n] ≈ dt·h(n dt). For the servo at dt = 1 ms, h[0] = 0.01980 against dt·h(0) = 0.02000; at the 20 ms of a 50 Hz flight computer, 0.330 against 0.400 (Problem 5).

### 2.4 Convolution by FFT, and why it wraps around

The DFT turns convolution into multiplication (`P02.10`), but the convolution it multiplies is **circular** of length L: (x ⊛ h)[n] = Σ_k x[k] h[(n − k) mod L]. If the linear result y has more than L samples, its tail folds back:

**c[n] = Σ_m y[n + mL],  n = 0, …, L − 1.**   (5)

With L = N (no padding) and M < N, the last M − 1 samples of the linear result are added onto the first M − 1; the rest is correct. Padding both sequences with zeros to **L ≥ N + M − 1** makes the folded copies zero, and then circular = linear. A power of two is chosen for the FFT. In `main()`, for N = 1000 and M = 200, the padded FFT result (L = 2048) equals the direct sum to 3.2 × 10⁻¹⁴, while the unpadded circular result has errors up to 31.5 in its first 199 samples and 3.9 × 10⁻¹⁴ after them.

**Cost.** The direct sum costs N·M multiplications. The FFT route costs three length-L transforms at about (L/2) log₂ L complex multiplications each, plus L for the product, about 4 real multiplications per complex one. For N = 100 000 and M = 1000 that is 1.0 × 10⁸ against 1.4 × 10⁷; for N = 100 and M = 10 the FFT is six times *more* expensive (5.9 × 10³ against 10³). For a long record, what matters is the shorter length: with N = 5500, the FFT wins only for M ≥ 123 (Problem 4). Short filters are applied directly, long ones by FFT (in blocks, extension 1).

### 2.5 Riemann sums for the integral

On a grid t_k = k dt, the rectangle rule gives y(t_k) ≈ dt Σ_{j=0}^{k} x_j h_{k−j} = dt (x * h)[k], the discrete convolution scaled by dt. It counts both end points of [0, t_k] with full weight, so its leading error is ½dt [x(0)h(t) + x(t)h(0)]: **first order in dt**. The trapezoid rule subtracts exactly that, leaving an O(dt²) error for smooth integrands:

y_trap(t_k) = dt[(x * h)[k] − ½(x₀ h_k + x_k h₀)].   (6)

On (e^{−t}u) * (e^{−3t}u) the code's maximum error halves with dt for the rectangle (0.100, 0.050, 0.025, 0.0125 at dt = 0.1, 0.05, 0.025, 0.0125 s) and quarters for the trapezoid (6.4 × 10⁻⁴, 1.6 × 10⁻⁴, 4.0 × 10⁻⁵, 1.0 × 10⁻⁵). The rectangle rule also biases the DC gain: dt Σ h_k = (dt/τ)/(1 − e^{−dt/τ}) ≈ 1 + dt/(2τ), 1.6 % high for the RC at dt = 50 μs. The same bias shows in the sampled centroid, τ − dt/2 = 1.575 ms instead of 1.6 ms.

The trapezoid rule is second order only if the integrand is smooth. A doublet **jumps**, and inside the strip that contains a jump neither rule knows where the jump is: the error is first order again, bounded by about A·dt·max|h| per jump, here 2A dt/τ. For the 2 deg doublet at dt = 1 ms the code's maximum error is 0.060 deg with the rectangle rule and 0.040 deg with the trapezoid, halving to 0.020 deg at 0.5 ms. For a piecewise-constant input there is an exact alternative: integrate h analytically over each strip, which turns the sum into one over increments of the step response (the **step-invariant** convolution, Problem 6).

### 2.6 Where the model stops being valid

Convolution with a fixed h is exact only for LTI systems. The servo has a rate limit and a deflection limit (typically ±25 to ±30 deg for an elevator), so its response to a large doublet is not the scaled response to a small one; above the rate limit it ramps instead of following an exponential. The RC filter's corner moves with component tolerance (±10–20 % for many capacitors) and is loaded by the converter's input. An impulse response measured on a test rig is noisy and truncated; a 2 s record of a system with τ = 0.5 s misses 1.8 % of the area, a 1.8 % DC-gain error. And the SYNTHETIC gyro record is a single tone plus white noise; real vibration is a set of harmonics on a broadband floor.

## 3. Build the computational model

Open `starter/main.py`. `next_power_of_two`, `operation_counts`, `rc_time_constant`, `doublet`, `doublet_response_closed_form`, `synthetic_gyro_record` and `main` are given.

### Step 1: the direct sum

Loop over the shorter sequence and add shifted, scaled copies of the longer one (2.1):

```python
def convolve_direct(x, h):
    if len(x) < len(h):
        x, h = h, x
    y = np.zeros(len(x) + len(h) - 1)
    for k, hk in enumerate(h):
        y[k:k + len(x)] += hk * x
    return y
```

`causal_lti_output(h, x)` keeps the first len(x) samples.

### Step 2: FFT multiplication and the circular convolution

`convolve_fft` pads to `L = next_power_of_two(N + M - 1)`, multiplies `np.fft.fft(x, L) * np.fft.fft(h, L)`, inverts, keeps N + M − 1 samples and returns the real part for real inputs. `circular_convolve(x, h, L)` implements the mod-L sum from its definition (add `np.roll(h_padded, k)` scaled by x[k]); it is the slow reference that shows what an unpadded FFT computes.

### Step 3: continuous convolution and the closed forms

`convolve_continuous(f, g, dt, rule)` is `dt * convolve_direct(f, g)[:N]` for the rectangle rule and (6) for the trapezoid. `exp_conv_exp` is (3) with the a = b limit; `rect_conv_rect` is the triangle.

### Step 4: first-order responses

`first_order_impulse_response(t, K, tau)` and `first_order_step_response` are the formulas of 2.3; `step_from_impulse` is the cumulative trapezoid integral and `impulse_from_step` the central difference; `impulse_response_centroid` is Σ t h / Σ h; `sampled_first_order` returns (1 − a)K aⁿ.

### Verify

Run `python -m pytest curriculum/phase-09/lesson-02-convolution-impulse-response`. The tests:

1. **Direct sum against FFT multiplication with zero padding** (written in the test with `numpy.fft`, independent of your code) and against `numpy.convolve`, for lengths (1, 1), (7, 3), (64, 64), (500, 37), to 10⁻¹¹.
2. **Your FFT convolution** equals your direct sum, real-valued, for three length pairs.
3. **Wrap-around**: unpadded circular convolution equals the DFT product, its first M − 1 samples equal y[n] + y[n + N] as (5) says, the rest is exact; padded to 61 or 64 it equals the linear result.
4. **Algebra**: commutative, associative, distributive, δ is the identity, a shifted δ delays, areas multiply.
5. **Causal output** equals the recursion computed by `scipy.signal.lfilter`.
6. **Convergence**: measured orders 1.0 ± 0.1 (rectangle) and 2.0 ± 0.1 (trapezoid) for three (a, b) pairs; trapezoid exact when a = b; the closed form is symmetric in (a, b) and continuous at a = b.
7. **Triangle**: rect * rect within 2 dt of the closed form, peak w at t = w, area w².
8. **First order**: the step by convolution matches 1 − e^{−t/τ} to 10⁻³K and `scipy.signal.step` to 10⁻⁹; area of h equals K; integral and derivative relations; centroid τ within dt; the sampled recursion has unit DC gain times K and ratio e^{−dt/τ}, and approaches dt·h(t) for dt ≪ τ.
9. **Doublet**: convolution against superposed steps and against `scipy.signal.lsim` within 2A dt/τ.
10. **main** writes the figure.

### Compare

`numpy.convolve` and `scipy.signal.fftconvolve` are the production tools; the first is the direct sum in compiled code, the second the padded FFT of Step 2 (with real-input transforms and a choice of fast lengths). Our direct sum agrees with both to rounding; our FFT route agrees with the direct sum to 10⁻¹⁰ and is faster only when the shorter sequence is long (2.4). `scipy.signal.lsim` does not convolve at all: it integrates the state-space form of the system, assuming the input varies linearly between samples, which is why it too differs from the exact doublet response by a first-order amount at the jumps (test 9). `scipy.signal.step` evaluates the exact step response and agrees with 1 − e^{−t/τ} to 10⁻⁹. The value of the from-scratch code is that every one of these differences can now be predicted before it is seen.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment** with a SYNTHETIC record. Write each prediction first.

1. **Direct, FFT, circular.** For seeded random x (N = 1000) and h (M = 200), compute the direct sum, the padded FFT product and the unpadded circular convolution (L = 1000). Predict where the circular result is wrong and by how much (5). Tabulate `operation_counts` for (N, M) = (100, 10), (1000, 100), (10⁴, 10³), (10⁵, 10³) and predict which method wins in each case.
2. **Convergence.** Compute (e^{−t}u) * (e^{−3t}u) on [0, 8] s by both rules for dt = 0.1, 0.05, 0.025, 0.0125 s; record the maximum error against (3) and the observed order. Predict the orders before running.
3. **RC front end on a SYNTHETIC record.** `synthetic_gyro_record()` returns 0.2 s at 20 kHz of a 10 deg/s rate step at 20 ms, a 3 deg/s line at 1200 Hz and 0.2 deg/s white noise (seed 9). Convolve it with h of the RC filter by the trapezoid rule. Predict the vibration amplitude after the filter (3 × 0.083) and the mean rate after the step. Then repeat with the rectangle rule and predict the bias in the mean.
4. **Servo doublet.** Convolve the ±2 deg doublet (0.5 s per half, starting at 0.1 s) with the servo's h at dt = 2, 1 and 0.5 ms; compare with the sum of shifted step responses. Predict the error and how it scales with dt.

Raw data format: CSV tables `N, M, direct_mults, fft_mults` (1); `dt_s, err_rect, err_trap` (2); `rule, vib_amp_before, vib_amp_after, mean_after_step` labelled SYNTHETIC in its header (3); `dt_s, err_rect_deg, err_trap_deg` (4).

What you should find. (1) The circular result is wrong only in its first 199 samples (errors up to 31.5) and exact to 4 × 10⁻¹⁴ after them; the FFT is 6 times more expensive at (100, 10), 1.4 times at (1000, 100), and 7 times cheaper at (10⁴, 10³) and (10⁵, 10³). (2) Orders 1.0 and 2.0 exactly as in 2.5. (3) On the SYNTHETIC record the vibration-plus-noise amplitude falls from 3.01 to 0.25 deg/s, as |H(1200 Hz)| = 0.083 predicts, and the mean rate after the step is 10.000 deg/s with the trapezoid rule; with the rectangle rule it reads 10.156 deg/s, the 1.6 % DC-gain bias of 2.5 appearing as a scale-factor error. (4) The doublet response peaks at 1.9999 deg at t = 0.600 s (the end of the first half, ten time constants after the start); the errors are first order, 0.080, 0.040 and 0.020 deg for the trapezoid rule at 2, 1 and 0.5 ms (0.060 deg for the rectangle at 1 ms), matching the bound 2A dt/τ.

## 5. Validate and improve

| quantity | theory (§2) | code | independent check | difference | cause |
|---|---|---|---|---|---|
| direct sum vs padded FFT | equal | 3.2 × 10⁻¹⁴ | `numpy.fft` in the test | — | rounding |
| unpadded circular, first M − 1 samples | y[n] + y[n + N] | max error 31.5 | (5) | exact as predicted | wrap-around |
| rectangle rule order | 1 | 1.00 | (3) | — | end-point weights |
| trapezoid rule order | 2 | 2.00 | (3) | — | — |
| RC step, dt = 50 μs | 1 − e^{−t/τ} | 8.1 × 10⁻⁵ max error | `scipy.signal.step` | — | trapezoid O((dt/τ)²) |
| mean delay of h | τ = 1.600 ms | 1.575 ms (sampled) | τ − dt/2 | 0 | rectangle centroid |
| 1200 Hz line after RC (SYNTHETIC) | 3 × 0.083 = 0.25 | 0.25 deg/s | — | < 0.01 | noise |
| doublet, trapezoid, 1 ms | error ≤ 2A dt/τ = 0.08 | 0.040 deg | superposed steps | within bound | jump inside a strip |
| sampled servo h[0], 1 ms | 1 − e^{−0.02} = 0.01980 | 0.01980 | dt h(0) = 0.0200 | 1 % | sampling |

Everything agrees with the prediction, including the two places where the obvious method is wrong (the unpadded FFT, the rectangle rule's DC bias). These are comparisons between derivations and computations; the SYNTHETIC record tests the pipeline, not a sensor. Improvements worth making: block convolution (overlap-add) for long records with medium filters, where neither a single FFT nor the direct sum is best; and a step-invariant convolution for piecewise-constant commands, which is exact at every grid point (Problem 6). Improvement *not* worth making: a higher-order quadrature rule for the doublet; the error is set by the jumps, not by the rule.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand convolution of Problem 1; the derivations of (1), (3), (4), (5) and of the rectangle rule's leading error; `convolution.png`; the four experiment tables (the gyro one labelled SYNTHETIC); and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending how a flight-test planning tool should compute the servo's response to a pilot-input sequence, with the accuracy you expect and the conditions under which the linear model fails.

### Adversarial self-checks

- **What if the record is ten times longer?** The direct sum costs ten times more; the padded FFT about ten times more plus a log factor; the circular error without padding is unchanged in size but now a smaller fraction of the record. Nothing about the physics changes.
- **What if τ → 0?** h → δ: the filter passes everything, the delay vanishes, and on a fixed grid the Riemann sums break down once dt is comparable to τ (the rectangle rule's DC bias dt/(2τ) grows without bound). What if τ → ∞? The low-pass becomes an integrator scaled by 1/τ: bounded only for a finite record.
- **Which assumption fails first in the real world?** The servo's rate limit: a 2 deg step in 50 ms needs about 40 deg/s initially; a 20 deg doublet needs 400 deg/s, more than many hobby and light-aircraft actuators deliver, and the convolution prediction is then wrong.
- **What measurement could falsify the model?** Command two doublets of different amplitude on a servo test rig (`P09.21`): if the responses do not scale, the system is not linear; if a delayed command gives a different shape, it is not time-invariant.
- **Does the result survive a change of grid or rule?** For smooth inputs, the trapezoid rule converges at second order and both rules agree as dt → 0; for the doublet, every rule is first order and the remedy is the step-invariant formula, not a smaller dt.
- **Which decision changes if the capacitor tolerance doubles?** With C 20 % high, τ = 1.92 ms: the 1200 Hz attenuation improves slightly to |H| = 0.069 but the mean delay grows to 1.92 ms. For the pitch loop the added delay matters more than the extra attenuation (`P09.10` will put a number on it).

## Check yourself

1. Convolve [1, 1, 1] with [1, 2] by hand. How many samples does the result have, and what is its sum?
2. Why does the response of an LTI system to any input follow from its impulse response? Name the two properties used.
3. What is the impulse response of an RC low-pass with τ = 2 ms, its step response, and the mean delay it adds to slow signals?
4. You multiply the 512-point FFTs of a 512-sample record and a 64-tap filter. Which output samples are wrong, and how do you fix it?
5. The rectangle-rule convolution of a first-order impulse response with a step settles at 1.0156 instead of 1. Which ratio sets the error, and what rule removes it?
6. Why is the trapezoid rule only first-order accurate for a doublet input?
7. A first-order servo with τ = 50 ms is simulated by y[n] = a y[n−1] + (1 − a)x[n] at 50 Hz. What is a, and what is the DC gain?
