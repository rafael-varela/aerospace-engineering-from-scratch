# P09.02  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Any simulated record you generate is SYNTHETIC: label it so in every table and plot.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: convolution by hand)

Let x = [1, 2, 1 + (s mod 4)] and h = [3, −1, 1], both starting at n = 0.

(a) Compute y = x * h by hand (sliding and adding) and check that Σy = (Σx)(Σh). (b) Compute the circular convolutions of length L = 3, 5 and 8 by hand and relate each to y with the wrap-around formula c[n] = Σ_m y[n + mL]. (c) Check (a) and (b) with `convolve_direct`, `convolve_fft` and `circular_convolve`.

Answer form: a five-sample sequence with its sum, three circular results with one sentence each.

## 2. (Derivation: convolution and first-order systems)

(a) Starting from the staircase approximation of an input, derive the convolution integral (2) for a causal LTI system, and say where linearity and time invariance enter. (b) Derive (3), (e^{−at}u) * (e^{−bt}u) = (e^{−at} − e^{−bt})/(b − a), and its limit as b → a, once by direct integration and once from the Laplace transform (`P02.08`). (c) From τẏ + y = x derive h(t) and s(t), and show that ∫₀^∞ t h dt / ∫₀^∞ h dt = τ. (d) Show that the rectangle-rule convolution dt Σ_{j=0}^{k} x_j h_{k−j} differs from the integral by ½dt[x(0)h(t_k) + x(t_k)h(0)] + O(dt²) for smooth signals, and hence that (6) is second order.

Answer form: four derivations.

## 3. (Unfamiliar scenario: pneumatic lag in a pitot-static line)

The pressure line from a pitot probe to an air-data sensor behaves as a first-order lag with τ = 0.2 + s/100 s. During a take-off roll the true airspeed rises at 2 m/s per second for 10 s, then stays constant.

(a) Write the indicated airspeed as a convolution and evaluate it in closed form for 0 ≤ t ≤ 10 s. (b) Show that once the transient has decayed the indicated speed lags the true speed by a constant amount, and give it. (c) Compute the same response with `convolve_continuous` (trapezoid, dt = 1 ms) and report the maximum difference from (a). (d) How large is the lag 3τ after the acceleration stops? (e) In one sentence: why does this matter for a take-off-speed callout?

Answer form: a formula, a speed lag in m/s, a numerical error, a speed, one sentence.

## 4. (Numerical, uses your code: when does the FFT pay?)

(a) For N = 4000 + 100s samples and filter lengths M = 8, 16, …, 1024, tabulate `operation_counts(N, M)` and find the smallest M for which the FFT route needs fewer multiplications. (b) Explain from the formulas why the FFT count does not change with M over this range. (c) With a seeded random record of 1000 samples and a random filter of M = 50 + s taps, compute the unpadded circular convolution of length 1000 and report the largest error against the linear result, in which samples it occurs, and the error elsewhere. (d) Time `convolve_direct` and `convolve_fft` on your machine for (N, M) = (10⁵, 10) and (10⁵, 2000). Report the timings as machine-specific and compare their ratio with the operation counts.

Answer form: a table and a crossover M, one explanation, two numbers with a location, two timings with a comment.

## 5. (Judgment: is convolution the right model?)

Answer each in a short paragraph with a number where one is asked for. (a) The elevator servo has a rate limit of 60 deg/s. Up to what doublet amplitude is the convolution prediction of the lesson (τ = 50 ms) valid, and what happens above it? (b) A colleague runs the servo model in a 50 Hz flight computer as y[n] = y[n−1] + (dt/τ)(x[n] − y[n−1]) "because h(0) dt = dt/τ". Compare its first impulse-response sample with the exact (1 − e^{−dt/τ}) and say which you would fly. (c) An impulse response of a system with τ = 0.5 s is recorded for only 2 s. What DC-gain error does the truncation cause? (d) The impulse response of a control-surface actuator is measured on a cold morning and used in a summer flight. What could change, and how would you detect it in flight data?

Answer form: four paragraphs, three with numbers.

## 6. (Hard: an exact convolution for piecewise-constant commands)

(a) Show that if x(t) is constant, equal to x_j, on each interval [t_j, t_{j+1}), then exactly y(t_k) = Σ_j x_j [s(t_k − t_j) − s(t_k − t_{j+1})], where s is the step response. Write this as a discrete convolution of x with the sequence d_m = s(m dt) − s((m − 1)dt). (b) Implement it with `convolve_direct` and `first_order_step_response`, and compute the doublet response of the lesson (τ = 50 ms, ±2 deg, 0.5 s halves starting at 0.1 s) at dt = 2, 1 and 0.5 ms. Compare the maximum error against `doublet_response_closed_form` with those of the rectangle and trapezoid rules. (c) Explain why the result is exact to rounding and why it is not a contradiction of 2.5. (d) For which inputs is the step-invariant formula *worse* than the trapezoid rule? Give an example.

Answer form: a derivation, a 3 × 3 error table, two short explanations.
