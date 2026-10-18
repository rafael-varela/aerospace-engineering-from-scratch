# P09.21  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and the probe scripts used to write this file, given to the figures shown. Results computed from `write_synthetic_csvs` / `make_synthetic_dataset` are SYNTHETIC: they show the analysis, not a measurement.

## 1. A P gain for a gear motor ($s = 15$: K = 550 deg/s per unit duty, τ = 45 ms, 12.5 % overshoot)

**Assumptions.** The dead zone is compensated, so the loop sees the linear plant K/(s(τs + 1)); continuous control (no sampling); no saturation.

(a) ζ = −ln 0.125/√(π² + ln²0.125) = **0.5519**. KpK = 1/(4τζ²) = 18.24 1/s, so **Kp = 0.03316 duty/deg**; ω_n = √(KpK/τ) = **20.13 rad/s**.

(b) `predict_closed_loop([550], [0.045, 1, 0], 0.03316, 0, 20, 1.5)`: **rise 0.0866 s, peak 0.1871 s, settling 0.2899 s**, overshoot 12.50 %. The estimate 4/(ζω_n) = 0.360 s is 24 % high: it bounds the decay envelope, while the response itself enters the 2 % band earlier.

(c) u₀ = 0.03316 × 20 + 0.09 = **0.753**, above the 0.6 limit: the duty saturates for the first samples and the prediction does not hold. The largest step with a linear start is (0.6 − 0.09)/Kp = **15.4°**.

(d) ±(0.10 − 0.09)/Kp = **±0.30°**.

**Check.** Rounding Kp to 0.033 for the firmware gives a predicted overshoot of 12.41 %: rounding down can only lower it.

**Common mistakes.** Using the apparent gain of one step test (which the dead zone has reduced) instead of K above the dead zone; forgetting the compensation term in the first-sample duty; quoting 4/(ζω_n) as the prediction.

## 2. The four results the lab is built on

(a) For a step of size u the Laplace transform of the output (no delay) is Ku/(s²(τs + 1)) = Ku[1/s² − τ/s + τ²/(τs + 1)]. Inverting: y(t) = Ku[t − τ + τe^{−t/τ}] = **Ku[t − τ(1 − e^{−t/τ})]**. A delay L replaces t by t − L for t ≥ L (and y = 0 before). As t → ∞ the exponential vanishes, leaving the line **Ku(t − τ − L)**: slope Ku, time-axis intercept τ + L. Subtracting, y − Ku(t − τ − L) = **Kuτ e^{−(t − L)/τ}** for t ≥ L. ∎

(b) The closed loop is KpK/(τs² + s + KpK) = (KpK/τ)/(s² + s/τ + KpK/τ). Matching s² + 2ζω_n s + ω_n²: **ω_n = √(KpK/τ)**, 2ζω_n = 1/τ, so **ζ = 1/(2√(τKpK))**. From M_p = exp(−πζ/√(1 − ζ²)), ζ = −ln M_p/√(π² + ln²M_p), and **Kp = 1/(4τζ²K)**. ∎

(c) With Ki = Kp/τ, C = Kp(s + 1/τ)/s and G = (K/τ)/(s + 1/τ): the factor (s + 1/τ) cancels and CG = KKp/(τs). Then T = CG/(1 + CG) = **1/(τ_cl s + 1), τ_cl = τ/(KKp)**. Its step is 1 − e^{−t/τ_cl}: it reaches a fraction f at t = −τ_cl ln(1 − f), so t₁₀ = τ_cl ln(10/9), t₉₀ = τ_cl ln 10 and **t_r = τ_cl ln 9**; it stays within 2 % after t = **τ_cl ln 50**. ∎ (Exact cancellation needs the true τ; with an identified τ the zero and pole are close, not equal, and a slow residual mode remains, small when they are close.)

(d) With the lag, L(s) = ω_c/(s(θs + 1)), which is the loop of (b) with τ → θ and KKp → ω_c: **ζ = 1/(2√(θω_c))**, ω_n = √(ω_c/θ). At θω_c = 1/4, ζ = 1 and ω_n = 2ω_c. The critically damped step is 1 − e^{−x}(1 + x), x = ω_n t; bisection gives x₁₀ = 0.5318 and x₉₀ = 3.8897, so t_r = 3.3579/ω_n = **1.679/ω_c**, against ln 9/ω_c = 2.197/ω_c for the first-order loop: **23.6 % shorter**. ∎ The delay slows the start (the S-shape) but steepens the middle, and the 10–90 % metric sees only the middle. The settling time does not improve the same way, and at larger θω_c the loop overshoots.

**Common mistakes.** In (a), forgetting that the tail intercept is τ + L, not τ (so L comes out as zero or τ is overestimated); in (c), quoting the rise time as τ_cl (that is the 63 % time from zero).

## 3. A camera-gimbal pointing loop ($s = 15$: τ = 22.5 ms; paper only)

**Assumptions.** The inner rate loop behaves as 1/(τs + 1); the outer loop delay is one sample plus half a sample of hold; no rate limit.

(a) `design_p_position(1, 0.0225, 5)`: ζ = 0.6901, **Kp = 23.33 (deg/s)/deg**, ω_n = 32.20 rad/s. Continuous prediction: **rise 0.0651 s, settling 0.1862 s**, overshoot 5.0 %.

(b) At 100 Hz, θ = 1.5 × 10 ms = **15 ms**. Crossover |L(jω)| = 1 at **ω_c = 21.08 rad/s**; phase margin 90° − atan(τω_c) = **64.6°** without the delay; the delay costs ω_cθ = 18.1°, leaving **46.5°**. With the delay as a lag: overshoot **19.9 %**, rise 0.0581 s. The sampled simulation (one-sample delay, 100 Hz): overshoot **22.2 %**, rise 0.0507 s, settling 0.287 s.

(c) At 400 Hz, θ = 3.75 ms, the delay costs 4.5° and the lag model gives **8.0 %** overshoot, the sampled simulation **7.95 %** (rise 0.0600 s, settling 0.189 s).

(d) `design_p_sampled` at 100 Hz for 5 %: **Kp = 13.75**, rise 0.0884 s, 36 % slower than the continuous design's 0.0651 s. At 400 Hz the same procedure gives Kp = 19.86 and a rise of 0.0708 s. **Recommendation:** run the outer loop at 400 Hz (or faster) if the processor allows it; at 100 Hz the delay, not the gimbal, sets the achievable speed.

**Check.** The lag model and the sampled simulation agree within 0.05 points at 400 Hz and differ by 2.3 points at 100 Hz, where ω_cθ = 0.32 rad is no longer small and the lag approximation weakens.

**Common mistakes.** Designing on the continuous model and assuming "100 Hz is fast" because it is several times ω_c/2π (3.4 Hz): the margin lost is ω_cθ, not a ratio of frequencies.

## 4. The full analysis on seed 15 (SYNTHETIC)

**Assumptions.** Everything here is SYNTHETIC. The synthetic closed-loop runs use the firmware gains fixed in the generator (Kp = 0.029 and 0.022; Kp = 0.33, Ki = 5.9), designed on seed 2026; your seed changes the noise, not the gains.

(a) Motor per-run least-squares fits (apparent gain deg/s per unit duty, τ, L): duty 0.30: 360.2 ± 0.5, 44.98 ± 1.29 ms, 0.00 ± 0.98 ms; duty 0.45: 419.6 ± 0.3, 44.03 ± 0.73 ms, 0.55 ± 0.55 ms; duty 0.60: 449.7 ± 0.2, 44.82 ± 0.46 ms, 0.06 ± 0.34 ms. Dead-zone line: **K = 539.3 deg/s per unit duty, u_dz = 0.0997**, mean **τ = 44.61 ms**. Servo: **K = 0.9294, τ = 55.59 ms, L = 12.90 ms**. Designs: Kp = 0.02974 (motor, 10 %); Kp = 0.3285, Ki = 5.910 1/s (servo, 0.4 s).

(b) (rise s / overshoot % / settling s)

| run | continuous | sampled model | SYNTHETIC |
|---|---|---|---|
| motor 15°, Kp 0.029 | 0.0988 / 9.56 / 0.317 | 0.0855 / 21.39 / 0.360 | 0.0867 / 20.49 / 0.390 |
| motor 60°, Kp 0.029 | 0.0988 / 9.56 / 0.317 | 0.1922 / 8.81 / 0.461 | 0.1922 / 9.19 / 0.469 |
| motor 15°, Kp 0.022 | 0.1280 / 5.12 / 0.368 | 0.1168 / 9.13 / 0.408 | 0.1260 / 7.13 / 0.430 |
| servo 20°, PI | 0.4012 / 0 / 0.715 | 0.2790 / 0 / 0.499 | 0.2664 / 0.69 / 0.514 |
| servo 45°, PI | 0.4012 / 0 / 0.715 | 0.2791 / 0.04 / 0.501 | 0.2669 / 0.27 / 0.462 |

(c) Over seeds 15–24: synthetic 15° overshoot **20.78 ± 0.33 %**, sampled-model prediction **21.98 ± 0.41 %**. The mean difference, 1.2 points, is about 3.6 run-to-run standard deviations: **significant**, so a systematic effect is missing from the model. It is the bench's asymmetric dead zone (0.10 forward, 0.11 back): running the bench alone with a symmetric 0.10 moves its overshoot from 20.9 ± 0.5 % to 21.8 ± 0.4 %, onto the model.

(d) All of the above is SYNTHETIC.

**Common mistakes.** Comparing one synthetic run with the model and calling a 1-point difference "noise" without measuring the run-to-run scatter; reporting the standard errors of the fits as the whole uncertainty (they ignore model error).

## 5. A faster servo loop, designed on the refined model (shipped SYNTHETIC dataset)

**Assumptions.** Identified model K = 0.9282, τ = 55.55 ms, L = 12.6 ms (SYNTHETIC, seed 2026); Ts = 20 ms; one-sample delay; ±60° limit; 0.066° quantisation; cancellation kept.

(a) Bisection gives Kp = 0.4129; rounded **Kp = 0.41, Ki = 7.4 1/s**. Sampled model: **rise 0.201 s, overshoot 0.33 %, settling 0.340 s**.

(b) The continuous model predicts **0.320 s** for the same gains: it would have told you these gains are too slow.

(c) The continuous design for 0.20 s is Kp = 0.6575, Ki = 11.84 1/s. On the sampled model it gives **rise 0.100 s, overshoot 6.27 %, settling 0.349 s**: twice as fast as intended, with overshoot. Put the sampled-model design (0.41, 7.4) on the board, and write its sampled-model prediction (0.201 s, 0.33 %, 0.340 s) in the dated prediction, with a note of which model produced it.

**Check.** Rounding changes the sampled rise time by under 0.001 s (0.2002 s before rounding, 0.2011 s after).

**Common mistakes.** Treating the sampled model's prediction as a measurement; re-tuning on the synthetic closed-loop data and then "predicting" it (a prediction must come before the data it predicts).

## 6. Reading a closed-loop test

(a) The integrator meeting the servo's internal dead band (or gearbox backlash): the error integrates until the servo jumps, overshoots the band, and integrates back: a slow limit cycle. Decide by setting Ki = 0: the hunting stops and a small steady error appears; or add an error dead band of one count to the integrator.
(b) Loop delay (the one-sample write, the hold and any dead time) eroding the phase margin. Decide by running the same gain at twice the loop rate: a delay-driven excess overshoot must shrink roughly in proportion to θ.
(c) Saturation: the large step meets the duty limit and the motor runs at its limited speed. Decide from the log: the applied duty sits at 0.6 during the rise.
(d) Asymmetric friction or dead zone (or a gravity load on the pointer). Decide by open-loop steps of the same duty in both directions: the slopes differ.
(e) The servo's internal slew-rate limit: larger steps hit it, and a first-order fit absorbs it as a longer τ and L. Decide by plotting the maximum measured rate against step size: it flattens at the slew limit.
(f) The dead zone (and an audible PWM frequency): the duty Kp|e| + compensation is below the breakaway duty. Decide by a slow open-loop duty ramp, reading the duty at which the pointer first moves.
(g) Battery sag: K is proportional to the pack voltage. Decide by measuring the pack before and after, and repeating one open-loop step at the end of the session.

## Check-yourself answers

1. y = Ku(1 − e^{−(t − L)/τ}) for t ≥ L. ln(1 − y/y_ss) = L/τ − t/τ is linear in t: the slope is −1/τ and the intercept L/τ.
2. Slope Ku; time-axis intercept τ + L. The remainder y − Ku(t − τ − L) is Kuτ e^{−(t − L)/τ}; a line through its logarithm gives τ, then L = intercept − τ.
3. ζ = 1/(2√(τKpK)), ω_n = √(KpK/τ). For 10 %, ζ = 0.5912 and Kp = 1/(4 × 0.045 × 0.5912² × 540) = 0.02944 duty/deg (ω_n = 18.80 rad/s).
4. The zero cancels the plant pole, leaving the integrator KKp/(τs) in the loop, whose unity-feedback closed loop is first order with τ_cl = τ/(KKp); rise τ_cl ln 9 = 2.197 τ_cl, 2 % settling τ_cl ln 50 = 3.912 τ_cl.
5. About 1.5 samples, 15 ms (plus any plant dead time); at 13.4 rad/s that costs 13.4 × 0.015 = 0.201 rad = 11.5° of phase margin (58.9° → 47.3° in the lesson's loop).
6. ±(0.10 − 0.09)/0.029 = ±0.34°. With 0.11 the net duty at a tiny error is ±0.01 beyond the dead zone and flips sign with the noise: the output chatters around the target. Without compensation the band would be ±0.10/0.029 = ±3.4°.
7. Stop test (the latched button zeroes the output within one sample: any runaway can be ended), limit test (the output is clipped at 0.6 or ±60°: bounded speed and current), sign test (positive duty gives positive angle: a wrong sign is positive feedback and drives the output into a stall), travel-limit test (the firmware stops before the end stops: no stalled motor, no broken potentiometer). In that order, all before the first closed-loop run, the stop test again at the start of every session.
8. You may write that on the SYNTHETIC dataset the continuous model underpredicted the overshoot (9.7 % against 19.8 %), that a refined model with the loop delay reproduces it, and what that implies for a real loop. You may not call the 19.8 % a measurement, omit SYNTHETIC from the title, summary or figures, or revise the dated prediction after seeing the synthetic closed-loop data.
