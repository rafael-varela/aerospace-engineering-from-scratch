# P09.11  Worked solutions

Worked for **s = 15**. Every number below is produced by `solution.py` (the functions named in each problem) or by the probe scripts that call it. All responses are simulations of a representative model, not flight measurements. Units: SI; angles in degrees where stated, radians inside the code.

## 1. The filtered PID

**Assumptions.** Kp = 2, Ki = 4 s⁻¹, Kd = 0.3 + 15/100 = 0.45 s, Tf = 0.02 s.

**Derivation.** C(s) = [(Kp Tf + Kd)s² + (Kp + Ki Tf)s + Ki]/[s(Tf s + 1)] = **(0.49 s² + 2.08 s + 4)/(0.02 s² + s)**.

**Result** (`pid_tf`):

| ω (rad/s) | 1 | 10 | 100 | 10⁴ |
|---|---|---|---|---|
| \|C\| | 4.0792 | 4.8612 | 21.9153 | 24.4997 |
| arg C (deg) | −60.50 | +53.88 | +24.13 | +0.26 |

High-frequency limit **Kp + Kd/Tf = 24.5**. The unfiltered PID at 100 rad/s has |Kp + Ki/(100j) + 45j| = **45.00**, about twice the filtered value, and grows without bound above that. At 1 rad/s the integral term (4/s, magnitude 4) dominates and pulls the phase toward −90°. At 10⁴ rad/s the filter has capped the derivative, and C is a real gain of 24.5. A rate gyro's noise at tens to hundreds of hertz would reach the servo multiplied by Kd ω without the filter: at 100 Hz (628 rad/s) that is 283 instead of 24.5.

**Independent check.** At ω = 1/Tf = 50 rad/s the derivative term alone is Kd ω/√2 ≈ 15.9 at +45°, which is consistent with the rise in |C| between 10 and 100 rad/s.

**Common mistakes.** Writing the high-frequency gain as Kd/Tf (forgetting Kp); evaluating Ki/(jω) with the wrong sign of phase.

## 2. Derivations

**(a)** With u = C_r r − C y, C_r = b Kp + Ki/s (plus D(s) when the derivative acts on the error) and C = Kp + Ki/s + D(s): y = P(C_r r − C y) ⇒ y(1 + PC) = P C_r r ⇒ **T = P C_r/(1 + PC)**. Clearing fractions, the poles are the roots of den_P den_C + num_P num_C. C_r does not appear, so b and the derivative input change only the zeros of T. The test `test_closed_loop_tf_is_P_Cr_over_1_plus_PC` checks four combinations pointwise and that they share a denominator.

**(b)** T = PC_r/(1 + PC) = P C (C_r/C)/(1 + PC). As s → 0, C ≈ Ki/s → ∞ and C_r ≈ Ki/s too, so C_r/C → 1. If P(0) ≠ 0 (finite or infinite), |PC| → ∞ and PC/(1 + PC) → 1. Hence **T(0) = 1**.

**(c)** χ(s) = den_P(s)(Tf s² + s) + num_P(s)[(Kp Tf + Kd)s² + (Kp + Ki Tf)s + Ki] = χ₀(s) + Kp[num_P(s)(Tf s² + s)] + Ki[num_P(s)(Tf s + 1)] + Kd[num_P(s) s²]. That is affine in (Kp, Ki, Kd). The three equations are **Re χ(λ) = 0, Im χ(λ) = 0, χ(μ) = 0**: a 3 × 3 real linear system.

**(d)** With Tf = 0 and P = β/(s(s + a)), χ = s²(s + a) + β(Kd s² + Kp s + Ki) = s³ + (a + βKd)s² + βKp s + βKi. Matching s³ + c₁s² + c₂s + c₃ gives **Kd = (c₁ − a)/β, Kp = c₂/β, Ki = c₃/β**. `test_three_pole_placement_on_a_design_model_places_every_pole` confirms the general solve reproduces these to 10⁻⁴ with Tf = 10⁻⁶.

**(e)** In steady state θ̇ = q = 0. Then α̇ = Z_α α = 0 ⇒ α = 0, and q̇ = M_δ(δ + d) = 0 ⇒ δ = −d. A PD gives δ = Kp(b r − θ), since the derivative of a constant is zero, so with b = 1, **e = r − θ = −d/Kp**. For the lesson (d = −2°, Kp = 2.417) this is 0.827°, which is what `main()` prints.

**Common mistakes.** Putting b into the characteristic polynomial; forgetting that a complex pole gives two real equations, so a pair and a real pole are already three.

## 3. Your servo (s = 15: ω_s = 25 rad/s)

**Result.** `ultimate_gain` gives **Ku = 6.0299, ω_u = 8.4184 rad/s, Pu = 0.7464 s**. The classic ZN PID is **Kp = 3.6179, Ki = 9.6948 s⁻¹, Kd = 0.3375 s**, with closed-loop poles −53.52, −18.25, −1.385 ± 6.711j, −2.808 and −1.355. Placing −4 ± 4j and −1.5 gives **Kp = 2.6166, Ki = 6.1528 s⁻¹, Kd = 0.4832 s**, with poles −54.71, −10.46, −4 ± 4j, −4.031 and −1.5.

| design | b | overshoot (%) | rise (s) | settle (s) | peak elevator (deg) |
|---|---|---|---|---|---|
| ZN | 1 | 66.8 | 0.173 | 2.89 | 38.8 |
| ZN | 0 | 9.4 | 0.346 | 2.62 | 13.3 |
| placed | 1 | 37.6 | 0.238 | 1.63 | 27.5 |
| placed | 0 | 8.3 | 0.516 | 2.32 | 8.3 |

**(d)** Fly the placed gains with b = 0. They have the least overshoot and stay well inside the ±15° limit (8.3° peak elevator). ZN with b = 0 has similar overshoot but its closed-loop pair −1.385 ± 6.711j has ζ = 0.20, so any disturbance rings. The deciding number is that damping ratio, together with the peak elevator.

**Independent check.** The faster servo raises Ku from 5.004 (ω_s = 20) to 6.030, as expected from less phase lag at the crossing.

**Common mistakes.** Judging ZN by its b = 0 step alone: the poles, and therefore the disturbance response, are the same as with b = 1.

## 4. Roll attitude (s = 15: L_p = −8 s⁻¹)

**(a)** p(s − L_p) = L_δ δ_a and φ = p/s ⇒ **φ/δ_a = 30/(s(s + 8))**, the design-model form with β = 30 and a = 8.

**(b)** `pid_place_three_poles([30], [1, 8, 0], [-3+3j, -3-3j, -3], 0.01)` → **Kp = 1.1882, Ki = 1.782 s⁻¹, Kd = 0.0305 s**. Closed-loop poles: −3 ± 3j, −3 and **−99** (the unplaced one, the derivative filter pole, moved slightly from −1/Tf = −100). The ideal-PID closed form of 2(d) gives 1.2, 1.8 and 0.0333 for the target s³ + 9s² + 36s + 54; the filter accounts for the small differences.

**(c)** u = C_r r − C φ and φ = P u ⇒ **U/R = C_r/(1 + PC)** = num_Cr den_P/χ.

**(d)** For a 30° bank command:

| b | overshoot (%) | rise (s) | settle (s) | peak aileron (deg) |
|---|---|---|---|---|
| 1 | 32.8 | 0.278 | 1.62 | 36.12 |
| 0 | 0.0 | 0.834 | 1.44 | 10.82 |

With b = 1 the initial proportional kick Kp × 30° = 35.6° plus the integral **saturates the ±20° aileron**, so the linear prediction does not hold. With b = 0 the peak is 10.8°: **no saturation**, and no overshoot because this plant has no zero.

**Independent check.** At t = 0⁺ with b = 1 the aileron is b Kp r = 1.1882 × 30° = 35.6°, close to the 36.1° peak; with b = 0 the command starts at zero.

**Common mistakes.** Using θ-loop intuition (the slow zero) for roll: φ/δ_a has no finite zero, so b = 0 removes all overshoot.

## 5. Windup and Tt (s = 15: 25° step)

**(a)** Like clamping: a small Tt drains the integrator as fast as the error fills it, so I barely grows while saturated.

**(b)** Placed gains, b = 1, ±15° limit; Ti = Kp/Ki = 0.4133 s, Td = Kd/Kp = 0.1914 s, √(Ti Td) = **0.2812 s**:

| scheme | Tt (s) | overshoot (%) | settle (s) |
|---|---|---|---|
| none | — | 69.6 | 2.79 |
| clamp | — | 10.8 | 2.28 |
| back-calc | 0.1 | 10.1 | 2.29 |
| back-calc | 0.2812 | 19.6 | 2.10 |
| back-calc | 0.3 | 20.7 | 2.08 |
| back-calc | 1.0 | 43.9 | 2.02 |

**(c)** The overshoot rises monotonically with Tt: a longer tracking time constant leaves more integrator charge at the moment the error reverses. At Tt = 0.1 s back-calculation is as good as clamping. At Tt = 1 s it recovers less than half of the windup penalty.

**Common mistakes.** Choosing Tt larger than Ti "to be safe": the anti-windup then does almost nothing.

## 6. Judgment

**(a)** Lesson plant (ω_s = 20 rad/s):

| design | GM | GM (dB) | at ω (rad/s) | PM (deg) | crossover (rad/s) |
|---|---|---|---|---|---|
| placed | 9.439 | 19.50 | 28.77 | 43.78 | 6.417 |
| ZN | 9.079 | 19.16 | 23.99 | 22.09 | 6.325 |

The gain margins are nearly equal. The ZN design's weakness is its phase margin (22°, a delay margin of 61 ms, against 119 ms for the placed design).

**(b)** A sustained oscillation at Ku is a closed loop on the stability boundary. Any increase in gain, from a gust or a speed change, makes it diverge, with loss of control and large structural loads. Flight test measures margins without crossing that boundary. A small frequency sweep or multisine is injected at the servo command with the loop closed. The broken-loop frequency response is computed from the measured signals, and gain and phase margins are read from it. This is done at a safe altitude, with a safety pilot and a disengage switch.

**(c)** In order of how soon each binds in a 20° pull-up: (1) the **servo rate limit**. The linear loop asks for 51° of elevator in a fraction of a second; a rate limit slows the response and increases the time saturated, so windup gets worse. (2) **Airspeed loss and the phugoid** over 5–20 s. Pitch-up bleeds speed, which lowers M_δ and the elevator effectiveness, so more elevator is needed and saturation lasts longer. (3) **Nonlinear lift near α = 14°**. The lift-curve slope falls, Z_α shrinks and stall approaches, so the response is no longer linear, and in the worst case the wing departs.

**(d)** Benefit: a derivative closer to Kd s up to higher frequency, a little more phase lead at crossover. Costs: the high-frequency gain grows from Kp + Kd/Tf = **25.5 to 233.7**, nine times more sensor noise to the servo. And Tf = 2 dt is barely resolved: the backward-Euler filter pole Tf/(Tf + dt) = **0.667** differs from the exact e^{−dt/Tf} = **0.607**, while at Tf = 0.02 s the two are 0.9524 and 0.9512.

## 7. Back-calculation

**(a)** With u = δ_max and v = P + I + D: İ = Ki e + (δ_max − P − D − I)/Tt = **Ki e − (I − I*)/Tt, I* = δ_max − P − D**.

**(b)** With e, P, D frozen, the equation is a first-order lag with fixed point **I∞ = I* + Ki Tt e**, reached with time constant Tt. Then v = P + D + I → **δ_max + Ki Tt e**: the unsaturated command parks a margin Ki Tt e above the limit.

**(c)** Tt → 0: v → δ_max, the integrator is held at the value that just saturates (clamping-like). Tt → ∞: the correction term vanishes, which is no anti-windup.

**(d)** Placed gains (Kp = 2.4174, Ki = 5.8484), r = 20° = 0.34907 rad, y = 0, δ_max = 0.26180 rad. P = Kp r = 0.8438 and D = 0, so **I* = −0.5820**:

| Tt (s) | predicted I∞ | final I | v − δ_max | Ki Tt e | fitted time constant (s) |
|---|---|---|---|---|---|
| 0.3 | 0.03041 | 0.03041 | 0.61244 | 0.61244 | 0.2991 |
| 0.1 | −0.37788 | −0.37788 | 0.20415 | 0.20415 | 0.0995 |

The fitted time constants are slightly below Tt because the discrete update multiplies the deviation by (1 − dt/Tt) per step, i.e. τ = −dt/ln(1 − dt/Tt).

**(e)** No, only roughly. In the closed-loop 20° step with Tt = 0.3 s, v − δ_max is 0.5717 against Ki Tt e = 0.6107 at 0.05 s, 0.4649 against 0.5822 at 0.15 s, and 0.0229 against 0.2504 at 0.497 s, just before the output leaves saturation. The error and the P and D terms change on the same 0.1–0.5 s time scale as Tt, so the integrator never reaches its quasi-steady value.

**Common mistakes.** Sign errors in (u − v): the correction must reduce I when v exceeds the upper limit.

## Check-yourself answers (lesson.md)

1. C(s) = Kp + Ki/s + Kd s/(Tf s + 1) = [(Kp Tf + Kd)s² + (Kp + Ki Tf)s + Ki]/[s(Tf s + 1)]; high-frequency gain Kp + Kd/Tf.
2. The poles are the roots of 1 + PC, and b appears only in C_r, which sets the zeros of T = P C_r/(1 + PC). With b = 0 the reference enters through the integrator only, which removes the PID zeros from T: in the lesson's design the overshoot falls from 41.1 % to 10.6 % with identical poles.
3. The smallest proportional gain at which the closed loop has roots on the imaginary axis, so the loop oscillates steadily. It is located where Im[den(jω) num(−jω)] = 0 with K = −den(jω)/num(jω) real and positive (Ku = 5.004 at 7.633 rad/s for the lesson's plant).
4. The characteristic polynomial is affine in three gains, so it can satisfy only three real root conditions. The other three roots (here −5.45 ± 1.79j and −53.3) must be computed and checked: they could be slow or poorly damped.
5. While the actuator is saturated, the integrator keeps accumulating an error the actuator cannot reduce; the stored integral then drives a large overshoot once the error reverses. Clamping stops integrating while saturated and pushing further; back-calculation keeps integrating but drains the integrator toward the limit with time constant Tt.
6. e = −d/Kp = 1°/2.5 = 0.4° (nose low for a nose-down moment); integral action removes it.
7. About 0.011°, half of 0.023°: the sampled loop is a first-order approximation, error ∝ dt (the test measures the ratio 2.0).
8. It deliberately drives the closed loop to the stability boundary; in flight any small change in gain (speed, gust) makes it diverge, with loss of control and structural loads. Margins are measured instead with small injected signals in a stable loop.
