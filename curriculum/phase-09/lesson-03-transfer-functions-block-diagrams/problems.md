# P09.03  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: three connections by hand)

Let G₁ = 2/(s + 1) and G₂ = (s + a)/(s² + 2s + 5) with a = s/5.

(a) By hand, write the series G₁G₂, the parallel G₁ + G₂ and the unity negative-feedback loop around G₁G₂ as ratios of polynomials. (b) Give the poles, zeros and DC gain of each. (c) Check everything with `series`, `parallel`, `feedback`, `poles`, `zeros` and `dc_gain`.

Answer form: three transfer functions, a 3 × 3 table of poles, zeros and DC gains.

## 2. (Derivation: the rules and the realisation)

(a) Derive (2), T = G/(1 + GH), from the signal equations of a negative-feedback loop, and the corresponding formula for positive feedback. (b) Show that the closed-loop zeros are the zeros of G together with the poles of H, and that T → 1/H as the gain of G grows without bound. (c) For n = 2, verify by direct computation that the controllable canonical realisation of (b₂s² + b₁s + b₀)/(s² + a₁s + a₀) satisfies C(sI − A)⁻¹B + D = G. (d) Derive the two block-diagram moves: moving a pick-off point from after a block G to before it, and moving a summing junction from before a block G to after it. What block must be inserted in each case?

Answer form: four derivations.

## 3. (Unfamiliar scenario: a reaction-wheel attitude loop)

A satellite's single-axis attitude obeys J θ̈ = τ_w with J = 10 kg m². The reaction wheel delivers the commanded torque through a lag 1/(0.1s + 1). The controller commands τ_c = K_p(θ_c − θ) − K_d ω, with K_p = 2 N m/rad and K_d = 8 N m s/rad, where ω = θ̇ is measured by a gyro (take it as ideal).

(a) Draw the diagram with three dynamic blocks (wheel, body rate ω/τ_w = 1/(Js), angle θ/ω = 1/s) and reduce it by hand to θ/θ_c. (b) Give the closed-loop poles and DC gain. (c) Simulate the unreduced diagram with `simulate_diagram` (write the `connect` function from your drawing) and compare with `simulate_tf` of your reduction for a 1 rad step over 40 s. (d) Report the overshoot, the 2 % settling time and the peak wheel torque per radian of command. What sets the settling time?

Answer form: a transfer function, three poles and a DC gain, a maximum difference, three numbers and one sentence.

## 4. (Numerical, uses your code: tuning the pitch damper)

Use the lesson's pitch model with K_θ = 1 + s/10.

(a) For K_q = 0, 0.1, …, 1.0, tabulate the largest real part of the closed-loop poles, the overshoot and the 2 % settling time of the unit-step response (simulate the reduced transfer function for 8 s). (b) Which K_q gives the smallest settling time? Which gives the fastest slowest pole? Are they the same, and why might they differ? (c) For your chosen K_q, confirm the reduction by simulating the unreduced diagram and report the largest difference.

Answer form: an 11-row table, two gains with an explanation, one number.

## 5. (Judgment: what the reduced transfer function hides)

(a) A plant 1/(s − 1) (an unstable mode, such as a statically unstable airframe) is driven by a controller (s − 1)/(s + 2). Show that the series connection reduces with `minreal` to 1/(s + 2). Simulate the unreduced diagram for a unit step, first exactly as drawn and then with a constant disturbance of 10⁻⁶ added at the plant input. Report the plant output at 5, 10 and 20 s in both cases and explain. (b) Compute the roots of the expanded (s + 1)ⁿ for n = 4, 10 and 16 with `np.roots` and report the largest error. What does this mean for reducing a diagram with many identical blocks? (c) The lesson's peak elevator is 2.295 rad per rad of attitude command. If the elevator saturates at 25 deg, up to what attitude step is the reduced transfer function valid? (d) Name two things in a real pitch loop that no transfer function of this lesson represents.

Answer form: a reduced TF, a 2 × 3 table with an explanation, three errors with one sentence, an angle, two items.

## 6. (Hard: altitude hold around the pitch loop)

An outer loop holds altitude h. The attitude loop of the lesson (K_q = 0.3, K_θ = 3) is the inner block T(s) = θ/θ_c. At airspeed V = 30 m/s the flight-path angle follows pitch with the lag of the short-period zero, γ/θ = (1/T_θ2)/(s + 1/T_θ2), and ḣ = Vγ. The altitude controller is θ_c = K_h(h_c − h).

(a) Write h/θ and the closed-loop h/h_c in terms of T(s); give the order of h/h_c and its DC gain. Explain why one closed-loop pole sits exactly at −1/T_θ2. (b) With `series` and `feedback`, find the largest real part of the closed-loop poles for K_h = 0.001, 0.005, 0.01, 0.02, 0.05 rad/m, and by bisection the K_h at which the loop becomes unstable, with the frequency of the crossing. (c) For K_h = 0.01 rad/m, simulate the unreduced diagram (the pitch-loop blocks plus a γ-lag block and an altitude integrator) for a 10 m step and compare with the reduced transfer function. Report the peak altitude, the 2 % settling time and the initial pitch command. (d) Is 0.57 deg of pitch per metre of altitude error a sensible gain for a small UAV? What limits would you add before flying it?

Answer form: two transfer functions with an explanation, a five-row table and a critical gain with a frequency, a difference and three numbers, a paragraph.
