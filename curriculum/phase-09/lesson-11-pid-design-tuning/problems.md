# P09.11  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which are from your code, and which are simulations of a representative model (none of them is a flight measurement).

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: the filtered PID in the frequency domain)

A PID has Kp = 2, Ki = 4 s⁻¹, Kd = 0.3 + s/100 s and Tf = 0.02 s.

(a) Write C(s) as one ratio of polynomials. (b) Evaluate |C(jω)| and arg C(jω) at ω = 1, 10, 100 and 10⁴ rad/s with `pid_tf`. (c) Give the high-frequency limit of |C| and compare |C(100j)| with the unfiltered PID Kp + Ki/s + Kd s at the same frequency. (d) In one sentence each, say which term dominates at 1 rad/s and at 10⁴ rad/s, and why a pitch-rate gyro's noise makes Tf necessary.

Answer form: one rational function, an eight-entry table, two numbers, two sentences.

## 2. (Derivation: the two-degree-of-freedom loop and pole placement)

(a) For u = Kp(b r − y) + (Ki/s)(r − y) − D(s) y and y = P u, derive T = y/r = P C_r/(1 + PC) and state C_r. Show that the characteristic polynomial does not depend on b or on whether D acts on y or on e. (b) Show that T(0) = 1 whenever Ki > 0 and P(0) ≠ 0 (P may have a pole at the origin). (c) Show that the characteristic polynomial den_P(Tf s² + s) + num_P num_C is affine in (Kp, Ki, Kd), and write the three real equations that make a complex pair λ, λ* and a real pole μ roots of it. (d) For P = β/(s(s + a)) and Tf = 0 solve them in closed form and show Kd = (c₁ − a)/β, Kp = c₂/β, Ki = c₃/β for the target s³ + c₁s² + c₂s + c₃. (e) Derive the steady pitch error of a PD loop (Ki = 0) under a constant disturbance moment equivalent to an elevator angle d, using the steady state of (1) in the lesson.

Answer form: five derivations.

## 3. (Numerical, uses your code: tune your own servo)

Use the lesson's pitch plant with the servo bandwidth changed to ω_s = 10 + s rad/s.

(a) Compute Ku, ω_u and Pu with `ultimate_gain`, and the classic Ziegler–Nichols PID gains (Tf = 0.02 s). (b) Compute the gains that place −4 ± 4j and −1.5 with `pid_place_three_poles`, and list all six closed-loop poles of both designs. (c) For a 10° step with dt = 1 ms, tabulate overshoot, 10–90 % rise time, 2 % settling time and peak elevator for both designs with b = 1 and b = 0. (d) Which design would you fly, and which number decides it?

Answer form: three numbers and three gains; three gains and two pole lists; an eight-row table; two sentences.

## 4. (Unfamiliar scenario: roll attitude of a small UAV)

Roll dynamics with the aileron servo neglected: ṗ = L_p p + L_δ δ_a, φ̇ = p, with L_p = −(5 + s/5) s⁻¹ (roll damping) and L_δ = 30 s⁻² (aileron power). The aileron is limited to ±20°.

(a) Write φ/δ_a and show it has the design-model form β/(s(s + a)). (b) Place −3 ± 3j and −3 with Tf = 0.01 s, derivative on the measurement; give the gains and all closed-loop poles. Which pole did you not place, and where is it? (c) Derive the transfer function from the bank command to the aileron, U/R = C_r/(1 + PC). (d) For a 30° bank command, give the overshoot, rise time, settling time and peak aileron with b = 1 and with b = 0, using `tf_step`. Does either design saturate the aileron?

Answer form: one TF; three gains and four poles with one sentence; one TF; an eight-entry table and a yes/no for each design.

## 5. (Numerical: windup and the tracking time constant)

With the lesson's plant, the placed gains (−4 ± 4j, −1.5, Tf = 0.02 s), b = 1 and the ±15° elevator limit, command a pitch step of (10 + s)°.

(a) Predict, before running, whether back-calculation with Tt → 0 behaves more like clamping or like no anti-windup, and why. (b) Tabulate overshoot and 2 % settling time for: no anti-windup, clamping, and back-calculation with Tt = 0.1, 0.3, 1.0 s and √(Ti Td) (give the value). (c) Explain the trend with Tt in terms of how fast the integrator is drained while saturated.

Answer form: one sentence; a six-row table with one number; three sentences.

## 6. (Judgment: margins, flight test and model limits)

(a) Compute the gain margin (with `ultimate_gain` applied to the loop L = PC) and the phase margin and crossover of L (a frequency sweep) for the ZN and placed designs on the lesson's plant. (b) The ZN experiment raises the proportional gain until the loop oscillates. Explain why that is unacceptable as a flight test of a pitch loop, and say what a flight-test engineer would do instead to find the margins. (c) The lesson's 20° step drives α to 14° with clamping anti-windup. List, in order, three effects missing from the model that would change the answer, with the direction each pushes the result. (d) A colleague proposes Tf = 0.002 s "for a more ideal derivative". Give one benefit and two costs, with numbers from `pid_tf` (the high-frequency gain) and from the sample time.

Answer form: a four-entry table; three short paragraphs and a list.

## 7. (Hard: what the back-calculation law converges to)

(a) While the output is saturated at u = δ_max, show that the back-calculation integrator obeys, in continuous time, İ = Ki e − (I − I*)/Tt, and find I* in terms of δ_max, the proportional term and the derivative term. (b) Show that if e and the other terms are frozen, I approaches I* + Ki Tt e with time constant Tt, and interpret this: what value does the unsaturated output v approach? (c) Show that as Tt → 0 the law keeps v at the limit, and as Tt → ∞ it becomes no anti-windup. (d) Verify (b) numerically with frozen signals: with the placed gains, b = 1, r = 20° and y = 0 held constant (the aircraft has not responded yet), call `pid_step` 3000 times at dt = 1 ms with Tt = 0.3 s and Tt = 0.1 s. Compare the final I and v − δ_max with your prediction, and fit the time constant of I's approach. (e) In the closed-loop 20° step with Tt = 0.3 s, is v − δ_max close to Ki Tt e during saturation? Explain why or why not.

Answer form: one ODE and one expression; one limit and one sentence; two limits; a two-row table with fitted time constants; two sentences.
