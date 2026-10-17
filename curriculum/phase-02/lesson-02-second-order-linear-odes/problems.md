# P02.02  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: an avionics shelf by hand)

An equipment shelf of mass m = 2 + s/10 kg hangs on isolators of total stiffness k = 8000 N/m with a damper c = 40 N s/m.

(a) Find ω_n (rad/s and Hz), ζ, ω_d, the two characteristic roots and the regime. (b) The shelf is displaced 5 mm and released from rest. Write x(t) and evaluate it at t = 0.1 s and 0.5 s. (c) Give the time constant of the envelope and the logarithmic decrement. (d) What damper coefficient would make the shelf critically damped? Confirm (a)–(d) with `free_response`, `characteristic_roots` and `damper_for_ratio`.

Answer form: six numbers with units, a formula, two evaluated values, one check.

## 2. (Derivation: the three regimes and the Wronskian)

(a) Starting from the guess x = e^{rt}, derive the characteristic equation of m ẍ + c ẋ + k x = 0 and show that its roots are ω_n(−ζ ± √(ζ² − 1)).
(b) For ζ < 1, combine the two complex exponential solutions with Euler's formula into the real form e^{−ζω_n t}(A cos ω_d t + B sin ω_d t) and determine A and B from x(0) = x₀, ẋ(0) = v₀.
(c) For ζ = 1, substitute x = u(t) e^{−ω_n t} into the equation and show that ü = 0, so the second solution is t e^{−ω_n t}.
(d) Prove Abel's identity: for the Wronskian W = x₁ẋ₂ − ẋ₁x₂ of any two solutions, dW/dt = −(c/m)W.
(e) Show that successive positive peaks of (5) are spaced by 2π/ω_d and derive δ = 2πζ/√(1 − ζ²) and its inversion ζ = δ/√(4π² + δ²).

Answer form: five short derivations.

## 3. (Unfamiliar scenario: the same equation in a sensor circuit)

The anti-aliasing stage of an accelerometer signal conditioner is a series RLC circuit, L q̈ + R q̇ + q/C = 0, with L = 10 mH, C = 4.7 μF and R = s Ω.

(a) Map (L, R, 1/C) onto (m, c, k) and find ω_n (rad/s and Hz) and ζ. (b) Classify the regime and find the value of R that makes the circuit critically damped. (c) The capacitor is charged and the circuit closed at t = 0 (q̇(0) = 0). Sketch q(t) qualitatively for your R and for the critical R, and say which one a filter designer would prefer and why. (d) Check (a) and (b) with `natural_frequency`, `damping_ratio`, `classify_damping` and `damper_for_ratio` by passing the electrical quantities as the mechanical arguments.

Answer form: a mapping, four numbers, a qualitative comparison, one check.

## 4. (Numerical, uses your code: your landing)

The lesson's strut (1200 kg, 150 kN/m) touches down at a sink rate of v₀ = 2 + s/10 m/s.

(a) For ζ = 0.2, 0.3, 0.4, 0.5, 0.7 and 1.0, tabulate the maximum stroke, the peak strut force, the bounce flag and the 2 % settling time (3 s record, 3001 points). (b) Find the ζ that minimises the peak force with `scipy.optimize.minimize_scalar` on [0.05, 1] and the bounce threshold by sweeping ζ from 0.1 to 0.6 in steps of 0.001. (c) Recommend a ζ and justify it with three numbers from your table. (d) Check the ζ = 0.4 response against `solve_ivp` (rtol 10⁻¹⁰) on the first-order form and report the maximum difference.

Answer form: a six-row table, two thresholds, a recommendation with numbers, one oracle difference.

## 5. (Judgment: what the strut model cannot say)

(a) The linear model predicts a force cv₀ on the airframe at the instant of touchdown. Explain why a real oleo strut (orifice damping F ∝ ẋ|ẋ|) does not, and what this does to the "optimum ζ" of part 4.
(b) The model lets the compression go negative and we call it a bounce. What actually happens physically at that instant, and what second spring in the real gear takes over?
(c) Which of the lesson's assumptions (linear damper, linear spring, rigid tyre, full weight at t = 0, one leg) would you remove first to predict the *maximum stroke* better, and which to predict the *peak load* better? One sentence each.
(d) Name one measurement from a drop test that would show, on its own, that the damper is not linear.

Answer form: four short paragraphs, each with a specific physical mechanism.

## 6. (Hard: the first peak and the least load)

(a) For the touchdown response x(t) = x_s + e^{−σt}[x₀ cos ω_d t + b sin ω_d t] with x₀ = −x_s, σ = ζω_n, b = (v₀ + σx₀)/ω_d, derive the time of the first velocity zero, tan(ω_d t*) = (−σx₀ + bω_d)/(σb + x₀ω_d), and evaluate it for ζ = 0.3 and 0.7 at v₀ = 3 m/s. Compare with `max_stroke`.
(b) Show that the strut force F = kX + cẊ can be written as m(g − Ẍ) and hence that its maximum over the response is where the deceleration is largest, not where the stroke is largest. Where is it for ζ = 0.3? For ζ = 1.0?
(c) The peak force has a minimum near ζ = 0.36 for v₀ = 3 m/s. Using (b), explain why the minimum exists, and predict (then verify with the code) whether the optimum ζ moves up or down when v₀ is halved to 1.5 m/s.

Answer form: a derivation with two checks, an argument with two locations, a prediction with a verification.
