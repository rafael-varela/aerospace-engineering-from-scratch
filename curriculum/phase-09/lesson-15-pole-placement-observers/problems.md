# P09.15  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which come from your code, and which are simulations (noise records are SYNTHETIC).

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: one axis by hand)

A satellite axis has J = 2 + s/10 kg m² and only a star tracker.

(a) By coefficient matching, find K for regulator poles with ζ = 0.8 and ω_n = 0.2 rad/s, and the reference gain N. (b) Find L for observer poles at −0.8 and −1.0 s⁻¹. (c) Check all three with `ackermann`, `observer_gain` and `reference_gain`, and list the four eigenvalues of `observer_closed_loop`. (d) What is the peak wheel torque for a 5° step with a perfect initial estimate?

Answer form: two gains and N; L; a code check with four eigenvalues; one torque.

## 2. (Derivation: Ackermann, duality, separation)

(a) For the controllable canonical form (`P09.13`), show that u = −Kx replaces the first row of A by −(a + K), so each gain entry sets one coefficient of the closed-loop polynomial. (b) Derive Ackermann's formula (2) of the lesson, using Cayley–Hamilton for A − BK, or by transforming the canonical-form result back with T = W_c W̃_c⁻¹. (c) Show that eig(A − LC) = eig(Aᵀ − CᵀLᵀ), and hence that L = ackermann(Aᵀ, Cᵀ, poles)ᵀ. (d) Derive the block-triangular form (5) and state the three consequences listed in §2.4. (e) Derive N, (3) of the lesson, and show N = k₁ for the satellite.

Answer form: five derivations.

## 3. (Numerical, uses your code: full-state pitch control)

Use the 4-state pitch model (`pitch_model`: α, q, θ, δ, input δ_cmd, output θ). Let a = 3 + s/10.

(a) Place the regulator poles at −a ± aj, −6 and −20 and give K and N. Check K against the eigenvalues of A − BK. (b) Design an observer from θ alone with every pole twice the corresponding regulator pole; give L and the eight eigenvalues of the loop. (c) Simulate a 5° pitch step for 4 s (dt = 1 ms, x̂(0) = 0) and give the overshoot, the peak elevator command and the peak angle of attack. (d) Compare with the ±15° elevator limit of `P09.11`. What in the structure u = Nr − Kx̂ causes the peak, and name two fixes.

Answer form: two vectors and a number with a check; a vector and eight eigenvalues; three numbers; two sentences and two fixes.

## 4. (Unfamiliar scenario: a thrust-vector-controlled rocket)

A small rocket's pitch dynamics in the dense atmosphere, aerodynamically unstable, are θ̈ = aθ + bδ, with a = 4 s⁻² and b = 10 s⁻² per rad of gimbal angle δ (representative numbers).

(a) Give the open-loop poles. (b) Find K for poles with ζ = 0.7, ω_n = 4 rad/s, and for a double pole at −2, and show by hand that k₁ = (ω_n² + a)/b. Why must k₁ exceed a/b whatever poles you choose? (c) Design an observer from θ with poles −10 and −12. (d) Now only a rate gyro works. Unlike the satellite, the rocket is observable from the rate alone. Show this with the rank test, give L for the same observer poles, and explain physically where the angle information comes from.

Answer form: two poles; two gains with a derivation and one sentence; one L; a rank, an L and one sentence.

## 5. (Judgment: what separation does not promise)

(a) The satellite observer and regulator of the lesson are built on J = 3 kg m², but the true inertia is 3.6 or 2.4 kg m². Compute the eigenvalues of the true loop (plant with the true J, observer and gains with the model's J) and the steady value of θ/r. Is the loop still stable, and why is the steady error still zero? (b) Explain why the separation principle is a statement about the model and not about the hardware. (c) A colleague proposes observer poles 20 times faster than the regulator "so that the estimate is always perfect". Use the lesson's SYNTHETIC noise results and (6) to argue against it, with numbers. (d) Which measurement in flight would tell you that the observer's model is wrong?

Answer form: two sets of four eigenvalues with one number and one sentence; three short paragraphs.

## 6. (Hard: measuring the decay rate of an oscillating error)

Use the lesson's satellite regulator and an observer with complex poles −0.3 ± 0.3j.

(a) Give L. (b) Simulate (dt = 0.05 s, 60 s, x(0) = (0, 0.001), x̂(0) = 0, r = 0.1 rad) and fit the decay of ‖e‖ with `fit_decay_rate` over 20–60 s and over 5–60 s. Why is the result not exactly 0.3, and why does the window matter? (c) Fit instead the norm of the modal coordinates z = V⁻¹e (V the eigenvectors of A − LC), and the local maxima of ‖e‖. Which estimate is unbiased, and why? (d) Show that for a complex pair σ ± jω the Euclidean norm of e(t) is e^{σt} times a bounded periodic function, and say what that function does to a least-squares fit over a window that is not a whole number of half-periods.

Answer form: one L; two rates and two sentences; two rates and one sentence; a derivation.
