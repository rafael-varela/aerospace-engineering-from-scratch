# P02.08  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a transform and an inverse by hand)

Let a = s/5 (your seed over five).

(a) Find F(s) = L{3 + 2e^{−at} + sin 4t} by the table and evaluate it at s = 2. Check with `laplace_numeric`.
(b) Invert Y(s) = (s + 2)/(s² + 4s + 13) by completing the square and the s-shift rule. Check with `residues` (give the poles and residues) and `inverse_laplace` at t = 0, 0.3, 1.0 s.
(c) Why does `laplace_numeric` refuse s = 0, and what happens to F(s) in (a) as s → 0?

Answer form: an expression and a number in (a), an expression, two poles, two residues and three values in (b), two sentences in (c).

## 2. (Derivation: the derivative rule, the step response, and the overshoot)

(a) Derive L{f′} = sF(s) − f(0) by integration by parts, stating the condition on f for the boundary term to vanish. Derive L{f″}.
(b) Transform y″ + 2ζω_n y′ + ω_n² y = ω_n² u(t) with y(0) = y₀, y′(0) = ẏ₀ and solve for Y(s).
(c) For a unit step from rest, split Y(s) = ω_n²/[s(s² + 2ζω_n s + ω_n²)] into partial fractions (take 0 < ζ < 1) and invert to obtain y(t) = 1 − e^{−ζω_n t}[cos ω_d t + (ζ/√(1 − ζ²)) sin ω_d t].
(d) Differentiate y(t), find the time of the first maximum, and derive the overshoot formula e^{−πζ/√(1 − ζ²)}. Evaluate both for the lesson's actuator.

Answer form: four derivations; two numbers in (d).

## 3. (Unfamiliar scenario: pressure lag in a pitot-static line)

The static-pressure line to an altimeter behaves as a first-order lag: τ ṗ_ind + p_ind = p_true, with τ = 0.5 + s/20 seconds. The aircraft climbs at a constant rate of (5 + s/4) m/s; near sea level the static pressure falls by 12 Pa per metre, so p_true(t) = R t with R = −12 × (5 + s/4) Pa/s, starting from p_ind(0) = 0.

(a) Write P_ind(s) and identify the poles. Why does `residues` refuse this function?
(b) Invert by hand (partial fractions with a repeated pole, or SymPy) and show that for t ≫ τ the indicated pressure lags the true one by a constant R τ.
(c) Convert the lag to an altitude error in metres. Is the altimeter reading high or low during the climb?
(d) Verify (b) by integrating the ODE with `solve_ivp` at t = 1, 5, 10 s.
(e) What does the final-value theorem, lim_{t→∞} f(t) = lim_{s→0} sF(s), say about the *error* e = p_true − p_ind? Apply it.

Answer form: a transform and two poles, a derivation and a constant, an altitude error with sign, three numbers, one limit.

## 4. (Numerical, uses your code: does a softer actuator meet the specification?)

A second actuator has ω_n = 25 rad/s and ζ = 0.3 + s/100.

(a) Give its poles, overshoot, peak time and 2 % settling time (from a 10 µs grid over 1 s), and compare the settling time with 4/(ζω_n).
(b) Does it meet "< 10 % overshoot and settled within 250 ms"? If not, which half fails?
(c) Compute the response to a 50 ms command pulse and report the peak fraction of the command that reaches the surface and when.
(d) Give |G(iω)| and the phase lag at 3, 10 and 25 rad/s. The aircraft's short-period mode is at about 3 rad/s: does the actuator distort commands at that frequency?

Answer form: two poles and four numbers in (a), a verdict in (b), two numbers in (c), six numbers and a sentence in (d).

## 5. (Judgment: what the transfer function cannot tell you)

(a) The bench step tests of the actuator at 1°, 3° and 10° give overshoots of 9.5 %, 9.4 % and 2 %, and the 10° response rises almost linearly at 45°/s. What is happening, and why does no transfer function describe it?
(b) The specification is checked with a 1° step. A pilot's largest command is 15°. Which of the four performance numbers (overshoot, peak time, settling time, bandwidth) become meaningless at 15°, and what would you measure instead?
(c) The poles of G are at −15 ± 20i. Sketch the region of the pole diagram in which a second-order actuator would meet "< 5 % overshoot, settled within 150 ms", and say whether the current actuator could be moved into it by changing ζ alone.
(d) Hydraulic stiffness drops by 30 % when the oil is cold. Which pole parameter moves, in which direction, and which specification item is at risk?

Answer form: four short paragraphs, each with a number.

## 6. (Hard: the shock-response curve of a half-sine pulse)

A half-sine command pulse u(t) = sin(πt/t_p) for 0 ≤ t ≤ t_p (zero after) is applied to the lesson's actuator (ω_n = 25 rad/s, ζ = 0.6).

(a) Show that U(s) = (ω_p/(s² + ω_p²))(1 + e^{−t_p s}) with ω_p = π/t_p, using the shift rule. (Hint: a half-sine pulse is a full sine plus the same sine delayed by t_p.)
(b) Hence write the response as y(t) = y_sin(t) + y_sin(t − t_p) u_{t_p}(t), where y_sin is the response to sin ω_p t for all t ≥ 0. Build y_sin with `response_with_initial_conditions` (U = ([ω_p], [1, 0, ω_p²])) and `inverse_laplace`, and evaluate y(t) for t_p = (10 + s) ms.
(c) Verify against `solve_ivp` with the pulse as input and against `convolve_causal` of the impulse response with the pulse.
(d) Sweep t_p = 25, 50, 100, 200, 400 ms and tabulate the peak of y. Explain the two limits: short pulses (peak ≈ pulse area × max g, the impulse limit) and long pulses (peak → 1 plus a small overshoot).

Answer form: a derivation, a formula, a peak value and time for your t_p, two agreement figures, a five-row table with two explanations.
