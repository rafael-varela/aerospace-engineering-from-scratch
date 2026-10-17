# P09.04  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. State the metric conventions (10–90 % rise, 2 % settling) with every answer.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: reading the poles)

A second-order system has ω_n = s/3 rad/s and ζ = 0.3.

(a) Give its poles, overshoot, peak time, envelope settling bound and 4/(ζω_n). (b) Find the exact 10–90 % rise time with `rise_time_exact`. (c) Measure all metrics from `second_order_step` on a fine grid with `step_metrics` and tabulate formula against measurement. Which number disagrees most, and why?

Answer form: two poles, a 6-row table, one sentence.

## 2. (Derivation: the formulas)

(a) By partial fractions of ω_n²/(s(s² + 2ζω_n s + ω_n²)), derive y(t) = 1 − e^{−σt} sin(ω_d t + φ)/√(1 − ζ²) with φ = arccos ζ. (b) Derive t_p = π/ω_d and the overshoot formula (1), and invert it for ζ. (c) Derive the envelope bound (2) and show it is an upper bound on the true settling time. (d) Derive (3) for the servo and prove that the real part of the closed-loop poles is −1/(2τ) for every K > 1/(4τ). (e) Derive (5) for a zero at s = z, and show that the initial slope of the step response is −ω_n²/z.

Answer form: five derivations.

## 3. (Unfamiliar scenario: a thrust-vector gimbal actuator)

A small launch vehicle's thrust-vector gimbal is positioned by an electromechanical actuator modelled as K/(s(τs + 1)) with τ = 0.02 + s/1000 s, under proportional feedback. The specification is overshoot ≤ 15 % and peak time ≤ 60 ms.

(a) Sketch the region of the s-plane where both requirements hold (a sector and a horizontal band). (b) Find the largest K that meets the overshoot requirement and the peak time it gives. Is the specification met? (c) For which actuator time constants τ could a proportional loop meet both requirements? (d) Confirm (b) by simulating the step response.

Answer form: a sketch with two boundaries, a gain and a time with a yes/no, a bound on τ, a measured overshoot and peak time.

## 4. (Numerical, uses your code: a filtered potentiometer)

The servo of the lesson (τ = 50 ms) measures its angle with a potentiometer behind an RC filter of time constant τ_f = s/1000 s, so the loop has the feedback block 1/(τ_f s + 1).

(a) Write the closed-loop transfer function and say where the filter pole appears. (b) With the lesson's K = 14.31 s⁻¹, compute the poles, the dominance ratio, and the overshoot, peak time and settling time from `residue_step_response`. (c) Explain why a dominance ratio near 10 does not mean "the filter does not matter". (d) Find, by bisection, the K that restores 10 % overshoot with the filter, and report its settling time.

Answer form: a transfer function, three poles with a ratio and three metrics, a paragraph, a gain and a time.

## 5. (Judgment: choosing the servo gain)

The shortest 2 % settling time of the proportional servo is 0.282 s at K = 8.2 s⁻¹ (§2.3). (a) Compute the settling time and overshoot at K = 8.2 s⁻¹ for τ = 40, 50 and 60 ms. (b) Compute the overshoot and settling at K = 14.31 s⁻¹ for τ = 40 and 60 ms. (c) The servo's angle is logged at 50 Hz. Measure the metrics of the K = 8.2 s⁻¹ response from samples every 20 ms and compare with a 0.5 ms grid. (d) Which gain would you recommend, and what would you say to a colleague who quotes "settling 0.28 s" in a specification sheet?

Answer form: a 3-row table, a 2-row table, a comparison, a paragraph.

## 6. (Hard: a right-half-plane zero)

Use ζ = 0.5, ω_n = 1 rad/s and a zero at z = s/5 (right half plane).

(a) From (5), show that the step response starts with slope −ω_n²/z, and that every right-half-plane zero produces undershoot. (b) Compute the undershoot, the time of the minimum, the overshoot and the 2 % settling time with `zero_step_closed_form` and `step_metrics`; check the initial slope numerically. (c) Using the feedback formula of `P09.03`, show that unity feedback around a plant with a right-half-plane zero leaves the zero in the closed loop for every gain. (d) Explain physically why an aircraft's altitude response to the elevator has such a zero, and what it implies for an altitude-hold autopilot that is too aggressive.

Answer form: a derivation, four numbers and a slope, a proof, a paragraph.
