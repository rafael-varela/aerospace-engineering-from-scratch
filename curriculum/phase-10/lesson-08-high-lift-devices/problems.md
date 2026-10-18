# P10.08  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. Angles in degrees in answers, radians in formulas. Say which results are exact (closed form), which numerical, and which rest on an assumption.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: an aileron)

An aileron has a chord ratio E = 0.20 + s/100 and deflects δ = 10° trailing edge down.

(a) By hand, find θ_h, τ, Δcl, Δα_L0 and Δcm_c/4 from equations (1) and (2) of the lesson.
(b) Compare τ with the small-flap approximation (4/π)√E. For which E is that approximation within 2 %?
(c) Check (a) with `flap_increments` and with `thin_airfoil(flap_slope(E, δ), 0)`.

Answer form: five numbers; a comparison and a range of E; a check.

## 2. (Derivation: flap theory)

(a) Derive A0, Aₙ for the kinked camber line (slope 0 ahead of the hinge, −δ behind) and from them equations (1) and (2) of the lesson.
(b) Show that |Δcm_c/4|/δ is largest at E = 0.25 and equals 3√3/8 there.
(c) Show that τ → (4/π)√E as E → 0 by expanding θ_h near π.
(d) Show that c_h tends to −(π/2)(α + δ) as E → 1, and explain why physically.

Answer form: four derivations.

## 3. (Unfamiliar scenario: the elevator's hinge moment)

The trainer's horizontal tail (2.0 m²) has an assumed span of 3.5 m, so a chord of 0.571 m, with a 35 % elevator. On approach at 30 m/s the tail is at α_t = −2° and the elevator at δ = −10° (trailing edge up). Treat the tail as a 2-D flat plate.

(a) Compute c_hα and c_hδ for E = 0.35 with `hinge_moment`, then c_h at the approach condition.
(b) Turn c_h into a hinge moment H in N·m. Which way does it push the elevator?
(c) Why do aircraft designers add aerodynamic balance (a horn or a set-back hinge) to control surfaces, and which sign of c_hδ would make the controls unstable?

Answer form: three numbers; a moment and a direction; a short argument.

## 4. (Numerical, uses your code: slat deflection)

Set the slat deflection to δ_s = 20 + s/3 degrees (round to the nearest 0.5°, and keep it between 20° and 30° so the elements do not overlap), keeping the slat's trailing edge where `slat_configuration` puts it.

(a) At α = 10°, report the main element's and the slat's leading-edge peaks and the total cl. Compare with δ_s = 20° and 30°.
(b) Find cl_max and α_max with the suction-peak criterion (cp_crit from the clean 2412 at cl_max = 1.6). Which element limits it?
(c) Explain the trend with δ_s in terms of the slat's own incidence.

Answer form: a table; two numbers and an element; an argument.

## 5. (Numerical, uses your code: your own trainer)

Your trainer has mass 700 + 5s kg, S = 12 m², a 30 % slotted flap at 30° over a fraction 0.5 + s/100 of the wing area, and the lesson's other assumptions.

(a) Compute the clean and flapped CL_max, V_s and the margin at a 30 m/s approach.
(b) What CL_max does a 1.3 margin need? Is it met?
(c) What approach speed would you recommend, and what is the single assumption you would test first?

Answer form: two sets of three numbers; a number and a yes/no; a speed and an assumption.

## 6. (Judgment: trusting the cl_max criterion)

(a) The criterion gives a 25 % plain flap Δcl_max = 0.975 at 30° and 1.279 at 40°. Why should you not believe these numbers, and what physical process does the model miss?
(b) For the slotted flap at 30°, the flap's own leading-edge peak is −16.5, more than twice the critical value, yet the lesson excludes the flap from the criterion. Defend or attack that choice.
(c) Propose a measurement or a dataset that would calibrate the criterion for slotted flaps, and say what you would compare.

Answer form: three short arguments.
