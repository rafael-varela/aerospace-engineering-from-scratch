# P01.02  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer, say which results are exact and which are numerical, and name the assumption behind every optimum you report.

Several problems use your personal seed **s**: the day of the month you were born (1–31). The lesson's glider (350 kg, 10.5 m², AR 20, e = 0.9, C_D0 = 0.015) is used wherever "the glider" is named.

## 1. (Warm-up: your own glider)

Your glider has mass m = 300 + 5s kg, wing area S = 10 + s/10 m², aspect ratio AR = 15 + s/4, Oswald factor e = 0.9 and C_D0 = 0.012 + s/5000. At sea level (ρ = 1.225 kg/m³, g = 9.80665 m/s²), find by hand: (a) K; (b) C_L* and (L/D)max; (c) the best-glide airspeed (small-angle) and the glide angle in degrees; (d) the sink rate; (e) the still-air range from 1000 m. Check (b) with `best_glide_cl_newton` and (c)–(e) with `airspeed`, `sink_rate`, `glide_range`. Answer form: six numbers with units, four significant figures, and the code's values beside them.

## 2. (Derivation: minimum sink and the sensitivity's error term)

(a) Starting from w = V sin γ with the small-angle approximations and (2), show that minimising the sink rate at fixed weight and density means minimising C_D/C_L^(3/2), and derive C_L,ms = √(3C_D0/K). Apply the second-derivative test (or the end-behaviour argument) to show it is a minimum.
(b) Show that the ratio of the sink rate at best glide to the sink rate at minimum sink is 3^(3/4)/2, independent of the aircraft.
(c) Derive (8), d(L/D)max/dC_D0 = −(L/D)max/(2C_D0), and the second derivative d²(L/D)max/dC_D0² = (3/4)(L/D)max/C_D0². Hence show the linear approximation's error for a relative change δ = ΔC_D0/C_D0 is (3/8)(L/D)max δ² to leading order, and check it against the code's −0.106 for δ = 0.1.

Answer form: derivation with every rule named; one line per check.

## 3. (Unfamiliar scenario: a jet's best-range lift coefficient)

For a jet aircraft in cruise, Phase 10 will show that range is proportional to C_L^(1/2)/C_D rather than to C_L/C_D, because thrust-specific fuel consumption is fixed and the aircraft should fly fast.

(a) Maximise C_L^(1/2)/(C_D0 + K C_L²) with respect to C_L and show the optimum is C_L,jet = √(C_D0/(3K)) = C_L*/√3.
(b) Show the glide ratio at that point is (√3/2)(L/D)max, independent of the polar, and compute it for the glider's polar.
(c) Compute the corresponding airspeed for the glider's weight and area and compare with V*. By what factor is the jet-optimum speed faster, and why does the ratio not depend on the aircraft?
(d) Confirm (a) with `newton` applied to the derivative of C_L^(1/2)/C_D (write the derivative by hand) from a start of your choice, and report the iteration count.

Answer form: derivation, three numbers with units, one sentence, one integer.

## 4. (Numerical, uses your code: wind, and Newton's basin)

(a) Run `headwind_best_speed` and `ground_glide_ratio` for V_w = −20, −15, …, +20 m/s for the glider. Tabulate the best airspeed, its C_L, and the ground glide ratio. Plot the speed increase (V_best − V*) against V_w and fit a straight line through the headwind points: what fraction of the wind does the rule "add k × headwind" give, and over which range is it within 1 m/s?
(b) Show numerically that in a tailwind the best speed never drops below the minimum-sink speed, and explain in one sentence why it cannot.
(c) Start `best_glide_cl_newton` from C_L = 0.1, 0.2, …, 3.0 and record the outcome of each: converged to C_L*, converged elsewhere, or failed. Mark the boundaries of the basin and compare them with the inflection point √(3C_D0/K).
(d) For a start inside the basin, print |C_L,n − C_L*| at every iteration and verify that the number of correct digits roughly doubles per step.

Answer form: one table, one plot with a fitted slope and range, a basin map, a convergence table.

## 5. (Related rates and a decision: can we reach the field?)

The glider's lift source fails at h = 1500 m with a field d = 26 km directly upwind. The forecast headwind is 10 ± 3 m/s. The pilot flies the headwind-optimum airspeed V_b for the forecast wind.

(a) The required glide ratio is R(t) = d(t)/h(t) with ḋ = −(V_b − V_w) and ḣ = −w(V_b). Use the quotient rule to write dR/dt and show that dR/dt < 0 at the moment of failure exactly when the ground glide ratio exceeds d/h.
(b) Evaluate dR/dt at t = 0 for V_w = 7, 10 and 13 m/s (use `headwind_best_speed` and `sink_rate_of_speed`). For which winds does the glider reach the field?
(c) For each wind give the height needed at the start to reach the field with the forecast optimum speed, and the time to get there.
(d) Judgment: the wind is known to ±3 m/s. Would you turn towards the field or towards a nearer downwind field with a ground range of 12 km? State the decision rule you used in one sentence.

Answer form: one formula, three signed rates with units, three heights and times, one decision with a rule.

## 6. (Judgment: where the model stops)

For each statement, say which assumption of the lesson fails, estimate the size of the effect with one number, and give a measurable sign in flight data that it is happening:

(a) The glider is flown at the computed minimum-sink C_L of 1.6 but the wing stalls at C_L = 1.4.
(b) The glide starts at cruise speed 40 m/s immediately after the failure.
(c) The headwind at 1000 m is 10 m/s but only 4 m/s at 100 m.
(d) The polar was measured at 350 kg but the glider is flown at 420 kg with water ballast.

Answer form: four short paragraphs, each with a number and an observable.

## 7. (Hard, uses your code: sinking air)

The glider descends through air that is itself sinking at a uniform u m/s (a common situation on the lee side of a ridge). The sink relative to the ground is w(V) + u, and the ground glide ratio is V/(w(V) + u) in still horizontal wind.

(a) Derive the condition for the best speed by setting dG/dV = 0, and show it has the same form as (11) with V_w = 0 and w replaced by w + u. Interpret it geometrically on the polar plot.
(b) Modify (not overwrite) `headwind_best_speed` into a function `sinking_air_best_speed(u, ...)` and solve with Newton for u = 0.5, 1 and 2 m/s. Report the best airspeed, the C_L, and the ground glide ratio.
(c) By how much does 1 m/s of sinking air reduce the range from 1000 m? Compare with the effect of a 10 m/s headwind.
(d) Explain in two sentences why the optimum speed in sinking air is higher than best glide, using the derivative condition rather than intuition.

Answer form: one derivation, nine numbers with units, two comparisons, two sentences.
