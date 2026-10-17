# P05.04  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. Draw a free-body diagram for every part you analyse, state units on every numerical answer, and say which direction each friction force or moment acts and why.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: an aircraft parked on a ramp)

An aircraft stands on a sloping apron with its brakes set. The tyre-to-surface friction coefficient is μ = 0.4 + s/100 (illustrative).

(a) Draw the free-body diagram of the aircraft on a slope of angle α with all wheels braked, and show that it stays put while tan α ≤ μ. Find the steepest slope.
(b) Now only the main wheels are braked and the nose wheel rolls freely. Assume the mains carry 90 % of the normal load whatever the slope. Find the steepest slope.
(c) Why does the answer not depend on the aircraft's mass? Name one effect the 90 % assumption ignores.

Answer form: two angles in degrees and two short explanations.

## 2. (Derivation: three results the code relies on)

(a) Derive the capstan law T₂/T₁ = e^{μθ} from the equilibrium of a cable element in contact with a fixed drum.
(b) Show that a pin of radius r turning in a bearing with friction coefficient μ exerts a friction moment r_f |R| with r_f = rμ/√(1 + μ²). Then derive the tension ratio k of a cable over a pulley of radius R on such a pin, for a wrap angle θ, as the root of a quadratic, and show that it reduces to (1 + ρ)/(1 − ρ) at θ = 180° (ρ = r_f/R).
(c) Use virtual work to show that the stick force of a frictionless one-degree-of-freedom linkage is F = −H (dδ/dφ)/L_g.

Answer form: three derivations, each with a sketch.

## 3. (Unfamiliar scenario: an over-centre toggle lock)

A landing-gear downlock uses a toggle: two equal links of length l = 0.1 m, pinned end to end, lie nearly in a straight line between a fixed pin and a slider. Each link makes a small angle θ with that line. An actuator pushes the centre pin with a force P = 20 N perpendicular to the line, which straightens the toggle and pushes the slider with a force Q along the line. Neglect friction.

(a) Draw the free-body diagrams of the centre pin and of the slider and show that Q/P = 1/(2 tan θ) and that each link carries P/(2 sin θ).
(b) Evaluate Q and the link force at θ = 10°, 5° and 1°.
(c) Check Q at θ = 5° by virtual work, differencing the slider and centre-pin positions numerically.
(d) Why is the lock made to go slightly *past* the straight position against a stop, and what would happen to the gearing of the elevator linkage in this lesson if one of its pushrods came into line with an arm?

Answer form: a derivation, a table of three rows, a check, a short explanation.

## 4. (Numerical, uses your code: change the gearing)

Change the bellcrank's second arm of the lesson's linkage to b₂ = 0.08 + s/1000 m (`Linkage(b2=...)`; the second pushrod's length is recomputed).

(a) Predict the small-angle gearing from the arm ratios, then check it with `gearing`.
(b) Compute δ at φ = ±20° and the frictionless stick force at φ = +20°. Compare with the lesson's linkage.
(c) Find, with `brentq`, the arm b₂ that gives exactly δ = −25° at φ = +20°. Report δ at −20°, the stick forces at ±20° and the neutral gearing. Why is the travel not symmetric?

Answer form: a gearing, three numbers, an arm length with four numbers and an explanation.

## 5. (Code: pushrods against cables)

Take the frictionless stick force of the lesson's linkage at φ = 20° as the ideal force. A cable run delivers the same ideal force path, but its friction multiplies the force by the run's ratio k when driving and divides by k when back-driven.

(a) Compute the ratio of three single plain-bushing pulleys (R = 38.1 mm, pin radius 5 mm, μ_b = 0.15) at wraps of 90°, 45° and 30°, and of a fairlead with a 3° wrap at μ = 0.2.
(b) A run has four pulleys (90°, 45°, 30°, 90°) and three 3° fairleads. Compute k, the driving and back-driven forces and the band, for plain-bushing pulleys and for ball-bearing pulleys (equivalent μ_b = 0.005). Compare with the pushrod linkage with μ = 0.2 pin friction.
(c) The fourth pulley (90°) seizes and the cable slides over it with μ = 0.3. What does the pilot feel?

Answer form: four ratios, a table of three systems, a short paragraph.

## 6. (Judgment: is friction a problem here?)

A designer proposes, to save weight, to replace the ball-bearing pulleys of a cable elevator system with nylon plain-bushing pulleys, and to replace the 6.35 mm stick-pivot bolt with a 12.7 mm one because the stick pivot sees the largest reaction.

(a) Estimate, with the lesson's numbers, how each change affects the friction band at full travel.
(b) What else besides the band should be considered: return to centre, breakout, the effect of airspeed, wear and environmental effects?
(c) What test would you ask for before accepting either change?

Answer form: two estimates and a short engineering argument.
