# P05.01  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and draw a free-body diagram before writing any equilibrium equation.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: moments in three dimensions)

A ground force F = (−1000, 3000, 100s) N acts at r = (0, −1.2, 0.3) m, measured from a trunnion pin at the origin (x forward, y up, z to starboard).

(a) Compute the moment vector M = r × F by hand and its magnitude. Check with `moment_about_point`.
(b) Compute the moment about the fore-aft (x) axis through the origin, and about the axis through the origin along (1, 0, 1). Check with `moment_about_axis`.
(c) Verify that M · F = 0 and M · r = 0. Why must both hold?

Answer form: a vector and a magnitude, two scalars, two dot products and one sentence.

## 2. (Derivation: equivalent force systems)

(a) Prove that the moment of a force about an axis does not depend on which point of the axis is used.
(b) Prove the transport rule M_B = M_A + (A − B) × R, and deduce that a couple has the same moment about every point.
(c) Show that a planar force system with R ≠ 0 reduces to a single force, and that the central-axis point p = R × M_O/|R|² lies on that force's line of action. What does a planar system with R = 0 reduce to?

Answer form: three short proofs.

## 3. (Unfamiliar scenario: a rocket on its launch rail)

A model rocket of mass m = 10 + s/5 kg sits on a straight launch rail inclined at 80° to the horizontal, before ignition. It rides on two rail buttons, at 0.3 m and 1.1 m from its aft end along the rail, which can push on the rocket only perpendicular to the rail. A stop at the aft end of the rail (on the rail line) holds it axially. The rocket's centre of gravity is 0.9 m up the rail from the aft end and 0.05 m off the rail line, on the side of the rocket body (above the rail when the rail is tilted).

(a) Draw the free-body diagram. How many unknowns, how many equations?
(b) Find the two button forces and the stop force by hand, using a moment equation about a well-chosen point.
(c) Check with `solve_fbd_2d` and `equilibrium_residual_2d`.
(d) Which button carries more, and why does the 0.05 m offset matter?

Answer form: an FBD, three forces with units, a residual, two sentences.

## 4. (Numerical, uses your code: a raked main gear)

Your aircraft has mass m = 600 + 20s kg; use `landing_load_case` with its defaults. The leg is raked: the tyre contact point is 0.10 m **ahead** of the trunnion, C = (0.10, 0); A, B and E are as in `GearGeometry`.

(a) Compute the trunnion reactions and the brace force for D/V = 0, 0.25 and 0.5. Check each brace force with `brace_force_by_moments`.
(b) At D/V = 0 the brace is loaded although there is no drag. Explain why, and say whether it is in tension or compression.
(c) Find, by hand, the drag ratio at which the brace carries nothing, and confirm it with the code.
(d) Compare the brace force at D/V = 0.25 with the vertical leg of the lesson for the same mass. Which geometry would you choose, and what else would you need to know?

Answer form: a 3 × 3 table, a sentence, a number, a comparison.

## 5. (Derivation and numerics: wing lift distributions)

(a) Derive the resultant and centroid of the elliptic load w(y) = w₀√(1 − (y/s)²) on 0 ≤ y ≤ s. Hint: substitute y = s sin u.
(b) A linearly tapered wing with taper ratio λ = 0.5 carries a lift proportional to its chord, c(y) = c_r(1 − (1 − λ)y/s). Derive its centroid ȳ = s(1 + 2λ)/(3(1 + λ)).
(c) Schrenk's approximation takes the average of the elliptic and the chord-proportional distributions, each scaled to the same total lift. For s = 5 m and 5000 N on the semispan, compute the resultant and its position with `distributed_load_resultant` on 41, 161, 641 and 2561 points, and compare with the average of the two closed-form centroids.
(d) What root bending moment does each distribution predict? Why does the designer of a wing spar care which is right?

Answer form: two derivations, a table, two moments and a sentence.

## 6. (Judgment: what the free-body diagram left out)

Answer each in two or three sentences.

(a) The analysis treated the touchdown as a static load. Name two dynamic effects that change the brace load and say which way each pushes it.
(b) A retraction actuator is attached to the drag brace halfway between its pins. What changes in the free-body diagram of the brace, and in the analysis of the leg?
(c) The trunnion was modelled as a frictionless pin. If it were a bolted joint that transmits some moment, would the brace force go up or down? Argue from the FBD.
(d) A colleague computes the brace force by taking moments about the contact point C instead of A, and gets a different answer. Without seeing their work, list the two most likely mistakes.

Answer form: four short paragraphs.
