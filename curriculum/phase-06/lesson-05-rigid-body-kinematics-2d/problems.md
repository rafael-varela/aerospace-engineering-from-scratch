# P06.05  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Draw the linkage to scale before writing any equation.

Several problems use your personal seed **s**: the day of the month you were born (1–31). All dimensions are illustrative.

## 1. (Warm-up: classify, assemble, move)

(a) Classify by Grashof's criterion (ground r₁, crank r₂, coupler r₃, rocker r₄, in mm): (100, 40, 120, 80); (60, 100, 90, 110); (100, 70, 30, 90); (100, 80, 90, 95); (50 + s, 20, 60, 40).
(b) For the first linkage at θ₂ = (10s)°, find θ₃ and θ₄ on both branches with `position_circle_intersection`, and the transmission angle on each.
(c) On the branch +1, with ω₂ = 10 rad/s, find ω₃ and ω₄, the speeds of A and B, and check v_B = v_A + ω₃ × r_BA by hand.

Answer form: five names; four angles and two transmission angles; two rates, two speeds and a check.

## 2. (Derivation: the four-bar equations)

(a) Derive the loop closure (3) and its Jacobian (4), and show det J = r₃r₄ sin(θ₃ − θ₄). Interpret det J = 0 geometrically.
(b) Differentiate (3) once and twice to get (5) and (6). Write out the known right-hand side of (6) in full.
(c) Prove Kennedy's theorem for three bodies in plane motion, and use it to derive (7).
(d) Prove that a parallelogram linkage (r₁ = r₃, r₂ = r₄, open branch) has ω₃ = 0 and ω₄ = ω₂ at every position.

Answer form: four derivations.

## 3. (Unfamiliar scenario: a flapping-wing drive)

A small ornithopter flaps its wings with a crank-rocker driven by a motor at a constant 5 Hz: ground 55 mm, crank 15 mm, coupler 60 mm, rocker (the wing root lever) 25 mm.

(a) Check that the crank can turn fully.
(b) The rocker's extreme positions occur when the crank and the coupler are in line (extended and folded). Find both rocker angles by the law of cosines, and the flapping amplitude. Check with a full-turn sweep.
(c) Find the crank angles at the two extremes. Through how many degrees does the crank turn on the downstroke and on the upstroke, and what is the time ratio?
(d) Find the largest wing-lever angular velocity and acceleration over a cycle at 5 Hz, and the smallest transmission angle. Is the drive acceptable?

Answer form: an inequality; two angles and an amplitude; two angles, two travels and a ratio; three numbers and a sentence.

## 4. (Numerical, uses your code: your own servo arm)

Repeat the lesson's servo sweep with a servo arm of r₂ = 10 + s/3 mm (horn 25 mm, pushrod and ground 120 mm, neutral with the arm vertical).

(a) Classify it. Find the neutral pushrod and horn angles.
(b) Over ±45°: the deflection range, the gearing at neutral and its range, the best straight-line gearing and the worst departure, the minimum transmission angle.
(c) Check ω₄ and α₄ at servo +20° (ω₂ = 2 rad/s, α₂ = 3 rad/s²) against central differences with h = 10, 5 and 2.5 ms, and report the observed orders.

Answer form: a class and two angles; six numbers; a table and two orders.

## 5. (Judgment: choosing the horn)

A UAV elevator needs ±25° or more of travel, a gearing that varies by no more than 25 % over the travel, a transmission angle above 40° everywhere, and the servo's rated torque is 0.35 N·m against a hinge moment up to 0.3 N·m.

(a) From the lesson's horn-length sweep and your own runs, which horns between 15 and 40 mm meet each requirement?
(b) The servo's gear train has 1° of backlash. How much surface free play does that give for each acceptable horn? Which way does that push the choice?
(c) Why might a designer avoid the exactly linear 20 mm horn, even though it meets the gearing requirement best?
(d) Recommend a horn length in three sentences.

Answer form: a table of pass/fail; free play per horn and a direction; a reason; a recommendation.
