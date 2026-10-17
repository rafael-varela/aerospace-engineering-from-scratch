# P06.02  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Draw the free-body diagram in path coordinates before writing any equation.

Several problems use your personal seed **s**: the day of the month you were born (1–31). All aircraft numbers are illustrative.

## 1. (Warm-up: a pull-up and a turn at the same g)

A trainer flies at V = 150 + s m/s.

(a) What radius gives 6 g at the bottom of a pull-up? If the pilot held that radius and speed over the top of a loop, what would the load factor be there?
(b) For a sustained level turn at 6 g and the same speed, find the bank angle, the radius and the turn rate (in °/s).
(c) The pilot and gear are 85 kg; the pilot's head and helmet are 6 kg. What force does the seat exert on the pilot, and what force does the neck carry, at 6 g?
(d) What is the slowest speed at which a loop of radius 1000 m can be flown with n ≥ 0 at the top?

Answer form: a radius and a load factor; three numbers; two forces; a speed.

## 2. (Derivation: path equations and the dive recovery)

(a) From the free-body diagram (weight, lift, thrust, drag), derive equations (1) of the lesson. State where you assume the thrust acts along the velocity.
(b) Show that an accelerometer along the aircraft's normal axis reads n = L/W in units of g. What does it read in level flight, at the top of a 0 g parabola, and in a 60° bank level turn?
(c) Derive the altitude loss (3) and the time (4) of a constant-speed, constant-n recovery. Show that for a small dive angle (3) reduces to the drop of a circular arc, and find the radius of that arc.

Answer form: two equations with a stated assumption; a derivation and three numbers; two derivations and a limit.

## 3. (Unfamiliar scenario: a zero-g research flight)

A research aircraft flies weightless parabolas. It pulls up at 1.8 g at 160 m/s, then pushes over to n = 0 at γ = +45° and holds n = 0 until γ = −45°. During the zero-g arc the engines are set so that thrust equals drag.

(a) What is the radius of the 1.8 g pull-up at the bottom (at 160 m/s)?
(b) With n = 0 and T = D, show that the path is a projectile parabola. Find the duration of the weightless arc, the height gained from entry to the apex, the horizontal distance covered and the exit speed, by hand.
(c) Check (b) with `flight_path_rhs` and `solve_ivp`, stopping with an event at γ = −45°.
(d) What is the radius of curvature of the path at the apex, and why is it not infinite even though n = 0?

Answer form: a radius; four numbers; a comparison; a radius and a sentence.

## 4. (Numerical, uses your code: your own dive recovery)

A fighter dives at γ₀ = −(30 + s)° at 250 m/s and pulls to 7 g, holding its speed.

(a) Compute the closed-form altitude loss and time.
(b) With `recovery_altitude_loss`, find the altitude loss for onset rates of 2 and 6 g/s. Compare the extra loss with the delay estimate (5).
(c) Repeat (b) at 6 g/s with thrust equal to drag instead of constant speed. How much more altitude, and why?
(d) Show with `kinematic_load_factor` that a simulated ramp at 6 g/s is reproduced from the positions alone.

Answer form: two numbers; two losses and two comparisons; a loss and a sentence; a plot or a largest error.

## 5. (Judgment: setting a recovery trigger)

An automatic recovery system must start the pull before the aircraft is lower than the predicted altitude loss plus a margin.

(a) For the lesson's 60° dive at 6 g/s, how much does the predicted loss change for a 10 % error in airspeed? For one g less than planned (6 g)?
(b) Which assumption of the model would make the real loss *larger* than predicted, and which smaller? Name at least two of each.
(c) At what speed would the illustrative wing (S = 27.9 m², C_L,max ≈ 1.4, ρ = 0.909 kg/m³) no longer be able to pull 7 g for the 12 000 kg aircraft? What does that do to the trigger at low speed?
(d) Propose a margin and justify it in two sentences.

Answer form: two percentages; two lists; a speed and a sentence; a margin with a justification.
