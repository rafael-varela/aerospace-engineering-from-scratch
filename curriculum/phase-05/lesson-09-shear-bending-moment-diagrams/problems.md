# P05.09  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer, state your sign convention, and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31). Use the lesson's convention throughout: forces and distributed loads positive up, couples positive counterclockwise, V the sum of upward forces to the left of the cut, M positive sagging, g = 9.80665 m/s². The aircraft numbers are illustrative of the light-twin class, not the data of any particular type.

## 1. (Warm-up: a simply supported equipment beam)

A beam in an avionics bay spans L = 3 + s/10 m between a pin at x = 0 and a roller at x = L. It carries a uniform downward load of 1.5 kN/m over its whole length and a downward point load of 4.0 kN at x = L/3.

(a) Find both reactions by hand and check that they sum to the total load.
(b) Write V(x) on 0 ≤ x < L/3 and on L/3 < x ≤ L. Where does V change sign? Find the maximum bending moment and where it occurs.
(c) Check (a) and (b) with `solve_reactions` and `extreme_moments`, and sketch V and M.

Answer form: two reactions, two expressions and a location, one moment with its position, a sketch.

## 2. (Derivation: the slice relations, a linear load and the elliptical wing)

(a) From the free-body diagram of a slice dx carrying a distributed load w, derive dV/dx = w and dM/dx = V in the lesson's convention. Then integrate across a point force P and across a counterclockwise couple M₀ to obtain the jump conditions.
(b) Derive the linear-load row of the singularity-function table (lesson 2.4): show that the load "w_a + k(x − a) starting at a, minus w_b + k(x − b) starting at b" vanishes beyond b, and integrate it twice.
(c) A simply supported span L carries a load rising linearly from 0 at x = 0 to w₀ (downward) at x = L. Show that V = 0 at x = L/√3 and that the maximum moment is w₀L²/(9√3).
(d) For the elliptical lift w(x) = w₀√(1 − (x/s)²), show that the total lift on the half-wing is πw₀s/4, the root bending moment is w₀s²/3, and the centre of lift is at 4s/(3π) from the root.

Answer form: four derivations.

## 3. (Unfamiliar scenario: a wing on trestles in the hangar)

For maintenance, the light twin's right wing is removed and laid on two trestles. Model it as a beam from x = 0 (the root rib) to x = 6.7 m (the tip): its structure weighs 80 kg, spread uniformly, and the engine and nacelle, left installed, weigh 150 kg at x = 1.9 m. The trestles are at x₁ = 1.0 m (a pin) and x₂ = 4.0 + s/20 m (a roller). Gravity only (n = 1); the fuel tanks are empty.

(a) Find the trestle reactions. Which trestle carries more, and why?
(b) Find V just to the left and right of each trestle and of the engine, and the bending moments there. Where are the largest sagging and hogging moments, and what are they?
(c) Keeping x₁ = 1.0 m, where should the outboard trestle go to make the largest |M| along the wing as small as possible? Use the code to sweep x₂.
(d) Compare the largest |M| on the trestles with the root bending moment in the 3.8 g pull-up of the lesson. Would you worry about the wing's strength on the trestles? What would you worry about instead?

Answer form: two reactions, a table of six values, two moments with positions, a position and a moment, a ratio and two sentences.

## 4. (Numerical, uses your code: your own pull-up)

Use `WingCase(n=3.0 + s/20, engine_x=1.5 + s/50)` with the other values of the lesson.

(a) Find w₀ and the exact root shear and root moment of the lift alone. Does the centre of lift depend on n? Why?
(b) Solve the lift alone with 20, 80 and 320 panels on a cantilever fixed at x = 0 and tabulate the relative error in the root moment. Estimate the order of convergence.
(c) With inertia relief and 400 panels, find the root shear and moment and the fraction of the lift-alone root moment that remains. Check the root moment with `sections_check` on the exact ellipse (a `DistributedFunction`) and report the difference.
(d) Report the shear just inboard and just outboard of the engine and the moment at the engine. Verify that the jump equals n times the engine's weight.

Answer form: three numbers and a sentence, a three-row table and an order, three numbers and a difference, three numbers and a check.

## 5. (Numerical and judgment: which tank to burn first)

Each wing of the light twin (aircraft mass 1900 kg with full fuel, n = 3.8, 400 panels) has two tanks: an inboard tank on 0.5–3.0 m and an outboard tank on 3.0–(4.5 + s/20) m, each holding 40 kg of fuel spread uniformly. Wing structure and engine are as in the lesson. Build the loads from `WingCase(mass=..., fuel_mass=0.0).loads(...)` plus your own `DistributedLoad`s for the tanks.

(a) Find the root bending moment with full tanks (1900 kg).
(b) After 80 kg of fuel has been burned (40 kg per wing; the aircraft mass is 1820 kg), find the root bending moment if the inboard tanks were emptied first, and if the outboard tanks were emptied first.
(c) Find the root bending moment with all wing fuel gone (1740 kg). Explain why it can exceed the full-fuel value even though the aircraft is lighter.
(d) The root moment with no wing fuel is linear in the aircraft mass. Find the aircraft mass with empty wing tanks that gives the same root moment as the full-fuel case. This is the idea behind a *maximum zero-fuel weight*: what does it limit, and why?
(e) Recommend a burn order and state its cost.

Answer form: four moments, one mass, two short paragraphs.

## 6. (Judgment: is the diagram telling the truth?)

Answer each in two or three sentences, with a number where one is asked for.

(a) The engine in the lesson is a point force of 5590 N at 1.9 m. In reality the pylon attaches with two fittings 0.6 m apart, each taking half. Where does the bending moment change, and by at most how much? Compare with the moment at the engine (about 39.7 kN·m).
(b) A colleague solves the lift alone with 10 panels and reports a root moment of 98.2 kN·m, calling it "conservative because it is a coarse model". Is it conservative? By how much does it differ from the exact value?
(c) The aircraft taxis over a bump at n = 1 + s/30 with no lift on the wing and the landing gear in the fuselage. Compute the root bending moment (wing structure, lesson fuel and engine only) and say whether it is sagging or hogging. Which spar cap is in compression?
(d) The real lift distribution of the tapered wing puts the centre of lift 5 % further outboard than the ellipse. By how much does the lift-alone root moment rise? Is that larger or smaller than the effect of moving the engine by 0.7 m?

Answer form: four short paragraphs, four numbers.
