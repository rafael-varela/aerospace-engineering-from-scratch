# P04.03  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a docking)

A cargo vehicle of (8000 + 100s) kg approaches a 420 t station at (0.05 + s/1000) m/s relative to it.

(a) For a perfectly inelastic docking, the common velocity after capture (mm/s) and the energy the docking mechanism must absorb.
(b) The mechanism is too stiff and the collision has e = 0.2. Velocities after the contact; energy absorbed. What does the vehicle's crew see?
(c) Show that the energy absorbed in (a) is within 3 % of the vehicle's whole kinetic energy, and give the condition on the masses for "nearly all".

Answer form: two numbers; three numbers and a sentence; one inequality.

## 2. (Derivation: from pellets to the rocket equation)

(a) Derive dv = u dm/m for one throw, stating which frame u is measured in and where momentum conservation is applied. Repeat with u measured relative to the cart *before* the throw and show the two agree to first order in dm.
(b) Show that the n-pellet velocity (5) is a left Riemann sum of u∫dm/m and derive the error (7) from the trapezoid-rule correction. Check it at n = 3 for m₀ = 10 t, m_f = 2.5 t, u = 2942 m/s against `pellet_convergence`, and explain the size of the difference.
(c) Show that F_ext = d(mv)/dt, applied to a rocket that leaks propellant with zero relative speed in a frame where it moves at v, predicts a velocity change, and that the prediction depends on the frame. Then derive the correct equation (8).
(d) Derive the burnout altitude (10).

Answer form: four derivations; one comparison of two numbers.

## 3. (Unfamiliar scenario: an astronaut without a tether)

An astronaut with her suit has a mass of 120 kg and holds three tools of (2 + s/10) kg each. She is at rest 20 m from the station and can throw anything at 5 m/s relative to herself after the throw.

(a) Her speed if she throws all three tools at once, and if she throws them one at a time. How long does each take her to reach the station?
(b) Compare both with Tsiolkovsky for the same masses. Which is the upper bound and why?
(c) If instead she throws all three at 5 m/s relative to her *initial* rest frame, what is her speed? Why is this larger than the "relative to herself after the throw" answer, and which one describes a real throw better?
(d) She arrives at the station at the speed of (a, one at a time). Her suit's arms can absorb 30 J safely. Is that enough to stop her? (Treat the station as infinitely massive.)

Answer form: two speeds and two times; three numbers and a sentence; a speed and two sentences; a yes or no with a number.

## 4. (Numerical, uses your code: choosing the thrust-to-weight ratio)

Use the lesson's vehicle (m₀ = 10 t, m_f = 2.5 t, u = 2942 m/s), with drag k = 0.144 (1 + s/15) kg/m and H = 8.5 km.

(a) For T/W₀ = 1.2, 1.5, 2, 2.5, 3, 4, 5, 6 tabulate the burnout speed, the gravity loss, the drag loss and the maximum dynamic pressure q = ½ρ₀e^{−z/H}v² (ρ₀ = 1.225 kg/m³).
(b) Find the T/W₀ that maximises the burnout speed (`minimize_scalar`). What limits a real vehicle from flying at it?
(c) The structure is rated for 40 kPa. Find the largest T/W₀ that respects it (`brentq` on max-q) and the burnout speed there. How much burnout speed does the structural limit cost?
(d) Why does the optimum of (b) lie so far above real lift-off thrust-to-weight ratios? Name two things the model leaves out that push the real choice down.

Answer form: an eight-row table; a number and a sentence; two numbers and a difference; two sentences.

## 5. (Staging with two different stages)

A two-stage vehicle has 30 t of propellant in total and a payload of (200 + 10s) kg. Stage 1: u₁ = 2800 m/s, ε₁ = 0.08. Stage 2: u₂ = 3400 m/s, ε₂ = 0.12.

(a) Write `stack_dv` for this vehicle as a function of the first stage's propellant fraction f and find the best f and Δv. Compare with f = 0.5.
(b) The single-stage limits u ln(1/ε) of each stage alone. Why can the two-stage vehicle exceed both?
(c) Find the largest payload for which the best split gives 9.4 km/s, and the vehicle's gross lift-off mass. What fraction of the lift-off mass reaches orbit?
(d) Without computing: would moving propellant to the stage with the higher exhaust speed always help? Argue from (11).

Answer form: a fraction and two speeds; two speeds and a sentence; two masses and a percentage; two sentences.

## 6. (Judgment: momentum in practice)

For each, give a short argued answer (one to three sentences).

(a) A satellite's cold-gas thruster has u = 700 m/s; an electric thruster has u = 30 000 m/s. For a 100 m/s manoeuvre of a 500 kg satellite, how much propellant does each use, and why is the electric thruster not always the better choice?
(b) A sounding rocket's designers quote "thrust = 2 kN, burn time 10 s, so the impulse is 20 kN·s and the velocity change is 20 000/m". What is wrong with dividing by m?
(c) A ball dropped on a floor bounces to 64 % of its drop height. What is e, and is the collision elastic?
(d) A spacecraft vents leftover propellant through a nozzle at the end of its mission. Does that change its orbit? What if the vent is non-propulsive (two opposite nozzles)?
(e) The lesson's vertical ascent "reaches 127 km at 2.6 km/s". Is the vehicle in orbit? What is missing?
(f) Your simulation shows the total momentum of rocket plus exhaust drifting by 1 % over the burn in free space. Which mistakes would produce that?

Answer form: six labelled items.
