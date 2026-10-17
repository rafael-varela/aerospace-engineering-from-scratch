# P09.14  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which come from your code, and which are model predictions.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: one axis by hand)

One satellite axis: x = (θ, ω), A = [[0, 1], [0, 0]].

(a) For B = (0, 1/J)ᵀ and for B = (1, 0)ᵀ (a hypothetical actuator that changes the angle directly), compute W_c and its rank by hand. (b) For C = [1, 0] (star tracker) and C = [0, 1] (gyro), compute W_o and its rank, and the unobservable direction where there is one. (c) Apply the PBH test at λ = 0 to each of the four cases and confirm the verdicts. (d) Check everything with `ctrb`, `obsv` and the two PBH functions.

Answer form: four ranks with matrices, one direction, four PBH verdicts, a code check.

## 2. (Derivation: the tests and the minimum-energy input)

(a) Using Cayley–Hamilton, show that the set of states reachable from x₀ = 0 is the column space of W_c. (b) Prove the PBH test: if wB = 0 and wA = λw for some w ≠ 0, then (A, B) is not controllable; show what w x(t) does. (c) Derive AW(T) + W(T)Aᵀ + BBᵀ = e^{AT}BBᵀe^{AᵀT} from (4) of the lesson, and the Lyapunov equation (5) for stable A. (d) Show that the input (6) of the lesson reaches x_f, and that any other admissible input has at least as much energy. (e) For the double integrator with B = (0, 1/J)ᵀ, derive W(T) in closed form and the rest-to-rest input (7) and its energy, and the bang-bang energy for the same slew.

Answer form: five derivations.

## 3. (Numerical, uses your code: choosing the pyramid angle)

The four-wheel pyramid of the lesson has its spin axes at β from +z. Take β = 40 + s/2 degrees, and also β = 54.74°, 30° and 70°.

(a) For each β, compute the eigenvalues of WWᵀ with all wheels and with wheel 0 failed, and the smallest eigenvalue of W(20 s) for both. (b) Sweep β from 5° to 85° and find the angle that maximises the smallest eigenvalue of WWᵀ after a single failure. Compare with cos²β = 1/3. (c) Explain why a shallow pyramid (small β) is good for z-axis slews and bad after a failure.

Answer form: a four-row table; one angle with a comparison; two sentences.

## 4. (Unfamiliar scenario: rendezvous with a failed thruster)

The in-plane Clohessy–Wiltshire model of `P02.05` (x radial, y along-track, n = 1.1314 × 10⁻³ rad/s) is A = [[0, 0, 1, 0], [0, 0, 0, 1], [3n², 0, 0, 2n], [0, 0, −2n, 0]].

(a) The chaser's radial thruster has failed, leaving only along-track thrust, B = (0, 0, 0, 1)ᵀ. Is the in-plane motion controllable? (b) If instead only the radial thruster works, B = (0, 0, 1, 0)ᵀ, what does the rank test say, and which eigenvalue does PBH report? Explain physically what cannot be done. (c) The relative-navigation sensor measures only one coordinate. With C = [1, 0, 0, 0] (radial only) and C = [0, 1, 0, 0] (along-track only), which is observable? Give the unobservable direction where there is one, and explain it with the CW physics.

Answer form: a verdict with a rank; a rank, an eigenvalue and an explanation; two verdicts, one direction and one explanation.

## 5. (Judgment: what the yes/no answer hides)

(a) The rank test uses a tolerance of 10⁻⁹ relative to the largest singular value. Construct the orthogonal-wheel satellite with the z wheel at a fraction ε of its nominal torque, and find the largest ε for which `is_controllable` says False. Compare with the smallest eigenvalue of W(20 s) at ε = 10⁻³ and 10⁻⁶. What should an operator watch? (b) A single star tracker measures roll about its boresight several times less accurately than the cross-boresight angles. The lesson's model says "unobservable". Is that too pessimistic or too optimistic for designing an estimator, and what would you compute instead? (c) The pyramid stays controllable after any single failure. List two things controllability does not guarantee about the surviving configuration, with the lesson's numbers where possible.

Answer form: one threshold, two eigenvalues and one recommendation; one paragraph; two items.

## 6. (Hard: a three-axis slew with a failed wheel)

The satellite must slew from rest at zero attitude to (5°, −3°, 8°) at rest in T = 60 s.

(a) Compute the minimum-energy input and its energy for the full pyramid, the pyramid with wheel 0 failed, the pyramid with wheel 2 failed, and the orthogonal set. (b) For each, give the peak torque of every wheel. Is the 10 mN m limit an issue? (c) Explain why losing wheel 2 costs almost nothing for this slew while losing wheel 0 costs three times the energy. Use the slew direction and the wheel axes. (d) Show that the minimum-energy distribution among redundant wheels at each instant is the minimum-norm (pseudo-inverse) split of the required body torque, and verify it numerically at t = 0 for the full pyramid.

Answer form: a four-row table of energies; a table of peak wheel torques with a yes/no; two sentences; one derivation and one numerical check.
