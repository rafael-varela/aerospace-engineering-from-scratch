# P06.03  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Sketch the force against payout before you compute it.

Several problems use your personal seed **s**: the day of the month you were born (1–31). All gear and aircraft numbers are illustrative.

## 1. (Warm-up: your aircraft on the lesson's gear)

An aircraft of mass m = 18 + s/5 t engages at v₀ = 55 + s/2 m/s.

(a) Its kinetic energy and momentum. The constant force, the deceleration in g and the time for an ideal 100 m runout.
(b) It lands on the lesson's constant-force gear (c = 100 N·s²/m², k = 4225 N/m, F₀ = 0), which was set for 20 t at 65 m/s. With your code find the runout, the stopping time, the peak force and where it acts, and the stroke efficiency.
(c) Compute the distance-average and the time-average force. Which theorem gives each?
(d) Check with `simulate_arrestment` and `energy_audit`: the impulse, and how the energy splits between the damper and the k x term.

Answer form: five numbers; four numbers and a location; two forces and two names; an impulse and two percentages.

## 2. (Derivation: the closed form and the monotonic force)

(a) From m dv/dt = −F and F = F₀ + kx + cv², derive the linear ODE (4) for u = v²(x), and solve it with an integrating factor to get (5). Check your answer by substitution.
(b) Prove that F(x) obeys F′ = k − (2c/m)F along the stroke, solve it, and conclude that the peak force is at an end of the stroke.
(c) Derive the constant-force design k = 2cF/m. Show that this k is half of the spring stiffness that would stop the aircraft alone in the same distance.
(d) Show that the stopping time ∫dx/v has an integrable singularity at x_s and that the substitution x = x_s(1 − s²) removes it.

Answer form: four derivations.

## 3. (Unfamiliar scenario: catching a drone on a bungee)

A 25 kg fixed-wing drone is recovered by flying its hook into a line attached to a bungee, at 25 m/s relative to the ground. The runout must be 15 m.

(a) With a bungee alone (F = kx), find k, the peak force and deceleration, and the stopping time. After the stop the bungee pulls the drone back: at what speed does it pass the engagement point going backwards?
(b) Add a damper whose force opposes the motion in both directions, F = kx + c v|v|, and halve k. Find the c that keeps the runout at 15 m. What design is this, and what is its peak deceleration?
(c) Integrate the full motion with `solve_ivp` (write your own right-hand side with the sign-correct damper) until the drone returns to x = 0. How fast is it going, and what fraction of its original kinetic energy comes back? Explain the number with the energy split.
(d) Recommend one change that stops the rebound.

Answer form: four numbers and a speed; a damper coefficient, a name and a deceleration; a speed, a fraction and a sentence; a design change.

## 4. (Numerical, uses your code: one gear, many aircraft)

Use the lesson's constant-force gear and your aircraft from Problem 1.

(a) Tabulate the runout, the peak g and where it acts for masses 0.8, 1.0 and 1.2 times yours and speeds v₀ − 5, v₀ and v₀ + 5 m/s.
(b) The deck allows 105 m. What is the heaviest aircraft that the gear can take at your v₀ (`brentq` on `stopping_distance`)?
(c) Design a constant-force gear for your own aircraft at 100 m. How different are c and k from the lesson's?

Answer form: a 3 × 3 table; a mass; two coefficients and a sentence.

## 5. (Judgment: choosing a force law)

You must recommend a force law of the form F₀ + kx + cv² for a gear that will take aircraft from 16 t to 24 t at 60–70 m/s.

(a) From the lesson's tables, which mix gives the lowest peak at design, which is most tolerant of speed errors, and which of weight errors?
(b) The constant-force design stores half the energy in its k x term. Why is that a problem, and how could a real gear avoid it?
(c) The damper-heavy design peaks at engagement. Which neglected effect of section 2.6 would make that peak worse in reality?
(d) Make a recommendation in three sentences, including what must be set before each landing.

Answer form: three designs; a reason and a mechanism; an effect; a recommendation.
