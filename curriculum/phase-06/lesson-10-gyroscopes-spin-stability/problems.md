# P06.10  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: nutation of an oblate satellite)

The lesson's oblate satellite (I_t = 800, I₃ = 1200 kg·m²) spins at ω₃ = (0.5 + s/20) rad/s with a transverse rate of 0.02 rad/s.

(a) By hand: H, the body-frame nutation rate λ, the inertial coning rate H/I_t and the half-cone angle θ.
(b) The periods of the two motions. Which would an observer on the ground see as the "wobble" of the satellite's axis?
(c) Check the coning rate with `simulate_rigid` and `axis_coning_rate` over 40 s.

Answer form: four numbers; two periods and a sentence; a number.

## 2. (Derivation: energy, dissipation and the damper)

(a) For an axisymmetric body with angular momentum H whose axis is at θ from H, derive (6). Then prove, for any rigid body, T ≥ H²/(2I_max), with equality only for spin about the axis of largest moment.
(b) Derive the energy-sink rate (7) and show that it is positive for a prolate body and negative for an oblate one when P > 0.
(c) From (8), show that H = Iω + Jω_d is constant in inertial axes and that d/dt(½ω·Iω + ½J|ω_d|²) = −c|ω_d − ω|².
(d) Linearise (8) about co-rotating spin Ω and derive (9). Find its eigenvalues for c = 0 and interpret both.
(e) Show that for c → ∞ the nutation eigenvalue tends to i(I₃ − I_t)Ω/(I_t + J), and explain why dissipation then stops again.

Answer form: five derivations.

## 3. (Unfamiliar scenario: a spinning upper stage)

A solid-motor upper stage with its payload is spun to 60 rpm before ignition. Its moments are I₃ = (100 + 5s) kg·m² about the spin axis and I_t = 3I₃. During the 80 s burn the 60 kN thrust passes 1 mm from the centre of mass, a torque fixed in the body.

(a) The torque and H. The largest wander of H by (5).
(b) Simulate the burn with `simulate_rigid` and a constant body torque. The largest wander of H and the largest tilt of the spin axis from its starting direction. Why is the second larger?
(c) How far would the same torque turn a non-spinning stage in 80 s (from rest, about a transverse axis)?
(d) The stage is prolate. If its propellant slosh gives a nutation e-folding time of 300 s, should the designer worry during the burn? After it?

Answer form: two numbers and a number; two numbers and a sentence; an angle; two sentences.

## 4. (Numerical, uses your code: a gyroscope)

The oblate satellite spins purely about its axis at 1 rad/s (H = 1200 N·m·s), tilted (20 + s)° from the inertial z axis. A torque of (1 + s/10) N·m is kept perpendicular to H and to z (`gyroscope_torque`).

(a) The precession rate from (4) and its period.
(b) Simulate 100 s and fit the precession rate. Confirm that |H| and the tilt stay constant.
(c) Compare with the estimate τ/H that ignores the tilt. Which is right, and why does the tilt matter?

Answer form: two numbers; a number and two drifts; a number and a sentence.

## 5. (Numerical, uses your code: tuning a damper)

The lesson's prolate spinner (I_t = 1500, I₃ = 600 kg·m², 3 rad/s) carries a Kane damper with J = 60 kg·m².

(a) The nutation growth rate from `kane_nutation_mode` for c = 10, 30, 60, 100, 200 and 300 N·m·s, and the c that maximises it.
(b) For c = 30 + 2s N·m·s, simulate 900 s from a transverse rate of 10⁻³ rad/s. Report the fitted early growth rate, the times at which the axis passes 45° and 89° from H, and the final energy against H²/(2(I_t + J)).
(c) Predict from the linear rate alone how long the transverse rate takes to grow from 10⁻³ to 0.2 rad/s. Compare with the simulation's time to 45°.
(d) The same damper is mounted on an oblate spinner as a nutation damper. What value of c would you choose, and why?

Answer form: six numbers and a number; four numbers; a time and a sentence; a number and a sentence.

## 6. (Judgment: spin stability in practice)

For each, give a short argued answer (one to three sentences).

(a) A student says Explorer 1 tumbled because "the rigid-body equations are wrong". Correct the statement.
(b) Why is a spin-stabilised spacecraft usually built oblate (wider than it is long)?
(c) A long, thin launcher stage is spun for a short burn even though it is prolate. When is that acceptable?
(d) A colleague adds a nutation damper to a prolate spinner "to make it more stable". What happens?
(e) Which result of this lesson does not depend on the dissipation mechanism, and which does?
(f) The lesson's Explorer-like cylinder took 42 s to reach a flat spin in simulation. Why must that number not be quoted as Explorer 1's real tumbling time?

Answer form: six labelled items.
