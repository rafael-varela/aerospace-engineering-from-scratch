# P10.20  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact (closed form), which are numerical (the lattice), and which are empirical (a fitted formula). Axes: x downstream, y to starboard, z up, ground plane z = 0.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: one horseshoe and its image)

A single horseshoe vortex of span b_v and circulation Γ flies at height h above the ground.

(a) Sketch the Trefftz plane (looking downstream): the two trailing vortices, their images, and the senses of rotation.
(b) Show that the downwash at the centre falls from Γ/(π b_v/2) to that value times (4h/b_v)²/(1 + (4h/b_v)²). Evaluate the ratio at h/b_v = 0.05, 0.1, 0.25, 0.5, 1.
(c) An elliptic wing of span b is often replaced by one horseshoe of span b_v = πb/4. Evaluate the ratio at h/b = 0.1 and compare with the classical empirical estimate of Wieselsberger, 1 − (1 − 1.32 h/b)/(1.05 + 7.4 h/b). Which is closer to the lattice result of the lesson (0.531), and why is the single horseshoe so pessimistic?

Answer form: a sketch; a derivation and five numbers; two numbers and a sentence.

## 2. (Derivation: why images work, and what they change)

(a) Show that a two-dimensional point vortex of strength Γ at (y₀, h) and one of strength −Γ at (y₀, −h) induce zero vertical velocity everywhere on z = 0. Extend the argument to a straight three-dimensional vortex segment parallel to the plane and to one that is not, and say why reflecting every point of a horseshoe and reversing its circulation satisfies the ground condition exactly.
(b) The Trefftz-plane induced drag is Di = (ρ/2) ∫ Γ(y) w(y) dy with w the far-wake downwash. Explain why the image wake reduces w at a given Γ(y), and why that reduces Di at a given lift.
(c) At a fixed angle of attack the lift rises in ground effect. Explain why, from the tangency condition at the control points. Would the lift rise if the wing were held at fixed circulation instead?

Answer form: a proof; two short arguments with a formula each.

## 3. (Unfamiliar scenario: an airliner's flare)

A large twin-engine airliner has a wing of about 60 m span and 7 m mean chord; at touchdown its wing sits about 5 m above the runway. Model the wing as a flat rectangle 60 × 7 m with `wing_in_ground_effect` (ns = 40, nc = 4) at α = 5°.

(a) Compute CL/CL_free and the induced-drag factor ratio K = (CDi/CL²)/(CDi/CL²)_free at h = 5, 10 and 20 m. Compare K with Wieselsberger's formula.
(b) The pilot holds the approach speed while descending through 20 m to 5 m. Describe what happens to the drag and to the angle of attack needed for the same lift, and why airliner pilots describe the aircraft as "floating" in the flare.
(c) Name two features of the real wing (not in your flat rectangle) that change the answer, and say in which direction.

Answer form: a table of six numbers and three empirical values; a paragraph; two features with directions.

## 4. (Numerical, uses your code: your own flare)

The reference trainer (750 kg, S = 12.0 m², b = 10.0 m, c = 1.2 m) approaches at 30 m/s. Its wing is at h = 0.5 + s/20 m above the runway at the flare.

(a) Compute CL/CL_free at fixed α and the induced-drag ratio K at your h, and the induced drag in free air and at your h for the lift that supports the aircraft.
(b) How much does the angle of attack needed for that lift fall (flat-plate lattice, angles from zero lift)?
(c) If the drag saving were not removed by reducing power, by how much would it change the aircraft's deceleration? Comment on the landing distance.
(d) Check that your lattice result has converged: repeat (a) with ns = 20 and 80 and report the spread.

Answer form: two ratios and two forces; two angles; an acceleration and a sentence; three values of K.

## 5. (Numerical, uses your code: the tail in ground effect)

The trainer's horizontal tail (2.0 m², assume aspect ratio 4, span 2.83 m) is 4.5 m behind the wing's quarter chord.

(a) Compute dε/dα at the tail in free air for tail heights 0 and 0.5 m above the wing plane, and compare with 2CL_α/(πAR). Why is the lattice value lower than the far-field formula even far downstream?
(b) Repeat in ground effect with the wing at h = 1, 2 and 5 m (tail at the wing's height).
(c) With a tail volume coefficient V_H = S_t l_t/(S c) and the tail's own lift slope from `wing_in_ground_effect` on a 2.83 m × 0.707 m rectangle, estimate the change in the aircraft's pitching-moment coefficient at α = 12° when the wing descends from free air to h = 1 m (tail efficiency 1). Which way does the nose go, and what must the pilot do?

Answer form: two gradients and a sentence; three gradients; a ΔCm, a direction and an action.

## 6. (Judgment: where the image model fails)

The lattice in this lesson is flat, thin and inviscid, the ground is flat and stationary relative to the air, and the wing does not move vertically.

(a) Measured ground effect on some thick or flapped wings very close to the ground shows a *loss* of lift at fixed α. Which physical mechanism, absent from a flat lattice, could cause it? (Recant, NACA TN-705, tested flapped wings near a ground board.)
(b) In a wind tunnel, a ground board has a boundary layer that the real runway does not have. What does it do to a measurement of ground effect, and how do tunnels avoid it?
(c) An aircraft descending at 3 m/s through ground effect sees a changing image system. Is the steady result of this lesson an upper or lower bound on the cushioning, and why?
(d) A WIG craft flies over waves 1 m high at h = 1.5 m. Which assumption fails first, and what is the hazard?

Answer form: four short arguments.
