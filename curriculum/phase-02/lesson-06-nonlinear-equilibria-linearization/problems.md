# P02.06  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: your own swing test)

A model aircraft of mass m = 3.0 + s/10 kg is hung from a pivot d = 0.30 m above its centre of mass. Its true pitch inertia about the centre of mass is I_cg = 0.10 + s/100 kg m².

(a) Compute I_pivot, the effective length L_eff, ω₀ and the small-angle period T₀.
(b) The technician releases it at 15°. Find the exact period, the inertia the small-angle formula infers from it, and the relative error, with `swing_test_error`.
(c) Compute the lever ratio I_pivot/I_cg and check that the error is close to (I_pivot/I_cg) θ₀²/8.
(d) What release amplitude keeps the inertia error below 1 %?

Answer form: four numbers with units in (a), three in (b), a ratio and a check in (c), an angle in (d).

## 2. (Derivation: the inverted equilibrium and the period series)

(a) Starting from θ̈ + 2ζω₀θ̇ + ω₀² sin θ = 0, write the first-order system, find its equilibria, and derive the Jacobian at a general θ.
(b) At the inverted equilibrium, show that the eigenvalues are λ = −ζω₀ ± ω₀√(1 + ζ²) for every ζ ≥ 0, that one is always positive, and give the time for a small tilt to double. Evaluate it for the UAV of the lesson with ζ = 0, 0.01 and 0.2.
(c) Starting from the energy integral, derive T = (4/ω₀) K(k) with k = sin(θ₀/2), stating the substitution.
(d) Using K(k) = (π/2)(1 + k²/4 + 9k⁴/64 + …) and the series of sin(θ₀/2), derive the θ₀² term of T/T₀ and verify the coefficient 1/16.

Answer form: a derivation with every step; three doubling times in seconds.

## 3. (Unfamiliar scenario: gravity-gradient libration of a satellite)

A satellite in a circular orbit with mean motion n (rad/s) and principal moments I_x, I_y, I_z pitches about its orbit-normal axis according to

θ̈ + 3 n² k sin θ cos θ = 0,   k = (I_x − I_z)/I_y,

where θ is the pitch angle of the body x axis from the local horizontal (body z points along the local vertical at θ = 0). Take an orbit period of 92.67 minutes (n = 2π/5560 s⁻¹) and k = 0.5 + s/100.

(a) Find all equilibria in [0, π] and classify each with the Jacobian. Which orientation of the satellite does the stable one correspond to?
(b) Give the small-amplitude libration frequency and period in minutes, and the ratio of libration period to orbit period. Why is it independent of the satellite's size?
(c) Show that with φ = 2θ the equation becomes a pendulum equation, and use `period_exact` to find by what factor the libration period at 30° amplitude exceeds the small-amplitude period.
(d) What would a negative k mean physically, and what happens to the classification?

Answer form: a table of equilibria with types, two periods, one factor, a short argument.

## 4. (Numerical, uses your code: identifying damping from a simulated swing)

Simulate the UAV swing test of the lesson with ζ = s/1000, released from rest at 20°, for twelve periods with `integrate_rk4` at h = 1 ms. Treat the output as a synthetic record of the test.

(a) Find the successive positive peaks of θ(t) and their times. From the first six peaks compute the logarithmic decrement δ = (1/5) ln(θ₁/θ₆) and the damping ratio estimate ζ = δ/√(4π² + δ²). Compare with the ζ you set.
(b) Measure the period with `period_from_trajectory`. Explain why it is longer than T₀ but shorter than the exact undamped period at 20°.
(c) Feed the measured period into `inertia_from_period`. What is the error in I_cg, and which effect dominates it: damping or amplitude?
(d) Propose and test a two-release procedure (release at 20° and at 10°) that recovers T₀ from the two measured periods using the θ₀² term of the series, and report the inertia error it leaves.

Answer form: six peak amplitudes, a decrement, a damping ratio, three periods, two inertia errors.

## 5. (Judgment: what limits a real swing test)

A real swing test times ten swings with a stopwatch (T uncertain by ±2 ms), measures d with a tape (±2 mm), and releases "by hand, about 10°" (amplitude between 5° and 15°).

(a) Differentiate the small-angle inertia formula to find ∂I_cg/∂T and ∂I_cg/∂d for the UAV of the lesson, and the relative error in I_cg each uncertainty causes.
(b) Using `swing_test_error`, give the inertia error at 5°, 10° and 15°.
(c) Rank the three error sources. Which one would you spend money to reduce, and how?
(d) Name one way the model of the lesson (rigid body, linear damping, planar swing) could be wrong in a way that none of these three uncertainties captures, and a measurement that would reveal it.

Answer form: two sensitivities with units, five percentages, a ranked list with a justification, one paragraph.

## 6. (Hard: the Van der Pol limit cycle and why the Jacobian cannot see it)

(a) For ẍ − μ(1 − x²)ẋ + x = 0 find the Jacobian at the origin and classify it for μ = 0.1, 1 and 3.
(b) With `limit_cycle_amplitude`, find the amplitude for μ = 0.1, 1 and 3 from a start inside (0.5, 0) and outside (3, 0) the cycle. (Use `t_end = 200` for μ = 0.1: the approach is slow.)
(c) Energy-balance argument for small μ: assume x ≈ A cos t over one cycle and compute the work done by the term μ(1 − x²)ẋ over the cycle, ∮ μ(1 − x²)ẋ² dt. Show that it is zero only for A = 2, positive for A < 2 and negative for A > 2, and explain why this makes the cycle attracting from both sides.
(d) Explain in two sentences why no linearisation at any point of the phase plane can predict the amplitude 2.

Answer form: three classifications, six amplitudes, one integral evaluated, two sentences.
