# P02.12  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: from a harmonic function to an analytic one)

Let u(x, y) = x³ − 3xy² + s·x.

(a) Show that u is harmonic.
(b) Find its conjugate harmonic v from the Cauchy–Riemann equations and write f = u + iv as a function of z.
(c) Evaluate f′(1 + i) by hand and with `complex_derivative`; check `cauchy_riemann_residual` there.
(d) Why is u + i(−v) *not* analytic?

Answer form: a two-line check, a function of z, one complex number, one sentence.

## 2. (Derivation: the Kutta condition and the Kutta–Joukowski theorem)

(a) Starting from the complex velocity (4) of the lesson, show that w = 0 at the trailing-edge point ζ = ζ₀ + Re^{−iβ} requires Γ = 4πUR sin(α + β). Keep the signs explicit: which sign convention for Γ makes this positive at positive α + β?
(b) Show that the two stagnation points on the circle sit at angles −β and π + 2α + β from the centre (measured from the +x axis). Where are they at α = −β?
(c) Derive Blasius' formula F_x − iF_y = (iρ/2)∮ w² dz from F = −∮ p **n** ds and Bernoulli, using the facts that n ds = −i dz on a counter-clockwise contour and that w dz is real on a streamline.
(d) Expand w at large |z| and use the residue theorem to obtain L′ = ρUΓ and D′ = 0. Explain in one sentence why the result does not depend on the shape of the body.

Answer form: four derivations.

## 3. (Unfamiliar scenario: a wake vortex near the ground)

A transport aircraft's trailing vortex, modelled as a line vortex of circulation Γ = 200 + 10s m²/s, is at height h = 10 + s m above a flat runway. The ground is a streamline, which the **method of images** achieves with a second vortex of opposite circulation at −ih.

(a) Write the complex potential and the complex velocity of the pair.
(b) Show that the vertical velocity on the ground is zero everywhere, find the horizontal velocity on the ground as a function of x, and give its maximum and where it occurs.
(c) Each vortex moves with the velocity induced by the other. Show the real vortex drifts horizontally at Γ/(4πh) and give the number. Which way does a vortex with clockwise circulation drift?
(d) Check with `circulation_and_flux` that a circle around the real vortex returns ∓Γ (state which sign and why) and that a circle enclosing both returns zero.

Answer form: two formulas, two numbers with a direction, two contour results.

## 4. (Numerical, uses your code: your own section)

Build the section with c = 0.25 m, x₀ = −c(0.05 + s/400), y₀ = c·s/300.

(a) Report R, β, chord, thickness and camber (as fractions of chord).
(b) At U = 50 m/s and α = 4°: Γ, L′ and C_l.
(c) The lift slope from α = 0 and 8°, against 2π(1 + 0.77 t/c); the zero-lift angle.
(d) Verify that the pressure-integrated lift equals ρUΓ to 0.1 % and that the drag is below 10⁻³ of the lift.

Answer form: a table of geometry and lift; four numbers for (c)–(d).

## 5. (Numerical, uses your code: how fast is too fast?)

The lesson section at α = s/4 degrees and U = 50 m/s.

(a) Γ, L′, C_l, and the lift from the pressure integral.
(b) The maximum surface speed, where it occurs, and the minimum C_p.
(c) The local Mach number at that point with a = 340 m/s, and the free-stream speed at which it would reach 0.3 (speeds in potential flow scale with U).
(d) Suppose the trailing-edge point in `airfoil_coordinates` were set to c(1 + 10⁻⁶) instead of exactly c. What would `airfoil_surface_velocity` return there, and what does that tell you about the Kutta condition as a numerical object?

Answer form: six numbers; a short paragraph for (d).

## 6. (Judgment: what the model may and may not say)

(a) List the four physical effects the potential-flow model omits, and for each say whether it makes the predicted lift too high, too low, or unaffected at α = 4°, with a rough percentage where you can.
(b) Which single feature of the measured C_p distribution of a real section at 4° would most clearly show the model's limit, and where on the chord is it?
(c) The model predicts C_l = 2.0 at 12°. Would you use that number for a wing-loading estimate? What would you use instead, and from where?
(d) A colleague argues that because the pressure integral and ρUΓ agree to five figures, "the lift is validated". Reply in three sentences distinguishing verification from validation.

Answer form: four short paragraphs, each with at least one number.
