# P10.07  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. All coefficients are dimensionless; give angles in degrees in your answers but compute in radians. Say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a two-term design by hand)

A designer wants cl = 0.5 and cm_c/4 = −0.08 at the ideal angle.

(a) Find A1, A2 and α_i by hand from section 2.4 of the lesson.
(b) Write dz/dx(θ) and integrate it to z(θ) by hand. Find the maximum camber and its position (numerically is fine).
(c) Check (a) and (b) with `camber_from_cl_cm` and `camber_from_coefficients`, then run `thin_airfoil` forward at α_i. What are cl, cm, A0 and the zero-lift angle?

Answer form: three numbers; a formula, a camber and a position; four numbers.

## 2. (Derivation: the inverse relation and the a = 1.0 line)

(a) Starting from equation (1) of the lesson, show that if Δcp(θ) is known and finite at the leading edge, the camber slope is dz/dx = α_i + Σ Aₙ cos nθ with Aₙ = (2/π)∫₀^π (Δcp/4) sin nθ dθ.
(b) Derive the closure condition α_i = −½ Σ Aₙ Iₙ and evaluate Iₙ for odd and even n.
(c) For a uniform load, show Aₙ = cl_i/(πn) for odd n. Using Σ_{n odd} cos nθ/n = −½ ln tan(θ/2), derive dz/dx = (cl_i/4π) ln((1 − x)/x) and integrate it to the a = 1.0 mean line.
(d) Show that the tail of the series makes an error in z of order 1/N², and compare with the lesson's observed orders.

Answer form: four derivations.

## 3. (Unfamiliar scenario: a section for a flying wing)

A tailless aircraft needs a section with cm_c/4 = 0 at its design cl of 0.3, so that it trims without a tail.

(a) Design the two-term camber line. What are A1, A2 and α_i? Sketch the camber line: what is unusual about its rear part, and why does it need that shape?
(b) Thicken it to 12 % with `thicken` and run the panel method at α_i and at 4°. Compare cm with the NACA 2412's panel cm at the same angles.
(c) What is the thin-airfoil zero-lift angle of your section, and what does the zero cm cost in lift at a given angle compared with a 2412?

Answer form: three numbers and a sentence; a table of cl and cm; an angle and a sentence.

## 4. (Numerical, uses your code: your own roof-top)

Take the roof-top load of the lesson with its break at a = 0.5 + s/100 (uniform to x = a, linear to zero at the trailing edge, total cl = 0.5).

(a) Run `inverse_design` with 200 terms. Report α_i, cm_c/4 and the maximum camber and its position.
(b) Check it forward with `thin_airfoil` (200 terms, m = 20000) and report the largest load error on 2–98 % chord.
(c) Thicken to 10 % and run the panel method at α_i. Report cl, cm and the mean load error on 10–90 % chord. Which of (b) and (c) measures the code, and which measures the model?

Answer form: four numbers; two numbers; three numbers and a sentence.

## 5. (Numerical, uses your code: is the selection robust?)

(a) Draw your own trainer weights with `np.random.default_rng(s).dirichlet(20 * WEIGHTS["trainer"])` and recompute the trainer's decision matrix (the raw criteria of `select_sections("trainer", 0.3309, 1.1122)`). Which section wins, by how much?
(b) Remove the NACA 0012 from the candidates and recompute the nominal totals. Which totals change, and why does removing a section that does not win change anything?
(c) Compute the flip points of your weights for the moment criterion. Write a two-sentence recommendation that a reviewer could check.

Answer form: a winner and a margin; a table and an explanation; a weight and two sentences.

## 6. (Judgment: what is the matrix missing?)

The real reference trainer uses the NACA 2412, which ranks fourth in the lesson's matrix.

(a) Name three criteria a real designer would add, and for each say which later lesson (or which NACA data) provides it.
(b) The thin-airfoil inverse design of the lesson produced a section (a = 1.0, cl_i = 0.4) that scores about the same as the 2412. List two things you would need to know before building a wing with it, and how you would find each out.
(c) A colleague proposes adjusting the weights until the 2412 wins, "to match practice". What is wrong with that, and what should be done instead?

Answer form: a list of three with sources; two items with methods; a short argument.
