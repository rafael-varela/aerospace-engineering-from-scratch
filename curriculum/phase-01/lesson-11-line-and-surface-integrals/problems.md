# P01.11  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical, with the number of segments or cells used.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: one vortex, three loops)

A point vortex of circulation Γ = 10 + s m²/s (counter-clockwise) sits at the origin.

(a) Write the velocity at (2, 0) m and at (0, −2) m as vectors, with magnitude and direction.
(b) By hand, give the circulation around a circle of radius 3 m centred on the origin, around a circle of radius 1 m centred at (5, 0) m, and the outward flux through each.
(c) Check all four numbers with `circulation` and `flux_2d` using n = 2000 segments, and report the numerical error against (b) together with the prediction of equation (9) in the lesson.

Answer form: two vectors with units; four numbers; a short table of numerical against exact.

## 2. (Derivation: the circulation of a gradient field on a loop that does, and does not, enclose the singularity)

(a) Show that **v** = (Γ/2π)(−y, x)/r² equals ∇φ with φ = Γθ/(2π), θ = atan2(y, x), at every point except the origin. Compute both partial derivatives explicitly.
(b) Using the fundamental theorem for line integrals, show that the circulation around any closed loop that does not enclose the origin is zero. State exactly where you use the fact that the loop does not enclose the origin.
(c) Show that around any counter-clockwise loop that encloses the origin once the circulation is Γ, by deforming the loop onto a small circle or by following θ once around.
(d) What is the circulation around a figure-of-eight loop whose left lobe encloses the vortex counter-clockwise and whose right lobe is traversed clockwise and encloses nothing? And around a loop that encloses the vortex twice?

Answer form: derivation with each step named; one sentence each for (d).

## 3. (Unfamiliar scenario: a hurricane as a Rankine vortex)

A mature hurricane has its peak sustained wind V_max = 50 m/s at the eyewall, radius a = 30 km, and the wind outside the eyewall falls roughly as 1/r. Model it as a Rankine vortex.

(a) Find Γ and the vorticity inside the core (2V_max/a).
(b) Find the circulation around circles of radius 15 km and 60 km centred on the eye, and the wind speed at 100 km.
(c) Compare Γ and the core vorticity with the lesson's aircraft wake vortex (Γ = 300 m²/s, a = 1.5 m, peak swirl 32 m/s). Which has the larger circulation and which the larger vorticity, and by what factors?
(d) Check (b) with `rankine_vortex_velocity` and `circulation`, n = 2000.

Answer form: six numbers with units; two ratios; one paragraph.

## 4. (Numerical, uses your code: how many segments do you need?)

(a) Run `circulation_convergence` for a point vortex with Γ = 1 and n = 64, 128, 256, 512. Tabulate the circle and square errors, compute the log₂ convergence rates, and compare the circle column with Γ[(n/π) tan(π/n) − 1] and with its leading term π²/(3n²).
(b) Derive the closed form (9) from the geometry of a chord: locate the midpoint of each chord, evaluate the vortex speed there, and dot it with the chord.
(c) Repeat (a) for the lifting cylinder (V∞ = 50 m/s, a = 0.5 m, Γ = −20 m²/s) on a square of half-width 2 m. Report the relative error at n = 64, 256 and 1024 and the n needed to reach 10⁻⁴.
(d) Why is the circle error for the cylinder flow on a circle *centred on the cylinder* identical to the error for the bare vortex? Test your explanation by moving the circle's centre to (0.3, 0.2) m.

Answer form: two tables, a derivation, two numbers, one paragraph.

## 5. (Surface integral: inlet mass flow with a boundary-layer profile)

The fan face of radius R = 1 m at 10 km (ρ = 0.4135 kg/m³) sees an axial velocity w(r) = V(1 − r/R)^{1/7} with V = 250 m/s, a crude boundary-layer profile that is zero at the wall.

(a) By hand, compute ṁ = ∬ ρw dA as a fraction of the uniform-profile value ρVπR². (Substitute u = 1 − r/R.)
(b) Compute ṁ with `surface_flux` on the disc for ns = nt = 20, 40, 80, 160, 320. Tabulate, estimate the convergence order from successive differences, and explain why it is not 2.
(c) Propose and test a change of variable that restores second-order convergence.

Answer form: a fraction and a number in kg/s; a table with an order; one paragraph and one verifying number.

## 6. (Judgment: where the models stop)

For each statement, say whether it is true for the lesson's model, give one measurable sign that the model has failed in the real flow, and name the assumption that fails:

(a) "The circulation around a wing section is the same on every loop that encloses it."
(b) "The mass flow through the fan face equals ρ V∞ A, with V∞ the flight speed."
(c) "The work done by Earth's gravity on a satellite between two orbits does not depend on the path."
(d) "The velocity at the centre of a wake vortex is infinite."

Answer form: four short paragraphs, each with one number or one measurable quantity.

## 7. (Hard, uses your code: the wake vortex pair)

A transport aircraft on approach (W = 686.7 kN, b = 34 m, V = 70 m/s, ρ = 1.225 kg/m³) sheds two counter-rotating vortices of circulation ±Γ, Γ = 4W/(πρVb), separated by b′ = πb/4. Place the vortex with circulation −Γ (clockwise) at (−b′/2, 0) and +Γ at (+b′/2, 0).

(a) Compute Γ and b′.
(b) Each vortex moves with the velocity induced by the other. Show that the pair descends at Γ/(2πb′) and compute that speed.
(c) Build the pair as a superposition of two `vortex_velocity` calls and compute the circulation around a loop enclosing only the left vortex, only the right vortex, and both. Explain the third result.
(d) The velocity at the midpoint between the vortices is the downwash a following aircraft meets if it flies exactly between them. Compute it and compare with the descent speed of the pair; explain the factor between them.
(e) A light aircraft of span 10 m crosses the wake 2 km behind. Using the circulation on a loop of 5 m radius around one core, estimate the rolling velocity difference across its span and say whether a 20 m/s swirl is survivable at its 40 m/s approach speed.

Answer form: five numbers with units; two short explanations; one judgment with a number.
