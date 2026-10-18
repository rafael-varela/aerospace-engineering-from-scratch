# P10.15  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which numbers are measured (Ladson's data), which are fitted, and which are computed. Angles in degrees in answers, radians in code.

## 1. (Warm-up: the section model)

With the lesson's fit at Re = 2 × 10⁶ (a₀ = 5.741 /rad):

(a) Compute cl at α = 10° for f = 1, 0.5 and 0. What fraction of the attached-flow lift does a half-separated section keep?
(b) Compute cl and cd of the flat-plate model at 45°. What is the resultant force coefficient, and in which direction does it act?
(c) Why is (1) of the lesson useless at 60° while the flat-plate model is reasonable?

Answer form: three numbers and a fraction; three numbers and a direction; a sentence.

## 2. (Derivation: fits and iterations)

(a) Show that (1) reduces to thin-airfoil theory for f = 1 when a₀ = 2π, and that a fully separated section keeps a quarter of it.
(b) For a model linear in its parameters, cl = Xp, derive the normal equations and the covariance σ²(XᵀX)⁻¹. Explain why the same formula, with X replaced by the Jacobian, gives approximate standard errors for the nonlinear fit.
(c) Show that a converged decambering iteration satisfies the nonlinear lifting-line equation (4) exactly, whatever slope a₀ the linear solver uses.
(d) For a straight-line section, write the textbook fixed point as A ← (1 − ω)A + ω(JA + b). Show that it converges iff |1 + ω(λ − 1)| < 1 for every eigenvalue λ of J, and hence, if all λ are real and negative, iff ω < 2/(1 − λ_min).

Answer form: four derivations.

## 3. (Unfamiliar scenario: a tilt-wing in transition)

A tilt-wing aircraft's 12 m² wing passes through 35°, 45° and 60° angle of attack at 20 m/s while converting from hover to forward flight.

(a) With the flat-plate model, compute cl, cd, the lift, the drag and the normal force at each angle.
(b) Which component helps the conversion and which resists it? Why do tilt-wing designs bathe the wing in propeller slipstream during this phase?
(c) Which assumption of the flat-plate model would you distrust most for this wing, and why?

Answer form: a table; two sentences; a sentence.

## 4. (Numerical, uses your code: washing out the tapered wing)

The 0.45 taper stalls first at 56 % of the semi-span, over the ailerons.

(a) Add linear washout of 0, 2, 4 and 6° to the tapered wing and find, with `stall_onset_alpha`, the angle, spanwise station and C_L of the first local stall.
(b) How much washout brings the first stall inboard of 30 % of the semi-span? What does it cost in C_L at stall onset?
(c) Compute e at α = 5° for each washout with `lifting_line_linear`. Why does e fall so much, and why is e at a single α not the right measure of a twisted wing's induced drag?

Answer form: a table; a number and a cost; four numbers and a paragraph.

## 5. (Numerical, uses your code and authentic data: a higher Reynolds number)

(a) Refit the Kirchhoff model to the Re = 8.86 × 10⁶ block. Find the section's stall angle and cl_max.
(b) Repeat the stall-pattern study for the three planforms with that section (bracket α between 0° and 35°).
(c) Compare with the lesson's Re = 2 × 10⁶ results. What changed and what did not? Is the trainer's wing ever at Re = 8.86 × 10⁶?

Answer form: two numbers; a table; a paragraph.

## 6. (Judgment: what can this model say about hysteresis?)

(a) Repeat the lesson's up-and-down sweep (12–22°) with 8, 16 and 32 Fourier terms, and tabulate C_L on the way down at 12, 14 and 16°.
(b) Why do the results depend on the number of terms past stall, but not below it?
(c) A colleague proposes to "predict the hysteresis loop of the trainer's wing" from this model. Write a two-sentence reply.
(d) What would a model need, beyond a steady 2-D section curve, to predict post-stall wing behaviour credibly?

Answer form: a table; a paragraph; two sentences; two features.
