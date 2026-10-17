# P10.10  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. Angles in degrees in answers, radians in the formulas. Say which results are exact, which numerical, and at what N.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: an elliptic wing by hand)

A flat elliptic wing of span 10 m has AR = 6 + s/5 and a0 = 2π.

(a) By hand, find CL_α, CL at α = 4°, the induced angle and CDi.
(b) Check with `elliptic_wing`, `solve_monoplane` and `wing_coefficients` (N = 10). What are A3 and δ?
(c) By what factor does CDi change if the span doubles at the same area and CL?

Answer form: four numbers; two numbers and a check; a factor.

## 2. (Derivation: the monoplane equation)

(a) Starting from (1) and (2) of the lesson, derive the monoplane equation (5), using Glauert's integral ∫₀^π cos nφ dφ/(cos φ − cos θ) = π sin nθ/sin θ.
(b) Derive CL = πAR A1 and CDi = πAR Σ nAₙ² by orthogonality.
(c) Solve (5) exactly for the elliptic chord and constant α, and show e = 1 and (8).
(d) Show that the rolling-moment coefficient (positive right wing down) is Cl = (πAR/4)A2.

Answer form: four derivations.

## 3. (Unfamiliar scenario: a partial-span flap)

The trainer's flat rectangular wing (b = 10 m, c = 1.2 m) gets plain flaps over |y| < 3 m, 25 % chord, 10° down. In lifting-line terms the flapped sections have α_L0 = −τδ (`P10.08`: τ = 0.609); the others α_L0 = 0.

(a) Solve at α = 0 with N = 40, 80, 160. What is ΔCL? What happens to the convergence, and why?
(b) Compare ΔCL with two quick rules: the section increment times the flapped-span fraction, and the same times the ratio of the 3-D to the 2-D lift slope.
(c) Report the local cl just inside and just outside the flap edge. What does a real wing do at that edge that the model cannot?

Answer form: three numbers and an explanation; two numbers and a judgment; two numbers and a sentence.

## 4. (Numerical, uses your code: a convergence study of your own)

Take the flat wing of span 10 m and area 12 m² with taper λ = s/40.

(a) Run `convergence_study` at α = 5° with N = 10, 20, 40, 80 (reference 1280). Report the errors and observed orders. Do they match section 2.6?
(b) Report δ and e at N = 80. Where is your λ relative to the minimum of δ?
(c) Without running it, predict the observed order if you added a linear washout to the rectangular wing, and say why.

Answer form: a table and a sentence; two numbers and a sentence; an order and a reason.

## 5. (Data-based: the lift slope of a measured wing)

NACA TN-1422 measured three aspect-ratio-9 wings with taper 0.4 at Reynolds number 4.4 × 10⁶ (rows in `curriculum/phase-10/data/finite_wing_cla.csv`; measured dCL/dα = 0.085–0.086 per degree).

(a) Model the wing (span 10 m, AR 9, taper 0.4, with and without 2° linear washout) and compute CL_α per degree for a0 = 2π, 0.95 × 2π and 0.9 × 2π.
(b) Which a0 reproduces the measurement? Is that plausible for a 10 %-thick section at this Reynolds number? Does washout change the slope, and why not?
(c) The report's own calculated values are 0.083–0.085 per degree. List three reasons a lifting-line slope and a measured slope can differ.

Answer form: a table; a value and two sentences; three reasons.

## 6. (Judgment: when not to trust the lifting line)

(a) A student applies the lesson's solver to a delta wing of AR 2. Explain, with section 2.7, what will go wrong.
(b) The rectangular trainer wing's lattice result (`P10.11`) has a lift slope 5 % lower and e = 0.970 against 0.934 here. Which would you use for the trainer's induced drag, and why?
(c) Why does the lifting-line equation still give a finite answer for a wing whose chord jumps (a step in planform), even though the load must then have a kink?

Answer form: three short arguments.
