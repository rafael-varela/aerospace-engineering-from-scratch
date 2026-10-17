# P05.16  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: dummy loads on a cantilever)

A cantilever of length L = 2.0 m and EI = 50 kN·m² carries a uniform load q = (1 + s/10) kN/m.

(a) By Castigliano with a dummy tip force Q, find the tip deflection. Check against qL⁴/(8EI).
(b) With a dummy tip *moment* C, find the tip slope. Check against qL³/(6EI).
(c) With the unit-load method, find the deflection at mid-length. Check against the closed form q x²(6L² − 4Lx + x²)/(24EI) at x = L/2.
(d) Reproduce (a) and (b) with `castigliano_deflection` on a function U(Q, C) built with `strain_energy_bending`.

Answer form: two deflections and a slope with units, each with its check.

## 2. (Derivations: the theorems and a redundant)

(a) Prove Castigliano's second theorem δ_j = ∂U/∂P_j for a linear elastic structure, starting from Clapeyron's U = ½ΣP_kδ_k and δ = FP.
(b) Prove that F is symmetric by applying two loads in either order and comparing the work done (Maxwell).
(c) Derive the prop force R = 3qL/8 of a propped cantilever under uniform load from ∂U/∂R = 0, and the deflection at mid-span of the propped beam by the unit-load method (hint: any statically admissible unit-load moment will do, for example that of the released cantilever).
(d) Derive the braced-spar strut force (8) of the lesson, including the spar's axial energy.

Answer form: four derivations.

## 3. (Unfamiliar scenario: a curved frame segment)

A curved bracket is a quarter circle of radius R = (0.4 + 0.01s) m and EI = 15 kN·m², built in at A = (R, 0) and free at B = (0, R) (centre at the origin). A vertical load P = 2 kN acts downward at B. Neglect axial and shear energy.

(a) Write the bending moment at the point at angle φ, show that U = πP²R³/(8EI), and find the vertical deflection of B.
(b) With a dummy horizontal load at B, find the horizontal deflection of B.
(c) Compare with a straight cantilever of length R under the same tip load. Why is the curved bracket more flexible?
(d) The lesson neglected shear deformation. For a cantilever with a tip load, the shear energy ∫αV²/(2GA) dx adds αPL/(GA) to the tip deflection. Show that the ratio of shear to bending deflection is 3αEI/(GAL²), and evaluate it for a solid 25 × 150 mm 2024-T3 rectangle (E = 73.1 GPa, G = 28 GPa, α = 1.2) at L = 0.3, 1.0 and 5.0 m. When is neglecting shear unsafe?

Answer form: two deflections, a ratio, a short explanation, three ratios and a sentence.

## 4. (Code: the engine-mount truss under your own load)

Use the three-bar engine-mount truss of the lesson (nodes (0, 0), (0, 0.3), (0, 0.6), (0.5, 0.25) m, bars from each wall node to the engine node, 4130 tubes 19.05 × 0.889 mm, E = 205 GPa). The engine node carries 6.0 kN at an angle of (−90 + 2s)° from the x axis.

(a) Find the bar forces with `truss_force_method` and confirm them with `truss_stiffness_method`.
(b) Find the displacement of the engine node in the direction of the load by the unit-load method with a unit load along that direction, and compare with the stiffness result projected onto it.
(c) Compute the three support reactions from the bar forces and check global equilibrium: ΣF = 0 and ΣM = 0 about the origin.
(d) Multiply the middle bar's EA by 0.1, 10 and 100 and tabulate the bar forces. How does the redundant bar's share depend on its stiffness, and what is its limit?

Answer form: bar forces, two displacements, reactions with two residuals, a table and a sentence.

## 5. (Code: where should the strut go?)

Use `braced_spar` with the lesson's spar and strut properties, but with the strut's fuselage point at h = (0.8 + 0.02s) m below the spar.

(a) For the elliptic lift, sweep the strut position a from 2.0 to 4.0 m in steps of 0.05 m and find the position that minimises the largest bending moment anywhere on the spar. Report S, R and the tip deflection there.
(b) Repeat for a uniform lift of the same total. How much does the best position move, and why?
(c) At a = 3.0 m, compare S and the root moment for the two distributions.

Answer form: two optimal positions, a table of S, R, max |M| and tip deflection, and a short explanation.

## 6. (Judgment: the finite-element model disagrees)

A colleague's finite-element model of the braced wing reports a strut force of 25 kN at the lesson's load case, against your 31.5 kN.

Write a half-page plan for finding out who is right. List at least four modelling differences that could explain a 20 % lower strut force (and the direction each would push S), say which two you would check first and how, and name one measurement on a real aircraft that would settle the question.

Answer form: a structured plan of about 250 words.
