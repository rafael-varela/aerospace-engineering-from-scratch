# P08.18  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. Say which results are exact, which are numerical, and on which grid. Every number you report from the solver comes with the grid it was computed on.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: the wall and the cell)

(a) For the 33 × 33 grid at Re = 100, compute h, the cell Reynolds number U h Re, and the lid term −2U/h of Thom's formula. Repeat the cell Reynolds number for the 17, 65 and 129 grids. Which grids are outside the range U h Re ≲ 2?
(b) The function ψ = x²(1 − x)²y²(1 − y)² vanishes with its normal derivative on every wall of the unit square. Find the exact vorticity on the bottom wall at x = 0.5. Evaluate Thom's estimate there with `wall_vorticity` on the 17, 33 and 65 grids, and the observed order of its error.
(c) Explain in two sentences why a first-order wall formula does not spoil the second-order convergence of the interior solution.

Answer form: three numbers and a list; an exact value, three estimates and an order; two sentences.

## 2. (Derivation: from Navier–Stokes to ψ and ω)

(a) Starting from the 2-D incompressible Navier–Stokes equations (1)–(2) of the lesson, take the curl and derive the vorticity transport equation (4). Show each step where continuity is used, and say why the pressure disappears.
(b) Show that u = ∂ψ/∂y, v = −∂ψ/∂x satisfy continuity identically and that ∇²ψ = −ω.
(c) Derive Thom's formula (5) for a wall moving at speed U_w, and its truncation error.
(d) Prove that the discrete divergence of (D_yψ, −D_xψ), with centred differences on a uniform grid, is zero at every interior node for any ψ that vanishes on the walls. Then show that the trapezoid-rule flux ∫₀¹ u dy through the vertical centreline equals (ψ_{N−2} − ψ_1)/2 + hU/2 and is O(h²).

Answer form: four derivations.

## 3. (Numerical, uses your code: a grid-convergence study at your Reynolds number)

Take Re = 50 + s.

(a) Solve on the 17, 33, 65 and 129 grids. Tabulate ψ_min (with `vortex_centre`), the vortex-centre position and u at the centre node.
(b) Compute the observed order, the Richardson-extrapolated value and the grid-convergence index of u(0.5, 0.5) and of ψ_min on 33/65/129, and the order on 17/33/65.
(c) How much does the vortex centre move between the 65 and 129 grids? Is the 65 × 65 grid adequate for a 1 % requirement on ψ_min?

Answer form: a table; two orders, two extrapolates and two GCIs per triple; a distance and a yes or no with the number.

## 4. (Unfamiliar scenario: a faster stream over the cavity)

At Re = 1000 the cavity is the 2-D model of a recess under a stream ten times faster (or ten times larger) than at Re = 100.

(a) Call `solve_cavity_newton(65, 1000.0)` from the Stokes start. What happens, and why can Newton's method fail from a poor first guess even when a solution exists?
(b) Use continuation: solve at Re = 100, then pass that solution as `initial` for Re = 400, then for Re = 1000. Report the iterations, ψ_min and the vortex centre at each step.
(c) Repeat (b) on the 33 and 129 grids. What happens on 33 × 33 at Re = 1000, and why (use the cell Reynolds number)?
(d) Compare ψ_min at Re = 1000 on the 65 and 129 grids. Is either grid converged? What does that say about using the Re = 100 grid study to choose a grid at Re = 1000?

Answer form: a sentence and a reason; a table; a sentence and a number; two values and two sentences.

## 5. (Numerical, uses your code: what each route costs)

(a) Time `solve_cavity_newton` and `march_cavity` on the 17, 33 and 65 grids at Re = 100; record Newton iterations, marching steps and dt.
(b) Explain from the two stability limits of 2.4 why the step count is the same on the 17 and 33 grids but three times larger on the 65 grid.
(c) Estimate the number of marching steps and the time a 129 × 129 march would take, and compare with Newton's 2 s. State the assumption behind your estimate.
(d) Show that the two routes reach the same discrete solution on the 33 grid, to a tolerance you justify.

Answer form: a table; an argument with the two limits; a number of steps, a time and an assumption; a number and a justification.

## 6. (Judgment: corners, published data and what to claim)

(a) Print the vorticity at the wall node next to the top-right corner, ω[N − 2, N − 1], and at the middle of the lid, ω[N//2, N − 1], on the 17, 33, 65 and 129 grids at Re = 100. Which converges and which grows? Explain with the lid's velocity jump.
(b) The maximum difference between your centreline u and Ghia's Table I is 0.0023 on the 65 grid and 0.0043 on the 129 grid. A colleague concludes the 129 solution is worse. Using your grid-convergence results, argue for or against.
(c) Your team wants to use this solver to estimate the recirculating flow in an open wheel well at Re = 10⁶ based on the well depth. List what the solver can and cannot tell them, and what you would do instead.

Answer form: two series and an explanation; an argument with numbers; two lists and a recommendation.
