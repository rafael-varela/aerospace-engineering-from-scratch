# P03.14  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which are numerical, and what the numerical error is.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: steepest descent on a quadratic by hand)

Let f(x) = ½(x₁² + κx₂²) with κ = 10 + s, starting from x₀ = (κ, 1).

(a) Compute the gradient, the exact line-search step α₀ = gᵀg/(gᵀQg) and the first iterate x₁ by hand. Show that x₁ = ((κ − 1)/(κ + 1)) × (κ, −1): the iterates shrink by the same factor and flip side every step.
(b) Hence show that f_{k+1}/f_k = ((κ − 1)/(κ + 1))² exactly for this start, so the bound (1) is attained, and predict the number of iterations needed to reduce f by 10⁻¹⁰.
(c) Run `steepest_descent(..., exact_quadratic=Q)` from x₀ and confirm the ratio and the count. Then run `bfgs(..., exact_quadratic=Q)` and report its iteration count and final H.

Answer form: the derivation, one predicted and one measured count, the BFGS H.

## 2. (Derivation: the BFGS update)

(a) Starting from the Taylor expansion of ∇f, derive the secant equation H_{k+1}y = s and explain why it does not determine H_{k+1} for n > 1.
(b) Verify algebraically that (3) satisfies the secant equation and is symmetric, and prove that it is positive definite whenever H_k is and yᵀs > 0.
(c) Show that the strong Wolfe curvature condition |φ'(α)| ≤ c₂|φ'(0)| with c₂ < 1 implies yᵀs > 0.
(d) Compute (3) by hand from H₀ = I with s = (1, 0) and y = (a, 1), a = 1 + s/10. Check the secant equation and the eigenvalues, and compare with the one-variable update H = s/y.

Answer form: three derivations and one 2 × 2 matrix with its check.

## 3. (Code: the step size of a finite-difference gradient, and scaling)

Use the cruise objective `fuel_per_km` with the given aircraft at x = (0.70 + s/1000, 1.0).

(a) Compute a reference gradient by Richardson extrapolation of two central differences (h = 2 × 10⁻⁴ and 10⁻⁴). Tabulate the central- and forward-difference errors for h = 10⁻¹ … 10⁻¹⁰ and identify the best h for each. Explain the two slopes and the turning point from the truncation and rounding terms.
(b) Re-parameterise with altitude in metres, x = (M, h_m), and compute the gradient and the finite-difference Hessian at the same physical point. Report the Hessian's condition number for both parameterisations.
(c) Run `bfgs` on both parameterisations from (0.70, 8 km) with tol = 10⁻⁷ and max_iter = 500. Report iterations, convergence flag and the final f for each, and explain the difference.
(d) Propose a stopping test that would have declared the unscaled run converged, and say what it costs in guarantees.

Answer form: two tables, two condition numbers, a comparison and a one-paragraph recommendation.

## 4. (Unfamiliar scenario: cruise-climb as fuel burns off)

The aircraft has burned fuel and now weighs W = 600 − 2s kN.

(a) Without running anything, predict from §2.7 the new optimum Mach, C_L and altitude, and the new fuel per km (ratio to the 600 kN value).
(b) Run `optimal_cruise` for the new weight and compare. Then tabulate the optimum for W = 600, 570, 540, 500 kN and fit h* against ln W.
(c) At the new weight, what is the penalty for flying "M 0.78 at 11 km"? For continuing to fly the 600 kN optimum point?
(d) Compute the finite-difference Hessian at the new optimum and give the 1 % window in Mach and altitude. Which of the two is the flight-planning constraint more likely to be binding on, and why?

Answer form: predictions with their justification, a table, two penalties, two half-widths and a paragraph.

## 5. (Code: the three methods on Rosenbrock)

Start from x₀ = (−1 − s/50, 1).

(a) Run steepest descent (tol 10⁻⁶), damped Newton and BFGS (tol 10⁻¹⁰), counting function and gradient evaluations with a wrapper. Tabulate iterations, evaluations and the final f.
(b) For each method tabulate the last five error ratios |e_{k+1}|/|e_k| and the order estimates from `convergence_order`. Classify each as linear, superlinear or quadratic and say why the order estimates at the very end are unreliable.
(c) Compare the final BFGS matrix H with the exact inverse Hessian at (1, 1). Why does BFGS converge to the right x without H converging to the right matrix?
(d) Run pure Newton (α = 1 always, no positive-definiteness check) from the same start. Report its iteration count and the largest f along its path, and say what this shows about the role of the line search.

Answer form: three tables, two explanations, one number pair.

## 6. (Judgment: what the optimiser cannot tell you)

(a) `optimal_cruise` returns a point and a positive-definite Hessian. List three things that are true of that point and three things an operator needs that it does not establish.
(b) The optimum fuel flow is 0.516 kg/s; the function is flat to 1 % over ±0.023 in Mach. An airline's cost index moves the cruise Mach up to 0.80. Estimate the fuel penalty from the Hessian (§2.7) and from a direct evaluation, and say which you would quote and why.
(c) Suppose the gradient came from flight-test data with 0.5 % noise in f instead of a smooth model. Which of the three methods would you use, with what step size, and what would you expect of BFGS's curvature condition? Justify with the rounding-error analysis of §2.6.
(d) The lighter aircraft of problem 4 has an unconstrained optimum above a typical ceiling. Write the operating problem as a constrained optimisation (variables, objective, constraints) and say which SciPy method you would reach for and what you would check in its output.

Answer form: two lists, two penalties, a method choice with reasoning, a problem statement.
