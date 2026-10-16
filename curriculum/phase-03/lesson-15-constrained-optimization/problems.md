# P03.15  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact (hand vertex) and which are numerical (iteration to a tolerance).

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: KKT by hand)

Minimise f(x, y) = (x − 3)² + (y − 2)² subject to x + y ≤ 2 + s/10.

(a) Write the constraint in the form g ≥ 0, the Lagrangian, and the four KKT conditions.
(b) Decide whether the constraint is active for your s and solve for (x*, y*, λ*). Interpret λ* as a shadow price: by how much does f* change if the right-hand side grows by 0.1?
(c) For which values of s is the constraint inactive? What are x*, y*, λ* then?

Answer form: three numbers with the KKT conditions written out; one sentence on (c).

## 2. (Derivation: prices, Newton, penalty)

(a) Let f*(ε) be the optimum of minimise f subject to g_i(x) ≥ −ε. Starting from the Lagrangian, show that d f*/dε = −λ_i* (equation (3) of the lesson). State where you use stationarity and where complementarity.
(b) Show that the KKT system of the quadratic program (6) with the active set 𝒲 is the Newton step (5) for the nonlinear KKT equations restricted to 𝒲. Which matrix must be nonsingular, and what does that require of the constraint gradients?
(c) For the quadratic penalty φ = f + (μ/2) Σ min(g_i, 0)², set ∇φ = 0 and compare with KKT stationarity to show that an active constraint is violated by about λ_i*/μ. Check the prediction against the lesson's run (λ = 0.464, μ = 100).
(d) Rewrite the strut problem in u = ln r, v = ln t, show it is a linear program, and derive the multipliers of the Euler + local vertex, λ_E = λ_l = m*/3, by converting the log-space multipliers to mass units.

Answer form: four short derivations; (c) and (d) end in numbers.

## 3. (Unfamiliar scenario: a propellant tank)

A cylindrical aluminium tank (ρ = 2700 kg/m³) with hemispherical ends must hold V = 2 + s/10 m³ at a pressure p = 0.3 MPa. Its radius R, cylinder length L and uniform wall thickness t must satisfy: hoop stress pR/t ≤ 180 MPa; t ≥ 1.2 mm (gauge); R ≤ 0.6 m (it must fit the fairing); L ≥ 0. Mass m = ρ t (2πRL + 4πR²).

(a) Eliminate L with the volume equation πR²L + (4/3)πR³ = V and write the problem in (R, t) with four normalised inequality constraints.
(b) With no radius limit, show that the fixed-gauge optimum is a sphere (L = 0), and find its radius for your V. Does it fit?
(c) Solve your problem by hand: which constraints are active, what are R, L, t and m, and what is the multiplier of the radius limit in kg per metre? What does a 10 mm wider fairing buy?
(d) Build the problem as a `Problem` and solve it with your `sqp`; confirm the design, the active set and the multipliers, and compare with SLSQP.

Answer form: a 2-variable problem statement, four numbers with units, the multiplier with its interpretation, one line of code output.

## 4. (Numerical, uses your code: your own strut)

Size the lesson's strut for P = (30 + s) kN and L = 2 + s/20 m, all other parameters unchanged.

(a) Predict from the vertex formulas which pair of constraints is active (check the two load thresholds of 2.6 for your L), then compute r*, t*, m* and the stress at the optimum by hand.
(b) Run `sqp` from (r, t) = (50 mm, 2 mm) with `hessian="bfgs"` and `"exact"`. Report iterations, the final KKT residual, the multipliers, and the observed order from the last three errors ‖x_k − x*‖ (use `analytic_optimum` for x*).
(c) Compute the shadow price d m*/dP from the multipliers and check it against `shadow_price_check`; then state the mass saved per kilonewton of load reduction.
(d) Compare with `scipy.optimize.minimize(method="SLSQP")`: iterations and distance to the vertex.

Answer form: a table (quantity, hand, SQP-BFGS, SQP-exact, SLSQP), the order estimates, two sentences.

## 5. (Numerical, uses your code: the price of a penalty)

Wrap the strut objective in a counter and run, from (50 mm, 2 mm): `penalty_minimize` for μ = 10, 10², …, 10⁶; `augmented_lagrangian` with μ = 1, 10 and 100; and `sqp`.

(a) Tabulate for the penalty method ‖x_μ − x*‖, the largest constraint violation and the number of objective evaluations. Fit the order of the error in 1/μ. Explain the growth of the evaluation count from the conditioning of the penalty Hessian.
(b) Tabulate outer iterations, final error and evaluations for the three augmented-Lagrangian runs. Why does μ = 100 need fewer outer iterations than μ = 10, and why is μ = 1 so slow?
(c) How many evaluations does SQP need? Give the smallest μ at which the penalty method matches SQP's constraint violation (10⁻⁸ or better), and the number of evaluations it would need, by extrapolating your table.

Answer form: two tables, a fitted order, three short explanations.

## 6. (Judgment: reading an optimiser's output)

A colleague sizes a wing box with 24 design variables (skin and spar thicknesses in metres) and 60 constraints (stresses in pascals, buckling loads in newtons) with SLSQP. It stops after 300 iterations reporting "Positive directional derivative for linesearch", and the returned point has one constraint violated by 3 × 10⁴ Pa.

(a) Name two scaling problems in the setup and say what each does to the BFGS Hessian and to the merit function.
(b) The violated constraint is a stress limit of 300 MPa. Is 3 × 10⁴ Pa a problem? How would you have known at a glance with normalised constraints?
(c) After rescaling, the run converges and reports a bound constraint with multiplier 2 × 10⁻⁹ kg/m. Is the bound active? Should it be relaxed?
(d) The two active constraints at the optimum are both buckling constraints whose knock-down factors are uncertain by ±30 %. The optimiser's answer is quoted as 412.3 kg. How many figures of that number would you defend, and what single calculation from this lesson supports your answer?

Answer form: four short paragraphs, each with at least one number.
