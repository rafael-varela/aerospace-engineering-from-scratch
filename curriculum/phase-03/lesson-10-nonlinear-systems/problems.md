# P03.10  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: one Newton step by hand)

F(x, y) = (x² + y² − (1 + s/10), xy − 0.3) = 0.

(a) Write the Jacobian. Starting from (x, y) = (1.5, 0.5), compute F, J, the Newton step and the new iterate by hand, and the residual norm before and after.
(b) Run `newton` from the same start to 10⁻¹⁴: report the iterates, the residual norms and the observed order from `convergence_order`.
(c) Solve the system exactly (hint: (x + y)² and (x − y)²) and compare. Why does Newton fail from (1, 1)?

Answer form: a matrix, two vectors, two norms, an iteration table, two exact roots, one sentence.

## 2. (Derivation: convergence, damping, Broyden, and the arctan threshold)

(a) For a scalar f with simple root x*, show from the Taylor expansion of f(x*) about x_k that e_{k+1} = −f''(ξ) e_k²/(2f'(x_k)) for some ξ between x_k and x*, and hence that |e_{k+1}| ≤ C e_k² with C = max|f''|/(2 min|f'|) near x*. What does the formula give when f'(x*) = 0 (a double root)?
(b) Show that ∇(½‖F‖²) = JᵀF and that the Newton direction δ = −J⁻¹F satisfies ∇m·δ = −‖F‖². Then show that for small λ, m(x + λδ) = (1 − 2λ) m(x) + O(λ²), so the Armijo condition (3) with c₁ < 1 can always be met.
(c) Show that among all matrices B satisfying Bs = y, the one minimising ‖B − B_k‖_F is B_k + (y − B_k s)sᵀ/(sᵀs). (Hint: write B = B_k + E; Es = y − B_k s; show the proposed E has the smallest Frobenius norm among all E with that action on s, using ‖E‖_F ≥ ‖Es‖/‖s‖.)
(d) For f(x) = arctan x, show that the Newton iterate from x₀ lands at −x₀ exactly when 2x₀ = (1 + x₀²) arctan x₀, solve this with your `newton` (one unknown), and argue why |x₀| above this threshold diverges and below it converges.

Answer form: four derivations; one number to five digits.

## 3. (Unfamiliar scenario: a two-node radiator with T⁴ nonlinearity)

An electronics box (node 1) dissipating Q₁ = (20 + s) W is bolted to a radiator panel (node 2) through a conductance G = 0.5 W/K. Both radiate to deep space (0 K) with emissivity ε = 0.85 from areas A₁ = 0.05 m² and A₂ = 0.20 m²; σ = 5.670 × 10⁻⁸ W/(m²·K⁴).

(a) Write the two energy balances F(T₁, T₂) = 0 and the analytic Jacobian. Check the Jacobian against `fd_jacobian` at (300, 300) K.
(b) Solve with `newton` from (300, 300) K; report the temperatures, the residual norms per iteration and the observed order. Verify the energy balance: total radiated power equals Q₁.
(c) Run `newton` from (1, 1) K. What happens in the first step and why? How many iterations does it take, and how does `damped_newton` behave from the same start?
(d) Solve with `broyden` from (300, 300) and compare evaluations with Newton. Then compute the temperature of a single isothermal node with area A₁ + A₂ and explain why T₁ is above it and T₂ below it.

Answer form: two equations, a 2 × 2 Jacobian, two temperatures with units, two iteration counts, one energy check, three sentences.

## 4. (Numerical, uses your code: a more stable aircraft and its ceiling)

Set C_mα = −1.2 − s/50 rad⁻¹ (a more forward centre of gravity), everything else as in the lesson.

(a) Trim at V = 200 m/s and h = 0, 5 and 10 km; tabulate α, δe, T and the Newton iteration count from the physical initial guess. Compare δe at 5 km with the baseline aircraft's 2.218° and explain the sign of the change.
(b) Thrust available is T_avail = T₀ ρ(h)/ρ₀ with T₀ = (70 + s) kN. Tabulate T_required and T_avail at h = 0, 5, 8, 10 km. Where do they cross?
(c) Find the absolute ceiling at 200 m/s as a root-finding problem in four unknowns (α, δe, T/W, h) by appending F₄ = (T − T_avail(h))/W to the trim residual and using `newton` with `fd_jacobian`, starting from the 5 km trim and h = 5 km. Report h, α, δe, T and the iteration count. Why is a differenced Jacobian the natural choice here?
(d) The ISA model in the code is valid to 11 km. What would you check before trusting a ceiling above it, and how would the code tell you?

Answer form: a 3-row table, a 4-row table, four numbers with units and a count, two sentences.

## 5. (Numerical, extends your code: the chord method)

The chord method is Newton with the Jacobian computed once at x₀ and reused; the "m-step chord" recomputes it every m iterations. Implement `chord(f, x0, m, tol, max_iter)` with `fd_jacobian`, counting residual evaluations. Apply it to the plate `HotPlate(n=8, beta=-0.004 (1 + s/30))` to 10⁻⁹.

(a) Run m = 1 (Newton), 2, 3, 5 and m = ∞ (pure chord). Tabulate iterations, residual evaluations and the observed order.
(b) Which m minimises the number of residual evaluations? Explain the trade-off in two sentences.
(c) Show that the pure chord method converges linearly with a factor equal to the spectral radius of I − J(x₀)⁻¹J(x*), and compare with the measured error ratios.
(d) Compare with `picard_plate` on the same problem: iterations and cost per iteration. Which would you use at n = 64, and why?

Answer form: a 5-row table, one m with two sentences, one derivation with a comparison, two sentences.

## 6. (Judgment: choosing the nonlinear solver)

For each case say which method you would use (Newton with analytic Jacobian, Newton with differenced Jacobian, damped Newton, Broyden, Picard, Newton–Krylov, continuation) and why, in one sentence each: (a) the 3-unknown trim of this lesson, solved once; (b) the same trim inside a flight simulator at 100 Hz, with the previous frame's state available; (c) trim of a CFD model where one residual evaluation costs 10 minutes and no Jacobian is available; (d) a 10⁶-unknown nonlinear heat-conduction finite-element model; (e) a structural problem with a snap-through (the load–deflection curve folds back); (f) a problem where Newton from the obvious start sends the iterate to negative temperatures; (g) a mildly nonlinear coefficient problem A(u)u = b where assembling A is cheap and the Jacobian of A(u)u is not available.

Answer form: seven labelled lines.
