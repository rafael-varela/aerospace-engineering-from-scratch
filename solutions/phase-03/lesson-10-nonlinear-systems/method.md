# P03.10  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. c = 1 + s/10; F(1.5, 0.5) = (2.5 − c, 0.45); J = [[3, 1], [0.5, 1.5]]. 2. Solve Jδ = −F (det J = 4); x₁ = x₀ + δ; evaluate F(x₁). 3. `newton(f, [1.5, 0.5], j, tol=1e-14)`; `iterate_errors(res.xs, res.x)`; `convergence_order`. 4. (x + y)² = c + 0.6, (x − y)² = c − 0.6, so x = (√(c + 0.6) + √(c − 0.6))/2, y = (√(c + 0.6) − √(c − 0.6))/2 (and the three sign variants). det J = 2x² − 2y² = 0 on x = ±y: the Jacobian is singular at (1, 1).

## 2.
(a) Taylor f(x*) = 0 about x_k with Lagrange remainder; e_{k+1} = x_{k+1} − x* = e_k − f(x_k)/f'(x_k) = f''(ξ)e_k²/(2f'(x_k)) up to sign. Double root: f'(x_k) ≈ f''e_k, so e_{k+1} ≈ e_k/2: linear. (b) ∂/∂x_j of ½ΣF_i² = ΣF_i ∂F_i/∂x_j = (JᵀF)_j. m(x + λδ) = ½(F + λJδ)ᵀ(F + λJδ) + O(λ²) = ½(1 − λ)²‖F‖² + O(λ²) = (1 − 2λ)m + O(λ²); for λ small enough this is below (1 − 2c₁λ)m when c₁ < 1. (c) Any E with Es = y − B_k s =: r has ‖E‖_F ≥ ‖r‖/‖s‖; E* = rsᵀ/(sᵀs) has E*s = r and ‖E*‖_F = ‖r‖‖s‖/(sᵀs) = ‖r‖/‖s‖: it attains the bound. (d) x₁ = x₀ − (1 + x₀²) arctan x₀ = −x₀ ⟺ 2x₀ = (1 + x₀²) arctan x₀; `newton(lambda v: [2 v - (1 + v²) arctan v], [1.5])`. Beyond the threshold |x₁| > |x₀| and the sequence grows; below it |x₁| < |x₀| and the next step is inside too.

## 3.
(a) F as in the hint; J = [[−G − 4εσA₁T₁³, G], [G, −G − 4εσA₂T₂³]]. `fd_jacobian(F, [300, 300])`. (b) `newton(F, [300, 300], J, tol=1e-10)`; `iterate_errors`, `convergence_order` (divide the errors by 300 so they are below 1); check εσ(A₁T₁⁴ + A₂T₂⁴) = Q₁ and G(T₁ − T₂) = εσA₂T₂⁴. (c) J(1, 1) ≈ G[[−1, 1], [1, −1]] + O(10⁻⁷): det ≈ 10⁻⁸, the step is 10⁸–10⁹ K; from there Newton on a T⁴ function shrinks the iterate by a factor 3/4 per step (f/f' = T/4 for a pure quartic), so it needs ≈ log(10⁹/300)/log(4/3) ≈ 50 steps. `damped_newton` shortens the first step and then converges. (d) `broyden(F, [300, 300])`, compare `nfev`. Single node: T = (Q₁/(εσ(A₁ + A₂)))^{1/4}; the conductance makes node 1 hotter and node 2 cooler than the isothermal case.

## 4.
(a) `ac = Aircraft(Cma=-1.2 - s/50)`; `solve_trim(ac, 200, h)` for h in (0, 5e3, 1e4). A larger |C_mα| means more nose-down moment per degree of α, so the elevator has to do less (or push the other way) for the same α: compare δe signs. (b) `isa_density(h)`; T_avail = T₀ρ/ρ₀. (c) `F4(z) = append(trim_residual(z[:3], ac, 200, isa_density(z[3]*1e3)), (z[2] W − T₀ρ(z[3])/ρ₀)/W)`; `newton(F4, z0, None, tol=1e-12)` with `z0 = [*trim5km.x, 5.0]` (h in km so that it is O(1)). (d) `isa_density` is the troposphere formula; above 11 km it is wrong (the real atmosphere is isothermal there, so ρ falls faster). Check the converged h against 11 km; add an assertion or extend the model.

## 5.
(a) `chord(f, x0, m, tol, max_iter)`: recompute `J = fd_jacobian(f, x)` when `k % m == 0`; count 1 + 64 per Jacobian and 1 per residual. (b) Evaluations = (number of Jacobians) × 64 + iterations; fewer Jacobians means more, cheaper iterations. (c) x_{k+1} = g(x_k) with g(x) = x − J₀⁻¹F(x); g'(x*) = I − J₀⁻¹J(x*); linear convergence with factor ρ(g'(x*)); compute with `np.linalg.eigvals`. (d) `picard_plate(pl, tol=1e-9)`: each iteration assembles A and solves a 64 × 64 system; no Jacobian ever. At n = 64 the differenced Jacobian costs 4 097 assemblies per Newton step, so Picard or a Newton–Krylov method with an analytic Jacobian-vector product wins.

## 6.
(a) Newton with analytic J. (b) Newton (or even one Newton step per frame) from the previous frame: continuation. (c) Broyden, started from a differenced Jacobian only if you can afford n + 1 runs, otherwise from a crude model Jacobian. (d) Newton–Krylov (preconditioned) or Picard if the nonlinearity is mild. (e) Pseudo-arclength continuation; plain parameter continuation fails at the fold. (f) Damped Newton, or a better start (continuation from a problem you can solve), or a change of variables (solve for log T). (g) Picard: assembling A(u) is the whole cost, and the Jacobian is not needed.
