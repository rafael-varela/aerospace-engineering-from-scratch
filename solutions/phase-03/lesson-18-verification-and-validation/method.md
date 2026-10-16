# P03.18  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. d₂₁ = f₂ − f₁, d₃₂ = f₃ − f₂; p = ln(d₃₂/d₂₁)/ln 2. 2. f_ex = f₁ − d₂₁/(2ᵖ − 1). 3. ε₂₁ = d₂₁/f₁; GCI = 1.25|ε₂₁|/(2ᵖ − 1). 4. Repeat 2–3 with p = 2; ratio = (ε₃₂/ε₂₁)/4 with ε₃₂ = d₃₂/f₂. 5. Factor = (GCI/0.001)^{1/p}.

## 2.
(a) u^{n+1} − uⁿ = Δt u_t + (Δt³/24)u_ttt + …; θu^{n+1} + (1 − θ)uⁿ = u + (θ − ½)Δt u_t + (Δt²/8)u_tt + …; D₂u = u_xx + (Δx²/12)u_xxxx + …. Insert into (1), subtract the PDE at (x_i, t_{n+½}). (b) Compute D₂[2t x(L − x)] = −4t exactly; the time difference of 2t x(L − x) over Δt is 2x(L − x) exactly; the CN average of α(−4t) + q = 2x(L − x) is exact because it is linear in t; the two sides agree identically; ends are zero. (c) f₂ − f₁ = Chᵖ(rᵖ − 1), f₃ − f₂ = Chᵖrᵖ(rᵖ − 1); divide; solve the first for Chᵖ and insert in f = f₁ − Chᵖ. GCI₃₂/GCI₂₁ = |ε₃₂/ε₂₁| = rᵖ f₁/f₂. (d) The θ-average of q at n, n + 1 equals q(t_{n+½}) + O(Δt²); q(t_n) = q(t_{n+½}) − (Δt/2)q_t + O(Δt²); the difference is a first-order truncation term; with Δt ∝ Δx the error is C₁Δx² + C₂Δx, order 1 as Δx → 0.

## 3.
(a) q_m(t) = −(T(0) − T_a)e^{−t/τ}/τ + k[(T_a + (T(0) − T_a)e^{−t/τ})⁴ − T_a⁴]; at t = 0: −(T(0) − T_a)/τ + k(T(0)⁴ − T_a⁴). (b) `rhs = lambda T, t: -k*(T**4 - Ta**4) + qm(t)`; Heun loop; error at 1000 s; `observed_order`. (c) Replace `t + h` by `t` in k₂; the increment ½h(k₁ + k₂) is wrong by ½h·hq_m' = O(h²) per step, O(h) globally. (d) The manufactured source is additive and exercises the same code path; the only thing not exercised is a bug that depends on q being zero (e.g. a source multiplied by the wrong factor would be caught; a wrong k would be caught because T_m⁴ enters through k).

## 4.
(a) `u_expr = (20+s)*sp.sin(sp.pi*x/L)*(1-sp.exp(-t/tau)) + s*x*(L-x)/L**2*sp.sin(2*sp.pi*t/tau)`; `manufactured_source(u_expr, x, t, alpha, beta)`; hand: differentiate term by term, q = u_t − αu_xx + βu at (0.1, 50). (b) `convergence_study(run, [25,50,100,200,400])` with `heat_theta(..., beta=beta, source=q_m)`; `observed_order`. (c) `richardson(f1, f2, f3)`, `required_refinement(gci, p, 1e-4)`, compare `f_extrap` with the 800-cell value.

## 5.
(a) `synthetic_measurements(P, seed=s)`; model on 200 cells, 600 steps; `validation_comparison(model(0), data, u_data, 0.012)`. (b) Hold-out fit; `validation_comparison(model(beta_ho), ...)["z"][-1]`; full fit; curvature by three χ² evaluations; z_truth = (β_true − β̂)/σ_β. (c) h = β̂ρcA/P; `heat_theta(..., 1200.0, 1200, ...)[-1][80]` for the three β's; report the β̂ prediction with ± (∂u/∂β)σ_β, and state the synthetic origin and the single-data-set calibration.

## 6.
(a) Verification ≠ validation ≠ "agrees with a correlation"; ask for: an order-of-accuracy study (manufactured or exact solution) for the discretisation, a grid-convergence study on a geometry with a bend, and the solver's iterative-convergence criterion and its effect on the result. (b) Pattern: "The simulation shall [quantity] [criterion] [condition], verified by [method]" — one for the observed order, one for the GCI of the reported coefficient, one for the comparison error against a referent with stated uncertainty. (c) p = ln((f₃ − f₂)/(f₂ − f₁))/ln 1.5; GCI = 1.25|ε₂₁|/(1.5ᵖ − 1); compare GCI with 10 % and 2 %. (d) The comparison error is 4 % with a validation uncertainty of √(8² + GCI²) %; a 4 % agreement at 8 % uncertainty is consistent with model errors up to about 12 %; say what further referent would tighten it.
