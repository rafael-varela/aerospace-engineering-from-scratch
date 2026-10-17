# P03.10  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. J = [[2x, 2y], [y, x]]; solve the 2 × 2 system by elimination. The exact roots come from x + y = √(c + 0.6) and x − y = √(c − 0.6). At (1, 1) look at det J.
2. (a) 0 = f(x*) = f(x_k) + f'(x_k)(x* − x_k) + ½f''(ξ)(x* − x_k)²; divide by f'(x_k) and use x_{k+1} = x_k − f/f'. (b) Chain rule on ½FᵀF; then m(x + λδ) = ½‖F + λJδ + O(λ²)‖² with Jδ = −F. (c) ‖E‖_F ≥ ‖Es‖₂/‖s‖₂ for any matrix (Cauchy–Schwarz on each row), with equality for a rank-one E = (Es)sᵀ/(sᵀs). (d) The Newton step is x − (1 + x²) arctan x; set it equal to −x.
3. F₁ = Q₁ − G(T₁ − T₂) − εσA₁T₁⁴, F₂ = G(T₁ − T₂) − εσA₂T₂⁴; ∂(T⁴)/∂T = 4T³. At (1, 1) K the radiative terms and their derivatives are negligible, so J ≈ G[[−1, 1], [1, −1]] is nearly singular and the step is enormous. Sum the two balances to get the energy check.
4. `Aircraft(Cma=-1.2 - s/50)`; `solve_trim`. For the ceiling, build `F4(z)` that calls `trim_residual` with `isa_density(z[3] * 1e3)` and appends the thrust condition, then `newton(F4, z0, None)`. ρ(h) is given by a formula whose derivative you would rather not write; that is what `fd_jacobian` is for.
5. Copy `newton`, compute `J = fd_jacobian(f, x)` only when `k % m == 0`. The chord iteration is x ← x − J₀⁻¹F(x), a fixed-point iteration whose derivative at x* is I − J₀⁻¹J(x*).
6. Count: cost of F, availability of J, distance of the start from the root, size n, and whether the solution curve folds.
