# P10.05  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. At each point compute β = atan2(η, ξ − L) − atan2(η, ξ), r₁, r₂. 2. u = β/2π, v = −ln(r₁/r₂)/2π (uniform clockwise sheet, local axes = global here). 3. Point vortex: u = Γη/(2πr²), v = −Γ(ξ − ½)/(2πr²). 4. `panel_geometry([0,1],[0,0])` and `linear_vortex_influence(px, py, g)`; add the a and b parts; compare the split at (0.5, 0.5) by symmetry about the midpoint.

## 2.
(a) Expand s = (s − ξ) + ξ; ∫(s − ξ)η/D ds = η · ½ ln D; ∫(s − ξ)(ξ − s)/D ds = −∫(1 − η²/D) ds. (b) u = [g_a(I₁ − I₃/L) + g_b I₃/L]/2π, v = −[g_a(I₂ − I₄/L) + g_b I₄/L]/2π; on the panel η → 0± gives β → ±π. (c) N + 1 unknowns, N tangency equations, one Kutta row. γ₀ and γ_N are the surface velocities along the node direction at the two trailing-edge nodes, whose directions are opposite. (d) Velocity jump = γ; inside zero ⇒ outside = ∓γ; Γ = ∮γ ds with the trapezoid rule on linear γ.

## 3.
1. `th = linspace(0, 2π, N+1)`, x = cos θ, y = sin θ. 2. `vortex_panel(x, y, α)` for α = 0, 5°; compare Γ with 4π sin α. 3. Relate to `P10.03` (2): u_θ(0) = 0 ⇒ Γ = 4πVa sin α.

## 4.
1. μ, τ from your seed; `kt_convergence(1.0, mu, tau, radians(3), [20, …, 320])`. 2. Orders by `observed_order`. 3. Relative error; with the observed order p, N needed ≈ N_k (e_k/target)^{1/p}. 4. Argue from the asymptotic range and the trailing-edge corner.

## 5.
1. `naca4(code, 80)`; `vortex_panel` at 0 and 4°; slope = Δc_l/Δα; α_L0 = −c_l(0)/slope; min of `cp` and its `xc`. 2. `thin_airfoil_naca4`. 3. Interpret c_m,c/4 ≠ 0 for a symmetric section via the a.c. position.

## 6.
(a) Recall the corner-flow exponent at a finite-angle edge; the error at the nearest panel is O(1) for every N. (b) Stop the boundary layer a little before the trailing edge (the lesson uses x/c = 0.99); estimate the effect with `section_drag`'s panel-count study. (c) Collect: N = 160 error ~0.01 % (§4), corner affects only local Cp, boundary layer ~2–3 % of c_l (`P10.06` displacement iteration) and ~12–17 % of the slope (`P10.03` Problem 6), finite wing: the elliptic-wing lift slope 2π AR/(AR + 2) is 19 % below 2π at the trainer's AR = 8.33 (`P10.09`, `P10.10`).
