# P03.18  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Form the two differences f₂ − f₁ and f₃ − f₂ first; everything (p, the extrapolation, the GCIs, the ratio) is built from them and from r. The asymptotic ratio with the formal order is (ε₃₂/ε₂₁)/r².
2. (a) Expand each of u^{n+1}, uⁿ about t_{n+½}; odd powers of Δt/2 cancel in the difference and in the average. (b) D₂ of a quadratic in x is exactly its second derivative; the average of a quantity linear in t over n and n + 1 is exactly its value at n + ½. (c) Subtract the three equations in pairs; the ratio of the two differences is rᵖ. (d) q(t_n) = q(t_{n+½}) − (Δt/2)q_t + …
3. q_m = T_m' + k(T_m⁴ − T_a⁴) with T_m' = −(T(0) − T_a)e^{−t/τ}/τ. The bug makes the second stage evaluate q_m one step early: a stage error of hq_m', which is a first-order perturbation of the increment.
4. `sp.symbols("x t")`, build the expression with `sp.sin`, `sp.exp`, `sp.pi`; pass `beta` to `manufactured_source` *and* to `heat_theta`. For (c) the heater node is `int(round(nx * 0.4))`.
5. `fit_loss_coefficient(lambda b: model(b)[:5], data[:5], u_data[:5])` fits on five rows. The curvature of χ² at the minimum: χ²'' ≈ (χ²(β̂ + δ) − 2χ²(β̂) + χ²(β̂ − δ))/δ² with δ = 2 × 10⁻⁵ s⁻¹. h = β ρ c A/P with A = 10⁻⁴ m², P = 0.04 m.
6. A plot against a correlation on a different geometry at a single grid is neither code verification (no order, no exact solution) nor validation for *your* intended use. For (c), p = ln((f₃ − f₂)/(f₂ − f₁))/ln 1.5 and GCI = 1.25|ε₂₁|/(rᵖ − 1); compare the GCI with the required accuracy *after* subtracting what the validation uncertainty will take.
