# P03.14  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. With Q = diag(1, κ) and x = (κ, 1) the gradient is (κ, κ): both components equal, so the step is the same in both coordinates and the new point is x₀ − α(κ, κ). Compute α from gᵀg/(gᵀQg) and factor out (κ − 1)/(κ + 1). Once x₁ is a multiple of (κ, −1), the same algebra repeats with a sign flip.
2. Expand ∇f(x_{k+1}) about x_k to first order. For positive definiteness write zᵀH₊z as a sum of two non-negative terms and ask when both can vanish. For (c), φ'(α) = ∇f(x + αp)ᵀp, so yᵀs = α(φ'(α) − φ'(0)). For (d), compute V = I − ρsyᵀ first; it is upper triangular here.
3. Truncation ∝ h², rounding ∝ ε|f|/h: the log–log plot has slopes +2 and −1 meeting near h ≈ (ε)^{1/3}. In metres, ∂f/∂h is 10⁴ times smaller than in tens of kilometres and ∂²f/∂h² is 10⁸ times smaller, which is what the condition number sees. A gradient tolerance is not scale-free; think about what is.
4. Above the tropopause fuel per km = c_T W/(a M (L/D)) with a and c_T constant, so the minimum over (M, C_L) does not involve W. Altitude comes from p = 2W/(γM²SC_L) and the isothermal p(h). The 1 % window follows from Δf = ½H_ii Δx_i².
5. Wrap f and grad in a counter. For the ratios use `error_sequence` and drop the rounding-floor entries; the order estimator divides two logarithms that are both near zero at the end. For (d), loop x ← x − H⁻¹g with `np.linalg.solve` and record f at every iterate.
6. A positive-definite Hessian and a zero gradient are local statements. The Hessian gives Δf ≈ ½H₁₁ΔM²; a direct evaluation includes the quartic wave-drag term the quadratic model does not. With noise δ in f the central-difference error is about δ|f|/h, which sets a much larger best h; think about what yᵀs looks like when the gradient error exceeds the gradient change over a step.
