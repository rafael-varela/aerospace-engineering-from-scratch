# P05.14  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Mohr's circle has its centre at (σ_x + σ_y)/2 and its radius √(((σ_x − σ_y)/2)² + τ_xy²). The out-of-plane principal stress of a plane state is zero, and it belongs in the sort before you take σ₁ − σ₃.
2. (a) Write ε in terms of σ, form σ : ε, and split σ into pI + s with tr(s) = 0; the cross term p tr(s) vanishes. (b) Put σ₂ = σ₃ + f(σ₁ − σ₃) with 0 ≤ f ≤ 1 and minimise 1 − f + f² over f. (c) Pure shear has σ₁ = τ, σ₃ = −τ.
3. `principal_stresses` takes the full tensor. A hydrostatic addition changes every principal stress by the same amount; ask which formulas depend on differences only.
4. Write a function of t that returns MS_yield(t) − 0.15 using `critical_point(M, T, P, 0.05, t)`, and bracket the root between 1 mm and 10 mm. For Tresca, evaluate `tresca(tube_stress_at(np.pi, ...))` at the compression side.
5. A compressive σ_x with shear still makes a small tensile σ₁. Coulomb–Mohr weights σ₁ by 1/S_t and σ₃ by 1/S_c.
6. Separate three questions: is the number right for what it claims to be; is it the right number to compute; and what decision does a zero margin support. Look up what "typical" means compared with an A- or B-basis allowable.
