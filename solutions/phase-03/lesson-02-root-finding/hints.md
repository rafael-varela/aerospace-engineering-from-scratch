# P03.02  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Work in radians throughout; f(E) = E − e sin E − M and f′(E) = 1 − e cos E. The error after a step is roughly the *next* step, and for quadratic convergence the error after the second step is about (second step)² × f″/2f′.
2. (a) Expand f(x*) about x_k and keep the second-order remainder; divide by f′(x_k). (b) Write f = (x − x*)² h and f′ = 2(x − x*)h + (x − x*)² h′; the Newton step is f/f′ ≈ (x − x*)/2. (c) g(x_k) − g(x*) = g′(ξ)(x_k − x*). (d) The secant line is the linear interpolant through (x_{k−1}, f) and (x_k, f); its error at x* is f″(ξ)(x* − x_{k−1})(x* − x_k)/2; the secant root is where the interpolant vanishes.
3. f′(H) = e cosh H − 1 ≥ e − 1 > 0, so f is strictly increasing. For a bracket use f(0) = −M_h < 0 and, since sinh H ≥ H for H ≥ 0, f(H) ≥ (e − 1)H − M_h, which is positive at H = M_h/(e − 1).
4. Reference root from `brentq` with xtol 10⁻¹⁵; errors via `errors_from_iterates`; fit via `fitted_convergence_order` (needs at least four errors above the floor, which Newton may not supply: say so). Digits per evaluation: compare log₁₀ of the error reduction per evaluation.
5. `orbit_state(a, e, t)[0]` is the position; its norm minus R is g(t). Each g evaluation is one Kepler solve. For (b), cos E = (1 − R/a)/e. Radius increases monotonically from perigee to apogee and decreases after, so there are two crossings per orbit.
6. Count evaluations: Newton costs a sin and a cos per step; the secant method one sin per step after the first two. The cap and fallback follow from "must never return a wrong E": a bracketing method cannot. For (d), relate |f| to the error through f′ ∈ [1 − e, 1 + e].
7. The mean value theorem on sin gives |sin E₁ − sin E₂| ≤ |E₁ − E₂|. Uniqueness: two fixed points would have to satisfy |E₁ − E₂| ≤ e|E₁ − E₂|. Existence: the sequence is Cauchy (the errors are a geometric series). The bound counts every step at the worst rate e; the actual rate is e|cos E*|.
