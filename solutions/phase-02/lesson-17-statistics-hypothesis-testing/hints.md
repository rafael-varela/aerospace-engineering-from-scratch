# P02.17  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Compute the deviations from the mean first and square them; divide by 5, not 6. The one-sided p-value is the upper tail of t₅ at your statistic; the two-sided one is twice it.
2. For (a), write xᵢ − x̄ = (xᵢ − μ) − (x̄ − μ), expand the sum of squares, and use E[(x̄ − μ)²] = σ²/n. For (c), both "p < α" and "lower bound > μ₀" reduce to the same inequality on t.
3. The two tests ask different questions of the same data; with n = 6 the t test rests on normality it cannot check, the permutation test on exchangeability only. Part (d) is about individual motors, not means: a tolerance bound.
4. Everything is a function call; the standard error of a simulated rate is √(p(1 − p)/N). For (d) write a ten-line loop with `np.random.default_rng(s)`.
5. `one_sided_tolerance_factor(12)` and `lower_tolerance_bound(batch)`. For (c) copy the body of `tolerance_bound_coverage_mc` and replace `k`. For (d) try n = 100, 300, 1000, 3000 and interpolate.
6. A p-value conditions on H0; the colleague's sentence conditions on the data. Repeated testing until significance is "optional stopping"; the fix is to fix n in advance from the power calculation. For (d), the rule involves n p₀.
