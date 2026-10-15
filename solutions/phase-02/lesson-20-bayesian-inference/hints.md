# P02.20  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The evidence is P(alarm) = P(alarm | stuck)P(stuck) + P(alarm | ok)P(ok). In odds form the Bayes factor is the ratio of the two alarm probabilities, and two independent alarms multiply it twice.
2. Every integral in this problem is ∫₀^∞ λ^{a−1}e^{−bλ}dλ = Γ(a)/b^a with a different a and b. For the Fisher information, differentiate log p twice in λ and take the expectation using E[k] = λT. For the χ² claim, compare the Gamma(k + 1, T) density of λ with the χ²_ν density of x = 2Tλ (a χ²_ν variable is Gamma(ν/2, 1/2)).
3. Multiply the two Gaussians and complete the square in b; the coefficient of b² gives the posterior precision and the coefficient of b gives the mean. For the grid, `scipy.stats.norm.logpdf(grid, mu0, s0)` is the log prior and `norm.logpdf(z, grid, sz)` the log likelihood as a function of the grid.
4. Everything is a function call; the thinking is in the scatter you attach to the chain's numbers (effective sample size) and in the one-sentence explanation of the gap in (d).
5. For (b), start from the Jeffreys prior and add 0 failures and 1.2 × 10⁶ h; the posterior is then Gamma(1.5, 1.26 × 10⁶). For (d), the probability of no failure in T hours under a Gamma(α, β) belief is (β/(β + T))^α.
6. The coverage experiment needs the true rate and the interval; when the prior is wrong the intervals are centred near the prior mean and the true rates are not. The diagnostic in (d) is the prior predictive probability of the observed count (or of a count at least as large).
