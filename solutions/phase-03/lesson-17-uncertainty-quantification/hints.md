# P03.17  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Work in relative terms: d ln Δv = d ln I_sp + d ln(ln MR), and d ln(ln MR) = d(ln MR)/ln MR = (dMR/MR)/ln MR. The relative sensitivities are 1 and 1/ln MR.
2. (a) Linearity of expectation and independence; the indicator 1[y < y_crit] is Bernoulli. (b) Multiply two generating functions with parameters s and t, take the expectation using E[e^{ξu}] = e^{u²/2}, and compare the coefficient of sᵐtⁿ on both sides. (c) Set t = a in the generating function. (d) Take the expectation of (5), then of its square, and use orthogonality; for an additive model every variance term belongs to exactly one input.
3. s_L = V²/(2a) is quadratic in V and inverse in a: relative sensitivities 2 and −1, second derivatives 2s₀ s_V² and 2s₀ s_a² in the standardised variables. Define `f = lambda xi: (V*(1 + 0.03*xi[:, 0]))**2 / (2*a*(1 + 0.08*xi[:, 1]))`.
4. `dataclasses.replace(RangeInputs(), LD_rel_std=..., R_required=1.4e6)` builds your case; everything else is the lesson's `main` with your parameters. For (c), N ≈ (1 − p)/(p ε²).
5. The pair average (f(ξ) + f(−ξ))/2 is the even part of f: its mean is the true mean, but its distribution is not the distribution of f, so quantiles of the pair averages are not quantiles of the range. For (b), a 95 % interval misses one time in twenty by construction.
6. Order the methods by cost and by what they can deliver: linear first (sanity), PCE regression or sparse grid for six inputs (moments, Sobol, surrogate), Monte Carlo on the surrogate for the tail, direct runs in the tail to check the surrogate. 10⁻⁴ at 90 % confidence needs, roughly, 2.3/10⁻⁴ surrogate samples for zero failures observed... think about what "confidence" means for a count.
