# P03.17  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. Δv₀ = g₀ I_sp ln MR. 2. Relative sensitivities: ∂ln Δv/∂ln I_sp = 1, ∂ln Δv/∂ln MR = 1/ln MR; u_rel² = (0.01)² + (0.02/ln MR)². 3. Shares: each term over u_rel². 4. z = (3700 − Δv₀)/(Δv₀ u_rel); P = Φ(z). 5. Halve the larger term's standard deviation and recompute u_rel.

## 2.
(a) E[μ̂] = (1/N)ΣE[y_n] = μ; Var μ̂ = (1/N²)ΣVar y_n = σ²/N by independence. Indicator: mean p, variance p(1 − p). Relative error √(p(1 − p)/N)/p = ε ⇒ N = (1 − p)/(pε²). (b) G(s)G(t) = e^{ξ(s+t) − (s² + t²)/2}; E[·] = e^{(s+t)²/2 − (s²+t²)/2} = e^{st} = Σ_n sⁿtⁿ/n!; the left side is Σ_{m,n} E[He_m He_n] sᵐtⁿ/(m! n!); compare. (c) e^{aξ} = e^{a²/2} e^{aξ − a²/2} = e^{a²/2} Σ He_n(ξ) aⁿ/n!. Then Σ c_n² n! = e^{a²} Σ a^{2n}/n! = e^{2a²} = E[e^{2aξ}]. (d) E[f] = Σ c_α E[He_α] = c₀; E[f²] = Σ c_α c_β E[He_α He_β] = Σ c_α² α!; subtract c₀². Additive: every α ≠ 0 has one nonzero entry, so Σ_i S_i = Σ_{α ≠ 0} c_α² α!/Var = 1.

## 3.
(a) s₀ = V²/(2a); relative sensitivities 2s_V and −s_a; u_rel² = (2s_V)² + s_a²; shares; second derivatives ∂²s/∂ξ_V² = 2s₀s_V², ∂²s/∂ξ_a² = 2s₀s_a² (from 1/(1 + s_aξ) ≈ 1 − s_aξ + s_a²ξ²); mean shift = s₀(s_V² + s_a²). (b)–(d) Code calls as in the lesson's `main` with d = 2; skewness with `scipy.stats.skew`; Sobol from `sobol_indices`.

## 4.
(a)–(c) `linear_propagation(f, 3)`, `pce_projection(f, 3, 3, 4)`, `pce_moments`, `sobol_indices`, `monte_carlo(f, 3, N, 17)`, `mc_probability(samples, 1.4e6)`; surrogate P from `pce_evaluate` on 4 × 10⁵ seeded normals; Gaussian P = Φ((1400 − mean)/std). N = (1 − p)/(0.01 p); PCE runs = 4³. (d) Mean ± std (k = 1), P with its error, dominant input with S_i, inputs assumed, runs used.

## 5.
(a) `monte_carlo(..., N, 17)["se_mean"]` and `antithetic_monte_carlo(..., N, 17)["se_mean"]`; equivalent N = N × (se_plain/se_anti)². (b) `bootstrap_ci(samples, lambda y: np.quantile(y, 0.05), 2000, 1)`; surrogate quantile with `np.quantile`; Gaussian mean + Φ⁻¹(0.05) std. (c) Compute the 5th percentile of the pair averages and compare with the true 5th percentile; the pair average has the variance of the even part only.

## 6.
(a) Linear (13 runs) for sanity; PCE by regression (degree 2 in d = 6 has 28 terms, ~60–100 runs) or sparse grid; Monte Carlo on the surrogate (10⁶ samples); direct runs concentrated in the tail region identified by the surrogate. (b) Correlation: transform to independent variables by a Cholesky factor of the covariance before building the PCE, and report Sobol indices only for the transformed variables; a bounded input: use a bounded distribution (and its orthogonal family, e.g. Legendre/Jacobi) rather than a truncated Gaussian pushed through Hermite polynomials. (c) Direct runs in the tail are the model; the surrogate is a polynomial extrapolating into a region where few quadrature points lie; raise the degree or fit locally and re-check. (d) State the number, its statistical error, the conditionality on elicited distributions, the independence and smoothness assumptions, and the measurements (flight data, arc-jet tests) that would move it.
