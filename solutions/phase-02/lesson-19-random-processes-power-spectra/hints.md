# P02.19  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Write Φ_w in terms of x = ω/b and differentiate the dimensionless part (1 + 3x²)/(1 + x²)² with the quotient rule; the numerator of the derivative factors. Everything else is substitution; the first zero of R_w is where 1 − aτ/2 = 0.
2. For (a) the substitution x = tan θ turns dx/(1 + x²)² into cos²θ dθ. For (b) differentiate the known transform of e^{−a|τ|} with respect to a to get the transform of |τ|e^{−a|τ|}. For (c) expand (1 − aτ/2)² e^{−2aτ} and integrate term by term; Var(σ̂²) = (2/T)∫_{−∞}^{∞}R² dτ. For (d) look at ω²Φ_w as ω → ∞.
3. An integrator has |H|² = 1/ω², and 1/ω² is not integrable at ω = 0; but at a fixed time t the integral ∫₀^t n dt′ is a sum of independent increments, and the variance of a sum of independent pieces adds. The one-sided G₀ is twice the two-sided S₀ by (3).
4. Everything is a function call; the only decision is what scatter to attach. Use (13) for the gust variance and remember that the RMS has half the relative scatter of the variance. For the peak frequency, write |H_a|²Φ_w at high ω/b and low ω/c and find where the two power laws cross.
5. Ask of each quantity: is the integral that defines it dominated by frequencies near the break frequency b, near the aircraft's c = 1/τ, or by the tail beyond both? Quantities that converge only because of a filter's roll-off are the ones to distrust.
6. The variance accuracy fixes T; the resolution fixes the segment length; K is the ratio of the two (times two for 50 % overlap, minus one). The three numbers do not all come out as you want; that is the point.
