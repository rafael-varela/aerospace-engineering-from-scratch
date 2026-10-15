# P02.18  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. With x = 2, 3, 4, 5 the mean is 3.5 and S_xx = 5 whatever s is; only S_xy and S_yy depend on your thrusts. The slope in N/MPa divided by 10⁶ is C_F A_t in m².
2. For (a) set ∂S/∂a = 0 and ∂S/∂b = 0; the first says the residuals sum to zero. For (c) write ŷ(x₀) = ȳ + b(x₀ − x̄) and note that ȳ and b are uncorrelated when x is centred. For (d), ∂/∂b of Γ(b + 1) is Γ(b + 1) ψ(b + 1) with ψ the digamma function (`scipy.special.digamma`).
3. A residual plot that bends (positive in the middle, negative at both ends, or the reverse) means the *model* is wrong, not a *point*. Fit only the range where the residuals look random.
4. Everything is a function call; for the comparison with the true value use `TRUE_THRUST_PARAMS`. LOO and Gauss–Newton answer different questions: "which polynomial predicts best" and "do these three physical parameters describe the data".
5. Leverage: a point far from x̄ pulls the line harder and its own residual is *smaller* than the fault that produced it, so it is harder to flag. For (d), P(|z| > 2) ≈ 0.046 per point.
6. (a) The band is for the *mean* response under the *linear* model; outside the data it is pure assumption. (b) Noise in x attenuates the slope by the factor S_xx/(S_xx + n σ_x²). (c) The band is pointwise: each x has its own 95 %.
