# P03.11  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Everything scales from α = k/(ρc_p) and the two time scales L²/α and Δx²/(2α). The FTCS step count grows as n² because Δt_max ∝ Δx²; the implicit count does not depend on n at all. Before trusting a large-r Crank–Nicolson run, look at g(π).
2. Expand every term about the same point (x_j, t^m for FTCS, t^{m+1} for BTCS), keep terms to Δt² and Δx⁴, and replace u_tt by α²u_xxxx using the PDE itself. For the factors, substitute u_j^m = gᵐ e^{ijθ} and use e^{−iθ} − 2 + e^{iθ} = −4 sin²(θ/2). For (c), expand e^{−z} and each factor to z³.
3. The ghost node makes the boundary row u_n^{m+1} = u_n + 2r(u_{n−1} − u_n). Gershgorin: every eigenvalue of the update matrix lies within a disc centred on a diagonal entry with radius the sum of the row's off-diagonal magnitudes; the interior rows give [1 − 4r, 1] and the boundary row gives the same interval. The mixed series of `P02.11` is exact and needs only a few terms after one hour.
4. `time_to_limit_exact` bisects the series; `time_to_limit` returns (time, steps). The BTCS error is first order in Δt, so halving the step halves it; the flux comes from `cold_face_flux_exact` at the exact time.
5. For a first-order scheme e ≈ CΔt, so the Δt that gives a target error is Δt₀ × (target/e₀). The r = 1/6 cancellation needs constant α, a uniform grid and the heat equation itself (so that u_tt = α²u_xxxx).
6. r = αΔt/Δx² and g(π) = (1 − 2r)/(1 + 2r); the decay of the sawtooth after N steps is |g(π)|ᴺ. The remedy is two BTCS half-steps (Rannacher start-up), which multiply the sawtooth by (1/(1 + 2r))². For (d), the time to a threshold scales as 1/α when the boundary and initial conditions are fixed.
