# P02.04  Hints

One hint per problem. Read only the one you need.

1. The sign of Δ answers "saddle or not"; the sign of τ answers "stable or not"; the sign of τ² − 4Δ answers "oscillatory or not". For the sense of rotation, put a trajectory on the positive x₁ axis and look at the sign of ẋ₂ = c x₁ there. Nullcline slopes are −a/b and −c/d.

2. (a) Write out det(A − λI) = (a − λ)(d − λ) − bc and collect powers of λ. (b) Substitute the complex solution into ẋ = Ax and separate real and imaginary parts; the equation is linear and real, so each part satisfies it alone. (c) Differentiate Q using ẋ₁ = a x₁ + b x₂ and ẋ₂ = c x₁ + d x₂ and collect; the x₁x₂ terms cancel on their own and the x₂² term carries a factor (a + d). (d) The two statements are the two rows of the matrix; then read ω_n² = Δ and 2ζω_n = −τ.

3. Add the two equations. Total mass is conserved, so one eigenvalue is zero and the other is −(k₁ + k₂). The line of fixed points is where ṁ₁ = 0, i.e. k₁m₁ = k₂m₂.

4. For s = 15 the matrix is the lesson's own, so your numbers should match §2.3; for other seeds only M_q changes and Δ = 0.96·(M_q/−1.2) − M_α … no: Δ = (Z_α/V₀)M_q − M_α = −0.8 M_q + 4. Use a fine time grid (10⁴ points per second) to find the peak and the 0.5° time; "stays below" means the *last* crossing.

5. Think about the ratio of time scales (a factor 30 or more means the fast mode finishes before the slow one has moved). An unstable mode with a 1.6 s doubling time is slower than a pilot's reaction, which is exactly why unaugmented relaxed-stability fighters are flyable for seconds and airliners are not certified that way. For the damping discrepancy, consider effects that add or remove pitch damping (downwash lag, flexible fuselage, Mach effects, the dα/dt derivative M_α̇).

6. With τ fixed at −2 the conditions τ² = 4Δ and Δ = 0 are single equations in M_α. For (c), write τ and Δ in terms of m = M_q + k and set τ² = 4ζ²Δ with ζ = 0.5; the quadratic in m has a negative discriminant. Then minimise ζ²(m) = τ²/(4Δ) over the allowed m to see how far from a spiral you stay. For (d), set Δ = τ²/(4ζ²) with τ = −2 and solve for M_α.
