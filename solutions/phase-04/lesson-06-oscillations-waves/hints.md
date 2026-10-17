# P04.06  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The frame term for one timing is 0.41/f; dividing the timing by N divides its uncertainty by N. Compare it with T₀ × (T/T₀ − 1) at 10°.
2. (a) Differentiate E and use the equation of motion. (b) Peaks are one damped period 2π/ω_d apart, and the envelope is e^{−ζω₀t}. (c) Set ω = ω₀(1 + ε) and keep first-order terms in ε and ζ. (d) Energy: ½L²θ̇² = gL(cos θ − cos θ₀); a quarter period is the time from θ₀ to 0. cos θ = 1 − 2 sin²(θ/2). (e) The vertical force on an element is the difference of T_s sin(slope) at its two ends.
3. Q = 1/2ζ. Solve X/X_static = 10 numerically on a fine frequency grid (or as a quadratic in ω²). The envelope is e^{−ζω₀t}.
4. Courant 1 is the special case where the scheme is exact; below it the scheme disperses but still conserves its discrete energy. For the convergence study, avoid ending at a whole number of periods.
5. The amplitude correction multiplies g by (T/T₀)² at θ₀ = 0.15 rad. Weighted mean: weights 1/u².
6. For each case ask: which term in the equation (restoring force, damping, mass distribution, support) is no longer what the model says?
