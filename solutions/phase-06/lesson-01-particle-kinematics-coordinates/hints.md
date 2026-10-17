# P06.01  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Put the numbers straight into (2) and (3); for the path description, the sign of a_t says faster or slower, and the sign of the vertical component of v × a says left (counterclockwise, positive) or right.
2. (b) Write a_θ = (1/r) d(r²θ̇)/dt. For the maxima, put the aircraft on the line y = d, x = Vt, so θ = atan(Vt/d), and differentiate. (c) For the helix, |v × a| and |v| are both constant; compute them once.
3. Elevation is atan(z/D). Maximise ė(t) by setting its derivative to zero; the condition becomes a relation between D and z. The rocket's acceleration is vertical, so its components along e_R and e_e are projections of one vertical vector.
4. The best step moves with the fourth root of the noise and the best error with its square root; compute the per-axis noise as √((σ_r² + (σ_θ r̄)²)/2) with r̄ the mean range.
5. Look up the error at h = 0.5 s in the lesson's noise table and express it in g. A bias that drifts with time adds to θ̇ itself; find which acceleration terms contain θ̇ multiplied by something large.
