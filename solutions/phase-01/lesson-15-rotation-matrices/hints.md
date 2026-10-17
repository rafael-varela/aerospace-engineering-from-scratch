# P01.15  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Build R₃ first (it only touches the top-left 2×2 block), then R₂, then R₁, and multiply from the right: R₂R₃ first, then R₁ times that. The (1,3) entry must come out as −sin θ; if it does not, you multiplied in the wrong order.
2. For (a), write (C₁C₂)ᵀ = C₂ᵀC₁ᵀ and use det(AB) = det A det B. For (c), set sin θ = 1, cos θ = 0 in the full matrix and use the angle-addition formulas on the products of sines and cosines.
3. "To *b* from *i*" is "to *b* from *t*" times "to *t* from *i*", and "to *b* from *t*" is the transpose of the mounting matrix you were given. The thrust direction in inertial axes is a column of C_b/i, not a row.
4. The code already has every function you need; the only new line is the exact step. Compute the misalignment yourself with `misalignment_angle` rather than reading angles back, because the total angle exceeds 180°.
5. In the new frame the axial component of **v** is unchanged, the component along the old perpendicular direction is scaled by cos *a*, and the component along **n** × **v** appears with a sign fixed by the passive convention; fix the sign by comparing with R₃(a) for **n** = e₃. For (c), differentiate arccos near 1.
6. For (a), compute δ = ω Δt first and ask whether "δ ≪ 1" holds. For (b), solve ωT δ²/3 = 0.1° for δ, then f = ω/δ. For (c), the matrix is fine everywhere; only the map from matrix to angles breaks.
