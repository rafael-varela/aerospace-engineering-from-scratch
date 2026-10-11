# P00.15  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. t ∝ h^½ g^−½, so each relative uncertainty is halved and they add in quadrature. "Contribute equally" means the two squared terms are equal.
2. Solve the quadratic in t and keep the positive root. For the derivative, differentiate under the square root and set v₀ = 0 before simplifying. For (c), write both shifts as functions of h and look at the exponents.
3. Everything on Mars follows from the same formulas with g and ρ changed; the terminal speed has ρ in the denominator under a square root. For (d), h = ½gt² and then look at the envelope table in `P00.13`.
4. The frame term enters through `u_single`; it is the part of u_meas that changes with fps. A constant bias moves every mean by the same amount but leaves every standard deviation unchanged.
5. v_t = √(2mg/(ρC_D A)) with A = πd²/4. The first-order correction is t₀·gh/(6v_t²); compare it with the closed form at the largest height. For the gap, subtract the two balls' drag-corrected times.
6. E_n is a ratio: a large denominator makes any numerator look small. A growing-with-h pattern and a shrinking-with-h pattern point at different assumptions; a constant pattern points at the release or the clock.
7. The bottom of the ball falls h_centre − r, not h_centre. The ratio of two square roots is the square root of the ratio, and a fixed error is a larger fraction of a small h.
