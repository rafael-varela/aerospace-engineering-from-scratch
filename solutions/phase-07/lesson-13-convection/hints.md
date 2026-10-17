# P07.13  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The film temperature is the mean of the surface and the free stream; at Re of order 10⁴ the high-Re bracket of Churchill–Bernstein is no longer close to one, which is why the two forms differ.
2. h_x ∝ x^(−1/2) integrates to 2x^(1/2); for the turbulent law, h_x ∝ x^(−0.2) integrates to x^0.8/0.8. Continuity of the mixed average needs only that the turbulent piece vanish when Re_L = Re_c.
3. Write ν = μ/ρ and α = k/(ρc_p) with ρ = p/(RT); Ra contains 1/(να). Radiation does not care about the pressure at all.
4. The wall limit is met only on the turbulent branch; scan the flow upward and bracket the first sign change before calling `brentq`. U = ṁ/(ρπD²/4).
5. In the average form the analogy reads St̄ Pr^(2/3) = C_F/2, with St̄ = h̄/(ρc_pU); the drag acts on both faces, so the reference area is twice the plate area.
6. For each step ask: is the correlation inside its stated range, is it the right face or regime, and is the omitted effect small compared with the correlation's own scatter?
