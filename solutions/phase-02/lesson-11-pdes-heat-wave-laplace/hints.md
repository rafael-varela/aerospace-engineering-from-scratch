# P02.11  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Everything is a function of Fo = αt/L². The one-term formula keeps only n = 0 of (4); the n = 1 term has the opposite sign at the back face.
2. For (b), apply X(0) = 0 first (kills the cosine), then X′(L) = 0 to what is left. For (c), ∫₀ᴸ sin λ_n x dx = (1 − cos λ_n L)/λ_n and cos λ_n L = 0. For (d), differentiate twice in each variable and add.
3. Neumann at both ends means X′(0) = X′(L) = 0: cosines, and n = 0 is allowed. A closed–open tube has one node and one antinode, so quarter-wavelengths fit, exactly as in the tile's eigenvalues (3).
4. The amplification factor g = 1 − 4r sin²(κΔx/2) is most negative at κΔx = π (the 2Δx zigzag). At r = 0.5, g = −1 there: no growth, no decay, so the zigzag never dies. Count steps as t/(rΔx²/α).
5. sin(nπa/L) vanishes when na/L is an integer. The series error sits at the kink (Gibbs); the leapfrog energy error comes from the central difference of u_x across the kink, which is a one-cell error.
6. Think about the *sign* each assumption pushes the answer: a fixed hot face from t = 0 adds heat early; an insulated back keeps it; constant cold-value α slows it. A ratio of two thermocouple readings removes α from the comparison.
