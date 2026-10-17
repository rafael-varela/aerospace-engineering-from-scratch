# P03.01  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Find the binade first: 2²² = 4.19 × 10⁶ ≤ r < 2²³ = 8.39 × 10⁶, so ⌊log₂r⌋ = 22 for every seed. Digits held = −log₁₀(ulp/r). The requirement is on the *separation*, the error is the ulp of the *absolute position*.
2. For (a) write fl(x̂ − ŷ) = (x̂ − ŷ)(1 + δ₃) and expand. For (b) the absolute error of cos θ is at most u × 1; the result has size θ²/2. For (c) both terms of E[x²] − m² are of size m² and carry absolute error ~ε m². For (d) write E(h) = h F₂/2 + 2ε F₀/h, set dE/dh = 0.
3. The ulp of a float32 at x is 2^(⌊log₂x⌋ − 23); 604 800 lies between 2¹⁹ and 2²⁰. Multiply a time spacing by c to get a distance. The representation fix separates the large integer part from the small fractional part.
4. With equal summands, once the partial sum S is large compared with 0.1 each addition rounds to the grid at S, whose spacing grows with S, and always in the same direction. The error per step is proportional to S ∝ n, so the total is proportional to n². Separate "the stored 0.1 is wrong" from "the additions round".
5. `numpy.spacing(numpy.float32(4e5))` is 0.03125 m, larger than σ. A number that cannot be stored cannot be averaged. Shifting the data moves them into a binade where the spacing is far below σ.
6. For each item compute ulp(magnitude) in binary32 and compare with the requirement; where the ulp is too large, ask whether the quantity can be stored relative to a reference (representation) or whether the loss is in a subtraction that a different formula avoids.
7. For Sterbenz, take x ≥ y without loss of generality; the hypothesis gives 0 ≤ x − y ≤ y ≤ x, so the exponent of x − y is at most the exponent of y. Both x and y are integer multiples of ulp(y), so x − y is too, and a multiple of ulp(y) that is at most y fits in 53 bits.
