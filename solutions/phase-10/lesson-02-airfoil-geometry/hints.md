# P10.02  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. At x = 0.25 < p = 0.4 you are on the forward parabola. The surface point is the camber point plus y_t along the *normal* to the camber line, (−sin θ, cos θ).
2. (a) Integrate the polynomial term by term; √x integrates to (2/3)x^{3/2}. (b) Keep only the √x term near x = 0 and square it. (c) Put x = 1 in the polynomial. (d) Differentiate the cubic and set x = p; the line is linear in k₁ because k₁ is an overall factor.
3. `area_between` and `surface_thickness_camber` do most of the work. Thickness here means vertical thickness at a station, because the battery is a flat slab.
4. Near the nose y_t ∝ √x, whose derivatives are infinite at 0; with x = (1 − cos β)/2 ≈ β²/4, √x ≈ β/2 is smooth in β.
5. A quoted pair (r, k₁) can be checked twice: does the camber peak at p, and does `naca5_ideal_cl(r, k1)` return 0.3?
6. Perpendicular and vertical thickness differ by a factor related to cos θ, and the surface points also move in x. Where is θ largest?
