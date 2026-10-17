# P01.11  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The vortex speed is Γ/(2πr), always tangential, counter-clockwise for Γ > 0. Circulation counts only what is *inside* the loop; flux counts only sources inside. Equation (9) predicts the numerical error of the circle before you run it.
2. Differentiate θ = atan2(y, x) with the chain rule: ∂θ/∂x = −y/r², ∂θ/∂y = x/r². The fundamental theorem needs a potential that is single-valued *along the loop*; ask where θ jumps by 2π and whether the loop crosses that cut a net number of times.
3. Peak wind at the core edge is Γ/(2πa), so Γ = 2πaV_max. Inside the core the circulation scales with the enclosed area; outside it is constant. Vorticity is the circulation per unit area of the core.
4. The midpoint of a chord of a circle of radius R subtending angle 2π/n lies at radius R cos(π/n); the chord has length 2R sin(π/n) and is parallel to the velocity there. For (d), pair each segment with the one diametrically opposite and look at the doublet's contribution.
5. With u = 1 − r/R the integrand becomes u^{1/7}(1 − u), a polynomial in u with one fractional power. The midpoint rule is second order only when the integrand has a bounded second derivative; (1 − r/R)^{1/7} does not at r = R. A substitution r = R(1 − w⁷) makes the integrand smooth in w.
6. For each statement ask which of the lesson's assumptions (two-dimensional, inviscid, steady, point singularity, uniform profile, point-mass gravity) it rests on, and what instrument would show the assumption breaking.
7. Each vortex sits in the field of the other: the induced velocity at distance b′ is Γ/(2πb′). At the midpoint, *both* vortices induce velocity in the same direction, and the distance is b′/2. Circulation is additive over the vortices a loop encloses.
