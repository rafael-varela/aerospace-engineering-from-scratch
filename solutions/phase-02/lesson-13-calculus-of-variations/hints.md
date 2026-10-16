# P02.13  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. (θ − sin θ)/(1 − cos θ) increases monotonically from 0 (at θ = 0) to ∞ (at θ = 2π) and equals π/2 at θ = π, so you know at once which side of π your root is on. k follows from y_f alone once θ_f is known.
2. The boundary term of the integration by parts is where η(a) = η(b) = 0 is used. For Beltrami, differentiate F − y′F_{y′} with the chain rule and watch the terms cancel against the Euler–Lagrange equation; an explicit x adds F_x, which does not cancel. For Snell, dx/ds = 1/√(1 + y′²) and v = 2√(gk) sin(θ/2) on the cycloid.
3. Snell: sin θ/c is the same at the ground and at the apex (where θ = 90°). A circle of radius R centred at height −c₀/a has, at height y, sin θ = (y + c₀/a)/R: proportional to c(y). The straight path has θ = 90° everywhere but c changes along it only if y changes, which it does not, so check what sin θ/c would need to be at a point just above the ground.
4. The Ritz bound is a minimum over a subset of admissible paths. The cycloid's slope at x = 0 is infinite; a finite sum of sines has a finite slope there. A term proportional to x^{2/3}(1 − x/x_f) has infinite slope at 0 and vanishes at both ends.
5. Exponent 2 means the first variation is zero (an extremal); exponent 1 means it is not. Dip criterion: θ_f > π, which by the monotonicity in hint 1 is x_f/y_f > π/2.
6. Friction does work μN ds along the whole track; the normal force vanishes on a vertical wall. A measured time ratio is independent of the ball's moment of inertia (both tracks see the same rescaled g), so the ratio is the clean test.
