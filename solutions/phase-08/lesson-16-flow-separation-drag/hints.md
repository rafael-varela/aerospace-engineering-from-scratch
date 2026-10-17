# P08.16  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The drag area of a cylinder in cross-flow is its diameter times its length. For (c), the shedding frequency is St V/d with St itself a function of Re: solve f(V) = 90 Hz with a root finder.
2. (a) Write the outward normal of the surface element at angle θ and take the component of −p n along the stream. (b) Split the integral at ±θ_s and use ∫cos θ sin²θ dθ = sin³θ/3. (c) Divide the force balance by μ²/(ρ_f d²) and look for Re = ρ_f V d/μ.
3. Find the first time at which V reaches 0.95 V_t in the output of `simulate_fall`; for (c), interpolate the speed and the time at z = 2 m.
4. For a parachute, D = W gives the area directly: no Reynolds number enters if c_d is taken constant. For the sphere, compute Re at the answer and compare with each correlation's range.
5. The form-factor correlation divides friction (which scales with chord) by thickness; the thickness term grows like (t/c)⁴.
6. Think of what the constant-pressure wake leaves out: friction, pressure recovery, unsteadiness, and the real c_p ahead of separation.
