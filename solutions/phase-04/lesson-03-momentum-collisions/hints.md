# P04.03  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. `collide_1d(m, v, 420e3, 0, e)` and `collision_energy_loss`. For (c) compare μ with m: μ/m = M/(m + M).
2. (a) Momentum before = momentum after, in the frame where the cart moved at v. (b) The trapezoid rule is exact for the average of left and right sums; their difference is Δm(f(m_f) − f(m₀)) for f = 1/m. (c) Galilean invariance: a velocity change cannot depend on the observer's frame. (d) ∫ln(m₀/(m₀ − ṁt))dt by substituting m = m₀ − ṁt.
3. Each throw is one momentum balance; `pellet_rocket_velocities(m0, ma, 5, 3)` does the one-at-a-time case. For (d) her kinetic energy is ½mv² with m = 120 kg.
4. `simulate_rocket(m0, mf, tw * m0 * G0 / u, u, k_drag=k)`; `velocity_losses`; q from `fl.z` and `fl.v`. For (c) `brentq` on max-q − 40 kPa between 1.2 and 4.
5. Build the stage list from f; `minimize_scalar` on −Δv. For (c) `brentq` on the payload. For (d) consider where the dead mass sits when each stage burns.
6. Use Δv = u ln R and its inverse; for (b) think of the mass during the burn; for (c) the bounce speed is √(2gh); for (d) momentum in the vented gas; for (e) the speed needed at 127 km and its direction.
