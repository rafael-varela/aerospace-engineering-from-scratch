# P06.03  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Everything in (b)–(d) is one function call each; for (c), divide the kinetic energy by a distance and the momentum by a time.
2. Multiply Newton's law by v to turn d/dt into d/dx. For (b), differentiate F₀ + kx + cu and replace u′ using (4).
3. A bungee alone stores all the energy and gives it back. Halving k at the same runout is exactly the constant-force condition of 2.3; the damper takes half, the bungee keeps half, and the return trip loses more to the damper.
4. Loop `stopping_distance` and `peak_force` over a grid; for (b), runout grows with mass, so the root of runout − 105 m is unique.
5. Read the peak-g columns of the two tables in part 4; then ask what happens to stored energy at the end of the stroke and to a damper force at the instant of impact.
