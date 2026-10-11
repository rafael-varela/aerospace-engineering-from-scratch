# P00.13  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The impact speed √(2 g h) does not depend on the mass. Compute it once; compute m g h four times. Then check each of the three envelope rows (mass, height, energy) separately.
2. For (a), put the two limits into m g h at the same time. For (b), solve ½ m v² = 40 J for v and notice how the answer moves when m gets small.
3. Elastic energy is ½ k x², and the speed limit in the envelope applies to a *launched* object. If the speed comes out too high, the substitute control is a smaller x: solve ½ k x² = ½ m v² for x with v = 5 m/s.
4. Run the checker first and copy its output verbatim. Then read the box as a reviewer would and ask of each cell "could a stranger act on this?". For (c), the parse step already gives you an empty dictionary; the question is only where to turn that into a finding.
5. 1 mA·h at 1 V is 3.6 J. Then separate *how much energy* from *how fast it can get out and through what*.
6. Tip speed is ω r with ω in rad/s, so convert rev/min with 2π/60. The checker takes an `envelope` argument; a copy of the default dictionary with one key changed is the tool for (b).
7. Read each NAR rule and ask "which of our nine questions does this answer?" Several rules answer two.
