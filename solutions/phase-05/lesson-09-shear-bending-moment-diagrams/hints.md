# P05.09  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Take moments about one support to get the other reaction directly. When you look for the maximum moment, remember that V can change sign by jumping across a point load as well as by passing smoothly through zero.
2. (a) Draw the slice with V and M on its left face and V + dV, M + dM on its right face, acting in the opposite senses. (b) Every term in the subtracted load has the same value and slope as the original at x = b, so the difference is zero there and stays zero. (c) The reaction at x = 0 is w₀L/6. (d) Substitute u = 1 − (x/s)² in the moment integral; for the force, the integral is a quarter of an ellipse's area.
3. Both trestles are vertical supports, so `solve_reactions([("pin", 1.0), ("roller", x2)], loads)` works even though the wing overhangs both of them. The engine sits near the inboard trestle. The overhangs are cantilevers: their moments are hogging.
4. `wc.loads(n, relief=False)` gives the panelled lift; `solve_reactions([("fixed", 0.0)], ...)` gives the root reactions. Errors from three panel counts that differ by a factor 4 give two observed orders, log(e₁/e₂)/log 4.
5. Fuel far from the root relieves more moment per kilogram than fuel near the root. When the fuel burns, the lift falls with the aircraft mass, but the lift's centre stays at 4s/3π.
6. Superpose: two half loads minus one whole load is a self-equilibrated set that only acts between the fittings. For (c), with no lift the only loads are the inertia loads, all pulling down.
