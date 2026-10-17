# P10.05  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. For a uniform sheet the u component is β/(2π) and the v component is −ln(r₁/r₂)/(2π). At (0.5, 0.5) the panel subtends a right angle.
2. (a) Write s = (s − ξ) + ξ and reuse I₁ and I₂. (c) N + 1 node values, N midpoints. (d) The jump in tangential velocity across a vortex sheet equals its strength.
3. With the first and last node at the rear point, the condition makes that point behave like a sharp trailing edge.
4. Use `kt_convergence` with your μ and τ; relative error = absolute error / exact c_l. Interpolate the panel count from the error table using the observed order.
5. `thin_airfoil_naca4(m, p, alpha)` gives the thin-airfoil values. Thickness changes the lift slope and the moment, not the camber's zero-lift angle much.
6. At a finite-angle trailing edge the exact velocity falls to zero over a tiny distance, as distance to a small power; no polygon resolves that.
