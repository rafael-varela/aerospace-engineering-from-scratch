# P05.16  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. With the dummy loads, M(x) = q(L − x)²/2 + Q(L − x) + C, measured from the root. Differentiate U under the integral sign before setting Q = C = 0; the integrand becomes M ∂M/∂Q / EI.
2. (a) U = ½PᵀFP; differentiate with respect to P_j and use the symmetry you prove in (b). (b) Apply P_a then P_b, and P_b then P_a; the final state and the energy are the same. (c) Release the prop and write M(x) with R; set ∂U/∂R = 0.
3. The moment at the point (R cos φ, R sin φ) is P times the horizontal distance from the load's line of action; the arc element is R dφ. For the horizontal dummy load, the lever arm is the vertical distance R − R sin φ.
4. Build F as a 2-vector from the angle; the unit vector along F is the direction of the unit load. Reactions at a support are the opposite of the sum of the bars' pulls on it: with the lesson's C, the reaction is (C N) at the support degrees of freedom (use all eight degrees of freedom when building C).
5. Loop over `np.arange(2.0, 4.0001, 0.05)`, keep the a with the smallest `np.max(np.abs(r["M"]))`. For a uniform lift pass `lambda x: np.full_like(x, 18.6e3 / 5.0)`.
6. Every flexibility the hand model left out or added changes the ratio of the strut's stiffness to the spar's. Ask which way each one moves that ratio.
