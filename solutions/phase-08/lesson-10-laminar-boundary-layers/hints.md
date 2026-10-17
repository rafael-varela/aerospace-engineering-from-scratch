# P08.10  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Everything scales with x/√Re_x = √(νx/U); the constants (4.91, 1.7208, 0.6641, 0.664) come from your Blasius solution. c_f ∝ x^(−1/2), so halving x multiplies it by √2.
2. (a) With ψ = √(νUx) f(η), compute ∂η/∂x = −η/(2x) and ∂η/∂y = √(U/(νx)); every term of the momentum equation then carries the same factor U²/x. (b) Substitute f = aF(aη) and count powers of a in each term. (c) Integrate the momentum equation from y = 0 to the edge, use continuity to replace v, and remember that U is constant.
3. The core is what the displacement thickness leaves: (W − 2δ*)². A wall that moves outward by dδ*/dx per metre exactly compensates; compare that slope with v/U at the edge.
4. Use the same η_max in the reference so that the only difference is the RK4 step. For (b), think about what f′(η_max) = 1 forces if the layer is not allowed its natural thickness. For (c), F′(ζ) has a limit; a = [F′(∞)]^(−1/2).
5. The march is first-order in x (backward Euler) and second-order in y; check which one changes the answer when you refine it. Richardson with order p: c_∞ ≈ c_fine + (c_fine − c_coarse)/(2^p − 1).
6. δ99/x = 4.91/√Re_x. For (b) think about the pressure distribution of an aerofoil, transition, and the surface; for (c), think of what a flush wall-shear or surface-temperature sensor would see.
