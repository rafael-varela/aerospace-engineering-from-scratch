# P06.12  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The area of a libration orbit is twice the integral of p dq between ±a; the largest p is at q = 0. On the separatrix H = ω₀², so p(0) = 2ω₀; the bob's speed is l times θ̇.
2. (a) Write dH in two ways and compare coefficients. (b) Differentiate ∮ p(q, E) dq with respect to E; ∂p/∂E = 1/p and dq/p = dt. (c) The divergence of (∂H/∂p, −∂H/∂q) involves mixed partial derivatives. (d) A 2 × 2 triangular matrix with ones on the diagonal has determinant 1.
3. L contains no t, so the energy function h = q̇ ∂L/∂q̇ − L is conserved, but it is not T + V here because the constraint (the spinning hoop) does work. Equilibria are where ∂H/∂θ = 0 with p = 0.
4. The exact flow's area change is the polygon's error, not physics. For RK4, a per-step area defect of order h^m accumulates over T/h steps.
5. δE = v_p δv at fixed r, and δa = 2a² δE/μ. The along-track law (13) needs δa, not δv.
6. Volume is not size in every direction. An adaptive RK45 is not symplectic. Over one orbit, accuracy matters more than structure. Drag is not Hamiltonian.
