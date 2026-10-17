# P02.12  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. u_x = v_y gives v up to a function of y; u_y = −v_x fixes that function. The result should look like the imaginary part of a polynomial in z = x + iy: expand (x + iy)³.
2. For (a), write the trailing-edge point as ζ₀ + Re^{−iβ} and multiply the velocity by e^{−iβ} so that the free-stream term becomes a sine. For (c), p = p_∞ + ½ρ(U² − |w|²) and, on a streamline, w dz = w̄ dz̄ (both real). For (d), only the 1/z term of w² has a residue.
3. A vortex at ih and its mirror image at −ih with opposite sign: by symmetry the vertical velocities on the real axis cancel. The velocity "induced on" a vortex is the velocity of all the *other* singularities at its position.
4. Everything is a function call; keep the chord from `airfoil_geometry` for C_l, and take the lift slope from two incidences, not from 2π.
5. Potential-flow speeds scale linearly with U, so M_local/U is a property of the section and incidence. For (d), both the numerator w_ζ and the denominator dz/dζ vanish linearly as the point approaches c; what happens to the ratio when Γ is slightly wrong?
6. Verification asks "did I solve the equations right?"; validation asks "are these the right equations?". Only a measurement answers the second.
