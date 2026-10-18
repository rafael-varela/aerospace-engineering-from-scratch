# P08.12  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. `flat_plate_transition(U, ν, N_crit)` does all of (a); for (b) the only change is N_crit = `mack_ncrit(0.003)`.
2. (a) Write the x- and y-momentum equations for the disturbance, take ∂/∂y of the first minus ∂/∂x of the second, and substitute the normal mode. (b) The imaginary part of ∫(U″/(U − c))|φ|²dy must vanish if c_i ≠ 0; what does that require of U″? (c) e^9 is a pure number.
3. Re = U_c h/ν = 5772.22 gives U_c directly; convert with U_b = 2U_c/3 and the length you are asked for.
4. Converged digits: compare each N with the next; stop when the change is below the digit you report. For the spectrum, compare the same eigenvalues at N = 60 and N = 100.
5. N_crit enters only through Re_θ,tr = Re_θ0 + N_crit/(dN/dRe_θ); H enters through both Re_θ0 and the slope, and steeply.
6. Think of what reaches the boundary layer without growing as a Tollmien–Schlichting wave: large disturbances, roughness, crossflow on a swept wing, contamination along the leading edge.
