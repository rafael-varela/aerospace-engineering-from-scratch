# P06.09  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. m(|r|²E − rrᵀ) has off-diagonal entries only where two components of r are non-zero; `principal_axes` returns the axes as columns.
2. Dot Euler's equations with ω and with Iω; subtract I₂ × (2T) from H²; for the elliptic solution, compare the coefficients of sn dn, cn dn and sn cn on both sides.
3. `elliptic_parameters` gives the type and period; the angle between a body axis and H is arccos of (row of C) · Ĥ_i.
4. `elliptic_parameters`, `simulate`, `flip_times`, `elliptic_solution`, `invariant_drift`; (9) says how errors in the state move the flip time.
5. K/λ is a quarter of the period; the linear estimate ignores the factor 4 inside ln(4/k′) and the shape of the turn-over.
6. Frames and conventions; sensitivity to the initial state; how close the two moments are; Poinsot; the elliptic solution covers only ω; energy dissipation.
