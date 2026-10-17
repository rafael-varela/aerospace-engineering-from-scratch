# P10.22  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Code verification asks "are the equations solved right?" with answers known exactly; solution verification asks "how big is the numerical error of this answer?"; validation asks "are these the right equations for this physics?", which needs measurements.
2. (a) Write f_i = f* + C h_iᵖ for three grids and eliminate f* and C. (c) The chord appears in the lifting-line equation only in the product c a₀; the reference area scales with the chord too.
3. The sensitivity ∂CL_α/∂κ is a central difference with `vlm_cla` at κ ± 0.01; u_input is that times u_κ.
4. `SCATTER_REL` sets the experimental scatter term; note that `validation_uncertainty` takes it as an argument, so pass it explicitly.
5. The GCI is meant to cover the difference between the fine-grid value and the grid-converged value with about 95 % confidence; compare it with |f(8 × 48) − f(16 × 96)|/f(16 × 96), which is itself an underestimate of the true error.
6. Think of the vortices that form along the side edges of a low-aspect-ratio wing at incidence; a linear lattice's wake stays flat.
