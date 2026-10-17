# P08.18  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. h = 1/(N − 1); the cell Reynolds number is U h Re in the lesson's scaling. On the bottom wall ψ_xx = 0, so the exact ω is −ψ_yy there; differentiate y²(1 − y)² twice and set y = 0.
2. (a) Differentiate the v-equation with respect to x and the u-equation with respect to y, and subtract; the terms (∂u/∂x + ∂v/∂y)ω vanish by continuity. (d) Write D_x and D_y as Kronecker products of a 1-D difference matrix with the identity; two such products acting on different indices commute.
3. Use `centre_values` and `vortex_centre` on each grid and `observed_order` on each triple; the GCI is 1.25 |f₃ − f₂|/(|f₃|(2^p − 1)).
4. Newton converges only inside a basin round the solution; the Stokes field at Re = 1000 is far from the Re = 1000 flow. Each continuation step must start inside the next basin.
5. Compare the diffusion limit h²Re/4 with the convection–diffusion limit 2/(Re U²) on each grid: whichever is smaller sets dt, and the physical time to reach the steady state is about the same on every grid.
6. The lid's velocity jumps from U to 0 at the corner over one cell, so the shear there is about U/h. For (b), compare the size of your code's own error (the GCI) with the size of the difference from the table.
