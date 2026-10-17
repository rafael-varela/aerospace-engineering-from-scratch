# P05.05  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Two parts: the flange (centroid h − 4 mm above the base) and the web (height h − 8 mm). Ix of each about its own centroid, then A d² to the section centroid.
2. (a) Parametrise the edge as x = x_i + τ(x_{i+1} − x_i), y likewise, τ from 0 to 1. (b) Write y = y_c + d and expand the square; one integral is a first moment about the centroid. (c) Substitute the rotated coordinates into ∫v² dA and ∫uv dA; use double-angle formulas.
3. The hole is the inner circle cut off at x = −2 mm: its arc runs between the angles where R_i cos t = −2 mm. Use the angle range from arccos(−2/27.5) to 2π − arccos(−2/27.5). Mass per length is ρA.
4. `z_section_polygon(h, b, t)`, `z_section_parts`, `z_section_thin`. Each eigenvector of C is the *normal* of a principal axis. For (d), minimise the negative ratio over b between 4 and 40 mm.
5. `wing_box_parts(width=0.40, height=0.12, t_skin=..., t_web=0.002, cap_area=300e-6, n_stringers=4, stringer_area=60e-6, idealisation=...)`. A skin is a thin rectangle at y = ±(h − t)/2.
6. Run `polygon_section(z_section_polygon(...))` at 1.5 and 1.7 mm. For the holes, add two negative `rect` parts to `i_section_parts`.
