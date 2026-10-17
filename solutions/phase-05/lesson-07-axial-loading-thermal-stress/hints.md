# P05.07  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The tube's area is π(D² − d²)/4 with d = D − 2t. A restrained bar's thermal stress does not depend on its geometry at all.
2. (a) Write A(x) = A₀ + (A₁ − A₀)x/L and substitute u = A(x). (b) Let the junction move by u; write each segment's elongation in terms of u, then its force from (3). (c) The annulus area is π(u − a)(u + a) with u = D(s)/2 and a = d/2. (d) Ask what is common to the bolt and members in each case: the same displacement (parallel) or the same force (series).
3. Restrained: σ = −EαΔT with ΔT negative. With a bellows, the line and the bellows are springs in series sharing the free shrinkage αLΔT.
4. Each joint has its own preload F_i = T/(K d). Add the thermal change (the same for every joint), then use `separation_load` and `joint_loads`. The window's lower end comes from F_m ≥ 0 at ultimate when cold; the upper from F_b ≤ Fty × A_t at limit when hot.
5. At radius r the cut carries the centrifugal force of everything outboard: ∫ρAω²r′ dr′ from r to r₀ + L. In the FE model, look at the tip displacement for every element count before you look at the root stress.
6. (a) Think about the joint diagram: what fraction of an alternating load the bolt sees, and when separation happens. (b) `bolt_stiffness(8, 1.25, 20, 0, E)`. (c) Edges, flanges near a corner, closely spaced bolts. (d) Bolt stretch, measured mechanically or ultrasonically, or a strain gauge.
