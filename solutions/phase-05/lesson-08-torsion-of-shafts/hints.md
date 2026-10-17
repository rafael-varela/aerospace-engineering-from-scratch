# P05.08  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. `shaft_torque` returns N·mm. The inner surface has r = d_o/2 − t. The Bredt–Batho stress uses the mid-line radius and is the average stress through the wall.
2. (a) Start from the rigid rotation of a section: how far does a surface line tilt? (b) Take moments of q ds about any point; p ds is twice a triangle's area. For the twist, the strain energy per unit volume in shear is τ²/2G. (c) Expand r_o⁴ − r_i⁴ with r_o = r_m + t/2, r_i = r_m − t/2. (d) The slit tube is one strip of width 2πr_m.
3. Build the nose as points on the semicircle from the top of the web round to its bottom; the closing edge from the last vertex back to the first is the web. Give that edge its own thickness.
4. Reuse the loop in `main`; the crossover is where the twist diameter equals the ultimate diameter: solve the twist formula for φ with d set to the ultimate value.
5. A re-entrant corner of a hole makes the stress function's gradient singular there; second-order accuracy needs a smooth solution.
6. (a) Think about what the engine does in a transient. (b) Thin walls and fast rotation. (c) Closed versus open, and how a doubler or a frame can close the load path. (d) Stress concentration and alternating stresses.
