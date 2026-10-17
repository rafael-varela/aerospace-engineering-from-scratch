# P08.02  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The field is linear in x and y, so its gradient is a constant matrix. The local acceleration ∂u/∂t is zero because nothing depends on t.
2. (a) A streakline point is a pathline evaluated at the present time, with the release time as the parameter. (b) Through each point a steady ODE has exactly one solution curve. (d) Look at how the two edges of a small rectangle move apart and turn.
3. The velocity is the same everywhere, so every particle moves along the same shape of path; only the starting time differs. Write the streakline point released at τ and look at the sum of squares of its coordinates measured from a suitable centre.
4. The worst phase is the one that maximises the error; sample a full period of the gust. Expand the sine in the closed-form streakline to second order in x/U.
5. Divergence-free means area-preserving; strain still stretches. In u = ax, v = 0 the area changes at the rate tr G.
6. Ask for each method: is the image a single instant or a time exposure, and are the marked particles released from one place or spread through the flow?
