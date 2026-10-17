# P05.03  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. At the roller only two forces are unknown: the inclined member and the horizontal one. Resolve vertically first, because the horizontal member has no vertical component.
2. (a) Resolve perpendicular to the line of the collinear pair. (b) Think of a girder with one diagonal moved to a panel that already has one. (c) The moment centre for each chord is the joint where the other two cut members meet.
3. Members and joints: m = 15, j = 9, r = 3. With a symmetric load, what is the shear at mid-span?
4. `fm, vol = depth_sweep(span, 6, 2000.0, ratios * span)`; `np.argmin(vol)`. Then `pratt_truss(6, span, 0.25 * span)` and `pratt_members(6)` to pick out the groups. The bay length is span/6.
5. `np.argmin` and `np.argmax` of the forces give the member indices; `member_geometry` gives their lengths. For (b), `np.vstack([em.members, [[4, 6]]])`.
6. For each, ask whether every member still has forces only at its two ends, and whether the structure is still determinate.
