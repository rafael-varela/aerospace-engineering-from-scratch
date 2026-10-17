# P05.02  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Count first, then look for geometry that wastes a reaction: all lines parallel (no resistance across them), or all lines through one point (no resistance to rotation about it). A link from (4, 0) to (6, 0) lies on the x axis, which passes through the pin.
2. (a) Take moments about the meeting point: every column has a zero there. (b) Rank–nullity applies to A (n columns) and to Aᵀ (n_eq columns). (c) Take moments about the root pin; for (c)'s second part, write the horizontal pin force as T times the strut's horizontal direction cosine.
3. Both hinge reactions cut the hinge line. The strut's moment about the hinge line is its force times the y·u_z − z·u_y combination at its door end.
4. Use `StrutWing(y_strut=..., lift_half=..., wing_weight=60*G0, n=3.8)`. At the minimum of the pin resultant one of its components vanishes; ask which one can.
5. The torque is carried only by the links' tangential components. For (b), `np.linalg.svd(A)[2][-1]` is a null vector of a 6 × 7 matrix of rank 6.
6. Think of what happens to each design when the structure around it deforms or heats up, and when one member breaks.
