# P05.17  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The two rectangles are 40 × 4 along the bottom and 4 × (H − 4) up the left side; their centroids are at their centres. An L-section has no axis of symmetry through its centroid parallel to x or y.
2. (a) Q = x, P = 0 gives the area; on an edge parametrise x = x_i + τ(x_{i+1} − x_i), y likewise, 0 ≤ τ ≤ 1. For ∫x dA use Q = x²/2. (b) Write J_zz, J_xx, J_yy as integrals of (x² + y²), (y² + z²), (x² + z²) dm and use ∫z² dm = m t²/12. (c) The variance of a sum of independent variables is the sum of the variances.
3. The two semicircles together make one full circle. By symmetry the slot's centroid is at its centre. Build the obround polygon from two half-circle arcs of 200 points each.
4. The angle constraint tuple is ("angle", i, j, k, a): the angle at point j between the lines to i and to k, signed, as `np.arctan2(cross, dot)`. Compute it from the nominal geometry so that the scheme stays consistent.
5. The gap is housing − bracket − shim. Worst-case smallest gap = nominal − Σt_i; RSS smallest gap = nominal − √Σt_i².
6. Count: a fully constrained outline needs 10 equations (including the datum and orientation) and the holes 9. Every extra number on the drawing is redundant; every location that is not from the datum is chained.
