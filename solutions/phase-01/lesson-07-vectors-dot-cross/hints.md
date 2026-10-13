# P01.07  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Convert to metres and metres per second before anything else, or keep kilometres consistently and convert $|\mathbf h|^2/\mu$ at the end ($\mu = 3.986 \times 10^5$ km³/s²). The cyclic pattern $h_1 = r_2 v_3 - r_3 v_2$ fixes the other two components.
2. (a) Expand all three squared components of the cross product and all cross terms of $(a_1 b_1 + a_2 b_2 + a_3 b_3)^2$; everything pairs up into $(a_1^2 + a_2^2 + a_3^2)(b_1^2 + b_2^2 + b_3^2)$. (b) Expand the determinant along its first row. (c) Differentiate each component of (3) with the ordinary product rule, then regroup.
3. The air velocity has magnitude TAS along the heading; the condition "track east" is that the north component of the ground velocity vanishes. In (c) the horizontal part of the air velocity is $\sqrt{\text{TAS}^2 - w^2}$ with $w$ the climb rate.
4. For the triangular case $\int_0^s y(1 - y/s)\,\mathrm{d}y = s^2/6$. The order is $\log_2$ of the ratio of successive errors; the elliptic tip singularity costs half an order on a uniform grid.
5. Each case is one of: rigid body, central force, measurement precision, axis convention. For (b), compute the perigee altitude and ask what the air does there.
6. Three points fix a plane; the normal is $(\mathbf p_2 - \mathbf p_1) \times (\mathbf p_3 - \mathbf p_1)$. The plane of an unperturbed orbit passes through the Earth's centre, so $d = \hat{\mathbf n} \cdot \mathbf p_1$ measures how far the fitted plane misses it. For (e), the lever arms are the distances between the fixes (thousands of kilometres) against $|\mathbf r|$ times a velocity error.
