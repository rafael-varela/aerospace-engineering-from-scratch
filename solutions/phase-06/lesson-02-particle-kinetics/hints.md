# P06.02  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The pull-up spends one g on the weight; the level turn spends it on the vertical component of a banked lift. The neck holds up only what is above it.
2. (b) An accelerometer cannot feel gravity, only the push of whatever holds it: it measures a − g_vec. (c) The numerator of dh/dγ is the derivative of its denominator.
3. With n = 0 the normal equation gives γ̇ = −g cos γ / V and, with T = D, the tangential one gives V̇ = −g sin γ: that is exactly the motion of a thrown stone. The radius of curvature at the apex is V²/a_n with a_n = g.
4. Call the given `recovery_altitude_loss` with `speed_hold=False` for (c); for (d) simulate with `simulate_pullout` and dense output, then pass the positions only.
5. The instant-g loss scales as V²; for the lift limit, set n W = ½ρV²S C_L,max and solve for V.
