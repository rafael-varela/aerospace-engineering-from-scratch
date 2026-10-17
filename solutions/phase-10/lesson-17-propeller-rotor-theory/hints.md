# P10.17  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The disk area is πD²/4. In (b), n = 100 rev/s, so C_T = T/(ρn²D⁴) and J = V/(nD); the two efficiencies must agree to round-off because (5) is (3) rewritten.
2. (a) Write the power two ways: the rate of work of the disk on the air passing through it, and the kinetic-energy flux. (b) Replace W² in the axial balance by (V + u_a)²/sin²φ and in the torque balance by (V + u_a)(Ωr − u_t)/(sin φ cos φ). (c) With a Lagrange multiplier on the thrust constraint, the stationarity condition on each annulus contains u_a but not r, so u_a must be the same everywhere.
3. In hover there is no useful work, so η = 0 and the figure of merit takes its place. Momentum theory: hover power scales as T^{3/2}, and a wing of span 1.5 m deflects far more air per second than four 10 in disks.
4. Build a new blade with `reference_propeller()` and overwrite its `beta` with `lambda r: np.arctan(p / (2 * np.pi * r))`. The thrust requirement is `uav_cruise_thrust()`.
5. In (8), u_a ≈ (dT/dr)/(4πrρF(V + u_a)) when u_a ≪ V: a small F multiplies every error in dT/dr by 1/F.
6. Look at what changes and what does not when cd0 moves: the thrust is set by the airframe, the rpm mostly by the lift, the power by both lift and drag.
