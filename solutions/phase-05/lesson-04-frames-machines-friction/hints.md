# P05.04  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Resolve the weight along and normal to the slope. The friction available is μ times the normal force on the *braked* wheels only.
2. (a) An element of angle dθ: the two tensions differ by dT, and their components toward the drum centre give the normal force. (b) When the joint turns, the reaction leans by the friction angle; its line misses the centre by r sin φ_f. For the pulley, square the moment equation. (c) The only admissible motion has one degree of freedom; write the work of the grip force and of the hinge moment for it.
3. The centre pin is held by two links and the actuator. By symmetry both links carry the same force. For (c), x_slider = 2l cos θ and y_pin = l sin θ.
4. The small-angle gearing is the product of two arm ratios. For (c), `brentq(lambda b: delta_deg(b) + 25, 0.06, 0.14)`.
5. `pulley_ratio(mu_b, r_b, R, theta)`, `capstan_ratio(mu, theta)`, `cable_run_ratio` on a list of `("pulley", mu_b, r_b, R, theta)` and `("fairlead", mu, theta)` tuples. A seized pulley is a fairlead with its full wrap.
6. Friction moments scale with the pin radius and with μ; the band scales with the load too.
