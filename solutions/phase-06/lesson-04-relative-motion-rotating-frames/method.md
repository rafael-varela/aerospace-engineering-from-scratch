# P06.04  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. r₀ = R U, v = 250 N (or 250 E); `five_term_acceleration(0, 0, Ω_vec, r₀, v, 0)`; the fictitious accelerations are −(coriolis) and −(centripetal); project on E, N, U. 2. Sideways: the component perpendicular to v in the horizontal plane. 3. atan(a_side/g).

## 2.
(a) For e fixed in B, (de/dt)_I = ω × e; for A = Σ Aᵢeᵢ use the product rule. (b) v_P = v_O + v_rel + ω × r; differentiate in I with (1): (dv_rel/dt)_I = a_rel + ω × v_rel; d(ω × r)/dt = α × r + ω × (v_rel + ω × r). (c) Dot (3) with v: v · (Ω × v) = 0; the other two terms are gradients. (d) r̈₁ = −2Ω × (v₀ + g t); integrate twice; at the equator Ω = (0, Ω, 0), v₀ = (0, v_n, v_u); evaluate the east component at T = 2v_u/g; check the vertical part for the time correction. (e) Substitute r = R_z(−Ω_z t)s; show s̈ = −(ω₀² + Ω_z²)s; initial conditions; the major axis turns with R_z(−Ω_z t).

## 3.
(a) ω = √(g/R); apparent weight (ωR ± v)²/R. (b) Start (r, 0) with velocity (0, ωr), r = 98 m; hit |position| = R at t = √(R² − r²)/(ωr); the spanner's angle atan(ωrt/r) against the floor's ωt; arc difference × R. (c) `solve_ivp` with RHS from `rotating_frame_acceleration(r, v, 0, (0, 0, ω))`, event radius = R. (d) √(2h/g).

## 4.
1. `launch_state`, `simulate_rotating`, `simulate_inertial`, `range_and_cross`. 2. `simulate_rotating(..., coriolis=False)`; differences; `first_order_impact_shift(v_enu, lat)` projected on the azimuth and its right-hand normal. 3. `simulate_rotating(..., Omega=0)`. 4. J = ½v² − μ/r − ½|Ω × r|² at both ends.

## 5.
(a) `first_order_impact_shift` with V = 330, 1000, 1600 m/s; range V²/g. (b) Feedback: the airliner's navigation is corrected by GPS, radio aids and the autopilot holding track; an inertial navigator integrates the uncorrected acceleration twice. (c) Air drag. (d) Compare flight times and required accuracy.
