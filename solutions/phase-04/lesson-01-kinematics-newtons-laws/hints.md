# P04.01  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. T = 2v₀ sin θ/g, R = v₀² sin 2θ/g, H = (v₀ sin θ)²/2g. For (b) read `Trajectory.range`, `.time_of_flight`, `.apex_height`, `.impact_speed` and compute the impact angle from `atan2(-vz[-1], vx[-1])`. For (c): a lower apex is reached and left sooner.
2. (a) Separate variables: dv/(1 − v²/v_t²) = g dt, and ∫dv/(1 − v²/v_t²) = v_t artanh(v/v_t). (b) Set s(t) = h with u = gt/v_t, write t = t₀(1 + ε) and keep first order in ε and in gh/v_t². (c) d(1/v)/dt = k. (d) ẋ = v₀cos θ − k∫|v|ẋ dt with a positive integrand. For the angle, use the velocity angle ψ as the variable: the normal component of the equation is vψ̇ = −g cos ψ.
3. k = ρC_D A/(2m) with A the canopy's projected area πd²/4; `terminal_speed`, `vertical_drop_time`, `vertical_drop_state`. For (c) `simulate_projectile(10.0, 0.0, k, z0=30.0, max_time=100)`. For (d) the package's velocity relative to the ground is its velocity relative to the air plus the wind.
4. `predict_fall_times`, `frame_timing(1, fps)[1]`, `normalized_error`. Loop seeds, `reduce_timings`, `fit_terminal_speed`. u_CD/C_D = 2u_v/v_t. For (d) a constant bias looks like extra drag at low heights and less at high ones; the bean bag's fit will return a v_t far below its real one.
5. `heading_for_track(25, track, W, 40)` on each leg. For (b) write the two one-way times and add; for the crosswind case the ground speed is √(V² − W²) both ways. The 270° leg has a crosswind component W sin(40° − 270°).
6. Think: Magnus force; the drag crisis at Re ≈ 2 × 10⁵; a C_D far from 0.5 and a Reynolds number below 10³; moving air; density; and for (f) compute k and the 2 m excess before deciding.
