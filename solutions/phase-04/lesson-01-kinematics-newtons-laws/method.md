# P04.01  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. v₀ = 5 + s/10; θ = 40°. T = 2v₀ sin θ/g, R = v₀² sin 2θ/g, H = v₀² sin²θ/(2g); `vacuum_state(T, v0, 40)` must return z = 0 and x = R. 2. `tr = simulate_projectile(v0, 40.0, 0.1396)`; read the four attributes; impact angle = degrees(atan2(−tr.vz[−1], tr.vx[−1])). 3. Compare each drag value with its vacuum value as a ratio; the time of flight is shorter because the apex is lower (less vertical speed survives the climb), and the descent from a lower apex, even slowed, is shorter in time.

## 2.
(a) dv/dt = g − kv² = g(1 − v²/v_t²); separate; integrate from v = 0 at t = 0 (this is "from rest": the constant of integration is zero); invert artanh; integrate v = ds/dt once more from s = 0. (b) h = (v_t²/g) ln cosh(gt/v_t); expand to u⁴; write t = t₀(1 + ε), expand u² and u⁴ to first order in ε, equate to gh/v_t² and solve for ε; evaluate at h = 2 m, v_t = 8.38 m/s and compare with `vertical_drop_time(2.0, 8.38)`. (c) dv/dt = −kv² ⇒ d(1/v)/dt = k ⇒ 1/v = 1/v₀ + kt; integrate v for x; set v = v₀/2 and read x. (d) Integrate ẍ = −k|v|ẋ: ẋ(t) = ẋ₀ − k∫₀ᵗ|v|ẋ dt′ with a positive integrand while ẋ > 0, so ẋ decreases. For the angle: change variable to the velocity angle ψ (vψ̇ = −g cos ψ); show du/dψ = ku³/(g cos³ψ) > 0 for u = ẋ; write dz in terms of dψ; impose Δz = 0; assume α ≤ θ and pair ψ with −ψ on [−α, α] to reach a contradiction.

## 3.
(a) A = π(1.2)²/4; k = 1.2 × 1.5 × A/(2m); v_t = √(g/k); t = `vertical_drop_time(30, v_t)`; impact speed from `vertical_drop_state(t, v_t)[1]`; vacuum t₀ = √(60/g), impact √(2g × 30). (b) A = π(0.4)²/4, C_D = 1.0; same steps; the time depends on the whole fall (most of it near v_t only when v_t is low), the impact speed on the end state. (c) `simulate_projectile(10, 0, k, z0=30, max_time=100)`: `.range`, `.vx[-1]`, `.vz[-1]`, `.time_of_flight`; vacuum 10√(2h/g); `horizontal_decay_state(T, 10, k)[0]`; order them: the pure decay ignores the vertical speed's contribution to |v| (less drag on x), the vacuum ignores drag altogether. (d) Add the wind vector times the flight time to the landing point, in the ground frame.

## 4.
(a) `t_vac, t_drag = predict_fall_times(hs, k)`; `u = frame_timing(1, f)[1]`; E_n = (t_drag − t_vac)/u. (b) For seed in range(20): `make_synthetic_drop_data(k, hs, f, n, seed)`; `reduce_timings(x, u)` per height; `fit_terminal_speed(hs, means, us)`; collect v_t and u_v; mean, standard deviation, mean u_v; u_CD/C_D = 2u_v/v_t. (c) Pick the cheapest combination with 2u_v/v_t ≤ 0.2; the frame term is Type B and does not average down, so n helps only the Type A part. (d) Add `release_bias_s=s/1000`; the fit absorbs the bias into v_t (a constant offset is partly fitted as a 1/v_t² term); the bean bag, which should fit a huge v_t, fits a small one too, which flags the bias.

## 5.
(a) Loop the four tracks through `heading_for_track(25, trk, W, 40)`; leg time = 2000/ground speed; sum. (b) Out: d/(V − W); back: d/(V + W); sum = 2dV/(V² − W²) = (2d/V)/(1 − W²/V²). Across: ground speed √(V² − W²) each way; sum = (2d/V)/√(1 − W²/V²); compare the two factors for W < V. (c) The 270° leg has crosswind W sin(40° − 270°) = W sin(−230°) = W sin 130° = 0.766 W; impossible when 0.766 W > V. (d) 1200 s divided by the lap time, floored.

## 6.
Name the assumption, the sign, the measurement: (a) no spin; (b) constant C_D (drag crisis); (c) a sphere with C_D = 0.5 and Re > 10³; (d) still air; (e) ρ; (f) compute k and (5) for the bearing and compare the excess with 0.41/f.
