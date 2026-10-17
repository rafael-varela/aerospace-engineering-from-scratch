# P06.02  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. R = V²/(g(n − 1)); n_top = V²/(gR) − 1. 2. `level_turn(V, 6)`; convert φ and ω to degrees. 3. Seat: n m g with m = 85 kg; neck: n m g with m = 6 kg. 4. n_top = V²/(gR) − 1 ≥ 0 ⇒ V ≥ √(gR).

## 2.
(a) Resolve W, L, T, D along e_t and e_n; a = V̇ e_t + V γ̇ e_n. (b) f = a − g_vec; f·e_n = Vγ̇ + g cos γ = L/m; evaluate for n = 1, 0 and 1/cos 60°. (c) dh/dγ = V² sin γ / (g(n − cos γ)); integrate from γ₀ to 0; time from dt = V dγ/(g(n − cos γ)) with u = tan(γ/2); expand cos γ₀ ≈ 1 − γ₀²/2 and ln(1 + ε) ≈ ε.

## 3.
(a) R = V²/(g(n − 1)) at n = 1.8. (b) Show (1) with n = 0, T = D is ẍ = 0, ḧ = −g; T = 2V sin 45°/g; Δh = (V sin 45°)²/(2g); x = V² sin 90°/g; exit speed by energy (same height as entry). (c) `solve_ivp(flight_path_rhs, ..., args=(lambda t, y: 0.0, False, 0.0))` with an event on γ + 45°. (d) ρ = V_apex²/g with V_apex = V cos 45°.

## 4.
1. `pullout_altitude_loss`, `pullout_time` at γ₀ = −(30 + s)°. 2. `recovery_altitude_loss(250, γ₀, 7, r)` for r = 2, 6; extra = loss − closed form; (5) = V sin|γ₀| (n − 1)/(2ṅ). 3. Same with `speed_hold=False`. 4. `simulate_pullout` with `ramp_load_factor(1, 7, 6)`; `sol.sol(t)`; `kinematic_load_factor(t, x, h)` against the ramp.

## 5.
(a) Run `recovery_altitude_loss` at 225, 250, 275 m/s and at 6 g; percentages. (b) Larger: lift-limited pull, slower real onset, thrust left on, pitch lag, tailwind gusts. Smaller: speed brakes slowing the aircraft, a higher n available, a shallower real dive. (c) V = √(2nW/(ρ S C_L,max)). (d) Combine the speed and onset sensitivities.
