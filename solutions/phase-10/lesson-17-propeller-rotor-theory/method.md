# P10.17  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. A = π(0.30)²/4; T = 1.5 + s/20. 2. `actuator_disk(T, 15.0, A)`: v, v_wake, P, η. 3. n = 100 rev/s: C_T = T/(ρn²D⁴), J = 15/(nD); `ideal_efficiency(C_T, J)`; compare. 4. `static_ideal_power(10.0, A)` and v_h = √(T/(2ρA)).

## 2.
(a) ṁ = ρA(V + v); T = ṁw; P = ½ṁw(2V + w) = T(V + v) ⇒ w = 2v. T = 2ρAV²a(1 + a) with a = v/V; substitute A = πD²/4, V = JnD; solve the quadratic in a; η = 1/(1 + a). (b) Axial: 4πrρF(V + u_a)u_a = ½ρ(V + u_a)²/sin²φ · Bc C_n ⇒ u_a/(V + u_a) = k ⇒ u_a = Vk/(1 − k). Torque: 4πr²ρF(V + u_a)u_t = ½ρ(V + u_a)(Ωr − u_t)/(sin φ cos φ) · Bc C_t r ⇒ u_t/(Ωr − u_t) = k′. Insert both in (V + u_a)cos φ = (Ωr − u_t) sin φ and multiply by (1 − k)(1 + k′)/(Ωr). (c) With cd = 0: dQ Ω/dT = Ωr tan φ = V + u_a. Minimise P = ∫(V + u_a)dT at fixed ∫dT, with dT = 4πrρ(V + u_a)u_a dr: stationarity gives (V + u_a)(V + 3u_a) = λ(V + 2u_a) on every annulus, an equation without r ⇒ u_a uniform.

## 3.
1. T = 1.2 × 9.81/4. 2. `rpm_for_thrust(blade, T, 0.0)`; `bem(blade, 0.0, n)` for P. 3. FM = `static_ideal_power(T, πR²)`/P. 4. P_elec = 4P/0.8; energy 11.1 × 2.2 × 3600 × 0.8 J; time = energy/P_elec. 5. Ratio of hover shaft power (4 rotors) to cruise shaft power; argue with P = T^{3/2}/√(2ρA) and the wing's induced power, which acts on a much larger "disk" moving at 12 m/s.

## 4.
1. New blade with p = 0.152 m. 2. `rpm_for_thrust`, then `bem` for J, η, P. 3. `sweep_advance_ratio(blade, n, J_grid)` and the maximum of η for both blades. 4. `bem(blade, 0.0, n_cruise)["T"]` for both. 5. Repeat with p = 0.100 + s/500.

## 5.
1. Double loop over n_el and relaxation; record `residual < 1e-10` and the iteration count. 2. `bem(..., n_el=160)["F"][-3:]`. 3. Argue from u_a ≈ (dT/dr)/(4πrρF(V + u_a)): the derivative of the update with respect to φ scales as 1/F. 4. `bem` for all n_el; errors against 1280; `observed_order`.

## 6.
(a) `reference_propeller(Section(cd0=...))` and the cruise solve for each value; tabulate. (b) Laminar separation bubbles and the Re-dependence of cl_max and cd (`P10.14`); judge the direction. (c) Interpolate α and cl at r = 0.75R from `bem(...)`; compare cl with the cl of the best cl/cd of the section model. (d) Thrust and torque against J on a test stand in a tunnel; a wake survey of u_a(r) behind the disk.
