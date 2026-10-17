# P06.03  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. KE = ½mv₀², p = mv₀; F = KE/100 m, a = F/m, t = 200 m/v₀. 2. `stopping_distance`, `stopping_time`, `peak_force`, `efficiency` with (F₀, k, c) = (0, 4225, 100). 3. `energy_average_force` (work–energy) and `momentum_average_force` (impulse–momentum). 4. `simulate_arrestment`, `energy_audit`: impulse/(m v₀), W_damper/KE0, E_k/KE0.

## 2.
(a) m v dv/dx = −F ⇒ d(½mu)/dx = −F; substitute F; integrating factor e^{αx}, α = 2c/m; particular solution a + bx; constant from u(0) = v₀². (b) F′ = k + cu′ = k − αF; solve; monotonic exponential approach to km/(2c). (c) F′ = 0 ⇔ k = αF = 2cF/m; with F = ½mv₀²/x_s, k = mv₀²/(2x_s²); a pure spring needs mv₀²/x_s². (d) u ≈ |u′(x_s)|(x_s − x) near the end, so 1/√u ~ (x_s − x)^{−1/2}; with x = x_s(1 − s²), dx = −2x_s s ds and √u ~ s√(|u′|x_s).

## 3.
(a) k = mv₀²/x², F = kx, t = (π/2)√(m/k); energy conservation for the return. (b) `brentq` on `stopping_distance(25, 25, 0, k/2, c) − 15`; compare with `constant_force_design(25, 25, 15)`; peak from `peak_force`. (c) RHS [v, −(kx + c v|v|)/m], event x = 0 falling; return speed; fraction (v_back/v₀)². (d) Make the k term not return energy (a one-way brake or ratchet), or add a pure brake.

## 4.
1. Nested loops over mass and speed with `stopping_distance`, `peak_force`. 2. `brentq(lambda m: stopping_distance(v0, m, 0, 4225, 100) − 105, 1e3, 1e5)`. 3. `constant_force_design(v0, m, 100)`.

## 5.
(a) Peak g at design (part 4 table 1); runout change with v₀ and with mass (table 2). (b) ½kx_s² is pushed back into the aircraft; a valve that holds the absorber at the stop, or a dissipative element instead of k. (c) The engagement force rise and cable dynamics. (d) Combine.
