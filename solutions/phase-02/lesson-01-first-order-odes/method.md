# P02.01  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. t_b = (m₀ − m_dry)/ṁ. 2. Δv = v_e ln(m₀/m_dry); gravity loss g t_b; v_b = Δv − g t_b. 3. h_b = v_e t_b − (v_e m_dry/ṁ) ln(m₀/m_dry) − ½ g t_b². 4. Apogee h_b + v_b²/(2g). 5. Build the parameter dict, call `vacuum_velocity(t_b, p)`, `vacuum_altitude(t_b, p)`, then `simulate_rocket(p, t_b, 400, "rk4", drag=False)` and subtract.

## 2.
(a) h = ∫₀ᵗ [v_e ln(m₀/(m₀ − ṁτ)) − gτ] dτ; substitute u = m₀ − ṁτ in the first integral, integrate ∫ ln(m₀/u) du by parts (u′ = 1, v = ln(m₀/u)), evaluate between u = m₀ and u = m(t), tidy. (b) dh = v dt = v dv/(dv/dt) = −v dv/(g + kv²); integrate from v_b to 0; the substitution w = g + kv² gives a logarithm. The formula holds while v > 0 because the drag term was written −kv², which assumes upward motion. (c) Put x = kv_b²/g, expand ln(1 + x), divide by 2k; the second term is the drag loss of altitude to first order in k; the series converges only for x < 1, so check which of the two k values is inside that range.

## 3.
(a) Order one; linear (T appears to the first power, coefficient constant); autonomous (no explicit t); also separable because it is autonomous. Standard form dT/dt + (G/C)T = (P + GT_a)/C. (b) Multiply by e^(Gt/C), integrate, apply T(0). Time constant τ = C/G; steady state T_∞ = T_a + P/G. (c) Set T(t) = 15 and solve for t: t = −τ ln((15 − T_∞)/(T₀ − T_∞)). (d) `integrate_fixed_step` with the right-hand side (P − G(T − T_a))/C, 200 RK4 steps to the time from (c); compare the final temperature with 15 °C. `equilibria(lambda T: (P - G*(T - Ta))/C, grid)` returns one stable root at T_∞.

## 4.
(a) `q = dict(sounding_rocket(), cd=0.3 + s/100)`; `ts, ys = simulate_rocket(q, 150.0, 3000, "rk4")`; `i_b = argmin|ts − 20|`, `i_a = argmax(ys[:, 0])`. (b) Repeat with the four step counts and tabulate `ys[:, 0].max()`. (c) `rho_b = air_density(ys[i_b, 0])`; `k = rho_b*cd*area/(2*m_dry)`; `quadratic_drag_coast(v_b, k, g)` and `v_b**2/(2g)`; compare with `ys[i_a, 0] − ys[i_b, 0]`. (d) `solve_ivp(lambda t, y: rocket_rhs(t, y, q, True, True), (0, 20), [0, 0, 250], rtol=1e-10, atol=1e-8)`.

## 5.
(a) `terminal_velocity(100, air_density(h), 0.5, area, 9.81)` for each h. (b) `simulate_rocket(sounding_rocket(), 400.0, 8000, "rk4")`; find the first index where h < 0 (impact) and the index where h is nearest 18 km on the way down; read v there. (c) τ(v*) = 1/(2k|v*|) at each altitude against the time to fall one scale height, 8500/|v|. Where τ ≫ that time the velocity cannot follow the equilibrium. (d) Think about the attitude of a falling rocket body (it does not fall nose-first by itself) and about the drag coefficient of a tumbling cylinder.

## 6.
(a) Euler: y_N − y₀ = h Σ_{n=0}^{N−1} f(t_n), the left Riemann sum. Trapezoid = left sum + ½h[f(b) − f(a)]; trapezoid error is O(h²); so Euler's error is −½h[f(b) − f(a)] + O(h²) (sign depends on which way you subtract). Compare ½(f(t_b) − f(0)) with the lesson's error/h. (b) Unroll: e_n ≤ ½h²M Σ_{j<n}(1 + hL)^j = ½h²M[(1 + hL)ⁿ − 1]/(hL); then (1 + hL)ⁿ ≤ e^{nhL}. (c) L = max|∂f/∂v| = 2k|v|_max with |v|_max the terminal speed; M = max|d²v/dt²| = max|f′(v) f(v)|; evaluate the bound for each h, run `integrate_fixed_step(..., "euler")` against `solve_ivp(rtol=1e-12)`, take the ratio. The explanation is the sign of ∂f/∂v.
