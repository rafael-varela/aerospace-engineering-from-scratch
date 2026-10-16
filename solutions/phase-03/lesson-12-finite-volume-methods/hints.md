# P03.12  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Δt = νΔx/a; steps = (distance/a)/Δt; D = ½aΔx(1 − ν); the variance grows by 2Dt. For (c), demand 2Dt ≤ (1.1² − 1)σ₀² and solve for Δx.
2. In (a) the only non-exact object is the time average of the interface flux. In (b) replace u_x and u_xx by centred differences and rearrange into a difference of fluxes. In (c) expand ū_{i−1} and the time level to second order and use u_tt = a²u_xx. In (d) use |1 − ν + νe^{−iθ}|² and expand.
3. Differentiate q = au in time using the conservation form, then multiply the advective form by a. Mass per metre is conserved; a ratio is carried unchanged. Conservation form with a periodic domain and non-periodic a is inconsistent at the wrap.
4. `sensor_history(u0, a, dx, dt, int(0.2/dt), scheme, i_s)` with `i_s = int(x_sensor/dx)`; the exact 50 % time at the sensor cell centre is (x_i − x_pulse − width/2)/a. Diffusion is symmetric about the front, so it shifts the 10 % point and not the 50 % point.
5. The upwind error on a smooth slug is proportional to D ∝ (1 − ν) once the grid is fine enough; at ν = 1 the stencil is a pure shift. Instability grows only the modes with |g| > 1, which for upwind at ν slightly above 1 are the shortest ones, and a smooth Gaussian contains almost none of them.
6. Set ½aΔx(1 − ν) = D_phys and solve for Δx; steps = t a/(νΔx). Boundedness is a property of monotone schemes only. Clipping adds mass; a limiter redistributes it.
