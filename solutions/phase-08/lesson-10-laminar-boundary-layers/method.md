# P08.10  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. U from s; Re_L = UL/ν. 2. `plate_boundary_layer([L/2, L], U, ν, ρ, consts)` for δ99, δ*, θ, c_f, τ_w. 3. Ratio of c_f at L/2 and L. 4. `friction_drag(U, L, b, ρ, ν, f″(0), sides=2)`. 5. Compare Re_L with the few-×10⁵ to few-×10⁶ range where flat-plate transition is usually observed.

## 2.
(a) u = ψ_y, v = −ψ_x; compute u_x, u_y, u_yy in terms of f and η; substitute; divide by U²/x. Boundary conditions: no slip (u = v = 0 at y = 0 for all x gives f′(0) = f(0) = 0) and u → U. (b) f′ = a²F′, f″ = a³F″, f‴ = a⁴F‴; f f″ = a⁴FF″; equation holds for any a. Choose a so that f′(∞) = 1. (c) Integrate u u_x + v u_y = ν u_yy from 0 to h; v(h) = −∫u_x dy; combine with U ∫u_x dy to obtain d/dx ∫u(U − u)dy = ν u_y(0) = τ_w/ρ. Then substitute θ = C√(νx/U) and c_f = 2f″(0)/√Re_x and solve for C.

## 3.
1. U; Re at x = 2 m; δ99 and δ* from the constants. 2. Core area (W − 2δ*)²; U_exit/U_in = W²/(W − 2δ*)². 3. Wall slope dδ*/dx = δ*/(2x); v/U from `blasius_integrals` v_edge/√Re_x. 4. Compare the growth of δ* for a turbulent layer qualitatively.

## 4.
(a) Loop over n; `shoot_blasius(η_max, n)`; reference `blasius_bvp(η_max, tol=1e-10)`; errors and `observed_order`. (b) Loop over η_max with n = 100 η_max; compare with the published value. (c) `rk4_blasius(1.0, 10.0, 2000)`; F′ at the end; f″(0) = F′^(−3/2); count RK4 integrations in your secant loop.

## 5.
1. y_max from the formula; four marches; `wall_cf_from_profile` on the last station; errors; orders. 2. Richardson with p from (1). 3. η = y√(U/(νx)) at two stations; interpolate f′ onto those η; maximum difference.

## 6.
(a) Solve 4.91/√Re_x = 0.1; x = Re_x ν/U. (b) Pressure gradient (`P08.11`), transition (`P08.12`, `P08.14`), roughness and insects (`P08.12`/`P08.14`), curvature and thickness. (c) A surface sensor along the chord: hot-film gauges, a china-clay or oil-flow pattern, or a surface temperature camera, read against the laminar prediction.
