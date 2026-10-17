# P05.01  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. F = (−1000, 3000, 100s) N. 2. M = (r_y F_z − r_z F_y, r_z F_x − r_x F_z, r_x F_y − r_y F_x); |M| = √(M·M). 3. M_x = x̂ · M; for the axis along (1, 0, 1), û = (1, 0, 1)/√2 and M_u = û · M. 4. `moment_about_point(r, F)`, `moment_about_axis(r, F, np.zeros(3), [1, 0, 1])`. 5. Dot products: by the cross-product definition both vanish.

## 2.
(a) û · ((r − P − λû) × F) = û · ((r − P) × F) − λ û · (û × F); the last term is a triple product with a repeated vector. (b) M_B = Σ (r_i − B) × F_i = Σ (r_i − A) × F_i + (A − B) × Σ F_i; for a couple R = 0. (c) Planar: R = (R_x, R_y, 0), M_O = (0, 0, M). Look for P with P × R = M_O; the solution set is a line parallel to R; p = R × M_O/|R|² satisfies p × R = M_O (expand the triple product). R = 0: a pure couple.

## 3.
1. m = 10 + s/5; W = m g. 2. Unit vectors a = (cos 80°, sin 80°) along the rail and n = (−sin 80°, cos 80°) across it. 3. Components of W: W_a = −W sin 80°, W_n = −W cos 80°. 4. Stop: S = −W_a. 5. Moment about the stop of W at (0.9 a + 0.05 n): 0.9 W_n − 0.05 W_a. 6. N₁ + N₂ = −W_n and 0.3 N₁ + 1.1 N₂ = −(moment of W). 7. Check with `solve_fbd_2d([cg], [W], [b1, b2, 0], [n, n, a])` and the residual. 8. Repeat with zero offset to see its effect.

## 4.
1. m = 600 + 20s; `V, D = landing_load_case(m)`. 2. `g = GearGeometry(contact=np.array([0.10, 0.0]))`. 3. `gear_reactions(g, V, q*V)` for q = 0, 0.25, 0.5; `brace_force_by_moments` for each. 4. At D = 0: M_A of V = 0.10 V (anticlockwise); T = −M_A / (moment of a unit brace force about A). 5. Zero brace: the line through C with direction (−D, V) passes through A when D/V = (C_x − A_x)/(A_y − C_y). 6. Compare with `gear_reactions(GearGeometry(), V, 0.25*V)`.

## 5.
(a) W = w₀ s ∫₀^{π/2} cos²u du; ∫ y w dy = w₀ s² ∫ sin u cos²u du. (b) Split the chord into a rectangle λc_r and a triangle (1 − λ)c_r; centroids s/2 and s/3. (c) Build w_e = `elliptic_lift(y, s, L)` and w_t = L c(y)/∫c dy; average; `distributed_load_resultant`. (d) Root bending moment = L ȳ for each.

## 6.
(a) Spin-up drag peaks then reverses (spring-back); dynamic overshoot of V. (b) The brace becomes a three-force member: its end forces are no longer along the pin line; the leg sees a force at B with two unknown components. (c) A moment at A takes part of D ℓ; ΣM_A gives T d = D ℓ − M_A. (d) A wrong lever for one force, a forgotten force that has a moment about C (the trunnion reactions do).
