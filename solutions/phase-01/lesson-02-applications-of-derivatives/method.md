# P01.02  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. W = mg. 2. K = 1/(π e AR). 3. C_L* = √(C_D0/K), (L/D)max = 1/(2√(C_D0 K)). 4. V* = √(2W/(ρ S C_L*)); γ = arctan(1/(L/D)max). 5. Sink = V*/(L/D)max (small angle) or `sink_rate` (exact). 6. Range = 1000 m × (L/D)max. 7. `best_glide_cl_newton(cd0, k)`, `airspeed(cl, W, 1.225, S)`, `sink_rate(...)`, `glide_range(1000.0, cd0, k)`.

## 2.
(a) w = V sin γ ≈ V C_D/C_L and V = √(2W/(ρS)) C_L^(−1/2), so w = √(2W/(ρS)) C_D C_L^(−3/2). Set g = (C_D0 + K C_L²) C_L^(−3/2); g′ = 2K C_L · C_L^(−3/2) − (3/2)(C_D0 + K C_L²) C_L^(−5/2) = 0 gives 2K C_L² = (3/2)(C_D0 + K C_L²). The end behaviour (g → ∞ at both ends) or g″ > 0 confirms a minimum. (b) Form w*/w_ms = [C_D(C_L*)/C_L*^(3/2)]/[C_D(C_L,ms)/C_L,ms^(3/2)] with C_D(C_L*) = 2C_D0 and C_D(C_L,ms) = 4C_D0; the C_L ratio is 3^(1/2), raised to 3/2. (c) (L/D)max = ½ K^(−1/2) C_D0^(−1/2): f′ = −¼ K^(−1/2) C_D0^(−3/2) = −f/(2C_D0); f″ = (3/8) K^(−1/2) C_D0^(−5/2) = (3/4) f/C_D0². Error ≈ ½ f″ (δ C_D0)² = (3/8) f δ².

## 3.
(a) h = C_L^(1/2)/(C_D0 + K C_L²); numerator of h′: ½C_L^(−1/2)(C_D0 + K C_L²) − 2K C_L^(3/2) = 0 ⟹ C_D0 + K C_L² = 4K C_L² ⟹ C_D0 = 3K C_L². (b) Substitute into C_L/C_D: C_D = C_D0 + C_D0/3 = (4/3)C_D0 and C_L = C_L*/√3; divide by (L/D)max = C_L*/(2C_D0). (c) V ∝ C_L^(−1/2): the speed ratio is 3^(1/4). (d) `newton(hprime, hprime_derivative, 0.5)` with the derivative written out, or with a central difference of h′ as the derivative.

## 4.
(a) Loop over winds, store (V_b, 2W/(ρ S V_b²), G). `np.polyfit(vw_pos, v_b_pos - v_star, 1)`. (b) Compare V_b at −20 m/s with `airspeed(min_sink_cl(...))`; the tangent from a point on the axis must touch where the polar slopes down towards the point, which is impossible left of the minimum. (c) Wrap each call in try/except and record the outcome; the basin edge sits where `ld_second_derivative` changes sign. (d) Re-implement the loop of `newton` with a print, or call `newton` with increasing `maxiter` and `tol = 0`.

## 5.
(a) Quotient rule: dR/dt = (ḋ h − d ḣ)/h² = [−(V_b − V_w) h + d w]/h² = (w/h)[d/h − (V_b − V_w)/w]. Negative iff (V_b − V_w)/w > d/h. (b) For each wind: V_b from `headwind_best_speed`, w from `sink_rate_of_speed(V_b)`, G = (V_b − V_w)/w; substitute d = 26 000 m, h = 1500 m. (c) Height needed = d/G; time = d/(V_b − V_w). (d) Compare the margin at the pessimistic wind (13 m/s) with the height available; state the rule (e.g. "commit only if the pessimistic case still reaches").

## 6.
For each: name the assumption; quantify (C_L margin, kinetic energy ½V² converted to height, change in G between 10 and 4 m/s, √(W₂/W₁) on speeds); name the observable (buffet/stall warning, airspeed decaying, groundspeed increasing during descent, sink rate at a given IAS).

## 7.
(a) G = V/(w + u); G′ = 0 ⟹ (w + u) − V w′ = 0: the line from the origin tangent to the polar shifted down by u. (b) Copy `headwind_best_speed`, replace `w` by `w + u` in `f` only (F′ is unchanged since u is constant), keep V_w = 0. (c) Range = 1000 × G; compare G(u = 1) with G(V_w = 10). (d) From w + u = V w′: with u > 0 the right side must be larger, and V w′ increases with V above best glide (w′ grows as 3aV² − b/V²), so V must increase.
