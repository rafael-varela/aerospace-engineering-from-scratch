# P05.17  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) A = 160 + 4(H − 4); x̄ = (160 × 20 + 4(H − 4) × 2)/A; ȳ = (160 × 2 + 4(H − 4)(4 + (H − 4)/2))/A. (b) `polygon_properties(P)`; centroid Qy/A, Qx/A; I_xx,c = I_xx − Aȳ², I_yy,c = I_yy − Ax̄², I_xy,c = I_xy − Ax̄ȳ. (c) `principal_axes(Ixx_c, Iyy_c, Ixy_c)`.

## 2.
(a) Green: ∬(∂Q/∂x − ∂P/∂y) dA = ∮(P dx + Q dy). Q = x: ∮x dy; on edge i, ∫₀¹ (x_i + τΔx)Δy dτ = (x_i + x_{i+1})Δy/2; summing and rearranging gives ½Σ(x_i y_{i+1} − x_{i+1} y_i). Q = x²/2 for ∫x dA: ∫₀¹ (x_i + τΔx)²Δy/2 dτ; expand and sum. (b) J_xx = ∫(y² + z²) dm = ρt I_xx,c + m t²/12, similarly J_yy; J_zz = ∫(x² + y²) dm = ρt(I_xx,c + I_yy,c). Subtract. (c) Var(Σa_ix_i) = Σa_i²σ_i²; normal σ_i = t_i/3 gives 3σ = √Σ(a_it_i)²; uniform σ_i = t_i/√3 gives 3σ = √3 √Σ(a_it_i)².

## 3.
(a) A_slot = c w + π(w/2)², c = 14, w = 6.6. Polygon: two arcs of radius w/2 about (50 ± 7, 15). (b) Plate properties by `section_properties` on the rectangle, minus the slot (area A_slot at (50, 15)); mass ρ A t. (c) Two arcs, two lines; tangency and equal-radius constraints; dimension the slot width (or radius), the centre-to-centre length (or overall length) and the location of one arc centre (or the slot's centre) from the datum.

## 4.
(a) a = signed angle at (60, 15) between (60, 0) − (60, 15) and (20, 40) − (60, 15); replace ("dx", 0, 3, 20.0) by ("angle", 1, 2, 3, a); `constraint_status`. (b) Replace ("dy", 0, 6, 30.0) by ("dy", 5, 6, 20.0); status; height of the upper hole above the datum = 10 + 20, so its tolerance is the sum of two ±0.1 tolerances. (c) Append ("dy", 2, 3, 26.0); status; `solve_sketch(z, cons, n)` and the maximum residual.

## 5.
(a) g = 49.5 − 48.0 − 1.0; t = (0.10, 0.05, 0.02 + 0.005s); worst case Σt; RSS √Σt²; g − Δ. (b) `monte_carlo_stack(nom, [1, -1, -1], tol, 200000, s, dist)`; 3 × std; mean(y < 0.30). (c) Worst case needs Σt ≤ g − 0.30; solve for the shim's t. Arguments: the requirement says "never"; RSS accepts a small fraction of failures and depends on the distribution assumption.

## 6.
Findings to check: sloping-edge height 25 redundant with 40 and 15; top-edge length 20 locates the corner correctly but is a chained dimension from the left edge (acceptable, but it is the same as 20 from the datum only because the left edge is the datum); the lightening hole chained from a bolt hole instead of the datum; no general tolerance; projection symbol contradicts the view arrangement; thickness given (good); no material or units note. Corrected scheme: the lesson's table in 2.1. Check with `bracket_sketch` plus or minus the offending constraints.
