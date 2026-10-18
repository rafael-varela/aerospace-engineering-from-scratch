# P10.09  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. b′ = πb/4; Γ = W/(ρVb′). 2. Centre downwash Γ/(πb′). 3. `horseshoe_velocity` at (4.5, 0, 0) and (10⁵, 0, 0); ε = w/V. 4. α_i = CL/(πAR); compare 2α_i with the far-field angle.

## 2.
(a) dl = r0 dt, r = P − A − t r0; the integral of 1/|r|³ along a line gives the difference of direction cosines; rewrite in terms of r1, r2. Symmetric case: r1 × r2 has size 2Lh, the bracket 2L/√(L² + h²)… (b) B = A + Ld, L → ∞: r2/|r2| → −d. (c) y = −(b/2) cos φ; dΓ/dy in φ; the integral with n = 1. (d) Minimise CD/CL and CD/CL^{3/2} in CL.

## 3.
1. Γ0 = 4W/(ρVπb); separation πb/4. 2. w = 2 × Γ0/(2π × separation/2). 3. v = Γ0/(2π × 5 m); compare with 55 m/s and with problem 1.

## 4.
1. `cosine_layout(10, n)`, φ = arccos(−2y/b), Γ, `lift_and_drag`, `span_efficiency`. 2. The ratio Di_near/Di_trefftz. 3. Sign of a₃ and the shape near the tips.

## 5.
1. `best_glide(CD0, e, 8.333)`; `speed_for_cl(750, 12, CL)`. 2. `drag_and_power(55, ...)`. 3. Plot; compare with 26.4 m/s.

## 6.
(a) e ≤ 1 for planar wings; the tip control point close to the tip vortex. (b) Downwash at the tail vs Trefftz drag. (c) Separate the computed wing e from the whole-aircraft Oswald factor.
