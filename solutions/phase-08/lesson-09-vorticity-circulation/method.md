# P08.09  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. r_max = 1.1209 r_c; `lamb_oseen_utheta(r_max, Gamma, rc**2)`. 2. `lamb_oseen_circulation(R, Gamma, rc**2)` for each R; percentage; `circulation_polygon(circle_loop((0, 0), R, 2000), lambda p: vortex_velocity(p, Gamma, rc**2))`. 3. `lamb_oseen_vorticity(0, Gamma, rc**2)`; expand u_θ for small r and compare u_θ/r with ω(0)/2. 4. `circulation_polygon` on the off-centre circle; `disk_area_integral` of the vorticity on the same disk.

## 2.
(a) Material derivative of ∮u·dl; substitute Euler; write −∮ dP − ∮ dΦ + ∮ d(|u|²/2) = 0. The viscous term adds ν∮∇²u·dl = −ν∮(∇×ω)·dl. (b) ∂ω/∂t through r_c²(t); ω_r and ω_rr; compare the two sides. u_θ = (1/r)∫₀^r ω r′ dr′; Γ(r) = 2πr u_θ. (c) Biot–Savart with r₀ = B − A, r₁ = P − A, r₂ = P − B; for the infinite line both cosines → ±1; for the semi-infinite, one cosine is 0. The horseshoe: two semi-infinite legs at b′/2 (bound vortex contributes nothing at its own centre); far downstream two infinite legs.

## 3.
(a) W = mg; b′ = πb/4; Γ₀ = W/(ρVb′); descent Γ₀/(2πb′). (b) w = Γ₀/(2π × 5 m), opposite signs at the two tips; angle = atan(w/V). (c) r_c² = 4νt; r_max = 1.1209 r_c; `lamb_oseen_utheta(r_max, Gamma0, rc2)`; halving: 1 − e^{−25/(4νt)} = 1/2. (d) Compare the laminar time with the minutes of separation rules; think about turbulence, instability and the ground.

## 4.
1. Turnover = 2π r/u_θ(r) at r = 0.5 + s/50. 2. Loop over t and ν: `material_loop_exact`, `circulation_polygon` with the velocity at the current core size, `material_circulation_by_area`, and the sum of segment lengths. 3. Repeat the inviscid four-turnover case with 1000, 2000, 4000 points; ratios of drifts. 4. Compare the differential rotation of the spreading vortex with that of the frozen one.

## 5.
1. b′ = π × 10/4; Γ₀ = W/(ρVb′). 2. `horseshoe_velocity` at (l_t, 0, 0) and (l_t, 0, 0.5); ε = atan(−w/V). 3. `biot_savart_segment` for the bound vortex alone; legs = total − bound. 4. The y sweep at z = 0.5 m. 5. Scale the downwash with Γ.

## 6.
Γ = πcV sin α; w(d) = Γ/(2πd); w/(Vα); solve for d at 10 %. For (c), list what the model contains (Kelvin, the induced downwash) and what it lacks (a continuous wake, the time the circulation takes to grow).
