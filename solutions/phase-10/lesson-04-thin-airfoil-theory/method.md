# P10.04  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. dz/dx = 4h(1 − 2x) = 4h cos θ. 2. A₀ = α − (1/π)∫4h cos θ dθ = α; A₁ = (2/π)∫4h cos² θ dθ = 4h; Aₙ = 0 for n ≥ 2. 3. α_L0 = −A₁/2 = −2h; c_l = π(2α + 4h); c_m,c/4 = −πh; ideal c_l = πA₁; x_cp = ¼(1 + π(A₁ − A₂)/c_l). 4. `thin_airfoil(lambda x: 4*h*(1-2*x), radians(3))`.

## 2.
(a) Substitute γ(θ) into the integral equation with dξ = (c/2) sin φ dφ; use the Glauert integral for n = 0, 1 and for each n; match terms with α − dz/dx; multiply by cos nθ and integrate over (0, π) to isolate Aₙ (orthogonality). (b) L′ = ρV∫γ dξ; M_LE = −ρV∫γ ξ dξ; transfer to c/4. (c) c_l = 0 ⇒ 2A₀ + A₁ = 0; substitute the integrals. (d) Integrate the piecewise slope; A₀ = α + δ(π − θ_h)/π, A₁ = (2δ/π) sin θ_h, A₂ = (δ/π) sin 2θ_h; c_l and c_m from (b); expand θ_h ≈ π − 2√E.

## 3.
1. `flap_cm_derivative(0.15)`; δ = −c_m,2412/(dc_m/dδ) (negative: trailing edge up). 2. `combine_slopes` of the 2412 line and `flap_camber_slope(0.15, δ)`; `thin_airfoil` for α_L0 and c_m. 3. α for c_l = 0.5 is 0.5/(2π) + α_L0 for each section. 4. Discuss the lost lift (higher α_L0) and the drag of a deflected tab.

## 4.
1. Double loop over m and p; `naca4_camber_slope`, `thin_airfoil(…, 0, breaks)`. 2. Ratios of columns. 3. `naca5_camber_slope(0.2025, 15.957)`; compare α_L0 and c_m with the 2412.

## 5.
1. Loop n; `lumped_vortex(slope, α, n, spacing)` for both spacings; errors against `thin_airfoil`; `observed_order`. 2. Locate x = 0.4 in each grid. 3. The parabola: `lumped_vortex` cm errors, order.

## 6.
Assumptions: small slopes (tan δ ≈ δ), attached flow over the hinge, thin section. Discuss viscous separation at the hinge as the first to fail. For (b), compare thickness effects on lift slope and on the flap's effectiveness, and add the moment change. For (c), weigh flap deflection, approach attitude and the stall margin.
