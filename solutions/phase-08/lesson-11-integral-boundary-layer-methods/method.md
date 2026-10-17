# P08.11  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. a = 2U/r. 2. θ = √(0.075ν/a); λ = 0.075; (l, H) = `thwaites_l_H(0.075)`; δ* = Hθ. 3. `falkner_skan(1.0)`: θ_η, δ*_η, f″(0); multiply by √(ν/a). 4. τ_w = μUe l/θ (Thwaites) and μ a x f″(0)√(a/ν) (exact) at x = 2 mm; sin(x/r)/(x/r) − 1.

## 2.
(a) Integrate each term from 0 to h; integrate v u_y by parts; continuity; collect into d/dx(Ue²θ) + Ue Ue′ δ* = τ_w/ρ; divide by Ue². (b) Multiply by 2Ueθ/ν; identify l and λ; replace the right side by 0.45 − 6λ; recognise d(θ²Ue⁶)/dx; integrate. (c) Substitute the two Ue(x) into the quadrature; solve λ = −0.09 for x. (d) Ue = 2U sin φ, dUe/dx = (2U/R)cos φ, dx = R dφ; λ contains only φ.

## 3.
(a) Quadrature with θ₀ = 0; λ = θ²Ue′/ν; solve λ = −0.09; fall = 1 − Ue/U₀. (b) θ₀ from the flat plate; `thwaites(ξ, Ue, Ue′, ν, theta0=θ₀)`; `separation_point`. (c) Composite Ue on one grid from 10⁻⁴ m; `march_boundary_layer`; `march_separation`; subtract l₀. (d) Shorter straight section, gentler angle, suction or a trip.

## 4.
1. For each m: β, `falkner_skan(β)`, `wedge_thwaites(m)`. 2. Convert both to θ√Re_x/x, c_f√Re_x, H. 3. Errors in per cent; plot. 4. Compare dλ/dx near separation for Howarth with the gap between −0.068 and −0.09.

## 5.
(a) Loop over nx; march; `march_separation`; extrapolate (differences shrink slowly, so estimate from the trend; compare with 0.1199). (b) Same for the cylinder; degrees = x_s/R in degrees. (c) The PDE is parabolic only while u > 0.

## 6.
(a) Compare λ just behind each peak. (b) Free shear layer → transition → reattachment: a laminar separation bubble; drag and, if it bursts, lift. (c) The boundary-layer approximation (thin layer, pressure imposed): measure the surface pressure distribution and look for a plateau.
