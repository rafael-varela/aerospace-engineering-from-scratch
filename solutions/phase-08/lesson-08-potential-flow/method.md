# P08.08  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. Γ = 0.3 × 4πUa. 2. `cylinder_stagnation_points(U, a, Gamma)` and their angles. 3. `cylinder_surface_cp` on a fine θ grid; max, min and their angles; Cp at 90° and 270°. 4. ρUΓ; `force_by_pressure(cylinder(U, a, Gamma), circle_contour(a, 128), 1.225, U)`; c_l = L′/(½ρU² 2a).

## 2.
(a) z = ae^{iθ}: W = 2Ua cos θ, real, so ψ = 0; log z on the circle is ln a + iθ, so the vortex adds a constant to ψ. (b) dW/dz on the surface; u_θ = −Im(e^{iθ} dW/dz)… or differentiate φ along θ; Cp = 1 − (u_θ/U)²; F = −∫₀^{2π} p (cos θ, sin θ) a dθ with the Bernoulli p; stagnation where u_θ = 0; merge when the two roots coincide at θ = 270°. (c) |dW/dz|² dz̄… use the surface as a streamline (dW real along it) to replace the conjugate; substitute the Laurent series, square, keep the 1/z term, apply ∮ dz/z = 2πi.

## 3.
(a) ψ = Ur sin θ + mθ/(2π); value at the stagnation point (θ = π); far downstream θ → 0 on the upper half: Uy → m/2. Stagnation where U + m/(2πx) = 0. (b) For each θ the surface radius r = m(π − θ)/(2πU sin θ); z = r e^{iθ}; `pressure_coefficient([uniform(U), source(m)], z, U)`; normal velocity from the tangent dz/dθ by differences. (c) argmin of Cp; the largest x where |Cp| exceeds the tolerance; convert to D from the nose. (d) Think about the stem, the fuselage and the angle of attack.

## 4.
1. L′ = mg/b; Γ = L′/(ρV); c_l = L′/(½ρV²c). 2. α = asin(c_l/(2π)); `cylinder_stagnation_points(V, c/4, Gamma)`. 3. `force_blasius(cylinder(V, c/4, Gamma), circle_contour(2a), rho)`; plate: `kutta_circulation`, `plate_complex_velocity` on `circle_contour(6a, 256)`, I = (iρ/2)∮(dW/dz)²dz, F_x = Re I, F_y = −Im I; rotate by α into lift and drag.

## 5.
For each strength: `rankine_oval_size`, `oval_contour(1, m, 1, 256)`, the maximum of |`complex_velocity`| on it, `force_by_pressure` on it. Compare the speed with 2 (cylinder) and 1 (no body).

## 6.
(a) Viscosity: the boundary layer and separation; the Kutta condition. (b) Γ = 2πa × 2U; c_l = Γ/(Ua). (c) Go through each prediction against the assumptions (inviscid, attached, thin).
