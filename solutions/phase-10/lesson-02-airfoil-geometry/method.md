# P10.02  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. y_t from the polynomial with t = TT/100. 2. y_c = (m/p²)(2px − x²) and dy_c/dx = (2m/p²)(p − x) with m = 0.04, p = 0.4. 3. θ = arctan(dy_c/dx). 4. x_u = x − y_t sin θ, y_u = y_c + y_t cos θ; x_l = x + y_t sin θ, y_l = y_c − y_t cos θ. 5. r = 1.1019 t²; TE thickness 2y_t(1) = 0.021 t. 6. Call `naca4_thickness` and `naca4_camber` at x = 0.25 and compare.

## 2.
(a) ∫₀¹ (a₀√x + a₁x + a₂x² + a₃x³ + a₄x⁴)dx = 2a₀/3 + a₁/2 + a₂/3 + a₃/4 + a₄/5; multiply by 10t. Camber: compare `polygon_area` of 2412 and 0012. (b) y_t ≈ 5t a₀ √x; y_t² = 25t²a₀² x = 2rx. (c) Sum the coefficients with a₄ = −0.1036, then with −0.1015 and double. (d) d/dx of the cubic at x = p equals zero gives 3p² − 6rp + r²(3 − r) = 0; rearrange; solve by bisection. k₁ is a factor of y_c, so A₁ ∝ k₁ and k₁ = 0.3/(πA₁ for k₁ = 1).

## 3.
1. `naca_coordinates(code, 400)`, `polygon_area` × c², × span, × density. 2. `surface_thickness_camber`, scale by c, find where thickness ≥ 0.020 m. 3. `area_between(code, 0.2, 0.6, c)` × span; ratio to the total.

## 4.
1. Reference area with n = 20 000. 2. Loop n; cosine via `naca_coordinates`; uniform by building the stations with `np.linspace` and the same formulas. 3. `observed_order`. 4. Argument: the error of the first cell scales like h^{3/2} for y ∝ √x; substitute x = (1 − cos β)/2. 5. First station after the nose: h for uniform, (1 − cos(π/n))/2 for cosine; compare with r = 1.1019 t².

## 5.
1. Loop p over 0.05–0.25; `naca5_constants(p)`; `naca5_ideal_cl(*quoted)`. 2. For the worst pair, compute the k₁ that the quoted r would need (k₁ = 0.3/`naca5_ideal_cl(r_quoted, 1)`) to separate the r and k₁ contributions. 3. `naca5_constants(0.15, 0.6)`.

## 6.
1. `surface_thickness_camber` on each section; maxima and positions. 2. Same for 23012; compare the mean-line peak with p. 3. Build a vertical-thickness version (y_c ± y_t at the same x) and compare with the perpendicular one along the upper surface; locate the largest difference; relate it to θ.
