# P10.03  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. Γ = g · 4πaV. 2. `stagnation_angles(V, 0, a, Γ)`. 3. Cp(θ) = 1 − (2 sin θ + Γ/(2πaV))²; minimum at θ = 90°, maximum 1 at the stagnation points. 4. L′ = ρVΓ; c_l = L′/(q · 2a) = 2Γ/(V · 2a). 5. Pressure check with `pressure_force_coefficients(1j*z, cp, 0, 2a, dθ)` on z = a e^{iθ}. 6. Sketch: two stagnation points merge at the bottom when g = 1, a single one leaves the surface for g > 1.

## 2.
(a) On ζ′ = ae^{iθ}, the velocity in polar components is u_r − iu_θ = (dw/dζ) e^{iθ}; expand. (b) w_z finite requires w_ζ(b) = 0; set u_θ(−β) = 0. (c) L′ = ρVΓ with Γ from (b); divide by ½ρV²c; for a thin symmetric section a = b, β = 0, c = 4b. 4a/c is the thickness effect (a > c/4 for a thick section), β the camber effect (α_L0 = −β). (d) p − p∞ = ½ρV²(1 − (u_θ/V)²); force per unit span ∮ (p − p∞)(−cos θ, −sin θ) a dθ; use ∮ sin θ = 0, ∮ sin² θ = π.

## 3.
1. Γ = 2πaU_s; g = Γ/(4πaV) = U_s/(2V). 2. Compare g with 1; if g > 1 the single stagnation point is in the flow below (or above) the rotor. 3. L′ = ρVΓ; L = L′ × height; c_l = L′/(½ρV² · 2a). 4. Effects: boundary-layer separation on the side moving against the stream (drag, reduced circulation) and finite height (tip vortices, induced drag, `P10.09`). The weakest assumption is that the fluid circulates at the full surface speed.

## 4.
1. t from your seed; `mu = fit_karman_trefftz(t, 0.02, radians(16))`; `section_shape(1, mu, tau)`; `circle_geometry(1, mu)`. 2. α_L0,chord = −β + φ; slope 8πa cos(α + β)/c; α_x from arcsin(c_l c/(8πa)) − β; α_chord = α_x + φ. 3. `airfoil_flow` + `pressure_force_coefficients`; `contour_in_z("ellipse", 0.8c, ...)`, `airfoil_velocity_at`, then `circulation_line_integral`, `blasius_force`, `momentum_force`; normalise forces by ½V²c.

## 5.
1. Fit (t, m) for τ = 0 and τ = 16°; exact c_l from the formula. 2. Loop n; `airfoil_flow` + `pressure_force_coefficients`. 3. Orders log₂(e_n/e_2n). 4. Argument: the trapezoid rule on a periodic analytic integrand converges exponentially; a branch-point singularity of the integrand on the contour (r^{n−1} with non-integer power) limits it to an algebraic rate.

## 6.
1. Read the CSV (`#` comments), select Table I, the Re block, and −4.5° ≤ α ≤ 6.5°; least-squares slope in radians; divide by 2π. 2. Compare with 1.094. 3. Argue from the displacement thickness near the trailing edge (it decambers and effectively thins the section) and its decrease with Re. 4. Explain the Kutta condition as the steady-state outcome of the starting vortex shed by viscous action at the sharp edge.
