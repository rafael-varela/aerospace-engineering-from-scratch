# P01.11  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. Velocity from equation (1): magnitude Γ/(2πr), direction (−y, x)/r. 2. Circulation = Γ if the loop encloses the origin, else 0; flux of a vortex is 0 on every loop. 3. Build the field as a lambda, call `circulation(field, circle_path(...))` and `flux_2d`, subtract the exact values. 4. Compare the circle error with Γ[(n/π) tan(π/n) − 1] ≈ Γπ²/(3n²).

## 2.
(a) Write θ = atan2(y, x); for x > 0, θ = arctan(y/x); differentiate with the chain rule and simplify with r² = x² + y²; argue the formulas hold on the other branches too. (b) On a loop that does not enclose the origin choose a ray from the origin that misses the loop; define θ continuously on the plane minus that ray; apply φ(B) − φ(A) with B = A. (c) Either follow θ: one counter-clockwise turn adds 2π to θ and Γ to φ; or deform the loop onto a small circle (the region between them contains no singularity, so the two circulations are equal by (b) applied to the boundary of that region) and use the §2.2 calculation. (d) Add the contributions of the lobes with their signs; winding number times Γ.

## 3.
(a) Γ = 2πaV_max; ω = 2V_max/a (vorticity = 2 × angular velocity of the solid-body core). (b) Inside: Γ(R/a)²; outside: Γ; speed at r > a: Γ/(2πr). (c) Form the two ratios. (d) `rankine_vortex_velocity` in a lambda, `circulation` with `circle_path(R, 2000)`; units are metres, so use 15 000 and 60 000.

## 4.
(a) `errs = circulation_convergence(field, 1.0, [64, 128, 256, 512])`; rates `np.log2(errs[:-1] / errs[1:])`; compare column 0 with the closed form. (b) Chord midpoint at radius R cos(π/n); speed there Γ/(2πR cos(π/n)); chord length 2R sin(π/n), parallel to the velocity; per-segment product; multiply by n; expand tan. (c) Field from `cylinder_flow_velocity` in a lambda; `rectangle_path(-2, 2, -2, 2, n // 4)`; relative error |Γ_n − Γ|/|Γ|; solve C/n² = 10⁻⁴ for n with C from the n = 64 point. (d) The doublet velocity at diametrically opposite midpoints is equal and the chords are opposite, so the doublet's contributions cancel in pairs; the vortex part is (9). Move the centre and the cancellation no longer applies to the doublet *exactly*, but the doublet is a gradient field with a single-valued potential, so its polygon circulation is still tiny; report what you observe.

## 5.
(a) ṁ = 2πρV∫₀^R (1 − r/R)^{1/7} r dr; substitute u = 1 − r/R, r = R(1 − u), dr = −R du; the integral becomes R²∫₀¹ u^{1/7}(1 − u) du = R²(7/8 − 7/15). (b) Build the field as a lambda returning (0, 0, ρV(1 − r/R)^{1/7}); loop over n; `order = log2(Δ_k / Δ_{k+1})` for successive differences. (c) Define a new surface `disc_w(w, t) = (R(1 − w⁷) cos t, R(1 − w⁷) sin t, 0)` on w ∈ [0, 1]; the normal flips sign (w runs inward), so take the absolute value or swap the parameter order; repeat (b).

## 6.
For each: identify the assumption; think of the instrument (PIV, pitot rake, tracking, anemometer); give the quantity it would show changing.

## 7.
(a) Γ = 4W/(πρVb); b′ = πb/4. (b) Velocity induced at one vortex by the other: Γ/(2πb′), direction from the right-hand rule; both move the same way. (c) Superpose two `vortex_velocity` calls with opposite signs and the two centres; loops: circles of radius 5 m around each centre and a rectangle around both. (d) Both vortices induce downward velocity Γ/(2π · b′/2) at the midpoint; sum; compare with (b). (e) Swirl at 5 m: Γ/(2π · 5 m); wingtips at ±5 m see opposite vertical velocities; angle-of-attack change atan(w/V); compare with a stall angle of about 15°.
