# P04.04  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. m = 2 + s/10. 2. `inertia_solid_cylinder(m, 0.1, 0.02)`: I₃ = ½mr², I₁ = m(3r² + ℓ²)/12. 3. ω = 6000 × 2π/60; h = I₃ω; T = ½I₃ω². 4. `parallel_axis(I, m, [0.5, 0, 0])`; d along a principal axis of the wheel ⇒ ddᵀ diagonal. 5. t = h/τ; body rate = h/135.4 by conservation of the axis-3 total.

## 2.
(a) dL/dt|i = dL/dt|b + ω × L; dL/dt|b = Iω̇. Write ω × (Iω) in components. (b) ω · (ω × Iω) = 0 ⇒ d(½ω·Iω)/dt = 0; (Iω) · (ω × Iω) = 0 ⇒ d|Iω|²/dt = 0. Inertial: dL_i/dt = R(dL/dt|b + ω × L) = Rτ = 0. (c) ω = (ε₁, ω₂ + ε₂, ε₃); drop products of ε; ε̈₁ = kε₁ with k = (I₂ − I₃)(I₁ − I₂)ω₂²/(I₁I₃); sign of k. (d) Body axes: ẑ = (0, sin θ, cos θ), ω = (0, Ω sin θ, ω₃), L = Iω; Ω ẑ × L = τ = mgl sin θ along axis 1 (sign as in the code). Quadratic in Ω; discriminant; the root 2c/(−b + √(b² − 4ac)) as ω₃ → ∞. (e) E = ½I₁(θ̇² + φ̇² sin²θ) + ½I₃ω₃² + mgl cos θ; L_z = I₁φ̇ sin²θ + L₃ cos θ = L₃ cos θ₀; set θ̇ = 0 and eliminate φ̇.

## 3.
(a) `inertia_solid_cylinder(m, 0.08, 2.0)`; ω₀ = 750 × 2π/60; L = I₃ω₀; T₀ = ½I₃ω₀². (b) Axis 3 has the smallest moment: an extreme axis. (c) T = L²/(2I) is smallest for the largest I; ω_f = L/I₁; T₀ − T_f. (d) State the major-axis rule.

## 4.
(a) `bang_bang_slew(135.4, radians(20 + s), 60)`. (b) Wheel torque −τ for t < 30 s, +τ after (the body feels the opposite); `simulate_rigid_body(I, 0, q_id, 60, wheel_torque=..., max_step=0.5)`; yaw = 2 atan2(q_z, q_w); body Iω + h. (c) Same with `h0=[hb, 0, 0]`; R from `quat_to_matrix`; tilt = arccos((R e₃)_z); max |ω₁|, |ω₂|; final |ω|. In (10), ω × h with h along axis 1 and ω along axis 3 gives a torque on axis 2. (d) Momentum: T ≥ 2IΘ/h_max; torque: T ≥ √(4IΘ/τ_max); check h at the torque-limited time.

## 5.
(a) ω₃ = 2π(40 + s); `steady_precession_rate`, `fast_top_precession`; `top_initial_state(th, w3, Om)`; simulate 1 s; `measured_precession_rate`. (b) ω_min = 2√(I₁mgl cos θ)/I₃. (c) `top_initial_state(th, w3, 0)`; simulate; `axis_tilt_and_azimuth`; `np.ptp(tilt)`; `nutation_turning_tilt(I1, I3*w3, m, l, th)`; peaks of the tilt; mean spacing ⇒ 2π/period. (d) Answer of (b) at 30°.

## 6.
(a) L accumulates as τt: t = 50/10⁻⁵. (b) L = Iω constant with I falling; T = L²/(2I) rises; her muscles do work. (c) Pairs of rotors spin opposite ways; speed one pair up, the other down. (d) The torque in the wrong frame, a sign error in ω × Iω, an unnormalised quaternion used for R, a loose solver tolerance. (e) The torque turns L; when released, the drop of the centre of mass feeds the precession and the nod. (f) Dissipation at fixed L drives the body to the major axis.
