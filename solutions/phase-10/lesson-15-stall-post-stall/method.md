# P10.15  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. `kirchhoff_cl(np.radians(10), f, 5.741)` for three f. 2. `flat_plate(np.radians(45))`; resultant √(cl² + cd²) = 2 sin α, normal to the chord. 3. Compare (1) at 60° with the plate's.

## 2.
(a) f = 1 ⇒ factor 1; f = 0 ⇒ ¼. (b) Minimise ‖Xp − y‖²: XᵀXp = Xᵀy; Cov(p) = (XᵀX)⁻¹Xᵀ Cov(y) X(XᵀX)⁻¹ = σ²(XᵀX)⁻¹; linearise the nonlinear model about the solution. (c) At convergence δ = α_eff − cl_sec/a₀, so a₀(α_eff − δ) = cl_sec: the linear solver's circulation equals the section's. (d) Error recursion and eigenvalues; with λ < 0 real, 1 + ω(λ − 1) > −1 ⇔ ω < 2/(1 − λ).

## 3.
1. q = ½ × 1.225 × 20². 2. `flat_plate` at 35°, 45°, 60°; L = qS cl, D = qS cd, N = qS·2 sin α. 3. Lift supports the aircraft, drag decelerates it; argue about attached flow from the slipstream.

## 4.
1. Loop over washout; `stall_onset_alpha(chord, twist, curve, a0, a_stall, tol_deg=0.01)`. 2. Find the smallest washout with η < 0.3. 3. `lifting_line_linear(np.radians(5), ...)["e"]`; note that a twisted wing's load shape changes with α.

## 5.
1. Fit, `section_stall_angle`. 2. For each planform, `stall_onset_alpha(..., lo_deg=0, hi_deg=35)` and `wing_cl_max` over 0–30°. 3. Compare stations and C_L/cl_max; compute the trainer's Re at its highest speed.

## 6.
1. Up-sweep then down-sweep with `lifting_line_nonlinear(..., n_terms=N, delta0=previous)`. 2. Uniqueness of solutions with a monotone section curve; resolution of stall cells. 3. State what the model can and cannot claim. 4. Spanwise and unsteady physics.
