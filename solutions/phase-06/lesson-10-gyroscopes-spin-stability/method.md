# P06.10  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. ω₃ = 0.5 + s/20. 2. H = √((I₃ω₃)² + (I_tω_t)²); λ = (I₃ − I_t)ω₃/I_t; Ω_n = H/I_t; θ = atan(I_tω_t/(I₃ω₃)). 3. Periods 2π/λ and 2π/Ω_n. 4. `simulate_rigid(I, [0.02, 0, w3], 40)`; `axis_coning_rate(res, H_i[0])`.

## 2.
(a) H along the axis: I₃ω₃ = H cos θ; across: I_tω_t = H sin θ; T = ½(I₃ω₃² + I_tω_t²). For a general body T = Σ Hₖ²/(2Iₖ) with Σ Hₖ² = H². (b) Differentiate (6) with respect to θ, set dT/dt = −P. (c) Sum: d(Iω + Jω_d)/dt + ω × (Iω + Jω_d) = 0, then transport to inertial axes; dot the first equation with ω, the second with ω_d, add. (d) ω = (z, Ω), ω_d = (z_d, Ω); keep first order; ω × Iω and ω × ω_d in complex form. (e) Add the two rows of (9) weighted by I_t and J; for large c, z_d → z.

## 3.
(a) τ = 60 × 10³ × 10⁻³; ω₃ = 2π rad/s; H = I₃ω₃; (5). (b) `simulate_rigid(diag(I_t, I_t, I₃), [0, 0, ω₃], 80, torque=lambda t, q, w: [τ, 0, 0])`; angle of H_i from H_i(0); angle of `spin_axis_inertial` from z. (c) ½(τ/I_t)t². (d) e^{80/300}; compare burn and coast times.

## 4.
(a) (4) with θ = 20 + s degrees; period 2π/Ω. (b) q₀ = (cos(θ/2), sin(θ/2), 0, 0); `simulate_rigid(I, [0, 0, 1], 100, torque=gyroscope_torque(I, tau), q0=q0)`; `inertial_momentum`; slope of the unwrapped atan2(H_y, H_x); |H| and the tilt from H_z/|H|. (c) τ/H.

## 5.
(a) Loop over c; optimise over log₁₀ c. (b) `simulate_with_damper(diag(1500, 1500, 600), 60, c, [1e-3, 0, 3], [0, 0, 3], 900, n_samples=9001)`; `fit_growth_rate(t, |ω_t|, 1e-2, 0.2)`; first times nutation exceeds 45° and 89°; `damper_energy` and H from `inertial_momentum(…, J)`. (c) ln(200)/σ. (d) The c that maximises the decay rate for the oblate body (repeat (a) with the oblate inertias).

## 6.
(a) Rigid equations are right for a rigid body; the satellite was not rigid. (b) Major-axis rule. (c) Burn short compared with the growth time. (d) The sign of (7). (e) (7) and the end state need only energy and momentum; the rate needs the sink. (f) Uniform cylinder and an illustrative damper.
