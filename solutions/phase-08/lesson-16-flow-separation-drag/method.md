# P08.16  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. `isa_troposphere(3000)`; `strut_drag(d, 0.3, 60, rho, mu)` gives Re, c_d, drag, power, shedding frequency. 2. `brentq` on `shedding_frequency(V, d, rho, mu) − 90` between a few m/s and 60 m/s. 3. Compare that speed with the take-off and approach speeds; options: stiffen, add mass, change shape, add a helical strake or a fairing.

## 2.
(a) n = (−cos θ, sin θ) for the point (−R cos θ, R sin θ); drag per unit span = −∮(p − p∞) n_x R dθ; divide by qd. Symmetry: substitute θ → π − θ in the integral (cos changes sign). (b) Front: ½∫(1 − 4 sin²θ) cos θ dθ from −θ_s to θ_s; wake: ½c_pb∫cos θ dθ from θ_s to 2π − θ_s. (c) 3πμVd = (ρ_p − ρ_f) g πd³/6; then multiply the general balance by ρ_f d²/μ² to get c_d Re² = (4/3)Ar.

## 3.
(a) `terminal_velocity` with each correlation. (b) `simulate_fall` to 6 s; first index with V ≥ 0.95V_t. (c) `numpy.interp` of V and t at z = 2 m; vacuum time √(2z/g); discuss how large the drag is compared with the weight at 2 m and how well the timing would resolve it.

## 4.
(a) D = ½ρV²c_d πD_c²/4 = mg ⇒ D_c; V ∝ ρ^(−1/2) at constant canopy. (b) `terminal_velocity(0.12, ρ_p, 1.225, 1.789e-5, cd=...)` for both; Re at the answers; check validity ranges; look at Morrison's curve near the answer. (c) Choose with the consequence in mind.

## 5.
(a) `strut_drag(..., t_over_c=...)` in a loop. (b) Separate 2C_f(c/t) and the form terms. (c) Difference of drags at t/c = 0.167 and 0.25.

## 6.
(a) List the model's omissions. (b) Pressure taps round the model and a manometer or scanner; the base pressure is one tap at θ = 180°. (c) Re of a golf ball at 60 m/s; where it sits relative to the smooth-sphere crisis; dimples trip the boundary layer; then the strut: compare a fairing's saving with a crisis's.
