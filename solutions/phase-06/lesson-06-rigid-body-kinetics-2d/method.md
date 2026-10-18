# P06.06  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. I = m(0.7R)²; θ = 4 + s/2 degrees. 2. N = mg cos θ; f_roll = mg sin θ · I/(I + mR²); compare with μ_s N. 3. If it rolls, a from (3), α = a/R. 4. (b) is (4). 5. (c) tan θ_max = μ_s(1 + mR²/I). At 30°: sliding, a = g(sin θ − μ_k cos θ), α = μ_k N R/I, slip acceleration a − αR.

## 2.
(a) dv_s/dt = dV/dt − R dω/dt; substitute (5) with T_b = 0. (b) Integrate each equation of (5) over the slide; add the rolling condition; solve the 3 × 3 linear system for J, V_f, ω_f. (c) E_f = ∫F v_s dt, F dt = −M_eff dv_s, integrate from v_s0 to 0. (d) Expand ½M(V₀² − V_f²) = J(V₀ − J/(2M)) and ½I(ω_f² − ω₀²) = JR(ω₀ + JR/(2I)); subtract; use J = M_eff v_s0. (e) ∫N dt piecewise; set μ_k ∫N dt = J; during the ramp v_s = v_s0 − ct², integrate. (f) M_eff → I/R²; v_s0 = (1 − f)V₀.

## 3.
(a) M = m, N = mg, T_b = 0. (b) V_f = V₀ − J/m with J = M_eff(V₀ − ω₀R), M_eff = m k²/(1 + k²); V_f < 0 ⇔ k²|ω₀|R > V₀. (c) t_s = J/(μ_k m g); x = V₀ t_s − ½ μ_k g t_s². (d) ½M_eff v_s0² and KE₀ − KE_f; ratio to KE₀. (e) `Wheel(m_wheel=3, R=0.15, k_ratio=sqrt(0.8), M=3)`, `simulate_touchdown(w, V0, 0.3, 0.5, ramp_load(m g, 0), 3.0, omega0=-(30 + s))`.

## 4.
(a) M_eff from (6); J = M_eff V₀; compare J with ½ μ_k N_max t_r; (8); V_f, ω_f from (7); slip distance from the matching branch; E_f from (9); peak drag μ_k N(t_s). (b) `simulate_touchdown(Wheel(), V0, mu_k, 0.8, ramp_load(N, 0.08), 0.3)`; `energy_audit`. (c) Pass rtol and atol. (d) Plot `res.t`, `res.slip`; a parabola during the ramp, a straight line after it.

## 5.
(a) `max_brake_torque(mu_s, 120e3, I, M, R)`; compare with μ_s N R; the brake must also decelerate the wheel's spin. (b) Two runs; `res.transitions`; deceleration from the last two samples; lock time from the transition list against ω(0.3) I/(T_b − μ_k N R). (c) `res.E_brake[-1]`, `res.E_friction[-1]`, `energy_audit`. (d) d = V₀²/(2a) for each deceleration.

## 6.
(a) Tread energy per unit area: same energy over a longer slip. (b) (1 − f)² = 0.01; overspin reverses v_s0 and the friction pushes the aircraft forward. (c) Moment equation about an accelerating point that is not G needs the term r × m a_G. (d) A non-terminal or missing slip event: Coulomb friction reverses sign every step. (e) A load removed suddenly from a spring–mass system excites oscillation (spring-back). (f) μ_k and N(t) set the time; the energy and final state need neither.
