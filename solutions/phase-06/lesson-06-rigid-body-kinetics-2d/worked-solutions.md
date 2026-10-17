# P06.06  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for **s = 15**; numerical results are from `solution.py` (and the probe scripts that call it) and are given to the figures shown. All "simulated" values are outputs of the model, not measurements.

## 1. A wheel on a ramp (s = 15: θ = 11.5°)

**Assumptions.** Rigid wheel; Coulomb friction; g = 9.80665 m/s².

(a) I = 8 × (0.175)² = 0.2450 kg·m²; N = mg cos θ = 76.878 N. The friction needed to roll is f = mg sin θ · I/(I + mR²) = **5.1437 N**, below μ_s N = 11.53 N, so **it rolls**: a = g sin θ/(1 + I/(mR²)) = **1.3122 m/s²**, α = a/R = **5.2487 rad/s²**.

(b) μ_min = tan θ/(1 + mR²/I) = **0.0669**.

(c) tan θ_max = 0.15 × (1 + 1/0.49) = 0.456: **θ_max = 24.52°**. At 30° it slides: a = g(sin 30° − 0.12 cos 30°) = **3.8842 m/s²**, α = μ_k N R/I = **8.3195 rad/s²**; the slip velocity grows at a − αR = **1.8043 m/s²**.

**Check.** At 30° the friction needed to roll would be 12.90 N, above μ_s N = 10.19 N, consistent with sliding. Dimensions: α = (N · m)/(kg · m²) = s⁻².

**Common mistakes.** Using μ_k to decide whether it rolls (the test is static friction); taking I = ½mR² for a wheel whose radius of gyration is given.

## 2. The spin-up in closed form

(a) v_s = V − ωR, so dv_s/dt = −F/M − R(F R)/I = −F(1/M + R²/I) = **−F/M_eff, 1/M_eff = 1/M + R²/I**.

(b) Over the slide, ∫M dV = −∫F dt gives **M(V₀ − V_f) = J**; ∫I dω = ∫F R dt gives **I(ω_f − ω₀) = JR**; at the end **V_f = ω_fR**. Substituting the first two into the third: V₀ − J/M = ω₀R + JR²/I, so **J = (V₀ − ω₀R)/(1/M + R²/I) = M_eff v_s0**. μ_k does not appear because these are three momentum statements and one kinematic condition: how the impulse is delivered in time does not matter, only its total, which the rolling condition fixes.

(c) E_f = ∫F v_s dt. From (a), F dt = −M_eff dv_s, so E_f = −M_eff ∫_{v_s0}^{0} v_s dv_s = **½ M_eff v_s0²**.

(d) ½M(V₀² − V_f²) = ½M(V₀ − V_f)(V₀ + V_f) = ½J(2V₀ − J/M) = J V₀ − J²/(2M). ½I(ω_f² − ω₀²) = ½JR(2ω₀ + JR/I) = J ω₀R + J²R²/(2I). The loss is the first minus the second: J(V₀ − ω₀R) − ½J²(1/M + R²/I) = J v_s0 − J²/(2M_eff). With J = M_eff v_s0: M_eff v_s0² − ½M_eff v_s0² = **½ M_eff v_s0²**, which is (9).

(e) With N = N_max t/t_r for t ≤ t_r: μ_k ∫₀^{t_s} N dt = μ_k N_max t_s²/(2t_r) = J gives **t_s = √(2Jt_r/(μ_kN_max))**, valid if t_s ≤ t_r, i.e. J ≤ ½μ_kN_max t_r. Otherwise μ_kN_max(t_s − t_r/2) = J gives **t_s = J/(μ_kN_max) + t_r/2**. In the first branch dv_s/dt = −μ_kN_max t/(M_eff t_r), so v_s = v_s0 − ct² with c = μ_kN_max/(2M_eff t_r), and s = v_s0 t_s − c t_s³/3.

(f) M → ∞: M_eff → I/R², V_f → V₀, and E_f → ½(I/R²)(V₀ − ω₀R)² = ½I(ω_f − ω₀)², which for ω₀ = 0 is **½Iω_f²**, the wheel's final kinetic energy. Pre-spin ω₀ = fV₀/R gives v_s0 = (1 − f)V₀ and **E_f ∝ (1 − f)²**.

**Checks.** At the branch point both expressions in (e) give t_s = t_r. In (d) the loss is positive for any v_s0 ≠ 0, as dissipation must be.

**Common mistakes.** Conserving kinetic energy across the slide; writing J = M V₀ (forgetting that the wheel must spin up, not the aircraft stop); dropping the R² in 1/M_eff.

## 3. A canister thrown with backspin (s = 15: V₀ = 3.5 m/s, ω₀ = −45 rad/s)

**Assumptions.** The canister rolls on its curved surface on a level floor; Coulomb friction; I = 0.8mR² = 0.0540 kg·m².

(a) The canister is both the "axle" and the "wheel": **M = m = 3 kg** in the translation equation, **N = mg = 29.42 N**, and T_b = 0. Equation (5) is m dV/dt = −F, I dω/dt = FR.

(b) M_eff = m k²/(1 + k²) and J = M_eff(V₀ − ω₀R), so V_f = V₀ − J/m = **(V₀ + k²ω₀R)/(1 + k²)**. It returns (V_f < 0) if **k²|ω₀|R > V₀** for backspin ω₀ < 0. Here k²|ω₀|R = 5.40 m/s > 3.5 m/s: **yes, it comes back**.

(c) v_s0 = 3.5 + 45 × 0.15 = 10.25 m/s; M_eff = 1.3333 kg; J = **13.667 N·s**; **V_f = −1.0556 m/s** (towards you), ω_f = −7.037 rad/s (rolling backwards, = V_f/R). The slide lasts t_s = J/(μ_k m g) = **1.5485 s**. The centre decelerates at μ_k g = 2.942 m/s² throughout the slide, so it stops 2.082 m out at 1.19 s, comes back, and starts to roll at x = V₀t_s − ½μ_kg t_s² = **1.893 m** from the release point.

(d) E_f = ½ M_eff v_s0² = **70.04 J**. KE₀ = ½ × 3 × 3.5² + ½ × 0.054 × 45² = 73.05 J; KE_f = 1.671 + 1.337 = 3.008 J; loss **70.04 J**, the same. The fraction is **0.959**: nearly all the energy of a backspin throw goes into heat.

(e) `spin_up_closed_form(3.5, 0.15, 0.054, 3, 0.3, 29.42, 0, omega0=-45)` gives the same J, t_s, V_f and E_f; the simulation gives t_s = 1.548458 s, V_f = −1.05556 m/s and the centre at 1.8926 m when rolling starts, equal to the hand values to the digits shown.

**Check.** Angular momentum about any fixed point on the floor is conserved during the slide (friction and N act through floor-level points, and N balances weight): m V₀R + Iω₀ = 3 × 3.5 × 0.15 − 0.054 × 45 = −0.855 kg·m²/s, and after: mV_fR + Iω_f = −0.475 − 0.380 = −0.855 kg·m²/s. ✓ This is a third route to V_f.

**Common mistakes.** Using the backspin's sign wrongly (v_s0 = V₀ + |ω₀|R); taking M_eff = I/R² as if the canister were attached to an aircraft.

## 4. Your own touchdown (s = 15: V₀ = 68 m/s, μ_k = 0.50, N_max = 115 kN, t_r = 0.08 s)

**Assumptions.** The lesson's wheel (I = 15.90 kg·m², R = 0.56 m, M = 15 000 kg, M_eff = 50.529 kg).

(a) J = M_eff V₀ = **3436.0 N·s**. ½μ_kN_max t_r = 2300 N·s < J, so the **second branch** of (8): t_s = 3436.0/57 500 + 0.04 = **99.76 ms**. V_f = 68 − 3436.0/15 000 = **67.7709 m/s**; ω_f = JR/I = **121.020 rad/s**. Slip distance: during the ramp v_s = v_s0 − ct² with c = 57 500/(2 × 50.529 × 0.08) = 7112 m/s³; s₁ = 68 × 0.08 − 7112 × 0.08³/3 = 4.226 m; v_s(t_r) = 22.48 m/s; then a linear fall over 19.76 ms adds 0.222 m: **s = 4.448 m**. E_f = ½M_eff V₀² = **116.82 kJ**; peak drag μ_kN_max = **57.5 kN**.

(b) At rtol = atol = 10⁻¹⁰: t_s differs by **7.4 × 10⁻¹² s**, E_f by **2.1 × 10⁻⁵ J**, the slip distance by **5.3 × 10⁻¹⁰ m**; the audit residual is **6.0 × 10⁻¹³**.

(c) At 10⁻⁶: **−3.4 × 10⁻⁷ s, −0.44 J, −1.5 × 10⁻⁵ m, audit 1.6 × 10⁻⁸**. The friction energy integral degrades most in absolute terms, yet it is still right to 4 parts per million; the kinetic-energy route is unaffected because the final state is set by the event.

(d) The slip velocity is a downward parabola (v_s0 − ct²) during the ramp, because the friction grows linearly; after 80 ms it falls along a straight line, because the friction is constant; it reaches zero at 99.76 ms with a kink (the event), after which it stays at zero.

**Check.** The wheel (and so M_eff) is the lesson's, and E_f depends on nothing but M_eff and v_s0, so E_f must be the lesson's 123.797 kJ times (68/70)² = 0.94367: 116.82 kJ. ✓ The different μ_k and load do not enter.

**Common mistakes.** Using the first branch without checking J against ½μ_kN_max t_r; forgetting the ramp's t_r/2 in the second branch.

## 5. Braking at the skid limit (s = 15: μ_s = 0.85)

**Assumptions.** As in the lesson; the brake torque is constant once applied.

(a) (11): **57.31 kN·m**; μ_sNR = **57.12 kN·m**. (11) is larger because part of the brake torque goes into decelerating the wheel's own rotation (the I term in I + MR²), which does not need runway friction.

(b) At 0.95 ×: phases **sliding → rolling** (spin-up at 103.59 ms, no skid); deceleration **6.460 m/s²** (= 0.95μ_sN/M). At 1.05 ×: **sliding → rolling → sliding (0.300 s) → locked (0.3853 s)**; the wheel locks **85.3 ms** after the skid starts, matching ω I/(T_b − μ_kNR) for the wheel's rate at 0.3 s; then the deceleration is **4.400 m/s²** (μ_kN/M).

(c) At 0.95 ×: brake **2006.6 kJ**, friction heat **123.8 kJ** (the spin-up only), audit residual 4.5 × 10⁻¹⁴. At 1.05 ×: brake **319.8 kJ**, friction heat **1295.7 kJ**, residual 4.5 × 10⁻¹⁴. Rolling braking puts the energy into the brake's heat sink, which is designed for it; a locked wheel puts it into the tyre and the runway.

(d) 70²/(2a): **379 m** rolling at 0.95 × against **557 m** locked. Exceeding the limit by 10 % costs 47 % more distance and a flat-spotted tyre; anti-skid aims to hold the wheel just below the static limit, which it does not know in advance, so it senses the incipient skid (the wheel decelerating faster than the aircraft) and releases.

**Check.** Energy audit closed to 10⁻¹³ in both runs: brake + friction heat + remaining kinetic energy = initial kinetic energy.

**Common mistakes.** Expecting a skidding wheel to brake harder; forgetting that once locked the brake does no work (ω = 0).

## 6. Judgment

(a) No. The energy is the same, but on a wet runway (low μ) it is spread over a longer slip distance (Experiment 1: 6.6 m at μ = 0.3 against 4.0 m at 0.8), so each part of the tread is heated less; on a dry runway the heat is more concentrated, and wear also depends on the tread temperature reached. The energy result alone does not settle the wear question.

(b) At 90 % pre-spin, (1 − f)² = 0.01: the heat falls from 124 kJ to 1.2 kJ. If the wheel is spun faster than V₀/R, v_s0 < 0, friction acts forward on the axle and decelerates the wheel: the energy is again ½M_eff v_s0², now from the wheel's excess spin, and the gear sees a forward load. The benefit depends on matching the speed, which changes with every approach.

(c) The contact point is not the centre of mass, and here the axle and the contact point accelerate. ΣM = I α holds about G (or a fixed point, or the instantaneous contact point only under special conditions); about any other point the moment of m a_G must be added. Taken about the contact, the equation silently drops the friction torque that spins the wheel.

(d) Coulomb friction was left switched on across v_s = 0: friction reverses with the slip's sign, so the integrator chatters around rolling, with the rim speed oscillating about the axle speed and steps collapsing. The fix is a terminal event at v_s = 0 and a switch to the rolling equations, as in `simulate_touchdown`.

(e) The leg is a spring with a mass at its end. A drag load that builds over 0.1 s and is then removed in an instant (at spin-up) leaves the leg deflected and lets it spring forward at its natural frequency; the peak loads in the leg and its attachment can exceed the static drag. It is a design case for gear structures.

(f) The friction coefficient and the normal-load history set the spin-up time, slip distance and peak drag, so they must be measured (or bounded) before those outputs are trusted. The final rolling state and the friction energy follow from momentum alone and can be trusted without them.

## Check-yourself answers

1. ΣF = m a_G and ΣM_G = I_G α. About the centre of mass G, because there it has no extra term even when G accelerates; about any other moving point an m a_G moment term appears.
2. The disk: its accelerations are g sin θ/(1 + ½) and g sin θ/(1 + 1), a ratio of 4/3 in favour of the disk.
3. 1/M_eff = 1/M + R²/I; the mass a friction force at the rim "feels" when it changes the slip velocity. For a very heavy aircraft M_eff → I/R².
4. Because the final state follows from the linear and angular impulses (both equal to J) and the rolling condition; μ_k sets only how fast J is delivered.
5. E_f = ½M_eff v_s0². For M → ∞ the aircraft loses 2 × ½Iω_f², so the heat is one half of its lost kinetic energy.
6. (1 − 0.7)² = 0.09: the energy falls to 9 %, by a factor of about 11.
7. Rolling, the tyre can supply up to μ_s N; locked, it supplies only μ_k N < μ_s N.
8. Because Coulomb friction reverses at v_s = 0; without the event the integrator steps across zero repeatedly and chatters, and the rolling constraint, which needs a different equation (friction as an unknown), is never applied.
