# P04.03  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. m = 8000 + 100s, v = 0.05 + s/1000. 2. `collide_1d(m, v, 420e3, 0.0, 0.0)`; `collision_energy_loss(..., 0.0)`. 3. Same with e = 0.2; the vehicle's velocity is negative: it rebounds. 4. ΔK/K₀ = μ/m = M/(m + M) > 0.97 ⇔ M > 32.3 m.

## 2.
(a) Before: mv. After: (m − dm)(v + dv) + dm(v + dv − u). Solve for dv. With u relative to the cart before the throw: pellet at v − u, cart at v + dv′: dv′ = u dm/(m − dm); compare by expanding 1/(m − dm). (b) v_n = uΣΔm/mᵢ with mᵢ the left endpoints of intervals of the decreasing variable m; trapezoid = (left + right)/2 is exact to O(Δm²) per interval; left − trapezoid = (Δm/2)Σ(f(mᵢ) − f(mᵢ₊₁)) telescopes; multiply by u. At n = 3, compare `pellet_convergence(..., [3])` with the formula; the remainder is O(1/n²). (c) Leak: dm leaves at velocity v (zero relative speed): the rocket's velocity cannot change (the pellet balance gives dv = 0). d(mv)/dt = m dv/dt + v dm/dt = 0 gives dv/dt = −(v/m)dm/dt ≠ 0, frame-dependent through v. Correct: apply (4) over dt with ṁ, add F_ext dt. (d) z = ∫v dt with v from (9); substitution m = m₀ − ṁt; ∫ln(m₀/m)dm = m ln(m₀/m) + m.

## 3.
(a) All at once: (120)v + 3m_t(v − 5) = 0. One at a time: `pellet_rocket_velocities(120 + 3m_t, 120, 5, 3)[-1]`. Times = 20/v. (b) `tsiolkovsky_dv(5, 120 + 3m_t, 120)`; Tsiolkovsky is the n → ∞ limit of the increasing sequence v_n. (c) Initial frame: 120v = 3m_t × 5. (d) ½ × 120 × v² against 30 J.

## 4.
(a) For each tw: `mdot = tw * 10000 * G0 / u`; `fl = simulate_rocket(..., k_drag=k)`; `velocity_losses(fl, u)`; `q = 0.5 * 1.225 * np.exp(-fl.z / 8500) * fl.v ** 2`; `q.max()`. (b) `minimize_scalar(lambda tw: -simulate_rocket(...).v[-1], bounds=(1.2, 10), method="bounded")`. (c) `brentq(lambda tw: maxq(tw) - 40e3, 1.2, 4)`; burnout there; difference from (b). (d) Think about engine mass growing with thrust, acceleration limits on crew and payload, and the fact that drag grows with C_D at transonic speed.

## 5.
(a) Stage list `[(f*mp, structure_mass(f*mp, 0.08), 2800), ((1-f)*mp, structure_mass((1-f)*mp, 0.12), 3400)]`; `minimize_scalar` on −`stack_dv(...)[0]`. (b) u ln(1/ε) for each. (c) Wrap (a) in a function of the payload; `brentq(lambda pl: best(pl) - 9400, 1, 5000)`; gross mass = payload + propellant + structures. (d) In (11), stage 2's mass ratio is limited by its own ε; moving propellant to it also increases the dead weight stage 1 must carry.

## 6.
(a) m_p = M(1 − e^{−Δv/u}) for each u; think about thrust level (time) and power. (b) m changes during the burn: use u ln R; impulse/m with m = initial mass underestimates. (c) e = v_up/v_down = √(h_up/h_down). (d) Momentum is carried away if the vent has a net direction; opposite nozzles cancel. (e) Orbital speed at 127 km is about 7.8 km/s horizontal. (f) The thrust term, the exhaust-velocity term (v − u vs −u), or the d(mv)/dt form.
