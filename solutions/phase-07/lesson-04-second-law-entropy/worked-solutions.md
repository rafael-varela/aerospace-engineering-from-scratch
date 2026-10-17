# P07.04  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numbers come from `solution.py` and the probe scripts that call it. Air: R = 287.045 J/(kg·K), NASA variable cp; T₀ = 288.15 K.

## 1. Three generators (s = 15)

(a) Q̇ = 250 W: Ṡ_gen = 250(1/290 − 1/340) = **0.1268 W/K**.

(b) T_f = (330 + 300)/2 = **315 K**; S_gen = 900 ln(315/330) + 900 ln(315/300) = **2.0431 J/K**. Small-ΔT formula: 900 × 30²/(4 × 315²) = **2.0408 J/K**, 0.11 % low, the size of the neglected fourth-order term.

(c) p₁ = 2.5 MPa: s_gen = 287.045 ln 25 = **923.96 J/(kg·K)**; lost work T₀ s_gen = **266.2 kJ/kg**.

**Check.** Each is positive, and each vanishes when its driving difference vanishes.

**Common mistakes.** Using °C in 1/T; putting the cold temperature in the negative term in (a); computing the throttle's T change from cp (an ideal gas's T does not change).

## 2. Derivations

(a) Go from 1 to 2 by the actual process and return by a reversible one. Clausius: ∫₁²δQ/T + ∫₂¹δQ_rev/T ≤ 0, and ∫₂¹δQ_rev/T = S₁ − S₂, so ∫₁²δQ/T ≤ S₂ − S₁. Define S_gen = S₂ − S₁ − ∫₁²δQ/T ≥ 0. For a control volume, mass carries entropy ṁs across the boundary: dS_cv/dt = ΣQ̇_k/T_k + Σ_in ṁs − Σ_out ṁs + Ṡ_gen; in steady flow the left side is zero.

(b) In a reversible process δq = T ds and δw = p dv, so du = T ds − p dv: **T ds = du + p dv**. With dh = du + p dv + v dp: **T ds = dh − v dp**. They relate properties, so they hold between any two states. Ideal gas: dh = cp(T) dT, v/T = R/p, so **ds = cp dT/T − R dp/p** and s₂ − s₁ = ∫cp dT/T − R ln(p₂/p₁).

(c) S_gen = −C₁ ln(T₁/T_f) − C₂ ln(T₂/T_f). Since ln y ≤ y − 1, −ln(T_i/T_f) ≥ 1 − T_i/T_f, so S_gen ≥ C₁ + C₂ − (C₁T₁ + C₂T₂)/T_f = 0, because C₁T₁ + C₂T₂ = (C₁ + C₂)T_f by the energy balance. Equality only if T₁ = T₂ = T_f.

(d) S_gen = Q_C/T_C − Q_H/T_H with Q_C = Q_H − W = Q_H(1 − η): S_gen ≥ 0 ⇔ (1 − η)/T_C ≥ 1/T_H ⇔ **η ≤ 1 − T_C/T_H**.

**Common mistakes.** Writing ∫δQ/T = ΔS for an irreversible process; evaluating δQ/T at the system's temperature instead of the boundary's.

## 3. A satellite (W = 200 W)

(a) Dissipation of work at 350 K: 200/350 = **0.5714 W/K**. Conduction 350 → 290 K: 200(1/290 − 1/350) = **0.1182 W/K**. Radiation 290 → 3 K: 200(1/3 − 1/290) = **65.977 W/K**. Total **66.667 W/K**.

(b) The sum telescopes: W/350 + (W/290 − W/350) + (W/3 − W/290) = **W/3**. Every joule of work ends up as heat in the 3 K reservoir; how it gets there does not change the total, only its distribution among the steps.

(c) Radiation to space dominates (99 %). That does not make the radiator bad: no radiator can reject heat near 3 K, and a radiator's job is to get rid of heat, so a large temperature drop to space is the price of a compact, light radiator (its area scales as 1/T⁴, lesson 14). It shows that the "lost work" ruler needs a sensible T₀. Against T₀ = 3 K almost all of the 200 W is "lost work", which is true but not useful; for designing the electronics' cooling, the meaningful irreversibility is the 0.118 W/K between chip and radiator.

**Common mistakes.** Forgetting the first step (work dissipated into heat generates entropy W/T even with no temperature difference).

## 4. Entropy is a property (s = 15: T₂ = 650 K, p₂ = 5 bar)

(a) Δs = **329.9897 J/(kg·K)**.

(b) Isobaric then isothermal: **329.98968**; isothermal then isobaric: **329.98968**; the wiggly path at n = 20 001 agrees with (a) to the trapezoid error (its error at n = 801 is already 4.4 × 10⁻⁴).

(c) Errors: **2.785 × 10⁻², 6.967 × 10⁻³, 1.742 × 10⁻³, 4.355 × 10⁻⁴** at n = 101, 201, 401, 801; observed orders **1.999, 2.000, 2.000**.

(d) cp = 1005: **315.08 J/(kg·K)** (−4.5 %). Mean cp = Δh/ΔT: **332.84 J/(kg·K)** (+0.86 %). The mean cp is better because it is the right average for h, though Δs weights cp by 1/T, so even it is not exact: the cold end, where cp is smaller, counts more in ∫cp dT/T than in ∫cp dT.

**Check.** The two piecewise paths agree with each other to 10⁻¹³ and with (a) to 2.8 × 10⁻⁶, the trapezoid error of the isobaric segment.

**Common mistakes.** Integrating R dp/p with p in bar and Pa mixed (the ratio is what matters, but both ends must use the same unit); thinking a path with a loop "adds" entropy (it does not: Δs is a property).

## 5. Redesigning the bleed chain (s = 15: p_prv = 0.35 MPa, precooler outlet 465 K)

(a)

| process | Ṡ_gen (W/K) | share | T₀Ṡ_gen (kW) |
|---|---|---|---|
| valve | 301.35 | 56.4 % | 86.83 |
| precooler (Q̇ = 246.2 kW, fan air out at 521.65 K) | 168.18 | 31.5 % | 48.46 |
| duct (0.35 → 0.28 MPa) | 64.05 | 12.0 % | 18.46 |
| mixer | 1.14 | 0.2 % | 0.33 |
| total | 534.71 | | 154.08 |

(b) **534.71 against 535.81 W/K.** The valve and the duct together throttle the air from 1.0 to 0.28 MPa whatever the intermediate pressure; their sum is R ln(1.0/0.28) = 365.4 W/K per kg/s either way. Moving p_prv only moves entropy generation from the valve to the duct. The small change in the total is the precooler's (465 K instead of 470 K).

(c) Fan air 0.8 kg/s: **113.61 W/K**; 1.2 kg/s: **168.18 W/K**; 2.0 kg/s: **224.02 W/K**. Less fan air is better for the second law (the cold stream runs hotter, closer to the bleed air), but it needs a larger exchanger for the same duty, and the cold outlet (619.9 K at 0.8 kg/s) must stay below the 700 K hot inlet: the temperature cross sets the limit. Choose the smallest flow the exchanger size and weight allow.

(d) Replace the valve by an expansion turbine (the air-cycle-machine idea). An upper bound on the gain is removing the valve entry: the total falls from 534.71 to **233.36 W/K (−56 %)**. The real gain is smaller: a turbine has its own losses, and by cooling the air it also changes the precooler's duty (Extension 1 of the lesson does it properly).

**Common mistakes.** Changing `BLEED` itself instead of a copy; reading a smaller valve entry as a smaller total.

## 6. Judgment

(a) `check_cycle([300, -90], [900, 300])` **rejects** it: S_gen = 90/300 − 300/900 = −0.0333 J/K < 0. Its claimed η = 0.70 exceeds the Carnot 0.667: impossible, whatever its mechanism.

(b) **Possible**: W = 250 J, Q_in = 700 J, η = **0.357** against Carnot **0.708** between 1200 and 350 K; S_gen = 0.619 J/K.

(c) Mixing is cheap only when the streams are close in temperature, pressure and composition: S_gen scales with the square of the temperature difference. Mixing 270 K with 297 K air costs 1.1 W/K; mixing bleed air at 700 K with cabin air would cost hundreds, and mixing different gases always costs −R_U Σ n ln x.

(d) Every lost-work figure T₀S_gen scales by 220/288.15 (falls by 24 %); S_gen itself does not depend on T₀, so the ranking and the decision (attack the valve) are unchanged.

(e) Throttling a liquid near saturation (lesson 03: water or ammonia through a valve): the liquid flashes to a two-phase mixture, h = constant does not mean T = constant, and s must come from the table (Extension 3). Also helium or nitrogen at 30 MPa (lesson 02), where real-gas effects make T change on throttling.

## Check-yourself answers

1. Kelvin–Planck: no cycle can turn heat from a single reservoir entirely into work. Clausius: heat does not flow by itself from cold to hot.
2. η_Carnot = 1 − T_C/T_H. An engine exceeding it could drive a reversed Carnot engine to produce net work from one reservoir, violating Kelvin–Planck.
3. It relates properties only, so it holds between any two equilibrium states, whatever process connects them.
4. S_gen = Q(1/T_C − 1/T_H) = QΔT/(T_H T_C); it vanishes linearly in ΔT.
5. Throttle: h₂ = h₁; ideal gas: T₂ = T₁; s₂ − s₁ = −R ln(p₂/p₁) = R ln(p₁/p₂), and with no heat crossing the boundary all of it is generated.
6. More cooling air for the same duty warms less, so the heat is received at lower temperatures, across a larger mean temperature difference.
7. It raises `ValueError`: a positive Clausius integral means S_gen < 0, a second-law violation.
8. Same: Δs, Δu, Δh (end-state properties). Different: Q, W and S_gen (zero for the reversible path; equal to Δs for the adiabatic free expansion).
