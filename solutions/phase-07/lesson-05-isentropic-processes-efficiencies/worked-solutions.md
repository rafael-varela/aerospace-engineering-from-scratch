# P07.05  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and probe scripts that call it, given to the figures shown. Machine efficiencies are typical, illustrative values.

## 1. One axial stage (s = 15: π = 1.35)

**Assumptions.** Steady, adiabatic, ideal gas, kinetic energy changes neglected (stagnation states).

(a) Constant cp: T₂ₛ = 288.15 × 1.35^0.2857 = **313.947 K**, T₂ = 288.15 + 25.797/0.90 = **316.814 K**. NASA air: T₂ₛ = **313.941 K**, T₂ = **316.805 K**.

(b) w = **28.797 kJ/kg** (both gases, to the figures shown); s_gen = **9.131 J/(kg·K)**; T₀ s_gen = **2.631 kJ/kg**, 9.1 % of the work.

(c) η_p = **0.9042**. With a = m ln π small, (8) of the lesson expands as η_c ≈ η_p[1 − a(1 − η_p)/(2η_p) + …]: the difference is first order in a = 0.0857, multiplied by (1 − η_p)/(2η_p) ≈ 0.056, about 0.5 % (0.9000 against 0.9042).

**Independent check.** At this small π the variable-cp correction to T₂ is 0.009 K: γ = 1.4 is essentially exact near room temperature, as test 3 asserts.

**Common mistakes.** Applying η to temperatures with variable cp (here harmless, at π = 40 not); quoting s_gen per mole.

## 2. Derivations

(a) T ds = dh − v dp with dh = cp(T) dT and v = RT/p gives ds = cp dT/T − R dp/p. Integrate from 1 to 2: s₂ − s₁ = [s°(T₂) − s°(T₁)] − R ln(p₂/p₁), with s°(T) = ∫cp dT/T. Setting s₂ = s₁: **s°(T₂ₛ) − s°(T₁) = R ln π**. With constant cp, s°(T₂) − s°(T₁) = cp ln(T₂/T₁), so T₂ₛ/T₁ = π^(R/cp); with cp − cv = R and γ = cp/cv, R/cp = (γ − 1)/γ: **T₂ₛ/T₁ = π^((γ−1)/γ)**.

(b) dh = v dp/η_p ⇒ cp dT = RT dp/(η_p p) ⇒ cp dT/T = (R/η_p) d ln p ⇒ **s°(T₂) − s°(T₁) = (R/η_p) ln π**. With constant cp, T₂/T₁ = π^(m/η_p), m = (γ − 1)/γ. Then η_c = (h₂ₛ − h₁)/(h₂ − h₁) = (T₂ₛ − T₁)/(T₂ − T₁) = **(π^m − 1)/(π^(m/η_p) − 1)**. Turbine: dh = η_p v dp gives s°(T₂) − s°(T₁) = η_p R ln(p₂/p₁); with π_t = p₁/p₂, T₂/T₁ = π_t^(−η_p m) and η_t = (T₁ − T₂)/(T₁ − T₂ₛ) = **(1 − π_t^(−η_p m))/(1 − π_t^(−m))**.

(c) Along the path, T ds = dh − v dp and v dp = η_p dh, so **T ds = (1 − η_p) dh** and ∫T ds = (1 − η_p) w. The entropy change s₂ − s₁ = [s°(T₂) − s°(T₁)] − R ln π = (R/η_p) ln π − R ln π, so **s_gen = R ln π (1/η_p − 1)**; cp never appears. (Adiabatic, so s_gen = s₂ − s₁.)

(d) Let x = π^m > 1 and k = 1/η_p > 1. Then η_c = (x − 1)/(x^k − 1). The function x^k is strictly convex for k > 1, so it lies above its tangent at x = 1: x^k > 1 + k(x − 1) for x ≠ 1. Hence x^k − 1 > k(x − 1) and **η_c < 1/k = η_p**. ∎

**Common mistakes.** Forgetting the −R ln(p/p°) term; using the polytropic exponent n of p vⁿ = const without relating it to η_p (they are related: (n − 1)/n = m/η_p for constant cp).

## 3. An electric cabin-air compressor

**Assumptions.** Steady, adiabatic, NASA air, inlet at the static outside state (ram ignored), dead state the outside air.

(a) π = 75 000/22 632 = 3.3139. T₂ₛ = 305.22 K, **T₂ = 327.31 K**; **w = 111.10 kJ/kg** (isentropic 88.88 kJ/kg); **power = 133.3 kW** at 1.2 kg/s.

(b) s_gen = **70.28 J/(kg·K)**; exergy destroyed T₀ s_gen = 216.65 × 70.28 = **15.23 kJ/kg** (13.7 % of the work).

(c) η_p = **0.8303**, above η_s = 0.80 as for any compressor. The air leaves at 327 K and must be cooled to about 290 K before the cabin: 37.5 kJ/kg is rejected in a heat exchanger to ram air. Part of the compressor work was spent on heating air that is then cooled again. That is not a design error: the cabin needs pressure, and adiabatic compression delivers pressure and temperature together.

**Independent check (Gouy–Stodola, `P07.08`).** The flow exergy of the delivered air relative to the outside air is ψ₂ = (h₂ − h₀) − T₀(s₂ − s₀) = 95.87 kJ/kg, which equals w − T₀ s_gen = 111.10 − 15.23 = 95.87 kJ/kg. The exergy put in as work either reaches the air or is destroyed.

**Common mistakes.** Using T₀ = 288.15 K; forgetting that 75 kPa is an absolute pressure.

## 4. A multi-stage compressor (s = 15: n = 9 stages, π = 35)

(a) T_exit = **871.48 K**; overall η_s = **0.8244**; overall η_p = **0.8863**.

(b) Polytropic limit at η_p = 0.88: **877.92 K**; the nine-stage machine is **6.45 K** cooler than it, because a stage of finite pressure ratio with η_s = 0.88 has η_p slightly above 0.88 (as Problem 1 showed for 0.90).

(c) **61 stages** bring the exit within 1 K of the limit: the convergence is first order in 1/n.

(d) First stage 39.31 kJ/kg, last stage **105.34 kJ/kg**: every stage has the same pressure ratio and efficiency, but the last compresses air that is three times hotter, and isentropic work scales with inlet temperature (the preheat effect).

**Independent check.** The constant-cp closed form at η_p = 0.88 and π = 35 gives η_s = 0.8111; the nine-stage constant-cp machine has η_p = 0.8865, close to the NASA machine's 0.8863: the stage-count effect hardly depends on cp.

**Common mistakes.** Splitting the pressure ratio as π/n instead of π^(1/n).

## 5. A high-pressure turbine

(a) NASA air: **T₂ = 1186.39 K, η_t = 0.9144, w = 495.76 kJ/kg** (out).

(b) γ = 1.4: T₂ = 1120.22 K, w = 482.01 kJ/kg (η_t = 0.9169, the closed form). γ = 1.33: T₂ = 1174.02 K, w = 492.81 kJ/kg. **γ = 1.33 is much closer** (0.6 % low in work, 12 K low in T₂, against 2.8 % and 66 K for 1.4), because NASA air has γ = **1.307** at 1600 K and an average near 1.32 over the expansion. Hot-section calculations with γ = 1.4 understate the work and the exit temperature.

(c) Reheat: friction in the early part of the expansion leaves the gas hotter than isentropic, and the later part expands that hotter gas, recovering part of the loss. The integrated machine is better than its infinitesimal stages, η_t > η_p.

**Common mistakes.** Passing the expansion ratio 4 instead of the pressure ratio 1/4 to a function that expects p₂/p₁; quoting turbine work as negative and then dividing by a positive ideal work.

## 6. What limits the answer?

(a) NASA air, π = 40:

| η_p | T₂ (K) | ΔT₂ from 0.90 (K) | w (kJ/kg) | Δw | η_s |
|---|---|---|---|---|---|
| 0.88 | 912.74 | +21.51 | 659.14 | +3.80 % | 0.8125 |
| 0.89 | 901.81 | +10.58 | 646.87 | +1.87 % | 0.8279 |
| 0.90 | 891.23 | 0 | 635.02 | 0 | 0.8433 |
| 0.91 | 880.99 | −10.25 | 623.57 | −1.80 % | 0.8588 |
| 0.92 | 871.06 | −20.18 | 612.49 | −3.55 % | 0.8743 |

A ±0.01 uncertainty in η_p moves the work by about ±1.8 %, more than the entire constant-cp error in work (1.45 %); it moves T₂ by about ±10 K, less than the constant-cp error in T₂ (38 K).

(b) Adequate: compressor work and cycle efficiency estimates up to π ≈ 20 (work error below 0.74 %), trends and trade studies. Not adequate: delivery temperature at high π (17 K error at π = 20, 38 K at 40), anything in the hot section (Problem 5), and fuel-flow calculations, which depend on the combustor inlet temperature.

(c) (1) **Heat loss through the casing** lowers the delivery temperature without lowering the shaft work, so the temperature-rise efficiency reads **high**. (2) **A bleed** taken from an intermediate stage leaves the delivery temperature unchanged but means the shaft work compresses less air than is delivered at the end only partly; an efficiency computed from shaft work and the delivery flow reads **high** unless the bleed is accounted for, while the temperature-rise efficiency is unaffected. (Others: probes reading closer to static than total temperature read low and push η high; non-uniform exit temperature profiles averaged by area instead of mass flow.)

**Common mistakes.** Treating an efficiency as a measured fact rather than a derived quantity with its own measurement errors.

## Check-yourself answers

1. s°(T₂ₛ) − s°(T₁) = R ln(p₂/p₁), with s°(T) = ∫cp(T)/T dT the entropy at the reference pressure p° = 1 bar, tabulated by the NASA polynomials.
2. Because the efficiency is meant to compare works, and the work of an adiabatic machine is the enthalpy change; with variable cp, temperature differences are not proportional to enthalpy differences.
3. The turbine: the reheat effect makes η_t > η_p, while the preheat effect makes η_c < η_p. With constant cp at ratio 12 and η_p = 0.90, η_t = 0.929 and η_c = 0.861.
4. Zero at η = 1 (reversible and adiabatic); positive at η = 0.85, s_gen = s₂ − s₁ > 0.
5. The γ = 1.4 work (644.25 against 635.02 kJ/kg): the real cp rises with T, so the real isentrope and the real compression end cooler, and the real gas needs less work for the same pressure ratio.
6. It falls (preheat), approaching the isentropic efficiency of a compressor with polytropic efficiency equal to the stage efficiency.
7. Friction heat is generated at temperatures above T₀, and heat at temperature T can still do work (T − T₀)/T per joule; only T₀ s_gen is permanently lost.
