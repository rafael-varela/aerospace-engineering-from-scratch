# P07.08  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and probe scripts that call it, given to the figures shown. Engine values are typical, illustrative values. Everything is computed; nothing is burned.

## 1. The exergy of a stream and of heat (s = 15: T = 750 K)

**Assumptions.** Dry-air ideal-gas mixture (NASA), dead state 288.15 K and 101 325 Pa, no kinetic energy.

(a) ψ(750 K, 1 MPa) = **384.37 kJ/kg**; ψ(750 K, p₀) = **195.01 kJ/kg**; ψ(T₀, 1 MPa) = **189.36 kJ/kg**; 195.01 + 189.36 = 384.37 ✓. For an ideal gas h depends on T only and s = s°(T) − R ln(p/p°), so ψ = [(h − h₀) − T₀(s° − s°₀)] + RT₀ln(p/p₀): a thermal part and a pressure part that add. The pressure part checks against RT₀ln(1 MPa/p₀) = 189.36 kJ/kg.

(b) ψ(T₀, 0.5 bar) = **−58.42 kJ/kg** = RT₀ln(0.5/1.01325). Negative: the stream is below ambient pressure and work must be spent to discharge it into the atmosphere (a vacuum is a resource only if the atmosphere does the work of filling it).

(c) q = H(750) − H(T₀) = **479.19 kJ/kg**; q(1 − T₀/T) = **295.09 kJ/kg**; the thermal flow exergy is **195.01 kJ/kg**, smaller. q(1 − T₀/T) would be right only if all the heat were available at 750 K; as the gas cools, its heat is released at falling temperatures down to T₀, and each increment dq is worth only dq(1 − T₀/T_current). Integrating gives ψ.

**Common mistakes.** Using q(1 − T₀/T_max) as the exergy of a cooling stream; forgetting the pressure part.

## 2. Derivations

(a) Energy: 0 = Σ_in ṁ(h + V²/2) − Σ_out ṁ(h + V²/2) + ΣQ̇_k − Ẇ. Entropy: 0 = Σ_in ṁs − Σ_out ṁs + ΣQ̇_k/T_k + Ṡ_gen. Multiply the second by T₀ and subtract: Ẇ = Σ_in ṁ(h + V²/2 − T₀s) − Σ_out ṁ(h + V²/2 − T₀s) + ΣQ̇_k(1 − T₀/T_k) − T₀Ṡ_gen. Since Σ_in ṁ = Σ_out ṁ in steady flow of a single stream (and per stream where species are tracked), add Σ_out ṁ(h₀ − T₀s₀) − Σ_in ṁ(h₀ − T₀s₀) = 0 to turn each bracket into ψ: **Ẇ = Σ_in ṁψ − Σ_out ṁψ + ΣQ̇_k(1 − T₀/T_k) − T₀Ṡ_gen**. ∎

(b) Inlet: D_i = ψ₀ + V₀²/2 − ψ₂. Compressor: D_c = ψ₂ + w_c − ψ₃. Combustor: D_cc = ψ₃ + Ex_fuel − m_pψ₄ (written per kilogram of air; the change of composition is carried by the fuel exergy). Turbine: D_t = m_pψ₄ − w_c − m_pψ₅. Nozzle: D_n = m_pψ₅ − m_pψ₉ − m_pV₉²/2. Add them: ψ₂, ψ₃, ψ₄, ψ₅ and w_c cancel, and ψ₀ = 0, leaving ΣD = Ex_fuel + V₀²/2 − m_pV₉²/2 − m_pψ₉, which rearranges to (5). ∎

(c) η_II = (ψ₃ − ψ₂)/w = (w − T₀s_gen)/w. η_s = w_s/w = 1 − (h₃ − h₃ₛ)/w. States 3s and 3 lie on the same isobar, where dh = T ds, so h₃ − h₃ₛ = ∫_{3s}^{3} T ds. On that isobar T ≥ T₃ₛ > T₀ (when the exit is above ambient), so ∫T ds ≥ T₀(s₃ − s₃ₛ) = T₀(s₃ − s₂) = T₀s_gen. Hence **η_II ≥ η_s**. ∎ (Sea level: 0.9286 against 0.85.)

(d) Each D_k = T₀(S_out − S_in) for the component (with the fuel's entropy at the combustor). Their sum is T₀(m_ps₉ − s_air,0 − n_f s_f). The fuel exergy (4) is B_air(T₀) + n_f b_f − B_prod(T₀) and the exhaust exergy is (H₉ − H_prod,0) − T₀(S₉ − S_prod,0). Then Ex_fuel − product − exhaust − ΣD = H_air(T₀) + V₀²/2 + n_f h_f − H₉ − m_pV₉²/2, the energy in minus the energy out of the whole engine, which is the sum of the component energy residuals. It vanishes only if every component's first law holds. ∎

**Common mistakes.** Treating the fuel's chemical exergy as heat at the flame temperature; forgetting the ram kinetic energy in flight.

## 3. What bleed air costs (1 kg/s from station 3)

**Assumptions.** Bleed at the sea-level engine's station 3 (630.58 K, 1191.58 kPa), throttle adiabatic, cooler rejecting to the atmosphere.

(a) ψ₃ = **324.36 kJ/kg**; after the throttle (630.58 K, 2 bar) **176.74 kJ/kg**; after the cooler (288.15 K, 2 bar) **56.24 kJ/kg**.

(b) The throttle destroys **147.62 kW** (per kg/s). Check: T₀R_mix ln(1191.58/200) = 147.62 kJ/kg ✓. A throttle is adiabatic with no work, so h is conserved, and an ideal gas's h depends on T only: T is unchanged and only the entropy rises, by R ln(p_in/p_out). The cooler then throws away a further 120.50 kJ/kg of thermal exergy.

(c) The compressor spent **351.10 kJ/kg** on this air; **56.24 kJ/kg** is delivered, 16 %. An ideal isothermal compressor taking ambient air to 2 bar needs exactly RT₀ln(2/1.01325) = **56.24 kJ/kg**, the delivered exergy, by definition. Over-compressing air and then throttling and cooling it is an expensive way to make cabin air, which is the thermodynamic case for an electric cabin compressor (`P07.05` Problem 3) that compresses only as far as the cabin needs.

**Common mistakes.** Expecting the throttle to cool an ideal gas; forgetting that the bleed air's compressor work was paid with fuel.

## 4. Your own engine (s = 15: T₄ = 1450 K, π_c = 15.5)

(a) Sea-level static, exergy in 995.66 kJ/kg (f = 0.02206, V₉ = 874.9 m/s):

| item | kJ/kg air | % of exergy in |
|---|---|---|
| product (jet kinetic energy) | 391.13 | 39.3 |
| exhaust exergy | 244.33 | 24.5 |
| combustor | 313.09 | 31.4 |
| compressor | 26.78 | 2.7 |
| turbine | 14.30 | 1.4 |
| nozzle | 4.36 | 0.4 |
| inlet | 1.67 | 0.2 |

Residual **2.9 × 10⁻¹⁰ J/kg**.

(b) **η_II = 0.3928, η_th = 0.4089**; compressor 0.9336, turbine 0.9657, nozzle 0.9890. Against the lesson's engine (0.3670), the higher pressure ratio and turbine temperature cut the combustor's share from 32.9 % to 31.4 % and the exhaust's from 26.0 % to 24.5 %.

(c) Cruise: **η_II = 0.4603**; combustor **24.7 %**, exhaust **23.8 %**.

**Common mistakes.** Dividing by the fuel exergy alone at cruise, leaving out the ram exergy.

## 5. How much does the fuel-exergy convention matter?

(a) Product 360.45 kJ/kg in all three. (i) Ex_fuel = 982.20 kJ/kg: **η_II = 0.3670**. (ii) f × (−ΔG°) = 961.71 kJ/kg: **0.3748**. (iii) f × LHV = 943.59 kJ/kg (this is the thermal efficiency): **0.3820**. Spread **1.5 points**.

(b) n_f R T₀ ln(p₃/p₀) = **0.77 kJ/kg of air**, 0.08 % of the fuel exergy: negligible.

(c) (1) **Condensation of the product water** at the dead state: liquid water has lower entropy and enthalpy than vapour, so B_prod falls and the fuel exergy rises (towards a "higher heating value" basis). (2) **Diffusion of CO₂ and H₂O to their atmospheric concentrations** (an unrestricted dead state): the products could do further work in spreading out, which also raises the fuel exergy. A third: the fuel is really a liquid, whose exergy differs from the vapour's by about its vaporisation terms.

(d) No. All three conventions change only the denominator; every destruction, the exhaust exergy and the product are computed from the same states. The ranking (combustor ≫ exhaust > compressor > turbine > nozzle > inlet) and their ratios are unchanged, so the engineering conclusion is the same: it is the cycle (pressure ratio, turbine temperature, exhaust recovery), not the components, that has most to give.

## 6. Attacking the biggest loss

(a) (1400 K, π = 12): combustor **32.9 %** (0.3285 of the fuel exergy), **η_II = 0.3670**, compressor delivery 630.58 K. (1600 K, π = 30): combustor **28.1 %**, **η_II = 0.4512**, delivery 819.24 K. The combustor's share falls with both.

(b) In the combustor, chemical energy with a work potential close to its heating value is released as random molecular energy in gas at 630–1400 K. Had that energy been delivered as heat at those temperatures, it would be worth only (1 − T₀/T) of itself, and the mixing of fuel and air and the reaction itself are spontaneous processes far from equilibrium, which generate entropy without any friction or pressure loss. The hotter the gas the energy is released into (hotter compressor delivery air, hotter products), the larger (1 − T₀/T) and the less is destroyed.

(c) The turbine inlet temperature is limited by blade materials and the cooling air that costs exergy of its own; the pressure ratio by compressor delivery temperature (discs, combustor inlet), the number of stages, and the falling specific work (`P07.06`). The exhaust's exergy can be partly recovered by a bypass stream that turns some of the core's hot-gas exergy into a slower, larger, more propulsively efficient cold jet (a turbofan; Extension 2), or on the ground by a recuperator or a steam bottoming cycle (`P07.07`).

## Check-yourself answers

1. The state of the surroundings (T₀, p₀, at rest, and for chemical exergy a reference composition) with which a system comes to equilibrium. Exergy is the work obtainable on the way there, so it depends on where "there" is; energy needs no reference beyond an arbitrary zero.
2. ψ = (h − h₀) − T₀(s − s₀) + V²/2. Yes: air at T₀ below ambient pressure has ψ = RT₀ln(p/p₀) < 0 (−58.42 kJ/kg at 0.5 bar).
3. Exergy destroyed = T₀Ṡ_gen. It follows from combining the first law (energy balance) and the second law (entropy balance) of a control volume.
4. Because it converts chemical energy, almost all of which is work potential, into thermal energy of gas at a few hundred to 1400 K, through spontaneous mixing and reaction; energy is conserved but the entropy rises steeply (1119.89 J/(kg·K) at sea level).
5. Higher (0.9286 at sea level): the extra work a real compressor needs is partly still present as exergy of its hotter exit air; only T₀s_gen is lost.
6. The combustion reaction increases the number of gas moles (18.75 → 23.5), so its entropy change is positive and −ΔG° = −ΔH° + T₀ΔS° exceeds −ΔH° (the LHV); diluting the products in excess air raises it a little more.
7. Not by the engine: it leaves the engine as exergy of the hot jet and is destroyed outside, as the jet mixes with the atmosphere. A bypass stream (turbofan), a recuperator or a bottoming cycle could recover part of it.
8. The combustor's destruction (33.4 %) and the exergy left in the hot exhaust (21.1 %).
