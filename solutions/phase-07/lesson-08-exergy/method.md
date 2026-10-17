# P07.08  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. T = 600 + 10s. 2. `flow_exergy(N_AIR, T, 1e6, T0, p0)`, `flow_exergy(N_AIR, T, p0, T0, p0)`, `flow_exergy(N_AIR, T0, 1e6, T0, p0)`; add the last two. 3. `flow_exergy(N_AIR, T0, 0.5e5, T0, p0)`; compare with R_mix T₀ ln(0.5/1.01325). 4. q = H(T) − H(T₀) with `mixture_enthalpy`; q(1 − T₀/T); compare.

## 2.
(a) Energy and entropy balances; multiply the entropy balance by T₀; subtract; use Σ_in ṁ = Σ_out ṁ to insert h₀ − T₀s₀; identify ψ and Q(1 − T₀/T). (b) Write (3) for inlet, compressor (work in), combustor (fuel exergy in), turbine (work out), nozzle (kinetic energy out); add; ψ₂, ψ₃, ψ₄, ψ₅ and w_c cancel. (c) η_II = 1 − T₀s_gen/w, η_s = 1 − (h₃ − h₃ₛ)/w; h₃ − h₃ₛ = ∫T ds on the isobar; bound T below by T₀. (d) Write each destruction as T₀Δs and the fuel exergy from B = H − T₀S; the sum reduces to a combination of enthalpies that is zero only if the energy balances hold.

## 3.
1. `st = engine_states(EngineSpec())`; T₃, p₃. 2. ψ at (T₃, p₃), (T₃, 2 bar), (T₀, 2 bar). 3. Throttle destruction = ψ difference; compare with T₀ R_mix ln(p₃/2 bar). 4. Compare ψ(T₀, 2 bar) with `st["w_c"]` and with R_mix T₀ ln(2 bar/p₀).

## 4.
1. `EngineSpec(T_tit=1300 + 10*s, pi_c=8 + s/2)`; `exergy_audit`. 2. Table of product, exhaust and `D`; `residual`. 3. `eta_II`, `eta_th`, `eta_II_comp`. 4. Same spec with T0=216.65, p0=22632, V0=236.1.

## 5.
(a) η_II = product/Ex with the three denominators: `Ex_fuel`, f × `standard_reaction_exergy()`, f × `lower_heating_value()`. (b) n_f = `states["n_fuel"]`; R T₀ ln(p₃/p₀). (c) Condensation of product water; diffusion of CO₂ and H₂O to atmospheric concentrations. (d) Compare the ranking of destructions, which does not depend on the denominator.

## 6.
(a) Read the sweep table; `engine_states(...)["T"][3]` for each case. (b) Think of combustion as heat released by chemistry into the gas, then compare with (1 − T₀/T) weighting at the gas temperatures. (c) Materials and cooling for T₄; delivery temperature, stages and size for π; bypass (turbofan) or heat recovery for the exhaust.
