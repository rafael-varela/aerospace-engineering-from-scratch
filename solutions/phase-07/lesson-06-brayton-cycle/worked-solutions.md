# P07.06  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and probe scripts that call it, given to the figures shown. Engine numbers are typical, illustrative values.

## 1. An ideal cycle by hand (s = 15: π = 9)

**Assumptions.** Air standard, isentropic compressor and turbine, no pressure losses, constant cp = 1004.7 J/(kg·K), R = 287.05 J/(kg·K).

(a) π^m = 9^0.2857 = 1.8734.

| state | T (K) | p (kPa) | v (m³/kg) |
|---|---|---|---|
| 1 | 288.15 | 101.33 | 0.81630 |
| 2 | 539.83 | 911.93 | 0.16992 |
| 3 | 1200.00 | 911.93 | 0.37772 |
| 4 | 640.53 | 101.33 | 1.81457 |

w_c = **252.86 kJ/kg**, w_t = **562.07 kJ/kg**, w_net = **309.22 kJ/kg**, q_in = **663.24 kJ/kg**, q_out = **354.02 kJ/kg**, **η = 0.4662**, back-work ratio **0.4499**. Check: 1 − 1/1.8734 = 0.4662 ✓, and q_in − q_out = 309.22 kJ/kg = w_net ✓.

(b) η_Carnot = 1 − 288.15/1200 = **0.7599**; η/η_Carnot = **0.614**.

(c) NASA air: w_net = **329.20 kJ/kg** (+6.5 %), **η = 0.4465** (−2.0 points). The real air's cp rises with temperature: the expansion delivers more work, but more heat is needed to reach 1200 K, and the hot gas's lower γ leaves the turbine exit hotter, so more heat is rejected.

**Common mistakes.** Using p in kPa with R in J/(kg·K) for v (a factor of 1000); computing T₄ = T₃ − (T₂ − T₁).

## 2. Derivations

(a) q_in = cp(T₃ − T₂), q_out = cp(T₄ − T₁). The isentropes give T₂/T₁ = π^m = T₃/T₄, so T₄/T₁ = T₃/T₂ = r and q_out/q_in = T₁(r − 1)/(T₂(r − 1)) = T₁/T₂; **η = 1 − T₁/T₂ = 1 − π^(−m)**. For any ideal gas, s₃ − s₂ = s₄ − s₁ = Δs because s₂ = s₁ and s₄ = s₃. Each heat transfer is internally reversible at constant pressure, so q_in = ∫T ds over 2→3 and q_out = ∫T ds over 4→1 (reversed). Define T̄_in = q_in/Δs and T̄_out = q_out/Δs; then **η = 1 − T̄_out/T̄_in**. T̄_in is an average of T over 2→3, so T̄_in < T₃; T̄_out is an average over 1→4, so T̄_out > T₁. Hence η < 1 − T₁/T₃. ∎

(b) With x = π^m and τ = T₃/T₁: w_t = cpT₃(1 − 1/x), w_c = cpT₁(x − 1), so **w_net = cpT₁[τ(1 − 1/x) − (x − 1)]**. dw/dx = cpT₁(τ/x² − 1) = 0 ⇒ **x = √τ**; d²w/dx² = −2cpT₁τ/x³ < 0, a maximum. Then T₂ = xT₁ = √(T₁T₃) and T₄ = T₃/x = √(T₁T₃).

(c) w_net = h₃ − h₄ − h₂ + h₁ with h₁, h₃ fixed. On the compressor isentrope, dh = v dp, so dh₂/d ln p₂ = v₂p₂ = RT₂, and p₂ = πp₁ gives dh₂/d ln π = RT₂. On the turbine isentrope the exit pressure is p₁ and the inlet pressure πp₁: s°(T₄) = s°(T₃) − R ln π, so cp(T₄) dT₄/T₄ = −R d ln π and dh₄ = cp dT₄ = −RT₄ d ln π. Therefore **dw_net/d ln π = RT₄ − RT₂**. ∎ No assumption on cp was made.

(d) ∮d(pv) = 0 for any closed path because pv is a function of state, and d(pv) = p dv + v dp, so ∮p dv = −∮v dp. For a reversible adiabatic steady-flow device, w_shaft,out = −∫v dp (from dh = v dp and the steady-flow energy equation); the isobars do no shaft work and have dp = 0. Summing the four devices, w_net = −∮v dp = ∮p dv.

**Common mistakes.** Treating T̄ as the arithmetic mean of the end temperatures; forgetting in (c) that T₄ moves when π changes.

## 3. Closed-cycle Brayton power in space

**Assumptions.** Ideal components, monatomic ideal gas with constant cp, R = 8.314/0.040 = 207.86 J/(kg·K), cp = 519.65 J/(kg·K).

(a) τ = 1150/400 = 2.875; π_opt = τ^(γ/(2(γ−1))) = 2.875^1.25 = **3.744** (the numerical optimiser gives 3.74367). There: **w_net = 100.57 kJ/kg, η = 0.4102**, back-work ratio **0.590**.

(b) At x = √τ, η = 1 − 1/x = **1 − √(T₁/T₃)**, which contains no γ. Air between the same temperatures: π_opt = 2.875^1.75 = **6.348**, w_net = **194.44 kJ/kg**, the same η = 0.4102. The gas changes the pressure ratio and the work per kilogram (through cp), not the efficiency at the work optimum. Monatomic gases are used in space Brayton systems for other reasons: a low pressure ratio means few compressor stages, and the inert gas does not attack hot metals.

(c) At T₁ = 400 K: q_out/w_net = 1.4376, index (q_out/w_net)/T₁⁴ = 5.616 × 10⁻¹¹ K⁻⁴. At T₁ = 500 K (π_opt = 2.832, η = 0.3406, w_net = 69.33 kJ/kg): q_out/w_net = 1.9358, index 3.097 × 10⁻¹¹ K⁻⁴. **Ratio 0.55: the hotter radiator is 45 % smaller per watt**, even though the cycle is 7 points less efficient. For a spacecraft, radiator mass usually dominates, so a higher T₁ is chosen; a power station rejects heat to a river or the air at no extra mass cost and is optimised for efficiency. (Using the mean rejection temperature, between T₄ and T₁, would make both radiators smaller but leaves the conclusion unchanged.)

**Common mistakes.** Using γ/(γ − 1) instead of γ/(2(γ − 1)); scaling the radiator with q_out alone instead of per unit net power.

## 4. Optimum and convergence (s = 15: T₃ = 1450 K)

(a) Constant cp: **π_opt = 16.907, w_net = 447.45 kJ/kg, η = 0.5542**, matching the closed form 16.907. NASA air: **π_opt = 20.828, w_net = 492.29 kJ/kg, η = 0.5532**; the closed form (16.907) does not apply, but `brayton_states` at 20.828 gives T₂ = T₄ to the optimiser's tolerance, as (c) of Problem 2 requires.

(b) NASA air at π = 12 (w_net = 476.47 kJ/kg, η = 0.4818):

| n | \|∮p dv − w_net\| (J/kg) | \|∮T ds − w_net\| (J/kg) |
|---|---|---|
| 20 | 832.15 | 47.69 |
| 40 | 207.95 | 11.92 |
| 80 | 51.98 | 2.982 |
| 160 | 13.00 | 0.746 |

Observed orders **2.00, 2.00, 2.00** (p–v) and **2.00, 2.00, 2.00** (T–s), as the trapezoid rule predicts.

(c) In the T–s plane the isentropes are vertical: ds = 0 there and they contribute exactly zero, so all the error comes from the gently curved isobars. In the p–v plane every leg is curved (the isentropes strongly, as v^(−γ)) and all four contribute error.

**Common mistakes.** Forgetting to close the loop (the last point must repeat the first); traversing the loop counter-clockwise and getting −w_net.

## 5. What the ideal cycle can and cannot say

(a) The ideal cycle has **no friction at all**: every component is isentropic. The 64 % (0.5083/0.7942 = 0.640) is not a loss inside the machine but the consequence of where it exchanges heat: it takes heat in at an entropic mean temperature of 934.71 K, not 1400 K, and rejects it at 459.56 K, not 288.15 K. A Carnot engine between 288 and 1400 K would need isothermal heat addition at 1400 K and isothermal rejection at 288 K, which a flowing gas heated by combustion cannot do; approaching it would need infinitely many intercoolers and reheaters with a perfect recuperator (the Ericsson cycle; `P07.07` takes the first step, one intercooler). The comparison with Carnot says why the cycle's efficiency is what it is; it does not say the gas turbine is wasteful.

(b)

| omitted | effect on η | where |
|---|---|---|
| compressor and turbine losses (η_c ≈ 0.85, η_t ≈ 0.88) | lowers: with a 4 % combustor pressure loss as well, from 0.483 to 0.341 at π = 12 (NASA air) | `P07.07` |
| pressure losses in combustor and ducts; cooling air | lowers | `P07.07` |
| real combustion (fuel mass, products' cp, the irreversibility of burning) | lowers the second-law efficiency; changes η by a few points | `P07.08` |

(Also: a recuperator or intercooler, which the ideal simple cycle lacks, can raise η; `P07.07`.)

(c) (1) Near π_max the work per kilogram tends to zero, so the engine would need an enormous air flow for any power. (2) The compressor delivery temperature approaches T₃, which the last compressor stages and the combustor inlet cannot tolerate, and a very high π needs many stages whose losses (in a real machine) erase the ideal gain: with real components the efficiency peaks at a finite π (`P07.07`).

## Check-yourself answers

1. Isentropic compression (compressor), constant-pressure heat addition (combustor), isentropic expansion (turbine), constant-pressure heat rejection (the exhaust and the atmosphere, closing the loop).
2. η = 1 − π^(−(γ−1)/γ): it depends on the pressure ratio (and γ), not on the turbine inlet temperature.
3. The fraction of the turbine work used to drive the compressor, w_c/w_t (0.42 in the example). Net work is the difference of two large numbers, so a small fractional loss in either becomes a large fractional loss in the net.
4. They are equal: T₂ = T₄ (= √(T₁T₃) for constant cp).
5. Every process is internally reversible, so the heat of each leg is ∫T ds and the cycle's net heat, ∮T ds, equals its net work. With an irreversible turbine, T ds on the turbine leg is no longer heat (the turbine is adiabatic but its entropy rises), and the enclosed area overstates the net work.
6. The Brayton cycle adds heat at a mean temperature below the maximum and rejects it at a mean temperature above the minimum, so 1 − T̄_out/T̄_in < 1 − T₁/T₃.
7. Only in the limit π → (T₃/T₁)^(γ/(γ−1)), where the heat added, the heat rejected and the net work all go to zero.
