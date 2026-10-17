# P07.08  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which are closed forms and which are numerical. Name the dead state before computing any exergy. All work in this problem set is computed; nothing is burned.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Dead state unless stated: T₀ = 288.15 K, p₀ = 101 325 Pa. Engine values are the typical, illustrative defaults of `EngineSpec`.

## 1. (Warm-up: the exergy of a stream and of heat)

Air (the lesson's dry-air mixture) is at T = 600 + 10s K.

(a) Compute its flow exergy at 1 MPa and at p₀ with `flow_exergy`, and the exergy of air at T₀ and 1 MPa. Show that the first equals the sum of the other two, and explain why.
(b) Compute the flow exergy of air at T₀ and 0.5 bar. Interpret the sign.
(c) The air at T and p₀ is cooled to T₀ at constant pressure. Compute the heat removed, q, and compare q(1 − T₀/T) with the flow exergy from (a). Which is larger, and why?

Answer form: three numbers and a sentence; a number and a sentence; three numbers and a sentence.

## 2. (Derivation: exergy balance, audit identity, compressor efficiency)

(a) From the steady-flow energy and entropy balances of a control volume with heat Q̇_k at boundary temperatures T_k, derive Ẇ = Σ_in ṁψ − Σ_out ṁψ + ΣQ̇_k(1 − T₀/T_k) − T₀Ṡ_gen.
(b) Add the exergy balances of the five components of the turbojet core and derive the audit identity, equation (5) of the lesson. Which internal terms cancel?
(c) Prove that an adiabatic compressor's second-law efficiency (ψ₃ − ψ₂)/w is at least its isentropic efficiency, when its exit temperature is above T₀.
(d) Show that, for the lesson's fuel-exergy definition, Ex_fuel − (product + exhaust exergy + ΣD) = 0 requires the first law to close in every component.

Answer form: four derivations.

## 3. (Unfamiliar scenario: what bleed air costs)

An aircraft takes 1 kg/s of air from the core compressor exit of the sea-level engine (station 3) for its air-conditioning packs. The air is throttled to 2 bar and then cooled to 288.15 K at 2 bar before use.

(a) Compute the flow exergy of the bleed air at station 3, after the throttle, and after the cooler.
(b) How much exergy does the throttle destroy? Check your answer against T₀R ln(p₃/p) for an ideal gas, and say why throttling an ideal gas leaves its temperature unchanged.
(c) Compare the exergy finally delivered with the compressor work spent on that air, and with the work of an ideal isothermal compressor taking ambient air to 2 bar. What does this say about bleed air, and about the electric cabin compressor of `P07.05` Problem 3?

Answer form: three numbers; a number and a check; two comparisons and two sentences.

## 4. (Numerical, uses your code: your own engine)

Take a turbine inlet temperature of 1300 + 10s K and a compressor pressure ratio of 8 + s/2, other inputs at their defaults.

(a) Run `exergy_audit` at sea-level static. Tabulate the product, the exhaust exergy and every destruction in kJ/kg and as a percentage of the exergy in, and confirm the closure.
(b) Give η_II, η_th and the component second-law efficiencies. Compare with the lesson's engine.
(c) Repeat at cruise (11 km, 216.65 K, 22 632 Pa, 236.1 m/s) and report η_II and the shares of the combustor and the exhaust.

Answer form: a table and a residual; six numbers and a comparison; three numbers.

## 5. (Judgment: how much does the fuel-exergy convention matter?)

(a) For the lesson's sea-level engine, compute η_II with the fuel exergy taken (i) by the lesson's equation (4), (ii) as f × (−ΔG°), (iii) as f × LHV. How far apart are they?
(b) Estimate the pressure term n_f R T₀ ln(p₃/p₀) in the fuel's exergy, which the lesson's assumption (vapour fuel supplied at combustor pressure) includes. Is it significant?
(c) Name two further physical effects that a less restricted dead state would include, and say in which direction each would move the fuel exergy.
(d) Does any of this change which component an engineer should work on? Explain.

Answer form: three efficiencies and a spread; a number and a sentence; two effects with directions; a paragraph.

## 6. (Judgment: attacking the biggest loss)

A colleague proposes to "fix the combustor, since it destroys a third of the fuel's exergy".

(a) Using the lesson's sweep, how do the combustor's share and η_II change with turbine inlet temperature and compressor pressure ratio? Give numbers for (1400 K, π = 12) and (1600 K, π = 30), and the compressor delivery temperature at each.
(b) Explain physically why a combustor destroys exergy even when it is adiabatic, perfectly mixed and loses no pressure, and why that destruction falls as the air entering it gets hotter.
(c) What limits how far the turbine inlet temperature and pressure ratio can be raised? What other route recovers part of the exhaust's 26 %?

Answer form: four numbers and two temperatures; a paragraph; two limits and one route.
