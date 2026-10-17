# P07.07  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and probe scripts that call it, given to the figures shown. Component values are typical, illustrative values.

## 1. The real simple cycle by hand (s = 15: π = 11)

**Assumptions.** Air standard, steady flow, adiabatic compressor and turbine, no pressure losses, cp = 1004.7 J/(kg·K).

(a) x = 11^0.2857 = 1.9840, τ = 4.8586. T₂ = 288.15[1 + 0.9840/0.85] = **621.72 K**; T₇ = 1400[1 − 0.88(1 − 1/1.9840)] = **788.97 K**. w_c = **335.13 kJ/kg**, w_t = **613.88 kJ/kg**, w_net = **278.75 kJ/kg**, q_in = cp(1400 − 621.72) = **781.90 kJ/kg**, **η = 0.3565**. The ideal cycle at π = 11: η = 0.4960, w_net = 412.73 kJ/kg. Two component efficiencies of 0.85 and 0.88 cost 14 points of efficiency and a third of the work.

(b) The code gives the same numbers for constant cp; NASA air gives **η = 0.3432**, w_net = 306.49 kJ/kg (q_in = 893.06 kJ/kg).

(c) w_net = 0 when η_tτ(1 − 1/x) = (x − 1)/η_c, i.e. **x = η_cη_tτ = 3.634** (π ≈ 91), against x = τ = 4.859 (π = 252.8) for the ideal cycle. The real cycle runs out of work at a much lower pressure ratio, which is why its efficiency also peaks at a finite π.

**Independent check.** First law: q_in − q_out = 781.90 − cp(788.97 − 288.15) = 781.90 − 503.15 = 278.75 kJ/kg = w_net ✓.

**Common mistakes.** Multiplying the compressor work by η_c instead of dividing; applying η_t to the turbine exit temperature directly instead of to the temperature drop.

## 2. Derivations

(a) ε = 1 and isentropic components: T₅ = T₇ = T₃/x, so q_in = cp(T₃ − T₃/x) = cpT₃(1 − 1/x) = w_t. Then η = (w_t − w_c)/w_t = 1 − cpT₁(x − 1)/(cpT₃(1 − 1/x)) = 1 − (T₁/T₃)(x − 1)x/(x − 1) = **1 − (T₁/T₃)x**. Equal to the simple cycle's 1 − 1/x when (T₁/T₃)x = 1/x, **x² = T₃/T₁**, π = (T₃/T₁)^(γ/(2(γ−1))). ∎

(b) w_c = (cpT₁/η_c)[a^m − 1 + (π/a)^m − 1]. d/da: m a^(m−1) − m π^m a^(−m−1) = 0 ⇒ a^(2m) = π^m ⇒ **a = √π**; the second derivative m(m − 1)a^(m−2) + m(m + 1)π^m a^(−m−2) is positive at a = √π (it equals 2m²π^((m−2)/2) > 0), so it is a minimum. Requirements: equal η_c, cooling back to the inlet temperature, no intercooler pressure loss, constant cp.

(c) With a = b = √π each compressor has temperature ratio x^½, so w_c = 2cpT₁(x^½ − 1); the recuperator still delivers T₅ = T₇, so q_in = w_t = cpT₃(1 − 1/x). **η = 1 − 2(T₁/T₃)(x^½ − 1)/(1 − 1/x).**

(d) Per component: s_gen,i = Δs_stream,i + (heat to reservoirs)/(their T) − (heat from reservoirs)/(their T). Summing over the closed loop, the stream entropy changes add to zero (the air returns to state 1), leaving Σs_gen = (q_ic + q_ex)/T₀ − q_in/T_H. The first law of the cycle gives q_ic + q_ex = q_in − w_net. Multiply by T₀: T₀Σs_gen = q_in − w_net − T₀q_in/T_H, so **w_net = q_in(1 − T₀/T_H) − T₀Σs_gen**. ∎

**Common mistakes.** In (d), forgetting that the combustor's s_gen includes −q_in/T_H (heat from a reservoir); in (b), splitting π arithmetically.

## 3. Should an APU have a recuperator? (π = 4, 1250 K, NASA air)

(a)

| ε | η | w_net (kJ/kg) | q_in (kJ/kg) | T₅ (K) |
|---|---|---|---|---|
| 0 (no recuperator) | 0.2220 | 196.03 | 883.17 | 451.79 |
| 0.60 | 0.3208 | 179.79 | 560.40 | 758.28 |
| 0.75 | 0.3748 | 179.79 | 479.71 | 831.80 |
| 0.85 | 0.4221 | 179.79 | 425.91 | 880.20 |
| 0.95 | 0.4831 | 179.79 | 372.12 | 928.15 |

w_net is the same for every ε because the turbine's expansion ratio is set by the fixed pressure losses (the cold-side 3 % and hot-side 4 %), and nothing the recuperator does downstream changes the turbine inlet or exit state. The 16.24 kJ/kg drop from 196.03 is the cost of those losses.

(b) Heat per unit work: 883.17/196.03 = 4.505 without, 425.91/179.79 = 2.369 with ε = 0.85: **47 % less fuel per kilowatt-hour**, at 8 % less power per kilogram of air.

(c) A ground APU runs for long hours and its weight costs nothing; an aircraft APU runs briefly and its weight is carried for the whole flight, so the fuel saved may be smaller than the fuel spent carrying a recuperator. (Also: volume in the tail cone.)

**Common mistakes.** Recomputing the turbine with a different exit pressure for each ε; quoting efficiency gain without the power loss.

## 4. The design study at your own temperature (s = 15: 1375 K)

Grid of step 0.5 in π, NASA air:

| configuration | η_max at π | w_max (kJ/kg) at π | η at π = 12 |
|---|---|---|---|
| simple | 0.3676 at 26.0 | 286.8 at 11.0 | 0.3385 |
| intercooled | 0.3711 at 40.0 (still rising) | 362.4 at 25.0 | 0.3283 |
| recuperated | 0.4557 at 5.5 | 273.0 at 11.5 | 0.4152 |
| intercooled + recuperated | 0.4883 at 7.5 | 351.3 at 26.5 | 0.4813 |

(b) (i) A small, light turboprop core: no heat exchangers (their weight dominates), a simple cycle near π ≈ 11–15, close to its work maximum (286.8 kJ/kg) where the engine is smallest for its power, accepting η ≈ 0.34. If the mission is long enough that fuel outweighs hardware, an intercooled-recuperated core at π ≈ 8 (η = 0.488) is the alternative. (ii) A continuous ground generator: intercooled-recuperated at π ≈ 7.5–10, for 0.48–0.49 against 0.34, with a few compressor stages; weight is irrelevant and fuel is the running cost.

**Common mistakes.** Choosing the efficiency peak of the intercooled cycle at π = 40, where it is only 0.37 and needs many stages.

## 5. Where the work goes (π = 12, 1400 K, intercooled + recuperated, NASA air)

(a) T₀ = 288.15 K:

| rank | component | s_gen (J/(kg·K)) | T₀s_gen (kJ/kg) |
|---|---|---|---|
| 1 | exhaust (to the atmosphere) | 186.69 | 53.80 |
| 2 | combustor (heat from 1400 K) | 165.34 | 47.64 |
| 3 | turbine | 106.74 | 30.76 |
| 4 | intercooler | 100.98 | 29.10 |
| 5 | recuperator | 82.15 | 23.67 |
| 6 | HP compressor | 52.37 | 15.09 |
| 7 | LP compressor | 51.69 | 14.89 |
| | total | 745.97 | 214.95 |

(b) q_in(1 − T₀/T_H) = 702.38 × (1 − 288.15/1400) = **557.82 kJ/kg**; minus 214.95 = **342.87 kJ/kg** = w_net ✓ (the code's residual is at round-off).

(c) Four of the five largest terms are heat transfers across finite temperature differences: the combustor takes heat from 1400 K into air at 792–1400 K, the exhaust dumps 500 K gas into 288 K air, the intercooler dumps 432 K air into a 288 K sink, and the recuperator moves heat across the gap between 852→500 K and 435→792 K. The turbomachines total 60.74 kJ/kg. A higher recuperator effectiveness shrinks both the exhaust and recuperator terms (smaller temperature differences); cooling the intercooler's air to a sink nearer its own temperature cannot help, because its heat is rejected anyway; the combustor term is largely the air-standard stand-in for combustion, whose real irreversibility `P07.08` computes.

## 6. Recuperators on high-pressure-ratio aero engines

(a) NASA air, 1400 K, lesson component set:

| π | simple | recuperated | intercooled + recuperated |
|---|---|---|---|
| 25 | 0.3726 | 0.3371 | 0.4542 |
| 30 | 0.3726 | 0.3087 | 0.4422 |
| 40 | 0.3632 | 0.2584 | 0.4207 |

At π = 40 the HP compressor delivers air at 886.8 K without intercooling but 523.7 K with it, while the turbine exhaust is at 667.5 K. Without intercooling T₇ < T₄ and the "recuperator" cools the compressed air; with intercooling T₇ > T₄ again and the recuperator has heat to give. That is the argument behind intercooled-recuperated aero-engine concepts.

(b) Doubled losses (4 %, 6 %, 8 %, 8 %): η = **0.4644** (−2.37 points) and w_net = **316.86 kJ/kg** (−7.6 % from 342.87).

(c) (1) **Weight and volume** of two heat exchangers and their ducts: against, because the weight is carried for every flight. (2) **Frontal area and drag** of the intercooler's cooling-air path: against. (3) **Reliability and maintenance** of thin-walled exchangers under thermal cycling: against. On the other side, (4) lower specific fuel consumption on long missions argues for, and is the reason the concepts keep being studied.

## Check-yourself answers

1. ε = (h₅ − h₄)/(h(T₇) − h₄): the fraction of the largest possible cold-side enthalpy rise. For ε = 1, T₅ = T₇.
2. Heat flows from the exhaust to the compressed air only if the turbine exit is hotter than the compressor exit (T₇ > T₄). For an ideal cycle that holds below π = (T₃/T₁)^(γ/(2(γ−1))), the pressure ratio of maximum work (15.9 at 1400 K).
3. η = 1 − (T₁/T₃)π^((γ−1)/γ). It falls with π, because a higher π raises the compressor exit temperature and lowers the turbine exit temperature, leaving the recuperator less to recover, while the simple cycle gains from compression directly.
4. Intercooling lowers the compressor delivery temperature, so the combustor must add more heat to reach the same turbine inlet temperature; the extra heat is more than the work saved is worth, because the heat removed in the intercooler is thrown away.
5. √π for each, when both compressors have the same efficiency, the intercooler returns the air to the inlet temperature, there is no intercooler pressure loss and cp is constant.
6. They reduce the turbine's expansion ratio below the compressor's pressure ratio: the compressor does the full work of compressing by π, but the turbine can only expand through π × (1 − Δ₁)(1 − Δ₂)…
7. w_net = q_in(1 − T₀/T_H) − T₀Σs_gen: the Carnot work of the heat minus T₀ times the entropy generated in every component (Gouy–Stodola).
8. The turbine efficiency (−3.15 points for −0.05): the turbine's work is the larger of the two works, so a fractional loss in it removes more from the net.
