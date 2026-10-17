# P07.07  Real Brayton Cycles

> Real components pull the efficiency down and the best pressure ratio with it. A recuperator gives back what the exhaust would have thrown away, but only while the exhaust is hotter than the compressed air; an intercooler saves compressor work, but only a recuperator turns that saving into efficiency.

| | |
|---|---|
| **Learning outcomes** | 1. Build the station table of a gas-turbine cycle with component isentropic efficiencies, duct and combustor pressure losses, an optional intercooler and an optional recuperator, and close the first law for every component and for the cycle. 2. Derive the ideal recuperated-cycle efficiency 1 − (T₁/T₃)π^((γ−1)/γ), the ideal intercooled-recuperated closed form, the best intercooler pressure split √π, and the crossover pressure ratio above which a recuperator hurts. 3. Implement and test the model so that its ideal limits recover `P07.06` and the closed forms, every component generates non-negative entropy, the total equals the Clausius integral, and the efficiency stays below Carnot. 4. Run a design study for a turboprop or APU core: compare simple, intercooled, recuperated and intercooled-recuperated cycles across pressure ratios, and rank the sensitivities. 5. Decide which configuration a given application should use and say what the model leaves out. |
| **Prerequisites** | `P07.06` (ideal Brayton cycle, optimum pressure ratio, T₂ = T₄ condition), `P07.05` (isentropic efficiency, variable-cp isentropes), `P07.04` (entropy generation and the Clausius integral), `P07.02` (steady-flow energy equation), `P03.14` (one-dimensional optimisation). Skill check: you can compute the exit temperature of a compressor of isentropic efficiency 0.85 at π = 12. |
| **Engineering question** | A turboprop or APU core of 5 kg/s (typical, illustrative) runs at π = 12 and 1400 K turbine inlet temperature, with compressors of isentropic efficiency 0.85, a turbine of 0.88, a 4 % combustor pressure loss, an optional intercooler (2 % loss) and an optional recuperator (effectiveness 0.85, 3 % and 4 % losses): what efficiency and power does each configuration give, and at which pressure ratio would each be best? The answers (NASA air): simple 0.341 and 1.49 MW; intercooled 0.330 and 1.78 MW; recuperated 0.424 and 1.42 MW; intercooled and recuperated 0.488 and 1.71 MW. The recuperated cycle peaks at π ≈ 5 (η = 0.461) and the intercooled-recuperated one at π ≈ 8 (η = 0.494), while the simple cycle's efficiency keeps rising to π ≈ 27 (η = 0.373). |
| **Assumptions and limits** | Air-standard cycle (air throughout; heat addition in place of combustion; heat rejection to close the loop), steady flow, adiabatic compressors and turbine of constant isentropic efficiency, pressure losses as fixed fractions, an intercooler that returns the air to a set temperature, and a recuperator described by an enthalpy-based effectiveness with equal flows on both sides. Constant-cp air (γ = 1.4) for the closed-form checks, the NASA dry-air model for the design numbers. Not modelled: turbine cooling air, fuel mass and combustion products (`P07.08`), mechanical losses, heat-exchanger size and weight, off-design operation, and the intercooler's coolant system. |
| **Analytical work** | The real simple-cycle work and efficiency with η_c and η_t; recuperator effectiveness and energy balance; the ideal recuperated and intercooled-recuperated efficiencies; the crossover pressure ratio; the best intercooler split; per-component entropy generation and its sum, the Clausius integral; the Gouy–Stodola form w_net = q_in(1 − T₀/T_H) − T₀ Σs_gen. |
| **Computational deliverable** | `starter/main.py`: `cycle_states`, `cycle_performance`, `component_balances`, `entropy_generation`, `ideal_recuperated_efficiency`, `ideal_icr_efficiency`, `simple_cycle_efficiency`, `best_split` (the NASA air model, the `P07.05` component functions, `CycleSpec`, `sweep` and `main` are given). Tests (18) prove: no recuperator and ideal components recover `P07.06` (efficiency and work); an ideal recuperator recovers 1 − (T₁/T₃)π^m, and with an intercooler the intercooled-recuperated closed form; the recuperator helps only below the crossover pressure ratio, where the turbine exit is hotter than the compressor exit; the first law closes for every component and overall to round-off; every component's entropy generation is non-negative and the total equals the Clausius integral; η < η_Carnot for every configuration; every loss lowers efficiency; intercooling alone raises work and lowers efficiency; √π is the best split; station pressures follow the losses; and the sweep has the predicted shape. |
| **Practical activity** | **computational experiment**: the four configurations at the design point; a pressure-ratio sweep of efficiency and specific work for each; a sensitivity study of the intercooled-recuperated cycle to component efficiencies, effectiveness, pressure losses, turbine inlet temperature and intercooler temperature; and the ideal limits checked against closed forms. |
| **Equipment** | None. A computer with the course's Python environment. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `real_brayton.png`, the four experiment tables, the derivations of the recuperated and intercooled-recuperated closed forms and of the √π split, and a half-page design note recommending a configuration for the core with the evidence. |
| **Assessment** | `assessments/phase-07/grading-rubric.md`; competencies 1–4. Pass: the closed forms are derived; all public tests pass; the simple real cycle at π = 12 is reproduced by hand with constant cp; the note recommends a configuration and a pressure ratio, names the largest sensitivity, and states what the model omits. |
| **References** | `@mit-ocw-16-050-thermal-energy` Part 2A 2.A.7 "Effect of Departures from Ideal Behavior - Real Cycle behavior", https://ocw.mit.edu/courses/16-050-thermal-energy-fall-2002/resources/07_part2a/ (required-free; compressor and turbine losses and the combustor pressure drop in the cycle). `@mit-16-unified-thermodynamics-propulsion-notes` 11.6 Performance of Jet Engines (11.6.5 Effect of Departures from Ideal Behavior), https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node131.html, and 18.5 Heat Exchangers (18.5.3 Efficiency of a Counterflow Heat Exchanger) (required-free; read online, link only). `@lienhard-heat-transfer-textbook` Ch. 3 "Heat exchanger design": 3.1 Function and configuration of heat exchangers; 3.3 Heat exchanger effectiveness, https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf (required-free; the effectiveness of a recuperator). `@libretexts-yan-engineering-thermodynamics` 6.11 "Applications of the second law of thermodynamics in open systems" (the two-stage compressor with an intercooler) (supplementary-free). |
| **Extensions** | ★ Add reheat: split the turbine, reheat to T₃ between the two halves, and show that intercooling plus reheat plus an ideal recuperator approaches Carnot as the number of stages grows (the Ericsson cycle). ★ Size the recuperator: with an overall heat transfer coefficient of your choice, use the ε–NTU relation for a counterflow exchanger to find the area for ε = 0.85, and how it grows as ε → 1. ★ Add turbine cooling: bleed a fraction of the compressor delivery air around the turbine inlet and find the bleed at which raising T₃ from 1400 to 1600 K stops paying. |
| **Time** | unmeasured |

## 1. The engineering problem

`P07.06` gave an ideal core at π = 12 and 1400 K an efficiency of 48 % (NASA air) and showed that the efficiency rises with pressure ratio forever. Real turboprops and APUs run far below that: the compressor and turbine are not isentropic, the combustor and ducts lose pressure, and because the compressor takes back about 40 % of the turbine's work, a few percent lost in each component is a large fraction of the net. The first job of this lesson is to put those losses into the cycle and see what is left.

The second job is a design decision. Two heat exchangers can be added. A **recuperator** (also called a regenerator) passes the hot turbine exhaust over the compressed air on its way to the combustor, so the combustor needs less heat for the same turbine inlet temperature. An **intercooler** cools the air between two compressors, so the second compressor works on denser, cooler air and needs less work. Some marine and armoured-vehicle gas turbines have used them; airliner engines use neither, because heat exchangers are heavy and bulky and because at their pressure ratios the exhaust is cooler than the compressed air. Which pays for a small core, and at what pressure ratio, is a question the model can answer for any set of component numbers.

The numbers are the illustrative core of `P07.06` (sea-level air, π = 12, turbine inlet temperature 1400 K) with typical, illustrative component values: isentropic efficiency 0.85 for each compressor and 0.88 for the turbine, a 4 % combustor pressure loss, an intercooler returning the air to 288.15 K with a 2 % loss, and a recuperator of effectiveness 0.85 losing 3 % on the cold side and 4 % on the hot side, at 5 kg/s of air. Success criteria, written first: the model recovers `P07.06` and the two closed forms of §2 in its ideal limits to 10⁻¹⁰; the first law closes for every component and the cycle to round-off; every component generates non-negative entropy and the total equals the Clausius integral to 10⁻¹⁰; efficiency stays below Carnot; every loss lowers efficiency.

## 2. Understand and derive

### 2.1 Stations and components

```
          intercooler (q_ic out)               recuperator (q_rec, inside the cycle)
    1  +-----+  2   +---+  3  +-----+  4     +-----------+  5   +------+  6  +-----+  7
 air ->| LPC |----->| IC|---->| HPC |------->| cold side |----->|  CC  |---->|  T  |---+
       +-----+      +---+     +-----+        |           |      +------+     +-----+   |
                                          8 <| hot side  |<--------------------------+
                                exhaust (q_ex out)  +-----------+
```

Station 1 is the LP compressor inlet; 2 the LP exit; 3 the intercooler exit; 4 the HP exit; 5 the recuperator's cold-side exit (combustor inlet); 6 the turbine inlet; 7 the turbine exit; 8 the exhaust. Without an intercooler the LP compressor does the whole pressure ratio and 3 = 4 = 2; without a recuperator 5 = 4 and 8 = 7. The **overall pressure ratio** is π = p₄/p₁. Each duct loses a fraction of its inlet pressure: p₃ = p₂(1 − Δ_ic), p₅ = p₄(1 − Δ_rc), p₆ = p₅(1 − Δ_cc). The turbine expands to whatever pressure lets the exhaust leave at p₁ after the hot side and the exhaust duct: p₇ = p₁/[(1 − Δ_rh)(1 − Δ_ex)]. So the turbine's expansion ratio p₆/p₇ is smaller than the compressor's π by the product of all the losses: at π = 12 with every loss present it is 10.73 instead of 12. That is how pressure losses cost work.

The components use `P07.05`: compressors h_out = h_in + (h_out,s − h_in)/η_c, the turbine h₇ = h₆ − η_t(h₆ − h₇ₛ), with the isentropic states from s°(T).

### 2.2 The real simple cycle

With constant cp, x = π^m, τ = T₃/T₁ (here T₃ means the turbine inlet temperature) and no pressure losses, the works are w_c = cpT₁(x − 1)/η_c and w_t = η_t cpT₃(1 − 1/x), so

**w_net = cpT₁[η_t τ(1 − 1/x) − (x − 1)/η_c],    q_in = cpT₁[τ − 1 − (x − 1)/η_c].**   (1)

The work now vanishes at x = η_cη_tτ rather than at x = τ, and the efficiency no longer rises without limit: at some π the extra compressor work outweighs the gain. With the illustrative components and NASA air, the simple cycle's efficiency peaks near π = 27 at 0.373 and its work near π = 12 at 298.5 kJ/kg (part 4). At the design point π = 12 it is η = 0.3407, against 0.4830 for the ideal cycle: 14 points lost to two component efficiencies and one pressure loss.

### 2.3 The recuperator

A recuperator transfers heat from the turbine exhaust (hot side, 7 → 8) to the compressed air (cold side, 4 → 5). The flows are equal (air standard), so the energy balance is

h₇ − h₈ = h₅ − h₄ = q_rec.   (2)

The most the cold side could gain is to be heated all the way to the hot inlet temperature T₇. The **effectiveness** is the fraction achieved,

**ε = (h₅ − h₄)/(h(T₇) − h₄),**   (3)

which for constant cp and equal flows is the familiar (T₅ − T₄)/(T₇ − T₄) (Lienhard's §3.3, `@lienhard-heat-transfer-textbook`, defines it for any exchanger through the capacity rates). With variable cp the enthalpy form keeps (2) exact. The combustor now only needs q_in = h₆ − h₅. Heat flows the right way only if **T₇ > T₄**: the turbine exhaust must be hotter than the compressor delivery air. At high pressure ratio T₄ is high and T₇ low; then (3) with ε > 0 *cools* the compressed air, and the "recuperator" makes things worse. `P07.06` showed that T₂ = T₄ (there, compressor exit equals turbine exit) is exactly the condition for maximum work. So:

**A recuperator helps only below the pressure ratio of maximum work of the corresponding simple cycle.**

**The ideal recuperated cycle.** With ε = 1, isentropic components, no losses and constant cp, T₅ = T₇ = T₃/x: the combustor adds q_in = cp(T₃ − T₃/x) = w_t. The efficiency is

**η = 1 − w_c/w_t = 1 − cpT₁(x − 1)/(cpT₃(1 − 1/x)) = 1 − (T₁/T₃)x = 1 − (T₁/T₃)π^((γ−1)/γ).**   (4)

It *falls* with π, the opposite of the simple cycle, and it equals the simple cycle's 1 − 1/x where x² = τ: the crossover pressure ratio is (T₃/T₁)^(γ/(2(γ−1))) = 15.90 at 1400 K. Below 15.9 the ideal recuperator helps; above it, it hurts. At π = 12, (4) gives 0.5814 against 0.5083 for the simple ideal cycle. With real components the crossover moves: the probe finds the real recuperated and real simple cycles equally efficient at π = 20.05 (NASA air), while T₇ = T₄ in the real simple cycle at π = 21.45; the recuperator's own pressure losses account for the difference.

### 2.4 The intercooler

Split the compression into an LP stage of ratio a and an HP stage of ratio b, with ab = π, and cool to T₁ between them. With constant cp and equal efficiencies, w_c = (cpT₁/η_c)[(a^m − 1) + (b^m − 1)]. Minimising a^m + (π/a)^m over a gives m a^(m−1) − m π^m a^(−m−1) = 0, so

**a = b = √π.**   (5)

Equal pressure ratios, equal work in each compressor. The saving at π = 12 (NASA air, real components, 2 % intercooler loss): the total compressor work falls from 351.10 to 292.98 kJ/kg, and the net work rises from 298.54 to 356.66 kJ/kg, +19 %. But the HP compressor now delivers air at 434.59 K instead of 630.58 K, so the combustor must add 1079.41 instead of 876.20 kJ/kg, and the **efficiency falls**, from 0.3407 to 0.3304. Intercooling buys work (a smaller engine for the same power), not efficiency, because the heat it rejects is extra heat that the combustor must replace. With a recuperator, the exhaust replaces it: the intercooled air is cold, T₇ − T₄ is large, and the recuperator has more to give. That is why intercooling and recuperation belong together.

**The ideal intercooled-recuperated cycle.** One intercooler to T₁ at √π, ε = 1, isentropic components, constant cp: w_c = 2cpT₁(x^½ − 1), w_t = cpT₃(1 − 1/x) = q_in, so

**η = 1 − 2(T₁/T₃)(x^½ − 1)/(1 − 1/x).**   (6)

At π = 12 and 1400 K: 0.6549, against 0.5814 recuperated and 0.5083 simple, and Carnot 0.7942. More intercoolers, and reheaters between turbine stages, push this towards Carnot (the Ericsson cycle, Extension 1).

### 2.5 The first and second laws for the whole cycle

Per kilogram of air: w_c = (h₂ − h₁) + (h₄ − h₃); w_t = h₆ − h₇; q_in = h₆ − h₅; heat rejected q_ic = h₂ − h₃ in the intercooler and q_ex = h₈ − h₁ to the atmosphere (the air-standard closure); and w_net = w_t − w_c. The first law for the cycle is

**q_in − q_ic − q_ex = w_net,**   (7)

and each component has its own balance (the recuperator's is (2)); the code checks all of them.

The second law: take the heat from a reservoir at T_H = T₆ (the hottest the gas gets) and reject all of it to the atmosphere at T₀ = T₁. Each component's entropy generation is its entropy rise plus the entropy it gives to its reservoirs: for the combustor s₆ − s₅ − q_in/T_H, for the intercooler s₃ − s₂ + q_ic/T₀, for the adiabatic machines s_out − s_in, for the recuperator the sum of its two sides. Summed around the cycle the stream entropies cancel, leaving

**Σ s_gen = (q_ic + q_ex)/T₀ − q_in/T_H ≥ 0,**   (8)

which is the Clausius integral of `P07.04` with its sign reversed. Multiply (8) by T₀ and use (7):

**w_net = q_in(1 − T₀/T_H) − T₀ Σ s_gen.**   (9)

The first term is the Carnot work of the heat; every component's T₀ s_gen is subtracted from it. That is the Gouy–Stodola theorem, and `P07.08` builds on it. For the intercooled-recuperated design point: q_in(1 − T₀/T_H) = 557.82 kJ/kg, T₀ Σ s_gen = 214.95 kJ/kg, w_net = 342.87 kJ/kg, and the sum closes exactly.

### 2.6 Where the model stops

- **Heat exchangers are heavy.** A recuperator of ε = 0.85 for 5 kg/s is a large, heavy component; its size grows steeply as ε → 1 (Extension 2). For an aircraft that weight often costs more fuel than the recuperator saves; for a tank, a ship or a ground APU it may not.
- **Turbine cooling.** Real turbines at 1400 K are cooled by compressor air that bypasses the combustor; the model has none.
- **Air standard.** The fuel adds mass and the products have a higher cp; the combustor is not a heat exchanger (`P07.08`).
- **Constant efficiencies.** Component efficiencies depend on size (small compressors are less efficient) and on operating point.
- **The intercooler needs a sink.** On an aircraft, cooling air must be ducted to it, with its own drag.

## 3. Build the computational model

Open `starter/main.py`. The NASA dry-air model, `Gas`, `constant_cp_gas`, the inversions, `isentropic_temperature`, `compressor_exit_temperature`, `turbine_exit_temperature` (from `P07.05`), the `CycleSpec` dataclass holding every input, `sweep`, the illustrative component sets `REAL`, `IC`, `REC`, `CONFIGS` and `main` are given.

### Step 1: the station table

`cycle_states(gas, spec)` fills T, p and h at the eight stations in the order the physics allows: compressors and intercooler first; then the pressures of the cold side, combustor, turbine exit and hot side from the losses; then the turbine (which does not depend on the recuperator); and only then the recuperator, from (3) and (2), because it needs T₇:

```python
h5 = h[3] + s.eps * (h[6] - h[3])        # (3): cold side gains a fraction of what is possible
h8 = h[6] - (h5 - h[3])                  # (2): the hot side loses exactly that
```

### Step 2: performance and the balances

`cycle_performance` returns the works, heats, η, the Carnot efficiency and the residual of (7). `component_balances` returns each component's energy residual. `entropy_generation(gas, spec)` returns each component's s_gen and their total, per §2.5.

### Step 3: closed forms and the split

`simple_cycle_efficiency` (`P07.06`), `ideal_recuperated_efficiency` (4) and `ideal_icr_efficiency` (6) are one line each. `best_split` minimises total compressor work over ln a with `minimize_scalar`.

Running `main()` prints

```
NASA air, pi = 12:
  configuration                 eta     w_net (kJ/kg)  power (MW)  q_in (kJ/kg)  q_rec   T4 (K)   T7 (K)   T8 (K)  residual (J/kg)
  simple                        0.3407     298.54       1.493       876.20     0.00   630.58   839.72   839.72   0.0e+00
  intercooled                   0.3304     356.66       1.783      1079.41     0.00   434.59   839.72   839.72   0.0e+00
  recuperated                   0.4238     284.74       1.424       671.90   204.30   630.58   852.16   664.52   0.0e+00
  intercooled + recuperated     0.4881     342.87       1.714       702.38   377.03   434.59   852.16   499.59   0.0e+00
```

then the station table of the intercooled-recuperated cycle, its component balances (all below 3 × 10⁻¹¹ J/kg) and entropy generation, the ideal limits, the sweep and the sensitivity study, and writes `real_brayton.png` (efficiency and work against π for the four configurations, with the Carnot line).

### Verify

`python -m pytest curriculum/phase-07/lesson-07-real-brayton-cycles` runs 18 tests:

1. **Lesson 06 recovered** (π = 4, 12, 30): ε = 0, no intercooler, ideal components: η = 1 − π^(−m) and w_net = cpT₁[τ(1 − 1/x) − (x − 1)] to 10⁻¹⁰.
2. **Recuperated closed forms** (π = 3, 8, 12): ε = 1 gives (4) to 10⁻¹⁰, with T₅ = T₇ and T₈ = T₄; adding the intercooler gives (6), which is higher.
3. **Crossover**: (4) equals the simple-cycle efficiency at (T₃/T₁)^(γ/(2(γ−1))); with ε = 0.8 the recuperator helps at π = 8 and hurts at π = 25, exactly as T₇ − T₄ changes sign.
4. **First law** (both gases): every component balance and (7) below 10⁻³ J/kg; q_rec equals the hot side's enthalpy drop.
5. **Second law** (both gases): every component's s_gen ≥ 0; the total equals (8) to 10⁻¹⁰; an ideal cycle generates entropy only in the combustor and the exhaust.
6. **Carnot**: 0 < η < 1 − T₁/T₃ for four configurations, three π, two gases.
7. **Losses**: each of nine degradations lowers η; real components lower both η and w_net.
8. **Intercooling**: raises work and lowers compressor work; lowers η without a recuperator, raises it with one.
9. **Split**: `best_split` returns √π for constant cp; an intercooler pressure loss moves it above √π.
10. **Pressures** follow the losses; **sweep** shapes: the recuperated cycle peaks at lower π, and the simple cycle's work has an interior maximum. **main** writes the figure.

### Compare

The ideal limits are the independent check: the model never contains (4) or (6), yet with ε = 1 and ideal components its state table reproduces them to 10⁻¹⁰ (0.581372 and 0.654905 at π = 12, printed by `main`). The entropy bookkeeping is a second, independent check on the energy bookkeeping: (9) holds only if both are right, and it closes to round-off.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Four configurations at π = 12.** The table in part 3. *Predict* the ranking by efficiency and by power before running. *Finding:* by efficiency, intercooled-recuperated (0.488) > recuperated (0.424) > simple (0.341) > intercooled (0.330); by power at 5 kg/s, intercooled (1.783 MW) > intercooled-recuperated (1.714) > simple (1.493) > recuperated (1.424). The recuperator costs power through its pressure losses; the intercooler buys power and costs efficiency, unless a recuperator is present.

2. **Pressure-ratio sweep.** π = 2 to 40 for each configuration (NASA air). *Predict* where each efficiency peaks. *Finding (selected rows; η, then w_net in kJ/kg):*

   | π | simple | intercooled | recuperated | intercooled + recuperated |
   |---|---|---|---|---|
   | 3 | 0.1856, 204.5 | 0.1816, 212.0 | 0.4308, 184.8 | 0.4382, 192.3 |
   | 5 | 0.2551, 262.4 | 0.2486, 282.6 | 0.4612, 245.2 | 0.4843, 265.3 |
   | 8 | 0.3066, 291.7 | 0.2976, 329.4 | 0.4515, 276.4 | 0.4942, 314.1 |
   | 12 | 0.3407, 298.5 | 0.3304, 356.7 | 0.4238, 284.7 | 0.4881, 342.9 |
   | 20 | 0.3677, 282.7 | 0.3594, 374.6 | 0.3682, 270.6 | 0.4671, 362.6 |
   | 30 | 0.3726, 250.5 | 0.3724, 376.6 | 0.3087, 239.7 | 0.4422, 365.8 |

   On a grid of step 1 in π: the simple cycle's efficiency peaks at 0.3731 near π = 27 and its work at 298.5 kJ/kg near π = 12; intercooled 0.3761 still rising at π = 40, work 377.1 kJ/kg near π = 26; recuperated 0.4612 near π = 5, work 284.7 near 12; intercooled-recuperated 0.4942 near π = 8, work 366.0 near π = 28. The recuperated cycles want low pressure ratios, which suits a small, few-stage core; the simple cycle wants a high one. The recuperated and simple cycles cross at π = 20.05.

3. **Sensitivity of the intercooled-recuperated cycle at π = 12.** Change one input at a time. *Predict* the largest effect. *Finding (base η = 0.4881):*

   | change | η | change in points |
   |---|---|---|
   | η_c 0.85 → 0.80 (both compressors) | 0.4630 | −2.52 |
   | η_t 0.88 → 0.83 | 0.4567 | −3.15 |
   | ε 0.85 → 0.75 | 0.4592 | −2.90 |
   | all pressure losses doubled | 0.4644 | −2.37 |
   | TIT 1400 → 1300 K | 0.4588 | −2.93 |
   | intercooler exit T₁ + 20 K | 0.4767 | −1.15 |

   The turbine efficiency matters most per five points, as the back-work argument predicts (the turbine's work is the larger of the two); recuperator effectiveness and turbine inlet temperature come next.

4. **The ideal limits.** Constant cp, π = 12, 1400 K. *Finding:* ε = 0 and ideal components: 0.508343 = `P07.06`'s 0.508343; ε = 1: 0.581372 = (4); ε = 1 with intercooling: 0.654905 = (6); Carnot 0.794179.

Raw data format: CSV tables `config, eta, w_net_kJkg, power_MW, q_in, q_rec, T4, T7, T8` (1); `pi, config, eta, w_net_kJkg` (2); `change, eta, delta_points` (3); `limit, eta_model, eta_closed_form` (4); plus `real_brayton.png`.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| ideal simple cycle | `P07.06`: 0.508343 | 0.508343 | < 10⁻¹⁰ | identity |
| ideal recuperated | (4): 0.581372 | 0.581372 | < 10⁻¹⁰ | identity |
| ideal intercooled-recuperated | (6): 0.654905 | 0.654905 | < 10⁻¹⁰ | identity |
| cycle first law (7) | 0 | 0.0 J/kg (all four configurations) | round-off | — |
| component balances | 0 | ≤ 2.9 × 10⁻¹¹ J/kg | round-off | Newton inversions to 10⁻¹⁰ K |
| Σ s_gen against (8) | 745.97 J/(kg·K) | 745.97 J/(kg·K) | < 10⁻¹⁰ relative | identity |
| (9): q_in(1 − T₀/T_H) − T₀Σs_gen | 557.82 − 214.95 kJ/kg | w_net = 342.87 kJ/kg | round-off | Gouy–Stodola |
| best split, constant cp, no loss | √π | `best_split` | < 10⁻⁵ | optimiser tolerance |
| best split, NASA air, 2 % loss | above √π = 3.4641 | 3.4993 | — | pressure loss in the intercooler |

The engineering answer: **at π = 12 and 1400 K the intercooled-recuperated core is 48.8 % efficient and gives 1.71 MW from 5 kg/s; the simple core gives 34.1 % and 1.49 MW. For this component set the best efficiency is 49.4 %, at π ≈ 8, where fewer compressor stages are needed; the recuperator alone gives 46.1 % at π ≈ 5. A recuperator stops paying above π ≈ 20, which is why high-pressure-ratio aero engines do not carry one.** The largest uncertainty is the turbine efficiency (3 points per 5 points of η_t), then the recuperator effectiveness and the turbine inlet temperature.

What to improve: add turbine cooling air (Extension 3), the fuel mass and products' properties (`P07.08`), and the recuperator's size and weight (Extension 2), which decide whether the 14-point efficiency gain is worth carrying.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1), (4), (5), (6), (8) and (9); `real_brayton.png`; the four experiment tables; and a half-page design note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6, 7) recommending a configuration and pressure ratio for the core, with the evidence and the omissions.

### Adversarial self-checks

- **What if ε → 0?** The recuperator does nothing but still loses pressure, so the cycle is slightly worse than the simple one; with its losses set to zero it is exactly the simple cycle (test 1).
- **What if π rises to 30 with the recuperator fitted?** T₇ < T₄ and the recuperator heats the exhaust with the compressed air: η = 0.3087 against 0.3726 without it.
- **What if the components become ideal?** With ε = 1 the closed forms (4) and (6) are recovered. For the intercooled-recuperated cycle with NASA air, removing the pressure losses alone lifts η from 0.488 to 0.509; ideal compressors and turbine as well, with ε still 0.85, give 0.625; ε = 1 gives 0.670.
- **Which assumption fails first in the real world?** Turbine inlet temperature without cooling air: at 1400 K a real turbine bleeds several percent of the compressor flow for cooling, which costs work and efficiency the model does not see.
- **What measurement could falsify the model?** Measured station temperatures and pressures of a recuperated engine: the recuperator's energy balance (2) must close with the measured flows, and its measured effectiveness, entered into the model with the measured component efficiencies, must reproduce the measured shaft power within the uncertainty of the measurements.
- **Does the result survive a change of gas model?** The ranking of the four configurations at π = 12 is the same with constant cp (0.477, 0.407, 0.353, 0.340) and with NASA air (0.488, 0.424, 0.341, 0.330); the efficiencies differ by one to two points.
- **Which decision changes if the recuperator effectiveness is 0.75 instead of 0.85?** At π = 12, η falls 2.9 points to 0.459; the best intercooled-recuperated efficiency falls from 0.494 at π = 8 to 0.459 at π = 11.5, and the recuperated cycle's from 0.462 at π = 5.5 to 0.425 at π = 7. Whether the recuperator is worth its weight on an aircraft becomes doubtful.

## Check yourself

1. Define recuperator effectiveness. What is the cold-side exit temperature for ε = 1?
2. Why does a recuperator help only below a certain pressure ratio? What is that pressure ratio for an ideal cycle?
3. Write the ideal recuperated-cycle efficiency. How does it change with π, and why is that the opposite of the simple cycle?
4. Why does intercooling without recuperation lower the efficiency, even though it saves compressor work?
5. What pressure split between two compressors minimises the total work, and under what conditions?
6. How do duct pressure losses reduce the work of a gas turbine, given that they involve no friction in the compressor or turbine?
7. State the relation between net work, the Carnot work of the heat supplied, and the entropy generated in the components.
8. Which component efficiency does the intercooled-recuperated cycle's efficiency depend on most, and why?
