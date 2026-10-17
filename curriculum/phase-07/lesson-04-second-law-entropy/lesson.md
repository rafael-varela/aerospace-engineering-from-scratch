# P07.04  The Second Law & Entropy

> Energy is conserved; its usefulness is not. Every real process generates entropy, and the entropy generated, times the temperature of the surroundings, is work you can never get back. Count it, process by process, and you know where to spend design effort.

| | |
|---|---|
| **Learning outcomes** | 1. State the second law in the Kelvin–Planck and Clausius forms and as S_gen ≥ 0, and derive the Clausius inequality ∮δQ/T ≤ 0 and the Carnot bound η ≤ 1 − T_C/T_H. 2. Derive the Gibbs relations T ds = du + p dv = dh − v dp, show that entropy is a property, and compute Δs of an ideal gas with variable cp. 3. Write entropy balances for closed systems and control volumes and compute S_gen for heat across a finite temperature difference, two bodies equilibrating, throttling, duct friction, a two-stream heat exchanger, mixing and liquid-line friction. 4. Implement and test an entropy-generation tracker whose every entry is non-negative, that rejects impossible cycles, and whose entropy changes are path-independent. 5. Rank the irreversibilities of an aircraft bleed-air system and say which one to attack first. |
| **Prerequisites** | `P07.01` (ideal-gas properties, s°(T), the NASA air model), `P07.02` (first law for closed systems and control volumes), `P07.03` (h_fg = T s_fg), `P03.04` (numerical integration, trapezoid rule and its order), `P01.05` (integration). Skill check: you can integrate cp(T)/T numerically and explain why the trapezoid rule's error falls as h². |
| **Engineering question** | An airliner's bleed-air system takes 1 kg/s of air from the compressor at 1.0 MPa and 700 K, throttles it in a pressure-regulating valve to 0.30 MPa, cools it to 470 K in a precooler against 1.2 kg/s of 320 K fan air, loses 20 kPa in the ducts, and finally mixes 0.5 kg/s of conditioned air at 270 K with 0.5 kg/s of recirculated cabin air at 297 K (typical, illustrative numbers). Where is entropy generated, and how much work does each step throw away? The answer: 535.8 W/K in all, which at T₀ = 288.15 K is 154.4 kW of lost work; the valve generates 345.6 W/K (64.5 %), the precooler 169.3 W/K (31.6 %), the duct 19.8 W/K (3.7 %) and the mixer 1.1 W/K (0.2 %). The throttle, not the heat exchanger, is the place to start. |
| **Assumptions and limits** | Dry air as an ideal gas with the NASA variable cp of `P07.01`, valid 200–1000 K; steady flow; each device adiabatic to the outside except where heat crosses its boundary at a stated temperature; heat exchanger and mixer without pressure loss; reservoirs at uniform, fixed temperatures; incompressible liquid with constant c for the fuel line. The lost-work figure T₀ S_gen needs a dead-state temperature T₀ (here 288.15 K); its proper development is lesson 08 (exergy). The model fails where the gas leaves its validity range, where kinetic energy matters (a throttle at high Mach), and when heat-transfer temperatures are not uniform. |
| **Analytical work** | The Clausius inequality and the definition dS = δQ_rev/T; the Gibbs relations; Δs of an ideal gas; entropy balances; S_gen of the seven processes in closed form, including its small-ΔT limit (linear for fixed Q, quadratic for two bodies); the Carnot bound; the ranking of the bleed-air chain. |
| **Computational deliverable** | `starter/main.py`: `ds_ideal_gas`, `ds_path_integral`, `heat_transfer_sgen`, `bodies_equilibrate`, `throttle_sgen`, `heat_exchanger`, `mix_streams`, `mix_gases_sgen`, `liquid_line_friction`, `clausius_integral`, `carnot_efficiency`, `check_cycle`, and the `EntropyLedger` methods `add`, `total`, `lost_work`, `ranked`. The NASA dry-air helpers (`air_cp`, `air_h`, `air_s0`, `air_T_from_h`), the bleed-air case and `bleed_air_ledger`, and `main` are given. Tests (15) prove S_gen ≥ 0 for 200 seeded random cases of every process, S_gen → 0 as ΔT → 0 with slopes 1 and 2, path independence along three paths and a closed loop, second-order convergence of the path integral, the equality of the irreversible and reversible routes, the exchanger and mixer energy balances, the refusal of a temperature cross, the ledger's ranking and refusal of negative entries, and that every one of 2000 random cycles the tracker accepts obeys the Carnot bound while every cycle that beats it is rejected. |
| **Practical activity** | **computational experiment**: the small-ΔT limit, path-independence and convergence study, the bleed-air ledger with a sensitivity sweep of the valve outlet pressure and the fan-air flow, and a cycle-claim audit. |
| **Equipment** | None (a computer with the course's Python environment). |
| **Safety** | None beyond normal computer use. |
| **Evidence** | `entropy_generation.png`; the printed checks; the experiment tables; a half-page note answering the engineering question with the ranked ledger, the lost work and the change you would make first. |
| **Assessment** | `assessments/phase-07/grading-rubric.md`; competencies 1–4. Pass: the Clausius inequality, the Gibbs relations and the Carnot bound derived; all public tests pass; the valve's S_gen and the two-block S_gen reproduced by hand; the note ranks the irreversibilities and defends the first design change. |
| **References** | `@libretexts-yan-engineering-thermodynamics` Ch. 6 "Entropy and the Second Law of Thermodynamics": 6.2 Heat engine; 6.4 The second law of thermodynamics- Kelvin-Planck and Clausius statements; 6.5 Carnot cycles; 6.6 Entropy and entropy generation (6.5.1 The inequality of Clausius; 6.5.2 Definition of entropy; 6.5.3 Entropy generation, Sgen); 6.7 The second law of thermodynamics for closed systems; 6.10 The second law of thermodynamics for open systems; 6.11 Applications of the second law of thermodynamics in open systems (required-free; https://eng.libretexts.org/Bookshelves/Mechanical_Engineering/Introduction_to_Engineering_Thermodynamics_(Yan)). `@mit-16-unified-thermodynamics-propulsion-notes` 4.1 Reversibility and Irreversibility in Natural Processes; 5.1 Concept and Statements of the Second Law; 5.3 Combined First and Second Law Expressions; 5.4 Entropy Changes in an Ideal Gas; 5.5 Calculation of Entropy Change in Some Basic Processes (required-free, read online, https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node39.html). `@mit-ocw-16-050-thermal-energy` Part 1C "Applications of the Second Law", 1.C.7 Examples of Lost Work in Engineering Processes (required-free; https://ocw.mit.edu/courses/16-050-thermal-energy-fall-2002/resources/05_part1c/). `@openstax-university-physics-vol2` 4.1 Reversible and Irreversible Processes; 4.4 Statements of the Second Law of Thermodynamics; 4.5 The Carnot Cycle; 4.6 Entropy (required-free). |
| **Extensions** | ★ Replace the pressure-regulating valve by an ideal expansion turbine driving the fan (the air-cycle machine idea) and compute the work recovered and the entropy it no longer generates. ★ Add a pressure loss to the heat exchanger and find the flow arrangement that minimises total S_gen for a fixed duty (the trade between heat-transfer and friction irreversibility). ★ Use the water table of `P07.03` to compute S_gen of throttling saturated liquid water from 1 MPa to 100 kPa (flashing), and compare with an ideal gas throttled over the same ratio. |
| **Time** | unmeasured |

## 1. The engineering problem

A turbofan does not only make thrust: it also supplies hot, high-pressure air to the aircraft for cabin pressurisation, air conditioning and wing anti-icing. This **bleed air** is taken from the compressor at about 1 MPa and 700 K (typical, illustrative values), far hotter and at higher pressure than the cabin needs. So the system immediately destroys what the engine worked to create: a valve throttles the pressure down, a heat exchanger (the precooler) dumps heat into fan air, ducts lose more pressure, and a mixer blends the conditioned air with recirculated cabin air. The first law balances every step; energy is not lost. Yet the fuel spent compressing that air is real, and the engine's designers want to know where its value goes.

The second law answers that. Every real process generates **entropy**, and the entropy generated, S_gen, is the bookkeeping quantity that measures irreversibility: zero for an ideal reversible process, positive for every real one, never negative. Multiplied by the temperature of the environment, T₀ S_gen is the work that a perfect machine could have extracted and that is now gone (lesson 08 develops this as exergy). This lesson builds a tracker that computes S_gen for each kind of irreversible process, refuses any negative result (a second-law violation is a bug or an impossible claim), and ranks the entries.

Success criteria, written first: S_gen ≥ 0 for every process over hundreds of random inputs; S_gen → 0 as the temperature difference → 0; the entropy change of the gas between two states is the same along every path, to the integration error, with the error falling at second order; the exchanger and mixer energy balances close; and every cycle the tracker accepts has an efficiency below the Carnot limit for its reservoirs.

## 2. Understand and derive

### 2.1 Statements of the second law

Experience shows two asymmetries that the first law does not forbid:

- **Kelvin–Planck:** no cyclic device can take heat from a single reservoir and convert all of it into work. Every heat engine must reject some heat to a colder reservoir.
- **Clausius:** heat does not flow on its own from a colder body to a hotter one. A refrigerator needs work.

The two statements are equivalent (a violation of one lets you build a violation of the other). A **reservoir** is a body so large that heat exchange does not change its temperature. A **reversible** process is one that can be undone leaving no change anywhere; friction, heat across a finite temperature difference, unrestrained expansion and mixing are the classic **irreversibilities**, and they are exactly the processes this lesson tracks.

### 2.2 The Carnot bound

Consider any heat engine that takes Q_H from a reservoir at T_H and rejects Q_C to one at T_C. A reversible engine (a Carnot engine) can be run backwards as a refrigerator. If any engine were more efficient than a reversible one between the same reservoirs, it could drive the reversed reversible engine and produce net work from the cold reservoir alone, violating Kelvin–Planck. So all reversible engines between T_H and T_C have the same efficiency, which can depend only on the two temperatures, and the absolute (thermodynamic) temperature scale is defined so that for them Q_C/Q_H = T_C/T_H. Hence

**η = W/Q_H ≤ η_Carnot = 1 − T_C/T_H,**

with equality only for a reversible engine.

### 2.3 The Clausius inequality and the definition of entropy

For a reversible two-reservoir cycle, Q_H/T_H − Q_C/T_C = 0; for an irreversible one, Q_C is larger, so Q_H/T_H − Q_C/T_C < 0. Writing heat into the cycle as positive and generalising to any number of exchanges at boundary temperatures T_k (each piece of a cycle can be built from small Carnot cycles),

**∮ δQ/T ≤ 0  (Clausius inequality),**

with equality for a reversible cycle. A quantity whose integral around every reversible cycle is zero is a property: so define

**dS = δQ_rev/T.**

This is entropy. Between two states 1 and 2 joined by a reversible path R and returning by any path, Clausius gives ∫₁² δQ/T ≤ S₂ − S₁, with equality if the path is reversible. The difference is the **entropy generated** inside the system:

**S₂ − S₁ = ∫₁² δQ/T + S_gen,  S_gen ≥ 0.**

That is the second law as a balance: entropy changes because it is carried in with heat (δQ/T, evaluated at the boundary temperature) and because irreversibilities inside create it. It is never destroyed. For a control volume in steady flow, mass carries entropy too:

**0 = Σ_k Q̇_k/T_k + Σ_in ṁ s − Σ_out ṁ s + Ṡ_gen.**

### 2.4 The Gibbs relations and Δs of an ideal gas

For a reversible process of a simple compressible system, δQ_rev = T ds and δW_rev = p dv; the first law du = δq − δw then gives

**T ds = du + p dv,  and with h = u + pv,  T ds = dh − v dp.**

Both relate properties only, so they hold for any process between the same states, reversible or not; that is how entropy changes are computed for irreversible processes. For an ideal gas, dh = cp dT and v = RT/p, so

**ds = cp(T) dT/T − R dp/p,  s₂ − s₁ = s°(T₂) − s°(T₁) − R ln(p₂/p₁),**

with s°(T) the NASA function of `P07.01`. For air from (300 K, 100 kPa) to (800 K, 1 MPa) the code gives **355.342 J/(kg·K)**. Because s is a property, integrating cp dT/T − R dp/p along any path between those states gives the same answer; Section 3 checks it along three very different paths.

### 2.5 The seven entropy generators

Each formula follows from the entropy balance with the first law.

**Heat across a finite temperature difference.** Heat Q leaves a reservoir at T_H (its entropy falls by Q/T_H) and enters one at T_C (rises by Q/T_C): **S_gen = Q(1/T_C − 1/T_H) = Q ΔT/(T_H T_C)**. For fixed Q it vanishes linearly as ΔT → 0. With Q = 1 kJ into 300 K: 0.833 J/K at ΔT = 100 K, 0.0111 J/K at 1 K, 1.11 × 10⁻³ J/K at 0.1 K.

**Two bodies equilibrating.** Masses with heat capacities C₁ = m₁c₁, C₂ = m₂c₂ at T₁, T₂ in an insulated box reach T_f = (C₁T₁ + C₂T₂)/(C₁ + C₂), and **S_gen = C₁ ln(T_f/T₁) + C₂ ln(T_f/T₂)**. For equal bodies at T ± ΔT/2 this is −C ln(1 − (ΔT/2T)²) ≈ C ΔT²/(4T²): **quadratic** in ΔT. Two 1 kg aluminium blocks (c = 900 J/(kg·K), a typical value) at 300 ± 5 K generate 0.250 J/K; at 300 ± 0.05 K, 2.5 × 10⁻⁵ J/K. Here the heat exchanged also shrinks with ΔT, which is why the order is two.

**Throttling an ideal gas.** A valve is adiabatic with no work, so h₂ = h₁ (`P07.02`); for an ideal gas h depends on T only, so T₂ = T₁, and the entropy balance gives **s_gen = s₂ − s₁ = R ln(p₁/p₂)**. The bleed valve, 1.0 → 0.30 MPa: 287.045 × ln(3.333) = **345.6 J/(kg·K)**. Adiabatic duct friction at low Mach number is the same formula: 0.30 → 0.28 MPa gives **19.8 J/(kg·K)**.

**A two-stream heat exchanger.** Adiabatic outside, no pressure loss, steady: the energy balance ṁ_h(h_h,in − h_h,out) = ṁ_c(h_c,out − h_c,in) fixes the cold outlet, and **Ṡ_gen = ṁ_h(s_h,out − s_h,in) + ṁ_c(s_c,out − s_c,in)**: the hot stream loses entropy, the cold one gains more, because it receives the same heat at lower temperatures. The precooler: Q̇ = 241.1 kW, fan air leaves at 517.52 K, **Ṡ_gen = 169.3 W/K**. A counterflow exchanger cannot make the hot outlet colder than the cold inlet or the cold outlet hotter than the hot inlet; the code refuses such a "temperature cross".

**Adiabatic mixing of streams.** The energy balance Σṁ_i h_i = (Σṁ_i) h_out fixes T_out, and **Ṡ_gen = Σ ṁ_i [s(T_out, p_out) − s(T_i, p_i)]**. Mixing 0.5 kg/s at 270 K with 0.5 kg/s at 297 K: T_out = 283.50 K, Ṡ_gen = **1.14 W/K**. Small, because the temperatures are close: it scales like the two-body case, quadratically.

**Mixing different gases.** Different ideal gases initially separate at the same T and p each expand into the whole volume: **S_gen = −R_U Σ n_i ln x_i**. 0.7808 mol N₂, 0.2095 mol O₂ and 0.0093 mol Ar: **4.687 J/K**. This is the mixing entropy of `P07.01`, now recognised as generated.

**Friction in a liquid line.** An incompressible liquid (h = cT + p/ρ) flowing adiabatically through a filter or line with pressure drop Δp keeps h constant, so it warms by Δp/(ρc) and **Ṡ_gen = ṁ c ln(T_out/T_in) ≈ ṁ Δp/(ρ T)**. A fuel line with 2 kg/s, Δp = 150 kPa, ρ = 800 kg/m³ and c = 2000 J/(kg·K) (typical, illustrative values for a kerosene-type fuel): ΔT = 0.094 K, Ṡ_gen = **1.250 W/K**, a lost work T₀Ṡ_gen of 360 W against a pumping power ṁΔp/ρ of 375 W. Almost all of the pump's work into the line is destroyed by friction.

### 2.6 Cycles and the tracker's acceptance rule

For a cycle exchanging heat Q_k (positive into the cycle) with reservoirs at T_k, the first law gives the net work W = ΣQ_k and the entropy balance over a cycle (the working fluid returns to its state) gives **S_gen = −Σ Q_k/T_k ≥ 0**: the Clausius inequality. The tracker's `check_cycle` computes S_gen and **rejects** the cycle if it is negative. For two reservoirs, S_gen ≥ 0 is algebraically the same as η ≤ 1 − T_C/T_H; for several, η is still below the Carnot efficiency between the hottest and coldest reservoirs used. A claimed Stirling converter taking 250 W at 1050 K and rejecting at 350 K (typical, illustrative; Carnot 0.6667): 100 W out is accepted (η = 0.400, Ṡ_gen = 0.1905 W/K); 166 W is accepted barely (η = 0.664, 0.0019 W/K); 180 W is rejected (Ṡ_gen would be −0.0381 W/K).

### 2.7 Lost work

Lesson 08 shows that an irreversibility generating S_gen in surroundings at T₀ destroys exactly T₀ S_gen of work potential (the Gouy–Stodola theorem). Here it is used as a ruler: at T₀ = 288.15 K the bleed valve's 345.6 W/K is 99.6 kW. For scale, the bleed air's enthalpy above T₀ is 425.1 kW. T₀ is a choice: at cruise the outside air is colder, which changes the lost-work figures but not the ranking.

### 2.8 Where the model fails

- **The air model's range.** The NASA lower interval covers 200–1000 K; hotter bleed or combustor states need the upper interval (`P07.01`).
- **Kinetic energy.** In a valve with a sonic jet, h₂ = h₁ holds between stagnation states, not across the jet; the formula s_gen = R ln(p₁/p₂) applies between the upstream and fully mixed downstream states.
- **Non-uniform boundary temperatures.** Q/T assumes heat crosses at one temperature; a real heat-exchanger wall has a temperature distribution, so S_gen inside the wall and boundary layers must be integrated (lessons 12–13).
- **Real fluids.** Throttling a liquid near saturation flashes it; the ideal-gas formula does not apply (Extension 3).

## 3. Build the computational model

Open `starter/main.py`. The NASA dry-air helpers from `P07.01` (lower interval), the bleed-air case, `bleed_air_ledger` and `main` are given.

### Step 1: entropy as a property

`ds_ideal_gas(T1, p1, T2, p2)` is the end-state formula of 2.4. `ds_path_integral(T_path, p_path, n)` samples a path λ ↦ (T(λ), p(λ)) on n points and applies the trapezoid rule to cp(T)/T dT and −R dp/p:

```python
lam = np.linspace(0.0, 1.0, n)
T, p = T_path(lam), p_path(lam)
f_T = air_cp(T) / T
return np.sum(0.5 * (f_T[1:] + f_T[:-1]) * np.diff(T) - R_AIR * 0.5 * (1/p[1:] + 1/p[:-1]) * np.diff(p))
```

### Step 2: the generators

`heat_transfer_sgen`, `bodies_equilibrate`, `throttle_sgen` (raise `ValueError` if p₂ > p₁), `heat_exchanger` (energy balance with `air_T_from_h`, then the two streams' Δs, then the temperature-cross check), `mix_streams`, `mix_gases_sgen` and `liquid_line_friction` are the formulas of 2.5.

### Step 3: cycles and the ledger

`clausius_integral`, `carnot_efficiency` and `check_cycle` are 2.6. `EntropyLedger.add` appends an entry and raises `ValueError` for S_gen < −tol (and does not record it); `total`, `lost_work` and `ranked` (largest first, with fractions of the total) are one line each.

Running `main()` prints

```
entropy is a property: air from (300 K, 100 kPa) to (800 K, 1 MPa), J/(kg K)
  end states:                      355.342061
  isobaric then isothermal         355.342044
  isothermal then isobaric         355.342044
  wiggly path                      355.342063
  wiggly-path trapezoid error, n = 101..801: 9.78e-02, 2.45e-02, 6.12e-03, 1.53e-03; observed orders [2. 2. 2.]
bleed-air chain (typical, illustrative): precooler Q = 241.1 kW, fan air leaves at 517.52 K; mix manifold T_out = 283.50 K
process                                     S_gen (W/K)   share  T0 S_gen (kW)
pressure-regulating valve (throttling)           345.59   64.5%          99.58
precooler (heat across a finite dT)              169.27   31.6%          48.78
duct friction                                     19.80    3.7%           5.71
mix manifold (mixing)                              1.14    0.2%           0.33
total                                            535.81                 154.39
```

(abridged) followed by the small-ΔT table, the fuel line, the gas mixing and the cycle claims, and writes `entropy_generation.png`: the ranked ledger as bars; S_gen against ΔT on log–log axes with slopes 1 and 2; and the three paths in the (T, p) plane.

### Verify

`python -m pytest curriculum/phase-07/lesson-04-second-law-entropy` runs 15 tests:

1. **Property**: the end-state formula reduces to −R ln(p₂/p₁) and cp ln(T₂/T₁) in its limits and is antisymmetric; three paths agree with it to 10⁻⁶ and a closed loop gives zero; the trapezoid error converges at order 2.0 ± 0.2; the throttle's Δs equals the heat-over-temperature of a reversible isothermal expansion between the same states.
2. **Second law**: 200 seeded random cases of each generator give S_gen ≥ 0; heat across ΔT and two bodies give slopes 1 and 2 on log–log axes, with the two-body value C ΔT²/(4T²) at small ΔT; zero difference gives zero; a throttle refuses a pressure rise.
3. **Balances**: the exchanger's energy balance to 10⁻¹⁰, a warmer cold inlet lowering S_gen and more cold flow raising it, a temperature cross refused; the mixer's energy balance, equal-temperature mixing generating nothing, R_U ln 2 per mole for equal halves of two gases; the liquid line against ṁΔp/(ρT).
4. **Ledger**: totals, lost work, ranking, refusal of a negative entry without recording it; the bleed-air ledger complete, positive and led by the valve.
5. **Cycles**: a Carnot cycle has a zero Clausius integral and the Carnot efficiency; 2000 random two-reservoir cycles (seeded) are accepted exactly when η is below Carnot, and every accepted one obeys the bound; a four-reservoir cycle is bounded by Carnot between its extremes.

### Compare

Path independence is the independent check of the property calculation: the end-state formula (NASA s° differences) and numerical integration along three paths agree to 2 × 10⁻⁵ J/(kg·K) at 20 001 points, and the wiggly path's error falls by exactly a factor of 4 per halving of the step (observed orders 2.00). The throttle provides a second: its entropy generation, computed from the balance, equals the entropy that a reversible isothermal expansion between the same states would receive as Q/T. The cycle checker is tested against the Carnot formula on 2000 random cycles. No result depends on one formula alone.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before you run it.

1. **The small-ΔT limit.** For Q = 1 kJ into 300 K and for two 1 kg aluminium blocks centred on 300 K, compute S_gen at ΔT = 100, 10, 1 and 0.1 K. *Predict* the slopes. *Finding:* heat across ΔT: 0.833, 0.1075, 0.01107, 1.111 × 10⁻³ J/K (slope 1); two blocks: 25.35, 0.250, 2.50 × 10⁻³, 2.50 × 10⁻⁵ J/K (slope 2). Reversible heat transfer is the ΔT → 0 limit, which is why it would take an infinitely large heat exchanger.
2. **Paths and convergence.** Integrate Δs from (300 K, 100 kPa) to (800 K, 1 MPa) along the three paths of `main()` and a closed loop. *Predict* the error at n = 101 and how it scales. *Finding:* 355.342061 (end states) against 355.342044, 355.342044 and 355.342063 at n = 20 001; the wiggly path's error is 9.78 × 10⁻², 2.45 × 10⁻², 6.12 × 10⁻³, 1.53 × 10⁻³ at n = 101, 201, 401, 801.
3. **The valve outlet pressure.** Sweep the valve outlet from 0.20 to 0.60 MPa with the rest of the chain fixed. Keep the 20 kPa duct loss. *Predict* how the valve's share of the total changes. *Finding:* at 0.20, 0.30, 0.40 and 0.60 MPa the valve generates 461.98, 345.59, 263.02 and 146.63 W/K (69.7, 64.5, 58.7 and 44.9 % of totals of 662.64, 535.81, 448.15 and 326.77 W/K; lost work 190.9, 154.4, 129.1 and 94.2 kW). At 0.60 MPa the precooler (169.3 W/K, unchanged) overtakes the valve: once the throttling ratio is small, heat transfer across the precooler's large temperature difference becomes the next target. (Problem 5 asks you to also move the precooler outlet.)
4. **The fan-air flow.** With the bleed cooled from 700 to 465 K, vary the fan-air flow. *Predict* whether more fan air raises or lowers the precooler's S_gen. *Finding:* 0.8 kg/s: fan air leaves at 619.87 K, Ṡ_gen = 113.6 W/K; 1.2 kg/s: 521.65 K, 168.2 W/K; 2.0 kg/s: 441.68 K, 224.0 W/K. More cooling air means a colder cold stream, a larger temperature difference and more entropy: the opposite of the intuition that "more coolant is better". The limit is the temperature cross.
5. **Cycle claims.** Audit the four Stirling-converter claims of 2.6 and two more: Q = (300, −180) J at (900, 300) K, and Q = (500, 200, −450) J at (1200, 800, 350) K. *Finding:* the first is accepted with η = 0.400 against a Carnot 0.667 (S_gen = 0.267 J/K); the second, η = 0.357 against 0.708, S_gen = 0.619 J/K; Q = (300, −90) J at (900, 300) K is rejected (S_gen = −0.0333 J/K).

Raw data format: CSV tables `dT_K, sgen_heat_JK, sgen_blocks_JK` (1); `path, n, ds_JkgK, error` (2); `p_prv_MPa, process, sgen_WK, share` (3); `mdot_fan_kgs, T_c_out_K, sgen_WK` (4); `claim, W, eta, eta_carnot, sgen, verdict` (5); plus `entropy_generation.png`.

## 5. Validate and improve

| quantity | independent value | code | difference | cause |
|---|---|---|---|---|
| Δs, air 300 K/100 kPa → 800 K/1 MPa | path integrals | end states 355.342061 J/(kg·K) | ≤ 1.7 × 10⁻⁵ | trapezoid at n = 20 001 |
| path-integral convergence | order 2 (trapezoid) | 2.00, 2.00, 2.00 | — | — |
| two blocks, ΔT = 10 K | C ΔT²/(4T²) = 0.2500 J/K | 0.2500 J/K | < 10⁻³ relative | higher-order terms |
| valve | R ln(p₁/p₂) | 345.59 W/K | — | identity |
| fuel line | ṁΔp/(ρT) = 1.2500 W/K | 1.2498 W/K | 0.02 % | ln expansion |
| Carnot cycle | Clausius integral 0 | < 10⁻¹² | — | — |
| 2000 random cycles | η ≤ η_Carnot if accepted | all | — | — |

The engineering answer: **the bleed-air chain generates 535.8 W/K, 154.4 kW of lost work at T₀ = 288.15 K. The pressure-regulating valve is responsible for 64.5 % (345.6 W/K), the precooler for 31.6 % (169.3 W/K), the ducts for 3.7 % and the mixer for 0.2 %.** Attack the valve first: either take bleed air from a lower compressor stage, closer to the pressure needed, or expand it through a turbine that does work (the air-cycle machine of real environmental-control packs, and the electric "bleedless" architecture of some recent airliners avoids the bleed altogether). Polishing the mixer would gain nothing. What to improve in the model: the precooler's own pressure losses and the choice of T₀ for an aircraft at cruise (lesson 08).

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of the Clausius inequality, the Gibbs relations, the Carnot bound and the seven S_gen formulas; `entropy_generation.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question.

### Adversarial self-checks

- **What if the pressure ratio of the valve doubles?** Its S_gen rises by R ln 2 = 199 J/(kg·K) per kilogram: logarithmic, so halving the throttling ratio is worth less than it sounds.
- **What if ΔT goes to zero?** Every heat-transfer entry goes to zero (linearly or quadratically), and the heat exchanger needs infinite area: reversibility is a limit, not a design.
- **Which assumption fails first?** The uniform boundary temperature in the precooler, and the dead-state temperature, which changes the lost-work figures with altitude (but not the ranking).
- **What measurement could falsify the model?** Measured inlet and outlet p and T of each component give s changes directly; a component whose measured S_gen came out negative would mean a wrong measurement or an unaccounted heat path, never a real second-law violation.
- **Does the result survive a change of method?** Path integrals against end states; the tracker's acceptance against the Carnot formula; S_gen against its small-ΔT asymptotes. Yes.
- **Which decision changes if the precooler data are twice as uncertain?** None: the valve leads by a factor of two.

## Check yourself

1. State the Kelvin–Planck and Clausius forms of the second law.
2. Write the Carnot efficiency. Why can no engine between the same reservoirs exceed it?
3. Why can T ds = dh − v dp be used for an irreversible process?
4. What is S_gen for heat Q crossing from T_H to T_C? How does it behave as T_H → T_C?
5. Show that throttling an ideal gas from p₁ to p₂ generates R ln(p₁/p₂) per kilogram.
6. Why does doubling the cooling-air flow of the precooler increase its entropy generation?
7. What does the tracker do with a cycle whose Clausius integral is positive, and why?
8. Air goes from state 1 to state 2 by a free expansion and, separately, by a reversible path. Which quantities are the same in both cases and which differ?
