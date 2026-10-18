# P07.08  Exergy

> Energy is never lost; work potential is. Every component destroys exergy at the rate T₀Ṡ_gen, and the fuel's exergy, the jet's power, the hot exhaust and those destructions must add up exactly.

| | |
|---|---|
| **Learning outcomes** | 1. Define the dead state, the exergy of heat Q(1 − T₀/T) and the flow exergy ψ = (h − h₀) − T₀(s − s₀) + V²/2, and derive the control-volume exergy balance and the Gouy–Stodola theorem (exergy destroyed = T₀S_gen) from the first and second laws. 2. Compute the exergy of a hydrocarbon fuel from NASA enthalpies and entropies, with the assumptions stated, and compare it with the heating value and the Gibbs energy of reaction. 3. Implement and test a station-by-station exergy audit of a turbojet core with real combustion products, in which each component's destruction is computed two independent ways and the whole audit closes to round-off. 4. Define and compute second-law efficiencies for the compressor, turbine, nozzle and the whole core, and show that they reach 1 only where the process is reversible. 5. Decide, from the audit, where improving an engine is worth the effort and where it is not. |
| **Prerequisites** | `P07.07` (w_net = q_in(1 − T₀/T_H) − T₀Σs_gen), `P07.06` (the Brayton cycle), `P07.05` (isentropic efficiency, variable-cp isentropes), `P07.04` (entropy generation), `P07.01` (NASA polynomials, enthalpy of formation, ideal-gas mixtures), `P03.02` (Newton's method). Skill check: you can write the steady-flow energy and entropy balances for a control volume with one inlet and one exit. |
| **Engineering question** | A turbojet core (typical, illustrative: sea-level static, inlet recovery 0.98, compressor pressure ratio 12 at η_c = 0.85, turbine inlet temperature 1400 K, combustor loss 4 %, η_t = 0.88, a fully expanded nozzle of efficiency 0.97, burning Jet-A) runs at 20 kg/s of air: of the work its fuel could do, how much becomes jet power, how much is thrown away in the hot exhaust, and how much does each component destroy? The answers: of 982.2 kJ per kg of air (45.13 MJ per kg of fuel), 36.7 % becomes jet kinetic energy (7.21 MW), 26.0 % leaves in the hot exhaust, 32.9 % is destroyed in the combustor, 2.6 % in the compressor, 1.3 % in the turbine, 0.4 % in the nozzle and 0.2 % in the inlet; the audit closes to 4 × 10⁻¹⁶ of the input. Even with reversible turbomachinery the second-law efficiency would only rise from 0.367 to 0.456. |
| **Assumptions and limits** | Steady flow; adiabatic components; ideal-gas mixtures from the NASA polynomials; complete combustion of Jet-A vapour (C₁₂H₂₃) to CO₂ and H₂O (no dissociation, CO, unburnt fuel or NOₓ); fuel supplied as vapour at T₀ and the combustor pressure (the pump work for liquid fuel and its heat of vaporisation are neglected); turbine power exactly equal to compressor power (no mechanical losses, no cooling air, no bleeds); a fully expanded nozzle. The **dead state** is the ambient air at (T₀, p₀) and the products are valued against themselves at (T₀, p₀), with no condensation of water and no further diffusion into the atmosphere (a "restricted" dead state). The audit is a first model: real cores lose exergy to turbine cooling, leakage, mechanical friction and incomplete mixing that it does not represent. |
| **Analytical work** | The work available from heat and from a flowing stream; flow exergy; the control-volume exergy balance; Gouy–Stodola; the reaction's −ΔG° and the fuel exergy under the stated dead state; the turbojet's station equations (inlet, compressor, combustor energy balance linear in the fuel flow, turbine–compressor power balance, nozzle); the audit identity; second-law efficiencies and the proof that a compressor's exceeds its isentropic efficiency. |
| **Computational deliverable** | `starter/main.py`: `mixture_entropy`, `lower_heating_value`, `standard_reaction_exergy`, `products`, `engine_states`, `first_law_residuals`, `flow_exergy`, `fuel_exergy`, `exergy_audit` (the NASA data for N₂, O₂, Ar, CO₂, H₂O and Jet-A(g), species and mixture enthalpy, cp and s°, the temperature inversions, the fuel's h and s, `element_counts`, `EngineSpec`, the printing helper and `main` are given). Tests (17) prove: combustion conserves C, H, O, N and Ar atoms; the heating value follows from the database's formation enthalpies; the mixing entropy is positive and exact; flow exergy vanishes at the dead state, is positive for any temperature other than T₀ at p₀, and equals RT₀ln(p/p₀) at T₀; the first law closes for every component at three operating points; the fuel-air ratio rises with turbine inlet temperature; T₀S_gen equals the exergy-balance destruction for every component and the audit closes to 10⁻¹² of the input at three operating points; every destruction is non-negative for twelve seeded random engines; second-law efficiencies are below 1 and equal 1 in the reversible limit; the fuel exergy lies between the LHV and 1.06 × LHV and within 3 % of −ΔG°; each loss raises its own destruction and lowers η_II; and the combustor dominates and shrinks as the turbine inlet temperature rises. |
| **Practical activity** | **computational experiment**: the exergy audit of the core at sea-level static and at cruise; the reversible-turbomachinery limit; a sweep of turbine inlet temperature and compressor pressure ratio for second-law efficiency and combustor destruction; and a seeded random test of the second law. |
| **Equipment** | None. A computer with the course's Python environment. |
| **Safety** | None beyond normal computer use. All results are computed; nothing is burned. |
| **Evidence** | Completed `main.py`, passing tests, `exergy_audit.png`, the audit tables at sea level, at cruise and in the reversible limit, the sweep table, the derivations of ψ, the exergy balance and the audit identity, and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-07/grading-rubric.md`; competencies 1–4. Pass: the exergy balance and Gouy–Stodola derived from the two laws; all public tests pass; the compressor's destruction reproduced by hand from its inlet and exit states; the note gives the audit, the fuel-exergy assumption and its effect, and says which improvement would gain most. |
| **References** | `@mit-16-unified-thermodynamics-propulsion-notes` Ch. 6 "Applications of the Second Law": 6.5 Irreversibility, Entropy Changes, and "Lost Work", https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node48.html; 6.6 Entropy and Unavailable Energy (Lost Work by Another Name), .../node49.html; 6.7 Examples of Lost Work in Engineering Processes (subsections "Lost work in Adiabatic Throttling: Entropy and Stagnation Pressure Changes" and "Propulsive Power and Entropy Flux"), .../node50.html (required-free; read online, link only: lost work = T₀ΔS_total). `@mit-ocw-16-050-thermal-energy` Part 1C: 1.C.5 Irreversibility, Entropy Changes, and Lost Work; 1.C.6 Entropy and Unavailable Energy (Lost Work by Another Name); 1.C.7 Examples of Lost Work in Engineering Processes, https://ocw.mit.edu/courses/16-050-thermal-energy-fall-2002/resources/05_part1c/ (supplementary-free; the openly licensed earlier version of the same notes). `@nasa-tp-2002-211556-glenn-coefficients` "Empirical Equations for Fitting Thermodynamic Functions", https://ntrs.nasa.gov/citations/20020085330 (required-free; h(T) and s°(T) of every species, including Jet-A(g), at the dead state and each station). |
| **Extensions** | ★ Replace complete combustion by equilibrium products from `P07.11` at 1400 K and find how much dissociation changes the combustor's destruction. ★ Add a bypass stream (a turbofan): split the turbine work between the core compressor and a fan of pressure ratio 1.6, and show where the exhaust's exergy goes when part of it drives a cold jet. ★ Add a propulsive-efficiency stage to the audit at cruise: separate the jet's kinetic energy into thrust power F V₀ and the kinetic energy left in the wake, and compare the wake loss with the thermal exhaust loss. |
| **Time** | unmeasured |

## 1. The engineering problem

`P07.07` ended with a striking identity: a cycle's net work is the Carnot work of the heat it takes in, minus T₀ times all the entropy its components generate. Each component's T₀s_gen is a quantity of work potential destroyed, and the destructions add up. That makes the second law an accounting tool. An engineer who knows how much work potential enters a machine, how much leaves as useful output, and how much each component destroys knows where the losses are, how large each is, and which are worth attacking.

The first law cannot do this. By energy alone, a turbojet "loses" its energy almost entirely in the hot exhaust, since the combustor is adiabatic and loses nothing. That is misleading twice over: a kilojoule of hot exhaust at 820 K is worth much less than a kilojoule of shaft work, and the combustor, which loses no energy, is where a third of the fuel's work potential disappears. The quantity that tracks work potential is **exergy** (also called availability): the most work a system could deliver in reaching equilibrium with its surroundings, the **dead state**.

This lesson builds the exergy audit of a turbojet core, using the illustrative engine of `P07.06` and `P07.07` (π = 12, 1400 K turbine inlet temperature, typical component efficiencies) but with the air-standard idealisation removed: real air, real Jet-A, real combustion products, the turbine driving the compressor and a nozzle making the jet. All results are computed; there is no fuel, flame or combustion of any kind outside the computer. Success criteria, written first: atoms conserved; the first law closing for every component to round-off; each component's destruction computed from its entropy generation and, independently, from its exergy balance, agreeing to round-off; the whole audit (fuel exergy = jet power + exhaust exergy + Σ destructions) closing to 10⁻¹²; every destruction non-negative for any reasonable engine; second-law efficiencies ≤ 1, equal to 1 exactly where the process is reversible.

## 2. Understand and derive

### 2.1 The dead state and the exergy of heat

Exergy is measured against the surroundings: here, still air at T₀ and p₀ (288.15 K and 101 325 Pa at sea level; 216.65 K and 22 632 Pa at 11 km). A system at the dead state can do no work, because nothing is left to exploit. Heat Q available at a temperature T can at best drive a Carnot engine rejecting to T₀ (`P07.04`), so its exergy is

**Ex_Q = Q(1 − T₀/T).**   (1)

Heat at T₀ is worth nothing; heat at 1400 K, with T₀ = 288.15 K, is worth 79 % of itself.

### 2.2 Flow exergy

A stream of fluid at (h, s, V) enters a device that may exchange heat only with the surroundings at T₀ and leaves at the dead state, at rest. Per unit mass, the first law gives w = (h + V²/2) − h₀ − q₀, where q₀ is the heat given to the surroundings; the second law gives s_gen = s₀ − s + q₀/T₀ ≥ 0. Eliminating q₀:

w = (h − h₀) − T₀(s − s₀) + V²/2 − T₀s_gen.

The work is largest when s_gen = 0. That maximum is the **flow exergy**:

**ψ = (h − h₀) − T₀(s − s₀) + V²/2.**   (2)

Three properties follow. ψ = 0 at the dead state. At p₀, any temperature other than T₀ gives ψ > 0: hot or cold, a temperature difference from the surroundings can drive an engine. At T₀, ψ = RT₀ln(p/p₀), which is **negative below p₀**: a stream below ambient pressure needs work to be pushed out into the atmosphere. (The tests check all three.) For air at 750 K and 1 MPa, ψ = 384.37 kJ/kg: 195.01 kJ/kg from its temperature and 189.36 kJ/kg from its pressure.

### 2.3 The exergy balance and the Gouy–Stodola theorem

For a steady control volume with inlets and exits, heat Q̇_k crossing its boundary where the temperature is T_k, and shaft power Ẇ out, write the energy balance and the entropy balance:

0 = Σ_in ṁ(h + V²/2) − Σ_out ṁ(h + V²/2) + Σ_k Q̇_k − Ẇ,
0 = Σ_in ṁs − Σ_out ṁs + Σ_k Q̇_k/T_k + Ṡ_gen.

Multiply the second by T₀ and subtract it from the first; mass conservation lets h₀ and s₀ be added and subtracted freely. The result is the **exergy balance**:

**Ẇ = Σ_in ṁψ − Σ_out ṁψ + Σ_k Q̇_k(1 − T₀/T_k) − T₀Ṡ_gen.**   (3)

The work delivered is the exergy brought in by flows and heat, minus the exergy carried out by flows, minus **T₀Ṡ_gen**. Since Ṡ_gen ≥ 0, this last term is always a loss: it is the **exergy destroyed**, and the statement that destruction equals T₀Ṡ_gen is the **Gouy–Stodola theorem**. `P07.07`'s w_net = q_in(1 − T₀/T_H) − T₀Σs_gen is (3) applied to a whole cycle.

The code uses (3) as an independent check. For each component, it computes the destruction twice: as T₀ × (entropy out − entropy in), and as (exergy in − exergy out − work out) from (2). The two use the same states but different combinations of h and s; they agree only if the component's energy balance holds.

### 2.4 The exergy of a fuel

The fuel's work potential is chemical. Burning it with air and bringing the products to the dead state releases the heating value as heat, but the reversible work of the same change is set by the Gibbs function, not the enthalpy. For one kilogram of air and n_f moles of fuel, with B = H − T₀S of each stream:

**Ex_fuel = B_air(T₀, p₀) + n_f b_fuel(T₀, p₃) − B_prod(T₀, p₀).**   (4)

It is the most work any device could extract from the air and fuel entering and the products leaving, all at the dead state. Assumptions, stated: the fuel is Jet-A **vapour** (the NASA species Jet-A(g), C₁₂H₂₃, M = 167.31 g/mol, h_f° = −249 657 J/mol), supplied at T₀ and the combustor pressure p₃ (the pump work for liquid fuel and the heat of vaporisation are neglected); the products are an ideal-gas mixture at (T₀, p₀) with **no condensation** of their water; and they are valued against themselves, not after diffusing into the atmosphere (a **restricted dead state**). Mixtures use s = Σn_i[s_i°(T) − R ln(y_i p/p°)], with the mixing term.

Three numbers for Jet-A from the NASA data:

- the **lower heating value** at 298.15 K (water as vapour), from the formation enthalpies: **43.353 MJ/kg**;
- **−ΔG°** of C₁₂H₂₃ + 17.75 O₂ → 12 CO₂ + 11.5 H₂O at 298.15 K, every species pure at 1 bar: **44.185 MJ/kg**, 1.019 × LHV, larger because the reaction turns 18.75 moles of gas into 23.5 and its entropy change is positive;
- (4) for the sea-level engine: 982.20 kJ per kg of air, **45.126 MJ per kg of fuel, 1.041 × LHV**: larger still, because the CO₂ and H₂O leave diluted in excess air, where their entropy is higher (their partial pressures are low), so B_prod is lower.

Different dead-state choices move the fuel exergy by a few percent, and with it every percentage in the audit (Problem 5). What does not move is the closure: the destructions and the products are computed from the same states whatever the fuel exergy is called.

### 2.5 The engine

![Turbojet core for the exergy audit: ambient air at station 0 (T₀, p₀, V₀) passes the inlet to 2, the compressor to 3, the combustor (where fuel is added) to 4, the turbine to 5 and the nozzle to the jet at 9 (T₉, p₀, V₉). A shaft connects the turbine to the compressor.](figures/turbojet-stations.svg)
*Figure 1. Stations of the turbojet core used in the exergy audit. The shaft returns the turbine's work to the compressor, so the product that leaves the core is the jet's kinetic-energy gain.*

Per kilogram of air (NASA mixtures throughout; the product mixture follows from the fuel flow):

- **Inlet** (0 → 2): adiabatic, h₂ = h₀ + V₀²/2; the stagnation pressure is the isentropic ram value times the recovery π_d.
- **Compressor** (2 → 3): p₃ = π_c p₂, isentropic efficiency η_c, work w_c = H₃ − H₂.
- **Combustor** (3 → 4): p₄ = p₃(1 − Δ_cc), T₄ = turbine inlet temperature. The products are N_air + n_f ν (with ν the stoichiometric changes: −17.75 O₂, +12 CO₂, +11.5 H₂O per mole of fuel), so the energy balance H_air(T₃) + n_f h_fuel(T₀) = Σ(N_air,i + n_f ν_i)h_i(T₄) is **linear in n_f** and is solved exactly, with no iteration. At sea level f = 0.02177 kg of fuel per kg of air, about a third of stoichiometric (0.06818).
- **Turbine** (4 → 5): its power equals the compressor's, H₄ − H₅ = w_c; the exit pressure comes from the isentropic efficiency η_t.
- **Nozzle** (5 → 9): expands to p₀ with isentropic efficiency η_n; the jet velocity V₉ = √(2(H₅ − H₉)/m_p), with m_p = 1 + f the mass of products.

The **product** of a turbojet core is the kinetic-energy gain of the flow, m_pV₉²/2 − V₀²/2: what the propulsion system then turns into thrust work (`P08` and beyond). The **exhaust exergy** ψ₉ (thermal and pressure part, relative to the products at the dead state) leaves with the jet and is destroyed as the jet mixes with the atmosphere: a loss, but outside the engine.

### 2.6 The audit identity and second-law efficiencies

Add the exergy balances (3) of the five components. The flow exergies of the internal stations cancel, the compressor's work input cancels the turbine's output, and what remains is

**Ex_fuel + V₀²/2 = (m_pV₉²/2 − V₀²/2) + V₀²/2 + ψ₉ + Σ T₀s_gen,i**,   (5)

that is, exergy in (fuel, plus the ram kinetic energy of the incoming air in the engine frame) = product + exhaust exergy + destructions. Because each destruction is computed from entropy and the fuel exergy from (4), (5) closes only if the first law closes in every component. The **second-law efficiency** of the core is the useful fraction,

**η_II = (m_pV₉²/2 − V₀²/2)/(Ex_fuel + V₀²/2),**   (6)

and of a component, the ratio of exergy delivered to exergy spent: compressor (ψ₃ − ψ₂)/w_c, turbine w_c/(ψ₄ − ψ₅), nozzle (kinetic energy out)/(ψ₅ − ψ₉). Each is 1 − (destruction)/(exergy spent), so each is ≤ 1 with equality exactly when s_gen = 0.

**A compressor's second-law efficiency exceeds its isentropic efficiency.** η_II = 1 − T₀s_gen/w_c and η_s = 1 − (h₃ − h₃ₛ)/w_c. Along the exit isobar from 3s to 3, h₃ − h₃ₛ = ∫T ds ≥ T₀(s₃ − s₃ₛ) = T₀s_gen, because T ≥ T₀ there. So η_II ≥ η_s: the extra work a real compressor needs is not all lost; part of it is still in the hotter air as exergy. At sea level, η_s = 0.85 and η_II = 0.9286.

### 2.7 Where the model stops

- **Combustion chemistry.** Complete combustion is assumed; at a lean 1400 K the equilibrium solver of `P07.11` can check how small dissociation is (Extension 1). Unburnt fuel and CO are not modelled.
- **Cooling air, leakage, bleeds and mechanical losses** destroy exergy in real cores and are absent.
- **The dead state.** The products' water would condense at 288 K, and the CO₂ and H₂O would eventually diffuse into the atmosphere; including either changes the fuel exergy by a few percent.
- **The product.** At cruise, the jet's kinetic energy is not all useful: only the thrust power counts, and the wake's kinetic energy is another loss (Extension 3).

## 3. Build the computational model

Open `starter/main.py`. The NASA coefficients for N₂, O₂, Ar, CO₂, H₂O and Jet-A(g), `species_props`, `mixture_enthalpy`, `mixture_cp`, `mixture_s0`, the two Newton inversions for mixtures, `fuel_enthalpy`, `fuel_entropy`, `element_counts`, `EngineSpec`, the printing helper and `main` are given. Amounts are moles per kilogram of air, so every energy is per kilogram of air.

### Step 1: mixtures and the fuel

`mixture_entropy(n, T, p)` adds the pressure and mixing term −RΣn_i ln(y_i p/p°) to `mixture_s0`, skipping absent species. `products(n_fuel)` is N_AIR + n_f ν (raise `ValueError` if O₂ would go negative). `lower_heating_value` and `standard_reaction_exergy` are the reaction's −ΔH and −ΔG at 298.15 K per kilogram of fuel.

### Step 2: the engine

`engine_states(spec)` follows §2.5 station by station. The combustor is one line once the enthalpies at T₄ are known:

```python
h4 = np.array([species_props(k, s.T_tit)[1] for k in SPECIES])
n_fuel = (na @ h4 - H3) / (fuel_enthalpy(s.T0) - NU_FUEL @ h4)
```

`first_law_residuals` returns each component's energy-in minus energy-out.

### Step 3: exergy

`flow_exergy` is (2) for a mixture of fixed composition; `fuel_exergy` is (4); `exergy_audit` computes s_gen and T₀s_gen for each component, the same destructions from exergy balances, the fuel and ram exergy, the product, the exhaust exergy, the residual of (5), η_II, the thermal efficiency (product over fuel heating value) and the component second-law efficiencies.

Running `main()` prints

```
Jet-A(g) C12H23 (NASA): LHV(298.15 K) = 43.353 MJ/kg, -dG°(298.15 K) = 44.185 MJ/kg, ratio 1.0192
turbojet core, sea-level static (typical, illustrative)
  f = 0.02177 (fuel-air ratio), w_c = 351.10 kJ/kg, V9 = 840.0 m/s, specific thrust (m_p V9 - V0) = 858.2 N s/kg
  fuel exergy 982.20 kJ/kg air = 45.126 MJ/kg fuel; ram exergy 0.00 kJ/kg
    product: jet KE gain       360.45   36.7 %
    exhaust (lost)             255.83   26.0 %
    destroyed: inlet             1.67    0.2 %   (balance     1.67, S_gen 5.80 J/(kg K))
    destroyed: compressor       25.07    2.6 %   (balance    25.07, S_gen 87.01 J/(kg K))
    destroyed: combustor       322.70   32.9 %   (balance   322.70, S_gen 1119.89 J/(kg K))
    destroyed: turbine          12.53    1.3 %   (balance    12.53, S_gen 43.50 J/(kg K))
    destroyed: nozzle            3.94    0.4 %   (balance     3.94, S_gen 13.68 J/(kg K))
  Gouy-Stodola closure: in - product - exhaust - sum D = -3.49e-10 J/kg (relative -3.6e-16); max |T0 S_gen - balance| = 5.1e-10 J/kg
  eta_II = 0.3670, eta_th = 0.3820; component eta_II: compressor 0.9286, turbine 0.9655, nozzle 0.9892
  at 20.0 kg/s of air: fuel 0.4353 kg/s, jet power 7.209 MW, combustor destruction 6.454 MW
```

(abridged), then the station tables and the reversible and cruise audits, the sweep, and writes `exergy_audit.png` (the audit as stacked bars; η_II and the combustor's share against pressure ratio).

### Verify

`python -m pytest curriculum/phase-07/lesson-08-exergy` runs 17 tests:

1. **Atoms**: C, H, O, N and Ar balance to 10⁻¹⁴ for three fuel amounts; too much fuel raises an error.
2. **Heating value** from the database's formation enthalpies to 10⁻⁵ (the database used R = 8.314510 J/(mol·K), we use the CODATA 2018 value, so each species' h(298.15 K) reproduces its Hf298 to about 2 J/mol); LHV < −ΔG° < 1.05 LHV.
3. **Mixing entropy** positive and equal to −RΣn_i ln y_i.
4. **Flow exergy**: zero at the dead state; positive at p₀ for T ≠ T₀; RT₀ln(p/p₀) at T₀; kinetic energy adds one for one.
5. **First law** (static, cruise, reversible): every component below 10⁻⁸ J/kg; lean of stoichiometric; p₅ > p₀ and V₉ > V₀.
6. **Fuel-air ratio** rises with turbine inlet temperature.
7. **Gouy–Stodola** (static, cruise, reversible): (5) closes to 10⁻¹² of the input; T₀s_gen equals the exergy-balance destruction for every component.
8. **Second law, random engines**: twelve seeded engines (rng seed 7) with T₀, V₀, π_d, π_c, η_c, T₄, Δ_cc, η_t and η_n drawn over plausible ranges; every destruction ≥ 0 and the audit closes.
9. **Second-law efficiencies**: 0 < η_II < η_th < 1 and component values in (0, 1); all exactly 1 in the reversible limit, where only the combustor and the exhaust lose exergy.
10. **Fuel exergy** between LHV and 1.06 LHV, within 3 % of −ΔG°.
11. **Each loss** raises its own component's destruction and lowers η_II; **the combustor** dominates and its share falls as T₄ rises. **main** writes the figure.

### Compare

Three routes, one answer. The destruction of each component by T₀s_gen and by its exergy balance agree to 5 × 10⁻¹⁰ J/kg. The audit identity (5) closes to 3.5 × 10⁻¹⁰ J/kg out of 982 200 J/kg. And the fuel exergy (4), computed from mixture B functions, lies within 2.1 % of −ΔG°, computed from pure species, with the difference explained by the mixing and dilution terms of §2.4.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All results are computed; there is no combustion, flame or fuel of any kind in this lesson. Write each prediction before running it.

1. **The sea-level audit.** *Predict* the largest destruction and roughly what share of the fuel exergy reaches the jet. *Finding:* the table in part 3. The combustor destroys 32.9 %, more than the inlet, compressor, turbine and nozzle together (43.21 kJ/kg, 4.4 %) by a factor of 7.5; the exhaust carries away 26.0 %; 36.7 % becomes jet kinetic energy (η_II = 0.367, thermal efficiency 0.382).

2. **The reversible limit.** Set π_d = 1, η_c = η_t = η_n = 1 and no combustor pressure loss. *Predict* η_II. *Finding:* the four turbomachine destructions are zero and their second-law efficiencies exactly 1; η_II rises to 0.456 (product 474.14 kJ/kg, 45.6 %); the combustor still destroys 347.36 kJ/kg (33.4 %, slightly more than before because the colder compressor delivery air, 580.54 instead of 630.58 K, needs more fuel: f = 0.02307) and the exhaust carries 219.22 kJ/kg (21.1 %). Perfect turbomachinery buys 9 points; the remaining 55 % is lost to combustion and the exhaust, which no component efficiency can recover.

3. **Cruise.** 11 km ISA (216.65 K, 22 632 Pa), Mach 0.8 (236.1 m/s), the same engine. *Predict* whether η_II rises or falls. *Finding:* ram compression raises the compressor inlet to 244.43 K and 33.81 kPa; the fuel exergy is 1083.99 kJ/kg plus 27.86 kJ/kg of ram exergy; the product is 484.39 kJ/kg (43.6 %), the exhaust 280.76 (25.3 %), the combustor 286.43 (25.8 %), compressor 18.78 (1.7 %), turbine 7.68 (0.7 %), nozzle 4.69 (0.4 %), inlet 1.26 (0.1 %). η_II = 0.436, η_th = 0.461. The cold dead state makes the 1400 K combustion less irreversible relative to its surroundings, and the ram compression adds pressure ratio for free.

4. **Turbine inlet temperature and pressure ratio.** Sea-level static, T₄ = 1200, 1400, 1600 K and π_c = 6, 12, 20, 30. *Predict* how the combustor's share moves. *Finding:*

   | T₄ (K) | π_c | η_II | η_th | combustor D / exergy in | exhaust / exergy in |
   |---|---|---|---|---|---|
   | 1200 | 6 | 0.2813 | 0.2931 | 0.3840 | 0.2955 |
   | 1200 | 12 | 0.3457 | 0.3609 | 0.3578 | 0.2321 |
   | 1200 | 30 | 0.3497 | 0.3662 | 0.3285 | 0.1834 |
   | 1400 | 6 | 0.2901 | 0.3016 | 0.3517 | 0.3297 |
   | 1400 | 12 | 0.3670 | 0.3820 | 0.3285 | 0.2605 |
   | 1400 | 30 | 0.4225 | 0.4408 | 0.3017 | 0.1957 |
   | 1600 | 6 | 0.2934 | 0.3044 | 0.3254 | 0.3592 |
   | 1600 | 12 | 0.3760 | 0.3905 | 0.3048 | 0.2864 |
   | 1600 | 30 | 0.4512 | 0.4695 | 0.2806 | 0.2135 |

   The combustor's share falls with both: higher π_c delivers hotter air to the combustor, and higher T₄ makes the products hotter, so the chemical energy is released into gas nearer its own "temperature" and less work potential is destroyed. The exhaust share falls with π_c (a cooler exhaust) and rises with T₄ (a hotter one). At 1200 K the second-law efficiency peaks between π = 12 and 30 (0.3652 at 20); at 1600 K it is still rising at 30.

5. **The second law, at random.** Twelve engines drawn with `numpy.random.default_rng(7)` over the ranges of test 8. *Predict* whether any destruction can come out negative. *Finding:* none does, and every audit closes to 10⁻¹² of its input (test 8).

Raw data format: CSV tables `item, kJ_per_kg_air, percent` for each audit (1–3); `T4_K, pi_c, eta_II, eta_th, D_comb_share, exhaust_share` (4); `seed, case, component, D_Jkg` (5); plus `exergy_audit.png`.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| atoms in = atoms out | exact | 10⁻¹⁴ | round-off | — |
| LHV from Hf298 | 43.353 MJ/kg | 43.353 MJ/kg | < 10⁻⁵ | R values differ (2 J/mol per species) |
| component first law | 0 | < 10⁻⁸ J/kg | round-off | Newton to 10⁻¹⁰ K |
| T₀s_gen against exergy balance | equal (Gouy–Stodola) | max difference 5.1 × 10⁻¹⁰ J/kg | round-off | — |
| audit (5) | 0 | −3.49 × 10⁻¹⁰ J/kg of 982 200 | 3.6 × 10⁻¹⁶ relative | round-off |
| reversible turbomachinery | η_II,comp = 1, D = 0 | 1.0000, 0.00 | < 10⁻¹² | — |
| fuel exergy (4) against −ΔG° | within a few % | 45.126 against 44.185 MJ/kg | 2.1 % | dilution and mixing of the products |
| η_II of the compressor against η_s | η_II ≥ η_s | 0.9286 against 0.85 | — | §2.6 |

The engineering answer: **of the 45.1 MJ of work potential in each kilogram of Jet-A, the sea-level core turns 36.7 % into jet power (7.21 MW at 20 kg/s), throws 26.0 % away in the hot exhaust and destroys 32.9 % in the combustor; the compressor, turbine, nozzle and inlet together destroy 4.4 %. Perfect turbomachinery would lift the second-law efficiency only from 0.367 to 0.456.** The two big items are attacked by cycle design, not component polish: a higher pressure ratio and turbine inlet temperature cut the combustor's share (to 28 % at π = 30, 1600 K), and a bypass stream or a recuperator recovers part of the exhaust's exergy.

What to improve: the dead-state convention, which moves every percentage by up to two points (Problem 5), and the missing losses (cooling air, mechanical), which would add perhaps a few percent of destruction in the turbine. What not to improve: the property model or the solver; the audit already closes to round-off.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2), (3), (4), (5) and of η_II ≥ η_s for a compressor; `exergy_audit.png`; the experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6, 7) answering the engineering question, stating the fuel-exergy assumption and recommending the one change that would gain most.

### Adversarial self-checks

- **What if every component is ideal?** The turbomachine destructions vanish and their second-law efficiencies are exactly 1 (test 9), but η_II stays at 0.456: combustion and the exhaust are not component inefficiencies.
- **What if T₀ changes?** Every destruction scales with T₀ at fixed s_gen, and the fuel exergy hardly changes; the cruise audit shows the combustor's share falling from 32.9 % to 25.8 % with the colder dead state.
- **Which assumption fails first in the real world?** The absence of turbine cooling air: real cores at 1400 K cool their turbines with several percent of the compressor flow, and the mixing of cooling air with hot gas destroys exergy the model does not see.
- **What measurement could falsify the model?** Fuel flow, station temperatures and pressures and thrust on a test bed: the measured fuel flow must match the combustor energy balance, and the measured thrust the jet velocity, within their uncertainties. A measured fuel flow higher than predicted points to combustion inefficiency or cooling flows.
- **Does the result survive a change of method?** Two routes to every destruction agree to 5 × 10⁻¹⁰ J/kg; the audit closes to 10⁻¹⁶; twelve random engines all obey the second law.
- **Which decision changes if the fuel exergy is defined differently?** Taking the fuel exergy as the LHV instead of (4) raises η_II from 0.367 to 0.382 (Problem 5); the ranking of the losses, and therefore what to work on, does not change.

## Check yourself

1. What is the dead state, and why does exergy need one while energy does not?
2. Write the flow exergy of a stream. Can it be negative? Give an example.
3. State the Gouy–Stodola theorem and the two balance laws it combines.
4. Why is the combustor the largest exergy destroyer in a gas turbine, even though it loses no energy?
5. A compressor has η_s = 0.85. Is its second-law efficiency higher or lower than 0.85, and why?
6. Why is the fuel exergy of Jet-A slightly larger than its lower heating value?
7. The hot exhaust carries 26 % of the fuel exergy out of the core. Is that exergy "destroyed" by the engine? Where is it destroyed, and how could some of it be recovered?
8. In the reversible-turbomachinery limit, which two items still prevent η_II = 1?
