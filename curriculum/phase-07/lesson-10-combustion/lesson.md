# P07.10  Chemical Reactions & Combustion

> Count the atoms, then count the energy from the elements up. A flame temperature is nothing more than the temperature at which the products hold all the enthalpy the reactants brought in.

| | |
|---|---|
| **Learning outcomes** | 1. Balance the complete combustion of any CxHyOz fuel with air or oxygen at any equivalence ratio, and convert between equivalence ratio, fuel–air ratio and air–fuel ratio. 2. Compute the lower and higher heating values of a fuel from enthalpies of formation, and check them against standard enthalpies of combustion and against h_fg from the steam table. 3. Write the first law for a steady, adiabatic combustor with enthalpies measured from the elements, and solve it for the adiabatic flame temperature with NASA polynomials and a bracketed root finder. 4. Implement and test a combustion module for Jet-A/air at combustor-inlet conditions, CH₄/O₂ and H₂/O₂, and find the fuel–air ratio that gives a turbine inlet temperature. 5. Judge where complete combustion stops being a usable model. |
| **Prerequisites** | `P07.01` (NASA polynomials, h including the enthalpy of formation, dry air), `P07.02` (steady-flow energy balance), `P07.03` (h_fg from the steam table), `P07.04` (entropy generation), `P03.02` (root finding with `brentq`), Bridge B2 or equivalent general chemistry (moles, balancing equations). Skill check: you can balance CH₄ + O₂ → CO₂ + H₂O by hand and find the root of a monotone function on a bracket. |
| **Engineering question** | A gas-turbine combustor receives compressor air at 800 K and must deliver 1700 K to the turbine (typical, illustrative values), burning Jet-A. What fuel–air ratio does that take, and how hot would a stoichiometric zone of the same combustor get? The answer, with complete combustion and the fuel as vapour at 298.15 K: φ = 0.396, f = 0.0270 (1 kg of fuel per 37.0 kg of air), and 2755 K at stoichiometric, the peak of the T_ad(φ) curve. |
| **Assumptions and limits** | Steady flow, adiabatic, no work, kinetic energy neglected; ideal-gas mixtures (so h, and the complete-combustion flame temperature, do not depend on pressure); **complete combustion**: lean and stoichiometric products are CO₂, H₂O, spare O₂ and the inerts; rich products follow a fixed order (each C takes one O to make CO, the H then makes H₂O, any O left turns CO into CO₂); no dissociation, no NO, no soot. Jet-A is the NASA CEA surrogate C₁₂H₂₃ as a vapour (its heat of vaporization is not in this course's data). The dry-air model of Phase 07. Fails: above about 2000 K, where dissociation absorbs energy (lesson `P07.11`: 154 K lower at stoichiometric for Jet-A/air at 30 bar, more than 1800 K lower for the oxygen flames); for very rich mixtures (soot, unburned fuel). **Computational only: there is no at-home combustion, flame, burning, fuel handling or fuel of any kind in this lesson.** |
| **Analytical work** | Atom balances and the stoichiometric oxygen x + y/4 − z/2; φ, f and AFR; enthalpy of formation and the "absolute" enthalpy; Hess's law for LHV and HHV and HHV − LHV = (water per kg fuel) × h_fg; the adiabatic-combustor first law; a hand estimate of T_ad with a mean cp; why T_ad peaks at stoichiometric; the entropy generated in a combustor. |
| **Computational deliverable** | `starter/main.py`: `o2_stoich`, `reactants`, `fuel_oxidizer_ratio`, `equivalence_ratio`, `complete_products`, `mixture_enthalpy`, `reactant_enthalpy`, `heating_values`, `adiabatic_flame_temperature`, `phi_for_temperature` (the NASA functions, `molar_mass`, `mixture_mass`, `element_totals` and `main` are given). Tests (24) prove atoms and mass balance exactly for five fuel–oxidizer pairs at five equivalence ratios, the lean/stoichiometric/rich product regimes, f ↔ φ inversion, Hess's law, HHV against the standard enthalpies of combustion of methane and hydrogen, HHV − LHV against the steam-table h_fg, the energy balance at T_ad to round-off, the peak at φ ≈ 1, the φ → 0 limit, the inlet-temperature trend, a two-step first-law identity, and the turbine-inlet solve. |
| **Practical activity** | **computational experiment** (computational only: no combustion, flame, burning, fuel handling or fuel of any kind): sweeps of T_ad against φ for Jet-A/air at 298 K and at an 800 K combustor inlet, CH₄/air, H₂/air, CH₄/O₂ and H₂/O₂; a turbine-inlet fuel–air ratio solve; an inlet-temperature sensitivity study. |
| **Equipment** | None (a computer with the course's Python environment). There is no fuel of any kind in this lesson; nothing is burned. |
| **Safety** | None beyond normal computer use. This lesson is computational only: no at-home combustion, flame, burning, fuel handling or fuel of any kind. |
| **Evidence** | `flame_temperature.png`; the heating-value table with both checks; the T_ad(φ) table; the turbine-inlet solve with its energy-balance residual; a half-page note answering the engineering question and saying which of its numbers lesson `P07.11` will change and by roughly how much. |
| **Assessment** | `assessments/phase-07/grading-rubric.md`; competencies 1–4. Pass: the Jet-A stoichiometry and the CH₄ LHV/HHV done by hand; all public tests pass; the CH₄/air T_ad estimated by hand with a mean cp to within 3 %; the note explains why the stoichiometric value is not credible and the turbine-inlet one nearly is. |
| **References** | `@mit-16-unified-thermodynamics-propulsion-notes` Ch. 15 "Generating Heat: Thermochemistry": 15.1 Fuels; 15.2 Fuel-Air Ratio; 15.3 Enthalpy of formation; 15.4 First Law Analysis of Reacting Systems; 15.5 Adiabatic Flame Temperature (required-free, read online; https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node108.html). `@mit-ocw-16-050-thermal-energy` Part 2C "Introduction to Thermochemistry": 2.C.1–2.C.5 (supplementary-free; https://ocw.mit.edu/courses/16-050-thermal-energy-fall-2002/resources/09_part2c/). `@openstax-chemistry-2e` 4.1 Writing and Balancing Chemical Equations; 4.4 Reaction Yields (Limiting Reactant); 5.3 Enthalpy (Standard Enthalpy of Combustion, Table 5.2; Standard Enthalpy of Formation; Hess's Law) (required-free; https://openstax.org/books/chemistry-2e/pages/5-3-enthalpy). `@nasa-tp-2002-211556-glenn-coefficients` "Assigned Enthalpy Values"; Table B1 (required-free; https://ntrs.nasa.gov/citations/20020085330). `@nasa-rp-1311-cea-analysis` Part I Ch. 1 Introduction (the hp problem); Ch. 9 "Input Calculations" (equivalence ratios) (supplementary-free; https://ntrs.nasa.gov/citations/19950013764). `@nasa-cea-github` `data/thermo.inp` (supplementary-free; https://github.com/nasa/cea). |
| **Extensions** | ★ Give the fuel a liquid enthalpy: subtract an assumed heat of vaporization (label it illustrative) from the Jet-A(g) enthalpy and find how much it changes the turbine-inlet φ. ★ Add a combustor efficiency η_b (the fraction of the LHV released) and find η_b that matches a measured exit temperature. ★ Add propane from `nasa9_p07` (its fit starts at 300 K) and compare its stoichiometric T_ad in air with methane's. |
| **Time** | unmeasured |

## 1. The engineering problem

A turbofan's combustor is a short annulus where fuel is sprayed into air that the compressor has squeezed to 20–40 bar and heated to 700–900 K. Its job is to raise the gas to the turbine inlet temperature, about 1700 K in a modern engine at take-off (typical, illustrative), and no hotter: the first turbine blades already run above their alloy's melting point and survive only by internal cooling. So the combustor designer needs two numbers: how much fuel per kilogram of air gives exactly that temperature, and how hot the gas gets where fuel and air meet near the chemically correct (stoichiometric) proportion, because that is where liners burn through and nitrogen oxides form.

Both numbers come from the same two conservation laws. Atoms are neither made nor destroyed, which fixes the products. Energy is conserved, which, once every enthalpy is measured from the chemical elements, says the products leave with exactly the enthalpy the reactants brought. This lesson builds that model for Jet-A, methane and hydrogen, burning in air and, for rocket chambers, in pure oxygen. It is the simplest useful model, "complete combustion", and it is wrong in a way the next lesson fixes: hot products partly fall apart (dissociate), and that costs energy.

**This lesson is computational only. There is no at-home combustion, flame, burning, fuel handling or fuel of any kind; all results are computed.**

Success criteria, written first: atoms balance exactly; the heating values reproduce standard enthalpies of combustion to 0.1 % and the steam table's h_fg to 0.1 %; the energy balance at the flame temperature closes to round-off; T_ad peaks at φ ≈ 1; the turbine-inlet solve reproduces 1700 K.

## 2. Understand and derive

### 2.1 Atoms: stoichiometry and the equivalence ratio

Burning a fuel C_xH_yO_z completely means every carbon atom ends in CO₂ and every hydrogen atom in H₂O. Balancing oxygen atoms,

C_xH_yO_z + a O₂ → x CO₂ + (y/2) H₂O,  **a = x + y/4 − z/2**.   (1)

For methane a = 2, for hydrogen a = ½, for the Jet-A surrogate C₁₂H₂₃ a = 12 + 23/4 = 17.75. Air brings inerts with the oxygen: in the phase's dry air each mole of O₂ comes with 0.7808/0.2095 = 3.727 mol of N₂, 0.0444 mol of Ar and 0.0019 mol of CO₂. So stoichiometric methane–air is CH₄ + 2 O₂ + 7.454 N₂ + 0.0888 Ar + 0.0038 CO₂ (methane is 9.48 % of the mixture by mole).

The **fuel–air ratio** f is the mass of fuel per mass of air; its stoichiometric value f_st is 0.05802 for CH₄, 0.02916 for H₂ and 0.06818 for Jet-A (air–fuel ratios 17.24, 34.29 and 14.67). The **equivalence ratio** compares a mixture with stoichiometric,

**φ = f/f_st**,  O₂ supplied = a/φ per mole of fuel.   (2)

φ < 1 is **lean** (spare oxygen in the products), φ > 1 **rich** (not enough oxygen). Lean products are simply x CO₂, (y/2) H₂O, a(1/φ − 1) O₂ and the inerts. Rich products need a rule, because atom balances alone leave more unknowns (CO, CO₂, H₂, H₂O) than equations. This lesson uses a fixed order: each carbon first takes one oxygen atom (CO), the hydrogen then takes what it can (H₂O, the rest stays H₂), and any oxygen still left turns CO into CO₂. It is a convention, not physics; lesson `P07.11` replaces it with equilibrium. If there are fewer O atoms than C atoms, solid carbon (soot) would form and the model refuses.

### 2.2 Energy measured from the elements

Lesson `P07.01` built the NASA polynomials so that h(298.15 K) equals the **enthalpy of formation** Hf: the enthalpy change of making one mole of the species from its elements in their reference states (O₂, H₂, N₂ gases, graphite) at 298.15 K and 1 bar. Elements have Hf = 0. Water vapour has Hf = −241.826 kJ/mol, liquid water −285.830 kJ/mol, CO₂ −393.510 kJ/mol, methane −74.600 kJ/mol, the Jet-A vapour surrogate −249.657 kJ/mol. With every species measured from the same elements, chemistry needs no separate bookkeeping: the enthalpy of any mixture at any T is Σ n_j h_j(T).

The first law for a steady combustor with no work and kinetic energy changes neglected is

Σ_reactants n_i h_i(T_i) + Q = Σ_products n_j h_j(T_out).   (3)

**Heating values.** Burn one mole of fuel with its stoichiometric oxygen, all at 298.15 K, and cool the products back to 298.15 K. The heat removed per kilogram of fuel is the heating value. If the water leaves as vapour it is the **lower heating value (LHV)**; if it condenses, the **higher heating value (HHV)**:

LHV = [Hf(fuel) − x Hf(CO₂) − (y/2) Hf(H₂O, g)]/M_fuel,  HHV with Hf(H₂O, l).   (4)

For methane, by hand: −74 600 + 393 510 + 2 × 241 826 = 802 562 J/mol, or 50.03 MJ/kg; with liquid water 890 570 J/mol, 55.51 MJ/kg. The code (using the polynomials, which reproduce Hf to about 2 J/mol) gives:

| fuel | LHV (MJ/kg) | HHV (MJ/kg) | HHV − LHV | (y/2)(M_w/M_f) h_fg(25 °C) from the steam table |
|---|---|---|---|---|
| CH₄ | 50.027 | 55.513 | 5.4859 | 5.4839 (+0.037 %) |
| H₂ | 119.960 | 141.788 | 21.8286 | 21.8205 (+0.037 %) |
| Jet-A(g) | 43.353 | 46.377 | 3.0246 | 3.0235 (+0.037 %) |

Two independent checks. Per mole, the HHVs are 890.57 kJ/mol for methane and 285.83 kJ/mol for hydrogen; OpenStax Chemistry 2e, Table 5.2, lists the standard enthalpies of combustion as −890.8 and −285.8 kJ/mol. And HHV − LHV must be the latent heat of the water formed: the NASA difference Hf(g) − Hf(l) = 44.004 kJ/mol against the steam table's h_fg at 25 °C, 2441.68 kJ/kg (NIST SRD 69 via Yan), agree to 0.04 % for every fuel because each fuel's difference is the same water. A gas turbine never condenses its exhaust water, so engines are rated on the LHV.

### 2.3 The adiabatic flame temperature

With Q = 0 in (3), the products' temperature T_ad satisfies

**H_products(T_ad) = H_reactants(T_fuel, T_air)**,   (5)

one equation in one unknown, monotone because every cp is positive: `brentq` on [200 K, 6000 K] (the NASA fits' range) finds it. Because ideal-gas enthalpies do not depend on pressure, neither does the complete-combustion T_ad.

**A worked estimate.** Stoichiometric CH₄/air from 298.15 K. Products per mole of fuel: 1.0038 CO₂, 2 H₂O, 7.454 N₂, 0.0888 Ar (10.547 mol, the same 0.29256 kg as the reactants). Splitting (5) into two steps (react at 298.15 K, releasing the LHV, then heat the products),

T_ad = 298.15 + 802 557 J / Σ n_j c̄_p,j,

where c̄_p is the mean cp of the products over [298.15 K, T_ad]. With cp at 298.15 K, Σ n cp = 323.4 J/K and the estimate is 2780 K, far too hot; with cp at 1300 K, 403.9 J/K and 2285 K; the code's exact mean is 395.7 J/K, giving **T_ad = 2326.4 K**. Products' cp rises steeply with temperature (CO₂ and H₂O vibrations wake up), so a room-temperature cp overestimates flame temperatures by hundreds of kelvin (16.Unified 15.5.1 and 15.5.2 make the same comparison).

**Why T_ad peaks at stoichiometric.** Lean of φ = 1 there is less fuel per kilogram of products, and the spare O₂ and N₂ soak up the heat. Rich of φ = 1 there is not enough oxygen to release the fuel's full heating value, so the extra fuel is a diluent too. In this model the peak is a sharp corner exactly at φ = 1 (the code finds the peak at φ = 1.000); with equilibrium chemistry (`P07.11`) the corner rounds off and the peak moves slightly rich.

**The turbine-inlet problem** inverts (5): find φ such that T_ad(φ) = 1700 K, another bracketed root. Because the turbine limit is far below stoichiometric T_ad, gas turbines run lean overall (φ ≈ 0.4) and dilute the hot primary zone with more air downstream.

### 2.4 What generates entropy

A combustor is one of the most irreversible devices in an engine: chemical reaction far from equilibrium, mixing of fuel and air, and heat flowing from the flame to cooler gas. For the turbine-inlet case (Jet-A at 298.15 K, air at 800 K, products at 1700 K, all at 30 bar), a probe using the NASA entropies with each species at its partial pressure gives S_gen = 1141 J/K per kilogram of air. At a dead-state temperature of 288.15 K that is 329 kJ/kg of exergy destroyed, against f × LHV = 1170 kJ/kg of fuel energy released: more than a quarter of the fuel's work potential is lost before the turbine. Lesson `P07.08`'s exergy audit puts this number in its place in the engine.

### 2.5 Where the model stops

- **Dissociation.** Above about 2000 K, CO₂ and H₂O partly fall apart into CO, H₂, OH, H and O, absorbing energy. The 2755 K stoichiometric Jet-A flame is about 150 K too hot; the 4931 K hydrogen–oxygen flame is about 1860 K too hot (`P07.11`). The 1700 K turbine-inlet solve is nearly unaffected (about 6 K).
- **Rich mixtures.** The fixed CO/H₂ order is a convention. Very rich hydrocarbon flames form soot and leave unburned fuel.
- **Fuel state.** Jet-A is a liquid; vaporizing it takes energy this model does not subtract.
- **Kinetics.** Real combustors burn in milliseconds; completeness (combustion efficiency) is usually above 99.5 % at cruise but lower at idle and altitude relight.

## 3. Build the computational model

Open `starter/main.py`. The NASA coefficients for N₂, O₂, Ar, CO₂, H₂O, CO, H₂, CH₄, Jet-A(g) and H₂O(L) (both intervals, 200–6000 K), the steam-table h_f and h_g at 25 °C, and `cp_species`, `h_species`, `molar_mass`, `mixture_mass`, `element_totals` and `main` are given.

### Step 1: stoichiometry

`o2_stoich(fuel)` is (1) from `NASA9[fuel]["elements"]`. `reactants(fuel, phi, oxidizer)` returns a dict per mole of fuel: O₂ = a/φ and, for air, each other component in proportion `AIR_X[k] / AIR_X["O2"]`. `fuel_oxidizer_ratio` is the fuel's molar mass over the oxidizer's mass; `equivalence_ratio` is (2).

### Step 2: products

`complete_products(fuel, phi, oxidizer)` counts the atoms of the reactants with `element_totals` and hands them out: C to CO first, then H to H₂O, then CO to CO₂, the rest of the O as O₂, N as N₂, Ar as Ar. Writing it on atoms, not on the fuel formula, makes the air's own CO₂ and any oxygenated fuel come out right automatically.

```python
el = element_totals(reactants(fuel, phi, oxidizer))
C, H, O = el.get("C", 0.0), el.get("H", 0.0), el.get("O", 0.0)
o_left = O - C                                 # every C to CO first
h2o = min(H / 2.0, o_left); o_left -= h2o
co2 = min(C, o_left)                           # then CO -> CO2
```

### Step 3: energy

`mixture_enthalpy(n, T)` is Σ n_j h_j(T). `reactant_enthalpy` evaluates the fuel at T_fuel and the oxidizer at T_ox. `heating_values` is (4) at 298.15 K. `adiabatic_flame_temperature` is `brentq` on (5) between 200 and 6000 K with tight tolerances; `phi_for_temperature` is `brentq` on φ in (10⁻⁴, 1).

Running `main()` prints (abridged)

```
COMPUTATIONAL ONLY: no combustion, flame or fuel of any kind; every number below is computed.

fuel       O2/mol fuel   f_stoich (air)   AFR_stoich   LHV (MJ/kg)   HHV (MJ/kg)   HHV-LHV vs table h_fg
CH4            2.00      0.05802           17.24       50.027       55.513       5.4859 vs 5.4839 (+0.037 %)
H2             0.50      0.02916           34.29      119.960      141.788       21.8286 vs 21.8205 (+0.037 %)
Jet-A(g)      17.75      0.06818           14.67       43.353       46.377       3.0246 vs 3.0235 (+0.037 %)

complete-combustion adiabatic flame temperatures, reactants at 298.15 K, phi = 1:
  CH4/air: 2326.4 K
  H2/air: 2520.4 K
  Jet-A(g)/air: 2418.1 K
  CH4/O2: 5166.5 K
  H2/O2: 4930.6 K
...
  peak 2754.7 K at phi = 1.000
  turbine inlet 1700 K needs phi = 0.3960, f = 0.02700 (1 kg fuel per 37.0 kg air); energy check f LHV = 1170.4 kJ/kg air
  energy balance residual at that point: -1.16e-07 J per mol fuel (|H_r| = 2.9628e+06)
```

and writes `flame_temperature.png`: T_ad against φ for Jet-A/air with 800 K and 298 K air, CH₄/air and H₂/air, with the turbine-inlet point marked; and the complete-combustion products of Jet-A/air against φ.

### Verify

Run `python -m pytest curriculum/phase-07/lesson-10-combustion`. The tests:

1. **Stoichiometry**: a for CH₄, H₂, Jet-A; air components in the fixed proportions; atoms and mass balanced to 10⁻¹³ for CH₄/air, H₂/air, Jet-A/air, CH₄/O₂, H₂/O₂ at φ = 0.3, 0.7, 1, 1.2, 1.5; the lean, stoichiometric and rich product sets; refusal when soot would form; f_st by hand and f ↔ φ round trips.
2. **Heating values**: Hess's law against the tabulated Hf to 10 J/mol; HHV of methane and hydrogen within 0.1 % of OpenStax Table 5.2; HHV − LHV against the steam-table h_fg within 0.1 % for all three fuels.
3. **Flame temperature**: the energy balance closes to 10⁻⁹ relative for five pairs at three conditions; the peak lies in 0.95 ≤ φ ≤ 1.1; φ → 0 returns the air temperature; hotter inlet air gives a hotter flame by less than the inlet rise; oxygen flames are more than 2000 K hotter than air flames; the two-step identity T_ad − T₀ = m_fuel LHV/(m_products c̄_p) holds to 10⁻⁹; the turbine-inlet φ reproduces 1700 K.
4. **main** writes the figure.

### Compare

NASA CEA (`@nasa-cea-github`) solves the same problem with equilibrium products (its "hp" problem). Do not call it from your tested code; lesson `P07.11` builds the equilibrium solver itself and shows the difference: 102 K for CH₄/air at 1 atm, 154 K for stoichiometric Jet-A/air at 30 bar, 6 K at the turbine-inlet φ.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment. There is no at-home combustion, flame, burning, fuel handling or fuel of any kind in this lesson; every result below is computed.** Write each prediction before running it.

1. **T_ad against φ, Jet-A/air with 800 K air.** φ = 0.3 to 1.5. Predict the peak and its location. *Finding:* 1501.9, 1708.0, 1903.5, 2089.5, 2266.9, 2436.5, 2598.9, **2754.7** K at φ = 0.3 … 1.0, then 2655.2, 2557.9, 2462.8, 2369.6, 2278.3 K at φ = 1.1 … 1.5. The peak is a corner at φ = 1.000. Rich of stoichiometric the CO rises from 3.23 to 11.83 mol per mol of fuel as the CO₂ falls to 0.19.
2. **The turbine inlet.** Solve for 1700 K. *Finding:* φ = 0.3960, f = 0.02700, 37.0 kg of air per kg of fuel; residual −1.16 × 10⁻⁷ J per mole of fuel on 2.96 × 10⁶ J. The same air could burn 2.5 times more fuel before reaching stoichiometric: most of a combustor's air is there to cool, not to burn.
3. **Inlet temperature.** T3 = 298.15, 600, 800, 900 K. Predict whether each 100 K of inlet heat adds 100 K of flame temperature. *Finding:* at φ = 1, 2418.1, 2617.0, 2754.7, 2825.5 K (about 0.69 K per K between 600 and 800 K); at φ = 0.396, 1299.0, 1537.2, 1700.0, 1783.3 K (about 0.81 K per K). Less than one-for-one because the products' cp at flame temperature exceeds the air's cp at inlet temperature, and more so the hotter the flame. A higher-pressure-ratio engine (hotter T3) needs less fuel for the same turbine inlet.
4. **Fuel and oxidizer.** Stoichiometric from 298.15 K: CH₄/air 2326.4 K, H₂/air 2520.4 K, Jet-A/air 2418.1 K, CH₄/O₂ 5166.5 K, H₂/O₂ 4930.6 K. Predict which is hottest. *Finding:* removing nitrogen more than doubles the temperature rise. These oxygen numbers are not believable (they are above where the products would stay molecules); they are the motivation for lesson `P07.11`.
5. **H₂/O₂ and CH₄/O₂ against φ** (0.5, 0.8, 1.0, 1.2, 1.5, 2.0): H₂/O₂ 3831, 4593, 4931, 4452, 3905, 3269 K; CH₄/O₂ 3763, 4713, 5166, 4500, 3704, 2753 K. Rocket engines run fuel-rich partly for this reason: the flame is cooler and the products lighter.

Raw data format: CSV tables `phi, f, T_ad_K, n_CO2, n_H2O, n_CO, n_H2, n_O2` (1); the turbine-inlet line (2); `T3_K, T_ad_phi1_K, T_ad_phi0396_K` (3); `fuel, oxidizer, T_ad_K` (4, 5); plus `flame_temperature.png`.

## 5. Validate and improve

| quantity | independent route | code | difference | cause |
|---|---|---|---|---|
| HHV of CH₄ | OpenStax Table 5.2: 890.8 kJ/mol | 890.57 kJ/mol | 0.03 % | data sets' rounding |
| HHV of H₂ | OpenStax Table 5.2: 285.8 kJ/mol | 285.83 kJ/mol | 0.01 % | — |
| HHV − LHV | steam-table h_fg at 25 °C | NASA Hf(g) − Hf(l) | 0.037 % | two databases |
| CH₄/air T_ad | hand, mean cp 395.7 J/K: 2326.4 K | 2326.4 K | 0 | (same identity) |
| energy balance at T_ad | H_r = H_p | residual 10⁻⁷ J on 3 × 10⁶ J | round-off | — |

The engineering answer: **to take 800 K air to 1700 K the combustor burns Jet-A at φ = 0.396, f = 0.0270 (1 kg per 37.0 kg of air); a stoichiometric zone would reach 2755 K in this model.** The first number is trustworthy within the model's assumptions (dissociation is negligible at 1700 K, and the liquid fuel's heat of vaporization would raise φ slightly); the second is an upper bound that equilibrium lowers by about 150 K.

What to improve: equilibrium products (lesson `P07.11`), then a liquid-fuel enthalpy. What not to improve: the polynomials, which already agree with independent compilations to 0.04 %.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand stoichiometry of Jet-A and the hand LHV/HHV of methane; the hand T_ad estimate; `flame_temperature.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question and saying which number `P07.11` will change.

### Adversarial self-checks

- **What if the air is twice as hot?** The flame does not rise twice as much: 0.69–0.81 K per kelvin of inlet (Experiment 3).
- **What if φ goes to zero?** The products are air at the inlet temperature (tested). What if φ is very rich? The model hands oxygen out in a fixed order and refuses when there is less O than C; reality makes soot first.
- **Which assumption fails first?** Complete combustion, above about 2000 K: stoichiometric flames are 100–2000 K cooler than this model says. Then the gaseous fuel.
- **What measurement could falsify the model?** A combustor-exit temperature traverse with a measured fuel and air flow: at φ = 0.4 the model predicts 1700 K to within a few kelvin of the equilibrium value; a measured exit 50 K cooler would point at heat loss to the liner or incomplete combustion, not at the model.
- **Does the result survive a change of method?** The heating values agree with an independent compilation and the steam table; the flame temperature solves a monotone equation to round-off.
- **Which decision changes if the uncertainty doubles?** At the turbine inlet, 1 % more fuel (Δφ = 0.004) adds about 8 K (dT/dφ ≈ 2010 K): fuel-flow metering, not the thermodynamic model, dominates the uncertainty in turbine inlet temperature.

## Check yourself

1. How many moles of O₂ does C₁₂H₂₃ need for complete combustion? What is its stoichiometric air–fuel ratio in the phase's dry air?
2. Define the equivalence ratio. Is a gas-turbine combustor overall lean or rich, and why?
3. Why can every enthalpy in a reacting mixture be added directly once the NASA polynomials include the enthalpy of formation?
4. What is the difference between the LHV and the HHV, and which one rates a jet engine?
5. Why does the complete-combustion adiabatic flame temperature not depend on pressure?
6. Why does a cp evaluated at room temperature overestimate a flame temperature?
7. Why does T_ad peak near φ = 1?
8. Which of this lesson's results would dissociation change most, and which hardly at all?
