# P07.01  Systems, Properties & State

> Draw the boundary first, then name two independent properties: everything else about the gas follows. For an ideal gas the "everything else" is one function of temperature, and its shape at 1700 K is not the shape it has at 300 K.

| | |
|---|---|
| **Learning outcomes** | 1. Choose a system or a control volume for an aerospace device, classify its properties as intensive or extensive, and say which processes and cycles it undergoes. 2. State the state postulate for a simple compressible substance and decide whether a given pair of properties fixes the state (including why (T, h) does not for an ideal gas). 3. Derive the ideal-gas relations cp − cv = R and s = s°(T) − R ln(p/p°), connect p v = R T to the kinetic-theory pressure of `P04.07`, and derive h and s° from the NASA cp polynomial. 4. Implement and test an ideal-gas property calculator from the NASA 9-coefficient polynomials for single species and for dry air, and check it against the enthalpy of formation, interval continuity, equipartition and the Lienhard air table. 5. Decide when a constant cp is good enough, with numbers. |
| **Prerequisites** | `P04.07` (ideal gas, kinetic theory, equipartition and γ), `P04.02` (work and energy), `P01.05` (integration), `P03.02` (root finding with a bracket). Skill check: you can integrate a polynomial term by term and find a root with `scipy.optimize.brentq`. |
| **Engineering question** | A turbofan combustor heats air from 850 K at its inlet to 1700 K at the turbine inlet (typical, illustrative values). How wrong is the textbook constant cp = 1005 J/(kg·K) for that heat? The answer: the variable-cp enthalpy rise is 1002.9 kJ/kg; constant cp gives 854.3 kJ/kg, 14.8 % low, and the same heat added with constant cp predicts a turbine inlet temperature of 1847.9 K, 148 K too hot. Near room temperature the constant is good to 0.5 %; above about 620 K cp is more than 5 % above 1005. |
| **Assumptions and limits** | Ideal gas (p v = R T; u and h depend on T only); thermally perfect but calorically imperfect (cp varies with T); fixed dry-air composition N₂ 0.7808, O₂ 0.2095, Ar 0.0093, CO₂ 0.0004 by mole (M = 28.9657 g/mol); no dissociation or reaction; NASA fits valid from 200 K (300 K for He) to 6000 K. Fails: at high pressure and low temperature (real-gas effects, lesson 03 measures Z), near condensation, above about 2000 K in air where O₂ starts to dissociate (composition is no longer fixed), and for combustion products, which are not air (lesson 10). |
| **Analytical work** | System and control-volume choices; the state postulate; cp − cv = R; integrating the NASA cp/R polynomial to h/(RT) and s°/R and fixing the constants; the entropy of an ideal-gas mixture with its mixing and pressure terms; the mean cp between two temperatures; the constant-cp error by hand. |
| **Computational deliverable** | `starter/main.py`: the three NASA forms (`cp_over_R`, `h_over_RT`, `s0_over_R`), the species wrappers (`cp_molar`, `h_molar`, `s0_molar`), `fit_jump`, the mixture functions (`mixture_molar_mass`, `mixture_cp`, `mixture_h`, `mixture_s`), dry air (`air_cp`, `air_h`, `air_u`, `air_s`, `air_gamma`), the equation of state (`ideal_gas_pressure`, `specific_volume`, `kinetic_pressure`), the inversions (`temperature_from_h`, `temperature_from_s`), the state postulate (`air_state`), and the comparison (`delta_h_variable`, `delta_h_constant`, `exit_temperature`). Tests (25) prove cp = dh/dT and T ds°/dT = cp for six species, h(298.15 K) = Hf298, continuity at 1000 K, the equipartition limits, dry-air cp within 0.5 % of Lienhard Table A.6, p v = R T equal to ρ v_rms²/3, round trips of the state from eight property pairs, refusal of dependent pairs, and the size of the constant-cp error. |
| **Practical activity** | **computational experiment**: a sweep of the constant-cp error from 200 to 2500 K, a turbine-inlet sensitivity study, a composition study (dry air against an N₂/O₂-only model and moist air), and a continuity audit of the NASA fits. |
| **Equipment** | None (a computer with the course's Python environment). |
| **Safety** | None beyond normal computer use. |
| **Evidence** | `ideal_gas_properties.png`; the printed verification tables; the experiment tables; a half-page note answering the engineering question with the variable-cp and constant-cp numbers and the temperature above which you would no longer use a constant cp. |
| **Assessment** | `assessments/phase-07/grading-rubric.md`; competencies 1–4. Pass: cp − cv = R and the NASA h and s° forms derived; all public tests pass; the combustor Δh reproduced by hand with a mean cp to 1 %; the note gives the constant-cp error at the turbine inlet and the threshold temperature with the evidence. |
| **References** | `@nasa-tp-2002-211556-glenn-coefficients` McBride, Zehe & Gordon, NASA/TP-2002-211556: "Standard States, Reference States, and Fundamental Constants"; "Empirical Equations for Fitting Thermodynamic Functions"; "Assigned Enthalpy Values"; "Using the Coefficients to Obtain Thermodynamic Functions" (required-free; https://ntrs.nasa.gov/citations/20020085330): the fit forms, the 1 bar standard state and the formation-enthalpy datum. `@nasa-cea-github` `data/thermo.inp` (supplementary-free; https://github.com/nasa/cea): the machine-readable file the coefficients were copied from. `@libretexts-yan-engineering-thermodynamics` Ch. 1 "Basic Concepts and Definitions": 1.3 System and surroundings; 1.4 Extensive and intensive properties; 1.5 State, process, and cycle; 3.2 Ideal gas and ideal gas equation of state; 4.2 Internal energy in a system; 5.2 Enthalpy (required-free; https://eng.libretexts.org/Bookshelves/Mechanical_Engineering/Introduction_to_Engineering_Thermodynamics_(Yan)). `@mit-16-unified-thermodynamics-propulsion-notes` 1.2 Definitions and Fundamental Ideas of Thermodynamics; 2.4 Specific Heats (required-free, read online, https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node11.html). `@openstax-university-physics-vol2` 2.3 Heat Capacity and Equipartition of Energy (required-free; the `P04.07` link from molecules to cv). `@lienhard-heat-transfer-textbook` Appendix A, Table A.6 "Thermophysical properties of gases at atmospheric pressure (101325 Pa)" (required-free; https://ahtt.mit.edu/): the independent air cp data. |
| **Extensions** | ★ Add the upper NASA interval of every species in `nasa9_p07` and compute air properties to 6000 K; mark where dissociation (lesson 11) would make the fixed-composition model wrong. ★ Fit your own 7-term cp/R polynomial to the Lienhard air table by least squares (`P02.18`) and compare its h(T) with the NASA mixture. ★ Build a moist-air calculator (mole fraction of H₂O from a dew point) and find how much 1 % water vapour changes cp and R. |
| **Time** | unmeasured |

## 1. The engineering problem

A turbofan's combustor takes air from the compressor at about 850 K and delivers gas to the turbine at about 1700 K (typical, illustrative values for a modern engine; real combustor-exit gas is combustion products, not air, which is lesson 10's business). The fuel flow needed to do that, the turbine work available afterwards and the cooling air the first turbine blades need all depend on one number: how much enthalpy one kilogram of gas gains between those two temperatures. The textbook answer, Δh = cp ΔT with cp = 1005 J/(kg·K), gives 854 kJ/kg. Is that good enough?

It is not, and the reason is molecular. At room temperature an N₂ or O₂ molecule stores energy in three translations and two rotations (`P04.07`: cp = 7/2 R per mole, which is 1005 J/(kg·K) for air). As the gas heats, the molecules begin to vibrate, a new way of storing energy, and cp climbs: by 1700 K it is 1229 J/(kg·K), 22 % above the textbook constant. To do this properly we need cp as a function of temperature, and from it h(T) and s(T, p). NASA has fitted exactly these functions for hundreds of species (McBride, Zehe & Gordon, NASA/TP-2002-211556), and the fits are the property data of this whole phase.

This lesson first sets up the language the rest of thermodynamics uses (system, boundary, property, state, process, cycle, the state postulate) and then builds the property calculator. Success criteria, written first: cp equals dh/dT of the polynomial to 10⁻⁶; h(298.15 K) reproduces each species' tabulated enthalpy of formation to 3 J/mol; the two NASA temperature intervals agree at 1000 K to 10⁻⁵; dry-air cp agrees with an independent table (Lienhard Table A.6) within its stated 0.2–2 % uncertainty; p v = R T equals the kinetic-theory pressure; and any two independent properties reproduce the same state.

## 2. Understand and derive

### 2.1 Systems, boundaries and control volumes

Thermodynamics starts by drawing a line. Everything inside is the **system**; everything outside is the **surroundings**; the line is the **boundary**. What crosses the boundary is what the laws account for.

- A **closed system** (a control mass) exchanges energy but no mass. A sealed pressurant tank with its valve shut is one; so is the gas in a piston-cylinder.
- An **open system** or **control volume** has mass crossing its boundary. A combustor, a compressor stage, a nozzle and a tank that is being vented are control volumes: you draw the box around the hardware and count what flows in and out.
- An **isolated system** exchanges neither mass nor energy. The universe is one; so is "the tank plus everything it exchanges heat with", which lesson 04 uses to count entropy.

```
        closed system                         control volume (combustor)
   +-------------------+                 fuel ──┐
   |  gas, mass fixed  |  <- Q, W            ┌──▼──────────────┐
   |  (tank, valve     |     cross the       │                 │
   |   shut)           |     boundary   air ─►   burning gas    ├─► gas to turbine
   +-------------------+                     │                 │    (mass crosses)
                                             └─────────────────┘
```

The choice is yours, and a good choice makes the problem easy. For the combustor question, the control volume is the combustor, and in steady flow the energy crossing its boundary as heat (from the fuel's chemistry, lesson 10) equals the enthalpy rise of the flow through it: q = h₂ − h₁. Lesson 02 derives that; here we only need h(T).

### 2.2 Properties, states and processes

A **property** is any quantity that has a definite value in a given equilibrium condition, independent of how the system got there: pressure p, temperature T, volume V, mass m, internal energy U, enthalpy H = U + pV, entropy S. Heat Q and work W are not properties; they are energy in transit, and their amounts depend on the path (lesson 02).

- **Extensive** properties scale with the amount of matter (V, U, H, S, m). Divide by mass to get **specific** properties: v = V/m, u = U/m, h = H/m, s = S/m. Lower case means per kilogram in this phase; per mole is written with a "molar" label.
- **Intensive** properties do not scale (p, T, v, ρ = 1/v, and every specific property).

A **state** is the set of all property values; **equilibrium** means they are uniform and unchanging when the system is isolated. A **process** takes the system from one state to another along a **path**; a process slow enough that the system passes through equilibrium states is **quasi-static**. A **cycle** returns to the starting state, so every property change over a cycle is zero (`ΔU_cycle = 0`); this is how one recognises a property.

### 2.3 The state postulate

How many properties must you specify before the rest are fixed? For a **simple compressible system** (one where the only relevant work is p dV: no electric, magnetic, surface-tension or gravitational effects worth counting), the answer is two, provided they are **independent**:

> The state of a simple compressible substance in equilibrium is fixed by two independent intensive properties.

"Independent" is the catch. For an ideal gas, u and h depend on T alone (2.5), so (T, h), (T, u) and (h, u) are each one property counted twice: they fix T but say nothing about p or v. In a two-phase mixture (lesson 03) p and T are tied together by the saturation curve and do not fix the state either. The code's `air_state(**two)` accepts any two of p, T, v, h, u, s and refuses the dependent pairs with a `ValueError`. Starting from (101 325 Pa, 288.15 K) it gives v = 0.816304 m³/kg, h = −15 479.2 J/kg, u = −98 191.2 J/kg and s = 6826.23 J/(kg·K) (the NASA datum makes h negative near room temperature: 2.6), and the pairs (p, v), (p, h), (v, s), (u, s) and (p, s) all return T = 288.150000 K and p = 101 325.0000 Pa.

### 2.4 The ideal-gas equation of state

`P04.07` derived the pressure of a gas of point molecules from the momentum they carry to a wall, p = (1/3) n m v_rms², and found v_rms = √(3RT/M). Per unit mass, ρ = n m, so

**p = ρ v_rms²/3 = ρ R T/M,  i.e.  p v = R T  with  R = R_U/M.**

R_U = 8.314 462 618 J/(mol·K) is the universal gas constant (CODATA 2018); R is the gas constant of the particular gas. For the fixed dry-air model of this phase, M = 28.9657 g/mol and **R = 287.045 J/(kg·K)**. The code checks the equality both ways: at ρ = 1.225 kg/m³ and 288.15 K, ρ v_rms²/3 and ρ R T are both 101 322.15 Pa. (ρ = 1.225 is the ISA's rounded sea-level density; the unrounded value is 1.22503.) The ideal gas assumes point molecules with no forces between them; lesson 03 measures how far real steam and ammonia depart from it with the compressibility factor Z = p v/(R T).

### 2.5 Internal energy, enthalpy and the specific heats

Two specific heats are defined as partial derivatives of properties:

**cv = (∂u/∂T)_v,  cp = (∂h/∂T)_p.**

For an ideal gas u depends on T only (the molecules' energy does not depend on how far apart they are, because they do not interact), so du = cv dT on any path. Then h = u + p v = u + R T also depends on T only, and differentiating,

**cp = cv + R,  γ = cp/cv = cp/(cp − R) > 1.**

Equipartition (`P04.07`) predicts the values: a monatomic gas has cv = 3/2 R and cp = 5/2 R at every temperature; a diatomic gas with rotation active and vibration frozen has cp = 7/2 R (γ = 1.4); as vibration wakes up cp rises towards 9/2 R. The NASA data show all three regimes: cp/R of helium is 2.5000 at 500 K; N₂ has 3.5029 at 300 K and 4.4533 at 3000 K. That is why cp of air is not constant.

### 2.6 The NASA 9-coefficient polynomials

NASA fits cp for each species in two temperature intervals (200–1000 K and 1000–6000 K for the species here) as

**cp/R = a₁T⁻² + a₂T⁻¹ + a₃ + a₄T + a₅T² + a₆T³ + a₇T⁴.**

Enthalpy and entropy follow by integration, because dh = cp dT and, at constant pressure, ds = cp dT/T:

- h(T) = ∫ cp dT + constant. Term by term, ∫a₁T⁻² dT = −a₁T⁻¹, ∫a₂T⁻¹ dT = a₂ ln T, ∫a₃ dT = a₃T, ∫a₄T dT = a₄T²/2, and so on. Dividing by RT gives the tabulated form **h/(RT) = −a₁T⁻² + a₂ ln(T)/T + a₃ + a₄T/2 + a₅T²/3 + a₆T³/4 + a₇T⁴/5 + b₁/T**, with the integration constant b₁.
- s°(T) = ∫ cp/T dT + constant: **s°/R = −a₁T⁻²/2 − a₂T⁻¹ + a₃ ln T + a₄T + a₅T²/2 + a₆T³/3 + a₇T⁴/4 + b₂**, at the standard pressure p° = 1 bar.

The constants fix the **datum**. NASA chooses b₁ so that h(298.15 K) equals the species' enthalpy of formation ΔH_f°: zero for N₂, O₂, Ar, H₂ and He (elements in their reference state), −393 510 J/mol for CO₂, −241 826 J/mol for H₂O vapour. That choice makes reacting mixtures (lesson 10) consistent automatically, and it is why air's h is slightly negative at room temperature: its 0.04 % of CO₂ contributes 0.0004 × (−393 510) = −157 J/mol. Only differences of h matter for non-reacting air. b₂ is fixed by the third law (lesson 04).

Two checks of the code follow from this structure. The fits were made with R = 8.314510 J/(mol·K); with the CODATA value used here, h(298.15 K) of CO₂ is −393 507.76 J/mol, 2.24 J/mol above the tabulated value, which is exactly −393 510 × (8.314 4626/8.314 510 − 1). H₂O is 1.38 J/mol off for the same reason. And the two intervals were fitted to join at 1000 K: the code finds jumps in cp/R, h/(RT) and s°/R of order 10⁻⁸ (largest 4 × 10⁻⁸, for CO₂ and H₂O), far below anything an engineering result depends on.

### 2.7 Entropy of an ideal gas and of a mixture

Lesson 04 defines entropy properly. Here we need its value. From the Gibbs relation T ds = dh − v dp (lesson 04 derives it), with dh = cp dT and v = RT/p,

**ds = cp dT/T − R dp/p,  s(T, p) = s°(T) − R ln(p/p°).**

For a mixture of ideal gases (Dalton's model, lesson 09): each species i behaves as if it alone filled the volume, at its partial pressure p_i = x_i p (x_i is its mole fraction). The mixture's molar properties are mole-fraction averages, and its entropy picks up a term for each species' partial pressure:

**h_mix = Σ x_i h_i(T),  cp,mix = Σ x_i cp,i(T),  s_mix = Σ x_i [s°_i(T) − R_U ln(x_i p/p°)],**

per mole of mixture; divide by M_mix = Σ x_i M_i for per-kilogram values. The −R_U Σ x_i ln x_i part is the **entropy of mixing**, 162.8 J/(kg·K) for dry air: positive, constant at fixed composition, and so cancelling from every difference this phase takes for air. It matters when the composition changes (lessons 04, 09, 11).

### 2.8 Dry air from 200 to 2500 K

The dry-air composition used throughout Phase 07 (N₂ 0.7808, O₂ 0.2095, Ar 0.0093, CO₂ 0.0004 by mole) is this course's own fixed convention, chosen to match standard dry air to four figures so that every lesson uses the same air. It gives:

| T (K) | 200 | 300 | 500 | 1000 | 1500 | 2000 | 2500 |
|---|---|---|---|---|---|---|---|
| cp (J/(kg·K)) | 1002.4 | 1004.8 | 1029.6 | 1141.1 | 1211.0 | 1250.4 | 1276.1 |
| γ | 1.4013 | 1.3999 | 1.3866 | 1.3361 | 1.3107 | 1.2980 | 1.2902 |

Below about 400 K cp is flat and 1005 is right to 1 %. It passes 1005 × 1.01 at 413 K, 1.05 × at 618 K, 1.10 × at 830 K and 1.20 × at 1452 K. (Above about 2000 K real air also begins to dissociate, which no fixed-composition model captures; the table's last two columns are for the frozen mixture.)

### 2.9 The combustor by hand

For steady flow through the combustor, q = h(1700 K) − h(850 K) (lesson 02). With the NASA mixture:

**Δh = 1002.91 kJ/kg,  mean cp = Δh/ΔT = 1002.91/850 = 1179.9 J/(kg·K).**

The mean cp is close to cp at the mean temperature 1275 K (1185.1 J/(kg·K)) because cp is nearly linear over this range: a by-hand estimate with cp(T_mean) is 0.4 % high. The constant 1005 gives 854.25 kJ/kg, **14.8 % low**. Turned around: if the combustor really adds 1002.91 kJ/kg and the designer predicts its exit temperature with constant cp, the prediction is 850 + 1002 910/1005 = **1847.9 K instead of 1700 K**, a 148 K error at the place in the engine where the turbine blades are most sensitive to temperature. For the compression from 288.15 to 850 K the error is −4.1 % (589.06 against 564.66 kJ/kg); for a turbine expansion from 1700 to 1000 K, −15.6 %.

### 2.10 Where the model fails

- **High pressure, low temperature.** Molecules attract and occupy volume; Z departs from 1. Steam at 4 MPa and 300 °C has Z = 0.89 (lesson 03). A 30 MPa pressurant tank (lesson 02) is outside the ideal-gas model's comfort zone.
- **Dissociation and ionisation.** Above about 2000 K O₂ begins to split; air behind a re-entry shock is a reacting mixture whose composition depends on p and T (lesson 11). The fits for O, N and NO exist in the same database; the fixed composition does not.
- **Combustion products.** Gas leaving a combustor contains CO₂ and H₂O, whose cp is larger than air's; lesson 10 computes it species by species.
- **Range of the fits.** 200–6000 K here (helium from 300 K). The code refuses temperatures outside them rather than extrapolate.
- **Non-equilibrium.** In a fast expansion (a rocket nozzle) vibration can lag the temperature; cp(T) assumes the molecules are in equilibrium at each T.

## 3. Build the computational model

Open `starter/main.py`. The NASA coefficients (copied unrounded from the NASA CEA `thermo.inp` file), the Phase 07 dry-air composition, the Lienhard Table A.6 air cp values, `coefficients` (which picks the interval), `_per_T` (which applies a formula element by element), `molar_mass`, `rms_speed` and `main` are given.

### Step 1: the three NASA forms

`cp_over_R(a, T)`, `h_over_RT(a, b1, T)` and `s0_over_R(a, b2, T)` are the three formulas of 2.6, one line each. Watch the signs of the a₁ and a₂ terms and the powers' divisors.

### Step 2: species properties

```python
def cp_molar(species, T, interval=None):
    return _per_T(lambda a, b, t: R_U * cp_over_R(a, t), species, T, interval)
```

`h_molar` multiplies `h_over_RT` by R_U T; `s0_molar` multiplies `s0_over_R` by R_U. `fit_jump(species)` evaluates all three at 1000 K with `interval=0` and `interval=1` and returns the differences made dimensionless.

### Step 3: mixtures and dry air

`mixture_molar_mass`, `mixture_cp`, `mixture_h` are mole-fraction sums divided by M_mix. `mixture_s` adds −R_U ln(x_i p/p°) for each species (skip species with x_i = 0). `air_cp`, `air_h`, `air_s` call them with the `AIR` dictionary; `air_u = air_h − R_air T`; `air_gamma = cp/(cp − R_air)`.

### Step 4: the equation of state

`ideal_gas_pressure(ρ, T, M)`, `specific_volume(p, T, M)` and `kinetic_pressure(ρ, T, M) = ρ v_rms²/3` with the given `rms_speed`.

### Step 5: inversions and the state postulate

h(T) and s(T, p) are strictly increasing in T (cp > 0), so `brentq` on 200–6000 K finds T from h, or from s at a given p, without any risk of a wrong root. `air_state(**given)` first rejects anything but exactly two properties and the three dependent pairs, then finds T (directly, from h or u by root finding, from p v/R, from (p, s) by root finding, or from (v, s) by root finding on s(T, RT/v)), then p, then everything else.

### Step 6: the comparison

`delta_h_variable(T1, T2) = h(T2) − h(T1)`, `delta_h_constant = cp ΔT`, and `exit_temperature(T1, q)` inverts h(T1) + q.

Running `main()` prints

```
dry air: M = 28.9657 g/mol, R = 287.045 J/(kg K)
h(298.15 K) against the tabulated enthalpy of formation (J/mol):
   CO2: h =   -393507.76, Hf298 =   -393510.00, difference +2.24
   H2O: h =   -241824.62, Hf298 =   -241826.00, difference +1.38
fit continuity at 1000 K (upper minus lower interval): dcp/R, dh/(RT), ds/R
    N2: +8.47e-09 +4.66e-10 +7.31e-09
   CO2: +2.46e-08 +1.44e-08 +1.51e-08
air cp: NASA mixture against Lienhard Table A.6 (J/(kg K))
     200 K: NASA  1002.4, Lienhard   1007, deviation -0.46 %
     300 K: NASA  1004.8, Lienhard   1006, deviation -0.12 %
    1000 K: NASA  1141.1, Lienhard   1141, deviation +0.01 %
    1500 K: NASA  1211.0, Lienhard   1210, deviation +0.09 %
  largest |deviation| over all 22 rows: 0.46 %
kinetic theory: rho v_rms^2 / 3 = 101322.15 Pa; rho R T = 101322.15 Pa
constant cp = 1005 J/(kg K) against variable cp
  cp(T) of air: 300 K 1004.8, 850 K 1110.1, 1700 K 1229.1 J/(kg K); gamma: 1.3999, 1.3487, 1.3047
  compression 288.15 -> 850 K: variable   589.06 kJ/kg, constant   564.66 kJ/kg, constant-cp error -4.14 %
  heat addition 850 -> 1700 K: variable  1002.91 kJ/kg, constant   854.25 kJ/kg, constant-cp error -14.82 %
  turbine expansion 1700 -> 1000 K: variable  -834.01 kJ/kg, constant  -703.50 kJ/kg, constant-cp error -15.65 %
  the same heat 1002.91 kJ/kg added at 850 K with constant cp gives 1847.9 K instead of 1700.0 K
```

(abridged) and writes `ideal_gas_properties.png`: cp of air against T with the Lienhard points and the constant 1005; γ(T); and the constant-cp enthalpy error from 298.15 K.

### Verify

`python -m pytest curriculum/phase-07/lesson-01-systems-properties-state` runs 25 tests:

1. **Internal consistency of the fits**: for N₂, O₂, Ar, CO₂, H₂O and H₂ at 350, 800, 1500 and 4000 K, a central difference of h equals cp, and T times a central difference of s° equals cp, to 10⁻⁶ (this catches a wrong sign or divisor in any of the three forms).
2. **The datum**: h(298.15 K) = Hf298 within 3 J/mol (the stated R mismatch).
3. **Continuity**: the two intervals agree at 1000 K to 10⁻⁵ in cp/R, h/(RT) and s°/R.
4. **Equipartition**: cp = 5/2 R for He and Ar; 7/2 R for N₂ at 300 K within 0.2 %; cp of N₂ rising monotonically towards 9/2 R.
5. **Dry air**: M and R; cp within 0.5 % of all 22 Lienhard rows (0.2 % above 400 K); h − u = RT; the pressure term of s; a positive mixing entropy; a one-species mixture reducing to the species.
6. **Equation of state**: p v = R T and ρ v_rms²/3 agree to 10⁻¹²; sea-level density gives 101 325 Pa.
7. **State postulate**: eight independent pairs reproduce the (2.5 MPa, 820 K) state to 10⁻⁸; dependent pairs and wrong counts raise `ValueError`; h and s inversions round-trip on both sides of 1000 K.
8. **The comparison**: constant cp within 0.5 % near 300 K and 10–20 % low between 850 and 1700 K; Δh equals the trapezoid integral of cp; `exit_temperature` inverts the heat addition and constant cp overpredicts it.

### Compare

The independent check is a different data set: Lienhard & Lienhard, *A Heat Transfer Textbook*, Table A.6, which tabulates cp of air at 1 atm from separate measurements and correlations, with a stated cp uncertainty of 0.2–2 %. The NASA mixture agrees within 0.12 % from 300 K up and within 0.46 % everywhere. The largest difference is at 200 K, where Lienhard's 1007 is real air at 1 atm and the NASA value is the zero-pressure ideal gas: at low temperature the attraction between molecules adds a little to cp, the same real-gas effect that lesson 03 quantifies for steam. At 1000–1500 K the two agree to 0.1 %. Two sources built in different ways agree within the smaller one's uncertainty; the calculator is right.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before you run it.

1. **Where does constant cp fail?** Sweep T from 200 to 2500 K and tabulate cp/1005 − 1 and the error of 1005 (T − 298.15) against h(T) − h(298.15). *Predict* the temperature at which each error first reaches 5 %. *Finding:* cp passes 1 %, 5 %, 10 % and 20 % above 1005 at 413, 618, 830 and 1452 K; the integrated Δh error from 298.15 K reaches 1 % at 510 K and 5 % at 930 K (the integral lags the local value because the cold part of the range is accurate).
2. **Turbine-inlet sensitivity.** For a combustor inlet of 850 K and turbine inlet temperatures of 1400, 1500, 1600, 1700 and 1800 K, compute q with variable cp, then the exit temperature constant cp would predict for that q. *Predict* whether the temperature error grows linearly with the turbine inlet temperature. *Finding:* q = 638.24, 758.83, 880.42, 1002.91 and 1126.22 kJ/kg; constant cp then predicts 1485.1, 1605.1, 1726.0, 1847.9 and 1970.6 K, errors of 85.1, 105.1, 126.0, 147.9 and 170.6 K. The error grows faster than linearly, because each extra kelvin of turbine inlet temperature is added where cp is largest.
3. **Composition.** Replace the Phase 07 dry air by an N₂/O₂-only model (0.79/0.21) and by moist air with 1 % water vapour by mole (the dry fractions scaled by 0.99). *Predict* which property changes most. *Finding:* the N₂/O₂ model has M = 28.8503 g/mol and R = 288.19 J/(kg·K) (0.4 % high, because the heavy argon is gone) and cp(300 K) = 1011.4; moist air has M = 28.8562 g/mol and cp(300 K) = 1010.2; Δh from 850 to 1700 K is 1010.3 and 1009.8 kJ/kg against 1002.9. Composition shifts properties by under 1 %; temperature shifts cp by 22 %.
4. **Continuity audit.** Run `fit_jump` for all seven species. *Predict* the largest jump. *Finding:* every jump is below 5 × 10⁻⁸ in cp/R, h/(RT) and s°/R.

Raw data format: CSV tables `T_K, cp_JkgK, cp_err, dh_err` (1); `T4_K, q_kJkg, T4_const_K` (2); `model, M_gmol, R_JkgK, cp300, dh_850_1700_kJkg` (3); `species, dcp, dh, ds` (4); plus `ideal_gas_properties.png`.

## 5. Validate and improve

| quantity | independent value | code | difference | cause |
|---|---|---|---|---|
| cp = dh/dT (N₂, 600 K, ΔT = ±10 K) | 30.109058 J/(mol·K) | 30.109254 | 6.5 × 10⁻⁶ | finite-difference truncation |
| h(298.15 K), CO₂ | Hf298 = −393 510 J/mol | −393 507.76 | +2.24 J/mol | fit R 8.314510 against CODATA |
| interval jump at 1000 K, cp/R | 0 | ≤ 2.5 × 10⁻⁸ | — | fit rounding |
| cp/R of He, N₂ at 300 K | 5/2, 7/2 (equipartition) | 2.5000, 3.5029 | 0, 0.08 % | N₂ vibration just beginning |
| dry-air cp, 300–1500 K | Lienhard Table A.6 | within 0.12 % | — | inside the table's 0.2–2 % |
| dry-air cp, 200 K | 1007 (Lienhard, 1 atm) | 1002.4 | −0.46 % | real-gas cp at 1 atm, 200 K |
| p, sea level | ρ v_rms²/3 = 101 322.15 Pa | ρRT = 101 322.15 Pa | < 10⁻¹² | identity |

The engineering answer: **between 850 K and 1700 K the variable-cp enthalpy rise of air is 1002.9 kJ/kg; constant cp = 1005 J/(kg·K) gives 854.3 kJ/kg, 14.8 % low, and would put the turbine inlet 148 K too hot for the same heat. Use variable cp above about 600 K, where cp is more than 5 % above 1005.** The model's own uncertainty (0.1 % against Lienhard) is two orders of magnitude smaller than the error it removes. What to improve next is not the data but the gas: combustion products (lesson 10) and, above 2000 K, dissociation (lesson 11).

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of cp − cv = R, the NASA h and s° forms and the mixture entropy; `ideal_gas_properties.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question.

### Adversarial self-checks

- **What if the temperature doubles?** From 300 to 600 K, cp of air rises from 1004.8 to about 1051 J/(kg·K) and γ falls from 1.400 to about 1.38; from 850 to 1700 K, cp rises by 11 % (1110.1 to 1229.1). Doubling T matters more the hotter you start, because vibration is the mechanism.
- **What if the pressure doubles?** For the ideal gas, h and cp do not change; s falls by R ln 2 = 199.0 J/(kg·K); v halves. Whether the ideal gas itself survives depends on how close you are to condensation (lesson 03).
- **What are the limiting cases?** Helium: cp = 5/2 R exactly at every T in the fit. Air near 300 K: 7/2 R. Air as T → ∞ in the frozen model: 9/2 R plus a little, but in reality dissociation intervenes.
- **Which assumption fails first?** For a turbine, the composition (products, not air); for a pressurant tank, the ideal gas (Z ≠ 1); for re-entry, the fixed composition (dissociation).
- **What measurement could falsify the model?** A calorimeter measurement of air's cp at 1000 K differing from 1141 J/(kg·K) by more than about 1 %; a speed of sound in hot air at 1500 K that does not match √(γRT) with γ = 1.311.
- **Does the result survive a change of method?** cp from the polynomial and from differentiating h; h from the closed form and from integrating cp with the trapezoid rule (agreeing to 10⁻⁷); Lienhard's independent table. Yes.

## Check yourself

1. Is a rocket engine's combustion chamber in steady firing better treated as a closed system or a control volume? Why?
2. Classify as intensive or extensive: V, v, T, U, h, m, p.
3. State the state postulate. Why do T and h not fix the state of an ideal gas?
4. Show that cp − cv = R for an ideal gas, starting from h = u + p v.
5. Why does cp of nitrogen rise between 300 K and 3000 K? Towards what value?
6. What does the constant b₁ of the NASA enthalpy fit represent, and why is h of air negative at 288 K?
7. Write s(T, p) of an ideal gas in terms of s°(T). By how much does s change when p doubles at constant T?
8. Above roughly what temperature is cp of air more than 5 % above 1005 J/(kg·K)?
