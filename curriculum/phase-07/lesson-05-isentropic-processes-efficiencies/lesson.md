# P07.05  Isentropic Processes & Efficiencies

> The isentropic process is the yardstick, not the machine. Efficiency says how far the machine falls short, and the entropy it generates says where the shortfall went.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the isentropic relation s°(T₂) − s°(T₁) = R ln(p₂/p₁) for an ideal gas from the Gibbs equation, and reduce it to T₂ = T₁ π^((γ−1)/γ) for constant cp. 2. Define the isentropic efficiency of a compressor and a turbine through enthalpies, and compute exit temperatures, work and entropy generation with variable cp from the NASA polynomials. 3. Derive the polytropic (small-stage) efficiency, its closed form for constant cp, and the result that a compressor's isentropic efficiency is below its polytropic efficiency while a turbine's is above. 4. Implement and test a compressor efficiency calculator that inverts s°(T) and h(T) by Newton's method, builds a multi-stage compressor and shows that many small stages reproduce the polytropic closed form at first order. 5. Decide when constant cp is good enough for a compressor or turbine calculation. |
| **Prerequisites** | `P07.01` (NASA 9-coefficient cp, h and s° for dry air), `P07.02` (steady-flow energy equation), `P07.04` (entropy generation, S_gen ≥ 0), `P03.02` (Newton's method), `P03.04` (numerical integration), `P01.05` (integration). Skill check: you can write the steady-flow energy equation for an adiabatic device and solve f(T) = c for T with Newton's method. |
| **Engineering question** | Sea-level air (288.15 K) enters the core compressor of a modern turbofan, which raises its pressure 40 times with a polytropic efficiency of 0.90 (typical, illustrative): how hot does the air leave, how much work does each kilogram take, what isentropic efficiency does that correspond to, and how wrong is a constant-cp (γ = 1.4) calculation? The answers: 891.2 K, 635.0 kJ/kg, η_s = 0.843; constant cp predicts 929.4 K (38 K too hot) and 1.45 % too much work, while both models agree exactly on the entropy generated, 117.65 J/(kg·K). |
| **Assumptions and limits** | Steady flow; adiabatic machines (no heat through the casing); kinetic and potential energy changes neglected, or all states read as stagnation states; ideal gas with frozen composition (the Phase 07 dry-air model, NASA polynomials, 200–6000 K); efficiencies constant along the machine. Real compressors depart from these through casing heat transfer, cooling-air bleeds, tip leakage, humidity, efficiencies that vary across the operating map, and surge and stall, none of which this model represents. |
| **Analytical work** | The Gibbs equation and the ideal-gas isentrope with variable and constant cp; isentropic efficiency of compressors and turbines; entropy generation in an adiabatic machine; the polytropic efficiency from dh = v dp/η_p and its closed form; the entropy generated along a polytropic path, R ln π (1/η_p − 1); the preheat and reheat effects; the integral of T ds as friction heat. |
| **Computational deliverable** | `starter/main.py`: `temperature_from_s0`, `temperature_from_h`, `isentropic_temperature_closed_form`, `isentropic_temperature`, `compressor_exit_temperature`, `turbine_exit_temperature`, `polytropic_exit_temperature`, `isentropic_efficiency`, `polytropic_efficiency`, `isentropic_from_polytropic`, `isentropic_from_polytropic_turbine`, `specific_work`, `entropy_generation`, `multistage_compressor`, `first_law_integrals` (the NASA data, `cp_air`, `h_air`, `s0_air`, the `Gas` container, `constant_cp_gas`, `_newton_invert`, `observed_order` and `main` are given). Tests (15) prove: the inversions are exact inverses; the variable-cp isentrope reduces to T₁π^((γ−1)/γ) for a constant-cp gas; η = 1 generates no entropy and η < 1 always generates some, more as η falls; the efficiency definitions round-trip; the polytropic closed forms and the compressor–turbine asymmetry; s_gen along a polytropic path is independent of cp; many small stages converge to the polytropic closed form at first order; and the first law closes three ways. |
| **Practical activity** | **computational experiment**: constant against variable cp across pressure ratios 1.5–60; a stage-count convergence study of a multi-stage compressor towards the polytropic limit; a compressor and a turbine at the same polytropic efficiency; and the first law closed three independent ways. |
| **Equipment** | None. A computer with the course's Python environment. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `compressor_efficiency.png`, the four experiment tables, the derivations of the isentrope, the polytropic closed form and s_gen = R ln π (1/η_p − 1), and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-07/grading-rubric.md`; competencies 1–4. Pass: the isentropic and polytropic relations are derived from the Gibbs equation; all public tests pass; the core-compressor exit temperature, work and isentropic efficiency are reproduced by hand with constant cp to three figures and explained against the variable-cp result; the note states which assumption limits the answer. |
| **References** | `@libretexts-yan-engineering-thermodynamics` 6.8 "Specific entropy of a state" (6.7.3 Determining the specific entropy of ideal gases; 6.7.4 Isentropic relations for an ideal gas), https://eng.libretexts.org/Bookshelves/Mechanical_Engineering/Introduction_to_Engineering_Thermodynamics_(Yan)/06%3A_Entropy_and_the_Second_Law_of_Thermodynamics/6.08%3A_Specific_entropy_of_a_state, and 5.4 "Applications of the mass and energy conservation equations in steady flow devices" (its opening subsection on turbines and compressors) (required-free; the isentrope with s°(T) and with constant cp, and the steady-flow energy equation for a compressor). `@mit-16-unified-thermodynamics-propulsion-notes` 2.4 "Specific Heats" (2.4.2 Reversible adiabatic processes for an ideal gas) and 6.7 "Examples of Lost Work in Engineering Processes", subsection "Adiabatic Efficiency of a Propulsion System Component (Turbine)", https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node50.html (required-free; read online, link only). `@mit-ocw-16-050-thermal-energy` Part 1C, 1.C.7 b) "Adiabatic Efficiency of a Propulsion System Component (Turbine)", https://ocw.mit.edu/courses/16-050-thermal-energy-fall-2002/resources/05_part1c/ (supplementary-free; the openly licensed earlier version of the same notes). `@nasa-tp-2002-211556-glenn-coefficients` "Empirical Equations for Fitting Thermodynamic Functions", https://ntrs.nasa.gov/citations/20020085330 (required-free; the s°(T)/R fit used for the variable-cp isentrope). `@openstax-university-physics-vol2` 3.6 "Adiabatic Processes for an Ideal Gas" (supplementary-free; pV^γ = constant from the physics side). `@lienhard-heat-transfer-textbook` Appendix A, Table A.6 "Thermophysical properties of gases at atmospheric pressure (101325 Pa)" (supplementary-free; the cp check in part 3). |
| **Extensions** | ★ Add a turbine cooling-air bleed: take a fraction of the compressor flow from stage k and return it to the turbine, and compute what it costs in compressor work. ★ Add water vapour to the air (`P07.09`) and find how much a humid day changes the core compressor's exit temperature at the same pressure ratio. ★ Fit a compressor "map": let the stage efficiency fall quadratically away from a design flow coefficient and compute the overall efficiency along a working line. |
| **Time** | unmeasured |

## 1. The engineering problem

Every gas turbine starts by squeezing air. In a modern turbofan the fan, the booster and the high-pressure compressor together raise the pressure of the core air by a factor of about 40 before it reaches the combustor; in a small turboprop or an auxiliary power unit (APU) the factor is 8 to 15. That compression is expensive: it absorbs well over half of what the turbine produces (`P07.06` puts a number on it), so a few percent of compressor efficiency is worth a few percent of fuel burned for the life of the engine. The air also arrives at the combustor hot, and how hot matters twice: it sets how much fuel is needed to reach the turbine inlet temperature, and it sets what material the last compressor stages and the discs that carry them must be made of.

So the designer needs three numbers for a given pressure ratio: the exit temperature, the work per kilogram, and the efficiency, with a definition of efficiency that means the same thing for a one-stage APU compressor and a fourteen-stage core. The reference process is the **isentropic** one (adiabatic and reversible): no real machine achieves it, but it is the best any adiabatic machine could do, and it is computable from properties alone. Efficiency compares the real machine with it. There are two common definitions, isentropic and polytropic, and they differ by several points at high pressure ratio; using the wrong one, or mixing them, is a classic source of a 30 K error in a cycle calculation.

Success criteria, written first: the variable-cp isentrope reduces to T₁π^((γ−1)/γ) for constant cp to 10⁻¹¹; η = 1 gives s_gen = 0 and any η < 1 gives s_gen > 0; the closed-form relation between isentropic and polytropic efficiency is reproduced by the solver; a compressor of many small stages converges to the polytropic closed form at first order in 1/n; and the work equals h₂ − h₁, the integral of cp dT, and the integral of v dp + T ds along the path.

## 2. Understand and derive

### 2.1 The first law for a compressor

Take a control volume around the compressor, in steady flow, with one inlet and one exit (`P07.02`). **Assumptions:** no heat crosses the casing (adiabatic), and the kinetic and potential energy changes between inlet and exit are negligible, or equivalently every state is a stagnation state, which is how turbomachinery engineers read them. The steady-flow energy equation per unit mass then says the shaft work put in is the enthalpy rise:

**w = h₂ − h₁.**   (1)

Nothing about efficiency or entropy has been used yet. The whole difficulty is predicting h₂.

### 2.2 The isentropic process

The **Gibbs equation**, which combines the first and second laws for a simple compressible substance, is T ds = dh − v dp. For an ideal gas, dh = cp(T) dT and v = RT/p, so

ds = cp(T) dT/T − R dp/p.

Define the **standard-state entropy** s°(T) = ∫ cp(T) dT/T from a reference temperature, at the reference pressure p° = 1 bar. That is what the NASA polynomials tabulate (`P07.01`). Then s(T, p) = s°(T) − R ln(p/p°), and a process with s₂ = s₁ satisfies

**s°(T₂ₛ) − s°(T₁) = R ln(p₂/p₁) = R ln π.**   (2)

The subscript s marks the isentropic exit state. Equation (2) holds for any ideal gas, whatever cp(T) does; it is solved for T₂ₛ by inverting s°(T), which rises monotonically because ds°/dT = cp/T > 0. If cp is constant, s°(T) = cp ln(T/T_ref) and (2) becomes

**T₂ₛ/T₁ = π^(R/cp) = π^((γ−1)/γ),**   (3)

using cp = γR/(γ − 1). For air with γ = 1.4 the exponent is 0.2857. Units check: R/cp is dimensionless, as an exponent must be.

The dry-air model has cp = 1004.8 J/(kg·K) and γ = 1.3999 at 300 K, but cp = 1120.9 J/(kg·K) and γ = 1.3442 at 900 K: vibration of the N₂ and O₂ molecules wakes up as the air heats (`P04.07`). A rising cp means each kelvin of temperature rise "costs" more entropy, so the variable-cp isentrope runs cooler than (3). At π = 40 from 288.15 K, (3) gives T₂ₛ = 826.70 K and (2) with NASA air gives 801.53 K: 25 K apart before any inefficiency enters.

### 2.3 Isentropic efficiency

A real adiabatic compressor generates entropy, through viscous friction in the blade boundary layers, shocks, tip leakage and mixing, so its exit entropy is higher than its inlet's and, at the same exit pressure, its exit temperature is higher than T₂ₛ. The **isentropic efficiency** compares the ideal work with the actual work for the same pressure ratio:

**η_c = (h₂ₛ − h₁)/(h₂ − h₁)  (compressor),    η_t = (h₁ − h₂)/(h₁ − h₂ₛ)  (turbine).**   (4)

The ratio is upside down for the turbine because a turbine produces less work than the ideal, while a compressor needs more. Both are defined through **enthalpy**, not temperature: with variable cp the two differ, and the work is what the efficiency is meant to describe. To find the exit state from η_c: compute h₂ₛ from (2), set h₂ = h₁ + (h₂ₛ − h₁)/η_c, and invert h(T).

```
   h                                   2 (actual)
   ^                                  /
   |         p2 isobar  ............ 2s ......... 
   |                    :          /|
   |                    :  w_s   /  |  w = h2 - h1 > w_s = h2s - h1
   |                    :      /    |
   |   p1 isobar .......:.... 1     |
   |                          |---->|  s2 - s1 = s_gen > 0
   +---------------------------------------------> s
```

### 2.4 Entropy generation

For an adiabatic control volume in steady flow, the entropy balance (`P07.04`) has no heat term, so all the entropy rise is generated inside:

**s_gen = s₂ − s₁ = s°(T₂) − s°(T₁) − R ln π ≥ 0.**   (5)

η = 1 makes T₂ = T₂ₛ and, by (2), s_gen = 0; any η < 1 makes h₂ > h₂ₛ, so T₂ > T₂ₛ and s_gen > 0. The tests check both. What the generated entropy costs in work is the subject of `P07.08`: multiplied by the ambient temperature T₀ it is the **exergy destroyed**, the work that could have been obtained and was not.

### 2.5 Polytropic efficiency

A 14-stage compressor at π = 40 is fourteen small compressors in series, each with a pressure ratio near 1.3. Their stage efficiencies may all be 0.90, but the overall isentropic efficiency is not 0.90. To get a measure of aerodynamic quality that does not depend on the pressure ratio, take the limit of a stage so small that its pressure rises by dp. Define the **polytropic efficiency** η_p by

**dh = v dp / η_p  (compressor),    dh = η_p v dp  (turbine).**   (6)

v dp is the reversible work of the infinitesimal step (from the Gibbs equation with ds = 0); the actual enthalpy rise is larger by 1/η_p. For an ideal gas, cp dT = (RT/p) dp/η_p, so cp dT/T = (R/η_p) d ln p, which integrates to

**s°(T₂) − s°(T₁) = (R/η_p) ln π  (compressor),    s°(T₂) − s°(T₁) = η_p R ln π  (turbine, π < 1).**   (7)

With constant cp, T₂/T₁ = π^(m/η_p) with m = (γ − 1)/γ, and substituting into (4):

**η_c = (π^m − 1)/(π^(m/η_p) − 1),    η_t = (1 − π_t^(−η_p m))/(1 − π_t^(−m)),**   (8)

with π_t = p₁/p₂ the turbine's expansion ratio. As π → 1 both tend to η_p. As π grows, the compressor's isentropic efficiency falls below η_p and the turbine's rises above it: at η_p = 0.90 and π = 40 the compressor's η_c is 0.840 (constant cp) and at an expansion ratio of 12 the turbine's η_t is 0.929.

**Why the asymmetry: preheat and reheat.** In a multi-stage compressor the friction of the early stages heats the air, and the later stages must then compress hotter air, which needs more work for the same pressure ratio (the isentropic work of a stage is roughly cpT_in(π_st^m − 1), proportional to its inlet temperature). The losses compound, so the overall machine is worse than its stages: the **preheat** effect. In a turbine the friction heat of the early stages leaves hotter gas for the later stages to expand, and some of the lost work is recovered: the **reheat** effect. The overall turbine is better than its stages.

### 2.6 What the friction did: two integrals

Along the polytropic path, the Gibbs equation splits the work into two parts: dh = v dp + T ds. With (6), T ds = dh − v dp = (1 − η_p) dh. So

∫ T ds = (1 − η_p) w,    ∫ v dp = η_p w.   (9)

∫T ds is the heat that friction generated inside the gas, 10 % of the work at η_p = 0.90. And (5) with (7) gives the entropy generated along a polytropic compression:

**s_gen = R ln π (1/η_p − 1),**   (10)

which does **not** depend on cp at all. At π = 40 and η_p = 0.90 it is 117.65 J/(kg·K) for constant cp and for NASA air alike, even though their exit temperatures differ by 38 K. It is also why, in (9), the friction heat ∫T ds = 63.50 kJ/kg is larger than T₀ s_gen = 33.90 kJ/kg: the friction heat is generated at temperatures well above T₀, and heat at high temperature still has some capacity to do work. Only T₀ s_gen is truly lost (`P07.08`).

### 2.7 Where the model stops

- **Heat through the casing.** Small compressors lose heat to their surroundings, which lowers T₂ and makes a temperature-based efficiency measurement read high. Large ones are close to adiabatic.
- **Bleeds and cooling air.** Real cores extract air from intermediate stages for turbine cooling and the cabin; the flow is not constant through the machine.
- **Humidity.** Water vapour changes R and cp a little (Extension 2).
- **Off design.** η depends on the flow and speed; at low flow the compressor surges. A constant η is a design-point statement.
- **Frozen ideal gas.** Excellent for air in a compressor (pressures to 4 MPa at 900 K are far from the critical point; nothing dissociates). In a turbine burning fuel the gas is combustion products, with a different cp (`P07.08` handles that).

## 3. Build the computational model

Open `starter/main.py`. The NASA coefficients for N₂, O₂, Ar and CO₂, the dry-air functions `cp_air`, `h_air`, `s0_air`, the `Gas` container (R plus cp, h and s° as functions of T), `constant_cp_gas(gamma)`, the vectorised Newton solver `_newton_invert`, `observed_order` and `main` are given. Every function you write takes a `Gas`, so the same code runs with NASA air and with a constant-cp gas; that is what makes test 2 possible.

### Step 1: invert the properties

`temperature_from_s0(gas, s0_target)` and `temperature_from_h(gas, h_target)` call `_newton_invert` with the right derivative: ds°/dT = cp/T and dh/dT = cp. Newton converges in a handful of iterations because both functions are smooth and monotonic.

### Step 2: isentropes and efficiencies

`isentropic_temperature(gas, T1, pi)` is (2): `temperature_from_s0(gas, gas.s0(T1) + gas.R * np.log(pi))`. It works for compression (π > 1) and expansion (π < 1). `compressor_exit_temperature` and `turbine_exit_temperature` apply (4) through enthalpy. `polytropic_exit_temperature` is (7), choosing 1/η_p or η_p by the sign of ln π. `isentropic_efficiency` and `polytropic_efficiency` invert those definitions from measured inlet and exit temperatures; the closed forms (3) and (8) are one line each.

### Step 3: work, entropy and stages

`specific_work` is (1) and `entropy_generation` is (5). `multistage_compressor(gas, T1, pi_total, n, eta_stage)` divides the pressure ratio into n equal stage ratios π^(1/n) and calls `compressor_exit_temperature` n times, each stage starting from the previous stage's exit temperature:

```python
pi_st = pi_total ** (1.0 / n_stages)
T = [float(T1)]
for _ in range(n_stages):
    T.append(compressor_exit_temperature(gas, T[-1], pi_st, eta_stage))
```

### Step 4: the first law, three ways

`first_law_integrals(gas, T1, pi, eta_p)` walks the polytropic path on 2001 points in ln p, computes T at each from (7) (vectorised), and returns h₂ − h₁, ∫cp dT by `scipy.integrate.quad`, ∫v dp = ∫RT d ln p and ∫T ds by the trapezoid rule. The first two must agree to quadrature precision and the last two must sum to the first (Gibbs).

Running `main()` prints

```
dry-air model: M = 28.9657 g/mol, R = 287.045 J/(kg K), cp(300 K) = 1004.8, cp(900 K) = 1120.9 J/(kg K), gamma(300 K) = 1.3999, gamma(900 K) = 1.3442

aero-engine core compressor (typical, illustrative): T1 = 288.15 K, overall pi = 40, eta_p = 0.9
  constant cp, gamma = 1.4    : T2s =  826.70 K, T2 =  929.41 K, eta_s = 0.8398, w =  644.25 kJ/kg (isentropic  541.06), s_gen = 117.65 J/(kg K), T0 s_gen =  33.90 kJ/kg
  dry air, NASA 9-coefficient : T2s =  801.53 K, T2 =  891.23 K, eta_s = 0.8433, w =  635.02 kJ/kg (isentropic  535.52), s_gen = 117.65 J/(kg K), T0 s_gen =  33.90 kJ/kg
  closed form, gamma 1.4: T1 pi^0.2857 = 826.70 K, eta_s(eta_p) = 0.8398
  first law, NASA air: h2 - h1 = 635.0215 kJ/kg; int cp dT = 635.0215; int v dp = 571.5194; int T ds = 63.5022; v dp + T ds = 635.0215 kJ/kg
```

followed by the stage table, a single stage, a turbine, and the sweep of part 4, and writes `compressor_efficiency.png` (isentropic against polytropic efficiency for compressors and turbines; the stage-count convergence plot).

### Verify

`python -m pytest curriculum/phase-07/lesson-05-isentropic-processes-efficiencies` runs 15 tests:

1. **Inversions**: `temperature_from_s0(gas, gas.s0(T))` and `temperature_from_h(gas, gas.h(T))` return T to 10⁻¹¹ on both sides of the 1000 K interval boundary.
2. **Reduction to the closed form**: the s°-based isentrope given a constant-cp gas (γ = 1.4, 5/3, 1.3; compression and expansion) equals T₁π^((γ−1)/γ) to 10⁻¹¹.
3. **NASA isentrope**: s_gen = 0 on it, it runs cooler than γ = 1.4, and agrees with γ = 1.4 within 0.02 K at π = 1.2 near room temperature.
4. **Second law** (both gases): η = 1 gives s_gen = 0 for a compressor and a turbine; η = 0.95, 0.9, 0.8 give s_gen > 0, increasing as η falls.
5. **Round trips**: isentropic and polytropic efficiencies recovered from computed exit states.
6. **Closed forms**: (7) and (8) for constant cp; η_c < η_p < η_t; both → η_p as π → 1; η_c falls with π.
7. **s_gen = R ln π (1/η_p − 1)** for both gases.
8. **Stages** (both gases): 8 to 64 stages converge to the polytropic exit temperature with observed order 0.85–1.15, and 256 stages recover η_p within 5 × 10⁻⁴; one stage is the whole compressor.
9. **Preheat**: stage works increase through the machine, sum to the overall work, and the overall η_s is below the stage efficiency.
10. **First law** (both gases): h₂ − h₁ = ∫cp dT to 10⁻⁹ and ∫v dp + ∫T ds to 10⁻⁵; ∫T ds/w = 1 − η_p; a reversible path has ∫T ds = 0.
11. **Sanity** and **main** writes the figure.

### Compare

Two independent checks on the property model. Lienhard's Table A.6 (`@lienhard-heat-transfer-textbook`, air at 1 atm) gives cp = 1006, 1051, 1121, 1174 and 1200 J/(kg·K) at 300, 600, 900, 1200 and 1400 K; the NASA dry-air model gives 1004.8, 1051.0, 1120.9, 1174.7 and 1200.4, within 0.12 % everywhere and within 0.06 % above 600 K. And the stage-by-stage compressor, which knows nothing about polytropic efficiency, converges to the polytropic closed form (7) as the stages shrink: the analytic definition and the machine built from stages are the same thing in the limit.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Constant against variable cp.** At η_p = 0.90 from 288.15 K, sweep π from 1.5 to 60 and record the isentropic efficiency (NASA air and the γ = 1.4 closed form), the exit-temperature difference and the error in work of the constant-cp model. *Predict* the sign of the temperature error and the π at which the work error passes 1 %. *Finding:*

   | π | η_s NASA | η_s γ = 1.4 | T₂ NASA − T₂ const (K) | work error, const cp |
   |---|---|---|---|---|
   | 1.5 | 0.8941 | 0.8941 | −0.03 | +0.00 % |
   | 5 | 0.8757 | 0.8755 | −1.73 | +0.11 % |
   | 10 | 0.8648 | 0.8641 | −6.17 | +0.32 % |
   | 20 | 0.8540 | 0.8522 | −16.83 | +0.74 % |
   | 30 | 0.8477 | 0.8450 | −27.69 | +1.12 % |
   | 40 | 0.8433 | 0.8398 | −38.18 | +1.45 % |
   | 60 | 0.8372 | 0.8324 | −57.79 | +2.03 % |

   Constant cp overpredicts the temperature everywhere and the work by more than 1 % above π ≈ 30. The temperature error grows much faster than the work error: a 38 K error in compressor delivery temperature is serious for a combustor or disc designer, while 1.45 % in work is within the uncertainty of the efficiency itself.

2. **Stages to the polytropic limit.** Split π = 40 into n = 1, 2, 4, …, 64 stages of isentropic efficiency 0.90 each (NASA air). *Predict* whether the overall isentropic efficiency rises or falls with n, and the convergence order. *Finding:*

   | n | T_exit (K) | overall η_s | overall η_p | distance to the polytropic limit (K) |
   |---|---|---|---|---|
   | 1 | 855.37 | 0.9000 | 0.9365 | 35.867 |
   | 2 | 869.75 | 0.8764 | 0.9213 | 21.483 |
   | 4 | 879.36 | 0.8613 | 0.9116 | 11.872 |
   | 8 | 884.97 | 0.8527 | 0.9061 | 6.258 |
   | 16 | 888.02 | 0.8481 | 0.9031 | 3.215 |
   | 32 | 889.60 | 0.8457 | 0.9016 | 1.630 |
   | 64 | 890.41 | 0.8445 | 0.9008 | 0.821 |

   The polytropic limit is 891.23 K; the observed orders are 0.92, 0.96, 0.98, 0.99, tending to 1. More stages of the same quality make a *worse* overall machine (preheat): an eight-stage compressor of 0.90 stages is an 0.853 compressor. A single stage at π = 1.35 and η_s = 0.90 has η_p = 0.904, raising the air from 288.15 K to 316.80 K with 28.80 kJ/kg (typical of an axial stage, illustrative).

3. **Compressor against turbine.** Expand NASA air from 1400 K through an expansion ratio of 12 at η_p = 0.90. *Predict* whether η_t is above or below 0.90. *Finding:* T₂ = 801.59 K, η_t = 0.9259 (closed form with γ = 1.4: 0.9288), 691.73 kJ/kg of work out. The same aerodynamic quality is worth 0.926 in a turbine and 0.843 in a compressor.

4. **The first law three ways.** For NASA air at π = 40, η_p = 0.90: h₂ − h₁ = 635.0215 kJ/kg, ∫cp dT = 635.0215 kJ/kg, ∫v dp = 571.5194 and ∫T ds = 63.5022 kJ/kg, summing to 635.0215 kJ/kg. *Predict* ∫T ds/w before running it (0.10, from (9)).

Raw data format: CSV tables `pi, eta_s_nasa, eta_s_cf, dT2_K, dw_percent` (1); `n, T_exit_K, eta_s, eta_p, err_K` (2); `machine, T_in_K, ratio, eta_p, T_out_K, eta_s, w_kJkg` (3); `dh, cpdT, vdp, Tds` (4); plus `compressor_efficiency.png`.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| isentrope, constant cp, π = 40 | (3): 826.70 K | s°-solver with constant-cp gas: 826.70 K | < 10⁻¹¹ relative | identity |
| s_gen at η = 1 | 0 | ≤ 10⁻⁹ J/(kg·K) | round-off | Newton tolerance 10⁻¹⁰ K |
| η_c from η_p = 0.90, π = 40 | (8): 0.8398 | constant-cp solver: 0.8398 | < 10⁻¹⁰ | identity |
| s_gen along the polytropic path | (10): 117.65 J/(kg·K) | NASA 117.65, constant cp 117.65 | < 10⁻⁹ | identity |
| 64 stages against polytropic T₂ | 891.23 K | 890.41 K | 0.82 K | first order in 1/n, as derived |
| first law | h₂ − h₁ | ∫cp dT equal to 10⁻⁹; ∫v dp + ∫T ds equal to 10⁻⁵ | trapezoid error | 2001-point path |
| NASA cp against Lienhard Table A.6 | — | within 0.12 % from 300 to 1400 K | — | fits of different data |

Every result is reproduced by an independent route with the discrepancy the method predicts. The engineering answer: **the core compressor delivers air at 891 K for 635 kJ/kg, an isentropic efficiency of 0.843 for a polytropic 0.90; a constant-cp calculation is 38 K too hot and 1.45 % too expensive, and both agree that 117.65 J/(kg·K) of entropy is generated, which costs 33.9 kJ/kg of exergy at 288 K.** The uncertainty that dominates the work is not the property model but η_p itself: one point of polytropic efficiency (0.89 instead of 0.90) adds 1.87 % to the work and 10.6 K to the delivery temperature, more than the whole variable-cp correction to the work (1.45 %), though less than its 38 K correction to the temperature.

What to improve: for compressor work at π ≤ 10, nothing (constant cp is within 0.3 %). For delivery temperature at high π, or anything hot, use variable cp. For a turbine fed by combustion products, the gas itself changes (`P07.08`).

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2), (3), (7), (8) and (10); `compressor_efficiency.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6, 7) answering the engineering question, with the variable-cp correction and the sensitivity to η_p.

### Adversarial self-checks

- **What if the pressure ratio doubles to 80?** For constant cp the exit temperature rises by the factor 2^(0.2857/0.9) = 1.25, to 1158 K; with NASA air it reaches 1082 K. Delivery temperature, which the last stages, the discs and the combustor inlet must survive, rises quickly with π, and it is one of the limits on the pressure ratio of a core.
- **What if η → 1?** T₂ → T₂ₛ, s_gen → 0, η_s = η_p = 1 and the work falls to the isentropic value (535.52 kJ/kg for NASA air at π = 40). The tests check the limit.
- **What if π → 1?** Isentropic and polytropic efficiencies coincide, and constant and variable cp agree; the distinction matters only for large pressure ratios.
- **Which assumption fails first in the real world?** Constant efficiency: a real compressor's efficiency varies with flow and speed, and at low flow it surges. Next, the adiabatic assumption in small machines.
- **What measurement could falsify the model?** Inlet and delivery temperature and pressure on a running compressor give η_s and η_p by `isentropic_efficiency` and `polytropic_efficiency`; with a torque meter the work too. A measured work that differs from h₂ − h₁ by more than the casing heat loss means the temperature probes are wrong.
- **Does the result survive a change in method?** The stage-by-stage machine and the polytropic integral agree in the limit; the first law closes by three routes.
- **Which decision changes if η_p is uncertain by ±0.02 instead of ±0.01?** At π = 40 the delivery temperature spans 871–913 K instead of 881–902 K, and the work ±3.6 % instead of ±1.8 %: enough to move a combustor inlet or disc temperature margin, and the fuel-flow estimate of `P07.08`.

## Check yourself

1. Write the isentropic relation for an ideal gas with variable cp, and say what s°(T) is.
2. Why is isentropic efficiency defined with enthalpy differences rather than temperature differences?
3. A compressor and a turbine both have polytropic efficiency 0.90 and pressure ratio 12. Which has the higher isentropic efficiency, and why?
4. What is the entropy generated per kilogram in an adiabatic compressor with η = 1? With η = 0.85?
5. At π = 40, which is larger: the work predicted with γ = 1.4, or the work with NASA air? Why?
6. A compressor is made of more and more stages, each of the same isentropic efficiency. Does its overall isentropic efficiency rise or fall, and what does it approach?
7. Why is the friction heat ∫T ds larger than the exergy destroyed T₀ s_gen?
