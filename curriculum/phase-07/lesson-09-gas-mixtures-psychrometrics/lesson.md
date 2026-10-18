# P07.09  Gas Mixtures & Psychrometrics

> Air at 11 km is drier than any desert. Whatever water an airliner cabin holds, the passengers breathed into it, and only the fresh-air flow carries it out again.

| | |
|---|---|
| **Learning outcomes** | 1. Compute the molar mass, mole and mass fractions, partial pressures (Dalton) and partial volumes (Amagat) of an ideal-gas mixture, and its enthalpy and entropy including the entropy of mixing. 2. Define and compute the humidity ratio, relative humidity, dew point and enthalpy of moist air per kilogram of dry air, from the steam table of `P07.03`, and check the saturation pressure against an independent estimate from the NASA polynomials. 3. Write the dry-air, water and energy balances of mixing, cooling-with-condensation and a ventilated cabin, and show that adiabatic mixing lies on a straight line in (ω, h). 4. Implement and test a moist-air module and a steady model of an airliner environmental control system (ECS) from cruise outside air to the cabin. 5. Decide what sets cabin humidity in cruise, and judge which inputs the answer is sensitive to. |
| **Prerequisites** | `P07.01` (NASA polynomials, dry air as a mixture), `P07.02` (control-volume mass and energy balances), `P07.03` (the saturated-water table and interpolation), `P07.04` (entropy generation), `P04.07` (ideal gas and the ISA), `P03.02` (root finding with a bracket), `P03.03` (interpolation). Skill check: you can write a steady-flow energy balance for a box with two inlets and one outlet, and interpolate a table in ln p against 1/T. |
| **Engineering question** | An airliner cruises at 11 km (216.65 K, 22.6 kPa, almost no water). With 186 occupants each breathing out 60 g/h of water, 0.25 kg/min of fresh air per occupant, half the cabin supply recirculated and the cabin held at 22 °C at an 8000 ft cabin altitude (typical, illustrative inputs), what relative humidity and dew point does the cabin settle at, and what does recirculation do to them? The answer: 18.4 % RH, a dew point of about −3 °C (an estimate below the triple point), and nothing: recirculation changes the supply temperature (12.4 °C at 50 %) but not the humidity, which only the fresh-air flow sets. |
| **Assumptions and limits** | Every gas is ideal, including water vapour (good to 0.5 % up to 60 °C at saturation; the table checks it); no enthalpy of mixing; steady flow; the dry-air model of Phase 07 (N₂ 0.7808, O₂ 0.2095, Ar 0.0093, CO₂ 0.0004 by mole, a curriculum convention); liquid water is incompressible. Saturation pressure from the table between 0.01 and 60 °C; below the triple point it is a Clausius–Clapeyron extrapolation over supercooled liquid, an estimate not checked against data (the frost point over ice is a little higher). The occupants' moisture, the fresh-air rate, the cabin heat load, the outside humidity and the bleed temperature are typical, illustrative values, not data. Fails: at high pressure (real-gas mixtures), near the critical point, in fog or cloud (liquid droplets carried with the air), and wherever ice forms. |
| **Analytical work** | Dalton and Amagat for ideal gases; mole and mass fractions; the entropy of mixing −R Σ x ln x; ω = ε p_v/(p − p_v); relative humidity and dew point; moist-air enthalpy per kilogram of dry air; the Gibbs-equality estimate of p_sat; the mixing line and lever rule; the cabin water balance ω_cabin = ω_out + ṁ_w/ṁ_fresh; the ECS energy balances; the entropy generated in the mixing manifold and the pack. |
| **Computational deliverable** | `starter/main.py`: mixtures (`molar_mass`, `mass_fractions`, `mole_fractions`, `partial_pressures`, `partial_volumes`, `mixture_cp`, `mixture_h`, `mixture_s`, `entropy_of_mixing`), water (`p_sat`, `p_sat_gibbs`, `p_sat_extrapolated`), moist air (`humidity_ratio`, `vapour_pressure`, `relative_humidity`, `dew_point`, `dry_air_enthalpy`, `vapour_enthalpy`, `moist_air_enthalpy`), processes (`mix_streams`, `cool_dehumidify`, `cabin_humidity`) and the system (`ecs_cruise`); the NASA functions, `isa_state`, `p_sat_any`, `omega_saturated`, `liquid_enthalpy`, `temperature_from_enthalpy`, `ground_case` and `main` are given. Tests (16) prove Dalton = Amagat, the fraction round trips, the entropy-of-mixing limits, the table against the NASA Gibbs estimate to 1 %, the inverses (ω ↔ p_v, dew point ↔ p_sat), the enthalpy references, mixing on the straight line with all three balances closing, condensation closing water and energy, the cabin's recirculation independence, and the ECS balances closing to round-off. |
| **Practical activity** | **computational experiment**: sweeps of recirculation fraction, fresh-air flow per occupant and outside humidity through the cruise ECS model; a hot, humid ground case with condensation in the pack; and a method comparison of two independent saturation-pressure routes. |
| **Equipment** | None (a computer with the course's Python environment). |
| **Safety** | None beyond normal computer use. |
| **Evidence** | `psychrometrics.png`; the printed check tables; the three sweep tables and the ground case; a half-page note answering the engineering question, naming the input the answer is most sensitive to and the one it is not. |
| **Assessment** | `assessments/phase-07/grading-rubric.md`; competencies 1–4. Pass: ω = ε p_v/(p − p_v) and the cabin water balance derived; all public tests pass; the cabin humidity ratio and RH reproduced by hand to three figures; the note explains why recirculation does not change the humidity and labels the sub-zero dew point as an estimate. |
| **References** | `@openstax-chemistry-2e` 9.3 "Stoichiometry of Gaseous Substances, Mixtures, and Reactions" (The Pressure of a Mixture of Gases: Dalton's Law; Collection of Gases over Water) (required-free; https://openstax.org/books/chemistry-2e/pages/9-3-stoichiometry-of-gaseous-substances-mixtures-and-reactions). `@mit-16-unified-thermodynamics-propulsion-notes` 7.1 "Entropy Change in Mixing of Two Ideal Gases" (required-free, read online; https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node54.html). `@libretexts-yan-engineering-thermodynamics` 7.1 Appendix A (Table A1: Saturated Water), the saturation data used here, generated from NIST SRD 69; 5.4 Applications of the mass and energy conservation equations in steady flow devices (5.3.5 Mixing chambers) (required-free; https://eng.libretexts.org/Bookshelves/Mechanical_Engineering/Introduction_to_Engineering_Thermodynamics_(Yan)). `@lienhard-heat-transfer-textbook` §11.2 Mixture compositions and species fluxes, subsection "Ideal gases"; §11.7 Simultaneous heat and mass transfer at low rates, subsections "The wet bulb temperature" and "Dew-point temperature and dehumidification" (supplementary-free; https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf). `@nasa-tp-2002-211556-glenn-coefficients` "Using the Coefficients to Obtain Thermodynamic Functions" (required-free; https://ntrs.nasa.gov/citations/20020085330): h, s° and g° of H₂O gas and liquid. |
| **Extensions** | ★ Add a wet-bulb temperature function (an adiabatic-saturation energy balance solved with `brentq`) and compare it with the dew point and dry-bulb temperature for the cabin and the ground case. ★ Add ozone and CO₂ to the cabin balance: with 186 occupants each exhaling about 0.3 L/min of CO₂ (illustrative), find the cabin CO₂ mole fraction against fresh-air flow, and see which of humidity and CO₂ limits a reduction in fresh air. ★ Model the air-cycle machine inside the pack (compressor, heat exchanger, turbine, `P07.05`) and find the turbine outlet temperature that this lesson's pack outlet requires. |
| **Time** | unmeasured |

## 1. The engineering problem

Passengers on long flights complain of dry eyes and throats, and they are right: the cabin of an airliner in cruise is usually between 10 % and 20 % relative humidity, drier than most deserts. Why? The air that ventilates the cabin comes from outside, at 11 km, where it is −56.5 °C. Air that cold can hold almost no water: its partial pressure of water vapour is a few pascals. The engines compress it (bleed air) or electric compressors do, the air-conditioning packs cool it back down, and it is mixed with filtered, recirculated cabin air and blown into the cabin. Nobody adds water. The only source of humidity is the people on board, who breathe and perspire, and the only sink is the outflow valve that lets cabin air overboard to make room for fresh air.

This lesson builds the tools to make that quantitative: how ideal gases mix, how much water a kilogram of air holds at a given temperature and pressure, what "relative humidity" and "dew point" mean, and how to write the balances of dry air, water and energy for every box in an environmental control system (ECS). It then answers the cabin question, and a second one that matters on the ground: on a hot, humid day the same pack cools outside air far below its dew point, and the water that condenses has to go somewhere.

Success criteria, written first: Dalton's and Amagat's rules agree to round-off; the saturation pressure from the steam table agrees with an independent estimate from the NASA polynomials within 1 % between 0.01 and 60 °C; adiabatic mixing lands on the straight line between its inlet states in (ω, h); every dry-air, water and energy balance in the ECS closes to round-off; the cabin humidity is reproduced by hand.

## 2. Understand and derive

### 2.1 Ideal-gas mixtures: Dalton and Amagat

A mixture of N species has amounts n_i (mol) and total n = Σ n_i. The **mole fraction** is x_i = n_i/n; the **mass fraction** is y_i = m_i/m = x_i M_i/M, where the mixture's **molar mass** is M = Σ x_i M_i. For the phase's dry air (N₂ 0.7808, O₂ 0.2095, Ar 0.0093, CO₂ 0.0004) M = 28.9657 g/mol and the mass fractions are N₂ 0.7551, O₂ 0.2314, Ar 0.0128, CO₂ 0.0006: oxygen is heavier than the average molecule, so its mass fraction is larger than its mole fraction.

Two ways to split an ideal-gas mixture at T, p and volume V:

- **Dalton**: give each species the whole volume. Alone, it would exert the **partial pressure** p_i = n_i R T/V. Adding them, Σ p_i = n R T/V = p, so p_i = x_i p.
- **Amagat**: give each species the mixture's pressure. Alone, it would occupy the **partial volume** V_i = n_i R T/p. Adding them, Σ V_i = V, so V_i = x_i V.

For ideal gases the two are the same statement, x_i = p_i/p = V_i/V, because ideal-gas molecules ignore each other. (For real gases at high pressure they differ, and neither is exact.) For one mole of air at 300 K and 1 atm the code finds the Dalton pressures sum to 101 325.000 Pa and the Amagat volumes to 24.61721 L, the mixture's volume.

Because the molecules ignore each other, each species also keeps its own enthalpy: **h_mix = Σ x_i h_i(T)** per mole, with no enthalpy of mixing. Entropy is different. Each species in the mixture sits at its partial pressure x_i p, not at p, so

**s_mix(T, p) = Σ x_i [s°_i(T) − R ln(x_i p / p°)]**,

and mixing pure gases at the same T and p raises the entropy by

**Δs_mix = −R Σ x_i ln x_i ≥ 0** per mole of mixture.

It is zero for a pure gas, R ln 2 = 5.763 J/(mol·K) for a 50/50 mixture of two gases, and 4.717 J/(mol·K) (162.8 J/(kg·K)) for making dry air out of its four pure components. Mixing is irreversible: this term is entropy generated, and un-mixing air costs work (oxygen concentrators exist because of it). The 16.Unified notes, section 7.1, derive the two-gas case; the code checks the general one.

### 2.2 Moist air

**Moist air** is dry air plus water vapour. Its pressure is p = p_a + p_v (Dalton), with p_v the partial pressure of the vapour. Engineers count water per kilogram of dry air, because the dry-air flow does not change through a humidifier, a cooler or a cabin while the water flow does. The **humidity ratio** is

ω = m_v/m_a = (n_v M_w)/(n_a M_a) = (M_w/M_a) p_v/p_a,

**ω = ε p_v/(p − p_v)**,  with ε = M_w/M_a = 18.01528/28.9657 = 0.62195.   (1)

Units: kg of water per kg of dry air (often quoted in g/kg).

The vapour pressure cannot exceed the **saturation pressure** p_sat(T): at that point the vapour is in equilibrium with liquid water and any more condenses. The **relative humidity** is

**φ = p_v / p_sat(T)**,   (2)

and the **dew point** T_dp is the temperature at which air cooled at constant p and ω becomes saturated: p_sat(T_dp) = p_v. A surface colder than the dew point gets wet (or frosted, below 0 °C).

Two consequences of (1) and (2) worth seeing once. At a fixed φ and T, ω grows as p falls: 50 % RH at 22 °C is 8.23 g/kg at sea level and 10.91 g/kg at 76.7 kPa (Problem 1), but the dew point is the same, 11.09 °C, because it depends on p_v alone. And at fixed ω, φ falls as T rises: warming the outside air of 11 km to cabin temperature without adding water gives φ = 0.23 % (below).

The **enthalpy of moist air** per kilogram of dry air is the sum of its parts:

**h = h_a(T) + ω h_v(T)**,   (3)

with dry air measured from 0 °C (h_a = 0 at 273.15 K, from the NASA mixture) and water on the steam-table reference (liquid at the triple point = 0), so the vapour carries h_g(0.01 °C) = 2500.92 kJ/kg plus its ideal-gas sensible enthalpy from the NASA H₂O polynomial. At 22 °C the vapour enthalpy is 2541.85 kJ/kg, between the table's saturated-vapour values at 20 and 25 °C (2537.43 and 2546.51): the vapour is an ideal gas to within the table's precision. Dry air's mean cp from 0 to 25 °C is 1004.2 J/(kg·K).

### 2.3 Where the saturation pressure comes from, and an independent check

The saturation pressure here comes from the saturated-water table of `P07.03` (Yan, Appendix A, Table A1, generated from NIST SRD 69): rows every 5 °C from 0.01 to 60 °C, interpolated linearly in ln p against 1/T. That is the right variable pair because of **Clausius–Clapeyron**: for an ideal-gas vapour and a liquid of negligible volume, d ln p_sat/d(1/T) = −h_fg M_w/R, nearly constant.

The independent check uses no steam-table data at all. At saturation the vapour and liquid have equal Gibbs energy per mole, g_v(T, p_sat) = g_L(T). For an ideal-gas vapour g_v = g°_v(T) + RT ln(p_sat/p°), and a liquid's g hardly depends on pressure, so

**p_sat ≈ p° exp[−(g°_v(T) − g°_L(T)) / (R T)]**,   (4)

with g° = h − T s° from the NASA polynomials for H2O (gas) and H2O(L). The code's comparison:

| T (°C) | table (Pa) | NASA Gibbs (Pa) | difference |
|---|---|---|---|
| 0.01 | 610.0 | 612.1 | +0.35 % |
| 20 | 2340.0 | 2338.4 | −0.07 % |
| 40 | 7380.0 | 7369.0 | −0.15 % |
| 60 | 19 950.0 | 19 846.9 | −0.52 % |

The largest difference over 0.01–60 °C is 0.52 %, and both sides explain it: the table prints pressures to 10 Pa (0.3 % at 0 °C), and saturated vapour is not quite ideal (its compressibility factor p v_g/(R_w T) from the same table is 0.995 at 60 °C, Problem 5), which (4) ignores. The latent heats agree too: h_fg at 0.01 °C is 2500.92 kJ/kg from the table and 2500.76 kJ/kg from the NASA enthalpies.

Below the triple point the table has no rows. The cruise cabin's dew point falls there, so the module extrapolates (4)'s cousin, Clausius–Clapeyron with the triple-point h_fg,

ln(p_sat/p_tp) = −(h_fg M_w/R)(1/T − 1/T_tp),

over **supercooled liquid**. That is an estimate: it is not checked against data in this course, and the frost point over ice is a little higher than the dew point over supercooled water. Every sub-zero dew point below is labelled as such.

### 2.4 Mixing two streams

Two moist streams (dry-air flows ṁ₁, ṁ₂) mix adiabatically with no work. Three balances:

dry air  ṁ₃ = ṁ₁ + ṁ₂;  water  ṁ₃ω₃ = ṁ₁ω₁ + ṁ₂ω₂;  energy  ṁ₃h₃ = ṁ₁h₁ + ṁ₂h₂.

So ω₃ and h₃ are the same mass-weighted average of the inlets. Eliminating the weights,

(ω₃ − ω₁)/(ω₂ − ω₃) = (h₃ − h₁)/(h₂ − h₃) = ṁ₂/ṁ₁:   (5)

**the mixed state lies on the straight line joining the inlet states in the (ω, h) plane**, splitting it in inverse proportion to the flows (the lever rule). It is *not* on the straight line in (ω, T), because cp depends a little on ω and on T: mixing equal flows at 2 °C, 4.35 g/kg and 24 °C, 8.00 g/kg gives 13.039 °C, not 13.000 °C. Mixing generates entropy (two temperatures and two humidities become one); this is the irreversibility of every mixing manifold.

### 2.5 Cooling with condensation

A pack's heat exchanger cools air at constant p from T_in to T_out. If T_out is above the inlet dew point, ω is unchanged and Q = ṁ_a(h_in − h_out). If it is below, the air leaves saturated, ω_out = ω_sat(T_out, p), and the rest of the water leaves as liquid (the water separator):

water  ṁ_a ω_in = ṁ_a ω_out + ṁ_cond;  energy  Q = ṁ_a(h_in − h_out) − ṁ_cond h_f(T_out).   (6)

The condensing water gives up its latent heat, so cooling humid air costs far more than its temperature change suggests.

### 2.6 The cabin and its ECS

![Cabin air system: outside air at 11 km (216.65 K, 22.6 kPa, ω_out about 5 × 10⁻⁵) is taken as engine bleed at about 200 °C, cooled in the pack, and the fresh flow ṁ_f joins the recirculated flow ṁ_r in the mixing manifold. The supply ṁ_s = ṁ_f + ṁ_r enters the cabin (22 °C, 75.3 kPa), which adds Q_sens and water ṁ_w; ṁ_f leaves overboard through the outflow valve. A dashed control volume encloses the manifold and cabin.](figures/cabin-ecs.svg)
*Figure 1. The cabin environmental control system. The recirculation loop lies inside the dashed control volume, so it drops out of the water balance (7): only the fresh flow carries the occupants' water overboard.*

Draw the control volume around the cabin **and** the manifold. Fresh air enters (ṁ_f, ω_out), the occupants add ṁ_w of water, and cabin air leaves through the outflow valve (ṁ_f, ω_cabin); the recirculation loop is inside the boundary and cancels. The water balance is

ṁ_f ω_out + ṁ_w = ṁ_f ω_cabin  ⇒  **ω_cabin = ω_out + ṁ_w/ṁ_f**.   (7)

Recirculation does not appear. It cannot: at steady state, water only leaves with the exhaust, and the exhaust flow equals the fresh flow. Recirculation saves energy and bleed air, and it raises the supply temperature for the same cooling, but the humidity is set by the fresh-air flow per kilogram of water produced.

The energy balance on the cabin alone fixes the supply state: ṁ_s h_s + Q_sens + ṁ_w h_v(T_cabin) = ṁ_s h_cabin, with ṁ_s = ṁ_f + ṁ_r the supply flow, Q_sens the net sensible heat gain (people, electronics, sunshine, minus losses through the cold skin) and the occupants' water entering as vapour. The manifold balance then fixes the pack outlet: ṁ_f h_pack + ṁ_r h_cabin = ṁ_s h_s. Add the two and the recirculation cancels again: the pack outlet temperature depends on the fresh flow and the total cabin load, not on recirculation.

**The worked answer.** Inputs (typical, illustrative): 186 occupants, 0.25 kg/min of fresh air each (ṁ_f = 0.775 kg/s), 60 g/h of water each (ṁ_w = 11.16 kg/h = 3.100 × 10⁻³ kg/s), ω_out = 5 × 10⁻⁵, half the supply recirculated, Q_sens = 15 kW, cabin 22 °C at the ISA pressure of 8000 ft (2438.4 m), 75 262 Pa. By (7):

ω_cabin = 5 × 10⁻⁵ + 3.100 × 10⁻³/0.775 = 4.050 × 10⁻³ = **4.050 g/kg**.

By (1) inverted, p_v = ω p/(ε + ω) = 0.004050 × 75 262/0.62600 = **486.9 Pa**. The table gives p_sat(22 °C) = 2645.4 Pa, so **φ = 18.4 %**. p_v is below the triple-point pressure (610 Pa in the table), so the dew point is the Clausius–Clapeyron estimate, 1/T_dp = 1/273.16 − (R/(h_fg M_w)) ln(486.9/610.0) with R/(h_fg M_w) = 1.8454 × 10⁻⁴ K⁻¹: **T_dp ≈ 270.09 K = −3.1 °C**. The outside air itself, warmed to 22 °C at cabin pressure without any water added, would be at φ = 0.23 %.

### 2.7 What generates entropy

Three irreversibilities run in this system: **mixing** in the manifold (2.73 °C dry air into 22 °C cabin air; the code's mixture entropies give 1.81 W/K, of which the dry-air temperature mixing alone is about 1.09 W/K and the rest is the mixing of two humidities); **heat transfer across a finite temperature difference** in the pack (if it simply dumped its 154.7 kW to ram air at 216.65 K, the air's entropy would fall by 422.7 W/K and the ram air's would rise by 714.1 W/K, S_gen = 291.4 W/K; a real air-cycle machine recovers some work in its turbine); and the **throttling** at the outflow valve, where cabin air falls from 75 kPa to 22.6 kPa. None of them changes the humidity answer, which is pure mass conservation.

### 2.8 Where the model stops

- **Below 0 °C the dew point is an estimate.** The table stops at the triple point; the extrapolation is over supercooled liquid; real cold surfaces in the fuselage frost rather than dew.
- **The inputs are illustrative.** The occupants' moisture, the fresh-air rate and the outside humidity are not data from any aircraft. The experiment below shows which ones matter.
- **Steady state.** Boarding, descent and a door opening on a humid ramp are transients the model does not follow.
- **Water stays vapour.** Condensation and frost on the cold inner skin and in the insulation blankets (a real problem: "rain in the plane" on descent) remove water the model keeps in the air, so the real cabin is, if anything, slightly drier.

## 3. Build the computational model

Open `starter/main.py`. The NASA functions (`cp_species`, `h_species`, `s_species`, `g_species`), `isa_state`, `p_sat_any`, `omega_saturated`, `liquid_enthalpy`, `temperature_from_enthalpy`, `ground_case` and `main` are given, with the NASA coefficients for N₂, O₂, Ar, CO₂, H₂O and H₂O(L) and the 0.01–60 °C saturation rows.

### Step 1: mixtures

`molar_mass`, `mass_fractions`, `mole_fractions` are one-liners over dicts. `partial_pressures(n, T, V)` and `partial_volumes(n, T, p)` are Dalton and Amagat, each species alone. `mixture_h` and `mixture_cp` are mole-weighted sums; `mixture_s` evaluates each species at its partial pressure; `entropy_of_mixing` is −R Σ x ln x (skip zero fractions: x ln x → 0).

### Step 2: water

`p_sat(T)` interpolates ln p against 1/T between table rows with `np.interp` (reverse the arrays: 1/T must increase) and raises `ValueError` outside 273.16–333.15 K. `p_sat_gibbs(T)` is (4) with `g_species("H2O", T) - g_species("H2O(L)", T)`. `p_sat_extrapolated(T)` is Clausius–Clapeyron from the triple-point row.

```python
def p_sat(T):
    Tk = SAT_WATER[:, 0] + 273.15; Tk[0] = T_TRIPLE
    return float(np.exp(np.interp(1.0 / T, 1.0 / Tk[::-1], np.log(SAT_WATER[::-1, 1] * 1e6))))
```

### Step 3: moist-air state

`humidity_ratio` and `vapour_pressure` are (1) and its inverse; `relative_humidity` is (2). `dew_point(omega, p)` solves p_sat(T) = p_v with `brentq` when p_v is at or above the triple-point pressure, and otherwise inverts the Clausius–Clapeyron estimate in closed form, returning `(T_dp, True)` to flag it. `dry_air_enthalpy`, `vapour_enthalpy` and `moist_air_enthalpy` are (3) with the references of 2.2.

### Step 4: processes

`mix_streams` is the three balances of 2.4, with T₃ recovered from h₃ by `temperature_from_enthalpy`. `cool_dehumidify` is (6). `cabin_humidity` is (7) plus the manifold's water balance for the supply humidity.

### Step 5: the ECS

`ecs_cruise(...)` evaluates the ISA states, the flows, (7), the cabin energy balance for h_supply, the manifold balance for h_pack, and the pack duty Q_pack = ṁ_f (h(T_bleed) − h_pack); it returns the states and the three residuals (dry air, water, energy) so the balances can be checked rather than trusted.

Running `main()` prints

```
dry air: M = 28.9657 g/mol, eps = M_w/M_a = 0.62195; mass fractions N2 0.7551, O2 0.2314, Ar 0.0128, CO2 0.0006
1 mol of air at 300 K, 1 atm: sum of Dalton partial pressures 101325.000 Pa, sum of Amagat partial volumes 24.61721 L (V = 24.61721 L)
entropy of mixing of air from pure N2, O2, Ar, CO2: 4.7167 J/(mol K) = 162.84 J/(kg K); R ln 2 = 5.7631

saturation pressure: table (NIST SRD 69 via Yan) against the NASA Gibbs-equality estimate
  ...
  largest difference 0.01-60 C: 0.52 %
  h_fg at 0.01 C: table 2500.92 kJ/kg, NASA 2500.76 kJ/kg

ECS at cruise (typical, illustrative inputs): outside 216.65 K, 22632 Pa, w = 5.0e-05 (p_v = 1.82 Pa); cabin 22.0 C at 75262 Pa (8000 ft)
  fresh 0.775 kg/s, recirculated 0.775 kg/s, supply 1.550 kg/s, occupant water 11.16 kg/h
  outside air brought to the cabin unchanged: RH 0.23 %
  cabin w = 4.050 g/kg, RH 18.4 %, dew point -3.1 C (Clausius-Clapeyron estimate below 0.01 C)
  supply w = 2.050 g/kg at 12.40 C; pack outlet 2.73 C; pack duty 154.7 kW from 200 C bleed
  residuals: dry air 0.0e+00 kg/s, water 4.3e-19 kg/s, energy 0.0e+00 W
```

followed by the sweeps of part 4, and writes `psychrometrics.png`: the psychrometric chart at cabin pressure (ω against dry-bulb temperature with lines of constant φ, the sub-zero band shaded as an estimate) with the pack-outlet, supply and cabin states on one straight mixing line; and the percentage difference between the two saturation-pressure routes from 0.01 to 60 °C.

### Verify

Run `python -m pytest curriculum/phase-07/lesson-09-gas-mixtures-psychrometrics`. The tests:

1. **Mixtures**: dry-air M = 28.965 g/mol; fractions round-trip to 10⁻¹²; Dalton and Amagat both give x_i and add to p and V; Δs_mix is 0 for a pure gas, R ln 2 and R ln 3 for equal parts, positive and largest at 50/50; mixture h is mole-weighted, mixture s minus the pure components equals Δs_mix, and cp = dh/dT.
2. **Saturation**: interpolation reproduces every table row; the NASA Gibbs estimate agrees within 1 % over 0.01–60 °C, and its Clausius–Clapeyron slope gives h_fg within 0.5 % of the table; the extrapolation is continuous at the triple point and monotone.
3. **State**: ω ↔ p_v round trips; φ = 1 at saturation; the dew point of saturated air is its temperature, of half-saturated air is where p_sat halves, and below the triple point it is flagged; the enthalpy references hold and the ideal-gas vapour enthalpy is within 0.3 % of the table's h_g at every row.
4. **Processes**: mixing closes dry air, water and energy to 10⁻¹¹, is collinear in (ω, h) and obeys the lever rule; mixing a stream with itself changes nothing; cooling above the dew point condenses nothing, below it leaves saturated air and closes water and energy; recirculation does not change the cabin humidity, and doubling the fresh air halves the rise.
5. **ECS**: all three residuals vanish; the outside state is the ISA's; supply lies between pack and cabin; no duct is below the cabin dew point; with no recirculation the supply is pack air; the pack outlet does not depend on recirculation and falls when the load rises.
6. **main** writes the figure.

### Compare

There is no psychrometric library in the course environment, and none is needed: the model's check is the agreement of two independent property routes (2.3) and the closure of every balance to round-off. If you have a psychrometric chart for sea level, do not compare with it directly: the cabin is at 75 kPa, where the same φ and T mean a larger ω (Problem 1).

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Recirculation.** Same fresh air, recirculation fraction 0, 0.3, 0.5, 0.6. Predict what happens to the cabin RH, the supply temperature and the pack outlet. *Finding:* RH 18.41 % at every fraction; supply 2.73, 8.54, 12.40, 14.33 °C; pack outlet 2.73 °C at every fraction. Recirculation changes how cold the supply air is, which matters for comfort and for duct condensation, but neither humidity nor the pack's outlet temperature.
2. **Fresh air per occupant.** 0.15, 0.25, 0.35, 0.50 kg/min each, recirculation 0.5. Predict how RH scales. *Finding:* RH 30.40, 18.41, 13.24, 9.35 % and dew points 3.88, −3.06 (estimate), −7.43 (estimate), −11.88 °C (estimate). The rise in ω is inversely proportional to the fresh flow (7), so RH roughly halves when the fresh air doubles. This is the trade an ECS designer faces: less fresh air means a more comfortable humidity and less bleed (fuel), but more CO₂ and odour (Extension 2).
3. **Outside humidity.** ω_out from 0 to 3 × 10⁻⁴. Predict how much the answer depends on an input we do not know well. *Finding:* RH 18.18, 18.23, 18.41, 18.63, 19.53 % for ω_out = 0, 10⁻⁵, 5 × 10⁻⁵, 10⁻⁴, 3 × 10⁻⁴. Even six times the assumed value changes RH by about one point: the occupants dominate, and the uncertainty in the outside humidity at cruise does not matter.
4. **Hot, humid ground day.** Outside 35 °C, 60 % RH, sea level, the same fresh flow, pack outlet 2 °C (illustrative). Predict whether water condenses and how much. *Finding:* the outside air holds 21.45 g/kg with a dew point of 26.1 °C; leaving the pack saturated at 2 °C it holds 4.351 g/kg, so **47.7 kg/h** of water condenses and must be separated (otherwise it is blown into the cabin as fog), and the pack removes 59.8 kW, much of it latent. On the ground the ECS dries the air; in cruise it has nothing to dry.
5. **Two routes to p_sat.** Plot the table against (4) from 0.01 to 60 °C (the figure's right panel). Predict the sign of the difference at 60 °C. *Finding:* within 0.52 %, with the NASA estimate low at the warm end, where the real saturated vapour is least ideal; the scatter at 0–15 °C is the table's 10 Pa print resolution.

Raw data format: CSV tables `recirc_fraction, phi_cabin_pct, T_supply_C, T_pack_C` (1); `fresh_kg_min_per_occ, phi_cabin_pct, T_dew_C, dew_is_estimate` (2); `w_out, phi_cabin_pct` (3); the ground-case line (4); `T_C, p_table_Pa, p_nasa_Pa` (5); plus `psychrometrics.png`.

## 5. Validate and improve

| quantity | analysis (§2) | code | difference | cause |
|---|---|---|---|---|
| cabin ω | (7): 4.050 g/kg | 4.050 g/kg | round-off | — |
| cabin RH | by hand 486.9/2645.4 = 18.4 % | 18.41 % | — | — |
| p_sat, 0.01–60 °C | NASA Gibbs (4) | table | ≤ 0.52 % | table resolution; vapour non-ideality |
| h_fg at 0.01 °C | NASA 2500.76 kJ/kg | table 2500.92 kJ/kg | 0.006 % | — |
| mixing state | lever rule (5) | `mix_streams` | collinear to 10⁻⁹ | — |
| ECS balances | dry air, water, energy | residuals 0, 4 × 10⁻¹⁹ kg/s, 0 W | round-off | — |

The engineering answer: **with 0.25 kg/min of fresh air and 60 g/h of water per occupant, the cruise cabin settles at 4.05 g/kg, 18.4 % RH at 22 °C, with a dew point near −3 °C (an extrapolated estimate); recirculation leaves this unchanged.** The answer is most sensitive to the fresh-air flow per occupant and the moisture per occupant (it is proportional to their ratio), and insensitive to the outside humidity. The sub-zero dew point carries the extrapolation's uncertainty; a verified saturation table over ice would replace it.

What to improve: the occupants' moisture rate (it is an illustrative number; measure or cite one before designing), and a transient model for boarding and descent. What not to improve: the gas model. Ideal mixing is good to well under 1 % at these pressures.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1), (5) and (7); `psychrometrics.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question, saying which input dominates and which does not, and labelling the dew point's status.

### Adversarial self-checks

- **What if the number of passengers doubles with the same total fresh air?** The humidity rise doubles: 8.05 g/kg, 36.4 % RH (dew point 6.4 °C, now from the table). Certification sets fresh air per occupant, not per aircraft, which is why the rise per occupant is what matters.
- **What if the fresh air goes to zero?** (7) says ω_cabin grows without bound: the cabin approaches saturation and the windows fog. A closed spacecraft cabin has no fresh air at all, and must condense the water out instead (Problem 3).
- **Which assumption fails first?** The steady state: on descent cold-soaked structure warms through the dew point of the cabin air and frost in the insulation melts. Then the sub-zero dew point estimate.
- **What measurement could falsify the model?** A cabin hygrometer at cruise with a known passenger count and pack flow: if the measured ω_cabin − ω_out is not close to ṁ_w/ṁ_f, either the moisture per occupant is wrong or water is being lost to cold surfaces.
- **Does the result survive a change of method?** The humidity is algebra (7); the saturation pressure has two independent routes that agree to 0.52 %; the enthalpies are checked against the table's h_g to 0.3 %.
- **Which decision changes if the uncertainty doubles?** Halving or doubling the per-occupant moisture moves RH between 9.3 % and 36.4 %: a decision to fit a cabin humidifier would hinge on that number, so it needs a source, not an assumption.

## Check yourself

1. Why do Dalton's and Amagat's rules give the same mole fractions for ideal gases but not for real gases?
2. Write the humidity ratio in terms of p_v and p, and say why engineers count water per kilogram of *dry* air.
3. Air at 22 °C and 50 % RH is taken from sea level to a cabin at 75 kPa without adding water. Does its relative humidity change? Its dew point?
4. What two quantities does the entropy of mixing −R Σ x ln x depend on, and why is it zero for a "mixture" of a gas with itself?
5. Why does the mixed state of two moist streams lie on a straight line in (ω, h) but not exactly in (ω, T)?
6. Why does recirculating cabin air not change the cabin's humidity at steady state?
7. A pack cools 35 °C, 60 % RH air to 2 °C. Which balance tells you how much water condenses, and why does the heat removed exceed ṁ cp ΔT?
8. Why is the cruise cabin dew point in this lesson labelled an estimate, and in which direction would ice change it?
