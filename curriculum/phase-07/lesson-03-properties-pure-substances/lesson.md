# P07.03  Properties of Pure Substances

> Near the dome a substance is not an ideal gas, and its properties come from a table. Read the table with the right interpolation, and check it against the laws it must obey before you trust a number from it.

| | |
|---|---|
| **Learning outcomes** | 1. Describe the saturation dome of a pure substance (compressed liquid, saturated mixture, superheated vapour, critical point) and identify the phase from (T, p) or (T, v). 2. Compute quality and two-phase mixture properties, and look up superheated vapour by interpolation in temperature and pressure. 3. Derive the Clausius–Clapeyron equation and use it, the identity h = u + pv and the compressibility factor Z to test a property table. 4. Implement and test a table lookup for water and ammonia, including the complete state from (T, v) and a rigid-tank first-law problem. 5. Decide how much ullage a sealed ammonia propellant tank must keep so that it never becomes liquid-full. |
| **Prerequisites** | `P03.03` (linear interpolation, interpolation error, choosing the interpolating variable), `P07.01` (properties, state postulate, ideal gas, NASA fits), `P07.02` (first law for a closed system), `P03.02` (bracketed root finding). Skill check: you can interpolate linearly between two table rows by hand and explain why interpolating ln p against 1/T might beat interpolating p against T. |
| **Engineering question** | The X-15's XLR99 engine burned liquid oxygen with anhydrous ammonia. Suppose an ammonia tank is filled and sealed at 15 °C and may soak to 40 °C on a hot ramp (typical, illustrative temperatures). What pressure keeps the ammonia liquid at each temperature, and what is the largest liquid fill that never lets the tank become liquid-full? The answers: 728.5 kPa at 15 °C and 1555.4 kPa at 40 °C; the liquid fill must not exceed 93.74 % of the tank volume (91.03 % if the soak can reach 50 °C). A 90 % fill rises to 95.96 % liquid at 40 °C and absorbs 66.1 MJ per cubic metre of tank on the way. |
| **Assumptions and limits** | The tables are equilibrium properties (no superheating of liquid or supercooling of vapour); the lookup interpolates linearly in T (saturation, superheated) and between isobars with p v linear in p and s linear in ln p; p_sat is interpolated as ln p linear in 1/T. Only the subset of Yan's tables embedded in the code is available (water superheated at 10, 100, 1000 and 4000 kPa; ammonia at 100, 500, 1000 and 2000 kPa up to 200 °C); there is no compressed-liquid lookup. Interpolation error grows near the critical point and where the table is coarse; the rigid tank assumes a perfectly rigid wall and uniform temperature. |
| **Analytical work** | Phase identification; the lever rule for quality; the Clausius–Clapeyron equation from the equality of the Gibbs function across the dome, and its ideal-gas form d ln p/d(1/T) = −h_fg/R; the round-off bound of h = u + pv; Z and its low-pressure limit; the rigid-tank energy balance Q = m(u₂ − u₁); the liquid-full fill limit v = v_f(T_max). |
| **Computational deliverable** | `starter/main.py`: `saturation_pressure` (linear or logarithmic), `saturation_temperature`, `saturation_state`, `quality`, `mixture_property`, `phase_from_Tp`, `phase_from_Tv`, `superheated_state`, `state_from_Tv`, `rigid_tank_heating`, `max_liquid_fill`, and the table checks `clausius_clapeyron_hfg`, `enthalpy_identity_residual`, `compressibility_factor`, `superheated_dh_vs_ideal_gas`. The embedded tables, `sat_table`, `critical_point`, `gas_constant`, `interp`, `nasa_h2o_enthalpy`, `leave_one_out_psat` and `main` are given. Tests (23) prove node reproduction, log-over-linear interpolation, inversion of p_sat, phase and quality consistency, the superheated lookup, (T, v) round trips including supercritical states, Clausius–Clapeyron to 0.2 %, h = u + pv within the round-off bound for all 376 entries, Z → 1 with 1 − Z ∝ p, the NASA ideal-gas enthalpy at 10 kPa to 0.5 %, and the rigid-tank first law and fill limit. |
| **Practical activity** | **data-based investigation**: the NIST SRD 69 property values for water and ammonia, as published in Yan's appendices, are tested against Clausius–Clapeyron, h = u + pv, the ideal-gas limit and the NASA ideal-gas fit, and two interpolation methods are compared by leaving rows out. |
| **Equipment** | None (a computer with the course's Python environment). |
| **Safety** | None beyond normal computer use. (Ammonia is toxic and corrosive, and liquid ammonia is stored under pressure; it is only computed here, never handled.) |
| **Evidence** | `property_tables.png`; the investigation tables (Clausius–Clapeyron errors, identity residuals against their bound, Z against p, leave-one-out errors, NASA comparison); the tank calculation; a half-page note answering the engineering question and stating how far you trust the table and why. |
| **Assessment** | `assessments/phase-07/grading-rubric.md`; competencies 1–4. Pass: Clausius–Clapeyron derived; all public tests pass; a quality and a superheated lookup done by hand to four figures; the note gives the fill limit with the tank pressures and explains one table check that would have caught a corrupted row. |
| **References** | `@libretexts-yan-engineering-thermodynamics` Ch. 2 "Thermodynamic Properties of a Pure Substance": 2.2 Pure substance; 2.3 Thermodynamic properties; 2.4 Phase diagrams; 2.5 Thermodynamic tables; 3.3 Real gas and compressibility factor (required-free; https://eng.libretexts.org/Bookshelves/Mechanical_Engineering/Introduction_to_Engineering_Thermodynamics_(Yan)); Ch. 7 "Appendices": 7.1 Appendix A - Thermodynamic Properties of Water (Table A1: Saturated Water; Table A2: Superheated Vapor, Water) and 7.2 Appendix B - Thermodynamic Properties of Ammonia (Table B1: Saturated Ammonia; Table B2: Superheated Ammonia), the dataset, whose values Yan generated from the NIST Chemistry WebBook, SRD 69 (Lemmon, McLinden and Friend). `@mit-16-unified-thermodynamics-propulsion-notes` 8.1 Behavior of Two-Phase Systems; 8.4 The Clausius-Clapeyron Equation (required-free, read online, https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node64.html). `@openstax-university-physics-vol2` 1.5 Phase Changes (required-free). `@libretexts-ellgen-thermodynamics-chemical-equilibrium` Ch. 12 "Applications of the Thermodynamic Criteria for Change": 12.8 The Clapeyron Equation; 12.9 The Clausius-Clapeyron Equation (supplementary-free). |
| **Extensions** | ★ Add the CO₂ tables (Yan Appendix D, `co2_saturated_T.csv`, `co2_superheated.csv`) and repeat the investigation for CO₂, the working fluid of supercritical-CO₂ power cycles. ★ Replace linear interpolation in T by a cubic spline (`P03.03`) and repeat the leave-one-out test for v_g, h_g and p_sat; where does the spline win and where does it overshoot? ★ Add the compressed-liquid table (Yan Table A3) and measure how good the common approximation "compressed liquid ≈ saturated liquid at the same T" is for v, u and h. |
| **Time** | unmeasured |

## 1. The engineering problem

The North American X-15 flew to the edge of space on a Reaction Motors XLR99 engine that burned liquid oxygen with anhydrous ammonia. Ammonia is an unusual fuel, and an instructive one for a thermodynamicist: at room temperature and atmospheric pressure it is a gas, so the vehicle keeps it liquid by keeping it under pressure. The same is true of ammonia's other aerospace job, as the working fluid of large spacecraft thermal-control loops, and of water in a spacecraft's evaporative cooler, which works the other way: it boils water at a fraction of a kilopascal to throw heat away. In all of these the fluid sits near its saturation dome, where the ideal-gas model of `P07.01` is not good enough and properties come from tables.

A sealed tank of liquid ammonia raises a classic design question. The tank holds liquid plus a little vapour (the ullage). As it warms, the liquid expands; the vapour space shrinks; and if the liquid ever fills the tank, the pressure is no longer set by the saturation curve but by the stiffness of the liquid and the tank, and it rises very steeply with temperature. Tanks are therefore filled to a maximum liquid fraction set by the highest temperature they may see. Computing that limit needs v_f(T) from a table and nothing else; computing the pressure and the heat absorbed needs the rest of the table.

The tables in this lesson are the NIST Chemistry WebBook values (Standard Reference Database 69, by Lemmon, McLinden and Friend) as published in the appendices of Yan's *Introduction to Engineering Thermodynamics*. Tables are data, and data can be wrong (a typing error, a rounding, a wrong row). So the lesson does not just read the table: it **tests** it against three laws it must obey and against an independent fit. Success criteria, written first: interpolation reproduces every node exactly; Clausius–Clapeyron, computed from the table's own p_sat(T), reproduces its h_fg to 0.2 % wherever p_sat is printed to four figures; h = u + pv holds within the bound set by the table's rounding for every entry; Z → 1 as p → 0; and superheated-steam enthalpy differences at 10 kPa agree with the NASA ideal-gas fit to 0.5 %.

## 2. Understand and derive

### 2.1 The phases of a pure substance

A **pure substance** has one chemical composition throughout (it may be in more than one phase). Heat water at constant pressure from cold: it is a **compressed (subcooled) liquid** while its temperature is below the saturation temperature for that pressure; at T_sat(p) it begins to boil and becomes a **saturated mixture** of liquid and vapour at constant temperature; when the last liquid evaporates it is **saturated vapour**; heated further it is **superheated vapour**. On a T–v diagram the saturated-liquid states (v_f) and the saturated-vapour states (v_g) form a dome that closes at the **critical point** (373.95 °C, 22 064 kPa for water; 132.25 °C, 11 339 kPa for ammonia), above which there is no distinction between liquid and vapour.

![T–v diagram with v on a log scale: the saturated-liquid line v_f and the saturated-vapour line v_g form a dome that closes at the critical point. Compressed liquid lies to the left of the dome, saturated mixture inside it and superheated vapour to the right; an isobar crosses the dome horizontally at T_sat(p).](figures/t-v-saturation-dome.svg)
*Figure 1. Phases of a pure substance on a T–v diagram. Inside the dome p and T are tied together by saturation, so an isobar is flat there; v_f and v_g meet at the critical point, above which liquid and vapour cannot be told apart.*

Inside the dome **p and T are not independent**: at saturation, fixing T fixes p = p_sat(T). That is the state postulate of `P07.01` at work: (T, p) does not fix a two-phase state; (T, v), (T, h) or (p, x) does.

### 2.2 Phase identification and quality

From **(T, p)**: if T is above the critical temperature, the state is supercritical; otherwise compare p with p_sat(T). p > p_sat: compressed liquid. p < p_sat: superheated vapour. p = p_sat: saturated, and you need one more property.

From **(T, v)**: compare v with v_f(T) and v_g(T). Inside the dome the **quality** x is the mass fraction of vapour, and any specific property y (v, u, h, s) is the mass-weighted average of the two phases:

**y = y_f + x (y_g − y_f) = y_f + x y_fg,  x = (y − y_f)/y_fg.**

This is the lever rule. A worked example: water at 175 °C with v = 0.1 m³/kg. The table gives v_f = 0.001121, v_g = 0.216581 m³/kg, so v_f < v < v_g: a saturated mixture with **x = (0.1 − 0.001121)/0.215460 = 0.4589**, at p = p_sat = **892.60 kPa**, with h = h_f + x h_fg = **1673.40 kJ/kg**.

### 2.3 Interpolation, and which variable to interpolate

Tables give values at discrete rows; between them we interpolate (`P03.03`). Linear interpolation has an error proportional to the square of the row spacing times the second derivative of the function. Its quality therefore depends on the variable used: a function that is nearly linear in some transformed variable interpolates far better in that variable.

p_sat(T) is strongly convex (for water it rises by a factor of about 1.8 every 10 K near room temperature: 2.34 kPa at 20 °C, 4.25 kPa at 30 °C), so a straight chord in (T, p) lies above the curve and overestimates it. The Clausius–Clapeyron equation (2.4) predicts that ln p is nearly linear in 1/T. The investigation (Section 4) measures both choices.

For superheated vapour the code interpolates linearly in T along each tabulated isobar, then between the two isobars that bracket p, using the variables that are linear in an ideal gas: **p v linear in p** (p v = RT exactly for an ideal gas), **s linear in ln p** (s = s°(T) − R ln p), u and h linear in p (constant for an ideal gas). Water at 500 kPa and 300 °C, interpolated from the 100 and 1000 kPa isobars, gives v = 0.52253 m³/kg, h = 3064.35 kJ/kg and s = 7.4535 kJ/(kg·K); Yan's Table A2 (not embedded in the code) lists 0.52261, 3064.55 and 7.4614. The errors (−0.015 %, −0.007 %, −0.11 %) come from interpolating across a factor of 10 in pressure; s is worst because ln p interpolation of s is exact only for an ideal gas.

### 2.4 The Clausius–Clapeyron equation

Along the saturation curve the two phases are in equilibrium, which (lesson 04 and `@mit-16-unified-thermodynamics-propulsion-notes` 8.4 develop this) means their specific Gibbs functions g = h − Ts are equal: g_f = g_g. Move along the curve by dT and dp. For each phase dg = v dp − s dT (from T ds = dh − v dp), and the two changes must be equal to keep g_f = g_g:

v_f dp − s_f dT = v_g dp − s_g dT  ⇒  dp/dT = (s_g − s_f)/(v_g − v_f) = s_fg/v_fg.

At constant T and p, the heat to evaporate a kilogram reversibly is h_fg = T s_fg. Hence

**dp_sat/dT = h_fg/(T v_fg).**

This is a law the table must obey: p_sat(T), v_fg and h_fg are tabulated independently, and Clausius–Clapeyron ties them together. **If** the vapour is ideal (v_g = RT/p) and v_f ≪ v_g, the equation becomes d ln p/d(1/T) = −h_fg/R, so ln p is linear in 1/T wherever h_fg is nearly constant: the reason for logarithmic interpolation. That ideal-gas form is only approximate: at 100 °C, the slope of ln p against 1/T across the 95–105 °C rows gives h_fg ≈ −R × slope = 2293.2 kJ/kg against the tabulated 2256.4 (1.6 % high), because v_f is not zero and the vapour is not ideal. The full form, used by the code, has no such approximation.

The code's check: at each interior row, form dp/dT = −(p/T²) d ln p/d(1/T) from central differences of ln p against 1/T across the neighbouring rows, multiply by T v_fg, and compare with the tabulated h_fg. At 100 °C it gives 2256.29 kJ/kg against 2256.40 (−0.005 %).

### 2.5 The identity h = u + pv, and table round-off

Enthalpy is defined as h = u + pv, and Yan's tables list h, u, p and v separately. They must agree, but only to the precision printed. With p in kPa and v in m³/kg, pv is in kJ/kg directly. The bound on |h − u − pv| from rounding alone is half a unit in the last digit of each: 0.005 + 0.005 kJ/kg for h and u, plus δp·v plus p·δv. For saturated water vapour at 0.01 °C (p printed as 0.00061 MPa, so δp = 0.005 kPa, and v_g = 205.991 m³/kg) the bound is about 1 kJ/kg; for most entries it is 0.01–0.02 kJ/kg. A residual larger than its bound means a corrupted entry. All 222 water entries and 154 ammonia entries pass, the largest water residual being 0.380 kJ/kg at 95 % of its bound, the largest ammonia one 0.009 kJ/kg.

### 2.6 The compressibility factor and the ideal-gas limit

**Z = pv/(RT)** measures how far a gas is from ideal (Z = 1). R = R_U/M is 0.461 523 kJ/(kg·K) for water and 0.488 210 kJ/(kg·K) for ammonia. Molecules attract each other, which lowers the pressure for a given v and T, so Z < 1 for vapours near their dome. At low pressure the molecules are far apart and Z → 1; to first order 1 − Z is proportional to p (the second virial coefficient). The table shows it: steam at 300 °C has Z = 0.9998, 0.9976, 0.9753 and 0.8902 at 10, 100, 1000 and 4000 kPa; ammonia at 60 °C has 0.9933, 0.9655, 0.9286 and 0.8454 at 100, 500, 1000 and 2000 kPa. At low pressure, then, the table's superheated steam should behave like the NASA ideal-gas H₂O of `P07.01`, and Section 4 checks it.

### 2.7 The sealed tank

A sealed rigid tank is a closed system with fixed v = V/m (`P07.02`). As it warms from T₁ to T₂, both states are fixed by (T, v). No work is done (rigid), so the first law gives **Q = m(u₂ − u₁)**. The tank stays two-phase as long as v > v_f(T): the liquid volume fraction (1 − x) v_f/v grows as the liquid expands. It becomes **liquid-full** at the temperature where v_f(T) = v. To stay two-phase up to T_max, the tank's v must be at least v_f(T_max). At the filling temperature T_fill this sets the largest liquid volume fraction:

**x_fill = (v_f(T_max) − v_f(T_fill))/v_fg(T_fill),  fill_max = (1 − x_fill) v_f(T_fill)/v_f(T_max).**

For ammonia, v_f = 0.001619 m³/kg at 15 °C and 0.001726 m³/kg at 40 °C, so **fill_max = 93.74 %**. Because v_g ≫ v_f, this is almost exactly v_f(15)/v_f(40) = 0.9380: the ullage must absorb the liquid's 6.6 % thermal expansion. A worked tank: 1 m³ filled 90 % liquid by volume at 15 °C holds 556.47 kg, of which only 0.103 % by mass (x₁ = 0.00103) is vapour, though it fills 10 % of the volume. At 40 °C it is at 1555.4 kPa, x₂ = 0.00087, the liquid fills 95.96 % of the volume, and the first law says it absorbed Q = **66.1 MJ**.

### 2.8 Where the model fails

- **Near the critical point.** v_fg → 0, h_fg → 0, and every property varies steeply; linear interpolation and the Clausius–Clapeyron check both degrade (the ammonia check is −0.5 % at 130 °C, against 0.02 % elsewhere).
- **Coarse rows and few digits.** Water's p_sat below 20 °C is printed to two or three significant figures, which limits the Clausius–Clapeyron check there to about 1 % (Section 4).
- **Between widely spaced isobars.** Interpolating from 100 to 1000 kPa loses 0.1 % in s (2.3); a finer table would remove it.
- **Liquid-full tanks.** Once v < v_f the state is compressed liquid; this lookup has no compressed-liquid table and refuses the state. In reality the pressure then rises by megapascals per kelvin: that is the failure the fill limit prevents.
- **Metastable states.** Rapid depressurisation can superheat a liquid briefly; tables give equilibrium only.

## 3. Build the computational model

Open `starter/main.py`. The embedded tables (Yan Tables A1, A2, B1 and B2, values from NIST SRD 69; each as given constants with its source in a comment), the NASA H₂O coefficients, `sat_table` (named columns, p in kPa), `critical_point`, `gas_constant`, `interp` (linear, refusing to extrapolate), `nasa_h2o_enthalpy`, `leave_one_out_psat` and `main` are given.

### Step 1: saturation lookup

`saturation_pressure(fluid, T_C, method)`: for `"linear"` call `interp`; for `"log"` find the bracketing rows and interpolate ln p linearly in 1/T (in kelvin). `saturation_state` interpolates everything else linearly in T and adds the "fg" differences. `saturation_temperature(fluid, p)` inverts the log interpolation with `brentq` on the table's temperature range.

### Step 2: phase and quality

`quality` and `mixture_property` are the lever rule; `phase_from_Tp` compares p with p_sat (and T with T_c); `phase_from_Tv` compares v with v_f and v_g and returns the quality inside the dome.

### Step 3: superheated vapour and the (T, v) state

`superheated_state(fluid, p, T)` returns the tabulated isobar directly if p is one, and otherwise interpolates between the two bracketing isobars as in 2.3:

```python
f = (p - p0) / (p1 - p0)
pv = (1 - f) * p0 * a[0] + f * p1 * b[0]          # p v linear in p
g = np.log(p / p0) / np.log(p1 / p0)
s = (1 - g) * a[3] + g * b[3]                     # s linear in ln p
```

`state_from_Tv` uses the quality inside the dome and root-finds p on v(p, T) outside it.

### Step 4: the tank and the table checks

`rigid_tank_heating` and `max_liquid_fill` are 2.7. `clausius_clapeyron_hfg`, `enthalpy_identity_residual`, `compressibility_factor` and `superheated_dh_vs_ideal_gas` are 2.4–2.6.

Running `main()` prints the five investigation sections and the engineering case (abridged):

```
1. Clausius-Clapeyron: h_fg = T v_fg dp/dT from the table's own p_sat(T) against tabulated h_fg
  water: 74 interior rows; median |error| 0.005 %, largest |error| below Tc - 30 C 0.96 % (at 10 C), nearest the critical point -0.1 % (at 370.00 C)
  ammonia: 36 interior rows; median |error| 0.009 %, largest |error| below Tc - 30 C 0.02 % (at -35 C), nearest the critical point -0.5 % (at 130.00 C)
2. h = u + p v, residual against what table round-off allows
  water: 222 entries; largest |residual| 0.380 kJ/kg; largest |residual|/bound 0.95; entries exceeding the bound: 0
  ammonia: 154 entries; largest |residual| 0.009 kJ/kg; largest |residual|/bound 0.83; entries exceeding the bound: 0
4. Interpolating p_sat between rows: leave every other row out
  water, linear: largest |error| 5.711 % (at 5 C), median 0.495 %
  water,    log: largest |error| 0.233 % (at 15 C), median 0.015 %
Engineering case (typical, illustrative): a sealed ammonia propellant tank filled at 15 C
  p_sat(15 C) = 728.5 kPa: the tank must hold at least this to keep ammonia liquid; p_sat(40 C) = 1555.4 kPa
  largest liquid fill at 15 C that stays two-phase up to 40 C: 93.74 %; up to 50 C: 91.03 %
```

and writes `property_tables.png`: the two saturation domes on T–v axes; the Clausius–Clapeyron error against T; and the leave-one-out p_sat errors, linear against logarithmic.

### Verify

`python -m pytest curriculum/phase-07/lesson-03-properties-pure-substances` runs 23 tests: saturation lookups reproduce the table nodes to 10⁻¹² and refuse T above T_c; logarithmic p_sat interpolation is at least five times better than linear and the linear chord lies above the convex curve; p_sat inverts; phase from (T, p) and from (T, v), quality from v, u, h and s, and mixture properties are mutually consistent; superheated nodes are reproduced and interpolated values lie between neighbouring isobars; (T, v) states round-trip in the superheated, two-phase and supercritical regions; Clausius–Clapeyron within 0.2 % wherever p_sat ≥ 10 kPa and T < T_c − 30 °C, with median error below 0.05 %; h = u + pv within its round-off bound for every entry; Z < 1, decreasing with p, within 1 % of 1 at the lowest isobar, with 1 − Z roughly proportional to p; the NASA comparison at 10 kPa within 0.5 % and more than 5 % off at 4 MPa; the tank's first law, pressure and fill limit.

### Compare

The established "tool" for this lesson is the table itself, and the check is against laws and an independent model rather than another table. Clausius–Clapeyron, a consequence of the second law, holds to a median 0.005 % (water) and 0.009 % (ammonia) on data whose h_fg, v_fg and p_sat were tabulated separately. The NASA ideal-gas H₂O fit, built from spectroscopy and statistical mechanics with no reference to steam tables, agrees with the table's superheated steam at 10 kPa to 0.05–0.25 %. Two independent routes, consistent to better than a quarter of a per cent: the table and the code are right.

## 4. Experiment at home, or investigate real data

**Mode: data-based investigation.** **Dataset:** saturated (temperature) and superheated tables for water (Yan Appendix A, Tables A1 and A2) and ammonia (Yan Appendix B, Tables B1 and B2), shared as CSV for this phase and embedded in the code as given constants. **Source and licence:** Yan states that every value was generated from the NIST Chemistry WebBook, SRD 69 (Lemmon, McLinden and Friend, "Thermophysical Properties of Fluid Systems"); the appendix pages on LibreTexts carry no declared licence, so the values are used as data with a citation. **What was measured and how:** the numbers are not raw measurements but values computed by NIST from reference equations of state fitted to many experiments; the tables round them to 4–6 significant figures. Write each prediction before running it.

1. **Clausius–Clapeyron across the dome.** For every interior row of both saturation tables compute h_fg from T v_fg dp/dT and compare with the tabulated h_fg. *Predict* where the agreement is worst. *Finding:* median errors 0.005 % (water) and 0.009 % (ammonia). Water's worst is 0.96 % at 10 °C, where p_sat is printed as 0.00123 MPa: a half-unit rounding of the last digit is 0.4 % of p, and differencing two such values amplifies it. Ammonia, printed to five significant figures throughout, stays within 0.02 % below 102 °C. Near the critical point the errors grow (−0.1 % for water at 370 °C, −0.5 % for ammonia at 130 °C) as v_fg → 0.
2. **The identity h = u + pv.** Compute the residual for all 376 entries and its round-off bound. *Predict* whether any entry fails. *Finding:* none; the largest residual uses 95 % (water) and 83 % (ammonia) of its bound. A deliberately corrupted entry (change one h by 1 kJ/kg) fails at once: try it.
3. **The ideal-gas limit.** Tabulate Z on every superheated isobar at 300 and 500 °C (water) and 60 and 200 °C (ammonia). *Predict* the pressure below which steam at 300 °C is ideal to 1 %. *Finding:* the Z values of 2.6; with 1 − Z ∝ p, Z = 0.99 at about 405–415 kPa for steam at 300 °C and 145–149 kPa for ammonia at 60 °C (estimated from the 100 and 1000 kPa, and 100 and 500 kPa, isobars respectively).
4. **Interpolation method.** Leave every other row of each saturation table out, interpolate p_sat at the missing rows from the rest, linearly and logarithmically. *Predict* the ratio of the errors. *Finding:* linear: largest 5.711 % (water) and 3.257 % (ammonia), median 0.495 % and 0.729 %; logarithmic: largest 0.233 % and 0.084 %, median 0.015 % and 0.026 %. The logarithmic method is about 25 to 40 times better, as Clausius–Clapeyron predicts.
5. **Against the NASA fit.** Compare h(T₂) − h(T₁) of superheated steam at 10 kPa with the NASA ideal-gas H₂O fit. *Predict* the size of the difference. *Finding:* +0.252 % (100 → 300 °C), +0.106 % (100 → 600 °C), −0.050 % (300 → 1000 °C) at 10 kPa; at 4000 kPa, 300 → 600 °C, the table is 13.30 % above the ideal gas, because the 300 °C state at 4 MPa (Z = 0.89) has a large negative enthalpy departure.

Raw data format: CSV tables `fluid, T_C, hfg_cc, hfg_table, err` (1); `fluid, entry, residual_kJkg, bound_kJkg` (2); `fluid, T_C, p_kPa, Z` (3); `fluid, method, T_C, rel_err` (4); `p_kPa, T1_C, T2_C, dh_table, dh_nasa` (5); plus `property_tables.png`.

## 5. Validate and improve

| check | prediction | result from the table | verdict |
|---|---|---|---|
| Clausius–Clapeyron, water, p_sat ≥ 10 kPa, T < 344 °C | 0 | within 0.2 %, median 0.005 % | table self-consistent |
| Clausius–Clapeyron, water, 5–15 °C | 0 | up to 0.96 % | p_sat printed to 2–3 figures |
| Clausius–Clapeyron, ammonia | 0 | median 0.009 %; −0.5 % at 130 °C | near-critical steepness |
| h − u − pv, 376 entries | within round-off bound | largest 0.95 of bound | no corrupted entry |
| Z at the lowest isobar | → 1 | 0.9998 (steam, 10 kPa, 300 °C); 0.9933 (ammonia, 100 kPa, 60 °C) | ideal-gas limit reached |
| superheated steam Δh at 10 kPa | NASA ideal gas | within 0.25 % | independent route agrees |
| p_sat interpolation | log ≫ linear | about 25–40 times smaller error | use ln p against 1/T |
| 500 kPa, 300 °C steam | Yan Table A2 row | v −0.015 %, h −0.007 %, s −0.11 % | wide isobar spacing |

The engineering answer: **a sealed ammonia tank filled at 15 °C needs at least 728.5 kPa to stay liquid and reaches 1555.4 kPa at 40 °C; its liquid fill must not exceed 93.74 % of the volume if it may soak to 40 °C (91.03 % for 50 °C), or it becomes liquid-full and its pressure is no longer bounded by the saturation curve.** The answer depends on v_f only, and the table's v_f values pass every consistency check here. What to improve: the temperature limit itself (a tank in the sun can exceed the air temperature), then the tank's own expansion and the dissolved pressurant, both of which this model leaves out.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivation of Clausius–Clapeyron and of the fill limit; `property_tables.png`; the five investigation tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 4, 5, 6 and 7) answering the engineering question and stating which checks the data passed.

### Adversarial self-checks

- **What if the maximum temperature rises by 10 K?** From 40 to 50 °C the fill limit falls from 93.74 % to 91.03 %, and the pressure rises from 1555.4 to 2034.0 kPa: a tank's fill rule must name its temperature.
- **What if the fill temperature is lower?** The colder liquid is denser, so the same mass occupies less volume, and the allowed fill fraction (by volume at filling) falls further (Problem 4).
- **Which limiting cases does the lookup get right?** x = 0 and x = 1 return the saturated values; at the critical point v_f = v_g; at low pressure Z → 1 and the table agrees with the ideal-gas fit.
- **Which assumption fails first in the real world?** The equilibrium assumption during fast draining (flashing, superheated liquid), then the rigid, isothermal tank.
- **What would falsify the table?** A row whose h, u, p and v disagree beyond rounding, or a Clausius–Clapeyron error far above its neighbours': both point at a typing error in that row, and both checks are in your code.
- **Does the result survive a change of method?** Logarithmic and linear p_sat interpolation agree at the nodes and differ between them by up to 5.7 %; the fill limit uses only tabulated v_f at node temperatures (15 °C and 40 °C are rows), so it does not depend on the interpolation at all.

## Check yourself

1. Why do p and T not fix the state of a saturated mixture? Name two properties that do.
2. Water at 120 °C and 300 kPa: compressed liquid, saturated or superheated?
3. Define quality and write the lever rule for h.
4. State the Clausius–Clapeyron equation and the assumption that turns it into "ln p is linear in 1/T".
5. Why does linear interpolation of p_sat against T overestimate p_sat between rows?
6. What is the compressibility factor, and what value does it approach at low pressure?
7. A rigid tank of ammonia warms. Why does the liquid fraction by volume increase even though some liquid evaporates?
8. In the identity check h = u + pv, why is a residual of 0.4 kJ/kg acceptable for saturated water vapour at 0.01 °C but not at 200 °C?
