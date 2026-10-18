# P07.09  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numbers are from `solution.py` and the probe scripts that call it, given to the figures shown. Moisture rates, flows and the skin temperature are typical, illustrative values. Dew points below 0.01 °C are Clausius–Clapeyron estimates over supercooled water.

## 1. The cabin's air at altitude (s = 15: 7500 ft = 2286.0 m)

**Assumptions.** ISA cabin pressure; dry-air composition fixed; ideal gases.

(a) `isa_state(2286.0)`: T = 273.29 K, **p = 76 712.6 Pa**. p_O2 = 0.2095 p = **16 071 Pa** against **21 228 Pa** at sea level: the passenger breathes **75.7 %** of the sea-level oxygen partial pressure (the ratio of the pressures, since x_O2 does not change).

(b) Mass fractions: **N₂ 0.7551, O₂ 0.2314, Ar 0.0128, CO₂ 0.0006** (M = 28.9657 g/mol). Entropy of mixing **4.7167 J/(mol·K) = 162.84 J/(kg·K)**.

(c)

| | ω (g/kg) | dew point (°C) | h (kJ/kg dry air) |
|---|---|---|---|
| cabin, 76.7 kPa | 10.912 | 11.09 | 49.83 |
| sea level | 8.226 | 11.09 | 43.00 |

p_v = 0.5 p_sat(22 °C) = 1322.7 Pa in both cases, so the **dew point**, which depends on p_v only, is the same. ω = ε p_v/(p − p_v) grows as p falls (less dry air per cubic metre for the same vapour), and h grows with ω because each extra gram of vapour carries about 2.54 kJ.

**Independent check.** By hand, ω = 0.62195 × 1322.7/(76 712.6 − 1322.7) = 10.91 g/kg.

**Common mistakes.** Using the sea-level psychrometric chart at cabin pressure; expecting the dew point to change with total pressure.

## 2. Derivations

(a) For the same V and T: p_a V = n_a R T and p_v V = n_v R T, so n_v/n_a = p_v/p_a. Multiply by M_w/M_a: **ω = m_v/m_a = (M_w/M_a) p_v/(p − p_v) = ε p_v/(p − p_v)**, with p_a = p − p_v by Dalton. Inverse: **p_v = ω p/(ε + ω)**. ε = 18.01528/28.9657 = **0.62195**.

(b) Dry air, water and energy: ṁ₃ = ṁ₁ + ṁ₂, ṁ₃ω₃ = ṁ₁ω₁ + ṁ₂ω₂, ṁ₃h₃ = ṁ₁h₁ + ṁ₂h₂. So ω₃ = (ṁ₁ω₁ + ṁ₂ω₂)/ṁ₃ and h₃ = (ṁ₁h₁ + ṁ₂h₂)/ṁ₃, from which ṁ₁(ω₃ − ω₁) = ṁ₂(ω₂ − ω₃) and ṁ₁(h₃ − h₁) = ṁ₂(h₂ − h₃). Hence **(ω₃ − ω₁)/(ω₂ − ω₃) = (h₃ − h₁)/(h₂ − h₃) = ṁ₂/ṁ₁**: the point (ω₃, h₃) is on the segment joining (ω₁, h₁) and (ω₂, h₂), dividing it in inverse proportion to the flows. For temperature: h = h_a(T) + ω h_v(T) ≈ (c_pa + ω c_pv)(T − T₀) + ω h_g0 near room temperature. Subtracting the latent parts (which balance by the water equation), ṁ₁c₁(T₃ − T₁) = ṁ₂c₂(T₂ − T₃) with c_k = c_pa + ω_k c_pv (and c_pa itself varies a little with T): T₃ is a mean weighted by ṁ c, not by ṁ. For equal flows at 2 °C, 4.35 g/kg and 24 °C, 8.00 g/kg the code gives **13.039 °C**, not 13.000 °C: the moister, warmer stream has the larger c.

(c) Control volume: cabin + manifold. In: fresh air (ṁ_f, ω_out) and the occupants' ṁ_w. Out: outflow (ṁ_f, ω_cabin); the recirculated stream crosses no boundary. Water: ṁ_f ω_out + ṁ_w = ṁ_f ω_cabin ⇒ **ω_cabin = ω_out + ṁ_w/ṁ_f**. Energy on the same volume: ṁ_f h(T_pack, ω_out) + Q_sens + ṁ_w h_v(T_cabin) = ṁ_f h(T_cabin, ω_cabin). Every term is independent of ṁ_r, so h_pack and hence **T_pack do not depend on recirculation**. (The code confirms: 2.73 °C at fractions 0 to 0.6.)

**Common mistakes.** Drawing the control volume around the cabin alone and then wondering where the recirculated water goes; averaging temperatures instead of enthalpies.

## 3. A spacecraft cabin with no fresh air (s = 15: CHX at 5.5 °C)

**Assumptions.** Steady; the CHX returns saturated air at its outlet temperature; cabin well mixed; the 1.8 kg/day per crew member is illustrative.

(a) ṁ_w = 4 × 1.8/86 400 = **83.33 mg/s** (7.2 kg/day). ω_out = ω_sat(5.5 °C, 101.325 kPa) = **5.5811 g/kg**. A 10 °C dew point means p_v = p_sat(10 °C) = 1230 Pa and ω_cabin = **7.6427 g/kg**. Water balance on the cabin: the CHX removes ṁ_a (ω_cabin − ω_out) = ṁ_w, so **ṁ_a = 8.333 × 10⁻⁵/(7.6427 − 5.5811) × 10³ = 0.04042 kg/s = 2.425 kg/min** of dry air.

(b) φ = 1230/p_sat(22 °C) = **46.5 %**.

(c) `cool_dehumidify(0.04042, 295.15, 0.0076427, 278.65, 101325)`: condensate 7.200 kg/day (closes the balance), **Q = 886.6 W**, of which **669.8 W** is sensible (dry air from 22 to 5.5 °C) and 216.8 W is the latent heat of the condensed water and the cooling of the vapour.

(d) ω_cabin − ω_out doubles, to 4.123 g/kg: ω_cabin = 9.704 g/kg, φ = **58.8 %**, dew point 13.6 °C (the code gives 58.85 % and 13.56 °C).

**Independent check.** Latent heat alone: 8.333 × 10⁻⁵ kg/s × h_fg ≈ 2.48 MJ/kg (the table's 2489.04 kJ/kg at 5 °C) ≈ 207 W, close to the 216.8 W non-sensible part (the rest is cooling the vapour that stays in the air).

**Common mistakes.** Forgetting that in a closed cabin the water must be condensed out; using the cabin RH instead of the dew point to set ω.

## 4. Your own airliner (s = 15: 180 occupants, 0.35 kg/min each, 40 % recirculation)

(a) `ecs_cruise(n_occ=180, m_fresh_per_occ=0.35/60, recirc_fraction=0.4)`: fresh **1.050 kg/s**, recirculated **0.700 kg/s**, supply **1.750 kg/s**; **ω_cabin = 2.907 g/kg, φ = 13.24 %, dew point −7.43 °C (estimate)**; supply **13.49 °C**, pack outlet **7.78 °C**, pack duty **204.3 kW**.

(b) ṁ_w = 180 × 0.060/3600 = 3.000 × 10⁻³ kg/s; ω = 5 × 10⁻⁵ + 3.000 × 10⁻³/1.050 = 2.907 × 10⁻³. p_v = 0.0029071 × 75 262/(0.62195 + 0.0029071) = 350.2 Pa; φ = 350.2/2645.4 = **13.2 %**.

(c) `brentq` on the fresh flow per occupant with φ = 0.20 gives **0.2297 kg/min per occupant**, with the pack outlet at **0.33 °C** and the cabin dew point −1.94 °C (estimate). It is feasible thermally (the outside air is too dry for the pack to ice at cruise), but it cuts the fresh air by a third: CO₂ and odour rise in proportion, and certification sets a minimum fresh-air flow per occupant, so humidity is bought with air quality.

**Common mistakes.** Computing ṁ_w per minute and ṁ_f per second; bracketing `brentq` so wide that the pack outlet falls below the NASA fit's 200 K and the solver fails.

## 5. How ideal is water vapour?

(a)

| T (°C) | table (Pa) | NASA Gibbs (Pa) | difference |
|---|---|---|---|
| 5 | 870.0 | 873.0 | +0.35 % |
| 20 | 2340.0 | 2338.4 | −0.07 % |
| 40 | 7380.0 | 7369.0 | −0.15 % |
| 60 | 19 950.0 | 19 846.9 | −0.52 % |

(b) R_w = 461.52 J/(kg·K). Z_g = **0.99632, 0.99893, 0.99651, 0.99483** at 5, 20, 40, 60 °C.

(c) Clausius–Clapeyron between 20 and 40 °C: **2433.2 kJ/kg** from the table, **2431.5 kJ/kg** from the NASA route; the table's h_fg at 30 °C is **2429.82 kJ/kg**. Both are within 0.15 %.

(d) The NASA estimate assumes the vapour is an ideal gas (Z = 1). Real saturated vapour has Z slightly below 1, more so as T and p_sat rise: its molecules attract each other, so the real vapour's Gibbs energy at a given pressure is a little lower than the ideal gas's, and the ideal-gas estimate under-predicts p_sat by roughly 1 − Z. At 60 °C the difference is −0.52 % and 1 − Z = 0.52 %. At 5 °C neither the "+0.35 %" nor Z = 0.9963 is real-gas behaviour: the table prints p_sat = 0.00087 MPa, a resolution of 10 Pa on about 873 Pa, so both carry up to ±0.6 % of rounding.

**Common mistakes.** Using R for air instead of R_w; treating a two-significant-figure table entry as exact.

## 6. Should the airline humidify the cabin?

(a) Water exported through the outflow valve at steady state is ṁ_f(ω − ω_out): **12.15 kg/h** at 20 % RH and **18.35 kg/h** at 30 % RH. The occupants supply 11.16 kg/h, so a humidifier must add **0.99 kg/h** (20 %) or **7.19 kg/h** (30 %): **9.9 kg** and **71.9 kg** of water over 10 hours.

(b) Cabin dew points: **−3.1 °C** (18.4 %, estimate), **−1.9 °C** (20 %, estimate), **3.7 °C** (30 %, from the table). All are far above −30 °C: any cabin air that reaches the skin deposits frost in the insulation, and more at higher humidity. On descent the frost melts and the water runs along the frames and drips into the cabin or the bilge, with corrosion and wet, heavier insulation as the result.

(c) For: comfort (eyes, skin, mucous membranes), less static electricity, possibly fewer passenger complaints on long flights. Against: water mass carried (72 kg here), condensation and corrosion in the structure, microbial growth in humidifiers and wet blankets, energy to evaporate the water. The number to measure first is the occupants' moisture output per person (the cabin humidity is proportional to it, and this lesson's 60 g/h is illustrative); the second is how much water is already lost to cold surfaces.

**Common mistakes.** Treating the humidifier water as "free" because it is small per minute; forgetting the skin is far below the dew point at any useful humidity.

## Check-yourself answers

1. Ideal-gas molecules do not interact, so each component's pressure in the full volume and its volume at the full pressure are both proportional to its amount: p_i/p = V_i/V = x_i. Real gases interact, so their pressures and volumes do not simply add, and the two rules give different, approximate answers.
2. ω = ε p_v/(p − p_v), ε = 0.62195. Dry air passes unchanged through humidifiers, coolers and cabins while the water flow changes, so per kilogram of dry air is the invariant basis for the balances.
3. Yes, both. With no water added the mole fraction of vapour is fixed, so p_v = x_v p falls in proportion to p while p_sat(22 °C) does not: φ falls to 50 % × 75/101.325 ≈ 37 %, and the dew point falls because p_v fell. (Contrast Problem 1, where φ was held at 50 % at both pressures and ω had to differ.)
4. Only on the mole fractions (and R). A gas mixed with itself has a single x = 1, so −R Σ x ln x = 0: there is no new arrangement of distinguishable molecules, and no entropy is generated.
5. Mass and energy balances make ω and h exactly mass-weighted averages, so the point is on the straight line in (ω, h). T is not a linear function of h when cp depends on ω and T, so in (ω, T) the mixed point is slightly off the line.
6. At steady state water leaves only through the outflow valve, whose flow equals the fresh flow; the recirculation loop is internal to a control volume around the cabin and manifold, so it cannot appear in that balance.
7. The water balance (ṁ_cond = ṁ_a(ω_in − ω_sat(T_out))); the heat removed includes the latent heat of the condensed water, which for 47.7 kg/h is several tens of kilowatts on top of the sensible cooling.
8. Its vapour pressure (486.9 Pa) is below the triple-point pressure, where the table stops; the module extrapolates Clausius–Clapeyron over supercooled liquid, which is not checked against data. Over ice the saturation pressure is lower than over supercooled water at the same temperature, so the frost point would be somewhat higher (warmer) than the estimate.
