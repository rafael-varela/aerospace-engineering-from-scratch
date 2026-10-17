# P07.03  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numbers come from `solution.py` and the probe script that calls it. Table values are the NIST SRD 69 values as published in Yan's appendices. R = 0.461 523 kJ/(kg·K) for water and 0.488 210 kJ/(kg·K) for ammonia.

## 1. Three water states (s = 15: T = 175 °C)

(a) At 175 °C, v_f = 0.001121 and v_g = 0.216581 m³/kg (a table row). v_f < 0.1 < v_g: **saturated mixture**. p = p_sat = **892.60 kPa**; x = (0.1 − 0.001121)/0.215460 = **0.4589**; h = **1673.40 kJ/kg**, u = **1584.14 kJ/kg**, s = **4.1711 kJ/(kg·K)**. (For s ≥ 29, T ≥ 245 °C, v_g falls below 0.1 and the state is superheated vapour; then you need the superheated table at the pressure that gives v = 0.1, which `state_from_Tv` finds by root finding.)

(b) p_sat(175 °C) = 892.6 kPa < 2000 kPa: **compressed liquid**. The embedded tables contain no compressed-liquid data, so the lookup refuses it; Yan's Table A3 has it, or, approximately, take the saturated liquid at 175 °C for v and u.

(c) p_sat(300 °C) = 8587.9 kPa > 500 kPa: **superheated vapour**. Interpolated between the 100 and 1000 kPa isobars: v = **0.52253 m³/kg**, u = **2803.09 kJ/kg**, h = **3064.35 kJ/kg**, s = **7.4535 kJ/(kg·K)**; Z = **0.9877**. Yan's Table A2 at 500 kPa, 300 °C: 0.52261, 2803.25, 3064.55, 7.4614. Differences: **−0.015 %, −0.006 %, −0.007 %, −0.11 %**.

**Check.** h − u − pv from the interpolated values: 3064.35 − 2803.09 − 500 × 0.52253 = −0.002 kJ/kg, within rounding.

**Common mistakes.** Comparing v with v_g at the wrong temperature; interpolating v linearly in p between widely spaced isobars (v ∝ 1/p, so the error would be tens of per cent; pv linear in p avoids it).

## 2. Derivations

(a) Two phases in equilibrium at (T, p) have g_f = g_g. Moving along the saturation curve, dg_f = dg_g, and for each phase dg = v dp − s dT. So v_f dp − s_f dT = v_g dp − s_g dT, giving dp/dT = s_fg/v_fg. Evaporation at constant T and p is reversible with q = h_fg = T s_fg, so **dp_sat/dT = h_fg/(T v_fg)**.

(b) With v_fg ≈ v_g = RT/p: dp/dT = h_fg p/(RT²), i.e. d ln p = (h_fg/R) dT/T² = −(h_fg/R) d(1/T): **d ln p/d(1/T) = −h_fg/R**. From the 95 and 105 °C rows: h_fg ≈ **2293.2 kJ/kg** against the tabulated **2256.4 kJ/kg** at 100 °C, **1.6 % high**. The neglected v_f is small (0.1 % of v_g at 100 °C); the larger error is the ideal-vapour assumption: Z of saturated steam at 100 °C is below 1, so RT/p overestimates v_g, and the estimate of h_fg is high by about the same fraction.

(c) A mixture of mass m holds m_g = x m of vapour and m_f = (1 − x)m of liquid. Extensive properties add: Y = m_f y_f + m_g y_g. Dividing by m: **y = (1 − x) y_f + x y_g = y_f + x y_fg**.

(d) f is convex if f(λa + (1 − λ)b) ≤ λf(a) + (1 − λ)f(b) for λ in [0, 1]. The linear interpolant between rows a and b at T = λa + (1 − λ)b is exactly λf(a) + (1 − λ)f(b): the chord. So **the chord lies above a convex function** and linear interpolation of p_sat(T) overestimates p_sat between rows. The leave-one-out test confirms it: every linear error is positive.

**Common mistakes.** Using Celsius in 1/T; forgetting that h_fg = T s_fg holds only along the saturation curve.

## 3. A flash evaporator (s = 15: p = 1.0 kPa)

**Assumptions.** Steady flow; adiabatic flashing at constant h; saturated-liquid feed at 20 °C; equilibrium in the chamber.

(a) T_sat(1.0 kPa) = **6.99 °C**.

(b) h_f(20 °C) = 83.91 kJ/kg = h_f(6.99 °C) + x h_fg(6.99 °C) gives **x = 0.0220**: 2.2 % of the feed flashes, and its latent heat cools the other 97.8 % to 6.99 °C.

(c) Evaporated completely: h_g(6.99 °C) − h_f(20 °C) = **2429.79 kJ/kg**. A kilogram of water carries away about 2.43 MJ, which is why evaporators are compact heat sinks for short mission phases.

(d) 0.61 kPa is water's triple-point pressure (the table's first row, 0.01 °C). Below it there is no liquid phase: water sprayed in would freeze and then sublimate. That is the principle of a sublimator, a different device.

**Check.** Energy: the vapour fraction absorbs x h_fg = 0.02195 × 2484.3 = 54.5 kJ per kg of feed, equal to the drop in h_f from 20 °C to 6.99 °C (83.91 − 29.38 = 54.5 kJ/kg); h is conserved, as a throttle requires.

**Common mistakes.** Using h_g(20 °C) for the feed; using the table's p column in MPa as kPa.

## 4. Your ammonia tank (s = 15: T_fill = 12.5 °C, T_max = 40 °C)

(a) p_sat(12.5 °C) = **669.9 kPa** (log interpolation between the 10 and 15 °C rows); p_sat(40 °C) = **1555.4 kPa**.

(b) Fill limit **93.22 %**. (Below the 15 °C case's 93.74 %, because colder liquid is denser and has further to expand.)

(c) m = 0.9/v_f + 0.1/v_g = **559.53 kg**; x₁ = **0.000941**; x₂ = **0.000752**; liquid at 40 °C **96.50 %** of the volume; Q = m(u₂ − u₁) = **72.97 MJ**.

(d) 97 % fill: v = 0.0016594 m³/kg, which equals v_f at **25.09 °C**: the tank is liquid-full 15 K below T_max. Above that the liquid is compressed and the pressure is no longer the saturation pressure; it rises steeply with temperature.

**Check.** At the fill limit, V/m equals v_f(40 °C) = 0.001726 m³/kg to 10⁻¹⁰ (the test checks the same identity at 15 °C).

**Common mistakes.** Taking the fill limit as a mass fraction; forgetting the vapour mass (negligible here, but not in principle).

## 5. How ideal is the vapour?

(a)

| fluid, T | isobars (kPa) | Z |
|---|---|---|
| steam, 300 °C | 10, 100, 1000, 4000 | 0.9998, 0.9976, 0.9753, 0.8902 |
| ammonia, 60 °C | 100, 500, 1000, 2000 | 0.9933, 0.9655, 0.9286, 0.8454 |

(b) With 1 − Z = k p: steam at 300 °C, **415 kPa** from the 100 kPa isobar and **405 kPa** from the 1000 kPa isobar; ammonia at 60 °C, **149 kPa** from 100 kPa and **145 kPa** from 500 kPa. The estimates agree to 3 %, so 1 − Z ∝ p is a good model up to about 1 MPa for steam and 0.5 MPa for ammonia at these temperatures.

(c) At 100 kPa, 150 → 600 °C: table 928.98, NASA 922.67 kJ/kg, **+0.68 %**. At 1000 kPa, 200 → 600 °C: table 870.32, NASA 826.35 kJ/kg, **+5.3 %**. The departure grows with pressure and with closeness to the dome at the cold end (200 °C at 1 MPa is only 20 K above saturation), in line with the Z values. At 10 kPa it was 0.05–0.25 % (lesson, Section 4).

(d) ρ_ideal = p/(RT) = Z ρ_real, so the ideal-gas density is **7.1 % low** (Z = 0.9286). The vent line would be undersized for a given mass flow.

**Common mistakes.** Using R_U instead of R = R_U/M; computing Z with T in °C.

## 6. Judgment

(a) **The ideal gas is fine.** At 10 kPa and 400 °C steam has Z within about 10⁻⁴ of 1 (Z = 0.9998 already at 300 °C), and the NASA fit agrees with the table's enthalpy to 0.1–0.2 %; either source works.

(b) **Compressed liquid: neither embedded table.** p_sat(20 °C) = 857.5 kPa < 3 MPa, so the ammonia is subcooled; use Yan's compressed-liquid data if available, or the saturated liquid at 20 °C for v and u (liquids are nearly incompressible) with h ≈ h_f + v_f(p − p_sat).

(c) **Near-critical: not this lookup.** p_sat(370 °C) = 21 044 kPa, so at 21 MPa the state is just outside the dome, within 4 K of the critical point, where properties vary steeply; the embedded superheated table does not reach it, and linear interpolation would be unreliable even if it did. Use a reference equation of state (the NIST source itself).

(d) **The full Yan table, not the embedded subset.** p_sat(−10 °C) = 290.7 kPa > 250 kPa, so the vapour is slightly superheated, close to the dome where Z is a few per cent below 1; the embedded 500 kPa isobar starts at 10 °C, so interpolate between Yan's 200 and 300 kPa isobars instead.

(e) **Use the intermediate isobars.** Between 1 and 4 MPa at 320 °C, Z falls from about 0.975 to 0.89 and v is far from linear in p; Yan's table has isobars at 1.2, 1.4, 1.6, 1.8, 2, 2.5 and 3 MPa, which cut the interpolation error by an order of magnitude; at minimum interpolate pv, not v.

## Check-yourself answers

1. In the dome p = p_sat(T): fixing one fixes the other, so they count as one property. (T, v), (T, x), (p, h) or (p, s) fix the state.
2. p_sat(120 °C) = 198.67 kPa < 300 kPa: **compressed liquid**.
3. x is the mass fraction of vapour; h = h_f + x h_fg.
4. dp_sat/dT = h_fg/(T v_fg). With v_f neglected and the vapour ideal, d ln p/d(1/T) = −h_fg/R, linear wherever h_fg is nearly constant.
5. p_sat(T) is convex, and the chord of a convex function lies above it.
6. Z = pv/(RT); it tends to 1 as p → 0.
7. The liquid's specific volume v_f grows with T; the tiny mass that evaporates is outweighed by the expansion of the remaining liquid, because v_g ≫ v_f means very little mass is needed to fill the shrinking vapour space.
8. The round-off bound depends on the entry: at 0.01 °C p is printed as 0.00061 MPa and v_g = 205.991 m³/kg, so δp·v alone is about 1 kJ/kg; at 200 °C the bound is about 0.02 kJ/kg, and 0.4 kJ/kg would signal a corrupted entry.
