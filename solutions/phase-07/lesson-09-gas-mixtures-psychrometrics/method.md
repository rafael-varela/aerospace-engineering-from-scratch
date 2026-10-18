# P07.09  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. z = (6000 + 100 s) × 0.3048 m; `isa_state(z)` and `isa_state(0)`. 2. p_O2 = 0.2095 p at both; ratio. 3. `mass_fractions(AIR_X)`; `entropy_of_mixing(AIR_X)` and divide by M for per kilogram. 4. p_v = 0.5 p_sat(295.15 K); ω from `humidity_ratio` at each p; dew point from `dew_point`; h from `moist_air_enthalpy`. 5. Say which depend on p_v only.

## 2.
(a) Ideal gas for each component in the same V and T; ratio n_v/n_a = p_v/p_a; multiply by M_w/M_a; p_a = p − p_v; invert for p_v. (b) Three balances; write ω₃ and h₃ as weighted means; form (ω₃ − ω₁)/(ω₂ − ω₃) and (h₃ − h₁)/(h₂ − h₃); both equal ṁ₂/ṁ₁, hence collinear. Then expand h₃ = h_a(T₃) + ω₃ h_v(T₃) and show T₃ is a weighted mean of T with weights that depend on cp of each stream. (c) Water balance on cabin + manifold; energy balance on cabin + manifold; recirculated terms cancel.

## 3.
1. ω_out = `omega_saturated(T_CHX, p)`; ω_cabin = `humidity_ratio(p_sat(283.15), p)`. 2. ṁ_w = 4 × 1.8/86 400 kg/s. 3. ṁ_a = ṁ_w/(ω_cabin − ω_out). 4. φ = `relative_humidity(295.15, ω_cabin, p)`. 5. `cool_dehumidify(ṁ_a, 295.15, ω_cabin, T_CHX, p)`; sensible part ṁ_a (h_a(295.15) − h_a(T_CHX)). 6. Halving ṁ_a doubles ω_cabin − ω_out.

## 4.
1. `ecs_cruise(n_occ=150+2s, m_fresh_per_occ=(0.2+s/100)/60, recirc_fraction=0.4)`. 2. ω_c = ω_out + ṁ_w/ṁ_f; p_v = ω p/(ε + ω); φ = p_v/p_sat(295.15 K). 3. `brentq` on the fresh flow per occupant for φ = 0.20 with a bracket that keeps the pack outlet above 200 K; read T_pack.

## 5.
1. Two functions at four temperatures. 2. R_w = R/M_w; Z = p_sat v_g/(R_w T) with p_sat from the table. 3. h_fg = −(R/M_w) [ln p(T₂) − ln p(T₁)]/(1/T₂ − 1/T₁) for each route. 4. Relate the sign of (p_NASA/p_table − 1) to 1 − Z; identify table rounding at low T.

## 6.
1. ω for φ = 0.2 and 0.3 at 22 °C and cabin pressure; water exported = ṁ_f (ω − ω_out) × 3600; subtract the occupants' 11.16 kg/h; multiply by 10 h. 2. Dew points with `dew_point`; compare with −30 °C. 3. Reasons; identify the most uncertain input from the lesson's sensitivity study.
