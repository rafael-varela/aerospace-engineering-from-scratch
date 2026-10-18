# P07.03  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) `phase_from_Tv("water", T, 0.1)`; if a mixture, `saturation_state` for p and `mixture_property` for h, u, s with the returned x. (b) `phase_from_Tp("water", T, 2000)`; note there is no compressed-liquid table. (c) `phase_from_Tp`; `superheated_state("water", 500, 300)`; `compressibility_factor`; percentage differences against Yan's row.

## 2.
(a) Along the curve dg_f = dg_g; substitute dg = v dp − s dT; solve for dp/dT; replace s_fg by h_fg/T. (b) v_fg ≈ RT/p; dp/p = (h_fg/R) dT/T²; d(1/T) = −dT/T². Numerically: rows at 95 and 105 °C, slope of ln p against 1/T, times −R. (c) V = m_f v_f + m_g v_g, divide by m. (d) Convexity: f(λa + (1 − λ)b) ≤ λf(a) + (1 − λ)f(b).

## 3.
(a) `saturation_temperature("water", p)`. (b) h_feed = h_f(20 °C); `quality("water", T_sat, "h", h_feed)`. (c) h_g(T_sat) − h_f(20 °C). (d) Look at the first row of the table (the triple point).

## 4.
(a) `saturation_pressure` at both temperatures. (b) `max_liquid_fill`. (c) m = 0.9V/v_f + 0.1V/v_g at T_fill; `rigid_tank_heating`. (d) m for 97 %; `brentq(lambda T: saturation_state("ammonia", T)["vf"] - V/m, T_fill, 120)`.

## 5.
(a) Loop over isobars with `superheated_state` and `compressibility_factor`. (b) k = (1 − Z)/p from each isobar; p₀.₉₉ = 0.01/k. (c) `superheated_dh_vs_ideal_gas`. (d) ρ_ideal/ρ_real = Z.

## 6.
Locate each state relative to the dome, the critical point and the embedded table range; decide whether Z is near 1; name the source.
