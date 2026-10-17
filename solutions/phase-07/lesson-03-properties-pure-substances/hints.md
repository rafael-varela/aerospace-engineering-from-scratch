# P07.03  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Compare v with v_f(T) and v_g(T) first; compare p with p_sat(T) for (b). The 500 kPa isobar is not embedded, so the lookup interpolates between 100 and 1000 kPa.
2. (a) dg = v dp − s dT for each phase; equal changes keep g_f = g_g; h_fg = T s_fg. (b) Central difference of ln p against 1/T in kelvin. (d) Write the chord as a weighted average and use the definition of convexity.
3. Throttling to the chamber keeps h constant: h_f(20 °C) = h_f(T_sat) + x h_fg(T_sat).
4. v = V/m is fixed. The fill limit makes v equal v_f(T_max). For (d), `brentq` on v_f(T) − V/m.
5. Z = pv/(RT) with R = R_U/M. If 1 − Z = k p, one isobar gives k.
6. Ask: how close is the state to the dome or the critical point, and is that region in the embedded tables?
