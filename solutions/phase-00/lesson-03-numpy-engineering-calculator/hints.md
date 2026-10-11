# P00.03  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The geopotential correction at a few hundred metres is tiny (a few centimetres) but compute it anyway; the exponent 5.2559 is where your hand arithmetic will lose figures, so carry five.
2. ∫ dp/p = ln p. With T constant the right-hand side is just a constant times h. H has the units of a length; check that R T/g₀ does too.
3. Only `T0`, `P0`, `G0`, `R_AIR`, `GAMMA`, `H_BASE` and `LAPSE` describe the planet. If anything else needs changing, your Earth code has a constant hiding in a formula.
4. Density falls with altitude, so your grid of densities is *decreasing*; `np.interp` wants increasing x. Reverse both arrays with `[::-1]`.
5. Give the altitude array a trailing `None` index and the speed array a leading one: shapes (21, 1) and (1, 6) broadcast to (21, 6). For the Mach table, remember that a depends on temperature, not density.
6. For (b), the layer-7 formula is equation (6) with p_b = 3.9564 Pa, T_b = 214.65 K at h_b = 71 km; do not forget to convert 100 km geometric to geopotential first. For the argument, ask what a user of the number would do with it.
