# P07.07  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. With x = π^0.2857 and τ = 1400/288.15, the compressor work is cpT₁(x − 1)/η_c and the turbine work η_t cpT₃(1 − 1/x).
2. With ε = 1 and ideal components the air enters the combustor at the turbine exit temperature, so the heat added equals the turbine work. For (b), minimise a^m + (π/a)^m.
3. The recuperator's pressure losses are the same at every ε in this problem, and the turbine does not know what the recuperator does downstream. Fuel per unit work is q_in/w_net.
4. Use `np.argmax` on the arrays `sweep` returns, and index the π array with it.
5. Heat crossing a large temperature difference generates entropy; so does hot exhaust thrown away to the atmosphere.
6. Intercooling lowers the HP compressor delivery temperature: compare T₄ with T₇ in each configuration.
