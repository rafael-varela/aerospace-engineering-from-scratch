# P07.04  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Heat leaving a reservoir at T removes Q/T; heat entering at T adds Q/T. For the blocks, the final temperature is the heat-capacity-weighted mean. A throttle keeps T of an ideal gas.
2. (a) Apply Clausius to a cycle made of the actual path and a reversible return. (c) Write S_gen as a sum of C ln(T_f/T_i) and use ln y ≤ y − 1 on T_i/T_f. (d) Substitute Q_C = Q_H − W.
3. Each step: (heat or work)/(temperature it ends at) minus what it started with. Add them and watch the intermediate terms cancel.
4. Use `np.where` to build the two-segment paths. Mean cp means (h(T₂) − h(T₁))/(T₂ − T₁).
5. `case = dict(BLEED)`, change the keys, then `bleed_air_ledger(case)`. The total throttling from 1.0 to 0.28 MPa is what matters for the valve plus the duct.
6. S_gen = −ΣQ/T. For (d), T₀ multiplies S_gen; S_gen itself does not depend on T₀.
