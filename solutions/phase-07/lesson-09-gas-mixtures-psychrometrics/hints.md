# P07.09  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The oxygen mole fraction does not change with altitude, so its partial pressure scales with the cabin pressure. For (c), hold φ and T fixed: p_v is then fixed, and so is anything that depends only on p_v.
2. (a) Write p_a V = n_a R T and p_v V = n_v R T for the same V and T, divide, and multiply by the molar masses. (b) Both ω₃ and h₃ are the same weighted average; eliminate the weights. (c) Put the control volume around the cabin and the manifold together, so the recirculation loop is internal.
3. With no fresh air, the CHX is the only water sink: the water removed per second is the dry-air flow through it times the drop in ω across it.
4. The flows come straight from the inputs; ω_cabin from the cabin-plus-manifold water balance; RH needs p_v from ω and the table's p_sat at 22 °C.
5. Z = p v/(R_w T) is 1 for an ideal gas; the Gibbs estimate assumes Z = 1 exactly. Look at how many significant figures the table prints for p at 5 °C.
6. At steady state the water exported through the outflow valve per hour is ṁ_fresh (ω_cabin − ω_out); whatever the occupants do not supply, a humidifier must.
