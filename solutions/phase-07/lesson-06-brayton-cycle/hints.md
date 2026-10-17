# P07.06  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The two isentropic temperature ratios are equal, T₂/T₁ = T₃/T₄ = π^0.2857, so T₄ follows from T₃ in one line.
2. Heat added and heat rejected carry the same entropy change, because the isentropes join the two isobars. For (c), along an isentrope dh = v dp.
3. Put x = √τ into η = 1 − 1/x. The radiator needs area for q_out, and the net power is w_net; the ratio you need is (q_out/w_net)/T₁⁴.
4. `brayton_states(gas, ...)` at the numerical optimum gives both T₂ and T₄. The isentropic legs are vertical in the T–s plane.
5. Carnot needs all heat added at T₃ and all heat rejected at T₁. Where does the ideal Brayton cycle add and reject its heat?
