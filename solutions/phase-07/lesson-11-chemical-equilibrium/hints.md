# P07.11  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Count moles: (1 − α) H₂O, α H₂, α/2 O₂, total 1 + α/2; each mole fraction times p/p° goes into K_p.
2. (a) Differentiate n_j ln n_j − n_j ln n carefully: the ln n term contributes −Σ n_k/n = −1, which cancels the +1. (c) Multiply the first equation by a_kj n_j and sum over j.
3. NO comes from ½N₂ + ½O₂; at a lean mixture both are abundant, so x_NO ≈ K_p √(x_N₂ x_O₂). Its K_p has a large ΔH°.
4. `oxidizer_fuel_to_phi("H2", of)` then `equilibrium_flame_temperature("H2", phi, "O2", p)`; M̄ = Σ x_j M_j.
5. Compare relative differences, not absolute ones, and look at the trace species.
6. Pass your own `species` list to `equilibrium_tp` inside a `brentq` on the energy balance.
