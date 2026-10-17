# P07.08  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. For an ideal gas, s(T, p) splits into a part that depends only on T and a part −R ln(p/p°); h depends only on T. Heat removed from a cooling gas comes out over a range of temperatures, not all at T.
2. Multiply the entropy balance by T₀ and subtract it from the energy balance; add and subtract h₀ and T₀s₀ using mass conservation. For (c), integrate T ds along the exit isobar from the isentropic exit state to the actual one.
3. Use `engine_states(EngineSpec())` for station 3 and `flow_exergy(N_AIR, T, p, T0, p0)` at each step. A throttle keeps h constant.
4. `exergy_audit(EngineSpec(T_tit=..., pi_c=...))` returns everything; divide by `Ex_fuel + Ex_ram` for percentages.
5. The product is the same in all three; only the denominator changes. Water at 288 K and the partial pressure of its vapour in the products: would it stay vapour?
6. Chemical energy released as heat into gas at T is worth (1 − T₀/T) of itself at best, and the gas in the combustor is far colder than an ideal reversible reaction would need.
