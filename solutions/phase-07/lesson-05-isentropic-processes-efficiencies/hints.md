# P07.05  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The stage is adiabatic, so s_gen is simply s₂ − s₁. For a pressure ratio close to 1 the isentropic and polytropic efficiencies differ only at second order in ln π.
2. Divide T ds = dh − v dp by T and use dh = cp dT, v = RT/p; for the polytropic step, put dh = v dp/η_p into cp dT = dh. For (d), x^k is convex in x, so it lies above its tangent at x = 1.
3. The dead state is the outside air, so T₀ is 216.65 K, not 288.15 K. Whatever the compressor does to the air's temperature, the cabin needs it near 290 K.
4. `isentropic_efficiency(AIR, T1, T[-1], pi)` and `polytropic_efficiency(...)` work on any inlet and exit pair. Each stage's isentropic work is roughly cp T_in (π_st^m − 1): what changes from stage to stage?
5. `polytropic_exit_temperature` with π = 1/4 expands. At 1600 K the vibrational modes of N₂ and O₂ are well excited, so γ is well below 1.4.
6. Casing heat loss cools the delivery air without reducing the shaft work; a bleed removes air partway through the machine.
