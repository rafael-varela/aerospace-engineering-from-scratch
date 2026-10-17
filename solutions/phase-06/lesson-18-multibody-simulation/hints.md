# P06.18  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. dF/ds of the gas force is n A_a² p/V with p and V at the current stroke; the frequency is (1/2π)√(k/M).
2. (a) Write r₂ in terms of (y₁, s), differentiate, and collect T as ½q̇ᵀM_r q̇. A virtual δy₁ does not change the strut length. (d) Allowed velocities satisfy Φ_q q̇ = 0; compute the power q̇ᵀΦ_qᵀλ.
3. With a constant force the energy is F × stroke; the work of gravity over the stroke is part of what the legs absorb.
4. Each run is a few tenths of a second of simulation; reuse `drop_metrics` and `energy_terms`. Small orifices make a large force when the stroke rate is high; large ones leave the work to the gas.
5. Without lift, gravity keeps doing work on the whole aircraft through the stroke and the tyre deflection. The bottoming condition is a function of one variable, the sink rate.
6. The orifice force peaks when the stroke rate peaks; the tyre force rises as fast as the tyre lets it. Which of these sets the time of the first peak?
