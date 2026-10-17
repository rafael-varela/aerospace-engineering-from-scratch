# P08.19  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The volume is the bottle's area times the water height, including the 30 mm below the hole. The draining time goes as D²/(C_d d²), so each relative uncertainty enters with its exponent.
2. (a) Keep V_s = (C_d a/A)V_j in Bernoulli's equation and solve for V_j. (c) Differentiate V ∝ √(L sin θ) with respect to L and θ. (d) Compare the Laplace pressure of a cylinder (σ/r) and of a sphere (2σ/r) with ρgh.
3. Read the drain file with `read_csv`, split by run, and call `fit_discharge_coefficient`. For the truth, compute what C_d would look like through the nominal hole and bottle.
4. In a contraction, continuity gives V_sc = V_ts/CR; Bernoulli between the two wall taps then has both velocities in it. For the aeroplane, use the ISA density at 3000 m.
5. u_V/V = ½√((u_L/L)² + (u_θ/tan θ)²) with L = Δh/sin θ; the two terms pull in opposite directions as θ shrinks.
6. Compare each Reynolds number with the range over which the quantity you are measuring (C_d, or the pitot reading) is known to be independent of it.
