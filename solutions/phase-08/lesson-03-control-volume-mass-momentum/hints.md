# P08.03  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Draw the control volume round the whole engine and cut the thrust mount; no air enters, the propellant comes from tanks inside. The ambient pressure acts on every face except where the exit pressure acts.
2. (b) Mass in through the inlet must leave through the exit faces or the side; the jet core leaves faster than it entered, so it occupies less area than the air it came from. (d) Adding W to every u_x multiplies W by the net mass outflow.
3. The pressure far above and far below is ambient, so the net pressure force on the stream-tube control volume vanishes; the rotor thrust is the only force. The power is the rate of increase of kinetic energy of the air.
4. `engine_case(dict(CRUISE, mdot_a=...))` does (a). For (c), the hot gas expands further, so both V_e and T_e change and A_e must be recomputed from continuity at p_a.
5. β = J A/(ρ Q²): the power-law profile is nearly flat, the laminar one is not. Gravity loss is g t_b: the longer the burn, the longer gravity works against the thrust.
6. On the stand the engine does not move, so V₀ = 0 if the air is still; the load cell sees whatever the momentum balance says the cut solid must carry.
