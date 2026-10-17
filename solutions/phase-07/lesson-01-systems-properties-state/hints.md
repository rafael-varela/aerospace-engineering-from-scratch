# P07.01  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Per kilogram means divide the molar value by M in kg/mol. The pressure enters s only through −R ln(p/p°).
2. (a) Integrate cp dT for h and cp dT/T for s°, term by term; then divide h by RT. (b) With only a₃, h/(RT) = a₃ + b₁/T. (c) Differentiate h = u + RT with respect to T. (d) Each species sits at its partial pressure x_i p.
3. γ = cp/(cp − R). f = 2/(γ − 1). The bending vibration of CO₂ is soft; it begins to store energy at low temperature.
4. `exit_temperature(800, q)` and `exit_temperature(800, q, constant_cp=True)`. Mean cp is q divided by the temperature rise.
5. `air_s(T1, p1)` returns a NumPy float; pass `float(s1)` to `air_state`. For an ideal gas, h is a function of T alone.
6. For each case ask: is the gas dense, hot enough to vibrate strongly or dissociate, or not air at all?
