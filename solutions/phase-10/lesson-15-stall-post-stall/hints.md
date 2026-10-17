# P10.15  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. ((1 + √f)/2)² is 1, 0.729 and 0.25 for f = 1, 0.5, 0. The resultant of cl and cd of the plate is c_n = 2 sin α, perpendicular to the plate.
2. (d) The error e = A − A* obeys e ← [(1 − ω)I + ωJ]e; its eigenvalues are 1 + ω(λ − 1).
3. Forces = q S × coefficient with q = ½ρV²; lift is perpendicular to the flight path, drag along it.
4. Build the twist as `lambda y: -np.radians(w) * np.abs(y) / 5.0` and pair it with the tapered chord from `tapered_wing(0.45)`.
5. `read_ladson_block("I", 8.86e6)`, then `fit_kirchhoff` and `section_stall_angle`; the higher-Re section stalls later, so widen the α bracket.
6. Below stall the section curve has a positive slope everywhere and the solution is unique; past stall, more terms resolve smaller stalled regions.
