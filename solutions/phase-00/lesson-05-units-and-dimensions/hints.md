# P00.05  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Write each conversion factor as a fraction whose numerator and denominator are the same physical length (or pressure, or density); cancel units like algebraic symbols. The knot is exactly 1852 m per hour.
2. Both numbers in (a) are definitions, so the product is a definition. For (b), integrate `F = dp/dt` over the burn; what did you assume about `m` while integrating? For (c), divide the unit of impulse by the unit of mass and simplify.
3. The datasheet does not say which gallon. Compute both; the ratio of the two gallons is the ratio of the two masses.
4. The missed impulse per firing is `x (lbf·s) − x (N·s)`, both expressed in N·s. Position error from a constant velocity error is `δv · T`; convert 30 days to seconds and the answer to km.
5. `Q_(np.array([...]), "knot")` is a single Quantity holding an array; your `dynamic_pressure` already works on it. Knots and km/h have the same dimension.
6. `pint` checks dimension, not meaning. Ask of each case: is the mistake visible in the *dimension*? Temperatures with offsets are a special case documented under "non-multiplicative units".
7. Parse the header with a regular expression, build `ureg(unit_string)`, and let `.to("N*s")` do the dimension check for you. For the paragraph, look at how many of the Board's contributing causes are about process rather than data format.
