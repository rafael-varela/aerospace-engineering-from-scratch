# P07.03  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which values come straight from a table row, which are interpolated, and which are typical, illustrative inputs.

Several problems use your personal seed **s**: the day of the month you were born (1–31). Tables: the embedded subsets of Yan's Appendices A and B (NIST SRD 69 values). Where a problem asks for a value Yan tabulates but the code does not embed, read it from Yan's appendix page.

## 1. (Warm-up: three water states)

Let T = 100 + 5s °C.

(a) Water at T with v = 0.1 m³/kg. Identify the phase. If it is a mixture, find p, x, h, u and s; if not, say what more you would need.
(b) Water at T and 2000 kPa. Identify the phase. Why can your lookup not give its properties?
(c) Water at 500 kPa and 300 °C. Identify the phase, interpolate v, u, h and s with `superheated_state`, compute Z, and compare with Yan's Table A2 row at 500 kPa, 300 °C.

Answer form: a phase and five numbers (or a sentence); a phase and a sentence; a phase, five numbers and four percentage differences.

## 2. (Derivation: what the table must obey)

(a) Derive the Clapeyron equation dp_sat/dT = h_fg/(T v_fg) from g_f = g_g along the saturation curve.
(b) Show that with v_f ≪ v_g and v_g = RT/p it becomes d ln p/d(1/T) = −h_fg/R. Use the 95 °C and 105 °C rows of the water table to estimate h_fg at 100 °C this way, and explain why it differs from the tabulated value.
(c) Derive the lever rule y = y_f + x y_fg from the definition of quality and the additivity of extensive properties.
(d) Prove that the chord of a convex function lies above it, and so that linear interpolation of p_sat(T) overestimates p_sat between rows.

Answer form: three derivations, one number with a sentence, one proof.

## 3. (Unfamiliar scenario: a spacecraft flash evaporator)

A flash evaporator rejects heat by spraying water into a chamber held at low pressure; part of it evaporates and is vented. Take the chamber pressure as p = 0.7 + 0.02s kPa and the water feed at 20 °C as saturated liquid (typical, illustrative).

(a) Find the chamber's saturation temperature with `saturation_temperature`.
(b) If no heat is supplied, what fraction of the feed flashes to vapour as it drops to the chamber state (adiabatic, h constant)? Use `quality` with h.
(c) In operation the feed is evaporated completely by heat from a coolant loop. How much heat does each kilogram of water carry away?
(d) Why must the chamber pressure stay above about 0.61 kPa, and what happens below it?

Answer form: one temperature, one fraction, one heat, two sentences.

## 4. (Numerical, uses your code: your ammonia tank)

A rigid 1 m³ ammonia tank is filled at T_fill = 5 + s/2 °C and may soak to T_max = 35 + s/3 °C (typical, illustrative).

(a) Find p_sat at both temperatures.
(b) Find the largest liquid fill by volume at T_fill that stays two-phase up to T_max.
(c) For a 90 % fill, find the mass, x₁, x₂, the liquid volume fraction at T_max and the heat absorbed (`rigid_tank_heating`).
(d) For a 97 % fill, find the temperature at which the tank becomes liquid-full (solve v_f(T) = V/m).

Answer form: two pressures, one fraction, five numbers, one temperature.

## 5. (Data-based investigation: how ideal is the vapour?)

(a) From the superheated tables, tabulate Z for steam at 300 °C and ammonia at 60 °C on every embedded isobar.
(b) Assuming 1 − Z ∝ p at low pressure, estimate the pressure at which each is ideal to 1 % (Z = 0.99). Use two different isobars and compare the estimates.
(c) Compare superheated-steam enthalpy differences with the NASA ideal-gas fit: 150 → 600 °C at 100 kPa, and 200 → 600 °C at 1000 kPa. Relate the differences to your Z values.
(d) A colleague uses the ideal-gas law for ammonia vapour at 1 MPa and 60 °C to size a vent line. How wrong is the density?

Answer form: a table, two pairs of pressures, two percentages and a sentence, one percentage.

## 6. (Judgment: table, ideal gas, or neither?)

For each case say which property source you would use (this lesson's tables, the ideal gas of `P07.01`, or something else), and why. Two sentences each.

(a) Steam at 10 kPa and 400 °C in a vacuum test chamber.
(b) Ammonia liquid at 20 °C and 3 MPa in a pump inlet.
(c) Water at 370 °C and 21 MPa.
(d) Ammonia vapour leaving a spacecraft radiator loop's evaporator at −10 °C and 250 kPa.
(e) Interpolating the steam table between 1 and 4 MPa at 320 °C.

Answer form: five short paragraphs.
