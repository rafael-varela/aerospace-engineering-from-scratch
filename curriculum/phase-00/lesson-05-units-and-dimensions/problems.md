# P00.05  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer; where a conversion is exact, say so.

## 1. (Warm-up: four conversions, twice)

Convert by hand, writing each factor as a fraction equal to one (e.g. `1852 m / 1 nmi`), then check with `pint`:

(a) 35 000 ft to m  (b) 250 knots to m/s  (c) 29.92 inHg to Pa  (d) 1 slug/ft³ to kg/m³.

Answer form: four numbers with units, the hand and `pint` values side by side, and the relative difference.

## 2. (Derivation: the pound-force and the velocity change)

(a) From the definitions 1 lb = 0.453 592 37 kg and `g₀` = 9.806 65 m/s², derive 1 lbf in newtons to ten significant figures and state why no uncertainty attaches to it.
(b) Starting from Newton's second law, derive `Δv = I / m` for a short firing, stating the assumption that makes it exact and the condition under which it fails by more than 1 %.
(c) Show that pound-force-seconds per pound-mass has the dimension of velocity, evaluate it in m/s, and explain why `I_sp` in "seconds" is nevertheless meaningful.

Answer form: derivations with every unit written, two numbers.

## 3. (Unfamiliar scenario: the fuel order)

A datasheet gives jet fuel density as 6.7 lb/gal. A tanker is ordered to deliver 2 400 gal.
(a) What mass of fuel arrives, in kg, if "gal" is the US gallon (3.785 411 784 L)?
(b) What mass if the supplier reads it as the imperial gallon (4.546 09 L)?
(c) Write a `pint` expression that computes (a) and *raises* if someone passes the volume as a plain float.
(d) The aircraft's take-off calculation needs the fuel mass to ±2 %. Is the gallon ambiguity within that tolerance?

Answer form: two masses with units, one code snippet, a yes/no with the percentage.

## 4. (Numerical, uses your code: the Mars Climate Orbiter re-enactment)

Seed a generator with the last four digits of your phone number. Draw `n` firing impulses in lbf·s uniformly from 0.2 to 0.6, with `n` of your choice between 200 and 2 000, for a spacecraft of mass `m` of your choice between 300 and 700 kg. Using `accumulated_delta_v_error`, report the missed `Δv`. Then, treating the missed `Δv` as a velocity error sustained for 30 days, estimate the resulting position error in km and compare it with the ~170 km altitude error the Board reports. Say clearly which of your inputs are illustrative.

Answer form: seed, `n`, `m`, missed `Δv` in m/s, a position error in km, and a two-sentence comparison.

## 5. (Numerical, uses your code: a flight log with units)

Build a NumPy array of true airspeeds 180, 200, 220, 240, 260 knots and wrap it in a single `Quantity`. With ISA density at 11 km (0.3639 kg/m³), compute dynamic pressure with your `dynamic_pressure` in Pa and in psi, and plot it against airspeed using your P00.04 module with units on both axes. Then repeat with airspeeds mistakenly given in km/h and show the two curves on one figure.

Answer form: a table of five rows, a figure, and one sentence on how far apart the curves are and why `pint` did not object.

## 6. (Judgment: where the units library stops)

For each of the following, say whether `pint` catches the error, and if not, what would:
(a) adding 20 °C to 5 °C; (b) using indicated airspeed where true airspeed is required; (c) passing `I_sp` in seconds to a function expecting m/s; (d) a torque in N·m passed to a function expecting energy in J; (e) `Q_(1, "m") + 1.0`.

Answer form: five one-line verdicts with reasons.

## 7. (Hard: design the interface that would have saved the mission)

Write a function `read_small_forces(path) -> Quantity` that reads a text file of one impulse per line, requires a header line of the form `# impulse [<unit>]`, parses the unit with your registry, and returns the impulses as a single `Quantity` array in N·s. It must raise a clear error if the header is missing, if the unit is not an impulse, or if any line is not a number. Write three tests: a correct lbf·s file converts to N·s with the right factor; a file with no header raises; a file whose header says `[lbf]` raises `DimensionalityError`. Then argue in one paragraph whether this interface would have been sufficient to prevent the MCO loss, using the Board's contributing causes.

Answer form: code, passing tests (use `tmp_path`), one paragraph.
