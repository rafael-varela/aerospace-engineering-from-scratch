# P00.05  Units & Dimensional Consistency

> A number without a unit is not a measurement. A spacecraft was lost to one.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the conversion factor between any two units of the same dimension from the definitions of the base units, and check it by a round trip. 2. Implement unit-checked engineering calculations with `pint`, so that a dimensionally wrong operation raises an error instead of returning a number. 3. Decide, for a given function interface, whether the units are carried by the data or by the documentation, and name the failure that each choice invites. 4. Reconstruct the Mars Climate Orbiter small-forces error numerically and explain why no single engineer could have seen it in the numbers. |
| **Prerequisites** | `P00.01`, `P00.03`; skill check: you can rearrange `F = m a` for any variable and state the SI unit of each. |
| **Engineering question** | A ground program writes a thruster impulse as a bare number into a file; a navigation program reads it. How large is the trajectory error if the two programs disagree about the unit, and how would you make that disagreement impossible rather than unlikely? |
| **Assumptions and limits** | `pint` checks *dimensions* and converts *units*; it cannot know that a quantity labelled "velocity" is airspeed rather than ground speed, or that a temperature is absolute when you meant a difference. Conversion factors here are exact SI definitions; measured constants (e.g. `g` at a site) are not. The MCO re-enactment uses illustrative firing sizes and spacecraft mass, not the mission's data. |
| **Analytical work** | The pound-force in newtons from the definitions of the pound and standard gravity; dimensional analysis of impulse, velocity change and specific impulse; the order of magnitude of a navigation error from an accumulated impulse error. |
| **Computational deliverable** | `main.py` with a single `UnitRegistry`, `round_trip`, `magnitude_in`, `add_checked`, `dynamic_pressure`, `total_impulse`, `delta_v`, and the MCO functions `lbf_s_to_N_s_factor`, `mco_misread_impulse`, `accumulated_delta_v_error`. Tests prove: round trips return to within 10⁻¹² relative, a mismatched addition raises `DimensionalityError`, swapped arguments are rejected, and the lbf·s/N·s factor equals its exact definition. |
| **Practical activity** | **data-based investigation**: read the Mars Climate Orbiter Mishap Investigation Board report and reconstruct the chain from a unit mismatch to a 170 km altitude error. |
| **Equipment** | None. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Your `main.py`, passing tests, and a one-page causal chain for the MCO loss in your own words with the unit error placed correctly within it. |
| **Assessment** | `assessments/phase-00/grading-rubric.md`; competencies 2 (analytical problem solving), 3 and 4; pass: all public tests green, Problem 2 derivation correct to the exact factor, causal chain names the root cause and at least three contributing causes. |
| **References** | `@pint-tutorial` "Initializing a Registry", "Defining a Quantity", "Converting to different units", "String parsing", "String formatting" (required-free). `@nasa-mco-mib` Mars Climate Orbiter Mishap Investigation Board Phase I Report, 1999-11-10, sections on the root cause and contributing causes (required-free). |
| **Extensions** | ★ Wrap a NumPy array in a `Quantity` and vectorize `dynamic_pressure` over a flight log. ★ Write a `pint` *context* that converts between geometric altitude and pressure altitude using the ISA model. |
| **Time** | unmeasured |

## 1. The engineering problem

On 23 September 1999 the Mars Climate Orbiter fired its engine to enter orbit around Mars, passed behind the planet, and was never heard from again. The Mishap Investigation Board's Phase I report, published seven weeks later, found that the spacecraft had reached Mars about 170 km lower than planned, deep enough into the atmosphere to be destroyed. The root cause it names is a unit error: a ground-software file called "Small Forces" reported the impulse of each small thruster firing in pound-force-seconds, while the navigation software that consumed the file, following the written interface specification, treated the numbers as newton-seconds. Every firing was booked as 4.45 times smaller than it really was. Over nine months of cruise, the small discrepancy in the modelled trajectory accumulated into a fatal one.

Nothing in the numbers looked wrong. The file contained plausible values. The navigation solutions were self-consistent. The error was between two programs that each did exactly what they were written to do. This is the shape of most unit errors: they are not arithmetic mistakes, they are *interface* mistakes, and interfaces are where the units are usually stripped off.

In this lesson you learn to keep the unit attached to the number for as long as possible, using the `pint` library, so that an addition of metres to seconds is a runtime error, a function that wants a density refuses a speed, and a conversion that would silently lose a factor of 4.45 cannot be written. The success criteria: every public test passes, including one that *demands* an exception for a mismatched addition; and you can reconstruct, in your own words and with your own numbers, how a file of plausible-looking impulses lost a spacecraft.

## 2. Understand and derive

### 2.1 Dimensions, units, and numbers

A *dimension* is what kind of quantity something is: length, mass, time, force. A *unit* is an agreed amount of that dimension: metre, foot, kilogram, slug. A *quantity* is a number times a unit: 35 000 ft. Two quantities can be added only if they share a dimension; they can be *compared* numerically only if they are first expressed in the same unit. Multiplication and division are always allowed and produce a new dimension: force × time is impulse, with dimension `[mass][length]/[time]`.

Every equation in physics is dimensionally homogeneous: both sides have the same dimension, and every term in a sum has the same dimension. This is the oldest sanity check in engineering, and it catches a surprising fraction of algebra errors before any number is computed. In P00.12 you will use it on whole designs. Here you make a computer do it on every line.

### 2.2 The pound-force, derived

The international pound is *defined* as exactly 0.453 592 37 kg. Standard gravity `g₀` is defined as exactly 9.806 65 m/s². The pound-force is the weight of one pound under standard gravity:

```
1 lbf = (0.453 592 37 kg)(9.806 65 m/s²) = 4.448 221 615 260 5 N   (exact, to the digits shown)
```

So one pound-force-second is 4.448 221 6 N·s. There is no measurement in this number; it is a product of two definitions. Any software that converts between these units and gets a different factor is wrong, not "approximately right". The tests in this lesson check your factor against the product of the two definitions to one part in 10¹².

### 2.3 Impulse, velocity change, and why small errors accumulate

A thruster producing force `F` for time `Δt` delivers impulse `I = F Δt`, dimension `[force][time]`. Newton's second law integrated over the burn gives the change in momentum, so for a spacecraft of mass `m` whose mass barely changes during a short firing,

```
Δv = I / m        units: (N·s)/kg = (kg·m/s²·s)/kg = m/s.
```

The MCO's attitude control used small thrusters to dump angular momentum from its reaction wheels (the Board calls these Angular Momentum Desaturation events). Each firing produced a small, unintended translational impulse along with the intended torque. The navigation team modelled those impulses to predict the trajectory. If each firing's impulse is under-estimated by the factor `k = 4.448`, then the modelled `Δv` is `1/k` of the true one and the *unmodelled* part is

```
Δv_missed = (1 − 1/k) · Δv_true = 0.775 · Δv_true    per firing.
```

Each individual `Δv` was tiny, millimetres per second. But position error grows with time: a velocity error `δv` sustained over time `T` becomes a position error of order `δv · T`, and over months of cruise with firings every few hours the sum crossed from "navigation noise" to "below the atmosphere". The report notes that the spacecraft's single solar array made it require far more desaturation firings than a symmetric design would have, which multiplied the effect.

### 2.4 Where the unit got lost

Here is the interface, schematically:

![Three boxes in a row: the spacecraft team’s thruster model, computing in lbf·s, writes the number 0.42 into the Small Forces file, meaning lbf·s; the navigation team’s trajectory model reads the same 0.42 as N·s.](figures/mco-interface.svg)
*Figure 1. The Mars Climate Orbiter interface, schematically. The number crossed the interface intact; the unit did not, so every impulse was read 4.448 times too small.*

The file carried numbers. The *unit* lived in two places that disagreed: in the spacecraft team's code (which computed in lbf·s) and in the interface specification (which said N·s). Nobody compared them, and the data format gave no way to notice. The Board's contributing causes include that the interface specification was not followed and that navigation anomalies noticed during cruise were not escalated and resolved. The engineering lesson is structural: whenever a unit is implied rather than carried, you have an interface that can only be checked by a human reading two documents at once.

### 2.5 Carrying the unit in code

`pint` attaches a unit to a number and tracks it through arithmetic:

```python
from pint import UnitRegistry
ureg = UnitRegistry()
Q_ = ureg.Quantity

I = Q_(0.42, "lbf * s")          # a quantity, not a number
I.to("N * s")                    # 1.868 N·s
Q_(1, "m") + Q_(1, "s")          # raises DimensionalityError
```

Three habits make this work in practice. One registry per program, created once and imported everywhere, because quantities from different registries cannot be combined. Convert to a plain number only at the boundary, with the unit named in the call: `magnitude_in(q, "m/s")`. And check dimensions at function entry, so that `dynamic_pressure(v, rho)` with the arguments swapped fails immediately rather than returning a number with the unit `kg²/(m⁴·s²)`.

### 2.6 Worked example: a cruise condition

An airliner at 11 km, where the ISA density is 0.3639 kg/m³, cruises at 250 knots true airspeed. Dynamic pressure `q = ½ ρ V²`. By hand: 250 kt × (1852 m)/(3600 s) = 128.6 m/s; `q = 0.5 × 0.3639 × 128.6² = 3 010 Pa`. Units: `(kg/m³)(m²/s²) = kg/(m·s²) = N/m² = Pa`. With `pint`, the knot-to-metre-per-second step disappears from your code and cannot be done wrongly:

```python
q = (0.5 * Q_(0.3639, "kg/m**3") * Q_(250, "knot")**2).to("Pa")   # 3009.6 Pa
```

Had you typed `Q_(250, "km/h")` by mistake, the answer would be 877 Pa and nothing would stop you; `pint` guards dimensions, not intentions. That is the limit of the tool and the reason P00.12 exists.

### 2.7 Three ways to carry a unit across an interface

| Where the unit lives | Example | What catches a mismatch | What does not |
|---|---|---|---|
| In documentation only | `impulse: float` with "N·s" in a specification document | A person reading two documents at once | Everything automatic; this is the MCO case |
| In the name | `impulse_N_s: float`, a column header `impulse [N*s]` | A reviewer reading the code; a parser that *checks* the header | Arithmetic: `impulse_N_s + thrust_lbf` is legal Python |
| In the data | `Q_(0.42, "lbf*s")`; a file row `0.42 lbf*s` | The library, on every operation | Value errors and wrong-but-dimensionally-consistent quantities |

The third row is what this lesson builds. It is not free: `pint` quantities are slower than NumPy arrays by a large factor in tight loops, and some libraries will not accept them. The engineering compromise, used throughout this curriculum, is to carry units through the *interface* (file headers, function arguments, return values) and strip them, by name, only inside numerical kernels: `magnitude_in(q, "N*s")` at the door, plain arrays inside, a `Quantity` on the way out.

A file format that would have made the MCO error impossible needs one line:

```
# impulse [lbf*s]
0.42
0.37
```

A reader that parses the header into `Q_(value, header_unit)` and converts to the unit *it* expects will either get the right number or raise. A reader that ignores the header is the 1999 reader.

### 2.8 A unit trap that is not an error: specific impulse

Rocket motors are rated by *specific impulse*, impulse per unit of propellant. Two conventions coexist. Per unit propellant *mass*, `I_sp = I / m_p` has units N·s/kg = m/s and is the exhaust velocity. Per unit propellant *weight* (under standard gravity), `I_sp = I / (m_p g₀)` has units of seconds. A motor with 2 500 m/s exhaust velocity has an `I_sp` of 255 s. Both are legitimate; data sheets use the second because the number is the same in SI and in pound units (lbf·s per lb is also seconds, since `1 lbf·s / 1 lb = g₀ · 1 s`). A dimension check alone cannot tell you which convention a file uses, because seconds and metres per second are both valid units of "specific impulse" under *different definitions*. Here the unit is not enough; you must also carry the definition. P00.12 returns to this; for now, when you see `I_sp` in seconds, multiply by `g₀` before using it as a velocity.

## 3. Build the computational model

Open `starter/main.py`. Implement the functions in the order below; `python main.py` prints the report once they are done.

### Step 1: one registry, one `Quantity` constructor

```python
ureg = UnitRegistry()
Q_ = ureg.Quantity
```

The tests build their quantities with `impl.ureg.Quantity`, i.e. with *your* registry, which is why it must be a module-level name.

### Step 2: boundaries

```python
def magnitude_in(q, unit):
    return float(q.to(unit).magnitude)

def round_trip(value, unit_a, unit_b):
    return magnitude_in(Q_(value, unit_a).to(unit_b), unit_a)
```

`magnitude_in` is the only place a `Quantity` becomes a bare number. Asking for a pressure "in metres" raises `DimensionalityError` inside `.to`, which is what you want. `round_trip` is the test of a conversion table: `a → b → a` must return the starting value to within floating-point rounding (about 10⁻¹² relative for these magnitudes).

### Step 3: addition that refuses nonsense

```python
def add_checked(a, b):
    if not isinstance(a, Quantity) or not isinstance(b, Quantity):
        raise TypeError("plain numbers carry no units")
    return a + b
```

`pint` already raises on `1 m + 1 s`. The function adds one thing: it refuses a *bare float*, because `Q_(1, "m") + 1.0` is the way units leak back out of a program (pint treats a bare number as dimensionless and raises, but the error message is confusing; a `TypeError` at the interface is clearer).

### Step 4: formulas with dimension checks at the door

```python
def dynamic_pressure(rho, v):
    if not rho.check("[mass] / [length] ** 3"):
        raise DimensionalityError(rho.units, "kg / m ** 3")
    if not v.check("[length] / [time]"):
        raise DimensionalityError(v.units, "m / s")
    return (0.5 * rho * v**2).to("Pa")
```

`Quantity.check` compares dimensions, not units, so slugs per cubic foot and knots pass. `total_impulse(F, Δt)` and `delta_v(I, m)` follow the same pattern. The final `.to("Pa")` is a second check: if the arithmetic produced anything that is not a pressure, the conversion fails.

### Step 5: the Mars Climate Orbiter functions

```python
def lbf_s_to_N_s_factor():
    return magnitude_in(Q_(1.0, "lbf * s"), "N * s")

def mco_misread_impulse(x):
    true = Q_(x, "lbf * s").to("N * s")     # what the file meant
    misread = Q_(x, "N * s")                # what the reader attached to the bare number
    return true, misread
```

`accumulated_delta_v_error(impulses, mass)` sums the firings, takes the difference between the true and misread totals, and divides by the mass with `delta_v`. The ratio `true/misread` is the same for every firing: the error is a pure scale factor, which is exactly why it was invisible in any single number.

### Verify

- *Analytic case:* `lbf_s_to_N_s_factor()` equals `0.45359237 × 9.80665` to 10⁻¹² relative; `dynamic_pressure(1.225 kg/m³, 100 m/s)` equals `6125 Pa`.
- *Round trips:* ft↔m, knot↔m/s, psi↔Pa, lbf↔N, slug/ft³↔kg/m³ all return to the start.
- *Conservation under change of units:* `dynamic_pressure` with inputs in slug/ft³ and knots gives the same pascals as with SI inputs (physics does not depend on the unit system).
- *Negative cases:* `add_checked(1 m, 1 s)` raises `DimensionalityError`; `dynamic_pressure(v, rho)` with arguments swapped raises; `delta_v(force, mass)` raises because a force is not an impulse; `parse_quantity("3 furlongs_per_fortnite")` raises `UndefinedUnitError`.
- *Linearity:* the accumulated error over ten identical firings is ten times that of one; doubling the mass halves it.

### Compare

Do the cruise-condition example once by hand (part 2.6), once with `pint`, and once with "bare NumPy" where you convert knots yourself. All three agree to four figures. Now introduce a deliberate error in each: in the hand calculation forget to square the speed; in the NumPy version mistype the knot→m/s factor as 0.5114 instead of 1852/3600 = 0.5144; in the `pint` version pass `Q_(250, "knot")` as `rho`. Only the third one stops the program. Units libraries catch *dimension* errors loudly and *value* errors not at all; hand checks catch value errors and are slow. You need both.

## 4. Experiment at home, or investigate real data

**Mode: data-based investigation.** Read the Mars Climate Orbiter Mishap Investigation Board Phase I Report (`@nasa-mco-mib`, a public NASA document). Before reading, write down your prediction: in a chain of events from "a file is written in the wrong unit" to "the spacecraft is lost", how many independent opportunities to catch the error do you expect the Board to identify? Then read, and build the chain yourself as a numbered list with one sentence per link, marking each link as *root cause*, *contributing cause*, or *missed opportunity*, using the Board's own classification where it gives one.

Then run your code with the report in front of you. Pick an illustrative firing size (the report does not give individual impulses; 0.2–0.6 lbf·s per event is a reasonable order of magnitude for small attitude thrusters, and you should say so), an illustrative firing count for nine months of cruise, and an illustrative spacecraft mass of a few hundred kilograms. Compute `accumulated_delta_v_error`. You will get a number of order 1 m/s. That number is *your* model's output, not a mission figure: the real trajectory error depended on when each firing happened and in which direction, which a scalar sum ignores.

Raw-data format for your write-up: a table `link_number, event, classification, page_reference`.

## 5. Validate and improve

Compare your causal chain with the Board's list of contributing causes. Where you have a link the Board does not, ask whether it is a real cause or hindsight. Where the Board has one you missed, that is the useful finding: it is usually an *organisational* link (staffing, training, communication between teams) rather than a technical one. Unit errors are technically trivial and organisationally hard; the report's recommendations are almost all about process.

Compare the three cruise-condition calculations from part 3. If the `pint` result and your hand result differ in the fourth figure, find the rounding: the knot is exactly 1852 m/h, and `pint` uses that definition; a hand value of 0.514 m/s per knot is 0.05 % low.

Improve: make `accumulated_delta_v_error` accept a direction for each firing (a unit vector) and sum impulses as vectors. Observe that the missed `Δv` can now partially cancel, and ask what that does to your confidence in the scalar estimate.

## 6. Defend your work

Submit: `main.py` and its test output; the causal chain table; the three-way cruise comparison with the deliberate errors and which were caught; and a half-page answer to the engineering question in the metadata table. Use `projects/_template/report-template.md`.

### Adversarial self-checks

- What if the firing count doubles? The missed `Δv` doubles; the position error after a fixed cruise time doubles. What if the mass doubles? Halves. Which of these did the MCO design actually change by having one solar array instead of two?
- Which assumption fails first in the real world? Constant mass during a firing is fine; *scalar* summing is not, because firings in different directions partly cancel. Does that make the real error larger or smaller than your estimate? (It depends on the firing geometry; the report indicates the effect was systematic, not random.)
- What measurement could falsify the "unit error" explanation? A navigation residual history that did *not* show a consistent bias in the direction of the modelled small forces. The Board found residuals that did show it, and that were noticed but not resolved.
- Does the result survive a change of `pint` version or registry definitions? The lbf→N factor is a definition and will not change; the knot was redefined historically, and older registries may differ in the fourth figure. Pin the factor in a test, as this lesson does.
- Which decision changes if the uncertainty in each firing's impulse were twice as large? None, if the uncertainty is random; a factor-4.45 bias is not uncertainty, and no error bar would have contained it. This is the difference between precision and correctness.

## Check yourself

1. State the dimension of impulse in base dimensions, and give two different units for it.
2. Why is 1 lbf = 4.448 221 615 N exact rather than measured?
3. `Q_(20, "degC") + Q_(5, "degC")` raises an error in `pint`. Why is that the correct behaviour, and how do you express "5 °C warmer"?
4. A function takes `altitude_ft: float`. Name the two places the unit now lives and the failure that follows when they disagree.
5. By what factor were the MCO small-forces impulses under-estimated, and in which direction was the resulting trajectory error (too high or too low at Mars)?
6. `pint` accepted `dynamic_pressure(Q_(0.36, "kg/m**3"), Q_(250, "km/h"))` and returned 877 Pa. Was that an error `pint` should have caught? What would catch it?
7. Give one quantity for which a dimension check passes but the physics is still wrong, from aerospace practice.
8. Why does `round_trip` test to 10⁻¹² and not to exact equality?
