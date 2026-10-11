# P00.05  Worked solutions

Each solution states: assumptions, units, the derivation, the result with uncertainty where applicable, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes.

## 1. Four conversions, twice

**Assumptions.** International foot (0.3048 m exactly), international nautical mile (1852 m exactly), conventional inch of mercury (3386.389 Pa, the value `pint` defines), slug = 1 lbf·s²/ft.

**Derivation and result.**

| | Hand chain | Hand value | `pint` | Rel. diff. |
|---|---|---|---|---|
| (a) | 35 000 ft × 0.3048 m/ft | 10 668 m (exact) | 10 668.0 m | 0 |
| (b) | 250 nmi/h × 1852 m/nmi × 1 h/3600 s | 128.61 m/s | 128.611 m/s | < 10⁻⁴ (rounding) |
| (c) | 29.92 inHg × 3386.39 Pa/inHg | 101 321 Pa | 101 320.7 Pa | 3 × 10⁻⁶ |
| (d) | 14.5939 kg/slug ÷ (0.3048 m/ft)³ = 14.5939 / 0.0283168 | 515.38 kg/m³ | 515.379 kg/m³ | < 10⁻⁵ |

**Check.** (c) 29.92 inHg is the standard sea-level pressure in aviation altimeters and should come out within 0.01 % of 101 325 Pa: it does (the difference is the rounding of 29.92). (d) Sea-level ISA density is 0.002377 slug/ft³ = 1.225 kg/m³, so 1 slug/ft³ = 1.225/0.002377 = 515.4 kg/m³. Agrees.

**Common mistakes.** Using the US survey foot (differs in the sixth figure); forgetting to cube the length factor in (d); using 1 knot = 0.5 m/s "roughly" in a calculation that is later quoted to four figures.

## 2. The pound-force and the velocity change

**(a)** `1 lbf = 0.453 592 37 kg × 9.806 65 m/s² = 4.448 221 615 260 5 N`. To ten significant figures, **4.448 221 615 N**. Both factors are *definitions* agreed by international convention (the 1959 international yard and pound agreement; the standard-gravity value adopted by the CGPM in 1901), so the product is exact and no uncertainty attaches to it. Local `g` varies by ±0.3 %, but `lbf` is defined with `g₀`, not with local gravity.

**(b)** Newton's second law for a body of momentum `p = m v`: `F = dp/dt`. Over a firing from `t₁` to `t₂`,

```
∫ F dt = ∫ d(m v) = m₂ v₂ − m₁ v₁.
```

Assume `m` constant during the firing (`m₁ = m₂ = m`): then `I = m Δv`, i.e. `Δv = I / m`. Units: `(N·s)/kg = (kg·m/s²·s)/kg = m/s`. The assumption fails by more than 1 % when the propellant expelled during the firing exceeds about 1 % of the vehicle mass; for an attitude thruster delivering 0.4 lbf·s = 1.8 N·s at an exhaust speed of ~2 000 m/s, the propellant used is 0.9 mg against hundreds of kilograms, so the assumption is exact to one part in 10⁸. For a main-engine burn consuming tens of percent of the mass, use the rocket equation `Δv = v_e ln(m₁/m₂)` instead.

**(c)** `lbf·s / lb = (lb · g₀ · s) / lb = g₀ · s = 9.806 65 m/s`, a velocity: so 1 lbf·s per lb of propellant is an exhaust velocity of 9.807 m/s, and `I_sp` quoted in "seconds" is the exhaust velocity divided by `g₀`. It is meaningful because impulse *per unit weight* is a legitimate ratio, and its numerical value is the same whether impulse is in lbf·s and propellant in lb, or impulse in N·s and propellant in kgf: the unit system drops out. The price is that `I_sp[s] × g₀` must be done before using it as a velocity.

**Check.** A motor with `I_sp = 255 s` has `v_e = 255 × 9.807 = 2 501 m/s`, typical of a solid propellant. Dimensional check on the rocket equation: `v_e × ln(ratio)` is a velocity. Good.

**Common mistakes.** Using local `g` (9.81 or 9.78) in the lbf definition; dropping `g₀` and treating 255 s as 255 m/s (a factor of ten error that survives a dimension check only if someone also mislabels the unit).

## 3. The fuel order

**Assumptions.** Density exactly 6.7 lb/gal at the delivery temperature; 1 lb = 0.453 592 37 kg; US gallon = 3.785 411 784 L; imperial gallon = 4.546 09 L.

**Derivation.** `m = ρ V`. With US gallons the gallon cancels: `6.7 lb/gal × 2 400 gal = 16 080 lb = 7 294 kg`. If the supplier delivers 2 400 imperial gallons of fuel whose density is 6.7 lb per *US* gallon: `2 400 × 4.546 09 / 3.785 411 784 = 2 882.4 US gal`, so `m = 6.7 × 2 882.4 = 19 312 lb = 8 759 kg`.

**Result.** (a) **7 294 kg**; (b) **8 759 kg**; the ratio is 1.2009, the ratio of the two gallons.

(c)
```python
rho = Q_(6.7, "lb / gallon")
V = Q_(2400, "gallon")
mass = (rho * V).to("kg")        # 7293.8 kg
# (rho * 2400.0).to("kg") raises DimensionalityError: lb/gallon is not a mass
```

(d) **No.** The ambiguity is 20 %, ten times the 2 % tolerance. A take-off computed with 7 294 kg of fuel when 8 759 kg is aboard is 1 465 kg heavy: the rotation speed, field length and climb gradient are all wrong. Density in kg/m³ (803 kg/m³ here) and volume in litres remove the ambiguity.

**Check.** 803 kg/m³ is within the normal Jet A range (775–840 kg/m³). Dimensional: `(lb/gal)(gal) = lb`.

**Common mistakes.** Converting lb to kg with 2.2 (0.4545, 0.2 % low) when the fuel-load tolerance is 2 %: acceptable here but a habit that fails elsewhere; writing `6.7 * 2400 * 0.4536` with no record of which gallon was meant, which is precisely the ambiguity `pint` cannot fix either unless the unit string says `imperial_gallon`.

## 4. The Mars Climate Orbiter re-enactment

**Assumptions (all illustrative).** Seed 1234, `n = 800` firings with impulses uniform in 0.2–0.6 lbf·s, `m = 500 kg`, missed `Δv` treated as a constant velocity error for 30 days. None of these are mission values.

**Derivation.** Total true impulse `Σx = 317.94 lbf·s = 1 414.3 N·s`; booked `317.94 N·s`; missed `1 096.3 N·s`. `Δv_missed = 1 096.3 / 500 = 2.19 m/s`. Position error over 30 days: `2.19 m/s × 2.592 × 10⁶ s = 5 680 km`.

**Result.** Missed `Δv` = **2.19 m/s**, position error ≈ **5 700 km** for these inputs.

**Comparison.** This is thirty times the ~170 km the Board reports. The inputs or the model are too crude, and it is worth saying which: (i) the real impulses per event were likely smaller than 0.2–0.6 lbf·s and not all in the critical direction; (ii) a velocity error is not sustained for 30 days unopposed: the navigation team re-estimated the orbit from tracking data and performed trajectory correction manoeuvres through cruise, which absorbed most of the accumulated mismodelling, and the fatal error was the part that accumulated after the last correction and in the direction that tracking could least observe. Your estimate is an upper bound on an unmanaged error and demonstrates the mechanism and its scale; it is not a reconstruction of the mission's trajectory.

**Check.** Dimensional: `(N·s)/kg × s = m`. Linearity: doubling `n` doubles both numbers; doubling `m` halves them (the tests check exactly this).

**Common mistakes.** Reporting the 5 700 km as if it reproduced the mishap; forgetting that the *booked* impulse is `x` N·s, not zero (the model under-estimated, it did not omit); quoting the result to four figures when every input is a round guess.

## 5. A flight log with units

**Assumptions.** ISA density at 11 km 0.3639 kg/m³; speeds are true airspeeds.

**Result.**

| TAS [knot] | TAS [m/s] | q [Pa] | q [psi] |
|---|---|---|---|
| 180 | 92.6 | 1 560 | 0.226 |
| 200 | 102.9 | 1 926 | 0.279 |
| 220 | 113.2 | 2 331 | 0.338 |
| 240 | 123.5 | 2 774 | 0.402 |
| 260 | 133.8 | 3 255 | 0.472 |

With the same numbers read as km/h the dynamic pressures are 455, 562, 680, 809, 949 Pa: every value is lower by the factor `(1.852)² = 3.430`, the square of the ratio knot/(km/h). On one figure the two curves are both parabolas; nothing on the figure says which is right. `pint` did not object because knots and km/h have the same dimension `[length]/[time]`; the unit was *carried* correctly, it was the *input* that was wrong. Only an external check (a typical airliner cruise is ~250 kt, not 250 km/h, which is a propeller trainer's speed) catches this.

**Check.** `q` at 250 kt was 3 010 Pa in the lesson; 240 and 260 kt bracket it. Dimensional: Pa throughout.

**Common mistakes.** Wrapping each speed in its own `Quantity` inside a Python loop (slow and unnecessary; one `Quantity` holds the whole array); plotting the `Quantity` directly and getting a unit-less axis label.

## 6. Where the units library stops

(a) **Caught.** `pint` raises an "ambiguous operation with offset unit" error: 20 °C + 5 °C has no physical meaning as absolute temperatures. Express the intent as `Q_(20, "degC") + Q_(5, "delta_degC")` = 25 °C.
(b) **Not caught.** Both are `[length]/[time]`. The frame (indicated/true/ground) is part of the quantity's *definition*, not its unit. Carry it in the variable name or a data field and check it in code review; P00.12's order-of-magnitude checks catch gross cases.
(c) **Caught**, but only because the two conventions differ in dimension (`[time]` vs `[length]/[time]`): `Q_(255, "s").to("m/s")` raises. If the file had used m/s but meant something else, it would not.
(d) **Not caught.** N·m and J are the same dimension `[mass][length]²/[time]²`; `pint` converts 1 N·m to 1 J without complaint. Only the name of the quantity distinguishes torque from energy; keep them in separately named variables and never let a torque flow into an energy balance without a cross product having happened first.
(e) **Caught.** `pint` treats 1.0 as dimensionless and raises `DimensionalityError` ("cannot convert from meter to dimensionless"). The lesson's `add_checked` raises a clearer `TypeError` before `pint` gets the chance.

**Common mistakes.** Believing a units library makes a calculation "unit-safe" in every sense; the honest statement is "dimension-safe".

## 7. The interface that would have saved the mission

```python
import re
import numpy as np
from pint import Quantity

HEADER = re.compile(r"^#\s*impulse\s*\[(.+)\]\s*$")

def read_small_forces(path) -> Quantity:
    lines = Path(path).read_text().splitlines()
    m = HEADER.match(lines[0]) if lines else None
    if m is None:
        raise ValueError("missing header line '# impulse [unit]'")
    unit = ureg(m.group(1))
    values = []
    for number, line in enumerate(lines[1:], start=2):
        if not line.strip():
            continue
        try:
            values.append(float(line))
        except ValueError:
            raise ValueError(f"line {number}: {line!r} is not a number") from None
    return (np.array(values) * unit).to("N * s")   # DimensionalityError if the header unit is not an impulse
```

Tests (with `tmp_path`): a file `# impulse [lbf*s]` / `1.0` / `0.5` returns `[4.448, 2.224] N·s`; a file with no header raises `ValueError`; a header `[lbf]` raises `DimensionalityError`. All three pass against the code above.

**Argument.** The interface removes the *mechanism* of the root cause: a unit mismatch between writer and reader now fails at the first read, loudly, before any trajectory is computed. It is therefore sufficient to prevent *this* error. It is not sufficient to prevent the *class* of error the Board describes, because most of the contributing causes are not about data: the navigation team's unfamiliarity with the spacecraft, the thin operations staffing, the interface specification not being verified end to end, and anomalies in the tracking residuals that were seen and discussed but not escalated to a resolution. A team with those properties will find another seam (a sign convention, a reference frame, a time standard) that no unit header covers. The header is necessary engineering hygiene; the Board's process recommendations are what address the cause.

**Common mistakes.** Returning the array in the file's unit instead of converting (the reader's expectation is the point); catching the `DimensionalityError` and substituting a default; not testing the failure paths.

## Check-yourself answers

1. `[mass][length]/[time]`. N·s and kg·m/s (also lbf·s).
2. Because both the pound (0.453 592 37 kg) and standard gravity (9.806 65 m/s²) are defined constants, not measurements; their product is therefore defined too.
3. Adding two *absolute* temperatures is physically meaningless (the result depends on the zero point chosen); a temperature *difference* is a different kind of quantity. Write the difference as `Q_(5, "delta_degC")`, or work in kelvin where the ambiguity disappears for differences.
4. In the parameter name (`_ft`) and in the caller's head (or documentation). When the caller passes metres, the function silently computes with a value 3.28 times too small; nothing raises, and the error is found only if someone notices the answer is implausible.
5. By 4.45 (the lbf→N factor). The modelled small forces were too small, so the modelled trajectory was less perturbed than the real one; the spacecraft arrived about 170 km *lower* than planned, into the atmosphere.
6. No: km/h is a valid speed unit and the dimension check passed, as it should. A value check catches it: 250 km/h = 135 kt is far below airliner cruise speed; P00.12's order-of-magnitude checks, or a plot of `q` against the aircraft's known operating envelope, would flag it.
7. Torque in N·m used as energy in J (same dimension); indicated vs true airspeed; `I_sp` in m/s vs the same number mislabelled; pressure altitude vs geometric altitude (both lengths). Any one is acceptable.
8. Floating-point conversion factors are rounded to about 16 significant figures and the two divisions do not cancel exactly, so `a → b → a` returns a number that differs from the start in the last bit or two. 10⁻¹² relative is a hundred times looser than that and still a thousand times tighter than any engineering tolerance, so it catches a wrong factor without false alarms.
