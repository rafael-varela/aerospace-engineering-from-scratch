# P00.04  Worked solutions

Each solution states: assumptions, units, the derivation, the result with uncertainty where applicable, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes.

## 1. Reading a label

**Assumptions.** The figures are from an aircraft or motor test, so "k" is most likely kilofeet or kilometres, "CL" the lift coefficient, and "lbs" pounds-force mislabelled as pounds-mass.

**Result.**

| Original | House style | What the reader could get wrong |
|---|---|---|
| `Altitude (k)` | `Altitude [km]` or `Altitude [kft]` | "k" is a prefix, not a unit. 35 k could be 35 km (above any airliner) or 35 kft (≈ 10.7 km, a normal cruise altitude); the two differ by a factor 3.3. |
| `CL` | `Lift coefficient [-]` | Without the bracket the reader cannot tell whether a unit was forgotten. `CL` is also used for "centreline" and "confidence level" in other documents. |
| `Thrust, lbs` | `Thrust [lbf]`, or convert and write `Thrust [N]` | "lbs" is pounds-mass, a unit of mass; thrust is a force. A reader converting with 0.4536 kg/lb gets kilograms, not newtons; the correct factor to newtons is 4.448 N/lbf (see P00.05). |

**Check.** Each corrected label has a quantity, square brackets, and a unit symbol that exists in SI or has an exact SI definition.

**Common mistakes.** Writing `Altitude [k]` (copying the error into the new format); choosing `[kg]` for thrust because "the number is the same on Earth" (it is not: 1 lbf = 4.448 N, 1 lb = 0.4536 kg, and the plot needs one of them).

## 2. The straight line on a log axis

**Assumptions.** Hydrostatic equilibrium, ideal gas, constant `g`, constant `T` in the layer, `h` is geopotential altitude.

**Derivation.** `dp/dh = −ρ g` and `ρ = p/(R T)` give `dp/p = −(g/(R T)) dh`. Integrating from `(h₀, p₀)`:

```
ln p − ln p₀ = −(g/(R T)) (h − h₀)      ⇒      p = p₀ exp[−(g/(R T))(h − h₀)].
```

So `ln p` is linear in `h` with slope `−g/(R T)`. Units: `g/(R T)` has units `(m/s²) / ((J/(kg·K))·K) = (m/s²)/(m²/s²) = 1/m`. Good.

**Result.** `g/(R T) = 9.80665 / (287.05 × 216.65) = 1.577 × 10⁻⁴ m⁻¹ = 0.1577 km⁻¹`. In decades per kilometre, divide by `ln 10 = 2.3026`: **0.0685 decades/km**, i.e. pressure falls by a factor `e^{−0.1577} = 0.854` for every kilometre climbed; a decade takes 14.6 km.

Troposphere: `T = T₀ + L h` is not constant, so `dp/p = −(g/(R (T₀ + L h))) dh` integrates to `ln p = ln p₀ + (g/(R L)) ln(T/T₀)`, a logarithm of a linear function, **not** a straight line in `h`. On the log plot the troposphere curves (concave down for `L < 0`), and the curvature carries the lapse rate.

**Check.** Limiting case `L → 0` in the troposphere formula: `(g/(R L)) ln(1 + L h/T₀) → g h/(R T₀)`, recovering the isothermal straight line. Numerically, 11 km to 20 km in the isothermal layer: `22 632 × 0.854⁹ = 5 475 Pa`, which matches the ISA table.

**Common mistakes.** Forgetting to divide by `ln 10` when reporting decades; using `L = +6.5 K/km`; treating the troposphere curve as a plotting error and "fixing" it with a different axis.

## 3. Choosing axes for a drag measurement

**Assumptions.** `D = ½ ρ V² S CD`, `ρ = 1.20 kg/m³`, `S = 0.050 m²`, `CD` constant, force uncertainty ±5 % (1 σ) and independent between points.

**Derivation.** `log D = log(½ ρ S CD) + 2 log V`. On log–log axes constant `CD` is a straight line of **slope 2**. Semi-log would give a parabola; linear–linear squeezes the 5 m/s point against the axis (ratio of forces 254).

**Result.** `CD = D / (½ ρ V² S) = D / (0.030 V²)` with `V` in m/s and `D` in N:

| V [m/s] | D [N] | CD [-] |
|---|---|---|
| 5 | 0.12 | 0.160 |
| 10 | 0.47 | 0.157 |
| 20 | 1.9 | 0.158 |
| 40 | 7.6 | 0.158 |
| 80 | 30.5 | 0.159 |

Mean `CD = 0.1584`, sample standard deviation 0.0012, standard error of the mean 0.0005. The point-to-point scatter (0.8 %) is much smaller than the stated ±5 % on each force, so the data are consistent with constant `CD`: **CD = 0.158 ± 0.008** (the ±5 % of a single point dominates; if the 5 % is a *systematic* balance error it does not average down). On the figure, all five bars overlap the slope-2 line.

**Check.** Dimensional: `N / (kg/m³ · m²/s² · m²) = N / (kg·m/s²) = 1`. Slope check from the end points: `log₁₀(30.5/0.12) / log₁₀(80/5) = 2.405/1.204 = 1.998`.

**Common mistakes.** Fitting a straight line to `D` against `V` and declaring the relation non-linear; quoting `CD` to five figures when the balance is good to 5 %; plotting `CD` against `V` on a log axis "because the lesson said so" (the ratio of `CD` values is 1.02, far below 30).

## 4. The checker versus the wild

**Result.** With `seed = 11` (any seed gives the same catch/miss table):

| Flaw | Caught | Message |
|---|---|---|
| no y label | yes | `axes 0: missing y label` |
| label `Thrust (N)` | yes | `... has no [unit] tag` |
| two labelled series, no legend | yes | `... 2 labelled series but no legend` |
| line in `"magenta"` | yes | `... line color 'magenta' is not a house color` |
| measured points joined by a line | **no** | (none) |

The checker cannot know a series was *measured*; a `Line2D` stores `linestyle`, `marker` and data, nothing about provenance. A checkable proxy: "a series drawn with error bars (an `ErrorbarContainer`) must have its connecting line hidden". Pseudocode:

```
for c in ax.containers:
    if isinstance(c, ErrorbarContainer):
        line = c.lines[0]            # the data line of the container
        if line.get_linestyle() not in ("None", "none", ""):
            problems.append("measured series joined by a line")
```

`plot_series` already satisfies it (`fmt="o"`), so the rule costs compliant figures nothing.

**Check.** Build the compliant thrust figure and run the new rule: no message. Build `ax.errorbar(t, f, yerr=s, fmt="-o")`: one message.

**Common mistakes.** Testing the rule by looking at the picture instead of asserting on the message list; forgetting that `ax.get_lines()` does *not* include the error-bar container's data line, so the colour rule must also loop over `ax.containers` if colours in error-bar series are to be policed.

## 5. Reproducibility of noise

**Assumptions.** Bars are ±1 σ; noise is Gaussian with exactly that σ, by construction of `synthetic_thrust_curve`.

**Derivation.** For a standard normal `Z`, `P(|Z| > 1) = 0.317`. Expected count of bars that miss the model: `25 × 0.317 = 7.9`. The count is binomial with standard deviation `√(25 × 0.317 × 0.683) = 2.3`.

**Result.** Example counts for seeds 6 and 7: 12 and 6. Both lie within about two binomial standard deviations of 7.9 (12 is 1.8 σ high; with hundreds of learners each reporting two seeds, counts of 12 or 3 will appear routinely). Your own seeds will differ; the expected count does not.

**Check.** Run 2000 seeds and average the count: 7.86 ± 0.05, consistent with `25 × 0.317 = 7.93`.

**Common mistakes.** Expecting every bar to touch the model ("the bars are the uncertainty, so the model should be inside them"); concluding that seed 6 is "bad data" because 12 bars miss. A ±1 σ bar is *supposed* to miss about a third of the time. If no bars missed, the σ would be overstated.

## 6. When the house style is wrong

**Result (two of many acceptable answers).**

- *Raw telemetry.* A 1 kHz accelerometer log of a 300 s flight has 300 000 measured points. Markers with bars become an opaque smear; the reader wants to see trends, so a thin line (no markers) with a shaded ±σ band is the honest choice. The rule "no lines through measured points" protects the reader from invented data between sparse points; with dense data nothing is invented.
- *Quantities through zero.* Pitching-moment coefficient against angle of attack crosses zero; a log axis is impossible and symmetric-log (`symlog`) hides the crossing's slope, which is the stability derivative the reader cares about. Use a linear axis even though the ratio of extremes is large.
- (*Also acceptable:* a slide for a projector where a 300 dpi PNG displays faster and no one will zoom; a colour-blind reviewer who asks for line styles instead of colours.)

**Observable sign of misapplication.** A log axis whose data occupy less than a decade (all points fall between two adjacent major ticks), or a vector SVG that is tens of megabytes because it contains a hundred thousand markers.

**Common mistakes.** Treating the style as a law rather than as a default; breaking a rule without writing down which reading you traded away.

## 7. A checker rule with a false positive

**Result.** The colleague is right: `Time since launch [s] (UTC 12:00:00 origin)` marks the unit unambiguously; the checker's "ends in `]`" test is a proxy that failed. Replace it by a regular expression that requires a bracketed, non-empty, bracket-free unit anywhere in the label:

```python
import re
UNIT_TAG = re.compile(r"\[[^\[\]]+\]")

def _has_unit_tag(label: str) -> bool:
    return UNIT_TAG.search(label) is not None
```

`Thrust (N)` still fails (round brackets); `Thrust []` fails (empty); `Thrust [N] (sea level)` passes. Test:

```python
def test_unit_tag_anywhere():
    assert _has_unit_tag("Time since launch [s] (UTC 12:00:00 origin)")
    assert not _has_unit_tag("Thrust (N)")
    assert not _has_unit_tag("Thrust []")
```

**Argument.** A checker is a filter on the author's attention. Each false positive teaches the author to skim its output; once the author has learned to ignore it, the checker catches nothing at all. In practice a rule that is wrong more than about one time in ten gets silenced, so rules should be slightly permissive (accept the odd ugly-but-unambiguous label) rather than slightly strict. A missed defect costs one review comment; a silenced checker costs every defect thereafter.

**Common mistakes.** Making the regex accept `(N)` "to be helpful" (two conventions in one report is worse than one); writing the test only for the passing case.

## Check-yourself answers

1. `pyplot` picks a backend when it is first imported; `matplotlib.use("Agg")` after that import is ignored or warns. Selecting first guarantees no window is ever opened under tests or on a server.
2. The unit: m/s, km/h, knots or ft/s (a factor up to 3.3 between them); and the frame: true airspeed, indicated airspeed, or ground speed, which differ by tens of percent at altitude.
3. One decade spacing is a factor 10; with decades at 10⁴ and 10⁵ Pa the 20 km value (5.47 kPa) sits 0.26 decades below 10⁴ and can be read to a few percent. On a linear axis spanning 0–101 kPa the same value is 5 % of the axis height; a ±5 % reading would require resolving 0.27 kPa, about 0.3 % of the axis, finer than a line width.
4. A line asserts a value at every point between two measurements. Nothing was measured there; thrust could have spiked between samples. Markers say "we know these points", bars say "this well", and the model line makes the between-point claim explicitly.
5. The author never called `ax.legend()` (the `label=` arguments only *register* names). `finish_axes(ax)` or `ax.legend()` fixes it.
6. PNG stores pixels: zooming blurs text and lines, and a 300 dpi page-width figure is several megabytes. SVG/PDF store geometry and stay sharp at any size. Converting to PNG is legitimate as the *last* step for a venue that cannot embed vector graphics (a web page, some slide tools), keeping the vector original.
7. Loses: the parabola's shape, which is the physical content (`CD0` as the vertex, `k` as the curvature) and the ability to read `CD0` off the axis. Might gain: visibility of the small `CD` at low `CL` if `CD0` were tiny compared with the induced term, which is rare (ratio ≈ 7 here).
8. No: the colour rule only checks membership in the palette. A rule that would: within one `Axes`, no two labelled series may share a colour unless they also differ in line style; implement by collecting `(color, linestyle)` pairs over `ax.get_lines()` and flagging duplicates.
