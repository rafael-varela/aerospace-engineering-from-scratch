# P00.04  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every axis and every number you report.

## 1. (Warm-up: reading a label)

The three labels below appear on figures in a flight-test report. For each, write the house-style replacement and say what a reader could get wrong with the original.

(a) `Altitude (k)`  (b) `CL`  (c) `Thrust, lbs`

Answer form: three corrected labels and one sentence each.

## 2. (Derivation: the straight line on a log axis)

Starting from hydrostatic balance `dp/dh = −ρ g` and the ideal gas law `p = ρ R T`, derive `p(h)` for an isothermal layer at temperature `T` beginning at `(h₀, p₀)`. Show that `ln p` is linear in `h` and give the slope in units of "decades per kilometre" for `T = 216.65 K`, `g = 9.80665 m/s²`, `R = 287.05 J/(kg·K)`. Then state, with a one-line argument, whether the troposphere (where `T = T₀ + L h`, `L = −6.5 K/km`) gives a straight line on the same axes.

Answer form: derivation, a number with units, and a yes/no with reason.

## 3. (Unfamiliar scenario: choosing axes)

A wind-tunnel run records drag force on a model at airspeeds 5, 10, 20, 40 and 80 m/s: 0.12, 0.47, 1.9, 7.6 and 30.5 N, each ±5 %. You expect `D = ½ ρ V² S CD` with constant `CD`.

(a) Which axis scaling (linear–linear, log–log, or semi-log) makes "constant `CD`" appear as a straight line, and what slope do you expect?
(b) Plot the data that way with `plot_series` and error bars, add the expected line, and decide from the figure whether `CD` is constant across the range. Report `CD` with `ρ = 1.20 kg/m³`, `S = 0.050 m²`.

Answer form: a figure (SVG) built with your module, a slope, and `CD` with an uncertainty estimated from the scatter.

## 4. (Numerical, uses your code: the checker versus the wild)

Write a function in a scratch file that builds five deliberately flawed figures using plain `plt` calls: missing y label, unit-less label, two labelled series without a legend, a line in an off-palette colour, and a measured series joined by a line. Run `check_figure` on each. Report which flaws it catches and which it cannot, and propose a rule for one it misses. Parameterize: use a seed of your choice for `synthetic_thrust_curve` and report the seed.

Answer form: a table of flaw → caught / not caught → message, plus one proposed rule written as a sentence and as pseudocode.

## 5. (Numerical, uses your code: reproducibility of noise)

Call `synthetic_thrust_curve(seed=s, n=25)` for `s` = your birth month and for `s + 1`. Plot both on one Axes with error bars, each in a *different* role (`"measured"` and `"reference"`), with the generating model in `"model"`. For each seed count how many of the 25 points have error bars that do not reach the model curve. Compare with the number you expect for ±1 σ bars, and explain the difference in two sentences.

Answer form: a figure, two counts, an expected count, and the explanation.

## 6. (Judgment: when the house style is wrong)

Give two situations in aerospace work where one of the house rules (log axis for ratios > 30, no lines through measured points, vector output only, role-based colours) should be broken, and state what you would do instead and why. Then name one measurable sign on a finished figure that the style has been misapplied.

Answer form: two short arguments and one observable sign.

## 7. (Hard: a checker rule with a false positive)

`check_figure` requires every label to end in `[unit]`. A colleague labels a histogram's y axis `Count [-]` and a polar plot's angular axis `Angle of attack [deg]`, and the checker passes both. A third colleague labels a time axis `Time since launch [s] (UTC 12:00:00 origin)` and the checker fails it. Decide whether the checker or the colleague is right, modify the rule so it accepts the third label without accepting `Thrust (N)`, and write the test that proves the change. Explain in one paragraph what makes a checker rule good: what fraction of false positives would make you turn it off?

Answer form: modified code, a passing test, a paragraph.
