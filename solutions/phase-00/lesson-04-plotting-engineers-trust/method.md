# P00.04  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. Identify the quantity each label is trying to name. 2. Identify the unit, if any, and whether the written unit is a real unit. 3. Rewrite as `Quantity [unit]`, using `[-]` for dimensionless. 4. For each original, name the specific misreading (wrong unit, ambiguous symbol, unit of mass used for force).

## 2.
1. Write `dp/dh = −ρ g` and `ρ = p/(R T)`. 2. Combine to `dp/p = −(g/(R T)) dh`. 3. Integrate from `(h₀, p₀)` with `T` constant. 4. Take `log₁₀` and read off the slope; convert from natural log to decades by dividing by `ln 10`; convert from per metre to per kilometre. 5. For the troposphere, repeat step 2 with `T(h)`; notice the integrand is no longer constant, so `ln p` is not linear in `h`.

## 3.
1. Take logs of `D = ½ ρ V² S CD`. 2. With constant `CD`, `log D` is linear in `log V` with a fixed slope: this selects log–log axes. 3. Build the figure with `new_figure(..., "Airspeed", "m/s", "Drag", "N")`, `plot_series` with `yerr = 0.05 D` in the `"measured"` role, `set_log_axes(ax, x=True, y=True)`, and a `"model"` line `½ ρ V² S CD_mean`. 4. Compute `CD` at each point as `D / (½ ρ V² S)`; take the mean and the standard deviation; compare the spread with the ±5 % bars. 5. Judge constancy by whether every bar touches the model line.

## 4.
1. Build each flawed figure with `plt.subplots()` and plain calls. 2. Run `check_figure` and record the messages. 3. For the flaw with no message, explain which object property would need inspecting (a `Line2D` has a linestyle and markers; whether its data were "measured" is not stored anywhere). 4. Propose a rule that *is* checkable, e.g. "a series with error bars must have `linestyle='none'`", and write it as pseudocode over `ax.containers`.

## 5.
1. Generate the two seeded datasets. 2. Draw both with error bars in different roles and the model as a line. 3. For each point test `|F − model| > σ`; count. 4. Compare each count with `25 × P(|Z| > 1)` for a standard normal `Z`. 5. Explain the gap with the binomial standard deviation of a count of 25 trials.

## 6.
1. For each house rule write down the reading it protects. 2. Find a situation where that reading is irrelevant or impossible (e.g. a time series of raw telemetry with thousands of points; a figure for a slide shown on a projector; a quantity that passes through zero). 3. State the replacement rule and its cost. 4. Pick an observable sign: something a reader can point to on the finished figure.

## 7.
1. Decide what the rule is *for* (unambiguous unit marking). 2. Replace "ends with `]`" by a regular expression that finds `[...]` with a non-empty, bracket-free unit anywhere in the label. 3. Confirm `Thrust (N)` still fails (round brackets are not the house convention). 4. Write a pytest that passes the long label and fails `Thrust (N)`. 5. Argue the trade-off: a checker that people silence is worth nothing, so the false-positive rate must stay below what a tired author will tolerate.
