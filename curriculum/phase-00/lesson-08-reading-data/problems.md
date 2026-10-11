# P00.08  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units and uncertainty on every numerical answer.

Several problems use your personal seed **s**: the day of the month you were born (1–31). Pass it as `seed=s` to `generate_flight_log`.

## 1. (Warm-up: the header)

A log header reads `time [ms],altitude [ft],airspeed [kt],pitch [rad]` and its first data row is `1500,3500,120,0.05`.

(a) Convert each value to the lesson-standard unit by hand, stating the factor and whether it is exact or rounded.
(b) Suppose a careless parser ignored the units and treated the row as `s, m, m/s, deg`. For each channel give the error as a factor and say which error would be the first to be noticed by a pilot reading the result, and which could survive to a report.
Answer form: four converted values with units; four factors; two sentences.

## 2. (Derivation: the interpolation error bound)

Let $y(t)$ be twice differentiable with $|y''(t)| \le M$ on $[t_i, t_{i+1}]$, $h = t_{i+1} - t_i$, and let $\hat y$ be the straight line through the end points.

(a) Prove that $|y(t) - \hat y(t)| \le Mh^2/8$ for all $t$ in the interval.
(b) Apply it: a barometric altitude channel has a 2.0 s gap during a pull-up with vertical acceleration up to $4\ \mathrm{m/s^2}$. Bound the error of filling the gap by interpolation. Compare with the sensor noise σ = 1.5 m and say whether the filled values can be trusted.
(c) Now the gap is 20 s during a steady climb at constant rate. What is $M$, what is the bound, and why does the long gap not matter here?

Answer form: proof with each step named; two numbers with units; one paragraph.

## 3. (Unfamiliar scenario: the clock wraps)

A low-cost logger stamps each row with a 16-bit millisecond counter that overflows from 65535 back to 0. A 300 s flight therefore shows time stamps that climb to 65.535 s, drop to zero, and climb again, four times.

(a) What does `parse_flight_log` do with this file, at which row, and why is that the right behaviour for a parser that knows nothing about the counter?
(b) Design an *unwrapping* step to run before the monotonic check: state the rule in one sentence and write it as 5 to 10 lines of NumPy (hint: `np.diff`, a cumulative count of wraps, `np.cumsum`).
(c) Name one failure your rule still cannot detect, and one consistency check on the data that would catch it.

Answer form: short argument, code, two sentences.

## 4. (Numerical, uses your code: Study A)

Generate a log with `n=600, seed=s, missing_fraction=0.10, altitude_unit="ft", airspeed_unit="kt"`. Parse it and resample at $\Delta t \in \{0.25, 0.5, 1, 2, 5\}$ s. For each $\Delta t$ report: number of samples, maximum altitude in m, and climb rate in m/s from a straight-line fit (`np.polyfit(t, alt, 1)`) over $5\ \mathrm{s} < t < T/3 - 5\ \mathrm{s}$ where $T$ is the record length.

Plot both quantities against $\Delta t$. Then compare with the truth built into the generator (cruise altitude 300 m; climb rate $300/(T/3)$ m/s) and with the interpolation bound from problem 2. Answer form: a 5-row table, two plots, one paragraph stating which quantity is robust to $\Delta t$ and why.

## 5. (Numerical, uses your code: how to fill a gap)

Build a known signal $y(t) = 100 \sin(2\pi t/40)$ m on $t = 0, 0.5, \ldots, 59.5$ s (a gentle altitude oscillation). Remove a run of $g \in \{1, 5, 20\}$ consecutive samples starting at $t = 25$ s. Fill the gap three ways: your `fill_missing` (linear), "hold last valid value", and "drop the rows" (compute statistics on the remaining samples only).

For each $g$ and each method report the RMS error over the gap (for "drop", report instead the error in the maximum and the mean of the record). Compare the linear RMS error with the bound $Mh^2/8$ using $h = (g+1) \times 0.5$ s and $M = 100 (2\pi/40)^2$ m/s². Answer form: a 3×3 table, the bound column, and a two-sentence rule for when each method is acceptable.

## 6. (Judgment: choosing a format)

An IMU logs 10 channels at 1 kHz for a 2 h flight.

(a) Estimate the file size as CSV (assume about 9 characters per value plus a time stamp of 10 characters) and as HDF5 with `float64` and with `float32` datasets. Show the arithmetic.
(b) `float32` carries about 7 significant figures. For an accelerometer channel with a range of ±160 m/s² and noise of 0.02 m/s², is `float32` adequate? For a GPS time stamp in seconds since 1980 (about $1.4 \times 10^9$ s) is it?
(c) Give one situation in which you would keep the CSV even though it is 2 to 3 times larger and 100 times slower to read.

Answer form: three sizes with the arithmetic; two yes/no with a number each; one paragraph.

## 7. (Hard, uses your code: aliasing)

Sample a pitch oscillation $\theta(t) = 5\sin(2\pi f t)$ degrees with $f = 1$ Hz at $\Delta t = 0.1$ s for 20 s. Resample it with `resample_uniform` onto grids of $\Delta t \in \{0.4, 0.6, 0.9\}$ s.

(a) For each grid, estimate the apparent frequency from the number of zero crossings and report the apparent amplitude.
(b) Explain the 0.6 s and 0.9 s results with the Nyquist condition and derive the apparent frequency formula $f_a = |f - k/\Delta t|$ for the nearest integer $k$.
(c) Your flight log is resampled to 0.5 s. What is the fastest pitch oscillation that can be represented, and what should the parser do if the raw record contains faster content? Propose a test (in prose) that would catch an aliased channel.

Answer form: a 3-row table, a derivation, a short argument.
