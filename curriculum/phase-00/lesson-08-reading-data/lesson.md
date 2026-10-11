# P00.08  Reading Data: CSV, HDF5, Telemetry Logs

> The file is not the data. The data is what the file meant, and the file rarely says.

| | |
|---|---|
| **Learning outcomes** | 1. Parse a text telemetry log with comment lines, units in the header, irregular time stamps and missing values into SI-unit NumPy arrays, refusing malformed input with a clear error. 2. Fill gaps by interpolation and resample a channel onto a uniform time base, and bound the error of doing so. 3. Decide between CSV and HDF5 for a given dataset, and round-trip a log through HDF5 when the library is available. 4. Recognise the three classic telemetry failures (wrong unit, time going backwards, silent NaN propagation) in a file you have never seen. |
| **Prerequisites** | `P00.03`, `P00.05`, `P00.07`; skill check: can index a NumPy array with a boolean mask; can convert feet to metres and knots to m/s by hand. |
| **Engineering question** | From a flight log whose rows arrive late, early or not at all, what was the aircraft's maximum altitude and climb rate, and how far can we trust those numbers? |
| **Assumptions and limits** | Channels vary smoothly between samples (linear interpolation is valid when the sample spacing is short compared with the time scale of the motion); time stamps are correct even when irregular; missing values are missing at random, not during a fault; units in the header are truthful. The method fails for signals that oscillate faster than half the sample rate (aliasing), for gaps longer than the motion's time scale, and for logs whose clock wraps or resets. |
| **Analytical work** | Unit conversions by hand; the error bound of linear interpolation, $M h^2/8$; a size estimate for CSV versus binary storage. |
| **Computational deliverable** | `starter/main.py`: `parse_header`, `convert_to_si`, `parse_flight_log`, `read_flight_log`, `fill_missing`, `resample_uniform`, `resample_log`, `summarize`, `write_hdf5`, `read_hdf5`. Tests prove: units round-trip (the same flight in ft/kt and in SI parses to the same numbers), malformed logs raise, linear signals are filled and resampled exactly, interpolation error decays as the square of the source spacing, resampling never invents a higher peak, HDF5 round-trips bit-for-bit (or is skipped with a message). |
| **Practical activity** | **data-based investigation** (synthetic data generated in code; no network): study how the recovered maximum altitude and climb rate depend on the resampling step and on the fraction of missing samples. |
| **Equipment** | None. `h5py` is optional; the lesson says so where it matters. |
| **Safety** | None beyond normal computer use |
| **Evidence** | Completed `main.py`, passing tests, `flight_log.png`, the two sensitivity tables from part 4 with plots, and a short note stating the uncertainty you attach to the maximum altitude and the climb rate. |
| **Assessment** | `assessments/phase-00/grading-rubric.md`; competencies 3 (computational) and 4 (interpretation); pass when all public tests pass *and* the note gives a numerical uncertainty with a reason. |
| **References** | `@numpy-userguide` "NumPy fundamentals" > "Data types" (NaN, float64) and "NumPy: the absolute basics for beginners" (required-free); `@nasa-mco-mib` Mars Climate Orbiter Mishap Investigation Board Phase I Report, root cause section (required-free, case study for the unit rule). |
| **Extensions** | ★ Read a real log: export a track from a flight-tracking site or a phone GPS app (GPX or CSV), identify the units and gaps, and run it through your parser. ★ Add a `time [hh:mm:ss]` format to `convert_to_si`. ★ Compare `np.interp` with `scipy.interpolate.CubicSpline` on the climb segment. |
| **Time** | unmeasured |

## 1. The engineering problem

Every flight test, wind-tunnel run and rocket launch ends with a file. Most of the engineering happens after that file is opened, and most of the mistakes happen while opening it. A telemetry log from a small aircraft or a model rocket typically arrives as **CSV** (comma-separated values: plain text, one sample per line, one value per column). The header may say `altitude [ft]` or `alt_m` or just `alt`. Some rows are missing because the radio link dropped. The time stamps are not evenly spaced because the logger wrote a line "when it could". Somewhere there is a blank cell, a `nan`, or a `---`.

The question we want answered is simple: what was the maximum altitude, what was the climb rate, and how sure are we? The answers are only as good as the parsing. If `altitude [ft]` is read as metres, the peak is wrong by a factor of 3.28. If a `nan` becomes a zero, the mean altitude drops and the "climb rate" acquires a spike. If the time stamps are used as if they were uniform, every derivative is wrong by the jitter. The success criterion for this lesson is a parser that gets the physical quantities right in SI units, says so, refuses to guess, and attaches a number to "how sure".

We use a synthetic log, generated in code by `generate_flight_log`, so that the truth is known and the tests can check the parser against it. This is a model of a flight log, not a flight. Part 4's extension points you at a real one.

## 2. Understand and derive

### 2.1 Units live in the header, and the header lies by omission

A column header such as `airspeed [kt]` carries two facts: the name of the channel and its unit. Our convention is `name [unit]`, with a known set of units. Conversion is multiplication by a fixed factor:

| from | to | factor | note |
|---|---|---|---|
| ft | m | 0.3048 | exact by definition (international foot) |
| kt | m/s | 1852/3600 = 0.51444… | one nautical mile per hour |
| km/h | m/s | 1/3.6 | |
| ms | s | 0.001 | |
| rad | deg | 180/π | we keep angles in degrees, the flight-data convention |

A header with no unit tells you nothing. Our rule is to accept the number *unchanged and unlabelled* (unit `""`) so the gap is visible downstream, and to refuse any unit we do not recognise with `ValueError`. Refusing feels unfriendly; it is what the Mars Climate Orbiter team wished their software had done when a file of thruster impulses in pound-force seconds was read as newton seconds (`@nasa-mco-mib`).

### 2.2 Missing values and NaN

A missing cell becomes `NaN` ("not a number"), a special `float64` value. NaN has two properties you must know: every arithmetic operation involving NaN returns NaN (so one gap poisons a whole `np.mean`), and every comparison with NaN is `False` (so `x == np.nan` never finds it; use `np.isnan`). NumPy provides `np.nanmax`, `np.nanmean` and friends that ignore NaN. The time column may *not* contain NaN: a sample without a time stamp cannot be placed, so the parser rejects it.

### 2.3 Time must go forward

Telemetry time stamps are irregular but must be strictly increasing. A repeated stamp means two samples claim the same instant; a decreasing one means the clock was reset, wrapped around (a 16-bit millisecond counter overflows every 65.536 s), or the rows were reordered. Interpolation with non-monotonic time silently produces garbage, so `parse_flight_log` checks `np.all(np.diff(t) > 0)` and raises with the row number.

### 2.4 Filling gaps and resampling: linear interpolation and its error

Between two valid samples $(t_i, y_i)$ and $(t_{i+1}, y_{i+1})$ the straight line through them is

$$\hat y(t) = y_i + \frac{y_{i+1} - y_i}{t_{i+1} - t_i}\,(t - t_i).$$

`np.interp(t_new, t_valid, y_valid)` evaluates this for every point in `t_new`, holding the end values constant outside the range. We use it twice: once to replace NaNs (evaluate at the gap's own time stamps using only the valid neighbours) and once to move all channels onto a uniform grid $t_k = t_0 + k\,\Delta t$.

**How wrong is the straight line?** If the true signal $y(t)$ has a second derivative bounded by $|y''| \le M$ on an interval of length $h$, then the interpolation error anywhere in that interval satisfies

$$|y(t) - \hat y(t)| \le \frac{M h^2}{8}.$$

(Proof sketch: the error $e(t) = y - \hat y$ vanishes at both ends; by Rolle's theorem $e'$ vanishes somewhere inside; Taylor expansion of $e$ about that point with $e'' = y''$ gives $|e| \le \tfrac{M}{2}(h/2)^2$.) Two consequences for flight data. First, the error is quadratic in the spacing: halving the gap quarters the error, which the test `test_resample_uniform_is_exact_for_linear_and_second_order_otherwise` checks. Second, $M$ is a physical quantity. For altitude, $y''$ is vertical acceleration; a light aircraft rarely exceeds $2\ \mathrm{m/s^2}$ in smooth flight, so a 0.5 s gap costs at most $2 \times 0.25/8 = 0.06$ m, while a 20 s gap (the radio dropped out) costs up to 100 m, which is the whole climb. The bound tells you when interpolation is honest and when it is fiction.

Linear interpolation also never produces a value outside the range of its two neighbours, so a resampled peak is never *higher* than the raw peak; it can be lower if the true peak fell between samples. That is a bias, not noise, and part 4 measures it.

### 2.5 Aliasing: the limit of resampling

A sinusoid of frequency $f$ sampled at interval $\Delta t$ is recoverable only if $\Delta t < 1/(2f)$ (the Nyquist condition). Resample a 1 Hz pitch oscillation onto a 0.6 s grid and you get a convincing 0.68 Hz oscillation that does not exist. Resampling to a *coarser* grid than the source can therefore invent frequencies; resampling to a finer one cannot add information but is safe. Problem 7 makes you see it.

### 2.6 CSV versus HDF5

CSV is human-readable, universal, and inefficient: each `float64` becomes about 9 characters, units live in a header that may be missing, and a 2-hour, 10-channel, 1 kHz IMU log is about 0.7 GB of text that takes minutes to parse. **HDF5** (Hierarchical Data Format 5) is a binary container: named datasets (arrays) with attributes (metadata such as units) in a tree, readable by every scientific language, and typically read in seconds. The lesson's rule: CSV for small logs that people must read and diff; HDF5 (or similar) for anything over a few megabytes or anything whose units you want carried *inside* the file. The Python library is `h5py`; it is not in this course's base requirements, so `write_hdf5` must detect its absence and say so clearly instead of crashing.

## 3. Build the computational model

Open `starter/main.py`. `generate_flight_log` is given, complete; everything else is yours. Run the tests after each function.

### Step 1: `parse_header`

Split on commas; for each field, match `name [unit]` with an optional bracket group. A regular expression does it in one line:

```python
import re
m = re.fullmatch(r"\s*([^\[\]]+?)\s*(?:\[\s*([^\[\]]*?)\s*\])?\s*", raw)
name, unit = m.group(1).strip().lower(), (m.group(2) or "").strip()
```

Lower-casing names means `Altitude` and `altitude` are the same channel. Raise on an empty name.

### Step 2: `convert_to_si`

Look the unit up in `UNIT_FACTORS`; multiply; return the converted array and the target unit. Empty unit: pass through. Unknown unit: `ValueError`. Resist the temptation to "assume metres".

### Step 3: `parse_flight_log`

```python
lines = [ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
header = parse_header(lines[0])
```

Then, for each data line: split, check the field count equals the header length, convert each token with a helper that returns `np.nan` for anything in `MISSING_TOKENS` and `float(token)` otherwise (a token that is neither raises `ValueError`; wrap it so the message carries the line number). Build a 2-D array, pick out the four channels by name, convert each with `convert_to_si`, then check the time column: no NaN, strictly increasing. Store the header units in `source_units` and the converted units in `units`.

### Step 4: `fill_missing` and `resample_uniform`

```python
ok = ~np.isnan(y)
filled = np.interp(t, t[ok], y[ok])          # fill gaps at the original stamps
n = int(np.floor((t[-1] - t[0]) / dt + 1e-9)) + 1
t_u = t[0] + dt * np.arange(n)               # uniform grid, never past t[-1]
y_u = np.interp(t_u, t, filled)
```

The `1e-9` guards against `(t[-1]-t[0])/dt` evaluating to `2.9999999` when it should be 3 (P00.07 again). Validate: `dt > 0`, same shapes, strictly increasing `t`, no NaN left in `y`.

### Step 5: `summarize`, HDF5, and `main`

`summarize` uses `np.nanmax` and counts NaNs across the three data channels. `write_hdf5` checks `HAVE_H5PY` first; if present, one dataset per channel with a `units` attribute. `read_hdf5` reverses it. Run `python main.py`: it writes `flight_log.csv` in ft and kt, reads it back in SI, resamples, round-trips HDF5 if possible, and saves `flight_log.png`.

### Verify

- **Analytic case:** a linear signal is filled and resampled exactly (`test_fill_missing_recovers_a_linear_signal_exactly`, `test_resample_uniform_is_exact_for_linear...`).
- **Limiting cases:** zero missing fraction parses with no NaN; a uniform source grid resampled at its own spacing reproduces itself; `dt` larger than the record gives a single sample.
- **Conservation-like property:** the resampled peak never exceeds the raw peak, the mean airspeed is preserved to 2 %, and the units dictionary survives (`test_resample_log_preserves_physics_of_the_profile`).
- **Convergence:** interpolation error on $t^2$ decays with rate $\ge 1.8$ in $\log_2$ as the source spacing halves.
- **Units round trip:** the same flight written in ft/kt and in SI parses to the same metres and m/s within the file's 3-decimal rounding.

### Compare

Compare `np.interp` with `scipy.interpolate.interp1d(kind="linear")` on the whole record: the two agree exactly (maximum difference 0.0), as they should, since both evaluate the same straight lines. Then try `scipy.interpolate.CubicSpline`: it follows the noise more closely and reports a different peak (303.29 m against 302.96 m for seed 15), and because a cubic is free to bend between samples it can overshoot where the signal has a corner (the climb-to-cruise transition), where the second derivative is not bounded. The error bound of 2.4 assumed smoothness; a corner is where that assumption fails, and the simpler method, which never leaves the range of its neighbours, is the safer one for a noisy record.

## 4. Experiment at home, or investigate real data

**Mode: data-based investigation** on synthetic data with known truth (climb to 300 m at 3.0 m/s, cruise, descend; noise σ = 1.5 m on altitude).

**Study A: resampling step.** Generate a log with your seed, `missing_fraction=0.10`, imperial units. Parse, then resample at $\Delta t \in \{0.25, 0.5, 1, 2, 5\}$ s. For each, record the maximum altitude and the climb rate (slope of a straight-line fit to altitude over the first third of the record, excluding 5 s at each end). Write your prediction first: which way does the peak move as $\Delta t$ grows, and does the climb rate move at all?

**Study B: missing fraction.** Fix $\Delta t = 0.5$ s and vary `missing_fraction` in $\{0, 0.05, 0.10, 0.20, 0.40\}$. Record the same two numbers and the maximum gap length (in samples) the generator produced.

Raw data format: two CSV tables, `dt_s, max_alt_m, climb_rate_m_s` and `missing_fraction, max_gap_samples, max_alt_m, climb_rate_m_s`. Save both plots.

## 5. Validate and improve

Lay the truth (300 m, 3.0 m/s), the parser's answer and the error bound side by side. You should find the maximum altitude a few metres *above* 300 m at fine $\Delta t$ (the noise, σ = 1.5 m, means the highest of ~200 cruise samples sits about 2.5σ above the mean), decreasing slowly as $\Delta t$ grows (coarse grids skip the noise peaks), and never below the true cruise altitude by more than the bound $Mh^2/8$ with $M$ the largest vertical acceleration in the profile. The climb rate is robust: a straight-line fit over 90 s averages out both noise and resampling, and it stays within 0.3 % of 3.0 m/s for every $\Delta t$ up to 5 s. Report that asymmetry: **peaks are fragile, slopes are robust**, and say why in terms of 2.4.

For the missing-value study, interpolation hides gaps of a few samples with error below the noise: with seed 15 the peak stays at 302.96 m from 0 % to 20 % missing, and the climb rate wanders by less than 0.2 %. At 40 % missing the longest run of consecutive gaps grows to about 7 samples (3.5 s). Inside the cruise that still costs nothing measurable, but if such a gap straddled the climb-to-cruise corner the straight line would cut the corner by up to $(\Delta\text{slope})\,h/4 = 3 \times 3.5/4 \approx 2.6$ m, above the noise. The honest uncertainty on the peak is therefore not σ alone but the larger of σ and the corner-cutting error of the longest gap that touches the peak region. If you discover a gap longer than the motion's time scale, the right improvement is not a better interpolant: it is to report the gap and give the peak as a range.

## 6. Defend your work

Submit `main.py`, the test log, `flight_log.png`, the two study tables with plots, and a half-page note using `projects/_template/report-template.md` that states: maximum altitude with an uncertainty and the reason for it; climb rate with an uncertainty; the longest gap and what you did about it.

### Adversarial self-checks

- What if the sample rate doubles? Which numbers in your tables change, which stay, and does the HDF5 file double in size while the CSV more than doubles?
- Which assumption fails first in the real world? Probably "units in the header are truthful": a logger firmware update changes `[kt]` to `[m/s]` without changing the data. What single test on the *values* (not the header) would catch a 1.94× airspeed?
- What measurement could falsify your interpolation? A second, independent altitude source (GPS versus barometric) sampled during a gap.
- Does the result survive a change in the resampling step? Study A says the climb rate does and the peak does not; state the $\Delta t$ beyond which the peak error exceeds the noise.
- Which decision changes if the noise doubles? With σ = 3 m, is a 2 m difference between two flights' peaks still worth reporting?

## Check yourself

1. Convert 3500 ft, 120 kt and 0.05 rad to SI (or degrees) by hand; which conversion factor is exact by definition?
2. Why does `np.mean` of an array with one NaN return NaN, and which function should you use instead when a gap is present?
3. A log has two rows with the same time stamp. Name two possible causes and say what your parser does.
4. State the linear interpolation error bound and compute it for a 2 s gap in altitude with vertical acceleration up to 1 m/s².
5. Can resampling onto a uniform grid raise the maximum of a channel? Lower it? Explain with one sentence each.
6. A 1 Hz oscillation is resampled onto a 0.6 s grid. What is the fastest signal that grid can represent, and what will you see?
7. You receive a 3 GB CSV of IMU data. Name two concrete advantages of converting it to HDF5 and one reason you might keep the CSV anyway.
8. Why does the parser raise on an unknown unit instead of warning and continuing?
