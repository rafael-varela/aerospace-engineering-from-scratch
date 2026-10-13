# P00.08  Worked solutions

Each solution states: assumptions, units, the derivation, the result with uncertainty where applicable, at least one independent check, and the common mistakes. Seeded problems are worked for $s = 15$.

## 1. The header

**Derivation.** ms → s: $1500 \times 10^{-3} = 1.5$ s (exact, a decimal prefix). ft → m: $3500 \times 0.3048 = 1066.8$ m (exact: the international foot is *defined* as 0.3048 m). kt → m/s: $120 \times 1852/3600 = 61.733\ldots$ m/s (the nautical mile is defined as 1852 m, so the factor is an exact rational, but its decimal expansion repeats; keep it as a fraction in code). rad → deg: $0.05 \times 180/\pi = 2.8648$° (irrational factor; rounded).

(b) Treating the row as `s, m, m/s, deg`: time wrong by $1000\times$ (1500 s instead of 1.5 s), altitude by $1/0.3048 = 3.28\times$ (3500 m instead of 1067 m), airspeed by $1/0.5144 = 1.94\times$ (120 m/s instead of 61.7 m/s), pitch by $1/57.3 = 0.0175\times$ (0.05° instead of 2.86°). A pilot notices the time first: a 1500 s first sample on a 300 s flight is impossible. The pitch error survives: 0.05° is a perfectly plausible level-flight pitch, and nothing in the number itself says it is wrong. Small-looking factors applied to plausible ranges are the dangerous ones.

**Common mistakes.** Using 0.5144 to four figures and losing the exactness; converting rad to deg by multiplying by π/180 (the inverse).

## 2. Interpolation error bound

**Assumptions.** $y \in C^2[t_i, t_{i+1}]$, $|y''| \le M$, $h = t_{i+1} - t_i$.

**Derivation.** (a) Let $e = y - \hat y$. The line has $\hat y'' = 0$, so $e'' = y''$, and $e(t_i) = e(t_{i+1}) = 0$. Fix an interior point $t$ and define $w(x) = (x - t_i)(x - t_{i+1})$ and $\varphi(x) = e(x) - e(t)\,w(x)/w(t)$. Then $\varphi$ vanishes at the three points $t_i, t, t_{i+1}$. By Rolle's theorem $\varphi'$ vanishes at two points between them, and $\varphi''$ at one point $\xi$ between those. But $\varphi'' = e'' - 2e(t)/w(t)$, so $e(t) = \tfrac12 y''(\xi)\, w(t)$. Finally $|w(t)| \le h^2/4$, attained at the midpoint, giving $|e(t)| \le M h^2/8$. ∎

(b) $M = 4\ \mathrm{m/s^2}$, $h = 2$ s: bound $= 4 \times 4/8 = 2.0$ m. That is 1.3σ of the sensor noise, so the filled values are about as uncertain as a real sample, and acceptable, provided the pull-up really is bounded by 4 m/s².

(c) Steady climb means constant $\dot h$ and $M = 0$: the bound is zero and the straight line is exact, whatever the gap length. The long gap does not matter because the signal *is* a line there. The danger of long gaps is not their length but the curvature they may hide.

**Check.** Dimensions: $\mathrm{m/s^2 \times s^2 = m}$. Problem 5 measures RMS errors of 0.22, 1.8 and 23.8 m against bounds of 0.31, 2.8 and 34 m: always below, by a factor 1.4 to 1.6, as expected for a bound on the *maximum* compared with an RMS.

**Common mistakes.** Using the gap length instead of the distance between the two valid neighbours for $h$; forgetting that the bound needs a bound on $y''$ *during the gap*, which you cannot measure from the data you do not have.

## 3. The clock wraps

(a) The counter wraps at $65.536$ s, so around data row 132 (at 0.5 s nominal spacing) the stamp drops from about 65.5 to about 0.0. `parse_flight_log` converts the time column, finds `np.diff(t) <= 0` there, and raises `ValueError("time must be strictly increasing; first violation near data row 133")`. That is right: a parser that silently "fixed" the time would have to guess the counter width, and a wrong guess shifts every later sample by tens of seconds while looking perfectly healthy.

(b) Rule: every time the raw stamp decreases, add one more full counter period to all later stamps.

```python
PERIOD_S = 65.536
wraps = np.cumsum(np.concatenate([[0], (np.diff(t_raw) < 0).astype(int)]))
t = t_raw + PERIOD_S * wraps
assert np.all(np.diff(t) > 0.0)
```

Run this on the raw time column before the monotonic check (for example as an optional `unwrap_period` argument to `parse_flight_log`).

(c) It cannot detect a dropout longer than one counter period: two wraps happen with no sample in between, the rule adds one period instead of two, and the whole remainder of the flight is 65.5 s early. A consistency check: the nominal sample interval is known (0.5 s); any gap in the unwrapped time larger than, say, 50 nominal intervals is flagged, and a gap close to $65.536 - \epsilon$ s is reported as a probable missed wrap.

**Common mistakes.** Detecting wraps with `diff <= 0` after rounding to three decimals, which also catches repeated stamps; unwrapping *after* `convert_to_si` when the counter is in milliseconds and the period was written in seconds.

## 4. Study A (seed 15)

**Assumptions.** Generator truth: cruise 300 m, climb at $300/(T/3)$ with $T = 299.5$ s, so 3.005 m/s; noise σ = 1.5 m; 10 % of cells missing.

**Result.**

| $\Delta t$ [s] | samples | max altitude [m] | climb rate [m/s] |
|---|---|---|---|
| 0.25 | 1199 | 302.96 | 2.999 |
| 0.5 | 600 | 302.96 | 2.999 |
| 1.0 | 300 | 302.89 | 3.000 |
| 2.0 | 150 | 302.08 | 3.002 |
| 5.0 | 60 | 301.15 | 3.007 |

The raw record's own `nanmax` is 303.37 m.

**Interpretation.** The peak is $300 + (2.0$ to $3.4)$ m: the highest of about 200 noisy cruise samples sits about 2σ above the mean, and coarser grids skip more of those noise peaks, so the reported maximum falls monotonically with $\Delta t$ (302.96 → 301.15 m) while never going below 300 m. The climb rate is within 0.2 % of the truth at every $\Delta t$: a straight-line fit over 90 s of data averages noise and interpolation alike, and interpolation is exact on a linear segment ($M = 0$, problem 2c). **Peaks are fragile, slopes are robust.** Uncertainty to report: maximum altitude $303 \pm 2$ m (dominated by noise, not by the parser); climb rate $3.00 \pm 0.01$ m/s.

**Check.** The bound $Mh^2/8$ with $M \approx 0$ in cruise predicts no interpolation error there, so the $\Delta t$ dependence of the peak must come from *sampling* the noise, not from interpolation error; consistent with the peak never exceeding the raw maximum.

**Common mistakes.** Fitting the climb rate over the whole first third including the corner at $T/3$; comparing with 3.000 m/s instead of $300/(T/3)$; reading the peak off the raw record with `np.max` and getting NaN.

## 5. How to fill a gap

**Assumptions.** $y = 100\sin(2\pi t/40)$ m, $\Delta t = 0.5$ s, gap starting at $t = 25$ s; $M = 100(2\pi/40)^2 = 2.467\ \mathrm{m/s^2}$.

**Result** (RMS error over the gap, m):

| $g$ | $h = (g+1)\Delta t$ [s] | linear | hold | bound $Mh^2/8$ | drop: Δmax, Δmean |
|---|---|---|---|---|---|
| 1 | 1.0 | 0.22 | 5.8 | 0.31 | 0, +0.8 |
| 5 | 3.0 | 1.8 | 16.8 | 2.8 | 0, +4.4 |
| 20 | 10.5 | 23.8 | 26.6 | 34 | 0, +22 |

Δmax is zero only because the 60 s record holds two equal peaks (at 10 s and 50 s) and the gap touches neither. Δmean is large because the gap sits in the trough ($y < 0$ for $20 < t < 40$ s): removing low samples raises the mean of what remains. "Drop" is unbiased only when gaps are unrelated to the signal's value, which a radio dropout during a dive is not.

**Rule.** Linear filling is acceptable while $Mh^2/8$ is below the noise or the tolerance you need (here up to about $g = 3$ for a 1 m tolerance). Holding the last value is never better than linear for a smooth signal: its error grows with the slope times the gap. Dropping rows is the honest choice for statistics (mean, slope) when gaps are random, but it biases any extremum that falls in a gap and it breaks uniform-grid methods.

**Common mistakes.** Reporting the error of "drop" as zero because no filled values exist; using the gap length $g\Delta t$ instead of $(g+1)\Delta t$ for $h$.

## 6. Choosing a format

(a) Rows: $2 \times 3600 \times 1000 = 7.2 \times 10^6$. CSV: $7.2 \times 10^6 \times (10 + 10 \times 9) = 7.2 \times 10^8$ bytes ≈ 0.72 GB. HDF5 `float64`: $7.2 \times 10^6 \times 11 \times 8 = 0.63$ GB. `float32`: 0.32 GB. (The CSV is only modestly larger here because 9 characters per value is already tight; with typical 12-character values it is 1.0 GB. The bigger difference is parse time: text must be converted digit by digit.)

(b) Accelerometer: needed dynamic range $320/0.02 = 1.6 \times 10^4$ levels, far below `float32`'s $2^{24} = 1.7 \times 10^7$. Yes, adequate. GPS time in seconds since 1980: near $1.4 \times 10^9$ the `float32` spacing is $2^{30-23} = 128$ s. No: it cannot even resolve one sample. Store time as `float64`, or as an integer count of milliseconds, or as seconds since the start of the file.

(c) Keep the CSV when humans must read or diff it: a configuration or calibration log that goes into version control, a file that must be opened by a spreadsheet, or a legal record where the raw bytes are the evidence. Readability is a requirement, not a convenience, in those cases.

**Common mistakes.** Forgetting the time column in the HDF5 estimate; believing `float32` is "7 digits" everywhere without computing the spacing at the actual magnitude.

## 7. Aliasing

**Result.**

| grid $\Delta t$ [s] | $1/\Delta t$ [Hz] | Nyquist [Hz] | apparent $f$ [Hz] | apparent amplitude [°] |
|---|---|---|---|---|
| 0.4 | 2.5 | 1.25 | 1.02 | 4.76 |
| 0.6 | 1.667 | 0.833 | 0.68 | 4.76 |
| 0.9 | 1.111 | 0.556 | 0.13 | 4.76 |

(The amplitude is 4.76 rather than 5 for every grid because the resampled values are linear interpolations of the 0.1 s source samples, which already clip the peaks slightly.)

**Derivation.** On the grid $t_k = k\Delta t$, $\sin(2\pi f t_k) = \sin(2\pi f k\Delta t - 2\pi n k) = \sin(2\pi (f - n/\Delta t) t_k)$ for any integer $n$, since $2\pi nk$ is a multiple of $2\pi$. The grid cannot distinguish $f$ from $f - n/\Delta t$; the signal that *looks* right is the one with the smallest $|f - n/\Delta t|$, so $f_a = \min_n |f - n/\Delta t|$, which is below the Nyquist frequency $1/(2\Delta t)$. For 0.6 s: $|1 - 1.667| = 0.667$ Hz (measured 0.68). For 0.9 s: $|1 - 1.111| = 0.111$ Hz (measured 0.13; the zero-crossing count is coarse with only 2.5 cycles in 20 s). For 0.4 s the Nyquist frequency is 1.25 Hz > 1 Hz, so no aliasing and $f_a = f$.

(c) At $\Delta t = 0.5$ s the fastest representable pitch oscillation is 1 Hz, and anything faster folds back into the band. The parser cannot fix this (the information is gone once resampled), so it should warn when the source spacing is much finer than $\Delta t$ and the channel has content above the new Nyquist frequency. A test: compare the variance of the channel before and after resampling, or compute the fraction of the source signal's power above $1/(2\Delta t)$ with `np.fft.rfft`; flag the channel if it exceeds a few percent. The lesson's own flight log passes trivially: its pitch changes on a 100 s time scale.

**Common mistakes.** Treating $1/\Delta t$ as the highest representable frequency (it is twice that); concluding that resampling to a finer grid "adds resolution".

## Check-yourself answers

1. 3500 ft = 1066.8 m; 120 kt = 61.73 m/s; 0.05 rad = 2.865°. The foot (0.3048 m) is exact by definition; so is the nautical mile (1852 m), but the knot factor 1852/3600 has a repeating decimal.
2. Every arithmetic operation with NaN yields NaN, so the sum, and hence the mean, is NaN. Use `np.nanmean` (and `np.nanmax`, `np.nanmin`), or mask with `~np.isnan`.
3. Two samples logged within one clock tick of a coarse clock; or a buffered retransmission after a link dropout wrote the same row twice. The parser raises `ValueError` naming the row, because interpolation needs strictly increasing time.
4. $|y - \hat y| \le Mh^2/8$. With $M = 1\ \mathrm{m/s^2}$, $h = 2$ s: 0.5 m.
5. Raise: no, a linear interpolant never leaves the range of its two neighbours. Lower: yes, if the true peak fell between source samples or on a sample the coarser grid skips.
6. Nyquist frequency $1/(2 \times 0.6) = 0.833$ Hz; the 1 Hz signal aliases to $|1 - 1/0.6| = 0.67$ Hz and you see a slower oscillation of about the same amplitude.
7. Advantages: units and metadata stored inside the file as attributes; binary arrays load in seconds and can be read partially (one channel, one time slice) without parsing the whole file. Keep the CSV when the record must be human-readable or diffable, or when a consumer cannot open HDF5.
8. A warning scrolls past and the wrong number flows on into every downstream calculation (the Mars Climate Orbiter pattern). Raising stops the pipeline at the one place where the person who knows the unit can fix it.
