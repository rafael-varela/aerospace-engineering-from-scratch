# P00.08  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. Look up each unit's factor and target: ms→s, ft→m, kt→m/s, rad→deg. 2. Multiply. 3. Mark each factor exact (a definition) or rounded (an irrational or a repeating decimal). 4. For (b) divide the wrong reading by the right one to get the error factor, then ask which wrong value would look impossible in the cockpit and which would look plausible.

## 2.
(a) 1. Set $e = y - \hat y$; note $e(t_i) = e(t_{i+1}) = 0$ and $e'' = y''$ because the line has no curvature. 2. Fix an interior $t$ and let $w(x) = (x - t_i)(x - t_{i+1})$. Build $\varphi(x) = e(x) - e(t)\,w(x)/w(t)$, which vanishes at $t_i$, $t$ and $t_{i+1}$. 3. Rolle twice: $\varphi''$ vanishes at some $\xi$, giving $e(t) = \tfrac12 y''(\xi)\, w(t)$. 4. Bound $|w(t)| \le h^2/4$ (its maximum is at the midpoint). (b) Plug $M = 4$, $h = 2$. (c) Constant climb rate means $y'' = 0$.

## 3.
(a) Trace the parser: the monotonic check finds the first negative `np.diff`; it raises with that row. (b) Compute `wraps = np.cumsum(np.concatenate([[0], np.diff(t_raw) < 0]))`, then `t = t_raw + 65.536 * wraps`. (c) Consider a dropout longer than one full counter period, and a check based on the known nominal sample interval.

## 4.
1. `text = generate_flight_log(n=600, seed=s, missing_fraction=0.10, altitude_unit="ft", airspeed_unit="kt")`. 2. `log = parse_flight_log(text)`; `T = log.time[-1]`. 3. Loop over `dt`: `u = resample_log(log, dt)`; `sel = (u.time > 5) & (u.time < T/3 - 5)`; `np.polyfit(u.time[sel], u.altitude[sel], 1)[0]`; `np.max(u.altitude)`. 4. Plot both against `dt`. 5. Compare with 300 m and `300/(T/3)`; for the peak also compare with the noise σ and with $Mh^2/8$.

## 5.
1. Build `t`, `y`; copy `y` and set the gap to NaN. 2. Linear: `fill_missing(t, y_gap)`. Hold: forward-fill in a loop or with `np.maximum.accumulate` on valid indices. Drop: `y_gap[~np.isnan(y_gap)]`. 3. RMS over the gap indices for linear and hold; for drop, compare `max` and `mean` with those of the full record. 4. Bound: $M h^2 / 8$ with $M = 100(2\pi/40)^2$ and $h = (g+1) \cdot 0.5$.

## 6.
(a) rows $= 2 \times 3600 \times 1000$; CSV bytes $\approx$ rows $\times (10 + 10 \times 9)$; HDF5 $\approx$ rows $\times 11 \times 8$ or $\times 4$. (b) Dynamic range needed = full scale / resolution; compare with $2^{24}$. For the time stamp, compute the `float32` spacing near $1.4 \times 10^9$ (`np.spacing(np.float32(1.4e9))`). (c) Think of readability, diffing in version control, and tools that cannot open HDF5.

## 7.
1. `t = np.arange(0, 20, 0.1)`; `theta = 5*np.sin(2*np.pi*t)`. 2. `resample_uniform(t, theta, dt)` for each `dt`. 3. Zero crossings: `np.sum(np.diff(np.sign(theta_u)) != 0) / 2 / (t_u[-1] - t_u[0])`. 4. Amplitude: `np.max(np.abs(theta_u))`. 5. Derive: $\sin(2\pi f k\Delta t) = \sin(2\pi (f - n/\Delta t) k\Delta t)$; choose $n$ to bring the frequency into $[-1/(2\Delta t), 1/(2\Delta t)]$. 6. For (c) compute $1/(2 \times 0.5)$ and think about comparing the spectrum or the variance of the raw and resampled channel.
