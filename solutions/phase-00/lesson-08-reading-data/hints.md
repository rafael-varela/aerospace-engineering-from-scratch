# P00.08  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Every factor is in `UNIT_FACTORS`. Only one of the four is a definition rather than a derived or irrational number. For (b), think about which wrong value a pilot would see as absurd (a 1066 m aircraft at 61 m/s pitched 0.05° is not absurd; what is?).
2. Define $e(t) = y(t) - \hat y(t)$; it is zero at both ends, so Rolle gives a point $\tau$ where $e'(\tau) = 0$. Taylor-expand $e$ about $\tau$ and remember $\hat y'' = 0$. The widest the interval can be on either side of $\tau$ is $h$, but you can do better with the two-sided argument: the maximum of $|e|$ is at most $h/2$ from one of the ends.
3. `np.diff(t) < 0` marks the wraps. The number of wraps *so far* at each sample is a cumulative sum of that boolean. Add $65.536 \times$ (that count) to the time.
4. Build a mask for the fit window from the resampled time, not the raw one. The truth for the climb rate is `300 / (T/3)` with `T = log.time[-1]`.
5. The bound uses $h$ as the full width between the two valid neighbours, which is one sample more than the gap length. "Hold" is a constant extrapolation, so its error grows with the *slope*, not the curvature.
6. $2\,\mathrm{h} \times 3600\,\mathrm{s/h} \times 1000\,\mathrm{Hz}$ rows. For `float32` compare the needed dynamic range (full scale divided by the smallest meaningful change) with $2^{24} \approx 1.7 \times 10^7$.
7. Count sign changes of the resampled signal and divide by two to get cycles. The sampled values of $\sin(2\pi f k\Delta t)$ are identical to those of $\sin(2\pi (f - n/\Delta t) k \Delta t)$ for any integer $n$.
