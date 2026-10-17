# P09.07  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
N = 20 + s; zeros at k·1000/N, k = 1..7; test 180N/1000 for an integer. `b = moving_average(N)`; `abs(digital_freq_response(b, [1.0], f, 1000))` for f = 2, 10, 180; `group_delay(b, [1.0], 2.0, 1000)` = (N − 1)/2 samples.

## 2.
(a) Σ_{k=0}^{N−1} e^{−jωk} = (1 − e^{−jωN})/(1 − e^{−jω}); factor e^{−jωN/2}(e^{jωN/2} − e^{−jωN/2}) and e^{−jω/2}(e^{jω/2} − e^{−jω/2}); the ratio of sines; zeros where Nω/2 = kπ with ω ≠ 0 mod 2π. (b) h[k]e^{−jωk} + h[k]e^{−jω(N−1−k)} = h[k]e^{−jω(N−1)/2}(e^{jω((N−1)/2−k)} + e^{−jω(…)}) = 2h[k]cos(…)e^{−jω(N−1)/2}; sum over pairs (and the middle tap for odd N); where the real factor is negative the phase jumps by π. (c) y[n] − y[n−1] = (T/2)(x[n] + x[n−1]); transfer (T/2)(1 + z⁻¹)/(1 − z⁻¹) ↔ 1/s; invert to z = (K + s)/(K − s); \|K + s\|² − \|K − s\|² = 4K Re s; on s = jΩ, \|z\| = 1; z = e^{jω} ⇒ s = K(1 − e^{−jω})/(1 + e^{−jω}) = jK tan(ω/2); set Ω = 2πf_a, ω = 2πf_d/f_s; for prewarping require K tan(πf_p/f_s) = 2πf_p.

## 3.
ω₀ = 2π·180; `bilinear(np.array([1, 0, w0**2]), np.array([1, w0/Q, w0**2]), 1000, f_prewarp=180)`; evaluate `digital_freq_response` on a fine grid for the −3 dB band; `group_delay` at 2 Hz; `np.abs(np.roots(az))`. For (b): `t, true, vib, noise = synthetic_gyro_record(s)`, run `difference_equation` on vib alone (residual via `tone_amplitude` at 180 Hz over t ≥ 0.5), on noise alone (std), and on the sum (RMS against true).

## 4.
`gyro_filter_comparison(seed=s)`; `windowed_sinc_lowpass(N, 40, 1000, 'hamming')` for N = 51, 101, 201 and its gain at 180 Hz; delay (N − 1)/2 ms. Rank by phase at 5 Hz (or delay) first, residual line second.

## 5.
For each f_c: `butterworth_analog(2, 2*np.pi*fc)`; `bilinear` with and without `f_prewarp=fc`; \|H(f_c)\|; the −3 dB frequency by searching a fine grid for \|H\| = 1/√2; compare with (f_s/π) arctan(πf_c/f_s).

## 6.
(a) `butterworth_digital(4, fc, 1000)`; `np.abs(np.roots(a)).max()`; `ar = np.round(a, 4)`; same for `ar`; DC gain `b.sum()/ar.sum()`. (b) `np.angle(digital_freq_response(b, a, 5, 1000))` for each of `lesson_filters()`; compare with 45°. (c) Forward-backward (zero-phase) filtering of a stored record.
