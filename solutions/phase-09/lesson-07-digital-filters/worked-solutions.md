# P09.07  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown. Every gyro record here is SYNTHETIC.

## 1. A moving average on a 1 kHz gyro (s = 15: N = 35)

**Assumptions.** f_s = 1000 Hz; taps 1/35.

(a) Zeros at k·1000/35 Hz: **28.57, 57.14, 85.71, 114.29, 142.86, 171.43, 200.00 Hz**. 180 × 35/1000 = 6.3 is not an integer, so **180 Hz is not nulled**.

(b) **\|H(2 Hz)\| = 0.9920, \|H(10 Hz)\| = 0.8105, \|H(180 Hz)\| = 0.0431** (−27.3 dB, a sidelobe between the zeros at 171 and 200 Hz). Delay **(35 − 1)/2 = 17 samples = 17.0 ms**.

(c) `digital_freq_response` and `group_delay` return these values.

**Check.** (3) at f = 10 Hz: sin(35π·0.01)/(35 sin(π·0.01)) = sin(1.0996)/(35 × 0.031411) = 0.8910/1.0994 = 0.8105.

**Common mistakes.** Putting the first zero at f_s/(2N); quoting the delay as N samples.

## 2. Three properties

(a) Σ_{k=0}^{N−1} e^{−jωk} = (1 − e^{−jωN})/(1 − e^{−jω}) = [e^{−jωN/2}·2j sin(Nω/2)]/[e^{−jω/2}·2j sin(ω/2)], so **H = (1/N) e^{−jω(N−1)/2} sin(Nω/2)/sin(ω/2)**. It vanishes where sin(Nω/2) = 0 but sin(ω/2) ≠ 0: ω = 2πk/N, k = 1, …, N − 1, that is **f = k f_s/N**. At ω = 0 the limit is 1.

(b) Pair k with N − 1 − k (equal taps): h[k](e^{−jωk} + e^{−jω(N−1−k)}) = 2h[k] cos(ω((N−1)/2 − k)) e^{−jω(N−1)/2}. For odd N the middle tap is h[M]e^{−jωM} with M = (N − 1)/2. Summing, **H = e^{−jω(N−1)/2} R(ω)** with R real. Where R > 0 the phase is exactly −ω(N − 1)/2; where R < 0 it is that minus π. The jumps of π happen only where R changes sign (in the stopband, at its zeros), so the group delay is (N − 1)/2 everywhere it is defined.

(c) The trapezoidal integrator y[n] = y[n−1] + (T/2)(x[n] + x[n−1]) has Y/X = (T/2)(1 + z⁻¹)/(1 − z⁻¹); replacing 1/s by it gives **s = (2/T)(1 − z⁻¹)/(1 + z⁻¹) = K(1 − z⁻¹)/(1 + z⁻¹)** with K = 2f_s. Solving for z: z = (K + s)/(K − s). With s = σ + jΩ, \|K + s\|² − \|K − s\|² = 4Kσ, so **σ < 0 ⇔ \|z\| < 1** and σ = 0 ⇔ \|z\| = 1. On the circle, s = K(1 − e^{−jω})/(1 + e^{−jω}) = K(e^{jω/2} − e^{−jω/2})/(e^{jω/2} + e^{−jω/2}) = **jK tan(ω/2)**. The digital filter at ω = 2πf_d/f_s therefore reproduces the analog one at Ω = K tan(πf_d/f_s). With K = 2f_s, the analog corner Ω = 2πf_a appears where 2πf_a = 2f_s tan(πf_d/f_s): **f_d = (f_s/π) arctan(πf_a/f_s)**. Requiring K tan(πf_p/f_s) = 2πf_p gives **K = 2πf_p/tan(πf_p/f_s)**.

**Check.** For small f_a/f_s, arctan x ≈ x and f_d ≈ f_a; the prewarped K tends to 2f_s as f_p → 0.

**Common mistakes.** Writing s = K(z − 1)/(z + 1) and then mixing powers of z and z⁻¹ in the code; forgetting that prewarping fixes only one frequency.

## 3. A notch for the vibration line (s = 15)

**Assumptions.** Prewarped bilinear transform at 180 Hz; f_s = 1000 Hz.

(a)

| | Q = 2 | Q = 10 |
|---|---|---|
| b | [0.815523, −0.694466, 0.815523] | [0.956717, −0.8147, 0.956717] |
| a | [1, −0.694466, 0.631046] | [1, −0.8147, 0.913434] |
| \|H(180 Hz)\| | 4 × 10⁻¹⁶ (a zero) | 2 × 10⁻¹⁵ |
| \|H(2 Hz)\| | 0.999988 | 1.000000 |
| −3 dB band | 146.4–217.2 Hz | 172.9–187.3 Hz |
| \|H(170)\|, \|H(190)\| | 0.272, 0.264 | 0.816, 0.808 |
| delay at 2 Hz | 0.394 ms | 0.079 ms |
| pole radius | 0.7944 | 0.9557 |

(b) SYNTHETIC record, seed 15 (not a measurement): both notches remove the 180 Hz line completely (residual 0.00000 deg/s). The output noise is **0.2675 deg/s (Q = 2)** and **0.2905 deg/s (Q = 10)** against 0.3 in, and the RMS error against the true rate is the same, 0.2676 and 0.2905 deg/s: the notch does nothing about broadband noise.

(c) A notch is right when the interference is a narrow line at a known frequency and the loop cannot afford the delay of a low-pass: it costs a fraction of a millisecond. It is wrong when the line moves (a Q = 10 notch passes 82 % of the line 10 Hz away) or when the problem is broadband noise. In practice flight controllers combine a gentle low-pass with a notch whose centre follows the measured rotor speed.

**Check.** The zeros of b(z) lie on the unit circle at angle ±2π·180/1000 = ±1.131 rad: the roots of 0.815523z² − 0.694466z + 0.815523 have product 1 and real part 0.694466/(2 × 0.815523) = 0.4258 = cos 1.131.

**Common mistakes.** Not prewarping (the notch lands below 180 Hz); expecting a notch to reduce noise.

## 4. Your own SYNTHETIC comparison (s = 15)

(a) SYNTHETIC record, seed 15 (not a measurement):

| filter | residual line | noise out (Parseval) | delay at 2 Hz | RMS error |
|---|---|---|---|---|
| MA 50 | 0.0000 deg/s | 0.0432 (0.0424) | 24.5 ms | 0.510 deg/s |
| FIR 101 Hamming | 0.0007 | 0.0794 (0.0808) | 50.0 ms | 1.034 |
| Butterworth 2, 40 Hz | 0.1188 | 0.0881 (0.0891) | 5.6 ms | 0.168 |
| Butterworth 4, 40 Hz | 0.0047 | 0.0847 (0.0859) | 10.4 ms | 0.232 |

The residual line and the delay are the same as for seed 7 (they do not depend on the noise); the noise and the RMS error change by a few per cent.

(b) Hamming FIR at 40 Hz: **51 taps −68.8 dB, 25 ms; 101 taps −73.1 dB, 50 ms; 201 taps −78.8 dB, 100 ms.** Doubling the length buys about 5 dB at 180 Hz (the stopband floor of the window, not the length, sets the level) and doubles the delay.

(c) For a 5 Hz crossover, rank by phase lag first: Butterworth 2 (−10.1°), Butterworth 4 (−18.7°), MA 50 (−44.1°, and 36 % loss at 10 Hz), FIR 101 (−90°). The phase at the crossover decides it; the residual line (0.12 deg/s for Butterworth 2) can be handled by a notch, while lost phase cannot be bought back.

**Common mistakes.** Ranking by attenuation alone; forgetting that RMS error against the true rate includes the delay.

## 5. Frequency warping

| f_c | Tustin \|H(f_c)\| | Tustin −3 dB at | (f_s/π) arctan(πf_c/f_s) | prewarped \|H(f_c)\| |
|---|---|---|---|---|
| 50 Hz | 0.7012 | 49.60 Hz | 49.59 Hz | 0.7071 |
| 150 Hz | 0.6500 | 140.18 Hz | 140.18 Hz | 0.7071 |
| 300 Hz | 0.4245 | 240.58 Hz | 240.58 Hz | 0.7071 |
| 450 Hz | 0.0501 | 304.03 Hz | 304.03 Hz | 0.7071 |

(b) Below about f_s/20 the warping is under 1 % and prewarping is a refinement; above about f_s/10 it moves the corner by several per cent and is essential; near f_s/2 plain Tustin cannot put a corner there at all (an analog 450 Hz corner lands at 304 Hz).

**Check.** The measured −3 dB frequencies agree with the warping law to the 0.005 Hz grid.

**Common mistakes.** Prewarping at the wrong frequency (it must be the one whose response matters, usually the corner or the notch centre).

## 6. Arithmetic and delay

(a) Fourth-order Butterworth, f_s = 1000 Hz:

| f_c | largest pole radius | radius after rounding a to 4 decimals | DC gain after rounding |
|---|---|---|---|
| 40 Hz | 0.908956 | 0.911441, stable | 0.977 |
| 5 Hz | 0.988051 | 1.050757, **unstable** | 0.009 |
| 1 Hz | 0.997598 | 1.000000, **unstable** | −7 × 10⁶ |

With a low corner the poles crowd near z = 1, and the denominator coefficients (for 1 Hz: 1, −3.98358, 5.95088, −3.95101, 0.98372) almost cancel: A(1) = Σa_k is a tiny difference of numbers of order 6, and rounding in the fourth decimal changes it completely. The engineer factors the filter into cascaded second-order sections (each with only two poles, so far less sensitive), uses enough coefficient bits, and checks the poles of the implemented coefficients, not the designed ones. Often the better fix is to decimate first so f_c/f_s is not so small.

(b) Phase at 5 Hz: **MA 50 −44.1°, FIR 101 −90.0°, Butterworth 2 −10.1°, Butterworth 4 −18.7°.** With 45° of margin, the FIR destroys the loop (−45° of margin), the moving average leaves about 1°, the fourth-order Butterworth leaves 26°, the second-order 35°. Only the IIR filters are acceptable, the second-order comfortably.

(c) She filters the stored record forwards and then backwards in time (zero-phase filtering): the phase of the backward pass cancels that of the forward pass, the magnitude is squared, and there is no delay. The flight controller cannot, because the backward pass needs samples from the future.

**Common mistakes.** Trusting the designed coefficients on a fixed-point target; reading phase lag without converting to margin.

## Check-yourself answers

1. H(e^{jω}) = Σb_k e^{−jωk} / Σa_k e^{−jωk}. FIR: no feedback (a = [a₀]); its output is a finite weighted sum of inputs, so a bounded input gives a bounded output, and it has no poles except at z = 0.
2. At k·500/10 = **50, 100 and 150 Hz**.
3. Pairing tap k with tap N − 1 − k gives e^{−jω(N−1)/2} times a real function, so the phase is −ω(N − 1)/2 (plus π jumps where the real factor changes sign): a delay of (N − 1)/2 samples at every frequency.
4. The ideal impulse response is infinitely long and non-causal; truncating it abruptly leaves sidelobes near −35 dB whatever the length. A Blackman taper lowers the worst sidelobe to about −76 dB (101 taps at 40 Hz) at the price of a wider transition band.
5. s = K(1 − z⁻¹)/(1 + z⁻¹) with K = 2f_s (or prewarped). The left half plane maps to the inside of the unit circle and the imaginary axis onto the circle.
6. At (1000/π) arctan(π·250/1000) = **211.9 Hz**. Prewarping sets K = 2π·250/tan(π/4) so that 250 Hz maps exactly; the corner is then at 250 Hz.
7. 0.3/√50 = **0.0424 deg/s** (the SYNTHETIC record gives 0.0427).
