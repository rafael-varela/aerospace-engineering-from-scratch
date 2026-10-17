# P09.04  Worked solutions

Worked for **s = 15**. Every number below is produced by `solution.py` (its `main()` or the functions named) or by a probe script that calls them; the "independent check" line says how it was verified by a different route. Conventions: 10–90 % rise time, 2 % settling time, unit steps from rest.

## 1. Reading the poles (ω_n = 5 rad/s, ζ = 0.3)

**Poles** −ζω_n ± jω_n√(1 − ζ²) = **−1.5 ± 4.770j s⁻¹**.

| metric | formula | measured (`step_metrics`, 0.1 ms grid) |
|---|---|---|
| overshoot | 37.233 % | 37.233 % |
| peak time | π/ω_d = 0.6587 s | 0.6587 s |
| rise time | exact 0.26427 s | 0.26427 s |
| 2 % settling | envelope 2.639 s; 4/(ζω_n) = 2.667 s | **2.246 s** |

The settling time disagrees most (−15 % against 4/(ζω_n)): both formulas describe when the *envelope* enters the band, and the oscillating response itself enters (and stays) earlier, at a crossing that depends on the phase of the oscillation.

**Common mistakes.** Using ω_n instead of ω_d in the peak time; quoting 4/(ζω_n) as "the" settling time.

## 2. Derivations

**(a)** G/s = ω_n²/(s((s + σ)² + ω_d²)) = 1/s − (s + 2σ)/((s + σ)² + ω_d²) = 1/s − (s + σ)/((s + σ)² + ω_d²) − (σ/ω_d)·ω_d/((s + σ)² + ω_d²). Inverting: y = 1 − e^{−σt}[cos ω_d t + (σ/ω_d) sin ω_d t]. With σ/ω_d = ζ/√(1 − ζ²), cos ω_d t + (ζ/√(1 − ζ²)) sin ω_d t = sin(ω_d t + φ)/√(1 − ζ²) where sin φ = √(1 − ζ²), cos φ = ζ, i.e. **φ = arccos ζ**.
**(b)** ẏ = (ω_n/√(1 − ζ²))e^{−σt} sin ω_d t (the impulse response), zero first at **t_p = π/ω_d**. There sin(π + φ) = −sin φ = −√(1 − ζ²), so y(t_p) = 1 + e^{−σπ/ω_d} = 1 + exp(−πζ/√(1 − ζ²)). Inverting, ln M = −πζ/√(1 − ζ²) ⇒ ln²M (1 − ζ²) = π²ζ² ⇒ **ζ = −ln M/√(π² + ln²M)**.
**(c)** |y − 1| = e^{−σt}|sin(ω_d t + φ)|/√(1 − ζ²) ≤ e^{−σt}/√(1 − ζ²). The right side falls to 0.02 at t = −ln(0.02√(1 − ζ²))/σ; after that time the response is certainly inside the band, so the true settling time is no larger.
**(d)** T = K/(τs² + s + K) = (K/τ)/(s² + s/τ + K/τ): ω_n = √(K/τ), 2ζω_n = 1/τ, so ζ = 1/(2√(Kτ)) and **ζω_n = 1/(2τ)**. The poles are −1/(2τ) ± √(1/(4τ²) − K/τ); for K > 1/(4τ) the root is imaginary and the real part is −1/(2τ) whatever K.
**(e)** G_z = G − (s/z)G, so Y_z = Y − (1/z)sY, and sY is the transform of ẏ (y(0) = 0): **y_z = y − ẏ/z**. At t = 0⁺, ẏ = 0 and ÿ = ω_n² (from the ODE with y = ẏ = 0 and input 1), so ẏ_z(0⁺) = −ω_n²/z.

## 3. Thrust-vector gimbal (τ = 0.035 s)

**(a)** Overshoot ≤ 15 % ⇔ ζ ≥ ζ_min = **0.5169**: the poles must lie in the sector within arccos 0.5169 = 58.9° of the negative real axis. Peak time ≤ 60 ms ⇔ ω_d ≥ π/0.06 = **52.36 rad/s**: outside the band |Im s| < 52.36.
**(b)** K = 1/(4ζ_min²τ) = **26.73 s⁻¹**: ω_n = 27.64 rad/s, ω_d = 23.66 rad/s, real part −1/(2τ) = −14.29 s⁻¹, **t_p = 0.133 s > 0.06 s: not met**. Raising K raises ω_d but violates the overshoot; lowering K lowers ω_d. No proportional gain meets both.
**(c)** At ζ = ζ_min, ω_n = 1/(2ζ_min τ) and t_p = 2πζ_min τ/√(1 − ζ_min²). t_p ≤ 0.06 s needs **τ ≤ 15.8 ms**: the actuator itself must be more than twice as fast, or the loop needs velocity feedback or lead compensation.
**(d)** Simulated with RK4: overshoot 15.000 %, peak time 0.1328 s (2 % settling 0.284 s).

**Common mistakes.** Treating the two requirements separately and picking a K for each; forgetting that the real part is fixed.

## 4. A filtered potentiometer (τ_f = 15 ms)

**(a)** T = K(τ_f s + 1)/(s(τs + 1)(τ_f s + 1) + K): the filter pole becomes a **closed-loop zero at −66.7 s⁻¹**.
**(b)** With K = 14.31 s⁻¹: poles **−71.80 and −7.435 ± 14.508j**, zero −66.67, dominance ratio **9.66**. Overshoot **20.1 %** (designed for 10 %), peak time 0.215 s, settling 0.510 s.
**(c)** The ratio says the real pole at −71.8 is far away and almost cancelled by the zero at −66.7: that pair contributes little. But the filter has moved the *dominant pair itself*, from −10 ± 13.64j (ζ = 0.591) without the filter to −7.44 ± 14.51j (ζ = 0.456): the extra lag in the feedback path reduces the loop's damping. Dominance tells you whether you may ignore a pole once you know the pair, not whether the pair you calculated without that block is the right one.
**(d)** Bisection gives **K = 9.97 s⁻¹** for 10 % overshoot, with settling **0.430 s** and peak time 0.283 s. The price of the filter is a 30 % lower gain and a slower servo.

**Common mistakes.** Putting the filter in the forward path (it then adds a pole but no zero); computing the dominance ratio and stopping there.

## 5. Choosing the servo gain

**(a)** K = 8.2 s⁻¹:

| τ | 40 ms | 50 ms | 60 ms |
|---|---|---|---|
| overshoot | 0.36 % | 1.97 % | 4.10 % |
| 2 % settling | 0.308 s | 0.282 s | **0.509 s** |

**(b)** K = 14.31 s⁻¹: τ = 40 ms → overshoot 6.29 %, settling 0.318 s; τ = 60 ms → overshoot **13.35 %** (violates 10 %), settling 0.376 s.
**(c)** At 50 Hz (20 ms samples) against a 0.5 ms grid: rise time 0.1877 against 0.1871 s, overshoot 1.962 against 1.970 %, settling 0.2823 against 0.2820 s (interpolation recovers these well), but the peak time is 0.400 against 0.3925 s, quantised to the 20 ms grid.
**(d)** K = 8.2 s⁻¹ has the best nominal settling only because its first peak just touches the band; a 20 % slower motor nearly doubles the settling time. K = 10.5 s⁻¹ (the 5 % design) is a more honest choice: overshoot 2.2, 5.0 and 7.8 % and settling 0.336, 0.414 and 0.453 s for τ = 40, 50 and 60 ms, inside a 10 % overshoot specification over the whole range. To the colleague: "0.28 s holds only at τ = 50 ms exactly; quote the settling time over the τ range you expect, with the convention (2 %), and do not design to a saw-tooth minimum."

## 6. A right-half-plane zero (z = 3)

**(a)** From 2(e) the initial slope is −ω_n²/z, negative for z > 0, while y_z → 1 > 0: the response must first go negative and then turn round. Every RHP zero (real, z > 0) therefore produces undershoot.
**(b)** **Undershoot 4.48 %**, minimum at **t = 0.280 s**, **overshoot 17.0 %**, **2 % settling 8.42 s**; the numerical initial slope is −0.3333 = −ω_n²/z.
**(c)** With unity feedback T = n_G/(d_G + n_G): the numerator is n_G for every gain, so the zero stays in the closed loop. Feedback moves poles, not zeros; only a cancellation (which for an RHP zero means an unstable controller pole, an internally unstable loop, `P09.03` Problem 5) could remove it.
**(d)** Raising the elevator trailing edge produces a downward force on the tailplane before the aircraft has rotated; the total lift falls and the aircraft sinks briefly before the increased angle of attack makes it climb. An altitude loop that reacts strongly to the initial sink commands still more elevator, deepening the sink; too high a gain or bandwidth goes unstable. The RHP zero sets an upper limit on the achievable altitude-loop bandwidth (roughly below the zero's frequency), which is why altitude loops are slow outer loops around a pitch loop.

**Common mistakes.** Calling the initial dip "a transient of the closed loop that more gain will fix"; mixing up the sign convention of the zero location.

## Check-yourself answers (lesson.md)

1. Rise time: time from 10 % to 90 % of the final value. Peak time: time of the largest value. Overshoot: how far the peak exceeds the final value, in percent. Settling time: the time after which the response stays within ±2 % of the final value. The overshoot depends only on ζ.
2. ω_n = |−3 + 4j| = 5 rad/s, ζ = 3/5 = 0.6, overshoot 100 exp(−0.6π/0.8) = 9.48 %, peak time π/4 = 0.785 s.
3. ζ = −ln 0.05/√(π² + ln² 0.05) = 0.690.
4. T = K/(τs² + s + K): 2ζω_n = 1/τ, so Re s = −ζω_n = −1/(2τ) for every underdamped K. The envelope decay rate, hence roughly the settling time (4/(ζω_n) = 8τ), cannot be improved by gain.
5. Yes: ratio 20/2 = 10. With the real pole at −3 the ratio is 1.5 and the pair is not dominant; the real pole's slow mode shapes the response.
6. Its residue can be small (a nearby zero almost cancels it), or the opposite: a pole with a small dominance ratio can be hidden. Check the residues of G(s)/s, or simply compute the response.
7. It makes the response start in the wrong direction (undershoot) and limits how fast the loop can be. Altitude response to elevator: the aircraft sinks briefly before climbing.
8. It bounds the decay envelope (dropping the factor 1/√(1 − ζ²) and rounding ln 50 up), and the response enters the band before the envelope does. For the 10 % design it said 0.400 s against 0.350 s measured, 14 % high.
