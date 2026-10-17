# P09.10  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown.

## 1. Your servo gain (s = 15: K = 40 s⁻¹)

**Assumptions.** L = 40/(s(0.05s + 1)(0.01s + 1)); unity feedback.

(a) **GM = 120/40 = 3.00 (9.54 dB)** at **ω_pc = 1/√(0.0005) = 44.72 rad/s**.

(b) **PM = 25.39° at ω_gc = 24.54 rad/s; DM = 0.4431 rad/24.54 rad/s = 18.06 ms.**

(c) **N = 0, P = 0, Z = 0**: stable. **M_s = 2.72**: the curve comes within 0.37 of −1.

(d) **No**: 9.5 dB passes 6 dB, but 25° fails 45°.

**Check.** The Routh limit for 0.0005s³ + 0.06s² + s + K is K < 0.06/0.0005 = 120, so GM = 120/40.

**Common mistakes.** Quoting GM in dB as 20 log₁₀(1/3); forgetting that a stable loop can still fail its specification.

## 2. The criterion and the servo's margin

(a) F = 1 + L has zeros at the closed-loop poles and poles at the open-loop poles. Traverse the D-contour **clockwise**, enclosing the right half plane (indentations keep imaginary-axis poles outside). By the argument principle F(C) encircles the origin clockwise Z − P times, Z and P counted in the right half plane. F encircles 0 exactly when L = F − 1 encircles −1. So **N = Z − P, i.e. Z = N + P, with N the clockwise encirclements of −1 by L.**

(b) Near the origin L ≈ c/s^k with c = lim s^kL, real and positive for these loops. On s = εe^{jφ}, φ from −90° to +90°, L ≈ (c/ε^k)e^{−jkφ}: magnitude c/ε^k, huge, and angle running from +k·90° down to −k·90°, a **clockwise arc of kπ**. For k = 1 (the servo) it is a half circle through the positive real axis, joining the curve's ends at ω = 0⁻ and 0⁺. For k = 2 (the satellite) it is a full circle, starting and ending on the negative real axis side at −180°.

(c) L(jω) = K/D(jω) with D(jω) = jω(1 + jτ_mω)(1 + jτ_fω) = −(τ_m + τ_f)ω² + jω(1 − τ_mτ_fω²). L is real and negative where Im D = 0 with Re D < 0: **ω_pc² = 1/(τ_mτ_f)**. There \|L\| = K/((τ_m + τ_f)ω_pc²) = Kτ_mτ_f/(τ_m + τ_f), so **GM = (τ_m + τ_f)/(Kτ_mτ_f)**. Routh for τ_mτ_f s³ + (τ_m + τ_f)s² + s + K: stable while (τ_m + τ_f)·1 > τ_mτ_f K, i.e. **K < (τ_m + τ_f)/(τ_mτ_f)**. GM = K_crit/K: the same boundary.

(d) A delay multiplies L(jω) by e^{−jωT}: magnitude 1, phase −ωT. The curve reaches −1 when, at a frequency where \|L\| = 1, the extra lag equals the available phase margin: ω_gc T = PM, so **DM = PM/ω_gc**. With several gain crossovers each has its own PM_i/ω_gc,i and the loop fails at the smallest; a high-frequency crossover with a moderate PM can set the delay margin even if a low one has a smaller PM.

**Common mistakes.** Counting N counterclockwise; drawing the indentation's image on the wrong side (it closes through the right for k = 1).

## 3. An unstable airframe with a computer

(a) **PM = 103.49° at ω_gc = 17.51 rad/s; DM = 1.806 rad/17.51 rad/s = 0.1032 s.**

(b)

| delay | Nyquist (N, P, Z) | Padé(8) closed-loop RHP roots |
|---|---|---|
| 0.5 DM = 0.0516 s | (−1, 1, 0) | 0 |
| 0.9 DM = 0.0929 s | (−1, 1, 0) | 0 |
| 1.1 DM = 0.1135 s | (+1, 1, 2) | 2 |
| 1.5 DM = 0.1548 s | (+1, 1, 2) | 2 |

(c) N jumps from −1 to +1: the delay rotates the high-frequency part of the curve clockwise until it wraps around −1 the other way, so the counterclockwise encirclement that cancelled the open-loop unstable pole is replaced by a clockwise one, and Z = 1 + 1 = 2. With P = 1, stability needs N = −1, because each unstable open-loop pole must be "undone" by one counterclockwise turn of the curve about −1.

(d) For ω > 0 the curve never crosses the negative real axis, so there is no upper phase crossover: the gain can be raised without limit. But it starts on the negative real axis at L(0) = −3, and if the gain falls by a factor 3 (to K = 20/3) that point reaches −1 and the encirclement is lost. The **lower gain margin of 1/3 (−9.5 dB)** is the one that matters: a weak or failing actuator destabilises this airframe.

**Check.** The closed-loop constant term 3K − 20 vanishes at K = 20/3: a root at s = 0, matching L(0) = −1.

**Common mistakes.** Reporting GM = ∞ as "robust"; forgetting P when counting.

## 4. The servo sweep

(a)

| K [s⁻¹] | GM | PM | ω_gc [rad/s] | DM | M_s |
|---|---|---|---|---|---|
| 5 | 27.60 dB | 73.58° | 4.853 | 264.61 ms | 1.206 |
| 10 | 21.58 dB | 60.42° | 9.070 | 116.27 ms | 1.390 |
| 20 | 15.56 dB | 43.21° | 15.587 | 48.38 ms | 1.769 |
| 40 | 9.54 dB | 25.39° | 24.541 | 18.06 ms | 2.716 |
| 60 | 6.02 dB | 15.55° | 31.038 | 8.75 ms | 4.207 |
| 100 | 1.58 dB | 3.94° | 40.776 | 1.69 ms | 15.741 |

(b) DM falls fastest (a factor 157 from K = 5 to 100, against a factor 19 in PM), because DM = PM/ω_gc and ω_gc rises as PM falls: both move against it.

(c) GM ≥ 6 dB needs K ≤ 120/10^{0.3} = **60.1 s⁻¹**; PM ≥ 45° needs K ≤ **18.7 s⁻¹** (at 18.67, GM = 16.2 dB, ω_gc = 14.83 rad/s, DM = 52.9 ms). **The phase-margin requirement binds**: largest K ≈ 19 s⁻¹.

**Common mistakes.** Tuning to the gain margin alone; reading M_s from too coarse a grid near the resonance.

## 5. The satellite: where to put the crossover

(a)

| K_p [N m/rad] | PM | ω_gc [rad/s] | DM |
|---|---|---|---|
| 0.1 | 44.46° | 0.1265 | 6.135 s |
| 0.3162 | 54.90° | 0.3162 | 3.030 s |
| 1.0 | 44.46° | 0.7907 | 0.981 s |
| 3.0 | 28.48° | 1.5960 | 0.312 s |

(b) PM(ω) = 180° + ∠L = arctan(T_dω) − arctan(τω). Its derivative T_d/(1 + T_d²ω²) − τ/(1 + τ²ω²) vanishes at T_d(1 + τ²ω²) = τ(1 + T_d²ω²), i.e. (T_d − τ) = T_dτ(T_d − τ)ω², so **ω = 1/√(T_dτ) = 0.3162 rad/s**, where **PM = arctan(√10) − arctan(1/√10) = 54.90°**. K_p = 0.3162 places the crossover there.

(c) With a 0.5 s delay, PM falls by 0.5ω_gc rad: **40.84°, 45.84°, 21.81° and −17.24°** for the four gains (the last is unstable, Z = 2). **K_p = 0.1 and 0.3162** keep at least 30°.

**Check.** The two K_p values with PM 44.46° are reciprocal about the optimum: ω_gc = 0.1265 and 0.7907 rad/s multiply to 0.1 = 1/(T_dτ).

**Common mistakes.** Assuming more gain always means less phase margin (for this loop, too little gain also lowers it).

## 6. What the margins miss

(a) τ_f = 0.02 s: **K_crit = 0.07/0.001 = 70 s⁻¹**. At K = 20: **GM = 3.5 (10.88 dB), PM = 35.79° at 15.22 rad/s, DM = 41.0 ms.**

(b) With a 20 ms delay at K = 20: **GM = 1.87 (5.45 dB) at ω_pc = 23.59 rad/s, PM = 25.35°** at the unchanged ω_gc = 15.59 rad/s (DM left: 28.4 ms). The delay subtracts ωT of phase at every frequency, so the phase reaches −180° at a lower frequency, where \|L\| is larger; the gain margin falls although the gain did not change.

(c) M_s = 1.77 means the closest approach of the curve to −1, over all frequencies, is 1/1.77 = 0.57. GM and PM test only two points (on the negative real axis and on the unit circle); the curve can pass close to −1 between them, where a simultaneous small gain and phase change would destabilise the loop. M_s bounds both margins from below at once (`P09.18`).

(d) Rate saturation of the actuator (the surface cannot move faster than its maximum rate): in a large command the loop behaves as if it had an extra, amplitude-dependent delay, and a loop with good small-signal margins can enter a limit cycle. Backlash and friction in the linkage are others.

**Common mistakes.** Treating a sensor filter as "fast enough to ignore"; reading margins measured at one amplitude as valid at all amplitudes.

## Check-yourself answers

1. Z = N + P: Z closed-loop poles in the right half plane, P open-loop poles there, N clockwise encirclements of −1 by L(s) as s goes clockwise around the D-contour.
2. The contour may not pass through a pole of L. For a single integrator, the small semicircle maps to a large clockwise half circle (through the positive real axis) joining the ends of the curve at ω = 0⁻ and 0⁺.
3. GM = 1/\|L(jω_pc)\| where ∠L(jω_pc) = −180°; PM = 180° + ∠L(jω_gc) where \|L(jω_gc)\| = 1.
4. (π/4)/10 = **78.5 ms**.
5. At 1/√(0.05 × 0.01) = **44.72 rad/s**; GM = 120/20 = 6, **15.56 dB**.
6. **Once, counterclockwise** (N = −1).
7. A delay leaves \|L\| unchanged and adds phase lag ωT that grows without bound, so it eventually uses up the phase margin at the gain crossover even if the phase never crossed −180° before.
8. The curve can pass close to −1 between the phase and gain crossovers (a large peak sensitivity M_s), so a combined small change in gain and phase destabilises it; or the loop may have a second, nearby crossover that the single numbers hide.
