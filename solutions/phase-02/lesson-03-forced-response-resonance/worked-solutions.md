# P02.03  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown.

## 1. Sensor mount ($s = 15$: f = 28 Hz)

**Assumptions.** Linear mount, steady state, SI units.

**Derivation.** (a) ω_n = √(40 000/1.2) = 182.57 rad/s = 29.06 Hz; √(km) = 219.1 kg/s, ζ = 15/(2 × 219.1) = 0.0342; ω = 175.93 rad/s, r = 0.9636. (b) (1 − r²)² = (0.0714)² = 0.00510; (2ζr)² = (0.0660)² = 0.00435; M = 1/√(0.00945) = **10.28**; X = 10.28 × 5/40 000 = **1.285 mm**; φ = atan2(0.0660, 0.0714) = **42.7°**. (c) D = (k − mω²)² + (cω)² = (40 000 − 37 140)² + (2639)² = 8.18 × 10⁶ + 6.96 × 10⁶ = 1.514 × 10⁷; A = 5 × 2860/D = 0.944 mm, B = 5 × 2639/D = 0.872 mm; √(A² + B²) = 1.285 mm ✓, atan2(B, A) = 42.7° ✓. (d) 1/(2ζ) = **14.6**; f_peak = 29.06√(1 − 2 × 0.00117) = **29.02 Hz**; Δf ≈ 2ζf_n = 1.99 Hz (code: 1.992 Hz); Q = **14.6**.

**Check.** The code gives M = 10.281, X = 1285.2 μm, φ = 42.71°, (A, B) = (944.3, 871.7) μm, Q = 14.57 ✓. Limiting case: at r = 0 the same force would give 125 μm; the mount amplifies ten-fold because the drive sits 4 % below resonance with Q = 15.

**Common mistakes.** Using f instead of ω in cω; forgetting to square the (2ζr) term; reading φ from atan instead of atan2 (wrong quadrant above resonance).

## 2. Two routes to the steady state

**Derivation.** (a) x_p = A cos ωt + B sin ωt, ẋ_p = ω(−A sin ωt + B cos ωt), ẍ_p = −ω²x_p. m ẍ_p + c ẋ_p + k x_p = [(k − mω²)A + cωB] cos ωt + [(k − mω²)B − cωA] sin ωt = F₀ cos ωt. Hence (k − mω²)A + cωB = F₀ and −cωA + (k − mω²)B = 0. Cramer: D = (k − mω²)² + (cω)², A = F₀(k − mω²)/D, B = F₀cω/D.

(b) z = Z e^{iωt}: (−mω² + icω + k)Z = F₀, Z = F₀/(k − mω² + icω) = F₀(k − mω² − icω)/D. Re[Z e^{iωt}] = Re Z cos ωt − Im Z sin ωt = [F₀(k − mω²)/D] cos ωt + [F₀cω/D] sin ωt. So A = F₀ Re H, B = −F₀ Im H with H = 1/(k − mω² + icω) ✓.

(c) M⁻² = D(u) = (1 − u)² + 4ζ²u, u = r². dD/du = −2(1 − u) + 4ζ² = 0 ⟹ u = 1 − 2ζ². D_min = 4ζ⁴ + 4ζ²(1 − 2ζ²) = 4ζ²(1 − ζ²), so M_peak = 1/(2ζ√(1 − ζ²)). A peak requires u > 0: ζ < 1/√2. (D″ = 2 > 0 confirms a minimum of D.)

(d) x_p = t(A cos ω_n t + B sin ω_n t): ẋ_p = (A cos + B sin) + tω_n(−A sin + B cos); ẍ_p = 2ω_n(−A sin + B cos) − tω_n²(A cos + B sin). m ẍ_p + k x_p = 2mω_n(−A sin + B cos) + t(k − mω_n²)(A cos + B sin) = 2mω_n(B cos ω_n t − A sin ω_n t) since k = mω_n². Equal to F₀ cos ω_n t ⟹ A = 0, B = F₀/(2mω_n).

(e) x = [2F₀/(m(ω_n² − ω²))] sin((ω_n − ω)t/2) sin((ω_n + ω)t/2). As ω → ω_n: sin((ω_n − ω)t/2) → (ω_n − ω)t/2, ω_n² − ω² = (ω_n − ω)(ω_n + ω) → 2ω_n(ω_n − ω), (ω_n + ω)/2 → ω_n. x → [2F₀/(2mω_n(ω_n − ω))] × [(ω_n − ω)t/2] sin ω_n t = (F₀/(2mω_n)) t sin ω_n t ✓.

**Check.** (b) at ω = 0: H = 1/k, A = F₀/k, B = 0: the static deflection ✓.

**Common mistakes.** In (a), dropping the sin-equation and "solving" for A alone; in (c), differentiating with respect to r instead of r² and making a sign slip; in (d), forgetting the product rule on t.

## 3. Designing an accelerometer ($s = 15$: f_max = 250 Hz)

(a) An accelerometer's output is the relative displacement of its proof mass, (a/ω_n²)M(r, ζ); it reads the base acceleration faithfully only where M ≈ 1, the stiffness-controlled regime r ≪ 1, so f_n must sit well above every frequency to be measured. An isolator wants the opposite: transmissibility ≪ 1 is the mass-controlled regime r ≫ √2, so its f_n must sit well below the lowest disturbance frequency.

(b) ζ = 0.7: |M − 1| first exceeds 0.05 at r = **0.591** (M falls to 0.95 there; it never rises above 1.05 since M_peak = 1/(2 × 0.7 × 0.714) = 1.0002). Hence f_n ≥ 250/0.591 = **423 Hz**.

(c)

| ζ | r_max for ±5 % | f_n,min (Hz) |
|---|---|---|
| 0.50 | 0.322 | 776 |
| 0.60 | 0.844 | 296 |
| 0.65 | 0.720 | 347 |
| 0.70 | 0.591 | 423 |
| 0.75 | 0.476 | 525 |

M⁻² = 1 + (4ζ² − 2)r² + r⁴. For ζ = 1/√2 the r² term vanishes (maximally flat) and M = (1 + r⁴)^{−1/2} falls below 0.95 at r⁴ = 0.108, r = 0.573; for ζ a little below, the negative r² coefficient makes M rise above 1 first, and if the rise stays under +5 % (M_peak = 1/(2ζ√(1 − ζ²)) < 1.05 holds for ζ > 0.56) it *delays* the fall to 0.95, widening the band; at ζ = 0.5 the peak is 1.155 and the +5 % limit is hit early. The optimum for a ±5 % tolerance is near ζ = 0.6, allowing f_n = 296 Hz instead of 423 Hz.

(d) The sensitivity is 1/ω_n²: f_n = 2500 Hz would give 100 times less output than f_n = 250 Hz for the same acceleration, and the signal would drown in electronic noise. Accelerometers are designed with f_n as low as the bandwidth allows, which is exactly why the ζ trade in (c) matters.

**Check.** Dimensions: a/ω_n² is (m/s²)/(s⁻²) = m ✓.

**Common mistakes.** Taking the ±5 % band to end at the peak; assuming ζ = 1/√2 is optimal for every tolerance (it is optimal for the *flatness at r = 0*, not for a finite band).

## 4. The tolerance band ($s = 15$: X_lim = 0.6 mm)

**Result.** (a) ζ = 0.02: the steady amplitude exceeds 0.6 mm between **52.03 and 57.77 Hz** (width 5.74 Hz; peak 1.675 mm). (b) ζ = 0.05: **53.47–56.22 Hz** (width 2.75 Hz; peak 0.671 mm, barely over). ζ = 0.10: never exceeded (peak 0.337 mm). The limit stops being exceeded when M_peak F₀/k = X_lim, i.e. 1/(2ζ√(1 − ζ²)) = 8.96, **ζ = 0.0559** (the light-damping estimate F₀/(2kX_lim) = 0.0558). (c) Half-power widths are 2.20 Hz (ζ = 0.02) and 5.51 Hz (ζ = 0.05). At ζ = 0.02 the tolerance band (5.7 Hz) is 2.6 times the half-power band because 0.6 mm is far below the peak (M = 8.96 against 25); at ζ = 0.05 it is half the half-power band because 0.6 mm is just under the peak. The two agree only if X_lim happens to be the −3 dB level, M_peak/√2. (d) From rest at 55 Hz the amplitude passes 0.6 mm at **t = 0.067 s** (`full_response`), against 0.064 s from X_ss(1 − e^{−t/τ}) = X_lim with τ = 0.145 s; the 3 ms difference is the half-cycle discretisation of "amplitude" on a 55 Hz carrier.

**Check.** Both bands are symmetric about 55 Hz to within 0.2 Hz, as they must be for light damping (the peak is at 54.98 Hz).

**Common mistakes.** Bracketing `brentq` across the peak (two roots in one interval); forgetting that the band collapses to nothing, not to zero width, once the peak is below the limit.

## 5. What the sweep and the force model hide

(a) X = m_e e ω²/√((k − mω²)² + (cω)²). Divide by mω²: X = (m_e e/m)/√((r⁻² − 1)² + (2ζ/r)²) → **m_e e/m** as ω → ∞. The high-frequency floor is set by the ratio of unbalance to total mass and is independent of k: stiffening moves the resonance out of the band but cannot reduce the vibration below m_e e/m, and if the drive sits above the new f_n the amplitude is *near* that floor regardless. Stiffening helps only when it moves the drive into the stiffness-controlled regime, where X = m_e e ω²/k grows as ω² (so a drive band that extends upward is bad news for a stiff design).

(b) 20 Hz in 2 s is 10 Hz/s, between the lesson's 6.7 Hz/s (91 % of the peak) and 20 Hz/s (78 %). The measured peak is lower and broader than the steady one, so the reported Q = 18 is biased *low* (true 25) and the apparent f_n is biased *high* by about 1.5 Hz. Require a rate below ζ²f_n² ≈ 1.2 Hz/s; at 0.33 Hz/s the lesson's table shows 99.9 %.

(c) Coulomb friction in the bolted joints: the dissipated energy per cycle grows linearly with amplitude while the stored energy grows quadratically, so the equivalent ζ *falls* as amplitude rises; a two-level shaker test shows Q increasing with force level. Material (viscoelastic) damping in the aluminium and any elastomer washers: energy per cycle roughly ∝ amplitude², so ζ is amplitude-independent; Q unchanged between levels. Micro-slip in joints at high amplitude adds a third mechanism with the opposite trend (ζ rising with amplitude).

(d) The phase of the response relative to the force across the band. One mode gives a single smooth 0 → 180° swing through 90° at f_n; a second mode inside 50–60 Hz adds a second swing (or a dip and recovery), visible even when the two amplitude peaks merge into one broad hump.

**Common mistakes.** Treating F₀ as constant for rotating machinery; believing a fast sweep over-reads the peak.

## 6. A chirp through resonance ($s = 15$: T = 3, 15, 75, 300 s)

**Derivation.** (a) Phase Φ(t) = 2π(45t + 10t²/T) rad; f_inst = Φ′/(2π) = 45 + 20t/T Hz, from 45 Hz at t = 0 to 65 Hz at t = T; rate 20/T Hz/s.

(b) RK4, 10 000 steps/s, from rest; steady resonant amplitude F₀/(cω_n) = 1.675 mm:

| T (s) | rate (Hz/s) | peak (mm) | fraction | f at peak (Hz) |
|---|---|---|---|---|
| 3 | 6.67 | 1.523 | 0.909 | 56.24 |
| 15 | 1.33 | 1.654 | 0.988 | 55.33 |
| 75 | 0.267 | 1.673 | 0.999 | 55.05 |
| 300 | 0.067 | 1.675 | 1.000 | 54.99 |

(c) τ = 1/(ζω_n) = 0.145 s; time inside the half-power band Δt = 2ζf_n/rate = 2.2 Hz/rate: 0.33, 1.65, 8.25, 33 s. 1 − e^{−Δt/τ} = 0.898, 1.000, 1.000, 1.000. The argument gets the fast case right (0.898 against 0.909) and says all slower sweeps are complete, which they are to 1 %. It slightly *under*-predicts at T = 3 s because the amplitude keeps growing after the drive has left the band's centre; a sharper estimate uses the time to traverse half the band.

(d) The amplitude lags the drive: it is still building when f_inst passes f_n, so the maximum occurs after f_n has been passed, at a *higher* frequency on an upward sweep (56.24 Hz at 6.67 Hz/s). A downward sweep should show the mirror image. Verified: F = F₀ cos(2π(65t − 10t²/T)) with T = 3 s peaks at 1.537 mm at **53.69 Hz**, 1.3 Hz *below* f_n (the upward sweep peaked 1.24 Hz above).

**Check.** As T → ∞ the peak and its frequency converge to the steady values 1.675 mm and 54.98 Hz ✓ (T = 300 s: 1.675 mm at 54.99 Hz).

**Common mistakes.** Using the instantaneous frequency 45 + 20t/T as the *argument* of the cosine (the phase is its integral, so the quadratic term has the factor ½); too few steps per cycle at 65 Hz.

## Check-yourself answers

1. Linearity: L[x_h + x_p] = 0 + F. The initial conditions determine the constants in x_h (the transient); x_p is fixed by the forcing alone.
2. H = 1/(k − mω² + icω). At ω_n: |H| = 1/(cω_n), phase lag 90°.
3. M = 1/√((1 − r²)² + (2ζr)²). For ζ = 0.02: r = 0.5 → 1.333; r = 1 → 25.0; r = 2 → 0.333.
4. At r = √(1 − 2ζ²), slightly below ω_n; no peak for ζ ≥ 1/√2.
5. Q = 55/2.2 = 25; ζ = 1/(2Q) = 0.02.
6. The guess is itself a homogeneous solution (c = 0, ω = ω_n), so substituting it gives zero, not F₀ cos ω_n t. Multiply the guess by t: x_p = (F₀/(2mω_n)) t sin ω_n t.
7. Beat period 2π/|ω_n − ω| = 1/|55 − 52| = 1/3 s; the envelope sin((ω_n − ω)t/2) has frequency |f_n − f|/2 = 1.5 Hz (two envelope lobes per beat period).
8. At r = √2. Above it, more damping transmits *more* motion, because the damper is a stiff path for high-frequency force.
