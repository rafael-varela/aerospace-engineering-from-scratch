# P06.17  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and a probe script that calls it, given to the figures shown. Random time histories are SYNTHETIC. The input spectrum is NASA-STD-7001C Table 2; the masses, requirements and the Problem 5(b) spectrum are illustrative.

## 1. Reading an isolator (s = 15: m = 0.45 kg)

**Assumptions.** Rigid IMU, single axis, viscous ζ = 0.1, four identical mounts in parallel.

(a) k = 48 kN/m; **fₙ = 51.98 Hz**; c = 2ζ√(km) = **29.39 N·s/m**; sag **0.0919 mm** under 1 g and **0.460 mm** under 5 g. Check: (1/2π)√(9.80665/0.0919 × 10⁻³) = 51.99 Hz (the sag was rounded to three figures).

(b) T(1) = **5.099**, T(2) = **0.356**, and at 500 Hz (r = 9.62) T = **0.0237**: an isolation efficiency of **97.6 %**.

(c) Exact peak **5.123 at r = 0.9903** (51.48 Hz), against 1/(2ζ) = 5: the approximation is 2.4 % low, because the damper adds 4ζ² to the numerator.

**Common mistakes.** Using one mount's stiffness; quoting T at r = 1 as the peak.

## 2. Derivations

(a) With x = Xe^{iωt}, y = Ye^{iωt}: (k − mω² + icω)X = (k + icω)Y. Divide by k and use ωₙ² = k/m, c/k = 2ζ/ωₙ: X/Y = (1 + 2iζr)/(1 − r² + 2iζr), **T = √[(1 + (2ζr)²)/((1 − r²)² + (2ζr)²)]**, lag φ = atan2(2ζr, 1 − r²) − arctan(2ζr), which lies in [0, π). For a force F₀e^{iωt} on the mass over a fixed base, X = F₀/(k − mω² + icω) and the force through the mount is F_T = (k + icω)X, so **F_T/F₀ = (1 + 2iζr)/(1 − r² + 2iζr)**: the same ratio (isolating the source and isolating the receiver are one problem).

(b) T² = (1 + 4ζ²r²)/((1 − r²)² + 4ζ²r²). At r² = 2 the numerator is 1 + 8ζ² and the denominator 1 + 8ζ²: **T = 1** for every ζ. For r² > 2, (1 − r²)² > 1, so the denominator is larger than the numerator and T < 1.

(c) With u = r², write T² = N/D, N = 1 + 4ζ²u, D = (1 − u)² + 4ζ²u. dT²/du = 0 ⇔ N′D = ND′: 4ζ²[(1 − u)² + 4ζ²u] = (1 + 4ζ²u)[−2(1 − u) + 4ζ²]. Expanding, the cubic terms cancel and 2ζ²u² + u − 1 = 0 remains (after dividing by −2), so **u = r_p² = (√(1 + 8ζ²) − 1)/(4ζ²)**. For small ζ, √(1 + 8ζ²) ≈ 1 + 4ζ² − 8ζ⁴, so r_p² ≈ 1 − 2ζ² → 1, and T_peak ≈ √(1 + 4ζ²)/(2ζ) → **1/(2ζ)**.

(d) Replace k + icω by k(1 + iη): X/Y = (1 + iη)/(1 − r² + iη), **T_h = √[(1 + η²)/((1 − r²)² + η²)]**. For r ≫ 1, T_h ≈ √(1 + η²)/r²: −40 dB/decade. The viscous T ≈ √(1 + 4ζ²r²)/r² → 2ζ/r: −20 dB/decade, because the damper force icωY grows with frequency and is passed straight to the mass; the hysteretic force ikηY does not grow.

**Common mistakes.** Writing the force on the mass as −k x − c ẋ (that is a fixed base); losing the numerator 1 + 2iζr.

## 3. A camera on a helicopter

**Assumptions.** Single axis, ζ = 0.1 viscous, 1/rev = 5.6 Hz, 4/rev = 22.4 Hz.

(a) **fₙ ≤ 7.59 Hz** (r ≥ 2.95 at 22.4 Hz). k = 2.5 × (2π × 7.59)² = **5.68 kN/m**, sag **4.32 mm** under 1 g.

(b) T ≤ 1.1 needs r ≤ 0.302 at 1/rev, so **fₙ ≥ 18.5 Hz**. The two bounds (fₙ ≤ 7.6 Hz and fₙ ≥ 18.5 Hz) do not overlap: **no single isolator meets both**. The frequencies are only a factor of four apart, and an isolator needs fₙ below f/√2 at the frequency it rejects and well above the frequency it must pass.

(c) With fₙ = 7.59 Hz: T(1/rev) = **2.11** (the 1/rev is amplified by a factor of two, near resonance) and T(8/rev) = **0.046**.

(d) Hysteretic η = 0.2 gives **fₙ ≤ 8.02 Hz**, slightly better but the same conclusion. For a single dominant frequency, a **tuned vibration absorber** (a small mass on a spring tuned to 4/rev) or a notch device removes 4/rev without softening the mount at 1/rev; this is what helicopters use.

**Common mistakes.** Choosing fₙ between the two frequencies, which amplifies both.

## 4. Your own design (s = 15: sway limit 0.75 mm)

(a) fₙ = **64 Hz**, ζ = **0.22**, k = **56.6 kN/m**, c = **61.9 N·s/m**; transmitted **2.732 Grms**, sway **0.748 mm**, peak **2.53**, T(500 Hz) = **0.0595**. The **sway** (0.748 against 0.75 mm) and **T(500 Hz)** (0.0595 against 0.06) bind; fₙ ≥ 60 Hz and the peak are slack.

(b) SYNTHETIC simulation, seed 15, 32 s: transmitted **2.752 Grms (+0.80σ)**, sway RMS **0.1495 mm against 0.1484 mm (+0.81σ)**, with u = 0.94 %. Consistent.

(c) The tighter sway limit would push the lesson's 5 Hz-grid design (65 Hz, ζ = 0.2, 0.757 mm) out of bounds; on the 1 Hz, 0.01 grid the search finds a point with more damping and a slightly lower fₙ that meets both binding constraints, and its Grms (2.73) is 5 % lower than the lesson's coarse-grid 2.88. The optimum sits on the intersection of the two binding constraints, and a coarse grid cannot land on an intersection.

**Common mistakes.** Forgetting the 5 g static sag in the sway; reporting a grid optimum as if it were exact.

## 5. "Add more damping to be safe"

(a)

| ζ | transmitted Grms | peak T | T(500 Hz) |
|---|---|---|---|
| 0.1 | 3.97 | 5.12 | 0.032 |
| 0.2 | 2.88 | 2.73 | 0.056 |
| 0.4 | 2.38 | 1.66 | 0.107 |

ζ = 0.4 lowers the transmitted Grms and the peak but **breaks the T(500 Hz) ≤ 0.06 requirement** by almost a factor of two.

(b) The illustrative high-frequency spectrum has **12.46 Grms** overall. Transmitted at fₙ = 65 Hz: **1.012, 0.739, 0.632, 0.815 Grms** for ζ = 0.05, 0.1, 0.2, 0.4. The optimum is near ζ = 0.2; beyond it, damping makes things worse.

(c) Damping lowers the resonance, so it helps when much of the input's energy lies within a factor of about √2 of fₙ (the workmanship spectrum at 65 Hz). It hurts when the energy lies far above fₙ, where a viscous mount transmits T ≈ 2ζ/r, proportional to ζ. "More damping is safer" is true for the resonant peak and false for high-frequency isolation. Recommendation: keep ζ ≈ 0.2, or get the extra damping from a hysteretic elastomer, which lowers the peak without the −20 dB/decade penalty.

**Common mistakes.** Judging damping by the peak alone.

## 6. Viscous or hysteretic? (s = 15)

(a) Weighted least squares on the 25 SYNTHETIC points: viscous **fₙ = 68.5 ± 0.7 Hz, ζ = 0.038 ± 0.008, χ² = 1933 for 23 dof**; hysteretic **fₙ = 70.27 ± 0.20 Hz, η = 0.300 ± 0.008, χ² = 23.7 for 23 dof**.

(b) Slope of log T against log f above 300 Hz (8 points): **−2.04**. The viscous model predicts −1, the hysteretic −2.

(c) **The hysteretic model**: χ²/dof ≈ 1 and the right slope; the viscous model is rejected at enormous significance (χ²/dof ≈ 84), and its fitted ζ is distorted (it tries to match the high-frequency fall with very little damping). A designer who used ζ = η/2 = 0.15 in a viscous analysis predicts T(500 Hz) = 0.047 where the hysteretic mount gives 0.021: the viscous analysis is conservative at high frequency by a factor of 2.3, but its peak and sway predictions are close. Report which model was used, and fit the model to the high-frequency data before trusting it there.

**Common mistakes.** Fitting without weights (the high-frequency points then carry no weight); judging a fit by its parameter error bars rather than its χ².

## Check-yourself answers

1. r = √2. At r² = 2 the numerator 1 + 4ζ²r² and the denominator (1 − r²)² + 4ζ²r² are both 1 + 8ζ², whatever ζ is.
2. Approximately 1/(2ζ) = 5; exactly √(1 + 4ζ²)/(2ζ) = 5.10 at r = 1, and the true peak is 5.12 at r = 0.990.
3. The viscous damper transmits a force icωY that grows with frequency, so T → 2ζ/r; the hysteretic force ikηY does not grow, so T → √(1 + η²)/r².
4. fₙ = (1/2π)√(g/δ) = (1/2π)√(9.80665/10⁻⁴) = 49.8 Hz.
5. Lower fₙ (T ≈ 2ζ fₙ/f at fixed ζ, for f ≫ fₙ). The sway grows: the static sag as 1/fₙ² and the random sway roughly as fₙ^(−3/2).
6. The static sag under the 5 g launch acceleration (0.294 mm) and three standard deviations of the random relative displacement (3 × 0.154 = 0.463 mm); the random part dominates.
7. Because much of the workmanship spectrum's energy lies near the 65 Hz resonance, where damping reduces amplification; the energy far above fₙ is already attenuated so strongly that the extra transmission there adds little to the total.
8. Elastomer stiffness changes with frequency (dynamic stiffening), temperature and amplitude; and the mass is not what was assumed (harness, bracket, connectors). Also the mounts' rocking modes couple with translation.
