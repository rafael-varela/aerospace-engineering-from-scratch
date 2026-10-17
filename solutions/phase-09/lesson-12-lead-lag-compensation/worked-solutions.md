# P09.12  Worked solutions

Worked for **s = 15**. Every number below is produced by `solution.py` (the functions named in each problem) or by the probe scripts that call it; responses are simulations of representative models. Units: SI, angles in degrees, frequencies in rad/s.

## 1. A lead by hand (s = 15: α = 0.25)

**Derivation.** φ_max = arcsin(0.75/1.25) = arcsin 0.6 = **36.87°**. Gain at ω_m: 1/√0.25 = **2.000 (6.02 dB)**. T = 1/(40 × 0.5) = 0.05 s: zero **20 rad/s**, pole **80 rad/s**. High-frequency gain 1/α = **4.0**.

**Check.** `lead_max_phase(0.25)` returns (36.870, 2.0000); |`freq_response`| at 40 rad/s is 2.0000.

**(c)** Centred at 400 rad/s, the zero and pole move to 200 and 800 rad/s; φ_max, the gain at ω_m and the high-frequency gain do not change, because the lead depends on ω only through ωT.

**Common mistakes.** Putting ω_m at the arithmetic mean (50 rad/s) of the zero and pole instead of the geometric mean (40 rad/s).

## 2. Derivations

**(a)** dφ/dω = T/(1 + ω²T²) − αT/(1 + α²ω²T²) = 0 ⇒ 1 + α²ω²T² = α(1 + ω²T²) ⇒ ω²T²α(α − 1) = α − 1 ⇒ **ω_max = 1/(T√α)**. With x = 1/√α: tan φ = (x − αx)/(1 + αx²) = (1 − α)/(2√α), so sin φ = tan φ/√(1 + tan²φ) = (1 − α)/√(4α + (1 − α)²) = **(1 − α)/(1 + α)**. Gain: |1 + jx|/|1 + jαx| = √(1 + 1/α)/√(1 + α) = **1/√α**.

**(b)** At s = 0, (s + z)/(s + z/β) = β; as s → ∞ it tends to 1. At ω = 10z the phase is arctan 10 − arctan 10β. It is negative for β > 1 and decreases toward arctan 10 − 90° = **−5.71°** as β → ∞.

**(c)** s L(s) = Kc[(Ts + 1)/(αTs + 1)] · [(s + z)/(s + z/β)] · K/(τs + 1) → **Kc · β · K** as s → 0.

**(d)** |L(jω)| = ω_n²/(ω√(ω² + 4ζ²ω_n²)) = 1 ⇒ (ω²)² + 4ζ²ω_n²ω² − ω_n⁴ = 0 ⇒ ω_c² = ω_n²(√(1 + 4ζ⁴) − 2ζ²). The phase is −90° − arctan(ω_c/(2ζω_n)), so PM = arctan(2ζω_n/ω_c) = **arctan[2ζ/√(√(1 + 4ζ⁴) − 2ζ²)]**.

| ζ | PM (deg) | 100ζ | overshoot (%) |
|---|---|---|---|
| 0.3 | 33.27 | 30 | 37.23 |
| 0.5 | 51.83 | 50 | 16.30 |
| 0.55 | 55.68 | 55 | 12.63 |
| 0.7 | 65.16 | 70 | 4.60 |

The rule PM ≈ 100ζ is good to about 4° up to ζ = 0.6 and drifts above that.

**Common mistakes.** Forgetting that the overshoot formula belongs to the *closed* loop while PM belongs to the *open* loop; differentiating arctan(αωT) without the inner factor αT.

## 3. Heavier surface (s = 15: τ = 0.065 s)

**(a)** Plant phase at 30 rad/s: **−152.850°**, so the lead must add 27.850°. `design_lead` gives **Kc = 1.9810, α = 0.3632, zero 18.079 rad/s, pole 49.781 rad/s, φ_max = 27.850°**. Measured on the loop: **ω_c = 30.000 rad/s, PM = 55.000°, GM = ∞** (no phase crossover), Kv = **39.620 s⁻¹**.

**(b)** `design_lead_lag(..., 200)`: **β = 5.5757**, lag zero 3.000 rad/s, lag pole **0.5381 rad/s**, lead **α = 0.3006, Kc = 1.7935**. Measured: ω_c = 30.000 rad/s, PM = 55.000°, GM = ∞, **Kv = 200.00 s⁻¹**. Delay margin 32.0 ms, as for any 55° at 30 rad/s.

**(c)**

| design | overshoot (%) | rise (ms) | settle (ms) | ramp error (mrad per rad/s) |
|---|---|---|---|---|
| lead | 13.93 | 42.5 | 153.5 | 25.240 |
| lead-lag | 16.97 | 41.0 | 435.0 | 5.000 |

**(d)** Compared with τ = 0.05 s (lead α = 0.4669, lead-lag α = 0.3896, β = 5.9534), α moved most: the heavier surface costs 6.5° more phase at 30 rad/s, which only the lead can supply. β fell slightly, because the lead's larger Kc (1.98 instead of 1.85) already supplies part of Kv.

**Independent check.** 1/Kv of the lead alone is 1/39.620 = 25.24 mrad per rad/s, equal to the simulated ramp error.

## 4. Gimbal

**(a)** With C = 1: **ω_c = 7.793 rad/s, PM = 43.21°, phase crossover 22.361 rad/s, GM = 6.000 (15.56 dB)**. By hand: the phase is −180° when arctan(0.1ω) + arctan(0.02ω) = 90°, i.e. (0.1ω)(0.02ω) = 1, ω = √500 = 22.36 rad/s. There |G| = 10/(22.36 × √6 × √1.2) = 10/60, so GM = 6.

**(b)** At 15 rad/s the plant phase is −163.01°, so the lead must add 33.01°: **Kc = 1.5326, α = 0.2947**. Measured: **ω_c = 15.000 rad/s, PM = 50.000°, phase crossover 38.967 rad/s, GM = 4.586 (13.23 dB)**. Closed-loop step: **17.6 % overshoot, settling 418 ms**. The gain margin fell from 6.0 to 4.6. The lead raises |L| at all frequencies above its zero, including near the phase crossover (which it also pushes up to 39 rad/s). The lesson's servo has no phase crossover at all (its phase never reaches −180°), so its GM is infinite and stays infinite.

**(c)** `design_lead(num, den, 40, 50)` raises **"needs 74.6 deg of lead at wc, more than 65.0 deg from one stage"**. The plant phase at 40 rad/s is −204.62°.

**Common mistakes.** Reporting the PM at the old crossover; forgetting that a third-order plant has a finite GM that a lead can reduce.

## 5. Judgment

**(a)** Lead alone, ω_c = 30 rad/s:

| PM (deg) | α | Kc | Kv (s⁻¹) | Kc/α |
|---|---|---|---|---|
| 40 | 0.8020 | 2.4216 | 48.43 | 3.020 |
| 50 | 0.5615 | 2.0263 | 40.53 | 3.609 |
| 55 | 0.4669 | 1.8478 | 36.96 | 3.957 |
| 60 | 0.3858 | 1.6796 | 33.59 | 4.354 |
| 70 | 0.2562 | 1.3686 | 27.37 | 5.343 |

The PM = 70° design passes the most high-frequency noise: **5.343/3.020 = 1.77 times** the PM = 40° design (2 mV rms of feedback noise becomes 10.7 mV rms of command instead of 6.0). More phase needs a smaller α, and the high-frequency gain grows as 1/α faster than Kc falls.

**(b)** The lag puts a pole at z/β = 0.504 rad/s and a zero at 3 rad/s. The closed loop of the lesson's design has poles at −23.66 ± 33.11j, −17.99 and **−3.251**, and zeros at −18.72 and **−3.000**. The slow pole sits next to the lag zero, so its residue is small but its time constant is 1/3.251 = 0.31 s. The 10–90 % rise takes 41.5 ms, and the last few percent then creep in on that 0.3 s time scale. The 2 % settling time is set by that tail (471 ms against 147.5 ms for the lead alone). Accept it when the ramp-following or load-rejection requirement (Kv) matters more than the last 2 % of a step. For a surface servo inside a pitch loop it usually does, because the outer loop keeps commanding ramps.

**(c)** A PI is a lag with its pole moved to the origin. It adds an integrator and raises the system type by one: zero error to a ramp for this type-1 plant, and zero error to a constant load torque (an aerodynamic hinge moment). Use it when that error must be zero, not just small. The cost is more phase lag at low frequency and a slow closed-loop pole that the integrator creates.

**(d)** At 1 kHz the hold delay is dt/2 = 0.5 ms. Its phase cost at 30 rad/s is 30 × 0.0005 rad = **0.86°**, against a 55° margin, and 0.5 ms against the **32 ms** delay margin. The continuous design can be used, discretised with a method from `P09.17`, provided the computation delay is also small.

## 6. Two stages (gimbal, ω_c = 40 rad/s, PM = 50°)

**(a)** Total lead **74.62°** (plant phase −204.62°), so each stage supplies **37.31°**: **α = 0.2452**, zero **19.81 rad/s**, pole **80.78 rad/s**, **Kc = 5.1790**.

**(b)** Measured: **ω_c = 40.000 rad/s, PM = 50.000°, phase crossover 96.46 rad/s, GM = 3.852 (11.71 dB)**. Closed-loop step: **20.8 % overshoot, rise 29.0 ms, settling 205 ms**.

**(c)** High-frequency gain Kc/α²: **86.13** for the two stages at 40 rad/s. At 30 rad/s (62.53° needed), one stage gives Kc/α = 2.7043/0.05975 = **45.26**, and two stages give 3.5034/0.3167² = **34.94**. Splitting a large lead between two stages lowers the noise gain for the same phase (here by 23 %), but pushing the crossover from 30 to 40 rad/s costs a factor 2.5. Two stages are not a free lunch: the price of more bandwidth is noise and a smaller gain margin.

**(d)** The uncompensated loop crosses −180° at 22.36 rad/s with GM = 6.0. The leads add phase there and move the crossing up to 96.46 rad/s, where the two stages' combined gain is large. |L| at the new phase crossover is 1/3.852 = 0.26 instead of 0.17, so the GM falls from 15.6 dB to 11.7 dB.

**Common mistakes.** Centring the two stages at different frequencies without re-solving for Kc; quoting the GM at the uncompensated phase crossover.

## Check-yourself answers (lesson.md)

1. C_lead = Kc(Ts + 1)/(αTs + 1), 0 < α < 1. Its phase is largest at ω_max = 1/(T√α), the geometric mean of zero and pole, where sin φ_max = (1 − α)/(1 + α).
2. Kv = Kc K = 200 needs Kc = 10. The crossover then rises to 61.7 rad/s, where the plant's phase is −162°, leaving 18.0° of margin.
3. The lag multiplies Kv by its DC gain β. Placed a decade below the crossover, it has gain ≈ 1 and costs at most 5.7° of phase there, so it hardly disturbs the margin.
4. T = 1/(40 × 0.5) = 0.05 s: zero 20 rad/s, pole 80 rad/s, φ_max = 36.87°, gain 2 (6.02 dB).
5. The lag changes the phase (−4.75° for the lesson's design) and the gain at ω_c slightly. Designing the lead on G alone would miss the PM and the crossover by that amount.
6. PM ≈ 100ζ (degrees) for ζ up to about 0.6. For PM = 55°, ζ ≈ 0.55 and the overshoot is e^{−πζ/√(1−ζ²)} ≈ 12.6 %; the lesson's lead loop gives 13.3 %.
7. 0.960 rad/30 rad/s = 32.0 ms. In a digital servo it is consumed by the sample-and-hold (dt/2), the computation delay, and sensor and filter lag.
8. The formulas give PM = 55° only if they and the code that implements them are right. Measuring |L| = 1 and the phase on the actual compensated loop, and confirming with `scipy.signal.bode`, checks everything at once by an independent route.
