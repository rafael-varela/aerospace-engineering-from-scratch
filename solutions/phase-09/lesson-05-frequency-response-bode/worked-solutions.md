# P09.05  Worked solutions

Worked for **s = 15**. Every number below is produced by `solution.py` (its `main()` or the functions named) or by a probe script that calls them; the "independent check" line says how it was verified by a different route. Every sine test is a SIMULATION on the model, not a measurement.

## 1. A servo loop by hand (K = 15 s⁻¹, τ = 0.1 s)

**(a)** L = 15 · s⁻¹ · (1 + s/10)⁻¹: **K_B = 15, type 1, one real pole corner at 10 rad/s** (`bode_form` agrees).

**(b, c)**

| ω (rad/s) | straight-line \|L\| | exact \|L\| | straight-line ∠L | exact ∠L |
|---|---|---|---|---|
| 1 | 23.52 dB | 23.48 dB | −90.0° | −95.7° |
| 10 | 3.52 dB | 0.51 dB | −135.0° | −135.0° |
| 100 | −36.48 dB | −36.52 dB | −180.0° | −174.3° |

The magnitude differences are the real-pole corrections, −0.043 dB a decade away and −3.01 dB at the corner; the phase differences are the decade rule's 5.71° error a decade either side, and zero at the corner.
**(d)** ω_c²(1 + 0.01ω_c²) = 225 ⇒ **ω_c = 10.398 rad/s** (closed form and `gain_crossover` agree to 10⁻¹⁵), **phase −136.1°**. The straight line is on its −40 dB/decade segment at 0 dB, crossing at √(K/τ) = 12.25 rad/s: 18 % high, because the line ignores the 3 dB corner loss.

**Common mistakes.** Starting the low-frequency line at 20 log₁₀ K at ω = 0 instead of ω = 1; reading the crossover off the straight line near a corner.

## 2. Derivations

**(a)** With x = Re{Ae^{jωt}} for t ≥ 0, y(t) = Re{A∫₀ᵗ h(τ)e^{jω(t−τ)}dτ} = Re{Ae^{jωt}[∫₀^∞ h e^{−jωτ}dτ − ∫_t^∞ h e^{−jωτ}dτ]}. The first integral is G(jω); the second is bounded by ∫_t^∞|h|dτ, and for a stable system h is a sum of terms e^{p_k t} (times polynomials for repeated poles), so the bound decays like e^{Re p_slow · t}. What remains is Re{AG(jω)e^{jωt}} = A|G| cos(ωt + ∠G).
**(b)** G₁G₂ = |G₁||G₂|e^{j(∠G₁ + ∠G₂)}; 20 log₁₀ of a product is the sum of the logs.
**(c)** At the corner |1/(1 + j)| = 1/√2: 20 log₁₀(1/√2) = **−3.0103 dB**. A decade away, |1/(1 + 0.1j)| = 1/√1.01 → **−0.0432 dB** below (and the same above, relative to the −20 dB/decade line). The decade rule gives 0° at ω_b/10 where the true phase is −arctan 0.1 = −5.71°, and −90° at 10ω_b where it is −84.29°: **5.71°** error, the maximum (the true and straight curves are both antisymmetric about the corner, where they agree at −45°).
**(d)** |G|⁻² = (1 − u²)² + 4ζ²u², u = ω/ω_n. Its derivative in v = u² is −2(1 − v) + 4ζ² = 0 ⇒ v = 1 − 2ζ², so **ω_r = ω_n√(1 − 2ζ²)** (real only for ζ < 1/√2), and |G|⁻² = 4ζ⁴ + 4ζ²(1 − 2ζ²) = 4ζ²(1 − ζ²), so **M_r = 1/(2ζ√(1 − ζ²))**. At u = 1, |G| = 1/(2ζ) while both asymptotes meet at 0 dB: the error is **−20 log₁₀(2ζ)** dB.

## 3. Notch for a bending mode (f_m = 11 Hz)

**(a)** M_r = 1/(2·0.02·√(1 − 0.0004)) = 25.0, **+27.96 dB** at 10.996 Hz.

**(b)** Notch with ζ_n = 0.5:

| f | 1 Hz | 2 Hz | 10.5 Hz | 11.5 Hz | 22 Hz |
|---|---|---|---|---|---|
| \|N\| | −0.036 dB | −0.151 dB | −20.66 dB | −21.05 dB | −1.60 dB |
| ∠N | −5.24° | −10.65° | −84.7° | +84.9° | +33.7° |

At exactly 11 Hz the magnitude is zero (−∞ dB) and the phase jumps by 180°.

**(c)** Mode × notch peaks at **+1.24 dB at 7.72 Hz** over 1–40 Hz (against +27.96 dB for the mode alone): the resonance is removed.
**(d)** The notch costs 5.2° of phase at 1 Hz and 10.6° at 2 Hz, taken from the pitch loop's margin, so the loop must be re-checked with the notch in it. If the mode is 10 % lower (9.9 Hz) the notch attenuates it by only 13.7 dB there, and the combination peaks at **+14.4 dB at 9.86 Hz**: a narrow notch tuned to a modelled frequency is fragile. Remedies: a wider notch (more phase cost), a notch scheduled with fuel and stores, or a low-pass filter.

**Common mistakes.** Reading the notch's phase at its centre; forgetting that a notch also has phase lead above its centre.

## 4. How long must a sine test wait? (ω = 2.5 rad/s)

**(a)** Exact: **|G| = −1.6034 dB, ∠G = −26.964°**.

**(b)**

| settle | gain error | phase error |
|---|---|---|
| 0 s | −0.321 dB | −0.839° |
| 1 s | +0.0062 dB | +0.048° |
| 3 s | +0.0013 dB | +0.0099° |
| 10 s | < 10⁻⁴ dB | < 10⁻⁴ ° |

**(c)** The residual is the transient, dominated by the slowest closed-loop pole −1.068 s⁻¹, whose term decays as e^{−1.068 t}: from 0 to 1 to 3 s the errors fall by roughly the factors that exponential predicts once the fast poles have died.
**(d)** About 3 s of settling meets 0.01° at this frequency (the slow pole does not depend on the test frequency, so the same wait applies across the sweep), plus three periods per point: for 20 log-spaced frequencies from 0.1 to 10 rad/s the sweep lasts about **929 s (15.5 min)**, dominated by the long periods at the low end. This is why chirps and multisines are used for flight tests (extension 2).

## 5. Judgment

**(a)** Good enough for the shape of the loop: slopes, which corners matter, a first estimate of the crossover and the asymptotic gain. Not good enough for anything decided by a few degrees or a few dB: phase margin, gain margin, filter phase cost at a given frequency. 22.6° is the difference between a comfortable and a marginal loop.
**(b)** Identical magnitudes; the zero at −2 adds phase lead +arctan(ω/2), the zero at +2 adds lag −arctan(ω/2), a difference of 2 arctan(ω/2) (90° at ω = 2 rad/s). A controller designed from magnitude alone would assume the lead and find the RHP-zero plant far less stable than expected; this is why Bode's gain–phase relation applies only to minimum-phase systems.
**(c)** At ±0.2 deg a gear train's backlash and the motor's static friction (stiction) consume part of the motion: the output moves less than the linear model predicts at every frequency. A 3 dB loss is a factor 0.71. Repeat at ±0.5, ±1 and ±2 deg: if the deficit shrinks with amplitude it is a dead-zone or friction effect, not a wrong linear model; a deficit that stays constant points to a wrong gain (supply voltage, potentiometer scale).
**(d)** An open-loop transfer function with a right-half-plane pole has a transient that grows; no steady-state sine response exists, and the airframe would depart. The test is flown with a stabilising loop closed; the open-loop response is recovered from the measured closed-loop signals and the known controller (for example as the ratio of the plant output to the plant input spectra).

## 6. A time delay

**(a)** e^{−jωT} = cos ωT − j sin ωT: magnitude 1, phase −ωT (rad). |L_T| = |L|, so the 0 dB crossing does not move; only the phase there does.
**(b)** **ω_c = 12.213 rad/s, ∠L = −121.41°, 180° + ∠L = 58.59°.** At ω_c the delay adds −ω_cT:

| T | delay phase | 180° + phase |
|---|---|---|
| 10 ms | −7.00° | 51.59° |
| 20 ms | −14.00° | 44.59° |
| 50 ms | −34.99° | 23.60° |

**(c)** Padé phase −2 arctan(ωT/2) against −ωT: ωT = 0.5: −28.07° vs −28.65°; 1: −53.13° vs −57.30°; 2: −90.0° vs −114.6°; 3: −112.6° vs −171.9°. At ω_c the Padé error is 0.009°, 0.069° and 1.03° for T = 10, 20 and 50 ms: excellent while ωT ≲ 1, useless beyond about ωT = 2.
**(d)** 58.59° = 1.0226 rad; T = 1.0226/12.213 = **83.7 ms** brings 180° + phase to zero (`P09.10` calls this the delay margin).

**Common mistakes.** Letting the delay change the magnitude; mixing degrees and radians in ω_cT.

## Check-yourself answers (lesson.md)

1. y_ss = 1.5 cos(2t − 60°).
2. A cascade multiplies complex gains: magnitudes multiply (so their logarithms, the dB values, add) and angles add.
3. 1/√2 (−3.01 dB) and −45°; the straight line says 0 dB there, so it is 3.01 dB too high.
4. −20 log₁₀(2 × 0.05) = +20 dB above the straight line.
5. −20 dB/decade and −90°; K_B = lim_{s→0} sL(s), the value of the low-frequency asymptote at ω = 1 rad/s.
6. Its pole at s = 0 makes it marginally stable: a sine input produces an output that drifts (the integral of the input's mean and of the transient), so no steady state exists. It is measured inside a closed position loop and L is recovered from the closed-loop data.
7. 27/360 × 1 s = 75 ms.
8. e^{−1} ≈ 37 % of the slowest transient term remains.
